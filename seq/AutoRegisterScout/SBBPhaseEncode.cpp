// * -------------------------------------------------------------------------- *
// * Includes                                                                   *
// * -------------------------------------------------------------------------- *
#include "MrServers/MrImaging/seq/Kernels/SBBPhaseEncode.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/libRT.h"
#include "MrServers/MrImaging/libSBB/libSBBmsg.h"          // SBB error messages
#include "MrServers/MrProtSrv/MrProt/KSpace/MrKSpace.h"
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrSysSpec.h"
#include "MrServers/MrMeasSrv/SeqFW/libGSL/libGSL.h"


SBBPhaseEncode::SBBPhaseEncode( SBBList *pSBBList )
    : SeqBuildBlock                     (pSBBList)
    , m_dEpsilon                        (0.0000001)
    , m_lStartTime                      (-1)
    , m_lTimeToEcho                     (-1)
    , m_dDeltaMoment                    (0.0)
    , m_dOffsetMoment                   (0.0)
    , m_lFirstLinParToMeasure           (0)
    , m_bFirstLinParToMeasureDefined    (false)
    , m_lLastLinParToMeasure            (0)
    , m_bLastLinParToMeasureDefined     (false)
    , m_dFactorForPixelSize             (1.0)
    , m_dFactorForOffsetMoment          (1.0)
    , m_bFlowComp                       (false)
    , m_eGradientAxis                   (None)
    , m_eTableDirction                  (SEQ::DIR_ASCENDING)
    , m_eDimension                      (SEQ::DIM_2)
{
    m_sTB.setIdent ("phase encoding gradient");
    m_sTB_FC.setIdent ("phase encoding gradient flow compensated");
}



bool SBBPhaseEncode::prep (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo*)
{

    static const char *ptModule = "SBBPhaseEncode::prep";

    // * ---------------------------------------------------------------------- *
    // * Check settings                                                         *
    // * ---------------------------------------------------------------------- *
    if ( m_eGradientAxis == None )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s No valid gradient axis has been defined.", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }

    if ( m_bFlowComp && (m_lTimeToEcho == -1) )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Time to echo has not been defined.", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }

    if ( m_bFirstLinParToMeasureDefined == false )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s First line or partition has not been defined.", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }

    if ( m_bLastLinParToMeasureDefined == false )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Last line or partition has not been defined.", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Read the dimension setting from the protocol                           *
    // * ---------------------------------------------------------------------- *
    m_eDimension = pMrProt->kSpace().dimension();



    // * ---------------------------------------------------------------------- *
    // * Update of the gradient performance settings (FAST, NORMAL, WHISPER)    *
    // * ---------------------------------------------------------------------- *
    updateGradientPerformance (pMrProt->gradSpec().mode());



    // * ---------------------------------------------------------------------- *
    // * Incremental moment to get from one point in k-space to an adjacent     *
    // * point                                                                  *
    // * ---------------------------------------------------------------------- *
    switch ( m_eGradientAxis )  {
        case Phase:     m_dDeltaMoment = fGSLGetPEDeltaMoment( pMrProt );   break;
        case Slice:     m_dDeltaMoment = fGSLGet3DDeltaMoment( pMrProt );   break;
        default:
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "No valid gradient axis defined." ) )  return ( false );
    }

    if ( m_eTableDirction == SEQ::DIR_DESCENDING )  {
        m_dDeltaMoment *= -1.0;
    }


    // * ---------------------------------------------------------------------- *
    // * Maximum moments in k-space                                             *
    // * ---------------------------------------------------------------------- *
    double dMaxMomentFirstLine = m_dDeltaMoment * m_dFactorForPixelSize * m_lFirstLinParToMeasure + m_dOffsetMoment * m_dFactorForOffsetMoment;
    double dMaxMomentLastLine  = m_dDeltaMoment * m_dFactorForPixelSize * m_lLastLinParToMeasure  + m_dOffsetMoment * m_dFactorForOffsetMoment;



    if ( m_bFlowComp )  {

        // * ------------------------------------------------------------------ *
        // * Ramp up and ramp down time of the phase encoding gradients m_sTB   *
        // * and m_sTB_FC                                                       *
        // * ------------------------------------------------------------------ *
        m_lRampDuration = static_cast<long>(ceil (maximum (m_sTB.getMaxMagnitude()    * m_sTB.getMinRiseTime(),
                                                           m_sTB_FC.getMaxMagnitude() * m_sTB_FC.getMinRiseTime())));

        m_lRampDuration = RoundUpGRT (m_lRampDuration);



        // * ------------------------------------------------------------------ *
        // * Determine the maximum amplitude and check for a division by zero   *
        // * ------------------------------------------------------------------ *
        double dMaxAmplitude = minimum (m_sTB.getMaxMagnitude(), m_sTB_FC.getMaxMagnitude());
        if ( fabs (dMaxAmplitude) < m_dEpsilon  )  {
            if ( ! pSeqLim->isContextPrepForBinarySearch() ) {
                setNLSStatus (SBB_ERROR, ptModule, "Invalid maximum gradient amplitude encountered");
            } else {
                setNLSStatus (SBB_ERROR);
            }
            return ( false );
        }




        // * -------------------------------------------------------------- *
        // * Duration of m_sTB minus duration of m_sTB_FC must be negative, *
        // * i.e. the first gradient is always shorter than the second one. *
        // * This ensures that the sqrt can be calculated in any situation. *
        // * -------------------------------------------------------------- *
        double dDeltaDuration = -1.0 * maximum (fabs(dMaxMomentFirstLine), fabs(dMaxMomentLastLine)) / dMaxAmplitude;



        // * --------------------------------------------------------------- *
        // * Calculate the duration of m_sTB for the maximum moment          *
        // * dMaxMomentFirstLine. The gradient raster time is not considered *
        // * at this point.                                                  *
        // * --------------------------------------------------------------- *
        double dDuration1 = m_lRampDuration * m_lRampDuration + 2.0 * dDeltaDuration * dDeltaDuration;
        dDuration1 -= 4.0 * m_lTimeToEcho * dDeltaDuration    + 2.0 * m_lRampDuration * dDeltaDuration;

        if ( dDuration1 < 0.0 )  {
            if ( ! pSeqLim->isContextPrepForBinarySearch() ) {
                setNLSStatus (SBB_ERROR, ptModule, "Invalid arguennt for sqrt encountered");
            } else {
                setNLSStatus (SBB_ERROR);
            }
            return ( false );
        }

        dDuration1 =  0.5 * sqrt (dDuration1) - m_lRampDuration / 2.0;



        // * ------------------------------------------------------------------ *
        // * Calculate the duration of m_sTB_FC for the maximum moment          *
        // * dMaxMomentFirstLine. The gradient raster time is not considered at *
        // * this point.                                                        *
        // * ------------------------------------------------------------------ *
        double dTemp1     = (2.0 * dDuration1 - 2.0 * m_lTimeToEcho - m_lRampDuration) / 2.0;
        double dTemp2 = m_lRampDuration * m_lRampDuration + 4.0 * m_lTimeToEcho * m_lRampDuration;
        dTemp2 += 8.0 * m_lRampDuration * dDuration1      + 4.0 * m_lTimeToEcho * m_lTimeToEcho;
        dTemp2 += 8.0 * dDuration1 * dDuration1;
        dTemp2 =  0.5 * sqrt (dTemp2);

        double dDuration2 = 0.0;
        if ( dTemp1 > dTemp2 )  {
            dDuration2 = dTemp1 - dTemp2;
        } else {
            dDuration2 = dTemp1 + dTemp2;
        }


        // * ------------------------------------------------------------------ *
        // * Round gradient durations to the gradient raster time               *
        // * ------------------------------------------------------------------ *
        dDuration1 = RoundUpGRT (dDuration1);
        dDuration2 = RoundUpGRT (dDuration2);


        // * ------------------------------------------------------------------ *
        // * Check duration values:                                             *
        // *  - Duration must be longer than ramp up/down times                 *
        // *  - Exception: Both durations are 0 ==> set ramp time to 0, too     *
        // *    I.e. both gradients are not used and inherently "flow           *
        // *    compensated".                                                   *
        // * ------------------------------------------------------------------ *
        if ( dDuration1 < m_lRampDuration || dDuration2 < m_lRampDuration )  {

            if ( (fabs(dDuration1) < m_dEpsilon) && (fabs(dDuration2) < m_dEpsilon) )  {

                // * ---------------------------------------------------------- *
                // * Both gradient pulses have duration 0 and amplitude 0.0.    *
                // * I.e. both gradients are not used and automatically flow    *
                // * compensated.                                               *
                // *    ==> ramp time must be 0, too                            *
                // * ---------------------------------------------------------- *
                m_lRampDuration = 0;

            } else {

                // * ---------------------------------------------------------- *
                // * Not handled by the algorithm, yet                          *
                // * ---------------------------------------------------------- *
                if ( pSeqLim->isContextNormal() )  {
                    setNLSStatus(SEQU_ERROR, ptModule, "Pulse duration shorter than ramp time.\n");
                } else {
                    setNLSStatus(SEQU_ERROR);
                }
                return ( false );

            }

        }


        // * ------------------------------------------------------------------ *
        // * Prepare the phase encoding gradients for the maximum amplitude     *
        // * Amplituds valid for the actual phase encoding step have to be set  *
        // * prior to execution of the SBB.                                     *
        // * ------------------------------------------------------------------ *
        if ( ! m_sTB.prepAmplitude (m_lRampDuration, dDuration1, m_lRampDuration, 0.0) )  {
            if ( pSeqLim->isContextNormal() )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Preparation of phase encoding gradient failed.\n");
            } else {
                setNLSStatus(SEQU_ERROR);
            }
            return ( false );
        }


        if ( ! m_sTB_FC.prepAmplitude (m_lRampDuration, dDuration2, m_lRampDuration, 0.0) )  {
            if ( pSeqLim->isContextNormal() )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Preparation of phase encoding gradient failed.\n");
            } else {
                setNLSStatus(SEQU_ERROR);
            }
            return ( false );
        }

    } else {

        long lDuration = 0;

        // * ------------------------------------------------------------------ *
        // * Calculate a non-flow compensated encoding table                    *
        // * ------------------------------------------------------------------ *
        long lStatus = fGSLGetShortestTabTiming (
                         m_dDeltaMoment * m_dFactorForPixelSize,       // IMP: Moment between two steps [mT/m * us]
                         SEQ::DIR_ASCENDING,                           // IMP: true <=> (-) -> (+), i.e. table direction
                         m_dOffsetMoment * m_dFactorForOffsetMoment,   // IMP: the offset MOMENT of the table [us(mT/m)]
                         m_lFirstLinParToMeasure,                      // IMP: number of first line/partition (note: Echo=0)
                         m_lLastLinParToMeasure,                       // IMP: number of last line/partition (note: Echo=0)
                         m_sTB.getMaxMagnitude(),                      // IMP: maximum allowed gradient amplitude [mT/m]
                         m_sTB.getMinRiseTime(),                       // IMP: minimum allowed rise time [us/(mT/m)]
                        &m_lRampDuration,                              // EXP: minimum required ramp time [us] (up/down)
                        &lDuration                                     // EXP: minimum required duration  [us]
                      );

        if ( pSeqLim->isContextNormal() )  {
            if ( setNLSStatus (lStatus, ptModule, "Calculate of non-flow compensated encoding table failed") ) return ( false );
        } else {
            if ( setNLSStatus (lStatus) ) return ( false );
        }




        // * ---------------------------------------------------------------------- *
        // * Prepare the phase encoding gradients for the maximum amplitude         *
        // * Amplituds valid for the actual phase encoding step have to be set      *
        // * prior to execution of the SBB.                                         *
        // * ---------------------------------------------------------------------- *
        if ( ! m_sTB.prepAmplitude (m_lRampDuration, lDuration, m_lRampDuration, m_sTB.getMaxMagnitude()) )  {
            if ( pSeqLim->isContextNormal() )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Preparation of phase encoding gradient failed.\n");
            } else {
                setNLSStatus(SEQU_ERROR);
            }
            return ( false );
        }



        // * ------------------------------------------------------------------ *
        // * Flow compensation gradient m_sTB_FC is not required.               *
        // * ------------------------------------------------------------------ *
        if ( ! m_sTB_FC.prepAmplitude (0, 0, 0, 0.0) )  {
            if ( pSeqLim->isContextNormal() )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Preparation of phase encoding gradient failed.\n");
            } else {
                setNLSStatus(SEQU_ERROR);
            }
            return ( false );
        }


    }



    // * ---------------------------------------------------------------------- *
    // * Calculate the duration of the SBB                                      *
    // * ---------------------------------------------------------------------- *
    m_lSBBDurationPerRequest_us = m_sTB.getTotalTime() + m_sTB_FC.getTotalTime();



    setPrepared();


    return ( true );

}





bool SBBPhaseEncode::check (SeqLim* pSeqLim)
{

    static const char *ptModule = "SBBPhaseEncode::check";

    // * ---------------------------------------------------------------------- *
    // * Check first gradient table                                             *
    // * ---------------------------------------------------------------------- *
    if (! m_sTB.check() )  {

        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB.check failed.", ptModule);
        }

        setNLSStatus (m_sTB.getNLSStatus());
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Check second gradient table / used for flow compensation               *
    // * ---------------------------------------------------------------------- *
    if ( m_bFlowComp )  {

        if (! m_sTB_FC.check() )  {

            if ( pSeqLim->isContextNormal() )  {
                TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB_FC.check failed.", ptModule);
            }

            setNLSStatus ( m_sTB_FC.getNLSStatus() );
            return ( false );
        }

    }



    return ( true );

}




bool SBBPhaseEncode::disable ()
{

    static const char *ptModule = "SBBPhaseEncode::disable";

    if ( ! m_sTB.prepAmplitude (0, 0, 0, 0.0) )  {
        setNLSStatus (m_sTB.getNLSStatus(), ptModule, "Disabling gradient m_sTB failed");
        return (false);
    }

    if ( ! m_sTB_FC.prepAmplitude (0, 0, 0, 0.0) )  {
        setNLSStatus (m_sTB_FC.getNLSStatus(), ptModule, "Disabling gradient m_sTB_FC failed");
        return (false);
    }


    return (true);

}




bool SBBPhaseEncode::resize (long lDuration)
{

    // * ------------------------------------------------------------------ *
    // * Assure proper duration / gradient raster time                      *
    // * ------------------------------------------------------------------ *
    long lNewDuration = RoundUpGRT (lDuration);


    if ( !m_bFlowComp && (lNewDuration > m_sTB.getDuration()) )  {

        m_sTB.setDuration ( lNewDuration );

        // * ------------------------------------------------------------------ *
        // * Re-calculate the duration of the SBB                               *
        // * ------------------------------------------------------------------ *
        m_lSBBDurationPerRequest_us = m_sTB.getTotalTime() + m_sTB_FC.getTotalTime();

    }

    return ( true );

}






bool SBBPhaseEncode::updateAmplitude (SeqLim* pSeqLim, long lLinPar)
{

    static const char *ptModule = "SBBPhaseEncode::updateAmplitude";



    // * ---------------------------------------------------------------------- *
    // * Check dimension and encoding direction for concistency                 *
    // * ---------------------------------------------------------------------- *
    if ( (m_eDimension == SEQ::DIM_2) && (m_eGradientAxis == Slice) && (lLinPar != 0))  {
        setNLSStatus (SBB_ERROR, ptModule, "Partition index must be 0 for 2D imaging.");
        return ( false );        
    }

    

    // * ---------------------------------------------------------------------- *
    // * Check the actual line / partition index                                *
    // * ---------------------------------------------------------------------- *
    if ( (lLinPar < m_lFirstLinParToMeasure) || (lLinPar > m_lLastLinParToMeasure) )  {
        setNLSStatus (SBB_ERROR, ptModule, "Line / partition index out of range.");
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Actual gradient moment                                                 *
    // * ---------------------------------------------------------------------- *
    double dMoment = 0.0;

    if ( SysProperties::isPhaseEncodingEnabled() )  {
        dMoment = lLinPar * m_dDeltaMoment + m_dOffsetMoment;
    } else {
        dMoment = m_dOffsetMoment;
    }





    // * ---------------------------------------------------------------------- *
    // * Check whether the moment is larger than 0, calculation has to be       *
    // * performedonly only in this situation.                                  *
    // * ---------------------------------------------------------------------- *
    if ( fabs(dMoment) > m_dEpsilon )  {

        if ( m_bFlowComp )  {

            double dTemp;

            // * ------------------------------------------------------------------ *
            // * Amplitude of the first gradient (m_sTB)                            *
            // * ------------------------------------------------------------------ *
            double dAmplitude = dMoment * (6.0 * m_lTimeToEcho + 3.0 * m_lRampDuration + 3.0 * m_sTB_FC.getDuration());
            dTemp = -6.0 * m_lRampDuration * m_sTB.getDuration() - 3.0 * m_sTB.getDuration() * m_sTB_FC.getDuration() - 3.0 * m_sTB.getDuration() * m_sTB.getDuration();

            if ( fabs(dTemp) > m_dEpsilon )  {
                dAmplitude /= dTemp;
            } else {
                setNLSStatus (SBB_ERROR, ptModule, "Invalid divisor.");
                return ( false );
            }


            m_sTB.setAmplitude (dAmplitude);

            if (! m_sTB.check() )  {

                if ( pSeqLim->isContextNormal() )  {
                    TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB.check failed.", ptModule);
                }

                setNLSStatus (m_sTB.getNLSStatus());
                return ( false );
            }


            // * ------------------------------------------------------------------ *
            // * Amplitude of the second gradient (m_sTB_FC)                        *
            // * ------------------------------------------------------------------ *
            dAmplitude =  dMoment * (9.0 * m_lRampDuration + 6.0 * m_sTB_FC.getDuration() + 3.0 * m_sTB.getDuration() + 6.0 * m_lTimeToEcho);
            dTemp = 6.0 * m_lRampDuration * m_sTB_FC.getDuration() + 3.0 * m_sTB.getDuration() * m_sTB_FC.getDuration() + 3.0 * m_sTB_FC.getDuration() * m_sTB_FC.getDuration();

            if ( fabs(dTemp) > m_dEpsilon )  {
                dAmplitude /= dTemp;
            } else {
                setNLSStatus (SBB_ERROR, ptModule, "Invalid divisor.");
                return ( false );
            }


            m_sTB_FC.setAmplitude (dAmplitude);

            if (! m_sTB_FC.check() )  {

                if ( pSeqLim->isContextNormal() )  {
                    TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB_FC.check failed.", ptModule);
                }

                setNLSStatus (m_sTB_FC.getNLSStatus());
                return ( false );
            }

        } else {

            // * ------------------------------------------------------------------ *
            // * Prepare non-flow compensated gradient table for the current line   *
            // * partiton.                                                          *
            // * ------------------------------------------------------------------ *
            m_sTB.prepMomentumTOT (dMoment);


            if (! m_sTB.check() )  {

                if ( pSeqLim->isContextNormal() )  {
                    TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB.check failed.", ptModule);
                }

                setNLSStatus (m_sTB.getNLSStatus());
                return ( false );
            }

        }

    } else {

        // * ---------------------------------------------------------------------- *
        // * The moment is zero ==> the amplitudes of both gradients have to be     *
        // * zero, too.                                                             *
        // * ---------------------------------------------------------------------- *

        m_sTB.setAmplitude (0.0);

        if (! m_sTB.check() )  {

            if ( pSeqLim->isContextNormal() )  {
                TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB.check failed.", ptModule);
            }

            setNLSStatus (m_sTB.getNLSStatus());
            return ( false );
        }

        m_sTB_FC.setAmplitude (0.0);

        if (! m_sTB_FC.check() )  {

            if ( pSeqLim->isContextNormal() )  {
                TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_sTB_FC.check failed.", ptModule);
            }

            setNLSStatus (m_sTB_FC.getNLSStatus());
            return ( false );
        }

    }





    return ( true );

}






bool SBBPhaseEncode::run (MrProt*, SeqLim* pSeqLim, SeqExpo*, sSLICE_POS*)
{

    // * ---------------------------------------------------------------------- *
    // * Check successful preparation of SBBPhaseEncode                         *
    // * ---------------------------------------------------------------------- *
    if ( ! isPrepared() )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s SBBReadOut not prepared, cannot execute", ptModule);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Check the start time of the SBB                                        *
    // * ---------------------------------------------------------------------- *
    if ( m_lStartTime == -1 )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s No valid start time has been defined.", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }




    // * ---------------------------------------------------------------------- *
    // * Execute timing                                                         *
    // * ---------------------------------------------------------------------- *
    switch ( m_eGradientAxis )  {
        case Slice:
            fRTEI(m_lStartTime,         0,      0,      0,      0,       0,    &m_sTB,      0);

            if ( m_bFlowComp )  {
                fRTEI(m_lStartTime +
                      m_sTB.getTotalTime(), 0,  0,      0,      0,       0, &m_sTB_FC,      0);
            }
        break;

        case Phase:
            fRTEI(m_lStartTime,         0,      0,      0,    &m_sTB,       0,      0,      0);

            if ( m_bFlowComp )  {
                fRTEI(m_lStartTime +
                      m_sTB.getTotalTime(), 0,  0,      0, &m_sTB_FC,       0,      0,      0);
            }
        break;

        default:
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "No valid gradient axis defined." ) )  return ( false );
    }


    return ( true );

}





void SBBPhaseEncode::updateGradientPerformance (SEQ::Gradients eGradientMode)
{

    m_sTB.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBPhaseEncode::GradTable) );
    m_sTB.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBPhaseEncode::GradTable) );

    m_sTB_FC.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBPhaseEncode::GradTable) );
    m_sTB_FC.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBPhaseEncode::GradTable) );

}
