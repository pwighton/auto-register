// * -------------------------------------------------------------------------- *
// * Includes                                                                   *
// * -------------------------------------------------------------------------- *
#include "MrServers/MrImaging/seq/Kernels/SBBReadOut.h"
#include "MrServers/MrImaging/libSBB/SBBList.h"
#include "MrServers/MrImaging/libSeqUtil/KernelCalculationLimits.h"
#include "MrServers/MrImaging/libSBB/libSBBmsg.h"          // SBB error messages
#include "MrServers/MrMeasSrv/SeqIF/libRT/libRT.h"
#include "MrServers/MrImaging/ut/libSeqUT.h"
#include "MrServers/MrProtSrv/MrProt/KSpace/MrKSpace.h"
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrRXSpec.h"
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrSysSpec.h"
#include "MrServers/MrMeasSrv/SeqFW/libGSL/libGSL.h"




SeqBuildBlockReadOut::SeqBuildBlockReadOut( SBBList* pSBBList )
    : SeqBuildBlock                 (pSBBList)
    , m_dEpsilon                    (0.0000001)
    , m_lStartTime                  (-1)
    , m_lContrast                   (-1)
    , m_lReadOutColumns             (0)
    , m_lReadOutColumnsBefore       (0)
    , m_dReadOutAsymBefore          (0.5)
    , m_dExactReadOutAsymBefore     (0.5)
    , m_dReadOutAsymAfter           (0.5)
    , m_dExactReadOutAsymAfter      (0.5)
    , m_dEffectiveDW_us             (0.0)
    , m_dROMomentIn                 (0.0)
    , m_dROMomentOut                (0.0)
    , m_dROMomentToKSpaceCenter     (0.0)
    , m_dROMomentAfterKSpaceCenter  (0.0)
    , m_dROMomentFlatTop            (0.0)
    , m_bUseROP                     (false)
    , m_bRampUpInSBB                (true)
    , m_bRampDownInSBB              (true)
    , m_lRampUpTime                 (0)
    , m_lRequestedRampUpTime        (-1)
    , m_lRampDownTime               (0)
    , m_lRequestedRampDownTime      (-1)
    , m_bSymmetricTiming            (false)
    , m_lADCShift                   (0)
    , m_lTimeToEcho                 (0)
    , m_bEchoOnGRT                  (false)
    , m_eROPolarity                 (Positive)
    , m_dFactorForPixelSizeRO       (1.0)
    , m_dFactorForUpperPixelSizeRO  (1.0)
    , m_lStartTimeADC               (0)
    , m_lMinTimeBetweenADCs         (0)
    , m_lLeadTimeBetweenADCs        (0)
    , m_lTimeBetweenADCsPreviousRO  (0)
    , m_bFreqPhaseUpdated           (false)
    , m_lStartTimeRO                (0)
    , m_lStartTimeROP               (0)
    , m_bUseRampSampling            (false)
{
    m_ADC.setIdent ("ADC");
    m_RO.setIdent (RTEIDENT_PulseRO);
}




bool SeqBuildBlockReadOut::updateFreqPhase (long lLineIndex, long lPartitionIndex, sSLICE_POS* pSLC, double dAddPhase)
{

    // * ------------------------------------------------------------------ *
    // * Frequency and phase due to offset shift                            *
    // * ------------------------------------------------------------------ *
    m_ADCSet.prepSet( *pSLC, m_ADC, m_RO, lLineIndex, lPartitionIndex );    /*! EGA-05 !*/
    m_ADCNeg.prepNeg( *pSLC, m_ADC, m_RO, lLineIndex, lPartitionIndex );    /*! EGA-05 !*/



    // * ------------------------------------------------------------------ *
    // * Additional phase, e.g. rf spoiling                                 *
    // * ------------------------------------------------------------------ *
    m_ADCSet.increasePhase ( dAddPhase );    /*! EGA-05 !*/
    m_ADCNeg.decreasePhase ( dAddPhase );    /*! EGA-05 !*/


    // * ------------------------------------------------------------------ *
    // * Frequency / phase has been updated successfully                    *
    // * ------------------------------------------------------------------ *
    m_bFreqPhaseUpdated = true;

    return ( true );
}






bool SeqBuildBlockReadOut::init (MrProt*, SeqLim*, SeqExpo*)
{

    m_lStartTime                 = -1;
    m_lContrast                  = -1;
    m_lReadOutColumns            = 0;
    m_lReadOutColumnsBefore      = 0;
    m_dReadOutAsymBefore         = 0.5;
    m_dExactReadOutAsymBefore    = 0.5;
    m_dReadOutAsymAfter          = 0.5;
    m_dExactReadOutAsymAfter     = 0.5;
    m_dEffectiveDW_us            = 0.0;
    m_dROMomentIn                = 0.0;
    m_dROMomentOut               = 0.0;
    m_dROMomentToKSpaceCenter    = 0.0;
    m_dROMomentAfterKSpaceCenter = 0.0;
    m_dROMomentFlatTop           = 0.0;
    m_bUseROP                    = false;
    m_bRampUpInSBB               = true;
    m_bRampDownInSBB             = true;
    m_lRampUpTime                = 0;
    m_lRequestedRampUpTime       = -1;
    m_lRampDownTime              = 0;
    m_lRequestedRampDownTime     = -1;
    m_bSymmetricTiming           = false;
    m_lADCShift                  = 0;
    m_lTimeToEcho                = 0;
    m_bEchoOnGRT                 = false;
    m_eROPolarity                = Positive;
    m_dFactorForPixelSizeRO      = 1.0;
    m_dFactorForUpperPixelSizeRO = 1.0;
    m_lStartTimeADC              = 0;
    m_lLeadTimeBetweenADCs       = 0;
    m_lMinTimeBetweenADCs        = 0;
    m_lLeadTimeBetweenADCs       = 0;
    m_lTimeBetweenADCsPreviousRO = 0;
    m_bFreqPhaseUpdated          = false;
    m_lStartTimeRO               = 0;
    m_lStartTimeROP              = 0;


    return ( true );

}





bool SeqBuildBlockReadOut::prep (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo*)
{

    static const char *ptModule = "SeqBuildBlockReadOut::prep";

    // * ---------------------------------------------------------------------- *
    // * Update of the gradient performance settings (FAST, NORMAL, WHISPER)    *
    // * ---------------------------------------------------------------------- *
    if (! updateGradPerf ( pMrProt->gradSpec().mode() ) ) {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s Gradient performance values have not successfully been updated", ptModule);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }




    // * ---------------------------------------------------------------------- *
    // * Get scaling factors for gradient timing calculation                    *
    // * ---------------------------------------------------------------------- *
    if ( m_pCalcLimits )  {
        m_dFactorForPixelSizeRO      = m_pCalcLimits->getFactorForPixelSizeRO      ( pMrProt );
        m_dFactorForUpperPixelSizeRO = m_pCalcLimits->getFactorForUpperPixelSizeRO ( pMrProt );
    } else {
        m_dFactorForPixelSizeRO      = 1.0;
        m_dFactorForUpperPixelSizeRO = 1.0;
    }




    // * ---------------------------------------------------------------------- *
    // * Check contrast ID                                                      *
    // * ---------------------------------------------------------------------- *
    if ( m_lContrast == -1 )  {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s Invalid contrast, must be larger than or equal zero: %ld", ptModule, m_lContrast);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Both flags shift the ADC on the sub-GRT raster and cannot be used      *
    // * simultaneously                                                         *
    // * ---------------------------------------------------------------------- *
    if ( m_bSymmetricTiming && m_bEchoOnGRT )  {
        TRACE_PUT1 (TC_INFO, TF_SEQ, "%s Invalid combination of m_bSymmetricTiming and m_bEchoOnGRT", ptModule);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Calculatet effective dwell time (without oversampling) in us           *
    // * ---------------------------------------------------------------------- *
    m_dEffectiveDW_us = pMrProt->rxSpec().effDwellTime(pSeqLim->getReadoutOSFactor())[m_lContrast] / 1000.0;



    // * ---------------------------------------------------------------------- *
    // * Determin the minimum duration between consecutive ADCs                 *
    // * ---------------------------------------------------------------------- *
    m_lMinTimeBetweenADCs = static_cast<long>(ceil(SysProperties::getMinDurationBetweenReadoutAndReadout(pMrProt->rxSpec().realDwellTime()[m_lContrast] / 1000.0)));



    // * ---------------------------------------------------------------------- *
    // * Calculation of readout asymmetries:                                    *
    // *                                                                        *
    // *                 dAsymmetryBefore       dAsymmetryAfter                 *
    // *              <--------------------><-------------------->              *
    // *              0                   127                     255           *
    // *              ____________________________________________              *
    // *              |                    |                      |             *
    // *              |                    |                      |             *
    // *       _______|                    |                      |______       *
    // * ---------------------------------------------------------------------- *
    m_lReadOutColumnsBefore = static_cast<long>(pMrProt->kSpace().baseResolution() * m_dReadOutAsymBefore) / 2L * 2L;

    if ( m_lReadOutColumnsBefore )  {

        m_dExactReadOutAsymBefore = static_cast<double>(m_lReadOutColumnsBefore) / pMrProt->kSpace().baseResolution();

    } else {

        if (m_bUseRampSampling) 
        {
            m_dExactReadOutAsymBefore = static_cast<double>(m_lReadOutColumnsBefore) / pMrProt->kSpace().baseResolution();
        }
        else
        {
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s Invalid number of readout points prior to echo: %ld", ptModule, m_lReadOutColumnsBefore);
            setNLSStatus ( SBB_ERROR );
            return ( false );        
        }
    }


    long lReadOutColumnsAfter  = static_cast<long>(pMrProt->kSpace().baseResolution() * m_dReadOutAsymAfter) / 2L * 2L;

    if ( lReadOutColumnsAfter )  {

        m_dExactReadOutAsymAfter = static_cast<double>(lReadOutColumnsAfter) / pMrProt->kSpace().baseResolution();

    } else {

        TRACE_PUT2(TC_INFO, TF_SEQ, "%s Invalid number of readout points after the echo: %ld", ptModule, lReadOutColumnsAfter);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }

    m_lReadOutColumns = m_lReadOutColumnsBefore + lReadOutColumnsAfter;



    // * ---------------------------------------------------------------------- *
    // * Preparation of the ADC:                                                *
    // * ---------------------------------------------------------------------- *
    m_ADC.Mdh.setKSpaceCentreColumn ( (unsigned short)( m_dExactReadOutAsymBefore * pMrProt->kSpace().baseResolution() + 0.5 ));

    m_ADC.prep ( m_lReadOutColumns, (long)(pMrProt->rxSpec().effDwellTime(pSeqLim->getReadoutOSFactor())[m_lContrast]) ); 



    // * ---------------------------------------------------------------------- *
    // * Consider whether the echo should be shifted to the gradient raster     *
    // * time                                                                   *
    // * ---------------------------------------------------------------------- *
    if ( m_bEchoOnGRT )  {

        // * ------------------------------------------------------------------ *
        // * Artificial shift of the ADC in order to locate the echo maximum on *
        // * the gradient raster time.                                          *
        // * _----------------------------------------------------------------- *
        m_lADCShift = RoundUpGRT (static_cast<long>( round (m_dEffectiveDW_us * m_lReadOutColumnsBefore, 6) ) )
                                - static_cast<long>( round (m_dEffectiveDW_us * m_lReadOutColumnsBefore, 6) );
    } else {
        m_lADCShift = 0;
    }




    // * ---------------------------------------------------------------------- *
    // * Calculate readout gradient and export its moment to k-space center in  *
    // * order to select the ROP gradient outside of the SBB appropriately      *
    // * ---------------------------------------------------------------------- *
    double dMomentRO = fGSLGetROMoment ( pMrProt ) * (m_dExactReadOutAsymBefore + m_dExactReadOutAsymAfter);
    m_RO.setAmplitude ( m_eROPolarity * dMomentRO / m_ADC.getDuration() );

    if ( fabs(m_RO.getAmplitude()) > m_RO.getMaxMagnitude() )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Readout gradient amplitude exceeded", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Max. rise up/down time using the maximum allowed gradient amplitude.   *
    // * Used to limit ramp times to a reasonable length                        *
    // * ---------------------------------------------------------------------- *
    long lMaxRampTime = RoundUpGRT ( static_cast<long>( round (m_RO.getMaxMagnitude() * m_RO.getMinRiseTime(), 6)));



    // * ---------------------------------------------------------------------- *
    // * Calculate ramp up time or use a requested ramp up time if          *
    // * specified                                                          *
    // * ---------------------------------------------------------------------- *
    if ( m_lRequestedRampUpTime >= 0 )  {

        // * ------------------------------------------------------------------ *
        // * Use the externally defined ramp up time                            *
        // * ------------------------------------------------------------------ *
        m_RO.setRampUpTime ( m_lRequestedRampUpTime );
    } else {

        // * ------------------------------------------------------------------ *
        // * Internally calculte the ramp up time                               *
        // * ------------------------------------------------------------------ *
        long lActualRampUpTime = RoundUpGRT (static_cast<long>(fabs(m_RO.getAmplitude() * m_RO.getMinRiseTime() * 
                                                                    m_dFactorForPixelSizeRO)));    // * Rise time for actual amplitude * 

        m_RO.setRampUpTime ( minimum( lMaxRampTime, lActualRampUpTime ) );
    }

    m_lRampUpTime = m_RO.getRampUpTime();



    // * ---------------------------------------------------------------------- *
    // * Calculate ramp down time or use a requested ramp up time if specified  *
    // * ---------------------------------------------------------------------- *
    if ( m_lRequestedRampDownTime >= 0 )  {

        // * ------------------------------------------------------------------ *
        // * Use the externally defined ramp down time                          *
        // * ------------------------------------------------------------------ *
        m_RO.setRampDownTime ( m_lRequestedRampDownTime );
    } else {

        // * ------------------------------------------------------------------ *
        // * Internally calculte the ramp down time                             *
        // * ------------------------------------------------------------------ *
        long lActualRampDownTime = RoundUpGRT (static_cast<long>(fabs(m_RO.getAmplitude() * m_RO.getMinRiseTime() * 
                                                                      m_dFactorForPixelSizeRO)));  // * Rise time for actual amplitude *

        m_RO.setRampDownTime ( minimum( lMaxRampTime, lActualRampDownTime ) );
    }

    m_lRampDownTime = m_RO.getRampDownTime();


    // * ---------------------------------------------------------------------- *
    // * ADC shifts due to rounding to the gradient raster time have to be      *
    // * considered when calculating the duration of the readout gradient.      *
    // * ---------------------------------------------------------------------- *
    m_RO.setDuration ( RoundUpGRT ( static_cast<long>(m_RO.getRampUpTime() + m_lADCShift + m_ADC.getDuration()) ) );



    // * ---------------------------------------------------------------------- *
    // * Prepare the readout gradient                                           *
    // * ---------------------------------------------------------------------- *
    if ( ! m_RO.prep() )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Preparation of readout gradient failed.", ptModule);
        }
        setNLSStatus ( m_RO.getNLSStatus() );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Shift of the ADC related to the start time of the readout gradient.    *
    // * Sometimes necessary since the gradient raster time is 10 us and the    *
    // * "ADC raster time" is 1 us. Used e.g. for TrueFISP                      *
    // * ---------------------------------------------------------------------- *
    if ( m_bSymmetricTiming && (m_dReadOutAsymBefore == 0.5) )  {
        m_lADCShift += maximum ( (long)((m_RO.getFlatTopTime() - m_ADC.getRoundedDuration()) / 2) , 0L);
    }



    // * ---------------------------------------------------------------------- *
    // * Calculate the readout moment that has to be prephased by an ROP        *
    // * gradient                                                               *
    // * ---------------------------------------------------------------------- *
    m_dROMomentToKSpaceCenter = ( 0.5 * m_RO.getRampUpTime() + m_lADCShift + m_dEffectiveDW_us * m_lReadOutColumnsBefore ) * m_RO.getAmplitude();
    if ((m_bUseRampSampling && (m_lReadOutColumnsBefore == 0))) 
    {
        m_dROMomentToKSpaceCenter -= 0.5 * m_RO.getRampUpTime() * m_RO.getAmplitude();
    }

    // * ---------------------------------------------------------------------- *
    // * Calculate the readout moment that has to repahsed by an ROD gradient   *
    // * outside of the SBB                                                     *
    // * ---------------------------------------------------------------------- *
    m_dROMomentAfterKSpaceCenter = m_RO.getMomentumTOT() - m_dROMomentToKSpaceCenter;



    if ( m_lContrast == 0 )  {

        // * ------------------------------------------------------------------ *
        // * Preparation of the readout prephasing gradient, ROP is not used    *
        // * here                                                               *
        // * ------------------------------------------------------------------ *
        m_bUseROP = false;
        m_ROP.setAmplitude   (0.0);
        m_ROP.setRampUpTime  (0);
        m_ROP.setDuration    (0);
        m_ROP.setRampDownTime(0);

        if ( ! m_ROP.prep() )  {
            if ( pSeqLim->isContextNormal() )  {
                TRACE_PUT1(TC_INFO, TF_SEQ, "%s Preparation of readout gradient failed.", ptModule);
            }
            setNLSStatus ( m_RO.getNLSStatus() );
            return ( false );
        }

        m_dROMomentOut = m_dROMomentAfterKSpaceCenter;

    } else {

        if ( m_dROMomentToKSpaceCenter * m_dROMomentIn > m_dEpsilon )  {

            // * -------------------------------------------------------------- *
            // * Preceding and actual readout gradient have the same polarity   *
            // * -------------------------------------------------------------- *

            // * -------------------------------------------------------------- *
            // * Calcution of the shortest possible prephasing gradient in      *
            // * readout direction considering calculation limits               *
            // * -------------------------------------------------------------- *
            m_bUseROP = true;
            m_ROP.prepSymmetricTOTShortestTime( (-1.0) * (m_dROMomentToKSpaceCenter + m_dROMomentIn) * m_dFactorForPixelSizeRO );
            m_ROP.prepMomentumTOT ( (-1.0) * (m_dROMomentToKSpaceCenter + m_dROMomentIn) );

        } else {

            // * -------------------------------------------------------------- *
            // * Preceding and actual readout gradient have different           *
            // * polarities                                                     *
            // * -------------------------------------------------------------- *
            double dMomentToRefocus = (-1.0) * (m_dROMomentIn + m_dROMomentToKSpaceCenter);


            // * -------------------------------------------------------------- *
            // * The moment to refocus and readout amplitude have the same      *
            // * polarity                                                       *
            // *   ==> increase duration of the readout gradient                *
            // * -------------------------------------------------------------- *
            if ( m_eROPolarity * dMomentToRefocus > m_dEpsilon )  {

                dMomentToRefocus = round ( dMomentToRefocus, 6 );
                m_RO.setAmplitude ( round ( m_RO.getAmplitude(), 6 ) );
                long lDeltaRODuration = RoundDownGRT (fabs(dMomentToRefocus / m_RO.getAmplitude()));
                m_RO.setDuration (m_RO.getDuration() + lDeltaRODuration);


                // * ---------------------------------------------------------- *
                // * Start of ADC has be delayed                                *
                // * ---------------------------------------------------------- *
                m_lADCShift += lDeltaRODuration;

                if ( ! m_RO.prep() )  {
                    if ( pSeqLim->isContextNormal() )  {
                        TRACE_PUT1(TC_INFO, TF_SEQ, "%s Preparation of readout gradient failed.", ptModule);
                    }
                    setNLSStatus ( m_RO.getNLSStatus() );
                    return ( false );
                }

                // * ---------------------------------------------------------- *
                // * Re-evaluate readout moments before k-space center          *
                // * ---------------------------------------------------------- *
                m_dROMomentToKSpaceCenter += lDeltaRODuration * m_RO.getAmplitude();


                // * ---------------------------------------------------------- *
                // * Re-evaluating the moment that has to be prephased          *
                // * ---------------------------------------------------------- *
                dMomentToRefocus = (-1.0) * (m_dROMomentIn + m_dROMomentToKSpaceCenter);


                // * ---------------------------------------------------------- *
                // * Prephase residual readout moment                           *
                // * ---------------------------------------------------------- *
                m_bUseROP = true;
                m_ROP.prepSymmetricTOTShortestTime( dMomentToRefocus * m_dFactorForPixelSizeRO );
                m_ROP.prepMomentumTOT ( dMomentToRefocus );


            } else {

                // * ---------------------------------------------------------- *
                // * The moment to refocus and readout amplitude have opposite  *
                // * polarity                                                   *
                // *   ==> an extra prephasing gradient is required             *
                // * ---------------------------------------------------------- *
                if ( m_eROPolarity * dMomentToRefocus < (-1.0 * m_dEpsilon) )  {

                    // * ------------------------------------------------------ *
                    // * Use ROP gradient to balance the readout gradient       *
                    // * moment considering calculation limits                  *
                    // * ------------------------------------------------------ *
                    m_bUseROP = true;
                    m_ROP.prepSymmetricTOTShortestTime( dMomentToRefocus * m_dFactorForPixelSizeRO );
                    m_ROP.prepMomentumTOT ( dMomentToRefocus );

                } else {  // dMomentToRefocus == 0

                    // * ------------------------------------------------------ *
                    // * Preparation of the readout prephasing gradient, ROP is *
                    // * not used here                                          *
                    // * ------------------------------------------------------ *
                    m_bUseROP = false;
                    m_ROP.setAmplitude   (0.0);
                    m_ROP.setRampUpTime  (0);
                    m_ROP.setDuration    (0);
                    m_ROP.setRampDownTime(0);

                    if ( ! m_ROP.prep() )  {
                        if ( pSeqLim->isContextNormal() )  {
                            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Preparation of readout gradient failed.", ptModule);
                        }
                        setNLSStatus ( m_RO.getNLSStatus() );
                        return ( false );
                    }

                }

            }

        }

        m_dROMomentOut = m_dROMomentIn + m_ROP.getMomentumTOT() + m_RO.getMomentumTOT();
    }




    // * ---------------------------------------------------------------------- *
    // * A fixed ramp up time of the RO gradient has been requested but the ROP *
    // * gradient has to precede the readout blcok.                             *
    // * ---------------------------------------------------------------------- *
    if ( m_bUseROP && (m_lRequestedRampUpTime >= 0) )  {
        if ( pSeqLim->isContextNormal() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s ROP gradient must not be used in combination with a requested RO ramp up time", ptModule);
        }
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Calculate the readout moment considering the flat top, only            *
    // * ---------------------------------------------------------------------- *
    m_dROMomentFlatTop = m_RO.getFlatTopTime() * m_RO.getAmplitude();



    // * ---------------------------------------------------------------------- *
    // * Calculation of the start times                                         *
    // * ---------------------------------------------------------------------- *
    m_lStartTimeADC = m_lADCShift + m_ROP.getTotalTime() + m_RO.getRampUpTime();

    if ((pMrProt->kSpace().dimension() == SEQ::DIM_3) &&  
    pMrProt->kSpace().trajectory() == SEQ::TRAJECTORY_RADIAL && 
    pMrProt->kSpace().asymmetricEchoAllowed() &&
    (m_lReadOutColumnsBefore == 0)) 
    {
        m_lStartTimeADC -= m_RO.getRampUpTime();
    }

    if ( ! m_bRampUpInSBB )  {
        if ( m_bUseROP )  {
            m_lStartTimeADC -= m_ROP.getRampUpTime();
        } else {
            m_lStartTimeADC -= m_RO.getRampUpTime();
        }
    }

    m_lStartTimeRO  = m_ROP.getTotalTime(); 
    if ( ! m_bRampUpInSBB )  {
        if ( m_bUseROP )  {
            m_lStartTimeRO -= m_ROP.getRampUpTime();
        } else {
            m_lStartTimeRO -= m_RO.getRampUpTime();
        }
    }


    // * ------------------------------------------------------------------ *
    // * Calculate the duration of an additional delay that ensures the     *
    // * minimum delay between consecutive ADCs.                            *
    // * ------------------------------------------------------------------ *
    if ( m_lContrast > 0 )  {

        // * -------------------------------------------------------------- *
        // * Time between the end of the previos ADC and the current one    *
        // * -------------------------------------------------------------- *
        long lTimeBetweenADCs = m_lTimeBetweenADCsPreviousRO + m_lStartTimeADC;
        if ((m_bUseRampSampling && (m_lReadOutColumnsBefore == 0))) 
        {
            lTimeBetweenADCs += m_RO.getRampUpTime();
        }

        if ( lTimeBetweenADCs < m_lMinTimeBetweenADCs )  {
            m_lLeadTimeBetweenADCs = RoundUpGRT (m_lMinTimeBetweenADCs - lTimeBetweenADCs);
        } else {
            m_lLeadTimeBetweenADCs = 0;
        }

    } else {
        m_lLeadTimeBetweenADCs = 0;
    }


    // * ---------------------------------------------------------------------- *
    // * Start times with additional lead time                                  *
    // * ---------------------------------------------------------------------- *
    m_lStartTimeADC += m_lLeadTimeBetweenADCs;
    m_lStartTimeRO  += m_lLeadTimeBetweenADCs;
    m_lStartTimeROP =  m_lLeadTimeBetweenADCs;




    // * ------------------------------------------------------------------ *
    // * Calculate the time from the begin of the ADC to the echo maximum   *
    // * ------------------------------------------------------------------ *
    m_lTimeToEcho = m_lStartTimeADC + static_cast<long>( round (m_dEffectiveDW_us * m_lReadOutColumnsBefore, 6) );




    // * ------------------------------------------------------------------ *
    // * Calculate the duration per request                                 *
    // * ------------------------------------------------------------------ *
    m_lSBBDurationPerRequest_us = m_lLeadTimeBetweenADCs;
    if ( m_bUseROP )  {

        m_lSBBDurationPerRequest_us += m_ROP.getFlatTopTime() + m_ROP.getRampDownTime() + m_RO.getDuration();

        if ( m_bRampUpInSBB )  {  m_lSBBDurationPerRequest_us += m_ROP.getRampUpTime();  }

    } else {

        m_lSBBDurationPerRequest_us += m_RO.getFlatTopTime();

        if ( m_bRampUpInSBB )  {  m_lSBBDurationPerRequest_us += m_RO.getRampUpTime();  }

    }
    

    if ( m_bRampDownInSBB )  {
        m_lSBBDurationPerRequest_us += m_RO.getRampDownTime();
    }



    setPrepared();

    return ( true );

}



bool SeqBuildBlockReadOut::check (MrProt*, SeqLim* pSeqLim, SeqExpo*)
{

    // * ---------------------------------------------------------------------- *
    // * Check readout prephasing gradient                                      *
    // * ---------------------------------------------------------------------- *
    if (! m_ROP.check() )  {
        if ( pSeqLim->isContextNormal() )  {
            setNLSStatus (m_ROP.getNLSStatus(), "m_ROP.check");
        } else {
            setNLSStatus (m_ROP.getNLSStatus());
        }
        return ( false );
    }


    // * ---------------------------------------------------------------------- *
    // * Check readout gradient                                                 *
    // * ---------------------------------------------------------------------- *
    if (! m_RO.check() )  {
        if ( pSeqLim->isContextNormal() )  {
            setNLSStatus (m_RO.getNLSStatus(), "m_RO.check");
        } else {
            setNLSStatus (m_RO.getNLSStatus());
        }
        return ( false );
    }


    return ( true );

}



bool SeqBuildBlockReadOut::run (MrProt*, SeqLim*, SeqExpo*, sSLICE_POS*)
{

    // * ------------------------------------------------------------------ *
    // * Check successful preparation of SBBReadOut                         *
    // * ------------------------------------------------------------------ *
    if ( ! isPrepared() )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s SBBReadOut not prepared, cannot execute", ptModule);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }



    // * ------------------------------------------------------------------ *
    // * Check whether frequency and phase have been updated prior to       *
    // * execution of SBBReadOut::run(...).                                 *
    // * ------------------------------------------------------------------ *
    if ( ! m_bFreqPhaseUpdated )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s Frequency / phase has not been updated, cannot execute", ptModule);
        setNLSStatus ( SBB_ERROR );
        return ( false );
    }
    m_bFreqPhaseUpdated = false;



    // * ---------------------------------------------------------------------- *
    // * Set Mdh flags                                                          *
    // * ---------------------------------------------------------------------- *
    if ( m_eROPolarity == Positive )  {
        m_ADC.Mdh.deleteFromEvalInfoMask ( MDH_REFLECT );      /*! EGA-Any !*/
    } else {
        m_ADC.Mdh.addToEvalInfoMask      ( MDH_REFLECT );      /*! EGA-Any !*/
    }



    // * ---------------------------------------------------------------------- *
    // * Execute timing                                                         *
    // * ---------------------------------------------------------------------- *
    fRTEI(m_lStartTime + m_lStartTimeADC
                               , &m_ADCSet,         0,    &m_ADC,           0,          0,          0,          0);
    fRTEI(m_lStartTime + m_lStartTimeADC + m_ADC.getRoundedDuration() 
                               , &m_ADCNeg,         0,          0,          0,          0,          0,          0);

    if ( m_bUseROP )  {
        fRTEI(m_lStartTime + m_lStartTimeROP
                               ,         0,         0,          0,          0,     &m_ROP,          0,          0);
    }
    fRTEI(m_lStartTime + m_lStartTimeRO
                               ,         0,         0,          0,          0,      &m_RO,          0,          0);


    return ( true );
}





bool SeqBuildBlockReadOut::updateGradPerf (enum SEQ::Gradients eGradientMode)
{

    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for read out gradient                         *
    // * ---------------------------------------------------------------------- *
    m_ROP.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, GradGroupROP) );
    m_ROP.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, GradGroupROP) );

    m_RO.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, GradGroupRO) );
    m_RO.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, GradGroupRO) );


    return ( true );

}


bool SeqBuildBlockReadOut::getROAmplitudeRange (double *dMinAmpl, double *dMaxAmpl)
{

    if ( isPrepared() )  {

        *dMinAmpl = m_dFactorForUpperPixelSizeRO * m_RO.getAmplitude();
        *dMaxAmpl = m_dFactorForPixelSizeRO      * m_RO.getAmplitude();
        return ( true );

    } else {

        *dMinAmpl = 0.0;
        *dMaxAmpl = 0.0;
        return ( false );

    }

}


long SeqBuildBlockReadOut::getlROFlatToEcho (void) const
{
    return ( static_cast<long>(m_lADCShift + m_dEffectiveDW_us * m_lReadOutColumnsBefore ) );
}


void SeqBuildBlockReadOut::increaseADCPhase (double dPhaseOffset)
{
    m_ADCSet.increasePhase ( dPhaseOffset );    /*! EGA-05 !*/
    m_ADCNeg.decreasePhase ( dPhaseOffset );    /*! EGA-05 !*/
}




double SeqBuildBlockReadOut::getdPartialFourierRead (void) const
{
    return ( m_dExactReadOutAsymBefore + m_dExactReadOutAsymAfter );
}



void SeqBuildBlockReadOut::getdReadOutAsym (double & dAsymBefore, double & dAsymAfter) const
{
    dAsymBefore = m_dExactReadOutAsymBefore;
    dAsymAfter  = m_dExactReadOutAsymAfter;
}


double SeqBuildBlockReadOut::round (double dValue, unsigned long ulPrecision)
{
    double dPrec = 1.0;

    switch ( ulPrecision )  {
        case 0:     dPrec = 1.0e0;      break;
        case 1:     dPrec = 1.0e1;      break;
        case 2:     dPrec = 1.0e2;      break;
        case 3:     dPrec = 1.0e3;      break;
        case 4:     dPrec = 1.0e4;      break;
        case 5:     dPrec = 1.0e5;      break;
        case 6:     dPrec = 1.0e6;      break;
        case 7:     dPrec = 1.0e7;      break;
        case 8:     dPrec = 1.0e8;      break;

        default: for ( unsigned long lI = 0; lI < ulPrecision; lI++ )  {
            dPrec *= 10.0;
        }
    }

    return ( floor(dValue * dPrec) / dPrec );
}




