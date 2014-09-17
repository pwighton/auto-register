//   Read the documentation to learn more about C++ code generator
//   versioning.

//	-----------------------------------------------------------------------------
//	  Copyright (C) Siemens AG 1998  All Rights Reserved.
//	-----------------------------------------------------------------------------
//
//	 Project: NUMARIS/4
//	    File: \n4\pkg\MrServers\MrImaging\seq\Kernels\SBBGREKernel.cpp
//	 Version:
//	  Author: Clinical
//	    Date: n.a.
//	    Date: n.a.
//
//	    Lang: C++
//
//	 Descrip: MR::MrServers::MrImaging::seq::Kernels
//
//	 Classes:
//
//	-----------------------------------------------------------------------------




// * -------------------------------------------------------------------------- *
// * Includes                                                                   *
// * -------------------------------------------------------------------------- *
#include "MrServers/MrImaging/seq/Kernels/SBBGREKernel.h"
#include "MrServers/MrImaging/seq/Kernels/SBBReadOut.h"
#include "MrServers/MrImaging/libSeqUtil/libSeqUtil.h"

#include "MrServers/MrProtSrv/MrProt/KSpace/MrKSpace.h"
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrRXSpec.h"
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrSysSpec.h"
#include "MrServers/MrMeasSrv/SeqIF/SeqBuffer/SeqLim.h"
#include "MrServers/MrMeasSrv/SeqIF/SeqBuffer/SeqExpo.h"

#if defined SEND_TTC_DATA_IN_MDH
#include "MrServers/MrImaging/seq/common/AdvMRA/MdhDef.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/RTController.h"
#endif


#define DEBUG_ORIGIN                        DEBUG_SBB



SBBGREKernel::SBBGREKernel (SBBList* pSBBList)
      : SeqBuildBlock                       (pSBBList),
        m_dEpsilon                          (0.0000001),
        m_RF_WE121                          (NULL, FALSE, 180.0, 3, "Excit"),
        m_RF_WE121_short                    (NULL, FALSE,  90.0, 3, "Excit"),
        m_RF_WE11_short                     (NULL, FALSE,  45.0, 2, "Excit"),
        m_RF_WE11_medium                    (NULL, FALSE,  75.0, 2, "Excit"),
        m_RF_WE11_90                        (NULL, FALSE,  90.0, 2, "Excit"),
        m_RF_WE11                           (NULL, FALSE, 180.0, 2, "Excit"),
        m_pRF_WE                            (NULL),
        m_pRF_Exc                           (NULL),
        m_PhaseStabScan                     (NULL),
        m_aRO                               (NULL),
        m_dRFSpoilPhase                     (0.0),
        m_dMomentum_TB_3D                   (0.0),
        m_dMomentum_TB_3DR                  (0.0),
        m_dMomentROP                        (0.0),
        m_bRampsOutsideOfEventBlock         (true),
        m_bUsePERewinder                    (true),
        m_bBalancedGradMomentSliceSelect    (false),
        m_bSendOscBit                       (false),
        m_bSendWakeUpBit                    (false),
        m_lStartTimeWakeUp                  (500),
        m_bUpdateRotationMatrix             (false),
        m_bPerformNegativeTEFillCheck       (true),
        m_bEchoOnGRT                        (false),
        m_bSuppressTraces                   (false),
        m_eWEMode                           (WE_1_2_1_180_Deg),
        m_bUseWEOffsetFrequency             (false),
        m_dWEOffsetFrequency                (0.0),
        m_dWEBandwidthTimeProduct           (2.0),
        m_dRelRODephasingMoment             (0.5),
        m_dRelSSDephasingMoment             (0.0),
        m_eReadOutGradType                  (Spoiler),
        m_lNumberOfEchoes                   (1),
        m_eROPolarity                       (Bipolar),
        m_adROAsymmetryBefore               (NULL),
        m_adROAsymmetryAfter                (NULL),
        m_lFirstLineToMeasure               (0),
        m_lKSpaceCenterLine                 (0),
        m_lLastLineToMeasure                (0),
        m_lLineNumber                       (0),
        m_lFirstPartToMeasure               (0),
        m_lKSpaceCenterPartition            (0),
        m_lLastPartToMeasure                (0),
        m_lPartNumber                       (0),
        m_lExcitationAllTime                (0),
        m_lTEFromExcitation                 (0),
        m_lPrephaseTime                     (0),
        m_lRephaseTime                      (0),
        m_lReadoutTime                      (0),
        m_lOverTime                         (0),
        m_lPrephaseFillTimeRead             (0),
        m_lPrephaseFillTimeSlice            (0),
        m_lPrephaseFillTimePhase            (0),
        m_eLimitingPrephaseTime             (Read),
        m_eLimitingRephaseTime              (Read),
        m_lMinDelayBetweenADCs              (0),
        m_lMinDelayBetweenADCAndRF          (0),
        m_lMinDelayBetweenRFAndADC          (0),
        m_dFactorForPixelSizeRO             (1.0),
        m_dFactorForUpperPixelSizeRO        (1.0),
        m_dFactorForPixelSizePE             (1.0),
        m_dFactorForPixelSize3D             (1.0),
        m_dFactorForSliceThickness          (1.0),
        m_alTEMin                           (NULL),
        m_alTEFill                          (NULL),
        m_lTRMin                            (0),
        m_lTRFill                           (0),
        m_lTRFillEnd                        (0),
        m_bLastScanInConcat                 (false),
        m_bLastScanInMeas                   (false),
        m_bFlowCompRead                     (false),
        m_bFlowCompPhase                    (false),
        m_bFlowCompSliceSelect              (false),
        m_bFlowCompPartitionEncode          (false),
        m_bAvoidAcousticResonances          (false),
        m_ulUTIdent                         (0),
        m_bsuccessfulAllocation             (true),
        m_pRF_Pre                           (NULL),
        m_sGFTss                            ("GRE-GFTss"),
        m_lRampTimePreviousEB_us            (0),
        m_iAdjUpdate                        (RT_MDSUPDATE_ADJ_NONE),
        m_lTimeToFirstPreScanRF             (0),
        m_lPreScanFillTime                  (0),
        m_lDelayAdjustSSPosition            (0),
        m_padFlipAngleArray                 (NULL),
        m_uiNmbrOfFlipAngles                (2),


        m_bFoundFirstTimeStamp              (false),
        m_dLeadTimeMSec                     (0.0),
        m_dLastTimeToCenterPlayedOutSec     (-1.0),
        m_bUseRampSampling                  (false)
{
  
  long lI;    // * Loop counter *

  // * Set the ident string of the SBB *
  setIdent ("SBBKernel");


  // * Memory allocation for multiple contrasts *
  m_bsuccessfulAllocation = true;


    if ( (m_aRO = new SeqBuildBlockReadOut [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;

    if ( (m_alTEMin             = new long   [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
    if ( (m_alTEFill            = new long   [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
    if ( (m_abNegativeTE        = new bool   [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
    if ( (m_adROAsymmetryBefore = new double [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
    if ( (m_adROAsymmetryAfter  = new double [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;



    // * ---------------------------------------------------------------------- *
    // * Initialization of arrays                                               *
    // * ---------------------------------------------------------------------- *
    if ( m_bsuccessfulAllocation ) {
        for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_alTEMin[lI]        = 0;
        for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_alTEFill[lI]       = 0;
        for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_abNegativeTE[lI]   = true;
        for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_adROAsymmetryBefore[lI] = 0.5;
        for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_adROAsymmetryAfter[lI]  = 0.5;
    }


    m_sP_3D.setIdent ("Slice refocusing gradient");
    m_sP_3D_FC.setIdent ("Slice refocusing gradient with flow compensation");
    m_sTB_3DR.setIdent("3D phase encoding rewinder");
    m_sTB_PER.setIdent("Phase encoding rewinder");
    m_sP_ROP.setIdent ("Read out prephasing gradient");
    m_sP_ROD.setIdent ("Read out dephasing gradient");

}



SBBGREKernel::~SBBGREKernel()
{
    delete [] m_aRO;  
    delete [] m_alTEMin;
    delete [] m_alTEFill;
    delete [] m_abNegativeTE;
    delete [] m_adROAsymmetryBefore;
    delete [] m_adROAsymmetryAfter;
}



// **********************************************************************************
// * Name        : setRFPulseForExcitation                                          *
// *                                                                                *
// * Class       : SBBGREKernel                                                     *
// *                                                                                *
// * Description : Initialization of the member m_pRF_Exc.                          *
// *               m_pRF_Exc is an instance of sRF_PULSE.                           *
// *                                                                                *
// * Return      : true in case of successful initialization else false             *
// *                                                                                *
// **********************************************************************************
bool SBBGREKernel::setRFPulseForExcitation (sRF_PULSE* pSRF)
{

    static const char *ptModule = "SBBGREKernel::setRFPulseForExcitation";


    resetTimingCalculated();
    resetPrepared();    // SBBGREKernel is no longer prepared due the selection of a new rf-pulse

    if ( pSRF == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF == NULL");
        return ( false );
    }

    if ( pSRF->isTypeExcitation() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF->isTypeExcitation() == false");
        return ( false );
    }

    if ( pSRF->isPrepared() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF->isPrepared()==false");
        return ( false );
    }

    m_pRF_Exc = pSRF;

    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Sets a pointer to a sRF_PULSE instance of the     *
//	* excitation pulses of the gradient echo kernel and *
//	* the prescan (TrueFisp)                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC003D9
bool SBBGREKernel::setRFPulseForExcitation (sRF_PULSE* pSRF_Exc, sRF_PULSE* pSRF_Pre)
{
  // **********************************************************************************
  // * Name        : setRFPulseForExcitation                                          *
  // *                                                                                *
  // * Class       : SBBGREKernel                                                     *
  // *                                                                                *
  // * Description : Initialization of the member m_pRF_Exc and m_pRF_Pre.            *
  // *               m_pRF_Exc and m_pRF_Pre are instances of sRF_PULSE.              *
  // *               This function is used for TrueFisp sequences.                    *
  // *                                                                                *
  // * Return      : true in case of successful initialization else false             *
  // *                                                                                *
  // **********************************************************************************


    static const char *ptModule = "SBBGREKernel::setRFPulseForExcitation";


    resetTimingCalculated();
    resetPrepared();    // SBBGREKernel is no longer prepared due the selection of a new rf-pulse



    // * -------------------------------------------------------------------------- *
    // * Check pSRF_Exc                                                             *
    // * -------------------------------------------------------------------------- *
    if ( pSRF_Exc == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Exc == NULL");
        return ( false );
    }



    if ( pSRF_Exc->isTypeExcitation() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Exc->isTypeExcitation() == false");
        return ( false );
    }



    if ( pSRF_Exc->isPrepared() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Exc->isPrepared()==false");
        return ( false );
    }




    // * -------------------------------------------------------------------------- *
    // * Check pSRF_Pre                                                             *
    // * -------------------------------------------------------------------------- *
    if ( pSRF_Pre == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Pre == NULL");
        return ( false );
    }



    if ( pSRF_Pre->isTypeExcitation() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Pre->isTypeExcitation() == false");
        return ( false );
    }



    if ( pSRF_Pre->isPrepared() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Pre->isPrepared()==false");
        return ( false );
    }




    m_pRF_Exc = pSRF_Exc;
    m_pRF_Pre = pSRF_Pre;

    return ( true );

}





bool SBBGREKernel::updateGradPerf (enum SEQ::Gradients eGradientMode)
{

    static const char *ptModule = "SBBGREKernel::updateGradPerf";


    // * ---------------------------------------------------------------------- *
    // * Check pointer m_pRF_WE                                                 *
    // * ---------------------------------------------------------------------- *
    if ( m_pRF_WE == NULL )  {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL");
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for slice select gradient                     *
    // * ---------------------------------------------------------------------- *
    double adSSExternalPulseMinRiseTime[3]  = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_SS),
                                                getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_SS),
                                                getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_SS) };

    double adSSExternalPulseMaxMagnitude[3] = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_SS),
                                                getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_SS),
                                                getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_SS) };

    double adSSWEPulseMinRiseTime[3]        = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_WE),
                                                getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_WE),
                                                getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_WE) };


    double adSSWEPulseMaxMagnitude[3]       = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_WE),
                                                getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_WE),
                                                getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_WE) };


    m_pRF_WE->setMinRiseTimes  ( adSSExternalPulseMinRiseTime , SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF );
    m_pRF_WE->setMaxMagnitudes ( adSSExternalPulseMaxMagnitude, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF );
    m_pRF_WE->setMinRiseTimes  ( adSSWEPulseMinRiseTime       , SBBBinomialPulses_GRAD_PERF_BINOMIAL    );
    m_pRF_WE->setMaxMagnitudes ( adSSWEPulseMaxMagnitude      , SBBBinomialPulses_GRAD_PERF_BINOMIAL    );

    m_pRF_WE->setMinRiseTimeGSWD  (getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_SS), SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    m_pRF_WE->setMaxMagnitudeGSWD (getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_SS), SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    m_pRF_WE->setMinRiseTimeGSWD  (getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_WE), SBBBinomialPulses_GRAD_PERF_BINOMIAL);
    m_pRF_WE->setMaxMagnitudeGSWD (getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_WE), SBBBinomialPulses_GRAD_PERF_BINOMIAL);



    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for slice refocussing gradients               *
    // * ---------------------------------------------------------------------- *
    m_sP_3D .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_SS) );
    m_sP_3D .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_SS) );

    m_sP_3D_FC.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_SS) );
    m_sP_3D_FC.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_SS) );



    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for phase encoding gradient                   *
    // * ---------------------------------------------------------------------- *
    double adPEMinRiseTime[3]  = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_PE),
                                   getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_PE),
                                   getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_PE) };

    double adPEMaxMagnitude[3] = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_PE),
                                   getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_PE),
                                   getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_PE) };

    m_TB_PE.setMinRiseTimes  ( adPEMinRiseTime , SBBPhaseEncode::GradTable );
    m_TB_PE.setMaxMagnitudes ( adPEMaxMagnitude, SBBPhaseEncode::GradTable );

    m_TB_PE.setMinRiseTimeGSWD  (getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_PE), SBBPhaseEncode::GradTable);
    m_TB_PE.setMaxMagnitudeGSWD (getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_PE), SBBPhaseEncode::GradTable);


    m_sTB_PER.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_PE) );
    m_sTB_PER.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_PE) );



    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for partition encoding gradient               *
    // * ---------------------------------------------------------------------- *
    double ad3DMinRiseTime[3]  = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_3D),
                                   getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_3D),
                                   getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_3D) };

    double ad3DMaxMagnitude[3] = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_3D),
                                   getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_3D),
                                   getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_3D) };


    m_TB_3D.setMinRiseTimes  ( ad3DMinRiseTime , SBBPhaseEncode::GradTable );
    m_TB_3D.setMaxMagnitudes ( ad3DMaxMagnitude, SBBPhaseEncode::GradTable );

    m_TB_3D.setMinRiseTimeGSWD  (getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_3D), SBBPhaseEncode::GradTable);
    m_TB_3D.setMaxMagnitudeGSWD (getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_3D), SBBPhaseEncode::GradTable);


    m_sTB_3DR.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_3D) );
    m_sTB_3DR.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_3D) );



    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for read out gradient                         *
    // * ---------------------------------------------------------------------- *
    m_sP_ROP   .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );
    m_sP_ROP   .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );

    m_sP_ROP_FC.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );
    m_sP_ROP_FC.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );

    m_sP_ROD   .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROD) );
    m_sP_ROD   .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROD) );


    double adROMinRiseTime[3]  = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_RO),
                                   getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_RO),
                                   getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_RO) };

    double adROMaxMagnitude[3] = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_RO),
                                   getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_RO),
                                   getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_RO) };


    double adROPMinRiseTime[3]  = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_ROP),
                                    getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_ROP),
                                    getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_ROP) };

    double adROPMaxMagnitude[3] = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_ROP),
                                    getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_ROP),
                                    getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_ROP) };

    for ( long lI=0; lI < m_lNumberOfEchoes; lI++ )  {
        m_aRO[lI].setMinRiseTimes  ( adROMinRiseTime , SBBReadOut::GradGroupRO );
        m_aRO[lI].setMaxMagnitudes ( adROMaxMagnitude, SBBReadOut::GradGroupRO );

        m_aRO[lI].setMinRiseTimeGSWD  ( getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME, SBBGREKernel_GRAD_GROUP_RO), SBBReadOut::GradGroupRO );
        m_aRO[lI].setMaxMagnitudeGSWD ( getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME, SBBGREKernel_GRAD_GROUP_RO), SBBReadOut::GradGroupRO );


        m_aRO[lI].setMinRiseTimes  ( adROPMinRiseTime , SBBReadOut::GradGroupROP );
        m_aRO[lI].setMaxMagnitudes ( adROPMaxMagnitude, SBBReadOut::GradGroupROP );

        m_aRO[lI].setMinRiseTimeGSWD  ( getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME, SBBGREKernel_GRAD_GROUP_ROP), SBBReadOut::GradGroupROP );
        m_aRO[lI].setMaxMagnitudeGSWD ( getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME, SBBGREKernel_GRAD_GROUP_ROP), SBBReadOut::GradGroupROP );
    }




    return ( true );

}



bool SBBGREKernel::calculateTiming (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{

    // **************************************************************************
    // * Name        : getcalculateTiming                                       *
    // *                                                                        *
    // * Class       : SBBGREKernel                                             *
    // *                                                                        *
    // * Description : Calculation of the timing of the GRE kernel              *
    // *                                                                        *
    // * Return      : bool                                                     *
    // *                                                                        *
    // **************************************************************************

    static const char *ptModule = "SBBGREKernel::calculateTiming";

    static GPAProxy   theGPA;                   // * Information about the gradient switching time            *


    long   lI;                                  // * Loop counter                                             *
    double dROMomentBeforePhaseStab   = 0.0;



    if ( mIsTRrun )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s running\n", ptModule);
    }


    // * ---------------------------------------------------------------------- *
    // * Debug: running message and imported parameters                         *
    // * ---------------------------------------------------------------------- *
    printImports ( pMrProt, pSeqLim, ptModule );



    // * ---------------------------------------------------------------------- *
    // * Initialize variables                                                   *
    // * ---------------------------------------------------------------------- *
    if (! calculateInit(pMrProt, pSeqLim, pSeqExpo) ) return ( false );


    
    // * ---------------------------------------------------------------------- *
    // * Check kernel settings                                                  *
    // * ---------------------------------------------------------------------- *
    if (! calculateCheckSetting( pMrProt, pSeqLim ) ) return ( false );




    // * ---------------------------------------------------------------------- *
    // * Get scaling factors for gradient timing calculation                    *
    // * ---------------------------------------------------------------------- *
    m_dFactorForPixelSizeRO      = m_pCalcLimits->getFactorForPixelSizeRO      ( pMrProt );
    m_dFactorForUpperPixelSizeRO = m_pCalcLimits->getFactorForUpperPixelSizeRO ( pMrProt );
    m_dFactorForPixelSizePE      = m_pCalcLimits->getFactorForPixelSizePE      ( pMrProt );
    if ( pMrProt->kSpace().dimension() == SEQ::DIM_3 )  {
        m_dFactorForPixelSize3D  = m_pCalcLimits->getFactorForPixelSize3D      ( pMrProt );
    }
    m_dFactorForSliceThickness   = m_pCalcLimits->getFactorForSliceThickness   ( pMrProt );




    // *----------------------------------------------------------------------- *
    // * Selects the excitation pulse                                           *
    // *----------------------------------------------------------------------- *
    if (! calculateSelectWEPulse(pMrProt, pSeqLim, pSeqExpo) ) return ( false );




    // * ---------------------------------------------------------------------- *
    // * Update of the gradient performance settings (FAST, NORMAL, WHISPER)    *
    // *                                                                        *
    // * CAUTION: Function updateGradientPerformance MUST be called after       *
    // *          function calculateSelectWEPulse to initilize m_pRF_WE.        *
    // * ---------------------------------------------------------------------- *
    if (! updateGradPerf ( pMrProt->gradSpec().mode() ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Prepare readout SBBs                                                   *
    // * ---------------------------------------------------------------------- *
    if (! calculatePrepRO (pMrProt, pSeqLim, pSeqExpo) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Calcution of the readout prephasing gradient                           *
    // * Flow compensation in readout direction is taken into account.          *
    // * ---------------------------------------------------------------------- *
    if (! calculateROP ( pMrProt, pSeqLim ) ) return ( false );



    // *----------------------------------------------------------------------- *
    // * Prepares the excitation pulse                                          *
    // *----------------------------------------------------------------------- *
    if (! calculatePrepWEPulse (pMrProt, pSeqLim, pSeqExpo) ) return ( false );



    // *----------------------------------------------------------------------- *
    // *                                                                        *
    // * Calculate the gradient timing for phase encoding, kz                   *
    // *                                                                        *
    // *----------------------------------------------------------------------- *
    if (! calculate3DGradientTables ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );



    // *----------------------------------------------------------------------- *
    // *                                                                        *
    // * Calculate gradient timing for phase encoding, ky                       *
    // *                                                                        *
    // *----------------------------------------------------------------------- *
    if (! calculatePEGradientTables ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Which prephasing process is time limiting?                             *
    // *                                                                        *
    // * Prephasing in slice select, phase encode, read out direction           *
    // *                                                                        *
    // * or the m_lADCZeroBefore[0] duration                                    *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    if (! calculatePrephaseTime ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Initialization of phase stabilization scans                            *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->phaseStabilize() )  {
        m_PhaseStabScan.setlLastKernelRORampDownTimeInsideSBB ( m_aRO[m_lNumberOfEchoes-1].getlRampDownTime() );
        m_PhaseStabScan.setps_First_Kernel_ADC   ( &(m_aRO[0].ADC())       );
        m_PhaseStabScan.setdFactorForPixelSizeRO ( m_dFactorForPixelSizeRO );
    

        dROMomentBeforePhaseStab = m_aRO[m_lNumberOfEchoes-1].getdROMomentOut();
    

        m_PhaseStabScan.setlTimeOfLastKernelADCInSBB   ( 0 );
        m_PhaseStabScan.setdROMomentBeforePhaseStabADC ( dROMomentBeforePhaseStab );

	    switch ( m_aRO[pMrProt->contrasts()-1].geteROPolarity() )  {
	        case SeqBuildBlockReadOut::Positive:    m_PhaseStabScan.setePolarity ( SeqBuildBlockPhaseStabScan::Positive );  break;
	        case SeqBuildBlockReadOut::Negative:    m_PhaseStabScan.setePolarity ( SeqBuildBlockPhaseStabScan::Negative );  break;
	        default:                                m_PhaseStabScan.setePolarity ( SeqBuildBlockPhaseStabScan::Undefined );
	    }



        // * ---------------------------------------------------------------------- *
        // * Calculate the RO ramp down time of the phase stabilization SBB to      *
        // * fit to ROD gradient of the kernel                                      *
        // * ---------------------------------------------------------------------- *
        switch (m_eReadOutGradType)  {
            case Spoiler:
                m_PhaseStabScan.setRORampDownTime ( RoundUpGRT( maximum ( m_PhaseStabScan.getROMaxMagnitude() * m_PhaseStabScan.getROMinRiseTime(),
                                                                          m_sP_ROD.getMaxMagnitude() * m_sP_ROD.getMinRiseTime() ) ) );
            break;

            case Constant:
            case Rewinder:
            case Symmetric_Rewinder:
                // * No RO ramp down time is specified by SBBGREKernel               *
                // * ==> SBBPhaseStabScan calculates the ramp down time by its own.  *
            break;
    
            default: 
                if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." ) )  return ( false );
        }

    
        m_PhaseStabScan.setLastRampDownOutsideSBB ( true );     // * Last gradient ramp is outside of the SBB *

    }

    if (! m_PhaseStabScan.prep ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );




    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Calculation of the read out dephasing gradient                         *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    if (! calculateROD( pMrProt ) ) return ( false );





    // * ---------------------------------------------------------------------- *
    // * read out gradient ramp down times                                      *
    // * ---------------------------------------------------------------------- *
    if (! calculateRORampDownTime( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );
  


    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Which process is time limiting for rephasing?                          *
    // *                                                                        *
    // * Rephasing in slice select or phase encode direction, the end of the    *
    // *                                                                        *
    // * read out dephasing gradient or the lADCZeroAfter duration ?            *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    if (! calculateRephaseTime ( pMrProt ) ) return ( false );




    // * ---------------------------------------------------------------------- *
    // *                                                                        *
    // * Calculation of the duration that gradient ramps are allowed to dangle  *
    // * out of the event block                                                 *
    // *                                                                        *
    // * ---------------------------------------------------------------------- *
    if (! calculateOverTime( pMrProt, pSeqLim ) ) return ( false );




    // * ---------------------------------------------------------------------- *
    // * Additional calculations                                                *
    // * ---------------------------------------------------------------------- *
    if (! calculateAddOn ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );




    // * ---------------------------------------------------------------------- *
    // * Calculation of the minimum TE times                                    *
    // * ---------------------------------------------------------------------- *
    if (! calculateTEMin(pMrProt) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Calculate the TE fill times                                            *
    // * ---------------------------------------------------------------------- *
    if (! calculateTEFill ( pMrProt ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Gradient ramps that are dangling out of the event block                *
    // * (m_lOverTime) do not contribute to the minimum TR time                 *
    // * ---------------------------------------------------------------------- *
    if (!  calculateTRMin( pMrProt ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Calculation of start times of the readout SBBs                         *
    // * ---------------------------------------------------------------------- *
    if (! calculateROStartTime (pMrProt, pSeqLim, pSeqExpo) )  return ( false );



    // * ---------------------------------------------------------------------- *
    // * Additional calculations after TE and TR calculation                    *
    // * ---------------------------------------------------------------------- *
    if (! calculateAddOnPost ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Check whether m_lStartTimeWakeUp is shorter that the minimum TR time   *
    // * and correct m_lStartTimeWakeUp if necessary.                           *
    // * ---------------------------------------------------------------------- *
    if ( m_lStartTimeWakeUp > m_lTRMin )  {
        m_lStartTimeWakeUp = m_lTRMin;
    }



    // * ---------------------------------------------------------------------- *
    // * Check for negative TE fill time                                        *
    // * ---------------------------------------------------------------------- *
    bool bNegativeTEOccured = false;


    if ( m_bPerformNegativeTEFillCheck )  {

        // * TE fill times of all previous contrasts *
        long lTotalTEFill = 0;

        // * ------------------------------------------------------------------ *
        // * Check contrasts                                                    *
        // *                                                                    *
        // * One has to use pMrProt->contrasts() instead of m_lNumberOfEchoes   *
        // * since dess has one contrast but m_lNumberOfEchoes = 2              *
        // * ------------------------------------------------------------------ *
        for ( lI = 0; lI < pMrProt->contrasts() /*m_lNumberOfEchoes*/; lI++ )  {

            if ( (m_alTEMin[lI] + lTotalTEFill) > pMrProt->te()[lI] )  {

                // * ---------------------------------------------------------- *
                // * Negative TE fill time occured                              *
                // * ---------------------------------------------------------- *
                m_abNegativeTE[lI] = true;
                bNegativeTEOccured = true;

            } else {
                m_abNegativeTE[lI] = false;
            }

            lTotalTEFill += m_alTEFill[lI];

        }

    }


    // * ---------------------------------------------------------------------- *
    // * Return NLS status if negative TE fill times occured                    *
    // * ---------------------------------------------------------------------- *
    if ( bNegativeTEOccured )  {

        if ( !pSeqLim->isContextPrepForBinarySearch() && !m_bSuppressTraces)  {
            if ( setNLSStatus ( SEQU__NEGATIV_TEFILL, ptModule, "Negative TE fill time" ) ) return ( false );
        } else {
            if ( setNLSStatus ( SEQU__NEGATIV_TEFILL ) ) return ( false );
        }

    }


    // *----------------------------------------------------------------------- *
    // * Debug: intermediate results                                            *
    // *----------------------------------------------------------------------- *
    if ( mIsTRint )  {
        TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Intermediate Results                              *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s theGPA.getGradSwitchTime()    = %11.6f [us]\n\n", ptModule, theGPA.getGradSwitchTime()    );
        if ( m_pRF_WE )  {
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_WE->getGSData()->getlTimeToFirstFlatTop()  = %6ld [us]\n\n"  , ptModule, m_pRF_WE->getGSData().getlTimeToFirstFlatTop()  );
        } else {
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "Invalid pointer detected" ) ) return ( false );
        }
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForPixelSizeRO         = %11.6f\n"     , ptModule, m_dFactorForPixelSizeRO       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForUpperPixelSizeRO    = %11.6f\n"     , ptModule, m_dFactorForUpperPixelSizeRO  );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForPixelSizePE         = %11.6f\n"     , ptModule, m_dFactorForPixelSizePE       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForPixelSize3D         = %11.6f\n"     , ptModule, m_dFactorForPixelSize3D       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForSliceThickness      = %11.6f\n\n"   , ptModule, m_dFactorForSliceThickness    );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lTEFromExcitation             = %6ld [us]\n"  , ptModule, m_lTEFromExcitation           );

        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    }


    // *----------------------------------------------------------------------- *
    // * Debug: results                                                         *
    // *----------------------------------------------------------------------- *
    printResults ( pSeqLim, ptModule );
  

    // * ---------------------------------------------------------------------- *
    // * Timing has successfully been calculated                                *
    // * ---------------------------------------------------------------------- *
    setTimingCalculated();



    if ( mIsTRend )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s finished\n", ptModule);
    }


    return ( true );  // Returning without error

}


// ******************************************************************************
// * Name        : calculateCheckSetting                                        *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Checks whether the parameter settings are valid to preceed   *
// *               in kernel calculation.                                       *
// *               The function sets a NLS status and returns a false if the    *
// *               settings are not valid.                                      *
// *                                                                            *
// * Return      : BOOL                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateCheckSetting (MrProt* pMrProt, SeqLim*)
{
  
    static const char *ptModule = "SBBGREKernel::calculateCheckSetting";


    // * ---------------------------------------------------------------------- *
    // * Check whether dynamic memory allocation has been successfully carried  *
    // * out. Otherwise, quit the calculation of the kernel timing              *
    // * ---------------------------------------------------------------------- *
    if ( m_bsuccessfulAllocation == false )  {
        if ( setNLSStatus(SEQU_ERROR, ptModule, "Can't allocate dynamic memory in SBBGREKernel\n") )    return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Check whether a valid rf-pulse has been defined                        *
    // * ---------------------------------------------------------------------- *
    if ( m_pRF_Exc == NULL )  {
        setNLSStatus(SEQU_ERROR, ptModule, "No valid rf-pulse defined.\n");
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Check whether spoiler on slice select axis                             *
    // * ---------------------------------------------------------------------- *
    if ( m_bBalancedGradMomentSliceSelect && ( m_dRelSSDephasingMoment != 0.0 ) ) {
        setNLSStatus(SEQU_ERROR, ptModule, "Slice select axis: No spoiler possible in case of a balanced gradient scheme. \n");
        return ( false );
    }

    if ( ( m_bUsePERewinder == false ) && ( m_dRelSSDephasingMoment != 0.0 ) ) {
        setNLSStatus(SEQU_ERROR, ptModule, "Slice select axis: No spoiler possible if gradient pulse has been switched off.\n");
        return ( false );
    }


    if ( ( m_eReadOutGradType == Symmetric_Rewinder ) && ( m_adROAsymmetryBefore[0] != 0.5 ) && ( m_adROAsymmetryAfter[0] != 0.5 ) )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Error: No symmetric ADC timing.\n");
        return ( false );
    }


    // * ---------------------------------------------------------------------- *
    // * Check whether the number of contrast is inside the valid range         *
    // * Otherwise, quit the calculation of the kernel timing                   *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->contrasts() > lMaxNumberOfContrasts )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Number of contrasts to large\n");
        return ( false );
    }



    return ( true );

}



// ******************************************************************************
// * Name        : calculateSetRFWEPulse                                        *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Sets the member m_pRF_WE to the selected water selective     *
// *               rf-pulse                                                     *
// *                                                                            *
// * Return      : BOOL                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateSelectWEPulse (MrProt* pMrProt, SeqLim*, SeqExpo* )
{

    static const char *ptModule = "SBBGREKernel::calculateSelectWEPulse";
  


    switch ( m_eWEMode )  {
        case WE_1_1_45_Deg:
            m_pRF_WE = &m_RF_WE11_short;
        break;

        case WE_1_1_75_Deg:

            // * -------------------------------------------------------------- *
            // * Define parameters for sinc pulses used within WE_1_1_75_Deg    *
            // * WE pulse                                                       *
            // * -------------------------------------------------------------- *
            if ( pMrProt->txSpec().excitation() == SEQ::EXCITATION_SLICE_SELECTIVE  )  {

                m_sRFWE1.setDuration             (       2000 );
                m_sRFWE1.setBandwidthTimeProduct ( 12.7 * 0.8 );
                m_sRFWE1.setSamples              (        500 );
                m_sRFWE1.setInitialPhase         (        0.0 );
                m_sRFWE1.setThickness            ( pMrProt->sliceSeries().aFront().thickness() * 0.8 );

                if(! m_RF_WE11_medium.registerBinomialRFPulse( 0, &m_sRFWE1 ) )   {
                    setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_medium.registerBinomialRFPulse failed.");
                    return ( false );
                }

                m_sRFWE2.setDuration             ( m_sRFWE1.getDuration()             );
                m_sRFWE2.setBandwidthTimeProduct ( m_sRFWE1.getBandwidthTimeProduct() );
                m_sRFWE2.setSamples              ( m_sRFWE1.getSamples()              );
                m_sRFWE2.setInitialPhase         ( m_sRFWE1.getInitialPhase()         );
                m_sRFWE2.setThickness            ( m_sRFWE1.getThickness()            );

                if(! m_RF_WE11_medium.registerBinomialRFPulse( 1, &m_sRFWE2 ) )  {
                    setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_medium.registerBinomialRFPulse failed.");
                    return ( false );
                }

            } else {
                m_RF_WE11_medium.unregisterBinomialRFPulses();
            }


            m_pRF_WE = &m_RF_WE11_medium;

        break;

        case WE_1_1_90_Deg:

            // * -------------------------------------------------------------- *
            // * Define parameters for sinc pulses used within WE_1_1_90_Deg    *
            // * WE pulse                                                       *
            // * -------------------------------------------------------------- *
            m_sRFWE1.setDuration             (  600 );
            m_sRFWE1.setBandwidthTimeProduct (  6.0 );
            m_sRFWE1.setSamples              (  300 );
            m_sRFWE1.setInitialPhase         (  0.0 );
            m_sRFWE1.setThickness            ( pMrProt->sliceSeries().aFront().thickness() );

            if(! m_RF_WE11_90.registerBinomialRFPulse( 0, &m_sRFWE1 ) )  {
                setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_90.registerBinomialRFPulse failed.");
                return ( false );
            }



            m_sRFWE2.setDuration             ( m_sRFWE1.getDuration()             );
            m_sRFWE2.setBandwidthTimeProduct ( m_sRFWE1.getBandwidthTimeProduct() );
            m_sRFWE2.setSamples              ( m_sRFWE1.getSamples()              );
            m_sRFWE2.setInitialPhase         ( m_sRFWE1.getInitialPhase()         );
            m_sRFWE2.setThickness            ( m_sRFWE1.getThickness()            );

            if(! m_RF_WE11_90.registerBinomialRFPulse( 1, &m_sRFWE2 ) )  {
                setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_90.registerBinomialRFPulse failed.");
                return ( false );
            }

            m_pRF_WE = &m_RF_WE11_90;

        break;

        case WE_1_1_180_Deg:
            m_pRF_WE = &m_RF_WE11;
        break;
        
        case WE_1_2_1_90_Deg:
            m_pRF_WE = &m_RF_WE121_short;
        break;

        case WE_1_2_1_180_Deg:
            m_pRF_WE = &m_RF_WE121;
        break;
        
        default:  
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown water excitation mode." ) )  return ( false );
    }



    return ( true );
    
}


bool SBBGREKernel::calculatePrepWEPulse (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{

    static const char *ptModule = "SBBGREKernel::calculatePrepWEPulse";


    // *----------------------------------------------------------------------- *
    // * Prepares the excitation pulse                                          *
    // *----------------------------------------------------------------------- *
    if ( (m_pRF_WE == NULL) || (m_pRF_Exc == NULL) ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }
    
    m_pRF_WE->setUseWEOffsetFrequency     ( m_bUseWEOffsetFrequency, m_dWEOffsetFrequency );
    m_pRF_WE->setBandwidthTimeProduct     ( m_dWEBandwidthTimeProduct );
    m_pRF_WE->setNonSelectiveWERFDuration ( m_pRF_Exc->getDuration()  );
    m_pRF_WE->setLastRampDownOutsideSBB   (  true );
    m_pRF_WE->setUseOwnEventBlock         ( false );    // * Water selective pulse was used without its own event block *
    m_pRF_WE->setThickness                ( m_pRF_Exc->getThickness() );    // * Slice thickness of the WE pulse        *



    if (! m_pRF_WE->setPointerToCalculationLimits ( m_pCalcLimits ) )   {
        setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->setPointerToCalculationLimits");
        return ( false );
    }
  
    // hughtin6: we set the values in the flip angle array with these memember variables which will have been set by the sequence
    // So first test the pointer to see if it is NULL i.e. it has not been set in the sequence.  If this is the case then we need 
    // to ignore this
    if(m_padFlipAngleArray != NULL)
    {
        if (! m_pRF_WE->setFlipAngleArray ( m_padFlipAngleArray, m_uiNmbrOfFlipAngles ) )  {
            setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->setFlipAngleArray");
            return ( false );
        }
    }
    

      
  
    if (! m_pRF_WE->setExcitationRFPulse(m_pRF_Exc, pMrProt, pSeqExpo) )  {
        setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->setExcitationRFPulse");
        return ( false );
    }

    

    if (! m_pRF_WE->prep ( pMrProt, pSeqLim, pSeqExpo ) )   {
        setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->prep");
        return ( false );
    }




    m_lExcitationAllTime = m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getDurationPerRequest();
    m_lTEFromExcitation  = m_pRF_WE->getTEContribution();



    return ( true );

}






// * -------------------------------------------------------------------------- *
// * Name        : calculate3DGradientTables                                    *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the 3D gradient table and its rewinder        *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// * -------------------------------------------------------------------------- *
bool SBBGREKernel::calculate3DGradientTables (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
  
    static const char *ptModule = "SBBGREKernel::calculate3DGradientTables";


    long   lStatus             = SEQU__NORMAL;

    long   lRamp               = 0;
    long   lDuration           = 0;
  



    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }


    // * ---------------------------------------------------------------------- *
    // * Calculate gradient timing for the slice selection gradient und the     *
    // * phase encode rewinder, kz                                              *
    // * ---------------------------------------------------------------------- *
    if ( m_bFlowCompSliceSelect )   {

        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated slice refocussing gradients               *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->txSpec().excitation() == SEQ::EXCITATION_SLICE_SELECTIVE  )  
        {
        if ( ! m_pRF_WE->calculateShortestFlowCompGradients(pMrProt, pSeqLim, &m_sP_3D_FC, &m_sP_3D, 0.0, SEQ::AXIS_SLICE, m_dFactorForSliceThickness) )  {
            setNLSStatus (m_pRF_WE->getNLSStatus(),
                          ptModule,
                          "m_pRF_WE->calculateShortestFlowCompGradients(...) for slice selection",
                          pSeqLim->isContextPrepForBinarySearch());
            return ( false );
        }
        }
        else
        {
            // * -------------------------------------------------------------- *
            // * Slice refocussing gradient m_sP_3D is not required.            *
            // * -------------------------------------------------------------- *
            m_sP_3D.setAmplitude    (0.0);
            m_sP_3D.setRampUpTime   (0);
            m_sP_3D.setDuration     (0);
            m_sP_3D.setRampDownTime (0);
            // * -------------------------------------------------------------- *
            // * Flow compensation gradient m_sP_3D_FC is not required.         *
            // * -------------------------------------------------------------- *
            m_sP_3D_FC.setAmplitude    (0.0);
            m_sP_3D_FC.setRampUpTime   (0);
            m_sP_3D_FC.setDuration     (0);
            m_sP_3D_FC.setRampDownTime (0);        
        }
      
        m_dMomentum_TB_3D = m_sP_3D_FC.getMomentumTOT() + m_sP_3D.getMomentumTOT();


        if ( m_bFlowCompPartitionEncode )  {

            // * -------------------------------------------------------------- *
            // * Enable flow compensated phase encoding table                   *
            // * -------------------------------------------------------------- *
            m_TB_3D.setbFlowComp (true);

        } else {

            // * -------------------------------------------------------------- *
            // * Disable flow compensated phase encoding table                  *
            // * -------------------------------------------------------------- *
            m_TB_3D.setbFlowComp (false);

        }

        // * ------------------------------------------------------------------ *
        // * Initialize flow compensated phase encoding table                   *
        // * ------------------------------------------------------------------ *
        m_TB_3D.seteGradientAxis          (SBBPhaseEncode::Slice     );
        m_TB_3D.setlTimeToEcho            (m_aRO[0].getlTimeToEcho() );
        m_TB_3D.seteTableDirction         (SEQ::DIR_ASCENDING        );
        m_TB_3D.setdFactorForPixelSize    (m_dFactorForPixelSize3D   );
        m_TB_3D.setdFactorForOffsetMoment (m_dFactorForSliceThickness);
        m_TB_3D.setdOffsetMoment          (0.0                       );
        m_TB_3D.setlFirstLinParToMeasure  (m_lFirstPartToMeasure     );
        m_TB_3D.setlLastLinParToMeasure   (m_lLastPartToMeasure      );


        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated phase encoding table                      *
        // * ------------------------------------------------------------------ *
        if (! m_TB_3D.prep( pMrProt, pSeqLim, pSeqExpo) ) {
            setNLSStatus (m_TB_3D.getNLSStatus(),
                          ptModule,
                          "m_TB_3D.prep(...) failed",
                          pSeqLim->isContextPrepForBinarySearch());
            return ( false );
        }

    } else {

        // * ------------------------------------------------------------------ *
        // * Flow compensation gradient m_sP_3D_FC is not required.             *
        // * ------------------------------------------------------------------ *
        m_sP_3D_FC.setAmplitude    (0.0);
        m_sP_3D_FC.setRampUpTime   (0);
        m_sP_3D_FC.setDuration     (0);
        m_sP_3D_FC.setRampDownTime (0);


        if ( m_bFlowCompPartitionEncode )  {

            // * -------------------------------------------------------------- *
            // * Prepare slice refocussing gradient                             *
            // * -------------------------------------------------------------- *
            m_dMomentum_TB_3D = m_pRF_WE->getGSData().getdRequiredRefocusingMoment();

            // * Calculates the timing and amplitude for a moment m_dMomentum_TB_3D * m_dFactorForSliceThickness *
            m_sP_3D.prepSymmetricTOTShortestTime(m_dMomentum_TB_3D * m_dFactorForSliceThickness);
            // * Keeps the timing but sets the amplitude for a moment m_dMomentum_TB_3D *
            m_sP_3D.prepMomentumTOT (m_dMomentum_TB_3D);


            m_dMomentum_TB_3D = m_sP_3D_FC.getMomentumTOT() + m_sP_3D.getMomentumTOT();


            // * -------------------------------------------------------------- *
            // * Initialize flow compensated phase encoding table               *
            // * -------------------------------------------------------------- *
            m_TB_3D.setbFlowComp     (true);
            m_TB_3D.setdOffsetMoment (0.0 );

        } else {

            // * -------------------------------------------------------------- *
            // * Slice refocussing gradient m_sP_3D is not required.            *
            // * -------------------------------------------------------------- *
            m_sP_3D.setAmplitude    (0.0);
            m_sP_3D.setRampUpTime   (0);
            m_sP_3D.setDuration     (0);
            m_sP_3D.setRampDownTime (0);


            m_dMomentum_TB_3D = m_pRF_WE->getGSData().getdRequiredRefocusingMoment();


            // * -------------------------------------------------------------- *
            // * Initialize none-flow compensated phase encoding table          *
            // * -------------------------------------------------------------- *
            m_TB_3D.setbFlowComp     (false            );
            m_TB_3D.setdOffsetMoment (m_dMomentum_TB_3D);

        }

        // * ------------------------------------------------------------------ *
        // * Initialize flow compensated phase encoding table                   *
        // * ------------------------------------------------------------------ *
        m_TB_3D.seteGradientAxis          (SBBPhaseEncode::Slice     );
        m_TB_3D.seteTableDirction         (SEQ::DIR_ASCENDING        );
        m_TB_3D.setdFactorForPixelSize    (m_dFactorForPixelSize3D   );
        m_TB_3D.setdFactorForOffsetMoment (m_dFactorForSliceThickness);
        m_TB_3D.setlFirstLinParToMeasure  (m_lFirstPartToMeasure     );
        m_TB_3D.setlLastLinParToMeasure   (m_lLastPartToMeasure      );


        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated phase encoding table                      *
        // * ------------------------------------------------------------------ *
        if (! m_TB_3D.prep( pMrProt, pSeqLim, pSeqExpo) ) {
            setNLSStatus (m_TB_3D.getNLSStatus(),
                          ptModule,
                          "m_TB_3D.prep(...) failed",
                          pSeqLim->isContextPrepForBinarySearch());
            return ( false );
        }

    }

    

    // * ---------------------------------------------------------------------- *
    // * Calculate gradient timing for phase encode rewinder, kz                *
    // * ---------------------------------------------------------------------- *
    if ( m_bUsePERewinder == true ) {
  
        if ( m_bBalancedGradMomentSliceSelect )  {
            m_dMomentum_TB_3DR = m_dMomentum_TB_3D * m_dFactorForSliceThickness;  // * Balanced gradient moments for TrueFisp *
        } else  {
            m_dMomentum_TB_3DR =  -1.0 * m_dMomentum_TB_3D * m_dRelSSDephasingMoment * m_dFactorForSliceThickness;
        }

        double dDeltaMomentumKz = fGSLGet3DDeltaMoment( pMrProt );

        lStatus = fGSLGetShortestTabTiming (
                    dDeltaMomentumKz * m_dFactorForPixelSize3D,      // IMP: Moment between two steps [mT/m * us]
                    SEQ::DIR_DESCENDING,                             // IMP: true <=> (-) -> (+), i.e. table direction
                    m_dMomentum_TB_3DR,                              // IMP: the offset MOMENT of the table [us(mT/m)]
                    m_lFirstPartToMeasure,                           // IMP: number of first line/partition (note: Echo=0)
                    m_lLastPartToMeasure,                            // IMP: number of last line/partition (note: Echo=0)
                    m_sTB_3DR.getMaxMagnitude(),                     // IMP: maximum allowed gradient amplitude [mT/m]
                    m_sTB_3DR.getMinRiseTime(),                      // IMP: minimum allowed rise time [us/(mT/m)]
                   &lRamp,                                           // EXP: minimum required ramp time [us] (up/down)
                   &lDuration                                        // EXP: minimum required duration  [us]
                );

        if ( setNLSStatus (lStatus, ptModule, "fGSLGetShortestTabTiming for 3D", pSeqLim->isContextPrepForBinarySearch()) )  {
            return ( false );
        }

    } else {
        lRamp     = 0L;
        lDuration = 0L;
    }

    m_sTB_3DR.setRampUpTime   ( lRamp     );    // Set gradient timing
    m_sTB_3DR.setRampDownTime ( lRamp     );
    m_sTB_3DR.setDuration     ( lDuration );


    return ( true );

}





// * -------------------------------------------------------------------------- *
// *                                                                            *
// * Name        :  SBBGREKernel::calculatePEGradientTables                     *
// *                                                                            *
// * Description :  Calculation of the PE gradient table and its rewinder       *
// *                The phase encoding gradient is flow compensated when the    *
// *                member variable m_bFlowCompPhase is set to true.            *
// *                                                                            *
// * Parameter   :  pMrProt:                                                    *
// *                    Pointer to the protocol structure                       *
// *                                                                            *
// * Return      :  bool                                                        *
// *                                                                            *
// * -------------------------------------------------------------------------- *
bool SBBGREKernel::calculatePEGradientTables (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
  
    static const char *ptModule = "SBBGREKernel::calculatePEGradientTables";

    long   lStatus             = SEQU__NORMAL;

    long   lRamp               = 0;
    long   lDuration           = 0;
  

    // * ---------------------------------------------------------------------- *
    // * Initialize phase encoding table                                        *
    // * ---------------------------------------------------------------------- *
    m_TB_PE.seteGradientAxis         (SBBPhaseEncode::Phase    );
    m_TB_PE.setlTimeToEcho           (m_aRO[0].getlTimeToEcho());
    m_TB_PE.seteTableDirction        (SEQ::DIR_ASCENDING       );
    m_TB_PE.setdOffsetMoment         (0.0                      );
    m_TB_PE.setdFactorForPixelSize   (m_dFactorForPixelSizePE  );
    m_TB_PE.setlFirstLinParToMeasure (m_lFirstLineToMeasure    );
    m_TB_PE.setlLastLinParToMeasure  (m_lLastLineToMeasure     );


    // * ---------------------------------------------------------------------- *
    // * Calculate the phase encoding table                                     *
    // * ---------------------------------------------------------------------- *
    if ( m_bFlowCompPhase )  {

        // * ------------------------------------------------------------------ *
        // * Enable flow compensation in phase encoding direction               *
        // * ------------------------------------------------------------------ *
        m_TB_PE.setbFlowComp (true);


        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated phase encoding table                      *
        // * ------------------------------------------------------------------ *
        if (! m_TB_PE.prep( pMrProt, pSeqLim, pSeqExpo) ) {
            setNLSStatus (m_TB_PE.getNLSStatus(),
                          ptModule,
                          "m_TB_PE.prep(...) failed",
                          pSeqLim->isContextNormal());
            return ( false );
        }


        // * ------------------------------------------------------------------ *
        // * Calculation of the phase encoding rewinder gradient                *
        // * ------------------------------------------------------------------ *
        if ( m_bUsePERewinder == true) {

            lStatus = fGSLGetShortestTabTiming (
                        fGSLGetPEDeltaMoment(pMrProt) * m_dFactorForPixelSizePE,    // IMP: Moment between two steps [mT/m * us]
                        SEQ::DIR_DESCENDING,                                        // IMP: true <=> (-) -> (+), i.e. table direction
                        0.0,                                                        // IMP: the offset MOMENT of the table [us(mT/m)]
                        m_lFirstLineToMeasure,                                      // IMP: number of first line/partition (note: Echo=0)
                        m_lLastLineToMeasure,                                       // IMP: number of last line/partition (note: Echo=0)
                        m_sTB_PER.getMaxMagnitude(),                                // IMP: maximum allowed gradient amplitude [mT/m]
                        m_sTB_PER.getMinRiseTime(),                                 // IMP: minimum allowed rise time [us/(mT/m)]
                       &lRamp,                                                      // EXP: minimum required ramp time [us] (up/down)
                       &lDuration                                                   // EXP: minimum required duration  [us]
                      );

            if ( setNLSStatus (lStatus, ptModule, "fGSLGetShortestTabTiming for PE") ) return ( false );

            m_sTB_PER.setRampUpTime   ( lRamp     );
            m_sTB_PER.setRampDownTime ( lRamp     );
            m_sTB_PER.setDuration     ( lDuration );

        } else {

            m_sTB_PER.setRampUpTime   ( 0 );
            m_sTB_PER.setRampDownTime ( 0 );
            m_sTB_PER.setDuration     ( 0 );

        }

    } else {

        // * ------------------------------------------------------------------ *
        // * Disable flow compensation in phase encoding direction              *
        // * ------------------------------------------------------------------ *
        m_TB_PE.setbFlowComp (false);


        // * ------------------------------------------------------------------ *
        // * Prepare flow compensated phase encoding table                      *
        // * ------------------------------------------------------------------ *
        if (! m_TB_PE.prep( pMrProt, pSeqLim, pSeqExpo) ) {
            setNLSStatus (m_TB_PE.getNLSStatus(),
                          ptModule,
                          "m_TB_PE.prep(...) failed",
                          pSeqLim->isContextNormal());
            return ( false );
        }


        // * ------------------------------------------------------------------ *
        // * Set the gradient timing for phase encode rewinding                 *
        // * ------------------------------------------------------------------ *
        if ( m_bUsePERewinder == true) {
            m_sTB_PER.setRampUpTime   ( m_TB_PE.getlRampDuration() );
            m_sTB_PER.setRampDownTime ( m_TB_PE.getlRampDuration() );
            m_sTB_PER.setDuration     ( m_TB_PE.getDurationPerRequest() - m_TB_PE.getlRampDuration() );
        } else {
            m_sTB_PER.setRampUpTime   ( 0 );
            m_sTB_PER.setRampDownTime ( 0 );
            m_sTB_PER.setDuration     ( 0 );
        }

    }



    return ( true );

}








bool SBBGREKernel::calculatePrepRO (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{

    long lI = 0;


    // * ---------------------------------------------------------------------- *
    // * Specify settings for the various contrasts                             *
    // * ---------------------------------------------------------------------- *
    for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {

        // * ------------------------------------------------------------------ *
        // * Contrast ID, ranging from 0 to m_lNumberOfEchoes-1                 *
        // * ------------------------------------------------------------------ *
        m_aRO[lI].setlContrast ( lI );


        // * ------------------------------------------------------------------ *
        // * Readout polarity                                                   *
        // *    Monopolar  : Positive for all contrasts                         *
        // *    Bipolar    : Alternating positive - negative                    *
        // * ------------------------------------------------------------------ *
        switch ( m_eROPolarity )  {
            case Bipolar:
                if ( lI == 0 )  {
                    m_aRO[lI].seteROPolarity (SeqBuildBlockReadOut::Positive);
                } else {
                    if ( m_aRO[lI-1].geteROPolarity() == SeqBuildBlockReadOut::Positive )  {
                        m_aRO[lI].seteROPolarity (SeqBuildBlockReadOut::Negative);
                    } else {
                        m_aRO[lI].seteROPolarity (SeqBuildBlockReadOut::Positive);
                    }
                }
            break;

            case Monopolar:
                m_aRO[lI].seteROPolarity (SeqBuildBlockReadOut::Positive);
            break;

            default:
                setNLSStatus ( SEQU_ERROR, ptModule, "Unkown readout polarity." );
                return ( false );
        }

        
        // * ------------------------------------------------------------------ *
        // * Gradient ramps                                                     *
        // * ------------------------------------------------------------------ *
        if ( lI == 0 )  
        {
             if (m_bUseRampSampling)
             {
                m_aRO[lI].setbRampUpInSBB (true);
             }
             else
             {
            m_aRO[lI].setbRampUpInSBB (false);
             }
        } else {
            m_aRO[lI].setbRampUpInSBB (true );
        }

        if ( lI == m_lNumberOfEchoes-1 )  {
            m_aRO[lI].setbRampDownInSBB (false);
        } else {
            m_aRO[lI].setbRampDownInSBB (true);
        }


        // * ------------------------------------------------------------------ *
        // * Readout (a)symmetry                                                *
        // * ------------------------------------------------------------------ *
        if ( m_eReadOutGradType == Symmetric_Rewinder ) {
            m_aRO[lI].setbSymmetricTiming ( true );
        }
        m_aRO[lI].setdReadOutAsym ( m_adROAsymmetryBefore[lI], m_adROAsymmetryAfter[lI] );

        // * ------------------------------------------------------------------ *
        // * Ramp sampling measurement                                          *
        // * ------------------------------------------------------------------ *
        m_aRO[lI].setbUseRampSampling(m_bUseRampSampling);

        // * ------------------------------------------------------------------ *
        // * Calculation limit                                                  *
        // * ------------------------------------------------------------------ *
        m_aRO[lI].setPointerToCalculationLimits (m_pCalcLimits);


        // * ------------------------------------------------------------------ *
        // * Positioning of the echo time on the gradient raster                *
        // * ------------------------------------------------------------------ *
        m_aRO[lI].setbEchoOnGRT (m_bEchoOnGRT);


        // * ------------------------------------------------------------------ *
        // * Data from previous contrasts                                       *
        // * ------------------------------------------------------------------ *
        if ( lI > 0 )  {
            m_aRO[lI].setdROMomentIn                (m_aRO[lI-1].getdROMomentOut() );
            m_aRO[lI].setlTimeBetweenADCsPreviousRO (m_aRO[lI-1].getlTimeAfterADC());
        }



        // * ------------------------------------------------------------------ *
        // * Preparation of SBBReadOut                                          *
        // * ------------------------------------------------------------------ *
        if ( ! m_aRO[lI].prep (pMrProt, pSeqLim, pSeqExpo) )  {

            if ( pSeqLim->isContextNormal() )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Preparation of readout gradient failed.\n");
            } else {
                setNLSStatus(SEQU_ERROR);
            }
            return ( false );

        }

        m_lReadoutTime += m_aRO[lI].getDurationPerRequest();
    }


    return ( true );

}






// ******************************************************************************
// * Name        : calculateRORampDownTime                                      *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the read out ramp down times                  *
// *                                                                            *
// * Return      : BOOL                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateRORampDownTime (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo *pSeqExpo)
{

    // * ---------------------------------------------------------------------- *
    // * The ramp down time of the last readout has to be fitted to ramp up     *
    // * time of the spoiler if phase stabilization is disabled                 *
    // * ---------------------------------------------------------------------- *
    if ( ! pMrProt->phaseStabilize() )  {
        m_aRO[m_lNumberOfEchoes-1].setlRequestedRampDownTime( m_sP_ROD.getRampUpTime() );    
        m_aRO[m_lNumberOfEchoes-1].prep (pMrProt, pSeqLim, pSeqExpo);
    }    

  
    return ( true );
  
}



// ******************************************************************************
// * Name        : calculateROP                                                 *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the read out prephasing gradient              *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateROP (MrProt* pMrProt, SeqLim* pSeqLim)
{
  

    static const char *ptModule = "SBBGREKernel::calculateROP";
  
  
    if ( m_bFlowCompRead )  {

        // * ------------------------------------------------------------------ *
        // * Calculate flow compensation gradients                              *
        // * ------------------------------------------------------------------ *
        MeasNucleus MeasNucleus(pMrProt->txSpec().nucleusInfoArray()[0].nucleus());
        

        sGRAD_PULSE_TRAP cGradTemp;
        cGradTemp.setMaxMagnitude ( m_aRO[0].getdROMaxMagnitude() );
        cGradTemp.setRampUpTime   ( m_aRO[0].getlRampUpTime()     );
        cGradTemp.setAmplitude    ( m_aRO[0].getdROAmplitude()    );
        cGradTemp.setDuration     ( cGradTemp.getRampUpTime() + m_aRO[0].getlROFlatToEcho() );

        if (! fSUCalcShortestFlowComp ( 
                                        MeasNucleus,
                                        0.0,
                                       &cGradTemp,
                                        1.0,
                                        m_dFactorForPixelSizeRO,
                                       &m_sP_ROP_FC,
                                       &m_sP_ROP
                                      )
           )
        {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                if ( setNLSStatus(SEQU_ERROR, ptModule, "Error in calculation of flow compensation (read out). \n") ) return ( false );
            } else  {
                if ( setNLSStatus(SEQU_ERROR, ptModule) ) return ( false );
            }
        }

        m_dMomentROP = m_sP_ROP.getMomentumTOT();

    } else {
        
        // * ------------------------------------------------------------------ *
        // * Flow compensation gradient sP_ROP_FC is not required.              *
        // * ------------------------------------------------------------------ *
        m_sP_ROP_FC.setAmplitude    (0.0);
        m_sP_ROP_FC.setRampUpTime   (0);
        m_sP_ROP_FC.setDuration     (0);
        m_sP_ROP_FC.setRampDownTime (0);



        // * ------------------------------------------------------------------ *
        // * Moment of the prephasing gradient in read out direction            *
        // * considering asymmetry of k-space sampling                          *
        // * ------------------------------------------------------------------ *
         m_dMomentROP = (-1.0) * m_aRO[0].getdROMomentToKSpaceCenter();



        // * ------------------------------------------------------------------ *
        // * Calcution of the shortest possible prephasing gradient in          *
        // * readout direction                                                  *
        // * ------------------------------------------------------------------ *

        // * Calculates the timing and amplitude for a moment m_dMomentROP * m_dFactorForPixelSizeRO *
        m_sP_ROP.prepSymmetricTOTShortestTime(m_dMomentROP * m_dFactorForPixelSizeRO);
        // * Keeps the timing but sets the amplitude for a moment m_dMomentROP *
        m_sP_ROP.prepMomentumTOT (m_dMomentROP);

    }


    return ( true );

}


// ******************************************************************************
// * Name        : calculateROD                                                 *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the read out dephasing gradient               *
// *                                                                            *
// * Return      : BOOL                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateROD (MrProt* pMrProt)
{
  
    static const char *ptModule = "SBBGREKernel::calculateROD";
  
    long   lI;                         // * Loop counter *
    double dSign               = 1.0;
    double dMomentToRewind     = 0.0;

  
  
    // * ---------------------------------------------------------------------- *
    // * Sign of the read out dephasing gradient                                *
    // * ---------------------------------------------------------------------- *
    lI = m_lNumberOfEchoes - 1;                     // * Index of the last echo *

    dSign = m_aRO[lI].geteROPolarity();

    if ( pMrProt->phaseStabilize() ) {  // * One more ADC in case of phase stabilization *
        dSign *= -1.0;
    } 



    // * Moment of the readout gradient without asymmetric readout *
    double dMomentRO = fGSLGetROMoment ( pMrProt );


    // * ---------------------------------------------------------------------- *
    // * Duration, amplitude and ramp times of the read out dephasing gradient  *
    // * ---------------------------------------------------------------------- *
    switch ( m_eReadOutGradType )  {

        case Spoiler:

            m_sP_ROD.setAmplitude    ( dSign /* * m_eROPolarity*/ * m_sP_ROD.getMaxMagnitude() );
            if ( pMrProt->phaseStabilize() ) {
                m_sP_ROD.setRampUpTime ( m_PhaseStabScan.getRORampDownTime() );
            } else {

                double dMinROAmplitude = 0.0;
                double dMaxROAmplitude = 0.0;

                m_aRO[lI].getROAmplitudeRange (&dMinROAmplitude , &dMaxROAmplitude );

                m_sP_ROD.setRampUpTime ( RoundUpGRT( fabs(dSign * m_sP_ROD.getMaxMagnitude() - dMinROAmplitude) 
                                                     * m_sP_ROD.getMinRiseTime() ) );

            }

            if ( fabs(m_sP_ROD.getAmplitude()) < m_dEpsilon  )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Readout amplitude is zero. Cannot proceed.\n");
                return ( false );
            } else {
                // The momentum of read out dephasing gradient is m_dRelRODephasingMoment of the momentum of the read out gradient *
                m_sP_ROD.setDuration ( m_sP_ROD.getRampUpTime() + RoundUpGRT( fabs(dMomentRO * m_dFactorForPixelSizeRO * m_dRelRODephasingMoment / m_sP_ROD.getAmplitude() ) ) );
            }

            m_sP_ROD.setRampDownTime ( RoundUpGRT( fabs(m_sP_ROD.getAmplitude()) * m_sP_ROD.getMinRiseTime() ) );

        break;

        case Constant:
            if ( pMrProt->phaseStabilize() ) {

                m_sP_ROD.setAmplitude    ( m_PhaseStabScan.getROAmplitude()    );
                m_sP_ROD.setRampUpTime   ( m_PhaseStabScan.getRORampDownTime() );

                if ( fabs(m_PhaseStabScan.getROAmplitude()) < m_dEpsilon  )  {
                    setNLSStatus(SEQU_ERROR, ptModule, "Phase stabilization readout amplitude is zero. Cannot proceed.\n");
                    return ( false );
                } else {
                    m_sP_ROD.setDuration ( RoundUpGRT( fabs(dMomentRO * m_dFactorForPixelSizeRO * m_dRelRODephasingMoment / m_PhaseStabScan.getROAmplitude() ) ) );
                }


            } else {

                m_sP_ROD.setAmplitude    ( m_aRO[lI].getdROAmplitude()  );
                m_sP_ROD.setRampUpTime   ( m_aRO[lI].getlRampDownTime() );

                // * The gradient was switched on as long as the ADC was switched on.   *
                // * The read gradient is ramped down while the read out dephasing      *
                // * was ramped up at the same time, resulting in a constant gradient   *
                // * amplitude.                                                         *

                m_sP_ROD.setDuration     ( RoundUpGRT( (double)( m_aRO[lI].getlRampDownTime() /*m_sP_RO[lI].getRampDownTime()*/ /*+ m_lADCZeroAfter[lI]*/)) );
            }

            m_sP_ROD.setRampDownTime ( RoundUpGRT(fabs(m_sP_ROD.getAmplitude() * m_dFactorForPixelSizeRO * m_sP_ROD.getMinRiseTime()) ) );
        break;

        case Symmetric_Rewinder:
            m_sP_ROD.setRampUpTime   ( m_sP_ROP.getRampUpTime()   );
            m_sP_ROD.setDuration     ( m_sP_ROP.getDuration()     );
            m_sP_ROD.setRampDownTime ( m_sP_ROP.getRampDownTime() );
            m_sP_ROD.setAmplitude    ( m_sP_ROP.getAmplitude()    );
        break;

        case Rewinder:
            dMomentToRewind = m_aRO[m_lNumberOfEchoes-1].getdROMomentOut();

            // * Sets the amplitude, duration and ramp times to realize the gradient moment -1.0 * dMomentumToRewind *
            m_sP_ROD.prepSymmetricTOTShortestTime( -1.0 * dMomentToRewind * m_dFactorForPixelSizeRO);
            m_sP_ROD.prepMomentumTOT             ( -1.0 * dMomentToRewind                          );
        break;
    
        default: 
            setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." );
            return ( false );

    }
  
    return ( true );
  
}




// ******************************************************************************
// * Name        : calculatePrephaseTime                                        *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the prephase duration                         *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculatePrephaseTime (MrProt* , SeqLim* , SeqExpo*)
{
  
    static const char *ptModule = "SBBGREKernel::calculatePrephaseTime";

    long lTimeSlice   = 0;
    long lTimePhase   = 0;
    long lTimeRead    = 0;

  
  
    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }


    // * ---------------------------------------------------------------------- *
    // * Calculate the maximum prephasing time on all gradient axes             *
    // * ---------------------------------------------------------------------- *
    lTimeSlice =  m_pRF_WE->getRampTimeOutsideSBB();
    lTimeSlice += m_sP_3D.getTotalTime() + m_sP_3D_FC.getTotalTime() + m_TB_3D.getDurationPerRequest();

    lTimePhase =  m_TB_PE.getDurationPerRequest();

    lTimeRead  =  m_sP_ROP_FC.getTotalTime() + m_sP_ROP.getTotalTime() + m_aRO[0].getlRampUpTime();
    if (m_bUseRampSampling)
        lTimeRead  -= m_aRO[0].getlRampUpTime();


    m_lPrephaseTime = max (lTimeSlice     , lTimePhase);
    m_lPrephaseTime = max (m_lPrephaseTime, lTimeRead );
    m_lPrephaseTime = max (m_lPrephaseTime, m_lMinDelayBetweenRFAndADC);


    // * ---------------------------------------------------------------------- *
    // * Which gradient axis limits the prephasing time ?                       *
    // * ---------------------------------------------------------------------- *
    if ( lTimeSlice == m_lPrephaseTime ) {
        m_eLimitingPrephaseTime = Slice;
    } else if ( lTimePhase == m_lPrephaseTime ) {
        m_eLimitingPrephaseTime = Phase;
    } else if ( lTimeRead == m_lPrephaseTime )  {
        m_eLimitingPrephaseTime = Read;
    } else {
        m_eLimitingPrephaseTime = None;
    }


    // * ---------------------------------------------------------------------- *
    // * m_lADCZeroBefore[0] might not be a multiple of the gradient raster     *
    // * time -> m_lPrephaseTime[0] has be rounded to the gradient raster       *
    // * time.                                                                  *
    // * ---------------------------------------------------------------------- *
    m_lPrephaseTime = RoundUpGRT ( m_lPrephaseTime );



    // * ---------------------------------------------------------------------- *
    // * All prephasing gradients are stretched to fill the maximum             *
    // * prephasing time.                                                       *
    // * ---------------------------------------------------------------------- *
    m_lPrephaseFillTimeSlice = m_lPrephaseFillTimePhase = m_lPrephaseFillTimeRead = 0;

    if ( lTimeSlice < m_lPrephaseTime )  {
        if ( m_bFlowCompRead || m_bFlowCompSliceSelect || m_bFlowCompPartitionEncode || m_bFlowCompPhase ) {
            m_lPrephaseFillTimeSlice = m_lPrephaseTime - lTimeSlice;
        } else {
            m_TB_3D.resize (m_lPrephaseTime - m_TB_3D.getlRampDuration() - m_pRF_WE->getRampTimeOutsideSBB());
        }

    }

    if ( lTimePhase < m_lPrephaseTime )  {
        if ( m_bFlowCompRead || m_bFlowCompSliceSelect || m_bFlowCompPartitionEncode || m_bFlowCompPhase )  {
            m_lPrephaseFillTimePhase = m_lPrephaseTime - lTimePhase;
        } else {
            m_TB_PE.resize (m_lPrephaseTime - m_TB_PE.getlRampDuration());
        }
        
    }

    if ( lTimeRead < m_lPrephaseTime )  {
        if ( m_bFlowCompRead || m_bFlowCompSliceSelect || m_bFlowCompPartitionEncode || m_bFlowCompPhase ) {

            m_lPrephaseFillTimeRead = m_lPrephaseTime - lTimeRead;

        } else {
            // * -------------------------------------------------------------- *
            // * No flow compensation in read out direction                     *
            // * Use all the prephase time for m_sP_ROP. m_sP_ROP_FC is zero.   *
            // * -------------------------------------------------------------- *
            m_sP_ROP.setDuration ( m_sP_ROP.getDuration() + m_lPrephaseTime - lTimeRead );
        }
        
    }



    // * ---------------------------------------------------------------------- *
    // * Final calculatation of the relative prephasing gradient amplitude      *
    // * in read out direction                                                  *
    // * No recalculation of m_sP_ROP in case of flow compensation:             *
    // *    - flow compensation in read out direction: the timing must not be   *
    // *      changed to preserve flow                                          *
    // *    - slice compensation in slice direction: the read out prephasing    *
    // *      gradient should be stretched to reduce flow effects               *
    // * compensation.                                                          *
    // * ---------------------------------------------------------------------- *
    if ( ! (m_bFlowCompRead || m_bFlowCompSliceSelect || m_bFlowCompPartitionEncode || m_bFlowCompPhase) ) {
        if ( fabs(static_cast<double>(m_sP_ROP.getDuration())) < m_dEpsilon  )  {
            if (!m_bUseRampSampling) 
            {
                setNLSStatus(SEQU_ERROR, ptModule, "Readout prephasing duration is zero. Cannot proceed.\n");
                return ( false );
            }
        } else {
            m_sP_ROP.setAmplitude ( m_dMomentROP / static_cast<double>(m_sP_ROP.getDuration()) );
        }
    }



    return ( true );
  
}



// ******************************************************************************
// * Name        : calculateRephaseTime                                         *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the rephase duration                          *
// *               All rephasing gradients are stretched to fit into the        *
// *               rephasing interval                                           *
// *                                                                            *
// * Return      : void                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateRephaseTime (MrProt* pMrProt)
{

    static const char *ptModule = "SBBGREKernel::calculateRephaseTime";

    long  lTimeSlice    = 0;
    long  lTimePhase    = 0;
    long  lTimeRead     = 0;

    long  lOldGradTime  = 0;     // * Temp. variable     *
    long  lNewGradTime  = 0;



    if ( m_pRF_WE == NULL )  {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE = NULL.\n");
        return ( false );
    }



    // * -------------------------------------------------------------------------- *
    // * Calculate the maximum rephasing time on the slice select, phase encode     *
    // * or read out gradient axis                                                  *
    // * In case of phase stabilization ROD is carried out AFTER the phase          *
    // * stabilization scan and DOES NOT contribute to rephasing time.              *
    // * -------------------------------------------------------------------------- *

    lTimeSlice = m_sTB_3DR.getTotalTime();
    lTimePhase = m_sTB_PER.getTotalTime();
    lTimeRead  = m_sP_ROD.getTotalTime();


    if ( m_eReadOutGradType == Symmetric_Rewinder )  {
         lTimeSlice += m_pRF_WE->getRampTimeOutsideSBB();
    }

    if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) )  {   // * The read out and read out dephasing gradients do         *
         lTimeRead += m_aRO[m_lNumberOfEchoes-1].getlRampDownTime();                           // * NOT overlap in case of m_eReadOutGradType equal Rewinder *
    }

    m_lRephaseTime = maximum(lTimeSlice, lTimePhase);
    if ( ! pMrProt->phaseStabilize() )  {
        m_lRephaseTime = maximum(m_lRephaseTime, lTimeRead );
    }

    m_lRephaseTime = RoundUpGRT ( m_lRephaseTime );


    // * ---------------------------------------------------------------------- *
    // * The rephasing gradients in slice select and phase encode               *
    // * direction and the dephasing gradient in read out direction are         *
    // * stretched to fill the maximum rephasing time.                          *
    // *                                                                        *
    // * If this operation would be carried out in case of an gradient duration *
    // * that is equal to 0, dangling gradient ramps would be impossible        *
    // * ---------------------------------------------------------------------- *

    if (lTimeSlice < m_lRephaseTime) {
        if ( m_sTB_3DR.getRampDownTime() != 0 ) {
            m_sTB_3DR.setDuration( m_sTB_3DR.getDuration() + m_lRephaseTime - lTimeSlice );
        }
    }

    if (lTimePhase < m_lRephaseTime) {
        if ( m_sTB_PER.getRampDownTime() != 0 ) {
            m_sTB_PER.setDuration( m_sTB_PER.getDuration() + m_lRephaseTime - lTimePhase );
        }
    }

    if ( (lTimeRead < m_lRephaseTime) && ( ! pMrProt->phaseStabilize() ) ) {
        if ( (m_eReadOutGradType == Spoiler) || (m_eReadOutGradType == Rewinder) ) {
            // * The amplitude of the read out dephasing gradient can be reduced    *
            // * since its duration is increased:                                   *
            // * NewAmplitude = OldAmplitude * OldGradientTime / NewGradientTime    *


            lOldGradTime = static_cast<long>(m_sP_ROD.getFlatTopTime() + 0.5 * ( m_sP_ROD.getRampUpTime() + m_sP_ROD.getRampDownTime() ));
            lNewGradTime = static_cast<long>(lOldGradTime + m_lRephaseTime - lTimeRead );


            if ( fabs(lNewGradTime) < m_dEpsilon )  {
                setNLSStatus(SEQU_ERROR, ptModule, "Divisor is zero. Cannot proceed.\n");
                return ( false );
            } else {
                m_sP_ROD.setAmplitude( m_sP_ROD.getAmplitude() * lOldGradTime / lNewGradTime );
            }

        }
        m_sP_ROD.setDuration( m_sP_ROD.getDuration() + m_lRephaseTime - lTimeRead );
    }
  
    return ( true );  

}


// ******************************************************************************
// * Name        : calculateOverTime                                            *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the duration that gradient ramps are allowed  *
// *               to stick out of the event block                              *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateOverTime (MrProt* pMrProt, SeqLim*)
{

    static const char *ptModule = "SBBGREKernel::calculateOverTime";

    long              lMinDurationOfRephase = 0;


    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * lMinDurationOfRephase is the minimum duration of the rephase time      *
    // * that ensures that the time between the end of the ADC and the start    *
    // * of the next rf-pulse is greater than lMinDurationBetweenADCAndRF       *
    // * ---------------------------------------------------------------------- *
    if ( m_bRampsOutsideOfEventBlock == true )  {

        if ( pMrProt->phaseStabilize() )  {

            // * The last ADC is the phase stabilization ADC *
            lMinDurationOfRephase = maximum( m_lMinDelayBetweenADCAndRF - m_pRF_WE->getRequiredLeadTime() - m_pRF_WE->getGSData().getlTimeToFirstFlatTop(), 0L );
            lMinDurationOfRephase = RoundUpGRT ( lMinDurationOfRephase );

            m_lOverTime = minimum ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );

            if ( (m_sP_ROD.getTotalTime() - m_lOverTime) < lMinDurationOfRephase ) {
                m_lOverTime -= lMinDurationOfRephase - (m_sP_ROD.getTotalTime() - m_lOverTime);
            }

        } else {

            // * Last kernel ADC *
            lMinDurationOfRephase = maximum( m_lMinDelayBetweenADCAndRF - m_pRF_WE->getRequiredLeadTime() - m_pRF_WE->getGSData().getlTimeToFirstFlatTop(), 0L );
            lMinDurationOfRephase = RoundUpGRT ( lMinDurationOfRephase );

            // * ------------------------------------------------------------------ *
            // * Only the gradient with the longest FlatTopTime (end of the event   *
            // * block) can have a dangling gradient ramp                           *
            // * ------------------------------------------------------------------ *
            if ( m_sTB_3DR.getDuration() > m_sTB_PER.getDuration() )  {

                if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) )  {
      
                    if ( m_sTB_3DR.getDuration() > m_aRO[0].getlRampDownTime() + m_sP_ROD.getDuration() ) {
                        m_lOverTime = minimum ( m_sTB_3DR.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                    } else {
                        m_lOverTime = minimum ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                    }
        
                } else {
      
                    if ( m_sTB_3DR.getDuration() > m_sP_ROD.getDuration() ) {
                        m_lOverTime = minimum ( m_sTB_3DR.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                    } else {
                        m_lOverTime = minimum ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                    }     
        
                }

            } else {

                if ( m_sTB_PER.getDuration() > m_sP_ROD.getDuration() )  {
                    m_lOverTime = minimum ( m_sTB_PER.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                } else {
                    m_lOverTime = minimum ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                }
            }


            if ( (m_lRephaseTime - m_lOverTime) < lMinDurationOfRephase ) {
                m_lOverTime -= lMinDurationOfRephase - (m_lRephaseTime - m_lOverTime);
            }

        }


    } else {
  
        if ( pMrProt->phaseStabilize() )  {

            // * No gradient ramps outside of the event block *
            if ( m_sP_ROD.getTotalTime() < lMinDurationOfRephase ) {
                m_lOverTime -= lMinDurationOfRephase - m_sP_ROD.getTotalTime();
            } else {
                m_lOverTime =  0;
            }

        } else {

            // * No gradient ramps outside of the event block *
            if ( m_lRephaseTime < lMinDurationOfRephase ) {
                m_lOverTime -= lMinDurationOfRephase - m_lRephaseTime;
            } else {
                m_lOverTime =  0;
            }

        }

    }

    
    return ( true );

}


// ******************************************************************************
// * Name        : calculateTEMin                                               *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the minimum TE times                          *
// *                                                                            *
// * Return      : BOOL                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateTEMin (MrProt *pMrProt)
{
  
    long lTEMin = m_lTEFromExcitation + m_lPrephaseTime;

    // * ---------------------------------------------------------------------- *
    // * lTEMin is rounded to the gradient raster time in order to facilitate   *
    // * the evaluation of the TE fill times                                    *
    // * ---------------------------------------------------------------------- *
    for ( long lI = 0; lI < pMrProt->contrasts(); lI++ )  {

        m_alTEMin[lI] = RoundUpGRT ( lTEMin + m_aRO[lI].getlTimeToEcho() );

        lTEMin += m_aRO[lI].getDurationPerRequest();

    }


    return ( true );

}



// ******************************************************************************
// * Name        : calculateTEFill                                              *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the TE fill times                             *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateTEFill (MrProt* pMrProt)
{

    long lI = 0;

    // * ---------------------------------------------------------------------- *
    // * Calculate TE fill times                                                *
    // * ---------------------------------------------------------------------- *
    m_alTEFill[0] = pMrProt->te()[0] - m_alTEMin[0];

    long lTotalTEFill = 0;
    for ( lI = 1; lI < m_lNumberOfEchoes; lI++ )  {
        lTotalTEFill   += m_alTEFill[lI-1];
        m_alTEFill[lI] =  pMrProt->te()[lI] - m_alTEMin[lI] - lTotalTEFill;
    }


    // * ---------------------------------------------------------------------- *
    // * Reset unused TE fill times                                             *
    // * ---------------------------------------------------------------------- *
    for ( lI = m_lNumberOfEchoes; lI < lMaxNumberOfContrasts-1; lI++ )  {
        m_alTEFill[lI] = 0;
    }


    return ( true );

}




// ******************************************************************************
// * Name        : calculateTRMin                                               *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Calculation of the minimum TR time                           *
// *               Acoustic resonances could be avoided.                        *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateTRMin (MrProt* pMrProt)
{


    // * ---------------------------------------------------------------------- *
    // * Gradient ramps that are dangling out of the event block (m_lOverTime)  *
    // * do not contribute to the minimum TR time                               *
    // * ---------------------------------------------------------------------- *
    m_lTRMin = m_lExcitationAllTime  + m_lDelayAdjustSSPosition + m_lPrephaseTime +
               m_lReadoutTime        + m_lRephaseTime           - m_lOverTime ;


    if ( m_bAvoidAcousticResonances ) { 
        for ( long lI = 0; lI < pMrProt->contrasts(); lI++ )  {
            m_lTRMin += maximum (m_alTEFill[lI], 0L);
        }
    } else {
        for ( long lI = 0; lI < pMrProt->contrasts(); lI++ )  {
            m_lTRMin += m_alTEFill[lI];
        }
    }


    // * ---------------------------------------------------------------------- *
    // * TR contributions from phase stabilization                              *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->phaseStabilize() )  {
        m_lTRMin += m_PhaseStabScan.getDurationPerRequest() + m_sP_ROD.getTotalTime();
    }


    // * ---------------------------------------------------------------------- *
    // * Actual time between the last ADC of the current kernel call and the    *
    // * following kernel call.                                                 * 
    // * ---------------------------------------------------------------------- *
    long lActualDelayBetweenADC = m_lRephaseTime - m_lOverTime + m_lExcitationAllTime + m_lPrephaseTime + m_alTEFill[0];


    // * ---------------------------------------------------------------------- *
    // * Calculate the minimum delay between the last ADC of the current kernel *
    // * call and the following kernel call.                                    *
    // * ---------------------------------------------------------------------- *
    long lMinDelayBetweenADC_PCI = 0;
    
    if ( pMrProt->phaseStabilize() )  {
        lMinDelayBetweenADC_PCI = RoundUpGRT (m_PhaseStabScan.getlMinDurationToNextADC (m_aRO[0].ADC()));
    } else {
        lMinDelayBetweenADC_PCI = RoundUpGRT (m_aRO[m_lNumberOfEchoes-1].ADC().getMinDistanceBetweenReadout(m_aRO[0].ADC(), true));
    }


    if ( lActualDelayBetweenADC < lMinDelayBetweenADC_PCI )  {
        m_lTRMin += lMinDelayBetweenADC_PCI - lActualDelayBetweenADC;
    }



    // * ---------------------------------------------------------------------- *
    // * Increase TR if necessary in order to avoid acoustic resonances         *
    // * ---------------------------------------------------------------------- *
    if ( m_bAvoidAcousticResonances )  {

        long lTRTemp = getAllowedTR_us (m_lTRMin);


        // * ------------------------------------------------------------------ *
        // * Set TRFill time if TR has to be increased due to acoustic          *
        // * resonances                                                         *
        // * ------------------------------------------------------------------ *
        if ( lTRTemp > m_lTRMin )  {
            m_lTRFill = lTRTemp - m_lTRMin;
        } else {
            m_lTRFill = 0;
        }
        
    }


    return ( true );


}



// * -------------------------------------------------------------------------- *
// *                                                                            *
// * Name        :  SBBGREKernel::calculateROStartTime                          *
// *                                                                            *
// * Description :  Calculation of the start times of the readout SBBs          *
// *                                                                            *
// * Return      :  bool                                                        *
// *                                                                            *
// * -------------------------------------------------------------------------- *
bool SBBGREKernel::calculateROStartTime (MrProt*, SeqLim*, SeqExpo*)
{

    long lI         = 0;
    long lStartTime = m_lExcitationAllTime + m_lPrephaseTime;

    for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {

        lStartTime += m_alTEFill[lI];
        m_aRO[lI].setlStartTime (lStartTime);
        lStartTime += m_aRO[lI].getDurationPerRequest();

    }

    return ( true );

}






// ******************************************************************************
// * Name        : prep                                                         *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Additional calculations. Not used in GREKernel.              *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::calculateAddOn (MrProt* , SeqLim* , SeqExpo* )
{

    return ( true );
    
}





bool SBBGREKernel::calculateAddOnPost (MrProt* , SeqLim* , SeqExpo*)
{

    return ( true );
    
}




// ******************************************************************************
// * Name        : prep                                                         *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Prepare and check of the gardients, tables and ADC events    *
// *               Calculation of the energy per request and the duration per   *
// *               request.                                                     *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::prep (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{

    static const char *ptModule = "SBBGREKernel::prep";

    double            dAmplitudeOfLeftNeighbour = 0.0;  // * Used for gradient slew rate checks *



    // *------------------------------------------------------------------------*
    // * Debug: running message and imported parameters                         *
    // *------------------------------------------------------------------------*
    if ( mIsTRrun )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s running\n", ptModule);
    }



    if ( mIsTRinp )  {
        TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Imported Parameter                                *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);

        if ( isTimingCalculated() )  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Kernel was calculated.\n" ,ptModule);
        } else  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s Kernel has not been calculated, yet.\n" ,ptModule);
        }
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    }


    if ( ! isTimingCalculated() )  {

        if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
            setNLSStatus(SEQU_ERROR, ptModule, "Kernel timing has not been calculated, yet.\n");     return ( false );
        } else {
            setNLSStatus(SEQU_ERROR);     return ( false );
        }

    }



    // *------------------------------------------------------------------------*
    // * Prepare the osc bit                                                    *
    // *------------------------------------------------------------------------*
    if ( m_bSendOscBit == true )  {
        m_sOscBit.setIdent    ( RTEIDENT_Osc0 );
        m_sOscBit.prep        ( 0, 10 );          // * Channel: 0, duration: 10us *
    }




    // *------------------------------------------------------------------------*
    // * Prepare the wake up bit                                                *
    // * WakeUp event must have length zero. This is set already by the         *
    // * constructor!                                                           *
    // *------------------------------------------------------------------------*
  
    // WakeUp event must have length zero. This is set already by the constructor!
  
    m_sWakeUp.lCode       = SYNCCODE_WAKEUP;
    m_sWakeUp.setIdent    ( "WakeUp" );



    // * ---------------------------------------------------------------------- *
    // * Check the excitation pulse                                             *
    // * ---------------------------------------------------------------------- *
    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }
  
    if (! m_pRF_WE->checkGradients(pMrProt, pSeqLim) )  {
        if ( ! pSeqLim->isContextPrepForBinarySearch() ) {
            setNLSStatus ( m_pRF_WE->getNLSStatus(), "m_pRF_WE->checkGradients" );    return ( false );
        } else {
            setNLSStatus ( m_pRF_WE->getNLSStatus() );    return ( false );
        }
    }



    // * ---------------------------------------------------------------------- *
    // * Prepare and check gradient table rewinder in slice select direction    *
    // * ---------------------------------------------------------------------- *
    if (! prepAndCheck3D( pMrProt, pSeqLim ) ) return ( false );
  


    // * ---------------------------------------------------------------------- *
    // * Prepare and check phase encoding table                                 *
    // * ---------------------------------------------------------------------- *
    if (! m_TB_PE.updateAmplitude (pSeqLim, m_lFirstLineToMeasure) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            setNLSStatus (m_TB_PE.getNLSStatus(), "m_TB_PE.updateAmplitude(m_lFirstLineToMeasure)");
        } else  {
            setNLSStatus (m_TB_PE.getNLSStatus());
        }
        return ( false );
    }


    if (! m_TB_PE.updateAmplitude (pSeqLim, m_lLastLineToMeasure) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            setNLSStatus (m_TB_PE.getNLSStatus(), "m_TB_PE.updateAmplitude(m_lLastLineToMeasure)");
        } else  {
            setNLSStatus (m_TB_PE.getNLSStatus());
        }
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Prepare and check phase encoding rewinder table                        *
    // * ---------------------------------------------------------------------- *
    if ( m_bUsePERewinder == true) {
        if (! m_sTB_PER.prepPE(m_sTB_PER.getRampUpTime(), m_sTB_PER.getDuration(), m_sTB_PER.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, m_lFirstLineToMeasure))  {
            setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE");
            return ( false );
        }

        if (! m_sTB_PER.check() )  {
            setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.check");
            return ( false );
        }

        if (! m_sTB_PER.prepPE(m_sTB_PER.getRampUpTime(), m_sTB_PER.getDuration(), m_sTB_PER.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, m_lLastLineToMeasure))  {
            setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE");
            return ( false );
        }

        if (! m_sTB_PER.check() )  {
            setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.check");
            return ( false );
        }

    } else {
        if (! m_sTB_PER.prepPE(m_sTB_PER.getRampUpTime(), m_sTB_PER.getDuration(), m_sTB_PER.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, 0) )  {
            setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE");
            return ( false );
        }
    }


    // * ---------------------------------------------------------------------- *
    // * Prepare and check read out prephasing gradient                         *
    // * ---------------------------------------------------------------------- *
    if (! prepAndCheckROP( pMrProt, pSeqLim ) ) return ( false );



    // * ---------------------------------------------------------------------- *
    // * Prepare and check read out dephasing gradient                          *
    // * ---------------------------------------------------------------------- *
    if (! m_sP_ROD.prepAmplitude(m_sP_ROD.getRampUpTime(), m_sP_ROD.getDuration(), m_sP_ROD.getRampDownTime(), m_sP_ROD.getAmplitude() ))  {
        setNLSStatus (m_sP_ROD.getNLSStatus(), "SBBGREKernel::prep m_sP_ROD.prepAmplitude");
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Prepare and check additional gradients                                 *
    // * ---------------------------------------------------------------------- *
    if (! prepAndCheckAddOn ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );



    if ( pMrProt->phaseStabilize() )  {
        dAmplitudeOfLeftNeighbour = m_PhaseStabScan.getROAmplitude();            // * Amplitude of phase stabilize RO *
    } else {
        dAmplitudeOfLeftNeighbour = m_aRO[m_lNumberOfEchoes-1].getdROAmplitude(); // * Amplitude of last kernel RO     *
    }

    switch ( m_eReadOutGradType )  {
        case Spoiler:
            if (! m_sP_ROD.check( dAmplitudeOfLeftNeighbour, 0.0 ) )  {
                setNLSStatus (m_sP_ROD.getNLSStatus(), ptModule, "SBBGREKernel::prep m_sP_ROD.check", pSeqLim->isContextNormal());
                return ( false );
            }
        break;

        case Constant:
            if (! m_sP_ROD.check( dAmplitudeOfLeftNeighbour, 0.0 ) )  {
                setNLSStatus (m_sP_ROD.getNLSStatus(), ptModule, "SBBGREKernel::prep m_sP_ROD.check", pSeqLim->isContextNormal());
                return ( false );
            }
        break;

        case Rewinder:
        case Symmetric_Rewinder:
            if (! m_sP_ROD.check() )  {
                setNLSStatus (m_sP_ROD.getNLSStatus(), ptModule, "SBBGREKernel::prep m_sP_ROD.check", pSeqLim->isContextNormal());
                return ( false );
            }
        break;
    
        default: 
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." ) )  return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Calculate energy of one request                                        *
    // * ---------------------------------------------------------------------- *
    m_dEnergyPerRequest_Ws = m_pRF_WE->getEnergyPerRequest();


    // * ---------------------------------------------------------------------- *
    // * Calculate SBB duration of one request                                  *
    // * ---------------------------------------------------------------------- *
    m_lSBBDurationPerRequest_us = m_lTRMin + m_lTRFill + m_lTRFillEnd;




    if ( mIsTRend )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s finished\n", ptModule);
    }



    // * ---------------------------------------------------------------------- *
    // * Successful preparation of the GRE kernel                               *
    // * ---------------------------------------------------------------------- *
    setPrepared();
    
    return ( true );

}




// *****************************************************************************
// *                                                                            *
// * Name        : prepAndCheck3D                                               *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Prepares and checks the 3D gradient tables.                  *
// *               Returns true if prep and check has been successfully         *
// *               been finished otherwise false.                               *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::prepAndCheck3D (MrProt* pMrProt, SeqLim* pSeqLim)
{
  
    // * ---------------------------------------------------------------------- *
    // * Prepare and check gradient table in slice select direction             *
    // * ---------------------------------------------------------------------- *
    if (! m_TB_3D.updateAmplitude (pSeqLim, m_lFirstPartToMeasure) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            setNLSStatus (m_TB_3D.getNLSStatus(), "m_TB_3D.updateAmplitude(...)");
        } else  {
            setNLSStatus (m_TB_3D.getNLSStatus());
        }
        return ( false );
    }


    if (! m_TB_3D.updateAmplitude (pSeqLim, m_lLastPartToMeasure) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            setNLSStatus (m_TB_3D.getNLSStatus(), "m_TB_3D.updateAmplitude(...)");
        } else  {
            setNLSStatus (m_TB_3D.getNLSStatus());
        }
        return ( false );
    }


 
    if ( m_bFlowCompSliceSelect )  {

        if (! m_sP_3D.prepAmplitude(m_sP_3D.getRampUpTime(), m_sP_3D.getDuration(), m_sP_3D.getRampDownTime(), m_sP_3D.getAmplitude() ))  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                setNLSStatus (m_sP_3D.getNLSStatus(), "m_sP_3D.prepAmplitude");
            } else  {
                setNLSStatus (m_sP_3D.getNLSStatus());
            }
            return ( false );
        }

        if (! m_sP_3D.check())  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                setNLSStatus (m_sP_3D.getNLSStatus(), "m_sP_3D.check");
            } else  {
                setNLSStatus (m_sP_3D.getNLSStatus());
            }
            return ( false );
        }



        if (! m_sP_3D_FC.prepAmplitude(m_sP_3D_FC.getRampUpTime(), m_sP_3D_FC.getDuration(), m_sP_3D_FC.getRampDownTime(), m_sP_3D_FC.getAmplitude() ))  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                setNLSStatus (m_sP_3D_FC.getNLSStatus(), "m_sP_3D_FC.prepAmplitude");
            } else  {
                setNLSStatus (m_sP_3D_FC.getNLSStatus());
            }
            return ( false );
        }

        if (! m_sP_3D_FC.check())  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                setNLSStatus (m_sP_3D_FC.getNLSStatus(), "m_sP_3D_FC.check");
            } else  {
                setNLSStatus (m_sP_3D_FC.getNLSStatus());
            }
            return ( false );
        }

    } else {

        if ( m_bFlowCompPartitionEncode )  {

            if (! m_sP_3D.prepAmplitude(m_sP_3D.getRampUpTime(), m_sP_3D.getDuration(), m_sP_3D.getRampDownTime(), m_sP_3D.getAmplitude() ))  {
                if (! pSeqLim->isContextPrepForBinarySearch() )  {
                    setNLSStatus (m_sP_3D.getNLSStatus(), "m_sP_3D.prepAmplitude");
                } else  {
                    setNLSStatus (m_sP_3D.getNLSStatus());
                }
                return ( false );
            }

            if (! m_sP_3D.check())  {
                if (! pSeqLim->isContextPrepForBinarySearch() )  {
                    setNLSStatus (m_sP_3D.getNLSStatus(), "m_sP_3D.check");
                } else  {
                    setNLSStatus (m_sP_3D.getNLSStatus());
                }
                return ( false );
            }

        }

    }



    // * ---------------------------------------------------------------------- *
    // * Prepare and check gradient table rewinder in slice select direction    *
    // * ---------------------------------------------------------------------- *
    if ( m_bUsePERewinder == true) {
        if ( ! m_bBalancedGradMomentSliceSelect )  { 
            m_dMomentum_TB_3DR = -1.0 * m_dMomentum_TB_3D * m_dRelSSDephasingMoment;
        } else { 
            m_dMomentum_TB_3DR = m_dMomentum_TB_3D;    // * Balanced slice selecet gradients, used for True Fisp *
        }
    
        if (! m_sTB_3DR.prep3D(m_sTB_3DR.getRampUpTime(), m_sTB_3DR.getDuration(), m_sTB_3DR.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, m_dMomentum_TB_3DR, m_lFirstPartToMeasure))  {
            setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D");
            return ( false );
        }

        if (! m_sTB_3DR.check() )  {
            setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.check");
            return ( false );
        }

        if (! m_sTB_3DR.prep3D(m_sTB_3DR.getRampUpTime(), m_sTB_3DR.getDuration(), m_sTB_3DR.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, m_dMomentum_TB_3DR, m_lLastPartToMeasure))  {
            setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D");
            return ( false );
        }

        if (! m_sTB_3DR.check() )  {
            setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.check");
            return ( false );
        }

    } else {
        if (! m_sTB_3DR.prep3D(m_sTB_3DR.getRampUpTime(), m_sTB_3DR.getDuration(), m_sTB_3DR.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, 0) )  {
            setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D");
            return ( false );
        }
    }

    return ( true );

}




bool SBBGREKernel::prepAndCheckROP (MrProt* , SeqLim* )
{
    // **********************************************************************************
    // * Name        : prepAndCheckROP                                                  *
    // *                                                                                *
    // * Class       : SBBGREKernel                                                     *
    // *                                                                                *
    // * Description : Prepares and checks the read out prephasing gradients.           *
    // *               Returns true if prep and check has been successfully been        *
    // *               finished otherwise false.                                        *
    // *                                                                                *
    // * Return      : bool                                                             *
    // *                                                                                *
    // **********************************************************************************




    // * -------------------------------------------------------------------------- *
    // * Read out flow compensation gradient                                        *
    // * -------------------------------------------------------------------------- *
    if ( m_bFlowCompRead ) {

        if (! m_sP_ROP_FC.prepAmplitude(m_sP_ROP_FC.getRampUpTime(), m_sP_ROP_FC.getDuration(), m_sP_ROP_FC.getRampDownTime(), m_sP_ROP_FC.getAmplitude() ))  {
            setNLSStatus (m_sP_ROP_FC.getNLSStatus(), "m_sP_ROP_FC.prepAmplitude");
            return ( false );
        }
  
        if (! m_sP_ROP_FC.check())  {
            setNLSStatus (m_sP_ROP_FC.getNLSStatus(), "m_sP_ROP_FC.check");
            return ( false );
        }

    }



    // * -------------------------------------------------------------------------- *
    // * Read out prephasing gradient                                               *
    // * -------------------------------------------------------------------------- *
    if (! m_sP_ROP.prepAmplitude(m_sP_ROP.getRampUpTime(), m_sP_ROP.getDuration(), m_sP_ROP.getRampDownTime(), m_sP_ROP.getAmplitude() ))  {
        setNLSStatus (m_sP_ROP.getNLSStatus(), "m_sP_ROP.prepAmplitude");
        return ( false );
    }
  
    if (! m_sP_ROP.check())  {
        setNLSStatus (m_sP_ROP.getNLSStatus(), "m_sP_ROP.check");
        return ( false );
    }

    return ( true );


}



bool SBBGREKernel::prepAndCheckAddOn (MrProt* , SeqLim* , SeqExpo* )
{
  // **********************************************************************************
  // * Name        : run                                                              *
  // *                                                                                *
  // * Class       : SBBGREKernel                                           *
  // *                                                                                *
  // * Description : Additional preparations. Not used in GREKernel.                  *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

    return ( true );

}

// ******************************************************************************
// * Name        : run                                                          *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Executes the basic GRE kernel                                *
// *                                                                            *
// * Return      : bool                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::run (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, sSLICE_POS* pSLC)
{


    static const char *ptModule = "SBBGREKernel::run";


    long         lI;                         // * Loop counter *
  
    long         lGradLineNumber      = 0;   // * Line/partition number for gradient    *
    long         lGradPartNumber      = 0;   // * table amplitude calcilation           *
                                             // * lGradLineNumber = 0 -> k-space center *

    // * ---------------------------------------------------------------------- *
    // * Calculate line and partition number for gradient table amplitude       *
    // * calculation                                                            *
    // * ---------------------------------------------------------------------- *
    lGradLineNumber = m_lLineNumber - m_lKSpaceCenterLine;
    lGradPartNumber = m_lPartNumber - m_lKSpaceCenterPartition;



    // * ---------------------------------------------------------------------- *
    // * Set the Mdh of ADC[1] with data from ADC[0]                            *
    // * ---------------------------------------------------------------------- *
    for ( lI = 1; lI < m_lNumberOfEchoes; lI++ )  {
        m_aRO[lI].Mdh() = m_aRO[0].Mdh();
    }



    // * ---------------------------------------------------------------------- *
    // * Move LastScanInConcat and LastScanInMeas Mdh flags from the first to   *
    // * the last ADC.                                                          *
    // * ---------------------------------------------------------------------- *
    if ( m_lNumberOfEchoes > 1 )  {
        if ( m_aRO[0].Mdh().isLastScanInConcat() )  {
            for ( lI = 0; lI < m_lNumberOfEchoes-1; lI++ )  {
                m_aRO[lI].Mdh().deleteFromEvalInfoMask ( MDH_LASTSCANINCONCAT );
                m_aRO[lI+1].Mdh().addToEvalInfoMask ( MDH_LASTSCANINCONCAT );
            }
        }

      if ( m_aRO[0].Mdh().isLastScanInMeas() )  {
          for ( lI = 0; lI < m_lNumberOfEchoes-1; lI++ )  {
            m_aRO[lI].Mdh().deleteFromEvalInfoMask ( MDH_LASTSCANINMEAS );
            m_aRO[lI+1].Mdh().addToEvalInfoMask ( MDH_LASTSCANINMEAS );
          }
      }

    }




    // * ---------------------------------------------------------------------- *
    // * Remove the phase and/or slice FFT flags from the second contrast ADC   *
    // * if phase stabilization is selected                                     *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->phaseStabilize() && ( pMrProt->contrasts() == 2 ) )  {
        m_aRO[1].Mdh().deleteFromEvalInfoMask ( MDH_PHASEFFT );
        m_aRO[1].Mdh().deleteFromEvalInfoMask ( MDH_D3FFT    );
    }




    // * ---------------------------------------------------------------------- *
    // * If last scan in concat or meas, the corresponding flags have to set    *
    // * for the phase stabiliazation ADC.                                      *
    // * ---------------------------------------------------------------------- *
    m_PhaseStabScan.setbLastScanInConcat( m_bLastScanInConcat );
    m_PhaseStabScan.setbLastScanInMeas  ( m_bLastScanInMeas   );




    // * ---------------------------------------------------------------------- *
    // * General function to set Mdh-flags                                      *
    // * ---------------------------------------------------------------------- *
    runSetMdh();



    for ( lI=0; lI<m_lNumberOfEchoes; lI++ )  {


        // * ------------------------------------------------------------------ *
        // * Activate online ice-process                                        *
        // * ------------------------------------------------------------------ *
        m_aRO[lI].Mdh().addToEvalInfoMask ( MDH_ONLINE );
    
        // * ------------------------------------------------------------------ *
        // * Set the line, partition, echo and set parameter                    *
        // * ------------------------------------------------------------------ *
        m_aRO[lI].Mdh().setClin ( m_lLineNumber );
        m_aRO[lI].Mdh().setCpar ( m_lPartNumber );


        runSetMdhCeco ( lI );
        runSetMdhCset ( lI );
    


        // * ---------------------------------------------------------------------- *
        // * Set setKSpaceCentreColumn                                              *
        // * ---------------------------------------------------------------------- *
        m_aRO[lI].Mdh().setKSpaceCentreColumn ( (unsigned short)(m_aRO[lI].getdExactReadOutAsym() * pMrProt->kSpace().baseResolution() + 0.5) );




        // * ---------------------------------------------------------------------- *
        // * Specify k-space center line and partition numbers                      *
        // * ---------------------------------------------------------------------- *
        m_aRO[lI].Mdh().setKSpaceCentreLineNo      ( m_lKSpaceCenterLine      );
        m_aRO[lI].Mdh().setKSpaceCentrePartitionNo ( m_lKSpaceCenterPartition );




        // * -------------------------------------------------------------------------- *
        // * Define the time(s) between RF and ADCs for phase stabilization             *
        // * -------------------------------------------------------------------------- *
        if ( pMrProt->phaseStabilize() )  {
            m_aRO[lI].Mdh().setTimeSinceLastRF ( pMrProt->te()[lI] );
        }

    }




    // * ---------------------------------------------------------------------- *
    // * Set the frequency and phase properties of the rf pulse and ADC         *
    // * ---------------------------------------------------------------------- *
    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }

    double dPhaseCycle = runPhaseCycle();

    m_pRF_WE->setAdditionalPhase ( dPhaseCycle + m_dRFSpoilPhase );

    for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {
        m_aRO[lI].updateFreqPhase (lGradLineNumber, lGradPartNumber, pSLC);
        m_aRO[lI].increaseADCPhase (dPhaseCycle + m_dRFSpoilPhase);
    }




    // * ---------------------------------------------------------------------- *
    // * Calculate phase encode gradient amplitudes                             *
    // * ---------------------------------------------------------------------- *
    if (! m_TB_PE.updateAmplitude (pSeqLim, lGradLineNumber) )  { 
        if ( pSeqLim->isContextNormal() )  {
            setNLSStatus (m_TB_PE.getNLSStatus(), "Update of gradient amplitude for next phase encoding step failed");
        } else {
            setNLSStatus (m_TB_PE.getNLSStatus());
        }
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Update the gradient amplitude of the partition encoding gradient table *
    // * ---------------------------------------------------------------------- *    
    if ( ! m_TB_3D.updateAmplitude (pSeqLim, lGradPartNumber) )  {
        if ( pSeqLim->isContextNormal() )  {
            setNLSStatus (m_TB_3D.getNLSStatus(), ptModule, "Update of gradient amplitude for next partition encoding step failed");
        } else {
            setNLSStatus (m_TB_3D.getNLSStatus());
        }
        return (false);
    }



    if ( m_bUsePERewinder == true ) {     // * Prepare rewinder amplitudes             *
        if (! m_sTB_PER.prepPE(pMrProt, lGradLineNumber ) )    { setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE"); return ( false ); }
        if (! m_sTB_3DR.prep3D(pMrProt, lGradPartNumber ) )    { setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D"); return ( false ); }
    } else {                       // * No rewinder -> amplitude equal zero !!! *
        if (! m_sTB_PER.prepAmplitude (0.0) )    { setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE"); return ( false ); }
        if (! m_sTB_3DR.prepAmplitude (0.0) )    { setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D"); return ( false ); }
    }

#if defined SEND_TTC_DATA_IN_MDH
    // * ------------------------------------------------------------------------------------- *
    // * For display of timestamps in the image, we need to know the amount of time elapsed
    // * from the start of the measurement until the center of k-space is reached for the very
    // * first measurement-repeat.
    // * ------------------------------------------------------------------------------------- *
    if (fRTIsReadoutEnabled())
    {
        if (!m_bFoundFirstTimeStamp)
        {
            // Remember the time preceeding the first ADC.  This value is passed below
            // in a sixteen-bit integer, so we require that the lead time is less than
            // ~65 seconds.  It is highly unlikely that the first adc appears after
            // 65 seconds, but just in case it is, it is limited to a maximum of 65 seconds.
            m_bFoundFirstTimeStamp = true;
            m_dLeadTimeMSec = RTController::getInstance().getAbsTimeOfEventBlockMSec();
            double dMax = 65530.;
            if (m_dLeadTimeMSec > dMax)
            {
                TRACE_PUT3(TC_INFO, TF_SEQ, "%s : m_dLeadTimeMSec %f too large for 16-bit integer; truncated to %f",ptModule,m_dLeadTimeMSec,dMax);
                m_dLeadTimeMSec = dMax;
            }
        }

        // Deliver the lead time and the timestamp of each ADC to the ICE program.
        // The timestamp is delivered after it is split into two parts: a) whole seconds
        // (ushTime_sec) and b) the remaining milliseconds (ushTime_msec).  This is done
        // because the MDH parameter which stores the values is only 16 bit.  If we used
        // only one parameter we would have enough bits to store a maximum time of only
        // 2**16-1 msec ~ 65 seconds.
        double dTMSec    = RTController::getInstance().getAbsTimeOfEventBlockMSec();
        double dTimeMSec = dTMSec - m_dLeadTimeMSec;
        unsigned short ushLeadTime_msec = ((unsigned short) (m_dLeadTimeMSec+0.5));
        unsigned short ushTime_sec      = ((unsigned short) floor(dTimeMSec/1000.));
        unsigned short ushTime_msec     = ((unsigned short) (dTimeMSec - (1000.*ushTime_sec)));
        for ( lI=0; lI<m_lNumberOfEchoes; lI++ )
        {
            m_aRO[lI].Mdh().setIceProgramPara(MDHFREE_ADVMRA_LEADTIME_MS,ushLeadTime_msec);
            m_aRO[lI].Mdh().setIceProgramPara(MDHFREE_ADVMRA_CLOCK_SEC  ,ushTime_sec);
            m_aRO[lI].Mdh().setIceProgramPara(MDHFREE_ADVMRA_CLOCK_MSEC ,ushTime_msec);
        }

        // Set the time-to-center value
        if ((lGradLineNumber == 0) && (lGradPartNumber == 0))
        {
            m_dLastTimeToCenterPlayedOutSec = dTMSec/1000.;
        }
        
#ifdef DEBUG
        const bool isCenter=((lGradLineNumber == 0) && (lGradPartNumber == 0));

        // Display the time when the center of kSpace is played out.
        if ((lGradLineNumber == 0) && (lGradPartNumber == 0))
        {
            TRACE_PUT3(TC_ALWAYS, TF_SEQ, "%s : TTC %8.3f sec;%s",ptModule,dTMSec/1000.,(isCenter?" (*C*)":""));
        }
#endif

    }
#endif


    // * ---------------------------------------------------------------------- *
    // * Execute the kernel                                                     *
    // * ---------------------------------------------------------------------- *
    if (! runTiming ( pMrProt, pSeqLim, pSeqExpo, pSLC, lGradLineNumber, lGradPartNumber ) )  {
        return ( false );
    } else {
        return ( true  );
    }


}




// ******************************************************************************
// * Name        : runSetMdh                                                    *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : General function to set Mdh-flags                            *
// *                                                                            *
// * Return      : void                                                         *
// *                                                                            *
// ******************************************************************************
void SBBGREKernel::runSetMdh ()
{}




// ******************************************************************************
// * Name        : runSetMdhCeco                                                *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Sets the Mdh Ceco index                                      *
// *                                                                            *
// * Return      : void                                                         *
// *                                                                            *
// ******************************************************************************
void SBBGREKernel::runSetMdhCeco (long Ceco)
{
    MdhProxy& rMDH = this->m_aRO[Ceco].Mdh();
    rMDH.setCeco ( Ceco );
    if( Ceco > 0 )
    {
        //  ECO is pf type ICE_PP_FIX therefore only one Soda object for all echoes
        rMDH.setIsMdsReferencePosition(false);
        m_aRO[Ceco].ADC().setCoilSelectIndex(m_aRO[0].ADC().getCoilSelectIndex());
    }  
}




// ******************************************************************************
// * Name        : runSetMdhCset                                                *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : -                                                            *
// *                                                                            *
// * Return      : void                                                         *
// *                                                                            *
// ******************************************************************************
void SBBGREKernel::runSetMdhCset (long Cset)
{

    //lint -e{550}
    if( false ) { long lDummy; lDummy = Cset; }

}




double SBBGREKernel::runPhaseCycle ()
{
    return ( 0.0 );
}





// ******************************************************************************
// * Name        : runTiming                                                    *
// *                                                                            *
// * Class       : SBBGREKernel                                                 *
// *                                                                            *
// * Description : Executes the GRE kernel                                      *
// *                                                                            *
// * Return      : BOOL                                                         *
// *                                                                            *
// ******************************************************************************
bool SBBGREKernel::runTiming (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, sSLICE_POS* pSLC, long lGradLineNumber, long lGradPartNumber)
{

    static const char *ptModule = "SBBGREKernel::runTiming";

    long lT             = 0;
    long lI             = 0;
    long lStatus        = SEQU__NORMAL;
    long lOffsetTime    = 0;
    //  m_lTRFillEnd may be modified by runAdditionalTiming
    long lTRFillEnd     = m_lTRFillEnd;


    //lint -e{550}
    if( false ) { long lDummy; lDummy = lGradLineNumber = lGradPartNumber; }


    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL\n");
        return ( false );
    }



    // * ------------------------------------------------------------------------------------------------------------------- *
    // * |             Start Time        |     NCO     |    SRF   |    ADC   |          Gradient Events          | Sync    | *
    // * |               (usec)          |    Event    |   Event  |   Event  |   phase   |   read    |   slice   | Event   | *
    // * ------------------------------------------------------------------------------------------------------------------- *
    fRTEBInit( &pSLC->m_sROT_MATRIX,this->m_iAdjUpdate);
    this->m_iAdjUpdate = RT_MDSUPDATE_ADJ_NONE;
    m_bUpdateRotationMatrix = false;




    // * ---------------------------------------------------------------------- *
    // * Execute some additional timing elements                                *
    // * The duration of these timining elements has to be specified using      *
    // * lOffsetTime                                                            *
    // * ---------------------------------------------------------------------- *
    if ( ! runAdditionalTiming ( pMrProt, pSeqLim, pSeqExpo, pSLC, &lOffsetTime ) )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Error during execution of GREKernel::runAdditionalTiming().");
        return ( false );
    }

    lT += lOffsetTime;





    // * ---------------------------------------------------------------------- *
    // *     Execute the osc bit                                                *
    // * ---------------------------------------------------------------------- *
    if ( m_bSendOscBit )  {
        fRTEI(lT                         ,            0,         0,         0,          0,          0,          0,&m_sOscBit);
    }
  
  
    // * ---------------------------------------------------------------------- *
    // *     Execute the wake up bit                                            *
    // * ---------------------------------------------------------------------- *
    if ( m_bSendWakeUpBit )  {
        fRTEI(lOffsetTime+m_lStartTimeWakeUp ,        0,         0,         0,          0,          0,          0,&m_sWakeUp);
    }

  
    // * ---------------------------------------------------------------------- *
    // *     Execute the water selective excitation pulse                       *
    // * ---------------------------------------------------------------------- *
    // hughtin6: but first we must set the RunIndex so the correct flip angle is played out
    // we must first check to see if there is a flip angle array and if not then we should not set this
    if(m_padFlipAngleArray != NULL) m_pRF_WE->setRunIndex (m_lRunIndex);                                                            // sets the flip angle we are going to run
    m_pRF_WE->setStartTimeInEventBlock ( lOffsetTime+m_pRF_WE->getRequiredLeadTime() );
    if (! m_pRF_WE->run(pMrProt, pSeqLim, pSeqExpo, pSLC) )  {    
        setNLSStatus(m_pRF_WE->getNLSStatus(), ptModule, "Error in m_pRF_WE->run(...) at runtime");
        return ( false );
    }



    // * ---------------------------------------------------------------------- *
    // *   Phase encode and slice select events                                 *
    // * ---------------------------------------------------------------------- *
    lT = lOffsetTime + m_lExcitationAllTime + m_lDelayAdjustSSPosition;


    if ( m_bFlowCompSliceSelect )  {

        fRTEI(lT + m_pRF_WE->getRampTimeOutsideSBB(), 0,         0,         0,          0,          0,&m_sP_3D_FC,        0);
        fRTEI(lT + m_sP_3D_FC.getTotalTime() +
               m_pRF_WE->getRampTimeOutsideSBB(),     0,         0,         0,          0,          0,   &m_sP_3D,        0);

        if ( !m_bFlowCompSliceSelect && !m_bFlowCompPartitionEncode )  {
            m_TB_3D.setlStartTime (lT + m_alTEFill[0] + m_lPrephaseFillTimeSlice);
        } else {
            m_TB_3D.setlStartTime (lT + m_pRF_WE->getRampTimeOutsideSBB() + m_sP_3D.getTotalTime() + m_sP_3D_FC.getTotalTime() + m_alTEFill[0] + m_lPrephaseFillTimeSlice);
        }
        m_TB_3D.run (pMrProt, pSeqLim, pSeqExpo, pSLC);

    } else {

        if ( m_bFlowCompPartitionEncode )  {

            fRTEI(lT + m_pRF_WE->getRampTimeOutsideSBB(), 0,     0,         0,          0,          0,   &m_sP_3D,        0);

            m_TB_3D.setlStartTime (lT + m_sP_3D.getTotalTime() + m_alTEFill[0] + m_lPrephaseFillTimeSlice);
            m_TB_3D.run (pMrProt, pSeqLim, pSeqExpo, pSLC);

        } else {

            m_TB_3D.setlStartTime (lT + m_pRF_WE->getRampTimeOutsideSBB());
            m_TB_3D.run (pMrProt, pSeqLim, pSeqExpo, pSLC);

        }

    }


    m_TB_PE.setlStartTime (lT + m_alTEFill[0] + m_lPrephaseFillTimePhase);
    m_TB_PE.run (pMrProt, pSeqLim, pSeqExpo, pSLC);


    lT += m_lPrephaseTime + m_lReadoutTime;

    for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {  lT += m_alTEFill[lI];  }


    fRTEI(lT                             ,            0,         0,         0, &m_sTB_PER,          0, &m_sTB_3DR,        0);



    lT += m_lRephaseTime + m_lTRFill - m_lOverTime;
    if ( pMrProt->phaseStabilize()           )  {    lT+=m_PhaseStabScan.getDurationPerRequest() + m_sP_ROD.getTotalTime();    };

    fRTEI(lT                             ,            0,         0,         0,          0,          0,          0,        0);

    if (m_lTRFillEnd>0)  {
        fRTEI(lT+=m_lTRFillEnd           ,            0,         0,         0,          0,          0,          0,        0);
    }




    // * ---------------------------------------------------------------------- *
    // *     lT set to zero                                                     *
    // *     read out events                                                    *
    // * ---------------------------------------------------------------------- *
    lT = lOffsetTime + m_lExcitationAllTime + m_lDelayAdjustSSPosition;

    if ( m_bFlowCompRead )  {
        fRTEI(lT + m_alTEFill[0] +
              m_lPrephaseFillTimeRead ,               0,         0,         0,          0,&m_sP_ROP_FC,         0,        0);
    }

    fRTEI(lT + m_sP_ROP_FC.getTotalTime() +
               m_alTEFill[0] +
               m_lPrephaseFillTimeRead   ,            0,         0,         0,          0,  &m_sP_ROP,          0,        0);
    lT += m_lPrephaseTime;


    // * ---------------------------------------------------------------------- *
    // * Execute readout modules                                                *
    // * ---------------------------------------------------------------------- *

    if ( pMrProt->phaseStabilize() )  {

        lT += m_lReadoutTime + m_lRephaseTime;
        for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {  lT += m_alTEFill[lI];  }

        m_PhaseStabScan.setlStartTime ( lT );
        m_PhaseStabScan.setlTimeFromRFToStartSBB ( m_PhaseStabScan.getlStartTime() - m_lExcitationAllTime + m_lTEFromExcitation );
        m_PhaseStabScan.setdRFSpoilPhase ( m_dRFSpoilPhase );
        m_PhaseStabScan.run ( pMrProt, pSeqLim, pSeqExpo, pSLC );

        lT += m_PhaseStabScan.getDurationPerRequest();

    } else {

        lT += m_lReadoutTime;
        for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {  lT += m_alTEFill[lI];  }

        if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) ) {  lT += m_aRO[m_lNumberOfEchoes-1].getlRampDownTime();  }

    }

    fRTEI(lT                               ,          0,         0,         0,          0,  &m_sP_ROD,          0,        0);


    // * ---------------------------------------------------------------------- *
    // * Execute readout modules                                                *
    // *                                                                        *
    // *                        I M P O R T A N T                               *
    // *                        -----------------                               *
    // * m_PhaseStabScan.run(...) modifies the Mdh of m_aROs Mdh. For that      *
    // * reason the run function of the phase stabilization scan HAS to         *
    // * executed prior to the run function of m_aRO.                           *
    // * ---------------------------------------------------------------------- *
    for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {
        const long lStartTime_us = m_aRO[lI].getlStartTime();
        m_aRO[lI].setlStartTime(lStartTime_us+lOffsetTime);
        m_aRO[lI].run (pMrProt, pSeqLim, pSeqExpo, pSLC);
        m_aRO[lI].setlStartTime(lStartTime_us);
    }






    #ifndef VXWORKS
        mSEQTest(pMrProt, pSeqLim, pSeqExpo, m_ulUTIdent , 30, lGradLineNumber, pSLC->getSliceIndex(), 0, lGradPartNumber); /*! EGA-All !*/
    #endif

    lStatus = fRTEBFinish();

    if ( setNLSStatus( lStatus ) )  {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s : fRTEBFinish() [*0030*] returned with an error : NLS_CODE : 0x%lx"  , ptModule, lStatus);
        return ( false );
    }

    //  Restore TRFillEnd
    m_lTRFillEnd = lTRFillEnd;


    return ( true );

}

//	This function can be used to execute additional timing 
//prior to the kernel.
//	Duration of this timing has to be specified in us using 
//the argument lOffset.
//##ModelId=3D352EBB0132
bool SBBGREKernel::runAdditionalTiming (MrProt* , SeqLim* , SeqExpo* , sSLICE_POS* , long* plDuration)
{

    *plDuration = 0;
  //  GS gradient of SBBExcitation is ramped up at SBBExcitation::StartTimeInEventBlock
    //  GREKErnel set SBBExcitation::StartTimeInEventBlock to lOffsetTime+m_pRF_WE->getRequiredLeadTime()
    if( this->m_pRF_WE )
    {
        GradAxisDataForSBBExcitation sGSData = m_pRF_WE->getGSData();
        if( this->m_lRampTimePreviousEB_us > m_pRF_WE->getRequiredLeadTime()+sGSData.getlTimeToFirstFlatTop() )
        {
            *plDuration = this->m_lRampTimePreviousEB_us-m_pRF_WE->getRequiredLeadTime();
            m_lTRFillEnd += m_pRF_WE->getRequiredLeadTime();
            m_sGFTss.setAxis(SEQ::AXIS_SLICE);
            m_sGFTss.setStartTime(0);
            m_sGFTss.set(this->m_lRampTimePreviousEB_us,this->m_lRampTimePreviousEB_us,sGSData.getlTimeToFirstFlatTop(),sGSData.getdFirstAmplitude());
            if( !this->m_sGFTss.prep() )
            {
                TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at '%s'(%d)",__FILE__,__LINE__);
                setNLSStatus(this->m_sGFTss.getNLSStatus());
                return false;
            }
            this->m_sGFTss.run();
        }
    }

    return ( true );

}



bool SBBGREKernel::calculatePreScanTiming (MrProt* , SeqLim* , SeqExpo* )
{
  // **********************************************************************************
  // * Name        : calculatePreScanTiming                                           *
  // *                                                                                *
  // * Class       : SBBGREKernel                                           *
  // *                                                                                *
  // * Description : Can be used to calculate the timing of the prescan for TrueFisp  *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

  return ( true );  // Returning without error
  
}


bool SBBGREKernel::prepPreScan (SeqLim* )
{
  // **********************************************************************************
  // * Name        : prepPreScan                                                      *
  // *                                                                                *
  // * Class       : SBBGREKernel                                           *
  // *                                                                                *
  // * Description : Can be used to trepare and check the gardients used for the      *
  // *               TrueFisp prescan                                                 *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

  return ( true);
  
}


bool SBBGREKernel::runPreScan (MrProt* , SeqLim* , SeqExpo* , sSLICE_POS* )
{
  // **********************************************************************************
  // * Name        : run                                                              *
  // *                                                                                *
  // * Class       : SBBGREKernel                                           *
  // *                                                                                *
  // * Description : Can be used to execute the prescan for TrueFisp                  *
  // *                                                                                *
  // * Return      : An NLS status code                                               *
  // *                                                                                *
  // **********************************************************************************

  return ( true );

}






bool SBBGREKernel::printImports (MrProt* pMrProt, SeqLim* pSeqLim, const char* ptModule)
{
    // **************************************************************************
    // * Name        : printImports                                             *
    // *                                                                        *
    // * Class       : SBBGREKernel                                             *
    // *                                                                        *
    // * Description : Print all imported parameters                            *
    // *                                                                        *
    // * Return      : bool                                                     *
    // *                                                                        *
    // **************************************************************************

    long   lI;    // * Loop counter *
  
  


    if ( mIsTRinp )  {
        TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Imported Parameter                                *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);

        for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
            TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_adROAsymmetryBefore[lI][%1ld] = %11.6f   ", ptModule, lI, m_adROAsymmetryBefore[lI]);
            TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_adROAsymmetryAfter[lI][%1ld]  = %11.6f   ", ptModule, lI, m_adROAsymmetryAfter[lI]);
        }

        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_FirstLineToMeasure      = %8ld \n", ptModule, m_lFirstLineToMeasure  );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lLastLineToMeasure      = %8ld \n", ptModule, m_lLastLineToMeasure   );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lFirstPartToMeasure     = %8ld \n", ptModule, m_lFirstPartToMeasure  );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lLastPartToMeasure      = %8ld \n", ptModule, m_lLastPartToMeasure   );
        if (m_bBalancedGradMomentSliceSelect)  {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_bBalancedGradMomentSliceSelect = balanced \n", ptModule);
        }  else {
            TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_bBalancedGradMomentSliceSelect = not balanced \n", ptModule);
        }

        switch ( m_eReadOutGradType )  {
            case Spoiler:            TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Spoiler \n"            , ptModule);      break;
            case Constant:           TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Constant\n"            , ptModule);      break;
            case Rewinder:           TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Rewinder\n"            , ptModule);      break;
            case Symmetric_Rewinder: TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Symmetric_Rewinder\n"  , ptModule);      break;
            default:                 TRACE_PUT1(TC_INFO, TF_SEQ, "%s: Unknown read out gradient type \n"                , ptModule);
        }

        if ( m_pRF_Exc )  {
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_Exc->getDuration()  = %8ld \n", ptModule, m_pRF_Exc->getDuration());
        } else {
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "Invalid pointer detected" ) ) return ( false );
        }

        for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  { 
            TRACE_PUT3(TC_INFO, TF_SEQ, "%s MrProt->te()[%1ld]           = %8ld [us]\n", ptModule, lI, pMrProt->te()[lI]  );
        }
        TRACE_PUT2(TC_INFO, TF_SEQ,   "%s MrProt->tr()[0]           = %8ld [us]\n", ptModule, pMrProt->tr()[0]  );
        TRACE_PUT2(TC_INFO, TF_SEQ,   "%s m_lTRFill                 = %8ld [us]\n", ptModule, m_lTRFill         );
        TRACE_PUT1(TC_INFO, TF_SEQ,   "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    }
  
    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Print all results                                 *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBF0297
bool SBBGREKernel::printResults (SeqLim* pSeqLim, const char* ptModule)
{

    // **************************************************************************
    // * Name        : printResults                                             *
    // *                                                                        *
    // * Class       : SBBGREKernel                                             *
    // *                                                                        *
    // * Description : Print all results of the kernel calculation              *
    // *                                                                        *
    // * Return      : bool                                                     *
    // *                                                                        *
    // **************************************************************************

    long   lI;    // * Loop counter *
  
  
    if ( mIsTRres )  {
        TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Exported Parameter                                *\n" ,ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lExcitationAllTime          = %6ld [us]\n"        , ptModule, m_lExcitationAllTime         );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseTime               = %6ld [us]\n"        , ptModule, m_lPrephaseTime              );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lRephaseTime                = %6ld [us]\n"        , ptModule, m_lRephaseTime               );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lOverTime                   = %6ld [us]\n"        , ptModule, m_lOverTime                  );
        for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
            TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_alTEMin[%1ld]                  = %6ld [us]\n"   , ptModule, lI, m_alTEMin[lI]            );
        }
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lTRMin                      = %6ld [us]\n"        , ptModule, m_lTRMin                     );
        for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
            TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_alTEFill[%1ld]                 = %6ld [us]\n\n" , ptModule, lI, m_alTEFill[lI]           );
        }
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lStartTimeWakeUp            = %6ld [us]\n"        , ptModule, m_lStartTimeWakeUp           );

        if ( m_pRF_WE )  {
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_WE->getGSData().getdFirstAmplitude()   = %11.6f[mT/m]\n"    , ptModule, m_pRF_WE->getGSData().getdFirstAmplitude()   );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_WE->getGSData().getdRequiredRefocusingMoment() = %11.6f[mT/m]\n\n"  , ptModule, m_pRF_WE->getGSData().getdRequiredRefocusingMoment() );
        } else {
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "Invalid pointer detected" ) ) return ( false );
        }


        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseFillTimeSlice      = %6ld [us]\n"      , ptModule, m_lPrephaseFillTimeSlice     );
        if ( m_eLimitingPrephaseTime == Slice ) {
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampUpTime()      = %6ld [us] *\n"      , ptModule, m_sTB_3D.getRampUpTime()     );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getDuration()        = %6ld [us] *\n"      , ptModule, m_sTB_3D.getDuration()       );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampDownTime()    = %6ld [us] *\n"      , ptModule, m_sTB_3D.getRampDownTime()   );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dMomentum_TB_3D             = %11.6f[mT/m*us]\n\n", ptModule, m_dMomentum_TB_3D            );
        } else {
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampUpTime()      = %6ld [us]\n"        , ptModule, m_sTB_3D.getRampUpTime()     );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getDuration()        = %6ld [us]\n"        , ptModule, m_sTB_3D.getDuration()       );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampDownTime()    = %6ld [us]\n"        , ptModule, m_sTB_3D.getRampDownTime()   );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dMomentum_TB_3D             = %11.6f[mT/m*us]\n\n", ptModule, m_dMomentum_TB_3D            );
        }

        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3DR.getRampUpTime()     = %6ld [us]\n"        , ptModule, m_sTB_3DR.getRampUpTime()    );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3DR.getDuration()       = %6ld [us]\n"        , ptModule, m_sTB_3DR.getDuration()      );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3DR.getRampDownTime()   = %6ld [us]\n\n"      , ptModule, m_sTB_3DR.getRampDownTime()  );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dMomentum_TB_3DR            = %11.6f[mT/m*us]\n\n", ptModule, m_dMomentum_TB_3DR           );
    
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseFillTimePhase      = %6ld [us]\n"      , ptModule, m_lPrephaseFillTimePhase     );
//        if ( m_eLimitingPrephaseTime == Phase ) {
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampUpTime()      = %6ld [us] *\n"  , ptModule, m_sTB_PE.getRampUpTime()     );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getDuration()        = %6ld [us] *\n"  , ptModule, m_sTB_PE.getDuration()       );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampDownTime()    = %6ld [us] *\n\n", ptModule, m_sTB_PE.getRampDownTime()   );
//        } else {
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampUpTime()      = %6ld [us]\n"    , ptModule, m_sTB_PE.getRampUpTime()     );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getDuration()        = %6ld [us]\n"    , ptModule, m_sTB_PE.getDuration()       );
//            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampDownTime()    = %6ld [us]\n\n"  , ptModule, m_sTB_PE.getRampDownTime()   );
//        }

        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PER.getRampUpTime()     = %6ld [us]\n"        , ptModule, m_sTB_PER.getRampUpTime()    );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PER.getDuration()       = %6ld [us]\n"        , ptModule, m_sTB_PER.getDuration()      );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PER.getRampDownTime()   = %6ld [us]\n\n"      , ptModule, m_sTB_PER.getRampDownTime()  );
    
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseFillTimeRead       = %6ld [us]\n"      , ptModule, m_lPrephaseFillTimeRead      );
        if ( m_eLimitingPrephaseTime == Read ) {
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampUpTime()      = %6ld [us] *\n"      , ptModule, m_sP_ROP.getRampUpTime()     );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getDuration()        = %6ld [us] *\n"      , ptModule, m_sP_ROP.getDuration()       );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampDownTime()    = %6ld [us] *\n"      , ptModule, m_sP_ROP.getRampDownTime()   );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getAmplitude()       = %11.6f[mT/m]\n\n" , ptModule, m_sP_ROP.getAmplitude()      );

            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampUpTime()   = %6ld [us] *\n"     , ptModule, m_sP_ROP_FC.getRampUpTime()  );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getDuration()     = %6ld [us] *\n"     , ptModule, m_sP_ROP_FC.getDuration()    );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampDownTime() = %6ld [us] *\n"     , ptModule, m_sP_ROP_FC.getRampDownTime());
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getAmplitude()    = %11.6f[mT/m] *\n\n", ptModule, m_sP_ROP_FC.getAmplitude()   );
        } else {
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampUpTime()      = %6ld [us]\n"        , ptModule, m_sP_ROP.getRampUpTime()     );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getDuration()        = %6ld [us]\n"        , ptModule, m_sP_ROP.getDuration()       );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampDownTime()    = %6ld [us]\n"        , ptModule, m_sP_ROP.getRampDownTime()   );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getAmplitude()       = %11.6f[mT/m]\n\n"   , ptModule, m_sP_ROP.getAmplitude()      );

            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampUpTime()   = %6ld [us]\n"      , ptModule, m_sP_ROP_FC.getRampUpTime()  );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getDuration()     = %6ld [us]\n"      , ptModule, m_sP_ROP_FC.getDuration()    );
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampDownTime() = %6ld [us]\n"      , ptModule, m_sP_ROP_FC.getRampDownTime());
            TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getAmplitude()    = %11.6f[mT/m]\n\n" , ptModule, m_sP_ROP_FC.getAmplitude()   );
        }


        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getRampUpTime()      = %6ld [us]\n"       , ptModule, m_sP_ROD.getRampUpTime()      );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getDuration()        = %6ld [us]\n"       , ptModule, m_sP_ROD.getDuration()        );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getRampDownTime()    = %6ld [us]\n"       , ptModule, m_sP_ROD.getRampDownTime()    );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getAmplitude()       = %11.6f[mT/m]\n\n"  , ptModule, m_sP_ROD.getAmplitude()       );

        TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    }  

    return ( true );

}

//	The water excitation pulse may be used with an offset 
//frequency. This
//	function enables disables the use of the offset frequency 
//and specifies the
//	frequency in Hz. The offset frequency is an optional 
//parameter. Its default
//	value 0.0 Hz.
//##ModelId=3D352EBF0247
void SBBGREKernel::useWEOffsetFrequency (bool bUseOffsetFrequency, double dWEOffsetFrequency)
{
    m_bUseWEOffsetFrequency = bUseOffsetFrequency;
    
    if ( m_bUseWEOffsetFrequency ) {
    	m_dWEOffsetFrequency = dWEOffsetFrequency;
    } else {
    	m_dWEOffsetFrequency = 0.0;
    }


}


bool SBBGREKernel::calculateInit (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{


    resetTimingCalculated();

    m_lMinDelayBetweenADCs     = static_cast<long>(SysProperties::getMinDurationBetweenReadoutAndReadout(pMrProt->rxSpec().realDwellTime()[0] / 1000.0));
    m_lMinDelayBetweenADCAndRF = SysProperties::getMinDurationBetweenReadoutAndRFPulse();
    m_lMinDelayBetweenRFAndADC = SysProperties::getMinDurationBetweenRFPulseAndReadout();


    // * ---------------------------------------------------------------------- *
    // * Check whether the required number of contrasts exceeds the maximum     *
    // * number that can be handled by the kernel.                              *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->contrasts() > lMaxNumberOfContrasts )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            if ( setNLSStatus(SEQU_ERROR, ptModule, "Maximum number of contrasts exceeded. \n") ) return ( false );
        } else  {
            if ( setNLSStatus(SEQU_ERROR ) ) return ( false );
        }
    } else {
        m_lNumberOfEchoes = pMrProt->contrasts();
    }

    m_lReadoutTime = 0;


    // * ---------------------------------------------------------------------- *
    // * Initialize readout SBBs                                                *
    // * ---------------------------------------------------------------------- *
    long lI;
    for ( lI = 0; lI < m_lNumberOfEchoes; lI++ )  {
        m_aRO[lI].init(pMrProt, pSeqLim, pSeqExpo);
    }

    // * ---------------------------------------------------------------------- *
    // * reset basic timing parameters                                          *
    // * ---------------------------------------------------------------------- *
    m_lTRMin = m_lTRFill = m_lTRFillEnd = m_lPreScanFillTime = 0;

    for ( lI = 0; lI < lMaxNumberOfContrasts; lI++ )  {
        m_alTEMin[lI] = m_alTEFill[lI] = 0;
    }


    return ( true );

}




bool SBBGREKernel::initForbiddenTRHandling (long* alHarmOsc, long lNoOfHarmOsc = 1)
{
    return ( m_ForbiddenTR.init(alHarmOsc, lNoOfHarmOsc) );

}




long SBBGREKernel::getAllowedTR_us (long lTR_us) const
{

    long        lAllowedTR_us = lTR_us;
    double      dMinForbiddenTR, dMaxForbiddenTR;
    unsigned    lI = 0;

    while ( lI < m_ForbiddenTR.getSize() )  {

        dMinForbiddenTR =  m_ForbiddenTR.getInterval()[lI].getMin() * 1000.0;
        dMaxForbiddenTR =  m_ForbiddenTR.getInterval()[lI].getMax() * 1000.0;

        dMinForbiddenTR = fSDSRoundDownGRT (dMinForbiddenTR);
        dMaxForbiddenTR = fSDSRoundUpGRT   (dMaxForbiddenTR);

        lI++;

        if ( (lAllowedTR_us > dMinForbiddenTR) && (lAllowedTR_us < dMaxForbiddenTR) ) {
            lAllowedTR_us  = static_cast<long>( dMaxForbiddenTR );
        }

    }

    return ( fSDSRoundUpGRT(lAllowedTR_us) );

}



MdhProxy& SBBGREKernel::Mdh (long lContrast)
{
    return ( m_aRO[lContrast].Mdh() );
}



sREADOUT& SBBGREKernel::ADC(long lContrast)
{
    return ( m_aRO[lContrast].ADC() );
}




//long SBBGREKernel::RoundUpGRT (long lDuration)
//{
//    return ( ( (lDuration + eGradRasterTime - 1 ) / eGradRasterTime) * eGradRasterTime );
//}


long SBBGREKernel::RoundDownGRT (long lDuration)
{
    return ( (lDuration / eGradRasterTime) * eGradRasterTime );
}



double SBBGREKernel::fabs (double dValue)
{
    return ( dValue < 0.0 ?  -1.0 * dValue : dValue );
}


double SBBGREKernel::getadReadOutAsymUI (long lIndex) const
{
    return (m_aRO[lIndex].getdReadOutAsymUI());
}

// Sets the flip-angle array in the GRE Kernel with the given values.  Checks first of all if the input variables are meaningful
// The Kernel is set from the sequence.  First this function checks that the input paramters are emaningful.  If it is then it copies the pointer across
bool SBBGREKernel::setFlipAngleArray (double * padFlipAngleArray, unsigned int uiNmbrOfFlipAngles)
{
    if ((padFlipAngleArray == NULL) || (uiNmbrOfFlipAngles<1))  {
        TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: FlipAngleArray passed is NULL or size of array <1", ptModule);
        return (false);
    } else {
        m_padFlipAngleArray = padFlipAngleArray;
        m_uiNmbrOfFlipAngles = uiNmbrOfFlipAngles;
        return (true);
    }
}



double SBBGREKernel::getdPartialFourierRO (long lIndex) const
{

    if ( (-1 < lIndex) && (lIndex < lMaxNumberOfContrasts) )  {
        return (m_aRO[lIndex].getdPartialFourierRead());
    } else {
        return (-1.0);
    }

}



void SBBGREKernel::getdReadOutAsym (long lIndex, double & dAsymBefore, double & dAsymAfter) const
{
    if ( (-1 < lIndex) && (lIndex < lMaxNumberOfContrasts) )  {
        m_aRO[lIndex].getdReadOutAsym(dAsymBefore, dAsymAfter);
    } else {
        dAsymBefore = -1.0;
        dAsymAfter  = -1.0;
    }
}





void SBBGREKernel::setbEchoOnGRT (bool bFlag)
{
    m_bEchoOnGRT = bFlag;
}


long SBBGREKernel::getTimeToRFCenter_us() const
{
    if(m_pRF_WE) return m_pRF_WE->getRequiredLeadTime()+m_pRF_WE->getlRFCenterTime();
    return 0;
}

long SBBGREKernel::setRampTimeofPreviousEB(long lVal_us)
{
    return m_lRampTimePreviousEB_us = lVal_us;
}
long SBBGREKernel::getRampTimeofPreviousEB() const
{
    return m_lRampTimePreviousEB_us;
}

void SBBGREKernel::adjUpdate(int iAdjUpdate)
{
    m_iAdjUpdate = iAdjUpdate;
}


// * -------------------------------------------------------------------------- *
// * -------------------------------------------------------------------------- *
// *                                                                            *
// *                              OLD KERNEL DEFINITION                         *
// *                                                                            *
// * -------------------------------------------------------------------------- *
// * -------------------------------------------------------------------------- *

//	*****************************************************
//	*                                                   *
//	* Constructor                                       *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC1004C
SeqBuildBlockGREKernel::SeqBuildBlockGREKernel (SBBList* pSBBList)
      : m_RF_WE121(NULL,FALSE,180.0,3,"Excit"),
        m_RF_WE121_short(NULL,FALSE,90.0,3,"Excit"),
        m_RF_WE11_short(NULL, FALSE, 45.0,2,"Excit"),
        m_RF_WE11_90(NULL, FALSE, 90.0,2,"Excit"),
        m_RF_WE11_medium(NULL, FALSE, 75.0,2,"Excit"),
        m_RF_WE11(NULL,FALSE, 180.0,2,"Excit"),
        m_pRF_WE(NULL),
        m_pRF_Exc(NULL),
        m_pRF_Pre(NULL),
        m_lTimeToFirstPreScanRF(0),
        m_lPreScanFillTime(0),
        m_lFirstLineToMeasure(0),
        m_lLastLineToMeasure(0),
        m_lFirstPartToMeasure(0),
        m_lLastPartToMeasure(0),
        m_eReadOutGradType(Spoiler),
        m_bBalancedGradMomentSliceSelect(false),
        m_bSendOscBit(false),
        m_lExcitationAllTime(0),
        m_lTEFromExcitation(0),
        m_lPrephaseTime(0),
        m_lADCDatTime(NULL),
        m_lRephaseTime(0),
        m_lTimeToFirstRF(0),
        m_lSSFlatBeforeRF(0),
        m_dMomentum_TB_3D(0.0),
        m_dMomentum_TB_3DR(0.0),
        m_dMomentROP(0.0),
        m_dMomentRO(0.0),
        m_dRelRODephasingMoment(0.5),
        m_dRelSSDephasingMoment(0.0),
        m_lTimeToADC(NULL),
//        m_lRODuration(0),
        m_lADCZeroBefore(NULL),
        m_lADCZeroAfter(NULL),
        m_lROFlatToEcho(NULL),
        m_adRelWinWidthCol(NULL),
        m_adRelWinAsymCol(NULL),
        m_adRelCenterCol(NULL),
        m_dFactorForPixelSizeRO(1.0),
        m_dFactorForUpperPixelSizeRO(1.0),
        m_dFactorForPixelSizePE(1.0),
        m_dFactorForPixelSize3D(1.0),
        m_dFactorForSliceThickness(1.0),
        m_bUseRewinder(true),
        m_bRampsOutsideOfEventBlock(true),
        m_bPrepareTrueFispPreScan(false),
        m_bSendWakeUpBit(false),
        m_lStartTimeWakeUp(500),
        m_alTEMin(NULL),
        m_lTRMin(0),
        m_alTEFill(NULL),
        m_lTRFill(0),
        m_lTRFillEnd(0),
        m_bPerformNegativeTEFillCheck(true),
        m_lLineNumber(0),
        m_lKSpaceCenterLine(0),
        m_lPartNumber(0),
        m_lKSpaceCenterPartition(0),
        m_dRFSpoilPhase(0),
        m_dWEBandwidthTimeProduct(2.0),
        m_ulUTIdent(0),
        m_bCalculated(false),
        m_bsuccessfulAllocation(true),
        m_lTrueFispSSDephasingDuration(0),
        m_lNumberOfEchoes(1),
        m_eWEMode(_1_2_1_180_Deg),
        m_bUpdateRotationMatrix(false),
        m_eLimitingPrephaseTime(Read),
        m_eLimitingRephaseTime(Read),
        m_lDelayAdjustADCPosition(0),
        m_lDelayAdjustSSPosition(0),
        m_bUseWEOffsetFrequency(false),
        m_dWEOffsetFrequency(0.0),
        m_lPhaseEventDelay(0),
        m_bLastScanInConcat(false),
        m_bLastScanInMeas(false),
        m_bFlowCompRead(false),
        m_bFlowCompSlice(false),
        m_lPrephaseFillTimeRead(0),
        m_lPrephaseFillTimeSlice(0),
        m_lPrephaseFillTimePhase(0),
        m_lOverTime(0),
        m_PhaseStabScan(NULL)
,m_lMinDelayBetweenADCs(0)
,m_lMinDelayBetweenADCAndRF(0)
,m_lMinDelayBetweenRFAndADC(0)
  ,m_lDiffBetweenROAndADCTiming(0)
  ,m_lReadoutTime(0)
  ,m_eROPolarity(Positive)
  ,m_bSuppressTraces(false)
  , SeqBuildBlock ( pSBBList )
{
  
  long lI;    // * Loop counter *

  // * This is no longer a poor SBB *
  setIdent ("SBBKernel");


  // * Memory allocation for multiple contrasts *
  m_bsuccessfulAllocation = true;

  if ( (m_sP_RO            = new sGRAD_PULSE_RO [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;

  if ( (m_alTEMin          = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_alTEFill         = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_lADCDatTime      = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_lTimeToADC       = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_lADCZeroBefore   = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_lADCZeroAfter    = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_lROFlatToEcho    = new long           [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_adRelWinWidthCol = new double         [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_adRelWinAsymCol  = new double         [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;
  if ( (m_adRelCenterCol   = new double         [lMaxNumberOfContrasts]) == NULL ) m_bsuccessfulAllocation = false;








  // * ---------------------------------------------------------------------- *
  // * Initialization of arrays                                               *
  // * ---------------------------------------------------------------------- *
  if ( m_bsuccessfulAllocation ) {
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_alTEMin[lI]          = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_alTEFill[lI]         = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_lADCDatTime[lI]      = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_lTimeToADC[lI]       = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_lADCZeroBefore[lI]   = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_lADCZeroAfter[lI]    = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_lROFlatToEcho[lI]    = 0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_adRelWinWidthCol[lI] = 0.0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_adRelWinAsymCol[lI]  = 0.0;
    for (lI=0; lI<lMaxNumberOfContrasts; lI++  ) m_adRelCenterCol[lI]   = 0.0;
  }




  m_sTB_3D.setIdent ("3D phase encoding gradient");
  m_sTB_3DR.setIdent("3D phase encoding rewinder");
  m_sTB_PE.setIdent ("Phase encoding gradient");
  m_sTB_PER.setIdent("Phase encoding rewinder");
  m_sP_ROP.setIdent ("Read out prephasing gradient");
  for (lI=0; lI<lMaxNumberOfContrasts; lI++)  {
    m_sP_RO[lI].setIdent (RTEIDENT_PulseRO);
    m_ADC[lI].setIdent ("ADC");
  }
  m_sP_ROD.setIdent ("Read out dephasing gradient");

}


//##ModelId=3D352EC10041
SeqBuildBlockGREKernel::~SeqBuildBlockGREKernel()
{
  
  delete [] m_sP_RO;
  
  delete [] m_alTEMin;
  delete [] m_alTEFill;
  delete [] m_lADCDatTime;
  delete [] m_lTimeToADC;
  delete [] m_lADCZeroBefore;
  delete [] m_lADCZeroAfter;
  delete [] m_lROFlatToEcho;
  delete [] m_adRelWinWidthCol;
  delete [] m_adRelWinAsymCol;
  delete [] m_adRelCenterCol;

}



//	*****************************************************
//	*                                                   *
//	* Sets a pointer to a sRF_PULSE instance            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC003B1
bool SeqBuildBlockGREKernel::setRFPulseForExcitation (sRF_PULSE* pSRF)
{
  // **********************************************************************************
  // * Name        : setRFPulseForExcitation                                          *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Initialization of the member m_pRF_Exc.                          *
  // *               m_pRF_Exc is an instance of sRF_PULSE.                           *
  // *                                                                                *
  // * Return      : true in case of successful initialization else false             *
  // *                                                                                *
  // **********************************************************************************

    static const char *ptModule = "SeqBuildBlockGREKernel::setRFPulseForExcitation";


    m_bCalculated = false;
    resetPrepared();    // SeqBuildBlockGREKernel is no longer prepared due the selection of a new rf-pulse

    if ( pSRF == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF == NULL");
        return ( false );
    }

    if ( pSRF->isTypeExcitation() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF->isTypeExcitation() == false");
        return ( false );
    }

    if ( pSRF->isPrepared() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF->isPrepared()==false");
        return ( false );
    }

    m_pRF_Exc = pSRF;

    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Sets a pointer to a sRF_PULSE instance of the     *
//	* excitation pulses of the gradient echo kernel and *
//	* the prescan (TrueFisp)                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC003D9
bool SeqBuildBlockGREKernel::setRFPulseForExcitation (sRF_PULSE* pSRF_Exc, sRF_PULSE* pSRF_Pre)
{
  // **********************************************************************************
  // * Name        : setRFPulseForExcitation                                          *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Initialization of the member m_pRF_Exc and m_pRF_Pre.            *
  // *               m_pRF_Exc and m_pRF_Pre are instances of sRF_PULSE.              *
  // *               This function is used for TrueFisp sequences.                    *
  // *                                                                                *
  // * Return      : true in case of successful initialization else false             *
  // *                                                                                *
  // **********************************************************************************


    static const char *ptModule = "SeqBuildBlockGREKernel::setRFPulseForExcitation";


    m_bCalculated = false;
    resetPrepared();    // SeqBuildBlockGREKernel is no longer prepared due the selection of a new rf-pulse



    // * -------------------------------------------------------------------------- *
    // * Check pSRF_Exc                                                             *
    // * -------------------------------------------------------------------------- *
    if ( pSRF_Exc == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Exc == NULL");
        return ( false );
    }



    if ( pSRF_Exc->isTypeExcitation() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Exc->isTypeExcitation() == false");
        return ( false );
    }



    if ( pSRF_Exc->isPrepared() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Exc->isPrepared()==false");
        return ( false );
    }




    // * -------------------------------------------------------------------------- *
    // * Check pSRF_Pre                                                             *
    // * -------------------------------------------------------------------------- *
    if ( pSRF_Pre == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Pre == NULL");
        return ( false );
    }



    if ( pSRF_Pre->isTypeExcitation() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Pre->isTypeExcitation() == false");
        return ( false );
    }



    if ( pSRF_Pre->isPrepared() == false )  {
        setNLSStatus(SEQU_ERROR, ptModule, "pSRF_Pre->isPrepared()==false");
        return ( false );
    }




    m_pRF_Exc = pSRF_Exc;
    m_pRF_Pre = pSRF_Pre;

    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Sets the gradient performance due to the protocol *
//	* settings (FAST, NORMAL, WHISPER).                 *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC00389
bool SeqBuildBlockGREKernel::updateGradPerf (enum SEQ::Gradients eGradientMode)
{
  // **********************************************************************************
  // * Name        : updateGradientPerformance                                        *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Sets the gradient mode (FAST, NORMAL or WHISPER)                 *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

    static const char *ptModule = "SeqBuildBlockGREKernel::updateGradPerf";


    // * -------------------------------------------------------------------------- *
    // * Check pointer m_pRF_WE                                                     *
    // * -------------------------------------------------------------------------- *
    if ( m_pRF_WE == NULL )  {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL");
        return ( false );
    }



    // * -------------------------------------------------------------------------- *
    // * Set gradient performance for slice select gradient                         *
    // * -------------------------------------------------------------------------- *
    double adSSExternalPulseMinRiseTime[3]  = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_SS),
                                                getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_SS),
                                                getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_SS) };

    double adSSExternalPulseMaxMagnitude[3] = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_SS),
                                                getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_SS),
                                                getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_SS) };

    double adSSWEPulseMinRiseTime[3]        = { getMinRiseTime  (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_WE),
                                                getMinRiseTime  (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_WE),
                                                getMinRiseTime  (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_WE) };


    double adSSWEPulseMaxMagnitude[3]       = { getMaxMagnitude (SEQ::GRAD_FAST   , SBBGREKernel_GRAD_GROUP_WE),
                                                getMaxMagnitude (SEQ::GRAD_NORMAL , SBBGREKernel_GRAD_GROUP_WE),
                                                getMaxMagnitude (SEQ::GRAD_WHISPER, SBBGREKernel_GRAD_GROUP_WE) };


    m_pRF_WE->setMinRiseTimes  ( adSSExternalPulseMinRiseTime , SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF );
    m_pRF_WE->setMaxMagnitudes ( adSSExternalPulseMaxMagnitude, SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF );
    m_pRF_WE->setMinRiseTimes  ( adSSWEPulseMinRiseTime       , SBBBinomialPulses_GRAD_PERF_BINOMIAL    );
    m_pRF_WE->setMaxMagnitudes ( adSSWEPulseMaxMagnitude      , SBBBinomialPulses_GRAD_PERF_BINOMIAL    );

    m_pRF_WE->setMinRiseTimeGSWD  (getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_SS), SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    m_pRF_WE->setMaxMagnitudeGSWD (getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_SS), SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
    m_pRF_WE->setMinRiseTimeGSWD  (getMinRiseTime (SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_WE), SBBBinomialPulses_GRAD_PERF_BINOMIAL);
    m_pRF_WE->setMaxMagnitudeGSWD (getMaxMagnitude(SEQ::GRAD_FAST_GSWD_RISETIME,SBBGREKernel_GRAD_GROUP_WE), SBBBinomialPulses_GRAD_PERF_BINOMIAL);



    // * -------------------------------------------------------------------------- *
    // * Set gradient performance for slice refocussing gradients                   *
    // * -------------------------------------------------------------------------- *
    m_sP_3D_FC.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_SS) );
    m_sP_3D_FC.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_SS) );



    // * -------------------------------------------------------------------------- *
    // * Set gradient performance for phase encoding gradient                       *
    // * -------------------------------------------------------------------------- *
    m_sTB_PE .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_PE) );
    m_sTB_PE .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_PE) );

    m_sTB_PER.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_PE) );
    m_sTB_PER.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_PE) );



    // * -------------------------------------------------------------------------- *
    // * Set gradient performance for phase encoding gradient                       *
    // * -------------------------------------------------------------------------- *
    m_sTB_3D .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_3D) );
    m_sTB_3D .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_3D) );

    m_sTB_3DR.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_3D) );
    m_sTB_3DR.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_3D) );



    // * -------------------------------------------------------------------------- *
    // * Set gradient performance for read out gradient                             *
    // * -------------------------------------------------------------------------- *
    m_sP_ROP   .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );
    m_sP_ROP   .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );

    m_sP_ROP_FC.setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );
    m_sP_ROP_FC.setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROP) );

    m_sP_RO[0] .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_RO)  );
    m_sP_RO[0] .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_RO)  );

    m_sP_RO[1] .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_RO)  );
    m_sP_RO[1] .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_RO)  );

    m_sP_ROD   .setMinRiseTime  ( getMinRiseTime  ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROD) );
    m_sP_ROD   .setMaxMagnitude ( getMaxMagnitude ( eGradientMode, SBBGREKernel_GRAD_GROUP_ROD) );



    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the GRE kernel                     *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC002FC
bool SeqBuildBlockGREKernel::calculateTiming (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
  // **********************************************************************************
  // * Name        : getcalculateTiming                                               *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the timing of the GRE kernel                      *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

    static const char *ptModule = "SeqBuildBlockGREKernel::calculateTiming";

    static GPAProxy   theGPA;                   // * Information about the gradient switching time            *

    long   lI;                                  // * Loop counter                                             *
    double dROMomentBeforePhaseStab   = 0.0;


    m_bCalculated = false;
    m_lMinDelayBetweenADCs = SysProperties::getMinDurationBetweenReadoutAndReadout(pMrProt->rxSpec().realDwellTime()[0] / 1000.0) /*static_cast<long>(theRX4Proxy.getMinDurationBetweenReadoutAndReadout(pMrProt->rxSpec().realDwellTime()[0] / 1000.0))*/;


    if ( mIsTRrun )  {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s running\n", ptModule);
    }


    // * -------------------------------------------------------------------------- *
    // * Debug: running message and imported parameters                             *
    // * -------------------------------------------------------------------------- *
    printImports ( pMrProt, pSeqLim, ptModule );



    // * -------------------------------------------------------------------------- *
    // * Initialize variables                                                       *
    // * -------------------------------------------------------------------------- *
    calculateInit (pMrProt, pSeqLim, pSeqExpo);



    
    // * -------------------------------------------------------------------------- *
    // * Check kernel settings                                                      *
    // * -------------------------------------------------------------------------- *
    if (! calculateCheckSetting( pMrProt, pSeqLim ) ) return ( false );



    if (    pMrProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_YES
         || pMrProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_READOUT_ONLY ) {
        m_bFlowCompRead = true;
    } else {
        m_bFlowCompRead = false;
    }
  

    if (    pMrProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_YES
         || pMrProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_SLICESEL_ONLY ) {
        m_bFlowCompSlice = true;
    } else {
        m_bFlowCompSlice = false;
    }




    // * -------------------------------------------------------------------------- *
    // * Get scaling factors for gradient timing calculation                        *
    // * -------------------------------------------------------------------------- *
    m_dFactorForPixelSizeRO      = m_pCalcLimits->getFactorForPixelSizeRO      ( pMrProt );
    m_dFactorForUpperPixelSizeRO = m_pCalcLimits->getFactorForUpperPixelSizeRO ( pMrProt );
    m_dFactorForPixelSizePE      = m_pCalcLimits->getFactorForPixelSizePE      ( pMrProt );
    if ( pMrProt->kSpace().dimension() == SEQ::DIM_3 )  {
        m_dFactorForPixelSize3D    = m_pCalcLimits->getFactorForPixelSize3D      ( pMrProt );
    }
    m_dFactorForSliceThickness   = m_pCalcLimits->getFactorForSliceThickness   ( pMrProt );




    // *--------------------------------------------------------------------------- *
    // * Prepares the ADC                                                           *
    // *--------------------------------------------------------------------------- *
    calculatePrepADC ( pMrProt, pSeqLim );




    // *--------------------------------------------------------------------------- *
    // * Selects the excitation pulse                                               *
    // *--------------------------------------------------------------------------- *
    if (! calculateSelectWEPulse(pMrProt, pSeqLim, pSeqExpo) ) return ( false );




    // * -------------------------------------------------------------------------- *
    // * Update of the gradient performance settings (FAST, NORMAL, WHISPER)        *
    // *                                                                            *
    // * CAUTION: Function updateGradientPerformance MUST be called after           *
    // *          function calculateSelectWEPulse to initilize m_pRF_WE.            *
    // * -------------------------------------------------------------------------- *
    if (! updateGradPerf ( pMrProt->gradSpec().mode() ) ) return ( false );



    // *--------------------------------------------------------------------------- *
    // * Prepares the excitation pulse                                              *
    // *--------------------------------------------------------------------------- *
    if (! calculatePrepWEPulse (pMrProt, pSeqLim, pSeqExpo) ) return ( false );





    // ******************************************************************************
    // ******************************************************************************
    // *                                                                            *
    // * Calculate the gradient timing for phase encoding, kz                       *
    // *                                                                            *
    // ******************************************************************************
    // ******************************************************************************
    if (! calculate3DGradientTables ( pMrProt, pSeqLim ) ) return ( false );






    // ******************************************************************************
    // ******************************************************************************
    // *                                                                            *
    // * Calculate gradient timing for phase encoding, ky                           *
    // *                                                                            *
    // ******************************************************************************
    // ******************************************************************************
    if (! calculatePEGradientTables ( pMrProt ) ) return ( false );






    // ******************************************************************************
    // ******************************************************************************
    // *                                                                            *
    // * Calculation of the read out timing                                         *
    // *                                                                            *
    // ******************************************************************************
    // ******************************************************************************

    // * -------------------------------------------------------------------------- *
    // * Required area under read out gradient flat top [mT/m * us]                 *
    // * -------------------------------------------------------------------------- *
    m_dMomentRO = fGSLGetROMoment ( pMrProt );




    // * -------------------------------------------------------------------------- *
    // * Amplitude of the read out gradients [mT/m]                                 *
    // * -------------------------------------------------------------------------- *
    calculateROAmplitude();
  



    // * -------------------------------------------------------------------------- *
    // * read out gradient ramp up times                                            *
    // * -------------------------------------------------------------------------- *  
    if (! calculateRORampUpTime ( pMrProt )) return ( false );




    // * -------------------------------------------------------------------------- *
    // * read out gradient ramp down times                                          *
    // * -------------------------------------------------------------------------- *
    if (! calculateRORampDownTime( pMrProt ) ) return ( false );
  
  


    // * -------------------------------------------------------------------------- *
    // * Calcution of the read out prephasing gradient                              *
    // * Flow compensation in readout direction is taken into account.              *
    // * -------------------------------------------------------------------------- *
    if (! calculateROP ( pMrProt, pSeqLim ) ) return ( false );




    // * -------------------------------------------------------------------------- *
    // * Calculation of the length of the read out gradient flat tops               *
    // * -------------------------------------------------------------------------- *
    if (! calculateRODuration (pMrProt) ) return ( false );




  // * ---------------------------------------------------------------------- *
  // * Adjust the position of the ADC for "symmetric rewinder", i.e. shift    *
  // * the center of the ADC to the center of the RO gradient.                *
  // * RO duration is a multiple of the gradient raster time in contrast to   *
  // * the ADC duration                                                       *
  // * ---------------------------------------------------------------------- *
  if ( m_eReadOutGradType == Symmetric_Rewinder ) {
      m_lDelayAdjustADCPosition = abs ( m_ADC[0].getRoundedDuration() - m_sP_RO[0].getFlatTopTime() ) / 2;
  } else {
      m_lDelayAdjustADCPosition = 0;
  }



  // **************************************************************************
  // **************************************************************************
  // *                                                                        *
  // * Which prephasing process is time limiting?                             *
  // *                                                                        *
  // * Prephasing in slice select, phase encode, read out direction           *
  // *                                                                        *
  // * or the m_lADCZeroBefore[0] duration                                    *
  // *                                                                        *
  // **************************************************************************
  // **************************************************************************
  if (! calculatePrephaseTime ( pMrProt, pSeqLim ) ) return ( false );




  // **************************************************************************
  // **************************************************************************
  // *                                                                        *
  // * Initialization of phase stabilization scans                            *
  // *                                                                        *
  // **************************************************************************
  // **************************************************************************
  if ( pMrProt->phaseStabilize() )  {
    m_PhaseStabScan.setlLastKernelRORampDownTimeInsideSBB ( m_sP_RO[m_lNumberOfEchoes-1].getRampDownTime() );
    m_PhaseStabScan.setps_First_Kernel_ADC   ( &m_ADC[0]                     );
    m_PhaseStabScan.setdFactorForPixelSizeRO ( m_dFactorForPixelSizeRO       );
    

    dROMomentBeforePhaseStab   = m_sP_ROP_FC.getMomentumTOT() + m_sP_ROP.getMomentumTOT() + m_sP_RO[0].getMomentumTOT();
    if ( pMrProt->contrasts() == 2 )  {   
      dROMomentBeforePhaseStab   +=  m_sP_RO[1].getMomentumTOT();
    }
    

    m_PhaseStabScan.setlTimeOfLastKernelADCInSBB   ( m_lADCZeroAfter[m_lNumberOfEchoes-1] );
    m_PhaseStabScan.setdROMomentBeforePhaseStabADC ( dROMomentBeforePhaseStab );


    // * ---------------------------------------------------------------------- *
    // * Calculate the RO ramp down time of the phase stabilization SBB to      *
    // * fit to ROD gradient of the kernel                                      *
    // * ---------------------------------------------------------------------- *
    switch (m_eReadOutGradType)  {
      case Spoiler:
        m_PhaseStabScan.setRORampDownTime ( fSDSRoundUpGRT( maximum ( m_PhaseStabScan.getROMaxMagnitude() * m_PhaseStabScan.getROMinRiseTime(),
                                                                  m_sP_ROD.getMaxMagnitude() * m_sP_ROD.getMinRiseTime() ) ) );
      break;

      case Constant:
      case Rewinder:
      case Symmetric_Rewinder:
        // * No RO ramp down time is specified by SBBGREKernel               *
        // * ==> SBBPhaseStabScan calculates the ramp down time by its own.  *
      break;
    
      default: 
        if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." ) )  return ( false );
    }

    
    m_PhaseStabScan.setLastRampDownOutsideSBB ( true );     // * Last gradient ramp is outside of the SBB *

  }

  if (! m_PhaseStabScan.prep ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );


  // **************************************************************************
  // **************************************************************************
  // *                                                                        *
  // * Calculation of the read out dephasing gradient                         *
  // *                                                                        *
  // **************************************************************************
  // **************************************************************************
  if (! calculateROD( pMrProt ) ) return ( false );
  



  // **************************************************************************
  // **************************************************************************
  // *                                                                        *
  // * Which process is time limiting for rephasing?                          *
  // *                                                                        *
  // * Rephasing in slice select or phase encode direction, the end of the    *
  // *                                                                        *
  // * read out dephasing gradient or the lADCZeroAfter duration ?            *
  // *                                                                        *
  // **************************************************************************
  // **************************************************************************
  if (! calculateRephaseTime ( pMrProt ) ) return ( false );




  // **************************************************************************
  // **************************************************************************
  // *                                                                        *
  // * Calculation of the duration that gradient ramps are allowed to dangle  *
  // * out of the event block                                                 *
  // *                                                                        *
  // **************************************************************************
  // **************************************************************************
  if (! calculateOverTime( pMrProt, pSeqLim ) ) return ( false );




  // * ---------------------------------------------------------------------- *
  // * Additional calculations                                                *
  // * ---------------------------------------------------------------------- *
  if (! calculateAddOn ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );




  // * ---------------------------------------------------------------------- *
  // * Calculation of the minimum TE times                                    *
  // * ---------------------------------------------------------------------- *
  if (! calculateTEMin() ) return ( false );



  // * ---------------------------------------------------------------------- *
  // * Calculate the TE fill times                                            *
  // * ---------------------------------------------------------------------- *
  if (! calculateTEFill ( pMrProt ) ) return ( false );
  



  // * ---------------------------------------------------------------------- *
  // * Calculate start time of the phase stabilization scan                   *
  // * ---------------------------------------------------------------------- *
//  lStartTimeSBBPhaseStabScan = m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getDurationPerRequest() +
//                               m_lPrephaseFillTimeRead + m_sP_ROP_FC.getTotalTime() + m_sP_ROP.getTotalTime() + 2 * m_alTEFill[0] + m_sP_RO[0].getDuration() +
//                               m_lRephaseTime;


//  if ( pMrProt->contrasts() == 2 )  {
//    lStartTimeSBBPhaseStabScan +=  m_sP_RO[0].getRampDownTime() + 2 * m_alTEFill[1] + m_sP_RO[1].getDuration();

//  }

//  m_PhaseStabScan.setlStartTime            ( lStartTimeSBBPhaseStabScan );
//  m_PhaseStabScan.setlTimeFromRFToStartSBB ( lStartTimeSBBPhaseStabScan - m_lExcitationAllTime + m_lTEFromExcitation );




  // * ---------------------------------------------------------------------- *
  // * Calculate m_lTimeToADC[1]                                              *
  // * ---------------------------------------------------------------------- *
  if (! calculateTimeToADC1 ( pMrProt ) ) return ( false );
  



  // * ---------------------------------------------------------------------- *
  // * Gradient ramps that are dangling out of the event block                *
  // * (m_lOverTime) do not contribute to the minimum TR time                 *
  // * ---------------------------------------------------------------------- *
  if (!  calculateTRMin( pMrProt ) ) return ( false );




  // * ---------------------------------------------------------------------- *
  // * Check whether m_lStartTimeWakeUp is shorter that the minimum TR time   *
  // * and correct m_lStartTimeWakeUp if necessary.                           *
  // * ---------------------------------------------------------------------- *
  if ( m_lStartTimeWakeUp > m_lTRMin )  {
    m_lStartTimeWakeUp = m_lTRMin;
  }




  // * ---------------------------------------------------------------------- *
  // * Check for negative TE fill time                                        *
  // * ---------------------------------------------------------------------- *
  if ( m_bPerformNegativeTEFillCheck )  {
    if ( (m_alTEMin[0] > pMrProt->te()[0]) )  {

      if ( !pSeqLim->isContextPrepForBinarySearch() && !m_bSuppressTraces)  {
        if ( setNLSStatus ( SEQU__NEGATIV_TEFILL, ptModule, "Negative TE fill time" ) ) return ( false );
      } else {
        if ( setNLSStatus ( SEQU__NEGATIV_TEFILL ) ) return ( false );
      }

    }
  
  
    if ( pMrProt->contrasts() == 2 )  {
      if ( ((m_alTEMin[1] + 2 * m_alTEFill[0]) > pMrProt->te()[1]) )  {

        if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
          if ( setNLSStatus ( SEQU__NEGATIV_TEFILL, ptModule, "Negative TE fill time" ) ) return ( false );
        } else {
          if ( setNLSStatus ( SEQU__NEGATIV_TEFILL ) ) return ( false );
        }

      }
    }
  }



  // *----------------------------------------------------------------------- *
  // * Debug: intermediate results                                            *
  // *----------------------------------------------------------------------- *
  if ( mIsTRint )  {
    TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Intermediate Results                              *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s theGPA.getGradSwitchTime()    = %11.6f [us]\n\n", ptModule, theGPA.getGradSwitchTime()    );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_WE->getGSData()->getlTimeToFirstFlatTop()  = %6ld [us]\n\n"  , ptModule, m_pRF_WE->getGSData().getlTimeToFirstFlatTop()  );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForPixelSizeRO         = %11.6f\n"     , ptModule, m_dFactorForPixelSizeRO       );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForUpperPixelSizeRO    = %11.6f\n"     , ptModule, m_dFactorForUpperPixelSizeRO  );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForPixelSizePE         = %11.6f\n"     , ptModule, m_dFactorForPixelSizePE       );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForPixelSize3D         = %11.6f\n"     , ptModule, m_dFactorForPixelSize3D       );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dFactorForSliceThickness      = %11.6f\n\n"   , ptModule, m_dFactorForSliceThickness    );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lTEFromExcitation             = %6ld [us]\n"  , ptModule, m_lTEFromExcitation           );

    for (lI=0; lI < m_lNumberOfEchoes; lI++)  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_ADC[%1ld].getDuration()        = %11.6f [us]\n", ptModule, lI, m_ADC[lI].getDuration()   );
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_lROFlatToEcho[%1ld]            = %6ld [us]\n"  , ptModule, lI, m_lROFlatToEcho[lI]       );
    }

    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
  }


  // *----------------------------------------------------------------------- *
  // * Debug: results                                                         *
  // *----------------------------------------------------------------------- *
  printResults ( pSeqLim, ptModule );
  


  m_bCalculated = true;




  if ( mIsTRend )  {
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s finished\n"                    ,ptModule);
  }


  return ( true );  // Returning without error

}

//	*****************************************************
//	*                                                   *
//	* Checks whether the parameter settings are valid   *
//	* to preceed in kernel calculation                  *
//	* The function sets a NLs status and returns a      *
//	* FALSE if the settings are not valid               *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBE00AB
bool SeqBuildBlockGREKernel::calculateCheckSetting (MrProt* pMrProt, SeqLim* pSeqLim)
{
  // **********************************************************************************
  // * Name        : calculateCheckSetting                                            *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Checks whether the parameter settings are valid to preceed in    *
  // *               kernel calculation.                                              *
  // *               The function sets a NLS status and returns a false if the        *
  // *               settings are not valid.                                          *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
    static const char *ptModule = "SeqBuildBlockGREKernel::calculateCheckSetting";


    // * -------------------------------------------------------------------------- *
    // * Check whether dynamic memory allocation has been successfully carried out. *
    // * Otherwise, quit the calculation of the kernel timing                       *
    // * -------------------------------------------------------------------------- *
    if ( m_bsuccessfulAllocation == false )  {
        if ( setNLSStatus(SEQU_ERROR, ptModule, "Can't allocate dynamic memory in SBBGREKernel\n") )    return ( false );
    }



    // * -------------------------------------------------------------------------- *
    // * Check whether a valid rf-pulse has been defined                            *
    // * -------------------------------------------------------------------------- *
    if ( m_pRF_Exc == NULL )  {
        setNLSStatus(SEQU_ERROR, ptModule, "No valid rf-pulse defined.\n");
        return ( false );
    }



    // * -------------------------------------------------------------------------- *
    // * Check whether a valid pre pulse has been defined for TrueFisp              *
    // * -------------------------------------------------------------------------- *
    if ( m_bPrepareTrueFispPreScan )  {
        if ( m_pRF_Pre == NULL )  {
            setNLSStatus(SEQU_ERROR, ptModule, "No valid pre rf-pulse defined for TrueFisp.\n");
            return ( false );
        }
    }



    // * -------------------------------------------------------------------------- *
    // * Check restrictions for flow comp                                           *
    // *    - only 2D imaging supported in combination with flow compensation in    *
    // *      slice direction                                                       *
    // *      only flow comp in slice or read supported                             *
    // *    - only ReadOutGradType Constant or Spoiler supported                    *
    // * -------------------------------------------------------------------------- *
    if (    (pMrProt->kSpace().dimension() == SEQ::DIM_3) 
         && (pMrProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_SLICESEL_ONLY || pMrProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_YES) )   {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            if ( setNLSStatus(SEQU_ERROR, ptModule, "Flow compensation in slice direction supported for 3D imaging. \n") ) return ( false );
        } else  {
            if ( setNLSStatus(SEQU_ERROR ) ) return ( false );
        }

    }




    // * -------------------------------------------------------------------------- *
    // * Check whether spoiler on slice select axis                                 *
    // * -------------------------------------------------------------------------- *
    if ( m_bBalancedGradMomentSliceSelect && ( m_dRelSSDephasingMoment != 0.0 ) ) {
        setNLSStatus(SEQU_ERROR, ptModule, "Slice select axis: No spoiler possible in case of a balanced gradient scheme. \n");
        return ( false );
    }

    if ( ( m_bUseRewinder == false ) && ( m_dRelSSDephasingMoment != 0.0 ) ) {
        setNLSStatus(SEQU_ERROR, ptModule, "Slice select axis: No spoiler possible if gradient pulse has been switched off.\n");
        return ( false );
    }


    if ( ( m_eReadOutGradType == Symmetric_Rewinder ) &&
         ( (m_adRelWinWidthCol[0] != 1.0) || (m_adRelWinAsymCol[0] != 0.5) || (m_adRelCenterCol[0] != 0.5) ) )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Error: No symmetric ADC timing.\n");
        return ( false );
    }


    // * -------------------------------------------------------------------------- *
    // * Check whether the number of contrast is inside the valid range             *
    // * Otherwise, quit the calculation of the kernel timing                       *
    // * -------------------------------------------------------------------------- *
    if ( pMrProt->contrasts() > lMaxNumberOfContrasts )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Number of contrasts to large\n");
        return ( false );
    }



    // * -------------------------------------------------------------------------- *
    // * Phase stabilization does not support negative readout polarity             *
    // * -------------------------------------------------------------------------- *
    if ( (m_eROPolarity != Positive) && pMrProt->phaseStabilize() )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Phase stabilization does not support negative readout polarity\n");
        return ( false );
    }


    m_lNumberOfEchoes = pMrProt->contrasts();


    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Sets the member m_pRF_WE to the selected          *
//	* binominal pulse                                   *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBE0046
bool SeqBuildBlockGREKernel::calculateSelectWEPulse (MrProt* pMrProt, SeqLim*, SeqExpo* )
{
    // **********************************************************************************
    // * Name        : calculateSetRFWEPulse                                            *
    // *                                                                                *
    // * Class       : SeqBuildBlockGREKernel                                           *
    // *                                                                                *
    // * Description : Sets the member m_pRF_WE to the selected water selective         *
    // *               rf-pulse                                                         *
    // *                                                                                *
    // * Return      : BOOL                                                             *
    // *                                                                                *
    // **********************************************************************************

    static const char *ptModule = "SeqBuildBlockGREKernel::calculateSelectWEPulse";
    



    switch ( m_eWEMode )  {
        case _1_1_45_Deg:
            m_pRF_WE = &m_RF_WE11_short;
        break;

        case _1_1_75_Deg:

            // * -------------------------------------------------------------------------- *
            // * Define parameters for sinc pulses used within _1_1_75_Deg WE pulse         *
            // * -------------------------------------------------------------------------- *
            if ( pMrProt->txSpec().excitation() == SEQ::EXCITATION_SLICE_SELECTIVE  )  {

                m_sRFWE1.setDuration             (       2000 );
                m_sRFWE1.setBandwidthTimeProduct ( 12.7 * 0.8 );
                m_sRFWE1.setSamples              (        500 );
                m_sRFWE1.setInitialPhase         (        0.0 );
                m_sRFWE1.setThickness            ( pMrProt->sliceSeries().aFront().thickness() * 0.8 );

                if(! m_RF_WE11_medium.registerBinomialRFPulse( 0, &m_sRFWE1 ) )
                {
                    setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_medium.registerBinomialRFPulse failed.");
                    return ( false );
                }

                m_sRFWE2.setDuration             ( m_sRFWE1.getDuration()             );
                m_sRFWE2.setBandwidthTimeProduct ( m_sRFWE1.getBandwidthTimeProduct() );
                m_sRFWE2.setSamples              ( m_sRFWE1.getSamples()              );
                m_sRFWE2.setInitialPhase         ( m_sRFWE1.getInitialPhase()         );
                m_sRFWE2.setThickness            ( m_sRFWE1.getThickness()            );

                if(! m_RF_WE11_medium.registerBinomialRFPulse( 1, &m_sRFWE2 ) )
                {
                    setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_medium.registerBinomialRFPulse failed.");
                    return ( false );
                }

            }
            else
            {
                m_RF_WE11_medium.unregisterBinomialRFPulses();
            }


            m_pRF_WE = &m_RF_WE11_medium;

        break;

        case _1_1_90_Deg:

            // * -------------------------------------------------------------------------- *
            // * Define parameters for sinc pulses used within _1_1_90_Deg WE pulse         *
            // * -------------------------------------------------------------------------- *
            m_sRFWE1.setDuration             (  600 );
            m_sRFWE1.setBandwidthTimeProduct (  6.0 );
            m_sRFWE1.setSamples              (  300 );
            m_sRFWE1.setInitialPhase         (  0.0 );
            m_sRFWE1.setThickness            ( pMrProt->sliceSeries().aFront().thickness() );

            if(! m_RF_WE11_90.registerBinomialRFPulse( 0, &m_sRFWE1 ) )
            {
                setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_90.registerBinomialRFPulse failed.");
                return ( false );
            }



            m_sRFWE2.setDuration             ( m_sRFWE1.getDuration()             );
            m_sRFWE2.setBandwidthTimeProduct ( m_sRFWE1.getBandwidthTimeProduct() );
            m_sRFWE2.setSamples              ( m_sRFWE1.getSamples()              );
            m_sRFWE2.setInitialPhase         ( m_sRFWE1.getInitialPhase()         );
            m_sRFWE2.setThickness            ( m_sRFWE1.getThickness()            );

            if(! m_RF_WE11_90.registerBinomialRFPulse( 1, &m_sRFWE2 ) )
            {
                setNLSStatus ( SEQU_ERROR, ptModule, "m_RF_WE11_90.registerBinomialRFPulse failed.");
                return ( false );
            }

            m_pRF_WE = &m_RF_WE11_90;

        break;

        case _1_1_180_Deg:
            m_pRF_WE = &m_RF_WE11;
        break;
        
        case _1_2_1_90_Deg:
            m_pRF_WE = &m_RF_WE121_short;
        break;

        case _1_2_1_180_Deg:
            m_pRF_WE = &m_RF_WE121;
        break;
        
        default:  
            if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown water excitation mode." ) )  return ( false );
    }



    return ( true );
    
}


bool SeqBuildBlockGREKernel::calculatePrepWEPulse (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
  //## begin SeqBuildBlockGREKernel::calculatePrepWEPulse%3CC9612D02EF.body preserve=yes


    static const char *ptModule = "SeqBuildBlockGREKernel::calculatePrepWEPulse";


    // *--------------------------------------------------------------------------- *
    // * Prepares the excitation pulse                                              *
    // *--------------------------------------------------------------------------- *
    if ( (m_pRF_WE == NULL) || (m_pRF_Exc == NULL) ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }
    
    m_pRF_WE->setUseWEOffsetFrequency     ( m_bUseWEOffsetFrequency, m_dWEOffsetFrequency );
    m_pRF_WE->setBandwidthTimeProduct     ( m_dWEBandwidthTimeProduct );
    m_pRF_WE->setNonSelectiveWERFDuration ( m_pRF_Exc->getDuration()  );
    m_pRF_WE->setLastRampDownOutsideSBB   (  true );
    m_pRF_WE->setUseOwnEventBlock         ( false );    // * Water selective pulse was used without its own event block *
    m_pRF_WE->setThickness                ( m_pRF_Exc->getThickness() );    // * Slice thickness of the WE pulse        *



    if (! m_pRF_WE->setPointerToCalculationLimits ( m_pCalcLimits ) )   {
        setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->setPointerToCalculationLimits");
        return ( false );
    }
  
    if (! m_pRF_WE->setExcitationRFPulse(m_pRF_Exc, pMrProt, pSeqExpo) )  {
        setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->setExcitationRFPulse");
        return ( false );
    }

    if (! m_pRF_WE->prep ( pMrProt, pSeqLim, pSeqExpo ) )   {
        setNLSStatus (m_pRF_WE->getNLSStatus(), "m_pRF_WE->prep");
        return ( false );
    }




    m_lExcitationAllTime = m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getDurationPerRequest();
    m_lTEFromExcitation  = m_pRF_WE->getTEContribution();



    return ( true );

  //## end SeqBuildBlockGREKernel::calculatePrepWEPulse%3CC9612D02EF.body
}




//	*****************************************************
//	*                                                   *
//	* Prepares the ADCs                                 *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBD03DE
void SeqBuildBlockGREKernel::calculatePrepADC (MrProt* pMrProt, SeqLim* pSeqLim)
{
  // **********************************************************************************
  // * Name        : calculatePrepADC                                                 *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Prepares the ADCs                                                *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  long   lI                  = 0;     // * Loop counter         *
  double dEffectiveDwellTime = 0.0;   // * ADC dwell time in us *
  
  
  for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
    m_ADC[lI].prep(pMrProt->kSpace().baseResolution(), static_cast<long>(pMrProt->rxSpec().effDwellTime( pSeqLim->getReadoutOSFactor() )[lI]));
  
  
    dEffectiveDwellTime = pMrProt->rxSpec().effDwellTime( pSeqLim->getReadoutOSFactor() )[lI] / 1000.0;

    m_lADCDatTime[lI]    = fSDSRoundUpGRT ( m_ADC[lI].getDuration() * m_adRelWinWidthCol[lI] );   

    // * Round m_lROFlatToEcho[lI] to the dwell time raster *
    m_lROFlatToEcho[lI]  = static_cast<long>( ( m_ADC[lI].getDuration() * m_adRelWinWidthCol[lI] * m_adRelWinAsymCol[lI] ) 
                                                / dEffectiveDwellTime + 0.5 );
    m_lROFlatToEcho[lI]  = static_cast<long>( m_lROFlatToEcho[lI] * dEffectiveDwellTime );


    // * Round m_lADCZeroBefore[lI] to the dwell time raster *
    m_lADCZeroBefore[lI] = static_cast<long>(max (0.0, m_ADC[lI].getDuration() * m_adRelCenterCol[lI] - m_lROFlatToEcho[lI])
                                             / dEffectiveDwellTime + 0.5);
    m_lADCZeroBefore[lI] = static_cast<long>( m_lADCZeroBefore[lI] * dEffectiveDwellTime );
    

    // * Avoid negative timing parameter due rounding of m_lADCDatTime[lI] *
    if ( m_ADC[lI].getDuration() > (m_lADCDatTime[lI] + m_lADCZeroBefore[lI]) ) {
      m_lADCZeroAfter[lI]   = max (0L, static_cast<long>(m_ADC[lI].getDuration() - m_lADCZeroBefore[lI] - m_lADCDatTime[lI] + 0.5));
    } else {
      m_lADCZeroAfter[lI]   = 0;
    }
  }
  
  // * ---------------------------------------------------------------------- *
  // * m_lPhaseEventDelay is used to prevent the last frequency/phase event   *
  // * to placed directly at the end of the event table.                      *
  // * This situation is prevented if m_lADCZeroAfter of the last ADC is zero *
  // * since there is enough time between the freq/phase event and the and of *
  // * end of event block due to the readout spoiler or phase rewinder 
  // * ---------------------------------------------------------------------- *
  if ( (m_lADCZeroAfter[m_lNumberOfEchoes-1] == 0) && 
       ((m_eReadOutGradType == Spoiler) || (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder)) ) {
      m_lPhaseEventDelay = 0;
  } else {
      m_lPhaseEventDelay = 10;
  }
}

//	*****************************************************
//	*                                                   *
//	* Calculation of the 3D gradient table and its      *
//	* rewinder for the GRE kernel                       *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBD0384
bool SeqBuildBlockGREKernel::calculate3DGradientTables (MrProt* pMrProt, SeqLim* pSeqLim)
{
  // **********************************************************************************
  // * Name        : calculate3DGradientTables                                        *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the 3D gradient table and its rewinder            *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
    static const char *ptModule = "SeqBuildBlockGREKernel::calculate3DGradientTables";


    long   lStatus             = SEQU__NORMAL;

    long   lRamp               = 0;
    long   lDuration           = 0;
  
    double dDeltaMomentumKz    = 0.0;



    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }


    // * -------------------------------------------------------------------------- *
    // * Calculation of flow compensation gradients in slice selection gradients    *
    // *                                                                            *
    // *                  !!!! ONLY 2D IMAGING SUPPORTED !!!!                       *
    // *                                                                            *
    // * -------------------------------------------------------------------------- *
    if ( m_bFlowCompSlice )   {

        sGRAD_PULSE_TRAP sP_Temp;

        sP_Temp.setMinRiseTime  ( m_sTB_3D.getMinRiseTime()  );
        sP_Temp.setMaxMagnitude ( m_sTB_3D.getMaxMagnitude() );

        if ( ! m_pRF_WE->calculateShortestFlowCompGradients(pMrProt, pSeqLim, &m_sP_3D_FC, &sP_Temp, 0.0, SEQ::AXIS_SLICE, m_dFactorForSliceThickness) )  {
            if ( ! pSeqLim->isContextPrepForBinarySearch() ) {
                setNLSStatus (m_pRF_WE->getNLSStatus(), ptModule, "m_pRF_WE->calculateShortestFlowCompGradients(...) for slice selection");
            } else {
                setNLSStatus (m_pRF_WE->getNLSStatus());
            }
            return ( false );
        }
      

        // * Outer flow compensation gradient *
        m_sTB_3D.setRampUpTime    ( sP_Temp.getRampUpTime()   );
        m_sTB_3D.setRampDownTime  ( sP_Temp.getRampDownTime() );
        m_sTB_3D.setDuration      ( sP_Temp.getDuration()     );
        m_sTB_3D.setAmplitude     ( sP_Temp.getAmplitude()    );

        m_dMomentum_TB_3D = sP_Temp.getMomentumTOT();

    } else {

        // * ---------------------------------------------------------------------- *
        // * Flow compensation gradient m_sP_3D_FC is not required.                 *
        // * ---------------------------------------------------------------------- *
        m_sP_3D_FC.setAmplitude    (0.0);
        m_sP_3D_FC.setRampUpTime   (0);
        m_sP_3D_FC.setDuration     (0);
        m_sP_3D_FC.setRampDownTime (0);



        m_dMomentum_TB_3D = m_pRF_WE->getGSData().getdRequiredRefocusingMoment();



        // * ---------------------------------------------------------------------- *
        // * Calculate the gradient timing for 3D phase encoding, kz                *
        // * ---------------------------------------------------------------------- *
        dDeltaMomentumKz     = fGSLGet3DDeltaMoment( pMrProt );

        lStatus = fGSLGetShortestTabTiming (
                    dDeltaMomentumKz * m_dFactorForPixelSize3D,       // IMP: Moment between two steps [mT/m * us]
                    SEQ::DIR_ASCENDING,                               // IMP: true <=> (-) -> (+), i.e. table direction
                    m_dMomentum_TB_3D * m_dFactorForSliceThickness,   // IMP: the offset MOMENT of the table [us(mT/m)]
                    m_lFirstPartToMeasure,                            // IMP: number of first line/partition (note: Echo=0)
                    m_lLastPartToMeasure,                             // IMP: number of last line/partition (note: Echo=0)
                    m_sTB_3D.getMaxMagnitude(),                       // IMP: maximum allowed gradient amplitude [mT/m]
                    m_sTB_3D.getMinRiseTime(),                        // IMP: minimum allowed rise time [us/(mT/m)]
                   &lRamp,                                            // EXP: minimum required ramp time [us] (up/down)
                   &lDuration                                         // EXP: minimum required duration  [us]
                  );

        if ( setNLSStatus (lStatus, ptModule, "fGSLGetShortestTabTiming for 3D") ) return ( false );


        m_sTB_3D.setRampUpTime   ( lRamp     );          // * Save gradient timing *
        m_sTB_3D.setRampDownTime ( lRamp     );
        m_sTB_3D.setDuration     ( lDuration );

    }





    // * -------------------------------------------------------------------------- *
    // * Calculate gradient timing for phase encode rewinder, kz                    *
    // * -------------------------------------------------------------------------- *

    if ( m_bUseRewinder == true ) {
  
        if ( m_bBalancedGradMomentSliceSelect )  {
            m_dMomentum_TB_3DR = m_dMomentum_TB_3D * m_dFactorForSliceThickness;  // * Balanced gradient moments for TrueFisp *
        } else  {
            m_dMomentum_TB_3DR =  -1.0 * m_dMomentum_TB_3D * m_dRelSSDephasingMoment * m_dFactorForSliceThickness;
        }


        lStatus = fGSLGetShortestTabTiming (
                    dDeltaMomentumKz * m_dFactorForPixelSize3D,      // IMP: Moment between two steps [mT/m * us]
                    SEQ::DIR_DESCENDING,                             // IMP: true <=> (-) -> (+), i.e. table direction
                    m_dMomentum_TB_3DR,                                    // IMP: the offset MOMENT of the table [us(mT/m)]
                    m_lFirstPartToMeasure,                           // IMP: number of first line/partition (note: Echo=0)
                    m_lLastPartToMeasure,                            // IMP: number of last line/partition (note: Echo=0)
                    m_sTB_3DR.getMaxMagnitude(),                     // IMP: maximum allowed gradient amplitude [mT/m]
                    m_sTB_3DR.getMinRiseTime(),                      // IMP: minimum allowed rise time [us/(mT/m)]
                   &lRamp,                                           // EXP: minimum required ramp time [us] (up/down)
                   &lDuration                                        // EXP: minimum required duration  [us]
                );

        if ( setNLSStatus (lStatus, ptModule, "fGSLGetShortestTabTiming for 3D") ) return ( false );

    } else {
        lRamp     = 0L;
        lDuration = 0L;
    }

    m_sTB_3DR.setRampUpTime   ( lRamp     );           // Save gradient timing
    m_sTB_3DR.setRampDownTime ( lRamp     );
    m_sTB_3DR.setDuration     ( lDuration );


    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the PE gradient table and its      *
//	* rewinder for the GRE kernel                       *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBD0354
bool SeqBuildBlockGREKernel::calculatePEGradientTables (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculatePEGradientTables                                        *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the PE gradient table and its rewinder            *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  static const char *ptModule = "SeqBuildBlockGREKernel::calculatePEGradientTables";

  long   lStatus             = SEQU__NORMAL;

  long   lRamp               = 0;
  long   lDuration           = 0;
  
  double dDeltaMomentumKy    = 0.0;


  dDeltaMomentumKy     = fGSLGetPEDeltaMoment( pMrProt );
  
  lStatus = fGSLGetShortestTabTiming (
              dDeltaMomentumKy * m_dFactorForPixelSizePE,      // IMP: Moment between two steps [mT/m * us]
              SEQ::DIR_ASCENDING,                              // IMP: true <=> (-) -> (+), i.e. table direction
              0.0,                                             // IMP: the offset MOMENT of the table [us(mT/m)]
              m_lFirstLineToMeasure,                           // IMP: number of first line/partition (note: Echo=0)
              m_lLastLineToMeasure,                            // IMP: number of last line/partition (note: Echo=0)
              m_sTB_PE.getMaxMagnitude(),                      // IMP: maximum allowed gradient amplitude [mT/m]
              m_sTB_PE.getMinRiseTime(),                       // IMP: minimum allowed rise time [us/(mT/m)]
             &lRamp,                                           // EXP: minimum required ramp time [us] (up/down)
             &lDuration                                        // EXP: minimum required duration  [us]
            );

  if ( setNLSStatus (lStatus, ptModule, "fGSLGetShortestTabTiming for PE") ) return ( false );

  m_sTB_PE.setRampUpTime   ( lRamp     );
  m_sTB_PE.setRampDownTime ( lRamp     );
  m_sTB_PE.setDuration     ( lDuration );


  // * ---------------------------------------------------------------------- *
  // * Set the gradient timing for phase encode rewinding                     *
  // * ---------------------------------------------------------------------- *
  if ( m_bUseRewinder == true) {
    m_sTB_PER.setRampUpTime   ( m_sTB_PE.getRampUpTime   () );
    m_sTB_PER.setRampDownTime ( m_sTB_PE.getRampDownTime () );
    m_sTB_PER.setDuration     ( m_sTB_PE.getDuration ()     );
  } else {
    m_sTB_PER.setRampUpTime   ( 0 );
    m_sTB_PER.setRampDownTime ( 0 );
    m_sTB_PER.setDuration     ( 0 );

  }


  return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation the amplitude of the read out         *
//	* gradients                                         *
//	* The amplitude changes its sign from contrast to   *
//	* contrast.                                         *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBD0352
void SeqBuildBlockGREKernel::calculateROAmplitude ()
{
  // **********************************************************************************
  // * Name        : calculateROAmplitude                                             *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation the amplitudes of the read out gradients             *
  // *               The sign changes from contrast to contrast.                      *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  long    lI;
  double  dSign   = 1.0;
  
  
  for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
    m_sP_RO[lI].setAmplitude ( dSign * m_eROPolarity * m_dMomentRO / m_ADC[lI].getDuration() );
    
    dSign *= -1.0;             // * Toggle the sign of the RO gradient amplitude *
  }
  
}

//	*****************************************************
//	*                                                   *
//	* Calculation of the RO ramp up times for the GRE   *
//	* kernel                                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBD0334
bool SeqBuildBlockGREKernel::calculateRORampUpTime (MrProt* )
{
  // **********************************************************************************
  // * Name        : calculateRORampUpTime                                            *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the RO ramp up times                              *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************

  long lMaxRampUpTime    = 0;     // * Ramp up time for maximal gradient amplitude       *
  long lActualRampUpTime = 0;     // * Ramp up time due to the actual gradient amplitude *
  
  long lI;                        // * Loop counter *
  



  // * ---------------------------------------------------------------------- *
  // * read out gradient ramp up times                                        *
  // * ---------------------------------------------------------------------- *  
  for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
    lMaxRampUpTime    = fSDSRoundUpGRT (m_sP_RO[lI].getMaxMagnitude() * m_sP_RO[lI].getMinRiseTime());    // * Rise time for max. amplitude   *
    lActualRampUpTime = fSDSRoundUpGRT (fabs(m_sP_RO[lI].getAmplitude() * m_sP_RO[lI].getMinRiseTime() * m_dFactorForPixelSizeRO));    // * Rise time for actual amplitude *
        
    m_sP_RO[lI].setRampUpTime ( min( lMaxRampUpTime, lActualRampUpTime ) );
  }

  return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the RO ramp up times for the GRE   *
//	* kernel                                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC03C9
bool SeqBuildBlockGREKernel::calculateRORampDownTime (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculateRORampDownTime                                          *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the read out ramp down times                      *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  static const char *ptModule = "SeqBuildBlockGREKernel::calculateRORampDownTime";

  long    lI;               // * Loop counter *
  double  dSign = 1.0;

  
  for (lI=0; lI < m_lNumberOfEchoes-1; lI++  )  {
    m_sP_RO[lI].setRampDownTime ( m_sP_RO[lI].getRampUpTime() );
  }

  if ( pMrProt->phaseStabilize() )  {

    // * The last RO gradient has to be ramped down with its standard rampt down time *
    // * in case of phase stabilization scans                                         *
    m_sP_RO[m_lNumberOfEchoes-1].setRampDownTime ( m_sP_RO[m_lNumberOfEchoes-1].getRampUpTime() );

  } else {

    // * The last RO gradient overlaps with a RO dephasing gradient.            *
    // * -> Its ramp down time has to be calculated separately.                 *
    lI = m_lNumberOfEchoes - 1;
    switch (m_eReadOutGradType)  {
      case Spoiler:
        // * The read out dephasing gradient was ramped from m_sP_RO[lI].getAmplitude() *
        // * to dSign * m_sP_ROD[lI].getMaxMagnitude()                                  *
        // * m_sP_RO[0].getAmplitude() * dFactorForUpperPixelSizeRO is the minimum      *
        // * allowed RO amplitude in the actual calculation interval                    *

        dSign = ( (lI % 2) ? -1.0 : 1.0 );     // * dSign = +1.0 in case of odd echoes  *
                                               // * dSign = -1.0 in case of even echoes *

        m_sP_RO[lI].setRampDownTime ( fSDSRoundUpGRT( fabs( (dSign * m_sP_ROD.getMaxMagnitude() 
                                      - m_dFactorForUpperPixelSizeRO * m_sP_RO[lI].getAmplitude() ) )
                                      * max (m_sP_RO[lI].getMinRiseTime(), m_sP_ROD.getMinRiseTime() ) ) );
      break;

      case Constant:
        m_sP_RO[lI].setRampDownTime ( 0 );
      break;

      case Rewinder:
      case Symmetric_Rewinder:
        m_sP_RO[lI].setRampDownTime ( fSDSRoundUpGRT( fabs(m_sP_RO[lI].getAmplitude() * m_sP_RO[lI].getMinRiseTime() * m_dFactorForPixelSizeRO) ) );
      break;
    
      default: 
        if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." ) )  return ( false );
    }
  }

  
  return ( true );
  
}

//	*****************************************************
//	*                                                   *
//	* Calculation of the RO durations for the GRE       *
//	* kernel                                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC03AB
bool SeqBuildBlockGREKernel::calculateRODuration (MrProt* )
{
  // **********************************************************************************
  // * Name        : calculateRODuration                                              *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the read out durations                            *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  double   dFirstEchoMoment      = 0.0;
  double   dMomentToSecondEcho   = 0.0;
  long     lDuration1            = 0;


  
  
  
  // * ---------------------------------------------------------------------- *
  // * First echo                                                             *
  // * ---------------------------------------------------------------------- *
  m_sP_RO[0].setDuration( m_sP_RO[0].getRampUpTime() + m_lADCDatTime[0] );


  if ( m_lNumberOfEchoes == 2 )  {
    // * Momentum of the RO prephasing gradient + expected momentum of the RO gradient of the first echo *
    dFirstEchoMoment    =  m_dMomentROP + m_sP_ROP_FC.getMomentumTOT() + m_sP_RO[0].getMomentumTOT();

    // * Minimum momentum of the RO gradient of the second echo *
    dMomentToSecondEcho =  m_sP_RO[1].getAmplitude() * ( m_sP_RO[1].getRampUpTime() * 0.5 + m_lROFlatToEcho[1] );

    if ( fabs(dFirstEchoMoment) < fabs(dMomentToSecondEcho) )  {
      m_sP_RO[0].setDuration( fSDSRoundUpGRT (m_sP_RO[0].getDuration() + ( (fabs(dMomentToSecondEcho) - fabs(dFirstEchoMoment)) / m_sP_RO[0].getAmplitude() )) );
    }

    // * -------------------------------------------------------------------- *
    // * Second echo                                                          *
    // * -------------------------------------------------------------------- *
    m_sP_RO[1].setDuration( m_sP_RO[1].getRampUpTime() + m_lADCDatTime[1] );


    // * Momentum of the RO prephasing and first echo RO gradients *
    dFirstEchoMoment = m_dMomentROP + m_sP_ROP_FC.getMomentumTOT() + m_sP_RO[0].getMomentumTOT();

    if ( fabs(dFirstEchoMoment) > fabs(dMomentToSecondEcho) )  {
      lDuration1 = static_cast<long>((fabs(dFirstEchoMoment) - fabs(dMomentToSecondEcho)) / fabs(m_sP_RO[1].getAmplitude()));

      m_sP_RO[1].setDuration( fSDSRoundUpGRT (m_sP_RO[1].getDuration() + lDuration1) );
      m_lROFlatToEcho[1] += lDuration1;
    }
  }
  

    // * -------------------------------------------------------------------------- *
    // * Calculate duration of the read out interval                                *
    // * -------------------------------------------------------------------------- *
    m_lReadoutTime += m_sP_RO[0].getFlatTopTime();

    if ( m_lNumberOfEchoes == 2 )  {
        m_lReadoutTime += m_sP_RO[0].getRampDownTime() + m_sP_RO[1].getDuration();
    }



  
  return ( true );
  
}

//	*****************************************************
//	*                                                   *
//	* Calculation of the read out prephasing moment and *
//	* sets the ramp up/down times of m_sP_ROP as well
//	*
//	* as its amplitude.                                 *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC035B
bool SeqBuildBlockGREKernel::calculateROP (MrProt* pMrProt, SeqLim* pSeqLim)
{
    // **********************************************************************************
    // * Name        : calculateROP                                                     *
    // *                                                                                *
    // * Class       : SeqBuildBlockGREKernel                                           *
    // *                                                                                *
    // * Description : Calculation of the read out prephasing gradient                  *
    // *                                                                                *
    // * Return      : bool                                                             *
    // *                                                                                *
    // **********************************************************************************
  

    static const char *ptModule = "SeqBuildBlockGREKernel::calculateROP";
  
  
    if ( m_bFlowCompRead )  {

        // * ---------------------------------------------------------------------- *
        // * Calculate flow compensation gradients                                  *
        // * ---------------------------------------------------------------------- *
        MeasNucleus MeasNucleus(pMrProt->txSpec().nucleusInfoArray()[0].nucleus());
        
        m_sP_RO[0].setDuration(m_sP_RO[0].getRampUpTime() + m_lROFlatToEcho[0]);
        if (! fSUCalcShortestFlowComp ( 
                                        MeasNucleus,
                                        0.0,
                                       &m_sP_RO[0],
                                        1.0,
                                        m_dFactorForPixelSizeRO,
                                       &m_sP_ROP_FC,
                                       &m_sP_ROP
                                      )
           )
        {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                if ( setNLSStatus(SEQU_ERROR, ptModule, "Error in calculation of flow compensation (read out). \n") ) return ( false );
            } else  {
                if ( setNLSStatus(SEQU_ERROR, ptModule) ) return ( false );
            }
        }


        m_dMomentROP = m_sP_ROP.getMomentumTOT();

    } else {
        
        // * ---------------------------------------------------------------------- *
        // * Flow compensation gradient sP_ROP_FC is not required.                  *
        // * ---------------------------------------------------------------------- *
        m_sP_ROP_FC.setAmplitude    (0.0);
        m_sP_ROP_FC.setRampUpTime   (0);
        m_sP_ROP_FC.setDuration     (0);
        m_sP_ROP_FC.setRampDownTime (0);



        // * ---------------------------------------------------------------------- *
        // * Moment of the prephasing gradient in read out direction                *
        // * considering asymmetry of k-space sampling                              *
        // * ---------------------------------------------------------------------- *
        if ( m_eReadOutGradType == Symmetric_Rewinder )  {
             m_dMomentROP =  0.5 * (m_sP_RO[0].getRampUpTime() + m_lADCDatTime[0]);
        } else {
            m_dMomentROP =  m_lROFlatToEcho[0] + 0.5 * m_sP_RO[0].getRampUpTime();
        }
        m_dMomentROP *= -1.0 * m_sP_RO[0].getAmplitude();



        // * ---------------------------------------------------------------------- *
        // * Calcution of the shortest possible prephasing gradient in              *
        // * readout direction                                                      *
        // * ---------------------------------------------------------------------- *

        // * Calculates the timing and amplitude for a moment m_dMomentROP * m_dFactorForPixelSizeRO *
        m_sP_ROP.prepSymmetricTOTShortestTime(m_dMomentROP * m_dFactorForPixelSizeRO);
        // * Keeps the timing but sets the amplitude for a moment m_dMomentROP *
        m_sP_ROP.prepMomentumTOT (m_dMomentROP);

    }


    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the read out dephasing gradient    *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC0333
bool SeqBuildBlockGREKernel::calculateROD (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculateROD                                                     *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the read out dephasing gradient                   *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
    static const char *ptModule = "SeqBuildBlockGREKernel::calculateROD";
  
    long   lI;                         // * Loop counter *
    double dSign               = 1.0;
    double dMomentToRewind     = 0.0;

  
  
  
  
    // * -------------------------------------------------------------------------- *
    // * Sign of the read out dephasing gradient                                    *
    // * -------------------------------------------------------------------------- *
    lI = m_lNumberOfEchoes - 1;                       // * Index of the last echo *

    dSign = ( (m_lNumberOfEchoes % 2) ? 1.0 : -1.0 ); // * dSign = +1.0 in case of odd echoes  *
                                                      // * dSign = -1.0 in case of even echoes *

    if ( pMrProt->phaseStabilize() ) {  // * One more ADC in case of phase stabilization *
        dSign *= -1.0;
    } 




    // * -------------------------------------------------------------------------- *
    // * Duration, amplitude and ramp times of the read out dephasing gradient      *
    // * -------------------------------------------------------------------------- *
    switch ( m_eReadOutGradType )  {

        case Spoiler:

            m_sP_ROD.setAmplitude    ( dSign * m_eROPolarity * m_sP_ROD.getMaxMagnitude() );
            if ( pMrProt->phaseStabilize() ) {
                m_sP_ROD.setRampUpTime ( m_PhaseStabScan.getRORampDownTime() );
            } else {
                m_sP_ROD.setRampUpTime ( m_sP_RO[lI].getRampDownTime()         );
            }

            // The momentum of read out dephasing gradient is m_dRelRODephasingMoment of the momentum of the read out gradient *
            m_sP_ROD.setDuration     ( m_sP_ROD.getRampUpTime() + fSDSRoundUpGRT( fabs(m_dMomentRO * m_dFactorForPixelSizeRO * m_dRelRODephasingMoment / m_sP_ROD.getAmplitude() ) ) );

            m_sP_ROD.setRampDownTime ( fSDSRoundUpGRT( fabs(m_sP_ROD.getAmplitude()) * m_sP_ROD.getMinRiseTime() ) );

        break;

        case Constant:
            if ( pMrProt->phaseStabilize() ) {

                m_sP_ROD.setAmplitude    ( m_PhaseStabScan.getROAmplitude()    );
                m_sP_ROD.setRampUpTime   ( m_PhaseStabScan.getRORampDownTime() );
                m_sP_ROD.setDuration     ( fSDSRoundUpGRT( fabs(m_dMomentRO * m_dFactorForPixelSizeRO * m_dRelRODephasingMoment / m_PhaseStabScan.getROAmplitude() ) ) );

            } else {

                m_sP_ROD.setAmplitude    ( m_sP_RO[lI].getAmplitude()    );
                m_sP_ROD.setRampUpTime   ( m_sP_RO[lI].getRampDownTime() );

                // * The gradient was switched on as long as the ADC was switched on.   *
                // * The read gradient is ramped down while the read out dephasing      *
                // * was ramped up at the same time, resulting in a constant gradient   *
                // * amplitude.                                                         *

                m_sP_ROD.setDuration     ( fSDSRoundUpGRT( (double)(m_sP_RO[lI].getRampDownTime() + m_lADCZeroAfter[lI])) );
            }

            m_sP_ROD.setRampDownTime ( fSDSRoundUpGRT(fabs(m_sP_ROD.getAmplitude() * m_dFactorForPixelSizeRO * m_sP_ROD.getMinRiseTime()) ) );
        break;

        case Symmetric_Rewinder:
            m_sP_ROD.setRampUpTime   ( m_sP_ROP.getRampUpTime()   );
            m_sP_ROD.setDuration     ( m_sP_ROP.getDuration()     );
            m_sP_ROD.setRampDownTime ( m_sP_ROP.getRampDownTime() );
            m_sP_ROD.setAmplitude    ( m_sP_ROP.getAmplitude()    );
        break;

        case Rewinder:
            // * Total moment in RO direction (RO dephasing gradient + all RO gradients)  *
            dMomentToRewind = m_dMomentROP;

            for (lI=0; lI < m_lNumberOfEchoes; lI++  )  {
                dMomentToRewind += m_sP_RO[lI].getMomentumTOT();
            }

            // * Sets the amplitude, duration and ramp times to realize the gradient moment -1.0 * dMomentumToRewind *
            m_sP_ROD.prepSymmetricTOTShortestTime( -1.0 * dMomentToRewind * m_dFactorForPixelSizeRO);
            m_sP_ROD.prepMomentumTOT             ( -1.0 * dMomentToRewind                          );
        break;
    
        default: 
            setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." );
            return ( false );

    }
  
    return ( true );
  
}

//	*****************************************************
//	*                                                   *
//	* Calculation of the prephase duration              *
//	* All prephasing gradients are stretched to fit     *
//	* into the prephasing interval                      *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC02F6
bool SeqBuildBlockGREKernel::calculatePrephaseTime (MrProt* , SeqLim* )
{
    // **********************************************************************************
    // * Name        : calculatePrephaseTime                                            *
    // *                                                                                *
    // * Class       : SeqBuildBlockGREKernel                                           *
    // *                                                                                *
    // * Description : Calculation of the prephase duration                             *
    // *                                                                                *
    // * Return      : BOOL                                                             *
    // *                                                                                *
    // **********************************************************************************
  
    static const char *ptModule = "SeqBuildBlockGREKernel::calculatePrephaseTime";

    long lTimeSlice   = 0;
    long lTimePhase   = 0;
    long lTimeRead    = 0;

  
  
    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }


    // * ---------------------------------------------------------------------- *
    // * Calculate the maximum prephasing time on all gradient axes             *
    // * ---------------------------------------------------------------------- *
    lTimeSlice    = m_pRF_WE->getRampTimeOutsideSBB() + m_sP_3D_FC.getTotalTime() + m_sTB_3D.getTotalTime();
    lTimePhase    = m_sTB_PE.getTotalTime();
    lTimeRead     = m_sP_ROP_FC.getTotalTime() + m_sP_ROP.getTotalTime() + m_sP_RO[0].getRampUpTime();




    m_lPrephaseTime = max (lTimeSlice     , lTimePhase                                    );
    m_lPrephaseTime = max (m_lPrephaseTime, lTimeRead                                     );
    m_lPrephaseTime = max (m_lPrephaseTime, m_lADCZeroBefore[0] + m_lMinDelayBetweenRFAndADC);


    // * ---------------------------------------------------------------------- *
    // * Which gradient axis limits the prephasing time ?                       *
    // * ---------------------------------------------------------------------- *
    if ( lTimeSlice == m_lPrephaseTime ) {
        m_eLimitingPrephaseTime = Slice;
    } else if ( lTimePhase == m_lPrephaseTime ) {
        m_eLimitingPrephaseTime = Phase;
    } else if ( lTimeRead == m_lPrephaseTime )  {
        m_eLimitingPrephaseTime = Read;
    } else {
        m_eLimitingPrephaseTime = None;
    }


    // * ---------------------------------------------------------------------- *
    // * m_lADCZeroBefore[0] might not be a multiple of the gradient raster     *
    // * time -> m_lPrephaseTime[0] has be rounded to the gradient raster       *
    // * time.                                                                  *
    // * ---------------------------------------------------------------------- *
    m_lPrephaseTime = fSDSRoundUpGRT ( m_lPrephaseTime );
    m_lTimeToADC[0] = m_lPrephaseTime - m_lADCZeroBefore[0];



    // * ---------------------------------------------------------------------- *
    // * All prephasing gradients are stretched to fill the maximum             *
    // * prephasing time.                                                       *
    // * ---------------------------------------------------------------------- *
    m_lPrephaseFillTimeSlice = m_lPrephaseFillTimePhase = m_lPrephaseFillTimeRead = 0;

    if ( lTimeSlice < m_lPrephaseTime )  {
        if ( m_bFlowCompRead || m_bFlowCompSlice ) {
            m_lPrephaseFillTimeSlice = m_lPrephaseTime - lTimeSlice;
        } else {
            m_sTB_3D.setDuration ( m_sTB_3D.getDuration() + m_lPrephaseTime - lTimeSlice );
        }

    }

    if ( lTimePhase < m_lPrephaseTime )  {
        if ( m_bFlowCompRead || m_bFlowCompSlice )  {
            m_lPrephaseFillTimePhase = m_lPrephaseTime - lTimePhase;
        } else {
            m_sTB_PE.setDuration ( m_sTB_PE.getDuration() + m_lPrephaseTime - lTimePhase );
        }
        
    }

    if ( lTimeRead < m_lPrephaseTime )  {
        if ( m_bFlowCompRead || m_bFlowCompSlice ) {

            m_lPrephaseFillTimeRead = m_lPrephaseTime - lTimeRead;

        } else {
            // * ------------------------------------------------------------------ *
            // * No flow compensation in read out direction                         *
            // * Use all the prephase time for m_sP_ROP. m_sP_ROP_FC is zero.       *
            // * ------------------------------------------------------------------ *
            m_sP_ROP.setDuration ( m_sP_ROP.getDuration() + m_lPrephaseTime - lTimeRead );
        }
        
    }



    // * ---------------------------------------------------------------------- *
    // * Final calculatation of the relative prephasing gradient amplitude      *
    // * in read out direction                                                  *
    // * No recalculation of m_sP_ROP in case of flow compensation:             *
    // *    - flow compensation in read out direction: the timing must not be   *
    // *      changed to preserve flow                                          *
    // *    - slice compensation in slice direction: the read out prephasing    *
    // *      gradient should be stretched to reduce flow effects               *
    // * compensation.                                                          *
    // * ---------------------------------------------------------------------- *
    if ( ! (m_bFlowCompRead || m_bFlowCompSlice) ) {
        m_sP_ROP.setAmplitude ( m_dMomentROP / static_cast<double>(m_sP_ROP.getDuration()) );
    }



    return ( true );
  
}

//	*****************************************************
//	*                                                   *
//	* Calculation of the rephase duration               *
//	* All rephasing gradients are stretched to fit      *
//	* into the rephasing interval                       *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC02CE
bool SeqBuildBlockGREKernel::calculateRephaseTime (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculateRephaseTime                                             *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the rephase duration                              *
  // *               All rephasing gradients are stretched to fit into the rephasing  *
  // *               interval                                                         *
  // *                                                                                *
  // * Return      : void                                                             *
  // *                                                                                *
  // **********************************************************************************

    static const char *ptModule = "SeqBuildBlockGREKernel::calculateRephaseTime";

    long  lTimeSlice    = 0;
    long  lTimePhase    = 0;
    long  lTimeRead     = 0;

    long  lOldGradTime  = 0;     // * Temp. variable     *
    long  lNewGradTime  = 0;



    if ( m_pRF_WE == NULL )  {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE = NULL.\n");
        return ( false );
    }



    // * -------------------------------------------------------------------------- *
    // * Calculate the maximum rephasing time on the slice select, phase encode     *
    // * or read out gradient axis                                                  *
    // * In case of phase stabilization ROD is carried out AFTER the phase          *
    // * stabilization scan and DOES NOT contribute to rephasing time.              *
    // * -------------------------------------------------------------------------- *

    lTimeSlice = m_sTB_3DR.getTotalTime();
    lTimePhase = m_sTB_PER.getTotalTime();
    lTimeRead  = m_sP_ROD.getTotalTime();


    if ( m_eReadOutGradType == Symmetric_Rewinder )  {
         lTimeSlice += m_pRF_WE->getRampTimeOutsideSBB();
    }

    if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) )  {   // * The read out and read out dephasing gradients do         *
         lTimeRead += m_sP_RO[m_lNumberOfEchoes-1].getRampDownTime();                          // * NOT overlap in case of m_eReadOutGradType equal Rewinder *
    }

    m_lRephaseTime = max(lTimeSlice ,    lTimePhase);
    if ( ! pMrProt->phaseStabilize() )  {
        m_lRephaseTime = max(m_lRephaseTime, lTimeRead );
    }
    m_lRephaseTime = max(m_lRephaseTime, m_lADCZeroAfter[m_lNumberOfEchoes-1]);

    m_lRephaseTime = fSDSRoundUpGRT ( m_lRephaseTime );


    // * -------------------------------------------------------------------------- *
    // * The rephasing gradients in slice select and phase encode                   *
    // * direction and the dephasing gradient in read out direction are             *
    // * stretched to fill the maximum rephasing time.                              *
    // *                                                                            *
    // * If this operation would be carried out in case of an gradient duration     *
    // * that is equal to 0, dangling gradient ramps would be impossible            *
    // * -------------------------------------------------------------------------- *

    if (lTimeSlice < m_lRephaseTime) {
        if ( m_sTB_3DR.getRampDownTime() != 0 ) {
            m_sTB_3DR.setDuration( m_sTB_3DR.getDuration() + m_lRephaseTime - lTimeSlice );
        }
    }

    if (lTimePhase < m_lRephaseTime) {
        if ( m_sTB_PER.getRampDownTime() != 0 ) {
            m_sTB_PER.setDuration( m_sTB_PER.getDuration() + m_lRephaseTime - lTimePhase );
        }
    }

    if ( (lTimeRead < m_lRephaseTime) && ( ! pMrProt->phaseStabilize() ) ) {
        if ( (m_eReadOutGradType == Spoiler) || (m_eReadOutGradType == Rewinder) ) {
            // * The amplitude of the read out dephasing gradient can be reduced    *
            // * since its duration is increased:                                   *
            // * NewAmplitude = OldAmplitude * OldGradientTime / NewGradientTime    *


            lOldGradTime = static_cast<long>(m_sP_ROD.getFlatTopTime() + 0.5 * ( m_sP_ROD.getRampUpTime() + m_sP_ROD.getRampDownTime() ));
            lNewGradTime = static_cast<long>(lOldGradTime + m_lRephaseTime - lTimeRead );


            m_sP_ROD.setAmplitude( m_sP_ROD.getAmplitude() * lOldGradTime / lNewGradTime );
        }
        m_sP_ROD.setDuration( m_sP_ROD.getDuration() + m_lRephaseTime - lTimeRead );
    }
  
    return ( true );  

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the duration that gradient ramps   *
//	* are allowed to stick out of the event block       *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC027F
bool SeqBuildBlockGREKernel::calculateOverTime (MrProt* pMrProt, SeqLim*)
{
  // **********************************************************************************
  // * Name        : calculateOverTime                                                *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the duration that gradient ramps are allowed to   *
  // *               stick out of the event block                                     *
  // *                                                                                *
  // * Return      : void                                                             *
  // *                                                                                *
  // **********************************************************************************
  

    static const char *ptModule = "SeqBuildBlockGREKernel::calculateOverTime";

    long              lMinDurationOfRephase = 0;

//    static RX4Proxy   theRX4Proxy;
    long              lMinDurationBetweenADCAndRF = 0;



    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
    }




    // * -------------------------------------------------------------------------- *
    // *                                                                            *
    // * lMinDurationOfRephase is the minimum duration of the rephase time          *
    // * that ensures that the time between the end of the ADC and the start        *
    // * of the next rf-pulse is greater than lMinDurationBetweenADCAndRF           *
    // *                                                                            *
    // * -------------------------------------------------------------------------- *

    if ( m_bRampsOutsideOfEventBlock == true )  {

        if ( pMrProt->phaseStabilize() )  {

        // * The last ADC is the phase stabilization ADC *
        lMinDurationBetweenADCAndRF = SysProperties::getMinDurationBetweenReadoutAndRFPulse() /*static_cast<long>(theRX4Proxy.getMinDurationBetweenReadoutAndRFPulse ( (double)m_PhaseStabScan.getlEffectiveDwellTime_ns() / pSeqLim->getReadoutOSFactor() ) )*/;
        lMinDurationOfRephase = max( lMinDurationBetweenADCAndRF - m_lPhaseEventDelay - m_pRF_WE->getRequiredLeadTime() - m_pRF_WE->getGSData().getlTimeToFirstFlatTop(), 0 );
        lMinDurationOfRephase = fSDSRoundUpGRT ( lMinDurationOfRephase );

        m_lOverTime = min ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );

        if ( (m_sP_ROD.getTotalTime() - m_lOverTime) < lMinDurationOfRephase ) {
            m_lOverTime -= lMinDurationOfRephase - (m_sP_ROD.getTotalTime() - m_lOverTime);
        }

        } else {

        // * Last kernel ADC *
        lMinDurationBetweenADCAndRF = SysProperties::getMinDurationBetweenReadoutAndRFPulse() /*static_cast<long>(theRX4Proxy.getMinDurationBetweenReadoutAndRFPulse ( pMrProt->rxSpec().realDwellTime()[pMrProt->contrasts()-1] / 1000.0))*/;
        lMinDurationOfRephase = m_lADCZeroAfter[m_lNumberOfEchoes-1] + max( lMinDurationBetweenADCAndRF - m_lPhaseEventDelay - m_pRF_WE->getRequiredLeadTime() - m_pRF_WE->getGSData().getlTimeToFirstFlatTop(), 0 );
        lMinDurationOfRephase = fSDSRoundUpGRT ( lMinDurationOfRephase );

        // * ---------------------------------------------------------------------- *
        // *                                                                        *
        // * Only the gradient with the longest FlatTopTime (end of the event block)*
        // * can have a dangling gradient ramp                                      *
        // *                                                                        *
        // * ---------------------------------------------------------------------- *
    
        if ( m_sTB_3DR.getDuration() > m_sTB_PER.getDuration() )  {

            if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) )  {
      
                if ( m_sTB_3DR.getDuration() > m_sP_RO[m_lNumberOfEchoes-1].getRampDownTime() + m_sP_ROD.getDuration() ) {
                    m_lOverTime = min ( m_sTB_3DR.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                } else {
                    m_lOverTime = min ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                }
        
            } else {
      
                if ( m_sTB_3DR.getDuration() > m_sP_ROD.getDuration() ) {
                    m_lOverTime = min ( m_sTB_3DR.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                } else {
                    m_lOverTime = min ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
                }     
        
            }

        } else {
            if ( m_sTB_PER.getDuration() > m_sP_ROD.getDuration() )  {
                m_lOverTime = min ( m_sTB_PER.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
            } else {
                m_lOverTime = min ( m_sP_ROD.getRampDownTime(), m_pRF_WE->getRequiredLeadTime() + m_pRF_WE->getGSData().getlTimeToFirstFlatTop() );
            }
        }
    
  
    
        if ( (m_lRephaseTime - m_lOverTime) < lMinDurationOfRephase ) {
            m_lOverTime -= lMinDurationOfRephase - (m_lRephaseTime - m_lOverTime);
        }

        }


    } else {
  
        if ( pMrProt->phaseStabilize() )  {

            // * No gradient ramps outside of the event block *
            if ( m_sP_ROD.getTotalTime() < lMinDurationOfRephase ) {
                m_lOverTime -= lMinDurationOfRephase - m_sP_ROD.getTotalTime();
            } else {
                m_lOverTime =  0;
            }

        } else {

            // * No gradient ramps outside of the event block *
            if ( m_lRephaseTime < lMinDurationOfRephase ) {
                m_lOverTime -= lMinDurationOfRephase - m_lRephaseTime;
            } else {
                m_lOverTime =  0;
            }

        }

    }

    
    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the minimum TE times               *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC0274
bool SeqBuildBlockGREKernel::calculateTEMin ()
{
    // **********************************************************************************
    // * Name        : calculateTEMin                                                   *
    // *                                                                                *
    // * Class       : SeqBuildBlockGREKernel                                           *
    // *                                                                                *
    // * Description : Calculation of the minimum TE times                              *
    // *                                                                                *
    // * Return      : BOOL                                                             *
    // *                                                                                *
    // **********************************************************************************
  
    long lDuration1 = 0;   // * Temp. variables *
    long lDuration2 = 0;



    // * -------------------------------------------------------------------------- *
    // * lTEMin is rounded to the gradient raster time in order to facilitate the   *
    // * evaluation of the TE fill times                                            *
    // * -------------------------------------------------------------------------- *

    m_alTEMin[0] = m_lTEFromExcitation + m_lDelayAdjustSSPosition + m_lPrephaseTime + m_lROFlatToEcho[0];
    m_alTEMin[0] = fSDSRoundUpGRT ((double)m_alTEMin[0]);
  


    // * -------------------------------------------------------------------------- *
    // * Reset m_lDiffBetweenROAndADCTiming if there is only one contrast           *
    // * -------------------------------------------------------------------------- *
    if ( m_lNumberOfEchoes == 1 )  {
        m_lDiffBetweenROAndADCTiming = 0L;
    }



    if ( m_lNumberOfEchoes == 2 )  {
        // * m_alTEMin[1] may be limited by overlapping RO gradients or by          *
        // * overlapping ADC events                                                 *

        // * Minimum TE due to RO gradients                                         *
        lDuration1 = fSDSRoundUpGRT ( m_lTEFromExcitation + m_lPrephaseTime + m_sP_RO[0].getFlatTopTime() + m_sP_RO[0].getRampDownTime() +
                                      m_sP_RO[1].getRampUpTime() + m_lROFlatToEcho[1] );

        // * Minimum TE due to ADC events                                           *
        lDuration2   = fSDSRoundUpGRT ( m_lTEFromExcitation + m_lTimeToADC[0] + m_ADC[0].getDuration() +
                                        m_lMinDelayBetweenADCs +      // **** neu ****
                                        m_lADCZeroBefore[1] + static_cast<long>( m_lADCDatTime[1] * m_adRelWinAsymCol[1] + 0.5) );


        m_alTEMin[1] = max (lDuration1, lDuration2);


        // * ---------------------------------------------------------------------- *
        // * The calculation of TRmin as well as the timinig in runTiming are based *
        // * on the assumption that the RO gradients are limiting the minimum TE.   *
        // * If this assumption does not hold an additional delay has to be         *
        // * considered.                                                            *
        // * ---------------------------------------------------------------------- *
        if ( lDuration1 > lDuration2 )  {
            m_lDiffBetweenROAndADCTiming = 0L;
        } else {
            m_lDiffBetweenROAndADCTiming = lDuration2 - lDuration1;
        }

    }


    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the TE fill times                  *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC024C
bool SeqBuildBlockGREKernel::calculateTEFill (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculateTEFill                                                  *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the TE fill times                                 *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************


  m_alTEFill[0]   = static_cast<long>( ( pMrProt->te()[0] - m_alTEMin[0] ) / 2.0 );
  
  if ( m_lNumberOfEchoes == 2 )  {
      m_alTEFill[1] = static_cast<long>( ( pMrProt->te()[1] - m_alTEMin[1] - 2 * m_alTEFill[0]) / 2.0 );
  } else {
      m_alTEFill[1] = 0;
  }

  return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Calculation of the minimum TR time                *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC0224
bool SeqBuildBlockGREKernel::calculateTRMin (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculateTRMin                                                   *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the minimum TR time                               *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************


    // * -------------------------------------------------------------------------- *
    // * Gradient ramps that are dangling out of the event block (m_lOverTime) do   *
    // * not contribute to the minimum TR time                                      *
    // * -------------------------------------------------------------------------- *
    m_lTRMin = m_lExcitationAllTime  + m_lDelayAdjustSSPosition + m_lPrephaseTime + (2 * m_alTEFill[0]) +
               m_lReadoutTime        + m_lRephaseTime  + m_lPhaseEventDelay - m_lOverTime ;


    if ( pMrProt->phaseStabilize() )  {
        m_lTRMin += m_PhaseStabScan.getDurationPerRequest() + m_sP_ROD.getTotalTime();
    }


    if ( m_lNumberOfEchoes == 2 )  {
        m_lTRMin += (2 * m_alTEFill[1]) + m_lDiffBetweenROAndADCTiming;
    }



    return ( true );


}

//	*****************************************************
//	*                                                   *
//	* Calculation of the time between the both ADC      *
//	* intervals                                         *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC01F2
bool SeqBuildBlockGREKernel::calculateTimeToADC1 (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : calculateTimeToADC1                                              *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Calculation of the time between both ADC intervals               *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  if ( m_lNumberOfEchoes == 2 )  {
    m_lTimeToADC[1] = pMrProt->te()[1] - m_lTEFromExcitation - m_lTimeToADC[0] - 2 * m_alTEFill[0] - 
                      static_cast<long>(m_ADC[0].getDuration()) -
                      static_cast<long>(m_ADC[1].getDuration() * m_adRelCenterCol[1]);

  }


  
  return ( true );
  
}

//	Can be used to perform additional calculations
//##ModelId=3D352EBC01A2
bool SeqBuildBlockGREKernel::calculateAddOn (MrProt* , SeqLim* , SeqExpo* )
{
  // **********************************************************************************
  // * Name        : prep                                                             *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Additional calculations. Not used in GREKernel.                  *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

    return ( true );
    
}

//	*****************************************************
//	*                                                   *
//	* Preparation                                       *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC00284
bool SeqBuildBlockGREKernel::prep (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
  // **********************************************************************************
  // * Name        : prep                                                             *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Prepare and check of the gardients, tables and ADC events        *
  // *               Calculation of the energy per request and the duration per       *
  // *               request.                                                         *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

  static const char *ptModule = "SeqBuildBlockGREKernel::prep";

  long              lI;                               // * Loop counter *
  double            dAmplitudeOfLeftNeighbour = 0.0;  // * Used for gradient slew rate checks *




  // *------------------------------------------------------------------------*
  // * Debug: running message and imported parameters                         *
  // *------------------------------------------------------------------------*
  if ( mIsTRrun )  {
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s running\n"                                                                ,ptModule);
  }

  if ( mIsTRinp )  {
    TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Imported Parameter                                *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);

    if ( m_bCalculated )  {
      TRACE_PUT1(TC_INFO, TF_SEQ, "%s Kernel was calculated.\n" ,ptModule);
    } else  {
      TRACE_PUT1(TC_INFO, TF_SEQ, "%s Kernel has not been calculated, yet.\n" ,ptModule);
    }
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
  }


  if (m_bCalculated == false)  {

    if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
      setNLSStatus(SEQU_ERROR, ptModule, "Kernel timing has not been calculated, yet.\n");     return ( false );
    } else {
      setNLSStatus(SEQU_ERROR);     return ( false );
    }

  }



  // *------------------------------------------------------------------------*
  // * Prepare the osc bit                                                    *
  // *------------------------------------------------------------------------*
  if ( m_bSendOscBit == true )  {
    m_sOscBit.setIdent    ( RTEIDENT_Osc0 );
    m_sOscBit.prep        ( 0, 10 );          // * Channel: 0, duration: 10us *
  }




  // *------------------------------------------------------------------------*
  // * Prepare the wake up bit                                                *
  // * WakeUp event must have length zero. This is set already by the         *
  // * constructor!                                                           *
  // *------------------------------------------------------------------------*
  
  // WakeUp event must have length zero. This is set already by the constructor!
  
  m_sWakeUp.lCode       = SYNCCODE_WAKEUP;
  m_sWakeUp.setIdent    ( "WakeUp" );



  // * ---------------------------------------------------------------------- *
  // * Check the excitation pulse                                             *
  // * ---------------------------------------------------------------------- *
  if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
  }
  
  if (! m_pRF_WE->checkGradients(pMrProt, pSeqLim) )  {
      if ( ! pSeqLim->isContextPrepForBinarySearch() ) {
          setNLSStatus ( m_pRF_WE->getNLSStatus(), "m_pRF_WE->checkGradients" );    return ( false );
      } else {
          setNLSStatus ( m_pRF_WE->getNLSStatus() );    return ( false );
      }
  }




  // * ---------------------------------------------------------------------- *
  // * Prepare and check gradient table rewinder in slice select direction    *
  // * ---------------------------------------------------------------------- *
  if (! prepAndCheck3D( pMrProt ) ) return ( false );
  
  


  // * ---------------------------------------------------------------------- *
  // * Prepare and check phase encoding table                                 *
  // * ---------------------------------------------------------------------- *
  if (! m_sTB_PE.prepPE(m_sTB_PE.getRampUpTime(), m_sTB_PE.getDuration(), m_sTB_PE.getRampDownTime(), pMrProt, SEQ::DIR_ASCENDING, 0.0, m_lFirstLineToMeasure))  { /*! EGA-03; EGA-01 !*/
    setNLSStatus (m_sTB_PE.getNLSStatus(), "m_sTB_PE.prepPE");
    return ( false );
  }

  if (! m_sTB_PE.check() )  {
    setNLSStatus (m_sTB_PE.getNLSStatus(), "m_sTB_PE.check");
    return ( false );
  }

  if (! m_sTB_PE.prepPE(m_sTB_PE.getRampUpTime(), m_sTB_PE.getDuration(), m_sTB_PE.getRampDownTime(), pMrProt, SEQ::DIR_ASCENDING, 0.0, m_lLastLineToMeasure))  { /*! EGA-03; EGA-01 !*/
    setNLSStatus (m_sTB_PE.getNLSStatus(), "m_sTB_PE.prepPE");
    return ( false );
  }

  if (! m_sTB_PE.check() )  {
    setNLSStatus (m_sTB_PE.getNLSStatus(), "m_sTB_PE.check");
    return ( false );
  }


  // * ---------------------------------------------------------------------- *
  // * Prepare and check phase encoding rewinder table                        *
  // * ---------------------------------------------------------------------- *
  if ( m_bUseRewinder == true) {
    if (! m_sTB_PER.prepPE(m_sTB_PER.getRampUpTime(), m_sTB_PER.getDuration(), m_sTB_PER.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, m_lFirstLineToMeasure))  {
      setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE");
      return ( false );
    }

    if (! m_sTB_PER.check() )  {
      setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.check");
      return ( false );
    }

    if (! m_sTB_PER.prepPE(m_sTB_PER.getRampUpTime(), m_sTB_PER.getDuration(), m_sTB_PER.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, m_lLastLineToMeasure))  {
      setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE");
      return ( false );
    }

    if (! m_sTB_PER.check() )  {
      setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.check");
      return ( false );
    }

  } else {
    if (! m_sTB_PER.prepPE(m_sTB_PER.getRampUpTime(), m_sTB_PER.getDuration(), m_sTB_PER.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, 0) )  {
      setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE");
      return ( false );
    }
  }


  // * ---------------------------------------------------------------------- *
  // * Prepare and check read out prephasing gradient                         *
  // * ---------------------------------------------------------------------- *
  if (! prepAndCheckROP( pMrProt, pSeqLim ) ) return ( false );



  // * ---------------------------------------------------------------------- *
  // * Prepare and check read out gradients                                   *
  // * ---------------------------------------------------------------------- *
  if (! prepAndCheckRO( pMrProt, pSeqLim ) ) return ( false );



  // * ---------------------------------------------------------------------- *
  // * Prepare and check read out dephasing gradient                          *
  // * ---------------------------------------------------------------------- *
  if (! m_sP_ROD.prepAmplitude(m_sP_ROD.getRampUpTime(), m_sP_ROD.getDuration(), m_sP_ROD.getRampDownTime(), m_sP_ROD.getAmplitude() ))  {
    setNLSStatus (m_sP_ROD.getNLSStatus(), "SeqBuildBlockGREKernel::prep m_sP_ROD.prepAmplitude");
    return ( false );
  }




  // * ---------------------------------------------------------------------- *
  // * Prepare and check additional gradients                                 *
  // * ---------------------------------------------------------------------- *
  if (! prepAndCheckAddOn ( pMrProt, pSeqLim, pSeqExpo ) ) return ( false );




  if ( pMrProt->phaseStabilize() )  {
    dAmplitudeOfLeftNeighbour = m_PhaseStabScan.getROAmplitude();            // * Amplitude of phase stabilize RO *
  } else {
    dAmplitudeOfLeftNeighbour = m_sP_RO[m_lNumberOfEchoes-1].getAmplitude(); // * Amplitude of last kernel RO     *
  }

  switch ( m_eReadOutGradType )  {
    case Spoiler:
      if (! m_sP_ROD.check( dAmplitudeOfLeftNeighbour, 0.0 ) )  {
        setNLSStatus (m_sP_ROD.getNLSStatus(), "SeqBuildBlockGREKernel::prep m_sP_ROD.check");
        return ( false );
      }
    break;

    case Constant:
      if (! m_sP_ROD.check( dAmplitudeOfLeftNeighbour, 0.0 ) )  {
        setNLSStatus (m_sP_ROD.getNLSStatus(), "SeqBuildBlockGREKernel::prep m_sP_ROD.check");
        return ( false );
      }
     break;

    case Rewinder:
    case Symmetric_Rewinder:
      if (! m_sP_ROD.check() )  {
        setNLSStatus (m_sP_ROD.getNLSStatus(), "SeqBuildBlockGREKernel::prep m_sP_ROD.check");
        return ( false );
      }
    break;
    
    default: 
      if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." ) )  return ( false );
  }





  // * ---------------------------------------------------------------------- *
  // * Calculate temp. variables for the kernel timing                        *
  // * ---------------------------------------------------------------------- *

  lI = m_lNumberOfEchoes - 1;
//  m_lRODuration = m_sP_RO[lI].getDuration();                                               // * Gradients m_sP_RO[lI] and m_sP_ROD overlap in case of *
                                                                                           // * eReadOutType equal Constant or Spoiler                *
//  if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) ) {  // * m_sP_RO[lI] was ramped down before m_sP_ROD was       *
//    m_lRODuration += m_sP_RO[lI].getRampDownTime();                                        // * ramped up (no overlap) in case of eReadOutType        *
//  }                                                                                        // * equal Rewinder.or Symmetric_Rewinder                  *



  // * ---------------------------------------------------------------------- *
  // * Calculate energy of one request                                        *
  // * ---------------------------------------------------------------------- *
  m_dEnergyPerRequest_Ws = m_pRF_WE->getEnergyPerRequest();


  // * ---------------------------------------------------------------------- *
  // * Calculate SBB duration of one request                                  *
  // * ---------------------------------------------------------------------- *
  m_lSBBDurationPerRequest_us = m_lTRMin + m_lTRFill + m_lTRFillEnd;




  if ( mIsTRend )  {
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s finished\n", ptModule);
  }



  // * ---------------------------------------------------------------------- *
  // * Successful preparation of the GRE kernel                               *
  // * ---------------------------------------------------------------------- *
  setPrepared();
    
  return ( true );


}

//	*****************************************************
//	*                                                   *
//	* Prepares and checks the 3D gradient tables        *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC017A
bool SeqBuildBlockGREKernel::prepAndCheck3D (MrProt* pMrProt)
{
  // **********************************************************************************
  // * Name        : prepAndCheck3D                                                   *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Prepares and checks the 3D gradient tables.                      *
  // *               Returns true if prep and check has been successfully been        *
  // *               finished otherwise false.                                        *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  
  
  
  // * ---------------------------------------------------------------------- *
  // * Prepare and check gradient table in slice select direction             *
  // * ---------------------------------------------------------------------- *
  if (! m_sTB_3D.prep3D(m_sTB_3D.getRampUpTime(), m_sTB_3D.getDuration(), m_sTB_3D.getRampDownTime(), pMrProt, SEQ::DIR_ASCENDING, m_dMomentum_TB_3D, m_lFirstPartToMeasure))  {
    setNLSStatus (m_sTB_3D.getNLSStatus(), "m_sTB_3D.prep3D");
    return ( false );
  }

  if (! m_sTB_3D.check() )  {
    setNLSStatus (m_sTB_3D.getNLSStatus(), "m_sTB_3D.check");
    return ( false );
  }

  if (! m_sTB_3D.prep3D(m_sTB_3D.getRampUpTime(), m_sTB_3D.getDuration(), m_sTB_3D.getRampDownTime(), pMrProt, SEQ::DIR_ASCENDING, m_dMomentum_TB_3D, m_lLastPartToMeasure))  {
    setNLSStatus (m_sTB_3D.getNLSStatus(), "m_sTB_3D.prep3D");
    return ( false );
  }

  if (! m_sTB_3D.check() )  {
    setNLSStatus (m_sTB_3D.getNLSStatus(), "m_sTB_3D.check");
    return ( false );
  }


  if ( m_bFlowCompSlice ) {

      if (! m_sP_3D_FC.prepAmplitude(m_sP_3D_FC.getRampUpTime(), m_sP_3D_FC.getDuration(), m_sP_3D_FC.getRampDownTime(), m_sP_3D_FC.getAmplitude() ))  {
          setNLSStatus (m_sP_3D_FC.getNLSStatus(), "m_sP_3D_FC.prepAmplitude");
          return ( false );
      }

      if (! m_sP_3D_FC.check())  {
          setNLSStatus (m_sP_3D_FC.getNLSStatus(), "m_sP_3D_FC.check");
          return ( false );
      }

  }
  
  
  
  
  // * ---------------------------------------------------------------------- *
  // * Prepare and check gradient table rewinder in slice select direction    *
  // * ---------------------------------------------------------------------- *
  if ( m_bUseRewinder == true) {
    if ( ! m_bBalancedGradMomentSliceSelect )  { 
      m_dMomentum_TB_3DR = -1.0 * m_dMomentum_TB_3D * m_dRelSSDephasingMoment;
    } else { 
      m_dMomentum_TB_3DR = m_dMomentum_TB_3D;    // * Balanced slice selecet gradients, used for True Fisp *
    }
    
    if (! m_sTB_3DR.prep3D(m_sTB_3DR.getRampUpTime(), m_sTB_3DR.getDuration(), m_sTB_3DR.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, m_dMomentum_TB_3DR, m_lFirstPartToMeasure))  {
      setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D");
      return ( false );
    }

    if (! m_sTB_3DR.check() )  {
      setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.check");
      return ( false );
    }

    if (! m_sTB_3DR.prep3D(m_sTB_3DR.getRampUpTime(), m_sTB_3DR.getDuration(), m_sTB_3DR.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, m_dMomentum_TB_3DR, m_lLastPartToMeasure))  {
      setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D");
      return ( false );
    }

    if (! m_sTB_3DR.check() )  {
      setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.check");
      return ( false );
    }

  } else {
    if (! m_sTB_3DR.prep3D(m_sTB_3DR.getRampUpTime(), m_sTB_3DR.getDuration(), m_sTB_3DR.getRampDownTime(), pMrProt, SEQ::DIR_DESCENDING, 0.0, 0) )  {
      setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D");
      return ( false );
    }
  }

  return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Prepares and checks the read out dephasing        *
//	* gradient                                          *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC013E
bool SeqBuildBlockGREKernel::prepAndCheckROP (MrProt* , SeqLim* )
{
    // **********************************************************************************
    // * Name        : prepAndCheckROP                                                  *
    // *                                                                                *
    // * Class       : SeqBuildBlockGREKernel                                           *
    // *                                                                                *
    // * Description : Prepares and checks the read out prephasing gradients.           *
    // *               Returns true if prep and check has been successfully been        *
    // *               finished otherwise false.                                        *
    // *                                                                                *
    // * Return      : bool                                                             *
    // *                                                                                *
    // **********************************************************************************




    // * -------------------------------------------------------------------------- *
    // * Read out flow compensation gradient                                        *
    // * -------------------------------------------------------------------------- *
    if ( m_bFlowCompRead ) {

        if (! m_sP_ROP_FC.prepAmplitude(m_sP_ROP_FC.getRampUpTime(), m_sP_ROP_FC.getDuration(), m_sP_ROP_FC.getRampDownTime(), m_sP_ROP_FC.getAmplitude() ))  {
            setNLSStatus (m_sP_ROP_FC.getNLSStatus(), "m_sP_ROP_FC.prepAmplitude");
            return ( false );
        }
  
        if (! m_sP_ROP_FC.check())  {
            setNLSStatus (m_sP_ROP_FC.getNLSStatus(), "m_sP_ROP_FC.check");
            return ( false );
        }

    }



    // * -------------------------------------------------------------------------- *
    // * Read out prephasing gradient                                               *
    // * -------------------------------------------------------------------------- *
    if (! m_sP_ROP.prepAmplitude(m_sP_ROP.getRampUpTime(), m_sP_ROP.getDuration(), m_sP_ROP.getRampDownTime(), m_sP_ROP.getAmplitude() ))  {
        setNLSStatus (m_sP_ROP.getNLSStatus(), "m_sP_ROP.prepAmplitude");
        return ( false );
    }
  
    if (! m_sP_ROP.check())  {
        setNLSStatus (m_sP_ROP.getNLSStatus(), "m_sP_ROP.check");
        return ( false );
    }

    return ( true );


}

//	*****************************************************
//	*                                                   *
//	* Prepares and checks the read out gradients        *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC00EE
bool SeqBuildBlockGREKernel::prepAndCheckRO (MrProt* pMrProt, SeqLim* pSeqLim)
{
  // **********************************************************************************
  // * Name        : prepAndCheckRO                                                   *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Prepares and checks the read out gradient.                       *
  // *               Returns true if prep and check has been successfully been        *
  // *               finished otherwise false.                                        *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  static const char *ptModule = "SeqBuildBlockGREKernel::prepAndCheckRO";

  long lI = 0;
  

  
  // * ---------------------------------------------------------------------- *
  // * Prepare and check read out gradient                                    *
  // * ---------------------------------------------------------------------- *
  for ( lI=0; lI < m_lNumberOfEchoes-1; lI++ )  {
    if (! m_sP_RO[lI].prepAmplitude(m_sP_RO[lI].getRampUpTime(), m_sP_RO[lI].getDuration(), m_sP_RO[lI].getRampDownTime(), m_sP_RO[lI].getAmplitude() ))  { /*! EGA-03; EGA-01 !*/
      setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].prepAmplitude");
      return ( false );
    }

    if (! m_sP_RO[lI].check())  {
      if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
        setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].check");  return ( false );
      } else {
        setNLSStatus (m_sP_RO[lI].getNLSStatus());  return ( false );
      }
    }
  }


  lI = m_lNumberOfEchoes-1;  // * Last echo *
  if (! m_sP_RO[lI].prepAmplitude(m_sP_RO[lI].getRampUpTime(), m_sP_RO[lI].getDuration(), m_sP_RO[lI].getRampDownTime(), m_sP_RO[lI].getAmplitude() ))  {   /*! EGA-03; EGA-01 !*/
    setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].prepAmplitude");
    return ( false );
  }

  if ( pMrProt->phaseStabilize() ) {

    if (! m_sP_RO[lI].check() )  {
      if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
        setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].check");   return ( false );
      } else {
        setNLSStatus (m_sP_RO[lI].getNLSStatus());   return ( false );
      }
    }

  } else {

    // * Check for allowed ramp down time of the read out gradient pulse
    switch(m_eReadOutGradType)  {
      case Spoiler:
          if (! m_sP_RO[lI].check( 0.0, m_sP_ROD.getAmplitude() ) )  {
              if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
                  setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].check");   return ( false );
              } else {
                  setNLSStatus (m_sP_RO[lI].getNLSStatus());   return ( false );
              }
          }
      break;

      case Constant:
          if (! m_sP_RO[lI].check( 0.0, m_sP_ROD.getAmplitude() ) )  {
              if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
                  setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].check");  return ( false );
              } else {
                  setNLSStatus (m_sP_RO[lI].getNLSStatus());  return ( false );
              }
          }
      break;

      case Rewinder:
      case Symmetric_Rewinder:
          if (! m_sP_RO[lI].check() )  {
              if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
                  setNLSStatus (m_sP_RO[lI].getNLSStatus(), "m_sP_RO[].check");  return ( false );
              } else {
                  setNLSStatus (m_sP_RO[lI].getNLSStatus());  return ( false );
              }
          }
      break;
    
      default: 
        if ( setNLSStatus ( SEQU_ERROR, ptModule, "Unknown read out gradient type." ) )  return ( false );
    }

  }

  
  return ( true );

}

//	Can be used to perform additional preparations
//##ModelId=3D352EBC008C
bool SeqBuildBlockGREKernel::prepAndCheckAddOn (MrProt* , SeqLim* , SeqExpo* )
{
  // **********************************************************************************
  // * Name        : run                                                              *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Additional preparations. Not used in GREKernel.                  *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Run function                                      *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EC00194
bool SeqBuildBlockGREKernel::run (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, sSLICE_POS* pSLC)
{
  // **********************************************************************************
  // * Name        : run                                                              *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Executes the basic GRE kernel                                    *
  // *                                                                                *
  // * Return      : An NLS status code                                               *
  // *                                                                                *
  // **********************************************************************************


  static const char *ptModule = "SeqBuildBlockGREKernel::run";


  long         lI;                         // * Loop counter *
  
  long         lGradLineNumber      = 0;   // * Line/partition number for gradient    *
  long         lGradPartNumber      = 0;   // * table amplitude calcilation           *
                                           // * lGradLineNumber = 0 -> k-space center *

  // * ---------------------------------------------------------------------- *
  // * Calculate line and partition number for gradient table amplitude       *
  // * calculation                                                            *
  // * ---------------------------------------------------------------------- *
  lGradLineNumber = m_lLineNumber - m_lKSpaceCenterLine;
  lGradPartNumber = m_lPartNumber - m_lKSpaceCenterPartition;



  // * ---------------------------------------------------------------------- *
  // * Set the Mdh of ADC[1] with data from ADC[0]                            *
  // * ---------------------------------------------------------------------- *
  if ( m_lNumberOfEchoes == 2 )  {
    m_ADC[1].Mdh = m_ADC[0].Mdh;
  }



  // * ---------------------------------------------------------------------- *
  // * Move LastScanInConcat and LastScanInMeas Mdh flags from the first to   *
  // * the second ADC.                                                        *
  // * ---------------------------------------------------------------------- *
  if ( pMrProt->contrasts() == 2 )  {
      if ( m_ADC[0].Mdh.isLastScanInConcat() )  {
          m_ADC[0].Mdh.deleteFromEvalInfoMask ( MDH_LASTSCANINCONCAT );
          m_ADC[1].Mdh.addToEvalInfoMask ( MDH_LASTSCANINCONCAT );
      }


      if ( m_ADC[0].Mdh.isLastScanInMeas() )  {
          m_ADC[0].Mdh.deleteFromEvalInfoMask ( MDH_LASTSCANINMEAS );
          m_ADC[1].Mdh.addToEvalInfoMask ( MDH_LASTSCANINMEAS );
      }
  }


  // * ---------------------------------------------------------------------- *
  // * Remove the phase and/or slice FFT flags from the second contrast ADC   *
  // * if phase stabilization is selected                                     *
  // * ---------------------------------------------------------------------- *
  if ( pMrProt->phaseStabilize() && ( pMrProt->contrasts() == 2 ) )  {
    m_ADC[1].Mdh.deleteFromEvalInfoMask ( MDH_PHASEFFT );
    m_ADC[1].Mdh.deleteFromEvalInfoMask ( MDH_D3FFT    );
  }




  // * ---------------------------------------------------------------------- *
  // * If last scan in concat or meas, the corresponding flags have to set    *
  // * for the phase stabiliazation ADC.                                      *
  // * ---------------------------------------------------------------------- *
  m_PhaseStabScan.setbLastScanInConcat( m_bLastScanInConcat );
  m_PhaseStabScan.setbLastScanInMeas  ( m_bLastScanInMeas   );




  // * --------------------------------------------------------------------- *
  // * General function to set Mdh-flags                                     *
  // * --------------------------------------------------------------------- *
  runSetMdh();



  for ( lI=0; lI<m_lNumberOfEchoes; lI++ )  {


    // * --------------------------------------------------------------------- *
    // * Activate online ice-process                                           *
    // * --------------------------------------------------------------------- *
    m_ADC[lI].Mdh.addToEvalInfoMask ( MDH_ONLINE );
    
    // * ---------------------------------------------------------------------- *
    // * Set the line, partition, echo and set parameter                        *
    // * ---------------------------------------------------------------------- *
    m_ADC[lI].Mdh.setClin ( m_lLineNumber ); // * Clin and Cpar have to be set  *
    m_ADC[lI].Mdh.setCpar ( m_lPartNumber ); // * outside of SEQLoop to ensure  *
                                             // * correct Mdh values in case    *
                                             // * reordering and elliptical     *
                                             // * scanning.                     *


    runSetMdhCeco ( lI );
    runSetMdhCset ( lI );
    


    // * ---------------------------------------------------------------------- *
    // * Set setKSpaceCentreColumn                                              *
    // * ---------------------------------------------------------------------- *
    m_ADC[lI].Mdh.setKSpaceCentreColumn ( (unsigned short)( m_adRelCenterCol[lI] * pMrProt->kSpace().baseResolution() + 0.5 ));

    // * ---------------------------------------------------------------------- *
    // * Set the cut off information                                            *
    // * ---------------------------------------------------------------------- *
    runSetMdhCutOff ( pMrProt, pSeqLim, lI );
    


    // * ---------------------------------------------------------------------- *
    // * Specify k-space center line and partition numbers                      *
    // * ---------------------------------------------------------------------- *
    m_ADC[lI].Mdh.setKSpaceCentreLineNo      ( m_lKSpaceCenterLine      );
    m_ADC[lI].Mdh.setKSpaceCentrePartitionNo ( m_lKSpaceCenterPartition );




    // * -------------------------------------------------------------------------- *
    // * Define the time(s) between RF and ADCs for phase stabilization             *
    // * -------------------------------------------------------------------------- *
    if ( pMrProt->phaseStabilize() )  {
        m_ADC[lI].Mdh.setTimeSinceLastRF ( pMrProt->te()[lI] );
    }




    // * ---------------------------------------------------------------------- *
    // * Set the frequency and phase properties of the rf pulse                 *
    // * Offcenter shift in read and slice directions are handeled              *
    // * by phase increments of the ADC                                         *
    // * ---------------------------------------------------------------------- *
    m_ADCSet[lI].prepSet( *pSLC, m_ADC[lI], m_sP_RO[lI], lGradLineNumber, lGradPartNumber );    /*! EGA-05 !*/
    m_ADCNeg[lI].prepNeg( *pSLC, m_ADC[lI], m_sP_RO[lI], lGradLineNumber, lGradPartNumber );    /*! EGA-05 !*/
  }



  // * ---------------------------------------------------------------------- *
  // * Set the frequency and phase properties of the rf pulse and ADC         *
  // * ---------------------------------------------------------------------- *
  if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL or m_pRF_Exc == NULL \n");
        return ( false );
  }
  m_pRF_WE->setAdditionalPhase ( m_dRFSpoilPhase );

  m_ADCSet[0].increasePhase ( m_dRFSpoilPhase );   // *   RF spoiling   *
  m_ADCNeg[0].decreasePhase ( m_dRFSpoilPhase );



  if ( m_lNumberOfEchoes == 2 )  {
    m_ADCSet[1].increasePhase( m_dRFSpoilPhase ); // *   RF spoiling   *
    m_ADCNeg[1].decreasePhase( m_dRFSpoilPhase );
  }




  // * ---------------------------------------------------------------------- *
  // * Perform phase cycle                                                    *
  // * ---------------------------------------------------------------------- *
  runPhaseCycle ();




  // * ---------------------------------------------------------------------- *
  // * Calculate phase encode gradient amplitudes                             *
  // * ---------------------------------------------------------------------- *
  if (! m_sTB_PE.prepPE (pMrProt, lGradLineNumber ) )   { setNLSStatus (m_sTB_PE.getNLSStatus(), "m_sTB_PE.prepPE"); return ( false ); }
  if (! m_sTB_3D.prep3D (pMrProt, lGradPartNumber ) )   { setNLSStatus (m_sTB_3D.getNLSStatus(), "m_sTB_3D.prep3D"); return ( false ); }


  if ( m_bUseRewinder == true ) {     // * Prepare rewinder amplitudes             *
    if (! m_sTB_PER.prepPE(pMrProt, lGradLineNumber ) )    { setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE"); return ( false ); }
    if (! m_sTB_3DR.prep3D(pMrProt, lGradPartNumber ) )    { setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D"); return ( false ); }
  } else {                       // * No rewinder -> amplitude equal zero !!! *
    if (! m_sTB_PER.prepAmplitude (0.0) )    { setNLSStatus (m_sTB_PER.getNLSStatus(), "m_sTB_PER.prepPE"); return ( false ); }
    if (! m_sTB_3DR.prepAmplitude (0.0) )    { setNLSStatus (m_sTB_3DR.getNLSStatus(), "m_sTB_3DR.prep3D"); return ( false ); }
  }




  // * ---------------------------------------------------------------------- *
  // * Execute the kernel                                                     *
  // * ---------------------------------------------------------------------- *
  if (! runTiming ( pMrProt, pSeqLim, pSeqExpo, pSLC, lGradLineNumber, lGradPartNumber ) )
    return ( false );
  else
    return ( true  );


}

//	Set the Mdh flags
//##ModelId=3D352EBC008A
void SeqBuildBlockGREKernel::runSetMdh ()
{
    // **********************************************************************************
    // * Name        : runSetMdh                                                        *
    // *                                                                                *
    // * Class       : SeqBuildBlockGREKernel                                           *
    // *                                                                                *
    // * Description : General function to set Mdh-flags                                *
    // *                                                                                *
    // * Return      : void                                                             *
    // *                                                                                *
    // **********************************************************************************


    // * -------------------------------------------------------------------------- *
    // * Sets the Mdh reflect flags depending on the contrast number                *
    // * -------------------------------------------------------------------------- *
    if ( m_eROPolarity == Positive )  {
        m_ADC[0].Mdh.deleteFromEvalInfoMask ( MDH_REFLECT );      /*! EGA-Any !*/
    } else {
        m_ADC[0].Mdh.addToEvalInfoMask ( MDH_REFLECT );           /*! EGA-Any !*/
    }
    

    if ( m_lNumberOfEchoes == 2 )  {
        if ( m_eROPolarity == Positive )  {
            m_ADC[1].Mdh.addToEvalInfoMask ( MDH_REFLECT );       /*! EGA-Any !*/
        } else {
            m_ADC[1].Mdh.deleteFromEvalInfoMask ( MDH_REFLECT );  /*! EGA-Any !*/
        }
    }


}

//	*****************************************************
//	*                                                   *
//	* Set the Mdh Ceco index                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC0062
void SeqBuildBlockGREKernel::runSetMdhCeco (long Ceco)
{
  // **********************************************************************************
  // * Name        : runSetMdhCeco                                                    *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Sets the Mdh Ceco index                                          *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  
  m_ADC[Ceco].Mdh.setCeco  ( Ceco );
  
}

//	*****************************************************
//	*                                                   *
//	* Set the Mdh Cset index                            *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBC002F
void SeqBuildBlockGREKernel::runSetMdhCset (long Cset)
{
  // **********************************************************************************
  // * Name        : runSetMdhCset                                                    *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : -                                                                *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************

  //lint -e{550}
  if( false ) { long lDummy; lDummy = Cset; }

}

//	*****************************************************
//	*                                                   *
//	* Set the Mdh pre and post cut off positions        *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBB03A1
void SeqBuildBlockGREKernel::runSetMdhCutOff (MrProt* pMrProt, SeqLim* pSeqLim, long lIndex)
{
  // **********************************************************************************
  // * Name        : runSetMdhCutOff                                                  *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Sets the Mdh pre and post cut off positions                      *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************
  

  m_ADC[lIndex].Mdh.setPreCutOff  ( static_cast<unsigned int>(m_lADCZeroBefore[lIndex] /                                    // * ADCZeroBefore is given in us *
                                    pMrProt->rxSpec().effDwellTime( pSeqLim->getReadoutOSFactor() )[lIndex] * 1000 ));      // * Dwell time is given in ns    *
  
  m_ADC[lIndex].Mdh.setPostCutOff ( static_cast<unsigned int>(m_lADCZeroAfter[lIndex] /                                     // * ADCZeroAfter is given in us  *
                                    pMrProt->rxSpec().effDwellTime( pSeqLim->getReadoutOSFactor() )[lIndex] * 1000 ) + 1 ); // * Dwell time is given in ns    *


}

//	This function can be used to perform phase cycles
//##ModelId=3D352EBB039F
void SeqBuildBlockGREKernel::runPhaseCycle ()
{
}

//##ModelId=3D352EBB024B
bool SeqBuildBlockGREKernel::runTiming (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, sSLICE_POS* pSLC, long lGradLineNumber, long lGradPartNumber)
{
  // **********************************************************************************
  // * Name        : runTiming                                                        *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Executes the GRE kernel                                          *
  // *                                                                                *
  // * Return      : BOOL                                                             *
  // *                                                                                *
  // **********************************************************************************

    static const char *ptModule = "SeqBuildBlockGREKernel::runTiming";

    long lT             = 0;
    long lStatus        = SEQU__NORMAL;
    long lOffsetTime    = 0;


    //lint -e{550}
    if( false ) { long lDummy; lDummy = lGradLineNumber = lGradPartNumber; }


    if ( m_pRF_WE == NULL ) {
        setNLSStatus(SEQU_ERROR, ptModule, "m_pRF_WE == NULL\n");
        return ( false );
    }



    // * ------------------------------------------------------------------------------------------------------------------- *
    // * |             Start Time        |     NCO     |    SRF   |    ADC   |          Gradient Events          | Sync    | *
    // * |               (usec)          |    Event    |   Event  |   Event  |   phase   |   read    |   slice   | Event   | *
    // * ------------------------------------------------------------------------------------------------------------------- *
    fRTEBInit( &pSLC->m_sROT_MATRIX);
    m_bUpdateRotationMatrix = false;




    // * ---------------------------------------------------------------------------- *
    // * Execute some additional timing elements                                      *
    // * The duration of these timining elements has to be specified using            *
    // * lOffsetTime                                                                  *
    // * ---------------------------------------------------------------------------- *
    if ( ! runAdditionalTiming ( pMrProt, pSeqLim, pSeqExpo, pSLC, &lOffsetTime ) )  {
        setNLSStatus(SEQU_ERROR, ptModule, "Error during execution of GREKernel::runAdditionalTiming().");
        return ( false );
    }

    lT += lOffsetTime;




    // * ---------------------------------------------------------------------------- *
    // *     Execute the osc bit                                                      *
    // * ---------------------------------------------------------------------------- *
    if ( m_bSendOscBit )  {
        fRTEI(lT                             ,          0,         0,         0,          0,          0,          0,&m_sOscBit);
    }
  
  
    // * ---------------------------------------------------------------------------- *
    // *     Execute the wake up bit                                                  *
    // * ---------------------------------------------------------------------------- *
    if ( m_bSendWakeUpBit )  {
        fRTEI(lOffsetTime+m_lStartTimeWakeUp ,          0,         0,         0,          0,          0,          0,&m_sWakeUp);
    }

  
    // * ---------------------------------------------------------------------------- *
    // *     Execute the water selective excitation pulse                             *
    // * ---------------------------------------------------------------------------- *
    m_pRF_WE->setStartTimeInEventBlock ( lOffsetTime+m_pRF_WE->getRequiredLeadTime() );
    if (! m_pRF_WE->run(pMrProt, pSeqLim, pSeqExpo, pSLC) )  {    
        setNLSStatus(m_pRF_WE->getNLSStatus(), ptModule, "Error in m_pRF_WE->run(...) at runtime");
        return ( false );
    }



    // * ---------------------------------------------------------------------------- *
    // *   Phase encode and slice select events                                       *
    // * ---------------------------------------------------------------------------- *
    lT = lOffsetTime + m_lExcitationAllTime + m_lDelayAdjustSSPosition;

    if ( m_bFlowCompSlice )  {
        fRTEI(lT + m_pRF_WE->getRampTimeOutsideSBB(), 0,         0,         0,          0,          0, &m_sP_3D_FC,        0);
    }

    fRTEI(lT + 2 * m_alTEFill[0] +
          m_lPrephaseFillTimePhase       ,            0,         0,         0,  &m_sTB_PE,          0,           0,        0);

    fRTEI(lT + m_sP_3D_FC.getTotalTime() +
               m_pRF_WE->getRampTimeOutsideSBB(),     0,         0,         0,          0,          0,   &m_sTB_3D,        0);

    lT += m_lPrephaseTime + 2 * m_alTEFill[0] + m_lReadoutTime + 2 * m_alTEFill[1] + m_lDiffBetweenROAndADCTiming;


    fRTEI(lT                             ,            0,         0,         0, &m_sTB_PER,          0,  &m_sTB_3DR,        0);



    lT += m_lRephaseTime + m_lPhaseEventDelay + m_lTRFill - m_lOverTime;
    if ( pMrProt->phaseStabilize()           )  {    lT+=m_PhaseStabScan.getDurationPerRequest() + m_sP_ROD.getTotalTime();    };

    fRTEI(lT                             ,            0,         0,         0,          0,          0,           0,        0);

    if (m_lTRFillEnd>0)  {
        fRTEI(lT+=m_lTRFillEnd           ,            0,         0,         0,          0,          0,           0,        0);
    }




    // * ---------------------------------------------------------------------------- *
    // *     lT set to zero                                                           *
    // *     read out events                                                          *
    // * ---------------------------------------------------------------------------- *
    lT = lOffsetTime + m_lExcitationAllTime + m_lDelayAdjustSSPosition;
    if ( m_bFlowCompRead )  {
        fRTEI(lT + 2 * m_alTEFill[0] +
              m_lPrephaseFillTimeRead ,               0,         0,         0,          0,&m_sP_ROP_FC,          0,        0);
    }

    fRTEI(lT + m_sP_ROP_FC.getTotalTime() +
               2 * m_alTEFill[0] +
               m_lPrephaseFillTimeRead   ,            0,         0,         0,          0,  &m_sP_ROP,           0,        0);
    lT += m_lPrephaseTime + 2 * m_alTEFill[0];
    fRTEI(lT - m_sP_RO[0].getRampUpTime(),            0,         0,         0,          0,&m_sP_RO[0],           0,        0);

    if ( m_lNumberOfEchoes == 2 )  {
        fRTEI(lT - m_sP_RO[0].getRampUpTime() 
                 + m_sP_RO[0].getTotalTime()
                 + m_lDiffBetweenROAndADCTiming
                 + 2 * m_alTEFill[1]       ,          0,         0,         0,          0,&m_sP_RO[1],           0,        0);
    }


    if ( pMrProt->phaseStabilize() )  {

        lT += m_lReadoutTime + 2 * m_alTEFill[1] + m_lDiffBetweenROAndADCTiming + m_lRephaseTime;
        m_PhaseStabScan.setlStartTime ( lT );
        m_PhaseStabScan.setlTimeFromRFToStartSBB ( m_PhaseStabScan.getlStartTime() - m_lExcitationAllTime + m_lTEFromExcitation );
        m_PhaseStabScan.setdRFSpoilPhase ( m_dRFSpoilPhase );
        m_PhaseStabScan.run ( pMrProt, pSeqLim, pSeqExpo, pSLC );

        lT += m_PhaseStabScan.getDurationPerRequest();

    } else {

        lT += m_lReadoutTime + 2 * m_alTEFill[1] + m_lDiffBetweenROAndADCTiming;
        if ( (m_eReadOutGradType == Rewinder) || (m_eReadOutGradType == Symmetric_Rewinder) ) {  lT += m_sP_RO[0].getRampDownTime();  }

    }



    fRTEI(lT                               ,          0,         0,         0,          0,  &m_sP_ROD,           0,        0);






    // * ---------------------------------------------------------------------------- *
    // *     lT set to zero                                                           *
    // *     ADC events                                                               *
    // * ---------------------------------------------------------------------------- *
    lT = lOffsetTime + m_lExcitationAllTime + m_lDelayAdjustSSPosition;

    if ( m_eReadOutGradType == Symmetric_Rewinder ) {
        fRTEI(lT+= m_lDelayAdjustADCPosition,         0,         0,         0,          0,          0,           0,        0);
    }

    fRTEI(lT+= m_lTimeToADC[0]   +
               2 * m_alTEFill[0]         , &m_ADCSet[0],         0, &m_ADC[0],          0,          0,           0,        0);

    fRTEI(lT+= m_ADC[0].getRoundedDuration(),&m_ADCNeg[0],       0,         0,          0,          0,           0,        0);

    if ( m_lNumberOfEchoes == 2 )  {
        fRTEI(lT+= m_lTimeToADC[1]         , &m_ADCSet[1],       0, &m_ADC[1],          0,          0,           0,        0);
        fRTEI(lT+= m_ADC[1].getRoundedDuration(),&m_ADCNeg[1],   0,         0,          0,          0,           0,        0);
    }


    #ifndef VXWORKS
        mSEQTest(pMrProt, pSeqLim, pSeqExpo, m_ulUTIdent , 30, lGradLineNumber, pSLC->getSliceIndex(), 0, lGradPartNumber); /*! EGA-All !*/
    #endif

    lStatus = fRTEBFinish();

    if ( setNLSStatus( lStatus ) )  {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s : fRTEBFinish() [*0030*] returned with an error : NLS_CODE : 0x%lx"  , ptModule, lStatus);
        return ( false );
    }


    return ( true );

}

//	This function can be used to execute additional timing 
//prior to the kernel.
//	Duration of this timing has to be specified in us using 
//the argument lOffset.
//##ModelId=3D352EBB0132
bool SeqBuildBlockGREKernel::runAdditionalTiming (MrProt* , SeqLim* , SeqExpo* , sSLICE_POS* , long* plDuration)
{

    *plDuration = 0;

    return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Just returning true                        *
//	* In SBBTRUFIKernel: Calculation of the *
//	* prescan kernel for TrueFisp    *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBF03E1
bool SeqBuildBlockGREKernel::calculatePreScanTiming (MrProt* , SeqLim* , SeqExpo* )
{
  // **********************************************************************************
  // * Name        : calculatePreScanTiming                                           *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Can be used to calculate the timing of the prescan for TrueFisp  *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

  return ( true );  // Returning without error
  
}

//	*****************************************************
//	*                                                   *
//	* Preparation of the prescan kernel for TrueFisp    *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBF03C3
bool SeqBuildBlockGREKernel::prepPreScan (SeqLim* )
{
  // **********************************************************************************
  // * Name        : prepPreScan                                                      *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Can be used to trepare and check the gardients used for the      *
  // *               TrueFisp prescan                                                 *
  // *                                                                                *
  // * Return      : bool                                                             *
  // *                                                                                *
  // **********************************************************************************

  return ( true);
  
}

//	*****************************************************
//	*                                                   *
//	* Run function of the prescan for TrueFisp          *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBF0355
bool SeqBuildBlockGREKernel::runPreScan (MrProt* , SeqLim* , SeqExpo* , sSLICE_POS* )
{
  // **********************************************************************************
  // * Name        : run                                                              *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Can be used to execute the prescan for TrueFisp                  *
  // *                                                                                *
  // * Return      : An NLS status code                                               *
  // *                                                                                *
  // **********************************************************************************

  return ( true );

}

//	*****************************************************
//	*                                                   *
//	* Print all imported parameters                     *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBF02E7
void SeqBuildBlockGREKernel::printImports (MrProt* pMrProt, SeqLim* pSeqLim, const char* ptModule)
{
  // **********************************************************************************
  // * Name        : printImports                                                     *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Print all imported parameters                                    *
  // *                                                                                *
  // * Return      : void                                                             *
  // *                                                                                *
  // **********************************************************************************

  long   lI;    // * Loop counter *
  
  


  if ( mIsTRinp )  {
    TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Imported Parameter                                *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);

    for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_adRelWinWidthCol[%1ld]     = %11.6f   ", ptModule, lI, m_adRelWinWidthCol[lI]);
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_adRelWinAsymCol[%1ld]      = %11.6f   ", ptModule, lI, m_adRelWinAsymCol[lI] );
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_adRelCenterCol[%1ld]       = %11.6f\n ", ptModule, lI, m_adRelCenterCol[lI]  );
    }

    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_FirstLineToMeasure      = %8ld \n", ptModule, m_lFirstLineToMeasure  );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lLastLineToMeasure      = %8ld \n", ptModule, m_lLastLineToMeasure   );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lFirstPartToMeasure     = %8ld \n", ptModule, m_lFirstPartToMeasure  );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lLastPartToMeasure      = %8ld \n", ptModule, m_lLastPartToMeasure   );
    if (m_bBalancedGradMomentSliceSelect)  {
      TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_bBalancedGradMomentSliceSelect = balanced \n", ptModule);
    }  else {
      TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_bBalancedGradMomentSliceSelect = not balanced \n", ptModule);
    }

    switch ( m_eReadOutGradType )  {
      case Spoiler:            TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Spoiler \n"            , ptModule);      break;
      case Constant:           TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Constant\n"            , ptModule);      break;
      case Rewinder:           TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Rewinder\n"            , ptModule);      break;
      case Symmetric_Rewinder: TRACE_PUT1(TC_INFO, TF_SEQ, "%s m_ReadOutGradType         = Symmetric_Rewinder\n"  , ptModule);      break;
      default:                 TRACE_PUT1(TC_INFO, TF_SEQ, "%s: Unknown read out gradient type \n"                , ptModule);
    }


    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_Exc->getDuration()  = %8ld \n", ptModule, m_pRF_Exc->getDuration());

    for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  { 
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s MrProt->te()[%1ld]           = %8ld [us]\n", ptModule, lI, pMrProt->te()[lI]  );
    }
    TRACE_PUT2(TC_INFO, TF_SEQ,   "%s MrProt->tr()[0]           = %8ld [us]\n", ptModule, pMrProt->tr()[0]  );
    TRACE_PUT2(TC_INFO, TF_SEQ,   "%s m_lTRFill                 = %8ld [us]\n", ptModule, m_lTRFill         );
    TRACE_PUT1(TC_INFO, TF_SEQ,   "%s * ------------------------------------------------------------------ *\n" ,ptModule);
  }
  
}

//	*****************************************************
//	*                                                   *
//	* Print all results                                 *
//	*                                                   *
//	*****************************************************
//##ModelId=3D352EBF0297
void SeqBuildBlockGREKernel::printResults (SeqLim* pSeqLim, const char* ptModule)
{
  // **********************************************************************************
  // * Name        : printResults                                                     *
  // *                                                                                *
  // * Class       : SeqBuildBlockGREKernel                                           *
  // *                                                                                *
  // * Description : Print all results of the kernel calculation                      *
  // *                                                                                *
  // * Return      : void                                                             *
  // *                                                                                *
  // **********************************************************************************

  long   lI;    // * Loop counter *
  
  
  if ( mIsTRres )  {
    TRACE_PUT0(TC_INFO, TF_SEQ, "\n\n");
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s *                  Exported Parameter                                *\n" ,ptModule);
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lExcitationAllTime          = %6ld [us]\n"        , ptModule, m_lExcitationAllTime         );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseTime               = %6ld [us]\n"        , ptModule, m_lPrephaseTime              );
    for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_lADCDatTime[%1ld]              = %6ld [us]\n"   , ptModule, lI, m_lADCDatTime[lI] );
    }
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lRephaseTime                = %6ld [us]\n"        , ptModule, m_lRephaseTime               );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lOverTime                   = %6ld [us]\n"        , ptModule, m_lOverTime                  );
    for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_alTEMin[%1ld]                  = %6ld [us]\n"   , ptModule, lI, m_alTEMin[lI]            );
    }
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lTRMin                      = %6ld [us]\n"        , ptModule, m_lTRMin                     );
    for ( lI=0; lI < m_lNumberOfEchoes; lI++ )  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_alTEFill[%1ld]                 = %6ld [us]\n\n" , ptModule, lI, m_alTEFill[lI]           );
    }
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lStartTimeWakeUp            = %6ld [us]\n"        , ptModule, m_lStartTimeWakeUp           );

    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_WE->getGSData().getdFirstAmplitude()   = %11.6f[mT/m]\n"    , ptModule, m_pRF_WE->getGSData().getdFirstAmplitude()   );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_pRF_WE->getGSData().getdRequiredRefocusingMoment() = %11.6f[mT/m]\n\n"  , ptModule, m_pRF_WE->getGSData().getdRequiredRefocusingMoment() );


    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseFillTimeSlice      = %6ld [us]\n"      , ptModule, m_lPrephaseFillTimeSlice     );
    if ( m_eLimitingPrephaseTime == Slice ) {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampUpTime()      = %6ld [us] *\n"      , ptModule, m_sTB_3D.getRampUpTime()     );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getDuration()        = %6ld [us] *\n"      , ptModule, m_sTB_3D.getDuration()       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampDownTime()    = %6ld [us] *\n"      , ptModule, m_sTB_3D.getRampDownTime()   );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dMomentum_TB_3D             = %11.6f[mT/m*us]\n\n", ptModule, m_dMomentum_TB_3D            );
    } else {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampUpTime()      = %6ld [us]\n"        , ptModule, m_sTB_3D.getRampUpTime()     );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getDuration()        = %6ld [us]\n"        , ptModule, m_sTB_3D.getDuration()       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3D.getRampDownTime()    = %6ld [us]\n"        , ptModule, m_sTB_3D.getRampDownTime()   );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dMomentum_TB_3D             = %11.6f[mT/m*us]\n\n", ptModule, m_dMomentum_TB_3D            );
    }

    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3DR.getRampUpTime()     = %6ld [us]\n"        , ptModule, m_sTB_3DR.getRampUpTime()    );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3DR.getDuration()       = %6ld [us]\n"        , ptModule, m_sTB_3DR.getDuration()      );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_3DR.getRampDownTime()   = %6ld [us]\n\n"      , ptModule, m_sTB_3DR.getRampDownTime()  );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_dMomentum_TB_3DR            = %11.6f[mT/m*us]\n\n", ptModule, m_dMomentum_TB_3DR           );
    
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseFillTimePhase      = %6ld [us]\n"      , ptModule, m_lPrephaseFillTimePhase     );
    if ( m_eLimitingPrephaseTime == Phase ) {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampUpTime()      = %6ld [us] *\n"  , ptModule, m_sTB_PE.getRampUpTime()     );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getDuration()        = %6ld [us] *\n"  , ptModule, m_sTB_PE.getDuration()       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampDownTime()    = %6ld [us] *\n\n", ptModule, m_sTB_PE.getRampDownTime()   );
    } else {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampUpTime()      = %6ld [us]\n"    , ptModule, m_sTB_PE.getRampUpTime()     );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getDuration()        = %6ld [us]\n"    , ptModule, m_sTB_PE.getDuration()       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PE.getRampDownTime()    = %6ld [us]\n\n"  , ptModule, m_sTB_PE.getRampDownTime()   );
    }

    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PER.getRampUpTime()     = %6ld [us]\n"        , ptModule, m_sTB_PER.getRampUpTime()    );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PER.getDuration()       = %6ld [us]\n"        , ptModule, m_sTB_PER.getDuration()      );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sTB_PER.getRampDownTime()   = %6ld [us]\n\n"      , ptModule, m_sTB_PER.getRampDownTime()  );
    
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_lPrephaseFillTimeRead       = %6ld [us]\n"      , ptModule, m_lPrephaseFillTimeRead      );
    if ( m_eLimitingPrephaseTime == Read ) {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampUpTime()      = %6ld [us] *\n"      , ptModule, m_sP_ROP.getRampUpTime()     );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getDuration()        = %6ld [us] *\n"      , ptModule, m_sP_ROP.getDuration()       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampDownTime()    = %6ld [us] *\n"      , ptModule, m_sP_ROP.getRampDownTime()   );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getAmplitude()       = %11.6f[mT/m]\n\n" , ptModule, m_sP_ROP.getAmplitude()      );

        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampUpTime()   = %6ld [us] *\n"     , ptModule, m_sP_ROP_FC.getRampUpTime()  );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getDuration()     = %6ld [us] *\n"     , ptModule, m_sP_ROP_FC.getDuration()    );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampDownTime() = %6ld [us] *\n"     , ptModule, m_sP_ROP_FC.getRampDownTime());
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getAmplitude()    = %11.6f[mT/m] *\n\n", ptModule, m_sP_ROP_FC.getAmplitude()   );
    } else {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampUpTime()      = %6ld [us]\n"        , ptModule, m_sP_ROP.getRampUpTime()     );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getDuration()        = %6ld [us]\n"        , ptModule, m_sP_ROP.getDuration()       );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getRampDownTime()    = %6ld [us]\n"        , ptModule, m_sP_ROP.getRampDownTime()   );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP.getAmplitude()       = %11.6f[mT/m]\n\n"   , ptModule, m_sP_ROP.getAmplitude()      );

        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampUpTime()   = %6ld [us]\n"      , ptModule, m_sP_ROP_FC.getRampUpTime()  );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getDuration()     = %6ld [us]\n"      , ptModule, m_sP_ROP_FC.getDuration()    );
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getRampDownTime() = %6ld [us]\n"      , ptModule, m_sP_ROP_FC.getRampDownTime());
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROP_FC.getAmplitude()    = %11.6f[mT/m]\n\n" , ptModule, m_sP_ROP_FC.getAmplitude()   );
    }

    for (lI=0; lI < m_lNumberOfEchoes; lI++ )  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_sP_RO[%1ld].getRampUpTime()    = %6ld [us]\n",      ptModule, lI, m_sP_RO[lI].getRampUpTime()  );
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_sP_RO[%1ld].getDuration()      = %6ld [us]\n",      ptModule, lI, m_sP_RO[lI].getDuration()    );
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_sP_RO[%1ld].getRampDownTime()  = %6ld [us]\n",      ptModule, lI, m_sP_RO[lI].getRampDownTime());
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_sP_RO[%1ld].getAmplitude()     = %11.6f[mT/m]\n\n", ptModule, lI, m_sP_RO[lI].getAmplitude()   );
    }

    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getRampUpTime()      = %6ld [us]\n"       , ptModule, m_sP_ROD.getRampUpTime()      );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getDuration()        = %6ld [us]\n"       , ptModule, m_sP_ROD.getDuration()        );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getRampDownTime()    = %6ld [us]\n"       , ptModule, m_sP_ROD.getRampDownTime()    );
    TRACE_PUT2(TC_INFO, TF_SEQ, "%s m_sP_ROD.getAmplitude()       = %11.6f[mT/m]\n\n"  , ptModule, m_sP_ROD.getAmplitude()       );

    for (lI=0; lI < m_lNumberOfEchoes; lI++ )  {
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_lTimeToADC[%1ld]               = %6ld [us]\n",       ptModule, lI, m_lTimeToADC[lI]      );
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_lADCZeroBefore[%1ld]           = %6ld [us]\n",       ptModule, lI, m_lADCZeroBefore[lI]  );
      TRACE_PUT3(TC_INFO, TF_SEQ, "%s m_lADCZeroAfter[%1ld]            = %6ld [us]\n",       ptModule, lI, m_lADCZeroAfter[lI]   );
    }
    TRACE_PUT1(TC_INFO, TF_SEQ, "%s * ------------------------------------------------------------------ *\n" ,ptModule);
  }  

}

//	The water excitation pulse may be used with an offset 
//frequency. This
//	function enables disables the use of the offset frequency 
//and specifies the
//	frequency in Hz. The offset frequency is an optional 
//parameter. Its default
//	value 0.0 Hz.
//##ModelId=3D352EBF0247
void SeqBuildBlockGREKernel::useWEOffsetFrequency (bool bUseOffsetFrequency, double dWEOffsetFrequency)
{
    m_bUseWEOffsetFrequency = bUseOffsetFrequency;
    
    if ( m_bUseWEOffsetFrequency ) {
    	m_dWEOffsetFrequency = dWEOffsetFrequency;
    } else {
    	m_dWEOffsetFrequency = 0.0;
    }


}

// Additional Declarations

bool SeqBuildBlockGREKernel::calculateInit (MrProt* pMrProt, SeqLim*, SeqExpo*)
{
    m_lReadoutTime = 0;

    m_lMinDelayBetweenADCs     = static_cast<long>(SysProperties::getMinDurationBetweenReadoutAndReadout(pMrProt->rxSpec().realDwellTime()[0] / 1000.0));
    m_lMinDelayBetweenADCAndRF = SysProperties::getMinDurationBetweenReadoutAndRFPulse();
    m_lMinDelayBetweenRFAndADC = SysProperties::getMinDurationBetweenRFPulseAndReadout();


    return ( true );

}







// * -------------------------------------------------------------------------- *
// *                   END OF OLD KERNEL DEFINITION                             *
// * -------------------------------------------------------------------------- *


