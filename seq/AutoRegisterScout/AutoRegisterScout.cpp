/*[ Compilation unit ********************************************************\
*
* Project: NUMARIS/4
*
* File        : \n4\pkg\MrServers\MrImaging\seq\a_gre.cpp
*
* Author      : Clinical
*
* Date        : n.a.
*
* Lang        : cpp
*
* EGA Requirement Key: As shown on the following lines:
*
*    Abbrev.   Translation                          Relevant for
*    -------   -----------                          ------------
*    EGA-All   All of the following keys:           All EGA requirements
*    EGA-01    000_EGA_BildOri_SW_SequenzROVz       GR/GP   polarity
*    EGA-02    000_EGA_BildPos_SW_SequenzSSelVz     GS      polarity
*    EGA-03    000_EGA_BildMass_SW_SequenzROPC      GR/GP   amplitude
*    EGA-04    000_EGA_BildPos_SW_SequenzSSel       GS      amplitude
*    EGA-05    000_EGA_BildPos_SW_NCOFrequenzSSel   SRF     frequency
*    EGA-06    000_EGA_BildPos_SW_NCOFrequenzRO     Readout frequency
*    EGA-07    000_EGA_BildOri_SW_OrientierungTest  Image orientation
*
*         Function(s): fSEQInit, fSEQPrep, fSEQCheck, fSEQRun, fSEQRunKernel
*
* Description : Variable gradient echo sequence (2D / 3D)
*
* Variants:
*    - gre
*    - gre_rt
*
\****************************************************************************/


/*] END: */

#define THICKNESSmm_RxGain      (10) /* Slice thick. limit for setting Rx gain */






#include "MrServers/MrMeasSrv/SeqFW/libSSL/libSSL.h"
#include "MrServers/MrProtSrv/MrProt/SeqDefines.h"
#include "MrServers/MrPerAns/PerProxies/GCProxy.h"
#include "MrServers/MrImaging/seq/Kernels/SBBGREKernel.h"
#include "MrServers/MrImaging/seq/AutoRegisterScout/LocalSeqLoop.h"
#include "MrServers/MrImaging/libSBB/SBBTSat.h"
#include "MrServers/MrImaging/libSBB/SBBLineTag.h"
#include "MrServers/MrImaging/libSBB/SBBGridTag.h"
#include "MrServers/MrImaging/libSeqUtil/RealTimeControl.h"
#include "MrServers/MrImaging/libSeqUtil/ReorderInfoGRE.h"
#include "MrServers/MrImaging/libSeqUtil/KernelCalculationLimits.h"
#include "MrServers/MrImaging/libUICtrl/UICtrl.h"
#include "MrServers/MrImaging/seq/SystemProperties.h "  // * Siemens system properties *
#include "MrServers/MrMeasSrv/SeqIF/Sequence/AutoPrepContext.h"
#include "MrServers\MrProtSrv\MrProt\ParametricMapping\ParametricMapping.h"


// Realtime includes

    #include "MrServers/MrPerAns/PerHostDevices/IntScan/IntScanMouseData/IntScanMouseData.h"
    #include "MrServers/MrImaging/libSeqUtil/fSUQuaternion.h"
    #include "MrServers/MrMeasSrv/SeqIF/Sequence/sequmsg.h"
    #include "MrServers/MrMeasSrv/SeqIF/libRT/SEQSemaphore.h"

#ifndef VXWORKS
    #include "MrServers/MrImaging/seq/AutoRegisterScout/AffineTransformReceiver.h"

    #include "MrServers/MrProtSrv/MrProtocol/libUILink/StdRoutines.h"
    #include "MrServers/MrProtSrv/MrProtocol/libUILink/UILinkLimited.h"
    #include "MrServers/MrProtSrv/MrProtocol/libUILink/UILinkSelection.h"
    #include "MrServers/MrProtSrv/MrProtocol/UILink/MrUILinkTR.h"
    #include "MrServers/MrProtSrv/MrProtocol/UILink/MrUILinkPhysio.h"
    #include "MrServers/MrProtSrv/MrProtocol/libUILink/UILinkArray.h"
    #include "MrServers/MrProtSrv/MrProtocol/UILink/StdProtRes/StdProtRes.h"
    #include "MrServers/MrProtSrv/MrProtocol/UILink/MrStdNameTags.h"
    #include "MrServers/MrImaging/libUICtrl/UICtrl.h"
    #include "MrServers/MrMeasSrv/SeqIF/Sequence/Sequence.h"
    #include <vector>
#endif



#ifdef SUPPORT_iPAT
    #include   "MrServers/MrImaging/seq/common/iPAT/iPAT.h"
#endif

#ifdef SUPPORT_CT
#include "MrServers/MrImaging/seq/common/CT/CT_UI.h"
#endif



// * ------------------------------------------------------------------------ *
// * Function prototypes                                                      *
// * ------------------------------------------------------------------------ *
static NLS_STATUS fSEQRunKernel                    (MrProt*, SeqLim*, SeqExpo*, long, long, long, long);
#ifdef SUPPORT_CT
//  SEQRunKernel used, if pProt->MDS().mdsModeMask() == SEQ::MDS_ONCO
static NLS_STATUS fSEQRunKernel_ct                 (MrProt*, SeqLim*, SeqExpo*, long, long, long, long);
#endif

NLS_STATUS        fCalculateReordering             (MrProt*, SeqLim*, long*, long*, long*, long*);
NLS_STATUS        fSelectAndPrepRFPulse            (MrProt*, SeqLim*, SeqExpo*, sRF_PULSE**);
NLS_STATUS        fSetKernelParameter              (MrProt*, SeqLim*, long, long, long, long);
NLS_STATUS        fCalculateGRETiming              (MrProt*, SeqLim*, SeqExpo*);
NLS_STATUS        fCalculateTEMin                  (MrProt*, SeqLim*, SeqExpo*, Sequence*, long*, long*);
NLS_STATUS        fCalculateNoOfQuickSatsPerConcat (MrProt*, long, long, long*, bool*);


#ifndef VXWORKS

    BOOL            pfGetMinTE                           (MrProt*, SeqLim*, SeqExpo*, Sequence*, long*);

    bool            segmentsGetLimitsHandler             (LINK_LONG_TYPE* const, std::vector<MrLimitLong>&, unsigned long&, long);
    bool            firstAcqWindowGetLimitsHandler       (LINK_LONG_TYPE* const, std::vector<MrLimitLong>&, unsigned long&, long);
    bool            firstAcqWindowInGetLimitsHandler     (LINK_LONG_TYPE* const, std::vector<MrLimitLong>&, unsigned long&, long);
    bool            trGetLimitsHandler                   (LINK_DOUBLE_TYPE* const, std::vector<MrLimitDouble>&, unsigned long&, long);

    bool            phasesIsAvailableHandler             (LINK_LONG_TYPE* const, long);

    unsigned        firstAcqWindowToolTipIdHandler       (LINK_LONG_TYPE* const, char**, long);

    long            contrastSetValueHandler              (LINK_LONG_TYPE* const, long, long);
    long            segmentsSetValueHandler              (LINK_LONG_TYPE* const, long, long);
    long            firstAcqWindowSetValueHandler        (LINK_LONG_TYPE* const, long, long);
    double          trSetValueHandler                    (LINK_DOUBLE_TYPE* const, double, long);
    unsigned        firstPhysioSignalModeSetValueHandler (LINK_SELECTION_TYPE* const, unsigned, long);

    unsigned        fSSolvePhaseStabilizeConflict        (LINK_BOOL_TYPE*      const, char**, const void*, const MrProt*, long);
    unsigned        fLSolveSegmentsConflict              (LINK_LONG_TYPE*      const, char**, const void*, const MrProt*, long);
    unsigned        fLSolveLongConflict                  (LINK_LONG_TYPE*      const, char**, const void*, const MrProt*, long);
    unsigned        fLSolveBaseResConflict               (LINK_LONG_TYPE*      const, char**, const void*, const MrProt*, long);
    unsigned        fSSolveSelectionConflict             (LINK_SELECTION_TYPE* const, char**, const void*, const MrProt*, long);
    unsigned        fSSolveInversionConflict             (LINK_SELECTION_TYPE* const, char**, const void*, const MrProt*, long);
    unsigned        fSSolveExcitationModeConflict        (LINK_SELECTION_TYPE* const, char**, const void*, const MrProt*, long);
    unsigned        fSSolveFirstSignalModeConflict       (LINK_SELECTION_TYPE* const, char**, const void*, const MrProt*, long);

    unsigned        fSSolveFlowCompConfilct              (LINK_SELECTION_TYPE* const, char**, const void*, const MrProt*, long);
    unsigned        fSSolveMSMModeConflict               (LINK_SELECTION_TYPE* const, char**, const void*, const MrProt*, long);

    template<class TYPE>
    bool _fLocalUICSetWithinLimits(MrUILinkBase* pThis, TYPE& rVal, const char* tTag, unsigned long ulMode, long lIndex, MrUILinkLimited<TYPE>::PFctSetValue pFctSet);
    bool fLocalUICSetWithinLimits(MrUILinkBase* const pThis, double& rVal, const char* tTag, unsigned long ulMode, long lIndex, MrUILinkLimited<double>::PFctSetValue pFctSet=NULL);
    bool fLocalUICSetWithinLimits(MrUILinkBase* const pThis, long& rVal, const char* tTag, unsigned long ulMode, long lIndex, MrUILinkLimited<long>::PFctSetValue pFctSet=NULL);

    template <class ParamType, class Param, class OrigFuncType>
    Param basicSetValueHandler ( ParamType, Param, long, OrigFuncType );
#endif



// * ------------------------------------------------------------------------ *
// * Control of variants                                                      *
// * ------------------------------------------------------------------------ *
enum SequenceType  {
        Gre, RT_Gre, Gre_HR
};




enum eEchoMode {
    Symmetric, Asymmetric
};



// * ------------------------------------------------------------------------ *
// * Actually compiled sequence                                               *
// * ------------------------------------------------------------------------ *
#ifdef GRE
    SequenceType  eSequence = Gre;
#endif

#ifdef GRE_RT
    SequenceType  eSequence = RT_Gre;
#endif

#ifdef GRE_HIGHRES
    SequenceType  eSequence = Gre_HR;
#endif



#ifdef SUPPORT_CT
enum
{
    CSAT_MODE_AUTO             = 0
,   CSAT_MODE_MANUAL           = 1
,   GRAD_SPOIL_SCHEME_PRODUCT  = 0
,   GRAD_SPOIL_SCHEME_SIMPLE   = 1
};

int32_t CSAT_MODE(const MrProt* pProt)
{
    return pProt->MDS().alFree(0);
}
int32_t CSAT_MODE(MrProt* pProt, int32_t i32NewVal)
{
    return pProt->MDS().alFree(i32NewVal,0);
}

int32_t GRAD_SPOIL_SCHEME(const MrProt* pProt)
{
    return pProt->MDS().alFree(1);
}
int32_t GRAD_SPOIL_SCHEME(MrProt* pProt, int32_t i32NewVal)
{
    return pProt->MDS().alFree(i32NewVal,1);
}

double CSAT_FLIP_ANGLE_deg(const MrProt* pProt)
{
    return CSATFLIPANGLEdeg+pProt->data()->adFlipAngleDegree[K_NO_FLIP_ANGLE-1];
}
double CSAT_FLIP_ANGLE_deg(MrProt* pProt, double dNewVal)
{
    return CSATFLIPANGLEdeg+(pProt->data()->adFlipAngleDegree[K_NO_FLIP_ANGLE-1]=dNewVal-CSATFLIPANGLEdeg);
}

int32_t CSAT_SHIFT(const MrProt* pProt)
{
    return pProt->MDS().alFree(2);
}
int32_t CSAT_SHIFT(MrProt* pProt, int32_t i32NewVal)
{
    return pProt->MDS().alFree(i32NewVal,2);
}

static double T1_FAT_us()
{
    return SysProperties::getNominalBZero() > 2.8
        ? 300000
        : 270000
        ;
}

#endif //  SUPPORT_CT



/*********************************/
/* Global variables and arrays */
/*********************************/

static double         dRFSpoilPhase              = 0.0;                  // * For control of RF spoiling            *
static double         dRFSpoilPhasePrevSlice     = 0.0;
static double         dRFSpoilIncrement          = 0.0;
static double         dRFSpoilIncrementPrevSlice = 0.0;
static double         dRFPrevSlicePosSag         = 999999.0;             // * Comparing slice position              *
static double         dRFPrevSlicePosCor         = 999999.0;
static double         dRFPrevSlicePosTra         = 999999.0;
static double         dRFPrevSliceNormalSag      = 999999.0;
static double         dRFPrevSliceNormalCor      = 999999.0;
static double         dRFPrevSliceNormalTra      = 999999.0;
static double         dRFSPOIL_INCREMENTdeg      = RFSPOIL_INCREMENTdeg; /* Fixed phase increment  */




static long           lKSpaceCenterLine          = 0;     // * Number of the echo (ky = 0)              //
static long           lLinesToMeasure            = 0;     // * Total number of measured lines with iPAT //
static long           lLinesToMeasureMax         = 0;     // * Total number of measured lines without   //
                                                          // * iPAT                                     //
static long           lKSpaceCenterPartition     = 0;     // * Number of the echo (kz = 0)              //
static long           lPartitionsToMeasure       = 0;     // * Total number of measured partitions with //
                                                          // * iPAT                                     //
static long           lPartitionsToMeasureMax    = 0;     // * Total number of measured partitions      //
                                                          // * withou iPAT                              //


static long           lNumberOfQuickPRSats       = 0;
static long           lNumberOfRegularPRSats     = 0;


static long           lNoOfKernelsBetweenQuickRSats     = 0;
static long           lNoOfKernelsBetweenQuickFatSats   = 0;
static long           lNoOfKernelsBetweenQuickWaterSats = 0;

static long           lNoOfQuickSpoilGradsPerConcat     = 0;


static bool           bIgnoreQuickRSat           = false;               // * If equal true, no quick sat *
static bool           bIgnoreQuickFatSat         = false;               // * scheme will be used although*
static bool           bIgnoreQuickWaterSat       = false;               // * selected in the UI          *
                                                                        // * This is necessary to ensure *
                                                                        // * at least two kernels between*
                                                                        // * two quick sats              *
static double         dMinSatSpoilMoment         = 0.0;
static double         dMaxSatSpoilMoment         = 0.0;                 // * Spoil moment for SBB SatSpoil*
static double         dConstSatSpoilMoment       = 0.0;
static long           lSpoilMomentIncr           = 3;
static double         dSatSpoilMomentStep        = 0.0;
static long           lSatSpoilDuration          = 0;                   // * Duration of SatSpoil                       *


const long   lNoOfContrasts = 16;              // * No. of ocntrasts supported by the sequence *

const double SS_SPOILMOMENT4MAMMOGRAPHY = 4.; /* spoil moment factor in SS direction, used only for mammography */
const double SS_SPOILMOMENT_DEF         = 0.; /* default spoil moment factor in SS direction, used if no "Breast Coil" is selected */

static SEQ::PhysioSignal PhysioSignalHigh          = SEQ::SIGNAL_NONE;   // * Selected trigger modes    *
static SEQ::PhysioSignal PhysioSignalLow           = SEQ::SIGNAL_NONE;   // * and signals               *
static SEQ::PhysioMethod PhysioMethodHigh          = SEQ::METHOD_NONE;
static SEQ::PhysioMethod PhysioMethodLow           = SEQ::METHOD_NONE;

#ifdef SUPPORT_CT
//  If zero a quick FatSat is excuted immediateley before the segment which contains
//  the center line is measured. Otherwise it is execute s_i32CSatShift kernel calls earlier.
static int32_t                 s_i32CSatShift = 0;
//  Segment which contains center line
static int32_t                 s_i32ZSeg = 0xffffffff;

//  If true fSEQRunKernel_ct is used even if TimCT is OFF
static bool                    s_bUseSEQRunKernel_ct = false;

#endif

#ifndef VXWORKS
    static bool bUseTEBinarySearch = /*true*/ false;     // Disable binary search for TE
                                            // Correct limits have to be calculated by the get-limits handler

    static AffineTransformReceiver transform_receiver("10.1.13.170", 15001);

    LINK_LONG_TYPE::PFctSetValue            pOriginalContrastSetFct                 = NULL;
    LINK_LONG_TYPE::PFctSetValue            pOriginalSegmentsSetFct                 = NULL;
    LINK_LONG_TYPE::PFctSetValue            pOriginalFirstAcqWindowSetFct           = NULL;
    LINK_LONG_TYPE::PFctSetValue            pOriginalFirstAcqWindowInSetFct         = NULL;
    LINK_DOUBLE_TYPE::PFctSetValue          pOriginalTRSetFct                       = NULL;
    LINK_SELECTION_TYPE::PFctSetValue       pOriginalFirstPhysioSignalModeSetFct    = NULL;

    LINK_LONG_TYPE::PFctSolve               pOrgSolve_baseRes                       = NULL;
    LINK_LONG_TYPE::PFctSolve               pOriginalSGSizeSolveHandler             = NULL;
    LINK_SELECTION_TYPE::PFctSolve          pOriginalDimensionSolveHandler          = NULL;

    LINK_LONG_TYPE::PFctGetLimits           pOriginalSegmentsGetLimitsFct           = NULL;
    LINK_LONG_TYPE::PFctGetLimits           pOriginalFirstAcqWindowGetLimitsFct     = NULL;
    LINK_LONG_TYPE::PFctGetLimits           pOriginalFirstAcqWindowInGetLimitsFct   = NULL;
    LINK_LONG_TYPE::PFctGetLimits           pOriginalImaPerSlabGetLimitsHandler     = NULL;
    LINK_DOUBLE_TYPE::PFctGetLimits         pOriginalTRGetLimtsFct                  = NULL;
    LINK_DOUBLE_TYPE::PFctGetLimits         pOriginalFlipAngleGetLimitsHandler      = NULL;
    LINK_DOUBLE_TYPE::PFctGetLimits         pOriginalTEGetLimtsFct                  = NULL;
    LINK_DOUBLE_TYPE::PFctGetLimits         pOriginalSliceOSGetLimtsHandler         = NULL;

    LINK_LONG_TYPE::PFctIsAvailable         pOriginalPhasesIsAvailableFct           = NULL;


    LINK_LONG_TYPE::PFctGetToolTipId        pOriginalFirstAcqWindowToolTipIdHandler = NULL;

#ifdef SUPPORT_CT
    LINK_SELECTION_TYPE::PFctSetValue       pOriginalDixonSetFct                    = 0;
    LINK_SELECTION_TYPE::PFctSetValue       pOriginalMDSModeSetFct                  = 0;
#endif
#endif


// ****************************************************************************
// * Real-time event structures                                               *
// ****************************************************************************

// * ------------------------------------------------------------------------ *
// * Slice position information (rotation matrices and shifts)                *
// * ------------------------------------------------------------------------ *
static sSLICE_POS     asSLC[K_NO_SLI_MAX];



// * ------------------------------------------------------------------------ *
// * RF Pulses                                                                *
// * ------------------------------------------------------------------------ *
static sRF_PULSE*     pSRF    = NULL;
static sRF_PULSE_SINC sSRF01  (RTEIDENT_SRFExcit); // * Used for 2D sequences *
static sRF_PULSE_RECT sSRF02  (RTEIDENT_SRFExcit); // * Used for 3D non-      *
                                                   // * selective excitation  *


// * ------------------------------------------------------------------------ *
// * Sync Bits                                                                *
// * ------------------------------------------------------------------------ *
static sSYNC          sHALT01("sHALT01");



// * ------------------------------------------------------------------------ *
// * Declare instances used class                                             *
// * ------------------------------------------------------------------------ *
static const long lNoOfRSats = 8;


static SBBList                 SBB;              // * List of SBBs that will  *
                                                 // * automatically prepared  *

static SeqBuildBlockRSat       SBBRSat[lNoOfRSats];                              // * RSats                   *
static SeqBuildBlockCSat       CSatFat   ( &SBB );                               // * CSat fat suppression    *
static SeqBuildBlockCSat       CSatWater ( &SBB );                               // * CSat water suppression  *
static SeqBuildBlockMSat       MSat      ( &SBB );                               // * MSat                    *
static SeqBuildBlockTSat       SBBTSat   ( &SBB );                               // * TSat                    *
static SeqBuildBlockSpoilGrad  SatSpoil  ( NULL );                               // * Spoiler gradient after  *
static SeqBuildBlockLineTag    SBBLineTag( &SBB );                               // * Line tagging            *
static SeqBuildBlockGridTag    SBBGridTag( &SBB );                               // * Grid tagging            *
                                                                                          // * all sats                *
static SBBGREKernel            GREKernel ( NULL );                               // * GRE kernel              *

static LocalSeqLoop            RunLoop;                                          // * Standard run loop       *
static ReorderInfoGRE          Reorder ( 131072 );                               // * Reordering              *
static KernelCalculationLimits CalcLimits;

static RealTimeControl       RTControl;

static long                    lMaxTimeBetweenQuickRSats_us   =   150000;        // * Maximum allowed time    *
static long                    lMaxTimeBetweenQuickWaterSats_us = 160000;

static const double            dWEFlipAngleLimit = 30;



// * ------------------------------------------------------------------------ *
// * Local helpers                                                            *
// * ------------------------------------------------------------------------ *
template< class TYPE >
static inline void UNUSED_ARG(TYPE arg)
{
    //lint -e{550}
    if( false ) { TYPE dummy; dummy = arg; }
}









NLS_STATUS fSEQReceive
(
  SeqLim  * pSeqLim,       // IMP: sequence limits
  SeqExpo * pSeqExpo,      // EXP: sequence exports
  SEQData & rSEQData       // IMP: data-obj to process
)
{


    static const char *  ptModule = "fSEQReceive:";
    //TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s called", ptModule );

#if defined(SUPPORT_PACE) | defined(SUPPORT_CT)
    if( !RunLoop.receive(pSeqLim,pSeqExpo,rSEQData) )
    {
        TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: SeqLoopBH::receive failed.", ptModule);
        return RunLoop.getNLSStatus();
    }
#else
    UNUSED_ARG(pSeqLim);
    UNUSED_ARG(pSeqExpo);
#endif

    if (!RTControl.receiveMouseData(rSEQData))
    {
        TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: RTControl.receiveMouseData failed.", ptModule);
        return SEQU_ERROR;
    }

    return (SEQU__NORMAL);

}






/*[ Function ****************************************************************\
*
* Name        : fCalculateReordering
*
* Description : Calculates the reordering scheme for steady state sequences
*
*               The function exports
*                        - lines to measure
*                        - k-space center line
*                        - partitions to measure     and
*                        - k-space center partition
*
* Return      : NLS status
*
\****************************************************************************/
NLS_STATUS fCalculateReordering (
       MrProt    *pMrProt,                      // * IMP: Measurement protocol                   *
       SeqLim    *pSeqLim,                      // * IMP: Sequence limits                        *
       ReorderInfoGRE* pReorder,                // * IMP: ReorderInfo
       long      *plLinesToMeasureMax,          // * EXP: Total number of measured lines without *
                                                // *      iPAT                                   *
       long      *plLinesToMeasure,             // * EXP: Total number of measured lines with    *
                                                // *      iPAT                                   *
       long      *plKSpaceCenterLine,           // * EXP: k-space center line                    *
       long      *plPartitionsToMeasureMax,     // * EXP: Total number of measured partitions    *
                                                // *      without iPAT                           *
       long      *plPartitionsToMeasure,        // * EXP: Total number of measured partitions    *
                                                // *      with iPAT                              *
       long      *plKSpaceCenterPartition       // * EXP: k-space center partition               *
)
{

    static const char *ptModule = { "fCalculateReordering" };
    long lLoopOrder = PAR_IN_LIN;

    // * -------------------------------------------------------------------------- *
    // * we need to know the selected coil element for the reordering later on      *
    // * -------------------------------------------------------------------------- *
    const MrCoilSelect &coilSelect = pMrProt->coilInfo().Meas();
    SelectedCoilElements mySelectedElements(coilSelect, pSeqLim->getCoilContext());



    // * -------------------------------------------------------------------------- *
    // * Check and update number of reference lines for segmented PAT measurements  *
    // * -------------------------------------------------------------------------- *
    #ifndef VXWORKS
        #ifdef SUPPORT_iPAT
            pReorder->setOmitLowerPartOfKSpaceIfPF ( true, true );
            if (pSeqLim->isContextPrepForMrProtUpdate())  {
                if (!fPATCheckAndUpdateReferenceLineNumber(pMrProt, pSeqLim, pReorder,lPATDefOptNoRefLines))  {
                    return ( SEQU_ERROR );
                }
            }
        #endif
    #endif




    // * -------------------------------------------------------------------------- *
    // * Check UI selcetion and calculate the reordering schemes                    *
    // *                                                                            *
    // * CAUTION: The reordering schemes HAVE to be calculated prior to the         *
    // *          calculation of the rf-energy and the measurement time             *
    // * -------------------------------------------------------------------------- *
    switch ( pMrProt->kSpace().dimension() )  {
        case SEQ::DIM_2:
            pReorder->setReorderMode ( LIN_IN_PAR, SEG_ASCENDING_IN_SEG_ASCENDING, LINEAR_DESCENDING );

            pReorder->setOmitLowerPartOfKSpaceIfPF ( true, true );
            pReorder->seteReorderMethod (ReorderInfoGRE::m_GRE);
            if (! pReorder->reorderGRE ( pMrProt, pSeqLim )) return ( pReorder->getNLSStatus() );

            pReorder->deletePartitionFTFlags();

        break;

        case SEQ::DIM_3:
            // * No segmentation for 3D sequences *
            if ( pMrProt->fastImaging().segments() > 1 )  {

                if ( ! pSeqLim->isContextPrepForBinarySearch() )  {
                    TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: No segmentation in cas of 3D imaging.", ptModule );
                }

                return ( SEQU_ERROR );
            }

            lLoopOrder = PAR_IN_LIN;


            // * ------------------------------------------------------------------ *
            // * if we detect the breast coil, we have to change the reordering to  *
            // * move breathing induced artifacts from inplane to through-plane     *
            // * ------------------------------------------------------------------ *
            if (mySelectedElements.isOneApplProp("BREAST"))  {

                VectorPat<double> NormalVector = pMrProt->sliceSeries().front().normal();

                double dNormalSag = fabs(NormalVector.sag());
                double dNormalCor = fabs(NormalVector.cor());
                double dNormalTra = fabs(NormalVector.tra());


                // * ---------------------------------------------------------- *
                // * Sagittal slice orientation                                 *
                // * ---------------------------------------------------------- *
                if ( (dNormalSag > dNormalCor) && (dNormalSag > dNormalTra) )  {

                    if (   (fabs(pMrProt->sliceSeries().front().rotationAngle()) <  45.0 )
                        || (fabs(pMrProt->sliceSeries().front().rotationAngle()) > 135.0 )  ) {

                        lLoopOrder = LIN_IN_PAR;

                    } else {

                        lLoopOrder = PAR_IN_LIN;

                    }

                } else {

                    // * ------------------------------------------------------ *
                    // * Coronal slice orientation                              *
                    // * ------------------------------------------------------ *
                    if ( (dNormalCor > dNormalSag) && (dNormalSag > dNormalTra) )  {

                        lLoopOrder = PAR_IN_LIN;

                    } else {

                        // * -------------------------------------------------- *
                        // * Transversal slice orientation                      *
                        // * -------------------------------------------------- *
                        if ( (dNormalTra > dNormalSag) && (dNormalTra > dNormalCor) )  {

                            lLoopOrder = LIN_IN_PAR;

                        } else {

                            // * ---------------------------------------------- *
                            // * Two or event three elements of the normal      *
                            // * vector are identical                           *
                            // *    ===> use standard loop structure            *
                            // * ---------------------------------------------- *
                            lLoopOrder = PAR_IN_LIN;

                        }
                    }
                }
            }


            pReorder->setReorderMode ( lLoopOrder, LINEAR_ASCENDING, LINEAR_ASCENDING );


            // ellipt.scanning and 3D-PAT require the reorder3dE mode
            if ( pMrProt->kSpace().ellipticalScanning() || (   (pMrProt->PAT().PATMode() != SEQ::PAT_MODE_NONE)
                                                            && (pMrProt->PAT().AccelFact3D() > 1)               )) {
                pReorder->setOmitLowerPartOfKSpaceIfPF ( true, true );
                pReorder->seteReorderMethod (ReorderInfo::m_3DE);
                if (! pReorder->reorder3dE ( pMrProt, pSeqLim )) return ( pReorder->getNLSStatus() );
            } else {
                pReorder->setOmitLowerPartOfKSpaceIfPF ( true, true );
                pReorder->seteReorderMethod (ReorderInfoGRE::m_GRE);
                if (! pReorder->reorderGRE ( pMrProt, pSeqLim )) return( pReorder->getNLSStatus() );
            }

            switch (lLoopOrder)  {
                case PAR_IN_LIN:
                    pReorder->deletePhaseFTFlags();
                break;

                case LIN_IN_PAR:
                    pReorder->deletePartitionFTFlags();
                break;

                default:
                    TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown loop order = %ld  ", ptModule, lLoopOrder);
                    return ( SEQU_ERROR );
            }

        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
            return ( SEQU_ERROR );
    }



    // * -------------------------------------------------------------------------- *
    // * Calculate the number of measured lines and partitions                      *
    // * -------------------------------------------------------------------------- *

    switch ( pMrProt->kSpace().dimension() ) {
        case SEQ::DIM_2:

            *plKSpaceCenterLine      = pReorder->getKSCenterLin();
            *plKSpaceCenterPartition = 0;

            if (pMrProt->PAT().PATMode() == SEQ::PAT_MODE_NONE)  {

                *plLinesToMeasure       = *plLinesToMeasureMax      = pReorder->getLinesToMeasure();
                *plPartitionsToMeasure  = *plPartitionsToMeasureMax = pMrProt->kSpace().partitions();

            } else {

                *plLinesToMeasure         = pReorder->getPATLinesToMeasure();
                *plLinesToMeasureMax      = pReorder->getLinesToMeasure();

                *plPartitionsToMeasure    = pReorder->getPATPartitionsToMeasure();
                *plPartitionsToMeasureMax = pMrProt->kSpace().partitions();

            }

        break;

        case SEQ::DIM_3:

            *plKSpaceCenterLine      = pReorder->getKSCenterLin();
            *plKSpaceCenterPartition = pReorder->getKSCenterPar();


            if (pMrProt->PAT().PATMode() == SEQ::PAT_MODE_NONE)  {

                *plLinesToMeasure           = *plLinesToMeasureMax      = pReorder->getLinesToMeasure();
                *plPartitionsToMeasure      = *plPartitionsToMeasureMax = pReorder->getPartitionsToMeasure();

            } else {

                *plLinesToMeasure           = pReorder->getPATLinesToMeasure();
                *plLinesToMeasureMax        = pReorder->getLinesToMeasure();

                *plPartitionsToMeasure      = pReorder->getPATPartitionsToMeasure();
                *plPartitionsToMeasureMax   = pReorder->getPartitionsToMeasure();

            }

        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
            return ( SEQU_ERROR );
        }

    return ( SEQU__NORMAL );

}

#ifdef SUPPORT_CT
//  TimePerShift_us ... is the time between succesive excitation pulses
//  TimeConst_us    ... is the time between the center of the FS-pulse and the center of the next
//                      excitation pulse
//  dTR_us          ... time between successive fat suppression pulses
static bool CT_QSAT_FLIP_ANGLE(double& rdCSATFlipAngle_deg, double dTR_us, double dT1_us, double dTnull_us)
{
    const double dE1 = exp(-dTR_us/dT1_us)
        ,        dE2 = exp(-dTnull_us/dT1_us)
        ;
    //  FlipAngle = acos(E2-1/(E2-E1))
    //  E1 = exp(-TR/T1)
    //  E2 = exp(-Tnull/T1)
    //  argument must be in the range -1, ..., 0 to be mapped to PI/2
        //  i)  E2 < 1 --> E2 > E1
        //  ii) (1-E2) < (E2-E1) -> 1+E1 < 2E2
    if( (dE2 > dE1) && ((1+dE1) < 2*dE2) )
    {
        rdCSATFlipAngle_deg = (180/M_PI)*acos((dE2-1)/(dE2-dE1));
        return true;
    }
    rdCSATFlipAngle_deg = 180.;
    return false;
}
static bool CT_QSAT_SHIFT_AND_FLIP_ANGLE(double& rdCSATFlipAngle_deg, int32_t& ri32Shift, int32_t i32ZPos, int32_t i32NSat, double dTR_us, double dT1_us, double dTimeConst_us, double dTimePerShift_us)
{
    if( i32NSat < 1 )
    {
        //  Nothing to do
        return true;
    }
    dTR_us /= i32NSat;
    //  Maximum Shift
    ri32Shift = i32ZPos;

    const double dE1 = exp(-dTR_us/dT1_us);
    double dTI_us, dE2;
    do
    {
        //  FlipAngle = acos(E2-1/(E2-E1))
        //  E1 = exp(-TR/T1)
        //  E2 = exp(-TI/T1)
        dTI_us = ri32Shift*dTimePerShift_us+dTimeConst_us;
        dE2    = exp(-dTI_us/dT1_us);

        //  argument must be in the range -1, ..., 0 to be mapped to PI/2
        //  i)  E2 < 1 --> E2 > E1
        //  ii) (1-E2) < (E2-E1) -> 1+E1 < 2E2
        if( (dE2 > dE1) && ((1+dE1) < 2*dE2) )
        {
            break;
        }
    }
    while(--ri32Shift >= 0);
    if( ri32Shift < 0 )
    {
        //  No solution found: i32TimeConst_us probably too large
        //  Set shift to zero and use flip angle 180
        rdCSATFlipAngle_deg = 180;
        ri32Shift           = 0;
    }
    else
    {
        rdCSATFlipAngle_deg = (180/M_PI)*acos((dE2-1)/(dE2-dE1));
    }
    return true;

}
#endif






/*[ Function ****************************************************************\
*
* Name        : fSelectAndPrepRFPulse
*
* Description : Declaration of the rf-pulse structure *pSRF depending on the
*               protocol settings and preparation
*
* Return      : NLS status
*
\****************************************************************************/
NLS_STATUS  fSelectAndPrepRFPulse (
       MrProt           *pMrProt,       // * IMP: Measurement protocol        *
       SeqLim           *pSeqLim,       // * IMP: Sequence Limits             *
       SeqExpo          *pSeqExpo,      // * IMP: Returned values             *
       sRF_PULSE       **pSRF           // * EXP: Pointer to excitation pulse *
)
{

  double dEmpiricalFactor = 1.0;             // * Empirical factor for sinc   *
                                             // * pulse calculation.          *
                                             // * The pulse profile and slice *
                                             // * thickness is modified for   *
                                             // * values <> 1.0               *



  static const char *ptModule = { "fSelectAndPrepRFPulse" };


  // * ---------------------------------------------------------------------- *
  // * Prepare the RF pulse structures                                        *
  // * ---------------------------------------------------------------------- *

  switch ( pMrProt->kSpace().dimension() ) {
    case SEQ::DIM_2:


      sSRF01.setTypeExcitation         ();
      sSRF01.setFlipAngle              ( pMrProt->flipAngle() );
      sSRF01.setInitialPhase           ( 0.0 );
      sSRF01.setThickness              ( pMrProt->sliceSeries().aFront().thickness() );

      switch ( pMrProt->txSpec().rfPulseType() ) {

            case  SEQ::RF_FAST:
                if (pMrProt->bandWidthPerPixel (pSeqLim->getReadoutOSFactor())[0] > 1560)    {
                    sSRF01.setSamples              (  100 );
                    sSRF01.setDuration             (  400 );
                }  else  {
                    sSRF01.setSamples              (  500 );
                    sSRF01.setDuration             ( 1000 );
                }
                sSRF01.setBandwidthTimeProduct ( 2.00 );
            break;

            case  SEQ::RF_NORMAL:
                sSRF01.setSamples              (  500 );
                sSRF01.setDuration             ( 2000 );
                sSRF01.setBandwidthTimeProduct ( 2.70 );
            break;

            case  SEQ::RF_LOW_SAR:
                sSRF01.setSamples              (  500 );
                sSRF01.setDuration             ( 4000 );
                sSRF01.setBandwidthTimeProduct ( 2.00 );
            break;

            default: TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown rf-pulse type = %d  ", ptModule, pMrProt->txSpec().rfPulseType());
                     return ( SEQU_ERROR );
      }

      // * Preparation of the rf-pulse *
      if (! sSRF01.prepSinc (pMrProt, pSeqExpo) ) return ( sSRF01.getNLSStatus() );

      *pSRF = &sSRF01;

    break;

    case SEQ::DIM_3:
      if ( pMrProt->txSpec().excitation() == SEQ::EXCITATION_SLICE_SELECTIVE )  {

        dEmpiricalFactor = 0.8;
        sSRF01.setTypeExcitation         ();
        sSRF01.setSamples                (  500 );
        sSRF01.setFlipAngle              ( pMrProt->flipAngle() );
        sSRF01.setInitialPhase           ( 0.0 );
        sSRF01.setThickness              ( pMrProt->sliceSeries().aFront().thickness() );
        sSRF01.setThickness              ( pMrProt->sliceSeries().aFront().thickness() * dEmpiricalFactor );

        switch ( pMrProt->txSpec().rfPulseType() ) {
          case  SEQ::RF_FAST:    sSRF01.setDuration             ( 1000 );
                                 sSRF01.setBandwidthTimeProduct (  6.4 * dEmpiricalFactor );
          break;

          case  SEQ::RF_NORMAL:  sSRF01.setDuration             ( 2000 );
                                 sSRF01.setBandwidthTimeProduct ( 12.7 * dEmpiricalFactor );
          break;

          case  SEQ::RF_LOW_SAR: sSRF01.setDuration             ( 4000 );
                                 sSRF01.setBandwidthTimeProduct ( 25.4 * dEmpiricalFactor );

          break;

          default: TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown rf-pulse type = %d  ", ptModule, pMrProt->txSpec().rfPulseType());
          return ( SEQU_ERROR );
        }

        #ifndef VXWORKS
           SeqUT.setRFThicknessInfo ( &sSRF01, dEmpiricalFactor * pMrProt->sliceSeries().aFront().thickness() );
        #endif

        // * Preparation of the rf-pulse *
        if (! sSRF01.prepSinc (pMrProt, pSeqExpo) ) return ( sSRF01.getNLSStatus() );

        *pSRF = &sSRF01;

      } else {

        sSRF02.setTypeExcitation         ();
        sSRF02.setDuration               (  100 );
        sSRF02.setSamples                (   50 );
        sSRF02.setInitialPhase           ( 0.00 );
        sSRF02.setFlipAngle              ( pMrProt->flipAngle() );
        sSRF02.setThickness              ( 10000.0 ); // * This is just a dummy slice thickness *

        if (! sSRF02.prepRect ( pMrProt, pSeqExpo) ) return ( sSRF02.getNLSStatus() );

        #ifndef VXWORKS
           SeqUT.setRFThicknessInfo ( &sSRF02, 10000.0 );
        #endif

        *pSRF = &sSRF02;

      }

    break;

    default:
      TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
      return ( SEQU_ERROR );
  }






  return ( SEQU__NORMAL );

}




/*[ Function ****************************************************************\
*
* Name        : fSetKernelParameter
*
* Description : Sets the parameter that are required to calculate the kernel
*               timing
*
* Return      : NLS status
*
\****************************************************************************/
NLS_STATUS fSetKernelParameter (MrProt  *pMrProt,
                                SeqLim  *pSeqLim,
                                long     lLinesToMeasureMax,          // * IMP: Total number of measured lines         *
                                long     lKSpaceCenterLine,           // * IMP: k-space center line                    *
                                long     lPartitionsToMeasureMax,     // * IMP: Total number of measurd partitions     *
                                long     lKSpaceCenterPartition       // * IMP: k-space center partition               *
)
{
    static const char *ptModule = { "fSetKernelParameter" };


    static GPAProxy  theGPA;
    static MSUProxy  theMSU;

    double adROMaxGradAmpl[3], adPEMaxGradAmpl[3], adSSMaxGradAmpl[3];


    // * ---------------------------------------------------------------------- *
    // * we need to know the selected coil element for an additional spoiler    *
    // * ---------------------------------------------------------------------- *
    const MrCoilSelect &coilSelect = pMrProt->coilInfo().Meas();
    SelectedCoilElements mySelectedElements(coilSelect, pSeqLim->getCoilContext());


    // * ---------------------------------------------------------------------- *
    // * Set the gradient performance parameters                                *
    // * ---------------------------------------------------------------------- *
    GREKernel.setDefaultGradientPerformance();



    // * ---------------------------------------------------------------------- *
    // * Reduce gradient amplitudes of the phase encode gradient tables and     *
    // * water excitation pulses on low field systems and water excitation      *
    // * to avoid hysteresis effects.                                           *
    // * ---------------------------------------------------------------------- *
    if ( ( theMSU.getNominalB0() < LOW_FIELD_B0 ) && ( pMrProt->preparationPulses().fatSuppression() == SEQ::WATER_EXCITATION ) )  {
        double adMaxGradAmpl[3] = { 4.0, 4.0, 4.0 };

        GREKernel.setMaxMagnitudes ( adMaxGradAmpl, SBBGREKernel_GRAD_GROUP_PE );
        GREKernel.setMaxMagnitudes ( adMaxGradAmpl, SBBGREKernel_GRAD_GROUP_3D );
        GREKernel.setMaxMagnitudes ( adMaxGradAmpl, SBBGREKernel_GRAD_GROUP_WE );
    }




    // * ---------------------------------------------------------------------- *
    // * Gradient amplitude adjustments for realtime gre  .                     *
    // * ---------------------------------------------------------------------- *
    if ( eSequence == RT_Gre )  {
        if ( fabs (theGPA.getGradMaxAmplAbsolute() - 20.0) < 0.1 )  {  // * Turbo and ultra gradients *

            adROMaxGradAmpl[0] = 0.875 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 14 mT/m *
            adROMaxGradAmpl[1] = 0.875 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 14 mT/m *
            adROMaxGradAmpl[2] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 10 mT/m *

            adPEMaxGradAmpl[0] = 0.625 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 10 mT/m *
            adPEMaxGradAmpl[1] = 0.625 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 10 mT/m *
            adPEMaxGradAmpl[2] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 10 mT/m *

            adSSMaxGradAmpl[0] = 0.625 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 10 mT/m *
            adSSMaxGradAmpl[1] = 0.625 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 10 mT/m *
            adSSMaxGradAmpl[2] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 10 mT/m *

        } else {

            if ( fabs (theGPA.getGradMaxAmplAbsolute() - 30.0) < 0.1 )  {  // * Qunatum gradients *

                adROMaxGradAmpl[0] = 0.863 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 19 mT/m *
                adROMaxGradAmpl[1] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 16 mT/m *
                adROMaxGradAmpl[2] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 10 mT/m *

                adPEMaxGradAmpl[0] = 0.727 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 16 mT/m *
                adPEMaxGradAmpl[1] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 16 mT/m *
                adPEMaxGradAmpl[2] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 10 mT/m *

                adSSMaxGradAmpl[0] = 0.727 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 16 mT/m *
                adSSMaxGradAmpl[1] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 16 mT/m *
                adSSMaxGradAmpl[2] = 1.000 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 10 mT/m *

            } else {

                if ( fabs (theGPA.getGradMaxAmplAbsolute() - 40.0) < 0.1 )  {  // * Sonata *

                    adROMaxGradAmpl[0] = 0.892 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 25 mT/m *
                    adROMaxGradAmpl[1] = 0.786 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 22 mT/m *
                    adROMaxGradAmpl[2] = 0.571 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 16 mT/m *

                    adPEMaxGradAmpl[0] = 0.714 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 20 mT/m *
                    adPEMaxGradAmpl[1] = 0.786 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 22 mT/m *
                    adPEMaxGradAmpl[2] = 0.571 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 16 mT/m *

                    adSSMaxGradAmpl[0] = 0.571 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );  // * 16 mT/m *
                    adSSMaxGradAmpl[1] = 0.571 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );  // * 16 mT/m *
                    adSSMaxGradAmpl[2] = 0.571 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );  // * 16 mT/m *

                } else {  // * unkown gradient system; 1/sqrt(3)*gradMaxAmpl is used on all axes *

                    adROMaxGradAmpl[0] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );
                    adROMaxGradAmpl[1] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );
                    adROMaxGradAmpl[2] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );

                    adPEMaxGradAmpl[0] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );
                    adPEMaxGradAmpl[1] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );
                    adPEMaxGradAmpl[2] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );

                    adSSMaxGradAmpl[0] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_FAST    );
                    adSSMaxGradAmpl[1] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_NORMAL  );
                    adSSMaxGradAmpl[2] = 0.577 * theGPA.getGradMaxAmpl (SEQ::GRAD_WHISPER );

                }

            }

        }

        GREKernel.setMaxMagnitudes ( adSSMaxGradAmpl, SBBGREKernel_GRAD_GROUP_SS  );
        GREKernel.setMaxMagnitudes ( adSSMaxGradAmpl, SBBGREKernel_GRAD_GROUP_WE  );
        GREKernel.setMaxMagnitudes ( adPEMaxGradAmpl, SBBGREKernel_GRAD_GROUP_PE  );
        GREKernel.setMaxMagnitudes ( adSSMaxGradAmpl, SBBGREKernel_GRAD_GROUP_3D  );
        GREKernel.setMaxMagnitudes ( adROMaxGradAmpl, SBBGREKernel_GRAD_GROUP_ROP );
        GREKernel.setMaxMagnitudes ( adROMaxGradAmpl, SBBGREKernel_GRAD_GROUP_RO  );
        GREKernel.setMaxMagnitudes ( adROMaxGradAmpl, SBBGREKernel_GRAD_GROUP_ROD );
    }



    // * ---------------------------------------------------------------------- *
    // * Gradient settings for Gre_HR                                           *
    // * ---------------------------------------------------------------------- *
    if ( eSequence == Gre_HR )  {
        adROMaxGradAmpl[0] = pMrProt->fastImaging().echoSpacing();
        adROMaxGradAmpl[1] = pMrProt->fastImaging().echoSpacing();
        adROMaxGradAmpl[2] = pMrProt->fastImaging().echoSpacing();

        GREKernel.setMaxMagnitudes ( adROMaxGradAmpl, SBBGREKernel_GRAD_GROUP_ROP );
        GREKernel.setMaxMagnitudes ( adROMaxGradAmpl, SBBGREKernel_GRAD_GROUP_RO  );
        GREKernel.setMaxMagnitudes ( adROMaxGradAmpl, SBBGREKernel_GRAD_GROUP_ROD );
    }



    // * ---------------------------------------------------------------------- *
    // * Set gradient performance for Kernel to support GSWD binary search for  *
    // * optimum risetime.                                                      *
    // *                                                                        *
    // * NOTE: This function should be called AFTER all gradient performance    *
    // *       data was set for fast normal whisper and BEFORE individual       *
    // *       settings for the GSWD mode are made.                             *
    // * ---------------------------------------------------------------------- *
    GREKernel.setGSWDGradientPerformance (pMrProt, pSeqLim);



    // * ---------------------------------------------------------------------- *
    // * Calculate timing of the GRE Kernel                                     *
    // * ---------------------------------------------------------------------- *
    if ( eSequence == Gre_HR  )  {
        CalcLimits.setLimitsForPixelSizeRO    ( 0.0 );
        CalcLimits.setLimitsForPixelSizePE    ( 0.0 );
        CalcLimits.setLimitsForPixelSize3D    ( 0.0 );
        CalcLimits.setLimitsForSliceThickness ( 0.0 );
    } else {
        CalcLimits.setDefaultLimits();
    }





    if (! GREKernel.setPointerToCalculationLimits ( &CalcLimits ))  return (GREKernel.getNLSStatus());

    GREKernel.setbUsePERewinder                 ( true     );
    GREKernel.setbSendOscBit                    ( true     );
    if ( eSequence == Gre_HR )  {
        GREKernel.seteReadOutGradType ( SBBGREKernel::Constant );
    } else {
        GREKernel.seteReadOutGradType ( SBBGREKernel::Spoiler  );
    }



    switch ( pMrProt->flowComp()[0] )  {
        case SEQ::FLOW_COMPENSATION_YES:
                    GREKernel.setbFlowCompPhase           ( true );
                    GREKernel.setbFlowCompRead            ( true );
                    GREKernel.setbFlowCompSlice           ( true );
        break;

        case SEQ::FLOW_COMPENSATION_SLICE_READ:
                    GREKernel.setbFlowCompPhase           ( false );
                    GREKernel.setbFlowCompRead            ( true  );
                    GREKernel.setbFlowCompSlice           ( true  );
        break;

        case SEQ::FLOW_COMPENSATION_SLICESEL_ONLY:
                    GREKernel.setbFlowCompPhase           ( false );
                    GREKernel.setbFlowCompRead            ( false );
                    GREKernel.setbFlowCompSlice           ( true  );
        break;

        case SEQ::FLOW_COMPENSATION_READOUT_ONLY:
                    GREKernel.setbFlowCompPhase           ( false );
                    GREKernel.setbFlowCompRead            ( true  );
                    GREKernel.setbFlowCompSlice           ( false );
        break;

        default:
                    GREKernel.setbFlowCompPhase           ( false );
                    GREKernel.setbFlowCompRead            ( false );
                    GREKernel.setbFlowCompSlice           ( false );

    }



    GREKernel.setbRampsOutsideOfEventBlock      ( false    );

    if ( pMrProt->readOutMode() == SEQ::READOUT_BIPOLAR )  {
        GREKernel.seteReadOutPolarity (SBBGREKernel::Bipolar  );
    } else {
        GREKernel.seteReadOutPolarity (SBBGREKernel::Monopolar);
    }

    if ( (eSequence == Gre) || (eSequence == Gre_HR)) {
        GREKernel.setlStartTimeWakeUp           ( 300      );
    }


    if ( theMSU.getNominalB0() < LOW_FIELD_B0 )  {

        switch ( pMrProt->kSpace().dimension() )  {       // * WE pulse used for OPEN system, 0.2 T *
            case SEQ::DIM_2:    GREKernel.seteWEMode ( SBBGREKernel::WE_1_2_1_90_Deg );     break;
            case SEQ::DIM_3:    GREKernel.seteWEMode ( SBBGREKernel::WE_1_1_75_Deg   );     break;
            default:            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
                                return ( SEQU_ERROR );
        }

        GREKernel.setdWEBandwidthTimeProduct    ( 2.9 );
        GREKernel.useWEOffsetFrequency          ( true, pMrProt->preparationPulses().binomialPulsesOffsetFreq() );

    } else {

        switch ( pMrProt->preparationPulses().fatSuppression() )  {
            case SEQ::WATER_EXCITATION:         GREKernel.seteWEMode ( SBBGREKernel::WE_1_2_1_180_Deg );    break;
            case SEQ::WATER_EXCITATION_FAST:    GREKernel.seteWEMode ( SBBGREKernel::WE_1_1_180_Deg );      break;
            default: ;
        }

        GREKernel.setdWEBandwidthTimeProduct    ( 2.0 );
        GREKernel.useWEOffsetFrequency          ( false );
    }

    // for 3D, use higher BTProduct
    if (pMrProt->kSpace().dimension() == SEQ::DIM_3) {

        // * ---------------------------------------------------------------------- *
        // * Reduce BTP in order to reduce B1 amplitude to less than 17uT           *
        // * ---------------------------------------------------------------------- *
        if ( SysProperties::isVerio() )  {

            if ( pMrProt->flipAngle() <= dWEFlipAngleLimit )  {
                GREKernel.setdWEBandwidthTimeProduct ( 5.0 );
            } else {
                GREKernel.setdWEBandwidthTimeProduct ( 3.0 );
            }

        } else {
            GREKernel.setdWEBandwidthTimeProduct ( 10.0 );
        }

    }



    GREKernel.setlFirstLineToMeasure    ( - lKSpaceCenterLine                        );
    GREKernel.setlLastLineToMeasure     ( lLinesToMeasureMax - lKSpaceCenterLine - 1 );
    GREKernel.setlKSpaceCenterLine      ( lKSpaceCenterLine );

    GREKernel.setlFirstPartToMeasure    ( -lKSpaceCenterPartition                              );
    GREKernel.setlLastPartToMeasure     ( lPartitionsToMeasureMax - lKSpaceCenterPartition - 1 );
    GREKernel.setlKSpaceCenterPartition ( lKSpaceCenterPartition );




    // * -------------------------------------------------------------------------- *
    // * Non-selective excitation requires a stronger RO spoiler since spoiling due *
    // * to the slice selection is missing.                                         *
    // * -------------------------------------------------------------------------- *
    if ( pMrProt->txSpec().excitation() == SEQ::EXCITATION_VOLUME_SELECTIVE )  {
        GREKernel.setdRelRODephasingMoment (1.0);
    } else {
        if ( fabs( theMSU.getNominalB0() - 2.9) <= 0.1 )  {  // 3T system
            GREKernel.setdRelRODephasingMoment (1.0);
        } else {
            GREKernel.setdRelRODephasingMoment (0.5);
        }
    }




    // * ---------------------------------------------------------------------- *
    // * if we detect the breast coil, we add a spoil moment in slice select    *
    // * direction to get rid of lines at edges of the object                   *
    // * ---------------------------------------------------------------------- *
    GREKernel.setdRelSSDephasingMoment(SS_SPOILMOMENT_DEF);
    if (mySelectedElements.isOneApplProp("BREAST"))  {
           GREKernel.setdRelSSDephasingMoment(SS_SPOILMOMENT4MAMMOGRAPHY);
    }



    return ( SEQU__NORMAL );

}



/*[ Function ****************************************************************\
*
* Name        : fSetAsymmetry
*
* Description : This function defines the echo asymmetry of the SBBGREKernel
*
*               The following modes are supported:
*                   Symmetric : The full number of data points prior to the
*                               echo maximum are acquired.
*                   Asymmetric: The number of data points acquired prior to
*                               the echo is reduced. The number of data points
*                               after the echo maximum is reduced in case of
*                               double contrast.
*
* Return        : NLS status
*
\****************************************************************************/
NLS_STATUS fSetAsymmetry ( MrProt* pMrProt, eEchoMode Mode )
{

    static const char *ptModule = { "fSetAsymmetry" };

    static GPAProxy  theGPA;
    static MSUProxy  theMSU;

    long lI;  // Loop counter


    switch ( Mode )  {
        case Symmetric:
            for ( lI = 0; lI < pMrProt->contrasts(); lI++ )  {
                GREKernel.setadReadOutAsym ( lI, 128.0 / 256.0, 128.0 / 256.0 );
            }
        break;

        case Asymmetric:
            // * The acquisition window of the first echo is reduced for Syntata gradient systems *
            // * in case of multiple contrasts to enable in-phase opposed-phase measurements      *
            if (    (fabs(theMSU.getNominalB0() - 1.5) <= 0.1) && ( theGPA.getGradMaxAmplAbsolute() <= 35.0 )
                 && (pMrProt->contrasts() > 1) )  {
                GREKernel.setadReadOutAsym ( 0, 86.0 / 256.0, 114.0 / 256.0 );
            } else {
                GREKernel.setadReadOutAsym ( 0, 86.0 / 256.0, 128.0 / 256.0 );
            }


            for ( lI = 1; lI < pMrProt->contrasts(); lI++ )  {
                if ( (fabs(theMSU.getNominalB0() - 1.5) <= 0.1) && ( theGPA.getGradMaxAmplAbsolute() <= 35.0 ) )  {
                    GREKernel.setadReadOutAsym ( lI, 80.0 / 256.0,  128.0 / 256.0 );
                } else {
                    GREKernel.setadReadOutAsym ( lI, 86.0 / 256.0, 128.0 / 256.0 );
                }
            }
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown echo mode 2 = %d  ", ptModule, Mode);
            return ( SEQU_ERROR );
    }



    return ( SEQU__NORMAL );

}











/*[ Function ****************************************************************\
*
* Name        : fCalculateGRETiming
*
* Description : Calculates the timing of the GRE kernel
*               Symmetric echoes are used at first. If this fails due a
*               negative TE fill time an asymmetric echo is used. If this
*               fails again, the function retrurns the NLS status
*
* Return      : NLS status
*
\****************************************************************************/
NLS_STATUS fCalculateGRETiming ( MrProt*   pMrProt,
                                 SeqLim*   pSeqLim,
                                 SeqExpo*  pSeqExpo
)
{

    static GPAProxy  theGPA;
    static MSUProxy  theMSU;





    // * ---------------------------------------------------------------------- *
    // * Asymmetric echo is used during PrepForBinarySearch in order to         *
    // * calculate the TE limits appropreately.                                 *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->kSpace().asymmetricEchoAllowed() && pSeqLim->isContextPrepForBinarySearch() )  {

        if ( SEQU__NORMAL != fSetAsymmetry (pMrProt, Asymmetric) )  {
            return ( SEQU_ERROR );
        }


        GREKernel.setbSuppressTraces (true);
        if (! GREKernel.calculateTiming (pMrProt, pSeqLim, pSeqExpo))  {
            GREKernel.setbSuppressTraces (false);
            return (GREKernel.getNLSStatus());
        }

        GREKernel.setbSuppressTraces (false);
        return ( SEQU__NORMAL );

    } else {

        // * ---------------------------------------------------------------------- *
        // * Calculate and prepare the GRE kernel                                   *
        // * using symmetric echoes at first                                        *
        // * ---------------------------------------------------------------------- *
        if ( SEQU__NORMAL != fSetAsymmetry (pMrProt, Symmetric ) )  {
            return ( SEQU_ERROR );
        }


        GREKernel.setbSuppressTraces (true);
        if (! GREKernel.calculateTiming (pMrProt, pSeqLim, pSeqExpo))  {

            // * ------------------------------------------------------------------ *
            // * If negative TE fill time try asymmetric echoes (if allowed)        *
            // * ------------------------------------------------------------------ *
            if ( (GREKernel.getNLSStatus() == SEQU__NEGATIV_TEFILL) && (pMrProt->kSpace().asymmetricEchoAllowed()) )  {


                if ( SEQU__NORMAL != fSetAsymmetry (pMrProt, Asymmetric ) )  {
                    return ( SEQU_ERROR );
                }


                if (! GREKernel.calculateTiming (pMrProt, pSeqLim, pSeqExpo))  {
                    GREKernel.setbSuppressTraces (false);
                    return (GREKernel.getNLSStatus());
                }

            } else {

                GREKernel.setbSuppressTraces (false);
                return (GREKernel.getNLSStatus());
            }
        }

        GREKernel.setbSuppressTraces (false);
        return ( SEQU__NORMAL );

    }

}








/*[ Function ****************************************************************\
*
* Name        : fCalculateNoOfQuickSatsPerConcat
*
* Description : This function calculates the number of quick sats that have
*               to be executed per concatenation.
*
* Return      : NLS status
*
\****************************************************************************/
NLS_STATUS fCalculateNoOfQuickSatsPerConcat (
                MrProt     *pMrProt,                       // * Protocol structure              *
                long        lMaxTimeBetweenQuickSats_us,   // * Maximum time that is allowed to *
                                                           // * pass between two quick sats     *
                                                           // * The time is give in us.         *
                long        lMaxNoOfSlicesPerConcat,       // * Slices per concatenatio         *
                long       *lNoOfQuickSatsPerConcat,       // * No. of quick sats that have to  *
                                                           // * executed in a concatenation     *
                long       *lNoOfKernelsBetweenQuickSats,  // * No. of kernels that will be     *
                                                           // * executed between two quick sats *
                bool       *bIgnoreQuickSat,               // * Ignore quick sat                *
                long       lSync = 1)                      // * Kernel between quick sats used  *
                                                           // * for syncronization              *
{

    static const char *ptModule = {"fCalculateNoOfQuickSatsPerConcat"};




    // * ------------------------------------------------------------------------- *
    // * Calculate the number of quick sats that have to be executed in one        *
    // * concatenation.                                                            *
    // * ------------------------------------------------------------------------- *
#ifdef SUPPORT_CT
    if( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct )
    {
         *lNoOfQuickSatsPerConcat      = pMrProt->MDS().mdsTableSpeedNumerator();
         *lNoOfKernelsBetweenQuickSats = (lMaxNoOfSlicesPerConcat + *lNoOfQuickSatsPerConcat - 1)/(*lNoOfQuickSatsPerConcat);
    }
    else
#endif
    if ( (lMaxTimeBetweenQuickSats_us > 0) &&  (lSync > 0) )  {

        // ceil: TR / lMaxTimeBetweenQuickSats_us
        *lNoOfQuickSatsPerConcat = (pMrProt->tr()[0] + lMaxTimeBetweenQuickSats_us - 1) / lMaxTimeBetweenQuickSats_us;


        if ( *lNoOfQuickSatsPerConcat == 0 )  {

            *lNoOfKernelsBetweenQuickSats = lMaxNoOfSlicesPerConcat;
        } else  {

            // ceil: lMaxNoOfSlicesPerConcat / (*lNoOfQuickSatsPerConcat)
            *lNoOfKernelsBetweenQuickSats = (lMaxNoOfSlicesPerConcat + *lNoOfQuickSatsPerConcat - 1) /
                                            (*lNoOfQuickSatsPerConcat);
        }


        // * Syncronization *
        if ( *lNoOfKernelsBetweenQuickSats % lSync )  {
            *lNoOfKernelsBetweenQuickSats = (*lNoOfKernelsBetweenQuickSats + lSync - 1) / lSync * lSync;
        }


        // ceil: lMaxNoOfSlicesPerConcat / (*lNoOfKernelsBetweenQuickSats)
        *lNoOfQuickSatsPerConcat = (lMaxNoOfSlicesPerConcat + *lNoOfKernelsBetweenQuickSats - 1) /
                                   (*lNoOfKernelsBetweenQuickSats);

    } else {

        TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: No valid value given for lMaxTimeBetweenQuickSats_us.", ptModule);
        return ( SEQU_ERROR );

    }



    // * ------------------------------------------------------------------------- *
    // * If the number of quick sats is to large, ignore quick sat, i.e. one sat   *
    // * per kernel                                                                *
    // * ------------------------------------------------------------------------- *
    if ( *lNoOfQuickSatsPerConcat >= lMaxNoOfSlicesPerConcat ) {

        *bIgnoreQuickSat = true;
        *lNoOfQuickSatsPerConcat = lMaxNoOfSlicesPerConcat;

    } else {

        *bIgnoreQuickSat = false;

    }



    return ( SEQU__NORMAL );

}







/*[ Function ****************************************************************\
*
* Name        : fCalculateTEMin
*
* Description : Calculates the minimum TE values of the GRE kernel
*               Following steps are performed:
*                 - calculate the reordering scheme in order to determine
*                   the lines/partitions to measure and the k-space center
*                   lines/partitions
*                 - define the rf-pulse structure
*                 - calculate the GRE timing
*
* Return      : NLS status
*
\****************************************************************************/
NLS_STATUS fCalculateTEMin (
       MrProt      *pMrProt,       // * IMP: Measurement protocol   *
       SeqLim      *pSeqLim,       // * IMP: Sequence limits        *
       SeqExpo     *pSeqExpo,      // * IMP: Sequence exports       *
       Sequence    *pSequence,     // * IMP: Sequence class         *
       long        *alTEMin,       // * EXP: Array of minimum TE    *
                                   // *      times                  *
       long        *alTEFill       // * EXP: Array of TE fill times *
)
{
    NLS_STATUS  lStatus = SEQU__NORMAL;



    // * ---------------------------------------------------------------------- *
    // * Calculate the reordering scheme                                        *
    // * -> lLinesToMeasureMax, lKSpaceCenterLine, ...                          *
    // * ---------------------------------------------------------------------- *
    lStatus = fCalculateReordering ( pMrProt,
                                     pSeqLim,
                                    &Reorder,
                                    &lLinesToMeasureMax,
                                    &lLinesToMeasure,
                                    &lKSpaceCenterLine,
                                    &lPartitionsToMeasureMax,
                                    &lPartitionsToMeasure,
                                    &lKSpaceCenterPartition );
    if ( lStatus != SEQU__NORMAL ) return ( lStatus );



    // * ---------------------------------------------------------------------- *
    // * Set and prepare the rf-pulse structure                                 *
    // * ---------------------------------------------------------------------- *
    {

        AutoPrepContext MyContext (pSequence, pMrProt);

        if ( SEQU__NORMAL != (lStatus = fSelectAndPrepRFPulse (pMrProt, pSeqLim, pSeqExpo, &pSRF) ) ) return ( lStatus );
        if (! GREKernel.setRFPulseForExcitation (pSRF) )  return ( GREKernel.getNLSStatus() );



        // * ------------------------------------------------------------------ *
        // * Set all other kernel parameters such as:                           *
        // *    calculation limits, echo asymmetry, ...                         *
        // * ------------------------------------------------------------------ *
        fSetKernelParameter( pMrProt,
                             pSeqLim,
                             lLinesToMeasureMax,
                             lKSpaceCenterLine,
                             lPartitionsToMeasureMax,
                             lKSpaceCenterPartition );



        // * ------------------------------------------------------------------ *
        // * Calculate and prepare the GRE kernel                               *
        // *                                                                    *
        // * using symmetric echoes at first                                    *
        // * ------------------------------------------------------------------ *
        lStatus = fCalculateGRETiming ( pMrProt, pSeqLim, pSeqExpo );

    }



    for ( long lI = 0; lI < pMrProt->contrasts(); lI++ )  {

        alTEMin[lI]  = GREKernel.getalTEMin()[lI];
        alTEFill[lI] = GREKernel.getalTEFill()[lI];

    }



    return ( lStatus );

}








#ifndef VXWORKS







/*[ Function ****************************************************************\
*
* Name        : pfGetMinTE
*
* Description : This function calculates the TE values
*               of the GRE kernel within a solve handler
*               fSEQPrep is called prior to TE calculation to make *pMrProt
*               and *pSeqExpo consistent
*
*               Following steps are performed:
*                 - fSEQPrep is called
*                 - calculate the reordering scheme in order to determine
*                   the lines/partitions to measure and the k-space center
*                   lines/partitions
*                 - define the rf-pulse structure
*                 - calculate the GRE timing
*
* Return      : BOOL
*
\****************************************************************************/
BOOL pfGetMinTE (
       MrProt      *pMrProt,       // * IMP: Measurement protocol  *
       SeqLim      *pSeqLim,       // * IMP: Sequence limits       *
       SeqExpo     *pSeqExpo,      // * IMP: Sequence exports      *
       Sequence    *pSeq,          // * IMP:                       *
       long        *alTEMin        // * EXP: Array of TE times     *
)
{
  NLS_STATUS lStatus = SEQU__NORMAL;

  long lI = 0;  // * Loop counter *
  long alTEFill[lNoOfContrasts];

  // * ---------------------------------------------------------------------- *
  // * Prior to the TE calculation *pMrProt and *pSeqExpo have to be made     *
  // * consistent by calling fSEQPrep                                         *
  // * ---------------------------------------------------------------------- *
  bool bOldPerformNegativeTEFillCheck = GREKernel.getbPerformNegativeTEFillCheck();
  GREKernel.setbPerformNegativeTEFillCheck (false);


  BOOL bStatus = pSeq->prepareForBinarySearch(pMrProt);
  if ( !bStatus )  {
    if ( bOldPerformNegativeTEFillCheck ) GREKernel.setbPerformNegativeTEFillCheck ( bOldPerformNegativeTEFillCheck );

    return ( false );    // * No valid protocol -> return with error *
  }

  GREKernel.setbPerformNegativeTEFillCheck ( bOldPerformNegativeTEFillCheck );
  lStatus = fCalculateTEMin ( pMrProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill );

  for ( lI=0; lI<pMrProt->contrasts(); lI++ )  {
    if ( pMrProt->te()[lI] > alTEMin[lI] ) alTEMin[lI] = pMrProt->te()[lI];
  }


  if ( (lStatus == SEQU__NORMAL) || (lStatus == SEQU__NEGATIV_TEFILL) )  {
    return ( true );
  } else {
    return ( false );
  }

}







// * -------------------------------------------------------------------------- *
// * -------------------------------------------------------------------------- *
// *                                                                            *
// *            Definition of special UI parameter for gre_highres              *
// *                                                                            *
// * -------------------------------------------------------------------------- *
// * -------------------------------------------------------------------------- *
bool fDIsAvailableMaxGradAmpl(LINK_DOUBLE_TYPE* const , long)
{
    return (true);
}




unsigned fDGetLabelIdGetMaxGradAmpl (LINK_DOUBLE_TYPE* const, char* arg_list[], long )
{
    static const char* const pszLabel0    = "RO Grad Ampl.";

    arg_list[0] = (char*)pszLabel0;

    return MRI_STD_STRING;
}



unsigned fDGetUnitId (LINK_DOUBLE_TYPE* const, char* arg_list[], long )
{
    static const char* const pszUnit = "mT/m";
    arg_list[0] = (char*)pszUnit;

    return MRI_STD_STRING;
}



double fDGetValueMaxGradAmpl (LINK_DOUBLE_TYPE* const pThis, long )
{
    return pThis->prot().fastImaging().echoSpacing();
}



double fDSetValueMaxGradAmpl (LINK_DOUBLE_TYPE* const pThis, double value, long )
{
    return (pThis->prot().fastImaging().echoSpacing(value));
}



bool fDGetLimitsHandler (LINK_DOUBLE_TYPE* const, std::vector<MrLimitDouble>& rLimitVector, unsigned long& rulVerify, long )
{
    rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;
    rLimitVector.resize(1);
    rLimitVector[0].setEqualSpaced (10.0, 60.0, 0.1); // min, max, incr

    return true;
}







/*[ Function ****************************************************************\
*
* Name        : segmentsGetLimitsHandler
*
* Description : Segments limits handler dimms the segment parameter if
*               respiratory triggering has been selected
*
* Return      : bool
*
\****************************************************************************/
bool segmentsGetLimitsHandler (LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    SEQ::PhysioSignal                           TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal                           TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod                           TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod                           TrigMethodLow   = SEQ::METHOD_NONE;


    pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

    if ( (TrigSignalHigh == SEQ::SIGNAL_RESPIRATION) && (TrigMethodHigh == SEQ::METHOD_TRIGGERING) ) {

        // dimm the segment parameter
        rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;

        return ( false );

    } else {

        return ( (*pOriginalSegmentsGetLimitsFct)(pThis, rLimitVector, rulVerify, lIndex) );

    }

}




/*[ Function ****************************************************************\
*
* Name        : firstAcqWindowGetLimitsHandler
*
* Description : needs to be equivalent to firstAcqWindowInGetLimitsHandler
*
* Return      : bool
*
\****************************************************************************/
bool firstAcqWindowGetLimitsHandler (LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    SEQ::PhysioSignal                           TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal                           TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod                           TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod                           TrigMethodLow   = SEQ::METHOD_NONE;


    pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

    if ( (TrigSignalHigh == SEQ::SIGNAL_RESPIRATION) && (TrigMethodHigh == SEQ::METHOD_TRIGGERING) ) {

        rLimitVector.resize(1);
        rulVerify = LINK_LONG_TYPE::VERIFY_OFF;

        //  The minimum acquisition window accepted by sequence is equal to
        //  tr_ms * Phases * AcqWindowFactor, where
        //  the AcqWindowFactor is exported by the sequence
        if( !pThis->sequence().prepareForBinarySearch(&pThis->prot()) )
        {
            //  should never happen
            return false;
        }

        long lMinScanWindow_ms = (pThis->prot().tr()[0] + 1000 - 1) / 1000;

            return rLimitVector[0].setEqualSpaced ( lMinScanWindow_ms, 10000, 1 );

    } else {

        return ( (*pOriginalFirstAcqWindowGetLimitsFct)(pThis, rLimitVector, rulVerify, lIndex) );

    }

}

/*[ Function ****************************************************************\
*
* Name        : firstAcqWindowInGetLimitsHandler
*
* Description : needs to be equivalent to firstAcqWindowGetLimitsHandler
*
* Return      : bool
*
\****************************************************************************/
bool firstAcqWindowInGetLimitsHandler (LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    SEQ::PhysioSignal                           TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal                           TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod                           TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod                           TrigMethodLow   = SEQ::METHOD_NONE;


    pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

    if ( (TrigSignalHigh == SEQ::SIGNAL_RESPIRATION) && (TrigMethodHigh == SEQ::METHOD_TRIGGERING) ) {

        rLimitVector.resize(1);
        rulVerify = LINK_LONG_TYPE::VERIFY_OFF;

        //  The minimum acquisition window accepted by sequence is equal to
        //  tr_ms * Phases * AcqWindowFactor, where
        //  the AcqWindowFactor is exported by the sequence
        if( !pThis->sequence().prepareForBinarySearch(&pThis->prot()) )
        {
            //  should never happen
            return false;
        }

        long lMinScanWindow_ms = (pThis->prot().tr()[0] + 1000 - 1) / 1000;

            return rLimitVector[0].setEqualSpaced ( lMinScanWindow_ms, 10000, 1 );

    } else {

        return ( (*pOriginalFirstAcqWindowInGetLimitsFct)(pThis, rLimitVector, rulVerify, lIndex) );

    }

}



/*[ Function ****************************************************************\
*
* Name        : trGetLimitsHandler
*
* Description : TR limits handler considering various multi slice modes:
*
*               sequential : segments are within TR
*                            ==> acquisition window length < TR * phases
*                            ==> can be handled by standard TR limits handler
*
*               interleaved: segments are NOT within TR
*                            ==> acquisition window length < TR * segments
*
* Return      : bool
*
\****************************************************************************/
bool trGetLimitsHandler (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    static const char *ptModule = {"trGetLimitsHandler"};

    const double                                dEpsilon        = 0.001;                // * Used to provide rounding errors *

    SEQ::PhysioSignal                           TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal                           TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod                           TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod                           TrigMethodLow   = SEQ::METHOD_NONE;

    long                                        lDelay_ms       = 0;
    long                                        lScanWindow_ms  = 0;
    long                                        lPhysioPhases   = pThis->prot().physiology().phases();

    double                                      dMaxTR_ms       = 0.0;

    MrLimitDouble                               TRInterval;
    std::vector< MrLimitDouble >                OrigTRLimits;
    std::vector< MrLimitDouble >::size_type     nCntr           = 0;


    // * -------------------------------------------------------------------------- *
    // * Calculation of TR limits                                                   *
    // * -------------------------------------------------------------------------- *
    switch ( pThis->prot().kSpace().multiSliceMode() ) {

        case SEQ::MSM_SEQUENTIAL:
                return ( (*pOriginalTRGetLimtsFct)(pThis, rLimitVector, rulVerify, lIndex) );
        break;

        case SEQ::MSM_INTERLEAVED:

                rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;

                // * -------------------------------------------------------------- *
                // * Determine triger delay time and acquisition window length      *
                // * -------------------------------------------------------------- *
                pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

                switch ( TrigSignalHigh )   {

                    case SEQ::SIGNAL_EKG:
                    case SEQ::SIGNAL_PULSE:
                    case SEQ::SIGNAL_EXT:         lDelay_ms      = pThis->prot().physiology().triggerDelay ( TrigSignalHigh ) / 1000;
                                                  lScanWindow_ms = pThis->prot().physiology().scanWindow   ( TrigSignalHigh );
                    break;

                    case SEQ::SIGNAL_RESPIRATION: lDelay_ms      = 0;
                                                  lScanWindow_ms = pThis->prot().physiology().scanWindow ( TrigSignalHigh );

                    break;

                    case SEQ::SIGNAL_NONE:
                        switch ( TrigSignalLow )  {
                            case SEQ::SIGNAL_EKG:
                            case SEQ::SIGNAL_PULSE:
                            case SEQ::SIGNAL_EXT:         lDelay_ms      = pThis->prot().physiology().triggerDelay ( TrigSignalLow ) / 1000;
                                                          lScanWindow_ms = pThis->prot().physiology().scanWindow   ( TrigSignalLow );
                            break;

                            case SEQ::SIGNAL_RESPIRATION: lDelay_ms      = 0;
                                                          lScanWindow_ms = pThis->prot().physiology().scanWindow ( TrigSignalLow );
                            break;

                            case SEQ:: SIGNAL_NONE:       lDelay_ms      = 0;
                                                          lScanWindow_ms = 0;
                            break;

                            default:
                                TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown physiologic signal 2 = %d  ", ptModule, TrigSignalLow);
                                return ( false );
                        }
                    break;

                    default:
                        TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown physiologic signal 1 = %d  ", ptModule, TrigSignalHigh);
                        return ( false );

                }



                // * -------------------------------------------------------------- *
                // * Determine the maximal TR value                                 *
                // * -------------------------------------------------------------- *
                dMaxTR_ms = pThis->seqLimits().getTR()[lIndex].getMax() / 1000.0;


                    if ( (lScanWindow_ms > 0) && (lPhysioPhases > 0) )    {
                    dMaxTR_ms = min(static_cast<double>(lScanWindow_ms - lDelay_ms) / lPhysioPhases, dMaxTR_ms);
                    }



                // * -------------------------------------------------------------- *
                // * Set the new TR limits                                          *
                // * -------------------------------------------------------------- *
                if(! (*pOriginalTRGetLimtsFct)(pThis, OrigTRLimits, rulVerify, lIndex) )  {  return ( false );  }
                rLimitVector.clear();


                for ( nCntr = 0; nCntr < OrigTRLimits.size(); ++nCntr )  {

                    double dMin  = OrigTRLimits[nCntr].minimum();
                    double dMax  = OrigTRLimits[nCntr].maximum();
                    double dIncr = OrigTRLimits[nCntr].incr();

                    if ( dMax <= dMaxTR_ms )  {     // copy original interval

                        if( TRInterval.setEqualSpaced( dMin, dMax + dEpsilon*dIncr, dIncr, OrigTRLimits[nCntr].color() ) )    {
                            rLimitVector.push_back ( TRInterval );
                        }

                        continue;
                    }


                    if ( (dMin <= dMaxTR_ms) && (dMax >= dMaxTR_ms) )  {  // dMaxTR_ms is within limit

                        if( TRInterval.setEqualSpaced( dMin, dMaxTR_ms + dEpsilon*dIncr, dIncr, OrigTRLimits[nCntr].color() ) )    {
                            rLimitVector.push_back ( TRInterval );
                        }

                        break;
                    }

                }

                return ( rLimitVector.size() > 0 );

        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence multi slice mode = %d  ", ptModule, pThis->prot().kSpace().multiSliceMode());
            return ( false );

    }


}










// * -------------------------------------------------------------------------- *
// *                                                                            *
// * Name        :  teGetLimitsHandler                                          *
// *                                                                            *
// * Description :  Calculates the TE limits                                    *
// *                                                                            *
// * Return      :  bool                                                        *
// *                                                                            *
// * -------------------------------------------------------------------------- *
bool teGetLimitsHandler (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    if ( pThis->prot().contrasts() == 1 )  {

        return (*pOriginalTEGetLimtsFct)(pThis, rLimitVector, rulVerify, lIndex);

    } else {

        if ( bUseTEBinarySearch )  {

            // * -------------------------------------------------------------- *
            // * Set TE limits and enable binary search                         *
            // * -------------------------------------------------------------- *
            rulVerify = LINK_DOUBLE_TYPE::VERIFY_BINARY_SEARCH;
            const ParLimArray<long>& _seqLimits = pThis->seqLimits().getTE();
            MrLimitDouble TEInterval;
            TEInterval.setEqualSpaced( _seqLimits[lIndex].getMin()/1000.0, _seqLimits[lIndex].getMax()/1000.0, _seqLimits[lIndex].getInc()/1000.0 );

            rLimitVector.clear();
            rLimitVector.push_back ( TEInterval );

            return ( rLimitVector.size() > 0 );

        } else {

            MrProt   *pMrProt  = (MrProt*)  &pThis->prot();
            SeqLim   *pSeqLim  = (SeqLim*)  &pThis->seqLimits();
            Sequence *pSeq     = &pThis->sequence();

            long lIncr = pSeqLim->getTE()[lIndex].getInc();
            long alTEFill[lNoOfContrasts];
            long alTEMin[lNoOfContrasts];
            long lTRFill = 0;
            long lI;



            // * -------------------------------------------------------------- *
            // * In principle there might be three TE intervalls:               *
            // *                                                                *
            // *               RedLow      green       RedHigh                  *
            // *            |___________|___________|___________|               *
            // *                             ^                                  *
            // *                             |                                  *
            // *                            TE                                  *
            // * -------------------------------------------------------------- *
            MrLimitDouble RedLow, Green, RedHigh;

            long lRedLowLeft, lRedLowRight;
            long lGreenLeft, lGreenRight;
            long lRedHighLeft, lRedHighRight;



            // * -------------------------------------------------------------- *
            // * Prepare the sequence to get the                                *
            // *    - minimum TE times                                          *
            // *    - TE fill times                                             *
            // *    - TR fill times                                             *
            // * -------------------------------------------------------------- *
            pSeq->prepareForBinarySearch( pMrProt );



            // * -------------------------------------------------------------- *
            // * Disable binary search in order to increase performance         *
            // * -------------------------------------------------------------- *
            rulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;



            // * -------------------------------------------------------------- *
            // * Reset limits vector                                            *
            // * -------------------------------------------------------------- *
            rLimitVector.clear();




            for ( lI = 0; lI < pMrProt->contrasts(); lI++ )  {
                alTEMin[lI]  = GREKernel.getalTEMin()[lI];
                alTEFill[lI] = GREKernel.getalTEFill()[lI];
            }

            // * -------------------------------------------------------------- *
            // * Calculate TRFill: minimum of all TRFill times from the various *
            // * concatenations                                                 *
            // * -------------------------------------------------------------- *
            lTRFill = RunLoop.getlTRFillInConcat(0);  // TRFill from first concat
            for ( lI = 1; lI < pMrProt->concatenations(); lI++ )  {
                lTRFill = minimum ( lTRFill, RunLoop.getlTRFillInConcat(lI) );
            }



            // * -------------------------------------------------------------- *
            // * Calculate lower red limits                                     *
            // * -------------------------------------------------------------- *
            if  ( lIndex >= 1 )  {

                lRedLowLeft = maximum (alTEMin[lIndex], pSeqLim->getTE()[lIndex].getMin());
                lRedLowLeft = minimum (lRedLowLeft    , pSeqLim->getTE()[lIndex].getMax());


                lRedLowRight = alTEMin[lIndex];
                for ( lI = 0; lI < lIndex; lI++ )  {
                    lRedLowRight += maximum (alTEFill[lI], 0L);
                }
                lRedLowRight -= lIncr;

                // TE has to be within hard limits
                lRedLowRight = maximum (lRedLowRight, pSeqLim->getTE()[lIndex].getMin());
                lRedLowRight = minimum (lRedLowRight, pSeqLim->getTE()[lIndex].getMax());

                if ( lRedLowRight > lRedLowLeft )  {
                    // * ---------------------------------------------------------- *
                    // * Add red intervall                                          *
                    // * ---------------------------------------------------------- *
                    RedLow.setEqualSpaced  (lRedLowLeft/1000.0, lRedLowRight/1000.0, lIncr/1000.0, MrLimit<double>::RED);
                    rLimitVector.push_back (RedLow);
                }

            }


            // * -------------------------------------------------------------- *
            // * Calculate green limits                                         *
            // * -------------------------------------------------------------- *
            // Left border of green limit
            lGreenLeft = alTEMin[lIndex];
            for ( lI = 0; lI < lIndex; lI++ )  {
                lGreenLeft += maximum (alTEFill[lI], 0L);
            }

            // TE has to be within hard limits
            lGreenLeft = maximum (lGreenLeft, pSeqLim->getTE()[lIndex].getMin());
            lGreenLeft = minimum (lGreenLeft, pSeqLim->getTE()[lIndex].getMax());



            // Left border of green limit
            lGreenRight = alTEMin[lIndex];
            for ( lI = 0; lI <= lIndex; lI++  )  {
                lGreenRight += maximum (alTEFill[lI], 0L);
            }

            // Consider fill time of succeeding contrast or TR fill time
            if ( lIndex < pMrProt->contrasts()-1 )  {
                lGreenRight += maximum (alTEFill[lIndex+1], 0L);
            } else {
                lGreenRight += maximum (lTRFill, 0L);
            }

            // TE has to be within hard limits
            lGreenRight = maximum (lGreenRight, pSeqLim->getTE()[lIndex].getMin());
            lGreenRight = minimum (lGreenRight, pSeqLim->getTE()[lIndex].getMax());



            // * -------------------------------------------------------------- *
            // * Add green intervall                                            *
            // * -------------------------------------------------------------- *
            Green.setEqualSpaced   (lGreenLeft/1000.0, lGreenRight/1000.0, lIncr/1000.0, MrLimit<double>::GREEN);
            rLimitVector.push_back (Green);


            // * -------------------------------------------------------------- *
            // * Calculate red high limits                                      *
            // * TR changes are not implemented up to now                       *
            // * -------------------------------------------------------------- *
            if ( lIndex < pMrProt->contrasts()-1 )  {

                // Left border of the upper red limit
                lRedHighLeft = lGreenRight + lIncr;

                // TE has to be within hard limits
                lRedHighLeft = maximum (lRedHighLeft, pSeqLim->getTE()[lIndex].getMin());
                lRedHighLeft = minimum (lRedHighLeft, pSeqLim->getTE()[lIndex].getMax());


                // Calculate maximum TE of contrast lIndex using all TE fill times
                // and the TR fill time
                lRedHighRight = alTEMin[pMrProt->contrasts()-1];
                for ( lI = 0; lI < pMrProt->contrasts()-1; lI++ )  {
                    lRedHighRight += maximum (alTEFill[lI], 0L);
                }
                lRedHighRight += maximum (lTRFill, 0L);

                // The last contrast must not exeed its TE hard limit
                lRedHighRight =  minimum (lRedHighRight, pSeqLim->getTE()[pMrProt->contrasts()-1].getMax());
                lRedHighRight -= (alTEMin[pMrProt->contrasts()-1] - alTEMin[lIndex]);



                // TE has to be within hard limits
                lRedHighRight = maximum (lRedHighRight, pSeqLim->getTE()[lIndex].getMin());
                lRedHighRight = minimum (lRedHighRight, pSeqLim->getTE()[lIndex].getMax());


                // * ---------------------------------------------------------- *
                // * Add upper red intervall                                    *
                // * ---------------------------------------------------------- *
                if ( lRedHighLeft <= lRedHighRight )  {
                    RedHigh.setEqualSpaced (lRedHighLeft/1000.0, lRedHighRight/1000.0, lIncr/1000.0, MrLimit<double>::RED);
                    rLimitVector.push_back (RedHigh);
                }

            }
#ifdef SUPPORT_CT
            MrProt* pProt = &pThis->prot();
            if( pProt->dixon() != SEQ::DIXON_NONE && rLimitVector.size() > 0 )
            {

                MeasNucleus sNucleus(pProt->txSpec().nucleusInfoArray()[0].nucleus());
                const double dConst  = -.5*1e6/(sNucleus.getLarmorConst()*SysProperties::getNominalB0()*CHEMICAL_SHIFT_FAT_PPM);

                std::vector<MrLimitDouble> sDixonLim;
                sDixonLim.reserve(1+static_cast<int>((rLimitVector.back().maximum()-rLimitVector.front().minimum())/(2*dConst)));

                //  Oppsed phase echo
                int iCntr = lIndex == 0 ? 1 : 2;

                double dTE_ms;

                std::vector<MrLimitDouble>::iterator it = rLimitVector.begin();
                MrLimitDouble sInterval;
                for(;it != rLimitVector.end();iCntr += 2)
                {
                    dTE_ms = iCntr*dConst/1000.;
                    while( dTE_ms > (*it).maximum() )
                    {
                        if( ++it == rLimitVector.end() ) break;
                    }
                    if( dTE_ms >= (*it).minimum() && dTE_ms <= (*it).maximum() )
                    {
                        dTE_ms = (*it).minimum()+(*it).incr()*static_cast<int>(0.5+(dTE_ms-(*it).minimum())/(*it).incr());

                        sInterval.setLonely(dTE_ms,(*it).color());
                        sDixonLim.push_back(sInterval);
                    }
                }
                rLimitVector = sDixonLim;
            }
#endif
            return ( rLimitVector.size() > 0 );

        }

    }

}




/*[ Function ****************************************************************\
*
* Name        : sliceOSGetLimitsHandler
*
* Description : This get limits handler just calls the original get limits
*               handler and sets the binary search to VERIFY_SCAN_ALL. This
*               is necessary since the minimum TE might slightly vary depending
*               on the OS value and cause internal errors.
*
* Return      : bool
*
\****************************************************************************/
bool sliceOSGetLimitsHandler (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    bool bResOrigHandler = (*pOriginalSliceOSGetLimtsHandler)(pThis, rLimitVector, rulVerify, lIndex);

    if ( pThis->prot().kSpace().phasePartialFourierFactor() != SEQ::PF_OFF )  {
        rulVerify = LINK_DOUBLE_TYPE::VERIFY_SCAN_ALL;  // Check every value
    }

    return ( bResOrigHandler );
}






/*[ Function ****************************************************************\
*
* Name        : imaPerSlabGetLimitsHandler
*
* Description : This get limits handler just calls the original get limits
*               handler and sets the binary search to VERIFY_SCAN_ALL. This
*               is necessary since the minimum TE might slightly vary depending
*               on the number of images per slab and cause internal errors.
*
* Return      : bool
*
\****************************************************************************/
bool imaPerSlabGetLimitsHandler (LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    bool bResOrigHandler = (*pOriginalImaPerSlabGetLimitsHandler)(pThis, rLimitVector, rulVerify, lIndex);

    if ( pThis->prot().kSpace().phasePartialFourierFactor() != SEQ::PF_OFF )  {
        rulVerify = LINK_LONG_TYPE::VERIFY_SCAN_ALL;    // Check every value
    }

    return ( bResOrigHandler );
}



/*[ Function ****************************************************************\
*
* Name        : flipAngleGetLimitsHandler
*
* Description :
*
* Return      : bool
*
\****************************************************************************/
bool flipAngleGetLimitsHandler (LINK_DOUBLE_TYPE* const pThis, std::vector<MrLimitDouble>& rLimitVector, unsigned long& rulVerify, long lIndex)
{

    if (    (pThis->prot().kSpace().dimension() == SEQ::DIM_3)
         && ((pThis->prot().preparationPulses().fatSuppression() == SEQ::WATER_EXCITATION) || (pThis->prot().preparationPulses().fatSuppression() == SEQ::WATER_EXCITATION_FAST) ) )  {

        rulVerify = LINK_LONG_TYPE::VERIFY_BINARY_SEARCH;

        rLimitVector.resize(2);

        rLimitVector[0].setEqualSpaced (pThis->seqLimits().getFlipAngle().getMin(), dWEFlipAngleLimit + pThis->seqLimits().getFlipAngle().getInc() / 2, pThis->seqLimits().getFlipAngle().getInc());
        rLimitVector[1].setEqualSpaced (dWEFlipAngleLimit + pThis->seqLimits().getFlipAngle().getInc(), pThis->seqLimits().getFlipAngle().getMax() + pThis->seqLimits().getFlipAngle().getInc() / 2, pThis->seqLimits().getFlipAngle().getInc());

        return rLimitVector.size() > 0;

    } else {
        return ((*pOriginalFlipAngleGetLimitsHandler)(pThis, rLimitVector, rulVerify, lIndex));
    }

}



/*[ Function ****************************************************************\
*
* Name        : phasesIsAvailableHandler
*
* Description : Is available handler for the phases parameter
*               In case of respiratory triggering phases will be hidden.
*
* Return      : bool
*
\****************************************************************************/
bool phasesIsAvailableHandler (LINK_LONG_TYPE* const pThis,long lIndex)
{

    SEQ::PhysioSignal                           TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal                           TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod                           TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod                           TrigMethodLow   = SEQ::METHOD_NONE;

    pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

    if ( (TrigSignalHigh == SEQ::SIGNAL_RESPIRATION) && (TrigMethodHigh == SEQ::METHOD_TRIGGERING) )  {

        return ( false );   // hide phases parameter

    } else {
        if ( pOriginalPhasesIsAvailableFct ) {
            return (*pOriginalPhasesIsAvailableFct)(pThis, lIndex);
        } else {
            return ( false );
        }
    }

}








/*[ Function ****************************************************************\
*
* Name        : firstAcqWindowToolTipIdHandler
*
* Description : Tool ID handler handler for the first acquisition window
*               respiratory triggering : Shows the number of segments that
*                                        result from the acq. window length
*               else                   : The standard tool tip ID handler
*                                        will be used
*
* Return      : unsigned
*
\****************************************************************************/
unsigned firstAcqWindowToolTipIdHandler (LINK_LONG_TYPE* const pThis, char** arg_list, long lIndex)
{

    SEQ::PhysioSignal                           TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal                           TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod                           TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod                           TrigMethodLow   = SEQ::METHOD_NONE;

    pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

    if ( (TrigSignalHigh == SEQ::SIGNAL_RESPIRATION) && (TrigMethodHigh == SEQ::METHOD_TRIGGERING) )  {

        static char tLine    [100];
        static char tToolTip[1100];
        const long  lMaxNoOfLines       = 11;
        const long  lTR_ms              = (pThis->prot().tr()[0] + 999 ) / 1000;
        const long  lActualAcqWindow    = (pThis->prot().fastImaging().segments() * pThis->prot().tr()[0] + 999) / 1000;
        const long  lAcqWindowMax_ms    = 10000;
        long        lMinAcqWindow       = 0;
        long        lMaxAcqWindow       = 0;
        long        lSegments           = 0;

        static const char* const tFormat = "\t%1!r! %2!r!\t\t%3!r! \n%4!S!";


        // * ----------------------------------------------------------------------- *
        // * Prepare the sequence                                                    *
        // * ----------------------------------------------------------------------- *
        if( !pThis->sequence().prepareForBinarySearch ( &pThis->prot() ) )  {    return ( 0 );    }



        arg_list[0] = (char*) tFormat;
        arg_list[1] = (char*) MRI_STD_ACQUISITION_WINDOW_LABEL;
        arg_list[2] = (char*) MRI_STD_UNIT_MS;
        arg_list[3] = (char*) MRI_STD_SEGMENTS_LABEL;

        tToolTip[0] = '\0';
        sprintf ( tLine, "----------------------------------------------------------------------------------\n" ); strcat ( tToolTip, tLine );


        lSegments     = maximum (pThis->seqLimits().getSegments().getMin(), pThis->prot().fastImaging().segments() - (lMaxNoOfLines/2) * pThis->seqLimits().getSegments().getInc());
        lMinAcqWindow = lSegments * lTR_ms;
        lMaxAcqWindow = lMinAcqWindow + pThis->seqLimits().getSegments().getInc() * lTR_ms - 1;


        if ( (lMinAcqWindow <= lActualAcqWindow) && (lMaxAcqWindow >= lActualAcqWindow) )  {
            sprintf ( tLine, ">>>\t%8ld    -- %7ld\t\t%8ld\n", lMinAcqWindow, lMaxAcqWindow, lSegments); strcat ( tToolTip, tLine );
        } else   {
            sprintf ( tLine, "\t%8ld    -- %7ld\t\t%8ld\n", lMinAcqWindow, lMaxAcqWindow, lSegments); strcat ( tToolTip, tLine );
        }

        for ( long lI = 0; lI < lMaxNoOfLines-1; lI++ )  {

            lSegments     += pThis->seqLimits().getSegments().getInc();
            lMinAcqWindow = lSegments * lTR_ms;
            lMaxAcqWindow = lMinAcqWindow + pThis->seqLimits().getSegments().getInc() * lTR_ms - 1;

            if ( lMinAcqWindow > lAcqWindowMax_ms )  {  break;  }
            if ( lMaxAcqWindow > lAcqWindowMax_ms )  {  lMaxAcqWindow = lAcqWindowMax_ms;  }

            if ( (lMinAcqWindow <= lActualAcqWindow) && (lMaxAcqWindow > lActualAcqWindow) )  {
                sprintf ( tLine, ">>>\t%8ld    -- %7ld\t\t%8ld\n", lMinAcqWindow, lMaxAcqWindow, lSegments); strcat ( tToolTip, tLine );
            } else   {
                sprintf ( tLine, "\t%8ld    -- %7ld\t\t%8ld\n", lMinAcqWindow, lMaxAcqWindow, lSegments); strcat ( tToolTip, tLine );
            }

        }


        arg_list[4] = tToolTip;

        return ( MRI_STD_STRING );   // hide phases parameter

    } else {
        if ( pOriginalFirstAcqWindowToolTipIdHandler ) {
            return (*pOriginalFirstAcqWindowToolTipIdHandler)(pThis, arg_list, lIndex);
        } else {
            return ( 0 );
        }
    }

}
// * ------------------------------------------------------------------------ *
// * Declaration of set functions                                             *
// * ------------------------------------------------------------------------ *




template<class TYPE>
bool _fLocalUICSetWithinLimits(MrUILinkBase* pThis, TYPE& rVal, const char* tTag, unsigned long ulMode, long lIndex, MrUILinkLimited<TYPE>::PFctSetValue pFctSet)
{
  const double dTolerance = 1.0e-6;

        MrUILinkLimited<TYPE>* pParam   = NULL;
        //  search the parameter
        if(lIndex < 0)
        {
                pParam = _search< MrUILinkLimited<TYPE> >(pThis,tTag);
        }
        else
        {
                pParam = _searchElm< MrUILinkLimited<TYPE> >(pThis,tTag);
        }
        if( !pParam || !pParam->isAvailable(lIndex) )
        {
                //  parameter not found or not available
                return false;
        }

        std::vector< MrLimit<TYPE> > _limits;

        unsigned long ulVerify = 0;
    if(pFctSet)
    {
        //  if no set-value-handler has been set, the get-limits-handler isn't
        //  invoked by the getLimits-Function
        //  so we retrieve the get-limits-handler and invoke it explicitly
        MrUILinkLimited<TYPE>::PFctGetLimits pLimitHandler = pParam->registerGetLimitsHandler(NULL);
        //  restore the get-limits-handler
        pParam->registerGetLimitsHandler(pLimitHandler);
        //  invoke the get-limits-handler
        if(pLimitHandler)
        {
            if(!(*pLimitHandler)(pParam,_limits, ulVerify,lIndex))
            {
                _limits.resize(1);
                _limits[0].setLonely(pParam->value(lIndex));
            }
        }
        else
        {
            pParam->getLimits(_limits, ulVerify,lIndex);
        }
    }
    else
    {
        pParam->getLimits(_limits, ulVerify,lIndex);
    }
        if(_limits.size() < 1)
        {
                //  should never happen
                return false;
        }
    TYPE inVal = rVal;
        for(std::vector< MrLimit<TYPE> >::size_type nCntr = 0; nCntr < _limits.size(); ++nCntr)
        {
                if(rVal < _limits[nCntr].minimum())
                {
                        if(!nCntr || (ulMode & UIC_SET_WITHIN_LIMITS_ROUNDUP))
                        {
                                rVal = _limits[nCntr].minimum();
                        }
                        else if(ulMode & UIC_SET_WITHIN_LIMITS_ROUNDDOWN)
                        {
                                rVal = _limits[nCntr-1].maximum();
                        }
                        else
                        {
                                rVal = rVal-_limits[nCntr-1].maximum() < _limits[nCntr].minimum() - rVal ? _limits[nCntr-1].maximum() : _limits[nCntr].minimum();
                        }
                        rVal = pFctSet ? (*pFctSet)(pParam,rVal,lIndex) : pParam->value(rVal,lIndex);
            break;
                }
                if(rVal <= _limits[nCntr].maximum())
                {
                        if(ulMode & UIC_SET_WITHIN_LIMITS_ROUNDUP)
                        {
                                rVal = _limits[nCntr].maximum() - _limits[nCntr].incr()*static_cast<int>(static_cast<double>(_limits[nCntr].maximum()-rVal)/_limits[nCntr].incr() + dTolerance );
                        }
                        else if(ulMode & UIC_SET_WITHIN_LIMITS_ROUNDDOWN)
                        {
                                rVal = _limits[nCntr].minimum() + _limits[nCntr].incr()*static_cast<int>(static_cast<double>(rVal-_limits[nCntr].minimum())/_limits[nCntr].incr() + dTolerance );
                        }
                        else
                        {
                                rVal = _limits[nCntr].minimum() + _limits[nCntr].incr()*static_cast<int>(static_cast<double>(rVal-_limits[nCntr].minimum()+_limits[nCntr].incr()/2.)/_limits[nCntr].incr());
                        }
                        rVal = pFctSet ? (*pFctSet)(pParam,rVal,lIndex) : pParam->value(rVal,lIndex);
                        break;
                }
        }
    if(nCntr == _limits.size())
    {
        //  The input value is greater than the maximum of the 'outermost' limit interval
        rVal = pFctSet ? (*pFctSet)(pParam,_limits.back().maximum(),lIndex) : pParam->value(_limits.back().maximum(),lIndex);
    }
    if(ulMode & UIC_SET_WITHIN_LIMITS_ROUNDUP)   return rVal >= inVal;
    if(ulMode & UIC_SET_WITHIN_LIMITS_ROUNDDOWN) return rVal <= inVal;
        return true;
}







/*[ Function ****************************************************************\
*
* Name        : switchOffInversionRecovery
*
* Description : This function switches off inversion recovery and
*               formats the output of a solve hander pop-up
*
*
* Return      : 0 if an error occurs
*               1 otherwise
*
\****************************************************************************/
unsigned switchOffInversionRecovery ( MrUILinkBase* const   pThis,
                                      char**                arg_list,  // receives confirmation message
                                      const MrProt*         pOrigProt, // Original protocol with old rf mode
                                      long&                 lCtrl
                                     )
{

    const long lOffset = 20;

    if ( arg_list && (lCtrl <= 3) )  {
        arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_LABEL;
        switch ( pOrigProt->preparationPulses().inversion() )  {
            case SEQ::SLICE_SELECTIVE:  arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_IR_SS;   break;
            case SEQ::VOLUME_SELECTIVE: arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_IR_NS;   break;
        }
        arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_NONE;
        arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
        arg_list[39 + lCtrl * lOffset] = (char*) '\n';
        lCtrl++;
    }


    // * ---------------------------------------------------------------------- *
    // * Switch off magnetization preparation                                   *
    // * ---------------------------------------------------------------------- *
    LINK_SELECTION_TYPE* pInversionMode = _search<LINK_SELECTION_TYPE> (pThis, MR_TAG_INVERSION );

    if(!pInversionMode || !pInversionMode->isEditable(0) || !pThis->seqLimits().getInversion().hasOption(SEQ::INVERSION_OFF))  {  return ( 0 );  }
    if( pInversionMode->value(MRI_STD_MAGN_PREP_NONE,0) != MRI_STD_MAGN_PREP_NONE)  {  return ( 0 );  }

    return ( 1 );

}










/*[ Function ****************************************************************\
*
* Name        : switchOffSaturationRecovery
*
* Description : This function switches off saturation recovery and
*               formats the output of a solve hander pop-up
*
*
* Return      : 0 if an error occurs
*               1 otherwise
*
\****************************************************************************/
unsigned switchOffSaturationRecovery ( MrUILinkBase* const   pThis,
                                       char**                arg_list,  // receives confirmation message
                                       const MrProt*         pOrigProt, // Original protocol with old rf mode
                                       long&                 lCtrl
                                     )
{

    const long lOffset = 20;

    if ( arg_list && (lCtrl <= 3) )  {

        arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_LABEL;
        switch ( pOrigProt->preparationPulses().satRecovery() )  {
            case SEQ::SATREC_SLICE_SELECTIVE:   arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_SR_SS;   break;
            case SEQ::SATREC_VOLUME_SELECTIVE:  arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_SR_NS;   break;
        }
        arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_MAGN_PREP_NONE;
        arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
        arg_list[39 + lCtrl * lOffset] = (char*) '\n';
        lCtrl++;

    }

    // * -------------------------------------------------------------- *
    // * Switch off saturation recovery                                 *
    // * -------------------------------------------------------------- *
    LINK_SELECTION_TYPE* pSatRecoveryMode = _search<LINK_SELECTION_TYPE> (pThis, MR_TAG_INVERSION );

    if(!pSatRecoveryMode || !pSatRecoveryMode->isEditable(0) || !pThis->seqLimits().getSaturationRecovery().hasOption(SEQ::SATREC_NONE))  {  return ( 0 );  }
    if( pSatRecoveryMode->value(MRI_STD_MAGN_PREP_NONE,0) != MRI_STD_MAGN_PREP_NONE)  {  return ( 0 );  }

    return ( 1 );

}






/*[ Function ****************************************************************\
*
* Name        : switchOffDarkBlood
*
* Description : This function switches off dark blood preparation and
*               formats the output of a solve hander pop-up
*
*
* Return      : 0 if an error occurs
*               1 otherwise
*
\****************************************************************************/
unsigned switchOffDarkBlood ( MrUILinkBase* const   pThis,
                              char**                arg_list,   // receives confirmation message
                              const MrProt*,                    // Original protocol with old rf mode
                              long&                 lCtrl
                            )
{

    const long lOffset = 20;

    if ( arg_list && (lCtrl <= 3) )  {
        arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_DARK_BLOOD_LABEL;
        arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_ON;
        arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_OFF;
        arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
        arg_list[39 + lCtrl * lOffset] = (char*) '\n';
        lCtrl++;
    }

    // * -------------------------------------------------------------- *
    // * Switch off darkblood preparation                               *
    // * -------------------------------------------------------------- *
    LINK_BOOL_TYPE*      pDarkBloodMode = _search<LINK_BOOL_TYPE> (pThis, MR_TAG_DARK_BLOOD);

    if(!pDarkBloodMode || !pDarkBloodMode->isEditable(0) || !pThis->seqLimits().getDarkBlood().hasOption(SEQ::OFF))  {  return ( 0 );  }
    if( pDarkBloodMode->value(false,0) != false)  {  return ( 0 );  }

    return ( 1 );
}






/*[ Function ****************************************************************\
*
* Name        : changeNoOfSegments
*
* Description : This function changes the number of segments and
*               formats the output of a solve hander pop-up
*
*
* Return      : 0 if an error occurs
*               1 otherwise
*
\****************************************************************************/
unsigned changeNoOfSegments ( MrUILinkBase* const   pThis,
                              char**                arg_list,  // receives confirmation message
                              const MrProt*         pOrigProt, // Original protocol with old rf mode
                              long                  lNoOfSegments,
                              long&                 lCtrl
                            )
{

    const long lOffset = 20;

    if ( arg_list && (lCtrl <= 3) )  {

        arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_SEGMENTS_LABEL;
        arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_INT;
        arg_list[25 + lCtrl * lOffset] = (char*) pOrigProt->fastImaging().segments();
        arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_INT;
        arg_list[31 + lCtrl * lOffset] = (char*) lNoOfSegments;
        arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
        arg_list[39 + lCtrl * lOffset] = (char*) '\n';
        lCtrl++;

    }


    // * ------------------------------------------------------------------ *
    // * Set new value                                                      *
    // * ------------------------------------------------------------------ *
    LINK_LONG_TYPE* pSegments = _search<LINK_LONG_TYPE>(pThis, MR_TAG_SEGMENTS);

    // * Check whether segmented can be selected *
    if(!pSegments || !pSegments->isEditable(0) || pThis->seqLimits().getSegments().getMax() < lNoOfSegments)   {  return ( 0 );  }

    // * Set sequetial slice mode *
    if(pSegments->value(lNoOfSegments,0) != lNoOfSegments) {  return ( 0 );  }

    return ( 1 );

}






/*[ Function ****************************************************************\
*
* Name        : basicSetValueHandler
*
* Description : This function is a template set handler that can be for
*               various UI parameters.
*
*               It calls the original UI set handler to set the relevant UI
*               parameter. Afterwards the number of segments is automatically
*               maximized if respiratory triggering has been selected by
*               the user
*
*
* Return      : New value of the relevant UI parameter
*
\****************************************************************************/
template <class ParamType, class Param, class OrigFuncType>
Param basicSetValueHandler ( ParamType         pThis,
                             Param             desiredValue,
                             long              lIndex,
                             OrigFuncType      pOriginalSetValueHandler)
{


    SEQ::PhysioSignal   TrigSignalHigh  = SEQ::SIGNAL_NONE;     // * Selected trigger modes    *
    SEQ::PhysioSignal   TrigSignalLow   = SEQ::SIGNAL_NONE;     // * and signals               *
    SEQ::PhysioMethod   TrigMethodHigh  = SEQ::METHOD_NONE;
    SEQ::PhysioMethod   TrigMethodLow   = SEQ::METHOD_NONE;


    // * -------------------------------------------------------------------------- *
    // * Set the UI parameter using the original set function                       *
    // * -------------------------------------------------------------------------- *
    if( pOriginalSetValueHandler )  {
        desiredValue = (*pOriginalSetValueHandler)( pThis, desiredValue, lIndex );
    }



        // * -------------------------------------------------------------------------- *
        // * Set number of segments to the maximum if respiratory triggering is         *
        // * enabled.                                                                   *
        // * -------------------------------------------------------------------------- *
    pThis->prot().physiology().getPhysioMode (TrigSignalHigh, TrigMethodHigh, TrigSignalLow, TrigMethodLow);

        if ( (TrigSignalHigh == SEQ::SIGNAL_RESPIRATION) && (TrigMethodHigh == SEQ::METHOD_TRIGGERING) )  {

                double dCurrentAcqWindow_ms = pThis->prot().physiology().scanWindow(TrigSignalHigh) - pThis->prot().physiology().triggerDelay(TrigSignalHigh);
                long   lSegments            = (long) (dCurrentAcqWindow_ms * 1000.0 / (double)pThis->prot().tr()[0]);  // 1000.0 is needed because TR is in us.


        lSegments = shift_to_grid_low (lSegments,
                                       pThis->seqLimits().getSegments().getMin(),
                                       pThis->seqLimits().getSegments().getMax(),
                                       pThis->seqLimits().getSegments().getInc());


            //  search the parameter
        LINK_LONG_TYPE      *pSegments = _search< LINK_LONG_TYPE > ( pThis, MR_TAG_SEGMENTS );
        if (pSegments) {  pSegments->value(lSegments,lIndex);  }

        }



    return ( desiredValue );

}








/*[ Function ****************************************************************\
*
* Name        : firstPhysioSignalModeSetValueHandler
*
* Description : It calls the original UI set handler to set the relevant UI
*               parameter. Afterwards the number of segments is automatically
*               maximized if respiratory triggering has been selected by
*               the user
*
* Return      : New value of the relevant UI parameter
*
\****************************************************************************/
unsigned firstPhysioSignalModeSetValueHandler ( LINK_SELECTION_TYPE* const pThis, unsigned lDesiredPhysioMode, long lIndex)
{
    return ( basicSetValueHandler ( pThis, lDesiredPhysioMode, lIndex, pOriginalFirstPhysioSignalModeSetFct ) );
}




/*[ Function ****************************************************************\
*
* Name        : firstAcqWindowSetValueHandler, firstAcqWindowInSetValueHandler
*
* Description : It calls the original UI set handler to set the relevant UI
*               parameter. Afterwards the number of segments is automatically
*               maximized if respiratory triggering has been selected by
*               the user
*               firstAcqWindowInSetValueHandler is needed for the CaptureCycle
*               and must do the same thing as firstAcqWindowSetValueHandler
*
* Return      : New value of the relevant UI parameter
*
\****************************************************************************/
long firstAcqWindowSetValueHandler ( LINK_LONG_TYPE* const pThis, long lDesiredAcqWindowLength, long lIndex)
{
    return ( basicSetValueHandler ( pThis, lDesiredAcqWindowLength, lIndex, pOriginalFirstAcqWindowSetFct ) );
}

long firstAcqWindowInSetValueHandler ( LINK_LONG_TYPE* const pThis, long lDesiredAcqWindowLength, long lIndex)
{
    return ( basicSetValueHandler ( pThis, lDesiredAcqWindowLength, lIndex, pOriginalFirstAcqWindowInSetFct ) );
}




/*[ Function ****************************************************************\
*
* Name        : trSetValueHandler
*
* Description : It calls the original UI set handler to set the relevant UI
*               parameter. Afterwards the number of segments is automatically
*               maximized if respiratory triggering has been selected by
*               the user
*
* Return      : New value of the relevant UI parameter
*
\****************************************************************************/
double trSetValueHandler ( LINK_DOUBLE_TYPE* const   pThis,
                           double                    dDesiredTRValue,
                           long                      lIndex
                         )
{
    return ( basicSetValueHandler ( pThis, dDesiredTRValue, lIndex, pOriginalTRSetFct ) );
}





/*[ Function ****************************************************************\
*
* Name        : teSetValueHandler
*
* Description : TE set value handler
*               Perfomance optimized for multiecho applications (no call
*               of getLimits handler in set handler)
*
* Return      : New value of the relevant UI parameter
*
\****************************************************************************/
double teSetValueHandler ( LINK_DOUBLE_TYPE* const   pThis,
                           double                    dDesiredTRValue,
                           long                      lIndex
                         )
{
    pThis->prot().te()[lIndex] = _flt2Int(dDesiredTRValue * 1000.0);

    return static_cast<double>(pThis->prot().te()[lIndex])/1000.0;
}



/*[ Function ****************************************************************\
*
* Name        : contrastSetValueHandler
*
* Description : This function calls the original UI method to set the number
*               of contrasts
*
* Return      :
*
\****************************************************************************/
long contrastSetValueHandler
(
    LINK_LONG_TYPE* const   pThis,
    long                    lDesiredContrasts,
    long                    lIndex
)
{

    MrProt      *pMrProt        = (MrProt*)  &pThis->prot();
    SeqLim      *pSeqLim        = (SeqLim*)  &pThis->seqLimits();
    SeqExpo     *pSeqExpo       = &pThis->sequence().getSEQ_BUFF()->getSeqExpo ();
    Sequence    *pSeq           = &pThis->sequence();

    long        alTEMin[lNoOfContrasts];
    long        alTEFill[lNoOfContrasts];
    long        lTotalTEFill    = 0;
    bool        bAnErrorOccured = false;
    double      dTEMinKernel_ms = 0.0;
    double      dTEMinUI_ms     = 0.0;



    // * ---------------------------------------------------------------------- *
    // * Contrasts before protocol modifications                                *
    // * ---------------------------------------------------------------------- *
    long lPrevContrasts = pThis->prot().contrasts();


    // * ---------------------------------------------------------------------- *
    // * Set the number of contrast using the original set function             *
    // * ---------------------------------------------------------------------- *
    if( pOriginalContrastSetFct )  {
        lDesiredContrasts = ( *pOriginalContrastSetFct )(pThis, lDesiredContrasts, lIndex);
    }



    // * ---------------------------------------------------------------------- *
    // * Use the bandwidth of the last contrast for all new contrasts           *
    // * ---------------------------------------------------------------------- *
    if ( lDesiredContrasts > lPrevContrasts )  {

        long lI = lPrevContrasts;

        while ( lI < lDesiredContrasts )  {
            pThis->prot().rxSpec().realDwellTime()[lI] = pThis->prot().rxSpec().realDwellTime()[lI-1];
            lI++;
        }

    }

    // If the T2-star map has been selected then put the emasurements to 1
    if(pMrProt->parametricMapping().parametricMapMode() == SEQ::PMAP_T2STAR_MAP) pMrProt->repetitions(0);





    // * ---------------------------------------------------------------------- *
    // * Update exports after protocol has been changed                         *
    // * ---------------------------------------------------------------------- *
    pSeq->prepareForBinarySearch( pMrProt );


    // * ---------------------------------------------------------------------- *
    // * Calculate the min. TE and TR for the selected base matrix              *
    // * ---------------------------------------------------------------------- *
    if ( SEQU__NORMAL !=  fCalculateTEMin (pMrProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill) ) bAnErrorOccured = true;


    // * ------------------------------------------------------------------ *
    // * Calculate TE fill times of existing contrasts                      *
    // * ------------------------------------------------------------------ *
    for ( long lI = 0; lI < lPrevContrasts; lI++ )  {
        lTotalTEFill += alTEFill[lI];
    }

    // * ------------------------------------------------------------------ *
    // * Set TE times of new contrasts, minimum TE times should be used     *
    // * ------------------------------------------------------------------ *
    for ( lI = lPrevContrasts; lI < lDesiredContrasts; lI++ )  {

        dTEMinKernel_ms = dTEMinUI_ms = static_cast<double>((alTEMin[lI] + lTotalTEFill) / 1000.0);
        fUICSetWithinLimits (pThis, dTEMinUI_ms, MR_TAG_TE, UIC_SET_WITHIN_LIMITS_ROUNDUP| UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, lI);

        // * -------------------------------------------------------------- *
        // * Consider additional TE fill times arising from the UI          *
        // * increment of TE.                                               *
        // * -------------------------------------------------------------- *
        lTotalTEFill += fSDSRoundUpGRT ((dTEMinUI_ms - dTEMinKernel_ms) * 1000);
    }



    return ( lDesiredContrasts );

}





/*[ Function ****************************************************************\
*
* Name        : segmentsSetValueHandler
*
* Description : This function calls the original UI method to set the number
*               of segments
*
* Return      : New number of segments
*
\****************************************************************************/
long segmentsSetValueHandler (LINK_LONG_TYPE* const pThis, long lDesiredSegments, long lIndex)
{

    // * ------------------------------------------------------------------------ *
    // * Set the number of segments using the original set function               *
    // * ------------------------------------------------------------------------ *
    if( pOriginalSegmentsSetFct )  {
        lDesiredSegments = ( *pOriginalSegmentsSetFct )(pThis, lDesiredSegments, lIndex);
    }

    return ( lDesiredSegments );

}








// * ------------------------------------------------------------------------ *
// * Declaration of solve handler functions                                   *
// * ------------------------------------------------------------------------ *

/*[ Function ****************************************************************\
*
* Name        : fSSolvePhaseStabilizeConflict
*
* Description : Solve handler for phase stabilization
*               TR will increasd if necessary
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolvePhaseStabilizeConflict
(
    LINK_BOOL_TYPE* const pThis,
    char**                arg_list,
    const void*,
    const MrProt*         pOrigProt,
    long
)
{
    return ( fUICSolveBoolParamConflict ( pThis, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );
}




/*[ Function ****************************************************************\
*
* Name        : fLSolveSegmentsConflict
*
* Description : Solve handler for segments
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fLSolveSegmentsConflict
(
    LINK_LONG_TYPE* const pThis,          // ...
    char**                arg_list,       // receives confirmation message
    const void*,                          // not needed
    const MrProt*         pOrigProt,      // Original protocol with old number of segments
    long                                  // Array index reserved
)
{

    return ( fUICSolveSegmentsConflict ( pThis, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );

}





/*[ Function ****************************************************************\
*
* Name        : fSSolveSelectionConflict



*
* Description : Solve handler for selection type parameters
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveSelectionConflict
(
    MrUILinkSelection<unsigned>* const _this,     // ...
    char**                             arg_list,  // receives confirmation message
    const void*,                                  // not needed
    const MrProt*                      pOrigProt, // Original protocol with old rf mode
    long                                          // Array index reserved
)
{
  return ( fUICSolveSelectionConflict ( _this, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );
}







/*[ Function ****************************************************************\
*
* Name        : fLSolveTEConflict
*
* Description : Solve handler for TE
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fLSolveTEConflict
(
    LINK_DOUBLE_TYPE* const pThis,          // ...
    char**                  arg_list,       // receives confirmation message
    const void*             pToAddMemory,   // not needed
    const MrProt*           pOrigProt,      // Original protocol with old number of segments
    long                    lIndex          // Array index reserved
)
{

    MrProt      *pMrProt        = (MrProt*)  &pThis->prot();
    SeqLim      *pSeqLim        = (SeqLim*)  &pThis->seqLimits();
    SeqExpo     *pSeqExpo       = &pThis->sequence().getSEQ_BUFF()->getSeqExpo();
    Sequence    *pSeq           = &pThis->sequence();

    long        alTEMin [lNoOfContrasts];
    long        alTEFill[lNoOfContrasts];
    long        lI              = 0;
    long        lTotalTEFill    = 0;
    bool        bAnErrorOccured = false;
    double      dTEMinKernel_ms = 0.0;
    double      dTEMinUI_ms     = 0.0;



    // * ---------------------------------------------------------------------- *
    // * Solve handler can fix TE problems occuring with multiple echoes        *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->contrasts() > 1 )  {


        // * -------------------------------------------------------------------*
        // * TE has been increased ==> increase TE times of higher contrasts if *
        // * necessary                                                          *
        // * -------------------------------------------------------------------*
        if ( pMrProt->te()[lIndex] > pOrigProt->te()[lIndex] )  {

            // * -------------------------------------------------------------- *
            // * Format the Confirmation Message                                *
            // * -------------------------------------------------------------- *
            if ( arg_list ) {
                pThis->formatConfirmationCause(arg_list, pOrigProt, lIndex);
            }


            // * -------------------------------------------------------------- *
            // * Update exports after protocol has been changed                 *
            // * -------------------------------------------------------------- *
            pSeq->prepareForBinarySearch( pMrProt );


            // * -------------------------------------------------------------- *
            // * Calculate the min. TE and TR for the selected base matrix      *
            // * -------------------------------------------------------------- *
            if ( SEQU__NORMAL !=  fCalculateTEMin (pMrProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill) ) bAnErrorOccured = true;



            // * -------------------------------------------------------------- *
            // * Calculate the accumulated TE fill times                        *
            // * -------------------------------------------------------------- *
            for ( lI = 0; lI <= lIndex; lI++ )  {
                lTotalTEFill += alTEFill[lI];
            }


            for ( lI = lIndex+1; lI < pMrProt->contrasts(); lI++ )  {

                // * ---------------------------------------------------------- *
                // * The solve handler has to increase TE when the minimum TE   *
                // * time from the kernel plus the accumulated TE fill times is *
                // * larger than the TE time in the original protocol.          *
                // * ---------------------------------------------------------- *
                if ( alTEMin[lI] + lTotalTEFill > pOrigProt->te()[lI] )  {
                    dTEMinKernel_ms = dTEMinUI_ms = static_cast<double>((alTEMin[lI] + lTotalTEFill) / 1000.0);
                    fUICSetWithinLimits (pThis, dTEMinUI_ms, MR_TAG_TE, UIC_SET_WITHIN_LIMITS_ROUNDUP| UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, lI);

                    // * ------------------------------------------------------ *
                    // * Consider additional TE fill times arising from the UI  *
                    // * increment of TE.                                       *
                    // * ------------------------------------------------------ *
                    lTotalTEFill += fSDSRoundUpGRT ((dTEMinUI_ms - dTEMinKernel_ms) * 1000);

                }

            }

            return ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) ? MRI_STD_CONFIRMATION_MSG : 0);

        }



        // * ------------------------------------------------------------------ *
        // * TE has been decreased ==> decrease TE times of lower contrasts if  *
        // * necessary                                                          *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->te()[lIndex] < pOrigProt->te()[lIndex] )  {

            long lTEDecrement = pOrigProt->te()[lIndex] - pMrProt->te()[lIndex];


            // * -------------------------------------------------------------- *
            // * Format the Confirmation Message                                *
            // * -------------------------------------------------------------- *
            if ( arg_list ) {
                pThis->formatConfirmationCause(arg_list, pOrigProt, lIndex);
            }


            // * -------------------------------------------------------------- *
            // * Update exports after protocol has been changed                 *
            // * -------------------------------------------------------------- *
            pSeq->prepareForBinarySearch( pMrProt );


            // * -------------------------------------------------------------- *
            // * Calculate the min. TE and TR for the selected base matrix      *
            // * -------------------------------------------------------------- *
            if ( SEQU__NORMAL !=  fCalculateTEMin ((MrProt*)pOrigProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill) ) bAnErrorOccured = true;



            // * -------------------------------------------------------------- *
            // * Decrease TE if necessary, starting from the end of the echo    *
            // * train                                                          *
            // * -------------------------------------------------------------- *
            lI = lIndex-1;

            while ( (lI >= 0) && (lTEDecrement > 0) )  {

                if ( lTEDecrement < alTEFill[lI+1] )  {

                    lTEDecrement =  0;

                } else {

                    // * ------------------------------------------------------ *
                    // * TE decrement is larger than TE fill of contrast (lI+1) *
                    // *    ==> TE of the previous contrast (lI) has to be      *
                    // *        adapted.                                        *
                    // * ------------------------------------------------------ *
                    lTEDecrement    -= alTEFill[lI+1];
                    dTEMinKernel_ms =  dTEMinUI_ms = static_cast<double>((pOrigProt->te()[lI] - lTEDecrement) / 1000.0);
                    fUICSetWithinLimits (pThis, dTEMinUI_ms, MR_TAG_TE, UIC_SET_WITHIN_LIMITS_ROUNDDOWN| UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, lI);


                    // * ------------------------------------------------------ *
                    // * Consider additional TE fill times arising from the UI  *
                    // * increment of TE.                                       *
                    // * ------------------------------------------------------ *
                    lTEDecrement += fSDSRoundUpGRT ((dTEMinKernel_ms - dTEMinUI_ms) * 1000);

                }

                lI--;

            }

            return ( pThis->tryProt((void*)pToAddMemory, pMrProt, lIndex) ? MRI_STD_CONFIRMATION_MSG : 0);

        }

    }


    // * ---------------------------------------------------------------------- *
    // * Can't solve the conflict                                               *
    // * ---------------------------------------------------------------------- *
    return ( 0 );

}





/*[ Function ****************************************************************\
*
* Name        : fLSolveSGSizeConflict
*
* Description : Solve handler for TE
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fLSolveSGSizeConflict
(
    LINK_LONG_TYPE* const pThis,          // ...
    char**                arg_list,       // receives confirmation message
    const void*           pToAddMemory,   // not needed
    const MrProt*         pOrigProt,      // Original protocol with old number of segments
    long                  lIndex          // Array index reserved
)
{

    MrProt      *pMrProt        = (MrProt*)  &pThis->prot();
    SeqLim      *pSeqLim        = (SeqLim*)  &pThis->seqLimits();


    // * ---------------------------------------------------------------------- *
    // * Check whether there is a quick SBBRSat                                 *
    // * ---------------------------------------------------------------------- *
    bool bQuickRSat = false;
    for (long lI=0; lI<pMrProt->satList().size(); lI++)  {
        if ( pMrProt->satList()[lI].quickMode() )  {
            bQuickRSat = true;
            break;
        }
    }



    // * ---------------------------------------------------------------------- *
    // * Binary search does not work properly in case of quick sat due to non   *
    // * convex parameter spaces. The minimum TR might not be found correctly,  *
    // * therefore this solve handler tests the minimum TR values explicitly.   *
    // * ---------------------------------------------------------------------- *
    if (    (pMrProt->sliceSeries().size() > pOrigProt->sliceSeries().size())
         && (    bQuickRSat || pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK
              || pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) )  {



        // * ------------------------------------------------------------------ *
        // * Estimated TR interval                                              *
        // * ------------------------------------------------------------------ *
        double  dTRMin_ms;
        double  dTRTemp_ms;
        long    lNewTR = pMrProt->tr()[0] + pSeqLim->getTR()[0].getInc()/2;
        long    lMaxTR = lNewTR + GREKernel.getlTRMin();


        // * ------------------------------------------------------------------ *
        // * Search interval between lNewTR and lMaxTR                          *
        // * ------------------------------------------------------------------ *
        while ( lNewTR <= lMaxTR )  {

            dTRTemp_ms = dTRMin_ms = static_cast<double>(lNewTR / 1000.0);     // * TR in ms *
            fUICSetWithinLimits ( pThis, dTRMin_ms, MR_TAG_TR, UIC_SET_WITHIN_LIMITS_ROUNDUP | UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, 0 );

            if ( pThis->tryProt((void*)pToAddMemory, pMrProt, lIndex) )  {
                return ( MRI_STD_CONFIRMATION_MSG );
            } else {
                if ( pSeqLim->getTR()[0].getInc() == 0 )  {
                    // Can't solve since TR increment is invalid
                    return ( 0 );
                }

                lNewTR = static_cast<long>(dTRTemp_ms * 1000.0 + pSeqLim->getTR()[0].getInc()/2);
            }

        }

        return ( (*pOriginalSGSizeSolveHandler)(pThis, arg_list, pToAddMemory, pOrigProt, lIndex) );

    } else {

        // * ------------------------------------------------------------------ *
        // * Call original solve handler in all other cases                     *
        // * ------------------------------------------------------------------ *
        return ( (*pOriginalSGSizeSolveHandler)(pThis, arg_list, pToAddMemory, pOrigProt, lIndex) );

    }

}






    /*[ Function ****************************************************************\
    *
    * Name        : fSSolveReadOutModeConflict
    *
    * Description : Solve handler for the readout mode
    *
    * Return      : confirmation message
    *
    \****************************************************************************/
    unsigned fSSolveReadOutModeConflict(LINK_SELECTION_TYPE* const pThis,
                                        char**                  arg_list,       // receives confirmation message
                                        const void*             pToAddMemory,   // not needed
                                        const MrProt*           pOrigProt,      // Original protocol with old number of segments
                                        long                    lIndex          // Array index reserved
    )
    {

        MrProt      *pMrProt        = (MrProt*)  &pThis->prot();
        SeqLim      *pSeqLim        = (SeqLim*)  &pThis->seqLimits();
        SeqExpo     *pSeqExpo       = &pThis->sequence().getSEQ_BUFF()->getSeqExpo();
        Sequence    *pSeq           = &pThis->sequence();

        long        alTEMin [lNoOfContrasts];
        long        alTEFill[lNoOfContrasts];
        long        lI              = 0;
        long        lTotalTEFill    = 0;
        bool        bAnErrorOccured = false;
        double      dTEMinKernel_ms = 0.0;
        double      dTEMinUI_ms     = 0.0;



        // * ------------------------------------------------------------------ *
        // * The user switches from bipolar to monopolar which might require    *
        // * longer TE times                                                    *
        // * ------------------------------------------------------------------ *
        if (    (pMrProt->contrasts() > 1)
             && (pOrigProt->readOutMode() == SEQ::READOUT_BIPOLAR  )
             && (pMrProt->readOutMode()   == SEQ::READOUT_MONOPOLAR) )  {


            // * -------------------------------------------------------------- *
            // * Format the Confirmation Message                                *
            // * -------------------------------------------------------------- *
            if ( arg_list ) {
                pThis->formatConfirmationCause(arg_list, pOrigProt, lIndex);
            }


            // * -------------------------------------------------------------- *
            // * Update exports after protocol has been changed                 *
            // * -------------------------------------------------------------- *
            pSeq->prepareForBinarySearch( pMrProt );


            // * -------------------------------------------------------------- *
            // * Calculate the min. TE and TR for the selected base matrix      *
            // * -------------------------------------------------------------- *
            if ( SEQU__NORMAL !=  fCalculateTEMin (pMrProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill) ) bAnErrorOccured = true;


            // * -------------------------------------------------------------- *
            // * Initialize lTotalTEFill                                        *
            // * -------------------------------------------------------------- *
            lTotalTEFill = alTEFill[0];


            for ( lI = 1; lI < pMrProt->contrasts(); lI++ )  {

                // * ---------------------------------------------------------- *
                // * The solve handler has to increase TE when the minimum TE   *
                // * time from the kernel plus the accumulated TE fill times is *
                // * larger than the TE time in the original protocol.          *
                // * ---------------------------------------------------------- *
                if ( alTEMin[lI] + lTotalTEFill > pOrigProt->te()[lI] )  {
                    dTEMinKernel_ms = dTEMinUI_ms = static_cast<double>((alTEMin[lI] + lTotalTEFill) / 1000.0);
                    fUICSetWithinLimits (pThis, dTEMinUI_ms, MR_TAG_TE, UIC_SET_WITHIN_LIMITS_ROUNDUP| UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, lI);

                    // * ------------------------------------------------------ *
                    // * Consider additional TE fill times arising from the UI  *
                    // * increment of TE.                                       *
                    // * ------------------------------------------------------ *
                    lTotalTEFill += fSDSRoundUpGRT ((dTEMinUI_ms - dTEMinKernel_ms) * 1000);

                }

            }

            return ( pThis->tryProt((void*)pToAddMemory, pMrProt, lIndex) ? MRI_STD_CONFIRMATION_MSG : 0);

        } else {

            // * -------------------------------------------------------------- *
            // * Can't solve the conflict                                       *
            // * -------------------------------------------------------------- *
            return ( 0 );

        }

    }






/*[ Function ****************************************************************\
*
* Name        : fSSolveBoolParamConflict
*
* Description : Solve handler for boolean type parameters
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveBoolParamConflict
(
    LINK_BOOL_TYPE* const              _this,     // ...
    char**                             arg_list,  // receives confirmation message
    const void*                        pAddMem,   // pointer to additional memory
    const MrProt*                      pOrigProt, // Original protocol with old rf mode
    long                               lPos       // Array index
)
{
  return ( fUICSolveBoolParamConflict ( _this, arg_list, pAddMem, pOrigProt, lPos, pfGetMinTE, NULL, NULL ) );
}




/*[ Function ****************************************************************\
*
* Name        : fSSolveLongConflict
*
* Description : Solve handler for long type parameters
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveLongConflict
(
    LINK_LONG_TYPE* const _this,     // ...
    char**                arg_list,  // receives confirmation message
    const void*,                     // not needed
    const MrProt*         pOrigProt, // Original protocol with old rf mode
    long                             // Array index reserved
)
{
  return ( fUICSolveLongParamConflict ( _this, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );
}




/*[ Function ****************************************************************\
*
* Name        : fSSolveInversionConflict
*
* Description : Solve handler for inversion pulse
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveInversionConflict
(
    MrUILinkSelection<unsigned>* const _this,     // ...
    char**                             arg_list,  // receives confirmation message
    const void*,                                  // not needed
    const MrProt*                      pOrigProt, // Original protocol with old rf mode
    long                                          // Array index reserved
)
{
  return ( fUICInversionConflict ( _this, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );
}





/*[ Function ****************************************************************\
*
* Name        : fSSolveDarkBloodConflict
*
* Description : Solve handler for darkblood preparation
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveDarkBloodConflict
(
    LINK_BOOL_TYPE* const       pThis,     // ...
    char**                      arg_list,  // receives confirmation message
    const void*,                           // not needed
    const MrProt*               pOrigProt, // Original protocol with old rf mode
    long                                   // Array index reserved
)
{
  return ( fUICDarkBloodConflict ( pThis, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );
}





/*[ Function ****************************************************************\
*
* Name        : fDSolveBaseResConflict
*
* Description : Solve handler for Base resolution conflicts
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fDSolveBaseResConflict
(
    LINK_LONG_TYPE* const pThis,          // ...
    char**                arg_list,       // receives confirmation message
    const void*           pToAddMemory,   // sometimes needed
    const MrProt*         pOrigProt,      // Original protocol
    long                  lIndex          // Array index reserved
)
{


    MrProt      *pMrProt        = (MrProt*)  &pThis->prot();
    SeqLim      *pSeqLim        = (SeqLim*)  &pThis->seqLimits();
    SeqExpo     *pSeqExpo       = &pThis->sequence().getSEQ_BUFF()->getSeqExpo();
    Sequence    *pSeq           = &pThis->sequence();

    if ( pOrigProt->contrasts() == 1 )  {

        return ( fUICSolveBaseResolutionConflict ( pThis, arg_list, pToAddMemory, pOrigProt, lIndex, pfGetMinTE, NULL, NULL, pOrgSolve_baseRes) );

    } else {

        long        alTEMin [lNoOfContrasts];
        long        alTEFill[lNoOfContrasts];
        long        lI              = 0;
        long        lTotalTEFill    = 0;
        bool        bAnErrorOccured = false;
        double      dTEMinKernel_ms = 0.0;
        double      dTEMinUI_ms     = 0.0;


        // * ------------------------------------------------------------------ *
        // * Format the Confirmation Message                                    *
        // * ------------------------------------------------------------------ *
        if ( arg_list ) {
            pThis->formatConfirmationCause (arg_list, pOrigProt, lIndex);
        }


        // * ------------------------------------------------------------------ *
        // * Update exports after protocol has been changed                     *
        // * ------------------------------------------------------------------ *
        pSeq->prepareForBinarySearch (pMrProt);


        // * ------------------------------------------------------------------ *
        // * Calculate the min. TE and TR for the selected base matrix          *
        // * ------------------------------------------------------------------ *
        if ( SEQU__NORMAL !=  fCalculateTEMin (pMrProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill) ) bAnErrorOccured = true;


        // * ------------------------------------------------------------------ *
        // * Calculate the accumulated TE fill times                            *
        // * ------------------------------------------------------------------ *
        for ( lI = 0; lI < lIndex; lI++ )  {
            lTotalTEFill += maximum (alTEFill[lI], 0L);
        }


        for ( lI = lIndex; lI < pMrProt->contrasts(); lI++ )  {

            // * -------------------------------------------------------------- *
            // * The solve handler has to increase TE when the minimum TE time  *
            // * from the kernel plus the accumulated TE fill times is larger   *
            // * than the TE time in the original protocol.                     *
            // * -------------------------------------------------------------- *
            if ( alTEMin[lI] + lTotalTEFill > pOrigProt->te()[lI] )  {
                dTEMinKernel_ms = dTEMinUI_ms = static_cast<double>((alTEMin[lI] + lTotalTEFill) / 1000.0);
                fUICSetWithinLimits (pThis, dTEMinUI_ms, MR_TAG_TE, UIC_SET_WITHIN_LIMITS_ROUNDUP| UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, lI);


                // * ---------------------------------------------------------- *
                // * Consider additional TE fill times arising from the UI      *
                // * increment of TE.                                           *
                // * ---------------------------------------------------------- *
                lTotalTEFill += fSDSRoundUpGRT ((dTEMinUI_ms - dTEMinKernel_ms) * 1000);

            }

            lTotalTEFill += maximum (alTEFill[lI], 0L);


        }

        return ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) ? MRI_STD_CONFIRMATION_MSG : 0);

    }

}







/*[ Function ****************************************************************\
*
* Name        : fSolveDoubleMultiTEConflict
*
* Description : Solve handler that changes TE for multiple contrasts
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSolveDoubleMultiTEConflict
(
    LINK_DOUBLE_TYPE* const pThis,     // ...
    char**                             arg_list,  // receives confirmation message
    const void*                        pToAddMemory, // not needed
    const MrProt*                      pOrigProt, // Original protocol with old rf mode
    long                               lIndex     // Array index reserved
)
{

    MrProt      *pMrProt        = (MrProt*)  &pThis->prot();
    SeqLim      *pSeqLim        = (SeqLim*)  &pThis->seqLimits();
    SeqExpo     *pSeqExpo       = &pThis->sequence().getSEQ_BUFF()->getSeqExpo();
    Sequence    *pSeq           = &pThis->sequence();

    long        alTEMin [lNoOfContrasts];
    long        alTEFill[lNoOfContrasts];
    long        lI              = 0;
    long        lTotalTEFill    = 0;
    bool        bAnErrorOccured = false;
    double      dTEMinKernel_ms = 0.0;
    double      dTEMinUI_ms     = 0.0;


    // * ---------------------------------------------------------------------- *
    // * Format the Confirmation Message                                        *
    // * ---------------------------------------------------------------------- *
    if ( arg_list ) {
        pThis->formatConfirmationCause (arg_list, pOrigProt, lIndex);
    }


    // * ---------------------------------------------------------------------- *
    // * Update exports after protocol has been changed                         *
    // * ---------------------------------------------------------------------- *
    pSeq->prepareForBinarySearch (pMrProt);


    // * ---------------------------------------------------------------------- *
    // * Calculate the min. TE and TR for the selected base matrix              *
    // * ---------------------------------------------------------------------- *
    if ( SEQU__NORMAL !=  fCalculateTEMin (pMrProt, pSeqLim, pSeqExpo, pSeq, alTEMin, alTEFill) ) bAnErrorOccured = true;


    // * ---------------------------------------------------------------------- *
    // * Calculate the accumulated TE fill times                                *
    // * ---------------------------------------------------------------------- *
    for ( lI = 0; lI < lIndex; lI++ )  {
        lTotalTEFill += maximum (alTEFill[lI], 0L);
    }


    for ( lI = lIndex; lI < pMrProt->contrasts(); lI++ )  {

        // * ---------------------------------------------------------- *
        // * The solve handler has to increase TE when the minimum TE   *
        // * time from the kernel plus the accumulated TE fill times is *
        // * larger than the TE time in the original protocol.          *
        // * ---------------------------------------------------------- *
        if ( alTEMin[lI] + lTotalTEFill > pOrigProt->te()[lI] )  {
            dTEMinKernel_ms = dTEMinUI_ms = static_cast<double>((alTEMin[lI] + lTotalTEFill) / 1000.0);
            fUICSetWithinLimits (pThis, dTEMinUI_ms, MR_TAG_TE, UIC_SET_WITHIN_LIMITS_ROUNDUP| UIC_SET_WITHIN_LIMITS_LIMITS_ARE_CONST, lI);


            // * ------------------------------------------------------ *
            // * Consider additional TE fill times arising from the UI  *
            // * increment of TE.                                       *
            // * ------------------------------------------------------ *
            lTotalTEFill += fSDSRoundUpGRT ((dTEMinUI_ms - dTEMinKernel_ms) * 1000);

        }

        lTotalTEFill += maximum (alTEFill[lI], 0L);


    }

    return ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) ? MRI_STD_CONFIRMATION_MSG : 0);

}






/*[ Function ****************************************************************\
*
* Name        : fDSolveSliceOSConflict
*
* Description : Solve handler for slice OS
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fDSolveSliceOSConflict
(
    LINK_DOUBLE_TYPE* const pThis,     // ...
    char**                             arg_list,  // receives confirmation message
    const void*,                                  // not needed
    const MrProt*                      pOrigProt, // Original protocol with old rf mode
    long                               lIndex     // Array index reserved
)
{
  return ( fUICSolveDoubleParamConflict ( pThis, arg_list, NULL, pOrigProt, lIndex, pfGetMinTE, NULL, NULL ) );
}




/*[ Function ****************************************************************\
*
* Name        : fDSolveFlipAngle
*
* Description : Solve handler for the flip angle
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fDSolveFlipAngle
(
    LINK_DOUBLE_TYPE* const pThis,     // ...
    char**                             ,        // receives confirmation message
    const void*                        pToAddMemory,    // not needed
    const MrProt*                      pOrigProt,       // Original protocol with old rf mode
    long                               lIndex           // Array index reserved
)
{


    // * -------------------------------------------------------------------------- *
    // * Solve handler is intended for 3D with water excitation                     *
    // * -------------------------------------------------------------------------- *
    if (    (pThis->prot().kSpace().dimension() == SEQ::DIM_3)
         && ((pThis->prot().preparationPulses().fatSuppression() == SEQ::WATER_EXCITATION) || (pThis->prot().preparationPulses().fatSuppression() == SEQ::WATER_EXCITATION_FAST) ) )  {

        long        lNoOfIter     = 0;             // Number of iterations
        const long  lMaxNoOfIter  = 10;

        double      dActualSliceThickness    = 0.0;
        double      dIncreasedSlcieThickness = 0.0;

        LINK_DOUBLE_TYPE* pSliceThickness = _search<LINK_DOUBLE_TYPE>(pThis, MR_TAG_SLICE_THICKNESS);

        // * ---------------------------------------------------------------------- *
        // * Determine and check the number of the images per slab                  *
        // * ---------------------------------------------------------------------- *
        const long  lImagesPerSlab = pThis->prot().kSpace().imagesPerSlab();
        if ( lImagesPerSlab == 0 )  {  return (0L);  }


        // * ---------------------------------------------------------------------- *
        // * Iteratively increase the slice thickness                               *
        // * ---------------------------------------------------------------------- *
        while ( (!pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex)) && (lNoOfIter < lMaxNoOfIter) )  {

            // * ------------------------------------------------------------------ *
            // * Calculate the actual slice thickness                               *
            // * ------------------------------------------------------------------ *
            dActualSliceThickness = pThis->prot().sliceSeries()[0].thickness() / lImagesPerSlab;


            // * ------------------------------------------------------------------ *
            // * Calculate the next slice thickness to be tested.                   *
            // * For better performance use 10 times the 3DPartitionThickness       *
            // * increment.                                                         *
            // * ------------------------------------------------------------------ *
            dIncreasedSlcieThickness = dActualSliceThickness + pThis->seqLimits().get3DPartThickness().getInc() * 10;
            if (pSliceThickness->value (dIncreasedSlcieThickness, 0) != dIncreasedSlcieThickness)  {  return ( 0 );  }


            // * ------------------------------------------------------------------ *
            // * Increment iteration counter                                        *
            // * ------------------------------------------------------------------ *
            lNoOfIter++;
        }


        // * ---------------------------------------------------------------------- *
        // * Check whether conflict has been solved by increasing the slice         *
        // * thickness.                                                             *
        // * ---------------------------------------------------------------------- *
        if( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {
            return ( MRI_STD_CONFIRMATION_MSG );
        } else {
            return ( 0L );
        }

    } else {

        // Sorry, can't solve the problem !!
        return ( 0L );

    }

}






/*[ Function ****************************************************************\
*
* Name        : fSSolveExcitationModeConflict
*
* Description : Solve handler for excitation mode
*
*               Non-selective excitation is only allowed for 3D imaging
*               The solve handler switches to 3D if the user selects
*               non-selective excitation in combination with 2D imaging,
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveExcitationModeConflict
(
    LINK_SELECTION_TYPE*  const pThis,     // ...
    char**                arg_list,        // receives confirmation message
    const void*           pToAddMemory,
    const MrProt*         pOrigProt, // Original protocol with old rf mode
    long                  lIndex     // Array index reserved
)
{

    MrProt*               pNewProt   = &pThis->prot();               // with new PhasePartialFourier value
    long                  lCtrl      = 0;


    // * -------------------------------------------------------------------- *
    // * The user has selected non-selective excitation, but there is still   *
    // * 2D imaging selected.                                                 *
    // * -------------------------------------------------------------------- *
    if ( (pNewProt->txSpec().excitation() == SEQ::EXCITATION_VOLUME_SELECTIVE) &&
         (pNewProt->kSpace().dimension() == SEQ::DIM_2) )                            {


        if ( arg_list ) {

            // * ------------------------------------------------------------------ *
            // * Format the Confirmation Message                                    *
            // *                                                                    *
            // * Your last change                                                   *
            // *                                                                    *
            // * \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!                              *
            // * ------------------------------------------------------------------ *
            arg_list[0]  = (char*) pThis->labelId(arg_list+1, lIndex);  // resurce ID -> label
            arg_list[10] = (char*) pThis->value(lIndex);                // new value
            arg_list[16] = (char*) pThis->unitId(arg_list+17, lIndex);

            MrProt _storage;
            void* pVoid     = pThis->storeInMemory(&_storage,lIndex);
            //    pThis->prot()   = *pOrigProt;
            pThis->recallMemory(pToAddMemory,pOrigProt,lIndex);

            arg_list[4]     = (char*) pThis->value(lIndex);            // old value

            pThis->recallMemory (pVoid,&_storage,lIndex);
            pThis->clearMemory  (pVoid,lIndex);



            // * -------------------------------------------------------------------- *
            // * has adapted the following parameters:                                *
            // *                                                                      *
            // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                        *
            // * -------------------------------------------------------------------- *
            arg_list[20] = (char*) MRI_STD_DIMENSION_LABEL;

            switch ( pOrigProt->kSpace().dimension() ) {
                case SEQ::DIM_1:    arg_list[24] = (char*) MRI_STD_DIMENSION_1D;     break;
                case SEQ::DIM_2:    arg_list[24] = (char*) MRI_STD_DIMENSION_2D;     break;
            }

            arg_list[30] = (char*) MRI_STD_DIMENSION_3D;
            arg_list[36] = (char*) MRI_STD_EMPTY;
            arg_list[39] = (char*) '\n';
            lCtrl++;

        }

        // * -------------------------------------------------------------------- *
        // * Action:                                                              *
        // *                                                                      *
        // * now we switch to 3D imaging                                          *
        // * -------------------------------------------------------------------- *

        LINK_SELECTION_TYPE* pDimensionMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIMENSION);

        // * Check whether sequential slice mode can be selected *
        if(!pDimensionMode || !pDimensionMode->isEditable(0) || !pThis->seqLimits().getDimension().hasOption(SEQ::DIM_3))  {  return ( 0 );  }

        // * Switch to 3D *
        if(pDimensionMode->value(MRI_STD_DIMENSION_3D,0) != MRI_STD_DIMENSION_3D)  {  return ( 0 );  }


        //  calls prepareForBinarySearch by default
        if( !pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

            return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

        } else {

            if(arg_list) {  arg_list[39] = (char*) '\0';  }
            return ( MRI_STD_CONFIRMATION_MSG );

        }


    } else {

        // Sorry, can't solve the problem !!
        return ( false );

    }

}





/*[ Function ****************************************************************\
*
* Name        : fSSolveFirstSignalModeConflict
*
* Description : Solve handler for the first physio signal mode
*
*               Respiratory triggering is not allowed in combination with
*               sequential slice
*               The solve handler switches to interleaved multi slice if
*               respiratory triggering has been selected by the user
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveFirstSignalModeConflict
(
    LINK_SELECTION_TYPE*  const pThis,     // ...
    char**                arg_list,        // receives confirmation message
    const void*           pToAddMemory,
    const MrProt*         pOrigProt, // Original protocol with old rf mode
    long                  lIndex     // Array index reserved
)
{

    MrProt*               pNewProt   = &pThis->prot();               // with new PhasePartialFourier value

    const long            lOffset    = 20;
    long                  lCtrl      = 0;

    const long lMinNoOfSegments = 1;
    long       lNewNoOfSegments = 1;

    // * ------------------------------------------------------------------------- *
    // * Determine trigger mode and signal                                         *
    // * ------------------------------------------------------------------------- *
    SEQ::PhysioSignal PhysioSignalHigh          = SEQ::SIGNAL_NONE;   // * Selected trigger modes    *
    SEQ::PhysioSignal PhysioSignalLow           = SEQ::SIGNAL_NONE;   // * and signals               *
    SEQ::PhysioMethod PhysioMethodHigh          = SEQ::METHOD_NONE;
    SEQ::PhysioMethod PhysioMethodLow           = SEQ::METHOD_NONE;

    pNewProt->physiology().getPhysioMode (PhysioSignalHigh, PhysioMethodHigh, PhysioSignalLow, PhysioMethodLow);



    // * ------------------------------------------------------------------------- *
    // * The user has selected respiratory triggering but sequetial multi slice    *
    // * mode is active.                                                           *
    // * ------------------------------------------------------------------------- *
    if ( ( PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION ) && ( PhysioMethodHigh == SEQ::METHOD_TRIGGERING )  &&
         ( pNewProt->kSpace().multiSliceMode() == SEQ::MSM_SEQUENTIAL ) ) {

        if ( arg_list ) {

            // * ------------------------------------------------------------------ *
            // * Format the Confirmation Message                                    *
            // *                                                                    *
            // * Your last change                                                   *
            // *                                                                    *
            // * \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!                              *
            // * ------------------------------------------------------------------ *
            arg_list[0]  = (char*) pThis->labelId(arg_list+1, lIndex);  // resurce ID -> label
            arg_list[10] = (char*) pThis->value(lIndex);                // new value
            arg_list[16] = (char*) pThis->unitId(arg_list+17, lIndex);

            MrProt _storage;
            void* pVoid     = pThis->storeInMemory(&_storage,lIndex);
            //    pThis->prot()   = *pOrigProt;
            pThis->recallMemory(pToAddMemory,pOrigProt,lIndex);

            arg_list[4]     = (char*) pThis->value(lIndex);            // old value

            pThis->recallMemory (pVoid,&_storage,lIndex);
            pThis->clearMemory  (pVoid,lIndex);



            // * -------------------------------------------------------------------- *
            // * has adapted the following parameters:                                *
            // *                                                                      *
            // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                        *
            // * -------------------------------------------------------------------- *
            arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_MSM_LABEL;

            switch ( pOrigProt->kSpace().multiSliceMode() ) {
                case SEQ::MSM_SEQUENTIAL:    arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_MSM_SEQUENTIAL;     break;
                case SEQ::MSM_SINGLESHOT:    arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_MSM_SINGLE_SHOT;    break;
            }

            arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_MSM_INTERLEAVED;
            arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
            arg_list[39 + lCtrl * lOffset] = (char*) '\n';
            lCtrl++;

        }

        // * -------------------------------------------------------------------- *
        // * Action:                                                              *
        // *                                                                      *
        // * now we switch to interleaved multi slice                             *
        // * -------------------------------------------------------------------- *

        LINK_SELECTION_TYPE* pMultiSliceMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_MULTI_SLICE_MODE);

        // * Check whether sequential slice mode can be selected *
        if(!pMultiSliceMode || !pMultiSliceMode->isEditable(0) || !pThis->seqLimits().getMultiSliceMode().hasOption(SEQ::MSM_INTERLEAVED))  {  return ( 0 );  }

        // * Switch to 3D *
        if(pMultiSliceMode->value(MRI_STD_MSM_INTERLEAVED,0) != MRI_STD_MSM_INTERLEAVED)  {  return ( 0 );  }




        // * ---------------------------------------------------------------------- *
        // * If interleaved multi slice mode has been enabled, inversion recovery   *
        // * has to be switched off.                                                *
        // * ---------------------------------------------------------------------- *
        if ( (   pOrigProt->preparationPulses().inversion() == SEQ::SLICE_SELECTIVE
              || pOrigProt->preparationPulses().inversion() == SEQ::VOLUME_SELECTIVE )
             && !switchOffInversionRecovery(pThis, arg_list, pOrigProt, lCtrl)  ) {
            return ( 0 );
        }



        // * ---------------------------------------------------------------------- *
        // * If interleaved multi slice mode has been enabled, saturation recovery  *
        // * has to be switched off.                                                *
        // * ---------------------------------------------------------------------- *
        if (    (   pOrigProt->preparationPulses().satRecovery() == SEQ::SATREC_SLICE_SELECTIVE
                 || pOrigProt->preparationPulses().satRecovery() == SEQ::SATREC_VOLUME_SELECTIVE )
             && ( !switchOffSaturationRecovery (pThis, arg_list, pOrigProt, lCtrl) ) )    {
            return ( 0 );
        }



        // * ------------------------------------------------------------------ *
        // * If interleaved multi slice mode has been enabled, dark blood       *
        // * preparation has to be switched off.                                *
        // * ------------------------------------------------------------------ *
        if ( pOrigProt->preparationPulses().darkBlood() && !switchOffDarkBlood (pThis, arg_list, pOrigProt, lCtrl) )  {
            return ( 0 );
        }



        // * ------------------------------------------------------------------ *
        // * Check whether the conflict has been solved                         *
        // * ------------------------------------------------------------------ *
        if ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

            if(arg_list) {  arg_list[39 + (lCtrl-1) * lOffset] = (char*) '\0';  }
            return ( MRI_STD_CONFIRMATION_MSG );

        } else {


            // * -------------------------------------------------------------- *
            // * Test whether there is a problem with the number of             *
            // * concatenations                                                 *
            // * Increment the number of concatenations in order to see whether *
            // * this solves the conflict.                                      *
            // * -------------------------------------------------------------- *
            long lNewNoOfConcats = pNewProt->concatenations() + 1;

            if (! fLocalUICSetWithinLimits (pThis, lNewNoOfConcats, MR_TAG_CONCATENATIONS, UIC_SET_WITHIN_LIMITS_ROUNDUP, -1) ) {  return ( 0 );  };


            if ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

                if ( arg_list && (pOrigProt->concatenations() != lNewNoOfConcats) ) {
                    arg_list[20 + lCtrl * lOffset] = (char*) MR_TAG_CONCATENATIONS;
                    arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_INT;
                    arg_list[25 + lCtrl * lOffset] = (char*) pOrigProt->concatenations();
                    arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_INT;
                    arg_list[31 + lCtrl * lOffset] = (char*) lNewNoOfConcats;
                    arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
                    arg_list[39 + lCtrl * lOffset] = (char*) '\0';
                }

                return ( MRI_STD_CONFIRMATION_MSG );

            } else {

                // * ---------------------------------------------------------- *
                // * Increment the number of concatenations by two in order to  *
                // * see whether this solves the conflict.                      *
                // * A larger increment in concatenations is not cosidered as a *
                // * solution due to measurement time constraints.              *
                // * ---------------------------------------------------------- *
                long lNewNoOfConcats = pNewProt->concatenations() + 1;

                if (! fLocalUICSetWithinLimits (pThis, lNewNoOfConcats, MR_TAG_CONCATENATIONS, UIC_SET_WITHIN_LIMITS_ROUNDUP, -1) ) {  return ( 0 );  };


                if ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

                    if ( arg_list && (pOrigProt->concatenations() != lNewNoOfConcats) ) {
                        arg_list[20 + lCtrl * lOffset] = (char*) MR_TAG_CONCATENATIONS;
                        arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_INT;
                        arg_list[25 + lCtrl * lOffset] = (char*) pOrigProt->concatenations();
                        arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_INT;
                        arg_list[31 + lCtrl * lOffset] = (char*) lNewNoOfConcats;
                        arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
                        arg_list[39 + lCtrl * lOffset] = (char*) '\0';
                    }

                    return ( MRI_STD_CONFIRMATION_MSG );

                } else {

                    //  calls prepareForBinarySearch by default
                    if ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {


                        if(arg_list) {  arg_list[39 + (lCtrl-1) * lOffset] = (char*) '\0';  }
                        return ( MRI_STD_CONFIRMATION_MSG );


                    } else {

                        if ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) )  {

                            if(arg_list) {  arg_list[39 + (lCtrl-1) * lOffset] = (char*) '\0';  }
                            return ( MRI_STD_CONFIRMATION_MSG );

                        } else {

                            if ( pNewProt->tr()[0] * pNewProt->fastImaging().segments() >  pNewProt->physiology().scanWindow(PhysioSignalHigh) * 1000 ) {

                                if ( pNewProt->tr()[0] > 0 ) {
                                    lNewNoOfSegments = pNewProt->physiology().scanWindow(PhysioSignalHigh) * 1000 / pNewProt->tr()[0];
                                } else {
                                    return ( 0 );
                                }

                                if ( arg_list ) {
                                    arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_SEGMENTS_LABEL;
                                    arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_INT;
                                    arg_list[25 + lCtrl * lOffset] = (char*) pOrigProt->fastImaging().segments();
                                    arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_INT;
                                    arg_list[31 + lCtrl * lOffset] = (char*) lNewNoOfSegments;
                                    arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
                                    arg_list[39 + lCtrl * lOffset] = (char*) '\n';
                                    lCtrl++;
                                }


                                LINK_LONG_TYPE* pSegments = _search< LINK_LONG_TYPE > ( pThis, MR_TAG_SEGMENTS );

                                // * Check whether phases dialog can be selected *
                                if(!pSegments || !pSegments->isEditable(0))  {  return ( 0 );  }

                                // * Set the new number of segments *
                                long lProposedNewNoOfSegments = lNewNoOfSegments;
                                if (! fLocalUICSetWithinLimits (pThis, lNewNoOfSegments, MR_TAG_SEGMENTS, UIC_SET_WITHIN_LIMITS_ROUNDDOWN, -1) ) {  return ( 0 );  };

                                // * The number of segments that has been set by fUICSetWithinLimits is *
                                // * larger than the proposed number of segments                        *
                                // * ==> TR will probably not fit into the acquisition window, no       *
                                // *     solution has been found !!!                                    *
                                if ( lNewNoOfSegments > lProposedNewNoOfSegments ) {  return ( 0 );  }



                                //  calls prepareForBinarySearch by default
                                if( !pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

                                    return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

                                } else {

                                    if(arg_list) {  arg_list[39 + (lCtrl-1) * lOffset] = (char*) '\0';  }
                                    return ( MRI_STD_CONFIRMATION_MSG );

                                }

                            } else {

                                // * ------------------------------------------------------------- *
                                // * The problem can not be fixed, by changing the number of       *
                                // * segments ==> the solve handler failed                         *
                                // * ------------------------------------------------------------- *
                                return ( 0 );

                            }

                        }

                    }

                }

            }

        }


    // * ------------------------------------------------------------------------- *
    // * The user has deselected respiratory triggering but interleaved multi      *
    // * slice mode active.                                                        *
    // * ------------------------------------------------------------------------- *
    }  else if ( ( PhysioSignalHigh != SEQ::SIGNAL_RESPIRATION ) || ( PhysioMethodHigh != SEQ::METHOD_TRIGGERING ) &&
                ( pNewProt->kSpace().multiSliceMode() == SEQ::MSM_INTERLEAVED ) ) {

        if ( arg_list ) {

            // * ------------------------------------------------------------------ *
            // * Format the Confirmation Message                                    *
            // *                                                                    *
            // * Your last change                                                   *
            // *                                                                    *
            // * \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!                              *
            // * ------------------------------------------------------------------ *
            arg_list[0]  = (char*) pThis->labelId(arg_list+1, lIndex);  // resurce ID -> label
            arg_list[10] = (char*) pThis->value(lIndex);                // new value
            arg_list[16] = (char*) pThis->unitId(arg_list+17, lIndex);

            MrProt _storage;
            void* pVoid     = pThis->storeInMemory(&_storage,lIndex);
            //    pThis->prot()   = *pOrigProt;
            pThis->recallMemory(pToAddMemory,pOrigProt,lIndex);

            arg_list[4]     = (char*) pThis->value(lIndex);            // old value

            pThis->recallMemory (pVoid,&_storage,lIndex);
            pThis->clearMemory  (pVoid,lIndex);



            // * -------------------------------------------------------------------- *
            // * has adapted the following parameters:                                *
            // *                                                                      *
            // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                        *
            // * -------------------------------------------------------------------- *
            if ( PhysioMethodHigh != SEQ::METHOD_NONE ) {
                arg_list[20] = (char*) MRI_STD_MSM_LABEL;
                arg_list[24] = (char*) MRI_STD_MSM_INTERLEAVED;
                arg_list[30] = (char*) MRI_STD_MSM_SEQUENTIAL;
            } else {
                arg_list[20] = (char*) MRI_STD_SEGMENTS_LABEL;
                arg_list[24] = (char*) MRI_STD_INT;
                arg_list[25] = (char*) pOrigProt->fastImaging().segments();
                arg_list[30] = (char*) MRI_STD_INT;
                arg_list[31] = (char*) lMinNoOfSegments;
            }
            arg_list[36] = (char*) MRI_STD_EMPTY;
            arg_list[39] = (char*) '\n';
            lCtrl++;

        }

        // * -------------------------------------------------------------------- *
        // * The behaviour of the solve handler depends on the selection of the   *
        // * user:                                                                *
        // *                                                                      *
        // * - If the user selects a different type of phyiologic triggering or   *
        // *   gating the multi slice mode will be switched from interleaved to   *
        // *   sequential multi slice.                                            *
        // *                                                                      *
        // * - If the user switches phyio mode to none, the interleaved multi     *
        // *   slice will not be changed, but the number of segments will be      *
        // *   changed to one.                                                    *
        // * -------------------------------------------------------------------- *

        if ( PhysioMethodHigh != SEQ::METHOD_NONE ) {

            // * switch interleaved multi slice mode to sequential *
            LINK_SELECTION_TYPE* pMultiSliceMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_MULTI_SLICE_MODE);

            // * Check whether sequential slice mode can be selected *
            if(!pMultiSliceMode || !pMultiSliceMode->isEditable(0) || !pThis->seqLimits().getMultiSliceMode().hasOption(SEQ::MSM_SEQUENTIAL))  {  return ( 0 );  }

            // * Switch to sequentiell multi slice mode *
            if(pMultiSliceMode->value(MRI_STD_MSM_SEQUENTIAL,0) != MRI_STD_MSM_SEQUENTIAL)  {  return ( 0 );  }


            // * ------------------------------------------------------------------ *
            // * It may be necessary to adjust the acquisition window length after  *
            // * changing the physio mode.                                          *
            // *    - Select the maximum acquisition window length. This is         *
            // *      necessary so that the solve handler fUICStandardSolveHandler  *
            // *      can successfully adjust TE, TI and TR                         *
            // *    - Adjust TE, TI and TR using fUICStandardSolveHandler           *
            // *    - decrease the acquisition window length to the minimum value   *
            // * ------------------------------------------------------------------ *

            // * Previous value of the acquisition window length *
            long lUIAcqWindow_ms = pThis->prot().physiology().scanWindow(PhysioSignalHigh);


            // * Set the acquisition window to its maximum value *
            pThis->prot().physiology().scanWindow(PhysioSignalHigh, 10000);


            // * Adjust TE, TI and TR *
            fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl );


            // * Set the acquisition window to its minimum value *
            long lMinAcqWindow_ms =  pThis->prot().tr()[0] / 1000;
            lMinAcqWindow_ms      *= pThis->prot().physiology().phases();
            lMinAcqWindow_ms      *= pThis->sequence().getSEQ_BUFF()->getSeqExpo().getAcqWindowFactor();
            lMinAcqWindow_ms      += pThis->prot().physiology().triggerDelay( PhysioSignalHigh ) / 1000;


            pThis->prot().physiology().scanWindow(PhysioSignalHigh, maximum(lMinAcqWindow_ms, lUIAcqWindow_ms) );


            // * Check whether the solve handler has been successful and return status *
            return ( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) ? MRI_STD_CONFIRMATION_MSG : 0 );

        } else {

            // * reduce the number of segments to one *
            LINK_LONG_TYPE* pSegments = _search<LINK_LONG_TYPE>( pThis, MR_TAG_SEGMENTS );

            // * Check whether segmented can be selected *
            if( !pSegments || !pSegments->isEditable(0) || pThis->seqLimits().getSegments().getMin() > lMinNoOfSegments )   {  return ( 0 );  }

            // * Reduce the number of segments to one *
            if( pSegments->value ( lMinNoOfSegments, 0 ) != lMinNoOfSegments )  {  return ( 0 );  }


            //  calls prepareForBinarySearch by default
            if( !pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

                return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

            } else {

                if(arg_list) {  arg_list[39 + (lCtrl-1) * lOffset] = (char*) '\0';  }
                return ( MRI_STD_CONFIRMATION_MSG );

            }

        }


    } else {

        // Sorry, can't solve the problem !!
        return ( 0 );

    }

}





/*[ Function ****************************************************************\
*
* Name        : fSSolveFlowCompConfilct
*
* Description : Solve handler for the flow compensation mode
*
*               Increases TE and/or TR if necessary
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveFlowCompConfilct
(
    LINK_SELECTION_TYPE*  const pThis,     // ...
    char**                arg_list,        // receives confirmation message
    const void*           pToAddMemory,
    const MrProt*         pOrigProt, // Original protocol with old rf mode
    long                  lIndex     // Array index reserved
)
{


    MrProt*               pNewProt   = &pThis->prot();               // with new PhasePartialFourier value

    const long            lOffset    = 20;
    long                  lCtrl      = 0;



    // * -------------------------------------------------------------------------- *
    // * The user has selected flow compensation in slice select direction          *
    // * There might be two different problems:                                     *
    // *    - combination with 3D imaging.                                          *
    // *    - TE and/or TR are to small                                             *
    // * -------------------------------------------------------------------------- *
    if (    pNewProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_YES
         || pNewProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_SLICESEL_ONLY )  {


        if ( arg_list ) {

            // * ------------------------------------------------------------------ *
            // * Format the Confirmation Message                                    *
            // *                                                                    *
            // * Your last change                                                   *
            // *                                                                    *
            // * \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!                              *
            // * ------------------------------------------------------------------ *
            arg_list[0]  = (char*) pThis->labelId(arg_list+1, lIndex);  // resurce ID -> label
            arg_list[10] = (char*) pThis->value(lIndex);                // new value
            arg_list[16] = (char*) pThis->unitId(arg_list+17, lIndex);

            MrProt _storage;
            void* pVoid     = pThis->storeInMemory(&_storage,lIndex);
            //    pThis->prot()   = *pOrigProt;
            pThis->recallMemory(pToAddMemory,pOrigProt,lIndex);

            arg_list[4]     = (char*) pThis->value(lIndex);            // old value

            pThis->recallMemory (pVoid,&_storage,lIndex);
            pThis->clearMemory  (pVoid,lIndex);

        }


        if ( pNewProt->kSpace().dimension() == SEQ::DIM_3 )  {

            // * ------------------------------------------------------------------ *
            // * Check whether switching to 2D solves the conflict                  *
            // * ------------------------------------------------------------------ *

            if ( arg_list ) {

                // * -------------------------------------------------------------- *
                // * has adapted the following parameters:                          *
                // *                                                                *
                // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                  *
                // * -------------------------------------------------------------- *
                arg_list[20 + lCtrl * lOffset] = (char*) MRI_STD_DIMENSION_LABEL;
                arg_list[24 + lCtrl * lOffset] = (char*) MRI_STD_DIMENSION_3D;
                arg_list[30 + lCtrl * lOffset] = (char*) MRI_STD_DIMENSION_2D;
                arg_list[36 + lCtrl * lOffset] = (char*) MRI_STD_EMPTY;
                arg_list[39 + lCtrl * lOffset] = (char*) '\n';
                lCtrl++;

            }


            // * -------------------------------------------------------------------- *
            // * Action:                                                              *
            // *                                                                      *
            // * now we switch to 2D imaging                                          *
            // * -------------------------------------------------------------------- *
            LINK_SELECTION_TYPE* pDimensionMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIMENSION);

            // * Check whether rf spoiling can be selected *
            if( !pDimensionMode || !pDimensionMode->isEditable(0) || !pThis->seqLimits().getDimension().hasOption(SEQ::DIM_2) ) {  return ( 0 );  }

            // * Switch rf spoiling on *
            if(pDimensionMode->value( MRI_STD_DIMENSION_2D, 0 ) != MRI_STD_DIMENSION_2D)  {  return ( 0 );  }



            //  calls prepareForBinarySearch by default
            if( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

                // * -------------------------------------------------------------- *
                // * Conflict has been solved                                       *
                // * -------------------------------------------------------------- *
                if(arg_list) {  arg_list[39] = (char*) '\0';  }
                return ( MRI_STD_CONFIRMATION_MSG );

            } else {

                // * -------------------------------------------------------------- *
                // * Maybe that TE and/or TR have to be increased, additionally     *
                // * -------------------------------------------------------------- *
                return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

            }

        } else {

            // * ------------------------------------------------------------------ *
            // * Increase TE and/or TR if necessary                                 *
            // * ------------------------------------------------------------------ *
            return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

        }

    } else {

        return ( fUICSolveSelectionConflict ( pThis, arg_list, NULL, pOrigProt, 0, pfGetMinTE, NULL, NULL ) );

    }

}






/*[ Function ****************************************************************\
*
* Name        : fSSolveDimensionConflict
*
* Description : Solve handler for the dimension
*
*               Increases TE and/or TR if necessary
*               Switch off flow compensation in slice directon
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveDimensionConflict
(
    LINK_SELECTION_TYPE*  const pThis,     // ...
    char**                arg_list,        // receives confirmation message
    const void*           pToAddMemory,
    const MrProt*         pOrigProt, // Original protocol with old rf mode
    long                  lIndex     // Array index reserved
)
{

    MrProt*     pNewProt    = &pThis->prot();   // with new dimension value
    const long  lOffset     = 20;
    long        lCtrl       = 0;


    // * -------------------------------------------------------------------------- *
    // * Call original solve handler                                                *
    // * -------------------------------------------------------------------------- *
    unsigned nOrgSolveRet = (pOriginalDimensionSolveHandler ? (*pOriginalDimensionSolveHandler)(pThis, arg_list, pToAddMemory, pOrigProt, lIndex) : 0);
    if ( nOrgSolveRet )  {

        // * ---------------------------------------------------------------------- *
        // * original solve handler exists and was been able to solve               *
        // * ---------------------------------------------------------------------- *
        return nOrgSolveRet;

    } else {

        // * ---------------------------------------------------------------------- *
        // * Format the Confirmation Message                                        *
        // *                                                                        *
        // * Your last change                                                       *
        // *                                                                        *
        // * \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!                                  *
        // * ---------------------------------------------------------------------- *
        if ( arg_list )  {
            arg_list[0]  = (char*) pThis->labelId(arg_list+1, lIndex);  // resurce ID -> label
            arg_list[10] = (char*) pThis->value(lIndex);                // new value
            arg_list[16] = (char*) pThis->unitId(arg_list+17, lIndex);

            MrProt _storage;
            void* pVoid     = pThis->storeInMemory(&_storage,lIndex);
            //    pThis->prot()   = *pOrigProt;
            pThis->recallMemory(pToAddMemory,pOrigProt,lIndex);

            arg_list[4]     = (char*) pThis->value(lIndex);            // old value

            pThis->recallMemory (pVoid,&_storage,lIndex);
            pThis->clearMemory  (pVoid,lIndex);

        }



        // * ---------------------------------------------------------------------- *
        // * Check whether 3D and flow compensation in slice direction is switched  *
        // * on at the same time.                                                   *
        // * ---------------------------------------------------------------------- *
        if ( pNewProt->kSpace().dimension() == SEQ::DIM_3  &&
             ( pNewProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_YES || pNewProt->flowComp()[0] == SEQ::FLOW_COMPENSATION_SLICESEL_ONLY ) ) {


            if ( arg_list ) {

                // * -------------------------------------------------------------- *
                // * has adapted the following parameters:                          *
                // *                                                                *
                // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                  *
                // * -------------------------------------------------------------- *
                arg_list[20] = (char*) MRI_STD_FLOW_COMP_LABEL;

                switch ( pOrigProt->flowComp()[0] ) {
                    case SEQ::FLOW_COMPENSATION_YES:
                        arg_list[24] = (char*) MRI_STD_FLOW_COMP_YES;
                        arg_list[30] = (char*) MRI_STD_FLOW_COMP_READ;
                    break;

                    case SEQ::FLOW_COMPENSATION_SLICESEL_ONLY:
                        arg_list[24] = (char*) MRI_STD_FLOW_COMP_SLICE;
                        arg_list[30] = (char*) MRI_STD_FLOW_COMP_NO;
                    break;
                }

                arg_list[36] = (char*) MRI_STD_EMPTY;
                arg_list[39] = (char*) '\n';
                lCtrl++;

            }



            // * ------------------------------------------------------------------ *
            // * Action:                                                            *
            // *                                                                    *
            // * Switch from "flow comp yes" to "flow comp read" or                 *
            // * from "flow comp slice" to "flow comp no".                          *
            // * ------------------------------------------------------------------ *
            LINK_SELECTION_TYPE* pFlowCompMode = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLOW_COMP);

            switch ( pOrigProt->flowComp()[0] ) {
                case SEQ::FLOW_COMPENSATION_YES:
                    if( ! pFlowCompMode || !pFlowCompMode->isEditable(0) || !pThis->seqLimits().getFlowCompensation()[0].hasOption(SEQ::FLOW_COMPENSATION_READOUT_ONLY) ) {  return ( 0 );  }
                    if(   pFlowCompMode->value( MRI_STD_FLOW_COMP_READ, 0 ) != MRI_STD_FLOW_COMP_READ)  {  return ( 0 );  }
                break;

                case SEQ::FLOW_COMPENSATION_SLICESEL_ONLY:
                    if( ! pFlowCompMode || !pFlowCompMode->isEditable(0) || !pThis->seqLimits().getFlowCompensation()[0].hasOption(SEQ::FLOW_COMPENSATION_NO) ) {  return ( 0 );  }
                    if(   pFlowCompMode->value( MRI_STD_FLOW_COMP_NO, 0 ) != MRI_STD_FLOW_COMP_NO)  {  return ( 0 );  }
                break;
            }

        }



        if ( (pOrigProt->fastImaging().segments() > 1) && (pNewProt->kSpace().dimension() == SEQ::DIM_3) )  {

            // * ------------------------------------------------------------------ *
            // * Reduce number of segments to one                                   *
            // * ------------------------------------------------------------------ *
            if ( !changeNoOfSegments ( pThis, arg_list, pOrigProt, 1, lCtrl ) )  {  return ( 0 );  }



            // * -------------------------------------------------------------- *
            // * If segments have been reduced to one, IR has to switched off.  *
            // * -------------------------------------------------------------- *
            if ( (    pOrigProt->preparationPulses().inversion() == SEQ::SLICE_SELECTIVE
                   || pOrigProt->preparationPulses().inversion() == SEQ::VOLUME_SELECTIVE )
                 && !switchOffInversionRecovery(pThis, arg_list, pOrigProt, lCtrl) ) {
                return ( 0 );
            }



            // * -------------------------------------------------------------- *
            // * If segments have been reduced to one, SR has to switched off.  *
            // * -------------------------------------------------------------- *
            if (     (    pOrigProt->preparationPulses().satRecovery() == SEQ::SATREC_SLICE_SELECTIVE
                       || pOrigProt->preparationPulses().satRecovery() == SEQ::SATREC_VOLUME_SELECTIVE)
                  && !switchOffSaturationRecovery (pThis, arg_list, pOrigProt, lCtrl) )   {
                return ( 0 );
            }



            // * -------------------------------------------------------------- *
            // * If segments have been reduced to one, dark blood preparation   *
            // * have to switched off.                                          *
            // * -------------------------------------------------------------- *
            if ( pOrigProt->preparationPulses().darkBlood() && !switchOffDarkBlood (pThis, arg_list, pOrigProt, lCtrl) )  {
                return ( 0 );
            }

        }



        // * ------------------------------------------------------------------ *
        // * Switch off SWI if the user selects 2D                              *
        // * ------------------------------------------------------------------ *
        if ( pOrigProt->SWI() && pNewProt->kSpace().dimension() == SEQ::DIM_2 )  {

            if ( arg_list ) {

                // * -------------------------------------------------------------- *
                // * has adapted the following parameters:                          *
                // *                                                                *
                // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                  *
                // * -------------------------------------------------------------- *
                arg_list[20] = (char*) MRI_STD_SWI_LABEL;
                arg_list[24] = (char*) MRI_STD_ON;
                arg_list[30] = (char*) MRI_STD_OFF;
                arg_list[36] = (char*) MRI_STD_EMPTY;
                arg_list[39] = (char*) '\n';
                lCtrl++;

            }

            // * ------------------------------------------------------------------ *
            // * Action:                                                            *
            // *                                                                    *
            // * Switch from "SWI on" to "SWI off" or                               *
            // * ------------------------------------------------------------------ *
            LINK_BOOL_TYPE* pSWI = _search<LINK_BOOL_TYPE>(pThis, MR_TAG_SWI);

            if( ! pSWI || !pSWI->isEditable(0) || !pThis->seqLimits().getSWI().hasOption(SEQ::OFF) ) {  return ( 0 );  }
            if(   pSWI->value( false, 0 ) != false)  {  return ( 0 );  }

        }




        // * ---------------------------------------------------------------------- *
        // * Check whether the conflict has been solved                             *
        // * ---------------------------------------------------------------------- *
        if( pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

            // * ------------------------------------------------------------------ *
            // * Conflict has been solved                                           *
            // * ------------------------------------------------------------------ *
            lCtrl--;
            if ( lCtrl <= 3 )  {
                if(arg_list) {  arg_list[39 + lCtrl * lOffset] = (char*) '\0';  }
            } else {
                if(arg_list) {  arg_list[39 + lCtrl * lOffset] = (char*) '\n';  }
            }
            return ( MRI_STD_CONFIRMATION_MSG );

        } else {

            // * ------------------------------------------------------------------ *
            // * The conflict has not been solved, yet.                             *
            // *                                                                    *
            // * Maybe that TE and/or TR have to be increased, additionally         *
            // * ------------------------------------------------------------------ *
            return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

        }

    }

}








unsigned fBSolveSWIConflict(LINK_BOOL_TYPE* const pThis,
                        char**                arg_list,       // receives confirmation message
                        const void*           pToAddMemory,   // not needed
                        const MrProt*         pOrigProt,      // Original protocol with old number of segments
                        long                  lCtrl          // Array index reserved
)
{

    MrProt      *pMrProt        = (MrProt*)  &pThis->prot();


    if ( pMrProt->SWI() )  {

        // * ------------------------------------------------------------------ *
        // * Check correct dimension                                            *
        // * ------------------------------------------------------------------ *
        if ( pThis->prot().kSpace().dimension() == SEQ::DIM_2 )  {

            LINK_SELECTION_TYPE* pDimensionMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_DIMENSION);

            // * Check whether rf spoiling can be selected *
            if( !pDimensionMode || !pDimensionMode->isEditable(0) || !pThis->seqLimits().getDimension().hasOption(SEQ::DIM_3) ) {  return ( 0 );  }

            // * Switch rf spoiling on *
            if(pDimensionMode->value( MRI_STD_DIMENSION_3D, 0 ) != MRI_STD_DIMENSION_3D)  {  return ( 0 );  }

        }



        // * ------------------------------------------------------------------ *
        // * Check correct flow compensation mode                               *
        // * ------------------------------------------------------------------ *
        if ( pThis->prot().flowComp()[0] != SEQ::FLOW_COMPENSATION_YES  )  {

            LINK_SELECTION_TYPE* pFlowComp = _searchElm<LINK_SELECTION_TYPE>(pThis, MR_TAG_FLOW_COMP);

            // * Check whether rf spoiling can be selected *
            if( !pFlowComp || !pFlowComp->isEditable(0) || !pThis->seqLimits().getFlowCompensation()[0].hasOption(SEQ::FLOW_COMPENSATION_YES) ) {  return ( 0 );  }

            // * Switch rf spoiling on *
            if(pFlowComp->value( MRI_STD_FLOW_COMP_YES, 0 ) != MRI_STD_FLOW_COMP_YES)  {  return ( 0 );  }

        }



        // * ------------------------------------------------------------------ *
        // * Check correct rf-spoiling mode                                     *
        // * ------------------------------------------------------------------ *
        if ( ! pMrProt->fastImaging().RFSpoiling() )  {

            LINK_BOOL_TYPE* pRFSpoiling = _search<LINK_BOOL_TYPE>(pThis, MR_TAG_RF_SPOILING);

            // * Check whether the UI parameter "rf-spoiling" exists and is *
            if(!pRFSpoiling || !pRFSpoiling->isEditable(0))  {  return ( 0 );  }
            // * Enable rf spoiling *
            if( pRFSpoiling->value(true,0) != true)  {  return ( 0 );  }
        }




        // * ------------------------------------------------------------------ *
        // * Check correct shim mode                                            *
        // * ------------------------------------------------------------------ *
        if ( pThis->prot().adjustment().adjShimMode() != SEQ::ADJSHIM_STANDARD )  {

            LINK_SELECTION_TYPE* pShimMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_SHIM_MODE);

            // * Check whether rf spoiling can be selected *
            if( !pShimMode || !pShimMode->isEditable(0) || !pThis->seqLimits().getAdjShim().hasOption(SEQ::ADJSHIM_STANDARD) ) {  return ( 0 );  }

            // * Switch rf spoiling on *
            if(pShimMode->value( MRI_STD_SHIM_MODE_STANDARD, 0 ) != MRI_STD_SHIM_MODE_STANDARD)  {  return ( 0 );  }

        }



        // * ------------------------------------------------------------------ *
        // * Check correct iPAT mode                                            *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->PAT().PATMode() == SEQ::PAT_MODE_SENSE )  {

            LINK_SELECTION_TYPE* pPATMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_PAT_MODE);

            // * Check whether rf spoiling can be selected *
            if( !pPATMode || !pPATMode->isEditable(0) || !pThis->seqLimits().getPATMode().hasOption(SEQ::PAT_MODE_GRAPPA) ) {  return ( 0 );  }

            // * Switch rf spoiling on *
            if(pPATMode->value( MRI_STD_PAT_MODE_GRAPPA, 0 ) != MRI_STD_PAT_MODE_GRAPPA)  {  return ( 0 );  }

        }


        // * ------------------------------------------------------------------ *
        // * Check correct uncombined images mode                               *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->uncombImages() )  {

            LINK_BOOL_TYPE* pUncombIma = _search<LINK_BOOL_TYPE>(pThis, MR_TAG_SAVE_UNCOMBINED);

            // * Check whether the UI parameter "elliptical scanning" exists, is editable and has the option "off" *
            if(!pUncombIma || !pUncombIma->isEditable(0))  {  return ( 0 );  }
            // * Switch off elliptical scanning *
            if( pUncombIma->value(false,0) != false)  {  return ( 0 );  }
        }



        // * ------------------------------------------------------------------ *
        // * Reconmode phase                                                    *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->reconstructionMode() == SEQ::RECONMODE_PHASE )  {
            LINK_SELECTION_TYPE* pReconMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_RECON_MODE);

            // * Check whether rf spoiling can be selected *
            if( !pReconMode || !pReconMode->isEditable(0) || !pThis->seqLimits().getReconstructionMode().hasOption(SEQ::RECONMODE_MAGN_PHASE) ) {  return ( 0 );  }

            // * Switch adaptive combine on *
            if(pReconMode->value( MRI_STD_RECON_MAGN_PHASE, 0 ) != MRI_STD_RECON_MAGN_PHASE)  {  return ( 0 );  }

        }



        // * ------------------------------------------------------------------ *
        // * Check correct coil combine mode                                    *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->coilCombineMode() == SEQ::COILCOMBINE_SUM_OF_SQUARES )  {

            LINK_SELECTION_TYPE* pCoilCombine = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_COIL_COMBINE_MODE);

            // * Check whether rf spoiling can be selected *
            if( !pCoilCombine || !pCoilCombine->isEditable(0) || !pThis->seqLimits().getCoilCombineMode().hasOption(SEQ::COILCOMBINE_ADAPTIVE_COMBINE) ) {  return ( 0 );  }

            // * Switch adaptive combine on *
            if(pCoilCombine->value( MRI_STD_COIL_COMBINE_ADAPTIVE_COMBINE, 0 ) != MRI_STD_COIL_COMBINE_ADAPTIVE_COMBINE)  {  return ( 0 );  }
        }



        // * ------------------------------------------------------------------ *
        // * Check iPAT in partitions direction                                 *
        // * ------------------------------------------------------------------ *
        if ( pMrProt->SWI() && (pMrProt->PAT().AccelFact3D() > 1) )  {

            LINK_LONG_TYPE* pAccelFact3D = _search<LINK_LONG_TYPE>(pThis, MR_TAG_PAT_ACC_3D);

            // * Check whether 3D acceleration factor can be set to 1 *
            if( !pAccelFact3D || !pAccelFact3D->isEditable(0) || (pThis->seqLimits().getAccelFactor3D().getMin() > 1) ) {  return ( 0 );  }

            // * Switch adaptive combine on *
            if(pAccelFact3D->value( 1, 0 ) != 1)  {  return ( 0 );  }


        }




        //  calls prepareForBinarySearch by default
        if( !pThis->tryProt((void*)pToAddMemory, pOrigProt, lCtrl) )  {

            return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

        } else {

            if(arg_list) {  arg_list[39] = (char*) '\0';  }
            return ( MRI_STD_CONFIRMATION_MSG );

        }

    } else {

        return ( 0 );

    }

}











bool fLocalUICSetWithinLimits(MrUILinkBase* const pThis, double& rVal, const char* tTag, unsigned long ulMode, long lIndex, MrUILinkLimited<double>::PFctSetValue pFctSet)
{
        return _fLocalUICSetWithinLimits(pThis, rVal, tTag, ulMode, lIndex, pFctSet);
}

bool fLocalUICSetWithinLimits(MrUILinkBase* const pThis, long& rVal, const char* tTag, unsigned long ulMode, long lIndex, MrUILinkLimited<long>::PFctSetValue pFctSet)
{
        return _fLocalUICSetWithinLimits(pThis, rVal, tTag, ulMode, lIndex, pFctSet);
}












/*[ Function ****************************************************************\
*
* Name        : fSSolveMSMModeConflict
*
* Description : Solve handler for the multi slice mode
*
*               Sequential slice mode is not allowed in combination with
*               respiratory triggering
*               The solve handler switches to triggering off if
*               sequential multi slice has been selected by the user
*
* Return      : confirmation message
*
\****************************************************************************/
unsigned fSSolveMSMModeConflict
(
    LINK_SELECTION_TYPE*  const pThis,     // ...
    char**                arg_list,        // receives confirmation message
    const void*           pToAddMemory,
    const MrProt*         pOrigProt, // Original protocol with old rf mode
    long                  lIndex     // Array index reserved
)
{

    MrProt*               pNewProt   = &pThis->prot();               // with new PhasePartialFourier value
    const long            lOffset    = 20;
    long                  lCtrl      = 0;



    // * -------------------------------------------------------------------------- *
    // * Format the Confirmation Message                                            *
    // *                                                                            *
    // * Your last change                                                           *
    // *                                                                            *
    // * \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!                                      *
    // * -------------------------------------------------------------------------- *
    if ( arg_list )  {
        arg_list[0]  = (char*) pThis->labelId(arg_list+1, lIndex);  // resurce ID -> label
        arg_list[10] = (char*) pThis->value(lIndex);                // new value
        arg_list[16] = (char*) pThis->unitId(arg_list+17, lIndex);

        MrProt _storage;
        void* pVoid     = pThis->storeInMemory(&_storage,lIndex);
        //    pThis->prot()   = *pOrigProt;
        pThis->recallMemory(pToAddMemory,pOrigProt,lIndex);

        arg_list[4]     = (char*) pThis->value(lIndex);            // old value

        pThis->recallMemory (pVoid,&_storage,lIndex);
        pThis->clearMemory  (pVoid,lIndex);
    }




    // * ------------------------------------------------------------------------- *
    // * Determine trigger mode and signal                                         *
    // * ------------------------------------------------------------------------- *
    SEQ::PhysioSignal PhysioSignalHigh          = SEQ::SIGNAL_NONE;   // * Selected trigger modes    *
    SEQ::PhysioSignal PhysioSignalLow           = SEQ::SIGNAL_NONE;   // * and signals               *
    SEQ::PhysioMethod PhysioMethodHigh          = SEQ::METHOD_NONE;
    SEQ::PhysioMethod PhysioMethodLow           = SEQ::METHOD_NONE;

    pNewProt->physiology().getPhysioMode (PhysioSignalHigh, PhysioMethodHigh, PhysioSignalLow, PhysioMethodLow);


    // * ------------------------------------------------------------------------- *
    // * The user has selected sequential slice but respiratory triggering is      *
    // * active.                                                                   *
    // * ------------------------------------------------------------------------- *
    if ( ( pNewProt->kSpace().multiSliceMode() == SEQ::MSM_SEQUENTIAL )  &&
         ( PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION ) && ( PhysioMethodHigh == SEQ::METHOD_TRIGGERING ) ) {

        if ( arg_list ) {

            // * ------------------------------------------------------------------ *
            // * has adapted the following parameters:                              *
            // *                                                                    *
            // * \t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!                      *
            // * ------------------------------------------------------------------ *
            arg_list[20] = (char*) MRI_STD_FIRST_SIGNAL_MODE_LABEL;
            arg_list[24] = (char*) MRI_STD_RESP_TRIGGER;
            arg_list[30] = (char*) MRI_STD_SIGNAL_NONE;
            arg_list[36] = (char*) MRI_STD_EMPTY;
            arg_list[39] = (char*) '\n';
            lCtrl++;

        }

        // * ---------------------------------------------------------------------- *
        // * Set new value                                                          *
        // * ---------------------------------------------------------------------- *
        LINK_SELECTION_TYPE* pFirstSignalMode = _search<LINK_SELECTION_TYPE>(pThis, MR_TAG_FIRST_SIGNAL_MODE);

        // * Check whether sequential slice mode can be selected *
        if(!pFirstSignalMode || !pFirstSignalMode->isEditable(0) || !pThis->seqLimits().getPhysioModes().hasOption(SEQ::SIGNAL_NONE, SEQ::METHOD_NONE, SEQ::SIGNAL_NONE, SEQ::METHOD_NONE) )  {  return ( 0 );  }

        // * Switch off physio mode *
        if(pFirstSignalMode->value ( MRI_STD_SIGNAL_NONE, 0 ) != MRI_STD_SIGNAL_NONE)  {  return ( 0 );  }


        //  calls prepareForBinarySearch by default
        if( !pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

            return ( fUICStandardSolveHandler ( pThis, arg_list, pOrigProt, pfGetMinTE, NULL, NULL, lCtrl ) );

        } else {

            if(arg_list) {  arg_list[39] = (char*) '\0';  }
            return ( MRI_STD_CONFIRMATION_MSG );

        }

    } else {

        // * ---------------------------------------------------------------------- *
        // * Switch off inversion recovery preparation                              *
        // * ---------------------------------------------------------------------- *
        if (    (pNewProt->kSpace().multiSliceMode()        == SEQ::MSM_INTERLEAVED)
             && (pOrigProt->preparationPulses().inversion() != SEQ::INVERSION_OFF  )
             && !switchOffInversionRecovery(pThis, arg_list, pOrigProt, lCtrl) ) {
            return ( 0 );
        }


        // * ---------------------------------------------------------------------- *
        // * Switch off saturation recovery preparation                             *
        // * ---------------------------------------------------------------------- *
        if (    (pNewProt->kSpace().multiSliceMode()          == SEQ::MSM_INTERLEAVED)
             && (pOrigProt->preparationPulses().satRecovery() != SEQ::SATREC_NONE)
             && !switchOffSaturationRecovery (pThis, arg_list, pOrigProt, lCtrl)      ) {
            return ( 0 );
        }



        // * ---------------------------------------------------------------------- *
        // * Switch off dark blood preparation                                      *
        // * ---------------------------------------------------------------------- *
        if (    (pNewProt->kSpace().multiSliceMode() == SEQ::MSM_INTERLEAVED)
             && pNewProt->preparationPulses().darkBlood()
             && !switchOffDarkBlood (pThis, arg_list, pOrigProt, lCtrl) ) {
            return ( 0 );
        }



        // * ---------------------------------------------------------------------- *
        // * Reduce number of segments to one                                       *
        // * ---------------------------------------------------------------------- *
        if (    ( pNewProt->fastImaging().segments() > 1 )
             && (!changeNoOfSegments ( pThis, arg_list, pOrigProt, 1, lCtrl ) ) ) {
            return ( 0 );
        }



        if( !pThis->tryProt((void*)pToAddMemory, pOrigProt, lIndex) )  {

            return ( 0 );

        } else {

            lCtrl--;
            if ( lCtrl <= 3 )  {
                if(arg_list) {  arg_list[39 + lCtrl * lOffset] = (char*) '\0';  }
            } else {
                if(arg_list) {  arg_list[39 + lCtrl * lOffset] = (char*) '\n';  }
            }
            return ( MRI_STD_CONFIRMATION_MSG );

        }

    }
}





/*[ Function ****************************************************************\
*
* Name        : fatSatModeGetOptionsHandler
*
* Description : Get-options handler for fat sat mode.
*
*               Fat mode weak / strong is allowed for quick fat, only
*
* Return      :
*
\****************************************************************************/
bool fatSatModeGetOptionsHandler (LINK_SELECTION_TYPE* const pThis, std::vector<unsigned>& rOptionVector, unsigned long& rulVerify, long )
{
    rOptionVector.clear();
    if (pThis->prot().preparationPulses().fatSuppression() != SEQ::FAT_SATURATION_QUICK) return ( false );

    rulVerify = LINK_SELECTION_TYPE::VERIFY_ON;

    const ParLimOption<SEQ::FatSatMode>& SeqOption = pThis->seqLimits().getFatSatMode();

    if( SeqOption.hasOption (SEQ::FAT_SAT_WEAK) )  {
        rOptionVector.push_back (MRI_STD_WEAK);
    }


    if( SeqOption.hasOption (SEQ::FAT_SAT_STRONG) )  {
        rOptionVector.push_back (MRI_STD_STRONG);
    }

    return ( rOptionVector.size() > 0 );
}





/*[ Function ****************************************************************\
*
* Name        : fatSatModeIsAvailable
*
* Description : is available handler for fat sat mode.
*
*               Fat mode weak / strong will only be displayed for quick fat sat
*
* Return      :
*
\****************************************************************************/
bool fatSatModeIsAvailable(LINK_SELECTION_TYPE* const pThis,long)
{
    if ( pThis->prot().preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK ) {
        return ( true );
    } else {
        return ( false );
    }
}

#ifdef SUPPORT_CT
//  The parameter is needed by the UILink layer in TimCT mode SEQ::MDS_ONCO
//  It returns the number of SMS segments per slice
static long fCTNSegGetVal(LINK_LONG_TYPE* const pThis, long)
{
    MrProt* pProt   = &pThis->prot();
    SeqLim* pSeqLim = &pThis->seqLimits();
    //  See
    static ReorderInfoGRE s_sReorder( 131072 );
    // Total number of measured lines without iPAT
    long lLinesToMeasureMax
        //  Total number of measured lines with iPAT
        ,  lLinesToMeasure
        //  k-space center line
        ,  lKSpaceCenterLine
        //  Total number of measured partitions without iPAT
        ,  lPartitionsToMeasureMax
        //  Total number of measured partitions with iPAT
        ,  lPartitionsToMeasure
        //  k-space center partition
        ,  lKSpaceCenterPartition
        ;

    NLS_STATUS lStatus = fCalculateReordering
        ( pProt
        , pSeqLim
        , &s_sReorder
        , &lLinesToMeasureMax
        , &lLinesToMeasure
        , &lKSpaceCenterLine
        , &lPartitionsToMeasureMax
        , &lPartitionsToMeasure
        , &lKSpaceCenterPartition
        );
    if( (lStatus & NLS_SEV) != NLS_SUCCESS )
    {
        return 0;
    }
    const int32_t i32NLinPerSeg = maximum(int32_t(1),pProt->MDS().mdsLinesPerSegment())
        ,         i32NSeg       = (lLinesToMeasure+i32NLinPerSeg-1)/i32NLinPerSeg
        ;
    return i32NSeg;
}

static bool fCTDixonIsAvailable(LINK_SELECTION_TYPE* const pThis, long)
{
    return pThis->seqLimits().getDixon().isAvailable()
        && pThis->prot().MDS().mdsModeMask() == SEQ::MDS_ONCO
        ;
}

static unsigned fCTDixonSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, long lPos)
{
    unsigned uRet = (*pOriginalDixonSetFct)(pThis,uNewVal,lPos);
    MrProt* pProt = &pThis->prot();
    if( pProt->dixon() != SEQ::DIXON_NONE )
    {
        //  Number of contrasts must be 2
        if( pProt->contrasts() != 2 )
        {
            if( LINK_LONG_TYPE* pContrasts = _search<LINK_LONG_TYPE>(pThis,MR_TAG_CONTRASTS) )
            {
                if( pContrasts->isAvailable(0) )
                {
                    pContrasts->value(2,0);
                }
            }
        }
        // First two in-phase echoes
        MeasNucleus sNucleus(pProt->txSpec().nucleusInfoArray()[0].nucleus());
        const double dConst  = -.5*1e6/(sNucleus.getLarmorConst()*SysProperties::getNominalB0()*CHEMICAL_SHIFT_FAT_PPM);
        double       dTE1_us = dConst
            ,        dTE2_us = 2*dConst
            ;


        Sequence*  pSeq   = &pThis->sequence();
        SeqLim*  pSeqLim  = &pSeq->getSEQ_BUFF()->getSeqLim();
        SeqExpo* pSeqExpo = &pSeq->getSEQ_BUFF()->getSeqExpo();

        //pSeq->prepareForBinarySearch(pProt);
        long alTEFill_us[lNoOfContrasts];
        long alTEMin_us [lNoOfContrasts];
        NLS_STATUS lStatus = fCalculateTEMin(pProt,pSeqLim,pSeqExpo,pSeq,alTEMin_us,alTEFill_us);
        if( (lStatus&NLS_SEV) == NLS_SUCCESS )
        {
            int iCntr = 1;
            while(dTE1_us < alTEMin_us[0] )
            {
                iCntr += 2;
                dTE1_us = iCntr*dConst;
                dTE2_us = (iCntr+1)*dConst;
            }
            if( fabs(pProt->te()[0]-dTE1_us) > 10 || fabs(pProt->te()[1]-dTE2_us) > 10 )
            {
                if( LINK_DOUBLE_TYPE* pTE = _searchElm<LINK_DOUBLE_TYPE>(pThis,MR_TAG_TE) )
                {
                    if( pTE->isAvailable(0) )
                    {
                        pTE->value(dTE1_us/1000,0);
                    }
                    if( pTE->isAvailable(1) )
                    {
                        pTE->value(dTE2_us/1000,1);
                    }
                }
            }
        }
    }
    return uRet;
}
static unsigned _formatFloat(char* arg_list[], double _val, int _prec)
{
    char _str[16];
    sprintf(_str,"%.*f",_prec,_val);
    int _x = 0;
    int _y = 0;
    sscanf(_str,"%d.%d",&_x,&_y);
    arg_list[0] = (char*) abs(_x);
    arg_list[1] = (char*)(_prec ? '.' : ' ');
    arg_list[2] = (char*)_prec;
    arg_list[3] = (char*)_y;
    return atof(_str) < 0. ? MRI_STD_NEG_FLOAT : MRI_STD_FLOAT;
}
static unsigned fCTSolveDixonConflict(LINK_SELECTION_TYPE* const pThis, char* arg_list[], const void* pVoid, const MrProt* pOrgProt, long lPos)
{
    MrProt* pProt = &pThis->prot();
    if( pProt->dixon() != SEQ::DIXON_NONE && pOrgProt->dixon() == SEQ::DIXON_NONE && pProt->contrasts() == 2 )
    {
        if( LINK_DOUBLE_TYPE* pBW = _searchElm<LINK_DOUBLE_TYPE>(pThis,MR_TAG_BANDWIDTH) )
        {
            if( pBW->isAvailable(0) )
            {
                //  Retrieve Limits
                std::vector<MrLimitDouble> sLimitVector;
                unsigned long ulVerify = LINK_DOUBLE_TYPE::VERIFY_OFF;
                pBW->getLimits(sLimitVector,ulVerify,0);
                double dBW1_Hz_pix =  pBW->value(0);
                double dBW2_Hz_pix =  pBW->value(1);
                std::vector<MrLimitDouble>::iterator it = sLimitVector.begin();
                while(it != sLimitVector.end() && (*it).minimum() < dBW1_Hz_pix )
                {
                    ++it;
                }
                if( it == sLimitVector.end() )
                {
                    //  Maximum already set
                    return 0;
                }
                unsigned uVal = pThis->value(lPos);
                unsigned uRet = 0;

                if( arg_list == 0 )
                {
                    //  Start with maximum and go down
                    std::vector<MrLimitDouble>::iterator it2 = sLimitVector.end();
                    while( --it2 != it )
                    {
                        //  Increase BW
                        pBW->value((*it2).minimum(),0);
                        pBW->value((*it2).minimum(),1);
                        //  Set Value
                        pThis->value(uVal,lPos);

                        if( pThis->tryProt(const_cast<void*>(pVoid),pOrgProt,lPos) )
                        {
                            //  Success
                            uRet = MRI_STD_CONFIRMATION_MSG;
                            break;
                        }
                    }
                    return uRet;
                }
                while( ++it != sLimitVector.end() )
                {
                    //  Increase BW
                    pBW->value((*it).minimum(),0);
                    pBW->value((*it).minimum(),1);
                    //  Set Value
                    pThis->value(uVal,lPos);

                    if( pThis->tryProt(const_cast<void*>(pVoid),pOrgProt,lPos) )
                    {
                        if( arg_list != 0 )
                        {
                            //  ////////////////////////////////////////////////////
                            //  Format the Confirmation Message
                            //  ////////////////////////////////////////////////////
                            //
                            //  Your last change
                            //
                            //  \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!
                            //
                            arg_list[0]  = (char*) pThis->labelId(arg_list+1,lPos);
                            arg_list[4]  = (char*) MRI_STD_DIXON_NONE;
                            arg_list[10] = (char*) pThis->value(lPos);
                            arg_list[16] = (char*) MRI_STD_EMPTY;
                            //
                            //  has adapted the following parameters
                            //
                            int _offset = 0;
                            int _total = maximum(2,int(pOrgProt->contrasts()));
                            for (int _i = 0; _i < _total;++_i,_offset += 20)
                            {

                                //  \t%21!r!:\t%24!r!\t--->\t%30!r!\t%36!r!%40!c!
                                //
                                //arg_list[20+_offset]  = (char*) (_storage->m_size > 1 ? MRI_STD_ELEMENT_LABEL : MRI_STD_BW_LABEL);
                                arg_list[21+_offset]  = (char*) pBW->labelId(arg_list+22+_offset,_i);
                                arg_list[24+_offset]  = (char*) _formatFloat(arg_list+25+_offset,_i == 0 ? dBW1_Hz_pix : dBW2_Hz_pix,1);
                                arg_list[30+_offset]  = (char*) _formatFloat(arg_list+31+_offset,pBW->value(_i),1);
                                arg_list[36+_offset]  = (char*) pBW->unitId(arg_list+37+_offset,_i);
                                arg_list[39+_offset]  = (char*) (_i == _total-1 ? '\0' : '\n');
                            }
                        }
                        //  Success
                        uRet = MRI_STD_CONFIRMATION_MSG;

                        break;
                    }
                }
                return uRet;
            }
        }
    }
    return 0;

}


static unsigned fCTMDSModeSetValue(LINK_SELECTION_TYPE* const pThis, unsigned uNewVal, long lPos)
{
    const MrProt* pProt = &pThis->prot();
    if( (uNewVal == MRI_STD_ON) && (pProt->kSpace().multiSliceMode() != SEQ::MSM_INTERLEAVED) )
    {
        if( LINK_SELECTION_TYPE* pSel = _search<LINK_SELECTION_TYPE>(pThis,MR_TAG_MULTI_SLICE_MODE) )
        {
            if( pSel->isAvailable(0) )
            {
                // Switch to interleaved mode. Forced mode is needed, since it may be necessary
                // to set the number of segments to one.
                pSel->value(MRI_STD_MSM_INTERLEAVED,0,MrUILinkBase::SET_FORCED);
            }
        }
    }
    if( pOriginalMDSModeSetFct != 0 )
    {
        return (*pOriginalMDSModeSetFct)(pThis,uNewVal,lPos);
    }
    return pThis->value(lPos);
}

static unsigned fCTSolveMDSModeConflict(LINK_SELECTION_TYPE * const pThis, char *arg_list[], const void *pVoid, const MrProt *pOrigProt, long lPos)
{
    MrProt* pProt = &pThis->prot();
    if( pProt->MDS().mdsModeMask() == SEQ::MDS_OFF )
    {
        if( pProt->dixon() != SEQ::DIXON_NONE )
        {
            //  Dixon not allowed
            if( LINK_SELECTION_TYPE* pDixon = _search<LINK_SELECTION_TYPE>(pThis,MR_TAG_DIXON) )
            {
                //  Turn Dixon-off
                //  This requires to turn TimCT-Onco temporary on.
                //  Otherwise the Dixon parameter is not available.
                pProt->MDS().mdsModeMask(SEQ::MDS_ONCO);
                if( pDixon->isAvailable(0) )
                {
                    pDixon->value(MRI_STD_DIXON_NONE,0);
                }
                pProt->MDS().mdsModeMask(SEQ::MDS_OFF);

                if( pProt->dixon() == SEQ::DIXON_NONE )
                {
                    if( arg_list )
                    {
                        //  Format header of confirmation message

                        //  Your last change
                        //
                        //  \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!

                        arg_list[0]  = (char*) pThis->labelId(arg_list+1,lPos);
                        arg_list[4]  = (char*) MRI_STD_ON;
                        arg_list[10] = (char*) MRI_STD_OFF;
                        arg_list[16] = (char*) MRI_STD_EMPTY;

                        //  has adapted the following parameters
                        //  \n\n\t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!
                        //
                        arg_list[20]  = (char*) pDixon->labelId(arg_list+21,0);
                        switch( pOrigProt->dixon() )
                        {
                        case SEQ::DIXON_WATER_IMAGE:
                            arg_list[24]  = (char*) MRI_STD_DIXON_WATER_IMAGE;
                            break;
                        case SEQ::DIXON_FAT_IMAGE:
                            arg_list[24]  = (char*) MRI_STD_DIXON_FAT_IMAGE;
                            break;
                        case SEQ::DIXON_WATER_FAT_IMAGES:
                            arg_list[24]  = (char*) MRI_STD_DIXON_WATER_FAT_IMAGES;
                            break;
                        default:
                            //  We should never get here
                            arg_list[24]  = (char*) MRI_STD_DIXON_NONE;
                        }
                        arg_list[30]  = (char*) MRI_STD_DIXON_NONE;
                        arg_list[36]  = (char*) pDixon->unitId(arg_list+37,0);
                        arg_list[39]  = (char*) '\0';
                    }
                    if( pThis->sequence().prepareForBinarySearch(pProt) )
                    {
                        //  Finished
                        return MRI_STD_CONFIRMATION_MSG;
                    }
                }
            }
        }  //  End dixon
        //  Try standard solve handler
        return fUICSolveSelectionConflict(pThis,arg_list,pVoid,pOrigProt,0,pfGetMinTE,0,0);
    }
    else  //  TimCT is turned on
    {
        //  In principle the following restrictions are not sequence specific
        //  (they are enforced via SeqLoop::ct_try_prot) and therefore the solve handler
        //  shouldn't be implemented in the sequence.
        //  However, multiple distributed solve-handlers make the code unclear.
        //  Therefore ...
        int  iOffset = 0;
        bool bSolved = false;
        if( pProt->preScanNormalizeFilter().on() )
        {
            if( LINK_BOOL_TYPE* pPreScanNormalize = _search<LINK_BOOL_TYPE>(pThis,MR_TAG_FLT_PRESCAN_NORMALIZE) )
            {
                //  Turn Prescan Normalize off
                //  We turn TimCT-Onco temporary off, so that we can disable
                //  the filter, if TimCT-Onco is used
                pProt->MDS().mdsModeMask(SEQ::MDS_OFF);
                if( pPreScanNormalize->isAvailable(0) )
                {
                    pPreScanNormalize->value(false,0);
                }
                pProt->MDS().mdsModeMask(SEQ::MDS_ONCO);

                if( !pProt->preScanNormalizeFilter().on() )
                {
                    bSolved = true;
                    if( arg_list )
                    {
                        //  Format header of confirmation message
                        if( iOffset ) arg_list[39+iOffset] = (char*) '\n';

                        if( iOffset < 60 )
                        {
                            //  has adapted the following parameters
                            //  \n\n\t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!
                            //
                            arg_list[20+iOffset]  = (char*) pPreScanNormalize->labelId(arg_list+21+iOffset,0);
                            arg_list[24+iOffset]  = (char*) MRI_STD_ON;
                            arg_list[30+iOffset]  = (char*) MRI_STD_OFF;
                            arg_list[36+iOffset]  = (char*) pPreScanNormalize->unitId(arg_list+37+iOffset,0);
                            arg_list[39+iOffset]  = (char*) '\0';
                            iOffset += 20;
                        }
                    }
                }
            }
        }  //  End: Pre-Scan Normalize
        if( pProt->BiFiCFilter().on() )
        {
            if( LINK_BOOL_TYPE* pB1Filter = _search<LINK_BOOL_TYPE>(pThis,MR_TAG_FLT_BIFIC) )
            {
                //  Turn Prescan Normalize off
                //  We turn TimCT-Onco temporary off, so that we can disable
                //  the filter, if TimCT-Onco is used
                pProt->MDS().mdsModeMask(SEQ::MDS_OFF);
                if( pB1Filter->isAvailable(0) )
                {
                    pB1Filter->value(false,0);
                }
                pProt->MDS().mdsModeMask(SEQ::MDS_ONCO);

                if( !pProt->BiFiCFilter().on() )
                {
                    bSolved = true;
                    if( arg_list )
                    {
                        //  Format header of confirmation message
                        if( iOffset ) arg_list[39+iOffset] = (char*) '\n';

                        if( iOffset < 60 )
                        {
                            //  has adapted the following parameters
                            //  \n\n\t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!
                            //
                            arg_list[20+iOffset]  = (char*) pB1Filter->labelId(arg_list+21+iOffset,0);
                            arg_list[24+iOffset]  = (char*) MRI_STD_ON;
                            arg_list[30+iOffset]  = (char*) MRI_STD_OFF;
                            arg_list[36+iOffset]  = (char*) pB1Filter->unitId(arg_list+37+iOffset,0);
                            arg_list[39+iOffset]  = (char*) '\0';
                            iOffset += 20;
                        }
                    }
                }
            }
        }  //  End: B1 Filter
        if( pProt->normalizeFilter().on() )
        {
            if( LINK_BOOL_TYPE* pNormalize = _search<LINK_BOOL_TYPE>(pThis,MR_TAG_FLT_NORM) )
            {
                //  Turn Normalize off
                //  We turn TimCT-Onco temporary off, so that we can disable
                //  the filter, if TimCT-Onco is used
                pProt->MDS().mdsModeMask(SEQ::MDS_OFF);
                if( pNormalize->isAvailable(0) )
                {
                    pNormalize->value(false,0);
                }
                pProt->MDS().mdsModeMask(SEQ::MDS_ONCO);

                if( !pProt->normalizeFilter().on() )
                {
                    bSolved = true;
                    if( arg_list )
                    {
                        //  Format header of confirmation message
                        if( iOffset ) arg_list[39+iOffset] = (char*) '\n';

                        if( iOffset < 60 )
                        {
                            //  has adapted the following parameters
                            //  \n\n\t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!
                            //
                            arg_list[20+iOffset]  = (char*) pNormalize->labelId(arg_list+21+iOffset,0);
                            arg_list[24+iOffset]  = (char*) MRI_STD_ON;
                            arg_list[30+iOffset]  = (char*) MRI_STD_OFF;
                            arg_list[36+iOffset]  = (char*) pNormalize->unitId(arg_list+37+iOffset,0);
                            arg_list[39+iOffset]  = (char*) '\0';
                            iOffset += 20;
                        }
                    }
                }
            }
        }  //  End: Normalize
        if( (pProt->distortionCorrFilter().mode() != SEQ::DISTCORR_NDIS) && (pProt->kSpace().multiSliceMode() != SEQ::MSM_SINGLESHOT) )
        {
            if( LINK_BOOL_TYPE* pDistCorr = _search<LINK_BOOL_TYPE>(pThis,MR_TAG_FLT_DISCOR) )
            {
                //  Turn Distortion correction off
                //  We turn TimCT-Onco temporary off, so that we can disable
                //  the filter, if TimCT-Onco is used
                pProt->MDS().mdsModeMask(SEQ::MDS_OFF);
                if( pDistCorr->isAvailable(0) )
                {
                    pDistCorr->value(false,0);
                }
                pProt->MDS().mdsModeMask(SEQ::MDS_ONCO);

                if( pProt->distortionCorrFilter().mode() == SEQ::DISTCORR_NDIS )
                {
                    bSolved = true;
                    if( arg_list )
                    {
                        //  Format header of confirmation message
                        if( iOffset ) arg_list[39+iOffset] = (char*) '\n';

                        if( iOffset < 60 )
                        {
                            //  has adapted the following parameters
                            //  \n\n\t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!
                            //
                            arg_list[20+iOffset]  = (char*) pDistCorr->labelId(arg_list+21+iOffset,0);
                            arg_list[24+iOffset]  = (char*) MRI_STD_ON;
                            arg_list[30+iOffset]  = (char*) MRI_STD_OFF;
                            arg_list[36+iOffset]  = (char*) pDistCorr->unitId(arg_list+37+iOffset,0);
                            arg_list[39+iOffset]  = (char*) '\0';
                            iOffset += 20;
                        }
                    }
                }
            }
        } // End: distortion correction
        if( pProt->intro() )
        {
            if( LINK_BOOL_TYPE* pIntro = _search<LINK_BOOL_TYPE>(pThis,MR_TAG_INTRO) )
            {
                pProt->MDS().mdsModeMask(SEQ::MDS_OFF);
                if( pIntro->isAvailable(0) )
                {
                    //  Turn Introduction off
                    pIntro->value(false,0);
                }
                pProt->MDS().mdsModeMask(SEQ::MDS_ONCO);

                if( !pProt->intro() )
                {
                    bSolved = true;
                    if( arg_list )
                    {
                        if( iOffset ) arg_list[39+iOffset] = (char*) '\n';

                        if( iOffset < 60 )
                        {
                            //  has adapted the following parameters
                            //  \n\n\t%21!r!:\t%25!r!\t--->\t%31!r!\t%37!r!%40!c!
                            //
                            arg_list[20+iOffset]  = (char*) pIntro->labelId(arg_list+21+iOffset,0);
                            arg_list[24+iOffset]  = (char*) MRI_STD_ON;
                            arg_list[30+iOffset]  = (char*) MRI_STD_OFF;
                            arg_list[36+iOffset]  = (char*) pIntro->unitId(arg_list+37+iOffset,0);
                            arg_list[39+iOffset]  = (char*) '\0';
                            iOffset += 20;
                        }
                    }
                }
            }
        }  //  End: Introduction
        if( !pThis->sequence().prepareForBinarySearch(pProt) )
        {
            //  try to increase TR
            char* arg_list2[MrUILinkBase::ARG_LIST_MAX_SIZE];
            memset(arg_list2,0,MrUILinkBase::ARG_LIST_MAX_SIZE*sizeof(char*));
            unsigned uResult = fUICSolveSelectionConflict(pThis,arg_list != 0 ? arg_list2 : 0,pVoid,pOrigProt,0,pfGetMinTE,0,0);
            if( uResult == MRI_STD_CONFIRMATION_MSG )
            {
                bSolved = true;
                if( arg_list )
                {
                    if( iOffset ) arg_list[39+iOffset] = (char*) '\n';
                    if( iOffset < 60 )
                    {
                        memcpy(arg_list+20+iOffset,arg_list2+20,19);
                        arg_list[39+iOffset]  = (char*) '\0';
                        iOffset += 20;
                    }
                }
            }
        }  // End standard solve handler
        if( bSolved )
        {
            if( arg_list )
            {
                //  Format header of confirmation message

                //  Your last change
                //
                //  \t%1!r!:\t%5!r!\t--->\t%11!r!\t%17!r!

                arg_list[0]  = (char*) pThis->labelId(arg_list+1,lPos);
                arg_list[4]  = (char*) MRI_STD_OFF;
                arg_list[10] = (char*) MRI_STD_ON;
                arg_list[16] = (char*) MRI_STD_EMPTY;
            }
            if( pThis->sequence().prepareForBinarySearch(pProt) )
            {
                //  Finished
                return MRI_STD_CONFIRMATION_MSG;
            }
        }
    }
    return 0;
}
#endif // SUPPORT_CT


#endif










/************************/
/* End of module header */
/************************/




#ifndef VXWORKS
// ------------------------------------------------------------------------------
// Functions   : fSEQConvProt
// ------------------------------------------------------------------------------
//
// Description : try to convert protocols from previous software versions
//
// Return      : SEQU_NORMAL for success
//               SEQU_ERROR  for error
//
// ------------------------------------------------------------------------------
NLS_STATUS fSEQConvProt (const MrProt &rMrProtSrc, MrProt &rMrProtDst)
{
    static const char * const ptModule = {"fSEQConvProt"};
    NLS_STATUS  lRet = SEQU_NORMAL;



    if( rMrProtSrc.getConvFromVersion() < 21310000 )  {

        // * ---------------------------------------------------------------------- *
        // * Flow compensation mode                                                 *
        // * ---------------------------------------------------------------------- *
        if ( rMrProtSrc.flowComp()[0] == SEQ::FLOW_COMPENSATION_YES )  {
             rMrProtDst.flowComp()[0] = SEQ::FLOW_COMPENSATION_SLICE_READ;
        }



        // * ---------------------------------------------------------------------- *
        // * Readout polarity                                                       *
        // * ---------------------------------------------------------------------- *
        const unsigned long ulROPolarity = 1UL;     // * Parameter index            *
        const unsigned long ulBipolar    = 1UL;     // * Selection: bipolar         *
        const unsigned long ulMonopolar  = 2UL;     // * Selection: monopolar       *


        if ( rMrProtSrc.wipMemBlock().alFree[ulROPolarity] != 0 )  {

            switch ( rMrProtSrc.wipMemBlock().alFree[ulROPolarity] )  {
                case ulBipolar   : rMrProtDst.readOutMode (SEQ::READOUT_BIPOLAR);
                case ulMonopolar : rMrProtDst.readOutMode (SEQ::READOUT_MONOPOLAR);
            }

            rMrProtDst.wipMemBlock().alFree[ulROPolarity] = 0L;

        }



        // * ---------------------------------------------------------------------- *
        // * SWI                                                                    *
        // * ---------------------------------------------------------------------- *
        const unsigned long ulSWI           = 2UL;      // * Parameter index                *
        const unsigned long ulSlideWin      = 8UL;
        const unsigned long ulMaskNo        = 9UL;
        const unsigned long ulSlices        = 10UL;
        const unsigned long ulVenCheckBox   = 11UL;
        const unsigned long ulArtCheckBox   = 12UL;
        const unsigned long ulOn            = 4UL;      // * Selection: On                  *
        const unsigned long ulOff           = 2UL;      // * Selection: Off                 *


        if ( rMrProtSrc.wipMemBlock().alFree[ulSWI] != 0 )  {

            switch ( rMrProtSrc.wipMemBlock().alFree[ulSWI] )  {
                case ulOn  : rMrProtDst.SWI(true );
                case ulOff : rMrProtDst.SWI(false);
            }

            rMrProtDst.wipMemBlock().alFree[ulSWI]          = 0L;
            rMrProtDst.wipMemBlock().alFree[ulSlideWin]     = 0L;
            rMrProtDst.wipMemBlock().alFree[ulMaskNo]       = 0L;
            rMrProtDst.wipMemBlock().alFree[ulSlices]       = 0L;
            rMrProtDst.wipMemBlock().alFree[ulVenCheckBox]  = 0L;
            rMrProtDst.wipMemBlock().alFree[ulArtCheckBox]  = 0L;

        }

    }




#ifdef SUPPORT_PACE
    if( !SeqLoopBH::convProt(rMrProtSrc,rMrProtDst) )
    {
        TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s SeqLoopBH::convProt(.,.) failed!", ptModule, SEQU_ERROR);
        lRet = SEQU_ERROR;
    }
#endif // SUPPORT_PACE
    return lRet;
}



#endif // ndef VXWORKS







/*[ Function ****************************************************************\
*
* Name        : fSEQInit
*
* Description : Defines the hard limits for the Seq/Change dialog.
*
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQInit
(
  SeqLim     *pSeqLim
)

{
    static const char *ptModule = {"fSEQInit"};
    NLS_STATUS  lStatus = SEQU__NORMAL;
    double      dMin, dMax, dInc, dDef;
    long        lMin, lMax, lInc, lDef;
    long        lI;                         // * Loop counter *


    UNUSED_ARG ( tString );



    // * ---------------------------------------------------------------------- *
    // * Declaration of required proxies                                        *
    // * ---------------------------------------------------------------------- *
    GCProxy     theGCProxy;
    MSUProxy    theMSUProxy;
    GPAProxy    theGPAProxy;

    #ifndef VXWORKS
        char      tMrTag[64];
    #endif




    // * ---------------------------------------------------------------------- *
    // * Initialization                                                         *
    // * ---------------------------------------------------------------------- *
    dMin = dMax = dInc = dDef = 0.;
    lMin = lMax = lInc = lDef = 0;



    // * --------------------------------------------------------------------------- *
    // * Add RSats to SBBList                                                        *
    // * --------------------------------------------------------------------------- *
    for ( lI=0; lI<lNoOfRSats; lI++ )  {
        SBBRSat[lI].addToSBBList(&SBB);
    }



    // * ---------------------------------------------------------------------- *
    // * Sequence imprint                                                       *
    // * ---------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
            pSeqLim->setSequenceHintText( (char *) "\n\
            Application: Gradient echo sequence for T1 and T2* weighted imaging.                  \n\
                 Basics: 2D and 3D imaging; up to 2 contrasts; segmentation; inversion pulse      \n\
                 Build : "__DATE__"                                                               \n");
        break;

        case RT_Gre:
            pSeqLim->setSequenceHintText( (char *) "\n\
            Application: Real time 2D gradient echo imaging                                       \n\
                         Only single slice.                                                       \n\
                 Basics: 2D imaging with Fisp or FLASH contrast                                   \n");
        break;

        case Gre_HR:
            pSeqLim->setSequenceHintText( (char *) "\n\
            Application: Gradient echo sequence for T1 and T2* weighted imaging.                  \n\
                 Basics: 2D and 3D imaging; up to 2 contrasts; segmentation; inversion pulse      \n\
                 Build : "__DATE__"                                                               \n");
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }


    pSeqLim->setMyOrigFilename                        (                    (char *) __FILE__);
    pSeqLim->setSequenceOwner                         (                    SEQ_OWNER_SIEMENS);


    // * ---------------------------------------------------------------------------------- *
    // * Specify sequence type                                                              *
    // * ---------------------------------------------------------------------------------- *
    pSeqLim->setSequenceType                                         (SEQ::SEQUENCE_TYPE_GRE);
    pSeqLim->getSequenceType().setDisplayMode                                   (SEQ::DM_OFF);



    if ( eSequence == RT_Gre ) {
        RTControl.setUseSpaceMouse(true);
    }
    RTControl.init(pSeqLim);


    if ( eSequence == RT_Gre ) {
        // * Allow the "do not save"-mode and make it the default for gre_rt
        pSeqLim->enableChangeStoreImages();
        pSeqLim->disableStoreImages();
        // * Allow endless measurement
        pSeqLim->setEndlessMeasurement                                      (SEQ::OFF, SEQ::ON);
        pSeqLim->getEndlessMeasurement().setDisplayMode                          (SEQ::DM_EDIT);
        // tell ICE that this is an interactive sequence where repetitions are not included in the data amount calculation
        pSeqLim->setInteractiveRealtime                                               (SEQ::ON);
    }


    // * ------------------------------------------------------------------------------------ *
    // * the system requirements: frequency, and gradient power                               *
    // * ------------------------------------------------------------------------------------ *
    pSeqLim->setAllowedFrequency                                            (8000000,500000000);
    pSeqLim->setRequiredGradAmpl                                                         ( 1.0);
    pSeqLim->setRequiredGradSlewRate                                                     ( 5.0);



    if (theGPAProxy.getGradMaxAmplNominal() < LOW_GRADIENT_SYSTEM)  {
        pSeqLim->setGradients                    (SEQ::GRAD_FAST, SEQ::GRAD_FAST_GSWD_RISETIME);
    } else {
      pSeqLim->setGradients               (SEQ::GRAD_FAST, SEQ::GRAD_NORMAL, SEQ::GRAD_WHISPER,
                                  SEQ::GRAD_FAST_GSWD_RISETIME, SEQ::GRAD_NORMAL_GSWD_RISETIME,
                                                               SEQ::GRAD_WHISPER_GSWD_RISETIME);
    }


    // * ------------------------------------------------------------------------------------- *
    // * Matrix                  (         min,         max,         inc,         def, psLIM );*
    // * ------------------------------------------------------------------------------------- *

    // * ------------------------------------------------------------------------------------- *
    // * Sequence Dimension                                                                    *
    // * ------------------------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            pSeqLim->setDimension                                    ( SEQ::DIM_2, SEQ::DIM_3 );
        break;

        case RT_Gre:
            pSeqLim->setDimension                                                ( SEQ::DIM_2 );
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }



    // * ------------------------------------------------------------------------------------- *
    // * Base matrix size of the image                                                         *
    // * ------------------------------------------------------------------------------------- *
    if (theMSUProxy.getNominalB0() < LOW_FIELD_B0)  {
        pSeqLim->setBaseResolution      (          64,         512,  SEQ::INC_64,          128);

        pSeqLim->setBandWidthPerPixel   (  0,         20,        130,         10,          130);
        pSeqLim->setBandWidthPerPixel   (  1,         20,        130,         10,          130);
    } else {
        pSeqLim->setBaseResolution      (          64,        1024,  SEQ::INC_64,          128);

        pSeqLim->setBandWidthPerPixel   (  0,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  1,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  2,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  3,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  4,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  5,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  6,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  7,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  8,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   (  9,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   ( 10,         30,       2000,         10,          260);
        pSeqLim->setBandWidthPerPixel   ( 11,         30,       2000,         10,          260);
    }


    pSeqLim->setAsymmetricEcho                                            ( SEQ::OFF, SEQ::ON );
    pSeqLim->setReadOutMode                      (SEQ::READOUT_BIPOLAR, SEQ::READOUT_MONOPOLAR);


    // * ------------------------------------------------------------------------------------- *
    // * Lines                                                                                 *
    // * ------------------------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            if (theMSUProxy.getNominalB0() < LOW_FIELD_B0)  {
                pSeqLim->setPELines       (        64,         512,SEQ::INC_GRE_SEGMENTS,  128);
            } else {
                pSeqLim->setPELines       (        64,        1024,SEQ::INC_GRE_SEGMENTS,  128);
            }

            pSeqLim->setSegments          (         1,         127,           2,             1);
        break;

        case RT_Gre:
            if (theMSUProxy.getNominalB0() < LOW_FIELD_B0)  {
                pSeqLim->setPELines       (        64,         512,           1,           128);
            } else {
                pSeqLim->setPELines       (        64,        1024,           1,           128);
            }

            pSeqLim->setSegments          (         1,           1,           1,             1);
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }

    pSeqLim->setPhaseOversampling         (     0.000,       1.000,       0.010,         0.000);
    pSeqLim->setPhasePartialFourierFactor               (SEQ::PF_OFF, SEQ::PF_7_8, SEQ::PF_6_8);


    // * ------------------------------------------------------------------------------------- *
    // * 3D-Partitions                                                                         *
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setPartition                 (        16,         256,           2,            32);
    pSeqLim->setImagesPerSlab             (        16,         512,           2,            32);
    pSeqLim->setSliceOversampling         (     0.000,       1.000,       0.010,         0.000);
    pSeqLim->setSlicePartialFourierFactor               (SEQ::PF_OFF, SEQ::PF_7_8, SEQ::PF_6_8);





    // * ------------------------------------------------------------------------------------- *
    // * Timing                        (  No.,      min,      max,      inc,      def, psLIM );*
    // * ------------------------------------------------------------------------------------- *

    // * ------------------------------------------------------------------------------------- *
    // * Set the TE, TI, TD, and TR limits and number of contrasts                             *
    // * ------------------------------------------------------------------------------------- *
    long  lTEMaxEasySite;

    if ( fSUIsEasySite (&lTEMaxEasySite) )  {

        pSeqLim->setTE            (    0,         500,       60000,          10,         10000);
        pSeqLim->setTE            (    1,         500,       80000,          10,         15000);
        pSeqLim->setTE            (    2,         500,       80000,          10,         20000);
        pSeqLim->setTE            (    3,         500,       80000,          10,         25000);
        pSeqLim->setTE            (    4,         500,       80000,          10,         30000);
        pSeqLim->setTE            (    5,         500,       80000,          10,         35000);
        pSeqLim->setTE            (    6,         500,       80000,          10,         40000);
        pSeqLim->setTE            (    7,         500,       80000,          10,         45000);
        pSeqLim->setTE            (    8,         500,       80000,          10,         50000);
        pSeqLim->setTE            (    9,         500,       80000,          10,         55000);
        pSeqLim->setTE            (   10,         500,       80000,          10,         60000);
        pSeqLim->setTE            (   11,         500,       80000,          10,         65000);



        if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
            pSeqLim->setContrasts    (              1,lNoOfContrasts,         1,             1);
        } else {
            pSeqLim->setContrasts    (              1,           1,           1,             1);
            pSeqLim->getContrasts().setDisplayMode                              ( SEQ::DM_OFF );
        }

    } else {

        if ( ( lMin =   500 ) > ( lMax = lTEMaxEasySite ) ) { lMin = lMax; }
        if ( ( lDef = 10000 ) >   lMax )                    { lDef = lMax; }
        pSeqLim->setTE            (    0,        lMin,        lMax,          10,          lDef);

        pSeqLim->setContrasts        (              1,           1,           1,             1);
        pSeqLim->getContrasts().setDisplayMode                                  ( SEQ::DM_OFF );

    }



    pSeqLim->setTI             (       0,        1000,     8000000,        1000,        300000);
    pSeqLim->setTR             (       0,        1000,    10000000,         100,        100000);
    pSeqLim->setTD             (       0,           0,    30000000,        1000,             0);


    // * ------------------------------------------------------------------------------------- *
    // * Slices/3D-Slabs         (        min,          max,         inc,         def, psLIM );*
    // * ------------------------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            pSeqLim->setSlices            (           1,K_NO_SLI_MAX,           1,           1);
        break;

        case RT_Gre:
            pSeqLim->setSlices            (           1,           1,           1,           1);
            pSeqLim->setTE       (       0,         500,       60000,          10,        5000);
            pSeqLim->setTR       (       0,        1000,    10000000,         100,       10000);
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }

    pSeqLim->setSliceThickness            (       1.000,     250.000,       0.500,       5.000);
    pSeqLim->setConcatenations            (           1,K_NO_SLI_MAX,           1,           1);
    if ( ( dMin = 5.0) > ( dMax = theGCProxy.getFoVMax() ) ) { dMin = dMax; }
    pSeqLim->setSlabThickness                                               (dMin,        dMax);
    pSeqLim->setSliceDistanceFactor       (       0.000,       8.000,       0.010,       0.200);
    pSeqLim->set3DPartThickness           (       0.050,       5.000,       0.010,       2.000);
    pSeqLim->setMinSliceResolution                                                    ( 0.500 );


    // * ------------------------------------------------------------------------------------- *
    // * In-plane FOV                                                                          *
    // * ------------------------------------------------------------------------------------- *
    if ( ( dMin = 20.0 ) > ( dMax = theGCProxy.getFoVMax() ) ) { dMin = dMax; }
    dInc =    1.0;
    dDef =  300.0;
    dMin = fSDSdRoundUpMinimum(dMin,dMax,dInc);
    if (dDef < dMin) dDef = dMin;
    if (dDef > dMax) dDef = dMax;
    pSeqLim->setReadoutFOV                (        dMin,        dMax,        dInc,        dDef);
    pSeqLim->setPhaseFOV                  (        dMin,        dMax,        dInc,        dDef);



    // * ------------------------------------------------------------------------------------- *
    // * Other Slice/Slab attributes                                                           *
    // * ------------------------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            pSeqLim->setMultiSliceMode              (SEQ::MSM_INTERLEAVED, SEQ::MSM_SEQUENTIAL);
        break;

        case RT_Gre:
            pSeqLim->setMultiSliceMode                                    (SEQ::MSM_SEQUENTIAL);
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }

    pSeqLim->setSliceSeriesMode             (SEQ::INTERLEAVED, SEQ::ASCENDING, SEQ::DESCENDING);
    pSeqLim->enableSliceShift                                                                ();
    pSeqLim->enableMSMA                                                                      ();
    pSeqLim->enableOffcenter                                                                 ();
    pSeqLim->setAllowedSliceOrientation                                   (SEQ::DOUBLE_OBLIQUE);

    // * ------------------------------------------------------------------------------------- *
    // * RF                                          (     min,       max,      inc,      def);*
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setFlipAngle                             (    1.000,   90.000,    1.000,   25.000);
    pSeqLim->setExtSrfFilename                                          ("%MEASDAT%/extrf.dat");
    pSeqLim->setRFSpoiling                                                  (SEQ::ON, SEQ::OFF);
    pSeqLim->setRFPulseType                                      (SEQ::RF_NORMAL, SEQ::RF_FAST);
    pSeqLim->setExcitationPulse                               (SEQ::EXCITATION_SLICE_SELECTIVE,
                                                            SEQ::EXCITATION_VOLUME_SELECTIVE);



    // * ------------------------------------------------------------------------------------- *
    // * Loop control                                                                          *
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setIntro                                                       (SEQ::ON, SEQ::OFF);
    pSeqLim->setAveragingMode                                (SEQ::INNER_LOOP, SEQ::OUTER_LOOP);
    pSeqLim->setEllipticalScanning                                          (SEQ::OFF, SEQ::ON);

    // * ------------------------------------------------------------------------------------- *
    // * Preparation Pulses                      (       min,       max,       inc,       def);*
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setRSatMode                                     ( SEQ::RSAT_REG, SEQ::RSAT_QUICK );
    pSeqLim->setRSats                             (         0,lNoOfRSats,         1,         0);
    pSeqLim->setRSatThickness                     (     3.000,   150.000,     1.000,    50.000);

    if ( eSequence == RT_Gre ) {
        pSeqLim->setTSats                         (         0,         0,         1,         0);
    } else {
        pSeqLim->setTSats                         (         0,         1,         1,         0);
        pSeqLim->setTSatThickness                 (     3.000,   150.000,     1.000,    40.000);
        pSeqLim->setTSatGapToSlice                (     0.000,    50.000,     1.000,    10.000);
        pSeqLim->setPSatMode      ( SEQ::PSAT_NONE, SEQ::PSAT_SINGLE_REG, SEQ::PSAT_DOUBLE_REG,
                                                SEQ::PSAT_SINGLE_QUICK, SEQ::PSAT_DOUBLE_QUICK);
        pSeqLim->setPSatThickness                 (     3.000,   150.000,     1.000,    50.000);
        pSeqLim->setPSatGapToSlice                (     5.000,    50.000,     0.100,    10.000);
    }

    if (theMSUProxy.getNominalB0() < LOW_FIELD_B0) {
        pSeqLim->setFatSuppression           ( SEQ::FAT_SUPPRESSION_OFF, SEQ::WATER_EXCITATION);
        pSeqLim->setWaterSuppression                               (SEQ::WATER_SUPPRESSION_OFF);
        pSeqLim->getWaterSuppression().setDisplayMode                           ( SEQ::DM_OFF );
    } else {
        pSeqLim->setFatSuppression             ( SEQ::FAT_SUPPRESSION_OFF, SEQ::FAT_SATURATION,
                                         SEQ::FAT_SATURATION_QUICK, SEQ::WATER_EXCITATION_FAST
                                                                       , SEQ::WATER_EXCITATION);
        pSeqLim->setWaterSuppression        (SEQ::WATER_SUPPRESSION_OFF, SEQ::WATER_SATURATION,
                                              SEQ::WATER_SATURATION_QUICK, SEQ::FAT_EXCITATION);
        pSeqLim->setFatSatMode                         (SEQ::FAT_SAT_STRONG, SEQ::FAT_SAT_WEAK);
    }



    pSeqLim->setMTC                                                         (SEQ::OFF, SEQ::ON);
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            pSeqLim->setInversion                   (SEQ::INVERSION_OFF, SEQ::VOLUME_SELECTIVE,
                                                                          SEQ::SLICE_SELECTIVE);

            pSeqLim->setSaturationRecovery      (SEQ::SATREC_NONE, SEQ::SATREC_SLICE_SELECTIVE,
                                                                  SEQ::SATREC_VOLUME_SELECTIVE);
            pSeqLim->setDarkBlood                                           (SEQ::OFF, SEQ::ON);
            pSeqLim->setTagging                      (SEQ::TAGGING_NONE, SEQ::TAGGING_GRID_TAG,
                                                                         SEQ::TAGGING_LINE_TAG);
            pSeqLim->setGridTagDistance           (     4.000,    20.000,     1.000,     8.000);
            pSeqLim->setLineTagDistance           (     4.000,   150.000,     1.000,    20.000);
        break;

        case RT_Gre:
            pSeqLim->setInversion                                          (SEQ::INVERSION_OFF);
            pSeqLim->setSaturationRecovery                                   (SEQ::SATREC_NONE);
            pSeqLim->setDarkBlood                                                    (SEQ::OFF);
            pSeqLim->setTagging                                             (SEQ::TAGGING_NONE);
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }


    // * ------------------------------------------------------------------------------------- *
    // * Acquisitions/Repetitions            (     min,       max,      inc,      def, psLIM );*
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setAverages                              (        1,       32,        1,        1);
    lInc = 100000;
    lMax = fSDSRoundToInc(2100000000-lInc+1,lInc);
    pSeqLim->setRepetitionsDelayTime                  (        0,     lMax,     lInc,        0);


    // * ------------------------------------------------------------------------------------- *
    // * Pause After Measurement             (     min,       max,      inc,      def, psLIM );*
    // * ------------------------------------------------------------------------------------- *
    if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
        // allow pause after measurement for TRIO
        if ( fabs( theMSUProxy.getNominalB0() - 2.9) <= 0.1 )  {  // * 3T system *
            pSeqLim->setMeasPause        (       0,   180000000,     1000000,           0);
        }

    }



    // * ------------------------------------------------------------------------------------- *
    // * Physiologic measurements (       min,          max,         inc,         def, psLIM );*
    // * ------------------------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            pSeqLim->setPhases            (          1, K_NO_SLI_MAX,           1,           1);

            pSeqLim->setCardiacScanWindowDialog                             (SEQ::ON, SEQ::OFF);
            pSeqLim->addPhysioMode            (SEQ::SIGNAL_NONE,        SEQ::METHOD_NONE      );
            pSeqLim->addPhysioMode            (SEQ::SIGNAL_CARDIAC,     SEQ::METHOD_TRIGGERING);
            pSeqLim->addPhysioMode            (SEQ::SIGNAL_RESPIRATION, SEQ::METHOD_TRIGGERING);
            pSeqLim->enableCaptureCycle                                                      ();
            pSeqLim->setRepetitions                   (        0,      511,        1,        0);
        break;

        case RT_Gre:
            pSeqLim->addPhysioMode                (SEQ::SIGNAL_NONE,    SEQ::METHOD_NONE      );
            // To speed up the sequence unit test, on debug versions only a small number of measuments
            //  is the default
            #ifdef DEBUG
                pSeqLim->setRepetitions               (        0,    31999,        1,       20);
            #else
                pSeqLim->setRepetitions               (        0,    31999,        1,     1999);
            #endif
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }



    // * ------------------------------------------------------------------------------------- *
    // * Data Receive & Image calculation                                                      *
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setReconstructionMode             (SEQ::RECONMODE_MAGNITUDE, SEQ::RECONMODE_PHASE,
                                                                     SEQ::RECONMODE_MAGN_PHASE);
    pSeqLim->set2DInterpolation                                             (SEQ::OFF, SEQ::ON);
    pSeqLim->setSWI                                                         (SEQ::OFF, SEQ::ON);

    pSeqLim->setDistortionCorrMode (SEQ::DISTCORR_NDIS, SEQ::DISTCORR_DIS2D, SEQ::DISTCORR_DIS3D);


    // * ------------------------------------------------------------------------------------- *
    // * Phase stabilization                                                                   *
    // * ------------------------------------------------------------------------------------- *
    switch ( eSequence ) {
        case Gre:
        case Gre_HR:
            if (theMSUProxy.getNominalB0() < LOW_FIELD_B0) {
                pSeqLim->setPhaseStabilisation                              (SEQ::ON, SEQ::OFF);
            } else {
                pSeqLim->setPhaseStabilisation                              (SEQ::OFF, SEQ::ON);
            }
        break;

        case RT_Gre:
            pSeqLim->setPhaseStabilisation                                           (SEQ::OFF);
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence type = %d  ", ptModule, eSequence);
            return ( SEQU_ERROR );
    }




    // * ------------------------------------------------------------------------------------- *
    // * Flow compensation                                                                     *
    // * ------------------------------------------------------------------------------------- *
    pSeqLim->setFlowCompensation            (  0,                    SEQ::FLOW_COMPENSATION_NO,
                                                           SEQ::FLOW_COMPENSATION_READOUT_ONLY,
                                                          SEQ::FLOW_COMPENSATION_SLICESEL_ONLY,
                                                             SEQ::FLOW_COMPENSATION_SLICE_READ,
                                                                    SEQ::FLOW_COMPENSATION_YES);
    pSeqLim->setFlowCompensation            (  1,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  2,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  3,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  4,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  5,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  6,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  7,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  8,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            (  9,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            ( 10,                    SEQ::FLOW_COMPENSATION_NO);
    pSeqLim->setFlowCompensation            ( 11,                    SEQ::FLOW_COMPENSATION_NO);




    // * ------------------------------------------------------------------------------------- *
    // * Use PAT Parameters                                                                    *
    // * ------------------------------------------------------------------------------------- *
#ifdef SUPPORT_iPAT
    if ( (eSequence == Gre) || (eSequence == Gre_HR) ) {

        // set default SeqLim for PAT: PATMode, RefScanMode, AccelFactorPE, RefLinesPE
        // if second argument is set to 'true', AccelFactor3D and RefLines3D will be initialized as well
        fPATSetDefaultSeqLim(pSeqLim,
                             true    ); // <- initialize 2D-PAT

        // allow EXTRA mode: ref.lines are scanned separately for each slice by SBBPATRefScan (within SeqLoop)
        pSeqLim->setRefScanMode      (SEQ::PAT_REF_SCAN_INPLACE,SEQ::PAT_REF_SCAN_EXTRA);
    }
#endif



    //* -------------------------------------------------------------------------------------- *
    //* Adjustment parameter                                                                   *
    //* -------------------------------------------------------------------------------------- *
    pSeqLim->setAdjShim                            (SEQ::ADJSHIM_TUNEUP, SEQ::ADJSHIM_STANDARD);


    //* -------------------------------------------------------------------------------------- *
    //* Disable the Consistency Checker that is performed after the measurement                *
    //* -------------------------------------------------------------------------------------- *
    if (eSequence == RT_Gre || eSequence == Gre_HR ) {
        pSeqLim->disableSAFEConsistencyCheck ();
    }


    // * ------------------------------------------------------------------------------------ *
    // * Database control                                                                     *
    // * ------------------------------------------------------------------------------------ *
    if ( eSequence == RT_Gre ) {
       pSeqLim->setMultipleSeriesMode (SEQ::MULTIPLE_SERIES_OFF, SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT);
    } else {
       pSeqLim->setMultipleSeriesMode (SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT, SEQ::MULTIPLE_SERIES_OFF);
    }


    //* -------------------------------------------------------------------------------------- *
    //* Disable the display of some UI parameter                                               *
    //* -------------------------------------------------------------------------------------- *
    pSeqLim->getEPIFactor().setDisplayMode                                      ( SEQ::DM_OFF );


#ifdef SUPPORT_CT
    pSeqLim->setAdjMDS(SEQ::ADJMDS_NONE,SEQ::ADJMDS_ADJUST);
    pSeqLim->getAdjMDS().setDisplayMode(SEQ::DM_OFF);
    pSeqLim->setMDSMode(SEQ::MDS_OFF, SEQ::MDS_ONCO);
    pSeqLim->setMdsRangeExtension    (0.0,         2000.0,       1.0,        500.0);
    pSeqLim->setApplicationDetails   (SEQ::APPLICATION_TIMCT);  // Disables inline composing
    pSeqLim->setMdsAllowedSliceOrientation (SEQ::MDS_TRA);
    pSeqLim->setTableSpeedNumerator(1,2,1,2);
    pSeqLim->setLinesPerSegment(1,127,1,15);
    pSeqLim->setDixon(SEQ::DIXON_NONE, SEQ::DIXON_WATER_FAT_IMAGES, SEQ::DIXON_WATER_IMAGE,SEQ::DIXON_FAT_IMAGE);
    //  Only default and availability of the following two parameters is used
    pSeqLim->setBreathHoldDuration_us(1000000,100000000,100000,15000000);
    pSeqLim->setBreathHoldRepeatedSegments(0,128,1,2);
#endif


    // * -------------------------------------------------------------------------- *
    // * Init registry entries (data rates, max. data size, ...)                    *
    // * -------------------------------------------------------------------------- *
    RunLoop.initRegistryEntries();



    #ifndef VXWORKS
        // * ---------------------------------------------------------------------- *
        // * Configure UI parameters                                                *
        // * ---------------------------------------------------------------------- *
        fStdImagingInitPost(pSeqLim);




        // * ---------------------------------------------------------------------- *
        // * Declaration of pointer to UI parameter classes                         *
        // * ---------------------------------------------------------------------- *
        LINK_BOOL_TYPE      *_phaseStabilize   = _search< LINK_BOOL_TYPE      > ( pSeqLim, MR_TAG_PHASE_STABILIZE          );
        LINK_BOOL_TYPE      *pDarkBlood        = _search< LINK_BOOL_TYPE      > ( pSeqLim, MR_TAG_DARK_BLOOD               );
        LINK_BOOL_TYPE      *pSWI              = _search< LINK_BOOL_TYPE      > ( pSeqLim, MR_TAG_SWI                      );

        LINK_LONG_TYPE      *_contrasts        = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_CONTRASTS                );
        LINK_LONG_TYPE      *_segments         = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_SEGMENTS                 );
        LINK_LONG_TYPE      *_baseResolution   = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_BASE_RESOLUTION          );
        LINK_LONG_TYPE      *_phase            = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_PHYSIO_PHASES            );
        LINK_LONG_TYPE      *_firstAcqWindow   = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_FIRST_ACQUISITION_WINDOW );
        LINK_LONG_TYPE      *_firstAcqWindowIn = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_FIRST_ACQUISITION_WINDOW_FOR_INTERNAL );
        LINK_LONG_TYPE      *pImagesPerSlab    = _search< LINK_LONG_TYPE      > ( pSeqLim, MR_TAG_IMAGES_PER_SLAB          );
        LINK_LONG_TYPE      *pSGSize           = _searchElm< LINK_LONG_TYPE   > ( pSeqLim, MR_TAG_SLICE_GROUP_LIST, MR_TAG_SG_SIZE );

        sprintf ( tMrTag, "%s.0\0", MR_TAG_TE );
        LINK_DOUBLE_TYPE    *_te               = _search< LINK_DOUBLE_TYPE    > ( pSeqLim, tMrTag                   );
        sprintf ( tMrTag, "%s.0\0", MR_TAG_TR );
        LINK_DOUBLE_TYPE    *_tr               = _search< LINK_DOUBLE_TYPE    > ( pSeqLim, tMrTag                   );
        sprintf ( tMrTag, "%s.0\0", MR_TAG_BANDWIDTH );
        LINK_DOUBLE_TYPE    *_bandwidth        = _search< LINK_DOUBLE_TYPE    > ( pSeqLim, tMrTag                   );
        LINK_DOUBLE_TYPE    *pSliceOS          = _search< LINK_DOUBLE_TYPE    > ( pSeqLim, MR_TAG_SLICE_OVERSAMPLING);
        LINK_DOUBLE_TYPE    *pFovRead          = _search< LINK_DOUBLE_TYPE    > ( pSeqLim, MR_TAG_READOUT_FOV       );
        LINK_DOUBLE_TYPE*    pSliceThickness   = _search< LINK_DOUBLE_TYPE    > (pSeqLim, MR_TAG_SLICE_THICKNESS    );
        LINK_DOUBLE_TYPE*    pFlipAngle        = _search< LINK_DOUBLE_TYPE    > ( pSeqLim, MR_TAG_FLIP_ANGLE        );

        LINK_SELECTION_TYPE *_fatSuppres       = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_FAT_SUPPRESSION   );
        LINK_SELECTION_TYPE *_waterSuppression = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_WATER_SUPPRESSION );
        LINK_SELECTION_TYPE *_gradMode         = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_GRADIENT_MODE     );
        LINK_SELECTION_TYPE *_rfMode           = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_RFPULSE_TYPE      );
        LINK_SELECTION_TYPE *_AsymEcho         = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_ASYMMETRIC_ECHO   );
        LINK_SELECTION_TYPE *_inversionMode    = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_INVERSION         );
        LINK_SELECTION_TYPE *_excitationMode   = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_EXCIT_PULSE       );
        LINK_SELECTION_TYPE *_triggerMode      = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_FIRST_SIGNAL_MODE );
        LINK_SELECTION_TYPE *_msmMode          = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_MULTI_SLICE_MODE  );
        LINK_SELECTION_TYPE *pFlowCompMode     = _searchElm< LINK_SELECTION_TYPE > (pSeqLim, MR_TAG_FLOW_COMP       );
        LINK_SELECTION_TYPE *pFatSatMode       = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_FAT_SAT_MODE      );
        LINK_SELECTION_TYPE *pDimension        = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_DIMENSION         );
        LINK_SELECTION_TYPE *pReadOutMode      = _search< LINK_SELECTION_TYPE > ( pSeqLim, MR_TAG_READOUT_MODE      );


        // * ---------------------------------------------------------------------- *
        // * Registration of set functions                                          *
        // * ---------------------------------------------------------------------- *
        if ( _tr               ) {  pOriginalTRSetFct                    = _tr->registerSetValueHandler               ( trSetValueHandler                    );  }
        if ( _contrasts        ) {  pOriginalContrastSetFct              = _contrasts->registerSetValueHandler        ( contrastSetValueHandler              );  }
        if ( _segments         ) {  pOriginalSegmentsSetFct              = _segments->registerSetValueHandler         ( segmentsSetValueHandler              );  }
        if ( _firstAcqWindow   ) {  pOriginalFirstAcqWindowSetFct        = _firstAcqWindow->registerSetValueHandler   ( firstAcqWindowSetValueHandler        );  }
        if ( _firstAcqWindowIn ) {  pOriginalFirstAcqWindowInSetFct      = _firstAcqWindowIn->registerSetValueHandler ( firstAcqWindowInSetValueHandler      );  }
        if ( _triggerMode      ) {  pOriginalFirstPhysioSignalModeSetFct = _triggerMode->registerSetValueHandler      ( firstPhysioSignalModeSetValueHandler );  }
        if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
            if ( _te           ) {  _te->registerSetValueHandler                                                      ( teSetValueHandler                    );  }
        }



        // * ---------------------------------------------------------------------- *
        // * Registration of get limits handler                                     *
        // * ---------------------------------------------------------------------- *
        if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
            if ( _tr             ) { pOriginalTRGetLimtsFct              = _tr->registerGetLimitsHandler             ( trGetLimitsHandler             ); }
            if ( _te             ) { pOriginalTEGetLimtsFct              = _te->registerGetLimitsHandler             ( teGetLimitsHandler             ); }
            if ( _firstAcqWindow ) { pOriginalFirstAcqWindowGetLimitsFct = _firstAcqWindow->registerGetLimitsHandler ( firstAcqWindowGetLimitsHandler ); }
            if ( _firstAcqWindowIn){ pOriginalFirstAcqWindowInGetLimitsFct = _firstAcqWindowIn->registerGetLimitsHandler ( firstAcqWindowInGetLimitsHandler ); }
            if ( _segments       ) { pOriginalSegmentsGetLimitsFct       = _segments->registerGetLimitsHandler       ( segmentsGetLimitsHandler       ); }
        }
        if ( pSliceOS            ) {  pOriginalSliceOSGetLimtsHandler    = pSliceOS->registerGetLimitsHandler        ( sliceOSGetLimitsHandler        ); }
        if ( pImagesPerSlab      ) {  pOriginalImaPerSlabGetLimitsHandler= pImagesPerSlab->registerGetLimitsHandler ( imaPerSlabGetLimitsHandler ); }
        if ( pFlipAngle          ) {  pOriginalFlipAngleGetLimitsHandler = pFlipAngle->registerGetLimitsHandler     ( flipAngleGetLimitsHandler  ); }



        // * ---------------------------------------------------------------------- *
        // * Registration of solve handler methods                                  *
        // * ---------------------------------------------------------------------- *
        if ( _segments         ) {  _segments->registerSolveHandler         ( fLSolveSegmentsConflict         ); }
        if ( _fatSuppres       ) {  _fatSuppres->registerSolveHandler       ( fSSolveSelectionConflict        ); }
        if ( _waterSuppression ) {  _waterSuppression->registerSolveHandler ( fSSolveSelectionConflict        ); }
        if ( _gradMode         ) {  _gradMode->registerSolveHandler         ( fSSolveSelectionConflict        ); }
        if ( _rfMode           ) {  _rfMode->registerSolveHandler           ( fSSolveSelectionConflict        ); }
        if ( _AsymEcho         ) {  _AsymEcho->registerSolveHandler         ( fSSolveSelectionConflict        ); }
        if ( _inversionMode    ) {  _inversionMode->registerSolveHandler    ( fSSolveInversionConflict        ); }
        if ( _excitationMode   ) {  _excitationMode->registerSolveHandler   ( fSSolveExcitationModeConflict   ); }
        if ( _baseResolution   ) {  pOrgSolve_baseRes =
                                   _baseResolution->registerSolveHandler    ( fDSolveBaseResConflict          ); }
        if ( _contrasts        ) {  _contrasts->registerSolveHandler        ( fSSolveLongConflict             ); }
        if ( _triggerMode      ) { _triggerMode->registerSolveHandler       ( fSSolveFirstSignalModeConflict  ); }
        if ( _msmMode          ) { _msmMode->registerSolveHandler           ( fSSolveMSMModeConflict          ); }
        if ( _phaseStabilize   ) { _phaseStabilize->registerSolveHandler    ( fSSolvePhaseStabilizeConflict   ); }
        if ( pFlowCompMode     ) { pFlowCompMode->registerSolveHandler      ( fSSolveFlowCompConfilct         ); }
        if ( pFatSatMode       ) { pFatSatMode->registerSolveHandler        ( fSSolveSelectionConflict        ); }
        if ( pDimension        ) { pOriginalDimensionSolveHandler =
                                   pDimension ->registerSolveHandler        ( fSSolveDimensionConflict        ); }
        if ( pDarkBlood        ) { pDarkBlood->registerSolveHandler         ( fSSolveDarkBloodConflict        ); }
        if ( pReadOutMode      ) { pReadOutMode->registerSolveHandler       ( fSSolveReadOutModeConflict      ); }
        if ( pSliceOS          ) { pSliceOS->registerSolveHandler           ( fDSolveSliceOSConflict          ); }
        if ( _te               ) { _te->registerSolveHandler                ( fLSolveTEConflict               ); }
        if ( _bandwidth        ) {  _bandwidth->registerSolveHandler        ( fSolveDoubleMultiTEConflict     ); }
        if ( pFovRead          ) { pFovRead->registerSolveHandler           ( fSolveDoubleMultiTEConflict     ); }
        if ( pSliceThickness   ) { pFovRead->registerSolveHandler           ( fSolveDoubleMultiTEConflict     ); }
        if ( pSGSize           ) { pOriginalSGSizeSolveHandler = pSGSize->registerSolveHandler (fLSolveSGSizeConflict); }
        if ( pSWI              ) { pSWI->registerSolveHandler               ( fBSolveSWIConflict              ); }
//        if ( pFlipAngle        ) { pFlipAngle->registerSolveHandler         ( fDSolveFlipAngle                );  }



        // * ---------------------------------------------------------------------- *
        // * Registration of get-options handler methods                            *
        // * ---------------------------------------------------------------------- *
        if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
            if ( pFatSatMode  ) {  pFatSatMode->registerGetOptionsHandler  ( fatSatModeGetOptionsHandler     ); }
        }






        // * ---------------------------------------------------------------------- *
        // * Registration of is available handler                                   *
        // * ---------------------------------------------------------------------- *
        if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
            if ( _phase       ) {  pOriginalPhasesIsAvailableFct = _phase->registerIsAvailableHandler ( phasesIsAvailableHandler );  }
            if ( pFatSatMode  ) {  pFatSatMode->registerIsAvailableHandler ( fatSatModeIsAvailable );  }
        }



        // * ---------------------------------------------------------------------- *
        // * Registration of get limits handler                                     *
        // * ---------------------------------------------------------------------- *
        if ( (eSequence == Gre) || (eSequence == Gre_HR) )  {
            if ( _firstAcqWindow ) {  pOriginalFirstAcqWindowToolTipIdHandler = _firstAcqWindow->registerGetToolTipIdHandler ( firstAcqWindowToolTipIdHandler ); }
        }



        if ( eSequence == Gre_HR )  {
            LINK_DOUBLE_TYPE* pMaxGradAmpl = _search< LINK_DOUBLE_TYPE    > (pSeqLim, MR_TAG_ECHO_SPACING    );

            if ( pMaxGradAmpl ) { pMaxGradAmpl->registerIsAvailableHandler ( fDIsAvailableMaxGradAmpl ); }
            if ( pMaxGradAmpl ) { pMaxGradAmpl->registerGetValueHandler    ( fDGetValueMaxGradAmpl    ); }
            if ( pMaxGradAmpl ) { pMaxGradAmpl->registerSetValueHandler    ( fDSetValueMaxGradAmpl    ); }
            if ( pMaxGradAmpl ) { pMaxGradAmpl->registerGetUnitIdHandler   ( fDGetUnitId              ); }
            if ( pMaxGradAmpl ) { pMaxGradAmpl->registerGetLimitsHandler   ( fDGetLimitsHandler       ); }
            if ( pMaxGradAmpl ) { pMaxGradAmpl->registerGetLabelIdHandler  ( fDGetLabelIdGetMaxGradAmpl ); }
        }

#ifdef SUPPORT_CT
        if( LINK_LONG_TYPE* pCT = _create<LINK_LONG_TYPE>(pSeqLim,"ct_sms") )
        {
            pCT->registerGetValueHandler(fCTNSegGetVal);
        }
        if( LINK_SELECTION_TYPE* pSel = _search<LINK_SELECTION_TYPE>(pSeqLim,MR_TAG_MDS_MODE) )
        {
            pSel->registerSolveHandler(fSSolveSelectionConflict);
        }
        if( LINK_SELECTION_TYPE* pSel = _search<LINK_SELECTION_TYPE>(pSeqLim,MR_TAG_DIXON) )
        {
            pOriginalDixonSetFct = pSel->registerSetValueHandler(fCTDixonSetValue);
            pSel->registerSolveHandler(fCTSolveDixonConflict);
            pSel->registerIsAvailableHandler(fCTDixonIsAvailable);
        }
        if( LINK_SELECTION_TYPE* pSel = _search<LINK_SELECTION_TYPE>(pSeqLim,MR_TAG_MDS_MODE) )
        {
            pSel->registerSolveHandler(fCTSolveMDSModeConflict);
            pOriginalMDSModeSetFct = pSel->registerSetValueHandler(fCTMDSModeSetValue);
        }
        lStatus = CT::fCTInitUI(pSeqLim);
        if( (lStatus & NLS_SEV) != NLS_SUCCESS )
        {
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
            return lStatus;
        }

#endif // SUPPORT_CT defined

    #endif // VXWORKS not defined


    // ---------------------------------------------------------------------------------
    // file containing the default postprocessing protocol (EVAProtocol)
    // also enables the parametric subcard to be used
    //----------------------------------------------------------------------------------
    #ifndef VXWORKS
        if ( eSequence == Gre) {
            pSeqLim->setParametricMapping(SEQ::PMAP_NONE, SEQ::PMAP_T2STAR_MAP);         //make the parametric sub-card available
            pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\Breast\\Breast.evp")); // this is the file which will include the EVA protocol
        } else {
            pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\Breast\\Breast.evp")); // this is the file which will include the EVA protocol
        }
    #endif

    #ifdef SUPPORT_PACE
        // * ------------------------------------------------------------------------ *
        // * Initialization of SeqLoopBH.                                             *
        // * ------------------------------------------------------------------------ *

        if( !RunLoop.init(pSeqLim) ) {
            TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: SeqLoopBH::init failed.", ptModule);
            return RunLoop.getNLSStatus();
        }

        // * ------------------------------------------------------------------------ *
        // * Enable SeqLoopBH to shift the sat positions at run time.                 *
        // * ------------------------------------------------------------------------ *

        if( !RunLoop.adaptRSatPos(SBBRSat, lNoOfRSats) )
        {
            TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: SeqLoopBH::adaptRSatPos failed.", ptModule);
            return RunLoop.getNLSStatus();
        }
    #endif

    #ifndef VXWORKS
        #ifdef SUPPORT_iPAT


            //-------------------------------------------------------------------------------------
            // register SetValue-handlers for segmented k-space sampling with PAT
            //-------------------------------------------------------------------------------------
            if (!fPATRegisterUILinkHandlersForSegmentedSequences(pSeqLim, &Reorder, SEQ::INC_GRE_SEGMENTS, lPATDefOptNoRefLines))
            {
                textTR("fPATRegisterUILinkHandlersForSegmentedSequences failed!");
                return SEQU_ERROR;
            }
        #endif
    #endif // VXWORKS

    #ifndef VXWORKS
    // prep for receiving transforms
    transform_receiver.initNetwork();
    #endif

    return(lStatus);
}

/*[ Function ****************************************************************\
*
* Name        : fSEQPrep
*
* Description : Prepares everything that the sequence needs at run time.
*
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQPrep
(
  MrProt       *pMrProt,       // * IMP: Measurement protocol  *
  SeqLim       *pSeqLim,       // * IMP: Sequence limits       *
  SeqExpo      *pSeqExpo       // * EXP: Returned values       *
)

{
    static const char *ptModule             = {"fSEQPrep"};


    const long   lMaxTimeBetweenQuickFatSatsWeak_us     = 75000;
    const long   lMaxTimeBetweenQuickFatSatsStrong_us   = 50000;


    NLS_STATUS   lStatus = SEQU__NORMAL;                                // * My return status            *
    long         lI;                                                    // * Helper variables            *
    #ifdef GSWD
        long     lTEInc, lTRInc, lTRNeeded;
        long     alTENeeded [lNoOfContrasts];
    #endif
    long         lPhasesToMeasure              = 1;                     // * Number of phases to measure *

    long         lTRMin                        = 0;                     // * Minimum TR                  *

    long         lNoOfMeasurements             = 0;                     // * Number of meas. repeats     *
    long         lSegments                     = 0;                     // * Number of segments from the *
                                                                        // * protocol                    *
    long         lSegmentsInTR                 = 0;                     // * Number of segments that are *
                                                                        // * executed with TR            *
    long         lScanTimeSatsEachSegment      = 0;                     // * Total scan time for sats    *
                                                                        // * for each segment            *
    long         lScanTimeSatsOnlyFirstSegment = 0;                     // * Scan time for sats that are *
                                                                        // * only executed for the first *
                                                                        // * segment                     *
    long         lMaxTimeBetweenQuickFatSats_us= 0;                     // * Max. allowed time between   *
                                                                        // * two quick fat sat pulses    *
    long         lQuickFatSatTimeInsideTR      = 0;                     // * Duration of quick CSats     *
    long         lQuickWaterSatTimeInsideTR    = 0;                     // * that are not included in TR *
    long         lQuickRSatTimeInsideTR        = 0;                     // * Duration of quick RSats     *
                                                                        // * that are not included in TR *
    long         lQuickSpoilGradTimeInsideTR   = 0;                     // * Duration of spoiler         *
                                                                        // * gradients that are not      *
                                                                        // * included in TR              *
    long         lNegativeFillTime             = 0;                     // * Binary search flag          *
    double       dRFEnergyInSBBs               = 0.0;                   // * RF energy in SBB calls      *
    double       dRFEnergyInKernel             = 0.0;                   // * RF energy of the GRE kernel *
    double       dRFEnergyInPrepPulse          = 0.0;                   // * RF energy in preparing      *
    double       dRFEnergyInDummyScans         = 0.0;
                                                                        // * pulses                      *
    double       dMeasureTimeUsec              = 0.0;                   // * Measurement time (usec)     *
    double       dTotalMeasureTimeUsec         = 0.0;                   // * Total measurement time      *


    long         lMaxNoOfSlicesInConcat        = 1;
    long         lNoOfQuickFatSatsPerConcat    = 0;
    long         lNoOfQuickWaterSatsPerConcat  = 0;
    long         lNoOfQuickRSatsPerConcat      = 0;


    long         lRequestsLineGridPulses       = 0;                     // * No. of line/grid tag pulses *
    long         lRequestsTSatPulses           = 0;                     // * No. of TSat pulses          *
    long         lRequestsCRMSatPulses         = 0;                     // * No. of C-, R- or MSat pulses*

    bool         bSuccess;


    char         ptIdentdummy[7];                                       // * Ident strings for sat pulses *


    int indi;
    int indj;


    // * ---------------------------------------------------------------------- *
    // * Required to get the actual B0 field strength                           *
    // * ---------------------------------------------------------------------- *
    static MSUProxy  theMSU;



    // * ---------------------------------------------------------------------- *
    // * For Gre_HR only:                                                       *
    // *                                                                        *
    // * Set default for maximum readout gradient amplitude                     *
    // * ---------------------------------------------------------------------- *
    if ( eSequence == Gre_HR )  {

        if ( ! pMrProt->fastImaging().echoSpacing() )  {

            if ( pSeqLim->isContextPrepForMrProtUpdate() )  {
                static GPAProxy theGPA;
                pMrProt->fastImaging().echoSpacing (theGPA.getGradMaxAmplAbsolute());
            } else {
                return ( SEQU_ERROR );
            }

        }

    }



    // * ---------------------------------------------------------------------- *
    // * we need to know the selected coil element for the reordering later on  *
    // * ---------------------------------------------------------------------- *
    const MrCoilSelect &coilSelect = pMrProt->coilInfo().Meas();
    SelectedCoilElements mySelectedElements(coilSelect, pSeqLim->getCoilContext());



    // * ---------------------------------------------------------------------- *
    // * make breast subcard available                                          *
    // * ---------------------------------------------------------------------- *
    if ( (eSequence == Gre) || (eSequence == Gre_HR) ) {
        pMrProt->applicationDetails(SEQ::APPLICATION_INLINE_BREAST);  // make breast subcard available
    }


    // * ---------------------------------------------------------------------- *
    // * Initialization                                                         *
    // * ---------------------------------------------------------------------- *
    lI = 0;
    lPhasesToMeasure = maximum (1L, pMrProt->physiology().phases()    );
    lSegments        = maximum (1L, pMrProt->fastImaging().segments() );

    // * Segments handling for various multi slice modes: *
    // *    sequential : segments are within TR           *
    // *    interleaved: segments are NOT within TR       *
    switch ( pMrProt->kSpace().multiSliceMode() ) {
        case SEQ::MSM_SEQUENTIAL:
                lSegmentsInTR = lSegments;
        break;

        case SEQ::MSM_INTERLEAVED:
                lSegmentsInTR = 1;
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence multi slice mode = %d  ", ptModule, pMrProt->kSpace().multiSliceMode());
            return ( SEQU_ERROR );

    }

    lNoOfMeasurements                 = pMrProt->repetitions() + 1;

    lNumberOfQuickPRSats              = 0;
    lNumberOfRegularPRSats            = 0;
    lNoOfQuickSpoilGradsPerConcat     = 0;



    // * ------------------------------------------------------------------------- *
    // * Determine trigger mode and signal                                         *
    // * ------------------------------------------------------------------------- *
    pMrProt->physiology().getPhysioMode (PhysioSignalHigh, PhysioMethodHigh, PhysioSignalLow, PhysioMethodLow);

#ifdef SUPPORT_CT
    lStatus = SeqLoop::ct_try_prot(pMrProt,pSeqLim->isContextNormal() != 0);
    if( (lStatus & NLS_SEV) != NLS_SUCCESS )
    {
        if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
        return lStatus;
    }
    if( (pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO) && (pMrProt->fastImaging().segments() > 1) )
    {
        if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
        return SEQU_ERROR;
    }
    if( pMrProt->dixon() != SEQ::DIXON_NONE )
    {
        //  2 contrasts
        if( pMrProt->contrasts() != 2 )
        {
            if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
            return SEQU_ERROR;
        }
        //  Adaptive combine
        if( pMrProt->coilInfo().Meas().getNumOfUsedRxChan() > 1 && pMrProt->coilCombineMode() != SEQ::COILCOMBINE_ADAPTIVE_COMBINE )
        {
            if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
            return SEQU_ERROR;
        }
        //  Only allowed, if TimCT-Onco is used
        if( pMrProt->MDS().mdsModeMask() == SEQ::MDS_OFF )
        {
            if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
            return SEQU_ERROR;
        }
    }
#endif


    // * ------------------------------------------------------------------------- *
    // * No PATMode 'mSENSE' with more than one contrast allowed (inplace)         *
    // * (second contrast is reconstructed with weights from first contrast,       *
    // *  which sometimes causes severe artefacts)                                 *
    // * ------------------------------------------------------------------------- *
    if ( (pMrProt->PAT().PATMode() == SEQ::PAT_MODE_SENSE) && (pMrProt->contrasts() > 1) && (pMrProt->PAT().RefScanMode()==SEQ::PAT_REF_SCAN_INPLACE) )  {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No PATMode 'mSENSE' with more than one contrast allowed.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }


    if ( (pMrProt->parametricMapping().parametricMapMode() == SEQ::PMAP_T2STAR_MAP) && (pMrProt->repetitions() > 0))
    {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No T2-star map with more than one measurement allowed.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }



    // * ------------------------------------------------------------------------- *
    // * SENSE not allowed in combination with ref.scan mode 'extra'               *
    // * (since we use the separate SBBPATRefScan, SENSE will not work)            *
    // * ------------------------------------------------------------------------- *
    if ( (pMrProt->PAT().PATMode() == SEQ::PAT_MODE_SENSE) && (pMrProt->PAT().RefScanMode()==SEQ::PAT_REF_SCAN_EXTRA) )  {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "PATMode 'mSENSE' not allowed in combination with ref.scan mode 'separate'.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }




    // * ------------------------------------------------------------------------- *
    // * No inversion or darkblood preparation without segmentation                *
    // * ------------------------------------------------------------------------- *
    if ( (pMrProt->fastImaging().segments() <= 1) &&
         (    (pMrProt->preparationPulses().inversion() != SEQ::INVERSION_OFF)
           || (pMrProt->preparationPulses().satRecovery() != SEQ::SATREC_NONE)
           || pMrProt->preparationPulses().darkBlood() ) )  {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No inversion or darkblood without segmentation.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }




    // * ------------------------------------------------------------------------- *
    // * Check consistency of preparation pulses                                   *
    // * ------------------------------------------------------------------------- *
    if ( (pMrProt->preparationPulses().satRecovery() != SEQ::SATREC_NONE) && pMrProt->preparationPulses().darkBlood() )  {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No saturation recovery pulse in combination with dark blood.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }




    // * ------------------------------------------------------------------------- *
    // * No real part reconstruction without inversion pulse                       *
    // * ------------------------------------------------------------------------- *
    if ( pMrProt->calcRealPartImages() && ( pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF ) ) {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No real part reconstruction without inversion pulse.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }





    // * ------------------------------------------------------------------------- *
    // * TSat is only allowed for sequential multi slice excitation                *
    // * ------------------------------------------------------------------------- *
    if ( (pMrProt->kSpace().multiSliceMode() == SEQ::MSM_INTERLEAVED) && (pMrProt->tSat().on()) )  {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (lStatus = SEQU_ERROR, "TSat is only allowed with sequential multi slice.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }


    // * ------------------------------------------------------------------------- *
    // * PACE respiratory triggering does not support 3D imaging        *
    // * ------------------------------------------------------------------------- *
    if ( (pMrProt->NavigatorParam().RespComp() & (SEQ::RESP_COMP_TRIGGER|SEQ::RESP_COMP_TRIGGER_AND_FOLLOW)) != 0
         && (pMrProt->kSpace().dimension() == SEQ::DIM_3) ) {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (lStatus = SEQU_ERROR, "Respiratory triggering is only allowed for 2D.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }




    // * ------------------------------------------------------------------------- *
    // * Segmentation is only allowed for sequential multi slice excitation        *
    // * This requirement does not hold for respiratory triggering.                *
    // * ------------------------------------------------------------------------- *
    if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) ) {

        if ( (pMrProt->kSpace().multiSliceMode() == SEQ::MSM_SEQUENTIAL)  )  {

            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                CheckStatusPB (lStatus = SEQU_ERROR, "Segmentation is only allowed with interleaved multi slice.");
            } else  {
                CheckStatusB  (SEQU_ERROR);
            }

        }

    } else {

        if ( (pMrProt->kSpace().multiSliceMode() == SEQ::MSM_INTERLEAVED) && (pMrProt->fastImaging().segments() > 1 ) )  {

            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                CheckStatusPB (lStatus = SEQU_ERROR, "Segmentation is only allowed with sequential multi slice.");
            } else  {
                CheckStatusB  (SEQU_ERROR);
            }

        }

    }



    // * ------------------------------------------------------------------------- *
    // * Respiratory triggering in steady-state does not support 3D imaging        *
    // * ------------------------------------------------------------------------- *
    if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) &&
         (pMrProt->kSpace().dimension() == SEQ::DIM_3) ) {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (lStatus = SEQU_ERROR, "Respiratory triggering is only allowed for 2D.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }




    // * ------------------------------------------------------------------------- *
    // * Respiratory triggering in steady-state does not supported multiple phases *
    // * ------------------------------------------------------------------------- *
    if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) &&
         (lPhasesToMeasure > 1) ) {

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (lStatus = SEQU_ERROR, "Respiratory triggering is only allowed with a single phase.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }

    }



    // * ------------------------------------------------------------------------- *
    // * With some preparation pulses only sequetial multislice mode is allowed    *
    // * ------------------------------------------------------------------------- *
    if ( pMrProt->kSpace().multiSliceMode() != SEQ::MSM_SEQUENTIAL )   {

        if( pMrProt->gridTag().on() )  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No GridTag without sequential multislice.");
            } else  {
            CheckStatusB  (SEQU_ERROR);
            }

        }

        if( pMrProt->lineTag().on() )  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No LineTag without sequential multislice.");
            } else  {
            CheckStatusB  (SEQU_ERROR);
            }
        }


        if( pMrProt->preparationPulses().inversion() == SEQ::VOLUME_SELECTIVE )  {
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SBB_ERROR, "No non-selective inversion without sequential multislice.");
            } else  {
            CheckStatusB  (SEQU_ERROR);
            }
        }

    }



    // * ---------------------------------------------------------------------- *
    // * Check consistency of sequence dimension and excitation type            *
    // * No non-selective excitation in case of2D                               *
    // * ---------------------------------------------------------------------- *
    if ( (pMrProt->kSpace().dimension() == SEQ::DIM_2) && (pMrProt->txSpec().excitation() == SEQ::EXCITATION_VOLUME_SELECTIVE) )  {
        CheckStatusB (SEQU_ERROR);
    }



    // * ---------------------------------------------------------------------- *
    // * Check consistency of the number of slices and excitation type          *
    // * ---------------------------------------------------------------------- *
    if ( (pMrProt->sliceSeries().size() != 1) && (pMrProt->txSpec().excitation() == SEQ::EXCITATION_VOLUME_SELECTIVE) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "No non-selective excitation with multiple slices/slabs.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }







    // * ---------------------------------------------------------------------- *
    // * Protocol checks for SWI:                                               *
    // *    - not supported for 2D                                              *
    // *    - requires flow compensation in all three spatial directions        *
    // *    - requires standard shim mode                                       *
    // *    - not supported in combination with SENSE                           *
    // *    - not supported in combination with "Save Uncombined"               *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->SWI() && pMrProt->kSpace().dimension() == SEQ::DIM_2 )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is only supported for 3D imaging.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && pMrProt->flowComp()[0] != SEQ::FLOW_COMPENSATION_YES )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported without flow compensation.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && pMrProt->adjustment().adjShimMode() != SEQ::ADJSHIM_STANDARD )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported without standard shim procedure.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && pMrProt->PAT().PATMode() == SEQ::PAT_MODE_SENSE )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported with SENSE.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && pMrProt->uncombImages() )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported with Save Uncombined Images.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && (pMrProt->coilCombineMode() == SEQ::COILCOMBINE_SUM_OF_SQUARES) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported with coil combine mode: sum of squares.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && (pMrProt->PAT().AccelFact3D() > 1) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported with iPAT in partitions direction.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && (! pMrProt->fastImaging().RFSpoiling()) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported without rf-spoiling.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }


    if ( pMrProt->SWI() && (pMrProt->reconstructionMode() == SEQ::RECONMODE_PHASE) )  {
        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (SEQU_ERROR, "SWI is not supported with phase image reconstruction.");
        } else  {
            CheckStatusB  (SEQU_ERROR);
        }
    }

    if (pSeqLim->isContextNormal() && transform_receiver.checkForTransform()) {
      pMrProt->sliceSeries().atPos(0).rotationAngle(
        transform_receiver.getTransformMatrixEl(0, 0));

      TRACE_PUT4(TC_ALWAYS, TF_SEQ, "%g %g %g %g",
                 transform_receiver.getTransformMatrixEl(0, 0),
                 transform_receiver.getTransformMatrixEl(0, 1),
                 transform_receiver.getTransformMatrixEl(0, 2),
                 transform_receiver.getTransformMatrixEl(0, 3)
                 );
      TRACE_PUT4(TC_ALWAYS, TF_SEQ, "%g %g %g %g",
                 transform_receiver.getTransformMatrixEl(1, 0),
                 transform_receiver.getTransformMatrixEl(1, 1),
                 transform_receiver.getTransformMatrixEl(1, 2),
                 transform_receiver.getTransformMatrixEl(1, 3)
                 );
      TRACE_PUT4(TC_ALWAYS, TF_SEQ, "%g %g %g %g",
                 transform_receiver.getTransformMatrixEl(2, 0),
                 transform_receiver.getTransformMatrixEl(2, 1),
                 transform_receiver.getTransformMatrixEl(2, 2),
                 transform_receiver.getTransformMatrixEl(2, 3)
                 );
    }

    // * ---------------------------------------------------------------------- *
    // * Calculate the rotation matrices and offsets for slices                 *
    // * ---------------------------------------------------------------------- *
    CheckStatusPB(lStatus = fSUPrepSlicePosArray (pMrProt, pSeqLim, asSLC),"fSUPrepSlicePosArray");

    if (false) {
      for (indi = 0; indi < 3; indi++) {
        for (indj = 0; indj < 3; indj++) {
          TRACE_PUT3(TC_ALWAYS, TF_SEQ, "m_sROT_MATRIX[%d][%d] = %f", indi, indj, asSLC[0].m_sROT_MATRIX.dMat[indi][indj]);
        }
      }
      //TRACE_PUT1(TC_ALWAYS, TF_SEQ, "slice array 0 angle: %g", pMrProt->sliceSeries().atPos(0).rotationAngle());
    }


    // * ---------------------------------------------------------------------- *
    // * Calculates the reordering scheme, lines/partitions to measure and the  *
    // * k-space center line/partition                                          *
    // *                                                                        *
    // * CAUTION: The reordering schemes HAVE to be calculated prior to the     *
    // *          calculation of the rf-energy and the measurement time         *
    // *          since these measures depend on the total number of scans      *
    // * ---------------------------------------------------------------------- *
    lStatus = fCalculateReordering ( pMrProt,
                                     pSeqLim,
                                    &Reorder,
                                    &lLinesToMeasureMax,
                                    &lLinesToMeasure,
                                    &lKSpaceCenterLine,
                                    &lPartitionsToMeasureMax,
                                    &lPartitionsToMeasure,
                                    &lKSpaceCenterPartition );
    CheckStatusB (lStatus);






    // * ---------------------------------------------------------------------- *
    // * Set the minimum distance factor                                        *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->preparationPulses().inversion() == SEQ::SLICE_SELECTIVE )  {
        RunLoop.setdDistFacMinIfConcYes ( 0.0 );
        RunLoop.setdDistFacMinIfConcNo  ( 0.2 );

    } else {
        RunLoop.setdDistFacMinIfConcYes ( 0.0 );
        RunLoop.setdDistFacMinIfConcNo  ( 0.0 );
    }




    // * ---------------------------------------------------------------------- *
    // * Set the number of concatenations                                       *
    // * ---------------------------------------------------------------------- *
    switch ( pMrProt->kSpace().multiSliceMode() )  {
        case SEQ::MSM_INTERLEAVED:
            RunLoop.setlConcMode ( Concatenations_NORMAL );
        break;

        case SEQ::MSM_SEQUENTIAL:
            RunLoop.setlConcMode ( Concatenations_SEQUENTIALSLICE );
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown multiSliceMode = %d  ", ptModule, pMrProt->kSpace().multiSliceMode());
            return ( SEQU_ERROR );
    }



    // * ---------------------------------------------------------------------- *
    // * Set the duration used for preparing scans                              *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->fastImaging().segments() > 1 ) {

        // !!! REMOVE WITH WHEN THE NEW KERNEL MODE FOR STEADY STATE TRIGGERING IS AVAILABLE !!!

        // * No preparation scans for steady-state respiratory triggering
        if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) ) {

            RunLoop.setlPreparingTime ( 0 );
        } else {
            RunLoop.setlPreparingTime ( 1 );
        }

    } else {
        RunLoop.setlPreparingTime ( 0 );
    }




    // * ---------------------------------------------------------------------- *
    // * Switch on phase correction scans for realpart reconstruction with      *
    // * partial fourier off, 7/8 and 6/8.                                      *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->calcRealPartImages() )  {
        RunLoop.setPhaseCorScans  ( true );
        RunLoop.setlPreparingTime (    0 );
    }






    // * ---------------------------------------------------------------------- *
    // * Select standard or steady-state triggering and define the delay that   *
    // * has to be used prior to the trigger halt event (necessary if gradient  *
    // * ramps are outside of the kernel event block).                          *
    // * ---------------------------------------------------------------------- *
    if ( PhysioMethodHigh == SEQ::METHOD_TRIGGERING ) {

        if ( PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION ) {

            RunLoop.seteTriggerMode ( SEQ::TRIGGER_STEADY_STATE );
            RunLoop.setlSteadyStateTriggerNoOfKernelsPerCheck ( 1 );
            RunLoop.setbSteadyStateTriggerDebugMode ( false );
            RunLoop.setbSteadyStateTriggerUnitTestActive ( mIsUnittestActive() != 0 );

            // * Segments are not within TR
            // *   ==> all lSegments have to be acquired without triggering
            RunLoop.setlLinesPerTrigger ( lSegments );


            // * Segments are not within TR
            // *  ==> Acq. Window length is TR * Phases * Segments
            pSeqExpo->setAcqWindowFactor  ( lSegments );

            if ( pMrProt->tr()[0] * pSeqExpo->getAcqWindowFactor() + pMrProt->physiology().triggerDelay(PhysioSignalHigh) >
                 pMrProt->physiology().scanWindow(PhysioSignalHigh) * 1000 )  {
                return ( SEQU_ERROR );
            }

        } else {

            RunLoop.seteTriggerMode ( SEQ::TRIGGER_STANDARD );

            // * Segments are within TR
            RunLoop.setlLinesPerTrigger  ( 1 );
            pSeqExpo->setAcqWindowFactor ( 1 );

        }

    } else {

        // * Reset data *
        RunLoop.setlLinesPerTrigger  ( 1 );
        pSeqExpo->setAcqWindowFactor ( 1 );

    }





    // * ---------------------------------------------------------------------- *
    // * Set the looping parameters for fSEQRunStd3dE                           *
    // * ---------------------------------------------------------------------- *
    RunLoop.setPhasesToMeasure  ( lPhasesToMeasure );

    switch ( pMrProt->kSpace().dimension() )  {
        case SEQ::DIM_2:
#ifdef SUPPORT_CT
            if( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct )
            {
                const int32_t i32NLinPerSeg = maximum(int32_t(1),pMrProt->MDS().mdsLinesPerSegment());
                const int32_t i32NSeg       = (lLinesToMeasure+i32NLinPerSeg-1)/i32NLinPerSeg;
                RunLoop.setLinesToMeasure( i32NSeg  );
                RunLoop.setPartitionsToMeasure( i32NLinPerSeg );
                RunLoop.setLoopLength3dE( 0 );
                break;
            }
#endif
            // * Segments loop is part of fSEQRunKernel ==> no. of calls is lLinesToMeasure / lSegments *
            // * Exception: respiratory triggering in steady-state                                      *
            if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) ) {
                RunLoop.setLinesToMeasure      ( lLinesToMeasure  );
            } else {
                RunLoop.setLinesToMeasure      ( lLinesToMeasure / lSegments);
            }
            RunLoop.setPartitionsToMeasure ( 1 );
            RunLoop.setLoopLength3dE       ( 0 );     // * Do not use 3dE loop structure of SEQLoop *
        break;

        case SEQ::DIM_3:
            // ellipt.scanning and 3D-PAT require the reorder3dE mode
            if ( pMrProt->kSpace().ellipticalScanning() || (   (pMrProt->PAT().PATMode() != SEQ::PAT_MODE_NONE)
                                                            && (pMrProt->PAT().AccelFact3D() > 1)               )) {
                RunLoop.setLoopLength3dE ( Reorder.getNoOfReorderIndices() );
            } else {
                RunLoop.setLinesToMeasure      ( lLinesToMeasure );
                RunLoop.setPartitionsToMeasure ( lPartitionsToMeasure / lSegments );
                RunLoop.setLoopLength3dE       ( 0 );   // * No elliptical scanning selected          *
            }
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
            return ( SEQU_ERROR );
    }




    // * ---------------------------------------------------------------------- *
    // * Disable effective TR calculation for 3D imaging or segmentation        *
    // * ---------------------------------------------------------------------- *
    if ( ( pMrProt->kSpace().dimension() == SEQ::DIM_3 ) || ( lSegments > 1 ) )  {
        RunLoop.setTReffective ( false );
    }



    // * ---------------------------------------------------------------------- *
    // * Perform TokTokTok only for the first measurement                       *
    // * ---------------------------------------------------------------------- *
    RunLoop.setPerformTokTokTok ( OnlyFirstRepetition );



    // * ---------------------------------------------------------------------- *
    // * The Osc Bits is not performed within RunLoop but within the sequence   *
    // * itself                                                                 *
    // * ---------------------------------------------------------------------- *
    RunLoop.setPerformOscBit ( false );



    // * ---------------------------------------------------------------------- *
    // * Disable TRFill in SeqLoop in case of unsegmented protocols to          *
    // * minimize the number of event blocks.                                   *
    // * The Kernel has to apply the TRFill by itself.                          *
    // * ---------------------------------------------------------------------- *
    if ( lSegmentsInTR == 1 ) {
        RunLoop.setPerformTRFill ( FALSE );
    } else {
        RunLoop.setPerformTRFill ( TRUE  );
    }



    // * ---------------------------------------------------------------------- *
    // * iPAT: switch on/off SBBPATRefScan for Ref.Scan Mode 'extra'            *
    // * ---------------------------------------------------------------------- *
    if (pMrProt->PAT().PATMode()!=SEQ::PAT_MODE_NONE)
    {
        // switch on/off SBBPATRefScan for RefScanMode 'extra'
        if (pMrProt->PAT().RefScanMode()==SEQ::PAT_REF_SCAN_EXTRA)
        {
            // advice SeqLoop to run SBBPATRefScan once (i.e. only for first measurement)
            RunLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_ONCE);
        }
        else
        {
            // inplace mode => no SBBPATRefScan from SeqLoop required
            RunLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_NEVER);
        }
    }
    else
    {
        RunLoop.setePATRefScanLoopMode(SeqLoop::PATRefScanLoopMode_NEVER);
    }



    // * ---------------------------------------------------------------------- *
    // * Prepare the run loop                                                   *
    // * ---------------------------------------------------------------------- *
    if (! RunLoop.prep(pMrProt, pSeqLim, pSeqExpo) ) return ( RunLoop.getNLSStatus() );






    // * ---------------------------------------------------------------------- *
    // * Prepare the RF pulse structures                                        *
    // * ---------------------------------------------------------------------- *
    CheckStatusB( lStatus = fSelectAndPrepRFPulse (pMrProt, pSeqLim, pSeqExpo, &pSRF ) );
    if (! GREKernel.setRFPulseForExcitation (pSRF) ) return ( GREKernel.getNLSStatus() );




    // * ---------------------------------------------------------------------- *
    // * Settings for the kernel calculation                                    *
    // * ---------------------------------------------------------------------- *
    CheckStatusB ( lStatus = fSetKernelParameter (pMrProt,
                                                  pSeqLim,
                                                  lLinesToMeasureMax,
                                                  lKSpaceCenterLine,
                                                  lPartitionsToMeasureMax,
                                                  lKSpaceCenterPartition) );



    // * ---------------------------------------------------------------------- *
    // * Calculate and prepare the GRE kernel                                   *
    // *                                                                        *
    // * using symmetric echoes at first                                        *
    // * ---------------------------------------------------------------------- *
    lStatus = fCalculateGRETiming (pMrProt, pSeqLim, pSeqExpo );
    if ( ((lStatus & NLS_SEV) == NLS_SEV) && !pSeqLim->isContextPrepForMrProtUpdate() )  {  return ( lStatus );  }


    // * ---------------------------------------------------------------------- *
    // * Define the number of requests of each RF pulse                         *
    // * ---------------------------------------------------------------------- *
    GREKernel.setRequestsPerMeasurement ( RunLoop.getKernelRequestsPerMeasurement (pSeqLim, lSegmentsInTR) );



    // * ---------------------------------------------------------------------- *
    // * Prepare the gradient echo kernel                                       *
    // * ---------------------------------------------------------------------- *
    if (! GREKernel.prep ( pMrProt, pSeqLim, pSeqExpo ) ) return ( GREKernel.getNLSStatus() );



    // * ---------------------------------------------------------------------- *
    // * Define the number of TSat pulses                                       *
    // * ---------------------------------------------------------------------- *
    lRequestsTSatPulses   = GREKernel.getRequestsPerMeasurement();




    // * ---------------------------------------------------------------------- *
    // * Define the number of C-, M-, or RSat pulses                            *
    // * Only one sat pulse for the first segment                               *
    // * At this point no quick sat scheme is considered                        *
    // * ---------------------------------------------------------------------- *
    lRequestsCRMSatPulses = GREKernel.getRequestsPerMeasurement() / lSegmentsInTR;



    // * ---------------------------------------------------------------------- *
    // * Calculate the maximum number of slices that will measured within       *
    // * one concatenation. This value is required for quick sat schemes.       *
    // * ---------------------------------------------------------------------- *
    lMaxNoOfSlicesInConcat = (long)ceil( (double)pMrProt->sliceSeries().size() / pMrProt->concatenations() );



    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    // * FatSat                                                                 *
    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK ) {


        switch ( pMrProt->preparationPulses().fatSatMode() )  {

            case SEQ::FAT_SAT_WEAK:
                lMaxTimeBetweenQuickFatSats_us = lMaxTimeBetweenQuickFatSatsWeak_us;
            break;

            case SEQ::FAT_SAT_STRONG:
                lMaxTimeBetweenQuickFatSats_us = lMaxTimeBetweenQuickFatSatsStrong_us;
            break;

            default:
                TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown fat sat mode = %d  ", ptModule, pMrProt->preparationPulses().fatSatMode() );
                return ( SEQU_ERROR );
        }


        CheckStatusB ( lStatus = fCalculateNoOfQuickSatsPerConcat (pMrProt,
                                                                  lMaxTimeBetweenQuickFatSats_us,
                                                                  lMaxNoOfSlicesInConcat,
                                                                 &lNoOfQuickFatSatsPerConcat,
                                                                 &lNoOfKernelsBetweenQuickFatSats,
                                                                 &bIgnoreQuickFatSat) );
    } else {

        lNoOfQuickFatSatsPerConcat      = 0;
        lNoOfKernelsBetweenQuickFatSats = 0;
        bIgnoreQuickFatSat              = false;

    }



    if ( ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK ) && !bIgnoreQuickFatSat )  {
        CSatFat.setRequestsPerMeasurement ( lRequestsCRMSatPulses / pMrProt->sliceSeries().size() * pMrProt->concatenations() * lNoOfQuickFatSatsPerConcat );
    } else {
        CSatFat.setRequestsPerMeasurement ( lRequestsCRMSatPulses );
    }
    CSatFat.setCSatMode                   ( SBBCSatCode_Fat       );
    CSatFat.setIdent                      ( "greFS"               );
    CSatFat.setGSWDGradientPerformance    ( pMrProt, pSeqLim      );



    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    // * WaterSat                                                               *
    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) {

        CheckStatusB ( lStatus = fCalculateNoOfQuickSatsPerConcat (pMrProt,
                                                                  lMaxTimeBetweenQuickWaterSats_us,
                                                                  lMaxNoOfSlicesInConcat,
                                                                 &lNoOfQuickWaterSatsPerConcat,
                                                                 &lNoOfKernelsBetweenQuickWaterSats,
                                                                 &bIgnoreQuickWaterSat) );
    } else {

        lNoOfQuickWaterSatsPerConcat      = 0;
        lNoOfKernelsBetweenQuickWaterSats = 0;
        bIgnoreQuickWaterSat              = false;

    }


    if ( ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) && !bIgnoreQuickWaterSat )  {
        CSatWater.setRequestsPerMeasurement ( lRequestsCRMSatPulses / pMrProt->sliceSeries().size() * pMrProt->concatenations() * lNoOfQuickWaterSatsPerConcat );
    } else {
        CSatWater.setRequestsPerMeasurement ( lRequestsCRMSatPulses );
    }

    CSatWater.setCSatMode                 ( SBBCSatCode_Water     );
    CSatWater.setIdent                    ( "greWS"               );
    CSatWater.setGSWDGradientPerformance  ( pMrProt, pSeqLim      );



    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    // * P / RSat                                                               *
    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *

    // * ---------------------------------------------------------------------- *
    // * Check how many regular sats are carried out using the quick mode       *
    // * ---------------------------------------------------------------------- *
    lNumberOfQuickPRSats   = lNumberOfRegularPRSats = 0;

    for (lI=0; lI<pMrProt->satList().size(); lI++)  {
        if ( pMrProt->satList()[lI].quickMode() )  {
            lNumberOfQuickPRSats++;
        } else {
            lNumberOfRegularPRSats++;
        }
    }



    // * ---------------------------------------------------------------------- *
    // * Calculate the number of quick RSats per concatenation                  *
    // * Quick RSats are syncronized with quick FatSats                         *
    // * ---------------------------------------------------------------------- *
    if ( lNumberOfQuickPRSats > 0 ) {

        CheckStatusB ( lStatus = fCalculateNoOfQuickSatsPerConcat (pMrProt,
                                                                   lMaxTimeBetweenQuickRSats_us,
                                                                   lMaxNoOfSlicesInConcat,
                                                                  &lNoOfQuickRSatsPerConcat,
                                                                  &lNoOfKernelsBetweenQuickRSats,
                                                                  &bIgnoreQuickRSat,
                                                                   maximum(1L, lNoOfKernelsBetweenQuickFatSats)) );
    } else {

        lNoOfQuickRSatsPerConcat      = 0;
        lNoOfKernelsBetweenQuickRSats = 0;
        bIgnoreQuickRSat              = false;

    }

    // * ---------------------------------------------------------------------- *
    // * Calculate the number of RSat requests                                  *
    // * ---------------------------------------------------------------------- *
    for (lI=0; lI<lNoOfRSats; lI++)  {
        if ( pMrProt->satList()[lI].quickMode() && !bIgnoreQuickRSat )  {
            SBBRSat[lI].setRequestsPerMeasurement( lRequestsCRMSatPulses / pMrProt->sliceSeries().size() * pMrProt->concatenations() * lNoOfQuickRSatsPerConcat );
        } else {
            SBBRSat[lI].setRequestsPerMeasurement( lRequestsCRMSatPulses );
        }
        sprintf ( ptIdentdummy, "greRS%1ld", lI+1 );
        SBBRSat[lI].setIdent                   ( ptIdentdummy     );
        SBBRSat[lI].setGSWDGradientPerformance ( pMrProt, pSeqLim );
    }

#ifdef SUPPORT_CT
    if( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct )
    {
        //  The spoiler SatSpoil is not used by fSEQRunKernel_ct
        SBBRSat[0].setGradSpoiling(pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION || pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK ? false : true);
        for (lI=1; lI<lNoOfRSats; lI++)
        {
            SBBRSat[lI].setLastRampDownOutsideSBB(true);
        }
        SBBRSat[0].setLastRampDownOutsideSBB(pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION || pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK );
    }
    else
    {
        //  Restore the default
        SBBRSat[0].setGradSpoiling(SBBRSat[0].getSATNumber() == 0 ? false : true);
        for (lI=0; lI<lNoOfRSats; lI++)
        {
            SBBRSat[lI].setLastRampDownOutsideSBB(false);
        }
    }
#endif








    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    // * MSat                                                                   *
    // * ---------------------------------------------------------------------- *
    // * ---------------------------------------------------------------------- *
    MSat.setRequestsPerMeasurement      ( lRequestsCRMSatPulses );
    MSat.setFrequencyOffset             ( MSATOFFSETFREQHz      );
    MSat.setIdent                       ( "greMS"               );
    MSat.setGSWDGradientPerformance     ( pMrProt, pSeqLim      );

    // * EffTrTime is the time between two MSats in ms *
    switch ( pMrProt->kSpace().multiSliceMode() )  {
        case SEQ::MSM_SEQUENTIAL:
            MSat.setEffTrTime               ( static_cast<double>(pMrProt->tr()[0]) / 1000.0 );
        break;

        case SEQ::MSM_INTERLEAVED:
            MSat.setEffTrTime               ( static_cast<double>( pMrProt->tr()[0]              ) /
                                              static_cast<double>( pMrProt->sliceSeries().size() ) *
                                              static_cast<double>( pMrProt->concatenations()     ) / 1000.0  );
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown multi slice mode = %d  ", ptModule, pMrProt->kSpace().multiSliceMode());
            return ( SEQU_ERROR );
    }



    // * ---------------------------------------------------------------------- *
    // * TSat                                                                   *
    // * ---------------------------------------------------------------------- *
    SBBTSat.setRequestsPerMeasurement      ( lRequestsTSatPulses   );
    SBBTSat.setIdent                       ( "greTS"               );
    SBBTSat.setGSWDGradientPerformance     ( pMrProt, pSeqLim      );


    // * ---------------------------------------------------------------------- *
    // * SatSpoil                                                               *
    // * ---------------------------------------------------------------------- *
    if (    SysProperties::isTrio()
         && (( lNumberOfQuickPRSats > 0 ) ||
             ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK )  ||
             ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK )) )  {
        dConstSatSpoilMoment = 31000.0;
    } else {
        dConstSatSpoilMoment = 0.0;
    }



    if ( fabs(theMSU.getNominalB0() - 3.0) <= 0.1 ) {

        if ( ( lNumberOfQuickPRSats > 0 ) ||
             ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK )  ||
             ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) )  {

            dMinSatSpoilMoment   = 135000.0;
            dMaxSatSpoilMoment   = 135000.0;

        } else {

            dMinSatSpoilMoment   = 65000.0;
            dMaxSatSpoilMoment   = 65000.0;

        }

    } else {

        if ( ( lNumberOfQuickPRSats > 0 ) ||
             ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK )  ||
             ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) )  {

            dMinSatSpoilMoment   = 31000.0;
            dMaxSatSpoilMoment   = 31000.0;

        } else {

            dMaxSatSpoilMoment   = 15000.0;
            dMinSatSpoilMoment   = 15000.0;

        }

    }



    SatSpoil.setPESpoilMomentConstant ( dMaxSatSpoilMoment   );
    SatSpoil.setROSpoilMomentConstant ( dConstSatSpoilMoment );
    SatSpoil.setSSSpoilMomentConstant ( dMaxSatSpoilMoment   );

    if ( pMrProt->tSat().on() )  {
        SatSpoil.setRequestsPerMeasurement ( lRequestsTSatPulses   );
    } else {
        SatSpoil.setRequestsPerMeasurement ( lRequestsCRMSatPulses );
    }
    SatSpoil.setGSWDGradientPerformance    ( pMrProt, pSeqLim      );
    SatSpoil.prep ( pMrProt, pSeqLim, pSeqExpo );           // * Prepare for the maximum gradient moment *
    lSatSpoilDuration = SatSpoil.getDurationPerRequest();

#ifdef SUPPORT_CT
    if( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct )
    {
        SatSpoil.resetPrepared();
        lSatSpoilDuration = 0;
    }
#endif





    // * ---------------------------------------------------------------------- *
    // * Define the number of line/grid tag pulses                              *
    // * Only one sat pulse for all slices and phases                           *
    // * ---------------------------------------------------------------------- *
    lRequestsLineGridPulses = RunLoop.getKernelRequestsPerMeasurement (pSeqLim) / lPhasesToMeasure;
    if ( pMrProt->phaseStabilize() && (RunLoop.getKernelRequestsPerMeasurement (pSeqLim) % lPhasesToMeasure != 0) ) {
        lRequestsLineGridPulses += pMrProt->sliceSeries().size();
    }


    // * ---------------------------------------------------------------------- *
    // * Grid tagging                                                           *
    // * ---------------------------------------------------------------------- *
    SBBGridTag.setRequestsPerMeasurement   ( lRequestsLineGridPulses );
    SBBGridTag.setIdent                    ( "greGT"                 );
    SBBGridTag.setGSWDGradientPerformance  ( pMrProt, pSeqLim        );


    // * ---------------------------------------------------------------------- *
    // * Line tagging                                                           *
    // * ---------------------------------------------------------------------- *
    SBBLineTag.setRequestsPerMeasurement   ( lRequestsLineGridPulses );
    SBBLineTag.setIdent                    ( "greLT"                 );
    SBBLineTag.setGSWDGradientPerformance  ( pMrProt, pSeqLim        );


    // * ---------------------------------------------------------------------- *
    // * Preparation of all SBBs                                                *
    // * ---------------------------------------------------------------------- *
    if(!SBB.prepSBBAll(pMrProt,pSeqLim,pSeqExpo,&dRFEnergyInSBBs)) return(SBB.getpSBBLastPrep ()->getNLSStatus());



    // * ---------------------------------------------------------------------- *
    // * Calculate scan times for sat-pulses                                    *
    // * ---------------------------------------------------------------------- *
    for (lI=0;lI<lNoOfRSats;lI++)  {
        if ( SBBRSat[lI].isQuickMode() && !bIgnoreQuickRSat )  {
            lQuickRSatTimeInsideTR        += SBBRSat[lI].getDurationPerRequest();
        } else {
            lScanTimeSatsOnlyFirstSegment += SBBRSat[lI].getDurationPerRequest();
        }
    }

    if ( lNumberOfQuickPRSats > 0 && !bIgnoreQuickRSat )  {
        lQuickRSatTimeInsideTR *= lNoOfQuickRSatsPerConcat;
    }


    if ( ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK  && !bIgnoreQuickFatSat ) ) {
        lQuickFatSatTimeInsideTR      += CSatFat.getDurationPerRequest() * lNoOfQuickFatSatsPerConcat;
    } else {
        lScanTimeSatsOnlyFirstSegment += CSatFat.getDurationPerRequest();
    }

    if ( ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK  && !bIgnoreQuickWaterSat ) ) {
        lQuickWaterSatTimeInsideTR    += CSatWater.getDurationPerRequest() * lNoOfQuickWaterSatsPerConcat;
    } else {
        lScanTimeSatsOnlyFirstSegment += CSatWater.getDurationPerRequest();
    }

    lScanTimeSatsOnlyFirstSegment += MSat.getDurationPerRequest();
    lScanTimeSatsEachSegment      += SBBTSat.getDurationPerRequest();



    if ( bIgnoreQuickRSat || bIgnoreQuickFatSat || bIgnoreQuickWaterSat || pMrProt->tSat().on() || (lNumberOfRegularPRSats > 0) ||
         (pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION)     ||
         (pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION) ||
         pMrProt->preparationPulses().MTC() )  {

        if ( pMrProt->tSat().on() )  {
            lScanTimeSatsEachSegment      += SatSpoil.getDurationPerRequest();
        } else {

            lScanTimeSatsOnlyFirstSegment += SatSpoil.getDurationPerRequest();
        }

    } else {

        if ( ( lNumberOfQuickPRSats > 0 )                                                        ||
             ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK )      ||
             ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) )  {

            for ( lI = 0; lI < lMaxNoOfSlicesInConcat; lI++ )  {

                if ( (lNumberOfQuickPRSats > 0) && !(lI % lNoOfKernelsBetweenQuickRSats) ) {

                    lNoOfQuickSpoilGradsPerConcat++;

                } else {

                    if ( (pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK) && !(lI % lNoOfKernelsBetweenQuickFatSats) ) {

                        lNoOfQuickSpoilGradsPerConcat++;

                    } else {

                        if ( (pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK) && !(lI % lNoOfKernelsBetweenQuickWaterSats) ) {

                            lNoOfQuickSpoilGradsPerConcat++;

                        }

                    }

                }

            }

            lQuickSpoilGradTimeInsideTR += SatSpoil.getDurationPerRequest() * lNoOfQuickSpoilGradsPerConcat;
                dSatSpoilMomentStep =  (dMaxSatSpoilMoment - dMinSatSpoilMoment) / lSpoilMomentIncr;

        }

    }

    RunLoop.setQuickSatTimeInsideTR ( lQuickRSatTimeInsideTR + lQuickFatSatTimeInsideTR + lQuickWaterSatTimeInsideTR +
                                      lQuickSpoilGradTimeInsideTR );

    // * ---------------------------------------------------------------------- *
    // * Calculation of the TR and TI fill times                                *
    // * ---------------------------------------------------------------------- *

        RunLoop.setbHandleTRTIConflict ( true );
        lTRMin = GREKernel.getlTRMin();
    #ifdef GSWD
        if ( GREKernel.getalTEFill()[0] < 0 )  {
            lTRMin -= 2 * GREKernel.getalTEFill()[0];
        }

        if ( ( GREKernel.getalTEFill()[1] < 0 ) && ( pMrProt->contrasts() == 2 ) )  {
            lTRMin -= 2 * GREKernel.getalTEFill()[1];
        }
    #endif

    if ( pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF )  {
        bSuccess = RunLoop.TrTiFillTimes ( pMrProt, pSeqLim, pSeqExpo,
                                           lScanTimeSatsOnlyFirstSegment + (lScanTimeSatsEachSegment + lTRMin) * lSegmentsInTR, 1,
                                           0,
                                           0,
                                           GREKernel.getlTRMin(),
                                           lScanTimeSatsOnlyFirstSegment + lScanTimeSatsEachSegment * lSegmentsInTR,
                                          &lNegativeFillTime );

        if ( !bSuccess && !pSeqLim->isContextPrepForMrProtUpdate() ) {  return ( RunLoop.getNLSStatus() );  }

    } else {
        bSuccess = RunLoop.TrTiFillTimes ( pMrProt, pSeqLim, pSeqExpo,
                                           lScanTimeSatsOnlyFirstSegment + (lScanTimeSatsEachSegment + lTRMin) * lSegmentsInTR,
                                           1,
                                           lScanTimeSatsOnlyFirstSegment + (lScanTimeSatsEachSegment + lTRMin) * (lSegmentsInTR / 2),
                                           lScanTimeSatsOnlyFirstSegment + (lScanTimeSatsEachSegment + lTRMin) * (lSegmentsInTR / 2),
                                           (lScanTimeSatsEachSegment + lTRMin) * (lSegmentsInTR / 2 + 1),
                                           lScanTimeSatsOnlyFirstSegment + (lScanTimeSatsEachSegment + lTRMin) * (lSegmentsInTR / 2),
                                          &lNegativeFillTime );

        if ( !bSuccess && !pSeqLim->isContextPrepForMrProtUpdate() )  {  return ( RunLoop.getNLSStatus() );  }

    }

    //  Since we need the TR-Fill CSat flip angle can't be calculated earlier
#ifdef SUPPORT_CT
    if( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct )
    {
        const int32_t i32NSeg              = RunLoop.getLinesToMeasure()
            ,         i32NLinPerSeg        = RunLoop.getPartitionsToMeasure()
            ,         i32NPrepScans1stSeg  = i32NSeg*i32NLinPerSeg-lLinesToMeasure
            ;
        //  First find the segment which contains the k-space center line
        int32_t i32CSeg, i32Index, i32CLinInSeg;
        for(i32CSeg = 0; i32CSeg != i32NSeg;++i32CSeg)
        {
            for(i32CLinInSeg = 0; i32CLinInSeg != i32NLinPerSeg;++i32CLinInSeg)
            {
                //i32Index = i32CSeg * i32NLinPerSeg + i32NLinPerSeg - i32NPrepScans1stSeg;
                i32Index = i32CSeg*i32NLinPerSeg+i32CLinInSeg-i32NPrepScans1stSeg;
                if( i32Index < 0 )
                {
                    continue;
                }
                if( Reorder.getLinNo(i32Index) == Reorder.getKSCenterLin() )
                {
                    s_i32ZSeg = i32CSeg;
                    //  Terminate outer loop
                    i32CSeg = i32NSeg-1;
                    break;
                }
            }
        }
        s_i32CSatShift = 0;
        if( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK )
        {
            if( CSAT_MODE(pMrProt) == CSAT_MODE_AUTO )
            {
                const int32_t i32NSlc  = pMrProt->sliceSeries().size()
                    ,         i32NCSat = minimum(pMrProt->MDS().mdsTableSpeedNumerator(),i32NSlc)
                    ;
                const double dTimePerShift_us = GREKernel.getlTRMin()+RunLoop.getlTRFillInConcat(0)
                    ,        dTimeConst_us    = CSatFat.getDurationPerRequest()-CSatFat.getdTimeToRFCenter_us()+GREKernel.getTimeToRFCenter_us()
                    ,        dT1_us           = T1_FAT_us()
                    ;
                double dCSatFlipAngle_deg = CSATFLIPANGLEdeg;
                if( !CT_QSAT_SHIFT_AND_FLIP_ANGLE(dCSatFlipAngle_deg,s_i32CSatShift,s_i32ZSeg,i32NCSat,double(pMrProt->tr()[0]),dT1_us,dTimeConst_us,dTimePerShift_us) )
                {
                    TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
                    return SEQU_ERROR;
                }
                dRFEnergyInSBBs -= CSatFat.getEnergyPerRequest() *  CSatFat.getRequestsPerMeasurement() * pMrProt->measurements();
                CSatFat.setFlipAngle(dCSatFlipAngle_deg);
                if( !CSatFat.prep(pMrProt,pSeqLim,pSeqExpo) )
                {
                    TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
                    return CSatFat.getNLSStatus();
                }
                dRFEnergyInSBBs += CSatFat.getEnergyPerRequest() *  CSatFat.getRequestsPerMeasurement() * pMrProt->measurements();
            }
            else
            {
                s_i32CSatShift = minimum(CSAT_SHIFT(pMrProt),s_i32ZSeg);
            }
        }
    }  //  END MDS ON
#endif
    // * ---------------------------------------------------------------------- *
    // * Calculate the total measurement time of one repetition                 *
    // * ---------------------------------------------------------------------- *
    // * Scan time without sats *
    dMeasureTimeUsec = RunLoop.getMeasurementTimeUsec ( pMrProt, pSeqLim );




    // * ---------------------------------------------------------------------- *
    // * Calculate the minimum additional measurement time of dummy scans in    *
    // * case of respiratory triggering.                                        *
    // * ---------------------------------------------------------------------- *
    if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) ) {
        dMeasureTimeUsec += (lLinesToMeasure / lSegments) * (double)(RunLoop.getlSteadyStateTriggerDebugDummyScans()) * (double)(pMrProt->tr()[0]) * (double)(pMrProt->averages()) * (double)(pMrProt->concatenations());

        if ( pMrProt->phaseStabilize() ) {
            dMeasureTimeUsec += pMrProt->tr()[0] * RunLoop.getlSteadyStateTriggerDebugDummyScans();
        }

    }

    // * ---------------------------------------------------------------------- *
    // * Calculate the measurement time of line and grid tag pulses             *
    // * ---------------------------------------------------------------------- *
    dMeasureTimeUsec += lRequestsLineGridPulses * SBBLineTag.getDurationPerRequest() +
                        lRequestsLineGridPulses * SBBGridTag.getDurationPerRequest();




    RunLoop.setdMeasureTimeUsec ( dMeasureTimeUsec );


    // * ---------------------------------------------------------------------- *
    // * Calculate the total measurement time, including measurement repeats    *
    // * ---------------------------------------------------------------------- *
    dTotalMeasureTimeUsec = RunLoop.getTotalMeasTimeUsec ( pMrProt, pSeqLim );




    // * ---------------------------------------------------------------------- *
    // * Calculate the energy of the GRE kernel                                 *
    // * ---------------------------------------------------------------------- *
    dRFEnergyInKernel = GREKernel.getEnergyPerRequest() * GREKernel.getRequestsPerMeasurement() *
                        lNoOfMeasurements;


    // * ---------------------------------------------------------------------- *
    // * Calculate the energy of the recovery pulses                            *
    // * ---------------------------------------------------------------------- *
    // * Sequence execution                                                     *
    dRFEnergyInPrepPulse =  ( RunLoop.getRfEnergyRecovery() + RunLoop.getRfEnergyDB() ) *
                            lNoOfMeasurements * pMrProt->sliceSeries().size() * pMrProt->averages() *
                            lPhasesToMeasure  * Reorder.getNoOfReorderIndices() / lSegmentsInTR;

    // * Preparation scans                                                      *
    // * getPreparingScansTotal() returns the number of prep. scans per         *
    // * measurement.                                                           *
    dRFEnergyInPrepPulse += ( RunLoop.getRfEnergyRecovery() + RunLoop.getRfEnergyDB() )
                            * RunLoop.getPreparingScansTotal() * lNoOfMeasurements;

    // * Phase stabilization reference scan                                     *
    // * IR-pulses are switched off during phase stabilize reference scans      *
    if ( pMrProt->phaseStabilize() )  {
        dRFEnergyInPrepPulse += RunLoop.getRfEnergyDB() * pMrProt->sliceSeries().size();
    }



    // * ---------------------------------------------------------------------- *
    // * Calculate the minimum additional energy of the kernel and sat pulses   *
    // * in dummy scans in case of respiratory triggering.                      *
    // * ---------------------------------------------------------------------- *
    if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) ) {
        dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->sliceSeries().size() * GREKernel.getEnergyPerRequest();
        dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->sliceSeries().size() * SBBTSat.getEnergyPerRequest();

        for ( lI = 0; lI < lNoOfRSats; lI++) {
            if ( lNumberOfQuickPRSats > 0 ) {
                dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->concatenations() * lNoOfQuickRSatsPerConcat * SBBRSat[lI].getEnergyPerRequest();
            } else {
                dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->sliceSeries().size() * SBBRSat[lI].getEnergyPerRequest();
            }
        }

        if ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK ) {
            dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->concatenations() * lNoOfQuickFatSatsPerConcat * CSatFat.getEnergyPerRequest();
        } else if ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION ) {
            dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->sliceSeries().size() * CSatFat.getEnergyPerRequest();
        }

        if ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) {
            dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->concatenations() * lNoOfQuickFatSatsPerConcat * CSatWater.getEnergyPerRequest();
        } else if ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION ) {
            dRFEnergyInDummyScans += (long)(lLinesToMeasure / lSegments) * pMrProt->averages() * pMrProt->sliceSeries().size() * CSatWater.getEnergyPerRequest();
        }

        if ( pMrProt->phaseStabilize() ) {
            dRFEnergyInDummyScans += pMrProt->sliceSeries().size() * GREKernel.getEnergyPerRequest();
            dRFEnergyInDummyScans += pMrProt->sliceSeries().size() * SBBTSat.getEnergyPerRequest();

            for ( lI = 0; lI < lNoOfRSats; lI++) {
                if ( lNumberOfQuickPRSats > 0 ) {
                    dRFEnergyInDummyScans += pMrProt->concatenations() * lNoOfQuickRSatsPerConcat * SBBRSat[lI].getEnergyPerRequest();
                } else {
                    dRFEnergyInDummyScans += pMrProt->sliceSeries().size() * SBBRSat[lI].getEnergyPerRequest();
                }
            }

            if ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK ) {
                dRFEnergyInDummyScans += pMrProt->concatenations() * lNoOfQuickFatSatsPerConcat * CSatFat.getEnergyPerRequest();
            } else if ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION ) {
                dRFEnergyInDummyScans += pMrProt->sliceSeries().size() * CSatFat.getEnergyPerRequest();
            }

            if ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) {
                dRFEnergyInDummyScans += pMrProt->concatenations() * lNoOfQuickFatSatsPerConcat * CSatWater.getEnergyPerRequest();
            } else if ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION ) {
                dRFEnergyInDummyScans += pMrProt->sliceSeries().size() * CSatWater.getEnergyPerRequest();
            }

        }

        dRFEnergyInDummyScans *= RunLoop.getlSteadyStateTriggerDebugDummyScans();
    }




    dRFEnergyInSBBs      += dRFEnergyInKernel;
    dRFEnergyInSBBs      += dRFEnergyInPrepPulse;
    dRFEnergyInSBBs      += dRFEnergyInDummyScans;

    #ifdef SUPPORT_PACE
    // * ---------------------------------------------------------------------- *
    // * Add energy emitted by the navigator                                    *
    // * ---------------------------------------------------------------------- *
        dRFEnergyInSBBs  += RunLoop.getScoutEnergy_Ws();
    #endif


    // * ---------------------------------------------------------------------- *
    // * Add energy from iPAT reference scan  (SBBPATRefScan)                   *
    // * ---------------------------------------------------------------------- *
    //
    if (pMrProt->PAT().RefScanMode()==SEQ::PAT_REF_SCAN_EXTRA)  {
        dRFEnergyInSBBs  += RunLoop.getdTotalRfEnergyPATRefScan();
    }



    // * ---------------------------------------------------------------------- *
    // * Set the gain of the receiver                                           *
    // * ---------------------------------------------------------------------- *
    lStatus = fSSLSetRxGain( fSURxGainCode(pMrProt, pSeqLim, THICKNESSmm_RxGain, GRE_Sequence) , pMrProt, pSeqLim);
    CheckStatusB(lStatus);



    // * ---------------------------------------------------------------------- *
    // * Sequence short string                                                  *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->SWI() )  {

        char String[256], ShortString[32], ContrastStr[2];

        String[0] = '\0';
        ShortString[0] = '\0';

        strcat (ShortString, "swi");

        switch ( pMrProt->kSpace().dimension() )  {
            case SEQ::DIM_2:   strcat (String, "swi2d");   break ;
            case SEQ::DIM_3:   strcat (String, "swi3d");   break ;
            default: TRACE_PUT1(TC_ALWAYS, TF_SEQ, "%s ERROR: Invalid dimension found in MrProt.", ptModule); return (false);
        }


        if ( pMrProt->contrasts() < 10 )  {
            sprintf (ContrastStr, "%1ld", pMrProt->contrasts());
        } else {
            sprintf (ContrastStr, "%2ld", pMrProt->contrasts());
        }

        strcat (String, ContrastStr);

        switch ( pMrProt->flowComp()[0] )  {
            case SEQ::FLOW_COMPENSATION_YES:   strcat (String, "r");   strcat(ShortString, "_r");   break;
            default: TRACE_PUT1(TC_ALWAYS, TF_SEQ, "%s ERROR: Invalid flow comp mode found in MrProt.", ptModule); return (false);
        }


        pSeqExpo->setSeqShortString (ShortString);
        pSeqExpo->setSequenceString (String);

    } else {

        if ( pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF )  {
            if ( pMrProt->fastImaging().RFSpoiling() )  {
                fSUSetSequenceString ("fl", pMrProt, pSeqExpo);
            } else {
                fSUSetSequenceString ("fi", pMrProt, pSeqExpo);
            }
        } else {
            fSUSetSequenceString ("tfl", pMrProt, pSeqExpo);
        }

    }






    // * ---------------------------------------------------------------------- *
    // * Prepare exports                                                        *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->calcRealPartImages() )  {
        pSeqExpo->setPhaseCorScans ( pMrProt->contrasts() * pMrProt->fastImaging().segments() );
    } else {
        pSeqExpo->setPhaseCorScans ( 0 );
    }
    #ifndef VXWORKS
        if (    ( pSeqLim->isContextPrepForScanTimeCalculation() || pSeqLim->isContextNormal() )
             && (RunLoop.getTrigHaltDuration()) )  {
            long lPhysioHalts = long (RunLoop.getdScanTimeTrigHalt (FirstMeas) + RunLoop.getdScanTimeTrigHalt (SecondMeas) * pMrProt->repetitions()) / RunLoop.getTrigHaltDuration();
            SeqUT.setExpectedPhysioHalts (lPhysioHalts);
            long lThreshold_us = pMrProt->tr()[0] * maximum (1L, pMrProt->physiology().phases()) * maximum (1L, pSeqExpo->getAcqWindowFactor())
                                  + pMrProt->physiology().triggerDelay(pMrProt->physiology().signal(1));
            pSeqExpo->setEstimatedMeasureTimeUsec (fSUGetEstimatedMeasTimeUsec(pMrProt, pSeqExpo, lPhysioHalts, lThreshold_us / 1000));
        } else {
            pSeqExpo->setEstimatedMeasureTimeUsec (0);
        }
    #endif


    pSeqExpo->setPreScans                      ( RunLoop.getlPreparingScans() );
    pSeqExpo->setMeasureTimeUsec               ( dMeasureTimeUsec );                  /*! EGA-08 !*/
    pSeqExpo->setTotalMeasureTimeUsec          ( dTotalMeasureTimeUsec );             /*! EGA-08 !*/
    pSeqExpo->setRFEnergyInSequence_Ws         ( dRFEnergyInSBBs );                   /*! EGA-08 !*/
    pSeqExpo->setMeasuredPELines               ( lLinesToMeasureMax );
    pSeqExpo->setMeasured3dPartitions          ( lPartitionsToMeasureMax );
    pSeqExpo->setRelevantReadoutsForMeasTime   ( RunLoop.getNumberOfRelevantADCs() );
    pSeqExpo->setApplicationCard               ( SEQ::APPLICATION_CARD_INLINE );
    pSeqExpo->setNoOfDisplayedAsymmetricEchoes (pMrProt->contrasts());
    for ( lI = 0; lI < pMrProt->contrasts(); lI++ ) {
        pSeqExpo->setAsymmetricEcho            ( lI, GREKernel.getadReadOutAsymUI(lI) );
    }

#ifdef SUPPORT_CT
#ifndef VXWORKS // Do Energy estimation & export only on the host!!!
    if( pSeqLim->isContextNormal() && pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO )
    {
        double dEnergyPerTR_Ws = dRFEnergyInSBBs
            / (lNoOfMeasurements*RunLoop.getLinesToMeasure())
            ;
        lStatus = RunLoop.ct_exportRFBlocks(pMrProt,pSeqLim,pSeqExpo,dEnergyPerTR_Ws);
        if( NLS_SEVERITY(lStatus) != NLS_SUCCESS )
        {
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
            return lStatus;
        }
    } // ContextNormal
#endif  //  VXWORKS not defined
#endif  //  SUPPORT_CT


    // * -------------------------------------------------------------------------- *
    // * Increase the FFT scale factor for 3D imaging                               *
    // * -------------------------------------------------------------------------- *
    if ( pMrProt->kSpace().dimension() == SEQ::DIM_3 )  {
        pSeqExpo->setAdditionalScaleFactor ( 2.0 );
    } else {
        pSeqExpo->setAdditionalScaleFactor ( 1.0 );
    }




    // * ---------------------------------------------------------------------- *
    // * Set phase correction parameter                                         *
    // * ---------------------------------------------------------------------- *
    pSeqExpo->setPCAlgorithm ( SEQ::PC_ALGORITHM_NONE );





    // * ---------------------------------------------------------------------- *
    // * Online fft                                                             *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->fastImaging().segments() > 1 )  {
        RunLoop.setbSegmentedLinesOrPartitions ( true  );
    } else {
        RunLoop.setbSegmentedLinesOrPartitions ( false );
    }


    // * ---------------------------------------------------------- *
    // * Restrict the number of relevant measurements used to       *
    // * calculate the amount of data, in order to allow multiple   *
    // * 3D-measurements even if they do exceed the memory limit    *
    // * (in particular if many receive coil-channels are used).    *
    // * Assumption: reconstruciton is faster than data limit is    *
    // *             reached during acquisition.                    *
    // * Since this assumption is difficult to assure for any       *
    // * protocol, we do this only for breast application!          *
    // * ---------------------------------------------------------- *
    if( (pMrProt->kSpace().dimension() == SEQ::DIM_3) && (pMrProt->contrasts() < 2) )  {    // ... furtheron restrict to 3d and single contrast

        const MrCoilSelect &coilSelect = pMrProt->coilInfo().Meas();
        SelectedCoilElements mySelectedElements(coilSelect, pSeqLim->getCoilContext());

        if (mySelectedElements.isOneApplProp("BREAST"))  {      // Breast application detected by checking the application properties of the coil

            pSeqExpo->setICEProgramParam (ICE_PROGRAM_PARA_RELEVANT_NO_REPS, 2);  // restrict data calculation to 2 repetitions
            if (! pSeqLim->isContextPrepForBinarySearch() )  {
                TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: Warning: restricted number of repetitions for memory calculation used", ptModule);
            }
        } else {
            pSeqExpo->setICEProgramParam (ICE_PROGRAM_PARA_RELEVANT_NO_REPS, 0);  // reset restriction
        }
    } else {
        pSeqExpo->setICEProgramParam (ICE_PROGRAM_PARA_RELEVANT_NO_REPS, 0);  // reset restriction
    }


    switch ( pMrProt->kSpace().dimension() )  {
        case SEQ::DIM_2:                                                                                           // * Contrasts and segments are handeled inside the kernel *
            if ( (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh == SEQ::METHOD_TRIGGERING) ) {
                // * The number of segments must not be passed as a additional      *
                // * multiplier since it is already considered by lLinesToMeasure.  *
                if (!RunLoop.setIceProgram (pMrProt, pSeqLim, pSeqExpo, RunLoop.getDataRateMBytePerSec(pMrProt, pSeqLim, pMrProt->contrasts()), SEQ::RS_LIN_IN_PAR ) ) return RunLoop.getNLSStatus();
            } else {
                if (!RunLoop.setIceProgram (pMrProt, pSeqLim, pSeqExpo, RunLoop.getDataRateMBytePerSec(pMrProt, pSeqLim, pMrProt->contrasts() * pMrProt->fastImaging().segments()), SEQ::RS_LIN_IN_PAR ) ) return RunLoop.getNLSStatus();
            }


            if (    (pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_OFF        )
                 || (pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_BREATH_HOLD) )  {
                pSeqExpo->setOnlineFFT (SEQ::ONLINE_FFT_PHASE);
            }

        break;

        case SEQ::DIM_3:
            // * ------------------------------------------------------------------ *
            // * Select ICE program                                                 *
            // * ------------------------------------------------------------------ *
            switch ( Reorder.getBasicScanningScheme() )  {
                case PAR_IN_LIN:
                    if (!RunLoop.setIceProgram (pMrProt, pSeqLim, pSeqExpo, RunLoop.getDataRateMBytePerSec(pMrProt, pSeqLim, pMrProt->contrasts()) * pMrProt->fastImaging().segments(), SEQ::RS_PAR_IN_LIN ) ) return RunLoop.getNLSStatus();
                break;

                case LIN_IN_PAR:
                    if (!RunLoop.setIceProgram (pMrProt, pSeqLim, pSeqExpo, RunLoop.getDataRateMBytePerSec(pMrProt, pSeqLim, pMrProt->contrasts()) * pMrProt->fastImaging().segments(), SEQ::RS_LIN_IN_PAR ) ) return RunLoop.getNLSStatus();
                break;

                default:
                    TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown loop order = %ld  ", ptModule, Reorder.getBasicScanningScheme());
                    return ( SEQU_ERROR );
            }

        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
            return ( SEQU_ERROR );
    }


    // Interactive

    if (! RTControl.prep(pMrProt, pSeqLim, pSeqExpo, asSLC)) {
        if (! pSeqLim->isContextPrepForBinarySearch() )
            TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: RTControl.prep failed.", ptModule);
        return(SEQU_ERROR) ;
    }

    // * ---------------------------------------------------------------------- *
    // * For Interactive Realtime, the online ICE program is needed             *
    // *  and the number of receiver channels must coincide with the setting    *
    // *  for switching on the OffIce-program in SeqLoop                        *
    // * ---------------------------------------------------------------------- *

    if ( eSequence == RT_Gre ) {
        pSeqExpo->setMaxReceiverChannels(RunLoop.getlChannelsLimitForOffICE());

        // * ---------------------------------------------------------------------- *
        // * For Interactive Realtime, the switch "endless measurement will set     *
        // *  the number of measurements to the highest possible value. To ensure   *
        // *  this behaviour also in CrDefProt and "repair", it is checked here.    *
        // * ---------------------------------------------------------------------- *

        if (pMrProt->endlessMeasurement() && (pSeqLim->getRepetitions().getMax() != pMrProt->repetitions())) {
            if (pSeqLim->isContextPrepForMrProtUpdate()) {
                pMrProt->repetitions(pSeqLim->getRepetitions().getMax());
            } else {
                if (! pSeqLim->isContextPrepForBinarySearch() )  {
                    TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: ERROR: Endless measurement needs maximum number of repetitions.", ptModule);
                }
                return SEQU_ERROR;
            }
        }
    }




    // * -------------------------------------------------------------------------- *
    // * Check of number of ICE blocks                                              *
    // * -------------------------------------------------------------------------- *



#ifdef SUPPORT_iPAT
        //---------------------------------------------------------------------------
        // call fPATPrepPost, which checks some PAT related restrictions
        //---------------------------------------------------------------------------
        lStatus = fPATPrepPost(pMrProt,pSeqLim,pSeqExpo,&Reorder);

        if (! pSeqLim->isContextPrepForBinarySearch() )  {
            CheckStatusPB (lStatus, "fPATPrepPost");
        } else  {
            CheckStatusB  (lStatus);
        }
#endif

    pSeqExpo->setICEProgramFilename  ("%CustomerIceProgs%\\IceVSend");


    /* --------------------------- */
    /* End of sequence preparation */
    /* --------------------------- */

    return ( lStatus );

FINISHED:

    pSeqExpo->setNegativeFillTime ( lNegativeFillTime );

    if ( NLS_SEVERITY(lStatus) == NLS_SUCCESS ) lStatus = SEQU_ERROR;
    return(lStatus);
}

/*[ Function ****************************************************************\
*
* Name        : fSEQCheck
*
* Description : Checks the real-time sequence for gradient overflows.
*
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQCheck
(
  MrProt         *pMrProt,     // * IMP: user choice parameters  *
  SeqLim         *pSeqLim,     // * IMP: limits from fSEQInit()  *
  SeqExpo        *pSeqExpo,    // * IMP: exports from fSEQPrep() *
  SEQCheckMode   * //pSeqCheckMode      // * IMP: Check mode            *
)

{


    static const char *ptModule = {"fSEQCheck"};
    NLS_STATUS   lStatus = SEQU__NORMAL;


    // * ---------------------------------------------------------------------- *
    // * Disable check for gre_highres                                          *
    // * ---------------------------------------------------------------------- *
    if ( eSequence == Gre_HR )  {
        return ( SEQU__NORMAL );
    }


  // * ---------------------------------------------------------------------- *
  // * Set the looping parameters                                             *
  // * ---------------------------------------------------------------------- *
  long lJ=0;
  long lNoOfRepetitions=0;



  if ( pMrProt->kSpace().ellipticalScanning() )  {

       lNoOfRepetitions = 2;

       for (lJ=0; lJ<lNoOfRepetitions; lJ++)  {

            RunLoop.setPairsOfLinesPartitionToCheck ( 0, lKSpaceCenterPartition );
            RunLoop.setPairsOfLinesPartitionToCheck ( lLinesToMeasureMax - 1, lKSpaceCenterPartition );

            RunLoop.setPairsOfLinesPartitionToCheck ( lKSpaceCenterLine, 0 );
            RunLoop.setPairsOfLinesPartitionToCheck ( lKSpaceCenterLine, lPartitionsToMeasureMax - 1 );

       }

  } else {

    switch ( pMrProt->kSpace().dimension() )  {
      case SEQ::DIM_2:

           lNoOfRepetitions = 4;

           for (lJ=0; lJ<lNoOfRepetitions; lJ++)  {
             RunLoop.setPairsOfLinesPartitionToCheck ( 0 );
             RunLoop.setPairsOfLinesPartitionToCheck ( lLinesToMeasureMax - 1 );
           }
      break;

      case SEQ::DIM_3:

           lNoOfRepetitions = 2;

           for (lJ=0; lJ<lNoOfRepetitions; lJ++)  {
             RunLoop.setPairsOfLinesPartitionToCheck ( 0, 0 );
             RunLoop.setPairsOfLinesPartitionToCheck ( 0, lPartitionsToMeasureMax - 1 );

             RunLoop.setPairsOfLinesPartitionToCheck ( lLinesToMeasureMax - 1, 0 );
             RunLoop.setPairsOfLinesPartitionToCheck ( lLinesToMeasureMax - 1, lPartitionsToMeasureMax - 1 );
           }
      break;

      default:
        TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
        return ( SEQU_ERROR );
    }

  }

  // * ---------------------------------------------------------------------- *
  // * Execute the check loops                                                *
  // * ---------------------------------------------------------------------- *
#ifdef SUPPORT_CT
  if (! RunLoop.check ( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct ? fSEQRunKernel_ct : fSEQRunKernel,
#else
  if (! RunLoop.check ( fSEQRunKernel,
#endif
                        pMrProt,
                        pSeqLim,
                        pSeqExpo,
                        asSLC,
                      &(GREKernel.ADC(0)) )
     ) return ( RunLoop.getNLSStatus() );



   #ifdef GRE_RT
       // This code checks for interactive scanning a higly stimulation-sensitive slice orientation
       //  Since a rotation matrix must be given, this corresponds to a different rotation matrix, if
       //  patient is lying on the side
       // Additionally, the original slice orientation is checked.
       // The original slice object has to be stored in a temporary variable, since fSeqRunKernel accesses the
       //  global object asSLC

    sSLICE_POS     asOrigSLC[1];
    sRTSlice_Pos   asCheckRTSLC[1];

    GradientPeRoSlVec  PRSOffVec;        // offset vector in Phase, Read, Slice, set to zero
    PRSOffVec.setPe(0);
    PRSOffVec.setRo(0);
    PRSOffVec.setSl(0);

    asOrigSLC[0] =  asSLC[0];                   // store original slice object
    double dRotMat[3][3] = {{0,0,0},{0,0,0},{0,0,0}} ;

    if (RTControl.getiOriginalPatPosition() < 3)  { // if the patient is lying supine/prone
        dRotMat[0][0]= 0.981365;dRotMat[0][1]=-0.118500;dRotMat[0][2]= 0.151261;
        dRotMat[1][0]= 0.192151;dRotMat[1][1]= 0.605209;dRotMat[1][2]=-0.772528;
        dRotMat[2][0]= 0.000000;dRotMat[2][1]= 0.787197;dRotMat[2][2]= 0.616701;
        } else  {                  // patient lies on the side
        dRotMat[0][0]= 0.222910;dRotMat[0][1]=-0.811064;dRotMat[0][2]=-0.540820;
        dRotMat[1][0]=-0.924546;dRotMat[1][1]=-0.000000;dRotMat[1][2]=-0.381070;
        dRotMat[2][0]= 0.309072;dRotMat[2][1]= 0.584958;dRotMat[2][2]=-0.749866;
    }

    GREKernel.updateRotationMatrix();
    asCheckRTSLC[0].prep (PRSOffVec, dRotMat, pMrProt);   // calculate RTSlc-object
    asSLC[0] = asCheckRTSLC[0];

    if (! RunLoop.check ( fSEQRunKernel,
                        pMrProt,
                        pSeqLim,
                        pSeqExpo,
                        asSLC,
                       &(GREKernel.ADC(0)) )
         ) lStatus = RunLoop.getNLSStatus();

    asSLC[0] = asOrigSLC[0];  // restore original slice object
    GREKernel.updateRotationMatrix();

    if ( NLS_SEVERITY(lStatus) != NLS_SUCCESS )  {
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s : SeqLoop.check of extreme slice returned an error : NLS_CODE : 0x%lx"  , ptModule, lStatus);
        return ( lStatus );
    }

  #endif


  /* --------------- */
  /* End of sequence */
  /* --------------- */


  return(lStatus);

}






/*[ Function ****************************************************************\
*
* Name        : fSEQRun
*
* Description : Executes the real-time sequence.
*
* Return      : An NLS status code.
*
\****************************************************************************/
NLS_STATUS fSEQRun
(
  MrProt       *pMrProt,      /* IMP: user choice parameters  */
  SeqLim       *pSeqLim,      /* IMP: limits from fSEQInit()  */
  SeqExpo      *pSeqExpo      /* IMP: exports from fSEQPrep() */
)

{
  static const char *ptModule = {"fSEQRun"};
  NLS_STATUS lStatus          = SEQU__NORMAL;


  mPrintTrace1 (DEBUG_RUN, DEBUG_CALL, "() <%s> started", pSeqLim->getLinkedSeqFilename() ) ;

  // * ---------------------------------------------------------------------- *
  // * Initialization of the unit test function                               *
  // * ---------------------------------------------------------------------- *
  mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ORIGIN_fSEQRunStart,0,0,0,0,0); /*! EGA-All !*/






  // * ---------------------------------------------------------------------- *
  // * Execute the measurement loops                                          *
  // * ---------------------------------------------------------------------- *
#ifdef SUPPORT_CT
    if (! RunLoop.run_new ( pMrProt->MDS().mdsModeMask() == SEQ::MDS_ONCO || s_bUseSEQRunKernel_ct ? fSEQRunKernel_ct : fSEQRunKernel,
#elif SUPPORT_PACE
    if (! RunLoop.run_new ( (NLS_STATUS (*)(TYPESFor_fSEQRunKernel)) fSEQRunKernel, // * IMP: pointer to kernel function *
#else
    if (! RunLoop.run ( (NLS_STATUS (*)(TYPESFor_fSEQRunKernel)) fSEQRunKernel,     // * IMP: pointer to kernel function *
#endif
                      pMrProt,                                                // * IMP: user choice parameters                 *
                      pSeqLim,                                                // * IMP: limits from fSEQInit()                 *
                      pSeqExpo,                                               // * IMP: exports from fSEQPrep()                *
                      asSLC,                                                  // * IMP: Rot.Matrices/Shifts                    *
                     &(GREKernel.ADC(0)) )                                    // * IMP: the ADC (for writing meas data header) *
        ) return ( RunLoop.getNLSStatus () );


  /* --------------- */
  /* End of sequence */
  /* --------------- */


  mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ORIGIN_fSEQRunFinish,0,0,0,0,0); /*! EGA-All !*/
  mPrintTrace1 (DEBUG_RUN, DEBUG_CALL | DEBUG_RETURN, "() <%s> finished", pSeqLim->getLinkedSeqFilename() ) ;
  return(lStatus);

}



/*[ Function ****************************************************************\
*
* Name        : fSEQRunKernel
*
* Description : Executes the basic timing of the real-time sequence.
*               This function is called by the function (libSBB)fSEQRunLoop.
*
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQRunKernel
(
  MrProt          *pMrProt,                // * IMP: user choice parameter    *
  SeqLim          *pSeqLim,                // * IMP: limits from fSEQInit()   *
  SeqExpo         *pSeqExpo,               // * IMP: exports of fSEQPrep()    *
  long            lKernelMode,             // * IMP: mode of kernel run       *
  long            lSlice,                  // * IMP: chronologic slice no.    *
  long            lPartition,              // * IMP: partition/segment number *
  long            lLine                    // * IMP: line/segment number      *
)

{
    static const char *ptModule         = {"fSEQRunKernel"};
    NLS_STATUS         lStatus          = SEQU__NORMAL;
    long               lI               = 0;   // * Loop counter *
    long               lSegments        = pMrProt->fastImaging().segments();
    long               lSegmentsInTR    = 0;
    long               lSegmentCounter  = 0;
    long               lIndex           = 0;

    bool               bPerformQuickRSat     = false;
    bool               bPerformQuickFatSat   = false;
    bool               bPerformQuickWaterSat = false;

    bool               bIsLastScanInMeas   = GREKernel.ADC(0).Mdh.isLastScanInMeas();     // * Store LastScanInMeas form SeqLoop *
    bool               bIsLastScanInConcat = GREKernel.ADC(0).Mdh.isLastScanInConcat();   // * Store LastScanInConcat form SeqLoop *

    bool               bToggleSpoiler   = false;    // * Spoiler gradients can be toggled to avoid artifacts    *
                                                    // * This option is used for quick sat techniques           *

    static long        lSatSpoilCounter = 0;

    static long        lPreviousIndex   = 0;  // * Reorder index of previous scan *
    static MSUProxy    theMSU;

    static long        lSignPE          = 1;
    long               lSignSS          = 1;
    long               lSpoilRO         = 0;

    if (eSequence == RT_Gre ) {

        // Tell the ADC whether we're swapped or not
        GREKernel.ADC(0).Mdh.setPRSwapped(RTControl.getbCurrentSwapState());

        // * ---------------------------------------------------------------------- *
        // * Activate WakeUp event near the end of a slice, if not called by SeqUT2 *
        // * ---------------------------------------------------------------------- *

        if ( ( lLine == maximum( lLinesToMeasure - 10L, 0L )) && (GREKernel.ADC(0).Mdh.getCacq() == pMrProt->averages()-1) && (! mIsUnittestActive()) ) {
            GREKernel.setbSendWakeUpBit ( true  );
        } else {
            GREKernel.setbSendWakeUpBit ( false );
        }

    }

    // * ---------------------------------------------------------------------- *
    // * Segments handling for various multi slice modes:                       *
    // *    sequential : segments are within TR                                 *
    // *    interleaved: segments are NOT within TR                             *
    // * ---------------------------------------------------------------------- *
    switch ( pMrProt->kSpace().multiSliceMode() ) {
        case SEQ::MSM_SEQUENTIAL:

                if ( (lKernelMode & KERNEL_STEADY_STATE_DUMMY_SCAN) == KERNEL_STEADY_STATE_DUMMY_SCAN     &&
// !!! REMOVE WITH WHEN THE NEW KERNEL MODE FOR STEADY STATE TRIGGERING IS AVAILABLE !!!
                     (PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION) && (PhysioMethodHigh != SEQ::METHOD_TRIGGERING) )  {
                    lSegmentsInTR = 1;
                } else {
                    lSegmentsInTR = lSegments;
                }
        break;

        case SEQ::MSM_INTERLEAVED:
                lSegmentsInTR = 1;
        break;

        default:
            TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence multi slice mode = %d  ", ptModule, pMrProt->kSpace().multiSliceMode());
            return ( SEQU_ERROR );

    }




    // * ---------------------------------------------------------------------- *
    // * Set origin code for unit test                                          *
    // * ---------------------------------------------------------------------- *
    if ( (lKernelMode & KERNEL_CHECK) == KERNEL_CHECK )  {
        GREKernel.setulUTIdent      ( RTEB_ORIGIN_fSEQCheck );
    } else {
        GREKernel.setulUTIdent      ( RTEB_ORIGIN_fSEQRunKernel );
    }





    // * UT exception for line / grid tag *
    #ifndef VXWORKS
        if( pMrProt->gridTag().on() || pMrProt->lineTag().on() )  {
            SeqUT.EnableTestCase(lTRClockErr, RTEB_ClockCheckTR);
        }
    #endif


    // * ---------------------------------------------------------------------- *
    // * Line / grid tag pulse                                                  *
    // * Executed only once for all segments, phases and slices                 *
    // * ---------------------------------------------------------------------- *
    if ( RunLoop.getFirstTimeInPhaseSliceLoop() )  {
        if(! SBBLineTag.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])  ) return (SBBLineTag.getNLSStatus());
        if(! SBBGridTag.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])  ) return (SBBGridTag.getNLSStatus());


        // * UT exception for line / grid tag *
        #ifndef VXWORKS
            if( pMrProt->gridTag().on() || pMrProt->lineTag().on() )  {
                SeqUT.DisableTestCase(lTRClockErr, RTEB_ClockCheckTR, "No common Tr with sequences with grid tag");
            }
        #endif

    }




    // * ------------------------------------------------------------------------- *
    // * Quick RSat: Perform a quick RSat after lNoOfKernelsBetweenQuickRSats      *
    // *             kernel calls                                                  *
    // * ------------------------------------------------------------------------- *
    if ( lNoOfKernelsBetweenQuickRSats == 0 ) {

        bPerformQuickRSat = true;

    } else {

        long lSliceNumber = RunLoop.getlSliceNumber( lKernelMode );

        if ( lSliceNumber % lNoOfKernelsBetweenQuickRSats ) {
             bPerformQuickRSat = false;
        } else {
             bPerformQuickRSat = true;
        }

    }





    // * ------------------------------------------------------------------------- *
    // * Quick FatSat: Perform a quick FatSat after                                *
    // *               lNoOfKernelsBetweenQuickFatSats kernel calls                *
    // * ------------------------------------------------------------------------- *
    if ( lNoOfKernelsBetweenQuickFatSats == 0 ) {

        bPerformQuickFatSat = true;

    } else {

        long lSliceNumber = RunLoop.getlSliceNumber( lKernelMode );

        if ( lSliceNumber % lNoOfKernelsBetweenQuickFatSats ) {
             bPerformQuickFatSat = false;
        } else {
             bPerformQuickFatSat = true;
        }

    }



    // * ------------------------------------------------------------------------- *
    // * Quick WatSater: Perform a quick WaterSat after                            *
    // *                 lNoOfKernelsBetweenQuickWaterSats kernel calls            *
    // * ------------------------------------------------------------------------- *
    if ( lNoOfKernelsBetweenQuickWaterSats == 0 ) {

        bPerformQuickWaterSat = true;

    } else {

        long lSliceNumber = RunLoop.getlSliceNumber( lKernelMode );

        if ( lSliceNumber % lNoOfKernelsBetweenQuickWaterSats ) {
             bPerformQuickWaterSat = false;
        } else {
             bPerformQuickWaterSat = true;
        }

    }



  // * --------------------------------------------------------------------------- *
  // * --------------------------------------------------------------------------- *
  // *                                                                             *
  // *           Execute the segments that are part of each TR intervall           *
  // *                                                                             *
  // * --------------------------------------------------------------------------- *
  // * --------------------------------------------------------------------------- *
  for ( lSegmentCounter=0; lSegmentCounter < lSegmentsInTR; lSegmentCounter++ )  {



    // * ---------------------------------------------------------------------- *
    // * WakeUp bit control for respiratory triggered Gre sequence              *
    // * ---------------------------------------------------------------------- *
    #ifdef VXWORKS
        if ( ( (eSequence == Gre) || (eSequence == Gre_HR) ) &&
             ( PhysioSignalHigh == SEQ::SIGNAL_RESPIRATION ) &&
             ( PhysioMethodHigh == SEQ::METHOD_TRIGGERING  )     ) {

            GREKernel.setbSendWakeUpBit ( false );

            if ( ( (lKernelMode & KERNEL_REAL_TIME) == KERNEL_REAL_TIME) && (lSegmentCounter == 0) )  {
                GREKernel.setbSendWakeUpBit ( true );
            }
        }
    #endif


    // * ---------------------------------------------------------------------- *
    // * RF spoiling control                                                    *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->fastImaging().RFSpoiling() )  {
      lStatus = fSUVerifyRFSpoil (
         pMrProt,                             // * IMP: user choice parameters  *
         pSeqLim,                             // * IMP: limits from fSEQInit()  *
         RunLoop.isbFirstSliceInConcat()?0:1, // * IMP: copy data   true/false  *
         lSlice,                              // * IMP: Anatomic slice number   *
        &dRFSpoilPhase,                       // * EXP: RF spoiling parameter   *
        &dRFSpoilIncrement,                   // * EXP: RF spoiling parameter   *
        &dRFSpoilPhasePrevSlice,              // * EXP: RF spoiling parameter   *
        &dRFSpoilIncrementPrevSlice,          // * EXP: RF spoiling parameter   *
        &dRFPrevSlicePosSag,                  // * EXP: RF spoiling parameter   *
        &dRFPrevSlicePosCor,                  // * EXP: RF spoiling parameter   *
        &dRFPrevSlicePosTra,                  // * EXP: RF spoiling parameter   *
        &dRFPrevSliceNormalSag,               // * EXP: RF spoiling parameter   *
        &dRFPrevSliceNormalCor,               // * EXP: RF spoiling parameter   *
        &dRFPrevSliceNormalTra                // * EXP: RF spoiling parameter   *
      );
      CheckStatusPB( lStatus,"fSUVerifyRFSpoiling");
    }




    // * ---------------------------------------------------------------------- *
    // * Set the phase for RF spoiling                                          *
    // * ---------------------------------------------------------------------- *
    if ( pMrProt->fastImaging().RFSpoiling() )  {
      dRFSpoilPhase     += (dRFSpoilIncrement += dRFSPOIL_INCREMENTdeg);
      dRFSpoilPhase      = fmod(dRFSpoilPhase    ,(double)RFMAXPHASEdeg);
      dRFSpoilIncrement  = fmod(dRFSpoilIncrement,(double)RFMAXPHASEdeg);
    }
    GREKernel.setdRFSpoilPhase ( dRFSpoilPhase );




    // * ---------------------------------------------------------------------- *
    // * Setting of the lines and partitions                                    *
    // * Setting of the phase and partition FFT flags                           *
    // * Setting of the first and last scan in slice flags                      *
    // * Setting of the segment counter (for phase correction)                  *
    // * ---------------------------------------------------------------------- *
    for ( lI=0; lI<pMrProt->contrasts(); lI++ )  {
        GREKernel.ADC(lI).Mdh.setCseg ( (short)lSegmentCounter );
    }
    Reorder.setIsLastAcquisition ( GREKernel.ADC(0).Mdh.getCacq() == pMrProt->averages()-1 );

    if ( (lKernelMode & KERNEL_PHASECOR) == KERNEL_PHASECOR)  {

        GREKernel.setlLineNumber ( lKSpaceCenterLine      );
        GREKernel.setlPartNumber ( lKSpaceCenterPartition );

    }  else  {

        if ( (lKernelMode & KERNEL_CHECK ) == KERNEL_CHECK ) {

            GREKernel.setlLineNumber ( lLine      );         // * No reordering in case of sequence check *
            GREKernel.setlPartNumber ( lPartition );         // * fSEQCheck sets the min. and max. values *

        } else {

            switch ( pMrProt->kSpace().dimension() ) {
                case SEQ::DIM_2:
                    lIndex = lLine * lSegmentsInTR + lSegmentCounter;
                break;

                case SEQ::DIM_3:
                    // ellipt.scanning and 3D-PAT require the reorder3dE mode
                    if ( pMrProt->kSpace().ellipticalScanning() || (   (pMrProt->PAT().PATMode() != SEQ::PAT_MODE_NONE)
                                                                    && (pMrProt->PAT().AccelFact3D() > 1)               )) {
                        lIndex = lLine;
                    } else {
                      lIndex = lLine * lPartitionsToMeasure + lPartition * lSegmentsInTR + lSegmentCounter;
                    }
                break;

                default:
                    TRACE_PUT2(TC_ALWAYS, TF_SEQ,"%s: Unknown sequence dimension = %d  ", ptModule, pMrProt->kSpace().dimension());
                    return ( SEQU_ERROR );
            }




            // * Use previous lIndex in case of staedy-state triggering *
            if ( (lKernelMode & KERNEL_STEADY_STATE_DUMMY_SCAN) == KERNEL_STEADY_STATE_DUMMY_SCAN)  {  // * Waiting for next trigger *
                lIndex = lPreviousIndex;
            } else {
                lPreviousIndex = lIndex;
            }

            GREKernel.setlLineNumber ( Reorder.getLinNo(lIndex) );
            GREKernel.setlPartNumber ( Reorder.getParNo(lIndex) );


            for ( lI=0; lI<pMrProt->contrasts(); lI++ )  {
                GREKernel.ADC(lI).Mdh.setPhaseFT     ( Reorder.isPhaseFT     ( lIndex ) );
                GREKernel.ADC(lI).Mdh.setPartitionFT ( Reorder.isPartitionFT ( lIndex ) );

                GREKernel.ADC(lI).Mdh.setFirstScanInSlice ( Reorder.isFirstScanInSlice ( lIndex ) );
                GREKernel.ADC(lI).Mdh.setLastScanInSlice  ( Reorder.isLastScanInSlice  ( lIndex ) );

                // * ---------------------------------------------------------------------- *
                // * Set PAT flags                                                          *
                // * ---------------------------------------------------------------------- *
                GREKernel.ADC(lI).Mdh.setPATRefScan       (Reorder.isPATRefScan(lIndex)      );
                GREKernel.ADC(lI).Mdh.setPATRefAndImaScan (Reorder.isPATRefAndImaScan(lIndex));

                // * -------------------------------------------------------------- *
                // * If LastScanInSlice is set by SeqLoop and last segment will     *
                // * be acquired ==> set LastScanInSlice flag in ADC                *
                // * -------------------------------------------------------------- *
                if ( bIsLastScanInMeas && (lSegmentCounter == lSegmentsInTR-1) ) {

                    if ( pMrProt->phaseStabilize() )  {
                        GREKernel.ADC(lI).Mdh.setLastScanInMeas( false );
                        GREKernel.setbLastScanInMeas ( true );
                    } else {
                        GREKernel.ADC(lI).Mdh.setLastScanInMeas( true );
                    }

                } else {
                    GREKernel.ADC(lI).Mdh.setLastScanInMeas( false );
                    GREKernel.setbLastScanInMeas ( false );
                }


                // * -------------------------------------------------------------- *
                // * If LastScanInConcat is set by SeqLoop and last segment will    *
                // * be acquired ==> set LastScanInConcat flag in ADC               *
                // * -------------------------------------------------------------- *
                if ( bIsLastScanInConcat && (lSegmentCounter == lSegmentsInTR-1)  ) {

                    if ( pMrProt->phaseStabilize() )  {
                        GREKernel.ADC(lI).Mdh.setLastScanInConcat( false );
                        GREKernel.setbLastScanInConcat ( true );
                    } else {
                        GREKernel.ADC(lI).Mdh.setLastScanInConcat( true  );
                    }

                } else {
                    GREKernel.ADC(lI).Mdh.setLastScanInConcat( false );
                    GREKernel.setbLastScanInConcat ( false );
                }

            }

        }

    }




    // * ---------------------------------------------------------------------- *
    // * Setting of the TRFill and TRFillEnd times                              *
    // * - TRFill is handled by the sequence if k-space is NOT segmented        *
    // * - TRFillEnd is handled by RunLoop                                      *
    // * ---------------------------------------------------------------------- *
    if ( lSegmentsInTR == 1 ) {
        GREKernel.setlTRFill ( RunLoop.getlTRFill() );
    } else {
        GREKernel.setlTRFill ( 0 );
    }
    GREKernel.setlTRFillEnd  ( 0 );





    // * ---------------------------------------------------------------------- *
    // * Execute sat pulses (only for the first segment)                        *
    // * ---------------------------------------------------------------------- *
    if ( lSegmentCounter == 0 )  {
        if(!MSat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])        ) return (MSat.getNLSStatus());

        for (lI= lNoOfRSats - 1; lI > -1; lI--)  { // * send RSats in inverse order *
            if ( SBBRSat[lI].isQuickMode() && !bIgnoreQuickRSat )  {
                if ( bPerformQuickRSat )  {
                    SBBRSat[lI].setPhaseInc ( dRFSpoilPhase );
                    if( !SBBRSat[lI].run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) ) return (SBBRSat[lI].getNLSStatus());
                }
            } else {
                SBBRSat[lI].setPhaseInc ( dRFSpoilPhase );
                if( !SBBRSat[lI].run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) ) return (SBBRSat[lI].getNLSStatus());
            }
        }
    }


    // * Perform one TSat pulse for each segment *
    SBBTSat.setPhaseInc ( dRFSpoilPhase );
    if ( !SBBTSat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])          ) return (SBBTSat.getNLSStatus()     );

    if ( lSegmentCounter == 0 )  {

        CSatFat.setPhaseInc   ( dRFSpoilPhase );
        CSatWater.setPhaseInc ( dRFSpoilPhase );

        if ( (pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK) && bPerformQuickFatSat ) {
            if ( !CSatFat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) ) return (CSatFat.getNLSStatus()  );
        } else {
            if ( pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION ) {
                if ( !CSatFat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) ) return (CSatFat.getNLSStatus()  );
            }
        }


        if ( ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK ) && bPerformQuickWaterSat )  {
            if ( !CSatWater.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) ) return (CSatWater.getNLSStatus());
        } else {
            if ( pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION ) {
                if ( !CSatWater.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) ) return (CSatWater.getNLSStatus());
            }
        }
    }


    // * ---------------------------------------------------------------------- *
    // * Determine whether the spoiler gradients should toggle.                 *
    // * Toggling is used for quick sat techniques.                             *
    // * ---------------------------------------------------------------------- *
    bToggleSpoiler =    ((pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK) && (! bIgnoreQuickFatSat) )
                     || ((pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK) && (! bIgnoreQuickWaterSat) )
                     || (( lNumberOfQuickPRSats > 0) && (! bIgnoreQuickRSat) );


    if ( bToggleSpoiler )  {
        if ( RunLoop.isbFirstSliceInConcat() )  {
            lSignPE  =  1;
            lSignSS  = -1;
            lSpoilRO =  1;
        } else {
            lSignSS  =  1;
            lSpoilRO =  0;
        }
    } else {
            lSignPE  =  1;
            lSignSS  =  1;
            lSpoilRO =  0;
    }

    // * If TSat is used perform one SatSpoil for each segment                                 *
    // * If Rsat, FatSat, WaterSat or MSat is used perform SatSpoil only for the first segment *
    // * If Quick PSat is used perform SatSpoil only for the first slice in  a concatenation   *
    if ( pMrProt->tSat().on() ||
         ( lSegmentCounter == 0 ) && ( bIgnoreQuickRSat || bIgnoreQuickFatSat || bIgnoreQuickWaterSat ||
         ( ( lNumberOfQuickPRSats > 0) && bPerformQuickRSat )  ||
         ( (pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK) && bPerformQuickFatSat ) ||
         ( (pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION_QUICK) && bPerformQuickWaterSat ) ||
         (lNumberOfRegularPRSats > 0) || (pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION) || (pMrProt->preparationPulses().waterSuppression() == SEQ::WATER_SATURATION) || pMrProt->preparationPulses().MTC() ) )  {


            SatSpoil.setPESpoilMomentConstant ( lSignPE * (dMaxSatSpoilMoment - dSatSpoilMomentStep * lSatSpoilCounter) );
            SatSpoil.setROSpoilMomentConstant ( lSpoilRO * dConstSatSpoilMoment );
            SatSpoil.setSSSpoilMomentConstant ( lSignSS * (dMinSatSpoilMoment + dSatSpoilMomentStep * lSatSpoilCounter) );


            lSatSpoilCounter += 2;

            if ( lSatSpoilCounter == lSpoilMomentIncr ) { lSatSpoilCounter = (lSpoilMomentIncr%2 ? 0 : 1); }
            if ( lSatSpoilCounter >  lSpoilMomentIncr ) { lSatSpoilCounter = (lSpoilMomentIncr%2 ? 1 : 0); }


        SatSpoil.setCalcMode      ( SPOIL_MODE_CONST  );
        SatSpoil.setAvailableTime ( lSatSpoilDuration );

        SatSpoil.prep ( pMrProt, pSeqLim, pSeqExpo );

        if ( !SatSpoil.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])  ) return (SatSpoil.getNLSStatus() );


        // * ------------------------------------------------------------------ *
        // * Toggle spoiler in phase encoding direction                         *
        // * ------------------------------------------------------------------ *
        if ( bToggleSpoiler )  {
            lSignPE *= -1;
        }

    }




    // * ---------------------------------------------------------------------- *
    // * Execute the GRE kernel                                                 *
    // * ---------------------------------------------------------------------- *
    if (! GREKernel.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )  {
      lStatus = GREKernel.getNLSStatus();

      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s GREKernel.run(...) return NLS error!", ptModule, lStatus);
      goto FINISHED;
    }



    if (eSequence == RT_Gre ) {
        if( GREKernel.ADC(0).Mdh.isLastScanInSlice() && (GREKernel.ADC(0).Mdh.getCacq() == pMrProt->averages()-1) )    {
            RTControl.calcAndSetNewOrientation(pMrProt, pSeqLim , asSLC , &(GREKernel.ADC(0)));
            GREKernel.updateRotationMatrix();
        }
    }

    #ifndef VXWORKS
      if ( lSegmentCounter == (pMrProt->fastImaging().segments() / 2 - 1) )  {
        mSEQTest(pMrProt, pSeqLim, pSeqExpo, RTEB_ClockCheckTI, 50, lLine, asSLC[lSlice].getSliceIndex(), 0, 0);
      }
    #endif

  }  // * End of segments loop *







  #ifndef VXWORKS
    mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ClockCheckTR,50,lLine,asSLC[lSlice].getSliceIndex(),0,0);

    if ( pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF )  {
      mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ClockCheckTI,50,lLine,asSLC[lSlice].getSliceIndex(),0,0);
    }
  #endif

  FINISHED:

  return(lStatus);
}

#ifdef SUPPORT_CT
//  Run kernel used if
NLS_STATUS fSEQRunKernel_ct
( MrProt          *pProt
, SeqLim          *pSeqLim
, SeqExpo         *pSeqExpo
, long            lKernelMode
, long            lSlice
, long            lClinInSeg
, long            lCSeg
)
{
    const int32_t  i32NSeg              = RunLoop.getLinesToMeasure()
        ,          i32NLinPerSeg        = RunLoop.getPartitionsToMeasure()
        ,          i32NPrepScans1stSeg  = i32NSeg*i32NLinPerSeg-lLinesToMeasure
        ,          i32NSlc              = pProt->sliceSeries().size()
        ;
    const int32_t   i32NEco              = pProt->contrasts();
    const bool bIsReadoutEnabled         = fRTIsReadoutEnabled();

    //  Fat-Saturation
    bool bExeCSat = pProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION;
    bool bExeRSat = lNumberOfRegularPRSats > 0;

    int iAdjUpdate = RunLoop.ct_adjUpdate(pProt);

    static uint32_t s_u32TimeSinceCSat = 0xffffffff; //  <limits> include doesn't work for i86

    if( pProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK )
    {
        if( pProt->MDS().mdsModeMask() == SEQ::MDS_ONCO )
        {
            bExeCSat = lCSeg == s_i32ZSeg-s_i32CSatShift;
        }
        else
        {

            //  Here we assume that the slice index is the excitation order
            //  This is not necessary the case, if multiple concatentions are used
            const int32_t i32NCSatPerTR = pProt->sliceSeries().mode() == SEQ::INTERLEAVED ? 2 : 1;
            bExeCSat = lSlice == 0
                ||    (i32NCSatPerTR == 2 && lSlice == (i32NSlc+1)/2)
                ;
        }
    }
    if( lNumberOfQuickPRSats > 0 )
    {
        if(  pProt->MDS().mdsModeMask() == SEQ::MDS_ONCO )
        {
            bExeRSat = lCSeg == s_i32ZSeg-s_i32CSatShift;
        }
        else
        {
            //  Here we assume that the slice index is the excitation order
            //  This is not necessary the case, if multiple concatentions are used
            const int32_t i32NCSatPerTR = pProt->sliceSeries().mode() == SEQ::INTERLEAVED ? 2 : 1;
            bExeRSat = lSlice == 0
                ||    (i32NCSatPerTR == 2 && lSlice == (i32NSlc+1)/2)
                ;
        }
    }
    if( bExeCSat )
    {
        s_u32TimeSinceCSat = 0;
    }
    else
    {
        ++s_u32TimeSinceCSat;
    }
    if( bExeRSat )
    {
        NLS_STATUS i32Status;
        SeqBuildBlockRSat* pRSat = SBBRSat+lNoOfRSats;
        while( pRSat-- != SBBRSat)
        {
            i32Status = SeqLoop::prepSLICE_POS(asSLC[lSlice],pProt,pSeqLim,pSeqExpo,RunLoop.getRTClock_us(),0,0);
            if( NLS_SEVERITY(i32Status) != NLS_SUCCESS )
            {
                TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
                return i32Status;
            }
            //  RF-Spoiling
            if( pProt->fastImaging().RFSpoiling() && pProt->MDS().mdsModeMask() == SEQ::MDS_OFF )
            {
                const int32_t i32j = RunLoop.getCTR();
                double dRFSpoilPhaseRSat = RFSPOIL_INCREMENTdeg*(i32j*i32j+i32j+2)/2.;
                dRFSpoilPhaseRSat = fmod(dRFSpoilPhaseRSat   ,(double)RFMAXPHASEdeg);
                pRSat->setPhaseInc( dRFSpoilPhaseRSat);
            }
            else
            {
                pRSat->setPhaseInc( 0 );
            }
            if( !bExeCSat && pRSat == SBBRSat )
            {
                pRSat->setRTLastRampDownOutsideSBB(true);
                GREKernel.setRampTimeofPreviousEB(pRSat->getRampTimeOutsideSBB());
            }
            pRSat->setAdjUpdate(iAdjUpdate);
            if( !pRSat->run(pProt,pSeqLim,pSeqExpo,&asSLC[lSlice]) )
            {
                TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
                return pRSat->getNLSStatus();
            }
            iAdjUpdate = pRSat->getAdjUpdate();
            if( !bExeCSat && pRSat == SBBRSat )
            {
                pRSat->setRTLastRampDownOutsideSBB(false);
            }
            RunLoop.addToRTClock_us(pRSat->getDurationPerRequest());
        }
    }  //  End ExeRSat
    if( bExeCSat )
    {
        if( GRAD_SPOIL_SCHEME(pProt) == GRAD_SPOIL_SCHEME_SIMPLE )
        {
            if( bExeRSat )
            {
                //  Add omitted spoil moment of last RSat
                if( !CSatFat.addRTSpoilMomFront((GP_RUT_RSAT+GPSpoil_DURATION_RSAT)*GPSpoil_AMPLITUDE_RSAT,(GP_RUT_RSAT+GPSpoil_DURATION_RSAT)*GPSpoil_AMPLITUDE_RSAT,(GP_RUT_RSAT+GPSpoil_DURATION_RSAT)*GPSpoil_AMPLITUDE_RSAT) )
                {
                    TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
                    return CSatFat.getNLSStatus();;
                }
            }
        }
        else
        {
            const bool bToggleSpoiler = pProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION_QUICK || lNumberOfQuickPRSats > 0;
            static long        lSatSpoilCounter = 0;
            static long        lSignPE          = 1;
            long               lSignSS          = 1;
            long               lSpoilRO         = 0;

            if ( bToggleSpoiler )
            {
                bool bFirstQuickSatInTR = RunLoop.isbFirstSliceInConcat();
                if(  pProt->MDS().mdsModeMask() == SEQ::MDS_ONCO )
                {
                    // RunLoop.isbFirstSliceInConcat means inner slice counter == 0
                    //  Here a saturation pulse is not necessarily played out before the first "Slice in Concat"
                    //  Therefore the condition is replaced by (First FS-Pulse in TR)
                    static int32_t s_i32CSatCntr = 0;
                    const int32_t i32NCSatPerTR = pProt->MDS().mdsTableSpeedNumerator();

                    bFirstQuickSatInTR = s_i32CSatCntr%i32NCSatPerTR == 0;
                    ++s_i32CSatCntr;
                }
                if ( RunLoop.isbFirstSliceInConcat() )
                {
                    lSignPE  =  1;
                    lSignSS  = -1;
                    lSpoilRO =  1;
                }
                else
                {
                    lSignSS  =  1;
                    lSpoilRO =  0;
                }
            }
            else
            {
                lSignPE  =  1;
                lSignSS  =  1;
                lSpoilRO =  0;
            }
            const double dMomPE = lSignPE  * (dMaxSatSpoilMoment - dSatSpoilMomentStep * lSatSpoilCounter);
            const double dMomRd = lSpoilRO * dConstSatSpoilMoment;
            const double dMomSl = lSignSS  * (dMinSatSpoilMoment + dSatSpoilMomentStep * lSatSpoilCounter);
            lSatSpoilCounter += 2;

            if ( lSatSpoilCounter == lSpoilMomentIncr ) { lSatSpoilCounter = (lSpoilMomentIncr%2 ? 0 : 1); }
            if ( lSatSpoilCounter >  lSpoilMomentIncr ) { lSatSpoilCounter = (lSpoilMomentIncr%2 ? 1 : 0); }

            if( !CSatFat.addRTSpoilMomBack(dMomPE,dMomRd,dMomSl) )
            {
                TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
                return CSatFat.getNLSStatus();
            }
            if ( bToggleSpoiler )
            {
                lSignPE *= -1;
            }
        }
        NLS_STATUS i32Status = SeqLoop::prepSLICE_POS(asSLC[lSlice],pProt,pSeqLim,pSeqExpo,RunLoop.getRTClock_us(),0,0);
        if( NLS_SEVERITY(i32Status) != NLS_SUCCESS )
        {
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
            return i32Status;
        }
        //  RF-Spoiling
        double dRFSpoilPhaseCSat = 0;
        if( pProt->fastImaging().RFSpoiling() && pProt->MDS().mdsModeMask() == SEQ::MDS_OFF  )
        {
            const int32_t i32j = RunLoop.getCTR();
            dRFSpoilPhaseCSat  = RFSPOIL_INCREMENTdeg*(i32j*i32j+i32j+2)/2.;
            dRFSpoilPhaseCSat  = fmod(dRFSpoilPhaseCSat,(double)RFMAXPHASEdeg);
        }
        CSatFat.setPhaseInc( dRFSpoilPhaseCSat );
        CSatFat.setRTLastRampDownOutsideSBB(true);
        CSatFat.setAdjUpdate(iAdjUpdate);
        if( !CSatFat.run(pProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
        {
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
            return CSatFat.getNLSStatus();
        }
        iAdjUpdate = CSatFat.getAdjUpdate();
        GREKernel.setRampTimeofPreviousEB(CSatFat.getRampTimeOutsideSBB());
        CSatFat.setRTLastRampDownOutsideSBB(false);

        RunLoop.addToRTClock_us(CSatFat.getDurationPerRequest());
    } //  END Exe CSat

    if( (lKernelMode & KERNEL_CHECK ) == KERNEL_CHECK )
    {
        GREKernel.setulUTIdent   ( RTEB_ORIGIN_fSEQCheck );
        GREKernel.setlLineNumber ( lCSeg      );
        GREKernel.setlPartNumber ( lClinInSeg );
    }
    else if( (lKernelMode & KERNEL_PHASECOR) == KERNEL_PHASECOR )
    {
        GREKernel.setulUTIdent   ( RTEB_ORIGIN_fSEQRunKernel );
        GREKernel.setlLineNumber ( lKSpaceCenterLine      );
        GREKernel.setlPartNumber ( lKSpaceCenterPartition );
    }
    else if( lCSeg == 0 && (lClinInSeg < i32NPrepScans1stSeg) )
    {
        GREKernel.setulUTIdent   ( RTEB_ORIGIN_fSEQRunKernel );
        //  Preparing scan 1st Line
        GREKernel.setlLineNumber ( lKSpaceCenterLine      );
        GREKernel.setlPartNumber ( lKSpaceCenterPartition );
        fRTSetReadoutEnable(0);
    }
    else
    {
        //  Imaging scan
        GREKernel.setulUTIdent   ( RTEB_ORIGIN_fSEQRunKernel );

        //  i) Retrieve reoder index
        long lIndex = 0;

        lIndex = lCSeg*i32NLinPerSeg+lClinInSeg-i32NPrepScans1stSeg;
        if( lIndex < 0 || lIndex >= lLinesToMeasure )
        {
            //  Something went wrong
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
            return SEQU_ERROR;
        }

        GREKernel.setlLineNumber ( Reorder.getLinNo(lIndex) );
        GREKernel.setlPartNumber ( Reorder.getParNo(lIndex) );

        //  ii) Complete MDH
        int32_t i32CEco = 0;
        for(;i32CEco != i32NEco;++i32CEco)
        {
            MdhProxy& rMDH = GREKernel.ADC(i32CEco).Mdh;
            rMDH.deleteFromEvalInfoMask(MDH_PHASESTABSCAN);
            rMDH.setPhaseFT     ( Reorder.isPhaseFT    ( lIndex ) );
            rMDH.setPartitionFT ( Reorder.isPartitionFT( lIndex ) );
            rMDH.setFirstScanInSlice ( lIndex == 0 );
            rMDH.setLastScanInSlice  ( lIndex == lLinesToMeasure-1 );
            rMDH.setPhaseFT     ( rMDH.isLastScanInSlice() );
            rMDH.setPartitionFT ( false ); // 2D

            rMDH.setIsMdsReferencePosition( Reorder.getLinNo(lIndex) == Reorder.getKSCenterLin() && i32CEco == 0 );

            if( pProt->MDS().mdsModeMask() == SEQ::MDS_ONCO )
            {
                rMDH.deleteFromEvalInfoMask(MDH_RTFEEDBACK);
                if( Reorder.getLinNo(lIndex) == Reorder.getKSCenterLin() )
                {
                    rMDH.addToEvalInfoMask(MDH_RTFEEDBACK);
                }
            }


            //  PAT flags
            rMDH.setPATRefScan      (Reorder.isPATRefScan(lIndex)      );
            rMDH.setPATRefAndImaScan(Reorder.isPATRefAndImaScan(lIndex));

            //  Last scan in meas/concat
            if( pProt->phaseStabilize() )
            {
                rMDH.setLastScanInMeas(false);
                GREKernel.setbLastScanInMeas( RunLoop.getbIsLastScanInMeas() );

                rMDH.setLastScanInConcat( false );
                GREKernel.setbLastScanInConcat( RunLoop.getbIsLastScanInConcat() );
            }
            else
            {
                rMDH.setLastScanInMeas  ( RunLoop.getbIsLastScanInMeas()   && i32CEco == i32NEco-1 );
                rMDH.setLastScanInConcat( RunLoop.getbIsLastScanInConcat() && i32CEco == i32NEco-1 );
            }
        }

    }

    //  TR-Fill
    GREKernel.setlTRFill ( RunLoop.getPerformTRFill() ? 0 : RunLoop.getlTRFill() );
    GREKernel.setlTRFillEnd( 0 );

    if( pProt->MDS().mdsModeMask() == SEQ::MDS_ONCO )
    {
        NLS_STATUS i32Status = RunLoop.prepSLICE_POS(asSLC[lSlice],pProt,pSeqLim,pSeqExpo,RunLoop.getRTClock_us(),0,&GREKernel.ADC(0));
        if( NLS_SEVERITY(i32Status) != NLS_SUCCESS )
        {
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
            return i32Status;
        }
    }

    //  RF-Spoiling
    if( pProt->fastImaging().RFSpoiling() && pProt->MDS().mdsModeMask() == SEQ::MDS_OFF  )
    {
        int32_t i32j = RunLoop.getCTR();
        dRFSpoilPhase = RFSPOIL_INCREMENTdeg*(i32j*i32j+i32j+2)/2.;
        dRFSpoilPhase = fmod(dRFSpoilPhase    ,(double)RFMAXPHASEdeg);
    }
    GREKernel.setdRFSpoilPhase ( dRFSpoilPhase );

    GREKernel.adjUpdate( iAdjUpdate );

    {
        //  Execute the Kernel
        if( !GREKernel.run(pProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
        {
            TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
            return GREKernel.getNLSStatus();
        }
        GREKernel.setRampTimeofPreviousEB(0);
    }
    iAdjUpdate = RT_MDSUPDATE_ADJ_NONE;

    RunLoop.addToRTClock_us(GREKernel.getDurationPerRequest()+GREKernel.getlTRFill());
    //  Restore readout state
    fRTSetReadoutEnable(bIsReadoutEnabled ? 1 : 0);

#ifndef VXWORKS
    mSEQTest(pProt,pSeqLim,pSeqExpo,RTEB_ClockCheckTR,50,lCSeg,asSLC[lSlice].getSliceIndex(),0,0);

    if ( pProt->preparationPulses().inversion() == SEQ::INVERSION_OFF )  {
        mSEQTest(pProt,pSeqLim,pSeqExpo,RTEB_ClockCheckTI,50,lCSeg,asSLC[lSlice].getSliceIndex(),0,0);
    }
#endif
    return SEQU_NORMAL;
}
#endif

/*---------------------------------------------------------------------------*/
/*  Copyright (C) Siemens AG 1998  All Rights Reserved.  Confidential        */
/*---------------------------------------------------------------------------*/
