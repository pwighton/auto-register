/*[ Compilation unit ********************************************************\
*
* Project: NUMARIS/4
*
* File        : \n4\pkg\MrServers\MrImaging\seq\a_gre\LocalSeqLoop.cpp
*
* Author      : koellner
*
* Date        : n.a. 
*
* Lang        : c++
*
*
* Description : Implementation of LocalSeqLoop member functions.
*
\****************************************************************************/


// * ------------------------------------------------------------------------ *
// * Class definition                                                         *
// * ------------------------------------------------------------------------ *

#include "MrServers/MrImaging/seq/a_gre/LocalSeqLoop.h"

// * ------------------------------------------------------------------------ *
// * Additional includes (also included by original SeqLoop implementation)   *
// * ------------------------------------------------------------------------ *

#include  "MrServers/MrImaging/libSeqUtil/libSeqUtil.h"
#include  "MrServers/MrMeasSrv/SeqIF/csequence.h"
#include  "MrServers/MrImaging/ut/libSeqUT.h"                // for mSEQTest


// * ------------------------------------------------------------------------ *
// * Implementation of LocalSeqLoop members                                   *
// * ------------------------------------------------------------------------ *

#define DEBUG_ORIGIN  DEBUG_SEQLOOP

LocalSeqLoop::LocalSeqLoop() : 
    m_lLinesPerTrigger(1),
    BASE_TYPE() {}






bool LocalSeqLoop::performStandardECG (MrProt* pProt, SeqLim* pSeqLim)
{

    if ( !(m_RawLinesCount % (m_lInnerSliceNumber * m_lLinesPerTrigger) ) ) {
        return BASE_TYPE::performStandardECG(pProt,pSeqLim);
    } else {
        return ( false );
    }

}


long LocalSeqLoop::getlSliceNumber (long lKernelMode)
{

    if ( lKernelMode & KERNEL_STEADY_STATE_DUMMY_SCAN ) {
        return ( m_SBBSteadyStateTrigger.getlSliceLoopCounter() );  // * Kernel is called for steady-state dummy scans
    } else {
        return ( m_lInnerSliceCounter );  // * Kernel is called for imaging *
    }

}






