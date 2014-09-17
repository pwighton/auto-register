/*[ Compilation unit ********************************************************\
*
* Project: NUMARIS/4
*
* File        : \n4\pkg\MrServers\MrImaging\seq\a_gre\LocalSeqLoop.h
*
* Author      : koellner
*
* Date        : n.a. 
*
* Lang        : c++
*
*
* Description : Definition of class LocalSeqLoop
*
\****************************************************************************/

#ifndef __GRE_LOCAL_SEQ_LOOP_H
#define __GRE_LOCAL_SEQ_LOOP_H

// * ------------------------------------------------------------------------ *
// * Base class definition                                                    *
// * ------------------------------------------------------------------------ *

#ifdef SUPPORT_PACE
    #include "MrServers/MrImaging/seq/common/libPACE/SLFB.h"
    typedef SLFB LOCAL_SEQ_LOOP_BASE_TYPE;
#else
    #include "MrServers/MrImaging/libSBB/SEQLoop.h"
    typedef SeqLoop LOCAL_SEQ_LOOP_BASE_TYPE;
#endif


// * ------------------------------------------------------------------------ *
// * Include for local use of SBBSteadyStateTrigger                           *
// * ------------------------------------------------------------------------ *
#include "MrServers/MrImaging/libSBB/SBBSteadyStateTrigger.h"

// * ------------------------------------------------------------------------ *
// * Definition of class LocalSeqLoop                                         *
// * ------------------------------------------------------------------------ *

class LocalSeqLoop : public LOCAL_SEQ_LOOP_BASE_TYPE
{
    public:

        typedef LOCAL_SEQ_LOOP_BASE_TYPE BASE_TYPE;

        LocalSeqLoop();

        long getlSliceNumber(long lKernelMode);

        virtual bool performStandardECG (MrProt*, SeqLim*);

        inline void setlLinesPerTrigger (long lLines) {  m_lLinesPerTrigger = lLines;  };

    protected:

        long m_lLinesPerTrigger;

};


#endif  //  __GRE_LOCAL_SEQ_LOOP_H
