/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2019 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_compgreedily.cpp
 * @ingroup DETECTORS
 * @brief  detector compgreedily (assigns the open cons and open vars of the partialdec greedily)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_compgreedily.h"
#include "cons_decomp.h"
#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "scip/clock.h"
#include <iostream>

/* constraint handler properties */
#define DEC_DETECTORNAME          "compgreedily"       /**< name of detector */
#define DEC_DESC                  "detector compgreedily" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'g'         /**< display character of detector */
#define DEC_ENABLED               FALSE       /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE       /**< should the finishing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */

/** parameter limits for emphasis default */

#define DEFAULT_LIMITHALFPERIMETERENABLEDFINISHING    20000   /** limit in terms of nrows + ncols for enabling finishing */
#define DEFAULT_LIMITHALFPERIMETERENABLEDORIGINAL     10000   /** limit in terms of nrows + ncols for enabling in detecting for unpresolved problem */


/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/**
 * @brief assigns all open constraints and open variables
 *
 *  strategy: assigns a cons (and related vars) to a new block if possible,
 *  if not to an existing block if possible (by means of prior var assignments)
 *  and finally to master, if there does not exist such a block
 */
static
void completeGreedily(
   gcg::PARTIALDECOMP* partialdec      /** partialdec to complete */
   )
{
   bool checkvar;
   bool isvarinblock;
   bool notassigned;

   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   /* tools to check if the openvars can still be found in a constraint yet*/
   std::vector<int> varinblocks; /* stores in which block the variable can be found */

   if( partialdec->getNBlocks() == 0 && partialdec->getNOpenconss() > 0 )
   {
      int block = partialdec->addBlock();
      std::vector<int>& openconss = partialdec->getOpenconssVec();
      partialdec->fixConsToBlock( openconss[0], block );
   }
   
   std::vector<int> del;

   /* check if the openvars can already be found in a constraint */
   std::vector<int> openvars = partialdec->getOpenvarsVec();
   for( int i = 0; i < partialdec->getNOpenvars(); ++ i )
   {
      varinblocks.clear();

      /* test if the variable can be found in blocks */
      for( int b = 0; b < partialdec->getNBlocks(); ++ b )
      {
         isvarinblock = false;
         std::vector<int>& conssforblock = partialdec->getConssForBlock(b);
         for( int k = 0; k < partialdec->getNConssForBlock(b) && !isvarinblock; ++ k )
         {
            for( int l = 0; l < detprobdata->getNVarsForCons( conssforblock[k] ); ++ l )
            {
               if( openvars[i] == detprobdata->getVarsForCons( conssforblock[k] )[l] )
               {
                  varinblocks.push_back( b );
                  isvarinblock = true;
                  break;
               }
            }
         }
      }
      if( varinblocks.size() == 1 ) /* if the variable can be found in one block set the variable to a variable of the block*/
      {
         partialdec->setVarToBlock(openvars[i], varinblocks[0]);
         del.push_back(openvars[i]);
         continue; /* the variable doesn't need to be checked any more */
      }
      else if( varinblocks.size() == 2 ) /* if the variable can be found in two blocks check if it is a linking var or a stairlinking var*/
      {
         if( varinblocks[0] + 1 == varinblocks[1] )
         {
            partialdec->setVarToStairlinking(openvars[i], varinblocks[0], varinblocks[1]);
            del.push_back(openvars[i]);
            continue; /* the variable doesn't need to be checked any more */
         }
         else
         {
            partialdec->setVarToLinking(openvars[i]);
            del.push_back(openvars[i]);
            continue; /* the variable doesn't need to be checked any more */
         }
      }
      else if( varinblocks.size() > 2 ) /* if the variable can be found in more than two blocks it is a linking var */
      {
         partialdec->setVarToLinking(openvars[i]);
            del.push_back(openvars[i]);
         continue; /* the variable doesn't need to be checked any more */
      }

      checkvar = true;

      /* if the variable can be found in an open constraint it is still an open var */
      for( int j = 0; j < partialdec->getNOpenconss(); ++ j )
      {
         checkvar = true;
         for( int k = 0; k < detprobdata->getNVarsForCons( j ); ++ k )
         {
            if( openvars[i] == detprobdata->getVarsForCons( j )[k] )
            {
               checkvar = false;
               break;
            }
         }
         if( ! checkvar )
         {
            break;
         }
      }

      /* test if the variable can be found in a master constraint yet */
        for( int k = 0; k < detprobdata->getNConssForVar( openvars[i] ) && checkvar; ++ k )
        {
           if( partialdec->isConsMastercons(detprobdata->getConssForVar(openvars[i])[k]) )
           {
              partialdec->setVarToMaster(openvars[i]);
              del.push_back(openvars[i]);
              checkvar = false; /* the variable does'nt need to be checked any more */
              break;
           }
        }
   }

   /* remove assigned vars from list of open vars */
   for(auto v : del)
      partialdec->deleteOpenvar(v);

   del.clear();
   partialdec->sort();

   std::vector<int> delconss;
   std::vector<int>& openconss = partialdec->getOpenconssVec();

   /* assign open conss greedily */
   for( int i = 0; i < partialdec->getNOpenconss(); ++ i )
   {
      std::vector<int> vecOpenvarsOfBlock; /* stores the open vars of the blocks */
      bool consGotBlockcons = false; /* if the constraint can be assigned to a block */

      /* check if the constraint can be assigned to a block */
      for( int j = 0; j < partialdec->getNBlocks(); ++ j )
      {
         /* check if all vars of the constraint are a block var of the current block, an open var, a linkingvar or a mastervar*/
         consGotBlockcons = true;
         for( int k = 0; k < detprobdata->getNVarsForCons( openconss[i] ); ++ k )
         {
            if( partialdec->isVarBlockvarOfBlock( detprobdata->getVarsForCons( openconss[i] )[k], j )
               || partialdec->isVarOpenvar( detprobdata->getVarsForCons( openconss[i] )[k] )
               || partialdec->isVarLinkingvar( detprobdata->getVarsForCons( openconss[i] )[k] )
               || partialdec->isVarStairlinkingvarOfBlock( detprobdata->getVarsForCons( openconss[i] )[k], j )
               || ( j != 0 && partialdec->isVarStairlinkingvarOfBlock( detprobdata->getVarsForCons( openconss[i] )[k], j - 1 ) ) )
            {
               if( partialdec->isVarOpenvar( detprobdata->getVarsForCons( openconss[i] )[k] ) )
               {
                  vecOpenvarsOfBlock.push_back( detprobdata->getVarsForCons( openconss[i] )[k] );
               }
            }
            else
            {
               vecOpenvarsOfBlock.clear(); /* the open vars don't get vars of the block */
               consGotBlockcons = false; /* the constraint can't be constraint of the block, check the next block */
               break;
            }
         }
         if( consGotBlockcons ) /* the constraint can be assigned to the current block */
         {
            partialdec->setConsToBlock( openconss[i], j );
            delconss.push_back(openconss[i]);
            for( size_t k = 0; k < vecOpenvarsOfBlock.size(); ++ k ) /* the openvars in the constraint get block vars */
            {
               partialdec->setVarToBlock( vecOpenvarsOfBlock[k], j );
               partialdec->deleteOpenvar( vecOpenvarsOfBlock[k] );
            }
            vecOpenvarsOfBlock.clear();

            break;
         }
      }

      if( !consGotBlockcons ) /* the constraint can not be assigned to a block, set it to master */
      {
         partialdec->setConsToMaster( openconss[i] );
         delconss.push_back(openconss[i]);
      }
   }

   /* remove assigned conss from list of open conss */
   for(auto c : delconss)
      partialdec->deleteOpencons(c);

   partialdec->sort();

   /* assign open vars greedily */
   for( int i = 0; i < partialdec->getNOpenvars(); ++ i )
   {
      notassigned = true;
      std::vector<int>& masterconss = partialdec->getMasterconss();
      for( int j = 0; j < partialdec->getNMasterconss() && notassigned; ++ j )
      {
         for( int k = 0; k < detprobdata->getNVarsForCons( masterconss[j] ); ++ k )
         {
            if( openvars[i] == detprobdata->getVarsForCons( masterconss[j] )[k] )
            {
               partialdec->setVarToMaster(openvars[i]);
               del.push_back(openvars[i]);
               notassigned = false;
               break;
            }
         }
      }
   }

   /* remove assigned vars from list of open vars */
   for(auto v : del)
      partialdec->deleteOpenvar(v);

   partialdec->sort();

   /* check if the open conss are all assigned */
   assert( partialdec->checkAllConssAssigned() );

   /* check if the open vars are all assigned */
   assert( partialdec->getNOpenvars() == 0 );

   assert( partialdec->checkConsistency( ) );
}


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */

#define freeCompgreedily NULL

/** destructor of detector to free detector data (called before the solving process begins) */

#define exitCompgreedily NULL

#define initCompgreedily NULL

static
DEC_DECL_PROPAGATEPARTIALDEC(propagatePartialdecCompgreedily)
{
   *result = SCIP_DIDNOTFIND;

   char decinfo[SCIP_MAXSTRLEN];
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;

   //assign open conss and vars greedily
   completeGreedily(partialdec);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock) );

   partialdecdetectiondata->detectiontime =  SCIPgetClockTime(scip, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;

   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(scip, temporaryClock));
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "compgreed");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

static
DEC_DECL_FINISHPARTIALDEC(finishPartialdecCompgreedily)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   gcg::PARTIALDECOMP* partialdec = new gcg::PARTIALDECOMP(partialdecdetectiondata->workonpartialdec);

   //assign open conss and vars greedily
   completeGreedily(partialdec);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock) );

   partialdecdetectiondata->detectiontime =  SCIPgetClockTime(scip, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "compgreed");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);

   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(scip, temporaryClock));

   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define detectorPostprocessPartialdecCompgreedily NULL


static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveCompgreedily)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   return SCIP_OKAY;
}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultCompgreedily)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDFINISHING ) );

   return SCIP_OKAY;
}

static
DEC_DECL_SETPARAMFAST(setParamFastCompgreedily)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   return SCIP_OKAY;

}


/*
 * detector specific interface methods
 */

/** creates the handler for compgreedily detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorCompgreedily(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create compgreedily detector data here*/
   detectordata = NULL;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeCompgreedily,initCompgreedily, exitCompgreedily, propagatePartialdecCompgreedily, finishPartialdecCompgreedily, detectorPostprocessPartialdecCompgreedily, setParamAggressiveCompgreedily, setParamDefaultCompgreedily, setParamFastCompgreedily) );

   /**@todo add compgreedily detector parameters */

   return SCIP_OKAY;
}
