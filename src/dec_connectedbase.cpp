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

/**@file   dec_connectedbase.cpp
 * @ingroup DETECTORS
 * @brief  detector connectedbase (completes the partialdec by bfs)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_connectedbase.h"
#include "cons_decomp.h"
#include "gcg.h"
#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

/* constraint handler properties */
#define DEC_DETECTORNAME          "connectedbase"       /**< name of detector */
#define DEC_DESC                  "detector connectedbase" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_DECCHAR               'C'         /**< display character of detector */
#define DEC_ENABLED               FALSE       /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      TRUE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE       /**< should the finishing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */
#define DEFAULT_USECONSSADJ       TRUE
/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
{
   SCIP_Bool useconssadj;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/**
 * @brief assigns all open constraints and open variables
 *
 *  strategy: assigns all conss and vars to the same block if they are connected,
 *  a cons and a var are adjacent if the var appears in the cons
 *
 *  @return scip return code
 */
static
SCIP_RETCODE completeByConnected(
   SCIP* scip,                      /**< scip data structure */
   gcg::PARTIALDECOMP* partialdec   /**< partialdecomp to complete */
   )
{
   int cons;
   int var;

   /* tools to check if the openvars can still be found in a constraint yet */
   std::vector<int> varinblocks; /* stores in which block the variable can be found */

   /* tools to update openvars */
   std::vector<int> openvarsToDelete( 0 );
   std::vector<int> oldOpenconss;

   int nconss = partialdec->getNConss();
   int nvars = partialdec->getNVars();
   
   std::vector<bool> isConsOpen( nconss, false );
   std::vector<bool> isConsVisited( nconss, false );

   std::vector<bool> isVarOpen( nvars, false );
   std::vector<bool> isVarVisited( nvars, false );

   std::queue<int> helpqueue = std::queue<int>();
   std::vector<int> neighborConss( 0 );
   std::vector<int> neighborVars( 0 );

   int nblocks = partialdec->getNBlocks();
   assert( (int) partialdec->getConssForBlocks().size() == nblocks );
   assert( partialdec->getNVarsForBlocks() == nblocks );
   assert( partialdec->getNTotalStairlinkingvars() == nblocks );

   SCIP_CALL( partialdec->refineToMaster( ) );

   if( nblocks < 0 )
   {
      nblocks = 0;
      partialdec->setNBlocks(0);
   }

   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();
   auto& openconss = partialdec->getOpenconssVec();
   auto& openvars = partialdec->getOpenvarsVec();

   /* initialize data structures */
   for( size_t c = 0; c < openconss.size(); ++ c )
   {
      cons = openconss[c];
      isConsOpen[cons] = true;
   }

   for( size_t v = 0; v < openvars.size(); ++ v )
   {
      var = openvars[v];
      isVarOpen[var] = true;
   }

   /* do breadth first search to find connected conss and vars */
   while( !openconss.empty() )
   {
      int newBlockNr;

      assert( helpqueue.empty() );
      helpqueue.push( openconss[0] );
      neighborConss.clear();
      neighborConss.push_back( openconss[0] );
      isConsVisited[openconss[0]] = true;
      neighborVars.clear();

      while( !helpqueue.empty() )
      {
         int nodeCons = helpqueue.front();
         assert( partialdec->isConsOpencons( nodeCons ) );
         helpqueue.pop();
         for( int v = 0; v < detprobdata->getNVarsForCons( nodeCons ); ++ v )
         {
            var = detprobdata->getVarsForCons( nodeCons )[v];
            assert( partialdec->isVarOpenvar( var ) || partialdec->isVarLinkingvar( var ) );

            if( isVarVisited[var] || partialdec->isVarLinkingvar( var ) )
               continue;

            for( int c = 0; c < detprobdata->getNConssForVar( var ); ++ c )
            {
               int otherNodeCons = detprobdata->getConssForVar( var )[c];
               if( !isConsOpen[otherNodeCons] || isConsVisited[otherNodeCons] )
               {
                  continue;
               }
               assert( partialdec->isConsOpencons( otherNodeCons ) );
               isConsVisited[otherNodeCons] = true;
               neighborConss.push_back( otherNodeCons );
               helpqueue.push( otherNodeCons );
            }
            isVarVisited[var] = true;
            neighborVars.push_back( var );
         }
      }

      /* assign found conss and vars to a new block */
      newBlockNr = partialdec->getNBlocks() + 1;
      partialdec->setNBlocks( newBlockNr );
      for( size_t i = 0; i < neighborConss.size(); ++ i )
      {
         cons = neighborConss[i];
         partialdec->setConsToBlock( cons, newBlockNr - 1 );
         if(partialdec->isConsOpencons(cons))
            partialdec->deleteOpencons( cons );
      }
      for( size_t i = 0; i < neighborVars.size(); ++ i )
      {
         var = neighborVars[i];
         partialdec->setVarToBlock( var, newBlockNr - 1 );
         if( partialdec->isVarOpenvar( var ) )
            partialdec->deleteOpenvar( var );
      }

      openconss = partialdec->getOpenconssVec();
   }

   /* assign left open vars to block 0, if it exists, and to master, otherwise */
   openvars = partialdec->getOpenvarsVec();
   for( size_t i = 0; i < openvars.size(); ++ i )
   {
      var = openvars[i];
      if( partialdec->getNBlocks() != 0 )
         partialdec->setVarToBlock( var, 0 );
      else
         partialdec->setVarToMaster( var );
      openvarsToDelete.push_back( var );
   }

   for( size_t i = 0; i < openvarsToDelete.size(); ++ i )
   {
      var = openvarsToDelete[i];
      if( partialdec->isVarOpenvar( var ) )
         partialdec->deleteOpenvar( var );
   }

   assert( partialdec->getNOpenconss() == 0 );
   assert( partialdec->getNOpenvars() == 0 );

   partialdec->prepare();

   assert( partialdec->checkConsistency( ) );

   return SCIP_OKAY;
}


/**
* @brief assigns all open constraints and open variables
*
*  strategy: assigns all conss and vars to the same block if they are connected
*  a cons and a var are adjacent if the var appears in the cons
*  \note this relies on the consadjacency structure of the detprobdata
*  hence it cannot be applied in presence of linking variables
*/
static
 void completeByConnectedConssAdjacency(
   SCIP* scip,                      /**< scip data structure */
   gcg::PARTIALDECOMP* partialdec   /**< partialdecomp to complete */
   )
{
    int cons;
    int var;

    /* tools to check if the openvars can still be found in a constraint yet */
    std::vector<int> varinblocks; /* stores in which block the variable can be found */

    /* tools to update openvars */
    std::vector<int> oldOpenconss;
    std::vector<int> openvarsToDelete;

   // note: this should not happen
    if( partialdec->getNLinkingvars() != 0 )
      completeByConnected(scip, partialdec);

   int nconss = partialdec->getNConss();
   int nvars = partialdec->getNVars();

    std::vector<bool> isConsOpen( nconss, false );
    std::vector<bool> isConsVisited( nconss, false );

    varinblocks = std::vector<int>(nvars, -1);

    std::queue<int> helpqueue = std::queue<int>();
    std::vector<int> neighborConss( 0 );

   int nblocks = partialdec->getNBlocks();
   assert( (int) partialdec->getConssForBlocks().size() == nblocks );
   assert( partialdec->getNVarsForBlocks() == nblocks );
   assert( partialdec->getNTotalStairlinkingvars() == nblocks );

   partialdec->refineToMaster();

   assert(partialdec->checkConsistency() );
   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   if( nblocks < 0 )
   {
      nblocks = 0;
      partialdec->setNBlocks(0);
   }

   auto& openconss = partialdec->getOpenconssVec();

   for( size_t c = 0; c < openconss.size(); ++ c )
    {
       cons = openconss[c];
       isConsOpen[cons] = true;
    }

    /* do breadth first search to find connected conss */
    while( !openconss.empty() )
    {
       int newBlockNr;

       assert( helpqueue.empty() );
       helpqueue.push( openconss[0] );
       neighborConss.clear();
       neighborConss.push_back( openconss[0] );
       isConsVisited[openconss[0]] = true;

       while( !helpqueue.empty() )
       {
          int nodeCons = helpqueue.front();
          assert( partialdec->isConsOpencons( nodeCons ) );
          helpqueue.pop();
          for( int c = 0; c < detprobdata->getNConssForCons( nodeCons ); ++c )
          {
             int othercons = detprobdata->getConssForCons( nodeCons )[c];

             if( isConsVisited[othercons] || partialdec->isConsMastercons( othercons ) || !isConsOpen[othercons] )
                continue;

             assert( partialdec->isConsOpencons( othercons ) );
             isConsVisited[othercons] = true;
             neighborConss.push_back( othercons );
             helpqueue.push( othercons );
          }
       }

       /* assign found conss and vars to a new block */
       newBlockNr = partialdec->getNBlocks() + 1;
       partialdec->setNBlocks( newBlockNr );
       for( size_t i = 0; i < neighborConss.size(); ++ i )
       {
          cons = neighborConss[i];

          partialdec->setConsToBlock( cons, newBlockNr - 1 );
          if(partialdec->isConsOpencons(cons))
            partialdec->deleteOpencons( cons );

          for( int j = 0; j < detprobdata->getNVarsForCons(cons); ++ j )
          {
             int newvar = detprobdata->getVarsForCons(cons)[j];

             if( partialdec->isVarLinkingvar(newvar) || varinblocks[newvar] != -1 )
                continue;

             assert(!partialdec->isVarMastervar( newvar) );
             partialdec->setVarToBlock( newvar, newBlockNr - 1 );
             varinblocks[newvar] = newBlockNr - 1;
             if( partialdec->isVarOpenvar(newvar) )
                partialdec->deleteOpenvar( newvar );
          }
       }
       openconss = partialdec->getOpenconssVec();
    }

    /* assign left open vars to block 0, if it exists, and to master, otherwise */
    auto openvars = partialdec->getOpenvarsVec();
    for( size_t i = 0; i < openvars.size(); ++ i )
    {
       var = openvars[i];
       if( partialdec->getNBlocks() != 0 )
          partialdec->setVarToBlock( var, 0 );
       else
          partialdec->setVarToMaster( var );
       openvarsToDelete.push_back( var );
    }

    for( size_t i = 0; i < openvarsToDelete.size(); ++ i )
    {
       var = openvarsToDelete[i];
      if( partialdec->isVarOpenvar( var ) )
         partialdec->deleteOpenvar( var );
    }

    assert( partialdec->getNOpenconss() == 0 );
    assert( partialdec->getNOpenvars() == 0 );

    partialdec->prepare();

    assert( partialdec->checkConsistency( ) );
 }


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
/** destructor of detector to free detector data (called when SCIP is exiting) */
static
DEC_DECL_FREEDETECTOR(freeConnectedbase)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}


/** destructor of detector to free detector data (called before the solving process begins) */
#define exitConnectedbase NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initConnectedbase NULL

#define propagatePartialdecConnectedbase NULL


static
DEC_DECL_FINISHPARTIALDEC(finishPartialdecConnectedbase)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
   char decinfo[SCIP_MAXSTRLEN];

   SCIP_Bool byconssadj;
   SCIP_Bool conssadjcalculated;

   gcg::PARTIALDECOMP* partialdec;
   partialdec = new gcg::PARTIALDECOMP(partialdecdetectiondata->workonpartialdec);

   SCIPgetBoolParam(scip, "detection/detectors/connectedbase/useconssadj", &byconssadj);
   conssadjcalculated = GCGconshdlrDecompGetConssAdjCalculated(scip);
   //complete the partialdec by bfs

   if( byconssadj && conssadjcalculated && partialdec->getNLinkingvars() == 0 )
      completeByConnectedConssAdjacency(scip, partialdec);
   else
      completeByConnected(scip, partialdec);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(scip, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "connected");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);

   partialdecdetectiondata->newpartialdecs[0]->addClockTime( SCIPgetClockTime(scip, temporaryClock)  );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


#define detectorPostprocessPartialdecConnectedbase NULL


static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveConnectedbase)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);


   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE ) );


   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultConnectedbase)
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
DEC_DECL_SETPARAMFAST(setParamFastConnectedbase)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE ) );

   return SCIP_OKAY;

}


/*
 * detector specific interface methods
 */

/** creates the handler for connectedbase detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorConnectedbase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create connectedbase detector data here*/
   detectordata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   detectordata->useconssadj = TRUE;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeConnectedbase, initConnectedbase, exitConnectedbase, propagatePartialdecConnectedbase, finishPartialdecConnectedbase, detectorPostprocessPartialdecConnectedbase, setParamAggressiveConnectedbase, setParamDefaultConnectedbase, setParamFastConnectedbase) );

   /* add consname detector parameters */
      /**@todo add connectedbase detector parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/connectedbase/useconssadj", "should the constraint adjacency be used", &detectordata->useconssadj, FALSE, DEFAULT_USECONSSADJ, NULL, NULL) );

   return SCIP_OKAY;
}
