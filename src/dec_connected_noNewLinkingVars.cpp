/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   dec_connected_noNewLinkingVars.c
 * @ingroup DETECTORS
 * @brief  detector connected_noNewLinkingVars (put your description here)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_connected_noNewLinkingVars.h"
#include "cons_decomp.h"
#include "gcg.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

/* constraint handler properties */
#define DEC_DETECTORNAME         "connected_noNewLinkingVars"       /**< name of detector */
#define DEC_DESC                 "detector connected_noNewLinkingVars" /**< description of detector*/
#define DEC_PRIORITY             0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              '?'         /**< display character of detector */
#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#define DEC_SKIP                 FALSE       /**< should detector be skipped if other detectors found decompositions */

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


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeConnected_noNewLinkingVars)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}


/** destructor of detector to free detector data (called before the solving process begins) */
#if 0
static
DEC_DECL_EXITDETECTOR(exitConnected_noNewLinkingVars)
{  /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitConnected_noNewLinkingVars NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initConnected_noNewLinkingVars)
{  /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initConnected_noNewLinkingVars NULL
#endif

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectConnected_noNewLinkingVars)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();  /*lint --e{527}*/

   return SCIP_OKAY;
}


static
bool haveConssCommonVars(
   int               firstCons,
   int               secondCons,
   gcg::Seeedpool*   seeedpool
   )
{
   for( int i = 0; i < seeedpool->getNVarsForCons(firstCons); ++i )
   {
      for( int j = 0; j < seeedpool->getNVarsForCons(secondCons); ++j )
      {
         if( seeedpool->getVarsForCons(firstCons)[i] == seeedpool->getVarsForCons(secondCons)[j])
         {
            return true;
         }
      }
   }
   return false;
}

/** Breadth First Search */
static
SCIP_RETCODE bfs(
   std::vector<int>      *visited,               /** vector to store the visited conss */
   std::vector<int>      *openConss,             /** vector with conss to be visited */
   gcg::Seeedpool*       seeedpool
   )
{
   std::queue<int> queue;
   std::vector<int> neighborNodes;
   int cons;

   std::vector<int>::iterator varIter;
   std::vector<int>::iterator varIterEnd;
   std::vector<int>::iterator it;

   queue.push( *(openConss->begin()) );
   openConss->erase(openConss->begin());
   while(!queue.empty())
   {
      cons = queue.front();
      visited->push_back(cons);
      queue.pop();
      varIter = openConss->begin();
      varIterEnd = openConss->end();

      for(; varIter != varIterEnd; ++varIter)
      {
         if(haveConssCommonVars(cons, *varIter, seeedpool))
         {
            queue.push(*varIter);
            neighborNodes.push_back(*varIter);
         }
      }

      for( size_t s = 0; s < neighborNodes.size(); ++s )
      {
         it = find(openConss->begin(), openConss->end(), neighborNodes[s]);
         assert( it != openConss->end() );
         openConss->erase(it);
      }

      neighborNodes.clear();
   }

   return SCIP_OKAY;
}

/** assigns openVars to Stairlinking if they can be found in two consecutive  blocks*/
static
SCIP_RETCODE considerImplicitsWithStairlinkingvars(
   gcg::Seeed*           seeed,
   gcg::Seeedpool*       seeedpool
   )
{
   std::vector<int> blocksOfOpenvar;
   std::vector<int> assignedOpenvars;
   bool foundInBlock;
   int var;
   int cons;

   seeed->considerImplicits(seeedpool);

   /** set vars to linking, if they can be found in more than one block */
   for(int i = 0; i < seeed->getNOpenvars(); ++i)
   {
      blocksOfOpenvar.clear();
      foundInBlock = false;
      var = seeed->getOpenvars()[i];
      for(int b = 0; b < seeed->getNBlocks(); ++b)
      {
         for(int c = 0; c < seeed->getNConssForBlock(b) && !foundInBlock; ++c)
         {
            cons = seeed->getConssForBlock(b)[c];
            for(int v = 0; v < seeedpool->getNVarsForCons(cons) && !foundInBlock; ++v)
            {
               if(seeedpool->getVarsForCons(cons)[v] == var)
               {
                  blocksOfOpenvar.push_back(b);
               }
            }
         }
      }
      if(blocksOfOpenvar.size() == 2 && blocksOfOpenvar[0] + 1 == blocksOfOpenvar[1])
      {
         seeed->setVarToStairlinking(var, blocksOfOpenvar[0], blocksOfOpenvar[1]);
         assignedOpenvars.push_back(var);
      }
   }

   for(size_t i = 0; i < assignedOpenvars.size(); ++i)
   {
      seeed->deleteOpenvar(assignedOpenvars[i]);
   }
   return SCIP_OKAY;
}

static
DEC_DECL_PROPAGATESEEED(propagateSeeedConnected_noNewLinkingVars)
{
   *result = SCIP_DIDNOTFIND;
   int cons;
   int var;
   bool assigned;
   std::vector<int> newOpenConss;
   std::vector<int> conssForBfs;
   std::vector<int> assignedConssForBfs;
   std::vector<int>::iterator it;
   std::vector<std::vector<int>> visitedConss = std::vector<std::vector<int>>(0); /** vector of vector with connected constraints */
   std::vector<int> emptyVector = std::vector<int>(0);
   int newBlocks = 0;
   int block;

   seeedPropagationData->seeedToPropagate->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
   gcg::Seeed* seeed;
   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);

   if(!seeed->areOpenVarsAndConssCalculated())
   {
      seeed->calcOpenconss();
      seeed->calcOpenvars();
      seeed->setOpenVarsAndConssCalculated(true);
   }

   considerImplicitsWithStairlinkingvars(seeed, seeedPropagationData->seeedpool);
   assert(seeedPropagationData->seeedToPropagate->checkConsistency());

   /** conss with stairlinkingvars are still openconss */
   for(int i = 0; i < seeed->getNOpenconss(); ++i)
   {
      assigned = false;
      cons = seeed->getOpenconss()[i];
      for(int v = 0; v < seeedPropagationData->seeedpool->getNVarsForCons(cons) && !assigned; ++v)
      {
         var = seeedPropagationData->seeedpool->getVarsForCons(cons)[v];
         for(int b = 0; b < seeed->getNBlocks() && !assigned; ++b)
         {
            if(seeed->isVarStairlinkingvarOfBlock(var, b))
            {
               newOpenConss.push_back(cons);
               assigned = true;
            }
         }
      }
   }

   for(int i = 0; i < seeed->getNOpenconss(); ++i)
   {
      cons = seeed->getOpenconss()[i];
      if(find(newOpenConss.begin(), newOpenConss.end(), cons) != newOpenConss.end())
         conssForBfs.push_back(cons);
   }

   /** assign the open conss which have common open vars with blockconss and assign the var */
   for(size_t i = 0; i < conssForBfs.size(); ++i)
   {
      assigned = false;
      cons = conssForBfs[i];
      for(int j = 0; j < seeedPropagationData->seeedpool->getNVarsForCons(cons) && !assigned; ++j)
      {
         var = seeedPropagationData->seeedpool->getVarsForCons(cons)[j];
         if(!seeed->isVarOpenvar(var))
            continue;
         for(int b = 0; b < seeed->getNBlocks() && !assigned; ++b)
         {
            for(int c = 0; c < seeed->getNConssForBlock(b) && !assigned; ++c)
            {
               for(int v = 0; v < seeedPropagationData->seeedpool->getNVarsForCons(seeed->getConssForBlock(b)[c]) && !assigned; ++v)
               {
                  if(var == seeedPropagationData->seeedpool->getVarsForCons(seeed->getConssForBlock(b)[c])[v])
                  {
                     seeed->setConsToBlock(cons, b);
                     seeed->deleteOpencons(cons);
                     seeed->setVarToBlock(var, b);
                     seeed->deleteOpenvar(var);
                     considerImplicitsWithStairlinkingvars(seeed, seeedPropagationData->seeedpool);
                     assignedConssForBfs.push_back(cons);
                     assigned = true;
                  }
               }
            }
         }
      }
   }


   for(size_t i = 0; i < assignedConssForBfs.size(); ++i)
   {
      it = find(conssForBfs.begin(), conssForBfs.end(), assignedConssForBfs[i]);
      conssForBfs.erase(it);
   }

   while(!conssForBfs.empty())
   {
      visitedConss.push_back(emptyVector);
      bfs(&visitedConss[newBlocks], &conssForBfs, seeedPropagationData->seeedpool); /** Breadth First Search */
      ++newBlocks;
   }

   if(newBlocks < 2)
   {
      seeedPropagationData->nNewSeeeds = 0;
      delete seeed;
   }
   else
   {
      for(size_t i = 0; i < visitedConss.size(); ++i)
      {
         block = seeed->addBlock();
         for(size_t c = 0; c < visitedConss[i].size(); ++c)
         {
            seeed->setConsToBlock(visitedConss[i][c], block);
            seeed->deleteOpencons(visitedConss[i][c]);
         }
      }

      seeed->considerImplicits(seeedPropagationData->seeedpool);
      seeedPropagationData->newSeeeds[0] = seeed;
      seeed->checkConsistency();
      seeedPropagationData->nNewSeeeds = 1;
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
/*
 * detector specific interface methods
 */

/** creates the handler for connected_noNewLinkingVars detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorConnected_noNewLinkingVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create connected_noNewLinkingVars detector data here*/
   detectordata = NULL;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, detectordata, detectConnected_noNewLinkingVars, freeConnected_noNewLinkingVars, initConnected_noNewLinkingVars, exitConnected_noNewLinkingVars, propagateSeeedConnected_noNewLinkingVars) );

   /**@todo add connected_noNewLinkingVars detector parameters */

   return SCIP_OKAY;
}
