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

/**@file   dec_connectedbase.c
 * @ingroup DETECTORS
 * @brief  detector connectedbase (put your description here)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_connectedbase.h"
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
#define DEC_DETECTORNAME         "connectedbase"       /**< name of detector */
#define DEC_DESC                 "detector connectedbase" /**< description of detector*/
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
DEC_DECL_FREEDETECTOR(freeConnectedbase)
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
DEC_DECL_EXITDETECTOR(exitConnectedbase)
{  /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitConnectedbase NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initConnectedbase)
{  /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initConnectedbase NULL
#endif

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectConnectedbase)
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

static
DEC_DECL_PROPAGATESEEED(propagateSeeedConnectedbase)
{
   *result = SCIP_DIDNOTFIND;

   int newBlocks = 0;
   std::vector<std::vector<int>> visitedConss = std::vector<std::vector<int>>(0); /** vector of vector with connected constraints */
   std::vector<int> emptyVector = std::vector<int>(0);
   std::vector<int> openConss;
   SCIP_HASHMAP* constoblock;
   int cons;

   //seeedPropagationData->seeedToPropagate->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));

   if(!seeedPropagationData->seeedToPropagate->areOpenVarsAndConssCalculated())
   {
      seeedPropagationData->seeedToPropagate->calcOpenconss();
      seeedPropagationData->seeedToPropagate->calcOpenvars();
      seeedPropagationData->seeedToPropagate->setOpenVarsAndConssCalculated(true);
   }
   assert(seeedPropagationData->seeedToPropagate->checkConsistency());

   for(int j = 0; j < seeedPropagationData->seeedToPropagate->getNOpenconss(); ++j)
   {
      openConss.push_back(seeedPropagationData->seeedToPropagate->getOpenconss()[j]);
   }

   while(!openConss.empty())
   {
      visitedConss.push_back(emptyVector);
      bfs(&visitedConss[newBlocks], &openConss, seeedPropagationData->seeedpool); /** Breadth First Search */
      ++newBlocks;
   }

   if(newBlocks < 2)
   {
      seeedPropagationData->nNewSeeeds = 0;
   }
   else
   {
      SCIP_CALL(SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), seeedPropagationData->seeedpool->getNConss()));

      for( int b = 0; b < seeedPropagationData->seeedToPropagate->getNBlocks(); ++b )
      {
         for( int c = 0; c < seeedPropagationData->seeedToPropagate->getNConssForBlock(b); ++c )
         {
            cons = seeedPropagationData->seeedToPropagate->getConssForBlock(b)[c];
            SCIP_CALL( SCIPhashmapSetImage( constoblock, (void*) (size_t)cons, (void*) (size_t) b+1) );
         }
      }

      for( int i = 0; i < newBlocks; ++i)
      {
         for( size_t c = 0; c < visitedConss[i].size(); ++c )
         {
            cons = visitedConss[i][c];
            SCIP_CALL( SCIPhashmapSetImage( constoblock, (void*) (size_t)cons, (void*) (size_t) seeedPropagationData->seeedToPropagate->getNBlocks() + i + 1));
         }
      }

      for( int j = 0; j < seeedPropagationData->seeedToPropagate->getNMasterconss(); ++j)
      {
         cons = seeedPropagationData->seeedToPropagate->getMasterconss()[j];
         SCIP_CALL( SCIPhashmapSetImage( constoblock, (void*) (size_t)cons, (void*) (size_t) seeedPropagationData->seeedToPropagate->getNBlocks() + newBlocks + 1));
      }

      gcg::Seeed* seeed;
      seeed = new gcg::Seeed(scip, seeedPropagationData->seeedpool->getNewIdForSeeed(), seeedPropagationData->seeedpool->getNDetectors(), seeedPropagationData->seeedpool->getNConss(), seeedPropagationData->seeedpool->getNVars());
      seeed->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
      seeed->filloutSeeedFromConstoblock(constoblock, seeedPropagationData->seeedToPropagate->getNBlocks() + newBlocks, seeedPropagationData->seeedpool);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
      seeedPropagationData->newSeeeds[0] = seeed;
      seeedPropagationData->nNewSeeeds = 1;
      seeedPropagationData->newSeeeds[0]->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
      SCIPhashmapFree(&constoblock);
   }

   *result = SCIP_SUCCESS;

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

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, detectordata, detectConnectedbase, freeConnectedbase, initConnectedbase, exitConnectedbase, propagateSeeedConnectedbase) );

   /**@todo add connectedbase detector parameters */

   return SCIP_OKAY;
}
