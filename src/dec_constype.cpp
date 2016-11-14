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

/**@file   dec_constype.cpp
 * @ingroup DETECTORS
 * @brief  detector constype (put your description here)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_constype.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip_misc.h"

#include <iostream>

/* constraint handler properties */
#define DEC_DETECTORNAME         "constype"       /**< name of detector */
#define DEC_DESC                 "detector constype" /**< description of detector*/
#define DEC_FREQCALLROUND        1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND         0           /** last round the detector gets called                              */
#define DEC_MINCALLROUND         0           /** first round the detector gets called                              */
#define DEC_PRIORITY             0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              '?'         /**< display character of detector */
#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#define DEC_SKIP                 FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL         FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
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

/** method to enumerate all subsets */
std::vector< std::vector<int> > getAllSubsets(std::vector<int> set)
{
    std::vector< std::vector<int> > subset;
    std::vector<int> empty;
    subset.push_back( empty );

    for (int i = 0; i < set.size(); i++)
    {
        std::vector< std::vector<int> > subsetTemp = subset;

        for (int j = 0; j < subsetTemp.size(); j++)
            subsetTemp[j].push_back( set[i] );

        for (int j = 0; j < subsetTemp.size(); j++)
            subset.push_back( subsetTemp[j] );
    }
    return subset;
}

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeConstype)
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
DEC_DECL_EXITDETECTOR(exitConstype)
{ /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitConstype NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initConstype)
{ /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initConstype NULL
#endif

/** detection function of detector */
static DEC_DECL_DETECTSTRUCTURE(detectConstype)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME)
;   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

static DEC_DECL_PROPAGATESEEED(propagateSeeedConstype)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CONS* cons;

   int seeedCounter = 0;
  gcg::Seeed* seeedOrig;
  gcg::Seeed* seeed;

  std::vector<consType> foundConstypes(0);
  std::vector<int> constypesIndices(0);

  seeedOrig = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);
  seeedOrig->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));


  for( size_t i = 0; i < seeedOrig->getNOpenconss(); ++i)
  {
      cons = seeedPropagationData->seeedpool->getConsForIndex(seeedOrig->getOpenconss()[i]);
      consType cT = GCGconsGetType(cons);

      /** find constype or not */
      std::vector<consType>::const_iterator constypeIter = foundConstypes.begin();
      for(; constypeIter != foundConstypes.end(); ++constypeIter)
      {
    	  if(*constypeIter == cT)
    		  break;
      }

      if( constypeIter  == foundConstypes.end()  )
      {
         foundConstypes.push_back(GCGconsGetType(cons) );
      }
  }

  for(int i = 0; i < foundConstypes.size(); ++i)
  {
      constypesIndices.push_back(i);
  }

  std::vector< std::vector<int> > subsetsOfConstypes = getAllSubsets(constypesIndices);

  SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), subsetsOfConstypes.size() - 1) );
  seeedPropagationData->nNewSeeeds = subsetsOfConstypes.size() - 1;


  if(!seeedOrig->areOpenVarsAndConssCalculated())
  {
      seeedOrig->calcOpenconss();
      seeedOrig->calcOpenvars();
      seeedOrig->setOpenVarsAndConssCalculated(true);
  }

  for(size_t subset = 0; subset < subsetsOfConstypes.size(); ++subset)
  {
      if(subsetsOfConstypes[subset].size() == 0)
          continue;

      seeed = new gcg::Seeed(seeedOrig, seeedPropagationData->seeedpool);
         /** set open cons that have type of the current subset to Master */
      for( size_t i = 0; i < seeed->getNOpenconss(); ++i)
      {
          for(size_t constypeId = 0; constypeId < subsetsOfConstypes[subset].size(); ++constypeId )
          {
              cons = seeedPropagationData->seeedpool->getConsForIndex(seeed->getOpenconss()[i]);
              if( GCGconsGetType   (cons) == foundConstypes[subsetsOfConstypes[subset][constypeId]] )
              {
                  seeed->bookAsMasterCons(seeed->getOpenconss()[i]);
              }
          }
      }
      seeed->flushBooked();
      seeedPropagationData->newSeeeds[seeedCounter] = seeed;
      seeedCounter++;

  }

  delete seeedOrig;
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
/*
 * detector specific interface methods
 */

/** creates the handler for constype detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorConstype(SCIP* scip /**< SCIP data structure */
)
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create constype detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, DEC_USEFULRECALL, detectordata, detectConstype, freeConstype, initConstype, exitConstype, propagateSeeedConstype));

   /**@todo add constype detector parameters */

   return SCIP_OKAY;
}
