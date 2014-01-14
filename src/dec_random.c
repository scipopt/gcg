/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

/**@file   dec_random.c
 * @ingroup DETECTORS
 * @brief  Random structure detection for testing purposes
 * @author Martin Bergner
 */

/**
 * This detector will partition the constraints of the problem randomly.
 * For each constraint, it will randomly pick a number between 0 and the
 * maxblocks parameter. Constraints assigned to maxblocks will be put in the
 * master problem. You can set the random seed via the seed parameter.
 *
 * If that maxblocks parameter is set to -1, it will default to number of
 * constraints.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_random.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "scip/clock.h"
#include <string.h>

/* constraint handler properties */
#define DEC_DETECTORNAME         "random"    /**< name of detector */
#define DEC_DESC                 "Random structure detection" /**< description of detector*/
#define DEC_PRIORITY             -10         /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'r'         /**< display character of detector */
#define DEC_ENABLED              FALSE       /**< should the detection be enabled */
#define DEC_SKIP                 FALSE       /**< should detector be skipped if others found detections */

#define DEFAULT_MAXBLOCKS        -1          /**< the maximal number of blocks, -1 defaults to average number of constraints */
#define DEFAULT_AVGCONSPERBLOCK  100         /**< average constraints per block to limit the maximal block number */
#define DEFAULT_SEED             -1          /**< random seed for the random number generator, -1 is the current time */

/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
{
   int                   seed;               /**< random seed for the random number generator */
   int                   maxblocks;          /**< the maximal number of blocks, -1 defaults to nconss/maxconsperblock */
   int                   avgconsperblock;    /**< the average number of constraints per block */
   SCIP_HASHMAP*         constoblock;        /**< hashmap to store partition */
   int                   nblocks;            /**< number of actual blocks found */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** randomly assings constraints to blocks or the master */
static
SCIP_RETCODE findRandomPartition(
   SCIP*              scip,                  /**< SCIP data structure */
   DEC_DETECTORDATA*  detectordata           /**< detector data structure */
   )
{
   unsigned int seed;
   int nconss;
   SCIP_CONS** conss;
   int i;
   int maxblocks;
   int* consblocks;
   int oldblock;
   int nblocks;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(detectordata->constoblock != NULL);

   nconss = SCIPgetNConss(scip);
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &conss, SCIPgetConss(scip), nconss) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &consblocks, nconss) );
   BMSclearMemoryArray(consblocks, nconss);
   
   if( detectordata->seed == -1 )
      seed = SCIPround(scip, SCIPclockGetTimeOfDay());
   else
      seed = detectordata->seed;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " (seed = %d)", seed);

   if( detectordata->maxblocks == -1 )
      maxblocks = nconss/detectordata->avgconsperblock;
   else
      maxblocks = detectordata->maxblocks;

   nblocks = 0;

   for( i = 0; i < nconss; ++i )
   {
      consblocks[i] = SCIPgetRandomInt(0, maxblocks, &seed);
   }

   SCIPsortIntPtr(consblocks, (void**)conss, nconss);
   oldblock = -1;
   for( i = 0; i < nconss; ++i )
   {
      assert(consblocks[i] >= oldblock);
      if( consblocks[i] != oldblock )
      {
         oldblock = consblocks[i];
         ++nblocks;
      }
      SCIPdebugMessage("Assigning cons <%s> to block %d.\n", SCIPconsGetName(conss[i]), nblocks);
      assert(nblocks > 0);
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, conss[i], (void*)(size_t) (nblocks)) ); /*lint !e866*/
   }

   detectordata->nblocks = nblocks;

   SCIPfreeMemoryArray(scip, &consblocks);
   SCIPfreeMemoryArray(scip, &conss);
   return SCIP_OKAY;
}


/*
 * detector callback methods
 */

/** destructor of detector to free detector data (called before the solving process begins) */
static
DEC_DECL_EXITDETECTOR(exitRandom)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   if( detectordata->constoblock  != NULL )
      SCIPhashmapFree(&detectordata->constoblock);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detection initialization function of detector (called before solving is about to begin) */
static
DEC_DECL_INITDETECTOR(initRandom)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->nblocks = 0;

   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), SCIPgetNConss(scip)) );

   return SCIP_OKAY;
}

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectRandom)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting random structure:");

   SCIP_CALL( findRandomPartition(scip, detectordata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d blocks.\n", detectordata->nblocks);

   if( detectordata->nblocks > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/
      SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[0])) );

      SCIP_CALL( DECfilloutDecdecompFromConstoblock(scip, (*decdecomps)[0], detectordata->constoblock, detectordata->nblocks, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );
      *ndecdecomps = 1;

      detectordata->constoblock = NULL;
      *result = SCIP_SUCCESS;
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
      if( detectordata->constoblock != NULL )
      {
         SCIPhashmapFree(&detectordata->constoblock);
      }
   }

   return SCIP_OKAY;
}


/*
 * detector specific interface methods
 */

/** creates the handler for random detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionRandom(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   detectordata->maxblocks = 0;
   detectordata->seed = -1;
   detectordata->constoblock = NULL;
   detectordata->nblocks = 0;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, detectordata, detectRandom, initRandom, exitRandom) );

   SCIP_CALL( SCIPaddIntParam(scip, "detectors/random/seed", "random seed for the random number generator, -1 is the current time", &detectordata->seed, FALSE, DEFAULT_SEED, -1, INT_MAX, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/random/maxblocks", "the maximal number of blocks, -1 defaults to avgconsperblock", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, -1, INT_MAX, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/random/avgconsperblock", "average constraints per block", &detectordata->avgconsperblock, FALSE, DEFAULT_AVGCONSPERBLOCK, 1, 10000, NULL, NULL ) );

   return SCIP_OKAY;
}
