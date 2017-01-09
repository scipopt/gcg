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

/**@file   dec_hrgpartition.cpp
 * @brief  arrowhead and bordered detector via graph partitioning (uses hmetis)
 * @ingroup DETECTORS
 * @author Martin Bergner
 *
 * Detects arrowhead (double bordered) decompositions as well as decompositions
 * with only linking variables or linking constraints.
 *
 * This detector needs hmetis and works only under Linux/MacOS, it further needs the Z-shell (zsh)
 * to enforce memory and time limits on hmetis as this is the only shell reliably doing that.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */
#include "dec_hrgpartition.h"

#if !defined(_WIN32) && !defined(_WIN64)
#include <cassert>
#include <cstring>
#include <cerrno>
#include <unistd.h>
#include <iostream>


#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "graph/matrixgraph.h"
#include "graph/hyperrowgraph.h"
#include "graph/graph_tclique.h"
#include "graph/weights.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "scip/clock.h"


#include <set>

using gcg::HyperrowGraph;
using gcg::MatrixGraph;
using gcg::Weights;

#define DEC_DETECTORNAME      "hrgpartition"    /**< name of the detector */
#define DEC_DESC              "enforces arrowhead structures using graph partitioning" /**< description of detector */
#define DEC_FREQCALLROUND        1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND         INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND         0           /** first round the detector gets called                              */
#define DEC_PRIORITY          1000           /**< priority of the detector */
#define DEC_DECCHAR           'r'            /**< display character of detector */
#define DEC_ENABLED           TRUE           /**< should detector be called by default */
#define DEC_ENABLEDFINISHING  FALSE          /**< should the finishing be enabled */
#define DEC_SKIP              FALSE          /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL      TRUE           /**< is it useful to call this detector on a descendant of the propagated seeed */


/* Default parameter settings */
#define DEFAULT_VARWEIGHT         1          /**< weight for variable nodes */
#define DEFAULT_VARWEIGHTBIN      2          /**< weight for binary variable nodes */
#define DEFAULT_VARWEIGHTINT      2          /**< weight for integer variable nodes */
#define DEFAULT_VARWEIGHTIMPL     2          /**< weight for implicit integer variable nodes */
#define DEFAULT_VARWEIGHTCONT     1          /**< weight for continous variable nodes */
#define DEFAULT_CONSWEIGHT        5          /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED          1          /**< random seed for the hmetis call */
#define DEFAULT_TIDY              TRUE       /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES        0.2        /**< percentage of dummy vertices*/
#define DEFAULT_CONSWEIGHT_SETPPC 5          /**< weight for constraint hyperedges that are setpartitioning or covering
                                                  constraints */
#define DEFAULT_MINBLOCKS         2          /**< value for the minimum number of blocks to be considered */
#define DEFAULT_MAXBLOCKS         20         /**< value for the maximum number of blocks to be considered */
#define DEFAULT_ALPHA             0.0        /**< factor for standard deviation of constraint weights */
#define DEFAULT_BETA              0.5        /**< factor of how the weight for equality and inequality constraints is
                                                  distributed (keep 1/2 for the same on both) */
#define DEFAULT_METIS_UBFACTOR    5.0        /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE     FALSE      /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB  TRUE       /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_REALNAME          FALSE      /**< whether the metis name should be real or temporary */
#define DEFAULT_TYPE              'r'        /**< type of the decomposition 'c' column hypergraph (single bordered, no
                                                  linking constraints), 'r' row hypergraph (single bordered, no linking
                                                  variables) and 'a' column-row hypergraph (arrowhead) */
/*
 * Data structures
 */

/** private detector data */
struct DEC_DetectorData
{
   /* Graph stuff for hmetis */
   MatrixGraph<gcg::GraphTclique>* graph;    /**< the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN];            /**< filename for the metis input file */

   /* weight parameters */
   int       varWeight;             /**< weight of a variable hyperedge */
   int       varWeightBinary;       /**< weight of a binary variable hyperedge */
   int       varWeightContinous;    /**< weight of a continuous variable hyperedge */
   int       varWeightInteger;      /**< weight of an integer variable hyperedge */
   int       varWeightImplint;      /**< weight of an implicit integer variable hyperedge */
   int       consWeight;            /**< weight of a constraint hyperedge */
   int       consWeightSetppc;      /**< weight of a setppc constraint hyperedge */
   SCIP_Real alpha;                 /**< factor for constraint coefficient value standard deviation */
   SCIP_Real beta;                  /**< factor for equality od inequality constraints */

   /* general parameters */
   SCIP_Real dummynodes;      /**< percent of dummy nodes */
   SCIP_Bool tidy;            /**< whether tempory metis files should be cleaned up */
   int       maxblocks;       /**< maximal number of blocks to test */
   int       minblocks;       /**< minimal number of blocks to test */

   /* metis parameters */
   int       randomseed;      /**< metis random seed */
   SCIP_Real metisubfactor;   /**< metis unbalance factor */
   SCIP_Bool metisverbose;    /**< should metis ouput be displayed */
   SCIP_Bool metisuseptyperb; /**< flag to indicate whether metis uses kway or rb partitioning */
   SCIP_Bool realname;        /**< flag to indicate real problem name or temporary filename for metis files */

   /* various data */
   SCIP_CLOCK* metisclock;    /**< clock to measure metis time */
   int         blocks;        /**< indicates the current block */
   SCIP_Bool   found;         /**< indicates whethere a decomposition has been found */
   char        type;          /**< type of the decomposition 'c' column hypergraph (single bordered, no linking
                                   constraints), 'r' row hypergraph (single bordered, no linking variables) and
                                   'a' column-row hypergraph (arrowhead) */
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */



/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeHrgpartition)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}



/** detector initialization method */
static
DEC_DECL_INITDETECTOR(initHrgpartition)
{
   int nconss;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);


   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);

   SCIP_CALL( SCIPcreateWallClock(scip, &detectordata->metisclock) );

   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DEC_DECL_EXITDETECTOR(exitHrgpartition)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);


   SCIP_CALL( SCIPfreeClock(scip, &detectordata->metisclock) );

   return SCIP_OKAY;
}

/** will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   SCIP_RESULT*          result              /**< result indicating whether the detection was successful */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];

   int status;

   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   remainingtime = DECgetRemainingTime(scip);

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   /* call metis via syscall as there is no library usable ... */
   if( !SCIPisInfinity(scip, DECgetRemainingTime(scip)) )
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"ulimit -t %.0f;hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               remainingtime,
               detectordata->tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }
   else
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               detectordata->tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }

   SCIP_CALL( SCIPresetClock(scip, detectordata->metisclock) );
   SCIP_CALL( SCIPstartClock(scip, detectordata->metisclock) );
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " %d", detectordata->blocks );
   status = system( metiscall );

   SCIP_CALL( SCIPstopClock(scip, detectordata->metisclock) );
   SCIPdebugMessage("time left before metis started: %f, time metis spend %f, remainingtime: %f\n", remainingtime, SCIPgetClockTime(scip, detectordata->metisclock),  remainingtime-SCIPgetClockTime(scip, detectordata->metisclock) );

   /* check error codes */
   if( status == -1 )
   {
      SCIPerrorMessage("System call did not succed: %s\n", strerror( errno ));
      SCIPerrorMessage("Call was %s\n", metiscall);
   }
   else if( status != 0 )
   {

      SCIPerrorMessage("Calling hmetis unsuccessful! See the above error message for more details.\n");
      SCIPerrorMessage("Call was %s\n", metiscall);
   }

   /* exit gracefully in case of errors */
   if( status != 0 )
   {
      return SCIP_ERROR;
   }

   (void) SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d", detectordata->tempfile, detectordata->blocks);
   SCIP_CALL( detectordata->graph->readPartition(metisout) );

   /* if desired delete the temoprary metis file */
   if( detectordata->tidy )
   {
      status = remove( metisout );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: %s\n", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Temporary file is in: %s\n", detectordata->tempfile);
   }
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** creates the temporary metis input file */
static
SCIP_RETCODE createMetisFile(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   int nvertices;
   int ndummyvertices;
   int fd;
   nvertices = detectordata->graph->getNNonzeroes();
   /*lint --e{524}*/
   ndummyvertices = SCIPceil(scip, detectordata->dummynodes*nvertices);
   detectordata->graph->setDummynodes(ndummyvertices);

   if( !detectordata->realname )
   {
      (void) SCIPsnprintf(detectordata->tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   }
   else
   {
      (void) SCIPsnprintf(detectordata->tempfile, SCIP_MAXSTRLEN, "gcg-%s-XXXXXX", SCIPgetProbName(scip));
   }

   fd = mkstemp(detectordata->tempfile);

   SCIP_CALL( detectordata->graph->writeToFile(fd, TRUE) );
   close(fd);
   return SCIP_OKAY;
}

/** are there conss and vars to be included by the graph */
static
bool graphCompletible(
   gcg::Seeedpool*  seeedpool,
   gcg::Seeed*      seeed
   )
{
   for(int c = 0; c < seeed->getNOpenconss(); ++c)
   {
      int cons = seeed->getOpenconss()[c];
      for(int v = 0; v < seeed->getNOpenvars(); ++v)
      {
         int var = seeed->getOpenvars()[v];
         for(int i = 0; i < seeedpool->getNVarsForCons(cons); ++i)
         {
            if(var == seeedpool->getVarsForCons(cons)[i])
            {
               return true;
            }
         }
      }
   }
   return false;
}

/** detection callback method */
static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildArrowhead)
{
   int i;
   int j;
   int ndecs;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   ndecs = detectordata->maxblocks-detectordata->minblocks+1;
   *ndecdecomps = 0;

   /* allocate space for output data */
   assert(detectordata->maxblocks >= detectordata->minblocks);
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, ndecs) );

   /* build the hypergraph structure from the original problem */

   Weights w(detectordata->varWeight, detectordata->varWeightBinary, detectordata->varWeightContinous,detectordata->varWeightInteger,detectordata->varWeightInteger,detectordata->consWeight);
   detectordata->graph = new HyperrowGraph<gcg::GraphTclique>(scip, w);

   SCIP_CALL( detectordata->graph->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( createMetisFile(scip, detectordata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");
   for( j = 0, i = detectordata->minblocks; i <= detectordata->maxblocks; ++i )
   {
      SCIP_RETCODE retcode;
      detectordata->blocks = i;
      /* get the partitions for the new variables from metis */
      retcode = callMetis(scip, detectordata, result);

      if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY )
      {
         continue;
      }

      SCIP_CALL( detectordata->graph->createDecompFromPartition(&(*decdecomps)[j]) );
      if( (*decdecomps)[j] != NULL )
      {
         *ndecdecomps += 1;
         ++j;
         detectordata->found = TRUE;
      }
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d decompositions found.\n",  *ndecdecomps);

   delete detectordata->graph;
   detectordata->graph = NULL;

   SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );

   if( detectordata->tidy )
   {
      int status = remove( detectordata->tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }


   *result = detectordata->found ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}
#endif

static
DEC_DECL_PROPAGATESEEED(propagateSeeedHrgpartition)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CLOCK* clock;
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   DEC_DETECTORDATA* detectordata = DECdetectorGetData(detector);
   int nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);


   SCIP_CALL( SCIPresetClock(scip, detectordata->metisclock) );

   /* add hrgpartition presolver parameters */

   std::vector<SCIP_Real> clockTimes;        /**< vector containing times in seconds  */
   int k;
   int j;
   int s;
   int nMaxSeeeds;
   int nNewSeeeds = 0;
   gcg::Seeed* seeed;
   gcg::Seeed** newSeeeds;
   std::vector<int> numberOfBlocks = seeedPropagationData->seeedpool->getCandidatesNBlocks();

   assert(scip != NULL);
   assert(detectordata != NULL);

   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);
   seeed->considerImplicits(seeedPropagationData->seeedpool);
   seeed->refineToMaster(seeedPropagationData->seeedpool);

   seeedPropagationData->seeedpool->decrementSeeedcount();

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   nMaxSeeeds = detectordata->maxblocks-detectordata->minblocks+1;

   /* allocate space for output data */
   assert(detectordata->maxblocks >= detectordata->minblocks);
   SCIP_CALL( SCIPallocBufferArray(scip, &(newSeeeds), 2 * nMaxSeeeds) );

   if(!graphCompletible(seeedPropagationData->seeedpool, seeed) || seeed->alreadyAssignedConssToBlocks() )
   {
      delete seeed;
      seeedPropagationData->nNewSeeeds = 0;
      SCIPfreeBufferArray(scip, &newSeeeds);
      SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
      SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* build the hypergraph structure from the original problem */

   Weights w(detectordata->varWeight, detectordata->varWeightBinary, detectordata->varWeightContinous,detectordata->varWeightInteger,detectordata->varWeightInteger,detectordata->consWeight);
   detectordata->graph = new HyperrowGraph<gcg::GraphTclique>(scip, w);

   SCIP_CALL( detectordata->graph->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed) );

   SCIP_CALL( createMetisFile(scip, detectordata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");
   SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &temporaryClock) );
   for( j = 0, k = 0; k < (int) numberOfBlocks.size(); ++k)
   {
      int nblocks = numberOfBlocks[k] - seeed->getNBlocks();
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      SCIP_RETCODE retcode;
      detectordata->blocks = nblocks;

      if(nblocks > seeedPropagationData->seeedToPropagate->getNOpenconss() || nblocks <= 0)
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
         SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
         continue;
      }

      retcode = callMetis(scip, detectordata, result);

      if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY)
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
         SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
         continue;
      }

      SCIP_CALL( detectordata->graph->createSeeedFromPartition(seeed, &newSeeeds[j], &newSeeeds[j+1], seeedPropagationData->seeedpool) );

      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      if( (newSeeeds)[j] != NULL )
      {
         nNewSeeeds = nNewSeeeds + 2;
         detectordata->found = TRUE;
         clockTimes.push_back(SCIPclockGetTime(temporaryClock));
         clockTimes.push_back(SCIPclockGetTime(temporaryClock)); // 2x because two seeeds where created
      }
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
      j = j + 2;
   }
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock ) );


   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d seeeds found.\n",  nNewSeeeds);

   delete detectordata->graph;
   detectordata->graph = NULL;
   delete seeed;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), nNewSeeeds) );
   seeedPropagationData->nNewSeeeds = nNewSeeeds;
   for(j = 0, s = 0; s < nNewSeeeds; ++j)
   {
      if(newSeeeds[j] != NULL)
      {
         seeedPropagationData->newSeeeds[s] = newSeeeds[j];
         seeedPropagationData->newSeeeds[s]->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
         ++s;
      }
   }
   SCIPfreeBufferArray(scip, &newSeeeds);

   if( detectordata->tidy )
   {
      int status = remove( detectordata->tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
         SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
         SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );
         return SCIP_WRITEERROR;
      }
   }

   SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );

   for( s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
      seeedPropagationData->newSeeeds[s]->addClockTime( SCIPclockGetTime(clock) + clockTimes[s] );
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );

   *result = detectordata->found ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   return SCIP_OKAY;

}

static
DEC_DECL_FINISHSEEED(finishSeeedHrgpartition)
{
   *result = SCIP_DIDNOTFIND;

   DEC_DETECTORDATA* detectordata = DECdetectorGetData(detector);
   int nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);


   SCIP_CALL( SCIPresetClock(scip, detectordata->metisclock) );

   /* add hrgpartition presolver parameters */

   std::vector<SCIP_Real> clockTimes;        /**< vector containing times in seconds  */
   int k;
   int j;
   int s;
   int nMaxSeeeds;
   int nNewSeeeds = 0;
   gcg::Seeed* seeed;
   gcg::Seeed** newSeeeds;
   std::vector<int> numberOfBlocks = seeedPropagationData->seeedpool->getCandidatesNBlocks();

   assert(scip != NULL);
   assert(detectordata != NULL);

   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);
   seeed->considerImplicits(seeedPropagationData->seeedpool);
   seeed->refineToMaster(seeedPropagationData->seeedpool);

   seeedPropagationData->seeedpool->decrementSeeedcount();

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   nMaxSeeeds = detectordata->maxblocks-detectordata->minblocks+1;

   /* allocate space for output data */
   assert(detectordata->maxblocks >= detectordata->minblocks);
   SCIP_CALL( SCIPallocBufferArray(scip, &(newSeeeds), 2 * nMaxSeeeds) );

   if(!graphCompletible(seeedPropagationData->seeedpool, seeed))
   {
      delete seeed;
      seeedPropagationData->nNewSeeeds = 0;
      SCIPfreeBufferArray(scip, &newSeeeds);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* build the hypergraph structure from the original problem */

   Weights w(detectordata->varWeight, detectordata->varWeightBinary, detectordata->varWeightContinous,detectordata->varWeightInteger,detectordata->varWeightInteger,detectordata->consWeight);
   detectordata->graph = new HyperrowGraph<gcg::GraphTclique>(scip, w);

   SCIP_CALL( detectordata->graph->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed) );

   SCIP_CALL( createMetisFile(scip, detectordata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");
   for( j = 0, k = 0; k < (int) numberOfBlocks.size(); ++k)
   {
      int nblocks = numberOfBlocks[k] - seeed->getNBlocks();
      SCIP_RETCODE retcode;
      detectordata->blocks = nblocks;

      if(nblocks > seeedPropagationData->seeedToPropagate->getNOpenconss() || nblocks <= 0)
      {
         continue;
      }

      retcode = callMetis(scip, detectordata, result);

      if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY)
      {
         continue;
      }

      SCIP_CALL( detectordata->graph->createSeeedFromPartition(seeed, &newSeeeds[j], &newSeeeds[j+1], seeedPropagationData->seeedpool) );

      if( (newSeeeds)[j] != NULL )
      {
         nNewSeeeds = nNewSeeeds + 2;
         detectordata->found = TRUE;
      }
      j = j + 2;
   }


   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d seeeds found.\n",  nNewSeeeds);

   delete detectordata->graph;
   detectordata->graph = NULL;
   delete seeed;
   assert(nNewSeeeds % 2 == 0);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), nNewSeeeds/2) );
   seeedPropagationData->nNewSeeeds = nNewSeeeds/2;
   for(j = 0, s = 0; s < nNewSeeeds/2; j+=2)
   {
      if((newSeeeds)[j] != NULL)
      {
         newSeeeds[j]->considerImplicits(seeedPropagationData->seeedpool);
         newSeeeds[j]->assignAllDependent(seeedPropagationData->seeedpool);
         seeedPropagationData->newSeeeds[s] = newSeeeds[j];
         assert(newSeeeds[j]->getNOpenconss() == 0);
         assert(newSeeeds[j]->getNOpenvars() == 0);
         ++s;
      }
   }
   SCIPfreeBufferArray(scip, &newSeeeds);

   if( detectordata->tidy )
   {
      int status = remove( detectordata->tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }


   *result = detectordata->found ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   return SCIP_OKAY;

}


/** creates the hrgpartition presolver and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorHrgpartition(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#if !defined(_WIN32) && !defined(_WIN64)
   DEC_DETECTORDATA *detectordata = NULL;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->blocks = -1;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_SKIP, DEC_USEFULRECALL, detectordata, detectAndBuildArrowhead, freeHrgpartition, initHrgpartition, exitHrgpartition, propagateSeeedHrgpartition, finishSeeedHrgpartition) );


   /* add hrgpartition presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/hrgpartition/beta", "factor on how heavy equality (beta) and inequality constraints are measured", &detectordata->beta, FALSE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/hrgpartition/alpha", "factor on how heavy the standard deviation of the coefficients is measured", &detectordata->alpha, FALSE, DEFAULT_ALPHA, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/varWeight", "Weight of a variable hyperedge", &detectordata->varWeight, FALSE, DEFAULT_VARWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/varWeightBinary", "Weight of a binary variable hyperedge", &detectordata->varWeightBinary, FALSE, DEFAULT_VARWEIGHTBIN, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/varWeightContinous", "Weight of a continuos variable hyperedge", &detectordata->varWeightContinous, FALSE, DEFAULT_VARWEIGHTCONT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/varWeightImplint", "Weight of a implicit integer variable hyperedge", &detectordata->varWeightImplint, FALSE, DEFAULT_VARWEIGHTIMPL, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/varWeightInteger", "Weight of a integer variable hyperedge", &detectordata->varWeightInteger, FALSE, DEFAULT_VARWEIGHTINT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/hrgpartition/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/hrgpartition/dummynodes", "percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/hrgpartition/consWeightSetppc", "Weight for constraint hyperedges that are setpartitioning or covering constraints", &detectordata->consWeightSetppc, FALSE, DEFAULT_CONSWEIGHT_SETPPC, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/hrgpartition/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/hrgpartition/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/hrgpartition/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/hrgpartition/realname", "Should the problem be used for metis files or a temporary name", &detectordata->realname, FALSE, DEFAULT_REALNAME, NULL, NULL) );
#endif
   return SCIP_OKAY;
}
