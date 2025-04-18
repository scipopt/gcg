/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_hrgpartition.cpp
 * 
 * @brief  arrowhead and bordered detector via graph partitioning (uses hmetis)
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

#include "gcg/dec_hrgpartition.h"

#if !defined(_WIN32) && !defined(_WIN64)
#include <cassert>
#include <cstring>
#include <cerrno>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>

#ifdef HMETIS_HEADER
#include "gcg/hmetis.h"
#else
#define HMETIS_EXECUTABLE "hmetis"
#endif

#include "gcg/cons_decomp.h"
#include "gcg/struct_decomp.h"
#include "gcg/pub_decomp.h"
#include "gcg/scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "graph/matrixgraph.h"
#include "graph/hyperrowgraph.h"
#include "graph/graph_tclique.h"
#include "graph/weights.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "scip/clock.h"


#include <set>

using gcg::HyperrowGraph;
using gcg::MatrixGraph;
using gcg::Weights;

#define DEC_NAME                  "hrgpartition"    /**< name of the detector */
#define DEC_DESC                  "enforces arrowhead structures using graph partitioning" /**< description of detector */
#define DEC_FREQCALLROUND         1           /**< frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0           /**< last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /**< first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /**< frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  0           /**< last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /**< first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              1000        /**< priority of the detector */
#define DEC_DECCHAR               'r'         /**< display character of detector */
#define DEC_ENABLED               FALSE       /**< should detector be called by default */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL          TRUE        /**< is it useful to call this detector on a descendant of the propagated partialdec */


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

#define DEFAULT_LIMITNCONSSNVARSDEFAULT  10000  /**< limit for sum of nvars and nconss for enabling this detector in default */
#define DEFAULT_ENABLEDFORLARGEPROBLEMS  FALSE  /** should this detector be enabled even the limit nconssnvars is exceeded */

#define DEFAULT_MAXNBLOCKCANDIDATES 3          /**< number of block number candidates to be considered */
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
#define DEFAULT_FALLBACK_NBLOCKS    8
#define SET_MULTIPLEFORSIZETRANSF 12500

/*
 * Data structures
 */

/** private detector data */
struct GCG_DetectorData
{

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
   int       limitnconssnvarsdefault;       /**<  limit for sum of nvars and nconss for enabling this detector in default*/
   SCIP_Bool enabledforlargeproblems;

   int       maxnblockcandidates;       /**< maximal number of block canddidates to test */

   /* metis parameters */
   int       randomseed;      /**< metis random seed */
   SCIP_Real metisubfactor;   /**< metis unbalance factor */
   SCIP_Bool metisverbose;    /**< should metis ouput be displayed */
   SCIP_Bool metisuseptyperb; /**< flag to indicate whether metis uses kway or rb partitioning */
   SCIP_Bool realname;        /**< flag to indicate real problem name or temporary filename for metis files */

   /* various data */
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
GCG_DECL_FREEDETECTOR(freeHrgpartition)
{
   GCG_DETECTORDATA* detectordata;

   assert(gcg != NULL);

   detectordata = GCGdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(GCGdetectorGetName(detector), DEC_NAME) == 0);

   SCIPfreeMemory(GCGgetOrigprob(gcg), &detectordata);

   return SCIP_OKAY;
}



/** detector initialization method */
static
GCG_DECL_INITDETECTOR(initHrgpartition)
{
   assert(gcg != NULL);

   assert(strcmp(GCGdetectorGetName(detector), DEC_NAME) == 0);

   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
GCG_DECL_EXITDETECTOR(exitHrgpartition)
{

   assert(gcg != NULL);

   assert(strcmp(GCGdetectorGetName(detector), DEC_NAME) == 0);


   return SCIP_OKAY;
}

/** will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
   GCG*                  gcg,                /**< GCG data struture */
   GCG_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   MatrixGraph<gcg::GraphTclique>* graph,    /**< the graph of the matrix */
   char                  tempfile[SCIP_MAXSTRLEN],  /**< filename for the metis input file */
   int                   nblocks,            /**< number of blocks */
   SCIP_RESULT*          result              /**< result indicating whether the detection was successful */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   SCIP_CLOCK* metisclock;
   SCIP* scip;

   int status;

   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   scip = GCGgetOrigprob(gcg);
   remainingtime = GCGgetRemainingTime(gcg);
   SCIP_CALL( SCIPcreateWallClock(scip, &metisclock) );

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   /* call metis via syscall as there is no library usable ... */
   if( !SCIPisInfinity(scip, GCGgetRemainingTime(gcg)) )
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"ulimit -t %.0f;" HMETIS_EXECUTABLE " %s %d -seed %d -ptype %s -ufactor %f %s\"",
               remainingtime,
               tempfile,
               nblocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }
   else
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"" HMETIS_EXECUTABLE " %s %d -seed %d -ptype %s -ufactor %f %s\"",
               tempfile,
               nblocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }

   SCIP_CALL( SCIPstartClock(scip, metisclock) );
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " %d", nblocks );
   status = system(metiscall);

   SCIP_CALL( SCIPstopClock(scip, metisclock) );
   SCIPdebugMessage("time left before metis started: %f, time metis spend %f, remainingtime: %f\n", remainingtime, SCIPgetClockTime(scip, metisclock),  remainingtime-SCIPgetClockTime(scip, metisclock) );

   SCIP_CALL( SCIPfreeClock(scip, &metisclock) );

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

   (void) SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d", tempfile, nblocks);
   SCIP_CALL( graph->readPartition(metisout) );

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
      SCIPinfoMessage(scip, NULL, "Temporary file is in: %s\n", tempfile);
   }
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** creates the temporary metis input file */
static
SCIP_RETCODE createMetisFile(
   GCG*                  gcg,                /**< GCG data struture */
   GCG_DETECTORDATA*     detectordata,       /**< detector data structure */
   int                   partialdecID,       /**< ID of partial decomposition to be used */
   MatrixGraph<gcg::GraphTclique>* graph,    /**< the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN]             /**< filename for the metis input file */

   )
{
   int nvertices;
   int ndummyvertices;
   int fd;
   SCIP* scip = GCGgetOrigprob(gcg);
   nvertices = graph->getNNonzeroes();
   /*lint --e{524}*/
   ndummyvertices = SCIPceil(scip, detectordata->dummynodes*nvertices);
   graph->setDummynodes(ndummyvertices);

   if( !detectordata->realname )
      {
         (void) SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-%c-%d.metis.XXXXXX", DEC_DECCHAR, partialdecID );
      }
      else
      {
         (void) SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-%s-%s-%c.metis.XXXXXX", SCIPgetProbName(scip), DEC_DECCHAR, partialdecID);
      }

   fd = mkstemp(tempfile);

   SCIP_CALL( graph->writeToFile(fd, TRUE) );
   close(fd);
   return SCIP_OKAY;
}

/** @returns whether the hyperrolgraph is connected */
static
bool connected(
   gcg::DETPROBDATA*  detprobdata,
   gcg::PARTIALDECOMP*      partialdec
   )
{
   std::vector<int> queue;
   std::vector<int> visited;
   std::vector<bool> inqueue(detprobdata->getNVars(), false);
   std::vector<bool> isvisited(detprobdata->getNVars(), false);
   int start = -1;

   if(partialdec->getNOpenvars() < 2)
      return false;

   start = partialdec->getOpenvars()[0];
   queue.push_back(start);
   inqueue[start] = true;
   do
   {
      int node = queue[0];
      queue.erase(queue.begin());
      inqueue[node] = false;
      visited.push_back(node);
      isvisited[node] = true;
      for( int c = 0; c < detprobdata->getNConssForVar(node); ++c )
      {
         int cons = detprobdata->getConssForVar(node)[c];
         if(!partialdec->isConsOpencons(cons))
            continue;
         for( int v = 0; v < detprobdata->getNVarsForCons(cons); ++v )
         {
            int var = detprobdata->getVarsForCons(cons)[v];
            if( !partialdec->isVarOpenvar(var) )
               continue;
            if( isvisited[var] )
               continue;
            if( inqueue[var] )
               continue;
            queue.push_back(var);
            inqueue[var] = true;
         }
      }
   } while( !queue.empty() );

   if((int)visited.size() != partialdec->getNOpenvars())
      return false;
   else
      return true;
}

/** detection function for partialdecs */
static
SCIP_RETCODE detection(
   GCG*                    gcg,                          /**< GCG data structure */
   GCG_DETECTORDATA*       detectordata,                 /**< detectordata of the detector */
   Partialdec_Detection_Data*   partialdecdetectiondata, /**< partialdecdetectiondata (including the detprobdata) where to store the new Partialdecs */
   gcg::PARTIALDECOMP*     partialdec,                   /**< partialdec to propagate */
   bool                    allowopenpartialdecs,         /**< whether new partialdecs should be stored in which this detector only assignes conss to master */
   SCIP_RESULT*            result                        /**< pointer where to store the result */
)
{
	/* add hrgpartition presolver parameters */
	char decinfo[SCIP_MAXSTRLEN];
	int maxnblockcandidates;
	int k;
   int j;
   int s;
   gcg::PARTIALDECOMP** newpartialdecs;
   SCIP_CLOCK* clock;
   SCIP_CLOCK* temporaryClock;
   std::vector<SCIP_Real> clockTimes;        /* vector containing times in seconds  */

   /* Graph stuff for hmetis */
   MatrixGraph<gcg::GraphTclique>* graph;    /* the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN];            /**< filename for the metis input file */
   SCIP* scip = GCGgetOrigprob(gcg);

   SCIP_CALL_ABORT( SCIPcreateClock(scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   *result = SCIP_DIDNOTFIND;

   std::vector<int> numberOfBlocks;
   partialdecdetectiondata->detprobdata->getSortedCandidatesNBlocks(numberOfBlocks);
   if( numberOfBlocks.empty() )
      numberOfBlocks.push_back(DEFAULT_FALLBACK_NBLOCKS);

   SCIP_CALL( SCIPgetIntParam(scip, "detection/detectors/hrgpartition/maxnblockcandidates", &maxnblockcandidates) );

   maxnblockcandidates = MIN(maxnblockcandidates, (int) numberOfBlocks.size() );

   SCIPdebugMessage("number of block numbers to test: %d\n", maxnblockcandidates);

   assert(scip != NULL);
   assert(detectordata != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_NAME);

   /* allocate space for output data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newpartialdecs), 2 * maxnblockcandidates));

    /* build the hypergraph structure from the original problem */

   Weights w(detectordata->varWeight, detectordata->varWeightBinary, detectordata->varWeightContinous,detectordata->varWeightInteger,detectordata->varWeightInteger,detectordata->consWeight);
   graph = new HyperrowGraph<gcg::GraphTclique>(gcg, w);

   SCIP_CALL( graph->createFromPartialMatrix(partialdecdetectiondata->detprobdata, partialdec) );
   SCIP_CALL( createMetisFile(gcg, detectordata, partialdec->getID(), graph, tempfile) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");

   SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &temporaryClock) );

   for( j = 0, k = 0; k < maxnblockcandidates; ++k )
   {
      int nblocks = numberOfBlocks[k] - partialdec->getNBlocks();
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      SCIP_RETCODE retcode;

      if( nblocks > partialdec->getNOpenvars() || nblocks <= 1 )
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
         SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
         continue;
      }

      retcode = callMetis(gcg, detectordata, graph, tempfile, nblocks, result);

      if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY)
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
         SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
         continue;
      }

      if( allowopenpartialdecs )
      {
         SCIP_CALL( graph->createPartialdecFromPartition(partialdec, &newpartialdecs[j], &newpartialdecs[j+1], partialdecdetectiondata->detprobdata));
      }
      else
      {
         SCIP_CALL( graph->createPartialdecFromPartition(partialdec, &newpartialdecs[j], NULL, partialdecdetectiondata->detprobdata));
      }

      if( newpartialdecs[j] != NULL )
      {
         if( !allowopenpartialdecs )
         {
            newpartialdecs[j]->considerImplicits();
            newpartialdecs[j]->refineToBlocks();
            assert(newpartialdecs[j]->getNOpenconss() == 0);
            assert(newpartialdecs[j]->getNOpenvars() == 0);
         }
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock) );

         detectordata->found = TRUE;
         (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "hr\\_%d", numberOfBlocks[k]);
         newpartialdecs[j]->addDetectorChainInfo(decinfo);

         if( allowopenpartialdecs )
         {
            clockTimes.push_back(SCIPgetClockTime(scip, temporaryClock) / 2);
            clockTimes.push_back(SCIPgetClockTime(scip, temporaryClock) / 2);
            newpartialdecs[j + 1]->addDetectorChainInfo(decinfo);
            j += 2;
         }
         else
         {
            clockTimes.push_back(SCIPgetClockTime(scip, temporaryClock));
            j++;
         }
      }
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }
   delete graph;

   int nnewpartialdecs = j;
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d partialdecs found.\n",  nnewpartialdecs);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(partialdecdetectiondata->newpartialdecs), nnewpartialdecs) );
   partialdecdetectiondata->nnewpartialdecs = nnewpartialdecs;
   for( s = 0; s < nnewpartialdecs; ++s )
   {
      partialdecdetectiondata->newpartialdecs[s] = newpartialdecs[s];
      partialdecdetectiondata->newpartialdecs[s]->addClockTime(clockTimes[s] + SCIPgetClockTime(scip, temporaryClock) / nnewpartialdecs);
   }

   SCIPfreeMemoryArray(scip, &newpartialdecs);
   SCIP_CALL_ABORT( SCIPfreeClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPfreeClock(scip, &clock) );

   if( detectordata->tidy )
   {
      int status = remove( tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: %s", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }

   *result = detectordata->found ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}


static
GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecHrgpartition)
{
   SCIP_CLOCK* temporaryClock;
   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;
   SCIP_Bool enabledforlarge;
   int limit;
   SCIP* scip = GCGgetOrigprob(gcg);

   SCIPgetBoolParam(scip, "detection/detectors/hrgpartition/enabledforlargeproblems", &enabledforlarge);
   SCIPgetIntParam(scip, "detection/detectors/hrgpartition/limitnconssnvarsdefault", &limit);

   if ( !enabledforlarge && (SCIPgetNConss(scip) + SCIPgetNVars(scip) > limit) )
   {
      partialdecdetectiondata->detectiontime = 0;
      partialdecdetectiondata->nnewpartialdecs = 0;
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Started propagate partialdec of detector %s and partial decomp %d \n", DEC_NAME, partialdec->getID() );

   SCIP_CALL_ABORT( SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );

   partialdec->considerImplicits();
   partialdec->refineToMaster();

   if( partialdec->alreadyAssignedConssToBlocks() )
   {
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock) );
      SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
      partialdecdetectiondata->detectiontime = 0;
      partialdecdetectiondata->nnewpartialdecs = 0;
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   if( !connected(partialdecdetectiondata->detprobdata, partialdec) )
   {
      partialdec->assignSmallestComponentsButOneConssAdjacency();
   }

   detection(gcg, GCGdetectorGetData(detector), partialdecdetectiondata, partialdec, true, result);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock) );
   partialdecdetectiondata->detectiontime = SCIPgetClockTime(scip, temporaryClock);
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );

   return SCIP_OKAY;
}


#define finishPartialdecHrgpartition NULL
#define detectorPostprocessPartialdecHrgpartition NULL


static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveHrgpartition)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   int newval;
   SCIP_Real modifier;
   SCIP* scip = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE ) );

     (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxcallround", name);
   SCIP_CALL( SCIPgetIntParam(scip, setstr, &newval) );
   ++newval;
   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "After Setting %s = %d\n", setstr, newval);


   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origmaxcallround", name);
   SCIP_CALL( SCIPgetIntParam(scip, setstr, &newval) );
   ++newval;
   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );

   /* check if no problem is read */
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
      SCIP_CALL( SCIPsetIntParam(scip, setstr, DEFAULT_MAXNBLOCKCANDIDATES ) );
      SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, DEFAULT_MAXNBLOCKCANDIDATES);
      return SCIP_OKAY;
   }

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;

   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);
   modifier += 1;

   newval = MAX( 0, DEFAULT_MAXNBLOCKCANDIDATES - modifier + 2 );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultHrgpartition)
{
   char setstr[SCIP_MAXSTRLEN];
   int newval;
   SCIP_Real modifier;
   const char* name = GCGdetectorGetName(detector);
   SCIP* scip = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDFINISHING ) );

   /* check if no problem is read */
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
      SCIP_CALL( SCIPsetIntParam(scip, setstr, DEFAULT_MAXNBLOCKCANDIDATES ) );
      SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, DEFAULT_MAXNBLOCKCANDIDATES);
      return SCIP_OKAY;
   }

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;

   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);
   modifier += 1;

   newval = MAX( 0, DEFAULT_MAXNBLOCKCANDIDATES - modifier );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}

static
GCG_DECL_SETPARAMFAST(setParamFastHrgpartition)
{
   char setstr[SCIP_MAXSTRLEN];
   int newval;
   SCIP_Real modifier;
   const char* name = GCGdetectorGetName(detector);
   SCIP* scip = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   if ( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNConss(scip) + SCIPgetNVars(scip) < 6000 )
      SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );
   else SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   /* check if no problem is read */
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
      SCIP_CALL( SCIPsetIntParam(scip, setstr, DEFAULT_MAXNBLOCKCANDIDATES ) );
      SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, DEFAULT_MAXNBLOCKCANDIDATES);
      return SCIP_OKAY;
   }

   modifier = ( (SCIP_Real)SCIPgetNConss(scip) + (SCIP_Real)SCIPgetNVars(scip) ) / SET_MULTIPLEFORSIZETRANSF;

   modifier = log(modifier) / log(2);

   if (!SCIPisFeasPositive(scip, modifier) )
      modifier = -1.;

   modifier = SCIPfloor(scip, modifier);
   modifier += 1;

   newval = MAX( 0, DEFAULT_MAXNBLOCKCANDIDATES - modifier - 2 );
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
   SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
   SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}


/** creates the hrgpartition presolver and includes it in SCIP */
extern "C"
SCIP_RETCODE GCGincludeDetectorHrgpartition(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_DETECTORDATA *detectordata = NULL;
   SCIP* scip = GCGgetOrigprob(gcg);
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;

   SCIP_CALL( GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeHrgpartition, initHrgpartition, exitHrgpartition, propagatePartialdecHrgpartition, finishPartialdecHrgpartition, detectorPostprocessPartialdecHrgpartition, setParamAggressiveHrgpartition, setParamDefaultHrgpartition, setParamFastHrgpartition) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/limitnconssnvarsdefault", "Limit for sum of nvars and nconss for enabling this detector in default", &detectordata->limitnconssnvarsdefault, TRUE, DEFAULT_LIMITNCONSSNVARSDEFAULT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hrgpartition/enabledforlargeproblems", "Should this detector be enabled even the limit nconssnvars is exceeded",&detectordata->enabledforlargeproblems, TRUE, DEFAULT_ENABLEDFORLARGEPROBLEMS, NULL, NULL) );

     /* add hrgpartition presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/maxnblockcandidates", "The maximal number of block number candidates", &detectordata->maxnblockcandidates, FALSE, DEFAULT_MAXNBLOCKCANDIDATES, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hrgpartition/beta", "Factor on how heavy equality (beta) and inequality constraints are measured", &detectordata->beta, FALSE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hrgpartition/alpha", "Factor on how heavy the standard deviation of the coefficients is measured", &detectordata->alpha, FALSE, DEFAULT_ALPHA, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/varWeight", "Weight of a variable hyperedge", &detectordata->varWeight, FALSE, DEFAULT_VARWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/varWeightBinary", "Weight of a binary variable hyperedge", &detectordata->varWeightBinary, FALSE, DEFAULT_VARWEIGHTBIN, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/varWeightContinous", "Weight of a continuos variable hyperedge", &detectordata->varWeightContinous, FALSE, DEFAULT_VARWEIGHTCONT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/varWeightImplint", "Weight of a implicit integer variable hyperedge", &detectordata->varWeightImplint, FALSE, DEFAULT_VARWEIGHTIMPL, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/varWeightInteger", "Weight of a integer variable hyperedge", &detectordata->varWeightInteger, FALSE, DEFAULT_VARWEIGHTINT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hrgpartition/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/randomseed", "Random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hrgpartition/dummynodes", "Percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hrgpartition/consWeightSetppc", "Weight for constraint hyperedges that are setpartitioning or covering constraints", &detectordata->consWeightSetppc, FALSE, DEFAULT_CONSWEIGHT_SETPPC, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hrgpartition/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hrgpartition/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hrgpartition/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hrgpartition/realname", "Should the problem be used for metis files or a temporary name", &detectordata->realname, FALSE, DEFAULT_REALNAME, NULL, NULL) );

   return SCIP_OKAY;
}

#endif
