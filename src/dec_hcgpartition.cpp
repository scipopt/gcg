/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   dec_hcgpartition.cpp
 * @brief  arrowhead and bordered detector via graph partitioning (uses hmetis)
 * @ingroup DETECTORS
 * @author Martin Bergner
 * @author Michael Bastubbe
 *
 * Detects arrowhead (double bordered) decompositions as well as decompositions
 * with only linking variables or linking constraints.
 *
 * This detector needs hmetis and works only under Linux/MacOS, it further needs the Z-shell (zsh)
 * to enforce memory and time limits on hmetis as this is the only shell reliably doing that.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


/* #define SCIP_DEBUG */
#include "dec_hcgpartition.h"


#if !defined(_WIN32) && !defined(_WIN64)
#include <cassert>
#include <cstring>
#include <cerrno>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <algorithm>


#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "graph/matrixgraph.h"
#include "graph/hypercolgraph.h"
#include "graph/graph_tclique.h"
#include "graph/weights.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "scip/clock.h"


#include <set>

using gcg::HypercolGraph;
using gcg::MatrixGraph;
using gcg::Weights;

#define DEC_DETECTORNAME         "hcgpartition"    /**< name of the detector */
#define DEC_DESC                 "enforces arrowhead structures using graph partitioning" /**< description of detector */
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  0     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              1000           /**< priority of the detector */
#define DEC_DECCHAR               'G'            /**< display character of detector */
#define DEC_ENABLED               FALSE           /**< should detector be called by default */
#define DEC_ENABLEDORIGINAL       FALSE        /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE          /**< should detector be called by default */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE          /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL          TRUE           /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */


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
#define DEFAULT_MAXNBLOCKCANDIDATES 1          /**< number of block number candidates to be considered */
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

#define SET_MULTIPLEFORSIZETRANSF 12500

/*
 * Data structures
 */

/** private detector data */
struct DEC_DetectorData
{
   /* Graph stuff for hmetis */
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
   int       maxnblockcandidates;       /**< maximal number of block canddidates to test */
   int       maxblocks;       /**< maximal number of blocks to test */
   int       minblocks;       /**< minimal number of blocks to test */

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

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeHcgpartition)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detector initialization method (called after problem was transformed) */
static

DEC_DECL_INITDETECTOR(initHcgpartition)
{
   int nconss;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata->found = FALSE;

   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);


   return SCIP_OKAY;
}

/** detector deinitialization method (called before the transformed problem is freed) */
static

DEC_DECL_EXITDETECTOR(exitHcgpartition)
{
   assert(scip != NULL);


   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);


   return SCIP_OKAY;
}

/** will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata,       /**< detector data data structure */
   MatrixGraph<gcg::GraphTclique>* graph,    /**< the graph of the matrix */
   char                  tempfile[SCIP_MAXSTRLEN],
   int                   nblocks,            /** number of blocks */
   SCIP_RESULT*          result              /**< result indicating whether the detection was successful */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];

   SCIP_CLOCK* metisclock;

   int status;

   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPcreateWallClock(scip, &metisclock);
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
               tempfile,
               nblocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }
   else
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
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
   status = system( metiscall );

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
      assert(false);
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
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata,        /**< detector data structure */
   int                   seeedID,             /**< used for speaking filenames */
   MatrixGraph<gcg::GraphTclique>* graph,    /**< the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN]
   )
{
   int nvertices;
   int ndummyvertices;
   int fd;
   nvertices = graph->getNNonzeroes();
   /*lint --e{524}*/
   ndummyvertices = SCIPceil(scip, detectordata->dummynodes*nvertices);
   graph->setDummynodes(ndummyvertices);

   if( !detectordata->realname )
   {
      (void) SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-%c-%d.metis.XXXXXX", DEC_DECCHAR, seeedID );
   }
   else
   {
      (void) SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-%s-%c-%d.metis.XXXXXX", SCIPgetProbName(scip), DEC_DECCHAR, seeedID);
   }

   fd = mkstemp(tempfile);

   SCIP_CALL( graph->writeToFile(fd, TRUE) );
   close(fd);
   return SCIP_OKAY;
}

/** returns, whether the hypercolgraph is connected */
static
bool connected(
   gcg::Seeedpool*  seeedpool,
   gcg::Seeed*      seeed
   )
{
   std::vector<int> queue;
   std::vector<int> visited;
   std::vector<bool> inqueue(seeedpool->getNConss(), false);
   std::vector<bool> isvisited(seeedpool->getNConss(), false);
   int start = -1;

   if(seeed->getNOpenconss() < 2)
      return false;

   start = seeed->getOpenconss()[0];

   queue.push_back(start);
   inqueue[start] = true;
   do
   {
      int node = queue[0];
      queue.erase(queue.begin());
      inqueue[node] = false;
      visited.push_back(node);
      isvisited[node] = true;
      for(int v = 0; v < seeedpool->getNVarsForCons(node); ++v)
      {
         int var = seeedpool->getVarsForCons(node)[v];
         if(!seeed->isVarOpenvar(var))
            continue;
         for(int c = 0; c < seeedpool->getNConssForVar(var); ++c)
         {
            int cons = seeedpool->getConssForVar(var)[c];
            if(!seeed->isConsOpencons(cons))
               continue;
            if( isvisited[cons] )
               continue;
            if( inqueue[cons] )
               continue;
            queue.push_back(cons);
            inqueue[cons] = true;
         }
      }
   } while(!queue.empty());

   if((int)visited.size() != seeed->getNOpenconss())
      return false;
   else
      return true;
}

/** detection function for seeeds */
static
SCIP_RETCODE detection(
   SCIP*                   scip,                         /**< SCIP data structure */
   DEC_DETECTORDATA*       detectordata,                 /**< detectordata of the detector */
   Seeed_Propagation_Data* seeedPropagationData,         /**< seeedPropagationData (including the seeedpool) where to store the new Seeeds */
   gcg::Seeed*             seeed,                        /**< seeed to propagate */
   bool                    border,                       /**< whether new seeeds should be stored in which this detector only assignes conss to master */
   SCIP_RESULT*            result                        /**< pointer where to store the result */
)
{

   /* add hcgpartition presolver parameters */
   char setstr[SCIP_MAXSTRLEN];
   char decinfo[SCIP_MAXSTRLEN];
   int maxnblockcandidates;
   int k;
   int j;
   int s;
   int nMaxSeeeds;
   int nNewSeeeds = 0;
   gcg::Seeed** newSeeeds;
   SCIP_CLOCK* clock;
   SCIP_CLOCK* temporaryClock;
   std::vector<SCIP_Real> clockTimes;        /**< vector containing times in seconds  */

   /* Graph stuff for hmetis */
   MatrixGraph<gcg::GraphTclique>* graph;    /**< the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN];            /**< filename for the metis input file */


   SCIP_CALL_ABORT( SCIPcreateClock(scip, &clock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock) );

   *result = SCIP_DIDNOTFIND;

   std::vector<int> numberOfBlocks = seeedPropagationData->seeedpool->getSortedCandidatesNBlocks();
   if(numberOfBlocks.empty())
      numberOfBlocks.push_back(8);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/hcgpartition/maxnblockcandidates");
   SCIP_CALL( SCIPgetIntParam(scip, setstr, &maxnblockcandidates) );

   maxnblockcandidates = MIN(maxnblockcandidates, (int) numberOfBlocks.size() );

   assert(scip != NULL);
   assert(detectordata != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   nMaxSeeeds = detectordata->maxblocks-detectordata->minblocks+1;

   /* allocate space for output data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newSeeeds), 2 * nMaxSeeeds) );

   /* build the hypergraph structure from the original problem */

   Weights w(detectordata->varWeight, detectordata->varWeightBinary, detectordata->varWeightContinous,detectordata->varWeightInteger,detectordata->varWeightInteger,detectordata->consWeight);
   graph = new HypercolGraph<gcg::GraphTclique>(scip, w);

   SCIP_CALL( graph->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed) );
   SCIP_CALL( createMetisFile(scip, detectordata, seeed->getID(), graph, tempfile) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");


   SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
   SCIP_CALL_ABORT( SCIPcreateClock(scip, &temporaryClock) );

   for( j = 0, k = 0; k < maxnblockcandidates; ++k)
   {
      int nblocks = numberOfBlocks[k] - seeed->getNBlocks();
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      SCIP_RETCODE retcode;
      if(nblocks > seeed->getNOpenconss() || nblocks <= 0)
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
         SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
         continue;
      }

      retcode = callMetis(scip, detectordata, graph, tempfile, nblocks, result);

      if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY)
      {
         SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
         SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
         continue;
      }

      SCIP_CALL( graph->createSeeedFromPartition(seeed, &newSeeeds[j], &newSeeeds[j+1], seeedPropagationData->seeedpool));
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      if( (newSeeeds)[j] != NULL )
      {
         nNewSeeeds = nNewSeeeds + 2;
         detectordata->found = TRUE;
         (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "hc\\_%d", numberOfBlocks[k]);
         newSeeeds[j]->addDetectorChainInfo(decinfo);
         newSeeeds[j+1]->addDetectorChainInfo(decinfo);

         clockTimes.push_back(SCIPclockGetTime(temporaryClock));
         clockTimes.push_back(SCIPclockGetTime(temporaryClock)); // 2x because two seeeds where created
      }
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
      j = j + 2;
   }

   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, clock ) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d seeeds found.\n",  nNewSeeeds);

   delete graph;
   graph = NULL;

   assert(nNewSeeeds % 2 == 0);
   if(border)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), nNewSeeeds) );
      seeedPropagationData->nNewSeeeds = nNewSeeeds;
      for(j = 0, s = 0; s < nNewSeeeds; ++j)
      {
         if(newSeeeds[j] != NULL)
         {
            seeedPropagationData->newSeeeds[s] = newSeeeds[j];
            ++s;
         }
      }
   }
   else
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), nNewSeeeds/2) );
      seeedPropagationData->nNewSeeeds = nNewSeeeds/2;
      for(j = 0, s = 0; s < nNewSeeeds/2; j+=2)
      {
         if(newSeeeds[j] != NULL)
         {
            seeedPropagationData->newSeeeds[s] = newSeeeds[j];
            ++s;
         }
      }
   }


   SCIPfreeMemoryArray(scip, &newSeeeds);

   if( detectordata->tidy )
   {
      int status = remove( tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
         SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
         SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );
         return SCIP_WRITEERROR;
      }
   }

   SCIP_CALL_ABORT( SCIPstopClock(scip, clock ) );
   if(border)
   {
      for( s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
         seeedPropagationData->newSeeeds[s]->addClockTime( SCIPclockGetTime(clock) + clockTimes[s] );
   }
   else
   {
      for( s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
         seeedPropagationData->newSeeeds[s]->addClockTime( SCIPclockGetTime(clock) + clockTimes[2*s] );
   }
   SCIP_CALL_ABORT(SCIPfreeClock(scip, &clock) );

   *result = detectordata->found ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}

/** Merge of propagation and finishing method from decomposition_toolbox */
static
SCIP_RETCODE fromToolbox(
   SCIP_Bool propagate,  /**< Shall the instance be propagated or finished */
   SCIP* scip, 
   DEC_DETECTOR* detector, 
   SEEED_PROPAGATION_DATA* seeedPropagationData, 
   SCIP_RESULT* result, 
   SCIP_DIALOGHDLR* dialoghdlr,
   SCIP_DIALOG* dialog
   )
{
   /* add hcgpartition presolver parameters */
   char decinfo[SCIP_MAXSTRLEN];
   gcg::Seeed** newSeeeds;
   DEC_DETECTORDATA* detectordata;
   gcg::Seeed* seeed;
   SCIP_RETCODE retcode;
   char* command;
   int commandlen;
   SCIP_Bool endoffile;
   int nblocks;
   /* Graph stuff for hmetis */
   MatrixGraph<gcg::GraphTclique>* graph;    /**< the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN];            /**< filename for the metis input file */

   seeed = seeedPropagationData->seeedToPropagate;
   detectordata = DECdetectorGetData(detector);

   *result = SCIP_DIDNOTFIND;

   std::vector<int> numberOfBlocks = seeedPropagationData->seeedpool->getSortedCandidatesNBlocks();
   if(numberOfBlocks.empty())
      numberOfBlocks.push_back(8);

   int nconss = seeedPropagationData->seeedpool->getNConss();
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);

   assert(scip != NULL);
   assert(detectordata != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);

   /* allocate space for output data */
   assert(detectordata->maxblocks >= detectordata->minblocks);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newSeeeds), 2 ) );

   /* build the hypergraph structure from the original problem */

   Weights w(detectordata->varWeight, detectordata->varWeightBinary, detectordata->varWeightContinous,detectordata->varWeightInteger,detectordata->varWeightInteger,detectordata->consWeight);
   graph = new HypercolGraph<gcg::GraphTclique>(scip, w);

   SCIP_CALL( graph->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed) );
   SCIP_CALL( createMetisFile(scip, detectordata, seeed->getID(), graph, tempfile) );

   detectordata->metisubfactor = DEFAULT_METIS_UBFACTOR; //@TODO: resolve s.t. this parameter does not have to be set manually here
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:\n");

   SCIPinfoMessage(scip, NULL, "Maximal number of blocks the decomposition can contain: %d\n", seeed->getNOpenconss());
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Type in the number of blocks that the decomposition should contain (e.g. \"5\") \nGCG/toolbox> ", &command, &endoffile) );
   commandlen = strlen(command);

   if( commandlen != 0 )
   {
      nblocks = atoi(command);
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Invalid input!\n");
      return SCIP_OKAY;
   }

   if(nblocks > seeed->getNOpenconss() || nblocks <= 0)
   {
      SCIPinfoMessage(scip, NULL, "Invalid number of blocks, choose at most %d\n", seeed->getNOpenconss());
      return SCIP_OKAY;
   }

   retcode = callMetis(scip, detectordata, graph, tempfile, nblocks, result);

   if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY)
   {
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   if( detectordata->tidy )
   {
      int status = remove( tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }

   SCIP_CALL( graph->createSeeedFromPartition(seeed, &newSeeeds[0], &newSeeeds[1], seeedPropagationData->seeedpool));
   delete graph;
   graph = NULL;
   if( (newSeeeds)[0] != NULL && !propagate ) //finishing successful
   {
      detectordata->found = TRUE;
      (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "hc\\_%d", numberOfBlocks[0]);
      newSeeeds[0]->addDetectorChainInfo(decinfo);
      seeedPropagationData->newSeeeds[0] = (newSeeeds)[0];
      ++(seeedPropagationData->nNewSeeeds);
      seeedPropagationData->newSeeeds[0]->setFinishingDetectorPropagated(detector);
      SCIPfreeMemoryArray(scip, &newSeeeds);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }
   else if( (newSeeeds)[1] != NULL && propagate ) //propagation successful
   {
      detectordata->found = TRUE;
      (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "hc\\_%d", numberOfBlocks[0]);
      newSeeeds[1]->addDetectorChainInfo(decinfo);
      seeedPropagationData->newSeeeds[0] = (newSeeeds)[1];
      ++(seeedPropagationData->nNewSeeeds);
      seeedPropagationData->newSeeeds[0]->setDetectorPropagated(detector);
      SCIPfreeMemoryArray(scip, &newSeeeds);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }
   else //propagation/finishing unsuccessful
   {
      SCIPfreeMemoryArray(scip, &newSeeeds);
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }
}

/** detector structure detection method, tries to detect a structure in the problem */
static
DEC_DECL_DETECTSTRUCTURE(detectHcgpartition)
{
   int i;
   int j;
   int ndecs;

   MatrixGraph<gcg::GraphTclique>* graph;    /**< the graph of the matrix */
   char tempfile[SCIP_MAXSTRLEN];            /**< filename for the metis input file */

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
   graph = new HypercolGraph<gcg::GraphTclique>(scip, w);

   SCIP_CALL( graph->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( createMetisFile(scip, detectordata, 0, graph, tempfile) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");
   for( j = 0, i = detectordata->minblocks; i <= detectordata->maxblocks; ++i )
   {
      SCIP_RETCODE retcode;
      /* get the partitions for the new variables from metis */
      retcode = callMetis(scip, detectordata, graph, tempfile, i, result);

      if( *result != SCIP_SUCCESS || retcode != SCIP_OKAY )
      {
         continue;
      }

      SCIP_CALL( graph->createDecompFromPartition(&(*decdecomps)[j]) );
      if( (*decdecomps)[j] != NULL )
      {
         *ndecdecomps += 1;
         ++j;
         detectordata->found = TRUE;
      }
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d decompositions found.\n",  *ndecdecomps);

   delete graph;
   graph = NULL;

   SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );

   if( detectordata->tidy )
   {
      int status = remove( tempfile );
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
DEC_DECL_PROPAGATESEEED(propagateSeeedHcgpartition)
{
   gcg::Seeed* seeed;
   seeed = seeedPropagationData->seeedToPropagate;

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );


  //assert( seeed->checkConsistency( seeedPropagationData->seeedpool ));

   seeed->considerImplicits();
  //assert( seeed->checkConsistency( seeedPropagationData->seeedpool ));

   seeed->refineToMaster();

   //assert( seeed->checkConsistency( seeedPropagationData->seeedpool ));

   if(!connected(seeedPropagationData->seeedpool, seeed) || seeed->alreadyAssignedConssToBlocks() )
   {
      seeed->assignSmallestComponentsButOneConssAdjacency();
   }

   detection(scip, DECdetectorGetData(detector), seeedPropagationData, seeed, TRUE, result);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );


   for( int s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
   {
      seeedPropagationData->newSeeeds[s]->addClockTime( SCIPclockGetTime(temporaryClock )  );
   }

   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );


   return SCIP_OKAY;
}

static
DEC_DECL_PROPAGATEFROMTOOLBOX(propagateFromToolboxHcgpartition)
{
   return fromToolbox(TRUE, scip, detector, seeedPropagationData, result, dialoghdlr, dialog );
}

static
DEC_DECL_FINISHFROMTOOLBOX(finishFromToolboxHcgpartition)
{
   return fromToolbox(FALSE, scip, detector, seeedPropagationData, result, dialoghdlr, dialog );
}

static
DEC_DECL_FINISHSEEED(finishSeeedHcgpartition)
{

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );


   gcg::Seeed* seeed = seeedPropagationData->seeedToPropagate;

   seeed->considerImplicits();
   seeed->refineToBlocks();

   if( !connected(seeedPropagationData->seeedpool, seeed ) )
   {
      seeed->assignSmallestComponentsButOneConssAdjacency();
   }

   detection(scip, DECdetectorGetData(detector), seeedPropagationData, seeed, FALSE, result);

   SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );


   for( int s = 0; s < seeedPropagationData->nNewSeeeds; ++s )
   {
      seeedPropagationData->newSeeeds[s]->considerImplicits();
      seeedPropagationData->newSeeeds[s]->refineToBlocks();
      seeedPropagationData->newSeeeds[s]->addClockTime( SCIPclockGetTime(temporaryClock )  );
      assert(seeedPropagationData->newSeeeds[s]->getNOpenconss() == 0);
      assert(seeedPropagationData->newSeeeds[s]->getNOpenvars() == 0);
   }

   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );


   return SCIP_OKAY;
}

#define detectorPostprocessSeeedHcgpartition NULL

static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveHcgpartition)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = DECdetectorGetName(detector);
   int newval;
   SCIP_Real modifier;

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
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
   SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, newval);

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxnblockcandidates", name);
      SCIP_CALL( SCIPsetIntParam(scip, setstr, newval ) );
      SCIPinfoMessage(scip, NULL, "%s = %d\n", setstr, newval);
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
DEC_DECL_SETPARAMDEFAULT(setParamDefaultHcgpartition)
{
   char setstr[SCIP_MAXSTRLEN];
   int newval;
   SCIP_Real modifier;

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDORIGINAL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, DEC_ENABLEDFINISHING ) );

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
     {
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
DEC_DECL_SETPARAMFAST(setParamFastHcgpartition)
{
   char setstr[SCIP_MAXSTRLEN];
   int newval;
   SCIP_Real modifier;


   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

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



/** creates the hcgpartition presolver and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorHcgpartition(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#if !defined(_WIN32) && !defined(_WIN64)
   DEC_DETECTORDATA *detectordata = NULL;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);


   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE, detectordata, detectHcgpartition, freeHcgpartition, initHcgpartition, exitHcgpartition, propagateSeeedHcgpartition, propagateFromToolboxHcgpartition, finishFromToolboxHcgpartition, finishSeeedHcgpartition, detectorPostprocessSeeedHcgpartition, setParamAggressiveHcgpartition, setParamDefaultHcgpartition, setParamFastHcgpartition) );


   /* add hcgpartition detector parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/maxnblockcandidates", "The maximal number of block number candidates", &detectordata->maxnblockcandidates, FALSE, DEFAULT_MAXNBLOCKCANDIDATES, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/maxblocks", "The maximal number of blocks (only used in legacy mode; detector is called for all block numbers in [minblocks,maxblocks])", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/minblocks", "The minimal number of blocks (only used in legacy mode; detector is called for all block numbers in [minblocks,maxblocks])", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hcgpartition/beta", "Factor on how heavy equality (beta) and inequality constraints are measured", &detectordata->beta, FALSE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hcgpartition/alpha", "Factor on how heavy the standard deviation of the coefficients is measured", &detectordata->alpha, FALSE, DEFAULT_ALPHA, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/varWeight", "Weight of a variable hyperedge", &detectordata->varWeight, FALSE, DEFAULT_VARWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/varWeightBinary", "Weight of a binary variable hyperedge", &detectordata->varWeightBinary, FALSE, DEFAULT_VARWEIGHTBIN, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/varWeightContinous", "Weight of a continuos variable hyperedge", &detectordata->varWeightContinous, FALSE, DEFAULT_VARWEIGHTCONT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/varWeightImplint", "Weight of a implicit integer variable hyperedge", &detectordata->varWeightImplint, FALSE, DEFAULT_VARWEIGHTIMPL, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/varWeightInteger", "Weight of a integer variable hyperedge", &detectordata->varWeightInteger, FALSE, DEFAULT_VARWEIGHTINT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hcgpartition/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/randomseed", "Random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hcgpartition/dummynodes", "Percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/hcgpartition/consWeightSetppc", "Weight for constraint hyperedges that are setpartitioning or covering constraints", &detectordata->consWeightSetppc, FALSE, DEFAULT_CONSWEIGHT_SETPPC, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detection/detectors/hcgpartition/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hcgpartition/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hcgpartition/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/hcgpartition/realname", "Should the problem be used for metis files or a temporary name", &detectordata->realname, FALSE, DEFAULT_REALNAME, NULL, NULL) );

#endif
   return SCIP_OKAY;
}
