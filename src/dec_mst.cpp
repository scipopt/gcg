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

/**@file   dec_mst.cpp
 * @ingroup DETECTORS
 * @brief  detector MST
 * @author Igor Pesic
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <time.h>     // for measuring time performance

#include "dec_mst.h"
#include "cons_decomp.h"
#include "graph/matrixgraph.h"
#include "graph/rowgraph_weighted.h"
#include "graph/graph_gcg.h"
#include "scip/clock.h"

using gcg::RowGraphWeighted;
using gcg::Weights;
using gcg::GraphGCG;


/* constraint handler properties */
#define DEC_DETECTORNAME          "mst"                               /**< name of detector */
#define DEC_DESC                  "detector based on MST clustering"  /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              910         /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'M'         /**< display character of detector */
#define DEC_ENABLED               FALSE        /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       FALSE        /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

/* Default parameter settings*/
#define DEFAULT_N_ITERATIONS              51
#define DEFAULT_JOHNSON_ENABLE            true
#define DEFAULT_INTERSECTION_ENABLE       false
#define DEFAULT_JACCARD_ENABLE            false
#define DEFAULT_COSINE_ENABLE             false
#define DEFAULT_SIMPSON_ENABLE            false
#define DEFAULT_POSTPROC_ENABLE           true
#define MAX_N_BLOCKS                      100

/*
 * Data structures
 */

/** detector handler data */
struct DEC_DetectorData
{
   std::vector< RowGraphWeighted<GraphGCG>*> *graphs;    /**< the graph of the matrix */
   SCIP_RESULT result;                                 /**< result pointer to indicate success or failure */
   SCIP_Bool found;
   int n_iterations;
   int n_similarities;                                  /**< number of active similarities */
   SCIP_Bool johnsonenable;
   SCIP_Bool intersectionenable;
   SCIP_Bool jaccardenable;
   SCIP_Bool cosineenable;
   SCIP_Bool simpsonenable;
   SCIP_Bool postprocenable;
};


/*
 * Local methods
 */

static std::vector<double> getEpsList(int length, double mid, bool isintersection)
{
   int n1, n2;
   if(isintersection)
   {
      n1 = (int) round((length+1) / 2.0);
      n2 = n1;
   }
   else
   {
      n2 = (int) round((length+1) / 4.0);
      n1 = abs(length - n2) + 1;
   }

   double s = mid;
   double end1 = mid + 0.9;      // lower boundary
   double end2 = mid + 0.4;      // upper boundary

   double q1 = pow(end1/s, 1.0/(double)(n1-1));
   double q2 = pow(end2/s, 1.0/(double)(n2-1));

   std::vector<double> geom_seq1(n1-1);
   std::vector<double> geom_seq2(n2);

   int j = 0;
   for(int i = n1 - 1; i > 0; i--)
   {
      geom_seq1[j] = 2*s-s*pow(q1, (double)i);
      j++;
   }
   for(int i = 0; i < n2; i++)
   {
      geom_seq2[i] = s*pow(q2, (double)i);
   }

   geom_seq1.insert( geom_seq1.end(), geom_seq2.begin(), geom_seq2.end() );

   assert((int)geom_seq1.size() == length);

   return geom_seq1;
}

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeMST)
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
static
DEC_DECL_EXITDETECTOR(exitMST)
{
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   delete detectordata->graphs;

   return SCIP_OKAY;
}


/** detection initialization function of detector (called before solving is about to begin) */
static
DEC_DECL_INITDETECTOR(initMST)
{  /*lint --e{715}*/

   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);


   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata->n_similarities = -1;
   detectordata->found = FALSE;
   detectordata->graphs = new std::vector< RowGraphWeighted<GraphGCG>*>();
   return SCIP_OKAY;
}

/** are there conss and vars to be included by the graph and have the conss common vars included by the graph */
static
bool graphCompletible(
   gcg::Seeedpool*  seeedpool,
   gcg::Seeed*      seeed
   )
{
   bool completible;

   //have the open conss open vars?
   for(int c = 0; c < seeed->getNOpenconss() && !completible; ++c)
   {
      int cons = seeed->getOpenconss()[c];
      for(int v = 0; v < seeed->getNOpenvars() && !completible; ++v)
      {
         int var = seeed->getOpenvars()[v];
         for(int i = 0; i < seeedpool->getNVarsForCons(cons) && !completible; ++i)
         {
            if(var == seeedpool->getVarsForCons(cons)[i])
            {
               completible = true;
            }
         }
      }
   }
   if(!completible)
      return false;

   //have the open conss common open vars?
   for(int c = 0; c < seeed->getNOpenconss(); ++c)
   {
      int cons1 = seeed->getOpenconss()[c];
      for(int d = c + 1; d < seeed->getNOpenconss(); ++d)
      {
         int cons2 = seeed->getOpenconss()[d];
         for(int v = 0; v < seeedpool->getNVarsForCons(cons1); ++v)
         {
            int var1 = seeedpool->getVarsForCons(cons1)[v];
            if(!seeed->isVarOpenvar(var1))
               continue;
            for(int w = 0; w < seeedpool->getNVarsForCons(cons2); ++w)
            {
               int var2 = seeedpool->getVarsForCons(cons2)[w];
               if(var1 == var2)
                  return true;
            }
         }
      }
   }
   return false;
}

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectMST)
{ /*lint --e{715}*/

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   *result = SCIP_DIDNOTFIND;

   Weights w(1, 1, 1, 1, 1, 1);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting MST structure:");

   time_t start, cp0, cp1, d_s, d_e;
   time(&start);

   std::vector<std::string> sim;

   if(detectordata->johnsonenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::JOHNSON, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Johnson");
   }
   if(detectordata->intersectionenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::INTERSECTION, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Intersection");
   }
   if(detectordata->jaccardenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::JACCARD, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Jaccard");
   }
   if(detectordata->cosineenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::COSINE, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Cosine");
   }
   if(detectordata->simpsonenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::SIMPSON, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Simspon");
   }
   time(&cp0);
   detectordata->n_similarities = (int) detectordata->graphs->size();

   double q = 10; // quantile to search for the percentile needed for the mid of the eps list
   std::vector<double> mids(detectordata->graphs->size());      // middle values for each eps list
   std::vector<std::vector<double> > epsLists(detectordata->graphs->size());
   for(int i = 0; i < (int)detectordata->graphs->size(); i++)
   {
      mids[i] = detectordata->graphs->at(i)->getEdgeWeightPercentile(q);
      if(i == 1 && detectordata->intersectionenable)
      {
         epsLists[i] = getEpsList(detectordata->n_iterations, mids[i], true); // case for intersection
      }
      else
      {
         epsLists[i] = getEpsList(detectordata->n_iterations, mids[i], false); // case for all except intersection
      }

   }

   time(&cp1);

   int max_ndecs = detectordata->n_iterations * detectordata->graphs->size();
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, max_ndecs) );

   const int max_blocks = std::min((int)round(0.3 * SCIPgetNConss(scip)), MAX_N_BLOCKS);

   *ndecdecomps = 0;
   time(&d_s);
   for(int i = 0; i < (int)detectordata->graphs->size(); i++)
   {
      RowGraphWeighted<GraphGCG>* graph = detectordata->graphs->at(i);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n  %s similarity: ", sim[i].c_str());
      int old_n_blocks = -1;
      int old_non_cl = -1;
      for(int j = 0; j < (int)epsLists[i].size() ; j++ )
      {
         double eps = epsLists[i][j];
         if(eps <= 0.0)
         {
            continue;
         }
         if(eps >= 1.0)
         {
            break;
         }

         // run MST with different eps
         SCIP_CALL( graph->computePartitionMST(eps, detectordata->postprocenable) );


         int n_blocks;
         SCIP_CALL( graph->getNBlocks(n_blocks) );
         int non_cl;

         SCIP_CALL( graph->nonClustered(non_cl) );

         // skip the case if we have too many blocks (it means we must increase eps) or if the clustering is the same as the last one
         if( n_blocks > max_blocks || n_blocks == 0 || (n_blocks == old_n_blocks && non_cl == old_non_cl) )
         {
            continue;
         }
         // stop. eps is already too big
         if( n_blocks == 1 && non_cl == 0)
         {
            break;
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n    Blocks: %d, Master Conss: %d/%d, ", n_blocks, non_cl, SCIPgetNConss(scip));
         old_n_blocks = n_blocks;
         old_non_cl = non_cl;


         SCIP_CALL( graph->createDecompFromPartition(&(*decdecomps)[*ndecdecomps]) );

         auto check = DECdecompGetNLinkingvars((*decdecomps)[*ndecdecomps]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Link Vars: %d. ", check);

         if( (*decdecomps)[*ndecdecomps] != NULL )
         {
            *ndecdecomps += 1;
            detectordata->found = TRUE;
         }
      }
      delete detectordata->graphs->at(i);
      detectordata->graphs->at(i) = NULL;
   }

   detectordata->graphs->clear();
   time(&d_e);
   double elapsed_graphs = difftime(cp0, start);
   double elapsed_mst = difftime(d_e, d_s);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d similarities used, %d decompositions found.\n", detectordata->n_similarities, *ndecdecomps);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "MST runtime: graphs: %.2lf, mst: %.2lf. \n", elapsed_graphs, elapsed_mst);

   SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );

   *result = *ndecdecomps > 0 ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   if( *ndecdecomps == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, decdecomps);
   }
   return SCIP_OKAY;
}

//#define propagateSeeedMST NULL

static
DEC_DECL_PROPAGATESEEED(propagateSeeedMST)
{ /*lint --e{715}*/

   int nNewSeeeds;
   gcg::Seeed* seeed;
   gcg::Seeed** newSeeeds;
   DEC_DETECTORDATA* detectordata = DECdetectorGetData(detector);
   std::vector<SCIP_Real> clockTimes1;        /**< vector containing times in seconds  */
   std::vector<SCIP_Real> clockTimes2;        /**< vector containing times in seconds  */
   std::vector<SCIP_Real> clockTimes3;        /**< vector containing times in seconds  */

   assert(scip != NULL);
   assert(detectordata != NULL);
   *result = SCIP_DIDNOTFIND;

   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate);
   seeed->refineToBlocks();

   if(!graphCompletible(seeedPropagationData->seeedpool, seeed))
   {
      delete seeed;
      seeedPropagationData->nNewSeeeds = 0;
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }
   Weights w(1, 1, 1, 1, 1, 1);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting MST structure:");

   time_t start, cp0, cp1, d_s, d_e;
   time(&start);

   std::vector<std::string> sim;

   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(scip, &temporaryClock) );

   if(detectordata->johnsonenable)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::JOHNSON, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Johnson");
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }
   if(detectordata->intersectionenable)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::INTERSECTION, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Intersection");
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }
   if(detectordata->jaccardenable)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::JACCARD, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Jaccard");
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }
   if(detectordata->cosineenable)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed,gcg::DISTANCE_MEASURE::COSINE, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Cosine");
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }
   if(detectordata->simpsonenable)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed,gcg::DISTANCE_MEASURE::SIMPSON, gcg::WEIGHT_TYPE::DIST));
      detectordata->graphs->push_back(g);
      sim.push_back("Simspon");
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes1.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }
   time(&cp0);
   detectordata->n_similarities = (int) detectordata->graphs->size();

   double q = 10; // quantile to search for the percentile needed for the mid of the eps list
   std::vector<double> mids(detectordata->graphs->size());      // middle values for each eps list
   std::vector<std::vector<double> > epsLists(detectordata->graphs->size());
   for(int i = 0; i < (int)detectordata->graphs->size(); i++)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      mids[i] = detectordata->graphs->at(i)->getEdgeWeightPercentile(q);
      if(i == 1 && detectordata->intersectionenable)
      {
         epsLists[i] = getEpsList(detectordata->n_iterations, mids[i], true); // case for intersection
      }
      else
      {
         epsLists[i] = getEpsList(detectordata->n_iterations, mids[i], false); // case for all except intersection
      }
      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes2.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes2.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }

   time(&cp1);

   int nMaxSeeeds = detectordata->n_iterations * detectordata->graphs->size();
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newSeeeds), 2 * nMaxSeeeds) );

   const int max_blocks = std::min((int)round(0.3 * SCIPgetNConss(scip)), MAX_N_BLOCKS);

   int n_seeeds_found = 0;
   nNewSeeeds = 0;
   time(&d_s);
   for(int i = 0; i < (int)detectordata->graphs->size(); i++)
   {
      SCIP_CALL_ABORT( SCIPstartClock(scip, temporaryClock) );
      RowGraphWeighted<GraphGCG>* graph = detectordata->graphs->at(i);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n  %s similarity: ", sim[i].c_str());
      int old_n_blocks = -1;
      int old_non_cl = -1;
      for(int j = 0; j < (int)epsLists[i].size() ; j++ )
      {
         double eps = epsLists[i][j];
         if(eps <= 0.0)
         {
            continue;
         }
         if(eps >= 1.0)
         {
            break;
         }

         // run MST with different eps
         SCIP_CALL( graph->computePartitionMSTForPartialGraph(seeedPropagationData->seeedpool, seeed, eps, detectordata->postprocenable) );


         int n_blocks;
         SCIP_CALL( graph->getNBlocks(n_blocks) );
         int non_cl;

         SCIP_CALL( graph->nonClustered(non_cl) );

         // skip the case if we have too many blocks (it means we must increase eps) or if the clustering is the same as the last one
         if( n_blocks > max_blocks || n_blocks == 0 || (n_blocks == old_n_blocks && non_cl == old_non_cl) )
         {
            continue;
         }
         // stop. eps is already too big
         if( n_blocks == 1 && non_cl == 0)
         {
            break;
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n    Blocks: %d, Master Conss: %d/%d, ", n_blocks, non_cl, SCIPgetNConss(scip));
         old_n_blocks = n_blocks;
         old_non_cl = non_cl;


         SCIP_CALL( graph->createSeeedFromPartition(seeed,&newSeeeds[n_seeeds_found], &newSeeeds[n_seeeds_found+1], seeedPropagationData->seeedpool));

         if((newSeeeds)[n_seeeds_found] != NULL)
         {
            nNewSeeeds += 2;
            detectordata->found = TRUE;
         }
         n_seeeds_found += 2;
      }
      delete detectordata->graphs->at(i);
      detectordata->graphs->at(i) = NULL;

      SCIP_CALL_ABORT( SCIPstopClock(scip, temporaryClock ) );
      clockTimes3.push_back(SCIPclockGetTime(temporaryClock));
      clockTimes3.push_back(SCIPclockGetTime(temporaryClock));
      SCIP_CALL_ABORT( SCIPresetClock(scip, temporaryClock ) );
   }

   delete seeed;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), nNewSeeeds) );
   seeedPropagationData->nNewSeeeds = nNewSeeeds;
   for(int j = 0, s = 0; s < nNewSeeeds; ++j)
   {
      if(newSeeeds[j] != NULL)
      {
         seeedPropagationData->newSeeeds[s] = newSeeeds[j];
         seeedPropagationData->newSeeeds[s]->addClockTime(clockTimes1[j] + clockTimes2[j] + clockTimes3[j]);
         seeedPropagationData->newSeeeds[s]->setDetectorPropagated(detector);
         ++s;
      }
   }
   SCIPfreeMemoryArray(scip, &newSeeeds);

   detectordata->graphs->clear();
   time(&d_e);
   double elapsed_graphs = difftime(cp0, start);
   double elapsed_mst = difftime(d_e, d_s);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d similarities used, %d decompositions found.\n", detectordata->n_similarities, nNewSeeeds);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "MST runtime: graphs: %.2lf, mst: %.2lf. \n", elapsed_graphs, elapsed_mst);

   *result = nNewSeeeds > 0 ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   if( nNewSeeeds == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, &(seeedPropagationData->newSeeeds));
   }

   SCIP_CALL_ABORT(SCIPfreeClock(scip, &temporaryClock) );
   return SCIP_OKAY;
}

#define finishSeeedMST NULL

#define detectorPostprocessSeeedMST NULL

#define setParamAggressiveMST NULL
#define setParamDefaultMST NULL
#define setParamFastMST NULL


/*
 * detector specific interface methods
 */

/** creates the handler for xyz detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorMST(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#if !defined(_WIN32) && !defined(_WIN64)
   DEC_DETECTORDATA *detectordata = NULL;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE,
      detectordata, detectMST, freeMST, initMST, exitMST, propagateSeeedMST, NULL, NULL, finishSeeedMST, detectorPostprocessSeeedMST, setParamAggressiveMST, setParamDefaultMST, setParamFastMST) );

   /* add arrowheur presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/mst/niterations", "Number of iterations to run mst with different eps.", &detectordata->n_iterations, FALSE, DEFAULT_N_ITERATIONS, 11, 1001, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/mst/johson", "Enable johson distance measure.", &detectordata->johnsonenable, FALSE, DEFAULT_JOHNSON_ENABLE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/mst/intersection", "Enable intersection distance measure.", &detectordata->intersectionenable, FALSE, DEFAULT_INTERSECTION_ENABLE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/mst/jaccard", "Enable jaccard distance measure.", &detectordata->jaccardenable, FALSE, DEFAULT_JACCARD_ENABLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/mst/cosine", "Enable cosine distance measure.", &detectordata->cosineenable, FALSE, DEFAULT_COSINE_ENABLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/mst/simpson", "Enable simpson distance measure.", &detectordata->simpsonenable, FALSE, DEFAULT_SIMPSON_ENABLE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/mst/postprocenable", "Enable post-processing step..", &detectordata->postprocenable, FALSE, DEFAULT_POSTPROC_ENABLE, NULL, NULL ) );

#endif
   return SCIP_OKAY;
}
