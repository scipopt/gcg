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

/**@file   dec_mcl.cpp
 * @ingroup DETECTORS
 * @brief  detector MCL
 * @author Igor Pesic
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <time.h>     // for measuring time performance
#include <string>
#include <iostream>
#include "dec_mcl.h"
#include "cons_decomp.h"
#include "graph/rowgraph_weighted.h"
#include "graph/graph_gcg.h"

using gcg::RowGraphWeighted;
using gcg::Weights;
using gcg::GraphGCG;


/* constraint handler properties */
#define DEC_DETECTORNAME         "mcl"                               /**< name of detector */
#define DEC_DESC                 "detector based on mcl clustering"  /**< description of detector*/
#define DEC_PRIORITY             910         /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'L'         /**< display character of detector */
#ifdef GSL
#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#else
#define DEC_ENABLED              FALSE        /**< should the detection be enabled */
#endif
#define DEC_SKIP                 FALSE       /**< should detector be skipped if other detectors found decompositions */

/* Default parameter settings*/
#define DEFAULT_N_ITERATIONS              13
#define DEFAULT_JOHNSON_ENABLE            true
#define DEFAULT_INTERSECTION_ENABLE       false
#define DEFAULT_JACCARD_ENABLE            false
#define DEFAULT_COSINE_ENABLE             false
#define DEFAULT_SIMPSON_ENABLE            false
#define DEFAULT_POSTPROC_ENABLE           true
#define MAX_N_ITERATIONS                  20
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
   int n_iterations;                                     /**< number of different inflate factors to use for MCL algorithm */
   int n_similarities;                                  /**< number of active similarities */
   SCIP_Bool johnsonenable;                           /**< use Johnson similarity measure */
   SCIP_Bool intersectionenable;                      /**< use Intersection similarity measure */
   SCIP_Bool jaccardenable;                           /**< use Jaccard similarity measure */
   SCIP_Bool cosineenable;                            /**< use Cosine similarity measure */
   SCIP_Bool simpsonenable;                           /**< use Simpson similarity measure */
   SCIP_Bool postprocenable;                          /**< enable post processing step after clustering */
};


/*
 * Local methods
 */


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeMCL)
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
DEC_DECL_EXITDETECTOR(exitMCL)
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
DEC_DECL_INITDETECTOR(initMCL)
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

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectMCL)
{ /*lint --e{715}*/

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   *result = SCIP_DIDNOTFIND;

   detectordata->n_iterations = std::min(detectordata->n_iterations, MAX_N_ITERATIONS);

   Weights w(1, 1, 1, 1, 1, 1);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting MCL structure:");

   time_t start, cp1, d_s, d_e;
   time(&start);

   std::vector<std::string> sim;

   if(detectordata->johnsonenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::JOHNSON, gcg::WEIGHT_TYPE::SIM));
      detectordata->graphs->push_back(g);
      sim.push_back("Johnson");
   }
   if(detectordata->intersectionenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::INTERSECTION, gcg::WEIGHT_TYPE::SIM));
      detectordata->graphs->push_back(g);
      sim.push_back("Intersection");
   }
   if(detectordata->jaccardenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::JACCARD, gcg::WEIGHT_TYPE::SIM));
      detectordata->graphs->push_back(g);
      sim.push_back("Jaccard");
   }
   if(detectordata->cosineenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::COSINE, gcg::WEIGHT_TYPE::SIM));
      detectordata->graphs->push_back(g);
      sim.push_back("Cosine");
   }
   if(detectordata->simpsonenable)
   {
      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
      SCIP_CALL( g->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip),
            gcg::DISTANCE_MEASURE::SIMPSON, gcg::WEIGHT_TYPE::SIM));
      detectordata->graphs->push_back(g);
      sim.push_back("Simspon");
   }
   detectordata->n_similarities = (int) detectordata->graphs->size();


   std::vector<double> inflatefactors(detectordata->n_iterations);
   double inflate = 1.1;
   for(int i = 0; i < detectordata->n_iterations; i++)
   {
      inflatefactors[i] = inflate;
      inflate += 0.05;
   }
   time(&cp1);

   int max_ndecs = detectordata->n_iterations * detectordata->graphs->size();
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, max_ndecs) );


   const int max_blocks = std::min((int)round(0.3 * SCIPgetNConss(scip)), MAX_N_BLOCKS);
   int n_decs_found = 0;

   *ndecdecomps = 0;
   time(&d_s);
   for(int i = 0; i < (int)detectordata->graphs->size(); i++)
   {
      RowGraphWeighted<GraphGCG>* graph = detectordata->graphs->at(i);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n    %s similarity:", sim[i].c_str());
      int old_n_blocks = -1;
      int old_non_cl = -1;
      for(int j = 0; j < (int)inflatefactors.size() ; j++ )
      {
         double inflatefactor = inflatefactors[j];

         // run MCL with different inflate factors
         int stoppedAfter;
         SCIP_CALL( graph->computePartitionMCL(stoppedAfter, inflatefactor, detectordata->postprocenable) );

         int n_blocks;
         SCIP_CALL( graph->getNBlocks(n_blocks) );
         int non_cl;
         SCIP_CALL( graph->nonClustered(non_cl) );


         // skip the case if we have 1 block (it means we must increase inflate factor) or if the clustering is the same as the last one
         if(n_blocks == 1 || (n_blocks == old_n_blocks && non_cl == old_non_cl) )
         {
            continue;
         }
         // stop. inflate factor is already too big
         if( n_blocks > max_blocks )
         {
            break;
         }

         old_n_blocks = n_blocks;
         old_non_cl = non_cl;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n      Inflate factor: %.2f,    ", inflatefactor);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " Stopped after: %d iters,    ", stoppedAfter);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " Blocks: %d, Master Conss: %d/%d, ", n_blocks, non_cl, SCIPgetNConss(scip));

         SCIP_CALL( graph->createDecompFromPartition(&(*decdecomps)[n_decs_found]) );

         auto check = DECdecompGetNLinkingvars((*decdecomps)[n_decs_found]);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Link Vars: %d. ", check);

         if( (*decdecomps)[n_decs_found] != NULL )
         {
            *ndecdecomps += 1;
            ++n_decs_found;
            detectordata->found = TRUE;
         }
      }
      delete detectordata->graphs->at(i);
      detectordata->graphs->at(i) = NULL;
   }

   detectordata->graphs->clear();
   time(&d_e);
   double elapsed_graphs = difftime(cp1, start);
   double elapsed_mcl = difftime(d_e, d_s);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d similarities used, %d decompositions found.\n", detectordata->n_similarities, *ndecdecomps);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "MCL runtime: graphs: %.2lf, mcl: %.2lf. \n", elapsed_graphs, elapsed_mcl);

   SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );

   *result = *ndecdecomps > 0 ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
   if( *ndecdecomps == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, decdecomps);
   }
   return SCIP_OKAY;
}

#define propagateSeeedMCL NULL
//static
//DEC_DECL_PROPAGATESEEED(propagateSeeedMCL)
//{ /*lint --e{715}*/
//
//   int nNewSeeeds;
//   gcg::Seeed* seeed;
//   gcg::Seeed** newSeeeds;
//   DEC_DETECTORDATA* detectordata = DECdetectorGetData(detector);
//
//   seeedPropagationData->seeedToPropagate->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
//   assert(scip != NULL);
//   assert(detectordata != NULL);
//   *result = SCIP_DIDNOTFIND;
//
//   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);
//   seeed->assignAllDependent(seeedPropagationData->seeedpool);
//   if(!graphCompletible(seeedPropagationData->seeedpool, seeed))
//   {
//      delete seeed;
//      seeedPropagationData->nNewSeeeds = 0;
//      *result = SCIP_SUCCESS;
//      return SCIP_OKAY;
//   }
//
//   detectordata->n_iterations = std::min(detectordata->n_iterations, MAX_N_ITERATIONS);
//
//   Weights w(1, 1, 1, 1, 1, 1);
//
//   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting MCL structure:");
//
//   time_t start, cp1, d_s, d_e;
//   time(&start);
//
//   std::vector<std::string> sim;
//
//   if(detectordata->johnsonenable)
//   {
//      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
//      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::JOHNSON, gcg::WEIGHT_TYPE::SIM));
//      detectordata->graphs->push_back(g);
//      sim.push_back("Johnson");
//   }
//   if(detectordata->intersectionenable)
//   {
//      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
//      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::INTERSECTION, gcg::WEIGHT_TYPE::SIM));
//      detectordata->graphs->push_back(g);
//      sim.push_back("Intersection");
//   }
//   if(detectordata->jaccardenable)
//   {
//      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
//      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::JACCARD, gcg::WEIGHT_TYPE::SIM));
//      detectordata->graphs->push_back(g);
//      sim.push_back("Jaccard");
//   }
//   if(detectordata->cosineenable)
//   {
//      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
//      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::COSINE, gcg::WEIGHT_TYPE::SIM));
//      detectordata->graphs->push_back(g);
//      sim.push_back("Cosine");
//   }
//   if(detectordata->simpsonenable)
//   {
//      RowGraphWeighted<GraphGCG>* g = new RowGraphWeighted<GraphGCG>(scip, w);
//      SCIP_CALL( g->createFromPartialMatrix(seeedPropagationData->seeedpool, seeed, gcg::DISTANCE_MEASURE::SIMPSON, gcg::WEIGHT_TYPE::SIM));
//      detectordata->graphs->push_back(g);
//      sim.push_back("Simspon");
//   }
//   detectordata->n_similarities = (int) detectordata->graphs->size();
//
//
//   std::vector<double> inflatefactors(detectordata->n_iterations);
//   double inflate = 1.1;
//   for(int i = 0; i < detectordata->n_iterations; i++)
//   {
//      inflatefactors[i] = inflate;
//      inflate += 0.05;
//   }
//   time(&cp1);
//
//   int nMaxSeeeds = detectordata->n_iterations * detectordata->graphs->size();
//   SCIP_CALL( SCIPallocBufferArray(scip, &(newSeeeds), 2 * nMaxSeeeds) );
//
//
//   const int max_blocks = std::min((int)round(0.3 * SCIPgetNConss(scip)), MAX_N_BLOCKS);
//   int n_seeeds_found = 0;
//
//   nNewSeeeds = 0;
//   time(&d_s);
//   for(int i = 0; i < (int)detectordata->graphs->size(); i++)
//   {
//      RowGraphWeighted<GraphGCG>* graph = detectordata->graphs->at(i);
//      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n    %s similarity:", sim[i].c_str());
//      int old_n_blocks = -1;
//      int old_non_cl = -1;
//      for(int j = 0; j < (int)inflatefactors.size() ; j++ )
//      {
//         double inflatefactor = inflatefactors[j];
//
//         // run MCL with different inflate factors
//         int stoppedAfter;
//         SCIP_CALL( graph->computePartitionMCLForPartialGraph(seeedPropagationData->seeedpool, seeed, stoppedAfter, inflatefactor, detectordata->postprocenable) );
//
//         int n_blocks;
//         SCIP_CALL( graph->getNBlocks(n_blocks) );
//         int non_cl;
//         SCIP_CALL( graph->nonClustered(non_cl) );
//
//
//         // skip the case if we have 1 block (it means we must increase inflate factor) or if the clustering is the same as the last one
//         if(n_blocks == 1 || (n_blocks == old_n_blocks && non_cl == old_non_cl) )
//         {
//            continue;
//         }
//         // stop. inflate factor is already too big
//         if( n_blocks > max_blocks )
//         {
//            break;
//         }
//
//         old_n_blocks = n_blocks;
//         old_non_cl = non_cl;
//         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n      Inflate factor: %.2f,    ", inflatefactor);
//         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " Stopped after: %d iters,    ", stoppedAfter);
//         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " Blocks: %d, Master Conss: %d/%d, ", n_blocks, non_cl, SCIPgetNConss(scip));
//
//         SCIP_CALL( graph->createSeeedFromPartition(seeed,&newSeeeds[n_seeeds_found], &newSeeeds[n_seeeds_found+1], seeedPropagationData->seeedpool));
//
//         if((newSeeeds)[n_seeeds_found] != NULL)
//         {
//            nNewSeeeds += 2;
//            detectordata->found = TRUE;
//         }
//         n_seeeds_found += 2;
//      }
//      delete detectordata->graphs->at(i);
//      detectordata->graphs->at(i) = NULL;
//   }
//
//   delete seeed;
//   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), nNewSeeeds) );
//   seeedPropagationData->nNewSeeeds = nNewSeeeds;
//   for(int j = 0, s = 0; s < n_seeeds_found; ++j)
//   {
//      if(newSeeeds[j] != NULL)
//      {
//         seeedPropagationData->newSeeeds[s] = newSeeeds[j];
//         seeedPropagationData->newSeeeds[s]->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));
//         ++s;
//      }
//   }
//   SCIPfreeBufferArray(scip, &newSeeeds);
//
//   detectordata->graphs->clear();
//   time(&d_e);
//   double elapsed_graphs = difftime(cp1, start);
//   double elapsed_mcl = difftime(d_e, d_s);
//
//   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d similarities used, %d decompositions found.\n", detectordata->n_similarities, nNewSeeeds);
//   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "MCL runtime: graphs: %.2lf, mcl: %.2lf. \n", elapsed_graphs, elapsed_mcl);
//
//   *result = nNewSeeeds > 0 ? SCIP_SUCCESS: SCIP_DIDNOTFIND;
//   if( nNewSeeeds == 0 )
//   {
//      SCIPfreeMemoryArrayNull(scip, &(seeedPropagationData->newSeeeds));
//   }
//   return SCIP_OKAY;
//}

/*
 * detector specific interface methods
 */

/** creates the handler for xyz detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorMCL(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#if !defined(_WIN32) && !defined(_WIN64)
   DEC_DETECTORDATA *detectordata = NULL;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP,
      detectordata, detectMCL, freeMCL, initMCL, exitMCL, propagateSeeedMCL) );

   /* add arrowheur presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/mcl/niterations", "Number of iterations to run MCL with different inflate factor (max=20).", &detectordata->n_iterations, FALSE, DEFAULT_N_ITERATIONS, 1, 20, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/mcl/johson", "Enable johson distance measure.", &detectordata->johnsonenable, FALSE, DEFAULT_JOHNSON_ENABLE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/mcl/intersection", "Enable intersection distance measure.", &detectordata->intersectionenable, FALSE, DEFAULT_INTERSECTION_ENABLE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/mcl/jaccard", "Enable jaccard distance measure.", &detectordata->jaccardenable, FALSE, DEFAULT_JACCARD_ENABLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/mcl/cosine", "Enable cosine distance measure.", &detectordata->cosineenable, FALSE, DEFAULT_COSINE_ENABLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/mcl/simpson", "Enable simpson distance measure.", &detectordata->simpsonenable, FALSE, DEFAULT_SIMPSON_ENABLE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/mcl/postprocenable", "Enable post-processing step..", &detectordata->postprocenable, FALSE, DEFAULT_POSTPROC_ENABLE, NULL, NULL ) );

#endif
   return SCIP_OKAY;
}
