/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file    bliss_automorph.cpp
 * @brief   automorphism recognition of SCIPs
 * @author  Daniel Peters
 * @author  Martin Bergner
 * @author  Jonas Witt
 * @author  Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */
#include "symmetry/automorphism.h"
#include "symmetry/pub_automorphism.h"
#include "gcg/scip_misc.h"
#include "scip/scip.h"
#include "gcg/gcg.h"
#include "scip/cons_linear.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/cons_decomp.hpp"

#include <cstring>

/** saves information of the permutation */
struct AUT_HOOK2
{
   SCIP_Bool aut;                            /**< true if there is an automorphism */
   unsigned int n;                           /**< number of permutations */
   SCIP_HASHMAP* varmap;                     /**< hashmap for permutated variables */
   SCIP_HASHMAP* consmap;                    /**< hashmap for permutated constraints */
   SCIP** scips;                             /**< array of scips to search for automorphisms */
   int* nodemap;                             /**< mapping of the nodes; filled generator-wise */
   int* conssperm;                           /**< mapping of constraints */
   gcg::DETPROBDATA* detprobdata;            /**< problem information the automorphism should be searched for */
   gcg::PARTIALDECOMP* partialdec;           /**< decomposition information */
   std::vector<int>* blocks;                 /**< array of blocks the automporphisms are searched for */
   SCIP* scip;
   int ncalls;
   int generatorlimit;
   AUT_GRAPH* graph;


   /** constructor for the hook struct*/
   AUT_HOOK2(
      SCIP_HASHMAP* varmap,                  /**< hashmap for permutated variables */
      SCIP_HASHMAP* consmap,                 /**< hashmap for permutated constraints */
      unsigned int n,                        /**< number of permutations */
      AUT_GRAPH* graph,                      /**< graph used to search for automorphisms */
      SCIP* scip,                            /**< SCIP data structure */
      gcg::DETPROBDATA* givendetprobdata,
      gcg::PARTIALDECOMP* givenpartialdec,
      std::vector<int>* givenblocks
      );

   /** destructor for hook struct */
   ~AUT_HOOK2();


   /** getter for the bool aut */
   SCIP_Bool getBool();

   /** setter for the bool aut */
   void setBool(SCIP_Bool aut);

   /** getter for the number of nodes */
   unsigned int getNNodes();

   /** getter for the variables hashmap */
   SCIP_HASHMAP* getVarHash();

   /** getter for the constraints hashmap */
   SCIP_HASHMAP* getConsHash();

   /** getter for the graph */
   AUT_GRAPH* getGraph();


};


void AUT_HOOK2::setBool( SCIP_Bool aut_ )
{
   aut = aut_;
}


AUT_HOOK2::~AUT_HOOK2()
{   /*lint -esym(1540,struct_hook::conssperm) */
   SCIPfreeMemoryArrayNull(scip, &nodemap);
   if( conssperm != NULL )
      SCIPfreeMemoryArrayNull(scip, &conssperm);
   conssperm  = NULL;
   detprobdata = NULL;
   partialdec = NULL;
}



SCIP_Bool AUT_HOOK2::getBool()
{
   return aut;
}

unsigned int AUT_HOOK2::getNNodes()
{
   return n;
}

SCIP_HASHMAP* AUT_HOOK2::getVarHash()
{
   return varmap;
}

SCIP_HASHMAP* AUT_HOOK2::getConsHash()
{
   return consmap;
}

AUT_GRAPH* AUT_HOOK2::getGraph()
{
   return graph;
}

/** constructor of the hook struct */
AUT_HOOK2::AUT_HOOK2(
   SCIP_HASHMAP*         varmap_,            /**< hashmap of permutated variables */
   SCIP_HASHMAP*         consmap_,           /**< hahsmap of permutated constraints */
   unsigned int          n_,                 /**< number of permutations */
   AUT_GRAPH*            graph_,             /**< graph used to search for automorphisms */
   SCIP*                 scip_,              /**< SCIP data structure */
   gcg::DETPROBDATA*     detprobdata_,
   gcg::PARTIALDECOMP*   partialdec_,
   std::vector<int>*     blocks_
   )
{
   size_t i;
   scip = scip_;
   aut = FALSE;
   n = n_;
   consmap = consmap_;
   varmap = varmap_;
   graph = graph_;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &nodemap, n_) ); /*lint !e666*/
   for (i = 0; i < n_; ++i)
      nodemap[i] = -1;

   conssperm = NULL;
   detprobdata = NULL;
   partialdec = NULL;
   blocks = NULL;

   ncalls = 0;
   generatorlimit = 0;

   detprobdata = detprobdata_;
   partialdec = partialdec_;
   blocks = blocks_;

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &(this->conssperm), detprobdata->getNConss() ) ); /*lint !e666*/
}

/** hook function to save the permutation of the graph; fhook() is called by metis for every generator,
 *  AUT_HOOK* hook  stores a combined mapping in nodemapping that is filled generator-wise */
static
void fhook(
   void*                 user_param,         /**< data structure to save hashmaps with permutation */
   unsigned int          N,                  /**< number of permutations */
   const unsigned int*   aut                 /**< array of permutations */
   )
{ /*lint -e715*/

   unsigned int i;
   unsigned int j;
   unsigned int n;
   int nvars;
   int nconss;
   SCIP_VAR** vars1 = NULL;
   SCIP_VAR** vars2 = NULL;
   SCIP_CONS** conss1 = NULL;
   SCIP_CONS** conss2 = NULL;
   SCIP* partialdecscip = NULL;
   gcg::DETPROBDATA* detprobdata = NULL;
   gcg::PARTIALDECOMP* partialdec = NULL;
   AUT_HOOK2* hook = (AUT_HOOK2*) user_param;

   j = 0;
   n = hook->getNNodes();

   /* new detection stuff */
   partialdec = hook->partialdec;
   detprobdata = hook->detprobdata;
   partialdecscip = NULL;

   ++hook->ncalls;

   if( hook->getBool() )
      return;

   // fallback check if library does not support termination
   if( hook->generatorlimit > 0 && hook->ncalls > hook->generatorlimit )
   {
      hook->setBool(false);
      return;
   }
   
  // SCIPdebugMessage("Looking for a permutation from [0,%u] bijective to [%u:%u] (N=%u) \n", n/2-1, n/2, n-1, N);
   for( i = 0; i < n / 2; i++ )
   {
      assert(aut[i] < INT_MAX);

      if( (aut[i]) >= n / 2 && hook->nodemap[i] == -1 )
      {
         assert(aut[i] < n);
//         SCIPdebugMessage("current generator: %u -> %u\n", i, aut[i]);
         hook->nodemap[i] = aut[i];
      }
   }

   for( i = 0; i < n / 2; i++ )
   {
//      SCIPdebugMessage("general mapping : %u -> %u\n", i, hook->nodemap[i]);
      if( hook->nodemap[i] >= (int) n / 2 )
         ++j;
   }

   if( j == n / 2 )
   {
      hook->setBool(TRUE);
   }

   for( i = n; i < N; ++i )
   {
      if( aut[i] != i )
      {
    //     SCIPdebugMessage("Master %u -> %u not the identity, no decomposition possible!\n", i, aut[i]);
         hook->setBool(false);
         break;
      }
   }

//   SCIPdebugMessage("Permutation %s found.\n", hook->getBool() ? "":"not");
//   SCIPdebugMessage("j = %u\n", j);

   if( !hook->getBool() )
      return;


   nvars = partialdec->getNVarsForBlock((*hook->blocks)[0]);

   assert(nvars == partialdec->getNVarsForBlock((*hook->blocks)[1]) );

   partialdecscip = detprobdata->getScip();

   SCIP_CALL_ABORT(SCIPallocBufferArray(partialdecscip, &vars1, nvars ));
   SCIP_CALL_ABORT(SCIPallocBufferArray(partialdecscip, &vars2, nvars ));
   nconss = partialdec->getNConssForBlock((*hook->blocks)[0]);
   assert(nconss == partialdec->getNConssForBlock((*hook->blocks)[1]));

   SCIP_CALL_ABORT( SCIPallocBufferArray(partialdecscip, &conss1, nconss ) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(partialdecscip, &conss2, nconss ) );

   for( int v = 0; v < nvars; ++v )
   {
      vars1[v] = detprobdata->getVar(partialdec->getVarsForBlock((*hook->blocks)[0])[v]);
      vars2[v] = detprobdata->getVar(partialdec->getVarsForBlock((*hook->blocks)[1])[v]);
   }

   for( int c = 0; c < nconss; ++c )
   {
      conss1[c] = detprobdata->getCons(partialdec->getConssForBlock((*hook->blocks)[0])[c]);
      conss2[c] = detprobdata->getCons(partialdec->getConssForBlock((*hook->blocks)[1])[c]);
   }

   for( i = 0; i < (unsigned int) nvars+nconss; i++ )
   {
      /* Assuming the following layout:
       *  0 ... nconss-1 = vertex ids for constraints
       *  nconss ... nconss+nvars-1 = vertex ids for variables
       *  nconss+nvars ... n-1 = nonzero entries (not relevant)
       */
      if( i < (unsigned int) nconss )
      {
         unsigned int consindex = i;
         unsigned int consindex2 = hook->nodemap[i]-n/2;
         assert( consindex2 < (unsigned int) nconss);
         SCIP_CONS* cons1 = conss1[consindex];
         SCIP_CONS* cons2 = conss2[consindex2];
         SCIP_CALL_ABORT( SCIPhashmapInsert(hook->getConsHash(), cons2, cons1) );
         SCIPdebugMessage("cons <%s> <-> cons <%s>\n", SCIPconsGetName(cons2), SCIPconsGetName(cons1));
      }
      else if( i < (unsigned int) nvars+nconss )
      {
         unsigned int varindex = i-nconss;
         unsigned int varindex2 = hook->nodemap[i]-nconss-n/2;
         assert( varindex2 < (unsigned int) nvars);
         SCIP_VAR* var1 = vars1[varindex];
         SCIP_VAR* var2 = vars2[varindex2];
         SCIP_CALL_ABORT( SCIPhashmapInsert(hook->getVarHash(), var2, var1) );
         SCIPdebugMessage("var <%s> <-> var <%s>\n", SCIPvarGetName(var2), SCIPvarGetName(var1));
      }
   }

   partialdecscip = detprobdata->getScip();
   SCIPfreeBufferArray(partialdecscip, &conss1);
   SCIPfreeBufferArray(partialdecscip, &conss2);
   SCIPfreeBufferArray(partialdecscip, &vars1);
   SCIPfreeBufferArray(partialdecscip, &vars2);

   hook->getGraph()->terminateSearch();
}


/** constructor for colorinfo arrays */
static SCIP_RETCODE allocMemoryNewDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   gcg::DETPROBDATA*     detprobdata,        /**< SCIP data structure */
   AUT_COLOR*            colorinfo,          /**< struct to save intermediate information */
   int                   nconss,             /**< number of constraints */
   int                   nvars,              /**< number of variables */
   int                   ncoeffs             /**< number of coefficients */
   )
{
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarraycoefs, ((size_t) ncoeffs )));
   colorinfo->alloccoefsarray = ncoeffs;
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarrayvars, (size_t) nvars));
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarrayconss, (size_t) nconss));
   return SCIP_OKAY;
}


/** destructor for colorinfoarrays */
static
SCIP_RETCODE freeMemory(
   SCIP*                 scip,               /**< SCIP data structure */
   AUT_COLOR*            colorinfo           /**< struct to save intermediate information */
   )
{
   int i;

   for(i = 0; i < colorinfo->lenvarsarray; i++ ){
      AUT_VAR* svar =  (AUT_VAR*) colorinfo->ptrarrayvars[i];
      delete svar;
   }
   for(i = 0; i < colorinfo->lenconssarray; i++ ){
      AUT_CONS* scons = (AUT_CONS*) colorinfo->ptrarrayconss[i];
      delete scons;
   }
   for(i = 0; i < colorinfo->lencoefsarray; i++ ){
      AUT_COEF* scoef = (AUT_COEF*) colorinfo->ptrarraycoefs[i];
      delete scoef;
   }

   SCIPfreeMemoryArray(scip, &colorinfo->ptrarraycoefs);
   SCIPfreeMemoryArray(scip, &colorinfo->ptrarrayconss);
   SCIPfreeMemoryArray(scip, &colorinfo->ptrarrayvars);
   return SCIP_OKAY;
}


/** set up a help structure for graph creation for new detection loop*/
static
SCIP_RETCODE setuparraysnewdetection(
   SCIP*                 scip,               /**< SCIP data structure */
   gcg::DETPROBDATA*     detprobdata,        /**< detprobdata corresponing to presolved or unpresolved problem */
   gcg::PARTIALDECOMP*   partialdec,         /**< partial decomp the  symmetry for two blocks is checked for */
   int                   nblocks,            /**< number of blocks the symmetry should be checked for */
   std::vector<int>      blocks,             /**< vectors of block indices the symmetry be checked for */
   AUT_COLOR*            colorinfo,          /**< data structure to save intermediate data */
   SCIP_RESULT*          result              /**< result pointer to indicate success or failure */
   )
{ /*lint -esym(593, scoef) */
   int i;
   int j;
   int b;
   int nconss;
   int nvars;
   SCIP_Bool added;

   added = FALSE;

   //allocate max n of coefarray, varsarray, and boundsarray in origscip
   nconss = partialdec->getNConssForBlock(blocks[0]) ;
   nvars = partialdec->getNVarsForBlock(blocks[0]) ;
   colorinfo->setOnlySign(FALSE);

   for( b = 0; b < nblocks && *result == SCIP_SUCCESS; ++b )
   {
      int block = blocks[b];

      assert( partialdec->getNVarsForBlock(blocks[b]) == nvars );
      assert( partialdec->getNConssForBlock(blocks[b]) == nconss );

      SCIPdebugMessage("Handling block %i (id %d %d x %d)\n", b, block, partialdec->getNConssForBlock(blocks[b]), partialdec->getNVarsForBlock(blocks[b]));
      //save the properties of variables in a struct array and in a sorted pointer array
      for( i = 0; i < nvars; i++ )
      {
         SCIP_VAR* var;
         AUT_VAR* svar;

         var = detprobdata->getVar(partialdec->getVarsForBlock(block)[i]);
         svar = new AUT_VAR(scip, var);
         //add to pointer array iff it doesn't exist
         SCIP_CALL( colorinfo->insert(svar, &added) );
         if( b > 0 && added)
         {
           *result = SCIP_DIDNOTFIND;
            break;
         }
         //otherwise free allocated memory
         if( !added )
            delete svar;
      }
      //save the properties of constraints in a struct array and in a sorted pointer array
      for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
      {
         int consid;
         SCIP_CONS* cons;

         consid = partialdec->getConssForBlock(block)[i];
         cons = detprobdata->getCons(consid);

         if( detprobdata->getNVarsForCons(consid) == 0 )
            continue;

         AUT_CONS* scons = new AUT_CONS(scip, cons);
         //add to pointer array iff it doesn't exist
         SCIP_CALL( colorinfo->insert(scons, &added) );
         if( b > 0 && added)
         {
           *result = SCIP_DIDNOTFIND;
           break;
         }
         //otherwise free allocated memory
         if( !added )
            delete scons;

         //save the properties of variables of the constraints in a struct array and in a sorted pointer array
         for( j = 0; j < detprobdata->getNVarsForCons(consid); j++ )
         {
            SCIP_Real val;
            AUT_COEF* scoef;
            val = detprobdata->getVal(consid, detprobdata->getVarsForCons(consid)[j]);
            scoef = new AUT_COEF(scip, val );
            //test, whether the coefficient is not zero
            if( !SCIPisZero(scip, scoef->getVal()) )
            {
               //add to pointer array iff it doesn't exist
               SCIP_CALL( colorinfo->insert(scoef, &added) );
               if( b > 0 && added)
               {
                  *result = SCIP_DIDNOTFIND;
                  break;
               }
            }
            //otherwise free allocated memory
            if( !added )
               delete scoef;
         }
      }
   }

   /* add color information for master constraints */
   for( i = 0; i < partialdec->getNMasterconss() && *result == SCIP_SUCCESS; ++i )
   {
      int masterconsid;
      SCIP_CONS* mastercons;

      masterconsid = partialdec->getMasterconss()[i];
      mastercons = detprobdata->getCons(masterconsid);

      /* add right color for master constraint */
      AUT_CONS* scons = new AUT_CONS(scip, mastercons);
      SCIP_CALL( colorinfo->insert(scons, &added) );

      /* if it hasn't been added, it is already present */
      if(!added)
         delete scons;

      for( j = 0; j < detprobdata->getNVarsForCons(masterconsid); ++j )
      {
         AUT_COEF* scoef;
         int varid;

         varid = detprobdata->getVarsForCons(masterconsid)[j];
         scoef= new AUT_COEF(scip, detprobdata->getVal(masterconsid, varid) );

         added = FALSE;

         if( !SCIPisZero(scip, scoef->getVal()) )
         {
            SCIP_CALL( colorinfo->insert(scoef, &added) );
         }

         if( !added )
            delete scoef;
      }
   }

   return SCIP_OKAY;
}

/** create a graph out of an array of scips */
static
SCIP_RETCODE createGraphNewDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   gcg::DETPROBDATA*     detprobdata,        /**< detection process information and data */
   gcg::PARTIALDECOMP*   partialdec,         /**< partial decomposition */
   std::vector<int>      blocks,             /**< vectors of block indices the symmetry be checked for */
   AUT_COLOR             colorinfo,          /**< data structure to save intermediate data  */
   AUT_GRAPH*            graph,              /**< graph needed for discovering isomorphism */
   int*                  pricingnodes,       /**< number of pricing nodes without master  */
   SCIP_RESULT*          result              /**< result pointer to indicate success or failure */
   )
{
   int i;
   int j;
   int b;
   int ncurvars;
   int* nnodesoffset;
   int color;
   int nconss;
   int nvars;
   int currentnode;
   int* pricingnonzeros;
   int* mastercoefindex;
   unsigned int nconsvarpairs;
   unsigned int nmasterconsnzs;
   unsigned int nnodes;
   unsigned int nnonemptyconss;
   int nblocks = (int)blocks.size();
   std::vector<bool> masterconssrelevant(partialdec->getNMasterconss(), false);

   if( *result != SCIP_SUCCESS )
      return SCIP_OKAY;

   pricingnonzeros = NULL;
   mastercoefindex = NULL;
   nnodesoffset = NULL;

   currentnode = 0;

   scip = detprobdata->getScip();

   SCIP_CALL( SCIPallocBufferArray(scip, &mastercoefindex, nblocks) );
   BMSclearMemoryArray(mastercoefindex, nblocks);

   nconss = partialdec->getNConssForBlock(blocks[0]);
   nvars = partialdec->getNVarsForBlock(blocks[0]);

   SCIP_CALL( SCIPallocBufferArray(scip, &nnodesoffset, nblocks) );
   BMSclearMemoryArray(nnodesoffset, nblocks);
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingnonzeros, nblocks) );
   BMSclearMemoryArray(pricingnonzeros, nblocks);

   nconsvarpairs = 0;
   nnonemptyconss = 0;
   for( i = 0; i < nconss; i++ )
   {
      ncurvars = detprobdata->getNVarsForCons(partialdec->getConssForBlock(blocks[0])[i]);
      if( ncurvars > 0 )
      {
         nconsvarpairs += ncurvars;
         nnonemptyconss++;
      }
   }

   nmasterconsnzs = 0;
   for( i = 0; i < partialdec->getNMasterconss(); ++i )
      nmasterconsnzs += partialdec->getNVarsOfBlockInMasterCons(i, blocks[0]);

   nnodes = nblocks * (nnonemptyconss + nvars + nconsvarpairs + nmasterconsnzs) + partialdec->getNMasterconss();
   graph->init(scip, nnodes);

   for( b = 0; b < nblocks && *result == SCIP_SUCCESS; ++b )
   {
      int block = blocks[b];
      int z = 0;

      SCIPdebugMessage("Pricing problem %d\n", block);
  
      nnodesoffset[b] = currentnode;

      //add a node for every constraint
      for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
      {
         int consid;
         SCIP_CONS* cons;

         consid = partialdec->getConssForBlock(block)[i];
         ncurvars = detprobdata->getNVarsForCons(consid);
         cons = detprobdata->getCons(consid);

         if( ncurvars == 0 )
            continue;

         color = colorinfo.get( AUT_CONS(scip, cons) );

         if(color == -1)
         {
            *result = SCIP_DIDNOTFIND;
            break;
         }

         SCIPdebugMessage("cons <%s> color %d\n", SCIPconsGetName(cons), color);
         graph->setColor(currentnode, color);
         currentnode++;
      }
      //add a node for every variable
      for( i = 0; i < nvars && *result == SCIP_SUCCESS; i++ )
      {
         int varid;
         SCIP_VAR* var;

         varid = partialdec->getVarsForBlock(block)[i];
         var = detprobdata->getVar(varid);

         color = colorinfo.get( AUT_VAR(scip, var) );

         if(color == -1) {
            *result = SCIP_DIDNOTFIND;
            break;
         }

         SCIPdebugMessage("var <%s> color %d\n", SCIPvarGetName(var), color);
         graph->setColor(currentnode, colorinfo.getLenCons() + color);
         currentnode++;
      }
      //connecting the nodes with an additional node in the middle
      //it is necessary, since only nodes have colors
      for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
      {
         int consid;
         SCIP_CONS* cons;
         int conscolor;

         consid = partialdec->getConssForBlock(block)[i];
         ncurvars = detprobdata->getNVarsForCons(consid);
         cons = detprobdata->getCons(consid);
         conscolor = colorinfo.get(AUT_CONS(scip, cons));

         if( ncurvars == 0 )
            continue;

         for( j = 0; j < ncurvars; j++ )
         {
            int varcolor;
            int varid;
            SCIP_VAR* var;
            SCIP_Real val;

            varid = detprobdata->getVarsForCons(consid)[j];
            var = detprobdata->getVar(varid);

            val = detprobdata->getVal(consid, varid);

            varcolor = colorinfo.get( AUT_VAR(scip, var )) + colorinfo.getLenCons(); /*lint !e864 */
            color = colorinfo.get( AUT_COEF(scip, val ));
            if( color == -1 )
            {
               *result = SCIP_DIDNOTFIND;
               break;
            }
            color += colorinfo.getLenCons() + colorinfo.getLenVar(); /*lint !e864 */
            graph->setColor(currentnode, color);
            currentnode++;
            graph->addEdge(nnodesoffset[b] + i, nnodesoffset[b] + nconss + nvars + z);
            graph->addEdge(nnodesoffset[b] + nconss + nvars + z, nnodesoffset[b]+nconss + partialdec->getVarProbindexForBlock(varid, block)     );
            SCIPdebugMessage("nz: c <%s> (id: %d, color: %d) -> nz (id: %d) (value: %f, color: %d) -> var <%s> (id: %d, color: %d) \n",
                              SCIPconsGetName(cons),
                              nnodesoffset[b] + i,
                              conscolor,
                              nnodesoffset[b] + nconss + nvars + z,
                              val,
                              color+colorinfo.getLenCons() + colorinfo.getLenVar(), /*lint !e864 */
                              SCIPvarGetName(var),
                              nnodesoffset[b]+nconss + partialdec->getVarProbindexForBlock(varid, block),
                              varcolor);
            z++;
         }
      }
      pricingnonzeros[b] = z;

      /* add coefficient nodes for nonzeros in the master */
      for( i = 0; i < partialdec->getNMasterconss() && *result == SCIP_SUCCESS; ++i )
      {
         int masterconsid;

         masterconsid = partialdec->getMasterconss()[i];
         ncurvars = detprobdata->getNVarsForCons(masterconsid);

         for( j = 0; j < ncurvars; ++j )
         {
            int varid;
            SCIP_VAR* var;
            SCIP_Real val;

            varid = detprobdata->getVarsForCons(masterconsid)[j];
            /* ignore if the variable belongs to a different block */
            if( !partialdec->isVarBlockvarOfBlock(varid, block) )
            {
//               SCIPdebugMessage("Var <%s> belongs to a different block (%d)\n", SCIPvarGetName(detprobdata->getVar(varid) ), block);
               continue;
            }

            var = detprobdata->getVar(varid);
            val = detprobdata->getVal(masterconsid, varid);
            color = colorinfo.get(AUT_COEF(scip, val));
            assert(color != -1);
            color += colorinfo.getLenCons() + colorinfo.getLenVar(); /*lint !e864 */

            masterconssrelevant[i] = true;

            /* add coefficent node for current coeff */
            graph->setColor(currentnode, color);
            assert(ABS(val) < SCIPinfinity(scip));
            SCIPdebugMessage("master nz for var <%s> (id: %d) (value: %f, color: %d)\n", SCIPvarGetName(var), currentnode, val, color);
            currentnode++;
         }
      }
      SCIPdebugMessage("Iteration %d: currentnode = %d\n", b, currentnode);
   }
   /* connect the created graphs with nodes for the master problem */

   SCIPdebugMessage( "handling %d masterconss\n", partialdec->getNMasterconss());
   *pricingnodes = currentnode;

   for( i = 0; i < partialdec->getNMasterconss() && *result == SCIP_SUCCESS; ++i )
   {
      int masterconsid;
      SCIP_CONS* mastercons;
      int masterconsnode;
      int conscolor;

      /**experimental */
      if( !masterconssrelevant[i] )
         continue;
      /*experimental end */

      masterconsid= partialdec->getMasterconss()[i];
      mastercons = detprobdata->getCons(masterconsid);
      ncurvars = detprobdata->getNVarsForCons(masterconsid);

      SCIPdebugMessage("Handling cons <%s>\n", SCIPconsGetName(mastercons));

      /* create node for masterconss and get right color */
      conscolor = colorinfo.get(AUT_CONS(scip, mastercons) );
      assert(conscolor != -1);
      graph->setColor(currentnode, conscolor);
      masterconsnode = currentnode;
      currentnode++;

      for( j = 0; j < ncurvars; ++j )
      {
         int varid;
         SCIP_VAR* var;
         SCIP_Real val;
         int blockid;
         int coefnodeindex;
         int bid;
         int varcolor;

         blockid = -1;
         bid = -1;
         varid = detprobdata->getVarsForCons(masterconsid)[j];

         var = detprobdata->getVar(varid);

         for( b = 0; b < nblocks; ++b )
         {
            if( partialdec->isVarBlockvarOfBlock(varid, blocks[b]) )
            {
               bid = b;
               blockid = blocks[b];
               break;
            }
         }

         /* ignore if the variable belongs to a different block */
         if( blockid == -1 )
         {
            //SCIPdebugMessage("Var <%s> belongs to a different block \n", SCIPvarGetName(var));
            continue;
         }
         val = detprobdata->getVal(masterconsid, varid);

         color = colorinfo.get(AUT_COEF(scip, val));
         assert(color != -1);
         color += colorinfo.getLenCons() + colorinfo.getLenVar(); /*lint !e864 */

         /* get coefficient node for current coefficient */
         coefnodeindex = nnodesoffset[bid] + nvars + nconss + pricingnonzeros[bid] + mastercoefindex[bid];
         ++(mastercoefindex[bid]);

         varcolor = colorinfo.get(AUT_VAR(scip, var));
         assert(varcolor != -1);
         varcolor += colorinfo.getLenCons();

         assert( (unsigned int) masterconsnode < graph->getNVertices());
         assert( (unsigned int) coefnodeindex < graph->getNVertices());
         /* master constraint and coefficient */
         graph->addEdge(masterconsnode, coefnodeindex);
         SCIPdebugMessage("ma: c <%s> (id: %d, color: %d) -> nz (id: %d) (value: <%.6f> , color: %d) -> pricingvar <%s> (id: %d, color: %d)\n",
            SCIPconsGetName(mastercons),
            masterconsnode, conscolor, coefnodeindex, val, color, SCIPvarGetName(var),
            nnodesoffset[bid] + nconss + varid, varcolor);

         /* get node index for pricing variable and connect masterconss, coeff and pricingvar nodes */
         graph->addEdge(coefnodeindex, nnodesoffset[bid] + nconss + partialdec->getVarProbindexForBlock(varid, blockid));
      }
   }

#ifndef NDEBUG
   if( *result == SCIP_SUCCESS )
      assert(currentnode == nnodes);
#endif

   //free all allocated memory
   SCIPfreeBufferArray(scip, &pricingnonzeros);
   SCIPfreeBufferArray(scip, &nnodesoffset);
   SCIPfreeBufferArray(scip, &mastercoefindex);

   return SCIP_OKAY;
}

/** compare two graphs w.r.t. automorphism */
SCIP_RETCODE cmpGraphPair(
   SCIP*                   scip,               /**< SCIP data structure */
   gcg::PARTIALDECOMP*     partialdec,         /**< partialdec the graphs should be compared for */
   int                     block1,             /**< index of first pricing prob */
   int                     block2,             /**< index of second pricing prob */
   SCIP_RESULT*            result,             /**< result pointer to indicate success or failure */
   SCIP_HASHMAP*           varmap,             /**< hashmap to save permutation of variables */
   SCIP_HASHMAP*           consmap,            /**< hashmap to save permutation of constraints */
   unsigned int            searchnodelimit,    /**< bliss search node limit (requires patched bliss version) */
   unsigned int            generatorlimit      /**< bliss generator limit (requires patched bliss version or version >=0.76) */
   )
{
   AUT_GRAPH graph;
   AUT_HOOK2 *ptrhook;
   AUT_COLOR colorinfo;
   std::vector<int> blocks(2, -1);
   gcg::DETPROBDATA* detprobdata;
   int nconss;
   int nvars;
   int ncoeffs;

   int pricingnodes;

   *result = SCIP_SUCCESS;

   assert(partialdec != NULL );

   blocks[0] = block1;
   blocks[1] = block2;
   pricingnodes = 0;

   if (partialdec->isAssignedToOrigProb() )
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(scip);
   else
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(scip);

   assert(detprobdata != NULL);

   //allocate max n of coefarray, varsarray, and boundsarray in origscip
   nconss = partialdec->getNConssForBlock(blocks[0]) ;
   nvars = partialdec->getNVarsForBlock(blocks[0]) ;
   ncoeffs = partialdec->getNCoeffsForBlock( blocks[0]);
   SCIP_CALL( allocMemoryNewDetection(scip, detprobdata, &colorinfo, nconss*2+partialdec->getNMasterconss(), nvars*2, ncoeffs*2 + partialdec->getNCoeffsForMaster() ) );
   SCIP_CALL( setuparraysnewdetection(scip, detprobdata, partialdec, 2, blocks, &colorinfo, result) );
   SCIPdebugMessage("finished setup array method.\n");
   if( *result == SCIP_SUCCESS )
   {
      SCIP_CALL( createGraphNewDetection(scip, detprobdata, partialdec, blocks, colorinfo, &graph,  &pricingnodes, result) );
      SCIPdebugMessage("finished create graph.\n");
      ptrhook = new AUT_HOOK2(varmap, consmap, (unsigned int) pricingnodes, &graph, detprobdata->getScip(), detprobdata, partialdec, &blocks);
      ptrhook->generatorlimit = generatorlimit;
      SCIPdebugMessage("finished creating aut hook.\n");

      graph.findAutomorphisms(ptrhook, &fhook, searchnodelimit, generatorlimit);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL , NULL, "finished calling bliss: number of reporting function calls (=number of generators): %d \n", ptrhook->ncalls);
      if( !ptrhook->getBool() )
         *result = SCIP_DIDNOTFIND;
      graph.destroy();
      delete ptrhook;
   }
   SCIP_CALL( freeMemory(scip, &colorinfo) );

   SCIPdebugMessage("finished find automorphisms.\n");

   return SCIP_OKAY;
}

static
int getSign(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value */
   )
{
   if( SCIPisNegative(scip, val) )
      return -1;
   if( SCIPisPositive(scip, val) )
      return 1;
   else
      return 0;
}

/** compare two values of two scips */
static
int comp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< value 1 to compare */
   SCIP_Real             val2,               /**< value 2 to compare */
   SCIP_Bool             onlysign            /**< use sign of values instead of values? */
   )
{
   SCIP_Real compval1;
   SCIP_Real compval2;

   if( onlysign )
   {
      compval1 = getSign(scip, val1);
      compval2 = getSign(scip, val2);
   }
   else
   {
      compval1 = val1;
      compval2 = val2;
   }


   if( SCIPisLT(scip, compval1, compval2) )
      return -1;
   if( SCIPisGT(scip, compval1, compval2) )
      return 1;
   else
      return 0;
}

/** compare two constraints of two scips */
static
int comp(
   SCIP*                 scip,               /**< SCIP data structure */
   AUT_CONS*             cons1,              /**< constraint 1 to compare */
   AUT_CONS*             cons2,              /**< constraint 2 to compare */
   SCIP_Bool             onlysign            /**< use sign of values instead of values? */
   )
{
   if( comp(scip, GCGconsGetRhs(scip, cons1->getCons()), GCGconsGetRhs(scip, cons2->getCons()), onlysign) != 0 )
      return comp(scip, GCGconsGetRhs(scip, cons1->getCons()), GCGconsGetRhs(scip, cons2->getCons()), onlysign);
   assert(SCIPisEQ(scip, GCGconsGetRhs(scip, cons1->getCons()), GCGconsGetRhs(scip, cons2->getCons())) || onlysign);

   if( comp(scip, GCGconsGetLhs(scip, cons1->getCons()), GCGconsGetLhs(scip, cons2->getCons()), onlysign) != 0 )
      return comp(scip, GCGconsGetLhs(scip, cons1->getCons()), GCGconsGetLhs(scip, cons2->getCons()), onlysign);
   assert(SCIPisEQ(scip, GCGconsGetLhs(scip, cons1->getCons()), GCGconsGetLhs(scip, cons2->getCons())) || onlysign);

   if( comp(scip, GCGconsGetNVars(scip, cons1->getCons()), GCGconsGetNVars(scip, cons2->getCons()), FALSE) != 0 )
      return comp(scip, GCGconsGetNVars(scip, cons1->getCons()), GCGconsGetNVars(scip, cons2->getCons()), FALSE);
   assert(SCIPisEQ(scip, GCGconsGetNVars(scip, cons1->getCons()), GCGconsGetNVars(scip, cons2->getCons())));

   return strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons1->getCons())), SCIPconshdlrGetName(SCIPconsGetHdlr(cons2->getCons())));
}


/** compare two variables of two scips */
static
int comp(
   SCIP*                 scip,               /**< SCIP data structure */
   AUT_VAR*              var1,               /**< variable 1 to compare */
   AUT_VAR*              var2,               /**< variable 2 to compare */
   SCIP_Bool             onlysign            /**< use sign of values instead of values? */
   )
{
   SCIP_VAR* origvar1;
   SCIP_VAR* origvar2;

   if( GCGvarIsPricing(var1->getVar()) )
         origvar1 = GCGpricingVarGetOriginalVar(var1->getVar());
      else
         origvar1 = var1->getVar();

   if( GCGvarIsPricing(var2->getVar()) )
         origvar2 = GCGpricingVarGetOriginalVar(var2->getVar());
      else
         origvar2 = var2->getVar();

   if( comp(scip, SCIPvarGetUbGlobal(origvar1), SCIPvarGetUbGlobal(origvar2), onlysign) != 0 )
      return comp(scip, SCIPvarGetUbGlobal(origvar1), SCIPvarGetUbGlobal(origvar2), onlysign);
   assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(origvar1), SCIPvarGetUbGlobal(origvar2)) || onlysign);

   if( comp(scip, SCIPvarGetLbGlobal(origvar1), SCIPvarGetLbGlobal(origvar2), onlysign) != 0 )
      return comp(scip, SCIPvarGetLbGlobal(origvar1), SCIPvarGetLbGlobal(origvar2), onlysign);
   assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(origvar1), SCIPvarGetLbGlobal(origvar2)) || onlysign);

   if( comp(scip, SCIPvarGetObj((origvar1)), SCIPvarGetObj(origvar2), onlysign) != 0 )
      return comp(scip, SCIPvarGetObj(origvar1), SCIPvarGetObj(origvar2), onlysign);
   assert(SCIPisEQ(scip, SCIPvarGetObj(origvar1), SCIPvarGetObj(origvar2)) || onlysign);

   if( SCIPvarGetType(origvar1) < SCIPvarGetType(origvar2) )
      return -1;
   if( SCIPvarGetType(origvar1) > SCIPvarGetType(origvar2) )
      return 1;
   return 0;
}

/** SCIP interface method for sorting the constraints */
static
SCIP_DECL_SORTPTRCOMP(sortptrcons)
{
   AUT_CONS* aut1 = (AUT_CONS*) elem1;
   AUT_CONS* aut2 = (AUT_CONS*) elem2;
   return comp(aut1->getScip(), aut1, aut2, FALSE);
}

/** SCIP interface method for sorting the constraints */
static
SCIP_DECL_SORTPTRCOMP(sortptrconssign)
{
   AUT_CONS* aut1 = (AUT_CONS*) elem1;
   AUT_CONS* aut2 = (AUT_CONS*) elem2;
   return comp(aut1->getScip(), aut1, aut2, TRUE);
}

/** SCIP interface method for sorting the variables */
static
SCIP_DECL_SORTPTRCOMP(sortptrvar)
{
   AUT_VAR* aut1 = (AUT_VAR*) elem1;
   AUT_VAR* aut2 = (AUT_VAR*) elem2;
   return comp(aut1->getScip(), aut1, aut2, FALSE);
}

/** SCIP interface method for sorting the variables */
static
SCIP_DECL_SORTPTRCOMP(sortptrvarsign)
{
   AUT_VAR* aut1 = (AUT_VAR*) elem1;
   AUT_VAR* aut2 = (AUT_VAR*) elem2;
   return comp(aut1->getScip(), aut1, aut2, TRUE);
}

/** SCIP interface method for sorting the constraint coefficients*/
static
SCIP_DECL_SORTPTRCOMP(sortptrval)
{
   AUT_COEF* aut1 = (AUT_COEF*) elem1;
   AUT_COEF* aut2 = (AUT_COEF*) elem2;
   return comp(aut1->getScip(), aut1->getVal(), aut2->getVal(), FALSE); /*lint !e864*/
}

/** SCIP interface method for sorting the constraint coefficients*/
static
SCIP_DECL_SORTPTRCOMP(sortptrvalsign)
{
   AUT_COEF* aut1 = (AUT_COEF*) elem1;
   AUT_COEF* aut2 = (AUT_COEF*) elem2;
   return comp(aut1->getScip(), aut1->getVal(), aut2->getVal(), TRUE); /*lint !e864*/
}


/** default constructor */
struct_colorinformation::struct_colorinformation()
 : color(0), lenconssarray(0), lenvarsarray(0), lencoefsarray(0), alloccoefsarray(0),
ptrarraycoefs(NULL), ptrarrayvars(NULL), ptrarrayconss(NULL), onlysign(FALSE)
{

}

/** inserts a variable to the pointer array of colorinformation */
SCIP_RETCODE struct_colorinformation::insert(
   AUT_VAR*              svar,               /**< variable which is to add */
   SCIP_Bool*            added               /**< true if a var was added */
   )
{
   int pos;

   if( !onlysign )
   {
      if( !SCIPsortedvecFindPtr(ptrarrayvars, sortptrvar, svar, lenvarsarray, &pos) )
      {
         SCIPsortedvecInsertPtr(ptrarrayvars, sortptrvar, svar, &lenvarsarray, NULL);
         *added = TRUE;
         color++;
      }
      else
         *added = FALSE;
   }
   else
   {
      if( !SCIPsortedvecFindPtr(ptrarrayvars, sortptrvarsign, svar, lenvarsarray, &pos) )
      {
         SCIPsortedvecInsertPtr(ptrarrayvars, sortptrvarsign, svar, &lenvarsarray, NULL);
         *added = TRUE;
         color++;
      }
      else
         *added = FALSE;
   }


   return SCIP_OKAY;
}

/** inserts a constraint to the pointer array of colorinformation */
SCIP_RETCODE struct_colorinformation::insert(
   AUT_CONS*             scons,              /**< constraint which is to add */
   SCIP_Bool*            added               /**< true if a constraint was added */
   )
{
   int pos;

   if( !onlysign )
   {
      if( !SCIPsortedvecFindPtr(ptrarrayconss, sortptrcons, scons,
            lenconssarray, &pos) )
      {
         SCIPsortedvecInsertPtr(ptrarrayconss, sortptrcons, scons,
               &lenconssarray, NULL);
         *added = TRUE;
         color++;
      }
      else
         *added = FALSE;
   }
   else
   {
      if( !SCIPsortedvecFindPtr(ptrarrayconss, sortptrconssign, scons,
            lenconssarray, &pos) )
      {
         SCIPsortedvecInsertPtr(ptrarrayconss, sortptrconssign, scons,
               &lenconssarray, NULL);
         *added = TRUE;
         color++;
      }
      else
         *added = FALSE;
   }

   return SCIP_OKAY;
}

/** inserts a coefficient to the pointer array of colorinformation */
SCIP_RETCODE struct_colorinformation::insert(
   AUT_COEF*             scoef,              /**< coefficient which is to add */
   SCIP_Bool*            added               /**< true if a coefficient was added */
   )
{
   int pos;

   if( !onlysign )
   {
      if( !SCIPsortedvecFindPtr(ptrarraycoefs, sortptrval, scoef, lencoefsarray, &pos) )
      {
         if( alloccoefsarray == 0 || alloccoefsarray < lencoefsarray + 1 )
         {
            int size = SCIPcalcMemGrowSize(scoef->getScip(), alloccoefsarray+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, &ptrarraycoefs, size) );
            alloccoefsarray = size;
         }

         SCIPsortedvecInsertPtr(ptrarraycoefs, sortptrval, scoef, &lencoefsarray, NULL);
         *added = TRUE;
         color++;
      }
      else
         *added = FALSE;
   }
   else
   {
      if( !SCIPsortedvecFindPtr(ptrarraycoefs, sortptrvalsign, scoef, lencoefsarray, &pos) )
      {
         if( alloccoefsarray == 0 || alloccoefsarray < lencoefsarray + 1 )
         {
            int size = SCIPcalcMemGrowSize(scoef->getScip(), alloccoefsarray+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, &ptrarraycoefs, size) );
            alloccoefsarray = size;
         }

         SCIPsortedvecInsertPtr(ptrarraycoefs, sortptrvalsign, scoef, &lencoefsarray, NULL);
         *added = TRUE;
         color++;
      }
      else
         *added = FALSE;
   }

   return SCIP_OKAY;
}

int struct_colorinformation::get(
   AUT_VAR               svar                /**< variable whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   if( !onlysign )
      found = SCIPsortedvecFindPtr(ptrarrayvars, sortptrvar, &svar, lenvarsarray, &pos);
   else
      found = SCIPsortedvecFindPtr(ptrarrayvars, sortptrvarsign, &svar, lenvarsarray, &pos);
   return found ? pos : -1;
}

int struct_colorinformation::get(
   AUT_CONS              scons               /**< constraint whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   if( !onlysign )
      found = SCIPsortedvecFindPtr(ptrarrayconss, sortptrcons, &scons, lenconssarray, &pos);
   else
      found = SCIPsortedvecFindPtr(ptrarrayconss, sortptrconssign, &scons, lenconssarray, &pos);
   return found ? pos : -1;
}

int struct_colorinformation::get(
   AUT_COEF              scoef               /**< coefficient whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   if( !onlysign )
      found = SCIPsortedvecFindPtr(ptrarraycoefs, sortptrval, &scoef, lencoefsarray, &pos);
   else
      found = SCIPsortedvecFindPtr(ptrarraycoefs, sortptrvalsign, &scoef, lencoefsarray, &pos);
   return found ? pos : -1;
}

SCIP_RETCODE struct_colorinformation::setOnlySign(
   SCIP_Bool            onlysign_            /**< new value for onlysign bool */
   )
{
   onlysign = onlysign_;

   return SCIP_OKAY;
}


SCIP_Bool struct_colorinformation::getOnlySign()
{
   return onlysign;
}


int struct_colorinformation::getLenVar()
{
   return lenvarsarray;
}

int struct_colorinformation::getLenCons()
{
   return lenconssarray;
}

SCIP_CONS* struct_cons::getCons()
{
   return cons;
}

SCIP* struct_cons::getScip()
{
   return scip;
}

SCIP_VAR* struct_var::getVar ()
{
   return var;
}

SCIP* struct_var::getScip()
{
   return scip;
}

SCIP* struct_coef::getScip()
{
   return scip;
}

SCIP_Real struct_coef::getVal()
{
   return val;
}

/** constructor of the variable struct */
struct_var::struct_var(
   SCIP*                 scip_,              /**< SCIP data structure */
   SCIP_VAR*             svar                /**< SCIP variable */
   )
{
   scip = scip_;
   var = svar;
}

/** constructor of the constraint struct */
struct_cons::struct_cons(
   SCIP*                 scip_,              /**< SCIP data structure */
   SCIP_CONS*            scons               /**< SCIP constraint */
   )
{
   scip = scip_;
   cons = scons;
}

/** constructor of the coefficient struct */
struct_coef::struct_coef(
   SCIP*                 scip_,              /**< SCIP data structure */
   SCIP_Real             val_                /**< SCIP value */
   )
{
   scip = scip_;
   val = val_;
}
