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

/**@file   hyperrowgraph_def.h
 * @brief  Column hypergraph
 * @author Martin Bergner
 * @author Annika Thome
 *
 * Hypergraph with a node for every variable and a hyperedge for every constraint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#ifndef GCG_HYPERROWGRAPH_DEF_H_
#define GCG_HYPERROWGRAPH_DEF_H_

#include "graph/hyperrowgraph.h"
#include "gcg/scip_misc.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include <set>
#include <algorithm>
#include <iostream>

namespace gcg
{
template <class T>
HyperrowGraph<T>::HyperrowGraph(
   GCG*                  gcgstruct,         /**< GCG data structure */
   Weights               w                  /**< weights for the given graph */
): MatrixGraph<T>(gcgstruct, w), graph(gcgstruct)
{
   this->graphiface = &graph;
   this->name = std::string("hyperrow");
}

template <class T>
HyperrowGraph<T>::~HyperrowGraph()
{
   // TODO Auto-generated destructor stub
}


/** writes the graph to the given file.
 *  The format is graph dependent
 */
template <class T>
SCIP_RETCODE HyperrowGraph<T>::writeToFile(
   int                fd,                    /**< filename where the graph should be written to */
   SCIP_Bool          edgeweights            /**< whether to write edgeweights */
   )
{
   FILE* file;
   SCIP* scip = GCGgetOrigprob(this->gcg);
   file = fdopen(fd, "w");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(scip, file, "%d %d %d\n", getNEdges(), getNNodes()+this->dummynodes, edgeweights ? 1 :0);

   for( int i = 0; i < getNEdges(); ++i )
   {
      std::vector<int> neighbors = getHyperedgeNodes(i);

      if( edgeweights )
      {
         SCIPinfoMessage(scip, file, "%d ", graph.getWeight(i+this->nvars));
      }
      for( size_t j = 0; j < neighbors.size(); ++j )
      {
         SCIPinfoMessage(scip, file, "%d ",neighbors[j]+1);
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   if( !fclose(file) )
      return SCIP_OKAY;
   else
      return SCIP_WRITEERROR;
}

template <class T>
int HyperrowGraph<T>::getNEdges()
{
   return this->nconss;
}

template <class T>
int HyperrowGraph<T>::getNNodes()
{
   return this->nvars;
}

template <class T>
int HyperrowGraph<T>::getNNeighbors(
   int i
   )
{
   assert(i >= 0);
   assert(i < getNNodes());

   return graph.getNNeighbors(i);
}

template <class T>
std::vector<int> HyperrowGraph<T>::getHyperedgeNodes(
   int i
   )
{
   assert(i >= 0);
   assert(i < getNEdges());

   std::vector<int> neighbors = graph.getHyperedgeNodes(i);
   return neighbors;
}

template <class T>
SCIP_RETCODE HyperrowGraph<T>::createDecompFromPartition(
   GCG_DECOMP**       decomp              /**< decomposition structure to generate */
   )
{
   int nblocks;
   SCIP_HASHMAP* constoblock = NULL;
   SCIP* scip = GCGgetOrigprob(this->gcg);
   int* nsubscipconss = NULL;
   int i;
   SCIP_CONS** conss = NULL;
   SCIP_Bool emptyblocks = FALSE;
   std::vector<int> partition = graph.getPartition();
   conss = SCIPgetConss(scip);
   nblocks = *(std::max_element(partition.begin(), partition.end()))+1;

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   BMSclearMemoryArray(nsubscipconss, nblocks);

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), this->nconss) );

   /* assign constraints to partition */
   for( i = 0; i < this->nconss; i++ )
   {

      std::set<int> blocks;
      std::vector<int> neighbors = getHyperedgeNodes(i);
      for( size_t k = 0; k < neighbors.size(); ++k )
      {
         if( partition[neighbors[k]] >= 0 )
            blocks.insert(partition[neighbors[k]]);
      }
      if( blocks.size() > 1 )
      {
         SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) (nblocks+1)) );
      }
      else
      {
         int block = *(blocks.begin());
         SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) (block +1)) );
         ++(nsubscipconss[block]);
      }
   }

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < nblocks; ++i )
   {
      if( nsubscipconss[i] == 0 )
      {
         SCIPdebugMessage("Block %d does not have any constraints!\n", i);
         emptyblocks = TRUE;
      }
   }

   if( !emptyblocks )
   {
      SCIP_CALL( GCGdecompCreate(this->gcg, decomp) );
      SCIP_CALL( GCGfilloutDecompFromConstoblock(this->gcg, *decomp, constoblock, nblocks, FALSE) );
   }
   else {
      SCIPhashmapFree(&constoblock);
      *decomp = NULL;
   }

   SCIPfreeBufferArray(scip, &nsubscipconss);
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE HyperrowGraph<T>::createPartialdecFromPartition(
   PARTIALDECOMP**     firstpartialdec,
   PARTIALDECOMP**     secondpartialdec,
   DETPROBDATA*        detprobdata
   )
{
   int nblocks;
   SCIP_HASHMAP* constoblock = NULL;
   SCIP* scip = GCGgetOrigprob(this->gcg);
   int* nsubscipconss = NULL;
   int i;
   SCIP_CONS** conss = NULL;
   SCIP_Bool emptyblocks = FALSE;
   std::vector<int> partition;

   if( firstpartialdec == NULL && secondpartialdec == NULL )
      return SCIP_INVALIDDATA;

   partition = graph.getPartition();
   conss = SCIPgetConss(scip);
   nblocks = *(std::max_element(partition.begin(), partition.end()))+1;

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   BMSclearMemoryArray(nsubscipconss, nblocks);

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), this->nconss) );

   /* assign constraints to partition */
   for( i = 0; i < this->nconss; i++ )
   {

      std::set<int> blocks;
      std::vector<int> neighbors = getHyperedgeNodes(i);
      for( size_t k = 0; k < neighbors.size(); ++k )
      {
         if( partition[neighbors[k]] >= 0 )
            blocks.insert(partition[neighbors[k]]);
      }
      if( blocks.size() > 1 )
      {
         SCIP_CALL( SCIPhashmapInsert(constoblock, (void*) (size_t)detprobdata->getIndexForCons(conss[i]), (void*) (size_t) (nblocks+1)) );
      }
      else
      {
         int block = *(blocks.begin());
         SCIP_CALL( SCIPhashmapInsert(constoblock, (void*) (size_t)detprobdata->getIndexForCons(conss[i]), (void*) (size_t) (block +1)) );
         ++(nsubscipconss[block]);
      }
   }

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < nblocks; ++i )
   {
      if( nsubscipconss[i] == 0 )
      {
         SCIPdebugMessage("Block %d does not have any constraints!\n", i);
         emptyblocks = TRUE;
      }
   }

   if( !emptyblocks )
   {
      bool original = detprobdata->isAssignedToOrigProb();
      if( firstpartialdec != NULL )
      {
         (*firstpartialdec) = new PARTIALDECOMP(this->gcg, original);
         SCIP_CALL((*firstpartialdec)->filloutPartialdecFromConstoblock(constoblock, nblocks));
      }
      if( secondpartialdec != NULL )
      {
         (*secondpartialdec) = new PARTIALDECOMP(this->gcg, original);
         SCIP_CALL((*secondpartialdec)->filloutBorderFromConstoblock(constoblock, nblocks));
      }
      SCIPhashmapFree(&constoblock);
   }
   else {
      SCIPhashmapFree(&constoblock);
      if( firstpartialdec != NULL )
      {
         *firstpartialdec = NULL;
      }
      if( secondpartialdec != NULL )
      {
         *secondpartialdec = NULL;
      }
   }

   SCIPfreeBufferArray(scip, &nsubscipconss);
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE HyperrowGraph<T>::createPartialdecFromPartition(
   PARTIALDECOMP*      oldpartialdec,
   PARTIALDECOMP**     firstpartialdec,
   PARTIALDECOMP**     secondpartialdec,
   DETPROBDATA*        detprobdata
   )
{
   int nblocks;
   SCIP_HASHMAP* constoblock = NULL;
   SCIP* scip = GCGgetOrigprob(this->gcg);
   int *nsubscipconss = NULL;
   int i;
   SCIP_Bool emptyblocks = FALSE;

   assert(oldpartialdec != NULL);

   if( firstpartialdec == NULL && secondpartialdec == NULL )
      return SCIP_INVALIDDATA;

   if( this->nconss == 0 )
   {
      if( firstpartialdec != NULL )
         (*firstpartialdec) = NULL;
      if( secondpartialdec != NULL )
         (*secondpartialdec) = NULL;
      return SCIP_OKAY;
   }

   std::vector<int> partition = graph.getPartition();
   nblocks = *(std::max_element(partition.begin(), partition.end()))+1;

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   BMSclearMemoryArray(nsubscipconss, nblocks);

   for( int b = 0; b < nblocks; ++b )
   {
       nsubscipconss[b] = 0;
   }

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), this->nconss) );

   //fillout conssForGraph
   vector<int> conssForGraph; /** stores the conss included by the graph */
   vector<bool> conssBool(oldpartialdec->getNConss(), false); /**< true, if the cons will be part of the graph */
   bool found;

   for( int c = 0; c < oldpartialdec->getNOpenconss(); ++c )
   {
      int cons = oldpartialdec->getOpenconss()[c];
      found = false;
      for( int v = 0; v < oldpartialdec->getNOpenvars() && !found; ++v )
      {
         int var = oldpartialdec->getOpenvars()[v];
         for( i = 0; i < detprobdata->getNVarsForCons(cons) && !found; ++i )
         {
            if( var == detprobdata->getVarsForCons(cons)[i] )
            {
               conssBool[cons] = true;
               found = true;
            }
         }
      }
   }

   for( int c = 0; c < oldpartialdec->getNOpenconss(); ++c )
   {
      int cons = oldpartialdec->getOpenconss()[c];
      if( conssBool[cons] )
         conssForGraph.push_back(cons);
   }

   /* assign constraints to partition */
   for( i = 0; i < this->nconss; i++ )
   {

      std::set<int> blocks;
      std::vector<int> neighbors = getHyperedgeNodes(i);
      for( size_t k = 0; k < neighbors.size(); ++k )
      {
         if( partition[neighbors[k]] >= 0 )
            blocks.insert(partition[neighbors[k]]);
      }
      if( blocks.size() > 1 )
      {
         SCIP_CALL( SCIPhashmapInsert(constoblock, (void*) (size_t) conssForGraph[i], (void*) (size_t) (nblocks+1)) );
      }
      else
      {
         int block = *(blocks.begin());
         SCIP_CALL( SCIPhashmapInsert(constoblock, (void*) (size_t) conssForGraph[i], (void*) (size_t) (block +1)) );
         ++(nsubscipconss[block]);
      }
   }

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < nblocks; ++i )
   {
      if( nsubscipconss[i] == 0 )
      {
         SCIPdebugMessage("Block %d does not have any constraints!\n", i);
         emptyblocks = TRUE;
      }
   }

   if( !emptyblocks )
   {
      if( firstpartialdec != NULL )
      {
         (*firstpartialdec) = new PARTIALDECOMP(oldpartialdec);
         SCIP_CALL( (*firstpartialdec)->assignPartialdecFromConstoblock(constoblock, nblocks) );
      }
      if( secondpartialdec != NULL )
      {
         (*secondpartialdec) = new PARTIALDECOMP(oldpartialdec);
         SCIP_CALL( (*secondpartialdec)->assignBorderFromConstoblock(constoblock, nblocks) );
      }
      SCIPhashmapFree(&constoblock);
   }
   else
   {
      SCIPhashmapFree(&constoblock);
      *firstpartialdec = NULL;
      *secondpartialdec = NULL;
   }

   SCIPfreeBufferArray(scip, &nsubscipconss);
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE HyperrowGraph<T>::createFromMatrix(
   SCIP_CONS**           conss,              /**< constraints for which graph should be created */
   SCIP_VAR**            vars,               /**< variables for which graph should be created */
   int                   nconss_,             /**< number of constraints */
   int                   nvars_               /**< number of variables */
   )
{
   int i;
   int j;
   SCIP_Bool success;
   SCIP* scip = GCGgetOrigprob(this->gcg);

   assert(conss != NULL);
   assert(vars != NULL);
   assert(nvars_ > 0);
   assert(nconss_ > 0);

   this->nvars = nvars_;
   this->nconss = nconss_;

   /* go through all variables */
   for( i = 0; i < this->nvars; ++i )
   {
      TCLIQUE_WEIGHT weight;

      /* calculate weight of node */
      weight = this->weights.calculate(vars[i]);

      this->graph.addNode(i, weight);
   }

   /* go through all constraints */
   for( i = 0; i < this->nconss; ++i )
   {
      SCIP_VAR** curvars = NULL;
      std::vector<int> hyperedge;
      TCLIQUE_WEIGHT weight;

      int ncurvars;
      SCIP_CALL( SCIPgetConsNVars(scip, conss[i], &ncurvars, &success) );
      assert(success);
      if( ncurvars == 0 )
         continue;

      /*
       * may work as is, as we are copying the constraint later regardless
       * if there are variables in it or not
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
      SCIP_CALL( SCIPgetConsVars(scip, conss[i], curvars, ncurvars, &success) );
      assert(success);

      /** @todo skip all variables that have a zero coeffient or where all coefficients add to zero */
      /** @todo Do more then one entry per variable actually work? */

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* var1 = NULL;
         int varIndex1;

         if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED)
            var1 = SCIPvarGetProbvar(curvars[j]);
         else
            var1 = curvars[j];

         if( !GCGisVarRelevant(var1) )
            continue;

         assert(var1 != NULL);
         varIndex1 = SCIPvarGetProbindex(var1);
         assert(varIndex1 >= 0);
         assert(varIndex1 < this->nvars);

         hyperedge.insert(hyperedge.end(), varIndex1);
      }
      /* calculate weight of hyperedge */
      weight = this->weights.calculate(conss[i]);

      this->graph.addHyperedge(hyperedge, weight);

      SCIPfreeBufferArray(scip, &curvars);
   }



   this->graph.flush();
   return SCIP_OKAY;
}



template <class T>
SCIP_RETCODE HyperrowGraph<T>::createFromPartialMatrix(
   DETPROBDATA* detprobdata,
   PARTIALDECOMP* partialdec
   )
{
     int i;
     int j;
     unordered_map<int, int> oldToNewVarIndex;
     TCLIQUE_WEIGHT weight;

     vector<bool> varsBool(partialdec->getNVars(), false); /**< true, if the var will be part of the graph */
     vector<bool> conssBool(partialdec->getNConss(), false); /**< true, if the cons will be part of the graph */
     vector<int> conssForGraph; /** stores the conss included by the graph */
     vector<int> varsForGraph; /** stores the vars included by the graph */

     //fillout conssForGraph and varsForGraph
     for(int c = 0; c < partialdec->getNOpenconss(); ++c)
     {
        int cons = partialdec->getOpenconss()[c];
        for(int v = 0; v < partialdec->getNOpenvars(); ++v)
        {
           int var = partialdec->getOpenvars()[v];
           for(i = 0; i < detprobdata->getNVarsForCons(cons); ++i)
           {
              if(var == detprobdata->getVarsForCons(cons)[i])
              {
                 varsBool[var] = true;
                 conssBool[cons] = true;
              }
           }
        }
     }

     for(int v = 0; v < partialdec->getNOpenvars(); ++v)
     {
        int var = partialdec->getOpenvars()[v];
        if(varsBool[var])
           varsForGraph.push_back(var);
     }
     for(int c = 0; c < partialdec->getNOpenconss(); ++c)
     {
        int cons = partialdec->getOpenconss()[c];
        if(conssBool[cons])
           conssForGraph.push_back(cons);
     }

     this->nconss = (int)conssForGraph.size();
     this->nvars = (int)varsForGraph.size();

     /* go through all variables */
     for( i = 0; i < this->nvars; ++i )
     {
        int oldVarId = varsForGraph[i];
        assert(varsBool[oldVarId]);

        /* calculate weight of node */
        weight = this->weights.calculate(detprobdata->getVar(oldVarId));

        oldToNewVarIndex.insert({oldVarId,i});
        this->graph.addNode(i, weight);
     }

     /* go through all open constraints */
     for( i = 0; i < this->nconss; ++i )
     {
        std::vector<int> hyperedge;
        int oldConsId = conssForGraph[i];

        assert(conssBool[oldConsId]);

        for( j = 0; j < detprobdata->getNVarsForCons(oldConsId); ++j )
        {
           int oldVarId = detprobdata->getVarsForCons(oldConsId)[j];
           if(!varsBool[oldVarId])
              continue;
           hyperedge.insert(hyperedge.end(), oldToNewVarIndex[oldVarId]);
        }
        /* calculate weight of hyperedge */
        weight = this->weights.calculate(detprobdata->getCons(oldConsId));
        this->graph.addHyperedge(hyperedge, weight);
     }


     this->graph.flush();
     return SCIP_OKAY;
}


} /* namespace gcg */

#endif
