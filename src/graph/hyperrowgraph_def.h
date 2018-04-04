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

#include "hyperrowgraph.h"
#include "scip_misc.h"
#include <set>
#include <algorithm>

namespace gcg
{
template <class T>
HyperrowGraph<T>::HyperrowGraph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               w                  /**< weights for the given graph */
): MatrixGraph<T>(scip, w), graph(scip)
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
   file = fdopen(fd, "w");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   SCIPinfoMessage(this->scip_, file, "%d %d %d\n", getNEdges(), getNNodes()+this->dummynodes, edgeweights ? 1 :0);

   for( int i = 0; i < getNEdges(); ++i )
   {
      std::vector<int> neighbors = getHyperedgeNodes(i);

      if( edgeweights )
      {
         SCIPinfoMessage(this->scip_, file, "%d ", graph.getWeight(i+this->nvars));
      }
      for( size_t j = 0; j < neighbors.size(); ++j )
      {
         SCIPinfoMessage(this->scip_, file, "%d ",neighbors[j]+1);
      }
      SCIPinfoMessage(this->scip_, file, "\n");
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
   DEC_DECOMP**       decomp              /**< decomposition structure to generate */
)
{
   int nblocks;
   SCIP_HASHMAP* constoblock;

   int *nsubscipconss;
   int i;
   SCIP_CONS **conss;
   SCIP_Bool emptyblocks = FALSE;
   std::vector<int> partition = graph.getPartition();
   conss = SCIPgetConss(this->scip_);
   nblocks = *(std::max_element(partition.begin(), partition.end()))+1;

   SCIP_CALL( SCIPallocBufferArray(this->scip_, &nsubscipconss, nblocks) );
   BMSclearMemoryArray(nsubscipconss, nblocks);

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(this->scip_), this->nconss) );

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
      SCIP_CALL( DECdecompCreate(this->scip_, decomp) );
      SCIP_CALL( DECfilloutDecompFromConstoblock(this->scip_, *decomp, constoblock, nblocks, FALSE) );
   }
   else {
      SCIPhashmapFree(&constoblock);
      *decomp = NULL;
   }

   SCIPfreeBufferArray(this->scip_, &nsubscipconss);
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
      SCIP_VAR **curvars;
      std::vector<int> hyperedge;
      TCLIQUE_WEIGHT weight;

      int ncurvars;
      SCIP_CALL( SCIPgetConsNVars(this->scip_, conss[i], &ncurvars, &success) );
      assert(success);
      if( ncurvars == 0 )
         continue;

      /*
       * may work as is, as we are copying the constraint later regardless
       * if there are variables in it or not
       */
      SCIP_CALL( SCIPallocBufferArray(this->scip_, &curvars, ncurvars) );
      SCIP_CALL( SCIPgetConsVars(this->scip_, conss[i], curvars, ncurvars, &success) );
      assert(success);

      /** @todo skip all variables that have a zero coeffient or where all coefficients add to zero */
      /** @todo Do more then one entry per variable actually work? */

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* var1;
         int varIndex1;

         if( SCIPgetStage(this->scip_) >= SCIP_STAGE_TRANSFORMED)
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

      SCIPfreeBufferArray(this->scip_, &curvars);
   }



   this->graph.flush();
   return SCIP_OKAY;
}


} /* namespace gcg */

#endif
