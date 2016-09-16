/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2016 Operations Research, RWTH Aachen University       */
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

/**@file   bipartitegraph_def.h
 * @brief  A bipartite graph
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_BIPARTITEGRAPH_DEF_H_
#define GCG_BIPARTITEGRAPH_DEF_H_

#include "bipartitegraph.h"
#include "scip_misc.h"

namespace gcg {

template <class T>
BipartiteGraph<T>::BipartiteGraph(
      SCIP*                 scip,              /**< SCIP data structure */
      Weights               w                 /**< weights for the given graph */
   ): MatrixGraph<T>(scip,w), graph(scip)
{
   this->graphiface = &graph;
   this->name = std::string("bipartite");
}

template <class T>
BipartiteGraph<T>::~BipartiteGraph()
{

}


/**
 * Builds a bipartite graph structure out of the matrix.
 *
 * The function will create an node for every constraint and every variable.
 * A constraint and a variable are adjacent if the variable appears in the constraint variable array.
 *
 * @todo The nonzeroness is not checked, all variables in the variable array are considered
 */
template <class T>
SCIP_RETCODE BipartiteGraph<T>::createFromMatrix(
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

   for( i = 0; i < this->nvars + this->nconss; ++i )
   {
      TCLIQUE_WEIGHT weight;

      /* note that the first nvars nodes correspond to variables */
      if( i < this->nvars )
         weight = this->weights.calculate(vars[i]);
      else
         weight = this->weights.calculate(conss[i-this->nvars]);

      this->graph.addNode(i, weight);
   }

   /* go through all constraints */
   for( i = 0; i < this->nconss; ++i )
   {
      SCIP_VAR **curvars;

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
         SCIP_VAR* var;
         int varIndex;

         if( SCIPgetStage(this->scip_) >= SCIP_STAGE_TRANSFORMED)
            var = SCIPvarGetProbvar(curvars[j]);
         else
            var = curvars[j];

         if( !GCGisVarRelevant(var) )
            continue;

         assert(var != NULL);
         varIndex = SCIPvarGetProbindex(var);
         assert(varIndex >= 0);
         assert(varIndex < this->nvars);

         SCIP_CALL( this->graph.addEdge(varIndex, this->nvars+i) );
      }
      SCIPfreeBufferArray(this->scip_, &curvars);
   }

   this->graph.flush();
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE BipartiteGraph<T>::createFromPartialMatrix(
                   std::vector<std::vector<int>>                                varsForConss,           /** stores for every constraint the indices of variables that are contained in the constraint */
                   std::vector<std::vector<int>>                                conssForVars,           /** stores for every variable the indices of constraints containing this variable */
                   std::vector<int>                                                     openVars,           /**< vector of variable indices that are not assigned yet */
                   std::vector<int>                                                     openConss,          /**< vector of constraint indices that are not assigned yet */
                   std::vector<SCIP_CONS*>                                              consToScipCons,     /** stores the corresponding scip constraints pointer */
                   std::vector<SCIP_VAR*>                                               varToScipVar,           /** stores the corresponding scip variable pointer */
                   int                                                                  nconss_,            /**< number of constraints */
                   int                                                                  nvars_              /**< number of variables */
     ){

         int i;
     int j;
     SCIP_Bool success;
     std::tr1::unordered_map<int, int> oldToNewVarIndex;
     std::tr1::unordered_map<int, int> oldToNewConsIndex;

     assert(varsForConss.size() == nconss_);
     assert(conssForVars.size() == nvars_);
     assert(nvars_ > 0);
     assert(nconss_ > 0);
     this->nvars = openVars.size();
     this->nconss = openConss.size;

     std::vector<int>::const_iterator varIter = openVars.begin();
     std::vector<int>::const_iterator varIterEnd = openVars.end();

     std::vector<int>::const_iterator consIter = openConss.begin();
     std::vector<int>::const_iterator consIterEnd = openConss.end();




     /** add node for every var */
     for( int i = 0 ; i < openVars.size(); ++i )
     {
         TCLIQUE_WEIGHT weight;

         /* note that the first nvars nodes correspond to variables */
         weight = this->weights.calculate(varToScipVar(openVars[i]) );
         oldToNewVarIndex.insert({ openVars[i],i});
         this->graph.addNode(i, weight);
     }


     /** add node for every cons */
     for( int j = 0 ; j < openConss.size(); ++j  )
     {
        TCLIQUE_WEIGHT weight;

        /* note that the first nvars nodes correspond to variables (legacy implementation) */
        weight = this->weights.calculate(consToScipCons(openConss[j] ) );
        oldToNewVarIndex.insert({ openConss[j], j});
        this->graph.addNode( this->nvars + j, weight);
     }

     /* go through all open constraints */
     for( i = 0; i < openConss.size(); ++i )
     {
        int oldConsId = openConss[i];

        std::vector<int>::const_iterator curVarIter = varsForConss[oldConsId].begin();
        std::vector<int>::const_iterator curVarIterEnd = varsForConss[oldConsId].end();

        for( ; curVarIter != curVarIterEnd; ++curVarIter )
        {
                if (oldToNewVarIndex.find(*curVarIter) == oldToNewVarIndex.end() )
                        continue;
                SCIP_CALL( this->graph.addEdge(oldToNewVarIndex[*curVarIter], this->nvars+i) );
        }
     }

     this->graph.flush();
     return SCIP_OKAY;
  }




template <class T>
int BipartiteGraph<T>::getNConsNodes()
{
   return this->nconss;
}

template <class T>
int BipartiteGraph<T>::getNVarNodes()
{
   return this->nvars;
}


} /* namespace gcg */

#endif
