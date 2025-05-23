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

/**@file   bipartitegraph_def.h
 * @brief  A bipartite graph
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_BIPARTITEGRAPH_DEF_H_
#define GCG_BIPARTITEGRAPH_DEF_H_

#include "graph/bipartitegraph.h"
#include "gcg/scip_misc.h"

namespace gcg {

template <class T>
BipartiteGraph<T>::BipartiteGraph(
      GCG*                  gcgstruct,        /**< GCG data structure */
      Weights               w                 /**< weights for the given graph */
   ): MatrixGraph<T>(gcgstruct,w), graph(gcgstruct)
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
   SCIP* scip = GCGgetOrigprob(this->gcg);

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
      SCIP_VAR** curvars = NULL;

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
         SCIP_VAR* var = NULL;
         int varIndex;

         if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED)
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
      SCIPfreeBufferArray(scip, &curvars);
   }

   this->graph.flush();
   return SCIP_OKAY;
}


template <class T>
SCIP_RETCODE BipartiteGraph<T>::createFromPartialMatrix(
                   DETPROBDATA*                                                   detprobdata,
                   PARTIALDECOMP*                                               partialdec
     ){

     int i;
     int j;
     unordered_map<int, int> oldToNewVarIndex;
     unordered_map<int, int> oldToNewConsIndex;
     std::vector<bool> varsBool(partialdec->getNVars(), false); /* true, if the var will be part of the graph */
     std::vector<bool> conssBool(partialdec->getNConss(), false); /* true, if the cons will be part of the graph */
     std::vector<int> conssForGraph; /* stores the conss included by the graph */
     std::vector<int> varsForGraph; /* stores the vars included by the graph */

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


     /* add node for every var */
     for( i = 0 ; i < this->nvars; ++i )
     {
         TCLIQUE_WEIGHT weight;
         int var = varsForGraph[i];

         /* note that the first nvars nodes correspond to variables */
         weight = this->weights.calculate(detprobdata->getVar(var));
         oldToNewVarIndex.insert({var,i});
         this->graph.addNode(i, weight);
     }


     /* add node for every cons */
     for(  j = 0 ; j < this->nconss; ++j  )
     {
        TCLIQUE_WEIGHT weight;
        int cons = conssForGraph[j];

        /* note that the first nvars nodes correspond to variables (legacy implementation) */
        weight = this->weights.calculate(detprobdata->getCons(cons));
        oldToNewConsIndex.insert({cons, j});
        this->graph.addNode( this->nvars + j, weight);
     }

     /* go through all open constraints */
     for( i = 0; i < this->nconss; ++i )
     {
        int oldConsId = conssForGraph[i];

        for( j = 0; j < detprobdata->getNVarsForCons(oldConsId); ++j )
        {
           int oldVarId = detprobdata->getVarsForCons(oldConsId)[j];
           if(!varsBool[oldVarId])
              continue;
           SCIP_CALL( this->graph.addEdge(oldToNewVarIndex[oldVarId], this->nvars+i) );
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
