/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   rowgraph.h
 * @brief  A row graph where each row is a node and rows are adjacent if they share a variable
 * @author Martin Bergner
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_ROWGRAPH_H_
#define GCG_ROWGRAPH_H_

#include "graph/graph.h"
#include "graph/bipartitegraph.h"
#include "graph/matrixgraph.h"

namespace gcg {
template <class T>
class RowGraph : public gcg::MatrixGraph<T>
{
protected:
   gcg::Graph<T> graph;

public:
   RowGraph(
      GCG*                  gcg,               /**< GCG data structure */
      Weights               w                  /**< weights for the given graph */
      );
   virtual ~RowGraph();

   virtual SCIP_RETCODE createDecompFromPartition(
      GCG_DECOMP**       decomp              /**< decomposition structure to generate */
      );

   /** amplifies a partialdec by dint of a graph created with open constraints and open variables of the partialdec */
   virtual SCIP_RETCODE createPartialdecFromPartition(
      PARTIALDECOMP*      oldpartialdec,            /**< partialdec which should be amplifies */
      PARTIALDECOMP**     firstpartialdec,          /**< pointer to buffer the new partialdec amplified by dint of the graph (can be NULL) */
      PARTIALDECOMP**     secondpartialdec,         /**< pointer to buffer the new partialdec whose border is amplified by dint of the graph (can be NULL) */
      DETPROBDATA*        detprobdata          /**< datprobdata the partialdecs correspond to */
      );

   virtual SCIP_RETCODE createFromMatrix(
      SCIP_CONS**           conss,              /**< constraints for which graph should be created */
      SCIP_VAR**            vars,               /**< variables for which graph should be created */
      int                   nconss_,            /**< number of constraints */
      int                   nvars_              /**< number of variables */
      );
};

} /* namespace gcg */
#endif /* GCG_ROWGRAPH_H_ */
