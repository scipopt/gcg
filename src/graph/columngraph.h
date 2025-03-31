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

/**@file   columngraph.h
 * @brief  A row graph where each column is a node and columns are adjacent if they appear in a row
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_COLUMNGRAPH_H_
#define GCG_COLUMNGRAPH_H_

#include "graph/graph.h"
#include "graph/bipartitegraph.h"
#include "graph/matrixgraph.h"

namespace gcg {
template <class T>
class ColumnGraph: public gcg::MatrixGraph<T>
{
private:
   gcg::Graph<T> graph;

public:
   ColumnGraph(
      GCG*                  gcg,               /**< GCG data structure */
      Weights               w                  /**< weights for the given graph */
      );
   virtual ~ColumnGraph();

   virtual SCIP_RETCODE createDecompFromPartition(
      GCG_DECOMP**      decomp                  /**< decomposition structure to generate */
      );

   virtual SCIP_RETCODE createFromMatrix(
      SCIP_CONS**           conss,              /**< constraints for which graph should be created */
      SCIP_VAR**            vars,               /**< variables for which graph should be created */
      int                   nconss_,             /**< number of constraints */
      int                   nvars_               /**< number of variables */
      );
};

} /* namespace gcg */
#endif /* GCG_COLUMNGRAPH_H_ */
