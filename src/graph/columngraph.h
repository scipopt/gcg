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

/**@file   columngraph.h
 * @brief  A row graph where each column is a node and columns are adjacent if they appear in a row
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_COLUMNGRAPH_H_
#define GCG_COLUMNGRAPH_H_

#include "graph.h"
#include "bipartitegraph.h"
#include "matrixgraph.h"

namespace gcg {
template <class T>
class ColumnGraph: public gcg::MatrixGraph<T>
{
private:
   gcg::Graph<T> graph;

public:
   ColumnGraph(
      SCIP*                 scip,              /**< SCIP data structure */
      Weights               w                  /**< weights for the given graph */
   );
   virtual ~ColumnGraph();
   //virtual SCIP_RETCODE writeToFile(
   //   const char*        filename,           /**< filename where the graph should be written to */
   //   SCIP_Bool          writeweights        /**< whether to write weights */
   //   );

   virtual SCIP_RETCODE createDecompFromPartition(
      DEC_DECOMP**      decomp                  /**< decomposition structure to generate */
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
