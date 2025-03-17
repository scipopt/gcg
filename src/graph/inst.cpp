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

/**@file   inst.cpp
 * @brief  Explicit instanciations for templates
 * @author Martin Bergner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
/*lint -e39*/
#include "graph/bipartitegraph_def.h"
#include "graph/columngraph_def.h"
#include "graph/rowgraph_def.h"
#include "graph/rowgraph_weighted_def.h"
#include "graph/hypercolgraph_def.h"
#include "graph/hyperrowgraph_def.h"
#include "graph/hyperrowcolgraph_def.h"
#include "graph/graph_def.h"
#include "graph/hypergraph_def.h"
#include "graph/graph_tclique.h"
#include "graph/graph_gcg.h"
#include "graph/matrixgraph_def.h"
#include "graph/graph_interface.h"
#include "graph/graphalgorithms_def.h"

namespace gcg {

/* graph instanciations for graphs using TCLIQUE graphs */
template class BipartiteGraph<GraphTclique>;
template class ColumnGraph<GraphTclique>;
template class RowGraph<GraphTclique>;
template class RowGraph<GraphGCG>;
template class RowGraphWeighted<GraphGCG>;
template class GraphAlgorithms<GraphGCG>;





template class HypercolGraph<GraphTclique>;
template class HyperrowGraph<GraphTclique>;
template class HyperrowcolGraph<GraphTclique>;
template class Graph<GraphTclique>;
template class Graph<GraphGCG>;
template class Hypergraph<GraphTclique>;
template class MatrixGraph<GraphTclique>;

}
