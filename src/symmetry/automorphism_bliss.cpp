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

#include "symmetry/pub_automorphism.h"
#include "bliss/graph.hh"
#include "scip/scip.h"

struct struct_graph_data
{
   SCIP* scip_;
   bliss::Graph* graph;
   bliss::Stats bstats;
   bool terminate;
};

SCIP_RETCODE struct_graph::init(
   SCIP* scip,
   int nvertices
   )
{
   graphdata = new struct_graph_data();
   graphdata->graph = new bliss::Graph(nvertices);
   assert(graphdata->graph->get_nof_vertices() == (unsigned int)nvertices);
   //graph = new bliss::Graph();
   graphdata->scip_ = scip;
   graphdata->terminate = false;
   return SCIP_OKAY;
}

SCIP_RETCODE struct_graph::destroy(
   )
{
   delete graphdata->graph;
   delete graphdata;
   return SCIP_OKAY;
}

void struct_graph::setColor(
   int vertex,
   int color
   )
{
   assert(vertex < graphdata->graph->get_nof_vertices());
   graphdata->graph->change_color((unsigned int)vertex, (unsigned int)color);
   //graph->add_vertex((unsigned int)color);
}

void struct_graph::addEdge(
   int v1,
   int v2
   )
{
   graphdata->graph->add_edge((unsigned int)v1, (unsigned int)v2);
}

unsigned int struct_graph::getNVertices(
   )
{
   return graphdata->graph->get_nof_vertices();
}

SCIP_RETCODE struct_graph::findAutomorphisms(
   void* userdata,
   void (*fhook)(void*, unsigned int, const unsigned int*),
   unsigned int searchnodelimit,
   unsigned int generatorlimit
   )
{
#ifdef BLISS_PATCH_PRESENT
   if( searchnodelimit > 0u || generatorlimit > 0u)
      graphdata->graph->set_search_limits(searchnodelimit, generatorlimit);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   auto report = [&](unsigned int n, const unsigned int* aut) {
      (*fhook)((void*)userdata, n, aut);
   };

   auto term = [&]() {
      return graphdata->terminate
            || (generatorlimit > 0 && graphdata->bstats.get_nof_generators() >= (long unsigned int) generatorlimit)
            || (searchnodelimit > 0 && graphdata->bstats.get_nof_nodes() >= (long unsigned int) searchnodelimit);
   };

   graphdata->graph->find_automorphisms(graphdata->bstats, report, term);
#else
   graphdata->graph->find_automorphisms(graphdata->bstats, fhook, userdata);
#endif
   return SCIP_OKAY;
}

void struct_graph::terminateSearch()
{
   graphdata->terminate = true;
}

extern "C"
void GCGgetBlissName(char* buffer, int len)
{
#ifdef BLISS_PATCH_PRESENT
   SCIPsnprintf(buffer, len, "bliss %sp", bliss::version);
#else
   SCIPsnprintf(buffer, len, "bliss %s", bliss::version);
#endif
}
