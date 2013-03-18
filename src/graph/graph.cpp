/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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

/**@file   graph.cpp
 * @brief  miscellaneous graph methods for structure detection
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "graph.h"
#include "tclique/tclique.h"
#include <fstream>

using std::ifstream;

namespace gcg {

Graph::Graph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               &w                 /**< weights for the given graph */
) : scip_(scip),tgraph(NULL),nconss(0),nvars(0),nnonzeroes(0),dummynodes(0),weights(w)
{
  TCLIQUE_CALL_EXC( tcliqueCreate(&tgraph) );
}

/** Destruktor */
Graph::~Graph()
{
   if(tgraph != NULL)
   {
      tcliqueFree(&tgraph);
      tgraph = NULL;
   }
}

int Graph::getNNodes() {
   return tcliqueGetNNodes(tgraph);
}

int Graph::getNEdges() {
   return tcliqueGetNEdges(tgraph);
}

int Graph::getNNeighbors(int i) {
   assert( i >= 0);
   return tcliqueGetLastAdjedge(tgraph,i)-tcliqueGetFirstAdjedge(tgraph, i)+1;
}

std::vector<int> Graph::getNeighbors(int i) {
   assert(i >= 0);
   std::vector<int> part(tcliqueGetFirstAdjedge(tgraph, i), tcliqueGetLastAdjedge(tgraph,i)+1);
   return part;
}

std::vector<int> Graph::getPartition()
{
   return partition;
}

SCIP_RETCODE Graph::writeToFile(
      const char* filename,
      SCIP_Bool writeweights
    )
{
   int nnodes;
   int nedges;
   FILE* file;
   assert(filename != NULL);
   file = fopen(filename, "wx");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   nnodes = Graph::getNNodes();
   nedges = Graph::getNEdges();

   SCIPinfoMessage(scip_, file, "%d %d\n", nnodes+dummynodes, nedges/2);

   for( int i = 0; i < nnodes; ++i )
   {
      int nneighbors = Graph::getNNeighbors(i);
      std::vector<int> neighbors = Graph::getNeighbors(i);

      if( writeweights )
      {
         SCIPinfoMessage(scip_, file, "%d ", Graph::getWeight(i));
      }
      for( int j = 0; j < nneighbors; ++j )
      {
         SCIPinfoMessage(scip_, file, "%d ", neighbors[j]+1);
      }
      SCIPinfoMessage(scip_, file, "\n");
   }

   for( int i = 0; i < dummynodes; ++i )
   {
      SCIPinfoMessage(scip_, file, "\n");
   }

   return SCIP_OKAY;
}

SCIP_RETCODE Graph::readPartition(
   const char* filename
)
{
   ifstream input(filename);
   if( !input.good() )
   {
      SCIPerrorMessage("Could not open file <%s> for reading\n", filename);
      return SCIP_READERROR;
   }
   partition.resize(getNNodes(), -1);
   for( int i = 0; i < getNNodes(); ++i )
   {
      int part = 0;
      if( !(input >> part) )
      {
         SCIPerrorMessage("Could not read from file <%s>. It may be in the wrong format\n", filename);
         return SCIP_READERROR;
      }
      partition[i] = part;
   }

   input.close();
   return SCIP_OKAY;
}

/** return the weight of given node */
int Graph::getWeight(
   int                i                   /**< the given node */
   )
{
   return tcliqueGetWeights(tgraph)[i];
}

}
