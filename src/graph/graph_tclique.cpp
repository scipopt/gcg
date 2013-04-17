/*
 * graph_tclique.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: thome
 */

#include <cassert>
#include "graph_tclique.h"

namespace gcg {

GraphTclique::GraphTclique()
{
   tcliqueCreate(&graph);
   }

GraphTclique::~GraphTclique()
{
   tcliqueFree(&graph);
}

int GraphTclique::getNNodes()
{
   return tcliqueGetNNodes(graph);
}

int GraphTclique::getNEdges()
{
   return tcliqueGetNEdges(graph);
}

SCIP_Bool GraphTclique::isEdge(int i, int j)
{
   assert(i >= 0);
   assert(j >= 0);

   return tcliqueIsEdge(graph, i, j);
}

int GraphTclique::getNNeighbors(int i)
{
   assert( i >= 0);
   return tcliqueGetLastAdjedge(graph,i)-tcliqueGetFirstAdjedge(graph, i)+1;
}

std::vector<int> GraphTclique::getNeighbors(int i)
{
   assert(i >= 0);
   std::vector<int> part(tcliqueGetFirstAdjedge(graph, i), tcliqueGetLastAdjedge(graph,i)+1);
   return part;
}

SCIP_RETCODE GraphTclique::addNode(int i, int weight)
{
   assert(i >= getNNodes());
   if(tcliqueAddNode(graph,i,weight))
      return SCIP_OKAY;
   else
      return SCIP_ERROR;
}

SCIP_RETCODE GraphTclique::deleteNode(int i)
{
   return SCIP_ERROR;
}

SCIP_RETCODE GraphTclique::addEdge(int i, int j)
{
   assert(i >=0);
   assert(i < getNNodes());
   assert(j >=0);
   assert(j < getNNodes());

   if(tcliqueAddEdge(graph,i,j))
      return SCIP_OKAY;
   else
      return SCIP_ERROR;
}

SCIP_RETCODE GraphTclique::deleteEdge(int i, int j)
{
   return SCIP_ERROR;
}

int* GraphTclique::graphGetFirstAdjedge(int i)
{
   return tcliqueGetFirstAdjedge(graph,i);
}

int* GraphTclique::graphGetLastAdjedge(int i)
{
   return tcliqueGetLastAdjedge(graph,i);
}

int GraphTclique::graphGetWeights(int i)
{
   return tcliqueGetWeights(graph)[i];
}

SCIP_RETCODE GraphTclique::graphFlush()
{
   if(tcliqueFlush(graph))
      return SCIP_OKAY;
   else
      return SCIP_ERROR;
}

} /* namespace gcg */
