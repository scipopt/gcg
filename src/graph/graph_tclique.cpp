/*
 * graph_tclique.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: thome
 */

#include <cassert>
#include "graph_tclique.h"

#define TCLIQUE_CALL_EXC(x)   do                                                                              \
                       {                                                                                      \
                          SCIP_Bool _restat_;                                                                 \
                          if( (_restat_ = (x)) != TRUE )                                                      \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                          \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

#define TCLIQUE_CALL(x)   do                                                                                  \
                       {                                                                                      \
                          SCIP_Bool _restat_;                                                                 \
                          if( (_restat_ = (x)) != TRUE )                                                      \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             return SCIP_ERROR;                                                               \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )



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


SCIP_RETCODE GraphTclique::flush()
{
   if(tcliqueFlush(graph))
      return SCIP_OKAY;
   else
      return SCIP_ERROR;
}

int GraphTclique::graphGetWeights(int i)
{
   assert( i >= 0);
   assert( i <= getNNodes());
   const TCLIQUE_WEIGHT* weights;
   weights = tcliqueGetWeights(graph);
   return weights[i];
}
} /* namespace gcg */
