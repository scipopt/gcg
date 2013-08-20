/*
 * graph_boost.cpp
 *
 *  Created on: Apr 17, 2013
 *      Author: thome
 */

#include "graph_boost.h"

namespace gcg {

GraphBoost::GraphBoost()
{
   // TODO Auto-generated constructor stub

}

GraphBoost::~GraphBoost()
{
   // TODO Auto-generated destructor stub
}

int GraphBoost::getNNodes()
{
   return 0;
}

int GraphBoost::getNEdges()
{
   return 0;
}

SCIP_Bool GraphBoost::isEdge(int i, int j)
{
   return true;
}

int GraphBoost::getNNeighbors(int i)
{
   return 0;
}

std::vector<int> GraphBoost::getNeighbors(int i)
{
   std::vector<int> test;

   return test;
}

SCIP_RETCODE GraphBoost::addNode(int i, int weight)
{
   return SCIP_OKAY;
}

SCIP_RETCODE GraphBoost::addEdge(int i, int j)
{
   return SCIP_OKAY;
}

//int* GraphBoost::graphGetFirstAdjedge(int i)
//int* GraphBoost::graphGetLastAdjedge(int i)
int GraphBoost::graphGetWeights(int i)
{
   return 0;
}

SCIP_RETCODE GraphBoost::graphFlush()
{
   return SCIP_OKAY;
}

} /* namespace gcg */
