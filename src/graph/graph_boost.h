/*
 * graph_boost.h
 *
 *  Created on: Apr 17, 2013
 *      Author: thome
 */

#ifndef GRAPH_BOOST_H_
#define GRAPH_BOOST_H_

#include "bridge.h"
#include "boost/graph/adjacency_list.hpp"

namespace gcg
{
class GraphBoost: public gcg::Bridge
{
public:
   GraphBoost();
   virtual ~GraphBoost();
   virtual int getNNodes();
   virtual int getNEdges();
   virtual SCIP_Bool isEdge(int i, int j);
   virtual int getNNeighbors(int i);
   virtual std::vector<int> getNeighbors(int i);
   virtual SCIP_RETCODE addNode(int i, int weight);
   virtual SCIP_RETCODE deleteNode(int i);
   virtual SCIP_RETCODE addEdge(int i, int j);
   virtual SCIP_RETCODE deleteEdge(int i, int j);
   virtual int* graphGetFirstAdjedge(int i);
   virtual int* graphGetLastAdjedge(int i);
   virtual int graphGetWeights(int i);

   virtual SCIP_RETCODE flush();
};

} /* namespace gcg */
#endif /* GRAPH_BOOST_H_ */
