/*
 * graph_tclique.h
 *
 *  Created on: Apr 16, 2013
 *      Author: thome
 */

#ifndef GRAPH_TCLIQUE_H_
#define GRAPH_TCLIQUE_H_

#include "bridge.h"
#include "tclique/tclique.h"

namespace gcg {

class GraphTclique: public gcg::Bridge
{
private:
   TCLIQUE_GRAPH* graph;

public:
   GraphTclique();
   virtual ~GraphTclique();
   virtual int getNNodes();
   virtual int getNEdges();
   virtual SCIP_Bool isEdge(int i, int j);
   virtual int getNNeighbors(int i);
   virtual std::vector<int> getNeighbors(int i);
   virtual SCIP_RETCODE addNode(int i, int weight);
   virtual SCIP_RETCODE deleteNode(int i);
   virtual SCIP_RETCODE addEdge(int i, int j);
   virtual SCIP_RETCODE deleteEdge(int i, int j);
   virtual int graphGetWeights(int i);

   virtual SCIP_RETCODE flush();
};














} /* namespace gcg */
#endif /* GRAPH_TCLIQUE_H_ */
