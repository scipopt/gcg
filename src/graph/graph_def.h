/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   graph_def.h
 * @brief  miscellaneous graph methods for structure detection
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_Debug

#ifndef GCG_GRAPH_DEF_H_
#define GCG_GRAPH_DEF_H_

#include "scip/scip.h"
#include "graph/graph.h"

namespace gcg {

template <class T>
Graph<T>::Graph(
   GCG*                  gcgstruct                /**< GCG data structure */
) : name("graph"),gcg(gcgstruct),graph(NULL),nconss(0),nvars(0),nnonzeroes(0),dummynodes(0)
{
   graph = new T();
}

template <class T>
Graph<T>::~Graph()
{
   if(graph != NULL)
      delete graph;

}

template <class T>
SCIP_RETCODE Graph<T>::addNNodes(int _n_nodes)
{
   return graph->addNNodes(_n_nodes);
}

template <class T>
SCIP_RETCODE Graph<T>::addNNodes(int _n_nodes, std::vector<int> weights)
{
   return graph->addNNodes(_n_nodes, weights);
}

template <class T>
int Graph<T>::getNNodes() {
   return graph->getNNodes();
}

template <class T>
int Graph<T>::getNEdges() {
   return graph->getNEdges();
}

template <class T>
SCIP_RETCODE Graph<T>::getEdges(std::vector<void*>& edges)
{
   return graph->getEdges(edges);
}

template <class T>
SCIP_RETCODE  Graph<T>::addNode(int i,int weight)
{
   SCIP_CALL( graph->addNode(i, weight) );
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE  Graph<T>::addNode()
{
   SCIP_CALL( graph->addNode() );
   return SCIP_OKAY;
}

/** adds the edge to the graph */
template <class T>
SCIP_RETCODE Graph<T>::addEdge(int i, int j)
{
   SCIP_CALL( graph->addEdge(i, j) );
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE Graph<T>::flush()
{
   SCIP_CALL( graph->flush() );
   return SCIP_OKAY;
}

template <class T>
SCIP_RETCODE Graph<T>::normalize()
{
   SCIP_CALL( graph->normalize() );
   return SCIP_OKAY;
}

template <class T>
int Graph<T>::edge(int i, int j) {
   assert( i>= 0);
   assert(j >= 0);

   int edge_ij=0;
   std::vector<int> Neighbors;

   Neighbors = getNeighbors(i);
   for(int k=0; k<(int)Neighbors.size(); k++)
   {
      if(Neighbors[k] == j)
      {
         edge_ij = 1;
         k = (int)Neighbors.size();
      }
   }
   return edge_ij;
}

template <class T>
int Graph<T>::getNNeighbors(int i) {
   assert( i >= 0);
   return graph->getNNeighbors(i);
}

template <class T>
std::vector<int> Graph<T>::getNeighbors(int i) {
   assert(i >= 0);

   return graph->getNeighbors(i);
}

template <class T>
void Graph<T>::setPartition(int i, int ID) {
   partition.resize(getNNodes(), -1);
   partition[i] = ID;
}

/** write the graph to a file */
template <class T>
SCIP_RETCODE Graph<T>::writeToFile(
      int                fd,
      SCIP_Bool writeweights
    )
{
   int nnodes;
   int nedges;
   FILE* file;
   file = fdopen(fd, "wx");
   SCIP* scip = GCGgetOrigprob(gcg);

   if( file == NULL )
      return SCIP_FILECREATEERROR;

   nnodes = Graph<T>::getNNodes();
   nedges = Graph<T>::getNEdges();

   SCIPinfoMessage(scip, file, "%d %d\n", nnodes+dummynodes, nedges/2);

   for( int i = 0; i < nnodes; ++i )
   {
      int nneighbors = Graph<T>::getNNeighbors(i);
      std::vector<int> neighbors = Graph<T>::getNeighbors(i);

      if( writeweights )
      {
         SCIPinfoMessage(scip, file, "%d ", Graph<T>::getWeight(i));
      }
      for( int j = 0; j < nneighbors; ++j )
      {
         SCIPinfoMessage(scip, file, "%d ", neighbors[j]+1);
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   for( int i = 0; i < dummynodes; ++i )
   {
      SCIPinfoMessage(scip, file, "\n");
   }

   return SCIP_OKAY;
}


/** read in the partition from a file */
template <class T>
SCIP_RETCODE Graph<T>::readPartition(
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
template <class T>
int Graph<T>::getWeight(
   int                i                   /**< the given node */
   )
{
   return graph->graphGetWeights(i);
}


/** adds the weighted edge to the graph */
template <class T>
SCIP_RETCODE Graph<T>::addEdge(int i, int j, double weight)
{
   SCIP_CALL( graph->addEdge(i, j, weight) );
   return SCIP_OKAY;
}

/** sets the weight of the edge in the graph */
template <class T>
SCIP_RETCODE Graph<T>::setEdge(int i, int j, double weight)
{
   SCIP_CALL( graph->setEdge(i, j, weight) );
   return SCIP_OKAY;
}

/** returns the weight of the edge in the graph */
template <class T>
double Graph<T>::getEdgeWeight(int i, int j)
{
   return graph->getEdgeWeight(i, j);
}

template <class T>
std::vector<std::pair<int, double> > Graph<T>::getNeighborWeights(int i)
{
   return graph->getNeighborWeights(i);
}


template <class T>
double Graph<T>::getEdgeWeightPercentile(double q)
{
   return graph->getEdgeWeightPercentile(q);
}



#ifdef WITH_GSL

template <class T>
void Graph<T>::expand(int factor)
{
   graph->expand(factor);
}

template <class T>
void Graph<T>::inflate(double factor)
{
   graph->inflate(factor);
}

template <class T>
void Graph<T>::colL1Norm()
{
   graph->colL1Norm();
}

template <class T>
void Graph<T>::prune()
{
   graph->prune();
}

template <class T>
bool Graph<T>::stopMCL(int iter)
{
   return graph->stopMCL(iter);
}

template <class T>
std::vector<int> Graph<T>::getClustersMCL()
{
   return graph->getClustersMCL();
}


template <class T>
void Graph<T>::initMCL()
{
   graph->initMCL();
}

template <class T>
void Graph<T>::clearMCL()
{
   graph->clearMCL();
}



#endif


} /* namespace gcg */

#endif
