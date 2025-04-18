/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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

/**@file   hyperrowcolgraph.h
 * @brief  A hypergraph with row and column nodes
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_HYPERROWCOLGRAPH_H_
#define GCG_HYPERROWCOLGRAPH_H_

#include "graph/matrixgraph.h"
#include "graph/hypergraph.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"

namespace gcg {
template <class T>
class HyperrowcolGraph: public MatrixGraph<T>
{
private:
   Graph<T> graph;
public:
   HyperrowcolGraph(
         GCG*                  gcg,               /**< GCG data structure */
         Weights               w                  /**< weights for the given graph */
      );
   virtual ~HyperrowcolGraph();

   SCIP_RETCODE createFromMatrix(
      SCIP_CONS**           conss,              /**< constraints for which graph should be created */
      SCIP_VAR**            vars,               /**< variables for which graph should be created */
      int                   nconss,             /**< number of constraints */
      int                   nvars               /**< number of variables */
      );

   /** creates a graph with open constraints and open variables of the partialdec */
   virtual SCIP_RETCODE createFromPartialMatrix(
      DETPROBDATA*          detprobdata,        /**< detection process information and data */
      PARTIALDECOMP*        partialdec          /**< partial decomposition to use for matrix */
      );

   /** writes the graph to the given file.
    *  The format is graph dependent
    */
   virtual SCIP_RETCODE writeToFile(
      int                fd,                  /**< file descriptor where the graph should be written to */
      SCIP_Bool          writeweights         /**< whether to write weights */
      );


   virtual SCIP_RETCODE createDecompFromPartition(
      GCG_DECOMP**       decomp              /**< decomposition structure to generate */
      );

   /** creates a new partialdec by dint of a graph created with all constraints and variables */
   virtual SCIP_RETCODE createPartialdecFromPartition(
      PARTIALDECOMP**      firstpartialdec,         /**< pointer to buffer the new partialdec created by dint of the graph (can be NULL) */
      PARTIALDECOMP**      secondpartialdec,        /**< pointer to buffer the new partialdec whose border is amplified by dint of the graph (can be NULL) */
      DETPROBDATA*         detprobdata              /**< detection process information and data */
      );

   /** amplifies a partialdec by dint of a graph created with open constraints and open variables of the partialdec */
   virtual SCIP_RETCODE createPartialdecFromPartition(
      PARTIALDECOMP*      oldpartialdec,            /**< partialdec which should be amplifies */
      PARTIALDECOMP**     firstpartialdec,          /**< pointer to buffer the new partialdec amplified by dint of the graph (can be NULL) */
      PARTIALDECOMP**     secondpartialdec,         /**< pinter to buffer the new partialdec whose border is amplified by dint of the graph (can be NULL) */
      DETPROBDATA*        detprobdata               /**< detection process information and data */
      );

   /**
    * reads the partition from the given file.
    * The format is graph dependent. The default is a file with one line for each node a
    */
   virtual SCIP_RETCODE readPartition(
      const char*        filename            /**< filename where the partition is stored */
      );

   virtual std::vector<int> getNeighbors(
      int i
      );

   virtual std::vector<int> getHyperedgeNodes(
      int i
      );

   std::vector<int> getConsNonzeroNodes(
      int i
      );

   std::vector<int> getVarNonzeroNodes(
      int i
      );

};

} /* namespace gcg */
#endif /* GCG_HYPERROWCOLGRAPH_H_ */
