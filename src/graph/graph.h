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

/**@file   graph.h
 * @brief  miscellaneous graph methods for structure detection
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef GCG_GRAPH_H_
#define GCG_GRAPH_H_
#include "objscip/objscip.h"
#include "tclique/tclique.h"
#include "weights.h"
#include "pub_decomp.h"

#include <exception>
#include <vector>
#include <string>

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

class Graph {
public:
   std::string name;
protected:
   SCIP* scip_;
   TCLIQUE_GRAPH* tgraph;
   int nconss;
   int nvars;
   int nnonzeroes;
   int dummynodes;
   Weights weights;
   std::vector<int> partition;

public:
   /** Constructor */
   Graph(
      SCIP*                 scip,              /**< SCIP data structure */
      Weights               w                  /**< weights for the given graph */
   );

   void swap(Graph & other) // the swap member function (should never fail!)
   {
      // swap all the members (and base subobject, if applicable) with other
      std::swap(partition, other.partition);
      std::swap(scip_ , other.scip_);
      std::swap(tgraph , other.tgraph);
      std::swap(nconss , other.nconss);
      std::swap(nvars , other.nvars);
      std::swap(nnonzeroes , other.nnonzeroes);
      std::swap(weights , other.weights);
      std::swap(dummynodes, other.dummynodes);

   }

   Graph& operator=(Graph other) // note: argument passed by value!
   {
      // swap this with other
      swap(other);

      return *this;
   }

   /** Destruktor */
   virtual ~Graph();

   /** return the number of nodes */
   virtual int getNNodes();

   /** return the number of edges (or hyperedges) */
   virtual int getNEdges();

   /** return the number of neighbor nodes of given node */
   virtual int getNNeighbors(
      int                i                   /**< the given node */
      );

   /** return the neighboring nodes of a given node */
   virtual std::vector<int> getNeighbors(
      int                i                   /**< the given node */
      );

   /** return a partition of the nodes */
   std::vector<int> getPartition();

   /** create graph from the matrix, to be overriden by the implementation*/
   virtual SCIP_RETCODE createFromMatrix(
      SCIP_CONS**        conss,              /**< constraints for which graph should be created */
      SCIP_VAR**         vars,               /**< variables for which graph should be created */
      int                nconss_,            /**< number of constraints */
      int                nvars_              /**< number of variables */
   ) { return SCIP_ERROR; };

   /** writes the graph to the given file.
    *  The format is graph dependent
    */
   virtual SCIP_RETCODE writeToFile(
      const char*        filename,           /**< filename where the graph should be written to */
      SCIP_Bool          writeweights = FALSE /**< whether to write weights */
    );

   /**
    * reads the partition from the given file.
    * The format is graph dependent. The default is a file with one line for each node a
    */
   virtual SCIP_RETCODE readPartition(
      const char*        filename            /**< filename where the partition is stored */
   );

   int getNNonzeroes() const
   {
      return nnonzeroes;
   }

   /** return the weight of given node */
   virtual int getWeight(
      int                i                   /**< the given node */
      );

   /** set the number of dummy nodes */
   void setDummynodes(int dummynodes_)
   {
      dummynodes = dummynodes_;
   };

   /** create decomposition based on the read in partition */
   virtual SCIP_RETCODE createDecompFromPartition(
      DEC_DECOMP**       decomp              /**< decomposition structure to generate */
      ) { return SCIP_ERROR;};

};

}

#endif
