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
#include <exception>

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

protected:
   SCIP* scip_;
   TCLIQUE_GRAPH* tgraph;
   int nconss;
   int nvars;
   int nnonzeroes;
   Weights weights;
   int* partition;

public:
/** Constructor */
   Graph(
      SCIP*                 scip,              /**< SCIP data structure */
      Weights               &w                 /**< weights for the given graph */
   ) : scip_(scip),tgraph(NULL),nconss(0),nvars(0),nnonzeroes(0),weights(w),partition(NULL)
   {
     TCLIQUE_CALL_EXC( tcliqueCreate(&tgraph) );
   }

   /** Destruktor */
   virtual ~Graph()
   {
      if(tgraph != NULL)
      {
         tcliqueFree(&tgraph);
         tgraph = NULL;
      }

      SCIPfreeMemoryArrayNull(scip_, &partition);
   }

   int getNNodes() {
      return tcliqueGetNNodes(tgraph);
   }

   int getNEdges() {
      return tcliqueGetNEdges(tgraph);
   }

   int getNNeighbors(int i) {
      assert( i >= 0);
      return tcliqueGetLastAdjedge(tgraph,i)-tcliqueGetFirstAdjedge(tgraph, i)+1;
   }

   int* getNeighbours(int i) {
      assert(i >= 0);
      return tcliqueGetFirstAdjedge(tgraph, i);
   }

   int* getPartition()
   {
      return partition;
   }

   /** create graph from the matrix, to be overriden by the implementation*/
   virtual SCIP_RETCODE createFromMatrix(
      SCIP_CONS**        conss,              /**< constraints for which graph should be created */
      SCIP_VAR**         vars,               /**< variables for which graph should be created */
      int                nconss_,            /**< number of constraints */
      int                nvars_              /**< number of variables */
   ) = 0;

   /** writes the graph to the given file.
    *  The format is graph dependent
    */
   virtual SCIP_RETCODE writeToFile(
      const char*        filename            /**< filename where the graph should be written to */
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
};
}

#endif
