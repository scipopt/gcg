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

/**@file   hyperrowcolgraph.h
 * @brief  Description
 * @author bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_HYPERROWCOLGRAPH_H_
#define GCG_HYPERROWCOLGRAPH_H_

#include "graph.h"

namespace gcg {

class HyperrowcolGraph: public Graph
{
public:
   HyperrowcolGraph(
         SCIP*                 scip,              /**< SCIP data structure */
         Weights               &w                 /**< weights for the given graph */
      );
   virtual ~HyperrowcolGraph();
   SCIP_RETCODE createFromMatrix(
      SCIP_CONS**           conss,              /**< constraints for which graph should be created */
      SCIP_VAR**            vars,               /**< variables for which graph should be created */
      int                   nconss,             /**< number of constraints */
      int                   nvars               /**< number of variables */
      );

   virtual SCIP_RETCODE writeToFile(
      const char* filename
    );
   SCIP_RETCODE readPartition(
      const char* filename
   );

   virtual int getNNodes();
   virtual int getNEdges();


   std::vector<int> getNeighbors(
         int i
      );
   std::vector<int> getHyperedgeNodes(
         int i
      );
};

} /* namespace gcg */
#endif /* GCG_HYPERROWCOLGRAPH_H_ */
