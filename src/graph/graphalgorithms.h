/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file    graphalgorithms.h
 * @brief   several metrics for graphs
 * @author  Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GRAPHALGORITHMS_H_
#define GRAPHALGORITHMS_H_

#include "hypergraph.h"

namespace gcg {

template<class T>
class GraphAlgorithms {

public:
   /** compute weighted sum of external degrees */
   static SCIP_Real computeSoed(
      Hypergraph<T>&     graph               /**< the hypergraph */
   );

   /** compute minimum hyperedge cut */
   static SCIP_Real computeMincut(
      Hypergraph<T>&     graph               /**< the hypergraph */
   );

   /** compute k-1 metric */
   static SCIP_Real computekMetric(
      Hypergraph<T>&     graph               /**< the hypergraph */
   );
};

} /* namespace gcg */



#endif /* GRAPHALGORITHMS_H_ */
