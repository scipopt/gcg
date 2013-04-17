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

/**@file   bridge.h
 * @brief  bridge
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef BRIDGE_H
#define BRIDGE_H
#include "objscip/objscip.h"
#include <vector>

namespace gcg
{

class Bridge
{


public:
   virtual ~Bridge() {};
   virtual int getNNodes() = 0;
   virtual int getNEdges() = 0;
   virtual SCIP_Bool isEdge(int i, int j);
   virtual int getNNeighbors(int i) = 0;
   virtual std::vector<int> getNeighbors(int i);
   virtual SCIP_RETCODE addNode(int i,int weight);
   virtual SCIP_RETCODE deleteNode(int i);
   virtual SCIP_RETCODE addEdge(int i, int j);
   virtual SCIP_RETCODE deleteEdge(int i, int j);
   virtual int* graphGetFirstAdjedge(int i);
   virtual int* graphGetLastAdjedge(int i);
   virtual int graphGetWeights(int i);


   virtual SCIP_RETCODE graphFlush();

};


}

#endif
