/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

/**@file   graph_boost.cpp
 * @brief  interface to boost graph library
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
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

SCIP_RETCODE GraphBoost::flush()
{
   return SCIP_OKAY;
}

} /* namespace gcg */
