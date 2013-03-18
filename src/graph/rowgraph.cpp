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

/**@file   rowgraph.cpp
 * @brief  A row graph where each row is a node and rows are adjacent if they share a variable
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
#include "rowgraph.h"

namespace gcg {

RowGraph::RowGraph(
   SCIP*                 scip,              /**< SCIP data structure */
   Weights               &w                 /**< weights for the given graph */
   ) : BipartiteGraph(scip, w)
{
   // TODO Auto-generated constructor stub

}

RowGraph::~RowGraph()
{
   // TODO Auto-generated destructor stub
}

/** writes row graph to file */
SCIP_RETCODE RowGraph::writeToFile(
   const char*        filename,           /**< filename where the graph should be written to */
   SCIP_Bool          writeweights         /**< whether to write weights */
   )
{
   int nedges;
   int* nrealneighbors;
   int** realneighbors;

   SCIP_Bool* handled;
   FILE* file;
   assert(filename != NULL);
   file = fopen(filename, "wx");
   if( file == NULL )
      return SCIP_FILECREATEERROR;

   nrealneighbors = 0;
   nedges = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip_, &handled, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip_, &realneighbors, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip_, &nrealneighbors, nconss) );

   SCIPdebug(tcliquePrintGraph(tgraph));
   for( int i = 0; i < nconss; ++i )
   {
      BMSclearMemoryArray(handled, nconss);
      handled[i] = TRUE;
      nrealneighbors[i] = 0;

      SCIP_CALL( SCIPallocMemoryArray(scip_, &realneighbors[i], nconss) );
      int nneighbors = getNNeighbors(nvars+i);

      SCIPdebugMessage("%d has %d neighbors\n", i+nvars, nneighbors);

      std::vector<int> neighbors = getNeighbors(i+nvars);
      for( int j = 0; j < nneighbors; ++j )
      {
         int neighbor = neighbors[j];
         int nneighborneighbors = getNNeighbors(neighbor);

         SCIPdebugMessage("\tneighbor %d has %d neighbors\n", neighbor, nneighborneighbors);
         std::vector<int> neighborneighbors = getNeighbors(neighbor);
         for( int k = 0; k < nneighborneighbors; ++k )
         {
            int neighborneighbor = neighborneighbors[k];

            SCIPdebugMessage("\t\t%d->%d->%d (", i+nvars, neighbor, neighborneighbor);
            if( !handled[neighborneighbor-nvars] )
            {
               SCIPdebugPrintf("x)\n");
               realneighbors[i][nrealneighbors[i]] = neighborneighbor-nvars;
               ++(nrealneighbors[i]);

               handled[neighborneighbor-nvars] = TRUE;
               ++nedges;
            }
            else
            {
               SCIPdebugPrintf("-)\n");
            }
         }
      }
   }

   SCIPinfoMessage(scip_, file, "%d %d\n", nconss, nedges);

   for( int i = 0; i < nconss; ++i)
   {
      for( int j = 0; j < nrealneighbors[i]; ++j )
      {
         SCIPinfoMessage(scip_, file, "%d ", realneighbors[i][j]+1);
      }
      SCIPinfoMessage(scip_, file, "\n");
      SCIPfreeMemoryArray(scip_, &realneighbors[i]);
   }

   for( int i = 0; i < dummynodes; ++i )
   {
      SCIPinfoMessage(scip_, file, "\n");
   }

   SCIPfreeMemoryArray(scip_, &handled);
   SCIPfreeMemoryArray(scip_, &realneighbors);
   SCIPfreeMemoryArray(scip_, &nrealneighbors);

   return SCIP_OKAY;
}

} /* namespace gcg */
