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

/**@file   graph_interface.h
 * @brief  miscellaneous graph interface methods
 * @author Martin Bergner
 * @author Annika Thome
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef GCG_GRAPHINTERFACE_H_
#define GCG_GRAPHINTERFACE_H_

#include "objscip/objscip.h"
#include "graph/weights.h"

#include <vector>
#include <fstream>

using std::ifstream;
namespace gcg {

class GraphInterface {
protected:
   std::vector<int> partition;
public:

   GraphInterface() {}

   virtual ~GraphInterface() {}

   /** return a partition of the nodes */
   virtual std::vector<int> getPartition() const { return partition; }

   /** assigns partition to a given node*/
   virtual void setPartition(int i, int nodeid) = 0;

   /** writes the graph to the given file.
    *  The format is graph dependent
    */
   virtual SCIP_RETCODE writeToFile(
      int                fd,                  /**< filename where the graph should be written to */
      SCIP_Bool          writeweights        /**< whether to write weights */
    ) = 0;


   /**
    * reads the partition from the given file.
    * The format is graph dependent. The default is a file with one line for each node a
    */
   virtual SCIP_RETCODE readPartition(
      const char*        filename            /**< filename where the partition is stored */
   ) = 0;


   /** create decomposition based on the read in partition */
   virtual SCIP_RETCODE createDecompFromPartition(
      GCG_DECOMP**       decomp              /**< decomposition structure to generate */
   )
   { /*lint -e715*/
      return SCIP_ERROR;
   }

   virtual SCIP_RETCODE flush() = 0;

};

}

#endif
