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

/**@file   class_miscvisualization.h
 * @brief  miscellaneous methods for visualizations
 * @author Hanna Franzen
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_MISCVISUALIZATION_H__
#define GCG_CLASS_MISCVISUALIZATION_H__

#include <iostream>
#include <string>
#include <fstream>

#include "class_seeedpool.h"

using namespace gcg;

class MiscVisualization
{

private:

public:

   /** constructor */
   MiscVisualization();


   /** destructor */
   ~MiscVisualization();


   /** gives a consistent filename for a (single) seeed visualization that includes the probname and seeedID
    *
    * @return standardized filename
    * */
   SCIP_RETCODE GCGgetVisualizationFilename(
      SCIP* scip,             /**< scip data structure */
      SeeedPtr seeed,         /**< seeed that is to be visualized */
      const char* extension,  /**< file extension (to be included in the name) */
      char* filename          /**< filename output */
      );


   /** gives the path of the file
    *
    * @return path of file
    * */
   char* GCGgetFilePath(
      SCIP* scip,       /**< scip data structure */
      FILE* file        /**< file */
      );


   /** compiles visualization files in gp or tex format and
    * opens the resulting pdf file with the
    *
    * @return path of file
    * */
   void GCGshowVisualization(
      SCIP* scip,       /**< scip data structure */
      char* filename    /**< filename (including path) */
      );


   /** checks in which seeedpool the seeed with given ID is stored and returns that seeedpool
    *
    * @returns pool: Seeedpool* where the Seeed was found
    */
   Seeedpool* GCGgetSeeedpoolForSeeed(
      SCIP* scip,       /**< SCIP data structure */
      int seeedid       /**< ID of Seeed */
      );

}; /* class MiscVisualization */

#endif /* SRC_CLASS_MISCVISUALIZATION_H_ */
