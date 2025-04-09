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

/**@file   miscvisualization.h
 * @brief  miscellaneous methods for visualizations
 * @author Hanna Franzen
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_MISCVISUALIZATION_H__
#define GCG_MISCVISUALIZATION_H__

#include <iostream>
#include <string>
#include <fstream>

#include "gcg/class_detprobdata.h"

using namespace gcg;

/** Gives a consistent filename for a (single) partialdec visualization that includes the probname and partialdecID.
 *
 */
void GCGgetVisualizationFilename(
   GCG* gcg,                        /**< GCG data structure */
   PARTIALDECOMP* partialdec,       /**< partialdec that is to be visualized */
   const char* extension,  /**< future file extension (to be included in the name) */
   char* filename          /**< filename output */
   );


/** Gives the path of the provided file.
 *
 */
void GCGgetFilePath(
   FILE* file,       /**< file */
   char* path        /**< buffer containing the path afterward, must be of length PATH_MAX! */
   );

#endif /* SRC_MISCVISUALIZATION_H_ */
