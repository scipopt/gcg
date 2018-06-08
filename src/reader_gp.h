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

/**@file   reader_gp.h
 * @brief  GP file reader writing gnuplot files
 * @author Martin Bergner
 * @author Hanna Franzen
 * @ingroup FILEREADERS
 *
 * This reader can write visualizations of seeeds to a .gp file.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_GP_H__
#define GCG_READER_GP_H__

#include "scip/scip.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif


/** includes the gp file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderGp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** writes a visualization for the given seeed */
extern SCIP_RETCODE GCGwriteGpVisualization(
   SCIP* scip,             /**< SCIP data structure */
   char* filename,         /**< filename (including path) to write to */
   char* outputname,       /**< filename for compiled output file */
   int seeedid             /**< id of seeed to visualize */
   );

/** makes a new makefile and readme for the given .gp file */
extern SCIP_RETCODE GCGgpWriteMakefileAndReadme(
   SCIP* scip,       /**< SCIP data structure */
   FILE* file        /**< File for which the makefile & readme are generated */
   );

#ifdef __cplusplus
}
#endif

#endif
