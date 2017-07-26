/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2016 Operations Research, RWTH Aachen University       */
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

/**@file   reader_tex.h
 * @brief  tex file reader for writing decomposition details to LaTeX files
 * @author Hanna Franzen
 * @ingroup FILEREADERS

 * This reader can write reports of decompositions to a tex file.
 * The gp reader is required for visualizations.

 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_TEX_H__
#define GCG_READER_TEX_H__

#include "scip/scip.h"
#include "type_decomp.h"

/** includes the tex file reader into SCIP */
extern SCIP_RETCODE SCIPincludeReaderTex(
   SCIP* scip /**< SCIP data structure */
   );

/** writes a visualization for the given seeed */
extern SCIP_RETCODE GCGwriteTexVisualization(
   char* filename,         /**< filename including path */
   int seeedid,            /**< id of seeed to visualize */
   SCIP_Bool statistics    /**< additionally to picture show statistics */
   );

/** writes a visualization of the family tree of the current seeedpool */
extern SCIP_RETCODE GCGwriteTexFamilyTree(
   char* filename   /**< filename including path */
   //@todo how to get/give the seeedpool?
   );

/*@todo more params? statistics, titlepage, toc, etc? */
/*@todo is a int* of seeedids the best option? */
/** writes a report for the given seeeds */
extern SCIP_RETCODE GCGwriteTexReport(
   char* filename,   /**< filename including path */
   int* seeedids     /**< ids of seeeds to visualize */
   );

/** makes a new makefile and readme for the given .tex file */
extern SCIP_RETCODE GCGtexWriteMakefileAndReadme(SCIP* scip, /**< SCIP data structure */
   FILE* file /**< File for which the makefile & readme are generated */
);

#endif /* GCG_READER_TEX_H__ */

