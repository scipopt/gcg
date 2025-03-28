/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
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

/**@file   reader_gp.h
 * @brief  GP file reader writing decompositions to gnuplot files
 * @author Martin Bergner
 * @author Hanna Franzen
 * @ingroup FILEREADERS-GCG
 *
 * This reader can write visualizations of partialdecs to a .gp file.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_GP_H__
#define GCG_READER_GP_H__

#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Output format of gnuplot. Specifies the output format that gnuplot will produce. */
enum GPOutputFormat
{
   GP_OUTPUT_FORMAT_PDF,
   GP_OUTPUT_FORMAT_PNG,
   GP_OUTPUT_FORMAT_SVG
};
typedef enum GPOutputFormat GP_OUTPUT_FORMAT;

/** Includes the gp file reader into SCIP
 * @returns SCIP status */
GCG_EXPORT
SCIP_RETCODE GCGincludeReaderGp(
   GCG*                  gcg                 /**< GCG data structure */
   );

/* Writes a visualization for the given partialdec */
GCG_EXPORT
SCIP_RETCODE GCGwriteGpVisualizationFormat(
   GCG* gcg,               /**< GCG data structure */
   char* filename,         /**< filename (including path) to write to */
   char* outputname,       /**< filename for compiled output file */
   int partialdecid,       /**< id of partialdec to visualize */
   GP_OUTPUT_FORMAT outputformat /**< the output format which gnuplot should emit */
   );

/** Writes a visualization as .pdf file for the given partialdec
 * @returns SCIP status */
GCG_EXPORT
SCIP_RETCODE GCGwriteGpVisualization(
   GCG* gcg,               /**< GCG data structure */
   char* filename,         /**< filename (including path), location of the output*/
   char* outputname,       /**< outputname is the name of the file for the compiled gnuplot output file */
   int partialdecid             /**< id of partialdec to visualize */
   );

/** Creates a block matrix and outputs its visualization as .pdf file
 * @returns SCIP return code
 * */
GCG_EXPORT
SCIP_RETCODE GCGWriteGpDecompMatrix(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           filename,           /**< filename the output should be written to (including directory) */
   const char*           workfolder,         /**< directory in which should be worked */
   SCIP_Bool             originalmatrix      /**< should the original (or transformed) matrix be written */
);

#ifdef __cplusplus
}
#endif

#endif
