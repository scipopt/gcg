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

#ifndef SRC_READER_TEX_TEMP_H_
#define SRC_READER_TEX_TEMP_H_

#include "scip/scip.h"
#include "type_decomp.h"

extern "C"{

/** includes the tex file reader into SCIP */
extern SCIP_RETCODE SCIPincludeReaderTex(SCIP* scip /**< SCIP data structure */
);

/** gets the path of the file */
extern SCIP_RETCODE GCGgetFilePath(SCIP* scip, /**< SCIP data structure */
FILE* file, /**< file */
char* pfile /**< return path of file */
);

/** write LaTeX code header & begin of document
 * The proper order in which a tex file is written goes as follows:
 *    -> GCGtexWriteHeaderCode         (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
extern SCIP_RETCODE GCGtexWriteHeaderCode(SCIP* scip, /**< SCIP data structure */
FILE* file /**< File pointer to write to */
);

/** write LaTeX code title page that includes general statistics about the problem
 *  * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    -> GCGtexWriteTitlepage          (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
extern SCIP_RETCODE GCGtexWriteTitlepage(SCIP* scip, /**< SCIP data structure */
FILE* file, /**< File pointer to write to */
int* npresenteddecomps /**< Number of decompositions to be shown in the file or NULL if unknown */
);

/** write LaTeX code for table of contents
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    -> GCGtexWriteTableOfContents    (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
extern SCIP_RETCODE GCGtexWriteTableOfContents(SCIP* scip, /**< SCIP data structure */
FILE* file /**< File pointer to write to */
);

/** write LaTeX code for one decomposition
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    -> GCGtexWriteDecompCode         (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
extern SCIP_RETCODE GCGtexWriteDecompCode(SCIP* scip, /**< SCIP data structure */
FILE* file, /**< File pointer to write to */
DEC_DECOMP* decomp /**< Decomposition pointer */
);

/** write LaTeX code for end of document
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    -> GCGtexWriteEndCode            (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
extern SCIP_RETCODE GCGtexWriteEndCode(SCIP* scip, /**< SCIP data structure */
FILE* file /**< File pointer to write to */
);

/** makes a new makefile and readme for the given .tex file
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    -> GCGtexWriteMakefileAndReadme  (optional but highly recommended)
 */
extern SCIP_RETCODE GCGtexWriteMakefileAndReadme(SCIP* scip, /**< SCIP data structure */
FILE* file /**< File for which the makefile & readme are generated */
);

/** Getter of parameter usegp */
extern
SCIP_Bool GCGtexGetUseGp(SCIP* scip /**< SCIP data structure */
);

/** Getter of parameter picturesonly */
extern
SCIP_Bool GCGtexGetPicturesonly(SCIP* scip /**< SCIP data structure */
);

/** Getter of parameter draftmode */
extern
SCIP_Bool GCGtexGetDraftmode(SCIP* scip /**< SCIP data structure */
);

/*extern C*/
}

#endif /* SRC_READER_TEX_TEMP_H_ */
