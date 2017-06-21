/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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
 * @ingroup FILEREADERS
 *
 * This file reader will write the decomposed or original matrix to a file usuable by gnuplot. You can use this writer
 * like and other writer.
 *
 * To display your decomposition, you can write the following in the interactive shell
 * \verbatim
GCG> write problem "prob.gp"
\endverbatim
 * If you use this command before reading in a decomposition, you will get a picture of the original matrix. If you
 * call the writer after a decomposition was specified or detected, it will write out a picture of the structured,
 * reordered matrix  where the structure is further highlighted indicated by boxes. You can create a PDF file by then
 * calling <code>gnuplot prob.gp</code> from your systems command line. This will create a pdf file in your current
 * directory which contains the image.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_GP_H__
#define GCG_READER_GP_H__

#define COLOR_WHITE     "#FFFFFF"
#define COLOR_BLUE      "#00549F"
#define COLOR_LBLUE     "#8EBAE5"
#define COLOR_PURPLE    "#7A6FAC"
#define COLOR_VIOLET    "#612158"
#define COLOR_CARMINE   "#A11035"
#define COLOR_RED       "#CC071E"
#define COLOR_MAGENTA   "#E30066"
#define COLOR_ORANGE    "#F6A800"
#define COLOR_YELLOW    "#FFED00"
#define COLOR_GRASS     "#BDAB27"
#define COLOR_GREEN     "#57AB27"
#define COLOR_CYAN      "#0098A1"
#define COLOR_TEAL      "#006165"
#define COLOR_BLACK     "#000000"

#define COLOR_MASTERVARS   COLOR_WHITE   /* for mastervars in block area */
#define COLOR_MASTERCONS   COLOR_BLUE    /* for mastercons */
#define COLOR_LINKING      COLOR_PURPLE
#define COLOR_STAIRLINKING COLOR_MAGENTA
#define COLOR_BLOCK        COLOR_TEAL
#define COLOR_OPEN         COLOR_GREEN
#define COLOR_NONZERO      COLOR_BLACK

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

/** writes the decomposition to the specific file */
SCIP_RETCODE SCIPwriteGp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition, /**< whether to write decomposed problem */
   SCIP_Bool             outputPDF           /**< if true give pdf file, if false give tex file instead */
   );

/** Getter of parameter draftmode */
SCIP_Bool GCGgpGetDraftmode(
   SCIP*                scip               /**< SCIP data structure */
   );

/** Setter of parameter draftmode */
void GCGgpSetDraftmode(
   SCIP*                scip,              /**< SCIP data structure */
   SCIP_Bool            usedraftmode       /**< new value for draftmode */
   );

#ifdef __cplusplus
}
#endif

#endif
