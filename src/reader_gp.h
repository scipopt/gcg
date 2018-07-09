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
 * @author Michael Feldmeier
 * @ingroup FILEREADERS
 *
 * This file reader will write the decomposed or original matrix to a file usuable by gnuplot. You can use this writer
 * like and other writer. The filereader is also able to generate a gnuplot plot of a quadratic programming instance.
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
extern
SCIP_RETCODE SCIPwriteGp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition  /**< whether to write decomposed problem */
   );

/** writes gnuplot code to the specific file to plot QP instance*/
extern
SCIP_RETCODE SCIPwriteQpGp(
   SCIP*                 scip,                     /** SCIP data structure */
   FILE*                 file,                     /** File pointer to write to */
   SCIP_CONS**           linearconstraints,        /** Array of linear constraints */
   int                   nLinearconstraints,       /** Number of linear constraints */
   SCIP_CONS**           quadraticconstraints,     /** Array of quadratic constraints */
   int                   nQuadraticconstraints,    /** Number of quadratic constraints */
   int*                  variableposition,         /** Array indicating at which position the variable at [i] is to be plot */
   int                   ncariables,               /** Number of variables */
   SCIP_Bool             writereorderedVariables   /** Use standard order of the variables or reposition as indicated in variablePosition? */
);

/** Compiles a given gnuplot file */
extern
SCIP_RETCODE GCGcompileGpFile(
   SCIP*                  scip,          /**< SCIP data structure */
   char*                  filename       /** path to .gp file */
);


#ifdef __cplusplus
}
#endif

#endif
