/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_ref_n.h
 * @brief  REF file reader
 * @author Lukas Kirchhart
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_REF_N_H__
#define __SCIP_READER_REF_N_H__


#include "scip/scip.h"
#include "type_decomp.h"


/** includes the ref file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderDec(
   SCIP*                 scip                /**< SCIP data structure */
   );


/* reads problem from file */
extern
SCIP_RETCODE SCIPreadDec(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   );


/* writes problem to file */
extern
SCIP_RETCODE SCIPwriteDec(
   SCIP*              scip,               /**< SCIP data structure */
   FILE*              file,               /**< output file, or NULL if standard output should be used */
   const char*        name,               /**< problem name */
   SCIP_Bool          transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE      objsense,           /**< objective sense */
   SCIP_Real          objscale,           /**< scalar applied to objective function; external objective value is
                                               extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real          objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**         vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                nvars,              /**< number of mutable variables in the problem */
   int                nbinvars,           /**< number of binary variables */
   int                nintvars,           /**< number of general integer variables */
   int                nimplvars,          /**< number of implicit integer variables */
   int                ncontvars,          /**< number of continuous variables */
   SCIP_CONS**        conss,              /**< array with constraints of the problem */
   int                nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*       result              /**< pointer to store the result of the file writing call */
   );

SCIP_RETCODE SCIPwriteDecomp(
   SCIP* scip,                                /**< SCIP data structure */
   FILE* file,                                /**< File pointer to write to */
   DECDECOMP* decdecomp,                      /**< Decomposition pointer */
   SCIP_Bool writeDecomposition                  /**< whether to write decomposed problem */
   );


SCIP_RETCODE SCIPReaderDecSetDecomp(
   SCIP* scip,
   DECDECOMP* decdecomp
   );

#endif
