/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_dec.h
 * @brief  DEC file reader
 * @author Lukas Kirchhart
 * @ingroup FILEREADERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_REF_N_H__
#define __SCIP_READER_REF_N_H__


#include "scip/scip.h"
#include "type_decomp.h"


/** includes the dec file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderDec(
   SCIP*                 scip             /**< SCIP data structure */
   );


/* reads problem from file */
extern
SCIP_RETCODE SCIPreadDec(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   );


/** write a DEC file for a given decomposition */
SCIP_RETCODE SCIPwriteDecomp(
   SCIP*      scip,              /**< SCIP data structure */
   FILE*      file,              /**< File pointer to write to */
   DECDECOMP* decdecomp,         /**< Decomposition pointer */
   SCIP_Bool  writeDecomposition /**< whether to write decomposed problem */
   );

#endif
