/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_ref.h
 * @brief  REF file reader for files .ref
 * @ingroup FILEREADERS
 * @author Gerald Gamrath
 * @author Christian Puchert
 * @author Martin Bergner
 *
 * This reader reads and writes a ref-file that defines the structur to be used for the decomposition.
 * The structure is defined constraint-wise, i.e., the number of blocks and the constraints belonging to each
 * block are defined. The constraints are numbered by the appearance in the problem.
 *
 * Constraints not mentioned in one of the blocks will remain in the master problem
 *
 * The format is the following
 * - first line: \#nblocks \#ncons_block_1 ... \#n_cons_block_n
 * - one line for each block with the indices of constraints to be put into that block separated by a comma
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_REF_H__
#define __SCIP_READER_REF_H__


#include "scip/type_scip.h"
#include "scip/type_reader.h"
#include "scip/type_result.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the ref file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderRef(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* reads problem from file */
extern
SCIP_RETCODE SCIPreadRef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   );

#ifdef __cplusplus
}
#endif

#endif
