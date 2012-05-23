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
 * @author Martin Bergner
 * @author Lukas Kirchhart
 * @ingroup FILEREADERS

 * This reader reads in a dec-file that defines the structur to be used for the decomposition.
 * The structure is defined constraint-wise, i.e., the number of blocks and the constraints belonging
 * to each block are  defined.  If needed, constraints can also be  forced into the master, even if
 * they could be transferred to one block.
 *
 * The keywords are:
 * - NBlocks: to be followed by a line giving the number of blocks
 * - Block i with 1 <= i <= nblocks: to be followed by the names of the constraints belonging to block i,
                  one per line.
 * - Masterconss: to be followed by names of constraints, one per line, that should go into the master,
 *                even if they only contain variables of one block and could thus be added to this block.

 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_DEC_H__
#define __SCIP_READER_DEC_H__


#include "scip/scip.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the dec file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderDec(
   SCIP*                 scip             /**< SCIP data structure */
   );


/* reads problem from file */
extern
SCIP_RETCODE SCIPreadDec(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   );


/** write a DEC file for a given decomposition */
SCIP_RETCODE SCIPwriteDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DECDECOMP*            decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition  /**< whether to write decomposed problem */
   );

#ifdef __cplusplus
}
#endif

#endif
