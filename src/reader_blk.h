/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_blk.h
 * @brief  BLK file reader
 * @author Gerald Gamrath
 * @ingroup FILEREADERS
 *
 * This reader reads in a blk-file that defines the structur to be used for the decomposition.
 * The structure is defined variable-wise, i.e., the number of blocks and the variables belonging to each block are
 * defined. Afterwards, each constraint that has only variables of one block is added to that block,
 * constraints having variables of more than one block go into the master. If needed, constraints can also be
 * forced into the master, even if they could be transferred to one block.
 *
 * The keywords are:
 * - Presolved: to be followed by either 0 or 1 indicating that the decomposition is for the presolved or unpresolved problem
 * - NBlocks: to be followed by a line giving the number of blocks
 * - Block i with 1 <= i <= nblocks: to be followed by the names of the variables belonging to block i, one per line.
 * - Masterconss: to be followed by names of constraints, one per line, that should go into the master,
 *                even if they only contain variables of one block and could thus be added to this block.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_BLK_H__
#define __SCIP_READER_BLK_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the blk file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderBlk(
   SCIP*                 scip                /**< SCIP data structure */
   );


/* reads problem from file */
extern
SCIP_RETCODE SCIPreadBlk(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   );

#ifdef __cplusplus
}
#endif

#endif
