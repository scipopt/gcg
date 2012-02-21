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
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_BLK_H__
#define __SCIP_READER_BLK_H__


#include "scip/scip.h"


/** includes the blk file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderBlk(
   SCIP*                 scip                /**< SCIP data structure */
   );


/* reads problem from file */
extern
SCIP_RETCODE SCIPreadBlk(
   SCIP*              scip,               /**< SCIP data structure */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   );


#endif
