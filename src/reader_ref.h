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
 * @brief  REF file reader for files *ref.txt
 * @ingroup FILEREADERS
 * @author Gerald Gamrath
 * @author Christian Puchert
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_REF_H__
#define __SCIP_READER_REF_H__


#include "scip/scip.h"
#include "type_decomp.h"

/** includes the ref file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderRef(
   SCIP* scip                 /**< SCIP data structure */
   );


/* reads problem from file */
extern
SCIP_RETCODE SCIPreadRef(
   SCIP*        scip,         /**< SCIP data structure */
   SCIP_READER* reader,       /**< the file reader itself */
   const char*  filename,     /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT* result        /**< pointer to store the result of the file reading call */
   );

#endif
