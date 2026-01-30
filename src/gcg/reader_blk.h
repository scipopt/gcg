/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_blk.h
 * @brief  BLK file reader for structure information
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @ingroup FILEREADERS-GCG
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

#ifndef GCG_READER_BLK_H__
#define GCG_READER_BLK_H__


#include "scip/scip.h"
#include "gcg/gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/** includes the blk file reader into SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeReaderBlk(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
