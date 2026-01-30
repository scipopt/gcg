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

/**@file   reader_dec.h
 * @brief  DEC file reader for structure information
 * @author Martin Bergner
 * @author Lukas Kirchhart
 * @author Micahel Bastubbe
 * @ingroup FILEREADERS-GCG

 * This reader reads and write files in .dec format. A data format to pass a (possibly partial) decomposition to GCG, prerequisite is a given MIP, whose constraints and variables are referred to by name
 * – everything behind a backslash (“\”) is a comment and is ignored
 * – information is given section-wise
 * – sections are started by key words (ignoring the case of the characters) and finished by starting a new section or reaching end of file
 * – each line in a section provides one value
 * – key words for sections are:
 *   – consdefaultmaster:
 *     – optional; followed by line with possible values: {0, 1}; default: 1; description: if set to 1 then (directly after file is read) each unassigned constraint is assigned to the master (needed for backward compatibility)
 * – presolved:
 *   – mandatory; followed by line with possible values: {0, 1}; description: if set to 0 (1) then the decomposition is considered for the unpresolved (presolved) problem
 * – nblocks
 *   – mandatory; possible values: N; description: number of (possibly empty) blocks this decomposition file has information for
 * – block (alternatives: blockconss or blockcons)
 *   – optional; directly followed by block index (starting with 1); each following line contains name of a constraint belonging to this block
 * – masterconss (alternative: mastercons)
 *   + optional; each following line contains name of a constraint belonging to the master
 * – blockvars
 *   + optional; directly followed by block index (starting with 1); each following line contains name of a variable belonging to this block
 * – mastervars (alternative: mastervar)
 *   + optional; each following line contains name of a master variable; (belongs explicitly only to master constraints)
 * – linkingvars (alternative: linkingvar)
 *   + optional; each following line contains name of a linking variable
 * – decomposition is rejected if there are any inconsistencies
 * – after reading (and and possibly assigning unassigned constraints because of consdefaultmaster, see above) implicit assignments are made:
 *   – unassigned constraints hitting at least two blocks -> assign to master;
 *   – unassigned variables hitting at least two blocks -> assign to linking ;
 * – all constraints of an unassigned variable are master constraints -> variable is master variable;
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_DEC_H__
#define GCG_READER_DEC_H__


#include "scip/scip.h"
#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the dec file reader into SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeReaderDec(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
