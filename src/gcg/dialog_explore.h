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

/**@file   dialog_explore.h
 * @ingroup DIALOGS
 * @brief  dialog menu for exploring decompositions
 * @author Michael Bastubbe
 * @author Hanna Franzen
 *
 * This file contains all dialog calls to build and use the explore menu.
 * The explore menu gives the user detailed information about all decompositions and a possibility to edit such.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_DIALOG_EXPLORE_H__
#define GCG_DIALOG_EXPLORE_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief method to handle user input for "explore" command
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGdialogExecExplore(
   GCG*                    gcg,                /**< GCG data structure */
   SCIP_DIALOGHDLR*        dialoghdlr,         /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog              /**< dialog for user input management */
   );


#ifdef __cplusplus
}
#endif

#endif
