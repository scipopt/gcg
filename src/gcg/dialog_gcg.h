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

/**@file   dialog_gcg.h
 * @ingroup DIALOGS
 * @brief  GCG user interface dialog
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_DIALOG_GCG_H__
#define GCG_DIALOG_GCG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** dialog execution method for the display additionalstatistics command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayAdditionalStatistics);

/** dialog execution method for the display statistics command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayStatistics);

/** dialog execution method print complete detection information */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecPrintDetectionInformation);

/** dialog execution method for adding block number candidate */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecChangeAddBlocknr);

/** dialog execution method for the display detectors command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayDetectors);

/** dialog execution method for the display constraint classifiers command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayConsClassifiers);

/** dialog execution method for the display variable classifiers command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayVarClassifiers);

/** dialog execution method for the display scores command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayScores);

/** dialog execution method for the display solvers command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplaySolvers);

/** dialog execution method for the display decomposition command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayDecomposition);

/** dialog execution method for the display nblockscandidates command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayNBlockcandidates);

/** dialog execution method for the presolve command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecPresolve);

/** dialog execution method for the master command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetMaster);

/** dialog execution method for the set loadmaster command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetLoadmaster);

/** dialog execution method for the detect command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecDetect);

/** dialog execution method for the select command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSelect);

/** dialog execution method for the transform command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecTransform);

/** dialog execution method for the optimize command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecOptimize);

/** dialog execution method for the set detectors fast command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsFast);

/** dialog execution method for the set detectors off command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsOff);

/** dialog execution method for the set detectors default command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsDefault);

/** dialog execution method for the set detectors aggressive command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsAggressive);

/** dialog execution method for the set heuristics fast command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetHeuristicsFast);

/** dialog execution method for the set heuristics off command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetHeuristicsOff);

/** dialog execution method for the set heuristics aggressive command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetHeuristicsAggressive);

/** dialog execution method for the set separators default command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetSeparatorsDefault);

/** dialog execution method for the set separators fast command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetSeparatorsFast);

/** dialog execution method for the set separators off command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetSeparatorsOff);

/** dialog execution method for the set separators aggressive command */
GCG_EXPORT
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetSeparatorsAggressive);

/** includes or updates the GCG dialog menus in SCIP */
GCG_EXPORT
SCIP_RETCODE GCGincludeDialogGcg(
   GCG*                  gcg                 /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
