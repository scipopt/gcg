/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog_explore.h
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

/*!
 * \brief help enum to avoid code duplication for the toolbox methods of the detectors
 */
enum toolboxtype {
   PROPAGATE,
   FINISH,
   POSTPROCESS
};

/**
 * @brief method too handle user input for "explore" command
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPdialogExecSelect(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );


#ifdef __cplusplus
}
#endif

#endif
