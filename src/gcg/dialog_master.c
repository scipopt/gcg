/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   dialog_master.c
 * @brief  user interface dialog for master problem
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/dialog_default.h"
#include "gcg/dialog_master.h"



/** dialog execution method telling that a command is not available */
SCIP_DECL_DIALOGEXEC(GCGmasterDialogExecNotAvailable)
{  /*lint --e{715}*/
   SCIPdialogMessage(scip, NULL, "Not available in the master problem\n");

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** creates a root dialog */
static
SCIP_RETCODE createRootMasterDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   )
{
   SCIP_CALL( SCIPincludeDialog(scip, root,
         NULL, SCIPdialogExecMenuLazy, NULL, NULL,
         "GCG (master)", "GCG's master main menu", TRUE, NULL) );
   
   SCIP_CALL( SCIPsetRootDialog(scip, *root) );
   SCIP_CALL( SCIPreleaseDialog(scip, root) );
   *root = SCIPgetRootDialog(scip);
   
   return SCIP_OKAY;
}


/** includes or updates the master dialog menus in GCG */
SCIP_RETCODE GCGincludeDialogMaster(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* masterprob;
   SCIP_DIALOG* root;
   SCIP_DIALOG* dialog;

   masterprob = GCGgetMasterprob(gcg);
   assert(masterprob != NULL);

   /* root menu */
   root = SCIPgetRootDialog(masterprob);
   if( root == NULL )
   {
      SCIP_CALL( createRootMasterDialog(masterprob, &root) );
   }
   
   /* change */
   if( !SCIPdialogHasEntry(root, "change") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         GCGmasterDialogExecNotAvailable, NULL, NULL,
         "change", "(not available in master problem)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

   /* free */
   if( !SCIPdialogHasEntry(root, "free") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         GCGmasterDialogExecNotAvailable, NULL, NULL,
         "free", "(not available in master problem)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

   /* newstart */
   if( !SCIPdialogHasEntry(root, "newstart") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         GCGmasterDialogExecNotAvailable, NULL, NULL,
         "newstart", "(not available in master problem)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

   /* optimize */
   if( !SCIPdialogHasEntry(root, "optimize") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         GCGmasterDialogExecNotAvailable, NULL, NULL,
         "optimize", "(not available in master problem)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

   /* presolve */
   if( !SCIPdialogHasEntry(root, "presolve") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         GCGmasterDialogExecNotAvailable, NULL, NULL,
         "presolve", "(not available in master problem)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

   /* read */
   if( !SCIPdialogHasEntry(root, "read") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         GCGmasterDialogExecNotAvailable, NULL, NULL,
         "read", "(not available in master problem)", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

   /* quit */
   if( !SCIPdialogHasEntry(root, "quit") )
   {
      SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
         NULL,
         SCIPdialogExecQuit, NULL, NULL,
         "quit", "switch back to the original problem's dialog", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
   }

    /* validatesolve */
    if( !SCIPdialogHasEntry(root, "validatesolve") )
    {
        SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
           NULL,
           GCGmasterDialogExecNotAvailable, NULL, NULL,
           "validatesolve", "(not available in master problem)", FALSE, NULL) );
        SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
        SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
    }

    /* concurrentopt */
    if( !SCIPdialogHasEntry(root, "concurrentopt") )
    {
        SCIP_CALL( SCIPincludeDialog(masterprob, &dialog,
           NULL,
           GCGmasterDialogExecNotAvailable, NULL, NULL,
           "concurrentopt", "(not available in master problem)", FALSE, NULL) );
        SCIP_CALL( SCIPaddDialogEntry(masterprob, root, dialog) );
        SCIP_CALL( SCIPreleaseDialog(masterprob, &dialog) );
    }

   SCIP_CALL( SCIPincludeDialogDefaultBasic(masterprob) );

   return SCIP_OKAY;
}
