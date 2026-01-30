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

/**@file   objdialog.cpp
 * @brief  C++ wrapper for dialogs
 * @author Kati Wolter
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -e64 disable useless NULL pointer errors */
#include <cassert>
#include "gcg/objdialog.h"


/*
 * Data structures
 */

/** dialog data */
struct SCIP_DialogData
{
   gcg::ObjDialog*       objdialog;          /**< dialog object */
   SCIP_Bool             deleteobject;       /**< should the dialog object be deleted when dialog is freed? */
};


/*
 * Callback methods of dialog
 */

extern "C"
{

#define dialogCopyObj NULL
/** destructor of dialog to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DIALOGFREE(dialogFreeObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_free(scip, dialog) );

   /* free dialog object */
   if( dialogdata->deleteobject )
      delete dialogdata->objdialog;

   /* free dialog data */
   delete dialogdata;
   SCIPdialogSetData(dialog, 0); /*lint !e64*/

   return SCIP_OKAY;
}


/** description output method of dialog */
static
SCIP_DECL_DIALOGDESC(dialogDescObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_desc(scip, dialog) );

   return SCIP_OKAY;
}

/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_exec(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}
}



/*
 * dialog specific interface methods
 */

/** creates the dialog for the given dialog object and includes it in SCIP */
SCIP_RETCODE GCGincludeObjDialog(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_DIALOG*          parentdialog,       /**< parent dialog */
   gcg::ObjDialog*       objdialog,          /**< dialog object */
   SCIP_Bool             deleteobject        /**< should the dialog object be deleted when dialog is freed? */
   )
{/*lint --e{429} */
   SCIP* scip = GCGgetOrigprob(gcg);
   assert(scip != 0);
   assert(objdialog != 0);
   assert(parentdialog != 0);

   /* create, include, and release dialog */
   if( !SCIPdialogHasEntry(parentdialog, objdialog->scip_name_) )
   {
      SCIP_DIALOGDATA* dialogdata; /*lint !e593*/
      SCIP_DIALOG* dialog;
      SCIP_RETCODE retcode;

      dialog = 0;

      /* create dialog data */
      dialogdata = new SCIP_DIALOGDATA;
      dialogdata->objdialog = objdialog;
      dialogdata->deleteobject = deleteobject;

      retcode = SCIPincludeDialog(scip, &dialog, dialogCopyObj, dialogExecObj, dialogDescObj, dialogFreeObj,
         objdialog->scip_name_, objdialog->scip_desc_, objdialog->scip_issubmenu_, dialogdata);
      if( retcode != SCIP_OKAY )
      {
          delete dialogdata;
          SCIP_CALL( retcode );
      }
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );  /*lint !e593*/
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );                /*lint !e593*/
   }  /*lint !e593*/

   return SCIP_OKAY;  /*lint !e593*/
}
