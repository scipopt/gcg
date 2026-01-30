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

/**@file   objdialog.h
 * @brief  C++ wrapper for dialogs
 * @author Kati Wolter
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_OBJDIALOG_H__
#define GCG_OBJDIALOG_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"
#include "gcg/gcg.h"


namespace gcg
{

/**
 *  @brief C++ wrapper for dialogs
 *
 *  This class defines the interface for dialogs implemented in C++. Note that there is a pure virtual function (this
 *  function has to be implemented). This function is: scip_exec().
 */
class ObjDialog : public scip::ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** GCG data structure */
   GCG* gcg;

   /** name of the dialog */
   char* scip_name_;

   /** description of the dialog */
   char* scip_desc_;

   /** default for whether the dialog is a menu */
   const SCIP_Bool scip_issubmenu_;

   /** default constructor */
   ObjDialog(
      GCG*               gcgstruct,          /**< GCG data structure */
      const char*        name,               /**< name of the dialog */
      const char*        desc,               /**< description of the dialog */
      SCIP_Bool          issubmenu           /**< default for whether the dialog is a menu */
      )
      : gcg(gcgstruct),
        scip_name_(0),
        scip_desc_(0),
        scip_issubmenu_(issubmenu)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjDialog()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
      gcg = NULL;
   }

   /** destructor of dialog to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_DIALOGFREE(x) in type_dialog.h
    */
   virtual SCIP_DECL_DIALOGFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** description output method of dialog
    *
    *  @see SCIP_DECL_DIALOGDESC(x) in type_dialog.h
    */
   virtual SCIP_DECL_DIALOGDESC(scip_desc)
   {  /*lint --e{715}*/
      SCIPdialogMessage(scip, 0, "%s", scip_desc_);
      return SCIP_OKAY;
   }

   /** execution method of dialog
    *
    *  @see SCIP_DECL_DIALOGEXEC(x) in type_dialog.h
    */
   virtual SCIP_DECL_DIALOGEXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the dialog for the given dialog object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyDialog* mydialog = new MyDialog(...);
 *       SCIP_CALL( GCGincludeObjDialog(gcg, &mydialog, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mydialog;    // delete dialog AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( GCGincludeObjDialog(gcg, new MyDialog(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyDialog is called here
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeObjDialog(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_DIALOG*          parentdialog,       /**< parent dialog */
   gcg::ObjDialog*       objdialog,          /**< dialog object */
   SCIP_Bool             deleteobject        /**< should the dialog object be deleted when dialog is freed? */
   );

#endif
