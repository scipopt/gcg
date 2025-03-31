/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   disp_master.c
 * 
 * @brief  master display columns
 * @author Gerald Gamrath
 * @author Christian Puchert
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/disp_master.h"
#include "scip/disp_default.h"
#include "gcg/gcg.h"

#define DISP_NAME_ORIGINAL         "original"
#define DISP_DESC_ORIGINAL         "display column printing a display line of the original SCIP instance"
#define DISP_HEAD_ORIGINAL         ""
#define DISP_WIDT_ORIGINAL         5
#define DISP_PRIO_ORIGINAL         80000
#define DISP_POSI_ORIGINAL         3550
#define DISP_STRI_ORIGINAL         TRUE

struct SCIP_DispData
{
   GCG*                    gcg;                /**< GCG data structure */
};

/*
 * Callback methods
 */

/** copy method for display plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DISPCOPY(dispCopyMaster)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(disp != NULL);

   /* call inclusion method of default SCIP display plugin */
   SCIP_CALL( SCIPincludeDispDefault(scip) );

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' printing a display column of the original SCIP instance */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputOriginal)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_ORIGINAL) == 0);
   assert(scip != NULL);
   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);

   SCIP_CALL( SCIPprintDisplayLine(GCGgetOrigprob(dispdata->gcg), file, SCIP_VERBLEVEL_HIGH, FALSE) );

   return SCIP_OKAY;
}

/** destructor method of display plugin */
static
SCIP_DECL_DISPFREE(SCIPdispFreeMaster)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;
   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);

   SCIPfreeBlockMemory(scip, &dispdata);

   return SCIP_OKAY;
}

/*
 * default display columns specific interface methods
 */

/** includes the default display columns in SCIP */
SCIP_RETCODE GCGincludeDispMaster(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP*                 masterprob          /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata;

   SCIP_CALL( SCIPallocBlockMemory(masterprob, &dispdata) );
   dispdata->gcg = gcg;

   SCIP_CALL( SCIPincludeDisp(masterprob, DISP_NAME_ORIGINAL, DISP_DESC_ORIGINAL, DISP_HEAD_ORIGINAL,
         SCIP_DISPSTATUS_AUTO, dispCopyMaster, SCIPdispFreeMaster, NULL, NULL, NULL, NULL, SCIPdispOutputOriginal, dispdata,
         DISP_WIDT_ORIGINAL, DISP_PRIO_ORIGINAL, DISP_POSI_ORIGINAL, DISP_STRI_ORIGINAL) );

   return SCIP_OKAY;
}
