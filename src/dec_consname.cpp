/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   dec_consname.c
 * @ingroup DETECTORS
 * @brief  detector for classical and blockdiagonal problems
 * @author Martin Bergner
 * @todo allow decompositions with only one pricing problem by just removing generalized covering and
 *       partitioning constraints
 * The detector will detect block diagonal matrix structures as wells as generalized
 * set partitioning or covering master problems.
 *
 * It works as follows:
 * - give regular expression
 * - all constraints whose names match the regular expression go into the master
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

#include <cstring>
#include <cassert>
#include <regex>

#include "dec_consname.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "consname"         /**< name of detector */
#define DEC_DESC                 "Build master constraints by name" /**< description of detector*/
#define DEC_PRIORITY             0              /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'N'            /**< display character of detector */

#define DEC_ENABLED              TRUE           /**< should the detection be enabled */
#define DEFAULT_REGEX            "(capacity)(.*)"         /**< default regular expression that is used to decide mastercons */
#define DEC_SKIP                 FALSE          /**< should detector be skipped if others found detections */

#define DEFAULT_SETPPCINMASTER   TRUE           /**< should the extended structure be detected */

/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
   char* regex;                              /**< regular expression that is used to decide mastercons */
   SCIP_Bool blockdiagonal;                  /**< flag to indicate whether the problem is block diagonal */
   SCIP_Bool setppcinmaster;                 /**< flag to indicate whether setppc constraints should always be in the master */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/* returns true if the constraint should be a master constraint and false otherwise */
static
SCIP_Bool isConsMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   char* consname;
   std::regex expr(detectordata->regex);
   consname = SCIPconsGetName(cons);

   if( std::regex_match(consname, expr) )
   {
      SCIPdebugPrintf("Name %s matches regular expression %s\n\n", consname, detectordata->regex);
      return TRUE;
   }

   return FALSE;
}

/** creates array for constraints in the master */
static
SCIP_RETCODE createMasterconssArray(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data */
   SCIP_CONS***          masterconss,        /**< pointer to return masterconss array */
   int*                  nmasterconss,       /**< pointer to return number of masterconss */
   SCIP_Bool*            masterisempty,      /**< indicates if the master is empty (all constraints in the subproblem) */
   SCIP_Bool*            pricingisempty      /**< indicates if the pricing is empty (all constraints in the master) */
   )
{
   int i;

   SCIP_CONS** conss;
   int nconss;

   assert(scip != NULL);
   assert(masterconss != NULL);
   assert(nmasterconss != NULL);
   assert(masterisempty != NULL);
   assert(pricingisempty != NULL);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, masterconss, nconss) );

   for (i = 0; i < nconss; ++i)
   {
      if( isConsMaster(scip, detectordata, conss[i]) )
      {
         SCIPdebugMessage("Constraint <%s> to be placed in master.\n", SCIPconsGetName(conss[i]));
         (*masterconss)[*nmasterconss] = conss[i];
         ++(*nmasterconss);
      }
   }

   if( *nmasterconss > 0 )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, masterconss, *nmasterconss) );
   }
   else
   {
      SCIPfreeMemoryArrayNull(scip, masterconss);
      *nmasterconss = 0;
      *masterconss = NULL;
   }

   *pricingisempty = (*nmasterconss == nconss);
   *masterisempty = (*nmasterconss == 0);

   return SCIP_OKAY;
}

/** looks for consname components in the constraints in detectordata */
static
SCIP_RETCODE findConnectedComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< detector data */
   DEC_DECOMP**          decomp,             /**< decomposition structure */
   SCIP_RESULT*          result              /**< result pointer to indicate success oder failure */
   )
{
   int nconss = 0;
   SCIP_CONS** conss = NULL;

   SCIP_Bool masterisempty;
   SCIP_Bool pricingisempty;


   assert(scip != NULL);
   assert(decomp != NULL);
   assert(result != NULL);

   SCIP_CALL( createMasterconssArray(scip, detectordata, &conss, &nconss, &masterisempty, &pricingisempty) );
   if( pricingisempty )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPfreeMemoryArrayNull(scip, &conss);
      return SCIP_OKAY;
   }

   SCIP_CALL( DECcreateDecompFromMasterconss(scip, decomp, conss, nconss) );

   if( DECdecompGetNBlocks(*decomp) > 1 )
      *result = SCIP_SUCCESS;
   else
   {
      SCIP_CALL( DECdecompFree(scip, decomp) );
      *decomp = NULL;
      *result = SCIP_DIDNOTFIND;
   }

   SCIPfreeMemoryArrayNull(scip, &conss);
   return SCIP_OKAY;
}



/** destructor of detector to free detector data (called when SCIP is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeConsname)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detection initialization function of detector (called before solving is about to begin) */
static
DEC_DECL_INITDETECTOR(detectorInitConsname)
{  /*lint --e{715}*/

   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->blockdiagonal = FALSE;

   return SCIP_OKAY;
}

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectConsname)
{
   *result = SCIP_DIDNOTFIND;
   *ndecdecomps = 0;
   *decdecomps = NULL;
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting structure by constraint name matching regular expression %s:", detectordata->regex);

   SCIP_CALL( findConnectedComponents(scip, detectordata, &((*decdecomps)[0]), result) );

   if( *result == SCIP_SUCCESS )
   {
      DEC_DECOMP *newdecomp;
      assert((*decdecomps)[0] != NULL);
      SCIP_CALL( DECcreatePolishedDecomp(scip, (*decdecomps)[0], &newdecomp) );
      if( newdecomp != NULL )
      {
         SCIP_CALL( DECdecompFree(scip, &((*decdecomps)[0])) );
         (*decdecomps)[0] = newdecomp;
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found with %d blocks.\n", DECdecompGetNBlocks((*decdecomps)[0]));
      detectordata->blockdiagonal = DECdecompGetType((*decdecomps)[0]) == DEC_DECTYPE_DIAGONAL;
      *ndecdecomps = 1;
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
   }

   if( *ndecdecomps == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, decdecomps);
   }
   return SCIP_OKAY;
}


/*
 * detector specific interface methods
 */

/** creates the handler for consname constraints and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorConsname(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create consname constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   detectordata->blockdiagonal = FALSE;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP,
      detectordata, detectorDetectConsname, detectorFreeConsname, detectorInitConsname, NULL) );

   /* add consname constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/consname/setppcinmaster", "Controls whether SETPPC constraints chould be ignored while detecting and be directly placed in the master", &detectordata->setppcinmaster, FALSE, DEFAULT_SETPPCINMASTER, NULL, NULL) );
   SCIP_CALL( SCIPaddStringParam(scip, "detectors/consname/regex", "All cons whose name match this regular expression will be mastercons", &detectordata->regex, FALSE, DEFAULT_REGEX, NULL, NULL) );

   return SCIP_OKAY;
}
