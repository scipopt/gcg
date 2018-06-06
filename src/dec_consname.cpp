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

/**@file   dec_consname.c
 * @ingroup DETECTORS
 * @brief  structure detection by constraint names (via regular expressions)
 * @author Jonas Witt
 *
 * The detector will detect a structure depending on the name of constraints
 *
 * It works as follows:
 * - given a regular expression,
 * - all constraints whose names match the regular expression will be master constraints,
 * - the pricing problems correspond to connected components in the remaining graph.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cstring>
#include <cassert>
#include <regex>

#include "dec_consname.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "consname"     /**< name of detector */
#define DEC_DESC                 "Build master constraints by name" /**< description of detector*/
#define DEC_PRIORITY             0              /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'N'            /**< display character of detector */

#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */

#define DEC_ENABLED              FALSE          /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL      FALSE        /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING     FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the finishing be enabled */
#define DEFAULT_REGEX            "(consname)(.*)" /**< default regular expression that is used to decide mastercons */
#define DEC_SKIP                 FALSE          /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL         FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE           FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */


/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
   char* regex;                              /**< regular expression that is used to decide mastercons */
};


/*
 * Local methods
 */

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
   const char* consname;
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

#define detectorInitConsname NULL

#define detectorExitConsname NULL

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
      assert((*decdecomps)[0] != NULL);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found with %d blocks.\n", DECdecompGetNBlocks((*decdecomps)[0]));
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


#define detectorPropagateSeeedConsname NULL
#define detectorFinishSeeedConsname NULL
#define detectorPostprocessSeeedConsname NULL

static
DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveConsname)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   return SCIP_OKAY;

}


static
DEC_DECL_SETPARAMDEFAULT(setParamDefaultConsname)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );

   return SCIP_OKAY;

}

static
DEC_DECL_SETPARAMFAST(setParamFastConsname)
{
   char setstr[SCIP_MAXSTRLEN];

   const char* name = DECdetectorGetName(detector);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(scip, setstr, FALSE ) );


   return SCIP_OKAY;

}


/*
 * detector specific interface methods
 */

/** creates the consmname detector and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorConsname(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create consname detector data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   detectordata->regex = NULL;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND,
      DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING,DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE, detectordata, detectorDetectConsname,
      detectorFreeConsname, detectorInitConsname, detectorExitConsname, detectorPropagateSeeedConsname, NULL, NULL, detectorFinishSeeedConsname, detectorPostprocessSeeedConsname, setParamAggressiveConsname, setParamDefaultConsname, setParamFastConsname) );

   /* add consname detector parameters */
   SCIP_CALL( SCIPaddStringParam(scip, "detection/detectors/consname/regex", "All cons whose name match this regular expression will be mastercons", &detectordata->regex, FALSE, DEFAULT_REGEX, NULL, NULL) );

   return SCIP_OKAY;
}
