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

/**@file   dec_connected.c
 * @ingroup DETECTORS
 * @brief  detector for classical and blockdiagonal problems
 * @author Martin Bergner
 * @todo allow decompositions with only one pricing problem by just removing generalized covering and
 *       partitioning constraints
 * The detector will detect block diagonal matrix structures as wells as generalized
 * set partitioning or covering master problems.
 *
 * It works as follows:
 * - It implicitly builds a graph with one vertex for every constraint and edges between constraints that
 *   share a node
 * - All vertices belonging to constraints of the form \f$\sum x_i = a \f$ for
 *   \f$x_i \in \mathbb Z, a\in \mathbb Z\f$ or of the form \f$\sum x_i \geq 1 \f$
 *   for \f$x_i \in \{0,1\} \f$ are removed
 * - The pricing problems correspond to connected components in the remaining graph
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "dec_connected.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"

/* constraint handler properties */
#define DEC_DETECTORNAME          "connected"    /**< name of detector */
#define DEC_DESC                  "Detector for classical and block diagonal problems" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0              /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               'C'            /**< display character of detector */

#define DEC_ENABLED               FALSE           /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       FALSE  /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE        /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the finishing be enabled */
#define DEFAULT_SETPPCINMASTER    TRUE           /**< should the extended structure be detected */
#define DEC_SKIP                  FALSE          /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            TRUE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
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
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int i;
   int nvars;
   SCIP_Bool relevant = TRUE;
   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("cons %s is ", SCIPconsGetName(cons));
   if( GCGconsGetType(cons) == setcovering || GCGconsGetType(cons) == setpartitioning || GCGconsGetType(cons) == logicor )
   {
      SCIPdebugPrintf("setcov, part or logicor.\n");
      return TRUE;
   }
   nvars = GCGconsGetNVars(scip, cons);
   vars = NULL;
   vals = NULL;
   if( nvars > 0 )
   {
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &vars, nvars) );
      SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &vals, nvars) );
      SCIP_CALL_ABORT( GCGconsGetVars(scip, cons, vars, nvars) );
      SCIP_CALL_ABORT( GCGconsGetVals(scip, cons, vals, nvars) );
   }

   /* check vars and vals for integrality */
   for( i = 0; i < nvars && relevant; ++i )
   {
      assert(vars != NULL);
      assert(vals != NULL);

      if( !SCIPvarIsIntegral(vars[i]) && !SCIPvarIsBinary(vars[i]) )
      {
         SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[i]) );
         relevant = FALSE;
      }
      if( !SCIPisEQ(scip, vals[i], 1.0) )
      {
         SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[i]), vals[i] );
         relevant = FALSE;
      }
   }

   /* free temporary data  */
   SCIPfreeMemoryArrayNull(scip, &vals);
   SCIPfreeMemoryArrayNull(scip, &vars);

   SCIPdebugPrintf("%s master\n", relevant ? "in" : "not in");
   return relevant;
}

/** creates array for constraints in the master */
static
SCIP_RETCODE createMasterconssArray(
   SCIP*                 scip,               /**< SCIP data structure */
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
      if( isConsMaster(scip, conss[i]) )
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

/** looks for connected components in the constraints in detectordata */
static
SCIP_RETCODE findConnectedComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP**          decomp,             /**< decomposition structure */
   SCIP_Bool             findextended,       /**< whether the classical structure should be detected */
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

   masterisempty = findextended;
   if( findextended )
   {
      SCIP_CALL( createMasterconssArray(scip, &conss, &nconss, &masterisempty, &pricingisempty) );
      if( pricingisempty )
      {
         *result = SCIP_DIDNOTFIND;
         SCIPfreeMemoryArrayNull(scip, &conss);
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( DECcreateDecompFromMasterconss(scip, decomp, conss, nconss) );

   if( DECdecompGetNBlocks(*decomp) > 1 )
      *result = SCIP_SUCCESS;
   else if( DECdecompGetNBlocks(*decomp) == 1 && !masterisempty && findextended )
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



/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeConnected)
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
DEC_DECL_INITDETECTOR(detectorInitConnected)
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

/** detector structure detection method, tries to detect a structure in the problem */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectConnected)
{
   int runs;
   int i;
   SCIP_Bool detectextended;
   *result = SCIP_DIDNOTFIND;
   runs = detectordata->setppcinmaster ? 2:1;
   detectextended = FALSE;

   *ndecdecomps = 0;
   *decdecomps = NULL;
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/

   for( i = 0; i < runs && *result != SCIP_SUCCESS; ++i )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting %s structure:", detectextended ? "set partitioning master":"purely block diagonal" );

      SCIP_CALL( findConnectedComponents(scip, &((*decdecomps)[0]), detectextended, result) );

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

      if( detectordata->setppcinmaster == TRUE && *result != SCIP_SUCCESS )
      {
         detectextended = TRUE;
      }
   }

   if( *ndecdecomps == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, decdecomps);
   }
   return SCIP_OKAY;
}

#define detectorExitConnected NULL
#define detectorPropagateSeeedConnected NULL
#define detectorFinishSeeedConnected NULL

#define detectorPostprocessSeeedConnected NULL
#define setParamAggressiveConnected NULL
#define setParamDefaultConnected NULL
#define setParamFastConnected NULL


/*
 * detector specific interface methods
 */

/** creates the handler for connected constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorConnected(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create connected constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);


   detectordata->blockdiagonal = FALSE;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING,DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE,
      detectordata, detectorDetectConnected, detectorFreeConnected, detectorInitConnected, detectorExitConnected, detectorPropagateSeeedConnected, NULL, NULL, detectorFinishSeeedConnected, detectorPostprocessSeeedConnected, setParamAggressiveConnected, setParamDefaultConnected, setParamFastConnected) );


   /* add connected constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/detectors/connected/setppcinmaster", "Controls whether SETPPC constraints chould be ignored while detecting and be directly placed in the master", &detectordata->setppcinmaster, FALSE, DEFAULT_SETPPCINMASTER, NULL, NULL) );

   return SCIP_OKAY;
}
