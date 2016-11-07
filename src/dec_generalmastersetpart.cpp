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

/**@file   dec_generalmastersetpart.cpp
 * @ingroup DETECTORS
 * @brief  detector generalmastersetpart (put your description here)
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_generalmastersetpart.h"
#include "cons_decomp.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip_misc.h"

#include <iostream>

/* constraint handler properties */
#define DEC_DETECTORNAME         "generalmastersetpart"       /**< name of detector */
#define DEC_DESC                 "detector generalmastersetpart" /**< description of detector*/
#define DEC_PRIORITY             0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              '?'         /**< display character of detector */
#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#define DEC_SKIP                 FALSE       /**< should detector be skipped if other detectors found decompositions */

/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
{
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(freeGeneralmastersetpart)
{
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** destructor of detector to free detector data (called before the solving process begins) */
#if 0
static
DEC_DECL_EXITDETECTOR(exitGeneralmastersetpart)
{ /*lint --e{715}*/

   SCIPerrorMessage("Exit function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define exitGeneralmastersetpart NULL
#endif

/** detection initialization function of detector (called before solving is about to begin) */
#if 0
static
DEC_DECL_INITDETECTOR(initGeneralmastersetpart)
{ /*lint --e{715}*/

   SCIPerrorMessage("Init function of detector <%s> not implemented!\n", DEC_DETECTORNAME);
   SCIPABORT();

   return SCIP_OKAY;
}
#else
#define initGeneralmastersetpart NULL
#endif

/** detection function of detector */
static DEC_DECL_DETECTSTRUCTURE(detectGeneralmastersetpart)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPerrorMessage("Detection function of detector <%s> not implemented!\n", DEC_DETECTORNAME)
;   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

static DEC_DECL_PROPAGATESEEED(propagateSeeedGeneralmastersetpart)
{
   *result = SCIP_DIDNOTFIND;

   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   bool relevant = true;
   int value;


   gcg::Seeed* seeed;
   seeed = new gcg::Seeed(seeedPropagationData->seeedToPropagate, seeedPropagationData->seeedpool);
   seeed->setDetectorPropagated(seeedPropagationData->seeedpool->getIndexForDetector(detector));

   if(!seeed->areOpenVarsAndConssCalculated())
   {
      seeed->calcOpenconss();
      seeed->calcOpenvars();
      seeed->setOpenVarsAndConssCalculated(true);
   }

   /** set open setpartitioning constraints to Master */
   for( int i = 0; i < seeed->getNOpenconss(); ++i)
   {
      cons = seeedPropagationData->seeedpool->getConsForIndex(seeed->getOpenconss()[i]);
      if( GCGconsGetType(cons) == setpartitioning )
      {
         seeed->setConsToMaster(seeed->getOpenconss()[i]);
         seeed->deleteOpencons(seeed->getOpenconss()[i]);
      }
      else if(GCGconsGetType(cons) != logicor && GCGconsGetType(cons) != setcovering && GCGconsGetType(cons) != setpacking )
      {
         nvars = GCGconsGetNVars(scip, cons);
         vars = NULL;
         vals = NULL;
         value = GCGconsGetLhs(scip, cons);

         if( SCIPisNegative(scip, value) || GCGconsGetRhs(scip,cons)!= value)
            relevant = false;
         if( nvars > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
            SCIP_CALL( GCGconsGetVars(scip, cons, vars, nvars) );
            SCIP_CALL( GCGconsGetVals(scip, cons, vals, nvars) );
         }
         for( int j = 0; j < nvars && relevant; ++j )
         {
            assert(vars != NULL);
            assert(vals != NULL);

            if( !SCIPvarIsIntegral(vars[j]) && !SCIPvarIsBinary(vars[j]) )
            {
               SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[j]) );
               relevant = FALSE;
            }
            if( !SCIPisEQ(scip, vals[j], 1.0) )
            {
               SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[j]), vals[j] );
               relevant = FALSE;
            }
         }
         SCIPfreeBufferArrayNull(scip, &vals);
         SCIPfreeBufferArrayNull(scip, &vars);

         if(relevant)
         {
            seeed->setConsToMaster(seeed->getOpenconss()[i]);
            seeed->deleteOpencons(seeed->getOpenconss()[i]);
         }
      }
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropagationData->newSeeeds), 1) );
   seeedPropagationData->newSeeeds[0] = seeed;
   seeedPropagationData->nNewSeeeds = 1;
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
/*
 * detector specific interface methods
 */

/** creates the handler for generalmastersetpart detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorGeneralmastersetpart(SCIP* scip /**< SCIP data structure */
)
{
   DEC_DETECTORDATA* detectordata;

   /**@todo create generalmastersetpart detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, detectordata, detectGeneralmastersetpart, freeGeneralmastersetpart, initGeneralmastersetpart, exitGeneralmastersetpart, propagateSeeedGeneralmastersetpart));

   /**@todo add generalmastersetpart detector parameters */

   return SCIP_OKAY;
}





