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

/**@file   gcg.cpp
 * @brief  methods for working with gcg structure
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "gcg/struct_gcg.h"
#include "gcg/gcgplugins.h"
#include "gcg/benders_gcg.h"
#include "gcg/pricer_gcg.h"
#include "gcg/objpricer_gcg.h"

/** create a gcg struct */
extern "C"
SCIP_RETCODE GCGcreate(
   GCG**                gcg               /**< pointer to store GCG data structure */
   )
{
   SCIP* origprob;
   *gcg = (GCG*) malloc(sizeof(GCG));
   SCIP_CALL( SCIPcreate(&origprob) );

   (*gcg)->origprob = origprob;
   (*gcg)->masterprob = NULL;
   (*gcg)->bendersmasterprob = NULL;
   (*gcg)->dwmasterprob = NULL;
   (*gcg)->relax = NULL;
   (*gcg)->pricer = NULL;
   (*gcg)->sepaorig = NULL;

   /* include gcg plugins */
   SCIP_CALL( GCGincludeGcgPlugins(*gcg) );
   return SCIP_OKAY;
}

/** free a gcg column */
extern "C"
SCIP_RETCODE GCGfree(
   GCG**                gcg               /**< pointer to gcg structure */
   )
{
   SCIP* origprob;

   if( gcg == NULL || (*gcg) == NULL )
      return SCIP_OKAY;

   origprob = (*gcg)->origprob;
   assert(origprob != NULL);
   SCIP_CALL( SCIPfree(&origprob) );
   free(*gcg);
   return SCIP_OKAY;
}

#ifndef NDEBUG
/** returns the original problem */
extern "C"
SCIP* GCGgetOrigprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->origprob;
}
#endif

#ifndef NDEBUG
/** returns the active master problem (may change until solving is initiated) */
extern "C"
SCIP* GCGgetMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->masterprob;
}
#endif

#ifndef NDEBUG
/** returns the benders master problem (also used to solve the original problem directly) */
extern "C"
SCIP* GCGgetBendersMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->bendersmasterprob;
}
#endif

#ifndef NDEBUG
/** returns the dw master problem */
extern "C"
SCIP* GCGgetDwMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->dwmasterprob;
}
#endif

/** returns the GCG data structure */
extern "C"
GCG* GCGmasterGetGcg(
   SCIP*                masterprob           /**< SCIP data structure */
   )
{
   SCIP_BENDERS* benders;
   SCIP_PRICER* pricer;

   assert(masterprob != NULL);

   /* retrieving the Benders' decomposition and the pricer plugins. There should only be one or the other for a given
    * master problem. If there are both, then an error is returned */
   benders = SCIPfindBenders(masterprob, "gcg");
   pricer = SCIPfindPricer(masterprob, "gcg");
   assert((benders != NULL && pricer == NULL) || (pricer != NULL && benders == NULL));

   if( benders != NULL && pricer == NULL )
   {
      return GCGbendersGetGcg(benders);
   }
   else if( pricer != NULL && benders == NULL )
   {
      return GCGpricerGetGcg(masterprob);
   }

   SCIPerrorMessage("There must exist either a pricer or a benders or neither, not both.\n");

   return NULL;
}

/** returns the GCG data structure */
extern "C"
GCG* GCGorigGetGcg(
   SCIP*                origprob             /**< SCIP data structure */
   )
{
   assert(origprob != NULL);
   return GCGrelaxGetGcg(origprob);
}

#ifndef NDEBUG
/** gets GCG's relaxator */
extern "C"
SCIP_RELAX* GCGgetRelax(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->relax;
}
#endif

#ifndef NDEBUG
/** gets the GCG pricer */
ObjPricerGcg* GCGgetObjPricer(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return reinterpret_cast<ObjPricerGcg*>(gcg->pricer);
}
#endif

#ifndef NDEBUG
/** gets the orig separator */
SCIP_SEPA* GCGgetSepaorig(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->sepaorig;
}
#endif
