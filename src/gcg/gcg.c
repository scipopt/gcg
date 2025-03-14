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

/**@file   gcg.c
 * @brief  methods for working with gcg structure
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "gcg/struct_gcg.h"
#include "gcg/gcgplugins.h"
#include "gcg/benders_gcg.h"
#include "gcg/pricer_gcg.h"

/** create a gcg struct */
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

   /* include gcg plugins */
   SCIP_CALL( GCGincludeGcgPlugins(*gcg) );
   return SCIP_OKAY;
}

/** free a gcg column */
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
SCIP* GCGgetOrigprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->origprob;
}
#endif

#ifndef NDEBUG
/** returns the master problem */
SCIP* GCGgetMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->masterprob;
}
#endif

#ifndef NDEBUG
/** returns the benders master problem */
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
SCIP* GCGgetDwMasterprob(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->dwmasterprob;
}
#endif

/** returns the GCG data structure */
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
GCG* GCGorigGetGcg(
   SCIP*                origprob             /**< SCIP data structure */
   )
{
   assert(origprob != NULL);
   return GCGrelaxGetGcg(origprob);
}

#ifndef NDEBUG
SCIP_RELAX* GCGgetRelax(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   assert(gcg != NULL);
   return gcg->relax;
}
#endif
