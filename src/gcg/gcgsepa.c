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

/**@file   gcgsepa.c
 * @brief  public methods for GCG separators
 * @author Christian Puchert
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/gcg.h"

#include "gcg/pub_gcgsepa.h"


/** resets the parameters to disable separators */
static
SCIP_RETCODE setSeparatorsOff(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /* set specific parameters for GCG basis separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "basis") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/basis/enable", FALSE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/enable = FALSE\n");
   }

   /* set specific parameters for GCG master separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "master") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/master/enable", FALSE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/enable = FALSE\n");
   }

   return SCIP_OKAY;
}

/** resets the parameters to their default value */
static
SCIP_RETCODE setSeparatorsDefault(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /* set specific parameters for GCG basis separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "basis") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/basis/enable", TRUE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/enable = TRUE\n");

      SCIP_CALL( SCIPsetIntParam(origprob, "sepa/basis/paramsetting", 0) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/paramsetting = %d\n", 0);
   }

   /* set specific parameters for GCG master separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "master") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/master/enable", TRUE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/enable = TRUE\n");

      SCIP_CALL( SCIPsetIntParam(origprob, "sepa/master/paramsetting", 0) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/paramsetting = %d\n", 0);
   }

   return SCIP_OKAY;
}

/** sets the parameters to aggressive values */
static
SCIP_RETCODE setSeparatorsAggressive(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /* set specific parameters for GCG basis separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "basis") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/basis/enable", TRUE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/enable = TRUE\n");

      SCIP_CALL( SCIPsetIntParam(origprob, "sepa/basis/paramsetting", 1) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/paramsetting = %d\n", 1);
   }

   /* set specific parameters for GCG master separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "master") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/master/enable", TRUE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/enable = TRUE\n");

      SCIP_CALL( SCIPsetIntParam(origprob, "sepa/master/paramsetting", 1) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/paramsetting = %d\n", 1);
   }

   return SCIP_OKAY;
}

/** sets the parameters to fast values */
static
SCIP_RETCODE setSeparatorsFast(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* origprob = GCGgetOrigprob(gcg);
   /* set specific parameters for GCG basis separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "basis") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/basis/enable", TRUE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/enable = TRUE\n");

      SCIP_CALL( SCIPsetIntParam(origprob, "sepa/basis/paramsetting", 2) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/basis/paramsetting = %d\n", 2);
   }

   /* set specific parameters for GCG master separator, if the separator is included */
#ifndef NDEBUG
   if( SCIPfindSepa(GCGgetMasterprob(gcg), "master") != NULL )
#endif
   {
      SCIP_CALL( SCIPsetBoolParam(origprob, "sepa/master/enable", TRUE) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/enable = TRUE\n");

      SCIP_CALL( SCIPsetIntParam(origprob, "sepa/master/paramsetting", 2) );
      SCIPverbMessage(origprob, SCIP_VERBLEVEL_NORMAL, NULL, "sepa/master/paramsetting = %d\n", 2);
   }

   return SCIP_OKAY;
}

/** sets separator parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separator parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separator is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separator are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turns off all separators
 */
SCIP_RETCODE GCGsetSeparators(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_PARAMSETTING     paramsetting        /**< parameter settings */
   )
{
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);

   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( setSeparatorsAggressive(gcg) );
      break;
   case SCIP_PARAMSETTING_OFF:
      SCIP_CALL(setSeparatorsOff(gcg));
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( setSeparatorsFast(gcg) );
      break;
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( setSeparatorsDefault(gcg) );
      break;
   default:
      SCIPerrorMessage("The given paramsetting is invalid!\n");
      break;
   }

   return SCIP_OKAY;
}
