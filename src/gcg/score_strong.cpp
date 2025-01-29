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

/**@file    score_strong.cpp
 * @brief   strong score
 * @author  Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <sstream>

#include <scip/scipdefplugins.h>

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.h"
#include "cons_decomp.hpp"
#include "scip_misc.h"
#include "score.h"
#include "score_strong.h"


/* score properties */
#define SCORE_NAME                "strong decomposition score"
#define SCORE_SHORTNAME           "strong"
#define SCORE_DESC                "strong decomposition score"

#define DEFAULT_STRONGTIMELIMIT                       30.         /**< timelimit for strong decompotition score calculation per partialdec */
#define DEFAULT_DUALVALRANDOMMETHOD                   1           /**< default value for method to dual initialization of dual values for strong decomposition: 1) naive, 2) expected equal, 3) expected overestimation */
#define DEFAULT_COEFFACTORORIGVSRANDOM                0.5         /**< default value for convex coefficient for orig dual val (1-this coef is factor for random dual value)  */

#define DEFAULT_SCORECOEF_FASTBENEFICIAL              1.          /**< coefficient for fast & beneficial in strong decomposition score computation */
#define DEFAULT_SCORECOEF_MEDIUMBENEFICIAL            0.75        /**< coefficient for not fast but beneficial in strong decomposition score computation */
#define DEFAULT_SCORECOEF_FASTNOTBENEFICIAL           0.3         /**< coefficient for fast & not beneficial in strong decomposition score computation */
#define DEFAULT_SCORECOEF_MEDIUMNOTBENEFICIAL         0.1         /**< coefficient for not & not beneficial in strong decomposition score computation */
#define DEFAULT_RANDPARTIALDEC                        23          /**< initial random partialdec */

/*
 * Data structures
 */
struct GCG_ScoreData
{
   std::vector<SCIP_Real>* dualvalsrandom;                        /**< vector of random dual values, used for strong detection scores */
   std::vector<SCIP_Real>* dualvalsoptimaloriglp;                 /**< vector of dual values of the optimal solved original lp, used for strong detection scores */
   int                   strongdetectiondualvalrandommethod;      /**< method to dual initialization of dual values for strong decomposition: 1) naive, 2) expected equal, 3) expected overestimation */
   SCIP_Real             coeffactororigvsrandom;                  /**< convex coefficient for orig dual val (1-this coef is factor for random dual value)  */
   SCIP_Bool             dualvalsoptimaloriglpcalculated;         /**< are the optimal dual values from original lp calulated? used for strong detection scores */
   SCIP_Bool             dualvalsrandomset;                       /**< are the random dual values set, used for strong detection scores */
   SCIP_Real             strongtimelimit;                         /**< timelimit for calculating strong decomposition score for one partialdec */
};

// TODO ref add comments! this is used in strong decomposition score calculation in @see shuffleDualvalsRandom()
enum GCG_Random_dual_methods
{
   GCG_RANDOM_DUAL_NAIVE                  =  0,         /**<  */
   GCG_RANDOM_DUAL_EXPECTED_EQUAL         =  1,         /**<  */
   GCG_RANDOM_DUAL_EXPECTED_OVERESTIMATE  =  2         /**<  */
};
typedef enum GCG_Random_dual_methods GCG_RANDOM_DUAL_METHOD;

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/**
 * @brief method that shuffles randomly and set dual variable values, used for strong detection score
 * @returns scip return code
 */
static
SCIP_RETCODE shuffleDualvalsRandom(
   SCIP* scip,             /**< SCIP data structure */
   GCG_SCORE* score,       /**< score */
   SCIP_Bool transformed   /**< whether the problem is tranformed yet */
   )
{
   GCG_SCOREDATA* scoredata = GCGscoreGetData(score);

   GCG_RANDOM_DUAL_METHOD usedmethod;
   SCIP_RANDNUMGEN* randnumgen;

   int method;
   int nconss = SCIPgetNConss(scip);

   SCIPgetIntParam(scip, "detection/scores/strong/dualvalrandommethod", &method);

   /* default method == 1 */
   usedmethod = GCG_RANDOM_DUAL_NAIVE;

   if( method == 2 )
      usedmethod = GCG_RANDOM_DUAL_EXPECTED_EQUAL;
   else if( method == 3 )
      usedmethod = GCG_RANDOM_DUAL_EXPECTED_OVERESTIMATE;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "set dual val random method to %d. \n", method );

   scoredata->dualvalsrandom->clear();
   scoredata->dualvalsrandom->resize(nconss, 0.);

   /* create random number generator */
   // TODO ref replace default for random partialdec with parameter
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, DEFAULT_RANDPARTIALDEC, TRUE) );

   gcg::DETPROBDATA* detprobdata = (transformed) ? GCGconshdlrDecompGetDetprobdataPresolved(scip) : GCGconshdlrDecompGetDetprobdataOrig(scip);

   /* shuffle dual multipliers of constraints*/

   /* 1) naive approach */
   if( usedmethod == GCG_RANDOM_DUAL_NAIVE )
   {
      for( int c = 0; c < nconss; ++c )
      {
         SCIP_Real dualval;
         SCIP_Real factor;
         SCIP_CONS* cons;

         cons = detprobdata->getCons(c);
         if( SCIPisInfinity( scip, -GCGconsGetLhs(scip, cons) ))
         {
            SCIP_Real modifier;
            modifier = 0.;
            if ( SCIPgetObjsense(scip) != SCIP_OBJSENSE_MAXIMIZE )
               modifier = -1.;

            factor = MAX(1., ABS(GCGconsGetRhs(scip, cons)  ) );
            dualval = SCIPrandomGetReal(randnumgen, 0.+modifier, 1.+modifier  ) * factor;
         }
         else if( SCIPisInfinity( scip, GCGconsGetRhs(scip, cons) ) )
         {
            SCIP_Real modifier;
            modifier = 0.;
            if ( SCIPgetObjsense(scip) != SCIP_OBJSENSE_MINIMIZE )
               modifier = -1.;

            factor = MAX(1., ABS(GCGconsGetLhs(scip, cons)  ) );
            dualval = SCIPrandomGetReal(randnumgen, 0.+modifier, 1.+modifier  ) * factor;
         }
         else
         {
            factor = MAX(1., ABS(GCGconsGetLhs(scip, cons)  ) );
            factor = MAX( factor, ABS(GCGconsGetRhs(scip, cons)   ) );
            dualval = SCIPrandomGetReal(randnumgen, -1., 1. ) * factor;
         }

         (*scoredata->dualvalsrandom)[c] = dualval;
      }

   } /* end naive approach */
   /* expected equal and expected overestimated approach */
   else if( usedmethod == GCG_RANDOM_DUAL_EXPECTED_EQUAL || usedmethod == GCG_RANDOM_DUAL_EXPECTED_OVERESTIMATE )
   {
      SCIP_Real largec  = 0.;
      for( int v = 0; v < SCIPgetNVars(scip); ++v )
         largec += ABS( SCIPvarGetObj(detprobdata->getVar(v)) );

      for( int c = 0; c < nconss; ++c )
      {
         SCIP_Real dualval;
         SCIP_CONS* cons;
         cons = detprobdata->getCons(c);
         double lambda;

         SCIP_Real divisor = 0.;

         SCIP_Real randomval;
         int nvarsincons = GCGconsGetNVars(scip, cons);
         SCIP_Real* valsincons = NULL;

         /* get values of variables in this constraint */
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &valsincons, nvarsincons) );
         GCGconsGetVals(scip, cons, valsincons, nvarsincons);

         /* add coefficients to divisor */
         for( int v = 0; v < nvarsincons; ++v )
         {
            divisor += ABS( valsincons[v] );
         }
         if ( usedmethod == GCG_RANDOM_DUAL_EXPECTED_EQUAL )
            divisor *= nconss;

         /* 1/lambda is the expected value of the distribution */
         lambda = divisor / largec;

         /* formula for exponential distribution requires a uniform random number in (0,1). */
         do
         {
            randomval = SCIPrandomGetReal(randnumgen, 0.0, 1.0);
         }
         while (randomval == 0.0 || randomval == 1.0);
         randomval = -std::log(randomval) / lambda;

         if( SCIPisInfinity(scip, -GCGconsGetLhs(scip, cons)) )
         {
            SCIP_Real modifier;
            modifier = 1.;
            if ( SCIPgetObjsense(scip) != SCIP_OBJSENSE_MAXIMIZE )
               modifier = -1.;

            dualval = modifier * randomval;
         }
         else if( SCIPisInfinity(scip, GCGconsGetRhs(scip, cons)) )
         {
            SCIP_Real modifier;
            modifier = 1.;
            if ( SCIPgetObjsense(scip) != SCIP_OBJSENSE_MINIMIZE )
               modifier = -1.;

            dualval = modifier * randomval;

         }
         else
         {
            SCIP_Real helpval = SCIPrandomGetReal(randnumgen, -1., 1. );

            if ( helpval < 0. )
               dualval =  -1. * randomval ;
            else dualval = randomval;
         }

         (*scoredata->dualvalsrandom)[c] = dualval;

         /* free storage for variables in cons */
         SCIPfreeBufferArrayNull(scip, &valsincons);
      }
   }

   return SCIP_OKAY;
}

/**
 * @brief return the a random value of the dual variable of the corresponding ; if it is not calculated yet it will be calculated
 * @returns the a random value of the dual variable of the corresponding
 */
static
SCIP_Real getDualvalRandom(
   SCIP* scip,             /**< SCIP data structure */
   GCG_SCORE* score,       /**< score */
   int  consindex,         /**< consindex  index of constraint the value is asked for */
   SCIP_Bool transformed   /**< is the problem transformed yet */
   )
{
   GCG_SCOREDATA* scoredata = GCGscoreGetData(score);

   if( !scoredata->dualvalsrandomset )
      shuffleDualvalsRandom(scip, score, transformed);
   scoredata->dualvalsrandomset = TRUE;

   return (*scoredata->dualvalsrandom)[consindex];
}

/**
 * @brief method to calculate and set the optimal dual values from original lp, used for strong detection score
 * @returns scip return code
 */
static
SCIP_RETCODE calculateDualvalsOptimalOrigLP(
   SCIP* scip,             /**< SCIP data structure */
   GCG_SCORE* score,       /**< score */
   SCIP_Bool transformed   /**< whether the problem is transormed yet */
   )
{
   GCG_SCOREDATA* scoredata = GCGscoreGetData(score);

   SCIP* scipcopy;
   SCIP_HASHMAP* origtocopiedconss;
   SCIP_Bool valid;
   int nvars;
   SCIP_VAR** copiedvars;
   int nconss;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "started calculating optimal dual values for original lp\n");

   SCIPhashmapCreate(&origtocopiedconss, SCIPblkmem(scip), SCIPgetNConss(scip));

   SCIPcreate(&scipcopy);

   SCIPcopy(scip, scipcopy, NULL, origtocopiedconss, "", FALSE, FALSE, FALSE, FALSE, &valid);

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scipcopy);
   copiedvars = SCIPgetVars(scipcopy);

   scoredata->dualvalsoptimaloriglp->clear();
   scoredata->dualvalsoptimaloriglp->resize(nconss, 0.);

   for( int var = 0; var < nvars ; ++var )
   {
      SCIP_VAR* copyvar = copiedvars[var];
      SCIP_Bool infeasible;

      if( SCIPvarGetType(copyvar) == SCIP_VARTYPE_BINARY )
         SCIPchgVarUbGlobal(scipcopy, copyvar, 1.);

      SCIPchgVarType(scipcopy, copyvar, SCIP_VARTYPE_CONTINUOUS, &infeasible);
   }

   copiedvars = SCIPgetVars(scipcopy);

   /* deactivate presolving */
   SCIPsetIntParam(scipcopy, "presolving/maxrounds", 0);

   /* deactivate separating */
   SCIPsetIntParam(scipcopy, "separating/maxrounds", 0);
   SCIPsetIntParam(scipcopy, "separating/maxroundsroot", 0);

   /* deactivate propagating */
   SCIPsetIntParam(scipcopy, "propagating/maxrounds", 0);
   SCIPsetIntParam(scipcopy, "propagating/maxroundsroot", 0);

   /* solve lp */
   SCIPsetIntParam(scipcopy, "lp/solvefreq", 1);

   /* only root node */
   SCIPsetLongintParam(scipcopy, "limits/nodes", 1);

   SCIPsetIntParam(scipcopy, "display/verblevel", SCIP_VERBLEVEL_FULL);

   SCIPtransformProb(scipcopy);

   SCIPsolve(scipcopy);

   gcg::DETPROBDATA* detprobdata = (transformed) ? GCGconshdlrDecompGetDetprobdataPresolved(scip) : GCGconshdlrDecompGetDetprobdataOrig(scip);

   for( int c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONS* copiedcons;
      SCIP_CONS* transcons = NULL;

      cons = detprobdata->getCons(c);
      if( !transformed )
      {
         SCIPgetTransformedCons(scip, cons, &transcons);
         if( transcons )
            cons = transcons;
         else
         {
            SCIPwarningMessage(scip, "Could not find constraint for random dual variable initialization when calculating strong decomposition score; skipping cons: %s \n", SCIPconsGetName(cons));
            continue;
         }
      }
      copiedcons = (SCIP_CONS*) SCIPhashmapGetImage(origtocopiedconss, (void*) cons);

      assert(copiedcons != NULL);
      assert( !SCIPconsIsTransformed(copiedcons) );

      transcons = NULL;
      SCIPgetTransformedCons(scipcopy, copiedcons, &transcons);
      if( transcons != NULL )
         copiedcons = transcons;

      (*scoredata->dualvalsoptimaloriglp)[c] = GCGconsGetDualsol(scipcopy, copiedcons);
      if( !SCIPisFeasEQ(scip, 0., (*scoredata->dualvalsoptimaloriglp)[c]) )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "optimal dual sol of constraint %s is %f \n", SCIPconsGetName(cons), (*scoredata->dualvalsoptimaloriglp)[c]);
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "finished calculating optimal dual values for original lp, start freeing\n");

   SCIPhashmapFree(&origtocopiedconss);
   SCIPfree(&scipcopy);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "finished freeing\n");

  return SCIP_OKAY;
}

/**
 * @brief returns the value of the optimal lp relaxation dual value of the given constrainr rid correspondoning problem of the detprobdata; if it is not calculated yet it will be calculated
 * @returns the value of the optimal lp relaxation dual value of the given constraint rid correspondoning problem of the detprobdata
 */
static
SCIP_Real getDualvalOptimalLP(
   SCIP* scip,             /**< SCIP data structure */
   GCG_SCORE* score,       /**< score */
   int  consindex,         /**< consindex index of constraint the value is asked for */
   SCIP_Bool transformed   /**< is the problem transformed yet */
   )
{
   GCG_SCOREDATA* scoredata = GCGscoreGetData(score);

   if( !scoredata->dualvalsoptimaloriglpcalculated )
      calculateDualvalsOptimalOrigLP(scip, score, transformed);
   scoredata->dualvalsoptimaloriglpcalculated = TRUE;

   return (*scoredata->dualvalsoptimaloriglp)[consindex];
}

/** @brief sets the pricing problem parameters 
 * @returns scip return code
*/
static
SCIP_RETCODE setTestpricingProblemParameters(
   SCIP*                 scip,               /**< SCIP data structure of the pricing problem */
   int                   clocktype,          /**< clocktype to use in the pricing problem */
   SCIP_Real             infinity,           /**< values larger than this are considered infinity in the pricing problem */
   SCIP_Real             epsilon,            /**< absolute values smaller than this are considered zero in the pricing problem */
   SCIP_Real             sumepsilon,         /**< absolute values of sums smaller than this are considered zero in the pricing problem */
   SCIP_Real             feastol,            /**< feasibility tolerance for constraints in the pricing problem */
   SCIP_Real             lpfeastolfactor,    /**< primal feasibility tolerance factor of LP solver in the pricing problem */
   SCIP_Real             dualfeastol,        /**< feasibility tolerance for reduced costs in LP solution in the pricing problem */
   SCIP_Bool             enableppcuts,       /**< should ppcuts be stored for sepa_basis */
   SCIP_Real             timelimit           /**< limit of time */
   )
{
   assert(scip != NULL);

   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(scip, "conflict/useprop", FALSE) );
   SCIP_CALL( SCIPsetCharParam(scip, "conflict/useinflp", 'o') );
   SCIP_CALL( SCIPsetCharParam(scip, "conflict/useboundlp", 'o') );
   SCIP_CALL( SCIPsetBoolParam(scip, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "conflict/usepseudo", FALSE) );

   /* reduce the effort spent for hash tables */
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/usevartable", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/useconstable", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/usesmalltables", TRUE) );

   /* disable expensive presolving */
   /* @todo test whether this really helps, perhaps set presolving emphasis to fast? */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/linear/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/setppc/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/logicor/presolpairwise", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/linear/presolusehashing", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/setppc/presolusehashing", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/logicor/presolusehashing", FALSE) );

   /* disable dual fixing presolver for the moment, because we want to avoid variables fixed to infinity */
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/dualfix/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/dualfix/maxprerounds", 0) );
   SCIP_CALL( SCIPfixParam(scip, "propagating/dualfix/freq") );
   SCIP_CALL( SCIPfixParam(scip, "propagating/dualfix/maxprerounds") );

   /* disable solution storage ! */
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", 0) );
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit ) );

   /* disable multiaggregation because of infinite values */
   SCIP_CALL( SCIPsetBoolParam(scip, "presolving/donotmultaggr", TRUE) );

   /* @todo enable presolving and propagation of xor constraints if bug is fixed */

   /* disable presolving and propagation of xor constraints as work-around for a SCIP bug */
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/xor/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/xor/propfreq", -1) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", (int)SCIP_VERBLEVEL_NORMAL) );
#if SCIP_VERSION > 210
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/printreason", FALSE) );
#endif
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", 0) );
   SCIP_CALL( SCIPfixParam(scip, "limits/maxorigsol") );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/catchctrlc", FALSE) );

   /* set clock type */
   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype", clocktype) );

   SCIP_CALL( SCIPsetBoolParam(scip, "misc/calcintegral", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/finitesolutionstore", TRUE) );

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/infinity", infinity) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", epsilon) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", sumepsilon) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", feastol) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/lpfeastolfactor", lpfeastolfactor) ); // canged from "numerics/lpfeastol"
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", dualfeastol) );

   /* jonas' stuff */
   if( enableppcuts )
   {
      int pscost;
      int prop;

      SCIP_CALL( SCIPgetIntParam(scip, "branching/pscost/priority", &pscost) );
      SCIP_CALL( SCIPgetIntParam(scip, "propagating/maxroundsroot", &prop) );
      SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", 11000) );
      SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );
      SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   }
   return SCIP_OKAY;
}

/** @brief creates the pricing problem constraints
 * @returns scip return code
 */
static
SCIP_RETCODE createTestPricingprobConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< the relaxator data data structure */
   gcg::PARTIALDECOMP*   partialdec,         /**< partialdec corresponding to the decomposition to test */
   int                   block,              /**< which pricing problem */
   SCIP_HASHMAP*         hashorig2pricingvar /**< hashmap mapping original to corresponding pricing variables */
   )
{
   SCIP_CONS* newcons;
   SCIP_HASHMAP* hashorig2pricingconstmp;
   int c;
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool success;

   assert(scip != NULL);

   gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();

   SCIP_CALL( SCIPhashmapCreate(&hashorig2pricingconstmp, SCIPblkmem(scip), detprobdata->getNConss() ) ); /*lint !e613*/

   assert(hashorig2pricingvar != NULL);

   for( c = 0; c < partialdec->getNConssForBlock(block); ++c )
   {
      SCIP_CONS* cons;

      cons = detprobdata->getCons(partialdec->getConssForBlock(block)[c]);

      SCIPdebugMessage("copying %s to pricing problem %d\n", SCIPconsGetName(cons), block);
      if( !SCIPconsIsActive(cons) )
      {
         SCIPdebugMessage("skipping, cons <%s> inactive\n", SCIPconsGetName(cons) );
         continue;
      }
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &cons) );
      assert(cons != NULL);

      /* copy the constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p%d_%s", block, SCIPconsGetName(cons));
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, cons, &newcons, SCIPconsGetHdlr(cons),
         hashorig2pricingvar, hashorig2pricingconstmp, name,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, &success) );

      /* constraint was successfully copied */
      assert(success);

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
#ifndef NDEBUG
      {
         SCIP_VAR** curvars;
         int ncurvars;

         ncurvars = GCGconsGetNVars(subscip, newcons);
         curvars = NULL;
         if( ncurvars > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
            SCIP_CALL( GCGconsGetVars(subscip, newcons, curvars, ncurvars) );

            SCIPfreeBufferArrayNull(scip, &curvars);
         }
      }
#endif
SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }

   SCIPhashmapFree(&hashorig2pricingconstmp);

   return SCIP_OKAY;
}




/*
 * score callback methods
 */

/** destructor of score to free user data (called when GCG is exiting) */
static
GCG_DECL_SCOREFREE(scoreFreeStrong)
{
   GCG_SCOREDATA* scoredata;

   assert(scip != NULL);

   scoredata = GCGscoreGetData(score);
   assert(scoredata != NULL);
   assert(strcmp(GCGscoreGetName(score), SCORE_NAME) == 0);

   delete scoredata->dualvalsoptimaloriglp;
   delete scoredata->dualvalsrandom;

   SCIPfreeMemory(scip, &scoredata);

   return SCIP_OKAY;
}

static
GCG_DECL_SCORECALC(scoreCalcStrong)
{
   /** @todo use and introduce scip parameter limit (for a pricing problem to be considered fractional solvable) of difference optimal value of LP-Relaxation and optimal value of artificial pricing problem */
   /* SCIP_Real gaplimitsolved; */

   /** @todo use and introduce scip parameter weighted limit (for a pricing problem to be considered fractional solvable) difference optimal value of LP-Relaxation and optimal value of artificial pricing problem */
   /* SCIP_Real gaplimitbeneficial; */

   GCG_SCOREDATA* scoredata = GCGscoreGetData(score);

   SCIP_Bool hittimelimit;
   SCIP_Bool errorpricing;
   int npricingconss = 0;

   SCIP_Real infinity;
   SCIP_Real epsilon;
   SCIP_Real sumepsilon;
   SCIP_Real feastol;
   SCIP_Real lpfeastolfactor;
   SCIP_Real dualfeastol;
   SCIP_Bool enableppcuts;

   SCIP_Bool benefical;
   SCIP_Bool fast;

   int clocktype;
   SCIP_Real dualvalmethodcoef;

   /* score works only on presolved  */
   gcg::PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);
   if( partialdec->isAssignedToOrigProb() )
   {
      *scorevalue = 0;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, " \n Attention! Strong decomposition score is not implemented for decomps belonging to the original problem \n\n");
      return SCIP_OKAY;
   }

   *scorevalue = 0.;

   /* ***** get all relevant parameters ***** */
   SCIPgetRealParam(scip, "detection/scores/strong/coeffactororigvsrandom", &dualvalmethodcoef);

   /* get numerical tolerances of the original SCIP instance in order to use the same numerical tolerances in master and pricing problems */
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/infinity", &infinity) );
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon) );
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/sumepsilon", &sumepsilon) );
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/feastol", &feastol) );
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/lpfeastolfactor", &lpfeastolfactor) ); // changed from "numerics/lpfeastol"
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/dualfeastol", &dualfeastol) );

   /* get clocktype of the original SCIP instance in order to use the same clocktype in master and pricing problems */
   SCIP_CALL( SCIPgetIntParam(scip, "timing/clocktype", &clocktype) );

   enableppcuts = FALSE;
   SCIP_CALL( SCIPgetBoolParam(scip, "sepa/basis/enableppcuts", &enableppcuts) );

   SCIP_Real timelimitfast  =  0.1 * scoredata->strongtimelimit;

   /* get number of pricingconss */
   for( int block = 0; block < partialdec->getNBlocks(); ++block )
   {
      npricingconss += partialdec->getNConssForBlock(block);
   }

   /* for every pricing problem calculate a corresponding score coeff and break if a pricing problem cannot be solved in the timelimit */
   for( int block = 0; block < partialdec->getNBlocks(); ++block )
   {
      SCIP* subscip;
      char name[SCIP_MAXSTRLEN];
      SCIP_HASHMAP* hashpricingvartoindex;
      SCIP_HASHMAP* hashorig2pricingvar;
      SCIP_Real score_coef;
      SCIP_Real weight_subproblem;
      std::stringstream subname;

      /* init all parameters */
      hittimelimit = FALSE;
      errorpricing = FALSE;

      subname << "temp_pp_" << block << ".lp";

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "started calculate strong decomposition subproblem for block %d \n", block );

      std::vector<SCIP_VAR*> indextopricingvar = std::vector<SCIP_VAR*>(SCIPgetNVars(scip), NULL);

      SCIP_CALL( SCIPhashmapCreate(&hashpricingvartoindex, SCIPblkmem(scip), SCIPgetNVars(scip)) ); /*lint !e613*/
      SCIP_CALL( SCIPhashmapCreate(&hashorig2pricingvar, SCIPblkmem(scip), SCIPgetNVars(scip)) ); /*lint !e613*/

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "testpricing_block_%d", block);

      benefical = FALSE;
      fast = FALSE;
      score_coef = 0.0;
      weight_subproblem = (SCIP_Real) partialdec->getNConssForBlock(block) / npricingconss;

      /* build subscip */
      SCIP_CALL( SCIPcreate(&subscip) );
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
      SCIP_CALL( setTestpricingProblemParameters(subscip, clocktype, infinity, epsilon, sumepsilon, feastol, lpfeastolfactor, dualfeastol, enableppcuts, scoredata->strongtimelimit) );
      SCIP_CALL( SCIPsetIntParam(subscip, "presolving/maxrounds", 0) );
      SCIP_CALL( SCIPsetIntParam(subscip, "lp/solvefreq", 1) );
      SCIP_CALL( SCIPcreateProb(subscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "started calculate strong decomposition, timelimit: %f  timelimitfast: %f \n",  scoredata->strongtimelimit, timelimitfast );

      /* copy variables */
      gcg::DETPROBDATA* detprobdata = partialdec->getDetprobdata();
      for( int var = 0; var < partialdec->getNVarsForBlock(block); ++var )
      {
         int varid;
         SCIP_VAR* origprobvar;
         SCIP_VAR* pricingprobvar;
         SCIP_Real obj;

         varid = partialdec->getVarsForBlock(block)[var];

         if ( partialdec->isAssignedToOrigProb() )
            origprobvar = detprobdata->getVar(varid);
         else
            origprobvar = SCIPvarGetProbvar(detprobdata->getVar(varid));

         /* calculate obj val from shuffled */
         obj = SCIPvarGetObj(origprobvar);
         for( int c = 0; c < detprobdata->getNConssForVar(varid); ++c )
         {
            int consid;
            SCIP_Real dualval;

            dualval = 0.;

            consid = detprobdata->getConssForVar(varid)[c];
            if ( partialdec->isConsMastercons(consid) )
            {
               if( SCIPisEQ( scip, dualvalmethodcoef, 0.0) )
                  dualval = getDualvalRandom(scip, score, consid, partialdec->isAssignedToOrigProb());
               else if( SCIPisEQ( scip, dualvalmethodcoef, 1.0) )
                  dualval = getDualvalOptimalLP(scip, score, consid, partialdec->isAssignedToOrigProb());
               else
                  dualval = dualvalmethodcoef * getDualvalOptimalLP(scip, score, consid, partialdec->isAssignedToOrigProb()) + (1. - dualvalmethodcoef) * getDualvalRandom(scip, score, consid, partialdec->isAssignedToOrigProb());
               obj -= dualval * detprobdata->getVal(consid, varid);
            }
         }

         /* round variable objective coeffs to decrease numerical troubles */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pr%d_%s", block, SCIPvarGetName(origprobvar));
         SCIP_CALL( SCIPcreateVar(subscip, &pricingprobvar, name, SCIPvarGetLbGlobal(origprobvar),
                  SCIPvarGetUbGlobal(origprobvar), obj, SCIPvarGetType(origprobvar),
                  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIPhashmapSetImage(hashorig2pricingvar, origprobvar, pricingprobvar);
         SCIPhashmapSetImage(hashpricingvartoindex, pricingprobvar, (void*) (size_t)varid);
         indextopricingvar[varid] = pricingprobvar;
         SCIP_CALL( SCIPaddVar(subscip, pricingprobvar) );
      }

      /* copy constraints */
      SCIP_CALL( createTestPricingprobConss(scip, subscip, partialdec, block, hashorig2pricingvar) );

      /* transform subscip */
      SCIP_CALL(SCIPtransformProb(subscip) );

      /* solve subscip */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "started solving subproblem for block %d \n", block );
      SCIP_CALL( SCIPsolve(subscip) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "finished solving subproblem in %f seconds \n", SCIPgetSolvingTime(subscip) );

      if( SCIPgetStatus(subscip) != SCIP_STATUS_OPTIMAL )
      {
         if( SCIPgetStatus(subscip) == SCIP_STATUS_TIMELIMIT )
            hittimelimit = TRUE;
         else
            errorpricing = TRUE;
      }

      if( errorpricing || hittimelimit )
      {
         if( hittimelimit )
            SCIPverbMessage(scip,  SCIP_VERBLEVEL_FULL, NULL, "Hit timelimit in pricing problem %d \n.", block);
         else
            SCIPverbMessage(scip,  SCIP_VERBLEVEL_FULL, NULL, "Error in pricing problem %d \n.", block);

         *scorevalue = 0.;
         SCIPhashmapFree(&hashpricingvartoindex);
         SCIPfree(&subscip);

         return SCIP_OKAY;
      }

      /* get coefficient */
      if( !SCIPisEQ( scip,  SCIPgetFirstLPLowerboundRoot(subscip), SCIPgetDualbound(subscip) ) )
         benefical = TRUE;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "first dual bound: %f ; dual bound end: %f \n",SCIPgetFirstLPLowerboundRoot(subscip), SCIPgetDualbound(subscip)   );

      if( SCIPisFeasLE( scip, SCIPgetSolvingTime(subscip), timelimitfast ) )
         fast = TRUE;

      if ( fast && benefical )
         score_coef = DEFAULT_SCORECOEF_FASTBENEFICIAL;

      if ( !fast && benefical )
         score_coef = DEFAULT_SCORECOEF_MEDIUMBENEFICIAL;

      if ( fast && !benefical )
         score_coef = DEFAULT_SCORECOEF_FASTNOTBENEFICIAL;

      if ( !fast && !benefical )
         score_coef = DEFAULT_SCORECOEF_MEDIUMNOTBENEFICIAL;

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "scorecoef for subproblem %d is %f with weighting factor %f\n", block, score_coef, weight_subproblem );

      *scorevalue += score_coef * weight_subproblem;

      /* free stuff */
      for( int var = 0; var < partialdec->getNVarsForBlock(block); ++var )
      {
         int varid = partialdec->getVarsForBlock(block)[var];
         SCIPreleaseVar(subscip, &indextopricingvar[varid]);
      }

      SCIPhashmapFree(&hashpricingvartoindex);

      SCIPfree(&subscip);
   }// end for blocks

   return SCIP_OKAY;
}




/*
 * score specific interface methods
 */

/** creates the strong decomposition score and includes it in SCIP */
SCIP_RETCODE GCGincludeScoreStrongDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   GCG_SCOREDATA* scoredata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &scoredata) );

   scoredata->dualvalsrandom = new std::vector<SCIP_Real>();
   scoredata->dualvalsoptimaloriglp = new std::vector<SCIP_Real>();
   scoredata->dualvalsrandomset = FALSE;
   scoredata->dualvalsoptimaloriglpcalculated = FALSE;

   assert(scoredata != NULL);

   SCIP_CALL(
      GCGincludeScore(scip, SCORE_NAME, SCORE_SHORTNAME, SCORE_DESC, scoredata,
         scoreFreeStrong, scoreCalcStrong) );

   SCIP_CALL( SCIPaddRealParam(scip, "detection/scores/strong/timelimit",
                               "Timelimit for strong decompositions score calculation per partialdec in seconds", &scoredata->strongtimelimit, FALSE,
                               DEFAULT_STRONGTIMELIMIT, 0., INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/scores/strong/dualvalrandommethod",
                              "Method for random dual values use for strong decomposition: 1: naive, 2: expected equality exponential distributed, 3: expected overestimation exponential distributed ",
                              &scoredata->strongdetectiondualvalrandommethod, FALSE,
                              DEFAULT_DUALVALRANDOMMETHOD, 1, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "detection/scores/strong/coeffactororigvsrandom",
                               "Convex coefficient for orig dual val, i.e. (1-this coef) is factor for random dual value", &scoredata->coeffactororigvsrandom, FALSE,
                               DEFAULT_COEFFACTORORIGVSRANDOM, 0., 1., NULL, NULL) );

   return SCIP_OKAY;
}
