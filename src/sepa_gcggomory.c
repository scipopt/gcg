/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
/**@file   sepa_gcggomory.c
 * @ingroup SEPARATORS
 * @brief  Gomory MIR Cuts
 * @author Tobias Achterberg
 */

/**@todo try k-Gomory-cuts (s. Cornuejols: K-Cuts: A Variation of Gomory Mixed Integer Cuts from the LP Tableau) */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "sepa_gcggomory.h"
#include "relax_gcg.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "gcggomory"
#define SEPA_DESC              "Gcggomory MIR cuts separator"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of gcggomory separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of gcggomory separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of gcggomory cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of gcggomory cuts separated per separation round in root node */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAXWEIGHTRANGE    1e+04 /**< maximal valid range max(|weights|)/min(|weights|) of row weights */

#define MAKECUTINTEGRAL        /* try to scale all cuts to integral coefficients */
/*#define MAKEINTCUTINTEGRAL*/     /* try to scale cuts without continuous variables to integral coefficients */
#define FORCECUTINTEGRAL       /* discard cut if conversion to integral coefficients failed */
//#define SEPARATEROWS           /* separate rows with integral slack */

#define BOUNDSWITCH              0.9999
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE
#define FIXINTEGRALRHS            FALSE
#define MAKECONTINTEGRAL          FALSE
#define MINFRAC                    0.05
#define MAXFRAC                    0.95

#define MAXAGGRLEN(nvars)          (0.1*(nvars)+1000) /**< maximal length of base inequality */


/** separator data */
struct SCIP_SepaData
{
   SCIP_Real             maxweightrange;     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   int                   maxrounds;          /**< maximal number of gcggomory separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of gcggomory separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of gcggomory cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of gcggomory cuts separated per separation round in root node */
   int                   lastncutsfound;     /**< total number of cuts found after last call of separator */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};




/*
 * local methods
 */

/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   SCIP_VAR**            cutvars,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real val;
   SCIP_Real absval;
   SCIP_Real cutsqrnorm;
   SCIP_Real act;
   SCIP_Real norm;
   int len;
   int v;

   assert(nvars == 0 || cutcoefs != NULL);
   assert(nvars == 0 || varsolvals != NULL);
   assert(cutvars != NULL);
   assert(cutvals != NULL);
   assert(cutlen != NULL);
   assert(cutact != NULL);
   assert(cutnorm != NULL);

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch( normtype )
   {
   case 'e':
      cutsqrnorm = 0.0;
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            absval = REALABS(val);
            norm = MAX(norm, absval);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 's':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 'd':
      for( v = 0; v < nvars; ++v )
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm = 1.0;
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}




/*
 * Callback methods
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeGcggomory)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called when problem solving starts) */
#define sepaInitGcggomory NULL


/** deinitialization method of separator (called when problem solving exits) */
#define sepaExitGcggomory NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolGcggomory NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolGcggomory NULL


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolGcggomory)
{  /*lint --e{715}*/
   SCIP* masterscip;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** mastervars;
   SCIP_COL** mastercols;
   SCIP_ROW** masterrows;
   SCIP_VAR** origvars;
   SCIP_COL** origcols;
   SCIP_ROW** origrows;
   SCIP_Real* varsolvals;
   SCIP_Real* binvrow;
   SCIP_Real* cutcoefs;
   SCIP_Real cutrhs;
   SCIP_Real cutact;
   SCIP_Real maxscale;
   SCIP_Longint maxdnom;
   int* basisind;
   int nmastervars;
   int nmastercols;
   int nmasterrows;
   int norigvars;
   int norigcols;
   int norigrows;
   int ncalls;
   int depth;
   int maxdepth;
   int maxsepacuts;
   int ncuts;
   int c;
   int i;
   SCIP_Bool success;
   SCIP_Bool cutislocal;
   char normtype;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   masterscip = GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   printf("sepa_gcggomory call %d at current node!\n", ncalls);

   /* only call the gcggomory cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(masterscip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if the LP solution is basic */
   if( !SCIPisLPSolBasic(masterscip) )
      return SCIP_OKAY;

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &origvars, &norigvars, NULL, NULL, NULL, NULL) );

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(masterscip, &mastercols, &nmastercols) );
   SCIP_CALL( SCIPgetLPRowsData(masterscip, &masterrows, &nmasterrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &origcols, &norigcols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &origrows, &norigrows) );

   if( nmastercols == 0 || nmasterrows == 0 || norigcols == 0 || norigrows == 0 )
      return SCIP_OKAY;


   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* set the maximal denominator in rational representation of gcggomory cut and the maximal scale factor to
    * scale resulting cut to integral values to avoid numerical instabilities
    */
   /**@todo find better but still stable gcggomory cut settings: look at dcmulti, gesa3, khb0525, misc06, p2756 */
   maxdepth = SCIPgetMaxDepth(scip);
   if( depth == 0 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/4 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/2 )
   {
      maxdnom = 100;
      maxscale = 100.0;
   }
   else
   {
      maxdnom = 10;
      maxscale = 10.0;
   }

   *result = SCIP_DIDNOTFIND;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nmasterrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nmasterrows) );
   varsolvals = NULL; /* allocate this later, if needed */

   /* get basis indices */
   SCIP_CALL( SCIPgetLPBasisInd(masterscip, basisind) );

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
      maxsepacuts = sepadata->maxsepacutsroot;
   else
      maxsepacuts = sepadata->maxsepacuts;

   SCIPdebugMessage("searching gcggomory cuts: %d cols, %d rows, maxdnom=%"SCIP_LONGINT_FORMAT", maxscale=%g, maxcuts=%d\n",
      nmastercols, nmasterrows, maxdnom, maxscale, maxsepacuts);

   /* for all basic columns belonging to integer variables, try to generate a gcggomory cut */
   ncuts = 0;
   for( i = 0; i < nmasterrows && ncuts < maxsepacuts && !SCIPisStopped(scip); ++i )
   {
      SCIP_Bool tryrow;

      tryrow = FALSE;
      c = basisind[i];
      if( c >= 0 )
      {
         SCIP_VAR* var;

         assert(c < nmastercols);
         var = SCIPcolGetVar(mastercols[c]);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_Real primsol;

            primsol = SCIPcolGetPrimsol(mastercols[c]);
            assert(SCIPgetVarSol(masterscip, var) == primsol); /*lint !e777*/

            if( SCIPfeasFrac(scip, primsol) >= MINFRAC )
            {
               SCIPdebugMessage("trying gcggomory cut for col <%s> [%g]\n", SCIPvarGetName(var), primsol);
               tryrow = TRUE;
            }
         }
      }
#ifdef SEPARATEROWS
      else
      {
         SCIP_ROW* row;

         assert(0 <= -c-1 && -c-1 < nrows); 
         row = rows[-c-1];
         if( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) )
         {
            SCIP_Real primsol;

            primsol = SCIPgetRowActivity(scip, row);
            if( SCIPfeasFrac(scip, primsol) >= MINFRAC )
            {
               SCIPdebugMessage("trying gcggomory cut for row <%s> [%g]\n", SCIProwGetName(row), primsol);
               tryrow = TRUE;
            }
         }
      }
#endif

      if( tryrow )
      {
         int j;
         int k;
         SCIP_Real* weights;
         char rowname[SCIP_MAXSTRLEN];
         SCIP_Bool allweights;

         SCIP_CALL( SCIPallocBufferArray(scip, &weights, norigrows));
         BMSclearMemoryArray(weights, norigrows);

         /* get the row of B^-1 for this basic integer variable with fractional solution value */
         SCIP_CALL( SCIPgetLPBInvRow(masterscip, i, binvrow) );

         allweights = TRUE;

         for ( j = 0; j < nmasterrows; j++ )
         {
            if ( SCIPisZero(scip, binvrow[j]) )
               continue;
            //printf("masterrow <%s>\n", SCIProwGetName(masterrows[j]));
            for ( k = 0; k < norigrows; k++ )
            {
               (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "m_%s", SCIProwGetName(origrows[k]));
               //printf("origrow <%s>, rowname <%s>\n", SCIProwGetName(origrows[k]), rowname);
               if( strcmp(SCIProwGetName(masterrows[j]), rowname) == 0 )
               {
                  //printf("masterrow <%s> =  origrow <%s>\n", SCIProwGetName(masterrows[j]), rowname);
                  weights[k] = binvrow[j];
                  break;
               }
            }
            if ( k == norigrows )
            {
               allweights = FALSE;
               printf("masterrow <%s> weight = %g\n", SCIProwGetName(masterrows[j]), binvrow[j]);
            }
         }

         /* create a MIR cut out of the weighted LP rows using the B^-1 row as weights */
         SCIP_CALL( SCIPcalcMIR(scip, sol, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, FIXINTEGRALRHS, NULL, NULL,
               (int) MAXAGGRLEN(nmastervars), sepadata->maxweightrange, MINFRAC, MAXFRAC,
               weights, 1.0, NULL, NULL, cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
         assert(ALLOWLOCAL || !cutislocal);
         SCIPdebugMessage(" -> success=%u: %g <= %g\n", success, cutact, cutrhs);

         SCIPfreeBufferArray(scip, &weights);

         /* if successful, convert dense cut into sparse row, and add the row as a cut */
         if( success /* && SCIPisFeasGT(scip, cutact, cutrhs) */ )
         {
            SCIP_VAR** cutvars;
            SCIP_Real* cutvals;
            SCIP_Real cutnorm;
            int cutlen;

            /* if this is the first successful cut, get the LP solution for all COLUMN variables */
            if( varsolvals == NULL )
            {
               int v;

               SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, norigvars) );
               for( v = 0; v < norigvars; ++v )
               {
                  if( SCIPvarGetStatus(origvars[v]) == SCIP_VARSTATUS_COLUMN )
                     varsolvals[v] = SCIPgetSolVal(scip, sol, origvars[v]);
               }
            }
            assert(varsolvals != NULL);

            /* get temporary memory for storing the cut as sparse row */
            SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, norigvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, norigvars) );

            /* store the cut as sparse row, calculate activity and norm of cut */
            SCIP_CALL( storeCutInArrays(scip, norigvars, origvars, cutcoefs, varsolvals, normtype,
                  cutvars, cutvals, &cutlen, &cutact, &cutnorm) );

            SCIPdebugMessage(" -> gcggomory cut for <%s>: act=%f, rhs=%f, norm=%f, eff=%f\n",
               c >= 0 ? SCIPvarGetName(SCIPcolGetVar(mastercols[c])) : SCIProwGetName(masterrows[-c-1]),
               cutact, cutrhs, cutnorm, (cutact - cutrhs)/cutnorm);

            /*if( SCIPisPositive(scip, cutnorm) && SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm) )*/
            {
               SCIP_ROW* cut;
               char cutname[SCIP_MAXSTRLEN];

               /* create the cut */
               if( c >= 0 )
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "gcggom%d_x%d", SCIPgetNLPs(masterscip), c);
               else
                  (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "gcggom%d_s%d", SCIPgetNLPs(masterscip), -c-1);
               SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), cutrhs, 
                                             cutislocal, FALSE, sepadata->dynamiccuts) );
               SCIP_CALL( SCIPaddVarsToRow(scip, cut, cutlen, cutvars, cutvals) );
               /*SCIPdebug(SCIPprintRow(scip, cut, NULL));*/

               assert(success);
#ifdef MAKECUTINTEGRAL
               /* try to scale the cut to integral values */
               SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                     maxdnom, maxscale, MAKECONTINTEGRAL, &success) );
#else
#ifdef MAKEINTCUTINTEGRAL
               /* try to scale the cut to integral values if there are no continuous variables
                *  -> leads to an integral slack variable that can later be used for other cuts
                */
               {
                  int k;
                  for( k = 0; k < cutlen && SCIPvarIsIntegral(cutvars[k]); ++k )
                  {}
                  if( k == cutlen )
                  {
                     SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                           maxdnom, maxscale, MAKECONTINTEGRAL, &success) );
                  }
               }
#endif
#endif

#ifndef FORCECUTINTEGRAL
               success = TRUE;
#endif

               if( success )
               {
#if 0
                  if( !SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIPdebugMessage(" -> gcggomory cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                        cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                        SCIPgetCutEfficacy(scip, NULL, cut));
                     /*SCIPdebug(SCIPprintRow(scip, cut, NULL));*/
                     success = FALSE;
                  }
                  else
#endif
                  {
                     SCIPdebugMessage(" -> found gcggomory cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                        cutname, SCIPgetRowSolActivity(scip, cut, sol), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                        SCIPgetCutEfficacy(scip, sol, cut),
                        SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                        SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
                     /*SCIPdebug(SCIPprintRow(scip, cut, NULL));*/
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, TRUE) );
                     if( !cutislocal )
                     {
                        SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                     }
                     *result = SCIP_SEPARATED;
                     ncuts++;
                  }
               }
               else
               {
                  SCIPdebugMessage(" -> gcggomory cut <%s> couldn't be scaled to integral coefficients: act=%f, rhs=%f, norm=%f, eff=%f\n",
                     cutname, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, sol, cut));
               }

               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }

            /* free temporary memory */
            SCIPfreeBufferArray(scip, &cutvals);
            SCIPfreeBufferArray(scip, &cutvars);
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArrayNull(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutcoefs);

   SCIPdebugMessage("end searching gcggomory cuts: found %d cuts\n", ncuts);

   sepadata->lastncutsfound = SCIPgetNCutsFound(scip);

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
#define sepaExeclpGcggomory NULL /* gcggomory cuts need a basic LP solution */




/*
 * separator specific interface methods
 */

/** creates the Gcggomory MIR cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaGcggomory(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
   sepadata->lastncutsfound = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeGcggomory, sepaInitGcggomory, sepaExitGcggomory,
         sepaInitsolGcggomory, sepaExitsolGcggomory, 
         sepaExeclpGcggomory, sepaExecsolGcggomory,
         sepadata) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gcggomory/maxrounds",
         "maximal number of gcggomory separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gcggomory/maxroundsroot",
         "maximal number of gcggomory separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gcggomory/maxsepacuts",
         "maximal number of gcggomory cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gcggomory/maxsepacutsroot",
         "maximal number of gcggomory cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gcggomory/maxweightrange",
         "maximal valid range max(|weights|)/min(|weights|) of row weights",
         &sepadata->maxweightrange, TRUE, DEFAULT_MAXWEIGHTRANGE, 1.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gcggomory/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

