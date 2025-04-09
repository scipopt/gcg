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

/**@file   clscons_miplibconstypes.cpp
 * 
 * @brief classifies constraints according to their miplib constypes
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/clscons_miplibconstypes.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include <vector>
#include <stdio.h>
#include <sstream>

#include "gcg/class_detprobdata.h"

#include "gcg/class_conspartition.h"
#include "gcg/scip_misc.h"

/* classifier properties */
#define CLSCONS_NAME                  "miplibconstype"       /**< name of classifier */
#define CLSCONS_DESC                  "miplib constypes"     /**< short description of classification*/
#define CLSCONS_PRIORITY              0

#define CLSCONS_ENABLED               TRUE


/*
 * Data structures
 */

/** classifier handler data */
struct GCG_ClassifierData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * classifier callback methods
 */

/** destructor of classifier to free user data (called when GCG is exiting) */
#define classifierFree NULL

static
GCG_DECL_CONSCLASSIFY(classifierClassify)
{
   gcg::DETPROBDATA* detprobdata;
   SCIP* origprob = GCGgetOrigprob(gcg);
   if( transformed )
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(gcg);
   }
   else
   {
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(gcg);
   }

   std::vector<int> nfoundconstypesrangedsinglecount( (int) SCIP_CONSTYPE_GENERAL + 1, 0 );
   std::vector<int> nfoundconstypesrangeddoublecount( (int) SCIP_CONSTYPE_GENERAL + 1, 0 );

   std::vector<int> classforcons = std::vector<int>( detprobdata->getNConss(), -1 );
   gcg::ConsPartition* classifier;

   /* firstly, assign all constraints to classindices */
   for( int c = 0; c < detprobdata->getNConss(); ++ c )
   {
      SCIP_CONS* cons;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real* vals = NULL;
      SCIP_VAR** vars = NULL;
      int nvars;
      int i;

      cons = detprobdata->getCons(c);

      nvars =  GCGconsGetNVars(origprob, cons );

      lhs = GCGconsGetLhs(origprob, cons);
      rhs = GCGconsGetRhs(origprob, cons);
      if( nvars != 0 )
      {
         SCIP_CALL_ABORT( SCIPallocBufferArray(origprob, &vals, nvars));
         SCIP_CALL_ABORT( SCIPallocBufferArray(origprob, &vars, nvars));
         SCIP_CALL_ABORT( GCGconsGetVals(origprob, cons, vals, nvars ) );
         SCIP_CALL_ABORT( GCGconsGetVars(origprob, cons, vars, nvars ) );
      }

      for( i = 0; i < nvars; i++ )
      {
         assert(!SCIPisZero(origprob, vals[i]) );
      }


      /* is constraint of type SCIP_CONSTYPE_EMPTY? */
      if( nvars == 0 )
      {
         SCIPdebugMsg(origprob, "classified as EMPTY: ");
         SCIPdebugPrintCons(origprob, cons, NULL);
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_EMPTY]++;
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_EMPTY]++;
         classforcons[c] = SCIP_CONSTYPE_EMPTY;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_FREE? */
      if( SCIPisInfinity(origprob, rhs) && SCIPisInfinity(origprob, -lhs) )
      {
         SCIPdebugMsg(origprob, "classified as FREE: ");
         SCIPdebugPrintCons(origprob, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_FREE]++;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_FREE]++;
         classforcons[c] = SCIP_CONSTYPE_FREE;
         SCIPfreeBufferArray(origprob, &vars);
         SCIPfreeBufferArray(origprob, &vals);
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_SINGLETON? */
      if( nvars == 1 )
      {
         SCIPdebugMsg(origprob, "classified as SINGLETON: ");
         SCIPdebugPrintCons(origprob, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SINGLETON] += 2 ;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SINGLETON]++;
         classforcons[c] = SCIP_CONSTYPE_SINGLETON;
         SCIPfreeBufferArray(origprob, &vars) ;
         SCIPfreeBufferArray(origprob, &vals) ;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_AGGREGATION? */
      if( nvars == 2 && SCIPisEQ(origprob, lhs, rhs) )
      {
         SCIPdebugMsg(origprob, "classified as AGGREGATION: ");
         SCIPdebugPrintCons(origprob, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_AGGREGATION]++;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_AGGREGATION]++;
         classforcons[c] = SCIP_CONSTYPE_AGGREGATION;
         SCIPfreeBufferArray(origprob, &vars) ;
         SCIPfreeBufferArray(origprob, &vals) ;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{VARBOUND}? */
      if( nvars == 2 )
      {
         SCIPdebugMsg(origprob, "classified as VARBOUND: ");
         SCIPdebugPrintCons(origprob, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_VARBOUND] += 2 ;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_VARBOUND]++;
         classforcons[c] = SCIP_CONSTYPE_VARBOUND;
         SCIPfreeBufferArray(origprob, &vars) ;
         SCIPfreeBufferArray(origprob, &vals) ;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{SETPARTITION, SETPACKING, SETCOVERING, CARDINALITY, INVKNAPSACK}? */
      {
         SCIP_Real scale;
         SCIP_Real b;
         SCIP_Bool unmatched;
         int nnegbinvars;

         unmatched = FALSE;
         nnegbinvars = 0;

         scale = REALABS(vals[0]);
         for( i = 0; i < nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(origprob, SCIPvarGetLbGlobal(vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(origprob, SCIPvarGetUbGlobal(vars[i]), 2.0);
            unmatched = unmatched || !SCIPisEQ(origprob, REALABS(vals[i]), scale);

            if( vals[i] < 0.0 )
               nnegbinvars++;
         }

         if( !unmatched )
         {
            if( SCIPisEQ(origprob, lhs, rhs) )
            {
               b = rhs/scale + nnegbinvars;
               if( SCIPisEQ(origprob, 1.0, b) )
               {
                  SCIPdebugMsg(origprob, "classified as SETPARTITION: ");
                  SCIPdebugPrintCons(origprob, cons, NULL);
                  nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SETPARTITION] += 1 ;
                  nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SETPARTITION]++;
                  classforcons[c] = SCIP_CONSTYPE_SETPARTITION;
                  SCIPfreeBufferArray(origprob, &vars) ;
                  SCIPfreeBufferArray(origprob, &vals) ;
                  continue;
               }
               else if( SCIPisIntegral(origprob, b) && !SCIPisNegative(origprob, b) )
               {
                  SCIPdebugMsg(origprob, "classified as CARDINALITY: ");
                  SCIPdebugPrintCons(origprob, cons, NULL);
                  nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_CARDINALITY] += 1 ;
                  nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_CARDINALITY]++;
                  classforcons[c] = SCIP_CONSTYPE_CARDINALITY;
                  SCIPfreeBufferArray(origprob, &vars);
                  SCIPfreeBufferArray(origprob, &vals);
                  continue;
               }
            }

            b = rhs/scale + nnegbinvars;
            if( SCIPisEQ(origprob, 1.0, b) )
            {
               SCIPdebugMsg(origprob, "classified as SETPACKING: ");
               SCIPdebugPrintCons(origprob, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SETPACKING] += 1 ;
               nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SETPACKING]++;
               classforcons[c] = SCIP_CONSTYPE_SETPACKING;
               rhs = SCIPinfinity(origprob);
            }
            else if( SCIPisIntegral(origprob, b) && !SCIPisNegative(origprob, b) )
            {
               SCIPdebugMsg(origprob, "classified as INVKNAPSACK: ");
               SCIPdebugPrintCons(origprob, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_INVKNAPSACK] += 1 ;
                nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_INVKNAPSACK]++;
                classforcons[c] = SCIP_CONSTYPE_INVKNAPSACK;
               rhs = SCIPinfinity(origprob);
            }

            b = lhs/scale + nnegbinvars;
            if( SCIPisEQ(origprob, 1.0, b) )
            {
               SCIPdebugMsg(origprob, "classified as SETCOVERING: ");
               SCIPdebugPrintCons(origprob, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SETCOVERING] += 1 ;
               nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SETCOVERING]++;
               classforcons[c] = SCIP_CONSTYPE_SETCOVERING;
               lhs = -SCIPinfinity(origprob);
            }

            if( SCIPisInfinity(origprob, -lhs) && SCIPisInfinity(origprob, rhs) )
            {
               SCIPfreeBufferArray(origprob, &vars);
               SCIPfreeBufferArray(origprob, &vals);
               continue;
            }
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{EQKNAPSACK, BINPACKING, KNAPSACK}? */
      /* @todo If coefficients or rhs are not integral, we currently do not check
       * if the constraint could be scaled (finitely), such that they are.
       */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         b = rhs;
         unmatched = FALSE;
         for( i = 0; i < nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(origprob, SCIPvarGetLbGlobal(vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(origprob, SCIPvarGetUbGlobal(vars[i]), 2.0);
            unmatched = unmatched || !SCIPisIntegral(origprob, vals[i]);

            if( SCIPisNegative(origprob, vals[i]) )
               b -= vals[i];
         }
         unmatched = unmatched || !detprobdata->isFiniteNonnegativeIntegral(b);

         if( !unmatched )
         {
            if( SCIPisEQ(origprob, lhs, rhs) )
            {
               SCIPdebugMsg(origprob, "classified as EQKNAPSACK: ");
               SCIPdebugPrintCons(origprob, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_EQKNAPSACK] += 1 ;
               nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_EQKNAPSACK]++;
               classforcons[c] = SCIP_CONSTYPE_EQKNAPSACK;
               SCIPfreeBufferArray(origprob, &vars);
               SCIPfreeBufferArray(origprob, &vals);
               continue;
            }
            else
            {
               SCIP_Bool matched;

               matched = FALSE;
               for( i = 0; i < nvars && !matched; i++ )
               {
                  matched = matched || SCIPisEQ(origprob, b, REALABS(vals[i]));
               }

               SCIPdebugMsg(origprob, "classified as %s: ", matched ? "BINPACKING" : "KNAPSACK");
               SCIPdebugPrintCons(origprob, cons, NULL);
               nfoundconstypesrangeddoublecount[matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK] += 1 ;
               nfoundconstypesrangedsinglecount[matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK]++;
               classforcons[c] = matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK;

            }

            if( SCIPisInfinity(origprob, -lhs) )
            {
               SCIPfreeBufferArray(origprob, &vars);
               SCIPfreeBufferArray(origprob, &vals);
               continue;
            }
            else
               rhs = SCIPinfinity(origprob);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{INTKNAPSACK}? */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         unmatched = FALSE;

         b = rhs;
         unmatched = unmatched || !detprobdata->isFiniteNonnegativeIntegral(b);

         for( i = 0; i < nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisNegative(origprob, SCIPvarGetLbGlobal(vars[i]));
            unmatched = unmatched || !SCIPisIntegral(origprob, vals[i]);
            unmatched = unmatched || SCIPisNegative(origprob, vals[i]);
         }

         if( !unmatched )
         {
            SCIPdebugMsg(origprob, "classified as INTKNAPSACK: ");
            SCIPdebugPrintCons(origprob, cons, NULL);
            nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_INTKNAPSACK] += 1 ;
            nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_INTKNAPSACK]++;
            classforcons[c] = SCIP_CONSTYPE_INTKNAPSACK;

            if( SCIPisInfinity(origprob, -lhs) )
            {
               SCIPfreeBufferArray(origprob, &vars);
               SCIPfreeBufferArray(origprob, &vals);
               continue;
            }
            else
               rhs = SCIPinfinity(origprob);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{MIXEDBINARY}? */
      {
         SCIP_Bool unmatched;

         unmatched = FALSE;
         for( i = 0; i < nvars && !unmatched; i++ )
         {
            if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS
               && (SCIPisLE(origprob, SCIPvarGetLbGlobal(vars[i]), -1.0)
                  || SCIPisGE(origprob, SCIPvarGetUbGlobal(vars[i]), 2.0)) )
               unmatched = TRUE;
         }

         if( !unmatched )
         {
            SCIPdebugMsg(origprob, "classified as MIXEDBINARY (%d): ", detprobdata->isRangedRow(lhs, rhs) ? 2 : 1);
            SCIPdebugPrintCons(origprob, cons, NULL);
            nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_MIXEDBINARY] += 1 ;
            nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_MIXEDBINARY]++;
            classforcons[c] = SCIP_CONSTYPE_MIXEDBINARY;
            SCIPfreeBufferArray(origprob, &vars) ;
            SCIPfreeBufferArray(origprob, &vals) ;
            continue;

         }
      }

      /* no special structure detected */
      SCIPdebugMsg(origprob, "classified as GENERAL: ");
      SCIPdebugPrintCons(origprob, cons, NULL);
      nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_GENERAL] += 1 ;
      nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_GENERAL]++;
      classforcons[c] = SCIP_CONSTYPE_GENERAL;
      SCIPfreeBufferArray(origprob, &vars);
      SCIPfreeBufferArray(origprob, &vals);
   }

   classifier = new gcg::ConsPartition(gcg, "constypes according to miplib", (int) SCIP_CONSTYPE_GENERAL + 1, detprobdata->getNConss() );

   /* set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::string name;
      std::stringstream text;
      switch( c )
      {
         case (int) SCIP_CONSTYPE_EMPTY:
            name = "empty";
            break;
         case SCIP_CONSTYPE_FREE:
            name = "free";
            break;
         case SCIP_CONSTYPE_SINGLETON:
            name = "singleton";
            break;
         case SCIP_CONSTYPE_AGGREGATION:
            name = "aggregation";
            break;
         case SCIP_CONSTYPE_VARBOUND:
            name = "varbound";
            break;
         case SCIP_CONSTYPE_SETPARTITION:
            name = "setpartition";
            break;
         case SCIP_CONSTYPE_SETPACKING:
            name = "setpacking";
            break;
         case SCIP_CONSTYPE_SETCOVERING:
            name = "setcovering";
            break;
         case SCIP_CONSTYPE_CARDINALITY:
            name = "cardinality";
            break;
         case SCIP_CONSTYPE_INVKNAPSACK:
            name = "invknapsack";
            break;
         case SCIP_CONSTYPE_EQKNAPSACK:
            name = "eqknapsack";
            break;
         case SCIP_CONSTYPE_BINPACKING:
            name = "binpacking";
            break;
         case SCIP_CONSTYPE_KNAPSACK:
            name = "knapsack";
            break;
         case SCIP_CONSTYPE_INTKNAPSACK:
            name = "intknapsack";
            break;
         case SCIP_CONSTYPE_MIXEDBINARY:
            name = "mixed binary";
            break;
         case SCIP_CONSTYPE_GENERAL:
            name = "general";
            break;
         default:
            name = "unknown";
            break;
      }


#ifdef WRITE_ORIG_CONSTYPES
         myfile << " " <<  nfoundconstypesrangeddoublecount[c] << ",";
#endif

      classifier->setClassName( c, name.c_str() );
      text << "This class contains all constraints that are of (miplib) constype \"" << name << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }

#ifdef WRITE_ORIG_CONSTYPES
      myfile << std::endl;
      myfile.close();
#endif



   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classforcons[i] );
   }



   classifier->removeEmptyClasses();
   SCIPverbMessage(origprob, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier \"%s\" yields a classification with %d different constraint classes \n", classifier->getName(), classifier->getNClasses() );

   detprobdata->addConsPartition(classifier);
   return SCIP_OKAY;
}

/*
 * classifier specific interface methods
 */

SCIP_RETCODE GCGincludeConsClassifierMiplibConstypes(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   GCG_CLASSIFIERDATA* classifierdata = NULL;

   SCIP_CALL(
      GCGincludeConsClassifier(gcg, CLSCONS_NAME, CLSCONS_DESC, CLSCONS_PRIORITY, CLSCONS_ENABLED, classifierdata,
                               classifierFree, classifierClassify));

   return SCIP_OKAY;
}
