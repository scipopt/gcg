/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_zerohalf.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  {0,1/2}-cuts separator
 * @author Leona Gottwald
 * @author Manuel Kutschka
 * @author Kati Wolter
 *
 * {0,1/2}-Chvátal-Gomory cuts separator. It solves the following separation problem:
 * Consider an integer program
 * \f[
 * \min \{ c^T x : Ax \leq b, x \geq 0, x \mbox{ integer} \}
 * \f]
 * and a fractional solution \f$x^*\f$ of its LP relaxation.  Find a weightvector \f$u\f$ whose entries \f$u_i\f$ are either 0 or
 * \f$\frac{1}{2}\f$ such that the following inequality is valid for all integral solutions and violated by \f$x^*\f$:
 * \f[
 * \lfloor(u^T A) x \rfloor \leq \lfloor u^T b\rfloor
 * \f]
 *
 * References:
 * - Alberto Caprara, Matteo Fischetti. {0,1/2}-Chvatal-Gomory cuts. Math. Programming, Volume 74, p221--235, 1996.
 * - Arie M. C. A. Koster, Adrian Zymolka and Manuel Kutschka. \n
 *   Algorithms to separate {0,1/2}-Chvatal-Gomory cuts.
 *   Algorithms - ESA 2007: 15th Annual European Symposium, Eilat, Israel, October 8-10, 2007, \n
 *   Proceedings. Lecture Notes in Computer Science, Volume 4698, p. 693--704, 2007.
 * - Arie M. C. A. Koster, Adrian Zymolka and Manuel Kutschka. \n
 *   Algorithms to separate {0,1/2}-Chvatal-Gomory cuts (Extended Version). \n
 *   ZIB Report 07-10, Zuse Institute Berlin, 2007. http://www.zib.de/Publications/Reports/ZR-07-10.pdf
 * - Manuel Kutschka. Algorithmen zur Separierung von {0,1/2}-Schnitten. Diplomarbeit. Technische Universitaet Berlin, 2007.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "string.h"
#include "scip/sepa_zerohalf.h"
#include "scip/cutsel_hybrid.h"
#include "zerohalf_selector.h"

#include "gcg.h"

#define SEPA_NAME              "zerohalf"
#define SEPA_DESC              "{0,1/2}-cuts separator"
#define SEPA_PRIORITY             -6000
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE
#define SEPA_DELAY                FALSE

#define DEFAULT_MAXROUNDS             5 /**< maximal number of zerohalf separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        20 /**< maximal number of zerohalf separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          20 /**< maximal number of zerohalf cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     100 /**< maximal number of zerohalf cuts separated per separation round in root node */
#define DEFAULT_MAXCUTCANDS        2000 /**< maximal number of zerohalf cuts considered per separation round */
#define DEFAULT_MAXSLACK            0.0 /**< maximal slack of rows to be used in aggregation */
#define DEFAULT_MAXSLACKROOT        0.0 /**< maximal slack of rows to be used in aggregation in the root node */
#define DEFAULT_GOODSCORE           1.0 /**< threshold for score of cut relative to best score to be considered good,
                                         *   so that less strict filtering is applied */
#define DEFAULT_BADSCORE            0.5 /**< threshold for score of cut relative to best score to be discarded */
#define DEFAULT_MINVIOL             0.1 /**< minimal violation to generate zerohalfcut for */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAXROWDENSITY      0.05 /**< maximal density of row to be used in aggregation */
#define DEFAULT_DENSITYOFFSET       100 /**< additional number of variables allowed in row on top of density */
#define DEFAULT_INITSEED         0x5EED /**< default initial seed used for random tie-breaking in cut selection */
#define DEFAULT_OBJPARALWEIGHT      0.0 /**< weight of objective parallelism in cut score calculation */
#define DEFAULT_EFFICACYWEIGHT      1.0 /**< weight of efficacy in cut score calculation */
#define DEFAULT_DIRCUTOFFDISTWEIGHT 0.0 /**< weight of directed cutoff distance in cut score calculation */
#define DEFAULT_GOODMAXPARALL       0.1 /**< maximum parallelism for good cuts */
#define DEFAULT_MAXPARALL           0.1 /**< maximum parallelism for non-good cuts */

/* SCIPcalcRowIntegralScalar parameters */
#define MAXDNOM                  1000LL
#define MAXSCALE                 1000.0

/* other defines */
#define MAXREDUCTIONROUNDS          100 /**< maximum number of rounds to perform reductions on the mod 2 system */
#define BOUNDSWITCH                 0.5 /**< threshold for bound switching */
#define MAXAGGRLEN(nvars)           ((int)(0.1*(nvars)+1000))

typedef struct Mod2Col MOD2_COL;
typedef struct Mod2Row MOD2_ROW;
typedef struct Mod2Matrix MOD2_MATRIX;
typedef struct TransIntRow TRANSINTROW;


/** enum for different types of row indices in ROWINDEX structure */

#define ROWIND_TYPE unsigned int
#define ORIG_RHS    0u
#define ORIG_LHS    1u
#define TRANSROW    2u

/* macro to get a unique index from the rowindex */
#define UNIQUE_INDEX(rowind) (3*(rowind).index + (rowind).type)



/** structure containing a transformed integral row obtained by relaxing an lp row */
struct TransIntRow
{
   SCIP_Real             slack;              /**< slack of row after transformation */
   SCIP_Real             rhs;                /**< right hand side value of integral row after transformation */
   SCIP_Real*            vals;               /**< values of row */
   int*                  varinds;            /**< problem variable indices of row */
   int                   size;               /**< alloc size of row */
   int                   len;                /**< length of row */
   int                   rank;               /**< rank of row */
   SCIP_Bool             local;              /**< is row local? */
};

/** structure representing a row in the mod 2 system */
struct Mod2Row
{
   ROWINDEX*             rowinds;            /**< index set of rows associated with the mod 2 row */
   MOD2_COL**            nonzcols;           /**< sorted array of non-zero mod 2 columns in this mod 2 row */
   SCIP_Real             slack;              /**< slack of mod 2 row */
   SCIP_Real             maxsolval;          /**< maximum solution value of columns in mod 2 row */
   int                   index;              /**< unique index of mod 2 row */
   int                   pos;                /**< position of mod 2 row in mod 2 matrix rows array */
   int                   rhs;                /**< rhs of row */
   int                   nrowinds;           /**< number of elements in rowinds */
   int                   rowindssize;        /**< size of rowinds array */
   int                   nnonzcols;          /**< number of columns in nonzcols */
   int                   nonzcolssize;       /**< size of nonzcols array */
};

/** structure representing a column in the mod 2 system */
struct Mod2Col
{
   SCIP_HASHSET*         nonzrows;           /**< the set of rows that contain this column */
   SCIP_Real             solval;             /**< solution value of the column */
   int                   pos;                /**< position of column in matrix */
   int                   index;              /**< index of SCIP column associated to this column */
};

/** matrix representing the modulo 2 system */
struct Mod2Matrix
{
   MOD2_COL**            cols;               /**< columns of the matrix */
   MOD2_ROW**            rows;               /**< rows of the matrix */
   TRANSINTROW*          transintrows;       /**< transformed integral rows obtained from non-integral lp rows */
   int                   ntransintrows;      /**< number of transformed integral rows obtained from non-integral lp rows */
   int                   nzeroslackrows;     /**< number of rows with zero slack */
   int                   nrows;              /**< number of rows of the matrix; number of elements in rows */
   int                   ncols;              /**< number of cols of the matrix; number of elements in cols */
   int                   rowssize;           /**< length of rows array */
   int                   colssize;           /**< length of cols array */
};


/** data of separator */

//typedef struct GCG_ZeroHalfData GCG_ZEROHALFDATA;



/** create an instance of cut indices */
SCIP_RETCODE GCGcreateCutIndicesFromRowindex(
   SCIP*                scip,          /**< SCIP data structure (master problem) */
   GCG_CUTINDICES**     cutindices,    /**< pointer to store cutindices */
   int                  nindices,      /**< number of indices */
   ROWINDEX*            rowindex       /**< indices */
)
{
   int i;

   assert(scip != NULL);
   assert(nindices > 0);
   assert(cutindices != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, cutindices) );
   (*cutindices)->nindices = nindices;
   (*cutindices)->indices = NULL;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*cutindices)->indices), nindices) );
   for( i = 0; i < nindices; i++ )
   {
      (*cutindices)->indices[i] = rowindex[i].index;
   }

   return SCIP_OKAY;
}

/** create an instance of cut indices */
SCIP_RETCODE GCGcreateCutIndicesFromArray(
   SCIP*                scip,          /**< SCIP data structure (master problem) */
   GCG_CUTINDICES**     cutindices,    /**< pointer to store cutindices */
   int                  nindices,      /**< number of indices */
   int*                 indices        /**< indices */
)
{
   assert(scip != NULL);
   assert(nindices > 0);
   assert(cutindices != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, cutindices) );
   (*cutindices)->nindices = nindices;
   (*cutindices)->indices = indices;

   return SCIP_OKAY;
}

/** free an instance of cut indices */
SCIP_RETCODE GCGfreeCutIndices(
   SCIP*             scip,               /**< SCIP data structure (master problem) */
   GCG_CUTINDICES**  cutindices          /**< pointer to cutindices */
)
{
   assert(scip != NULL);
   assert(cutindices != NULL);
   assert(*cutindices != NULL);

   if( (*cutindices)->nindices > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &((*cutindices)->indices), (*cutindices)->nindices);
      (*cutindices)->nindices = 0;
   }

   assert((*cutindices)->indices == NULL);
   assert((*cutindices)->nindices == 0);

   SCIPfreeBlockMemory(scip, cutindices);
   *cutindices = NULL;

   return SCIP_OKAY;
}


#define COLINFO_GET_MOD2COL(x) ((MOD2_COL*)  (((uintptr_t)(x)) & ~((uintptr_t)1)))
#define COLINFO_GET_RHSOFFSET(x) ((int)  (((uintptr_t)(x)) & ((uintptr_t)1)))
#define COLINFO_CREATE(mod2col, rhsoffset)  ((void*) (((uintptr_t)(mod2col)) | ((uintptr_t)(rhsoffset))))


#ifndef NDEBUG
static
void checkRow(MOD2_ROW* row)
{
   int i;
   SCIP_Real maxsolval = 0.0;

   for( i = 0; i < row->nnonzcols; ++i )
   {
      assert(row->nonzcols[i]->solval > 0.0);
      maxsolval = MAX(maxsolval, row->nonzcols[i]->solval);

      if( i + 1 < row->nnonzcols )
         assert(row->nonzcols[i]->index < row->nonzcols[i+1]->index);
   }

   assert(row->maxsolval == maxsolval); /*lint !e777*/
}
#else
#define checkRow(x)
#endif

/** compare to mod 2 columns by there index */
static
SCIP_DECL_SORTPTRCOMP(compareColIndex)
{
   MOD2_COL* col1;
   MOD2_COL* col2;

   col1 = (MOD2_COL*) elem1;
   col2 = (MOD2_COL*) elem2;

   if( col1->index < col2->index )
      return -1;
   if( col2->index < col1->index )
      return 1;

   return 0;
}

/** comparison function for slack of mod 2 rows */
static
SCIP_DECL_SORTPTRCOMP(compareRowSlack)
{
   MOD2_ROW* row1;
   MOD2_ROW* row2;
   SCIP_Bool slack1iszero;
   SCIP_Bool slack2iszero;

   row1 = (MOD2_ROW*) elem1;
   row2 = (MOD2_ROW*) elem2;

   slack1iszero = EPSZ(row1->slack, SCIP_DEFAULT_EPSILON);
   slack2iszero = EPSZ(row2->slack, SCIP_DEFAULT_EPSILON);

   /* zero slack comes first */
   if( slack1iszero && !slack2iszero )
      return -1;
   if( slack2iszero && !slack1iszero )
      return 1;
   if( !slack1iszero && !slack2iszero )
      return 0;

   /* prefer rows that contain columns with large solution value */
   if( row1->maxsolval > row2->maxsolval )
      return -1;
   if( row2->maxsolval > row1->maxsolval )
      return 1;

   /* rows with less non-zeros come first rows */
   if( row1->nnonzcols < row2->nnonzcols )
      return -1;
   if( row2->nnonzcols < row1->nnonzcols )
      return 1;

   return 0;
}

/** take integral real value modulo 2 */
static
int mod2(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_Real             val                 /**< value to take mod 2 */
)
{
   assert(SCIPisFeasIntegral(scip, val));
   val *= 0.5;
   return (REALABS(SCIPround(scip, val) - val) > 0.1) ? 1 : 0;
}

/** returns the integral value for the given scaling parameters, see SCIPcalcIntegralScalar() */
static
void getIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real*            sval,               /**< pointer to store the scaled value */
   SCIP_Real*            intval              /**< pointer to store the scaled integral value */
)
{
   SCIP_Real upviol;
   SCIP_Real downviol;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   *sval = val * scalar;
   downval = floor(*sval);
   upval = ceil(*sval);

   downviol = SCIPrelDiff(*sval, downval) - maxdelta;
   upviol = mindelta - SCIPrelDiff(*sval, upval);

   if( downviol < upviol )
      *intval = downval;
   else
      *intval = upval;
}

/** Tries to transform a non-integral row into an integral row that can be used in zerohalf separation */
static
SCIP_RETCODE transformNonIntegralRow(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_Real             maxslack,           /**< maximum slack allowed for transformed row */
   int                   sign,               /**< +1 or -1 scale to select the side of the row */
   SCIP_Bool             local,              /**< is the row only valid locally? */
   int                   rank,               /**< rank of row */
   int                   rowlen,             /**< length of row */
   SCIP_Real*            rowvals,            /**< coefficients of columns in row */
   SCIP_COL**            rowcols,            /**< columns of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   int*                  intvarpos,          /**< clean buffer array of size SCIPgetNVars that will be clean when the function returns */
   TRANSINTROW*          introw,             /**< pointer to return transformed row */
   SCIP_Bool*            success             /**< pointer to return whether the transformation succeeded */
)
{
   int i;
   int transrowlen;
   SCIP_Real transrowrhs;
   int* transrowvars;
   SCIP_Real* transrowvals;

   assert(scip != NULL);
   assert(sign == +1 || sign == -1);
   assert(rowvals != NULL || rowlen == 0);
   assert(rowcols != NULL || rowlen == 0);
   assert(intvarpos != NULL);
   assert(introw != NULL);
   assert(success != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &transrowvars, rowlen) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &transrowvals, rowlen) );
   transrowlen = 0;
   transrowrhs = rhs;

   /* first add all integral variables to the transformed row and remember their positions in the row */
   for( i = 0; i < rowlen; ++i )
   {
      int probindex;

      if( !SCIPcolIsIntegral(rowcols[i]) )  /*lint !e613*/
         continue;

      probindex = SCIPcolGetVarProbindex(rowcols[i]);
      transrowvars[transrowlen] = probindex;
      transrowvals[transrowlen] = sign * rowvals[i];
      intvarpos[probindex] = ++transrowlen;
   }

   /* now loop over the non-integral columns of the row and project them out using simple or variable bounds */
   *success = TRUE;

   for( i = 0; i < rowlen; ++i )
   {
      int closestvbdind;
      SCIP_Real closestbound;
      SCIP_VAR* vbdvar;
      SCIP_Real vbdcoef;
      SCIP_Real vbdconst;
      SCIP_VAR* colvar;
      SCIP_Real val;
      SCIP_Real closestvbd;
      SCIP_Bool localbound;

      if( SCIPcolIsIntegral(rowcols[i]) ) /*lint !e613*/
         continue;

      localbound = FALSE;

      colvar = SCIPcolGetVar(rowcols[i]); /*lint !e613*/

      val = sign * rowvals[i]; /*lint !e613*/

      /* if the value is positive we need to use a lower bound constraint */
      if( val > 0.0 )
      {
         /* retrieve simple variable bound */
         closestbound = SCIPvarGetLbGlobal(colvar);
         if( allowlocal && SCIPisSumGT(scip, SCIPvarGetLbLocal(colvar), closestbound) )
         {
            /* only use local bound if it is better thatn the global bound */
            closestbound = SCIPvarGetLbLocal(colvar);
            localbound = TRUE;
         }

         /* retrieve closest variable bound */
         SCIP_CALL( SCIPgetVarClosestVlb(scip, colvar, NULL, &closestvbd, &closestvbdind) );

         /* if a suitable variable bound exists which is at least as good as a local simple bound
          * or better than a global simple bound we use it
          */
         if( closestvbdind >= 0 && (SCIPisGT(scip, closestvbd, closestbound) || (localbound && SCIPisSumEQ(scip, closestvbd, closestbound))) )
         {
            vbdcoef = SCIPvarGetVlbCoefs(colvar)[closestvbdind];
            vbdvar = SCIPvarGetVlbVars(colvar)[closestvbdind];
            vbdconst = SCIPvarGetVlbConstants(colvar)[closestvbdind];
            closestbound = closestvbd;
         }
         else
         {
            closestvbdind = -1;
         }
      }
      else
      {
         /* retrieve simple variable bound */
         closestbound = SCIPvarGetUbGlobal(colvar);
         if( allowlocal && SCIPisSumLT(scip, SCIPvarGetUbLocal(colvar), closestbound) )
         {
            closestbound = SCIPvarGetUbLocal(colvar);
            localbound = TRUE;
         }

         /* retrieve closest variable bound */
         SCIP_CALL( SCIPgetVarClosestVub(scip, colvar, NULL, &closestvbd, &closestvbdind) );

         /* if a suitable variable bound exists which is at least as good as a local simple bound
          * or better than a global simple bound we use it
          */
         if( closestvbdind >= 0 && (SCIPisLT(scip, closestvbd, closestbound) || (localbound && SCIPisSumEQ(scip, closestvbd, closestbound))) )
         {
            vbdcoef = SCIPvarGetVubCoefs(colvar)[closestvbdind];
            vbdvar = SCIPvarGetVubVars(colvar)[closestvbdind];
            vbdconst = SCIPvarGetVubConstants(colvar)[closestvbdind];
            closestbound = closestvbd;
         }
         else
         {
            closestvbdind = -1;
         }
      }

      if( closestvbdind >= 0 )
      {
         SCIP_Real coef;
         int pos;

         coef = val * vbdcoef; /*lint !e644*/
         transrowrhs -= val * vbdconst; /*lint !e644*/

         pos = intvarpos[SCIPvarGetProbindex(vbdvar)] - 1; /*lint !e644*/
         if( pos >= 0 )
         {
            transrowvals[pos] += coef;
         }
         else
         {
            transrowvars[transrowlen] = SCIPvarGetProbindex(vbdvar);
            transrowvals[transrowlen] = coef;
            intvarpos[SCIPvarGetProbindex(vbdvar)] = ++transrowlen;
         }
      }
      else if( !SCIPisInfinity(scip, REALABS(closestbound)) )
      {
         local = local || localbound;
         transrowrhs -= val * closestbound;
      }
      else
      {
         *success = FALSE;
         break;
      }
   }

   for( i = 0; i < transrowlen;)
   {
      intvarpos[transrowvars[i]] = 0;
      if( SCIPisZero(scip, transrowvals[i]) )
      {
         --transrowlen;
         transrowvals[i] = transrowvals[transrowlen];
         transrowvars[i] = transrowvars[transrowlen];
      }
      else
         ++i;
   }

   if( transrowlen <= 1 )
      *success = FALSE;

   if( *success )
   {
      SCIP_Real mindelta;
      SCIP_Real maxdelta;
      SCIP_Real intscalar;
      int nchgcoefs;

      SCIP_VAR** vars = SCIPgetVars(scip);

      *success = ! SCIPcutsTightenCoefficients(scip, local, transrowvals, &transrowrhs, transrowvars, &transrowlen, &nchgcoefs);

      mindelta = -SCIPepsilon(scip);
      maxdelta = SCIPsumepsilon(scip);

      if( *success )
      {
         SCIP_CALL( SCIPcalcIntegralScalar(transrowvals, transrowlen, mindelta, maxdelta, MAXDNOM, MAXSCALE, &intscalar, success) );

         if( *success )
         {
            SCIP_Real floorrhs;
            SCIP_Real slack;

            transrowrhs *= intscalar; /*lint !e644*/

            /* slack is initialized to zero since the transrowrhs can still change due to bound usage in the loop below;
             * the floored right hand side is then added afterwards
             */
            slack = 0.0;
            for( i = 0; i < transrowlen; ++i )
            {
               SCIP_Real solval = SCIPgetSolVal(scip, sol, vars[transrowvars[i]]);
               SCIP_Real intval;
               SCIP_Real newval;

               getIntegralScalar(transrowvals[i], intscalar, mindelta, maxdelta, &newval, &intval);

               if( !SCIPisEQ(scip, intval, newval) )
               {
                  if( intval < newval )
                  {
                     SCIP_Real lb = local ? SCIPvarGetLbLocal(vars[transrowvars[i]]) : SCIPvarGetLbGlobal(vars[transrowvars[i]]);

                     if( SCIPisInfinity(scip, -lb) )
                     {
                        *success = FALSE;
                        break;
                     }

                     transrowrhs += (intval - newval) * lb;
                  }
                  else
                  {
                     SCIP_Real ub = local ? SCIPvarGetUbLocal(vars[transrowvars[i]]) : SCIPvarGetUbGlobal(vars[transrowvars[i]]);

                     if( SCIPisInfinity(scip, ub) )
                     {
                        *success = FALSE;
                        break;
                     }

                     transrowrhs += (intval - newval) * ub;
                  }
               }

               slack -= solval * intval;
               transrowvals[i] = intval;
            }

            if( *success )
            {
               floorrhs = SCIPfeasFloor(scip, transrowrhs);
               slack += floorrhs;

               if( slack <= maxslack )
               {
                  introw->rhs = floorrhs;
                  introw->slack = slack;
                  introw->vals = transrowvals;
                  introw->varinds = transrowvars;
                  introw->len = transrowlen;
                  introw->size = rowlen;
                  introw->local = local;
                  introw->rank = rank;

                  if( !SCIPisEQ(scip, floorrhs, transrowrhs) )
                     introw->rank += 1;
               }
               else
               {
                  *success = FALSE;
               }
            }
         }
      }
   }

   if( !(*success) )
   {
      SCIPfreeBlockMemoryArray(scip, &transrowvals, rowlen);
      SCIPfreeBlockMemoryArray(scip, &transrowvars, rowlen);
   }

   return SCIP_OKAY;
}


/** Tries to transform non-integral rows into an integral form by using simple and variable bounds */
static
SCIP_RETCODE mod2MatrixTransformContRows(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   GCG_ZEROHALFDATA *    zhdata,           /**< zerohalf separator data */
   MOD2_MATRIX*          mod2matrix,         /**< mod2 matrix structure */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_Real             maxslack            /**< maximum slack allowed for mod 2 rows */
)
{
   SCIP_ROW** rows;
   int nrows;
   int* intvarpos;
   int i;
   int maxnonzeros;
   SCIP_CONS** masterconss;
   SCIP* origscip;

   /* we only use the rows from the master constraints; Vladi used all */
   /* OG: SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );*/

   origscip = GCGgetOriginalprob(scip);
   masterconss = GCGgetMasterConss(origscip);
   nrows = GCGgetNMasterConss(origscip);

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, nrows) );

   for( i = 0; i < nrows; i++ )
   {
      SCIP_ROW* consrow;

      consrow = SCIPconsGetRow(scip, masterconss[i]);
      assert(consrow != NULL);
      rows[i] = consrow;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &mod2matrix->transintrows, 2*nrows) );
   mod2matrix->ntransintrows = 0;

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &intvarpos, SCIPgetNVars(scip)) );

   maxnonzeros = (int)(SCIPgetNLPCols(scip) * zhdata->maxrowdensity) + zhdata->densityoffset;

   for( i = 0; i < nrows; ++i )
   {
      int rowlen;
      SCIP_Real activity;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real lhsslack;
      SCIP_Real rhsslack;
      SCIP_Real* rowvals;
      SCIP_COL** rowcols;

      /* skip integral rows and rows not suitable for generating cuts */
      /* OG: if( SCIProwIsModifiable(rows[i]) || SCIProwIsIntegral(rows[i]) || (SCIProwIsLocal(rows[i]) && !allowlocal) || SCIProwGetNNonz(rows[i]) > maxnonzeros ) */
      if( SCIProwIsIntegral(rows[i]) || (SCIProwIsLocal(rows[i]) && !allowlocal) || SCIProwGetNNonz(rows[i]) > maxnonzeros )
         continue;

      lhs = SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]);
      rhs = SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]);
      activity = SCIPgetRowSolActivity(scip, rows[i], sol) - SCIProwGetConstant(rows[i]);

      /* compute lhsslack: activity - lhs */
      if( SCIPisInfinity(scip, -SCIProwGetLhs(rows[i])) )
         lhsslack = SCIPinfinity(scip);
      else
      {
         lhsslack = activity - lhs;
      }

      /* compute rhsslack: rhs - activity */
      if( SCIPisInfinity(scip, SCIProwGetRhs(rows[i])) )
         rhsslack = SCIPinfinity(scip);
      else
         rhsslack = rhs - activity;

      if( rhsslack > maxslack && lhsslack > maxslack )
         continue;

      rowlen = SCIProwGetNLPNonz(rows[i]);
      rowvals = SCIProwGetVals(rows[i]);
      rowcols = SCIProwGetCols(rows[i]);

      if( rhsslack <= maxslack )
      {
         SCIP_Bool success;
         TRANSINTROW* introw = &mod2matrix->transintrows[mod2matrix->ntransintrows];
         SCIP_CALL( transformNonIntegralRow(scip, sol, allowlocal, maxslack, 1, SCIProwIsLocal(rows[i]), SCIProwGetRank(rows[i]), \
               rowlen, rowvals, rowcols, rhs, intvarpos, introw, &success) );

         assert(success == 1 || success == 0);
         mod2matrix->ntransintrows += (int)success;
      }

      if( lhsslack <= maxslack )
      {
         SCIP_Bool success;
         TRANSINTROW* introw = &mod2matrix->transintrows[mod2matrix->ntransintrows];
         SCIP_CALL( transformNonIntegralRow(scip, sol, allowlocal, maxslack, -1, SCIProwIsLocal(rows[i]), SCIProwGetRank(rows[i]), \
               rowlen, rowvals, rowcols, -lhs, intvarpos, introw, &success) );

         assert(success == 1 || success == 0);
         mod2matrix->ntransintrows += (int)success;
      }
   }

   SCIPfreeCleanBufferArray(scip, &intvarpos);
   SCIPfreeBufferArray(scip, &rows);

   return SCIP_OKAY;
}


/** adds new column to the mod 2 matrix */
static
SCIP_RETCODE mod2MatrixAddCol(
   SCIP*                 scip,               /**< SCIP datastructure */
   MOD2_MATRIX*          mod2matrix,         /**< mod 2 matrix */
   SCIP_HASHMAP*         origvar2col,        /**< hash map for mapping of problem variables to mod 2 columns */
   SCIP_VAR*             origvar,            /**< problem variable to create mod 2 column for */
   SCIP_Real             solval,             /**< solution value of problem variable */
   int                   rhsoffset           /**< offset in right hand side due complementation (mod 2) */
)
{
   MOD2_COL* col;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, &col) );

   /* initialize fields */
   col->pos = mod2matrix->ncols++;
   col->index = SCIPvarGetProbindex(origvar);
   col->solval = solval;
   SCIP_CALL( SCIPhashsetCreate(&col->nonzrows, SCIPblkmem(scip), 1) );

   /* add column to mod 2 matrix */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &mod2matrix->cols, &mod2matrix->colssize, mod2matrix->ncols) );
   mod2matrix->cols[col->pos] = col;

   /* create mapping of problem variable to mod 2 column with its right hand side offset */
   assert(rhsoffset >= 0);
   SCIP_CALL( SCIPhashmapInsert(origvar2col, (void*) origvar, COLINFO_CREATE(col, rhsoffset)) ); /*lint !e571*/

   return SCIP_OKAY;
}

/** links row to mod 2 column */
static
SCIP_RETCODE mod2colLinkRow(
   BMS_BLKMEM*           blkmem,             /**< block memory shell */
   MOD2_COL*             col,                /**< mod 2 column */
   MOD2_ROW*             row                 /**< mod 2 row */
)
{
   SCIP_CALL( SCIPhashsetInsert(col->nonzrows, blkmem, (void*)row) );

   assert(SCIPhashsetExists(col->nonzrows, (void*)row));

   row->maxsolval = MAX(col->solval, row->maxsolval);

   return SCIP_OKAY;
}

/** unlinks row from mod 2 column */
static
SCIP_RETCODE mod2colUnlinkRow(
   MOD2_COL*             col,                /**< mod 2 column */
   MOD2_ROW*             row                 /**< mod 2 row */
)
{
   SCIP_CALL( SCIPhashsetRemove(col->nonzrows, (void*)row) );

   assert(!SCIPhashsetExists(col->nonzrows, (void*)row));
#ifndef NDEBUG
   {
      int nslots = SCIPhashsetGetNSlots(col->nonzrows);
      MOD2_ROW** rows = (MOD2_ROW**) SCIPhashsetGetSlots(col->nonzrows);
      int i;

      for( i = 0; i < nslots; ++i )
      {
         assert(rows[i] != row);
      }
   }
#endif

   return SCIP_OKAY;
}

/** unlinks row from mod 2 column */
static
void mod2rowUnlinkCol(
   MOD2_ROW*             row                 /**< mod 2 row */,
   MOD2_COL*             col                 /**< mod 2 column */
)
{
   int i;

   assert(row->nnonzcols == 0 || row->nonzcols != NULL);

   SCIP_UNUSED( SCIPsortedvecFindPtr((void**) row->nonzcols, compareColIndex, col, row->nnonzcols, &i) );
   assert(row->nonzcols[i] == col);

   --row->nnonzcols;
   BMSmoveMemoryArray(row->nonzcols + i, row->nonzcols + i + 1, row->nnonzcols - i); /*lint !e866*/

   if( col->solval >= row->maxsolval )
   {
      row->maxsolval = 0.0;
      for( i = 0; i < row->nnonzcols; ++i )
      {
         row->maxsolval = MAX(row->nonzcols[i]->solval, row->maxsolval);
      }
   }
}

/** adds a SCIP_ROW to the mod 2 matrix */
static
SCIP_RETCODE mod2MatrixAddOrigRow(
   SCIP*                 scip,               /**< scip data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory shell */
   MOD2_MATRIX*          mod2matrix,         /**< modulo 2 matrix */
   SCIP_HASHMAP*         origcol2col,        /**< hashmap to retrieve the mod 2 column from a SCIP_COL */
   SCIP_ROW*             origrow,            /**< original SCIP row */
   SCIP_Real             slack,              /**< slack of row */
   ROWIND_TYPE           side,               /**< side of row that is used for mod 2 row, must be ORIG_RHS or ORIG_LHS */
   int                   rhsmod2             /**< modulo 2 value of the row's right hand side */
)
{
   SCIP_Real* rowvals;
   SCIP_COL** rowcols;
   int rowlen;
   int i;
   MOD2_ROW* row;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &row) );

   row->index = mod2matrix->nrows++;
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &mod2matrix->rows, &mod2matrix->rowssize, mod2matrix->nrows) );
   mod2matrix->rows[row->index] = row;

   row->slack = MAX(0.0, slack);
   row->maxsolval = 0.0;
   row->rhs = rhsmod2;
   row->nrowinds = 1;
   row->rowinds = NULL;
   row->rowindssize = 0;

   if( SCIPisZero(scip, row->slack) )
      ++mod2matrix->nzeroslackrows;

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &row->rowinds, &row->rowindssize, row->nrowinds) );
   row->rowinds[0].type = side;
   row->rowinds[0].index = (unsigned int)SCIProwGetLPPos(origrow);

   row->nnonzcols = 0;
   row->nonzcolssize = 0;
   row->nonzcols = NULL;

   rowlen = SCIProwGetNNonz(origrow);
   rowvals = SCIProwGetVals(origrow);
   rowcols = SCIProwGetCols(origrow);

   for( i = 0; i < rowlen; ++i )
   {
      if( mod2(scip, rowvals[i]) == 1 )
      {
         void* colinfo;
         MOD2_COL* col;
         int rhsoffset;

         colinfo = SCIPhashmapGetImage(origcol2col, (void*)SCIPcolGetVar(rowcols[i]));

         /* extract the righthand side offset from the colinfo and update the righthand side */
         rhsoffset = COLINFO_GET_RHSOFFSET(colinfo);
         row->rhs = (row->rhs + rhsoffset) % 2;

         /* extract the column pointer from the colinfo */
         col = COLINFO_GET_MOD2COL(colinfo);

         if( col != NULL )
         {
            int k;

            k = row->nnonzcols++;

            SCIP_CALL( SCIPensureBlockMemoryArray(scip, &row->nonzcols, &row->nonzcolssize, row->nnonzcols) );
            row->nonzcols[k] = col;

            SCIP_CALL( mod2colLinkRow(blkmem, col, row) );
         }
      }
   }

   SCIPsortPtr((void**)row->nonzcols, compareColIndex, row->nnonzcols);

   checkRow(row);

   return SCIP_OKAY;
}

/** adds a transformed integral row to the mod 2 matrix */
static
SCIP_RETCODE mod2MatrixAddTransRow(
   SCIP*                 scip,               /**< scip data structure */
   MOD2_MATRIX*          mod2matrix,         /**< modulo 2 matrix */
   SCIP_HASHMAP*         origcol2col,        /**< hashmap to retrieve the mod 2 column from a SCIP_COL */
   int                   transrowind         /**< index to transformed int row */
)
{
   int i;
   SCIP_VAR** vars;
   BMS_BLKMEM* blkmem;
   MOD2_ROW* row;
   TRANSINTROW* introw;

   SCIP_CALL( SCIPallocBlockMemory(scip, &row) );

   vars = SCIPgetVars(scip);
   introw = &mod2matrix->transintrows[transrowind];

   blkmem = SCIPblkmem(scip);
   row->index = mod2matrix->nrows++;
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &mod2matrix->rows, &mod2matrix->rowssize, mod2matrix->nrows) );
   mod2matrix->rows[row->index] = row;

   row->slack = MAX(0.0, introw->slack);
   row->rhs = mod2(scip, introw->rhs);
   row->nrowinds = 1;
   row->rowinds = NULL;
   row->rowindssize = 0;
   row->maxsolval = 0.0;

   if( SCIPisZero(scip, row->slack) )
      ++mod2matrix->nzeroslackrows;

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &row->rowinds, &row->rowindssize, row->nrowinds) );
   row->rowinds[0].type = TRANSROW;
   row->rowinds[0].index = (unsigned int)transrowind;

   row->nnonzcols = 0;
   row->nonzcolssize = 0;
   row->nonzcols = NULL;

   for( i = 0; i < introw->len; ++i )
   {
      if( mod2(scip, introw->vals[i]) == 1 )
      {
         void* colinfo;
         MOD2_COL* col;
         int rhsoffset;

         colinfo = SCIPhashmapGetImage(origcol2col, (void*)vars[introw->varinds[i]]);

         /* extract the righthand side offset from the colinfo and update the righthand side */
         rhsoffset = COLINFO_GET_RHSOFFSET(colinfo);
         row->rhs = (row->rhs + rhsoffset) % 2;

         /* extract the column pointer from the colinfo */
         col = COLINFO_GET_MOD2COL(colinfo);

         if( col != NULL )
         {
            int k;

            k = row->nnonzcols++;

            SCIP_CALL( SCIPensureBlockMemoryArray(scip, &row->nonzcols, &row->nonzcolssize, row->nnonzcols) );
            row->nonzcols[k] = col;

            SCIP_CALL( mod2colLinkRow(blkmem, col, row) );
         }
      }
   }

   SCIPsortPtr((void**)row->nonzcols, compareColIndex, row->nnonzcols);

   checkRow(row);

   return SCIP_OKAY;
}

/** free all resources held by the mod 2 matrix */
static
void destroyMod2Matrix(
   SCIP*                 scip,               /**< scip data structure */
   MOD2_MATRIX*          mod2matrix          /**< pointer to mod2 matrix structure */
)
{
   int i;

   for( i = 0; i < mod2matrix->ncols; ++i )
   {
      SCIPhashsetFree(&mod2matrix->cols[i]->nonzrows, SCIPblkmem(scip));
      SCIPfreeBlockMemory(scip, &mod2matrix->cols[i]); /*lint !e866*/
   }

   for( i = 0; i < mod2matrix->nrows; ++i )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &mod2matrix->rows[i]->nonzcols, mod2matrix->rows[i]->nonzcolssize);
      SCIPfreeBlockMemoryArrayNull(scip, &mod2matrix->rows[i]->rowinds, mod2matrix->rows[i]->rowindssize);
      SCIPfreeBlockMemory(scip, &mod2matrix->rows[i]); /*lint !e866*/
   }

   for( i = 0; i < mod2matrix->ntransintrows; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &mod2matrix->transintrows[i].vals, mod2matrix->transintrows[i].size);
      SCIPfreeBlockMemoryArray(scip, &mod2matrix->transintrows[i].varinds, mod2matrix->transintrows[i].size);
   }
   SCIP* origscip = GCGgetOriginalprob(scip);

   SCIPfreeBlockMemoryArray(scip, &mod2matrix->transintrows, 2* GCGgetNMasterConss(origscip)); /*lint !e647*/

   SCIPfreeBlockMemoryArrayNull(scip, &mod2matrix->rows, mod2matrix->rowssize);
   SCIPfreeBlockMemoryArrayNull(scip, &mod2matrix->cols, mod2matrix->colssize);
}

/** build the modulo 2 matrix from all integral rows in the LP, and non-integral rows
 *  if the transformation to an integral row succeeds
 */
static
SCIP_RETCODE buildMod2Matrix(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   GCG_ZEROHALFDATA *    zhdata,           /**< zerohalf separator data */
   BMS_BLKMEM*           blkmem,             /**< block memory shell */
   MOD2_MATRIX*          mod2matrix,         /**< mod 2 matrix */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_Real             maxslack            /**< maximum slack allowed for mod 2 rows */
)
{
   SCIP_VAR** vars;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   SCIP_HASHMAP* origcol2col;
   int ncols;
   int nrows;
   int nintvars;
   int maxnonzeros;
   int i;
   SCIP* origscip;
   SCIP_CONS** mastercons;

   origscip = GCGgetOriginalprob(scip);
   mastercons = GCGgetMasterConss(origscip);
   nrows = GCGgetNMasterConss(origscip);

   SCIPallocBufferArray(scip, &rows, nrows);
   for( i = 0; i < nrows; i++ )
   {
      SCIP_ROW* consrow;

      consrow = SCIPconsGetRow(scip, mastercons[i]);
      assert(consrow != NULL);
      rows[i] = consrow;
   }

   /* OG: SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) ); */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   nintvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   vars = SCIPgetVars(scip);

   /* initialize fields */
   mod2matrix->cols = NULL;
   mod2matrix->colssize = 0;
   mod2matrix->ncols = 0;
   mod2matrix->rows = NULL;
   mod2matrix->rowssize = 0;
   mod2matrix->nrows = 0;
   mod2matrix->nzeroslackrows = 0;

   SCIP_CALL( SCIPhashmapCreate(&origcol2col, SCIPblkmem(scip), 1) );

   /* add all integral vars if they are not at their bound */
   for( i = 0; i < nintvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real lbsol;
      SCIP_Real ubsol;
      SCIP_Real primsol;
      SCIP_Bool useub;

      primsol = SCIPgetSolVal(scip, sol, vars[i]);

      lb = allowlocal ? SCIPvarGetLbLocal(vars[i]) : SCIPvarGetLbGlobal(vars[i]);
      lbsol = MAX(0.0, primsol - lb);
      if( SCIPisZero(scip, lbsol) )
      {
         SCIP_CALL( SCIPhashmapInsert(origcol2col, (void*) vars[i], COLINFO_CREATE(NULL, mod2(scip, lb))) ); /*lint !e571*/
         continue;
      }

      ub = allowlocal ? SCIPvarGetUbLocal(vars[i]) : SCIPvarGetUbGlobal(vars[i]);
      ubsol = MAX(0.0, ub - primsol);
      if( SCIPisZero(scip, ubsol) )
      {
         SCIP_CALL( SCIPhashmapInsert(origcol2col, (void*) vars[i], COLINFO_CREATE(NULL, mod2(scip, ub))) ); /*lint !e571*/
         continue;
      }

      if( SCIPisInfinity(scip, ub) ) /* if there is no ub, use lb */
         useub = FALSE;
      else if( SCIPisInfinity(scip, -lb) ) /* if there is no lb, use ub */
         useub = TRUE;
      else if( SCIPisLT(scip, primsol, (1.0 - BOUNDSWITCH) * lb + BOUNDSWITCH * ub) )
         useub = FALSE;
      else
         useub = TRUE;

      if( useub )
      {
         assert(ubsol > 0.0);

         /* coverity[var_deref_model] */
         SCIP_CALL( mod2MatrixAddCol(scip, mod2matrix, origcol2col, vars[i], ubsol, mod2(scip, ub)) );
      }
      else
      {
         assert(lbsol > 0.0);

         /* coverity[var_deref_model] */
         SCIP_CALL( mod2MatrixAddCol(scip, mod2matrix, origcol2col, vars[i], lbsol, mod2(scip, lb)) );
      }
   }

   maxnonzeros = (int)(SCIPgetNLPCols(scip) * zhdata->maxrowdensity) + zhdata->densityoffset;

   /* add all integral rows using the created columns */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real activity;
      SCIP_Real lhsslack;
      SCIP_Real rhsslack;
      int lhsmod2;
      int rhsmod2;

      /* skip non-integral rows and rows not suitable for generating cuts */
      /* OG: if( SCIProwIsModifiable(rows[i]) || !SCIProwIsIntegral(rows[i]) || (SCIProwIsLocal(rows[i]) && !allowlocal) || SCIProwGetNNonz(rows[i]) > maxnonzeros )*/
      if( !SCIProwIsIntegral(rows[i]) || (SCIProwIsLocal(rows[i]) && !allowlocal) || SCIProwGetNNonz(rows[i]) > maxnonzeros )
         continue;

      lhsmod2 = 0;
      rhsmod2 = 0;
      activity = SCIPgetRowSolActivity(scip, rows[i], sol) - SCIProwGetConstant(rows[i]);

      /* since row is integral we can ceil/floor the lhs/rhs after subtracting the constant */
      lhs = SCIPfeasCeil(scip, SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]));
      rhs = SCIPfeasFloor(scip, SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]));

      /* compute lhsslack: activity - lhs */
      if( SCIPisInfinity(scip, -SCIProwGetLhs(rows[i])) )
         lhsslack = SCIPinfinity(scip);
      else
      {
         lhsslack = activity - lhs;
         lhsmod2 = mod2(scip, lhs);
      }

      /* compute rhsslack: rhs - activity */
      if( SCIPisInfinity(scip, SCIProwGetRhs(rows[i])) )
         rhsslack = SCIPinfinity(scip);
      else
      {
         rhsslack = rhs - activity;
         rhsmod2 = mod2(scip, rhs);
      }

      if( rhsslack <= maxslack && lhsslack <= maxslack )
      {
         if( lhsmod2 == rhsmod2 )
         {
            /* maxslack < 1 implies rhs - lhs = rhsslack + lhsslack < 2. Therefore lhs = rhs (mod2) can only hold if they
             * are equal
             */
            assert(SCIPisEQ(scip, lhs, rhs));

            /* use rhs */
            /* coverity[var_deref_model] */
            SCIP_CALL( mod2MatrixAddOrigRow(scip, blkmem, mod2matrix, origcol2col, rows[i], rhsslack, ORIG_RHS, rhsmod2) );
         }
         else
         {
            /* use both */
            /* coverity[var_deref_model] */
            SCIP_CALL( mod2MatrixAddOrigRow(scip, blkmem, mod2matrix, origcol2col, rows[i], lhsslack, ORIG_LHS, lhsmod2) );
            SCIP_CALL( mod2MatrixAddOrigRow(scip, blkmem, mod2matrix, origcol2col, rows[i], rhsslack, ORIG_RHS, rhsmod2) );
         }
      }
      else if( rhsslack <= maxslack )
      {
         /* use rhs */
         /* coverity[var_deref_model] */
         SCIP_CALL( mod2MatrixAddOrigRow(scip, blkmem, mod2matrix, origcol2col, rows[i], rhsslack, ORIG_RHS, rhsmod2) );
      }
      else if( lhsslack <= maxslack )
      {
         /* use lhs */
         /* coverity[var_deref_model] */
         SCIP_CALL( mod2MatrixAddOrigRow(scip, blkmem, mod2matrix, origcol2col, rows[i], lhsslack, ORIG_LHS, lhsmod2) );
      }
   }

   /* transform non-integral rows */
   SCIP_CALL( mod2MatrixTransformContRows(scip, sol, zhdata, mod2matrix, allowlocal, maxslack) );

   /* add all transformed integral rows using the created columns */
   for( i = 0; i < mod2matrix->ntransintrows; ++i )
   {
      SCIP_CALL( mod2MatrixAddTransRow(scip, mod2matrix, origcol2col, i) );
   }

   SCIPhashmapFree(&origcol2col);
   SCIPfreeBufferArray(scip, &rows);

   return SCIP_OKAY;
}

/* compare two mod 2 columns for equality */
static
SCIP_DECL_HASHKEYEQ(columnsEqual)
{  /*lint --e{715}*/
   MOD2_COL* col1;
   MOD2_COL* col2;
   int nslotscol1;
   MOD2_ROW** col1rows;
   int i;

   col1 = (MOD2_COL*) key1;
   col2 = (MOD2_COL*) key2;

   if( SCIPhashsetGetNElements(col1->nonzrows) != SCIPhashsetGetNElements(col2->nonzrows) )
      return FALSE;

   nslotscol1 = SCIPhashsetGetNSlots(col1->nonzrows);
   col1rows = (MOD2_ROW**) SCIPhashsetGetSlots(col1->nonzrows);
   for( i = 0; i < nslotscol1; ++i )
   {
      if( col1rows[i] != NULL && !SCIPhashsetExists(col2->nonzrows, (void*)col1rows[i]) )
         return FALSE;
   }

   return TRUE;
}

/* compute a signature of the rows in a mod 2 matrix as hash value */
static
SCIP_DECL_HASHKEYVAL(columnGetSignature)
{  /*lint --e{715}*/
   MOD2_COL* col;
   MOD2_ROW** rows;
   uint64_t signature;
   int i;
   int nslots;

   col = (MOD2_COL*) key;

   nslots = SCIPhashsetGetNSlots(col->nonzrows);
   rows = (MOD2_ROW**) SCIPhashsetGetSlots(col->nonzrows);

   signature = 0;
   for( i = 0; i < nslots; ++i )
   {
      if( rows[i] != NULL )
         signature |= SCIPhashSignature64(rows[i]->index);
   }

   return signature;
}

/* compare two mod 2 rows for equality */
static
SCIP_DECL_HASHKEYEQ(rowsEqual)
{  /*lint --e{715}*/
   MOD2_ROW* row1;
   MOD2_ROW* row2;
   int i;

   row1 = (MOD2_ROW*) key1;
   row2 = (MOD2_ROW*) key2;

   assert(row1 != NULL);
   assert(row2 != NULL);
   assert(row1->nnonzcols == 0 || row1->nonzcols != NULL);
   assert(row2->nnonzcols == 0 || row2->nonzcols != NULL);

   if( row1->nnonzcols != row2->nnonzcols || row1->rhs != row2->rhs )
      return FALSE;

   for( i = 0; i < row1->nnonzcols; ++i )
   {
      if( row1->nonzcols[i] != row2->nonzcols[i] )
         return FALSE;
   }

   return TRUE;
}

/* compute a signature of a mod 2 row as hash value */
static
SCIP_DECL_HASHKEYVAL(rowGetSignature)
{  /*lint --e{715}*/
   MOD2_ROW* row;
   int i;
   uint64_t signature;

   row = (MOD2_ROW*) key;
   assert(row->nnonzcols == 0 || row->nonzcols != NULL);

   signature = row->rhs; /*lint !e732*/

   for( i = 0; i < row->nnonzcols; ++i )
      signature |= SCIPhashSignature64(row->nonzcols[i]->index);

   return signature;
}

/** removes a row from the mod 2 matrix */
static
SCIP_RETCODE mod2matrixRemoveRow(
   SCIP*                 scip,               /**< scip data structure */
   MOD2_MATRIX*          mod2matrix,         /**< the mod 2 matrix */
   MOD2_ROW*             row                 /**< mod 2 row */
)
{
   int i;
   int position = row->pos;

   checkRow(row);

   /* update counter for zero slack rows */
   if( SCIPisZero(scip, row->slack) )
      --mod2matrix->nzeroslackrows;

   /* remove the row from the array */
   --mod2matrix->nrows;
   mod2matrix->rows[position] = mod2matrix->rows[mod2matrix->nrows];
   mod2matrix->rows[position]->pos = position;

   /* unlink columns from row */
   for( i = 0; i < row->nnonzcols; ++i )
   {
      SCIP_CALL( mod2colUnlinkRow(row->nonzcols[i], row) );
   }

   /* free row */
   SCIPfreeBlockMemoryArrayNull(scip, &row->nonzcols, row->nonzcolssize);
   SCIPfreeBlockMemoryArray(scip, &row->rowinds, row->rowindssize);
   SCIPfreeBlockMemory(scip, &row);

   return SCIP_OKAY;
}

/** removes a column from the mod 2 matrix */
static
void mod2matrixRemoveCol(
   SCIP*                 scip,               /**< scip data structure */
   MOD2_MATRIX*          mod2matrix,         /**< the mod 2 matrix */
   MOD2_COL*             col                 /**< a column in the mod 2 matrix */
)
{
   int i;
   int nslots;
   MOD2_ROW** rows;
   int position;

   assert(col != NULL);

   /* cppcheck-suppress nullPointer */
   position = col->pos;

   /* remove column from arrays */
   --mod2matrix->ncols;
   mod2matrix->cols[position] = mod2matrix->cols[mod2matrix->ncols];
   mod2matrix->cols[position]->pos = position;

   /* cppcheck-suppress nullPointer */
   nslots = SCIPhashsetGetNSlots(col->nonzrows);
   /* cppcheck-suppress nullPointer */
   rows = (MOD2_ROW**) SCIPhashsetGetSlots(col->nonzrows);

   /* adjust rows of column */
   for( i = 0; i < nslots; ++i )
   {
      if( rows[i] != NULL )
         mod2rowUnlinkCol(rows[i], col);
   }

   /* free column */
   SCIPhashsetFree(&col->nonzrows, SCIPblkmem(scip));
   SCIPfreeBlockMemory(scip, &col);
}

/* remove columns that are (Prop3 iii) zero (Prop3 iv) identify indentical columns (Prop3 v) unit vector columns */
static
SCIP_RETCODE mod2matrixPreprocessColumns(
   SCIP*                 scip,               /**< scip data structure */
   MOD2_MATRIX*          mod2matrix,         /**< mod 2 matrix */
   GCG_ZEROHALFDATA*     zhdata
)
{
   int i;
   SCIP_HASHTABLE* columntable;

   SCIP_CALL( SCIPhashtableCreate(&columntable, SCIPblkmem(scip), mod2matrix->ncols,
                                  SCIPhashGetKeyStandard, columnsEqual, columnGetSignature, NULL) );

   for( i = 0; i < mod2matrix->ncols; )
   {
      MOD2_COL* col = mod2matrix->cols[i];
      int nnonzrows = SCIPhashsetGetNElements(col->nonzrows);
      if( nnonzrows == 0 )
      { /* Prop3 iii */
         mod2matrixRemoveCol(scip, mod2matrix, col);
      }
      else if( nnonzrows == 1 )
      { /* Prop3 v */
         MOD2_ROW* row;

         {
            int j = 0;
            MOD2_ROW** rows;
            rows = (MOD2_ROW**) SCIPhashsetGetSlots(col->nonzrows);
            while( rows[j] == NULL )
               ++j;

            row = rows[j];
         }

         checkRow(row);

         /* column is unit vector, so add its solution value to the rows slack and remove it */
         if( SCIPisZero(scip, row->slack) )
            --mod2matrix->nzeroslackrows;

         row->slack += col->solval;
         assert(!SCIPisZero(scip, row->slack));

         mod2matrixRemoveCol(scip, mod2matrix, col);
         ++zhdata->nreductions;

         checkRow(row);
      }
      else
      {
         MOD2_COL* identicalcol;
         identicalcol = (MOD2_COL*)SCIPhashtableRetrieve(columntable, col);
         if( identicalcol != NULL )
         {
            assert(identicalcol != col);

            /* column is identical to other column so add its solution value to the other one and then remove and free it */
            identicalcol->solval += col->solval;

            /* also adjust the solval of the removed column so that the maxsolval of each row is properly updated */
            col->solval = identicalcol->solval;

            mod2matrixRemoveCol(scip, mod2matrix, col);
         }
         else
         {
            SCIP_CALL( SCIPhashtableInsert(columntable, (void*)col) );
            ++i;
         }
      }
   }

   SCIPhashtableFree(&columntable);

   return SCIP_OKAY;
}

#define NONZERO(x)   (COPYSIGN(1e-100, (x)) + (x))

/** add original row to aggregation with weight 0.5 */
static
void addOrigRow(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            tmpcoefs,           /**< array to add coefficients to */
   SCIP_Real*            cutrhs,             /**< pointer to add right hand side */
   int*                  nonzeroinds,        /**< array of non-zeros in the aggregation */
   int*                  nnz,                /**< pointer to update number of non-zeros */
   int*                  cutrank,            /**< pointer to update cut rank */
   SCIP_Bool*            cutislocal,         /**< pointer to update local flag */
   SCIP_ROW*             row,                /**< row to add */
   int                   sign                /**< sign for weight, i.e. +1 to use right hand side or -1 to use left hand side */
)
{
   int i;
   SCIP_Real weight = 0.5 * sign;
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int rowlen;

   rowlen = SCIProwGetNNonz(row);
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   for( i = 0; i < rowlen; ++i )
   {
      SCIP_Real val;
      int probindex;

      probindex = SCIPcolGetVarProbindex(rowcols[i]);
      val = tmpcoefs[probindex];
      if( val == 0.0 )
      {
         nonzeroinds[(*nnz)++] = probindex;
      }

      val += weight * rowvals[i];
      tmpcoefs[probindex] = NONZERO(val);
   }

   if( sign == +1 )
   {
      *cutrhs += weight * SCIPfeasFloor(scip, SCIProwGetRhs(row) - SCIProwGetConstant(row));
   }
   else
   {
      assert(sign == -1);
      *cutrhs += weight * SCIPfeasCeil(scip, SCIProwGetLhs(row) - SCIProwGetConstant(row));
   }

   if( SCIProwGetRank(row) > *cutrank )
      *cutrank = SCIProwGetRank(row);
   *cutislocal = *cutislocal || SCIProwIsLocal(row);
}

/** add transformed integral row to aggregation with weight 0.5 */
static
void addTransRow(
   SCIP_Real*            tmpcoefs,           /**< array to add coefficients to */
   SCIP_Real*            cutrhs,             /**< pointer to add right hand side */
   int*                  nonzeroinds,        /**< array of non-zeros in the aggregation */
   int*                  nnz,                /**< pointer to update number of non-zeros */
   int*                  cutrank,            /**< pointer to update cut rank */
   SCIP_Bool*            cutislocal,         /**< pointer to update local flag */
   TRANSINTROW*          introw              /**< transformed integral row to add to the aggregation */
)
{
   int i;

   for( i = 0; i < introw->len; ++i )
   {
      int probindex = introw->varinds[i];
      SCIP_Real val = tmpcoefs[probindex];

      if( val == 0.0 )
      {
         nonzeroinds[(*nnz)++] = probindex;
      }

      val += 0.5 * introw->vals[i];
      tmpcoefs[probindex] = NONZERO(val);
   }

   *cutrhs += 0.5 * introw->rhs;

   *cutrank = MAX(*cutrank, introw->rank);
   *cutislocal = *cutislocal || introw->local;
}

/* calculates the cuts efficacy of cut */
static
SCIP_Real calcEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real             cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int                   cutnnz              /**< the number of non-zeros in the cut */
)
{
   SCIP_VAR** vars;
   SCIP_Real norm;
   SCIP_Real activity;
   int i;

   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);

   norm = SCIPgetVectorEfficacyNorm(scip, cutcoefs, cutnnz);
   vars = SCIPgetVars(scip);

   activity = 0.0;
   for( i = 0; i < cutnnz; ++i )
      activity += cutcoefs[i] * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);

   return (activity - cutrhs) / MAX(1e-6, norm);
}

/** computes maximal violation that can be achieved for zerohalf cuts where this row particiaptes */
static
SCIP_Real computeMaxViolation(
   MOD2_ROW*             row                 /**< mod 2 row */
)
{
   SCIP_Real viol;

   viol = 1.0 - row->slack;
   viol *= 0.5;

   return viol;
}

/** computes violation of zerohalf cut generated from given mod 2 row */
static
SCIP_Real computeViolation(
   MOD2_ROW*             row                 /**< mod 2 row */
)
{
   int i;
   SCIP_Real viol;

   viol = 1.0 - row->slack;

   for( i = 0; i < row->nnonzcols; ++i )
   {
      viol -= row->nonzcols[i]->solval;
   }

   viol *= 0.5;

   return viol;
}


/** remove rows that are (a) zero (b) identical to other rows (keep the one with smallest slack) (c) have slack greater
 * than 1 (d) for zero rows with 1 as rhs and slack less than 1, we can directly generate a cut and remove the row (Lemma 4)
 */
static
SCIP_RETCODE mod2matrixPreprocessRows(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   MOD2_MATRIX*          mod2matrix,         /**< the mod 2 matrix */
   GCG_ZEROHALFDATA*     zhdata,
   SCIP_Bool             allowlocal,          /**< should local cuts be allowed */
   GCG_CUTINDICES***     cutindices,
   int*                  ncutindices,
   int                   maxcuts
)
{
   int i;
   SCIP_HASHTABLE* rowtable;

   SCIP_CALL( SCIPhashtableCreate(&rowtable, SCIPblkmem(scip), mod2matrix->nrows,
                                  SCIPhashGetKeyStandard, rowsEqual, rowGetSignature, NULL) );

   for( i = 0; i < mod2matrix->nrows; )
   {
      MOD2_ROW* row = mod2matrix->rows[i];
      row->pos = i;

      checkRow(row);

      assert(row->nnonzcols == 0 || row->nonzcols != NULL);

      if( (row->nnonzcols == 0 && row->rhs == 0) || computeMaxViolation(row) < zhdata->minviol )
      { /* (a) and (c) */
         zhdata->nreductions += row->nnonzcols;
         SCIP_CALL( mod2matrixRemoveRow(scip, mod2matrix, row) );
      }
      else if( row->nnonzcols > 0 )
      { /* (b) */
         MOD2_ROW* identicalrow;
         identicalrow = (MOD2_ROW*)SCIPhashtableRetrieve(rowtable, (void*)row);
         if( identicalrow != NULL )
         {
            assert(identicalrow != row);
            assert(identicalrow->nnonzcols == 0 || identicalrow->nonzcols != NULL);

            checkRow(identicalrow);

            /* row is identical to other row; only keep the one with smaller slack */
            if( identicalrow->slack <= row->slack )
            {
               SCIP_CALL( mod2matrixRemoveRow(scip, mod2matrix, row) );
            }
            else
            {
               assert(SCIPhashtableExists(rowtable, (void*)identicalrow));

               SCIP_CALL( SCIPhashtableRemove(rowtable, (void*)identicalrow) );
               assert(!SCIPhashtableExists(rowtable, (void*)identicalrow));

               SCIP_CALL( SCIPhashtableInsert(rowtable, (void*)row) );

               SCIPswapPointers((void**) &mod2matrix->rows[row->pos], (void**) &mod2matrix->rows[identicalrow->pos]);
               SCIPswapInts(&row->pos, &identicalrow->pos);

               assert(mod2matrix->rows[row->pos] == row && mod2matrix->rows[identicalrow->pos] == identicalrow);
               assert(identicalrow->pos == i);
               assert(row->pos < i);

               SCIP_CALL( mod2matrixRemoveRow(scip, mod2matrix, identicalrow) );
            }
         }
         else
         {
            SCIP_CALL( SCIPhashtableInsert(rowtable, (void*)row) );
            ++i;
         }
      }
      else
      {
         GCG_CUTINDICES* cutindex;
         /* (d) */
         assert(row->nnonzcols == 0 && row->rhs == 1 && SCIPisLT(scip, row->slack, 1.0));

         /* OG: SCIP_CALL( generateZerohalfCut(scip, sol, mod2matrix, sepa, sepadata, allowlocal, row) ); */
         // transfer the indices associated with this cut to cutindices
         if( *ncutindices < maxcuts )
         {
            SCIP_CALL(GCGcreateCutIndicesFromRowindex(scip, &cutindex, row->nrowinds, row->rowinds) );
            (*cutindices)[*ncutindices] = cutindex;
            (*ncutindices)++;
         }


         if( zhdata->infeasible )
            goto TERMINATE;

         SCIP_CALL( mod2matrixRemoveRow(scip, mod2matrix, row) );
         ++i;
      }
   }
   TERMINATE:
   SCIPhashtableFree(&rowtable);

   return SCIP_OKAY;
}

/** add a mod2 row to another one */
static
SCIP_RETCODE mod2rowAddRow(
   SCIP*                 scip,               /**< scip data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory shell */
   MOD2_MATRIX*          mod2matrix,         /**< mod 2 matrix */
   MOD2_ROW*             row,                /**< mod 2 row */
   MOD2_ROW*             rowtoadd            /**< mod 2 row that is added to the other mod 2 row */
)
{
   SCIP_Shortbool* contained;
   int i;
   int j;
   int k;
   int nnewentries;
   int nlprows;
   MOD2_COL** newnonzcols;
   SCIP_Real newslack;
   SCIP* origscip;

   checkRow(row);
   checkRow(rowtoadd);

   assert(row->nnonzcols == 0 || row->nonzcols != NULL);
   assert(rowtoadd->nnonzcols == 0 || rowtoadd->nonzcols != NULL);
   origscip = GCGgetOriginalprob(scip);

   nlprows = GCGgetNMasterConss(origscip); //SCIPgetNLPRows(scip);
   row->rhs ^= rowtoadd->rhs;

   newslack = row->slack + rowtoadd->slack;
   blkmem = SCIPblkmem(scip);

   if( SCIPisZero(scip, row->slack) && !SCIPisZero(scip, newslack) )
      --mod2matrix->nzeroslackrows;

   row->slack = newslack;

   {
      /* the maximum index return by the UNIQUE_INDEX macro is 3 times
       * the maximum index value in the ROWINDEX struct. The index value could
       * be the lp position of an original row or the index of a transformed row.
       * Hence we need to allocate 3 times the maximum of these two possible
       * index types.
       */
      int allocsize = 3 * MAX(nlprows, mod2matrix->ntransintrows);
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &contained, allocsize) );
   }

   /* remember entries that are in the row to add */
   for( i = 0; i < rowtoadd->nrowinds; ++i )
   {
      contained[UNIQUE_INDEX(rowtoadd->rowinds[i])] = 1;
   }

   /* remove the entries that are in both rows from the row (1 + 1 = 0 (mod 2)) */
   nnewentries = rowtoadd->nrowinds;
   for( i = 0; i < row->nrowinds; )
   {
      if( contained[UNIQUE_INDEX(row->rowinds[i])] )
      {
         --nnewentries;
         contained[UNIQUE_INDEX(row->rowinds[i])] = 0;
         --row->nrowinds;
         row->rowinds[i] = row->rowinds[row->nrowinds];
      }
      else
      {
         ++i;
      }
   }

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &row->rowinds, &row->rowindssize, row->nrowinds + nnewentries) );

   /* add remaining entries of row to add */
   for ( i = 0; i < rowtoadd->nrowinds; ++i )
   {
      if( contained[UNIQUE_INDEX(rowtoadd->rowinds[i])] )
      {
         contained[UNIQUE_INDEX(rowtoadd->rowinds[i])] = 0;
         row->rowinds[row->nrowinds++] = rowtoadd->rowinds[i];
      }
   }

   SCIPfreeCleanBufferArray(scip, &contained);

   SCIP_CALL( SCIPallocBufferArray(scip, &newnonzcols, row->nnonzcols + rowtoadd->nnonzcols) );

   i = 0;
   j = 0;
   k = 0;
   row->maxsolval = 0.0;

   /* since columns are sorted we can merge them */
   while( i < row->nnonzcols && j < rowtoadd->nnonzcols )
   {
      if( row->nonzcols[i] == rowtoadd->nonzcols[j] )
      {
         SCIP_CALL( mod2colUnlinkRow(row->nonzcols[i], row) );
         ++i;
         ++j;
      }
      else if( row->nonzcols[i]->index < rowtoadd->nonzcols[j]->index )
      {
         row->maxsolval = MAX(row->maxsolval, row->nonzcols[i]->solval);
         newnonzcols[k++] = row->nonzcols[i++];
      }
      else
      {
         SCIP_CALL( mod2colLinkRow(blkmem, rowtoadd->nonzcols[j], row) );
         newnonzcols[k++] = rowtoadd->nonzcols[j++];
      }
   }

   while( i < row->nnonzcols )
   {
      row->maxsolval = MAX(row->maxsolval, row->nonzcols[i]->solval);
      newnonzcols[k++] = row->nonzcols[i++];
   }

   while( j <  rowtoadd->nnonzcols )
   {
      SCIP_CALL( mod2colLinkRow(blkmem, rowtoadd->nonzcols[j], row) );
      newnonzcols[k++] = rowtoadd->nonzcols[j++];
   }

   row->nnonzcols = k;
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &row->nonzcols, &row->nonzcolssize, row->nnonzcols) );
   BMScopyMemoryArray(row->nonzcols, newnonzcols, row->nnonzcols);

   SCIPfreeBufferArray(scip, &newnonzcols);

   assert(row->nnonzcols == 0 || row->nonzcols != NULL);
   checkRow(row);
   checkRow(rowtoadd);

   return SCIP_OKAY;
}

/* --------------------------------------------------------------------------------------------------------------------
 * callback methods of separator
 * -------------------------------------------------------------------------------------------------------------------- */

// SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, (unsigned int)sepadata->initseed, TRUE) );
// SCIPfreeRandom(scip, &sepadata->randnumgen);


/** perform the zerohalf cut separation */
SCIP_RETCODE GCGselectConstraintsZeroHalf(
   SCIP*                 scip,
   SCIP_SOL*             sol,
   SCIP_Bool             allowlocal,
   int                   depth,               /* current depth */
   GCG_ZEROHALFDATA*     zhdata,
   int                   ncalls,
   int                   maxcuts,
   GCG_CUTINDICES***     cutindices,
   int*                  ncutindices
)
{
   int i;
   int k;

   SCIP_Real maxslack;
   MOD2_MATRIX mod2matrix;
   MOD2_ROW** nonzrows;


   assert(zhdata != NULL);

   {

      /* only call the zerohalf cut separator a given number of times at each node */
      if( (depth == 0 && zhdata->maxroundsroot >= 0 && ncalls >= zhdata->maxroundsroot)
          || (depth > 0 && zhdata->maxrounds >= 0 && ncalls >= zhdata->maxrounds) )
         return SCIP_OKAY;

      maxslack = depth == 0 ? zhdata->maxslackroot : zhdata->maxslack;
      maxslack += 2 * SCIPfeastol(scip);
   }


   *ncutindices = 0;
   zhdata->infeasible = FALSE;

   SCIP_CALL( buildMod2Matrix(scip, sol, zhdata, SCIPblkmem(scip), &mod2matrix, allowlocal, maxslack) );

   SCIPdebugMsg(scip, "built mod2 matrix (%i rows, %i cols)\n", mod2matrix.nrows, mod2matrix.ncols);

   SCIP_CALL( SCIPallocBufferArray(scip, &nonzrows, mod2matrix.nrows) );

   for( k = 0; k < MAXREDUCTIONROUNDS; ++k )
   {
      int ncancel;

      zhdata->nreductions = 0;

      assert(mod2matrix.nzeroslackrows <= mod2matrix.nrows);
      SCIP_CALL( mod2matrixPreprocessRows(scip, sol, &mod2matrix, zhdata, allowlocal, cutindices, ncutindices, maxcuts) );
      assert(mod2matrix.nzeroslackrows <= mod2matrix.nrows);

      SCIPdebugMsg(scip, "preprocessed rows (%i rows, %i cols, %i cuts) \n", mod2matrix.nrows, mod2matrix.ncols,
                   *ncutindices);

      if( mod2matrix.nrows == 0 )
         break;

      if( *ncutindices >= maxcuts )
      {
         SCIPdebugMsg(scip, "enough cuts, stopping (%i rows, %i cols)\n", mod2matrix.nrows, mod2matrix.ncols);
         break;
      }

      SCIP_CALL( mod2matrixPreprocessColumns(scip, &mod2matrix, zhdata) );

      SCIPdebugMsg(scip, "preprocessed columns (%i rows, %i cols)\n", mod2matrix.nrows, mod2matrix.ncols);

      ncancel = mod2matrix.nrows;
      if( ncancel > 100 )
      {
         ncancel = 100;
         SCIPselectPtr((void**) mod2matrix.rows, compareRowSlack, ncancel, mod2matrix.nrows);
      }

      SCIPsortPtr((void**) mod2matrix.rows, compareRowSlack, ncancel);

      if( mod2matrix.ncols == 0 )
         break;

      assert(mod2matrix.nzeroslackrows <= mod2matrix.nrows);

      /* apply Prop5 */
      for( i = 0; i < ncancel; ++i )
      {
         int j;
         MOD2_COL* col = NULL;
         MOD2_ROW* row = mod2matrix.rows[i];

         if( SCIPisPositive(scip, row->slack) || row->nnonzcols == 0 )
            continue;

         SCIPdebugMsg(scip, "processing row %i/%i (%i/%i cuts)\n", i, mod2matrix.nrows, *ncutindices, maxcuts);

         for( j = 0; j < row->nnonzcols; ++j )
         {
            if( row->nonzcols[j]->solval == row->maxsolval ) /*lint !e777*/
            {
               col = row->nonzcols[j];
               break;
            }
         }

         assert( col != NULL );

         {
            int nslots;
            int nnonzrows;
            MOD2_ROW** rows;

            ++zhdata->nreductions;

            nnonzrows = 0;
            /* cppcheck-suppress nullPointer */
            nslots = SCIPhashsetGetNSlots(col->nonzrows);
            /* cppcheck-suppress nullPointer */
            rows = (MOD2_ROW**) SCIPhashsetGetSlots(col->nonzrows);

            for( j = 0; j < nslots; ++j )
            {
               if( rows[j] != NULL && rows[j] != row )
                  nonzrows[nnonzrows++] = rows[j];
            }

            for( j = 0; j < nnonzrows; ++j )
            {
               SCIP_CALL( mod2rowAddRow(scip, SCIPblkmem(scip), &mod2matrix, nonzrows[j], row) );
            }

            /* cppcheck-suppress nullPointer */
            row->slack = col->solval;
            --mod2matrix.nzeroslackrows;

            mod2matrixRemoveCol(scip, &mod2matrix, col);
         }
      }

      SCIPdebugMsg(scip, "applied proposition five (%i rows, %i cols)\n", mod2matrix.nrows, mod2matrix.ncols);

      if( zhdata->nreductions == 0 )
      {
         SCIPdebugMsg(scip, "no change, stopping (%i rows, %i cols)\n", mod2matrix.nrows, mod2matrix.ncols);
         break;
      }
   }

   for( i = 0; i < mod2matrix.nrows && *ncutindices < maxcuts; ++i )
   {
      GCG_CUTINDICES* cutindex;
      MOD2_ROW* row = mod2matrix.rows[i];

      if( computeMaxViolation(row) < zhdata->minviol ) // computeViolation(row)
         break;

      if( row->rhs == 0 )
         continue;

      /* OG: SCIP_CALL( generateZerohalfCut(scip, sol, &mod2matrix, sepa, sepadata, allowlocal, row) ); */
      // transfer the indices associated with this cut to cutindices
      SCIP_CALL( GCGcreateCutIndicesFromRowindex(scip, &cutindex, row->nrowinds, row->rowinds) );
      (*cutindices)[*ncutindices] = cutindex;
      (*ncutindices)++;

   }

   SCIPdebugMsg(scip, "total number of cuts found: %i\n", *ncutindices);

   /* If cuts where found we apply a filtering procedure using the scores and the orthogonalities,
    * similar to the sepastore. We only add the cuts that make it through this process and discard
    * the rest.
    */

   SCIPfreeBufferArray(scip, &nonzrows);

   destroyMod2Matrix(scip, &mod2matrix);

   return SCIP_OKAY;
}


/** creates the zerohalf separator and includes it in SCIP */
/*SCIP_RETCODE SCIPincludeSepaZerohalf(
   SCIP*                 scip
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;


   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
                                   SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpZerohalf, sepaExecsolZerohalf, sepadata) );

   assert(sepa != NULL);


   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyZerohalf) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeZerohalf) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolZerohalf) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolZerohalf) );


   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/maxrounds",
      "maximal number of zerohalf separation rounds per node (-1: unlimited)",
      &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/maxroundsroot",
      "maximal number of zerohalf separation rounds in the root node (-1: unlimited)",
      &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/maxsepacuts",
      "maximal number of zerohalf cuts separated per separation round",
      &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/initseed",
      "initial seed used for random tie-breaking in cut selection",
      &sepadata->initseed, FALSE, DEFAULT_INITSEED, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/maxsepacutsroot",
      "maximal number of zerohalf cuts separated per separation round in the root node",
      &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/maxcutcands",
      "maximal number of zerohalf cuts considered per separation round",
      &sepadata->maxcutcands, FALSE, DEFAULT_MAXCUTCANDS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/maxslack",
      "maximal slack of rows to be used in aggregation",
      &sepadata->maxslack, TRUE, DEFAULT_MAXSLACK, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/maxslackroot",
      "maximal slack of rows to be used in aggregation in the root node",
      &sepadata->maxslackroot, TRUE, DEFAULT_MAXSLACKROOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/goodscore",
      "threshold for score of cut relative to best score to be considered good, so that less strict filtering is applied",
      &sepadata->goodscore, TRUE, DEFAULT_GOODSCORE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/badscore",
      "threshold for score of cut relative to best score to be discarded",
      &sepadata->badscore, TRUE, DEFAULT_BADSCORE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/objparalweight",
      "weight of objective parallelism in cut score calculation",
      &sepadata->objparalweight, TRUE, DEFAULT_OBJPARALWEIGHT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/efficacyweight",
      "weight of efficacy in cut score calculation",
      &sepadata->efficacyweight, TRUE, DEFAULT_EFFICACYWEIGHT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/dircutoffdistweight",
      "weight of directed cutoff distance in cut score calculation",
      &sepadata->dircutoffdistweight, TRUE, DEFAULT_DIRCUTOFFDISTWEIGHT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/goodmaxparall",
      "maximum parallelism for good cuts",
      &sepadata->goodmaxparall, TRUE, DEFAULT_GOODMAXPARALL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/maxparall",
      "maximum parallelism for non-good cuts",
      &sepadata->maxparall, TRUE, DEFAULT_MAXPARALL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/minviol",
      "minimal violation to generate zerohalfcut for",
      &sepadata->minviol, TRUE, DEFAULT_MINVIOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
                               "separating/" SEPA_NAME "/dynamiccuts",
      "should generated cuts be removed from the LP if they are no longer tight?",
      &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
                               "separating/" SEPA_NAME "/maxrowdensity",
      "maximal density of row to be used in aggregation",
      &sepadata->maxrowdensity, TRUE, DEFAULT_MAXROWDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
                              "separating/" SEPA_NAME "/densityoffset",
      "additional number of variables allowed in row on top of density",
      &sepadata->densityoffset, TRUE, DEFAULT_DENSITYOFFSET, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}*/
