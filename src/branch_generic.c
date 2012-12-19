/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//#define SCIP_DEBUG
/**@file   branch_generic.c
 * @ingroup BRANCHINGRULES
 * @brief  branching rule based on vanderbeck's generic branching scheme
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "branch_generic.h"
#include "relax_gcg.h"
#include "cons_masterbranch.h"
#include "cons_origbranch.h"
#include "pricer_gcg.h"
#include "scip/cons_linear.h"
#include "type_branchgcg.h"
#include "pub_gcgvar.h"
#include "event_genericbranchvaradd.h"

#include <stdio.h>
#include <stdlib.h>

#define BRANCHRULE_NAME          "generic"
#define BRANCHRULE_DESC          "generic branching rule by Vanderbeck"
#define BRANCHRULE_PRIORITY      999999
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/** @todo enum and struct, extract from event handler, move to own file */
typedef SCIP_Real ComponentBoundSequence[3];  // [[comp], [sense], [bound]]    sense=1 means >=, sense=0 means <

struct GCG_BranchData
{
   ComponentBoundSequence**   C;             /**< S[k] bound sequence for block k */ //!!! sort of each C[i]=S[i] is important !!!
   int*                  sequencesizes;      /**< number of bounds in S[k] */
   int                Csize;
   ComponentBoundSequence*   S;              /**< component bound sequence which induce the child branching constraints */
   int                Ssize;
   int                blocknr;               /**< number of block branching was performed */
   int                childnr;
   SCIP_Real          lhs;
   int                nchildNodes;
   SCIP_Real*         childlhs;
   SCIP_CONS*         mastercons;            /**< constraint enforcing the branching restriction in the master problem */
   GCG_BRANCHDATA**   childbranchdatas;
   ComponentBoundSequence* consS;            /**< component bound sequence which induce the current branching constraint */
   int                consSsize;
   int                consblocknr;
};

struct GCG_Strip
{
   SCIP* scip;
   SCIP_VAR*          mastervar;             /**< master variable */
   SCIP_Real          mastervarValue;
   int                blocknr;               /**< number of the block in which the strip belong */
   SCIP_Real*         generator;             /**< corresponding generator to the mastervar */
   SCIP_Bool*         compisinteger;         /**< ? is comp with integer origvars? */
   int                generatorsize;
   ComponentBoundSequence**   C;             /**< often NULL, only needed for ptrilocomp */
   int                Csize;
   int*               sequencesizes;
};



/** set of component bounds in separate */
struct GCG_Record
{
   ComponentBoundSequence**   record;             /**< returnvalue of separte function */
   int                recordsize;
   int*               sequencesizes;
};
typedef struct GCG_Record GCG_RECORD;

/*
 * Callback methods
 */

/* define not used callback as NULL*/
#define branchCopyGeneric NULL
#define branchFreeGeneric NULL
#define branchExitGeneric NULL
#define branchInitsolGeneric NULL
#define branchExitsolGeneric NULL

/*
 * branching specific interface methods
 */

/** method for calculating the generator of mastervar*/
static
SCIP_RETCODE getGenerators(SCIP* scip, SCIP_Real** generator, int* generatorsize, SCIP_Bool** compisinteger, int blocknr, SCIP_VAR** mastervars, int nmastervars, SCIP_VAR* mastervar)
{
   int i;
   int j;
   int k;
   SCIP_VAR** origvarsunion;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;
   int nvarsinblock;

   i = 0;
   j = 0;
   k = 0;
   *generatorsize = 0;
   nvarsinblock = 0;
   origvarsunion = NULL;
   assert(mastervars != NULL);

   //	SCIPdebugMessage("get generator, block = %d, nvars = %d \n", blocknr, nmastervars);

   for(i=0; i<nmastervars; ++i)
   {
      origvars = GCGmasterVarGetOrigvars(mastervars[i]);
      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);

      if(blocknr != GCGvarGetBlock(mastervars[i]))
         continue;
      else
         ++nvarsinblock;
      if(*generatorsize == 0 && norigvars > 0)
      {
         *generatorsize = norigvars;
         SCIP_CALL( SCIPallocMemoryArray(scip, generator, *generatorsize) );
         SCIP_CALL( SCIPallocMemoryArray(scip, compisinteger, *generatorsize) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &origvarsunion, *generatorsize) );
         for(j=0; j<*generatorsize; ++j)
         {
            origvarsunion[j] = origvars[j];
            (*generator)[j] = 0;
            (*compisinteger)[j] = TRUE;
            if(SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_IMPLINT)
               (*compisinteger)[j] = FALSE;
         }
      }
      else
      {
         for(j=0; j<norigvars; ++j)
         {
            int oldgeneratorsize;

            oldgeneratorsize = *generatorsize;

            for(k=0; k<oldgeneratorsize; ++k)
            {
               if(origvarsunion[k] == origvars[j])
               {
                  break;
               }
               if(k == oldgeneratorsize-1) //norigvars-1)
               {
                  ++(*generatorsize);
                  SCIP_CALL( SCIPreallocMemoryArray(scip, generator, *generatorsize) );
                  SCIP_CALL( SCIPreallocMemoryArray(scip, compisinteger, *generatorsize) );
                  SCIP_CALL( SCIPreallocMemoryArray(scip, &origvarsunion, *generatorsize) );
                  origvarsunion[*generatorsize-1] = origvars[j];
                  (*generator)[*generatorsize-1] = 0;
                  (*compisinteger)[(*generatorsize)-1] = TRUE;
                  if(SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_IMPLINT)
                     (*compisinteger)[(*generatorsize)-1] = FALSE;
               }
            }
         }
      }
   }

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   for(i=0; i<norigvars; ++i)
   {
      for(j=0; j<*generatorsize; ++j)
      {
         if(origvarsunion[j]==origvars[i])
         {
            if(SCIPvarGetType(origvars[i]) == SCIP_VARTYPE_CONTINUOUS)
               (*compisinteger)[j] = FALSE;
            if(!SCIPisZero(scip, origvals[i]))
               (*generator)[j] = origvals[i];
            break;
         }
      }
   }

   SCIPfreeMemoryArray(scip, &origvarsunion);

   return SCIP_OKAY;
}

/** method for calculating the median over all fractional components values if its the minimum return ceil(arithm middle)*/
static
SCIP_Real GetMedian(SCIP* scip, SCIP_Real* array, int arraysize, SCIP_Real min)
{
   SCIP_Real Median;
   SCIP_Real swap;
   int l;
   int r;
   int i;
   int j;
   int MedianIndex;
   SCIP_Real arithmMiddle;

   r = arraysize -1;
   l = 0;
   arithmMiddle = 0;

   if( arraysize & 1)
      MedianIndex = arraysize/2;
   else
      MedianIndex = arraysize/2 -1;

   while(l < r-1)
   {
      Median = array[ MedianIndex ];
      i = l;
      j = r;
      do
      {
         while( SCIPisLT(scip, array[i], Median) )
            ++i;
         while( SCIPisGT(scip, array[j], Median) )
            --j;
         if( i <=j )
         {
            swap = array[i];
            array[i] = array[j];
            array[j] = swap;
            ++i;
            --j;
         }
      } while( i <=j );
      if( j < MedianIndex )
         l = i;
      if( i > MedianIndex )
         r = j;
   }
   Median = array[ MedianIndex ];

   if(  SCIPisEQ(scip, Median, min) )
   {
      for(i=0; i<arraysize; ++i)
         arithmMiddle += array[i];

      arithmMiddle /= arraysize;
      Median = SCIPceil(scip, arithmMiddle);
   }

   return Median;
}

// comparefunction for lexicographical sort
static
SCIP_DECL_SORTPTRCOMP(ptrcomp)
{
   struct GCG_Strip* strip1;
   struct GCG_Strip* strip2;
   int i;

   strip1 = (struct GCG_Strip*) elem1;
   strip2 = (struct GCG_Strip*) elem2;

   i = 0;

   assert(strip1->generatorsize == strip2->generatorsize);

   for( i=0; i< strip1->generatorsize; ++i)
   {
      if( strip1->generator[i] > strip2->generator[i] )
         return -1;
      if( strip1->generator[i] < strip2->generator[i] )
         return 1;
   }

   return 0;
}

// lexicographical sort using scipsort
// !!! changes the array
static
SCIP_RETCODE LexicographicSort( struct GCG_Strip** array, int arraysize)
{

   SCIPdebugMessage("Lexicographic sorting\n");

   SCIPsortPtr((void**)array, ptrcomp, arraysize );

   //change array
   return SCIP_OKAY;
}


// compare function for ILO: returns 1 if bd1 < bd2 else -1
static
int ILOcomp( SCIP* scip, struct GCG_Strip* strip1, struct GCG_Strip* strip2, ComponentBoundSequence** C, int NBoundsequences, int* sequencesizes, int p) // ComponentBoundSequence* S, int Ssize, int* IndexSet, int indexsetsize)
{
   int i;
   int ivalue;
   int j;
   int k;
   int l;
   int Nupper;
   int Nlower;
   SCIP_Bool returnvalue;
   ComponentBoundSequence** CopyC;
   int* newsequencesizes;

   i = -1;
   j = 0;
   k = 0;
   l = 0;
   Nupper = 0;
   Nlower = 0;

   //lexicographic Order ?
   if( C == NULL || NBoundsequences <= 1 )
      return (*ptrcomp)( strip1, strip2);// == -1);

   assert(C != NULL);
   assert(NBoundsequences > 0);

   //find i which is in all S in C on position p
   while( sequencesizes[k] < p )
   {
      ++k;
      assert(k < NBoundsequences);
   }
   i = C[k][p-1][0];
   assert(i < strip1->generatorsize && i < strip2->generatorsize);
   ivalue = C[k][p-1][2];

   assert(i >= 0);

   //calculate subset of C
   for(j=0; j< NBoundsequences; ++j)
   {
      if(sequencesizes[j] >= p)
      {
         assert(C[j][p-1][0] == i);
         if(C[j][p-1][1] == 1)
            ++Nupper;
         else
            ++Nlower;
      }
   }

   if( strip1->generator[i] >= ivalue && strip2->generator[i] >= ivalue )
   {
      k=0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nupper) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nupper) );
      for(j=0; j< NBoundsequences; ++j)
      {
         if( sequencesizes[j] >= p )
            assert(C[j][p-1][0] == i);

         if(sequencesizes[j] >= p && C[j][p-1][1] == 1)
         {
            CopyC[k] = NULL;
            SCIP_CALL( SCIPallocMemoryArray(scip, &(CopyC[k]), sequencesizes[j]) );
            for(l=0; l< sequencesizes[j]; ++l)
            {
               CopyC[k][l][0] = C[j][l][0];
               CopyC[k][l][1] = C[j][l][1];
               CopyC[k][l][2] = C[j][l][2];
            }
            newsequencesizes[k] = sequencesizes[j];
            ++k;
         }
      }

      if(k != Nupper)
      {
         SCIPdebugMessage("k = %d, Nupper+1 =%d\n", k, Nupper+1);
      }

      if(Nupper != 0)
         assert( k == Nupper );

      returnvalue = ILOcomp( scip, strip1, strip2, CopyC, Nupper, newsequencesizes, p+1);

      for(j=0; j< Nupper; ++j)
      {

         SCIPfreeMemoryArray(scip, &(CopyC[j]) );
      }
      SCIPfreeMemoryArray(scip, &newsequencesizes);
      SCIPfreeMemoryArray(scip, &CopyC);

      return returnvalue;
   }


   if( strip1->generator[i] < ivalue && strip2->generator[i] < ivalue )
   {
      k=0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nlower) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nlower) );
      for(j=0; j< NBoundsequences; ++j)
      {
         if( sequencesizes[j] >= p )
            assert(C[j][p-1][0] == i);

         if( sequencesizes[j] >= p && C[j][p-1][1] != 1)
         {
            CopyC[k] = NULL;
            SCIP_CALL( SCIPallocMemoryArray(scip, &(CopyC[k]), sequencesizes[j]) );
            for(l=0; l< sequencesizes[j]; ++l)
            {
               CopyC[k][l][0] = C[j][l][0];
               CopyC[k][l][1] = C[j][l][1];
               CopyC[k][l][2] = C[j][l][2];
            }
            newsequencesizes[k] = sequencesizes[j];
            ++k;
         }
      }

      if(k != Nlower)
      {
         SCIPdebugMessage("k = %d, Nlower+1 =%d\n", k, Nlower+1);
      }

      if(Nlower != 0)
         assert( k == Nlower);

      returnvalue = ILOcomp( scip, strip1, strip2, CopyC, Nlower, newsequencesizes, p+1);// S, Ssize, IndexSet, indexsetsize, p+1);

      for(j=0; j< Nlower; ++j)
      {

         SCIPfreeMemoryArray(scip, &(CopyC[j]) );
      }
      SCIPfreeMemoryArray(scip, &newsequencesizes);
      SCIPfreeMemoryArray(scip, &CopyC);

      return returnvalue;
   }
   if( strip1->generator[i] > strip2->generator[i])
      return 1;
   else
      return -1;
}

// comparefunction for induced lexicographical sort
static
SCIP_DECL_SORTPTRCOMP(ptrilocomp)
{
   struct GCG_Strip* strip1;
   struct GCG_Strip* strip2;
   int returnvalue;

   strip1 = (struct GCG_Strip*) elem1;
   strip2 = (struct GCG_Strip*) elem2;

   returnvalue = ILOcomp( strip1->scip, strip1, strip2, strip1->C, strip1->Csize, strip1->sequencesizes, 1); //NULL, 0, strip1->IndexSet, strip1->generatorsize, 1);

   return returnvalue;
}

// induced lexicographical sort
static
SCIP_RETCODE InducedLexicographicSort( SCIP* scip, struct GCG_Strip** array, int arraysize, ComponentBoundSequence** C, int NBoundsequences, int* sequencesizes )
{
   int i;

   SCIPdebugMessage("Induced Lexicographic sorting\n");

   if( NBoundsequences == 0 )
      return LexicographicSort( array, arraysize );
   assert( C!= NULL );

   assert(arraysize > 0);

   if(arraysize <= 1)
      return SCIP_OKAY;

   for( i=0; i<arraysize; ++i )
   {
      array[i]->scip = scip;
      array[i]->Csize = NBoundsequences;
      array[i]->sequencesizes = sequencesizes;
      array[i]->C = C;
   }

   SCIPsortPtr((void**)array, ptrilocomp, arraysize);

   return SCIP_OKAY;
}


// separation at the root node
static
SCIP_RETCODE Separate( SCIP* scip, struct GCG_Strip** F, int Fsize, int* IndexSet, int IndexSetSize, ComponentBoundSequence* S, int Ssize, struct GCG_Record** record )
{
   int i;
   int j;
   //int n;
   int k;
   int l;
   int Jsize;
   int* J;
   SCIP_Real median;
   SCIP_Real min;
   SCIP_Real max;
   int Fupper;
   int Flower;
   int* priority;
   struct GCG_Strip** copyF;
   ComponentBoundSequence* upperLowerS;
   SCIP_Real* alpha;
   SCIP_Real* compvalues;
   SCIP_Real  muF;
   SCIP_Real maxPriority;
   SCIP_Bool found;

   i = 0;
   j = 0;
   k = 0;
   l = 0;
   Jsize = 0;
   Fupper = 0;
   Flower = 0;
   max = 0;
   muF = 0;
   min = INT_MAX;
   maxPriority = INT_MIN;
   found = FALSE;
   priority = NULL;
   compvalues = NULL;
   J = NULL;
   copyF = NULL;
   upperLowerS = NULL;
   alpha = NULL;

   SCIPdebugMessage("Separate with ");

   if(Fsize == 0 || IndexSetSize == 0)
      return SCIP_OKAY;

   SCIPdebugMessage("Fsize = %d; Ssize = %d, IndexSetSize = %d\n", Fsize, Ssize, IndexSetSize);

   assert( F != NULL );
   assert( IndexSet != NULL );

   for(j=0; j<Fsize; ++j)
   {
      for( i=0; i<IndexSetSize; ++i)
      {
         assert(IndexSet[i] < F[j]->generatorsize);
         if(F[j]->generator[IndexSet[i]] > max )
         {
            max = F[j]->generator[IndexSet[i]];
         }
      }
   }

   if(max == 0)
      max = 1;

   SCIPdebugMessage("max = %g\n", max);

   for(j=0; j<Fsize; ++j)
      muF += max * F[j]->mastervarValue;
   SCIP_CALL( SCIPallocMemoryArray(scip, &alpha, IndexSetSize) );

   for(k=0; k<IndexSetSize; ++k)
   {
      ComponentBoundSequence* copyS;
      SCIP_Real mu_F;
      SCIP_Bool even;

      even = TRUE;
      mu_F = 0;
      i = IndexSet[k];
      copyS = NULL;
      alpha[k] = 0;

      assert(i < F[0]->generatorsize);

      if(!F[0]->compisinteger[i])
         continue;

      for(j=0; j<Fsize; ++j)
      {
         assert(i < F[j]->generatorsize);
         alpha[k] += F[j]->generator[i] * F[j]->mastervarValue;
         if(F[j]->generator[i] != 0 && F[j]->mastervarValue != 0)
         {
            //		SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
            //		SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
         }
      }
      if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
      {
         ++Jsize;
         //			SCIPdebugMessage("alpha[%d] = %g\n", k, alpha[k]);
      }
      if( SCIPisGT(scip, alpha[k] - SCIPfloor(scip, alpha[k]), 0) )
      {
         found = TRUE;

         //add to record
         ++Ssize;
         SCIP_CALL( SCIPallocMemoryArray(scip, &copyS, Ssize) );
         for(l=0; l < Ssize-1; ++l)
            for(j=0; j<3; ++j)
               copyS[l][j] = S[l][j];

         copyS[Ssize-1][0] = i;
         copyS[Ssize-1][1] = 1;

         SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
         for(l=0; l<Fsize; ++l)
         {
            compvalues[l] = F[l]->generator[i];
            if( SCIPisLT(scip, compvalues[l], min) )
               min = compvalues[l];
         }

         median = GetMedian(scip, compvalues, Fsize, min);
         j = 0;

         do
         {
            mu_F = 0;
            if(even)
            {
               median = median+j;
               even = FALSE;
            }
            else
            {
               median = median-j;
               even = TRUE;
            }

            for(l=0; l<Fsize; ++l)
            {
               if( SCIPisGE(scip, F[l]->generator[i], median) )
                  mu_F += F[l]->mastervarValue;
            }
            ++j;

         }while( SCIPisEQ(scip, mu_F - SCIPfloor(scip, mu_F), 0) );

         copyS[Ssize-1][2] = median;

         SCIPfreeBufferArray(scip, &compvalues);
         compvalues = NULL;

         (*record)->recordsize++;
         if((*record)->recordsize == 1 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
         }
         else
         {
            SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
         }
         (*record)->record[(*record)->recordsize-1] = copyS;
         (*record)->sequencesizes[(*record)->recordsize-1] = Ssize;
         --Ssize;
      }
   }

   if(found)
   {
      SCIPfreeMemoryArray(scip, &alpha);
      SCIPdebugMessage("one S found with size %d\n", (*record)->sequencesizes[(*record)->recordsize-1]);

      return SCIP_OKAY;
   }

   //discriminating components
   SCIP_CALL( SCIPallocMemoryArray(scip, &J, Jsize) );
   j=0;
   for(k=0; k<IndexSetSize; ++k)
   {
      if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
      {
         J[j] = IndexSet[k];
         ++j;
      }
   }

   //compute priority  (max-min)
   SCIP_CALL( SCIPallocMemoryArray(scip, &priority, Jsize) );
   for(j=0; j<Jsize; ++j)
   {
      int maxcomp;
      int mincomp;

      maxcomp = INT_MIN;
      mincomp = INT_MAX;

      for(l=0; l<Fsize; ++l)
      {
         assert(J[j] < F[l]->generatorsize);

         if(F[l]->generator[J[j]] > maxcomp)
            maxcomp = F[l]->generator[J[j]];

         if(F[l]->generator[J[j]] < mincomp)
            mincomp = F[l]->generator[J[j]];
      }
      //assert(J[j] < Jsize);
      priority[j] = maxcomp-mincomp;
   }

   //Partition
   SCIPdebugMessage("Partitioning\n");
   do{
      min = INT_MAX;
      maxPriority = INT_MIN;

      //max-min priority
      for(j=0; j<Jsize; ++j)
      {
         if(priority[j] > maxPriority && F[0]->compisinteger[J[j]])
         {
            maxPriority = priority[j];
            i = J[j];
         }
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
      for(l=0; l<Fsize; ++l)
      {
         compvalues[l] = F[l]->generator[i];
         if( SCIPisLT(scip, compvalues[l], min) )
            min = compvalues[l];
      }
      median = GetMedian(scip, compvalues, Fsize, min);
      SCIPfreeBufferArray(scip, &compvalues);

      if( !SCIPisEQ(scip, median, 0) )
      {
         SCIPdebugMessage("median = %g\n", median);
         SCIPdebugMessage("min = %g\n", min);
      }

      if( SCIPisEQ(scip, median, min) )
      {
         //here with max-min priority
         for(j=0; j<Jsize; ++j)
         {
            if( i == J[j])
            {
               J[j] = J[Jsize-1];
               break;
            }
         }
         --Jsize;

      }
      assert(Jsize>=0);
   }while( SCIPisEQ(scip, median, min) );

   ++Ssize;
   SCIP_CALL( SCIPallocMemoryArray(scip, &upperLowerS, Ssize) );
   for(l=0; l < Ssize-1; ++l)
      for(j=0; j<3; ++j)
         upperLowerS[l][j] = S[l][j];

   upperLowerS[Ssize-1][0] = i;
   upperLowerS[Ssize-1][1] = 1;
   upperLowerS[Ssize-1][2] = median;

   for(k=0; k<Fsize; ++k)
   {
      if( SCIPisGE(scip, F[k]->generator[i], median) )
         ++Fupper;
      else
         ++Flower;
   }

   //choose smallest partition
   if( (Flower <= Fupper && Flower > 0) || Fupper <= 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
      j = 0;
      for(k=0; k<Fsize; ++k)
      {
         if( SCIPisLT(scip, F[k]->generator[i], median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Flower;
   }
   else
   {
      upperLowerS[Ssize-1][1] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
      j = 0;
      for(k=0; k<Fsize; ++k)
      {
         if( SCIPisGE(scip, F[k]->generator[i], median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Fupper;
   }

   assert(j < Fsize+1);

   Separate( scip, copyF, Fsize, J, Jsize, upperLowerS, Ssize, record );


   SCIPfreeMemoryArray(scip, &copyF);
   SCIPfreeMemoryArray(scip, &upperLowerS);
   SCIPfreeMemoryArray(scip, &priority);
   SCIPfreeMemoryArray(scip, &J);
   SCIPfreeMemoryArray(scip, &alpha);

   //	return record;
   return SCIP_OKAY;
}

// choose a component bound sequence
static
SCIP_RETCODE ChoseS( SCIP* scip, struct GCG_Record** record, ComponentBoundSequence** S, int* Ssize )
{
   int minSizeOfMaxPriority;  //needed if the last comp priority is euqal to the one in other bound sequences
   int maxPriority;
   int i;
   int j;
   int Index;

   minSizeOfMaxPriority = INT_MAX;
   maxPriority = INT_MIN;
   i = 0;
   j = 0;
   Index = -1;

   SCIPdebugMessage("Chose S \n");

   assert((*record)->recordsize > 0);

   for( i=0; i< (*record)->recordsize; ++i )
   {
      assert((*record)->sequencesizes != NULL );
      assert((*record)->sequencesizes[i] > 0);
      if(maxPriority <= 1 || maxPriority == INT_MIN) // later by pseudocosts e.g.
      {
         if( maxPriority < 1 || maxPriority == INT_MIN)
         {
            maxPriority = 1; // only choose here first smallest S
            minSizeOfMaxPriority = (*record)->sequencesizes[i];
            Index = i;
         }
         else
            if( (*record)->sequencesizes[i] < minSizeOfMaxPriority )
            {
               minSizeOfMaxPriority = (*record)->sequencesizes[i];
               Index = i;
            }
      }
   }
   assert(maxPriority != INT_MIN);
   assert(minSizeOfMaxPriority != INT_MAX);
   assert(Index >= 0);

   *Ssize = minSizeOfMaxPriority;
   SCIP_CALL( SCIPallocMemoryArray(scip, S, *Ssize) );
   for(i=0; i< *Ssize;++i)
      for(j=0; j<3; ++j)
         (*S)[i][j] = (*record)->record[Index][i][j];

   assert(S!=NULL);
   assert(*S!=NULL);

   //free record
   for( i=0; i< (*record)->recordsize; ++i )
   {
      SCIPfreeMemoryArray(scip, &((*record)->record[i]) );
   }
   SCIPfreeMemoryArray(scip, &((*record)->record) );

   SCIPdebugMessage("with size %d \n", *Ssize);

   assert(*S!=NULL);

   return SCIP_OKAY;
}



// separation at a node other than the root node
static
SCIP_RETCODE Explore( SCIP* scip, ComponentBoundSequence** C, int Csize, int* sequencesizes, int p, struct GCG_Strip** F, int Fsize, int* IndexSet, int IndexSetSize, ComponentBoundSequence** S, int* Ssize, struct GCG_Record** record )
{
   int i;
   int j;
   int k;
   int l;
   SCIP_Real ivalue;
   int isense;
   SCIP_Real median;
   int Fupper;
   int Flower;
   int Cupper;
   int Clower;
   struct GCG_Strip** copyF;
   ComponentBoundSequence* copyS;
   ComponentBoundSequence** CopyC;
   int* newsequencesizes;
   SCIP_Real alpha_i;
   SCIP_Real  muF;
   SCIP_Bool found;
   SCIP_Real mu_F;
   SCIP_Real max;

   i = 0;
   j = 0;
   k = 0;
   l = 0;
   Fupper = 0;
   Flower = 0;
   Cupper = 0;
   Clower = 0;
   CopyC = NULL;
   newsequencesizes = NULL;
   copyF = NULL;
   CopyC = NULL;
   copyS = NULL;
   muF = 0;
   max = 0;
   found = FALSE;

   SCIPdebugMessage("Explore\n");

   SCIPdebugMessage("with Fsize = %d, Csize = %d, Ssize = %d\n", Fsize, Csize, *Ssize);

   //call separate?
   if( C == NULL || Fsize==0 || IndexSetSize==0 || Csize == 0)
   {
      SCIPdebugMessage("go to Separate\n");
      Separate( scip, F, Fsize, IndexSet, IndexSetSize, *S, *Ssize, record );

      if(S != NULL && *Ssize > 0 && *S != NULL)
      {
         SCIPfreeMemoryArray(scip, S);
         S = NULL;
         *Ssize = 0;
      }

      return SCIP_OKAY;
   }

   assert( C!=NULL );
   assert( Csize>0 );
   assert( F != NULL );
   assert( IndexSet != NULL );
   assert( sequencesizes != NULL );

   //find i which is in all S in C on position p
   while( sequencesizes[k] < p )
   {
      SCIPdebugMessage("sequencesizes[%d] = %d\n", k, sequencesizes[k]);
      ++k;
      if( k >= Csize )
      {
         SCIPdebugMessage("no %dth element bounded\n", p);
         Separate( scip, F, Fsize, IndexSet, IndexSetSize, *S, *Ssize, record );

         if(S != NULL && *Ssize > 0 && *S != NULL)
         {
            SCIPfreeMemoryArray(scip, S);
            S = NULL;
            *Ssize = 0;
         }

         return SCIP_OKAY;
      }
      assert( k < Csize );
   }
   i = C[k][p-1][0];
   isense = C[k][p-1][1];
   ivalue = C[k][p-1][2];

   assert(i < F[0]->generatorsize);

   SCIPdebugMessage("i = %d; ivalue = %g\n", i, ivalue);

   for(j=0; j<Fsize; ++j)
   {
      for( l=0; l<IndexSetSize; ++l)
      {
         assert(IndexSet[l] < F[j]->generatorsize);
         if(F[j]->generator[IndexSet[l]] > max )
         {
            max = F[j]->generator[IndexSet[l]];
         }
      }
   }

   if(max == 0)
      max = 1;

   for(j=0; j<Fsize; ++j)
   {
      muF += max * F[j]->mastervarValue;
      //muF += F[j]->generator[i] * F[j]->mastervarValue;
   }

   SCIPdebugMessage("muF = %g\n", muF);

   alpha_i = 0;
   for(j=0; j<Fsize; ++j)
   {
      if(isense == 1)
      {
         if( SCIPisGE(scip, F[j]->generator[i], ivalue) )
         {
            assert( i < F[j]->generatorsize);
            alpha_i += F[j]->generator[i] * F[j]->mastervarValue;
            if(F[j]->mastervarValue != 0 && F[j]->generator[i] != 0)
            {
               //					SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
               //					SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
            }
         }
      }
      else
      {
         if( SCIPisLT(scip, F[j]->generator[i], ivalue) ) //F[j]->generator[i] < ivalue)
         {
            alpha_i += F[j]->generator[i] * F[j]->mastervarValue;
            if(F[j]->mastervarValue != 0 && F[j]->generator[i] != 0)
            {
               SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
               SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
            }
         }

      }
   }
   if( alpha_i == 0 && isense != 1)
   {
      isense = 1;
      for(j=0; j<Fsize; ++j)
      {
         if(isense == 1)
         {
            if( SCIPisGE(scip, F[j]->generator[i], ivalue) )
            {
               alpha_i += F[j]->generator[i] * F[j]->mastervarValue;
               if(F[j]->mastervarValue != 0 && F[j]->generator[i] != 0)
               {
                  //					SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
                  //					SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
               }
            }
         }
         else
         {
            if( SCIPisLT(scip, F[j]->generator[i], ivalue) )
            {
               alpha_i += F[j]->generator[i] * F[j]->mastervarValue;
               if(F[j]->mastervarValue != 0 && F[j]->generator[i] != 0)
               {
                  SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
                  SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
               }
            }

         }
      }
   }
   SCIPdebugMessage("alpha[%d] = %g\n", i, alpha_i);

   if( SCIPisGT(scip, alpha_i - SCIPfloor(scip, alpha_i), 0) )
   {
      found = TRUE;
      SCIPdebugMessage("fractional alpha[%d] = %g\n", i, alpha_i);

      //add to record
      ++(*Ssize);
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyS, *Ssize) );
      for(l=0; l < *Ssize-1; ++l)
         for(j=0; j<3; ++j)
            copyS[l][j] = (*S)[l][j];

      copyS[*Ssize-1][0] = i;
      copyS[*Ssize-1][1] = isense;
      copyS[*Ssize-1][2] = ivalue;

      mu_F = 0;
      for(l=0; l<Fsize; ++l)
      {
         if(isense == 1 )
         {
            if( SCIPisGE(scip, F[l]->generator[i], ivalue) )
               mu_F += F[l]->mastervarValue;
         }
         else
         {
            if( SCIPisLT(scip, F[l]->generator[i], ivalue) )
               mu_F += F[l]->mastervarValue;
         }
      }

      if(SCIPisGT(scip, mu_F - SCIPfloor(scip, mu_F), 0))
      {
         (*record)->recordsize++;
         if((*record)->recordsize == 1 )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
         }
         else
         {
            SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
         }
         (*record)->record[(*record)->recordsize-1] = copyS;
         (*record)->sequencesizes[(*record)->recordsize-1] = *Ssize;
         --(*Ssize);
      }
      else
      {
         found = FALSE;
      }
   }

   if(found)
   {
      SCIPdebugMessage("found fractional alpha\n");
      return SCIP_OKAY;
   }

   ++(*Ssize);
   if(S==NULL || *S==NULL)
   {
      assert(*Ssize == 1);
      SCIP_CALL( SCIPallocMemoryArray(scip, S, *Ssize) );
   }
   else
   {
      assert(*Ssize >1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, S, *Ssize) );
   }
   median = ivalue;
   (*S)[*Ssize-1][0] = i;
   (*S)[*Ssize-1][1] = 1;
   (*S)[*Ssize-1][2] = median;

   for(k=0; k<Fsize; ++k)
   {
      if( SCIPisGE(scip, F[k]->generator[i], median) )
         ++Fupper;
      else
         ++Flower;
   }

   //calculate subset of C
   for(j=0; j< Csize; ++j)
   {
      if(sequencesizes[j] >= p)
      {
         if( SCIPisEQ(scip, C[j][p-1][1], 1) )
            ++Cupper;
         else
            ++Clower;
      }
   }

   if( SCIPisLE(scip, alpha_i, 0) && Fupper != 0 )
      Flower = INT_MAX; //Fupper
   if( SCIPisEQ(scip, alpha_i, muF) && Flower != 0)
      Fupper = INT_MAX; //Flower

   //choose smallest partition
   if( ((Fupper <= Flower && Fupper > 0 ) || Flower <= 0) && Fupper != INT_MAX )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
      j = 0;
      for(k=0; k<Fsize; ++k)
      {
         if( SCIPisGE(scip, F[k]->generator[i], median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Fupper;

      //new C
      k=0;
      if(Cupper > 0)
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Cupper) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Cupper) );
         for(j=0; j< Csize; ++j)
         {
            if(sequencesizes[j] >= p)
               assert(C[j][p-1][0] == i);

            if(sequencesizes[j] >= p &&  SCIPisEQ(scip, C[j][p-1][1], 1) )
            {
               CopyC[k] = C[j];
               newsequencesizes[k] = sequencesizes[j];
               ++k;
            }
         }
      }
      else
         CopyC = NULL;
      Csize = Cupper;
      SCIPdebugMessage("chose upper bound Fupper = %d, Cupper = %d\n", Fupper, Cupper);
   }
   else
   {
      SCIPdebugMessage("chose lower bound Flower = %d Clower = %d\n", Flower, Clower);
      (*S)[*Ssize-1][1] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
      j = 0;
      for(k=0; k<Fsize; ++k)
      {
         if( SCIPisLT(scip, F[k]->generator[i], median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Flower;

      //new C
      k=0;
      if(Clower > 0)
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Clower) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Clower) );
         for(j=0; j< Csize; ++j)
         {
            if(sequencesizes[j] >= p)
               assert(C[j][p-1][0] == i);

            if(sequencesizes[j] >= p && C[j][p-1][1] != 1)
            {
               CopyC[k] = C[j];
               newsequencesizes[k] = sequencesizes[j];
               ++k;
            }
         }
      }
      else
         CopyC = NULL;
      Csize = Clower;
   }
   assert( k <= Csize+1 );

   Explore( scip, CopyC, Csize, newsequencesizes, p+1, copyF, Fsize, IndexSet, IndexSetSize, S, Ssize, record );

   SCIPfreeMemoryArray(scip, &copyF);
   if(CopyC != NULL)
      SCIPfreeMemoryArray(scip, &CopyC);
   if(newsequencesizes != NULL)
      SCIPfreeMemoryArray(scip, &newsequencesizes);

   if(S != NULL && *Ssize > 0 && *S != NULL)
   {
      SCIPfreeMemoryArray(scip, S);
      S = NULL;
      *Ssize = 0;
   }

   return SCIP_OKAY;
}

// callup method for seperate
static
SCIP_RETCODE CallSeparate( SCIP* scip, struct GCG_Strip** F, int Fsize, ComponentBoundSequence** S, int* Ssize, ComponentBoundSequence** C, int Csize, int* CompSizes )
{
   int i;
   int* IndexSet;
   int IndexSetSize;
   GCG_RECORD* record;
   int exploreSsize;
   ComponentBoundSequence* exploreS;

   assert(Fsize > 0);
   assert(F!=NULL);
   exploreSsize = 0;
   exploreS = NULL;
   record = NULL;
   IndexSet = NULL;

   SCIPdebugMessage("Calling Separate\n");


   //record = (struct GCG_Record*) malloc(sizeof(struct GCG_Record));
   SCIP_CALL( SCIPallocBuffer(scip, &record) );
   record->recordsize = 0;
   record->record = NULL;
   record->sequencesizes = NULL;

   //calculate IndexSet
   IndexSetSize = F[0]->generatorsize;
   assert(IndexSetSize > 0);
   SCIP_CALL( SCIPallocMemoryArray(scip, &IndexSet, IndexSetSize) );
   for( i=0; i<IndexSetSize; ++i )
      IndexSet[i] = i;

   //rootnode?
   if( Csize<=0 )
      Separate( scip, F, Fsize, IndexSet, IndexSetSize, NULL, 0, &record );
   else
   {
      assert( C!=NULL );
      Explore( scip, C, Csize, CompSizes, 1, F, Fsize, IndexSet, IndexSetSize, &exploreS, &exploreSsize, &record);
      if(exploreS != NULL)
         SCIPfreeMemoryArray(scip, &exploreS);
   }

   ChoseS( scip, &record, S, Ssize );
   assert(*S!=NULL);

   #ifdef SCIP_DEBUG
   for( i=0; i<*Ssize; ++i)
   {
      SCIPdebugMessage("S[%d][0] = %g\n", i, (*S)[i][0]);
      SCIPdebugMessage("S[%d][1] = %g\n", i, (*S)[i][1]);
      SCIPdebugMessage("S[%d][2] = %g\n", i, (*S)[i][2]);
   }
   #endif
   SCIPfreeMemoryArray(scip, &IndexSet);
   SCIPfreeBuffer(scip, &record);

   return SCIP_OKAY;
}

/** callback deletion method for branching data*/
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteGeneric)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   SCIPdebugMessage("branchDataDeleteGeneric:\n");

   if(*branchdata == NULL)
   {
      SCIPdebugMessage("branchDataDeleteGeneric: cannot delete empty branchdata\n");

      return SCIP_OKAY;
   }

   if((*branchdata)->mastercons != NULL)
   {
      SCIPdebugMessage("branchDataDeleteGeneric: child blocknr %d, %s\n", (*branchdata)->consblocknr,
         SCIPconsGetName((*branchdata)->mastercons) );
   }
   else
   {
      SCIPdebugMessage("branchDataDeleteGeneric: child blocknr %d, empty mastercons\n", (*branchdata)->consblocknr);
   }

   /* release constraint that enforces the branching decision */
   if( (*branchdata)->mastercons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(GCGrelaxGetMasterprob(scip), &(*branchdata)->mastercons) );
      (*branchdata)->mastercons = NULL;
   }

   if((*branchdata)->S != NULL && (*branchdata)->Ssize > 0)
   {
      SCIPfreeMemoryArray(scip, &((*branchdata)->S));
      (*branchdata)->S = NULL;
   }

   if((*branchdata)->consS != NULL && (*branchdata)->consSsize > 0)
   {
      SCIPfreeMemoryArray(scip, &((*branchdata)->consS));
      (*branchdata)->consS = NULL;
   }

   if( (*branchdata)->nchildNodes > 0)
   {
      if((*branchdata)->childlhs != NULL)
      {
         SCIPfreeMemoryArray(scip, &((*branchdata)->childlhs));
         (*branchdata)->childlhs = NULL;
      }
      if((*branchdata)->childbranchdatas!= NULL) // this should not happen
      {
         SCIPfreeMemoryArray(scip, &((*branchdata)->childbranchdatas));
         (*branchdata)->childbranchdatas = NULL;
      }
   }

   SCIPfreeMemory(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

/** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_Bool pruneChildNodeByDominanceGeneric(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lhs,                /**< lhs for childnode which is checkes to be pruned */
   ComponentBoundSequence* childS,           /**< Component Bound Sequence defining the childnode */
   int                   childSsize,
   SCIP_CONS*            masterbranch,
   int                   childBlocknr             /**< number of the block for the childnode */
   )
{
   SCIP_CONS* cons;
   int i;

   i = 0;
   cons = NULL;

   SCIPdebugMessage("Prune by dominance\n");
   cons = GCGconsMasterbranchGetParentcons(masterbranch);

   if( cons == NULL)
   {
      SCIPdebugMessage("cons == NULL, not pruned\n");
      return FALSE;
   }
   while( GCGconsMasterbranchGetParentcons(cons) != NULL )
   {
      GCG_BRANCHDATA* parentdata;

      cons = GCGconsMasterbranchGetParentcons(cons);  //skip the first ancestor node
      parentdata = GCGconsMasterbranchGetBranchdata(cons);

      if(parentdata->blocknr == childBlocknr && parentdata->Ssize >= childSsize && parentdata->S != NULL)
      {
         for( i=0; i<childSsize; ++i)
         {
            if(parentdata->Ssize < childSsize)
               continue;

            if( !SCIPisEQ(scip, parentdata->childlhs[i], lhs) )
               continue;
            if( SCIPisEQ(scip, parentdata->S[i][0], childS[i][0]) && SCIPisEQ(scip, parentdata->S[i][2], childS[i][2]) )
            {
               if( !SCIPisEQ(scip, parentdata->S[i][1], childS[i][1]) && i < childSsize-1 )
                  break; //subset = FALSE;
               if( i == childSsize-1 )
               {
                  if(childSsize == parentdata->Ssize)
                  {
                     SCIPdebugMessage("child pruned\n");
                     return TRUE;
                  }
                  if( !SCIPisEQ(scip, parentdata->S[i][1], childS[i][1]) )  //child is induced by parantdata->S
                  {
                     SCIPdebugMessage("child pruned\n");
                     return TRUE;
                  }
               }
            }
            else
               break; //subset = FALSE; no subset or other lhs
         }
      }
   }

   SCIPdebugMessage("child not pruned\n");
   return FALSE;
}

/** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_RETCODE createChildNodesGeneric(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   ComponentBoundSequence* S,              /**< Component Bound Sequence defining the nodes */
   int                   Ssize,
   int                   blocknr,             /**< number of the block */
   struct GCG_Strip**    F,                   /**< strips with mu>0 */  //for rhs, will be small than
   int                   Fsize,
   SCIP_CONS*            parentcons,
   int*                  nmasternodes,
   SCIP_Bool             createorignodes
   )
{
   SCIP*  masterscip;
   int i;
   int p;
   int k;
   SCIP_Real pL;
   SCIP_Real L;
   SCIP_Real lhs;
   int nmastervars;
   int nmastervars2;
   int ncopymastervars;
   int nbranchcands;
   SCIP_Real mu;  // mu(S)
   SCIP_VAR** mastervars;
   SCIP_VAR** mastervars2;
   SCIP_VAR** branchcands;
   SCIP_VAR** copymastervars;
   GCG_BRANCHDATA* parentdata;

   SCIPdebugMessage("Vanderbeck branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, GCGrelaxGetNIdenticalBlocks(scip, blocknr));

   lhs = 0;
   p = 0;
   k = 0;
   i = 0;
   L = 0;
   pL = GCGrelaxGetNIdenticalBlocks(scip, blocknr);
   mu = 0;
   *nmasternodes = 0;
   parentdata = NULL;

   if(createorignodes)
   {
      parentdata = GCGconsMasterbranchGetBranchdata(parentcons);
      assert( parentdata != NULL);

      for(i=0; i<parentdata->nchildNodes; ++i)
      {
         SCIP_NODE* origchild;
         SCIP_CONS* origcons;
         char childname[SCIP_MAXSTRLEN];

         // define name for constraint
         (void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "child(%d, %g)", p, parentdata->childlhs[i]);
         SCIP_CALL( SCIPcreateChild(scip, &origchild, 0.0, SCIPgetLocalTransEstimate(scip)) );

         //create cons
         SCIP_CALL( GCGcreateConsOrigbranch(scip, &origcons, childname, origchild,
               GCGconsOrigbranchGetActiveCons(scip), branchrule, parentdata->childbranchdatas[i]) );
         SCIP_CALL( SCIPaddConsNode(scip, origchild, origcons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
      }

      return SCIP_OKAY;
   }


   if( parentcons != NULL)
   {
      parentdata = GCGconsMasterbranchGetBranchdata(parentcons);
      if(parentdata == NULL)
      {
         SCIP_CALL( SCIPallocMemory(scip, &parentdata) );
         parentdata->Ssize = 0;
         parentdata->consSsize = 0;
      }
      parentdata->nchildNodes = 0;
      parentdata->S = S;
      parentdata->Ssize = Ssize;
      parentdata->childlhs = NULL;
      parentdata->C = NULL;
      parentdata->sequencesizes = NULL;
      parentdata->blocknr = blocknr;

      SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childlhs), Ssize+1) );
   }
   else
      parentdata = NULL;

   // get variable data of the master problem
   masterscip = GCGrelaxGetMasterprob(scip);
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   nmastervars2 = nmastervars;
   assert(nmastervars >= 0);
   SCIP_CALL( SCIPallocMemoryArray(scip, &copymastervars, nmastervars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &mastervars2, nmastervars) );

   for(i=0; i<nmastervars; ++i)
   {
      mastervars2[i] = mastervars[i];
      copymastervars[i] = mastervars[i];
   }

   assert(scip != NULL);
   assert(Ssize > 0);
   assert(S != NULL);
   assert(F != NULL);
   assert(Fsize > 0);

   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

   SCIPdebugMessage("Vanderbeck branching rule: creating %d nodes\n", Ssize+1);

   for( p=0; p<Ssize+1; ++p )
   {
      GCG_BRANCHDATA* branchchilddata;

      mu = 0;
      branchchilddata = NULL;

      // allocate branchdata for same child and store information
      SCIP_CALL( SCIPallocMemory(scip, &branchchilddata) );
      branchchilddata->consblocknr = blocknr;
      branchchilddata->mastercons = NULL;
      branchchilddata->S = NULL;
      branchchilddata->consS = NULL;
      branchchilddata->C = NULL;
      branchchilddata->sequencesizes = NULL;
      branchchilddata->childlhs = NULL;
      branchchilddata->childbranchdatas = NULL;
      branchchilddata->childnr = p;
      branchchilddata->nchildNodes = 0;
      branchchilddata->Csize = 0;
      branchchilddata->Ssize = 0;
      branchchilddata->consSsize = 0;
      SCIPdebugMessage("p = %d \n", p);
      if( p == Ssize )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(branchchilddata->consS), Ssize) );
         branchchilddata->consSsize = Ssize;
      }
      else
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(branchchilddata->consS), p+1) );
         branchchilddata->consSsize = p+1;
      }
      for( k=0; k<=p; ++k)
      {
         ComponentBoundSequence compBound;

         if( k == Ssize )
         {
            assert( p == Ssize );
            compBound[0] = S[k-1][0];
            compBound[2] = S[k-1][2];
            compBound[1] = S[k-1][1];
            branchchilddata->consS[k-1][0] = compBound[0];
            branchchilddata->consS[k-1][1] = compBound[1];
            branchchilddata->consS[k-1][2] = compBound[2];
         }
         else
         {
            if( k < p )
            {
               compBound[0] = S[k][0];
               compBound[1] = S[k][1];
               compBound[2] = S[k][2];
            }
            else
            {
               compBound[0] = S[k][0];
               compBound[2] = S[k][2];
               if( S[p][1] == 1)
                  compBound[1] = 0;
               else
                  compBound[1] = 1;
            }
            branchchilddata->consS[k][0] = compBound[0];
            branchchilddata->consS[k][1] = compBound[1];
            branchchilddata->consS[k][2] = compBound[2];
         }
      }

      //last node?
      if( p == Ssize )
      {
         lhs = pL;
      }
      else
      {
         //calculate mu
         for(k=0;k<=p;++k)
         {
            ncopymastervars = nmastervars2;
            for( i=0; i<ncopymastervars; ++i)
            {
               SCIP_VAR* swap;
               SCIP_Real generator_i;
               SCIP_Real* generator;
               int generatorsize;
               SCIP_Bool* compisinteger;

               generator = NULL;
               compisinteger = NULL;

               if(GCGvarGetBlock(mastervars2[i]) == blocknr)
               {
                  getGenerators(scip, &generator, &generatorsize, &compisinteger, blocknr, mastervars, nmastervars, mastervars2[i]);
                  generator_i = generator[ (int) SCIPceil(scip, S[k][0]-0.5)];

                  if(S[k][1] == 1 )
                  {
                     if( SCIPisGE(scip, generator_i, S[k][2]) )
                     {
                        if( k == p )
                           mu += SCIPgetSolVal(masterscip, NULL, mastervars2[i]);
                     }
                     else
                     {
                        if(ncopymastervars > 0)
                        {
                           swap = mastervars2[i];
                           mastervars2[i] = mastervars2[ncopymastervars-1];
                           mastervars2[ncopymastervars-1] = swap;
                           --ncopymastervars;
                           --i;
                        }
                     }
                  }
                  else  //nested erasing
                  {
                     if( SCIPisLT(scip, generator_i, S[k][2]) )
                     {
                        if( k == p )
                           mu += SCIPgetSolVal(masterscip, NULL, mastervars2[i]);
                     }
                     else
                     {
                        if(ncopymastervars > 0)
                        {
                           swap = mastervars2[i];
                           mastervars2[i] = mastervars2[ncopymastervars-1];
                           mastervars2[ncopymastervars-1] = swap;
                           --ncopymastervars;
                           --i;
                        }
                     }
                  }
               }
               else
               {
                  if(ncopymastervars > 0)
                  {
                     swap = mastervars2[i];
                     mastervars2[i] = mastervars2[ncopymastervars-1];
                     mastervars2[ncopymastervars-1] = swap;
                     --ncopymastervars;
                     --i;
                  }
               }
            }
            nmastervars2 = ncopymastervars;
         }
         if( p == Ssize-1)
            L = SCIPceil(scip, mu);
         else
         {
            L = mu;
         }
         lhs = pL-L+1;
         SCIPdebugMessage("pL-L+1 = %g \n", pL-L+1);
         parentdata->childlhs[p] = lhs;
      }
      pL = L;

      branchchilddata->lhs = lhs;
      branchchilddata->childlhs = NULL;
      branchchilddata->C = NULL;
      branchchilddata->sequencesizes = NULL;
      SCIPdebugMessage("L = %g \n", L);
      SCIPdebugMessage("lhs set to %g \n", lhs);

      if( parentcons == NULL || !pruneChildNodeByDominanceGeneric(scip, lhs, branchchilddata->consS, branchchilddata->consSsize, parentcons, blocknr) )
      {
         ++(*nmasternodes);
         if( parentcons != NULL)
         {
            ++(parentdata->nchildNodes);
            if(parentdata->nchildNodes == 1)
            {
               //	SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childlhs), parentdata->nchildNodes) );
               SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childbranchdatas), parentdata->nchildNodes) );
            }
            else
            {
               //					SCIP_CALL( SCIPreallocMemoryArray(scip, &(parentdata->childlhs), parentdata->nchildNodes) );
               SCIP_CALL( SCIPreallocMemoryArray(scip, &(parentdata->childbranchdatas), parentdata->nchildNodes) );

            }
            SCIP_CALL( SCIPallocMemory(scip, &(parentdata->childbranchdatas[parentdata->nchildNodes-1])) );
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->S = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->C = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->sequencesizes = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->childbranchdatas = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->consS = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->mastercons = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->childlhs = NULL;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->consblocknr = branchchilddata->consblocknr;
            parentdata->childlhs[p] = lhs;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->lhs = lhs;
            parentdata->childbranchdatas[parentdata->nchildNodes-1]->consSsize = branchchilddata->consSsize;
            SCIPdebugMessage("branchchilddata->consSsize = %d\n", branchchilddata->consSsize);
            SCIPdebugMessage("branchchilddata->nchildnodes-1 = %d\n", parentdata->nchildNodes-1);
            SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childbranchdatas[parentdata->nchildNodes-1]->consS), branchchilddata->consSsize) );

            for(i=0; i<branchchilddata->consSsize; ++i)
            {
               SCIPdebugMessage("setting consS[%d][0] = %g\n", i, branchchilddata->consS[i][0]);
               parentdata->childbranchdatas[parentdata->nchildNodes-1]->consS[i][0] = branchchilddata->consS[i][0];
               SCIPdebugMessage("setting consS[%d][1] = %g\n", i, branchchilddata->consS[i][1]);
               parentdata->childbranchdatas[parentdata->nchildNodes-1]->consS[i][1] = branchchilddata->consS[i][1];
               SCIPdebugMessage("setting consS[%d][2] = %g\n", i, branchchilddata->consS[i][2]);
               parentdata->childbranchdatas[parentdata->nchildNodes-1]->consS[i][2] = branchchilddata->consS[i][2];
            }
         }
      }
      else
         parentdata->childlhs[p] = 0;
      SCIPfreeMemory(scip, &branchchilddata);
   }

   SCIPfreeMemoryArray(scip, &mastervars2);
   SCIPfreeMemoryArray(scip, &copymastervars);

   return SCIP_OKAY;
}

/** returns the number of successor nodes needed for branch_master while using the generic branching scheme */
/** @return SCIP_RETCODE, int as input, createorignodes auslagern */
int GCGbranchGenericGetNChildnodes(
   SCIP*            masterscip,
   SCIP_Bool        createorignodes
   )
{
   int nmasternodes;
   SCIP* scip;
   SCIP_Bool feasible;
   SCIP_Bool SinC;
   SCIP_VAR** branchcands;
   SCIP_VAR** allorigvars;
   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_CONS* masterbranchcons;
   SCIP_CONS* parentcons;
   int nbranchcands;
   GCG_BRANCHDATA* branchdata;
   SCIP_VAR* mastervar;
   SCIP_Real mastervarValue;
   ComponentBoundSequence* S;
   ComponentBoundSequence** C;
   struct GCG_Strip** F;
   int Ssize;
   int Csize;
   int Fsize;
   int* sequencesizes;
   int blocknr;
   int i;
   int c;
   int k;
   int allnorigvars;

   blocknr = -2;
   Ssize = 0;
   Csize = 0;
   Fsize = 0;
   i = 0;
   c = 0;
   feasible = TRUE;
   SinC = TRUE;
   S = NULL;
   F = NULL;
   nmasternodes = 0;
   branchdata = NULL;
   S = NULL;
   C = NULL;
   F = NULL;
   sequencesizes = NULL;

   assert(masterscip != NULL);

   SCIPdebugMessage("Get number of childnodes for Vanderbecks generic branching\n");


   /** @todo: Das ist pervers © Christian */
   /** die Methode sollte nur im Original aufgerufen werden */
   if(createorignodes)
   {
      scip = masterscip;
      masterscip = GCGrelaxGetMasterprob(scip);
   }
   else
   {
      scip = GCGpricerGetOrigprob(masterscip);
   }
   assert(scip != NULL);
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

   SCIP_CALL( SCIPgetVarsData(scip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   assert(nbranchcands > 0);

   for( i=0; i<nbranchcands; ++i )
   {
      mastervar = branchcands[i];
      assert(GCGvarIsMaster(mastervar));
      blocknr = GCGvarGetBlock(mastervar);
      if(blocknr >= -1)
         break;
   }
   if( blocknr < -1 )
   {
      feasible = TRUE;
      SCIPdebugMessage("Vanderbeck generic branching rule could not find variables to branch on!\n");

      return -1;
   }
   else
      feasible = FALSE;

   k = 0;
   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);

   //calculate F and the strips
   for( i=0; i<nbranchcands; ++i )
   {
      mastervar = branchcands[i];
      assert(GCGvarIsMaster(mastervar));

      if(blocknr == GCGvarGetBlock(mastervar))
      {
         mastervarValue = SCIPgetSolVal(masterscip, NULL, mastervar);
         if( SCIPisGT(scip, mastervarValue - SCIPfloor(scip, mastervarValue), 0) )
         {
            if(Fsize == 0)
            {
               ++Fsize;
               SCIP_CALL( SCIPallocMemoryArray(scip, &F, Fsize) );
            }
            else
            {
               ++Fsize;
               SCIP_CALL( SCIPreallocMemoryArray(scip, &F, Fsize) );
            }
            SCIP_CALL( SCIPallocBuffer(scip, &(F[Fsize-1]) ));
            F[Fsize-1]->blocknr = blocknr;
            F[Fsize-1]->mastervar = mastervar;
            F[Fsize-1]->mastervarValue= mastervarValue;

            getGenerators(scip, &(F[Fsize-1]->generator), &(F[Fsize-1]->generatorsize), &(F[Fsize-1]->compisinteger), blocknr, mastervars, nmastervars, mastervar);

            ++k;
         }
      }
   }

   //old data to regard?
   if( masterbranchcons != NULL )
   {
      //calculate C
      Csize = 1;
      SCIP_CALL( SCIPallocMemoryArray(scip, &C, Csize) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &sequencesizes, Csize) );

      parentcons = masterbranchcons;

      if( parentcons != NULL && GCGconsMasterbranchGetBranchdata(parentcons)->consS != NULL && GCGconsMasterbranchGetBranchdata(parentcons)->consblocknr == blocknr)
      {
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
         assert(branchdata != NULL);
         assert(branchdata->consSsize > 0);
         C[0] = NULL;
         SCIP_CALL( SCIPallocMemoryArray(scip, &(C[0]), branchdata->consSsize) );
         for(i=0; i<branchdata->consSsize; ++i)
         {
            C[0][i][0] = branchdata->consS[i][0];
            C[0][i][1] = branchdata->consS[i][1];
            C[0][i][2] = branchdata->consS[i][2];
         }
         sequencesizes[0] = branchdata->consSsize;

         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
      }
      else
         Csize = 0;
      while( parentcons != NULL )
      {
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
         assert(branchdata != NULL);
         if( branchdata->consS == NULL || branchdata->consSsize == 0)
            break;
         //S not yet in C ?
         SinC = FALSE;
         for( c=0; c<Csize && !SinC; ++c)
         {
            SinC = TRUE;
            if(branchdata->consSsize == sequencesizes[c])
            {
               for( i=0; i<branchdata->consSsize; ++i)
               {
                  if( !SCIPisEQ(scip, branchdata->consS[i][0], C[c][i][0]) || !SCIPisEQ(scip, branchdata->consS[i][1], C[c][i][1]) || !SCIPisEQ(scip, branchdata->consS[i][2], C[c][i][2]) )
                  {
                     SinC = FALSE;
                     break;
                  }
               }
            }
            else
               SinC = FALSE;
         }
         if(!SinC)
         {
            ++Csize;
            SCIP_CALL( SCIPreallocMemoryArray(scip, &C, Csize) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &sequencesizes, Csize) );
            C[Csize-1] = NULL;
            SCIP_CALL( SCIPallocMemoryArray(scip, &(C[Csize-1]), branchdata->consSsize) );

            /** @todo copy memory */
            for(i=0; i<branchdata->consSsize; ++i)
            {
               C[Csize-1][i][0] = branchdata->consS[i][0];
               C[Csize-1][i][1] = branchdata->consS[i][1];
               C[Csize-1][i][2] = branchdata->consS[i][2];
            }
            sequencesizes[Csize-1] = branchdata->consSsize;
         }
         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
      }

      if(C != NULL)
      {
         SCIPdebugMessage("Csize = %d\n", Csize);

         for(i=0; i<Csize; ++i)
         {
            for(c=0; c<sequencesizes[i]; ++c)
            {
               SCIPdebugMessage("C[%d][%d][0] = %g\n", i, c, C[i][c][0]);
               SCIPdebugMessage("C[%d][%d][1] = %g\n", i, c, C[i][c][1]);
               SCIPdebugMessage("C[%d][%d][2] = %g\n", i, c, C[i][c][2]);
            }
         }

         //		SCIP_CALL( InducedLexicographicSort(scip, F, Fsize, C, Csize, sequencesizes) );
         SCIP_CALL( CallSeparate(scip, F, Fsize, &S, &Ssize, C, Csize, sequencesizes) );
      }
      else
      {
         SCIPdebugMessage("C == NULL\n");
         //		SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
         SCIP_CALL( CallSeparate( scip, F, Fsize, &S, &Ssize, NULL, 0, NULL ) );
      }
      SCIPfreeMemoryArray(scip, &sequencesizes);
      for( i=0; i<Csize; ++i)
         SCIPfreeMemoryArray(scip, &(C[i]));
      SCIPfreeMemoryArray(scip, &C);
   }
   else
   {
      SCIPdebugMessage("masterbanchcons == NULL\n");
      //	SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
      SCIP_CALL( CallSeparate( scip, F, Fsize, &S, &Ssize, NULL, 0, NULL ) );
   }
   assert(S!=NULL);

   if( feasible )
   {
      SCIPdebugMessage("Vanderbeck generic branching rule could not find variables to branch on!\n");
      return -1;
   }

   /* create the |S|+1 child nodes in the branch-and-bound tree */
   SCIP_CALL( createChildNodesGeneric(scip, NULL, S, Ssize, blocknr, F, Fsize, masterbranchcons, &nmasternodes, createorignodes) );

   SCIPdebugMessage("free F\n");
   for(i=0; i<Fsize; ++i)
   {
      SCIPfreeMemoryArray(scip, &(F[i]->compisinteger));
      SCIPfreeMemoryArray(scip, &(F[i]->generator));
      SCIPfreeBuffer(scip, &(F[i]));
   }
   SCIPfreeMemoryArray(scip, &F);

   return nmasternodes;
}

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterGeneric)
{
   SCIP* masterscip;
   SCIP_VAR** mastervars;
   SCIP_VAR** copymastervars;
   SCIP_VAR** allorigvars;
   int allnorigvars;
   int nmastervars;
   int nvarsadded;
   int nnewmastervars;
   int oldnmastervars;
   int i;
   int p;
   char name[SCIP_MAXSTRLEN];

   i = 0;
   p = 0;
   nmastervars = 0;
   nvarsadded = 0;

   assert(scip != NULL);
   assert(branchdata != NULL);

   masterscip = scip;//GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &copymastervars, nmastervars) );

   for(i=0; i<nmastervars; ++i)
      copymastervars[i] = mastervars[i];

   oldnmastervars = nmastervars;

   SCIPdebugMessage("branchActiveMasterGeneric: Block %d, Ssize %d)\n", branchdata->consblocknr,
      branchdata->consSsize);

   if(branchdata->consS == NULL)
   {
      assert(branchdata->consSsize == 0);
      SCIPdebugMessage("root node:\n");
      return SCIP_OKAY;
   }

   #ifdef SCIP_DEBUG
   for( i=0; i<branchdata->consSsize; ++i)
   {
      SCIPdebugMessage("S[%d][0] = %g\n", i, branchdata->consS[i][0]);
      SCIPdebugMessage("S[%d][1] = %g\n", i, branchdata->consS[i][1]);
      SCIPdebugMessage("S[%d][2] = %g\n", i, branchdata->consS[i][2]);
   }
   #endif

   /* create corresponding constraint in the master problem, if not yet created */
   if( branchdata->mastercons == NULL && branchdata->consSsize > 0)
   {

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "child(%d, %g)", branchdata->consSsize, branchdata->lhs);

      // create constraint for child
      SCIP_CALL( SCIPcreateConsLinear(masterscip, &(branchdata->mastercons), name, 0, NULL, NULL,
            branchdata->lhs, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );

      //add mastervars
      for(p=0; p< branchdata->consSsize; ++p)
      {
         nnewmastervars = nmastervars;
         for(i=0; i<nnewmastervars; ++i)
         {
            SCIP_Real* generator;
            SCIP_Bool* compisinteger;
            int generatorsize;
            SCIP_Real generator_i;

            generatorsize = 0;
            generator = NULL;
            compisinteger = NULL;
            if( GCGvarGetBlock(copymastervars[i]) == branchdata->consblocknr )
            {
               getGenerators(scip, &generator, &generatorsize, &compisinteger, branchdata->consblocknr, mastervars, oldnmastervars, copymastervars[i]);

               assert((int) SCIPceil(scip, branchdata->consS[p][0]-0.5) < generatorsize);
               generator_i = generator[(int) SCIPceil(scip, branchdata->consS[p][0]-0.5)];

               if( SCIPisEQ(scip, branchdata->consS[p][1], 1) )
               {
                  if( SCIPisGE(scip, generator_i, branchdata->consS[p][2]) )
                  {
                     if( p == branchdata->consSsize-1 )
                     {
                        //add var to constraint
                        ++nvarsadded;
                        SCIP_CALL( SCIPaddCoefLinear(masterscip, branchdata->mastercons, copymastervars[i], 1.0) );
                     }
                  }
                  else
                  {
                     //small down array
                     copymastervars[i] = copymastervars[nnewmastervars-1];
                     --i;
                     --nnewmastervars;
                  }
               }
               else
               {
                  if( SCIPisLT(scip, generator_i, branchdata->consS[p][2]) )
                  {
                     if( p == branchdata->consSsize-1 )
                     {
                        //add var to constraint
                        ++nvarsadded;
                        SCIP_CALL( SCIPaddCoefLinear(masterscip, branchdata->mastercons, copymastervars[i], 1.0) );
                     }
                  }
                  else
                  {
                     //small down array
                     copymastervars[i] = copymastervars[nnewmastervars-1];
                     --i;
                     --nnewmastervars;
                  }
               }
            }
            if(generatorsize > 0)
            {
               SCIPfreeMemoryArray(scip, &generator);
               generator = NULL;
               SCIPfreeMemoryArray(scip, &compisinteger);
               compisinteger = NULL;
            }
            else
            {
               //small down array
               copymastervars[i] = copymastervars[nnewmastervars-1];
               --i;
               --nnewmastervars;
            }
         }
         nmastervars = nnewmastervars;
      }
   }
   /* add constraint to the master problem that enforces the branching decision */
   SCIP_CALL( SCIPaddCons(masterscip, branchdata->mastercons) );

   SCIPdebugMessage("%d vars added with lhs= %g\n", nvarsadded, branchdata->lhs);

   SCIPfreeMemoryArray(scip, &copymastervars);

   return SCIP_OKAY;
}

/** callback deactivation method */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterGeneric)
{
   SCIP* masterscip;

   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   masterscip = scip;//GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);

   SCIPdebugMessage("branchDeactiveMasterGeneric: Block %d, Ssize %d)\n", branchdata->consblocknr,
      branchdata->consSsize);

   /* remove constraint from the master problem that enforces the branching decision */
   assert(branchdata->mastercons != NULL);
   SCIP_CALL( SCIPdelCons(masterscip, branchdata->mastercons) );

   SCIP_CALL( SCIPreleaseCons(masterscip, &(branchdata->mastercons)));
   branchdata->mastercons = NULL;

   return SCIP_OKAY;
}



/** callback propagation method */
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterGeneric)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);
   assert(branchdata->consS != NULL);

   SCIPdebugMessage("branchPropMasterGeneric: Block %d ,Ssize %d)\n", branchdata->consblocknr,
      branchdata->consSsize);

   return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpGeneric)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execlp method of generic branching\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}



/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextGeneric)
{  /*lint --e{715}*/
   SCIP* masterscip;
   SCIP_Bool feasible;
   SCIP_Bool discretization;
   SCIP_CONS* masterbranchcons;
   GCG_BRANCHDATA* branchdata;
   int i;
   int j;
   i = 0;
   j = 0;
   feasible = TRUE;
   branchdata = NULL;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of Vanderbecks generic branching\n");

   *result = SCIP_DIDNOTRUN;

   /* the branching scheme only works for the discretization approach */
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
   {
      SCIPdebugMessage("Generic branching only for discretization approach\n");
      return SCIP_OKAY;
   }

   /* do not perform Ryan & Foster branching if we have neither a set partitioning nor a set covering structure */
   if( GCGrelaxIsMasterSetCovering(scip) || GCGrelaxIsMasterSetPartitioning(scip) )
   {
      SCIPdebugMessage("Generic branching executed on a set covering or set partitioning problem\n");
   }

   /* check whether the current original solution is integral */
   #ifdef SCIP_DEBUG
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
   #else
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), FALSE, TRUE, TRUE, TRUE, &feasible) );
   #endif
   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));

      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   masterscip = GCGrelaxGetMasterprob(scip);

   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   if(masterbranchcons != NULL)
      branchdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);

   if(branchdata!=NULL)
   {
      SCIPdebugMessage("branchdata->nchildNodes = %d\n", branchdata->nchildNodes);
      for(i=0; i<branchdata->nchildNodes; ++i)
      {
         SCIP_CONS* origcons;
         SCIP_NODE* origchild;
         GCG_BRANCHDATA* nodebranchdata;
         char childname[SCIP_MAXSTRLEN];

         //set branchdata
         SCIP_CALL( SCIPallocMemory(scip, &nodebranchdata) );
         nodebranchdata->lhs = branchdata->childbranchdatas[i]->lhs;
         nodebranchdata->nchildNodes = 0;
         nodebranchdata->consblocknr = branchdata->childbranchdatas[i]->consblocknr;
         nodebranchdata->mastercons = branchdata->childbranchdatas[i]->mastercons;
         nodebranchdata->S = NULL;
         nodebranchdata->consS = NULL;
         nodebranchdata->childbranchdatas = NULL;
         nodebranchdata->childlhs = NULL;
         nodebranchdata->C = NULL;
         nodebranchdata->sequencesizes = NULL;
         nodebranchdata->childnr = branchdata->childbranchdatas[i]->childnr;
         nodebranchdata->consSsize =	branchdata->childbranchdatas[i]->consSsize;
         SCIPdebugMessage("consSsize = %d, lhs = %g\n", nodebranchdata->consSsize, nodebranchdata->lhs);
         SCIP_CALL( SCIPallocMemoryArray(scip, &(nodebranchdata->consS), nodebranchdata->consSsize) );

         for(j=0; j< nodebranchdata->consSsize; ++j)
         {
            nodebranchdata->consS[j][0] = branchdata->childbranchdatas[i]->consS[j][0];
            nodebranchdata->consS[j][1] = branchdata->childbranchdatas[i]->consS[j][1];
            nodebranchdata->consS[j][2] = branchdata->childbranchdatas[i]->consS[j][2];
            SCIPdebugMessage("consS[%d][0] = %g\n", j, nodebranchdata->consS[j][0]);
            SCIPdebugMessage("consS[%d][1] = %g\n", j, nodebranchdata->consS[j][1]);
            SCIPdebugMessage("consS[%d][2] = %g\n", j, nodebranchdata->consS[j][2]);
         }

         // define name for constraint
         (void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "child(%d, %g)", i, nodebranchdata->lhs);
         SCIP_CALL( SCIPcreateChild(scip, &origchild, 0.0, SCIPgetLocalTransEstimate(scip)) );

         //create cons
         SCIP_CALL( GCGcreateConsOrigbranch(scip, &origcons, childname, origchild,
               GCGconsOrigbranchGetActiveCons(scip), branchrule, nodebranchdata) );
         SCIP_CALL( SCIPaddConsNode(scip, origchild, origcons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
      }
      if(branchdata->childbranchdatas != NULL)
      {
         for(i=0; i< branchdata->nchildNodes; ++i)
            SCIPfreeMemory(scip, &(branchdata->childbranchdatas[i]));
         SCIPfreeMemoryArray(scip, &(branchdata->childbranchdatas));
         branchdata->childbranchdatas = NULL;
      }
   }
   else
   {
      SCIPdebugMessage("Branchdata is NULL!\n");

      GCGbranchGenericGetNChildnodes(scip, TRUE);
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */ /*todo*/
static
SCIP_DECL_BRANCHEXECPS(branchExecpsGeneric)
{  /*lint --e{715}*/
   SCIP* masterscip;
   SCIP_CONS* masterbranchcons;
   GCG_BRANCHDATA* branchdata;
   int i;

   i = 0;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of Vanderbecks generic branching\n");

   return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   masterscip = GCGrelaxGetMasterprob(scip);

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of Vanderbecks generic branching\n");

   *result = SCIP_DIDNOTRUN;

   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   if(masterbranchcons!=NULL)
      branchdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);

   if(branchdata!=NULL)
   {
      for(i=0; i<branchdata->nchildNodes; ++i)
      {
         SCIP_CONS* origcons;
         SCIP_NODE* origchild;
         char childname[SCIP_MAXSTRLEN];
         // define name for constraint
         (void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "child(%d, %g)", i, branchdata->lhs);
         SCIP_CALL( SCIPcreateChild(scip, &origchild, 0.0, SCIPgetLocalTransEstimate(scip)) );

         //create cons
         SCIP_CALL( GCGcreateConsOrigbranch(scip, &origcons, childname, origchild,
               GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdata->childbranchdatas[i]) );
         SCIP_CALL( SCIPaddConsNode(scip, origchild, origcons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
      }
      SCIPfreeMemoryArray(scip, &(branchdata->childbranchdatas));
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitGeneric)
{
   assert(branchrule != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterGeneric,
         branchDeactiveMasterGeneric, branchPropMasterGeneric, NULL, branchDataDeleteGeneric) );

   return SCIP_OKAY;
}

/** creates the most infeasible LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleGeneric(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create branching rule data */
   branchruledata = NULL;

   SCIPdebugMessage("Include method of Vanderbecks generic branching\n");

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchCopyGeneric,
         branchFreeGeneric, branchInitGeneric, branchExitGeneric, branchInitsolGeneric,
         branchExitsolGeneric, branchExeclpGeneric, branchExecextGeneric, branchExecpsGeneric,
         branchruledata) );

   /* include event handler for adding generated mastervars to the branching constraints */
   SCIP_CALL( SCIPincludeEventHdlrGenericbranchvaradd(GCGrelaxGetMasterprob(scip)) );

   return SCIP_OKAY;
}
