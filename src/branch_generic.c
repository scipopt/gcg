/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_generic.c
 * @ingroup BRANCHINGRULES
 * @brief  branching rule based on vanderbeck's generic branching scheme
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG

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

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"

#include <stdio.h>
#include <stdlib.h>

#define BRANCHRULE_NAME          "generic"
#define BRANCHRULE_DESC          "generic branching rule by Vanderbeck"
#define BRANCHRULE_PRIORITY      99999
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

struct GCG_BranchData
{
   GCG_COMPSEQUENCE**   C;                  /**< S[k] bound sequence for block k */ //!!! sort of each C[i] = S is important !!!
   int*                 sequencesizes;      /**< number of bounds in S[k] */
   int                  Csize;
   SCIP_Real            lhs;
   SCIP_CONS*           mastercons;         /**< constraint enforcing the branching restriction in the master problem */
   GCG_COMPSEQUENCE*    consS;              /**< component bound sequence which induce the current branching constraint */
   int                  consSsize;
   int                  consblocknr;
};

/** set of component bounds in separate */
struct GCG_Record
{
   GCG_COMPSEQUENCE**   record;             /**< returnvalue of separate function */
   int                  recordsize;
   int*                 sequencesizes;
};
typedef struct GCG_Record GCG_RECORD;

/** an abstract strip, only for compare in ILO needed */
struct GCG_Strip
{
   SCIP*                scip;
   SCIP_VAR*            mastervar;
   GCG_COMPSEQUENCE**   C;             /**< current set of comp bound sequences */
   int                  Csize;
   int*                 sequencesizes;
};
typedef struct GCG_Strip GCG_STRIP;

/*
 * Callback methods
 */

/* define not used callback as NULL*/
#define branchFreeGeneric NULL
#define branchExitGeneric NULL
#define branchInitsolGeneric NULL
#define branchExitsolGeneric NULL

/*
 * branching specific interface methods
 */

/** computes the generator of mastervar for the entry in origvar
 * @return entry of the generator corresponding to origvar */
SCIP_Real getGeneratorEntry(
   SCIP_VAR*            mastervar,          /**< current mastervariable */
   SCIP_VAR*            origvar             /**< corresponding origvar */
   )
{
   int i;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   i = 0;
   assert(mastervar != NULL);
   assert(origvar != NULL);

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   for( i=0; i < norigvars; ++i)
   {
      if(origvars[i] == origvar )
      {
         //SCIPdebugMessage("origvar found \n");
         return origvals[i];
      }
   }

   return 0;
}

/** method for calculating the maximum over all generatorentries in F
 * @return maxentry */
static
SCIP_Real getMaxGeneratorEntry(
   SCIP*                scip,           /**< SCIP data structure */
   SCIP_VAR**           F,              /**< array of mastervars */
   int                  Fsize,          /**< number of mastervars */
   SCIP_VAR**           IndexSet,       /**< set of origvars to respect*/
   int                  IndexSetSize    /**< number of origvars to respect */
   )
{
   int i;
   int j;
   SCIP_Real maxentry;

   i = 0;
   j = 0;
   maxentry = 0;

   assert(F != NULL);
   assert(Fsize > 0);
   assert(IndexSet != NULL);
   assert(IndexSetSize > 0);

   for( i=0; i<Fsize; ++i)
   {
      for( j=0; j<IndexSetSize; ++j)
      {
         SCIP_Real generatorentry;

         generatorentry = getGeneratorEntry(F[i], IndexSet[j]);
         maxentry = MAX(generatorentry, maxentry);
      }
   }

   return maxentry;
}

/** method for initializing the set of respected indices */
static
SCIP_RETCODE InitIndexSet(
   SCIP*                scip,           /**< SCIP data structure */
   SCIP_VAR**           F,              /**< array of fractional mastervars */
   int                  Fsize,          /**< number of fractional mastervars */
   SCIP_VAR***           IndexSet,       /**< set to initialize */
   int*                 IndexSetSize   /**< size of the index set */
   )
{
   int i;
   int j;
   int k;
   SCIP_VAR** origvars;
   int norigvars;

   i = 0;
   j = 0;
   k = 0;
   norigvars = 0;
   *IndexSet = NULL;
   *IndexSetSize = 0;
   assert( F!= NULL);
   assert( Fsize > 0);

   for( i=0; i<Fsize; ++i )
   {
      origvars = GCGmasterVarGetOrigvars(F[i]);
      norigvars = GCGmasterVarGetNOrigvars(F[i]);

      if( *IndexSetSize == 0 && norigvars > 0 )
      {
         *IndexSetSize = norigvars;
         SCIP_CALL( SCIPallocMemoryArray(scip, IndexSet, *IndexSetSize) );
         for( j=0; j<*IndexSetSize; ++j )
         {
            (*IndexSet)[j] = origvars[j];
         }
      }
      else
      {
         for( j=0; j<norigvars; ++j )
         {
            int oldsize;

            oldsize = *IndexSetSize;

            for( k=0; k<oldsize; ++k )
            {
               // if variable already in union
               if( (*IndexSet)[k] == origvars[j] )
               {
                  break;
               }
               if( k == oldsize-1)
               {
                  // add variable to the end
                  ++(*IndexSetSize);
                  SCIP_CALL( SCIPreallocMemoryArray(scip, IndexSet, *IndexSetSize) );
                  (*IndexSet)[*IndexSetSize-1] = origvars[j];
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** method for calculating the median over all fractional components values using
 * the quickselect algorithm (or a variant of it)
 *
 * This method will change the array
 *
 * @return median or if the median is the mimimum return ceil(arithm middle)
 */
static
SCIP_Real GetMedian(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_Real*           array,              /**< array to find the median in (will be destroyed) */
   int                  arraysize,          /**< size of the array */
   SCIP_Real            min                 /**< minimum of array */
   )
{
   SCIP_Real Median;
   SCIP_Real swap;
   int l;
   int r;
   int i;
   int j;
   int MedianIndex;
   SCIP_Real arithmMiddle;

   assert(scip != NULL);
   assert(array != NULL);
   assert(arraysize > 0);

   r = arraysize -1;
   l = 0;
   arithmMiddle = 0;

   if( arraysize & 1 )
      MedianIndex = arraysize/2;
   else
      MedianIndex = arraysize/2 -1;

   while( l < r-1 )
   {
      Median = array[MedianIndex];
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
   Median = array[MedianIndex];

   if( SCIPisEQ(scip, Median, min) )
   {
//      SCIPdebugMessage("median = ceil(arithmMiddle), min = %f \n", min);

      for( i=0; i<arraysize; ++i )
         arithmMiddle += 1.0*array[i]/arraysize;

      Median = SCIPceil(scip, arithmMiddle);
   }
   //SCIPdebugMessage("median = %f \n", Median);

   return Median;
}

/** comparefunction for lexicographical sort */
static
SCIP_DECL_SORTPTRCOMP(ptrcomp)
{
   GCG_STRIP* strip1;
   GCG_STRIP* strip2;
   SCIP_VAR* mastervar1;
   SCIP_VAR* mastervar2;
   SCIP_VAR** origvars;
   int norigvars;
   int i;

   strip1 = (GCG_STRIP*) elem1;
   strip2 = (GCG_STRIP*) elem2;

   mastervar1 = strip1->mastervar;
   mastervar2 = strip2->mastervar;

   i = 0;
   assert(mastervar1 != NULL);
   assert(mastervar2 != NULL);

   if( GCGvarGetBlock(mastervar1) == -1 )
   {
      SCIPdebugMessage("linkingvar\n");
      assert(GCGvarIsLinking(mastervar1));
   }
   if( GCGvarGetBlock(mastervar2) == -1 )
   {
      SCIPdebugMessage("linkingvar\n");
      assert(GCGvarIsLinking(mastervar2));
   }

   origvars = GCGmasterVarGetOrigvars(mastervar1);
   norigvars = GCGmasterVarGetNOrigvars(mastervar1);

   //assert(strip1->generatorsize == strip2->generatorsize);

   for( i=0; i<norigvars; ++i )
   {
      if( getGeneratorEntry(mastervar1, origvars[i]) > getGeneratorEntry(mastervar2, origvars[i]) )
         return -1;
      if( getGeneratorEntry(mastervar1, origvars[i]) < getGeneratorEntry(mastervar2, origvars[i]) )
         return 1;
   }

   return 0;
}

/** lexicographical sort using scipsort
 * This method will change the array
 */
static
SCIP_RETCODE LexicographicSort(
   GCG_STRIP**           array,              /**< array to sort (will be changed) */
   int                  arraysize           /**< size of the array */
   )
{

   assert(array != NULL);
   assert(arraysize > 0);

   SCIPdebugMessage("Lexicographic sorting\n");

   SCIPsortPtr((void**)array, ptrcomp, arraysize );

   return SCIP_OKAY;
}


/** compare function for ILO: returns 1 if bd1 < bd2 else -1 with respect to bound sequence */
static
int ILOcomp(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_VAR*            mastervar1,         /**< first strip */
   SCIP_VAR*            mastervar2,         /**< second strip */
   GCG_COMPSEQUENCE**   C,                  /**< component bound sequence to compare with */
   int                  NBoundsequences,    /**< size of the bound sequence */
   int*                 sequencesizes,      /**< sizes of the bound sequences */
   int                  p                   /**< current depth in C*/
   )
{
   int ivalue;
   int j;
   int k;
   int l;
   int Nupper;
   int Nlower;
   SCIP_Bool returnvalue;
   GCG_COMPSEQUENCE** CopyC;
   SCIP_VAR* origvar;
   int* newsequencesizes;

   j = 0;
   k = 0;
   l = 0;
   Nupper = 0;
   Nlower = 0;
   origvar = NULL;
   newsequencesizes = NULL;

   /* lexicographic Order? */
   if( C == NULL || NBoundsequences <= 1 )
      return (*ptrcomp)( mastervar1, mastervar2);

   assert(C != NULL);
   assert(NBoundsequences > 0);

   /* find i which is in all S in C on position p */
   while( sequencesizes[k] < p )
   {
      ++k;
      assert(k < NBoundsequences);
   }
   origvar = C[k][p-1].component;
   assert(origvar != NULL);
   ivalue = C[k][p-1].bound;

   /* calculate subset of C */
   for( j=0; j< NBoundsequences; ++j )
   {
      if( sequencesizes[j] >= p )
      {
         assert(C[j][p-1].component == origvar);
         if( C[j][p-1].sense == GCG_COMPSENSE_GE )
            ++Nupper;
         else
            ++Nlower;
      }
   }

   if( SCIPisGE(scip, getGeneratorEntry(mastervar1, origvar), ivalue) && SCIPisGE(scip, getGeneratorEntry(mastervar2, origvar), ivalue) )
   {
      k=0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nupper) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nupper) );
      for( j=0; j< NBoundsequences; ++j )
      {
         if( sequencesizes[j] >= p )
            assert(C[j][p-1].component == origvar);

         if( sequencesizes[j] >= p && C[j][p-1].sense == GCG_COMPSENSE_GE )
         {
            CopyC[k] = NULL;
            SCIP_CALL( SCIPallocMemoryArray(scip, &(CopyC[k]), sequencesizes[j]) );
            for( l=0; l< sequencesizes[j]; ++l )
            {
               CopyC[k][l] = C[j][l];
            }
            newsequencesizes[k] = sequencesizes[j];
            ++k;
         }
      }

      if( k != Nupper )
      {
         SCIPdebugMessage("k = %d, Nupper+1 =%d\n", k, Nupper+1);
      }

      if( Nupper != 0 )
         assert( k == Nupper );

      returnvalue = ILOcomp( scip, mastervar1, mastervar2, CopyC, Nupper, newsequencesizes, p+1);

      for( j=0; j< Nupper; ++j )
      {
         SCIPfreeMemoryArray(scip, &(CopyC[j]) );
      }
      SCIPfreeMemoryArray(scip, &newsequencesizes);
      SCIPfreeMemoryArray(scip, &CopyC);

      return returnvalue;
   }


   if( SCIPisLT(scip, getGeneratorEntry(mastervar1, origvar), ivalue) && SCIPisLT(scip, getGeneratorEntry(mastervar2, origvar), ivalue) )
   {
      k=0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nlower) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nlower) );
      for( j=0; j< NBoundsequences; ++j )
      {
         if( sequencesizes[j] >= p )
            assert(C[j][p-1].component == origvar);

         if( sequencesizes[j] >= p && C[j][p-1].sense != GCG_COMPSENSE_GE )
         {
            CopyC[k] = NULL;
            SCIP_CALL( SCIPallocMemoryArray(scip, &(CopyC[k]), sequencesizes[j]) );
            for( l=0; l< sequencesizes[j]; ++l )
            {
               CopyC[k][l] = C[j][l];
            }
            newsequencesizes[k] = sequencesizes[j];
            ++k;
         }
      }

      if( k != Nlower )
      {
         SCIPdebugMessage("k = %d, Nlower+1 =%d\n", k, Nlower+1);
      }

      if( Nlower != 0 )
         assert( k == Nlower);

      returnvalue = ILOcomp( scip, mastervar1, mastervar2, CopyC, Nlower, newsequencesizes, p+1);

      for( j=0; j< Nlower; ++j )
      {

         SCIPfreeMemoryArray(scip, &(CopyC[j]) );
      }
      SCIPfreeMemoryArray(scip, &newsequencesizes);
      SCIPfreeMemoryArray(scip, &CopyC);

      return returnvalue;
   }
   if( SCIPisGT(scip, getGeneratorEntry(mastervar1, origvar), getGeneratorEntry(mastervar2, origvar)) )
      return 1;
   else
      return -1;
}

/** comparefunction for induced lexicographical sort */
static
SCIP_DECL_SORTPTRCOMP(ptrilocomp)
{
   GCG_STRIP* strip1;
   GCG_STRIP* strip2;
   int returnvalue;

   strip1 = (GCG_STRIP*) elem1;
   strip2 = (GCG_STRIP*) elem2;

   returnvalue = ILOcomp( strip1->scip, strip1->mastervar, strip2->mastervar, strip1->C, strip1->Csize, strip1->sequencesizes, 1);

   return returnvalue;
}

/** induced lexicographical sort */
static
SCIP_RETCODE InducedLexicographicSort(
   SCIP*                scip,               /**< SCIP ptr*/
   GCG_STRIP**          array,              /**< array of strips to sort in ILO*/
   int                  arraysize,          /**< size of the array */
   GCG_COMPSEQUENCE**   C,                  /**< current set o comp bound sequences*/
   int                  NBoundsequences,    /**< size of C */
   int*                 sequencesizes       /**< sizes of the sequences in C */
   )
{
   int i;

   SCIPdebugMessage("Induced Lexicographic sorting\n");

   if( NBoundsequences == 0 )
      return LexicographicSort( array, arraysize );
   assert( C!= NULL );

   assert(arraysize > 0);

   if( arraysize <= 1 )
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

/** partitions the strip according to the priority */
static
SCIP_RETCODE partition(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_VAR**           J,                  /**< */
   int*                 Jsize,              /**< */
   int*                 priority,           /**< branching priorities */
   SCIP_VAR**           F,                  /**< set of fractional solutions satisfying bounds */
   int                  Fsize,              /**< size of list of fractional solutions satisfying bounds */
   SCIP_VAR**           origvar,            /**< */
   SCIP_Real*           median              /**< */
   )
{
   int j;
   int l;
   SCIP_Real min;
   SCIP_Real maxPriority;
   SCIP_Real* compvalues;

   //SCIPdebugMessage("Jsize = %d\n", *Jsize);

   do
   {
      min = INT_MAX;
      maxPriority = INT_MIN;

      //max-min priority
      for ( j = 0; j < *Jsize; ++j )
      {
         if ( priority[j] > maxPriority && SCIPvarGetType(J[j]) != SCIP_VARTYPE_CONTINUOUS )
         {
            maxPriority = priority[j];
            *origvar = J[j];
         }
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
      for ( l = 0; l < Fsize; ++l )
      {
         compvalues[l] = getGeneratorEntry(F[l], *origvar);
         if ( SCIPisLT(scip, compvalues[l], min) )
            min = compvalues[l];
      }
      *median = GetMedian(scip, compvalues, Fsize, min);
      SCIPfreeBufferArray(scip, &compvalues);

      assert(min != INT_MAX);

      if ( !SCIPisEQ(scip, *median, 0) )
      {
         SCIPdebugMessage("median = %g\n", *median);
         SCIPdebugMessage("min = %g\n", min);
         SCIPdebugMessage("Jsize = %d\n", *Jsize);
      }

      if ( SCIPisEQ(scip, *median, min) )
      {
         //here with max-min priority
         for ( j = 0; j < *Jsize; ++j )
         {
            if ( *origvar == J[j] )
            {
               assert(priority[j] == 0);
               J[j] = J[*Jsize - 1];
               priority[j] = priority[*Jsize - 1];
               break;
            }
         }
         if( j < *Jsize )
            *Jsize = *Jsize-1;
      }
      assert(*Jsize >= 0);

   }while ( SCIPisEQ(scip, *median, min) && *Jsize > 0);

   return SCIP_OKAY;
}

/** add identified sequence to record */
static
SCIP_RETCODE addToRecord(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_RECORD*          record,             /**< record of identified sequences */
   GCG_COMPSEQUENCE*    S,                  /**< bound restriction sequence */
   int                  Ssize               /**< size of bound restriction sequence */
)
{
   int i;

   SCIPdebugMessage("recordsize=%d, Ssize=%d\n", record->recordsize, Ssize);

   if( record->recordsize == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(record->record), 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(record->sequencesizes), 1) );
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(record->record), record->recordsize+1) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(record->sequencesizes), record->recordsize+1) );
   }
   SCIP_CALL( SCIPallocMemoryArray(scip, &(record->record[record->recordsize]), Ssize) );
   for(i=0; i<Ssize;++i)
   {
      record->record[record->recordsize][i].component = S[i].component;
      record->record[record->recordsize][i].sense = S[i].sense;
      record->record[record->recordsize][i].bound = S[i].bound;
   }

   record->sequencesizes[record->recordsize] = Ssize; //+1 ?

   record->recordsize++;

//   SCIPdebugMessage("add to record with ");
//   for( i=0; i<Ssize; ++i)
//   {
//      SCIPdebugPrintf(" S[%d].component =%s", i, SCIPvarGetName(record->record[record->recordsize-1][i].component) );
//      SCIPdebugPrintf(" S[%d].sense =%d", i, record->record[record->recordsize-1][i].sense );
//      SCIPdebugPrintf(" S[%d].bound =%g\n", i, record->record[record->recordsize-1][i].bound );
//   }

   return SCIP_OKAY;
}


/** separation at the root node */
static
SCIP_RETCODE Separate(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_VAR**           F,                  /**< fractional strips respecting bound restrictions */
   int                  Fsize,              /**< size of the strips */
   SCIP_VAR**           IndexSet,           /**< index set */
   int                  IndexSetSize,       /**< size of index set */
   GCG_COMPSEQUENCE*    S,                  /**< ordered set of bound restrictions */
   int                  Ssize,              /**< size of the ordered set */
   GCG_RECORD*          record              /**< identified bound sequences */
   )
{
   int j;
   int k;
   int l;
   int Jsize;
   SCIP_VAR** J;
   SCIP_Real median;
   SCIP_Real min;
   SCIP_Real max;
   int Fupper;
   int Flower;
   int* priority;
   SCIP_VAR* origvar;
   SCIP_VAR** copyF;
   GCG_COMPSEQUENCE* upperLowerS;
   GCG_COMPSEQUENCE* upperS;
   SCIP_Real* alpha;
   SCIP_Real* compvalues;
   SCIP_Real  muF;
   SCIP_Bool found;

   assert(scip != NULL);
   assert((Fsize == 0) == (F == NULL));
   assert((IndexSetSize == 0) == (IndexSet == NULL));

   j = 0;
   k = 0;
   l = 0;
   Jsize = 0;
   Fupper = 0;
   Flower = 0;
   max = 0;
   muF = 0;
   min = INT_MAX;
   found = FALSE;
   priority = NULL;
   compvalues = NULL;
   J = NULL;
   origvar = NULL;
   copyF = NULL;
   upperLowerS = NULL;
   upperS = NULL;
   alpha = NULL;

   SCIPdebugMessage("Separate with ");

   /* if there are no fractional columns or potential columns, return */
   if( Fsize == 0 || IndexSetSize == 0 )
   {
      SCIPdebugPrintf("nothing, no fractional columns\n");
      return SCIP_OKAY;
   }

   SCIPdebugPrintf("Fsize = %d; Ssize = %d, IndexSetSize = %d\n", Fsize, Ssize, IndexSetSize);

   assert( F != NULL );
   assert( IndexSet != NULL );

   max = getMaxGeneratorEntry(scip, F, Fsize, IndexSet, IndexSetSize);

   if( max == 0 )
      max = 1;

   SCIPdebugMessage("max = %g\n", max);

   for( j=0; j<Fsize; ++j )
      muF += max * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);

   /* detect fractional alpha_i */
   SCIP_CALL( SCIPallocBufferArray(scip, &alpha, IndexSetSize) );

   for( k=0; k < IndexSetSize; ++k )
   {
      GCG_COMPSEQUENCE* copyS;
      SCIP_Real mu_F;
      SCIP_Bool even;
      SCIP_Real alphacontrol;
      SCIP_Real mucontrol;

      even = TRUE;
      mu_F = 0;
      origvar = IndexSet[k];
      copyS = NULL;
      alpha[k] = 0;
      alphacontrol = 0;
      mucontrol = 0;

      if( SCIPvarGetType(origvar) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
      for( l=0; l<Fsize; ++l )
      {
         compvalues[l] = getGeneratorEntry(F[l], origvar);
         if( SCIPisLT(scip, compvalues[l], min) )
            min = compvalues[l];
      }

      median = GetMedian(scip, compvalues, Fsize, min);
      SCIPfreeBufferArray(scip, &compvalues);
      compvalues = NULL;


      for( j = 0; j < Fsize; ++j )
      {
         SCIP_Real generatorentry;

         generatorentry = getGeneratorEntry(F[j], origvar);
         alpha[k] += generatorentry * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);

         if( SCIPisGE(scip, generatorentry, median) )
         {
            alphacontrol += generatorentry * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);
            mucontrol += SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);
         }
//         if( !SCIPisEQ(scip, generatorentry, 0) && !SCIPisEQ(scip, SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]), 0) )
//         {
//             SCIPdebugMessage("generatorentry(%s) = %g\n", SCIPvarGetName(origvar), generatorentry);
//             SCIPdebugMessage("mastervarvalue = %g\n", SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]));
//         }
      }
      if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
      {
         ++Jsize;
//          SCIPdebugMessage("alpha[%d] = %g\n", k, alpha[k]);
      }
      if( !SCIPisFeasIntegral(scip, alpha[k]) || !SCIPisFeasIntegral(scip, alphacontrol) //SCIPisGT(scip, alpha[k] - SCIPfloor(scip, alpha[k]), 0) || SCIPisGT(scip, alphacontrol - SCIPfloor(scip, alphacontrol), 0)
            || !SCIPisFeasIntegral(scip, mucontrol) )//  || SCIPisGT(scip, mucontrol - SCIPfloor(scip, mucontrol), 0))
      {
         SCIPdebugMessage("alpha[%d] = %g\n", k, alpha[k]);
         SCIPdebugMessage("alphacontrol = %g\n", alphacontrol);
         SCIPdebugMessage("mucontrol = %g\n", mucontrol);
         found = TRUE;
         /* ********************************** *
          *   add the current pair to record   *
          * ********************************** */

         /** @todo extract function */
         /* copy S */
         SCIP_CALL( SCIPallocMemoryArray(scip, &copyS, Ssize+1) );
         for( l=0; l < Ssize; ++l )
         {
            copyS[l] = S[l];
//            SCIPdebugMessage("copyS[%d].component = %s \n", l, SCIPvarGetName(copyS[l].component) );
//            SCIPdebugMessage("copyS[%d].sense = %d \n", l, copyS[l].sense );
//            SCIPdebugMessage("copyS[%d].bound = %g \n", l, copyS[l].bound );
         }

         /* create temporary array to compute median */
         SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
         for( l=0; l<Fsize; ++l )
         {
            compvalues[l] = getGeneratorEntry(F[l], origvar);
            if( SCIPisLT(scip, compvalues[l], min) )
               min = compvalues[l];
         }
         assert(median == GetMedian(scip, compvalues, Fsize, min));
         median = GetMedian(scip, compvalues, Fsize, min);
         SCIPfreeBufferArray(scip, &compvalues);
         compvalues = NULL;

         /** @todo mb: this is a fix for an issue that Marcel claims that Vanderbeck did wrong */
         j = 0;

         do
         {
            mu_F = 0;
            if( even )
            {
               median = median+j;
               even = FALSE;
            }
            else
            {
               median = median-j;
               even = TRUE;
            }

            for( l=0; l<Fsize; ++l )
            {
               if( SCIPisGE(scip, getGeneratorEntry(F[l], origvar), median) )
                  mu_F += SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[l]);
            }
            ++j;

         }while( SCIPisFeasIntegral(scip, mu_F) );

         SCIPdebugMessage("new median is %g, comp=%s, Ssize=%d\n", median, SCIPvarGetName(origvar), Ssize);

         /* add last bound change to the copy of S */
         copyS[Ssize].component = origvar;
         copyS[Ssize].sense = GCG_COMPSENSE_GE;
         copyS[Ssize].bound = median;

         /* add identified sequence to record */
         SCIP_CALL( addToRecord(scip, record, copyS, Ssize+1) );


         /* ********************************** *
          *  end adding to record              *
          * ********************************** */

      }
   }

   if( found )
   {
      if(alpha != NULL)
         SCIPfreeBufferArray(scip, &alpha);

      SCIPdebugMessage("one S found with size %d\n", record->sequencesizes[record->recordsize-1]);

      return SCIP_OKAY;
   }


   /* ********************************** *
    *  discriminating components         *
    * ********************************** */

   /** @todo mb: this is a filter */
   SCIP_CALL( SCIPallocMemoryArray(scip, &J, Jsize) );
   j=0;
   for( k=0; k<IndexSetSize; ++k )
   {
      if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
      {
         J[j] = IndexSet[k];
         ++j;
      }
   }
   assert( j == Jsize );

   /* ********************************** *
    *  compute priority  (max-min)       *
    * ********************************** */

   SCIP_CALL( SCIPallocMemoryArray(scip, &priority, Jsize) );
   for( j=0; j<Jsize; ++j )
   {
      int maxcomp;
      int mincomp;

      maxcomp = INT_MIN;
      mincomp = INT_MAX;

      origvar = J[j];

      for( l=0; l<Fsize; ++l )
      {
         SCIP_Real generatorentry;

         generatorentry = getGeneratorEntry(F[l], origvar);

         if( generatorentry > maxcomp )
            maxcomp = generatorentry;

         if( generatorentry < mincomp )
            mincomp = generatorentry;
      }
      priority[j] = maxcomp-mincomp;
   }

   //SCIPdebugMessage("Partitioning\n");
   SCIP_CALL( partition(scip, J, &Jsize, priority, F, Fsize, &origvar, &median) );

   /** @todo mb: this is a copy of S for the recursive call below */
   SCIP_CALL( SCIPallocMemoryArray(scip, &upperLowerS, Ssize+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &upperS, Ssize+1) );
   for( l=0; l < Ssize; ++l )
   {
      upperLowerS[l] = S[l];
      upperS[l] = S[l];
   }

   upperLowerS[Ssize].component = origvar;//i;
   upperS[Ssize].component = origvar;
   upperLowerS[Ssize].sense = GCG_COMPSENSE_LT; //GE?
   upperS[Ssize].sense = GCG_COMPSENSE_GE; //LT?
   upperLowerS[Ssize].bound = median;
   upperS[Ssize].bound = median;

   for( k=0; k<Fsize; ++k )
   {
      if( SCIPisGE(scip, getGeneratorEntry(F[k], origvar), median) )
         ++Fupper;
      else
         ++Flower;
   }

   /* ********************************** *
    *  choose smallest partition         *
    * ********************************** */

   SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fsize) );
   j = 0;

   //binary ? then only one part
   //if( max == 1 )
   /*
   if( (Flower <= Fupper && Flower > 0) || Fupper == 0 )
   {
      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisLT(scip, getGeneratorEntry(scip, F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Flower;
   }
   else
   {
      upperLowerS[Ssize].sense = GCG_COMPSENSE_LT;

      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisGE(scip, getGeneratorEntry(scip, F[k], origvar), median) )//F[k]->generator[i], median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Fupper;
   }
   */

   if( Flower > 0 )
   {
      j = 0;

      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisLT(scip, getGeneratorEntry(F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      //Fsize = Flower;
      assert(j < Fsize+1);
      if( Jsize == 0 && J != NULL )
      {
         SCIPfreeMemoryArray(scip, &J);
         J = NULL;
      }
      Separate( scip, copyF, Flower, J, Jsize, upperLowerS, Ssize+1, record );
   }

   if( Fupper > 0)
   {
      upperLowerS[Ssize].sense = GCG_COMPSENSE_GE; //LT?
      j = 0;

      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisGE(scip, getGeneratorEntry(F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      //Fsize = Fupper;
      assert(j < Fsize+1);
      if( Jsize == 0 && J != NULL )
      {
         SCIPfreeMemoryArray(scip, &J);
         J = NULL;
      }
      Separate( scip, copyF, Fupper, J, Jsize, upperS, Ssize+1, record );
   }

   //Separate( scip, copyF, Fsize, J, Jsize, upperLowerS, Ssize+1, record );

   SCIPfreeMemoryArray(scip, &copyF);

   if(upperLowerS != NULL)
      SCIPfreeMemoryArray(scip, &upperLowerS);
   if(upperS != NULL)
         SCIPfreeMemoryArray(scip, &upperS);
   SCIPfreeMemoryArray(scip, &priority);
   if( J != NULL )
      SCIPfreeMemoryArray(scip, &J);
   SCIPfreeBufferArray(scip, &alpha);

   return SCIP_OKAY;
}

/** choose a component bound sequence to create branching */
static
SCIP_RETCODE ChoseS(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_RECORD**         record,             /**< candidate of bound sequences */
   GCG_COMPSEQUENCE**   S,                  /**< pointer to return chosen bound sequence */
   int*                 Ssize               /**< size of the chosen bound sequence */
   )
{
   int minSizeOfMaxPriority;  //needed if the last comp priority is equal to the one in other bound sequences
   int maxPriority;
   int i;
   int Index;

   minSizeOfMaxPriority = INT_MAX;
   maxPriority = INT_MIN;
   i = 0;
   Index = -1;

   SCIPdebugMessage("Chose S \n");

   assert((*record)->recordsize > 0);

   SCIPdebugMessage("recordsize = %d \n", (*record)->recordsize);

   for( i=0; i< (*record)->recordsize; ++i )
   {
      assert((*record)->sequencesizes != NULL );
      assert((*record)->sequencesizes[i] > 0);
      if(maxPriority <= 1 || maxPriority == INT_MIN) // later by pseudocosts e.g.
      {
         if( maxPriority < 1 || maxPriority == INT_MIN )
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
   for( i=0; i< *Ssize;++i )
   {
      (*S)[i] =  (*record)->record[Index][i];
   }

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

/** updates the new set of sequences C in CopyC and the corresponding size array newsequencesizes
 *  returns the size of CopyC */
static
int computeNewSequence(
   int                   Csize,              /**< size of the sequence */
   int                   p,                  /**< index of node */
   SCIP_VAR*             origvar,            /**< another index generatorentry */
   int*                  sequencesizes,      /**< size of the sequences */
   GCG_COMPSEQUENCE**    C,                  /**< original sequence */
   GCG_COMPSEQUENCE**    CopyC,              /**< new sequence */
   int*                  newsequencesizes,   /**< output parameter for the new sequence */
   GCG_COMPSENSE         sense               /**< sense of the comparison */
   )
{
   int j;
   int k;
   for( k = 0, j = 0; j < Csize; ++j )
   {
      if ( sequencesizes[j] >= p )
         assert(C[j][p-1].component == origvar);

      if ( sequencesizes[j] >= p && C[j][p-1].sense == sense )
      {
         CopyC[k] = C[j];
         newsequencesizes[k] = sequencesizes[j];
         ++k;
      }
   }
   return k;
}

/** auxilary function to compute alpha for given index */
static
double computeAlpha(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   Fsize,              /**< size of F */
   GCG_COMPSENSE         isense,             /**< sense of the bound for origvar */
   double                ivalue,             /**< value of the bound for origvar */
   SCIP_VAR*             origvar,            /**< index of the variable */
   SCIP_VAR**            F                   /**< current fractional mastervars*/
   )
{
   int j;
   SCIP_Real alpha_i = 0;
   for ( j = 0; j < Fsize; ++j )
   {
      SCIP_Real generatorentry;

      generatorentry = getGeneratorEntry(F[j], origvar);

      if ( (isense == GCG_COMPSENSE_GE && SCIPisGE(scip, generatorentry, ivalue)) ||
           (isense == GCG_COMPSENSE_LT && SCIPisLT(scip, generatorentry, ivalue)) )
      {
         alpha_i += generatorentry * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);
         if ( !SCIPisEQ(scip, SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]), 0)
            && !SCIPisEQ(scip, generatorentry, 0) )
         {
//            SCIPdebugMessage("generator(%s) = %g\n", SCIPvarGetName(origvar), generatorentry);
//            SCIPdebugMessage("mastervarvalue = %g\n", SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]));
         }
      }
   }

   return alpha_i;
}

/** separation at a node other than the root node */
static
SCIP_RETCODE Explore(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COMPSEQUENCE**   C,                  /**< */
   int                  Csize,              /**< */
   int*                 sequencesizes,      /**< */
   int                  p,                  /**< */
   SCIP_VAR**           F,                  /**< Strip of fractional columns */
   int                  Fsize,              /**< size of the strips */
   SCIP_VAR**           IndexSet,           /**< */
   int                  IndexSetSize,       /**< */
   GCG_COMPSEQUENCE**   S,                  /**< component sequences */
   int*                 Ssize,              /**< length of component sequences */
   GCG_RECORD*          record              /**< */
   )
{
   int j;
   int k;
   int l;
   SCIP_Real ivalue;
   GCG_COMPSENSE isense;
   SCIP_Real median;
   int Fupper;
   int Flower;
   int Cupper;
   int Clower;
   int lowerSsize;
   SCIP_VAR** copyF;
   GCG_COMPSEQUENCE* copyS;
   GCG_COMPSEQUENCE* lowerS;
   GCG_COMPSEQUENCE** CopyC;
   int* newsequencesizes;
   SCIP_VAR* origvar;
   SCIP_Real alpha_i;
   SCIP_Real  muF;
   SCIP_Bool found;
   SCIP_Real nu_F;
   SCIP_Real max;
   SCIP_Real alphacontrol;
   SCIP_Real mucontrol;

   j = 0;
   k = 0;
   l = 0;
   alpha_i = 0;
   muF = 0;
   max = 0;
   Fupper = 0;
   Flower = 0;
   Cupper = 0;
   Clower = 0;
   lowerSsize = 0;
   alphacontrol = 0;
   mucontrol = 0;
   CopyC = NULL;
   newsequencesizes = NULL;
   copyF = NULL;
   CopyC = NULL;
   copyS = NULL;
   origvar = NULL;
   lowerS =  NULL;
   found = FALSE;

   SCIPdebugMessage("Explore\n");

   SCIPdebugMessage("with Fsize = %d, Csize = %d, Ssize = %d, p = %d\n", Fsize, Csize, *Ssize, p);

   /* *************************************** *
    *   if C=Ø, call separate and return that *
    * *************************************** */

   if( C == NULL || Fsize==0 || IndexSetSize==0 || Csize == 0 )
   {
      //SCIPdebugMessage("go to Separate\n");
      Separate( scip, F, Fsize, IndexSet, IndexSetSize, *S, *Ssize, record );

      if( S != NULL && *Ssize > 0 && *S != NULL )
      {
         SCIPfreeMemoryArray(scip, S);
         S = NULL;
         *Ssize = 0;
      }
      return SCIP_OKAY;
   }

   assert(C != NULL);
   assert(Csize > 0);
   assert(F != NULL);
   assert(IndexSet != NULL);
   assert(sequencesizes != NULL);

   /* ******************************************* *
    * find i which is in all S in C on position p *
    * ******************************************* */

   while( sequencesizes[k] < p )
   {
    //  SCIPdebugMessage("sequencesizes[%d] = %d\n", k, sequencesizes[k]);
      ++k;
      if( k >= Csize )
      {
         SCIPdebugMessage("no %dth element bounded\n", p);
         Separate( scip, F, Fsize, IndexSet, IndexSetSize, *S, *Ssize, record );

         if( S != NULL && *Ssize > 0 && *S != NULL )
         {
            SCIPfreeMemoryArray(scip, S);
            S = NULL;
            *Ssize = 0;
         }

         return SCIP_OKAY;
      }
      assert( k < Csize );
   }
   origvar = C[k][p-1].component;
   isense = C[k][p-1].sense;
   ivalue = C[k][p-1].bound;

   assert(origvar != NULL);
   //SCIPdebugMessage("orivar = %s; ivalue = %g\n", SCIPvarGetName(origvar), ivalue);

   max = getMaxGeneratorEntry(scip, F, Fsize, IndexSet, IndexSetSize);

   if( max == 0 )
      max = 1;

   for( j=0; j<Fsize; ++j )
   {
      muF += max * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);
   }

   //SCIPdebugMessage("muF = %g\n", muF);

   /* ******************************************* *
    * compute alpha_i                             *
    * ******************************************* */

   alpha_i = computeAlpha(scip, Fsize, isense, ivalue, origvar, F);

   if( alpha_i == 0 && isense != GCG_COMPSENSE_GE )
   {
      isense = GCG_COMPSENSE_GE;
      alpha_i = computeAlpha(scip, Fsize, isense, ivalue, origvar, F);
   }

//   SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
//   for( l=0; l<Fsize; ++l )
//   {
//      compvalues[l] = getGeneratorEntry(scip, F[l], origvar);
//      if( SCIPisLT(scip, compvalues[l], min) )
//         min = compvalues[l];
//   }

   median = ivalue;//GetMedian(scip, compvalues, Fsize, min);
//   SCIPfreeBufferArray(scip, &compvalues);
//   compvalues = NULL;

   for( j = 0; j < Fsize; ++j )
   {
      SCIP_Real generatorentry;

      generatorentry = getGeneratorEntry(F[j], origvar);

      if( SCIPisGE(scip, generatorentry, median) )
      {
         alphacontrol += generatorentry * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);
         mucontrol += SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);
      }
    }

   //SCIPdebugMessage("alpha(%s) = %g\n", SCIPvarGetName(origvar), alpha_i);

   /* ******************************************* *
    * if f > 0, add pair to record                *
    * ******************************************* */
   if( !SCIPisFeasIntegral(scip, alpha_i) || !SCIPisFeasIntegral(scip, alphacontrol) //SCIPisGT(scip, alpha_i - SCIPfloor(scip, alpha_i), 0)|| SCIPisGT(scip, alphacontrol - SCIPfloor(scip, alphacontrol), 0)
      || !SCIPisFeasIntegral(scip, mucontrol) )//SCIPisGT(scip, mucontrol - SCIPfloor(scip, mucontrol), 0) )
   {
      found = TRUE;
      //SCIPdebugMessage("fractional alpha(%s) = %g\n", SCIPvarGetName(origvar), alpha_i);

      /* ******************************************* *
       * compute nu_F                                *
       * ******************************************* */

      nu_F = 0;
      for( l = 0; l < Fsize; ++l )
      {
         if( (isense == GCG_COMPSENSE_GE && SCIPisGE(scip, getGeneratorEntry(F[l], origvar), ivalue))
            || (isense == GCG_COMPSENSE_LT && SCIPisLT(scip, getGeneratorEntry(F[l], origvar), ivalue)) )
         {
               nu_F += SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[l]);
         }
      }

      /* ******************************************* *
       * add to record                               *
       * ******************************************* */

      if( SCIPisGT(scip, nu_F - SCIPfloor(scip, nu_F), 0) )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &copyS, *Ssize+1) );
         for( l = 0; l < *Ssize; ++l )
         {
            copyS[l] = (*S)[l];
         }
         copyS[*Ssize].component = origvar;
         copyS[*Ssize].sense = isense;
         copyS[*Ssize].bound = ivalue;
         SCIP_CALL( addToRecord(scip, record, copyS, *Ssize+1) );
      }
      else
      {
         found = FALSE;
      }
   }

   if( found )
   {
      SCIPdebugMessage("found fractional alpha\n");
      return SCIP_OKAY;
   }

   /* add bound to the end of S */
   ++(*Ssize);
   if( S==NULL || *S==NULL )
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
   (*S)[*Ssize-1].component = origvar;
   (*S)[*Ssize-1].sense = GCG_COMPSENSE_GE;
   (*S)[*Ssize-1].bound = median;

   SCIP_CALL( SCIPallocMemoryArray(scip, &lowerS, *Ssize) );

   for( k=0; k<*Ssize-1; ++k)
   {
      lowerS[k].component = (*S)[k].component;
      lowerS[k].sense = (*S)[k].sense;
      lowerS[k].bound = (*S)[k].bound;
   }
   lowerSsize = *Ssize;
   lowerS[lowerSsize-1].component = origvar;
   lowerS[lowerSsize-1].sense = GCG_COMPSENSE_LT;
   lowerS[lowerSsize-1].bound = median;

   for( k=0; k<Fsize; ++k )
   {
      if( SCIPisGE(scip, getGeneratorEntry(F[k], origvar), median) )
         ++Fupper;
      else
         ++Flower;
   }

   //calculate subset of C
   for( j=0; j< Csize; ++j )
   {
      if( sequencesizes[j] >= p )
      {
         if( C[j][p-1].sense == GCG_COMPSENSE_GE )
         {
            ++Cupper;
         }
         else
         {
            ++Clower;
            assert( C[j][p-1].sense == GCG_COMPSENSE_LT );
         }
      }
   }

   SCIPdebugMessage("Cupper = %d, Clower = %d\n", Cupper, Clower);

   if( SCIPisLE(scip, alpha_i, 0) && Fupper != 0 )
      Flower = INT_MAX;
   if( SCIPisEQ(scip, alpha_i, muF) && Flower != 0 )
      Fupper = INT_MAX;

   /** @todo mb: one can easily merge these two branches as they are similar */
   /*
   //choose smallest partition
   if( ((Fupper <= Flower && Fupper > 0 ) || Flower <= 0) && Fupper != INT_MAX )
   {
      //SCIPdebugMessage("chose upper bound Fupper = %d, Cupper = %d\n", Fupper, Cupper);

      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
      for( j = 0, k = 0; k < Fsize; ++k )
      {
         if( SCIPisGE(scip, getGeneratorEntry(scip, F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Fupper;

      //new C
      if( Cupper > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Cupper) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Cupper) );
         k = computeNewSequence(Csize, p, origvar, sequencesizes, C, CopyC, newsequencesizes, GCG_COMPSENSE_GE);
      }
      else
      {
         CopyC = NULL;
         k = 0;
      }
      Csize = Cupper;
   }
   else
   {
     // SCIPdebugMessage("chose lower bound Flower = %d Clower = %d\n", Flower, Clower);
      (*S)[*Ssize-1].sense = GCG_COMPSENSE_LT;
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
      j = 0;
      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisLT(scip, getGeneratorEntry(scip, F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      Fsize = Flower;

      //new C
      if( Clower > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Clower) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Clower) );
         k = computeNewSequence(Csize, p, origvar, sequencesizes, C, CopyC, newsequencesizes, GCG_COMPSENSE_LT);
      }
      else
      {
         CopyC = NULL;
         k = 0;
      }
      Csize = Clower;
   }
   assert( k <= Csize+1 );
   */

   if(  Fupper > 0  && Fupper != INT_MAX )
   {
      SCIPdebugMessage("chose upper bound Fupper = %d, Cupper = %d\n", Fupper, Cupper);

      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
      for( j = 0, k = 0; k < Fsize; ++k )
      {
         if( SCIPisGE(scip, getGeneratorEntry(F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      //Fsize = Fupper;

      //new C
      if( Fupper > 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Cupper) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Cupper) );
         k = computeNewSequence(Csize, p, origvar, sequencesizes, C, CopyC, newsequencesizes, GCG_COMPSENSE_GE);
         if( k != Cupper )
         {
            SCIPdebugMessage("k = %d, p = %d\n", k, p);
         }
         assert(k == Cupper);
      }
      else
      {
         CopyC = NULL;
         k = 0;
      }
      //Csize = Cupper;
      Explore( scip, CopyC, Cupper, newsequencesizes, p+1, copyF, Fupper, IndexSet, IndexSetSize, S, Ssize, record );
   }

   if( Flower > 0 && Flower != INT_MIN )
   {
      SCIPdebugMessage("chose lower bound Flower = %d Clower = %d\n", Flower, Clower);
      //(*S)[*Ssize-1].sense = GCG_COMPSENSE_LT;
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
      j = 0;
      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisLT(scip, getGeneratorEntry(F[k], origvar), median) )
         {
            copyF[j] = F[k];
            ++j;
         }
      }
      //Fsize = Flower;

      //new C
      if( Flower > 0 )
      {
         if( CopyC != NULL )
         {
            SCIPfreeMemoryArray(scip, &CopyC);
            CopyC = NULL;
         }

         if( newsequencesizes != NULL )
         {
            SCIPfreeMemoryArray(scip, &newsequencesizes);
            newsequencesizes = NULL;
         }

         SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Clower) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Clower) );
         k = computeNewSequence(Csize, p, origvar, sequencesizes, C, CopyC, newsequencesizes, GCG_COMPSENSE_LT);
         if( k != Clower )
         {
            SCIPdebugMessage("k = %d, p = %d\n", k, p);
         }
         assert(k == Clower);
      }
      else
      {
         CopyC = NULL;
         k = 0;
      }
      //Csize = Clower;
      Explore( scip, CopyC, Clower, newsequencesizes, p+1, copyF, Flower, IndexSet, IndexSetSize, &lowerS, &lowerSsize, record );
   }



   //Explore( scip, CopyC, Csize, newsequencesizes, p+1, copyF, Fsize, IndexSet, IndexSetSize, S, Ssize, record );

   SCIPfreeMemoryArray(scip, &copyF);
   if( lowerS != NULL )
         SCIPfreeMemoryArray(scip, &lowerS);
   if( CopyC != NULL )
      SCIPfreeMemoryArray(scip, &CopyC);

   if( newsequencesizes != NULL )
      SCIPfreeMemoryArray(scip, &newsequencesizes);

   if( S != NULL && *Ssize > 0 && *S != NULL )
   {
      SCIPfreeMemoryArray(scip, S);
      S = NULL;
      *Ssize = 0;
   }

   return SCIP_OKAY;
}

/** callup method for seperate
 * decides whether Separate or Explore should be done */
static
SCIP_RETCODE ChooseSeparateMethod(
   SCIP*                scip,               /**< */
   SCIP_VAR**           F,
   int                  Fsize,              /**< */
   GCG_COMPSEQUENCE**   S,                  /**< */
   int*                 Ssize,              /**< */
   GCG_COMPSEQUENCE**   C,                  /**< */
   int                  Csize,              /**< */
   int*                 CompSizes,           /**< */
   int                  blocknr
   )
{
   SCIP* masterscip;
   int i;
   SCIP_VAR** IndexSet;
   SCIP_VAR** mastervars;
   int IndexSetSize;
   GCG_RECORD* record;
   int exploreSsize;
   GCG_COMPSEQUENCE* exploreS;
   GCG_STRIP** strips;
   int nstrips;
   int nmastervars;

   assert(Fsize > 0);
   assert(F != NULL);
   IndexSetSize = 0;
   exploreSsize = 0;
   exploreS = NULL;
   record = NULL;
   IndexSet = NULL;
   strips = NULL;
   nstrips = 0;

   SCIPdebugMessage("Calling Separate\n");

   SCIP_CALL( SCIPallocBuffer(scip, &record) );
   record->recordsize = 0;
   record->record = NULL;
   record->sequencesizes = NULL;

   //calculate IndexSet
   SCIP_CALL( InitIndexSet(scip, F, Fsize, &IndexSet, &IndexSetSize) );
   assert(IndexSetSize > 0);
   assert(IndexSet != NULL);

//   #ifdef SCIP_DEBUG
//   for( i=0; i<Csize; ++i )
//   {
//      SCIPdebugMessage("Compsizes[%d] = %d\n", i, CompSizes[i]);
//   }
//   #endif

   //rootnode?
   if( Csize<=0 )
      Separate( scip, F, Fsize, IndexSet, IndexSetSize, NULL, 0, record );
   else
   {
      assert( C!=NULL );
      Explore( scip, C, Csize, CompSizes, 1, F, Fsize, IndexSet, IndexSetSize, &exploreS, &exploreSsize, record);
      if( exploreS != NULL )
         SCIPfreeMemoryArray(scip, &exploreS);
   }

   assert(record != NULL);

   if(record->recordsize <= 0)
   {
      masterscip = GCGrelaxGetMasterprob(scip);
      assert(masterscip != NULL);
      SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

      for( i=0; i<nmastervars; ++i )
      {
         SCIP_Bool blockfound;
         SCIP_VAR** pricingvars;
         int u;

         if( GCGvarGetBlock(mastervars[i]) == -1 && GCGvarIsLinking(mastervars[i]) )
         {
            blockfound = FALSE;

            pricingvars = GCGlinkingVarGetPricingVars(mastervars[i]);
            assert(pricingvars != NULL );

            for( u=0; u<GCGlinkingVarGetNBlocks(mastervars[i]); ++u )
            {
               if( pricingvars[u] != NULL )
               {
                  if( GCGvarGetBlock(pricingvars[u]) == blocknr )
                  {
                     blockfound = TRUE;
                     break;
                  }
               }
            }
         }
         else
         {
            blockfound = (GCGvarGetBlock(mastervars[i]) == blocknr);
         }

         if( blockfound )
         {
            ++nstrips;

            if( nstrips == 1 )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &strips, nstrips) );
            }
            else
            {
               SCIP_CALL( SCIPreallocBufferArray(scip, &strips, nstrips) );
            }

            SCIP_CALL( SCIPallocBuffer(scip, &(strips[nstrips-1])) );

            strips[nstrips-1]->C = NULL;
            strips[nstrips-1]->mastervar = mastervars[i];
            strips[nstrips-1]->Csize = 0;
            strips[nstrips-1]->sequencesizes = NULL;
            strips[nstrips-1]->scip = NULL;
         }
      }

      InducedLexicographicSort(scip, strips, nstrips, C, Csize, CompSizes);

      for( i=0; i<nstrips; ++i)
      {
         SCIPfreeBuffer(scip, &(strips[i]));
         strips[i] = NULL;
      }

      if( strips != NULL )
      {
         SCIPfreeBufferArray(scip, &strips);
         strips = NULL;
      }

      //return SCIP_OKAY;
   }


   assert(record->recordsize > 0);

   ChoseS( scip, &record, S, Ssize );
   assert(*S!=NULL);

   #ifdef SCIP_DEBUG
   for( i=0; i<*Ssize; ++i )
   {
//      SCIPdebugMessage("S[%d].component = %s\n", i, SCIPvarGetName((*S)[i].component) );
//      SCIPdebugMessage("S[%d].sense = %d\n", i, (*S)[i].sense);
//      SCIPdebugMessage("S[%d].bound = %.6g\n", i, (*S)[i].bound);
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

   //SCIPdebugMessage("branchDataDeleteGeneric:\n");

   if( *branchdata == NULL )
   {
      SCIPdebugMessage("branchDataDeleteGeneric: cannot delete empty branchdata\n");

      return SCIP_OKAY;
   }

   if( (*branchdata)->mastercons != NULL )
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

   if( (*branchdata)->consS != NULL && (*branchdata)->consSsize > 0 )
   {
      SCIPfreeMemoryArray(scip, &((*branchdata)->consS));
      (*branchdata)->consS = NULL;
   }

   SCIPfreeMemory(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

/** check method for pruning ChildS directly on childnodes
 *  retrun TRUE if node is pruned */
static
SCIP_Bool checkchildconsS(
   SCIP*                 scip,
   SCIP_Real             lhs,               /**< lhs for childnode which is checkes to be pruned */
   GCG_COMPSEQUENCE*     childS,            /**< Component Bound Sequence defining the childnode */
   int                   childSsize,        /**< */
   SCIP_CONS*            parentcons,      /**< */
   int                   childBlocknr       /**< number of the block for the childnode */
   )
{
   int i;
   int nchildren;

   nchildren = GCGconsMasterbranchGetNChildcons(parentcons);
   assert(nchildren>0);

   //SCIPdebugMessage("nchildren = %d\n", nchildren);

   for(i=0; i<nchildren; ++i)
   {
      SCIP_CONS* childcons;
      GCG_BRANCHDATA* branchdata;
      SCIP_Bool same;
      int j;

      //SCIPdebugMessage("childnr = %d\n", i);

      same = TRUE;
      childcons = GCGconsMasterbranchGetChildcons(parentcons, i);
      if(childcons == NULL)
         continue;
      if(GCGconsMasterbranchGetbranchrule(childcons) != NULL)
         assert(strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetbranchrule(childcons)), "generic") == 0);

      branchdata = GCGconsMasterbranchGetBranchdata(childcons);
      assert(branchdata != NULL || GCGconsMasterbranchGetOrigbranchdata(childcons) != NULL);

      if( branchdata == NULL )
         branchdata = GCGconsMasterbranchGetOrigbranchdata(childcons);

      if( childBlocknr != branchdata->consblocknr || childSsize != branchdata->consSsize || !SCIPisEQ(scip, lhs, branchdata->lhs) )
         continue;

      assert(childSsize > 0 && branchdata->consSsize > 0);

      for(j=0; j< childSsize; ++j)
      {
         if( childS[j].component != branchdata->consS[j].component
            || childS[j].sense != branchdata->consS[j].sense
            || !SCIPisEQ(scip, childS[j].bound, branchdata->consS[j].bound) )
         {
            same = FALSE;
            break;
         }
      }

      if( same )
      {
         SCIPdebugMessage("child pruned \n");
         return TRUE;
      }
   }
//   SCIPdebugMessage("child not pruned \n");
   return FALSE;
}

/** check method for pruning ChildS indirectly by parentnodes
 *  retrun TRUE if node is pruned */
static
SCIP_Bool pruneChildNodeByDominanceGeneric(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_Real             lhs,               /**< lhs for childnode which is checkes to be pruned */
   GCG_COMPSEQUENCE*     childS,            /**< Component Bound Sequence defining the childnode */
   int                   childSsize,        /**< */
   SCIP_CONS*            masterbranchcons,      /**< */
   int                   childBlocknr       /**< number of the block for the childnode */
   )
{
   SCIP_CONS* cons;
   SCIP_Bool ispruned;

   ispruned = FALSE;
   cons = NULL;

   SCIPdebugMessage("Prune by dominance\n");
   cons = GCGconsMasterbranchGetParentcons(masterbranchcons);

   if( cons == NULL )
   {
      SCIPdebugMessage("cons == NULL, not pruned\n");
      return FALSE;
   }
   while( cons != NULL)
   {
      GCG_BRANCHDATA* parentdata;

      parentdata = GCGconsMasterbranchGetBranchdata(cons);
      if(parentdata == NULL)
      {
         //root node: check children for pruning
//         SCIPdebugMessage("parentdata == NULL, check children\n");
         return checkchildconsS(scip, lhs, childS, childSsize, cons, childBlocknr);
      }
      if(strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetbranchrule(cons)), "generic") != 0)
         return checkchildconsS(scip, lhs, childS, childSsize, cons, childBlocknr);

      ispruned = checkchildconsS(scip, lhs, childS, childSsize, cons, childBlocknr);

      if( ispruned )
      {
         //SCIPdebugMessage("child pruned 1\n");
         return TRUE;
      }

      cons = GCGconsMasterbranchGetParentcons(cons);
   }

   SCIPdebugMessage("child not pruned\n");
   return FALSE;
}

/** initialize branchdata at the node */
static
SCIP_RETCODE initNodeBranchdata(
   GCG_BRANCHDATA**      nodebranchdata,     /**< branching data to set */
   int                   blocknr             /**< block we are branching in */
   )
{
   SCIP_CALL( SCIPallocMemory(scip, nodebranchdata) );

   (*nodebranchdata)->consblocknr = blocknr;
   (*nodebranchdata)->mastercons = NULL;
   (*nodebranchdata)->consS = NULL;
   (*nodebranchdata)->C = NULL;
   (*nodebranchdata)->sequencesizes = NULL;
   (*nodebranchdata)->Csize = 0;
   (*nodebranchdata)->consSsize = 0;

   return SCIP_OKAY;
}

/** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_RETCODE createChildNodesGeneric(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,        /**< branching rule */
   GCG_COMPSEQUENCE*     S,                 /**< Component Bound Sequence defining the nodes */
   int                   Ssize,             /**< size of S*/
   int                   blocknr,           /**< number of the block */
   SCIP_CONS*            masterbranchcons,   /**< current masterbranchcons*/
   SCIP_RESULT*          result
   )
{
   SCIP*  masterscip;
   int i;
   int p;
   int k;
   SCIP_Real pL;
   SCIP_Real L;
   SCIP_Real lhs;
   SCIP_Real lhsSum;
   int nmastervars;
   int nmastervars2;
   int ncopymastervars;
   int nbranchcands;
   int nchildnodes;
   SCIP_Real mu;  // mu(S)
   SCIP_Real identicalcontrol;
   SCIP_VAR** mastervars;
   SCIP_VAR** mastervars2;
   SCIP_VAR** branchcands;
   SCIP_VAR** copymastervars;

   assert(scip != NULL);
   assert(Ssize > 0);
   assert(S != NULL);

   lhs = 0;
   lhsSum = 0;
   nchildnodes = 0;
   p = 0;
   k = 0;
   i = 0;
   L = 0;
   mu = 0;
   identicalcontrol = 0;

   SCIPdebugMessage("Vanderbeck branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, GCGrelaxGetNIdenticalBlocks(scip, blocknr));
   pL = GCGrelaxGetNIdenticalBlocks(scip, blocknr);

   // get variable data of the master problem
   masterscip = GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   nmastervars2 = nmastervars;
   assert(nmastervars >= 0);

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &copymastervars, mastervars, nmastervars) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &mastervars2, mastervars, nmastervars) );

   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

   SCIPdebugMessage("Vanderbeck branching rule: creating %d nodes\n", Ssize+1);

   for( p=0; p<Ssize+1; ++p )
   {
      GCG_BRANCHDATA* branchchilddata;
      SCIP_NODE* child;
      SCIP_CONS* childcons;
      char childname[SCIP_MAXSTRLEN];

      mu = 0;
      branchchilddata = NULL;

      // allocate branchdata for child and store information
      SCIP_CALL( initNodeBranchdata(&branchchilddata, blocknr) );

      //SCIPdebugMessage("p = %d \n", p);

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
      for( k=0; k<=p; ++k )
      {
         GCG_COMPSEQUENCE compBound;

         if( k == Ssize )
         {
            assert( p == Ssize );
            compBound = S[k-1];
            branchchilddata->consS[k-1] = compBound;
         }
         else
         {
            compBound = S[k];
            if( k >= p )
            {
               if( S[p].sense == GCG_COMPSENSE_GE )
                  compBound.sense = GCG_COMPSENSE_LT;
               else
                  compBound.sense = GCG_COMPSENSE_GE;
            }
            branchchilddata->consS[k] = compBound;
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
         ncopymastervars = nmastervars2;
         for( i=0; i<ncopymastervars; ++i )
         {
            SCIP_VAR* swap;
            SCIP_VAR** pricingvars;
            SCIP_Real generator_i;
            SCIP_Bool blockfound;
            int u;

            blockfound = TRUE;
            pricingvars = NULL;
            u = 0;

            if(i >= nmastervars2)
               break;

            if( GCGvarGetBlock(mastervars2[i]) == -1 )
            {
               assert( GCGvarIsLinking(mastervars2[i]) );
               blockfound = FALSE;

               pricingvars = GCGlinkingVarGetPricingVars(mastervars2[i]);
               assert(pricingvars != NULL );

               for( u=0; u<GCGlinkingVarGetNBlocks(mastervars2[i]); ++u )
               {
                  if( pricingvars[u] != NULL )
                  {
                     if( GCGvarGetBlock(pricingvars[u]) == blocknr )
                     {
                        blockfound = TRUE;
                        break;
                     }
                  }
               }
            }
            else
            {
               blockfound = (GCGvarGetBlock(mastervars2[i]) == blocknr);
            }

            if( blockfound )
            {
               generator_i = getGeneratorEntry(mastervars2[i], S[p].component);

               if( (S[p].sense == GCG_COMPSENSE_GE && SCIPisGE(scip, generator_i, S[p].bound)) ||
                  (S[p].sense == GCG_COMPSENSE_LT && SCIPisLT(scip, generator_i, S[p].bound) ) )
               {
                  mu += SCIPgetSolVal(masterscip, NULL, mastervars2[i]);
               }
               else if( ncopymastervars > 0 )
               {
                  swap = mastervars2[i];
                  mastervars2[i] = mastervars2[nmastervars2-1];//ncopymastervars-1];
                  mastervars2[nmastervars2-1] = swap;
                  --nmastervars2;//ncopymastervars;
                  --i;
               }
            }
            else if( nmastervars2 > 0)//ncopymastervars > 0 )
            {
               swap = mastervars2[i];
               mastervars2[i] = mastervars2[nmastervars2-1];//ncopymastervars-1];
               mastervars2[nmastervars2-1] = swap;
               --nmastervars2;//ncopymastervars;
               --i;
            }
         }
         //nmastervars2 = ncopymastervars;

         if( p == Ssize-1 )
         {
            L = SCIPceil(scip, mu);
            SCIPdebugMessage("mu = %g, \n", mu);
            assert(!SCIPisFeasIntegral(scip,mu));
         }
         else
         {
            L = mu;
            SCIPdebugMessage("mu = %g should be integer, \n", mu);
            assert(SCIPisFeasIntegral(scip,mu)); //SCIPisEQ(scip, mu - SCIPceil(scip, mu), 0));
         }
         lhs = pL-L+1;
         // SCIPdebugMessage("pL-L+1 = %g \n", pL-L+1);
      }
      SCIPdebugMessage("pL = %g \n", pL);
      pL = L;

      branchchilddata->lhs = lhs;
      SCIPdebugMessage("L = %g, \n", L);
      SCIPdebugMessage("lhs set to %g \n", lhs);
      assert(SCIPisFeasIntegral(scip, lhs));//SCIPisEQ(scip, lhs - SCIPceil(scip, lhs), 0) && lhs != 0);
      lhsSum += lhs;

      /* define names for origbranch constraints */
      (void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "node(%d,%d, %f) last comp=%s, sense %d, bound %g", p+1, blocknr, lhs,
         SCIPvarGetName(branchchilddata->consS[branchchilddata->consSsize-1].component),
         branchchilddata->consS[branchchilddata->consSsize-1].sense,
         branchchilddata->consS[branchchilddata->consSsize-1].bound);

      if( masterbranchcons == NULL || !pruneChildNodeByDominanceGeneric(scip, lhs, branchchilddata->consS, branchchilddata->consSsize, masterbranchcons, blocknr) )
      {
         if( masterbranchcons != NULL )
         {
            ++nchildnodes;
            for( i=0; i<branchchilddata->consSsize; ++i )
            {
//               SCIPdebugMessage("setting consS[%d].component = %s\n", i, SCIPvarGetName(branchchilddata->consS[i].component) );
//               SCIPdebugMessage("setting consS[%d].sense = %d\n", i, branchchilddata->consS[i].sense);
//               SCIPdebugMessage("setting consS[%d].bound = %g\n", i, branchchilddata->consS[i].bound);

            }
            SCIP_CALL( SCIPcreateChild(masterscip, &child, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
            SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &childcons, child, GCGconsMasterbranchGetActiveCons(masterscip)) );
            SCIP_CALL( SCIPaddConsNode(masterscip, child, childcons, NULL) );

            assert(branchchilddata != NULL);
            SCIP_CALL( GCGconsMasterbranchSetOrigConsData(masterscip, childcons, childname, branchrule, branchchilddata,
               NULL, 0, FALSE, FALSE, FALSE, NULL, 0, NULL, 0) );

            // release constraints
            SCIP_CALL( SCIPreleaseCons(masterscip, &childcons) );
            //SCIPdebugMessage("branchchilddata->consSsize = %d\n", branchchilddata->consSsize);
         }
      }
      else
         SCIPfreeMemory(scip, &branchchilddata);
   }
   SCIPdebugMessage("lhsSum = %g\n", lhsSum);

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

  for( i=0; i<nmastervars; ++i )
  {
     SCIP_VAR* mastervar;
     SCIP_VAR** pricingvars;
     SCIP_Bool blockfound;
     int u;

     mastervar = mastervars[i];
     pricingvars = NULL;
     blockfound = FALSE;
     u = 0;

     if( GCGvarGetBlock(mastervar) == -1 && GCGvarIsLinking(mastervar) )
     {
        assert( GCGvarIsLinking(mastervar) );
        blockfound = FALSE;

        pricingvars = GCGlinkingVarGetPricingVars(mastervar);
        assert(pricingvars != NULL );

        for( u=0; u<GCGlinkingVarGetNBlocks(mastervar); ++u )
        {
           if( pricingvars[u] != NULL )
           {
              if( GCGvarGetBlock(pricingvars[u]) == blocknr )
              {
                 blockfound = TRUE;
                 break;
              }
           }
        }
     }
     else
     {
        blockfound = (GCGvarGetBlock(mastervar) == blocknr);
     }

     if( blockfound )
     {
        identicalcontrol += SCIPgetSolVal(masterscip, NULL, mastervar);
     }

  }
  if(!SCIPisEQ(scip, identicalcontrol, GCGrelaxGetNIdenticalBlocks(scip, blocknr)))
  {
     SCIPdebugMessage("width of the block is only %g\n", identicalcontrol);
  }

  assert( SCIPisEQ(scip, identicalcontrol, GCGrelaxGetNIdenticalBlocks(scip, blocknr)) );
#endif

   assert( SCIPisEQ(scip, lhsSum, GCGrelaxGetNIdenticalBlocks(scip, blocknr) + Ssize) );

   SCIPfreeMemoryArray(scip, &mastervars2);
   SCIPfreeMemoryArray(scip, &copymastervars);

   if( nchildnodes <= 0 )
   {
      SCIPdebugMessage("node cut off, since all childnodes have been pruned\n");

      *result = SCIP_CUTOFF;
   }

   return SCIP_OKAY;
}

/** branching on copied origvar directly in master
 * @return SCIP_RETCODE */
static
SCIP_RETCODE branchDirectlyOnMastervar(
   SCIP*                scip,
   SCIP_VAR*            mastervar,
   SCIP_BRANCHRULE*     branchrule
   )
{
   SCIP* masterscip;
   GCG_BRANCHDATA* branchupchilddata;
   GCG_BRANCHDATA* branchdownchilddata;
   SCIP_NODE* upchild;
   SCIP_NODE* downchild;
   SCIP_CONS* upchildcons;
   SCIP_CONS* downchildcons;
   char upchildname[SCIP_MAXSTRLEN];
   char downchildname[SCIP_MAXSTRLEN];
   int bound;

   masterscip = GCGrelaxGetMasterprob(scip);
   assert(masterscip != NULL);

   bound = SCIPceil( scip, SCIPgetSolVal(masterscip, NULL, mastervar));

   // allocate branchdata for child and store information
   SCIP_CALL( initNodeBranchdata(&branchupchilddata, -3) );
   SCIP_CALL( initNodeBranchdata(&branchdownchilddata, -3) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(branchupchilddata->consS), 1) );
   branchupchilddata->consSsize = 1;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(branchdownchilddata->consS), 1) );
      branchdownchilddata->consSsize = 1;

   branchupchilddata->consS[0].component = mastervar;
   branchupchilddata->consS[0].sense = GCG_COMPSENSE_GE;
   branchupchilddata->consS[0].bound = bound;

   branchdownchilddata->consS[0].component = mastervar;
   branchdownchilddata->consS[0].sense = GCG_COMPSENSE_LT;
   branchdownchilddata->consS[0].bound = bound;


   (void) SCIPsnprintf(upchildname, SCIP_MAXSTRLEN, "node(1,-3, %f) direct up on comp=%s", branchupchilddata->consS[0].bound,
               SCIPvarGetName(branchupchilddata->consS[branchupchilddata->consSsize-1].component));
   (void) SCIPsnprintf(downchildname, SCIP_MAXSTRLEN, "node(1,-3, %f) direct up on comp=%s", branchdownchilddata->consS[0].bound,
               SCIPvarGetName(branchdownchilddata->consS[branchdownchilddata->consSsize-1].component));

   SCIP_CALL( SCIPcreateChild(masterscip, &upchild, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &upchildcons, upchild, GCGconsMasterbranchGetActiveCons(masterscip)) );
   SCIP_CALL( SCIPaddConsNode(masterscip, upchild, upchildcons, NULL) );

   SCIP_CALL( SCIPcreateChild(masterscip, &downchild, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &downchildcons, downchild, GCGconsMasterbranchGetActiveCons(masterscip)) );
   SCIP_CALL( SCIPaddConsNode(masterscip, downchild, downchildcons, NULL) );

   assert(branchupchilddata != NULL);
   SCIP_CALL( GCGconsMasterbranchSetOrigConsData(masterscip, upchildcons, upchildname, branchrule, branchupchilddata,
      NULL, 0, FALSE, FALSE, FALSE, NULL, 0, NULL, 0) );

   assert(branchdownchilddata != NULL);
   SCIP_CALL( GCGconsMasterbranchSetOrigConsData(masterscip, downchildcons, downchildname, branchrule, branchdownchilddata,
      NULL, 0, FALSE, FALSE, FALSE, NULL, 0, NULL, 0) );

   // release constraints
   SCIP_CALL( SCIPreleaseCons(masterscip, &upchildcons) );
   //SCIPdebugMessage("branchchilddata->consSsize = %d\n", branchupchilddata->consSsize);

   SCIP_CALL( SCIPreleaseCons(masterscip, &downchildcons) );
   //SCIPdebugMessage("branchchilddata->consSsize = %d\n", branchdownchilddata->consSsize);

   return SCIP_OKAY;
}

/** prepares informations for using the generic branching scheme
 * @return SCIP_RETCODE */
static
SCIP_RETCODE GCGbranchGenericInitbranch(
   SCIP*                masterscip,         /**< */
   SCIP_BRANCHRULE*     branchrule,
   SCIP_RESULT*         result
   )
{
   SCIP* origscip;
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
   GCG_COMPSEQUENCE* S;
   GCG_COMPSEQUENCE** C;
   SCIP_VAR** F;
   SCIP_VAR** pricingvars;
   int Ssize;
   int Csize;
   int Fsize;
   int* sequencesizes;
   int blocknr;
   int pricingblocknr;
   int nidenticalpricing;
   int i;
   int c;
   int allnorigvars;

   blocknr = -2;
   Ssize = 0;
   Csize = 0;
   Fsize = 0;
   pricingblocknr = -2;
   nidenticalpricing = 0;
   i = 0;
   c = 0;
   feasible = TRUE;
   SinC = TRUE;
   S = NULL;
   F = NULL;
   branchdata = NULL;
   S = NULL;
   C = NULL;
   F = NULL;
   sequencesizes = NULL;
   pricingvars = NULL;

   assert(masterscip != NULL);

   SCIPdebugMessage("get informations for Vanderbecks generic branching\n");

   origscip = GCGpricerGetOrigprob(masterscip);

   assert(origscip != NULL);
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   assert(nbranchcands > 0);

   for( i=0; i<nbranchcands; ++i )
   {
      mastervar = branchcands[i];
      assert(GCGvarIsMaster(mastervar));
      blocknr = GCGvarGetBlock(mastervar);
      if( blocknr >= -1 )
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

   //a special case; branch on copy of an origvar directly:- here say blocknr = -3
   if( blocknr == -1 && !GCGvarIsLinking(mastervar) )
      blocknr = -3;

   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   SCIPdebugMessage("branching in block %d \n", blocknr);

   if(blocknr == -1)
   {
     assert( GCGvarIsLinking(mastervar) );

      pricingvars = GCGlinkingVarGetPricingVars(mastervar);
      assert(pricingvars != NULL );

      for( i=0; i<GCGlinkingVarGetNBlocks(mastervar); ++i )
      {
         if(pricingvars[i] != NULL)
         {
            if(pricingblocknr == -2 )
            {
               pricingblocknr = GCGvarGetBlock(pricingvars[i]);
               nidenticalpricing = GCGrelaxGetNIdenticalBlocks(origscip, pricingblocknr);
            }
            else
            {
               assert( nidenticalpricing == GCGrelaxGetNIdenticalBlocks(origscip, GCGvarGetBlock(pricingvars[i])));
            }
         }
      }
      assert(pricingblocknr > -1);

      blocknr = pricingblocknr;
   }

   if( blocknr == -3 )
   {
      //direct branch on copied origvar
      SCIP_CALL( branchDirectlyOnMastervar(origscip, mastervar, branchrule) );

      return SCIP_OKAY;
   }

   //calculate F and the strips
   for( i=0; i<nbranchcands; ++i )
   {
      SCIP_Bool blockfound;
      int k;

      k = 0;
      mastervar = branchcands[i];
      assert(GCGvarIsMaster(mastervar));

      blockfound = TRUE;

      if( GCGvarGetBlock(mastervar) == -1 )
      {
         assert( GCGvarIsLinking(mastervar) );
         blockfound = FALSE;

         pricingvars = GCGlinkingVarGetPricingVars(mastervar);
         assert(pricingvars != NULL );

         for( k=0; k<GCGlinkingVarGetNBlocks(mastervar); ++k )
         {
            if( pricingvars[k] != NULL )
            {
               if( GCGvarGetBlock(pricingvars[k]) == GCGbranchGenericBranchdataGetConsblocknr(branchdata) )
               {
                  blockfound = TRUE;
                  break;
               }
            }
         }
      }
      else
      {
         blockfound = (blocknr == GCGvarGetBlock(mastervar));
      }

      if( blockfound )
      {
         mastervarValue = SCIPgetSolVal(masterscip, NULL, mastervar);
         if( SCIPisGT(origscip, mastervarValue - SCIPfloor(origscip, mastervarValue), 0) )
         {

            if( Fsize == 0 )
            {
               SCIP_CALL( SCIPallocMemoryArray(origscip, &F, Fsize+1) );
            }
            else
            {
               SCIP_CALL( SCIPreallocMemoryArray(origscip, &F, Fsize+1) );
            }
            F[Fsize] = mastervar;
            ++Fsize;
         }
      }
   }

   //old data to regard?
   if( masterbranchcons != NULL && GCGconsMasterbranchGetBranchdata(masterbranchcons) != NULL )
   {
      //calculate C
      //SCIPdebugMessage("calculating C\n");
      parentcons = masterbranchcons;
      Csize = 0;
      while( parentcons != NULL && GCGconsMasterbranchGetbranchrule(parentcons) != NULL
            && strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetbranchrule(parentcons)), "generic") == 0)
      {
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
         if( branchdata == NULL )
         {
            SCIPdebugMessage("branchdata is NULL\n");
            break;
         }
         if( branchdata->consS == NULL || branchdata->consSsize == 0 )
            break;
         if( branchdata->consblocknr != blocknr )
         {
            parentcons = GCGconsMasterbranchGetParentcons(parentcons);
            continue;
         }
         if( Csize == 0)
         {
            assert(branchdata != NULL);
            assert(branchdata->consSsize > 0);
            Csize = 1;
            SCIP_CALL( SCIPallocMemoryArray(origscip, &C, Csize) );
            SCIP_CALL( SCIPallocMemoryArray(origscip, &sequencesizes, Csize) );
            C[0] = NULL;
            SCIP_CALL( SCIPallocMemoryArray(origscip, &(C[0]), branchdata->consSsize) );
            for( i=0; i<branchdata->consSsize; ++i )
            {
               C[0][i] = branchdata->consS[i];
            }
            sequencesizes[0] = branchdata->consSsize;
            //SCIPdebugMessage("sequencesizes[0]= %d\n", sequencesizes[0]);

            parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         }
         else
         {
            //S not yet in C ?
            SinC = FALSE;
            for( c=0; c<Csize && !SinC; ++c )
            {
               SinC = TRUE;
               if( branchdata->consSsize == sequencesizes[c] )
               {
                  for( i=0; i<branchdata->consSsize; ++i )
                  {
                     if( branchdata->consS[i].component != C[c][i].component || branchdata->consS[i].sense != C[c][i].sense || !SCIPisEQ(origscip, branchdata->consS[i].bound, C[c][i].bound) )
                     {
                        SinC = FALSE;
                        break;
                     }
                  }
               }
               else
                  SinC = FALSE;
            }
            if( !SinC )
            {
               ++Csize;
               SCIP_CALL( SCIPreallocMemoryArray(origscip, &C, Csize) );
               SCIP_CALL( SCIPreallocMemoryArray(origscip, &sequencesizes, Csize) );
               C[Csize-1] = NULL;
               SCIP_CALL( SCIPallocMemoryArray(origscip, &(C[Csize-1]), branchdata->consSsize) );

               /** @todo copy memory */
               for( i=0; i<branchdata->consSsize; ++i )
               {
                  C[Csize-1][i] = branchdata->consS[i];
               }
               sequencesizes[Csize-1] = branchdata->consSsize;
               //SCIPdebugMessage("sequencesizes[%d]= %.d, consSsize = %d\n", Csize-1, sequencesizes[Csize-1], branchdata->consSsize);
            }
            parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         }
      }

      if( C != NULL )
      {
         SCIPdebugMessage("Csize = %d\n", Csize);

         for( i=0; i<Csize; ++i )
         {
            for( c=0; c<sequencesizes[i]; ++c )
            {
               SCIPdebugMessage("C[%d][%d].component = %s\n", i, c, SCIPvarGetName(C[i][c].component) );
               SCIPdebugMessage("C[%d][%d].sense = %d\n", i, c, C[i][c].sense);
               SCIPdebugMessage("C[%d][%d].bound = %.6g\n", i, c, C[i][c].bound);
            }
         }

         //SCIP_CALL( InducedLexicographicSort(scip, F, Fsize, C, Csize, sequencesizes) );
         SCIP_CALL( ChooseSeparateMethod(origscip, F, Fsize, &S, &Ssize, C, Csize, sequencesizes, blocknr) );
      }
      else
      {
         SCIPdebugMessage("C == NULL\n");
         //SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
         SCIP_CALL( ChooseSeparateMethod( origscip, F, Fsize, &S, &Ssize, NULL, 0, NULL, blocknr ) );
      }
      if( sequencesizes != NULL )
      {
         assert(Csize > 0);
         SCIPfreeMemoryArray(origscip, &sequencesizes);
      }
      for( i=0; i<Csize; ++i )
         if( C[i] != NULL )
            SCIPfreeMemoryArray(origscip, &(C[i]));
      if( C != NULL )
      {
         assert( Csize > 0);
         SCIPfreeMemoryArray(origscip, &C);
      }
   }
   else
   {
      SCIPdebugMessage("root node\n");
      //SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
      SCIP_CALL( ChooseSeparateMethod( origscip, F, Fsize, &S, &Ssize, NULL, 0, NULL, blocknr ) );
   }
   assert(S!=NULL);

   if( feasible )
   {
      SCIPdebugMessage("Vanderbeck generic branching rule could not find variables to branch on!\n");
      return -1;
   }

   /* create the |S|+1 child nodes in the branch-and-bound tree */
   SCIP_CALL( createChildNodesGeneric(origscip, branchrule, S, Ssize, blocknr, masterbranchcons, result) );

   SCIPdebugMessage("free F\n");
   SCIPfreeMemoryArray(origscip, &F);

   return SCIP_OKAY;
}

//from branch_master
static
SCIP_RETCODE GCGincludeMasterCopyPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );
   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleFullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );
   SCIP_CALL( SCIPincludeBranchruleMostinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );
   SCIP_CALL( SCIPincludeBranchruleRandom(scip) );
   SCIP_CALL( SCIPincludeBranchruleRelpscost(scip) );
   return SCIP_OKAY;
}
/** copy method for master branching rule */
static
SCIP_DECL_BRANCHCOPY(branchCopyGeneric)
{
   assert(branchrule != NULL);
   assert(scip != NULL);
   SCIP_CALL( GCGincludeMasterCopyPlugins(scip) );
   return SCIP_OKAY;
}

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterGeneric)
{
   SCIP* origscip;
   SCIP_VAR** mastervars;
   SCIP_VAR** copymastervars;
   SCIP_VAR** allorigvars;
   int allnorigvars;
   int nmastervars;
   int nvarsadded;
   int nnewmastervars;
   int i;
   int p;
   char name[SCIP_MAXSTRLEN];

   i = 0;
   p = 0;
   nmastervars = 0;
   nvarsadded = 0;
   copymastervars = NULL;

   assert(scip != NULL);
   assert(branchdata != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   if( branchdata->consblocknr == -3 )
   {
      assert(branchdata->consSsize == 1);
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "directchild(%d, %g) sense = %d",
            branchdata->consSsize, branchdata->consS[0].bound, branchdata->consS[0].sense);

      // create constraint for child
      if( branchdata->consS[0].sense == GCG_COMPSENSE_GE)
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &(branchdata->mastercons), name, 0, NULL, NULL,
            branchdata->consS[0].bound, SCIPinfinity(origscip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &(branchdata->mastercons), name, 0, NULL, NULL,
            -SCIPinfinity(origscip), branchdata->consS[0].bound-1, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );
      }

      SCIP_CALL( SCIPaddCoefLinear(scip, branchdata->mastercons, branchdata->consS[0].component, 1.0) );

      /* add constraint to the master problem that enforces the branching decision */
      SCIP_CALL( SCIPaddCons(scip, branchdata->mastercons) );

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPduplicateMemoryArray(origscip, &copymastervars, mastervars, nmastervars) );

   SCIPdebugMessage("branchActiveMasterGeneric: Block %d, Ssize %d)\n", branchdata->consblocknr,
      branchdata->consSsize);

   assert( (branchdata->consSsize == 0 ) == (branchdata->consS == NULL) );

   if( branchdata->consS == NULL )
   {
      assert(branchdata->consSsize == 0);
      SCIPdebugMessage("root node:\n");
      return SCIP_OKAY;
   }

   #ifdef SCIP_DEBUG
   for( i=0; i<branchdata->consSsize; ++i )
   {
//      SCIPdebugMessage("S[%d].component = %s\n", i, SCIPvarGetName(branchdata->consS[i].component) );
//      SCIPdebugMessage("S[%d].sense = %d\n", i, branchdata->consS[i].sense);
//      SCIPdebugMessage("S[%d].bound = %.6g\n", i, branchdata->consS[i].bound);
   }
   #endif

   /* create corresponding constraint in the master problem, if not yet created */
   if( branchdata->mastercons == NULL && branchdata->consSsize > 0 )
   {

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "child(%d, %g)", branchdata->consSsize, branchdata->lhs);

      // create constraint for child
      SCIP_CALL( SCIPcreateConsLinear(scip, &(branchdata->mastercons), name, 0, NULL, NULL,
            branchdata->lhs, SCIPinfinity(origscip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );

      //add mastervars
      for( p=0; p< branchdata->consSsize; ++p )
      {
         nnewmastervars = nmastervars;
         for( i=0; i<nnewmastervars; ++i )
         {
            SCIP_Real generator_i;
            SCIP_VAR** pricingvars;
            int k;
            SCIP_Bool blockfound;

            if( i >= nmastervars )
               break;


            if( GCGvarGetBlock(copymastervars[i]) == branchdata->consblocknr
               || ( GCGvarGetBlock(copymastervars[i]) == -1 && GCGvarIsLinking(copymastervars[i])) )
            {
               blockfound = TRUE;

               if( GCGvarGetBlock(copymastervars[i]) == -1 )
               {
                  assert( GCGvarIsLinking(copymastervars[i]) );
                  blockfound = FALSE;

                  pricingvars = GCGlinkingVarGetPricingVars(copymastervars[i]);
                  assert(pricingvars != NULL );

                  for( k=0; k<GCGlinkingVarGetNBlocks(copymastervars[i]); ++k )
                  {
                     if( pricingvars[k] != NULL )
                     {
                        if( GCGvarGetBlock(pricingvars[k]) == branchdata->consblocknr )
                        {
                           blockfound = TRUE;
                           break;
                        }
                     }
                  }
               }
               if( !blockfound )
               {
                  //small down array
                  copymastervars[i] = copymastervars[nmastervars-1];
                  --i;
                  --nmastervars;
                  continue;
               }

               generator_i = getGeneratorEntry(copymastervars[i], branchdata->consS[p].component);

               if( branchdata->consS[p].sense == GCG_COMPSENSE_GE )
               {
                  if( SCIPisGE(origscip, generator_i, branchdata->consS[p].bound) )
                  {
                     if( p == branchdata->consSsize-1 )
                     {
                        //add var to constraint
                        ++nvarsadded;
                        SCIP_CALL( SCIPaddCoefLinear(scip, branchdata->mastercons, copymastervars[i], 1.0) );
                     }
                  }
                  else
                  {
                     //small down array
                     copymastervars[i] = copymastervars[nmastervars-1];
                     --i;
                     --nmastervars;
                  }
               }
               else
               {
                  if( SCIPisLT(origscip, generator_i, branchdata->consS[p].bound) )
                  {
                     if( p == branchdata->consSsize-1 )
                     {
                        //add var to constraint
                        ++nvarsadded;
                        SCIP_CALL( SCIPaddCoefLinear(scip, branchdata->mastercons, copymastervars[i], 1.0) );
                     }
                  }
                  else
                  {
                     //small down array
                     copymastervars[i] = copymastervars[nmastervars-1];
                     --i;
                     --nmastervars;
                  }
               }
            }
            else
            {
               //small down array
               copymastervars[i] = copymastervars[nmastervars-1];
               --i;
               --nmastervars;
            }
         }
      }
   }
   /* add constraint to the master problem that enforces the branching decision */
   SCIP_CALL( SCIPaddCons(scip, branchdata->mastercons) );

   SCIPdebugMessage("%d vars added with lhs= %g\n", nvarsadded, branchdata->lhs);
   assert(nvarsadded > 0);

   if(copymastervars != NULL)
      SCIPfreeMemoryArray(origscip, &copymastervars);

   return SCIP_OKAY;
}

/** callback deactivation method */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterGeneric)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->mastercons != NULL);

   SCIPdebugMessage("branchDeactiveMasterGeneric: Block %d, Ssize %d)\n", branchdata->consblocknr,
      branchdata->consSsize);

   /* remove constraint from the master problem that enforces the branching decision */
   assert(branchdata->mastercons != NULL);
   SCIP_CALL( SCIPdelCons(scip, branchdata->mastercons) );

   SCIP_CALL( SCIPreleaseCons(scip, &(branchdata->mastercons)) );
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
   SCIP* origscip;
   SCIP_Bool feasible;
   SCIP_Bool discretization;

   feasible = TRUE;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("Execrel method of Vanderbecks generic branching\n");

   *result = SCIP_DIDNOTRUN;

   /* the branching scheme only works for the discretization approach */
   SCIP_CALL( SCIPgetBoolParam(origscip, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
   {
      SCIPdebugMessage("Generic branching only for discretization approach\n");
      return SCIP_OKAY;
   }

   /* do not perform Ryan & Foster branching if we have neither a set partitioning nor a set covering structure */
   if( GCGrelaxIsMasterSetCovering(origscip) || GCGrelaxIsMasterSetPartitioning(origscip) )
   {
      SCIPdebugMessage("Generic branching executed on a set covering or set partitioning problem\n");
   }

   /* check whether the current original solution is integral */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(origscip), TRUE, TRUE, TRUE, TRUE, &feasible) );
#else
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(origscip), FALSE, TRUE, TRUE, TRUE, &feasible) );
#endif

   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(origscip, GCGrelaxGetCurrentOrigSol(origscip)));

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   *result = SCIP_BRANCHED;

   GCGbranchGenericInitbranch(scip, branchrule, result);

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextGeneric)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execext method of generic branching\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */ /*todo*/
static
SCIP_DECL_BRANCHEXECPS(branchExecpsGeneric)
{  /*lint --e{715}*/
   SCIP_CONS* masterbranchcons;
   GCG_BRANCHDATA* branchdata;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of Vanderbecks generic branching\n");

   return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   masterbranchcons = GCGconsMasterbranchGetActiveCons(scip);

   if( masterbranchcons != NULL )
      branchdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);

   if( branchdata != NULL )
   {
      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitGeneric)
{
   SCIP* origscip;

   origscip = GCGpricerGetOrigprob(scip);
   assert(branchrule != NULL);
   assert(origscip != NULL);

   SCIPdebugMessage("Init method of Vanderbecks generic branching\n");

   SCIP_CALL( GCGrelaxIncludeBranchrule(origscip, branchrule, branchActiveMasterGeneric,
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
   SCIP_CALL( SCIPincludeEventHdlrGenericbranchvaradd(scip) );

   return SCIP_OKAY;
}

/** initializes branchdata */
SCIP_RETCODE GCGbranchGenericCreateBranchdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_BRANCHDATA**      branchdata          /**< branching data to initialize */
   )
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, branchdata) );
   (*branchdata)->consS = NULL;
   (*branchdata)->consSsize = 0;
   (*branchdata)->sequencesizes = 0;
   (*branchdata)->C = NULL;
   (*branchdata)->mastercons = NULL;
   (*branchdata)->consblocknr = -2;

   return SCIP_OKAY;
}

GCG_COMPSEQUENCE* GCGbranchGenericBranchdataGetConsS(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   )
{
   assert(branchdata != NULL);
   return branchdata->consS;
}

int GCGbranchGenericBranchdataGetConsSsize(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   )
{
   assert(branchdata != NULL);
   return branchdata->consSsize;
}

int GCGbranchGenericBranchdataGetConsblocknr(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   )
{
   assert(branchdata != NULL);
   return branchdata->consblocknr;
}

SCIP_CONS* GCGbranchGenericBranchdataGetMastercons(
   GCG_BRANCHDATA*      branchdata          /**< branching data to initialize */
   )
{
   assert(branchdata != NULL);
   return branchdata->mastercons;
}
