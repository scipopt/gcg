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
//#include "scip/cons_varbound.h"

#include <stdio.h>
#include <stdlib.h>

#define BRANCHRULE_NAME          "generic"
#define BRANCHRULE_DESC          "generic branching rule by Vanderbeck"
#define BRANCHRULE_PRIORITY      999999
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

struct GCG_BranchData
{
   GCG_COMPSEQUENCE**   C;                  /**< S[k] bound sequence for block k */ //!!! sort of each C[i]=S[i] is important !!!
   int*                 sequencesizes;      /**< number of bounds in S[k] */
   int                  Csize;
   GCG_COMPSEQUENCE*    S;                  /**< component bound sequence which induce the child branching constraints */
   int                  Ssize;
   int                  blocknr;            /**< number of block branching was performed */
   //int                  childnr;
   SCIP_Real            lhs;
   int                  nchildNodes;
   SCIP_Real*           childlhs;
   SCIP_CONS*           mastercons;         /**< constraint enforcing the branching restriction in the master problem */
   //GCG_BRANCHDATA**     childbranchdatas;
   GCG_COMPSEQUENCE*    consS;              /**< component bound sequence which induce the current branching constraint */
   int                  consSsize;
   int                  consblocknr;
};

struct GCG_Strip
{
   //SCIP* scip;
   SCIP_VAR*            mastervar;          /**< master variable */
   //SCIP_Real            mastervarValue;
   //int                  blocknr;            /**< number of the block in which the strip belong */
   //SCIP_Real*           generator;          /**< corresponding generator to the mastervar */
   //SCIP_Bool*           compisinteger;      /**< ? is comp with integer origvars? */
   //int                  generatorsize;
   //GCG_COMPSEQUENCE**   C;                  /**< default NULL, only needed for ptrilocomp */
   //int                  Csize;
   //int*                 sequencesizes;
};

typedef struct GCG_Strip GCG_STRIP;

/** set of component bounds in separate */
struct GCG_Record
{
   GCG_COMPSEQUENCE**   record;             /**< returnvalue of separte function */
   int                  recordsize;
   int*                 sequencesizes;
};
typedef struct GCG_Record GCG_RECORD;

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

/** method for calculating the generator of mastervar*/
/*
SCIP_RETCODE getGenerators(
   SCIP*                scip,               //
   SCIP_Real**          generator,          // /
   int*                 generatorsize,      // /
   SCIP_Bool**          compisinteger,      // /
   int                  blocknr,            // /
   SCIP_VAR**           mastervars,         // /
   int                  nmastervars,        // /
   SCIP_VAR*            mastervar           // /
   )
{
   int i;
   int j;
   int k;
   SCIP_VAR** origvarsunion;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   i = 0;
   j = 0;
   k = 0;
   *generatorsize = 0;
   norigvars = 0;
   origvarsunion = NULL;
   assert(mastervars != NULL);
   assert(mastervar != NULL);

   for( i=0; i<nmastervars; ++i )
   {
      origvars = GCGmasterVarGetOrigvars(mastervars[i]);

      if( blocknr != GCGvarGetBlock(mastervars[i]) )
         continue;

      norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);
      // if the
      if( *generatorsize == 0 && norigvars > 0 )
      {
         *generatorsize = norigvars;
         SCIP_CALL( SCIPallocMemoryArray(scip, generator, *generatorsize) );
         SCIP_CALL( SCIPallocMemoryArray(scip, compisinteger, *generatorsize) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &origvarsunion, *generatorsize) );
         for( j=0; j<*generatorsize; ++j )
         {
            origvarsunion[j] = origvars[j];
            (*generator)[j] = 0;
            (*compisinteger)[j] = SCIPvarGetType(origvars[j]) != SCIP_VARTYPE_CONTINUOUS;// && SCIPvarGetType(origvars[j]) != SCIP_VARTYPE_IMPLINT;
         }
      }
      else
      {
         for( j=0; j<norigvars; ++j )
         {
            int oldgeneratorsize;

            oldgeneratorsize = *generatorsize;

            for( k=0; k<oldgeneratorsize; ++k )
            {
               // if variable already in union
               if( origvarsunion[k] == origvars[j] )
               {
                  break;
               }
               if( k == oldgeneratorsize-1)
               {
                  // add variable to the end
                  ++(*generatorsize);
                  SCIP_CALL( SCIPreallocMemoryArray(scip, generator, *generatorsize) );
                  SCIP_CALL( SCIPreallocMemoryArray(scip, compisinteger, *generatorsize) );
                  SCIP_CALL( SCIPreallocMemoryArray(scip, &origvarsunion, *generatorsize) );

                  origvarsunion[*generatorsize-1] = origvars[j];
                  (*generator)[*generatorsize-1] = 0;
                  (*compisinteger)[(*generatorsize)-1] = SCIPvarGetType(origvars[j]) != SCIP_VARTYPE_CONTINUOUS ;//&& SCIPvarGetType(origvars[j]) != SCIP_VARTYPE_IMPLINT;
               }
            }
         }
      }
   }

   origvars = GCGmasterVarGetOrigvars(mastervar);
   norigvars = GCGmasterVarGetNOrigvars(mastervar);
   origvals = GCGmasterVarGetOrigvals(mastervar);

   // go through all original variables and if you find the original variable, add the value to the generator
   for( i=0; i<norigvars; ++i )
   {
      for( j=0; j<*generatorsize; ++j )
      {
         assert(origvars[i] != NULL);

         if( origvarsunion[j] == origvars[i] )
         {
            if( SCIPvarGetType(origvars[i]) == SCIP_VARTYPE_CONTINUOUS ) //|| SCIPvarGetType(origvars[i]) != SCIP_VARTYPE_IMPLINT )
               (*compisinteger)[j] = FALSE;
            else
               (*compisinteger)[j] = TRUE;
            if( !SCIPisZero(scip, origvals[i]) )
               (*generator)[j] = origvals[i];
            break;
         }
      }
   }

   SCIPfreeMemoryArray(scip, &origvarsunion);

   return SCIP_OKAY;
}
*/

/** computes the generator of mastervar for the entry in origvar
 * @return entry of the generator corresponding to origvar */
SCIP_Real getGeneratorEntry(
   SCIP*                scip,               /**< SCIP data structure */
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

         generatorentry = getGeneratorEntry(scip, F[i], IndexSet[j]);
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
      SCIPdebugMessage("median = ceil(arithmMiddle), min = %f \n", min);

      for( i=0; i<arraysize; ++i )
         arithmMiddle += 1.0*array[i]/arraysize;

      Median = SCIPceil(scip, arithmMiddle);
   }
   SCIPdebugMessage("median = %f \n", Median);

   return Median;
}

/** comparefunction for lexicographical sort */
//static
//SCIP_DECL_SORTPTRCOMP(ptrcomp)
//{
//   GCG_STRIP* strip1;
//   GCG_STRIP* strip2;
//   int i;
//
//   strip1 = (GCG_STRIP*) elem1;
//   strip2 = (GCG_STRIP*) elem2;
//
//   i = 0;
//
//   assert(strip1->generatorsize == strip2->generatorsize);
//
//   for( i=0; i< strip1->generatorsize; ++i )
//   {
//      if( strip1->generator[i] > strip2->generator[i] )
//         return -1;
//      if( strip1->generator[i] < strip2->generator[i] )
//         return 1;
//   }
//
//   return 0;
//}

/** lexicographical sort using scipsort
 * This method will change the array
 */
//static
//SCIP_RETCODE LexicographicSort(
//   GCG_STRIP**          array,              /**< array to sort (will be changed) */
//   int                  arraysize           /**< size of the array */
//   )
//{
//
//   assert(array != NULL);
//   assert(arraysize > 0);
//
//   SCIPdebugMessage("Lexicographic sorting\n");
//
//   SCIPsortPtr((void**)array, ptrcomp, arraysize );
//
//   return SCIP_OKAY;
//}


/** compare function for ILO: returns 1 if bd1 < bd2 else -1 with respect to bound sequence */
//static
//int ILOcomp(
//   SCIP*                scip,               /**< SCIP data structure */
//   GCG_STRIP*           strip1,             /**< first strip */
//   GCG_STRIP*           strip2,             /**< second strip */
//   GCG_COMPSEQUENCE**   C,                  /**< component bound sequence to compare with */
//   int                  NBoundsequences,    /**< size of the bound sequence */
//   int*                 sequencesizes,      /**< sizes of the bound sequences */
//   int                  p                   /**< */
//   )
//{
//   int i;
//   int ivalue;
//   int j;
//   int k;
//   int l;
//   int Nupper;
//   int Nlower;
//   SCIP_Bool returnvalue;
//   GCG_COMPSEQUENCE** CopyC;
//   int* newsequencesizes;
//
//   i = -1;
//   j = 0;
//   k = 0;
//   l = 0;
//   Nupper = 0;
//   Nlower = 0;
//
//   /* lexicographic Order? */
//   if( C == NULL || NBoundsequences <= 1 )
//      return (*ptrcomp)( strip1, strip2);// == -1);
//
//   assert(C != NULL);
//   assert(NBoundsequences > 0);
//
//   /* find i which is in all S in C on position p */
//   while( sequencesizes[k] < p )
//   {
//      ++k;
//      assert(k < NBoundsequences);
//   }
//   i = C[k][p-1].component;
//   assert(i < strip1->generatorsize && i < strip2->generatorsize);
//   ivalue = C[k][p-1].bound;
//
//   assert(i >= 0);
//
//   /* calculate subset of C */
//   for( j=0; j< NBoundsequences; ++j )
//   {
//      if( sequencesizes[j] >= p )
//      {
//         assert(C[j][p-1].component == i);
//         if( C[j][p-1].sense == GCG_COMPSENSE_GE )
//            ++Nupper;
//         else
//            ++Nlower;
//      }
//   }
//
//   if( strip1->generator[i] >= ivalue && strip2->generator[i] >= ivalue )
//   {
//      k=0;
//      SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nupper) );
//      SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nupper) );
//      for( j=0; j< NBoundsequences; ++j )
//      {
//         if( sequencesizes[j] >= p )
//            assert(C[j][p-1].component == i);
//
//         if( sequencesizes[j] >= p && C[j][p-1].sense == GCG_COMPSENSE_GE )
//         {
//            CopyC[k] = NULL;
//            SCIP_CALL( SCIPallocMemoryArray(scip, &(CopyC[k]), sequencesizes[j]) );
//            for( l=0; l< sequencesizes[j]; ++l )
//            {
//               CopyC[k][l] = C[j][l];
//            }
//            newsequencesizes[k] = sequencesizes[j];
//            ++k;
//         }
//      }
//
//      if( k != Nupper )
//      {
//         SCIPdebugMessage("k = %d, Nupper+1 =%d\n", k, Nupper+1);
//      }
//
//      if( Nupper != 0 )
//         assert( k == Nupper );
//
//      returnvalue = ILOcomp( scip, strip1, strip2, CopyC, Nupper, newsequencesizes, p+1);
//
//      for( j=0; j< Nupper; ++j )
//      {
//
//         SCIPfreeMemoryArray(scip, &(CopyC[j]) );
//      }
//      SCIPfreeMemoryArray(scip, &newsequencesizes);
//      SCIPfreeMemoryArray(scip, &CopyC);
//
//      return returnvalue;
//   }
//
//
//   if( strip1->generator[i] < ivalue && strip2->generator[i] < ivalue )
//   {
//      k=0;
//      SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nlower) );
//      SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nlower) );
//      for( j=0; j< NBoundsequences; ++j )
//      {
//         if( sequencesizes[j] >= p )
//            assert(C[j][p-1].component == i);
//
//         if( sequencesizes[j] >= p && C[j][p-1].sense != GCG_COMPSENSE_GE )
//         {
//            CopyC[k] = NULL;
//            SCIP_CALL( SCIPallocMemoryArray(scip, &(CopyC[k]), sequencesizes[j]) );
//            for( l=0; l< sequencesizes[j]; ++l )
//            {
//               CopyC[k][l] = C[j][l];
//            }
//            newsequencesizes[k] = sequencesizes[j];
//            ++k;
//         }
//      }
//
//      if( k != Nlower )
//      {
//         SCIPdebugMessage("k = %d, Nlower+1 =%d\n", k, Nlower+1);
//      }
//
//      if( Nlower != 0 )
//         assert( k == Nlower);
//
//      returnvalue = ILOcomp( scip, strip1, strip2, CopyC, Nlower, newsequencesizes, p+1);
//
//      for( j=0; j< Nlower; ++j )
//      {
//
//         SCIPfreeMemoryArray(scip, &(CopyC[j]) );
//      }
//      SCIPfreeMemoryArray(scip, &newsequencesizes);
//      SCIPfreeMemoryArray(scip, &CopyC);
//
//      return returnvalue;
//   }
//   if( strip1->generator[i] > strip2->generator[i] )
//      return 1;
//   else
//      return -1;
//}

/** comparefunction for induced lexicographical sort */
//static
//SCIP_DECL_SORTPTRCOMP(ptrilocomp)
//{
//   GCG_STRIP* strip1;
//   GCG_STRIP* strip2;
//   int returnvalue;
//
//   strip1 = (GCG_STRIP*) elem1;
//   strip2 = (GCG_STRIP*) elem2;
//
//   returnvalue = ILOcomp( strip1->scip, strip1, strip2, strip1->C, strip1->Csize, strip1->sequencesizes, 1);
//
//   return returnvalue;
//}

/** induced lexicographical sort */
//static
//SCIP_RETCODE InducedLexicographicSort(
//   SCIP*                scip,               /**< */
//   GCG_STRIP**          array,              /**< */
//   int                  arraysize,          /**< */
//   GCG_COMPSEQUENCE**   C,                  /**< */
//   int                  NBoundsequences,    /**< */
//   int*                 sequencesizes       /**< */
//   )
//{
//   int i;
//
//   SCIPdebugMessage("Induced Lexicographic sorting\n");
//
//   if( NBoundsequences == 0 )
//      return LexicographicSort( array, arraysize );
//   assert( C!= NULL );
//
//   assert(arraysize > 0);
//
//   if( arraysize <= 1 )
//      return SCIP_OKAY;
//
//   for( i=0; i<arraysize; ++i )
//   {
//      array[i]->scip = scip;
//      array[i]->Csize = NBoundsequences;
//      array[i]->sequencesizes = sequencesizes;
//      array[i]->C = C;
//   }
//
//   SCIPsortPtr((void**)array, ptrilocomp, arraysize);
//
//   return SCIP_OKAY;
//}

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
   double*              median              /**< */
   )
{
   int j;
   int l;
   double min;
   double maxPriority;
   double* compvalues;

   SCIPdebugMessage("Jsize = %d\n", *Jsize);

   do
   {
      min = INT_MAX;
      maxPriority = INT_MIN;

      //max-min priority
      for ( j = 0; j < *Jsize; ++j )
      {
         if ( priority[j] > maxPriority && SCIPvarGetType(J[j]) != SCIP_VARTYPE_CONTINUOUS )//F[0]->compisinteger[J[j]] )
         {
            maxPriority = priority[j];
            *origvar = J[j];
         }
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
      for ( l = 0; l < Fsize; ++l )
      {
         compvalues[l] = getGeneratorEntry(scip, F[l], *origvar);//F[l]->generator[*i];
         if ( SCIPisLT(scip, compvalues[l], min) )
            min = compvalues[l];
      }
      *median = GetMedian(scip, compvalues, Fsize, min);
      SCIPfreeBufferArray(scip, &compvalues);

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
               J[j] = J[*Jsize - 1];
               break;
            }
         }
         --(*Jsize);
      }
      assert(*Jsize >= 0);

   }while ( SCIPisEQ(scip, *median, min) );

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

   SCIPdebugMessage("add to record with ");
   for( i=0; i<Ssize; ++i)
   {
      SCIPdebugPrintf(" S[%d].component =%s", i, SCIPvarGetName(record->record[record->recordsize-1][i].component) );
      SCIPdebugPrintf(" S[%d].sense =%d", i, record->record[record->recordsize-1][i].sense );
      SCIPdebugPrintf(" S[%d].bound =%g\n", i, record->record[record->recordsize-1][i].bound );
   }

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
//   int i;
   int j;
   //int n;
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
   //GCG_STRIP** copyF;
   GCG_COMPSEQUENCE* upperLowerS;
   SCIP_Real* alpha;
   SCIP_Real* compvalues;
   SCIP_Real  muF;
   SCIP_Bool found;

   assert(scip != NULL);
   assert((Fsize == 0) == (F == NULL));
   assert((IndexSetSize == 0) == (IndexSet == NULL));

 //  i = 0;
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
   copyF = NULL;
   upperLowerS = NULL;
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

   //for( j = 0; j < Fsize; ++j )
   //{
      //for( i = 0; i < IndexSetSize; ++i )
      //{
         //assert(IndexSet[i] < F[j]->generatorsize);
         max = getMaxGeneratorEntry(scip, F, Fsize, IndexSet, IndexSetSize); //MAX(F[j]->generator[IndexSet[i]], max);
      //}
   //}

   if( max == 0 )
      max = 1;

   SCIPdebugMessage("max = %g\n", max);

   for( j=0; j<Fsize; ++j )
      muF += max * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);//F[j]->mastervarValue;

   /* detect fractional alpha_i */
   SCIP_CALL( SCIPallocBufferArray(scip, &alpha, IndexSetSize) );

   for( k=0; k < IndexSetSize; ++k )
   {
      GCG_COMPSEQUENCE* copyS;
      SCIP_Real mu_F;
      SCIP_Bool even;

      even = TRUE;
      mu_F = 0;
      //i = IndexSet[k];
      origvar = IndexSet[k];
      copyS = NULL;
      alpha[k] = 0;

      //assert(i < F[0]->generatorsize);
      //SCIPdebugMessage("origvar = %s ", SCIPvarGetName(origvar));

      if( SCIPvarGetType(origvar) == SCIP_VARTYPE_CONTINUOUS ) //!F[0]->compisinteger[i] )
         continue;

      //SCIPdebugMessage("i = %d, origvar = %s ", i, SCIPvarGetName(origvar));

      /** @todo mb: if we uses sorted master and origvars arrays and use the variable as index, we can get rid of generator here */
      for( j = 0; j < Fsize; ++j )
      {
         SCIP_Real generatorentry;

         generatorentry = getGeneratorEntry(scip, F[j], origvar);
         //assert(i < F[j]->generatorsize);
         alpha[k] += generatorentry * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);//F[j]->generator[i] * F[j]->mastervarValue;
         if( !SCIPisEQ(scip, generatorentry, 0) && !SCIPisEQ(scip, SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]), 0) ) //F[j]->generator[i], 0) && !SCIPisEQ(scip, F[j]->mastervarValue, 0) )
         {
             SCIPdebugMessage("generatorentry(%s) = %g\n", SCIPvarGetName(origvar), generatorentry);// F[j]->generator[i]);
             SCIPdebugMessage("mastervarvalue = %g\n", SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]));//F[j]->mastervarValue);
         }
      }
      if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
      {
         ++Jsize;
          SCIPdebugMessage("alpha[%d] = %g\n", k, alpha[k]);
      }
      if( SCIPisGT(scip, alpha[k] - SCIPfloor(scip, alpha[k]), 0) )
      {
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
            SCIPdebugMessage("copyS[%d].component =%s\n", l, SCIPvarGetName(copyS[l].component) );
            SCIPdebugMessage("copyS[%d].sense =%d\n", l, copyS[l].sense );
            SCIPdebugMessage("copyS[%d].bound =%g\n", l, copyS[l].bound );
         }

         /* create temporary array to compute median */
         SCIP_CALL( SCIPallocBufferArray(scip, &compvalues, Fsize) );
         for( l=0; l<Fsize; ++l )
         {
            compvalues[l] = getGeneratorEntry(scip, F[l], origvar);//F[l]->generator[i];
            if( SCIPisLT(scip, compvalues[l], min) )
               min = compvalues[l];
         }

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
               if( SCIPisGE(scip, getGeneratorEntry(scip, F[l], origvar), median) )//F[l]->generator[i], median) )
                  mu_F += SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);//F[l]->mastervarValue;
            }
            ++j;

         }while( SCIPisEQ(scip, mu_F - SCIPfloor(scip, mu_F), 0) );

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

         generatorentry = getGeneratorEntry(scip, F[l], origvar);
         //assert(J[j] < F[l]->generatorsize);

         if( generatorentry > maxcomp )//F[l]->generator[J[j]] > maxcomp )
            maxcomp = generatorentry;//F[l]->generator[J[j]];

         if( generatorentry < mincomp )//F[l]->generator[J[j]] < mincomp )
            mincomp = generatorentry;//F[l]->generator[J[j]];
      }

      priority[j] = maxcomp-mincomp;
   }

   SCIPdebugMessage("Partitioning\n");
   SCIP_CALL( partition(scip, J, &Jsize, priority, F, Fsize, &origvar, &median) );

   /** @todo mb: this is a copy of S for the recursive call below */
   SCIP_CALL( SCIPallocMemoryArray(scip, &upperLowerS, Ssize+1) );
   for( l=0; l < Ssize; ++l )
   {
      upperLowerS[l] = S[l];
   }

   upperLowerS[Ssize].component = origvar;//i;
   upperLowerS[Ssize].sense = GCG_COMPSENSE_GE;
   upperLowerS[Ssize].bound = median;

   for( k=0; k<Fsize; ++k )
   {
      if( SCIPisGE(scip, getGeneratorEntry(scip, F[k], origvar), median) )//F[k]->generator[i], median) )
         ++Fupper;
      else
         ++Flower;
   }

   /* ********************************** *
    *  choose smallest partition         *
    * ********************************** */

   SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fsize) );
   j = 0;
   if( (Flower <= Fupper && Flower > 0) || Fupper == 0 )
   {
      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisLT(scip, getGeneratorEntry(scip, F[k], origvar), median) )//F[k]->generator[i], median) )
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

   assert(j < Fsize+1);
   if( Jsize == 0 && J != NULL )
   {
      SCIPfreeMemoryArray(scip, &J);
      J = NULL;
   }


   /** @todo mb: This is a simplification according to the paper, only selecting one side */
   Separate( scip, copyF, Fsize, J, Jsize, upperLowerS, Ssize+1, record );

   SCIPfreeMemoryArray(scip, &copyF);
   SCIPfreeMemoryArray(scip, &upperLowerS);
   SCIPfreeMemoryArray(scip, &priority);
   if( J != NULL )
      SCIPfreeMemoryArray(scip, &J);
   SCIPfreeMemoryArray(scip, &alpha);

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

      generatorentry = getGeneratorEntry(scip, F[j], origvar);
      //assert( i < F[j]->generatorsize);
      if ( (isense == GCG_COMPSENSE_GE && SCIPisGE(scip, generatorentry, ivalue)) ||
           (isense == GCG_COMPSENSE_LT && SCIPisLT(scip, generatorentry, ivalue)) ) //F[j]->generator[i], ivalue)) )
      {
         alpha_i += generatorentry * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);//F[j]->mastervarValue;
         if ( !SCIPisEQ(scip, SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]), 0)
            && !SCIPisEQ(scip, generatorentry, 0) )//F[j]->mastervarValue != 0 && F[j]->generator[i] != 0 )
         {
            SCIPdebugMessage("generator(%s) = %g\n", SCIPvarGetName(origvar), generatorentry);//F[j]->generator[i]);
            SCIPdebugMessage("mastervarvalue = %g\n", SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]));
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
//   int i;
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
   SCIP_VAR** copyF;
   GCG_COMPSEQUENCE* copyS;
   GCG_COMPSEQUENCE** CopyC;
   int* newsequencesizes;
   SCIP_VAR* origvar;
   SCIP_Real alpha_i;
   SCIP_Real  muF;
   SCIP_Bool found;
   SCIP_Real nu_F;
   SCIP_Real max;

   //i = 0;
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
   CopyC = NULL;
   newsequencesizes = NULL;
   copyF = NULL;
   CopyC = NULL;
   copyS = NULL;
   origvar = NULL;
   found = FALSE;

   SCIPdebugMessage("Explore\n");

   SCIPdebugMessage("with Fsize = %d, Csize = %d, Ssize = %d\n", Fsize, Csize, *Ssize);

   /* *************************************** *
    *   if C=Ø, call separate and return that *
    * *************************************** */
   if( C == NULL || Fsize==0 || IndexSetSize==0 || Csize == 0 )
   {
      SCIPdebugMessage("go to Separate\n");
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
      SCIPdebugMessage("sequencesizes[%d] = %d\n", k, sequencesizes[k]);
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
   //i = C[k][p-1].component;
   origvar = C[k][p-1].component;
   isense = C[k][p-1].sense;
   ivalue = C[k][p-1].bound;

   //assert(i < F[0]->generatorsize);
   assert(origvar != NULL);// i>= 0);
   SCIPdebugMessage("orivar = %s; ivalue = %g\n", SCIPvarGetName(origvar), ivalue);

//   for( j=0; j<Fsize; ++j )
//   {
//      for( l=0; l<IndexSetSize; ++l )
//      {
//         assert(IndexSet[l] < F[j]->generatorsize);
//         max = MAX(max, F[j]->generator[IndexSet[l]]);
//      }
//   }
   max = getMaxGeneratorEntry(scip, F, Fsize, IndexSet, IndexSetSize); //MAX(F[j]->generator[IndexSet[i]], max);

   if( max == 0 )
      max = 1;

   for( j=0; j<Fsize; ++j )
   {
      muF += max * SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[j]);//F[j]->mastervarValue;
   }

   SCIPdebugMessage("muF = %g\n", muF);

   /* ******************************************* *
    * compute alpha_i                             *
    * ******************************************* */

   alpha_i = computeAlpha(scip, Fsize, isense, ivalue, origvar, F);

   if( alpha_i == 0 && isense != GCG_COMPSENSE_GE )
   {
      isense = GCG_COMPSENSE_GE;
      alpha_i = computeAlpha(scip, Fsize, isense, ivalue, origvar, F);
   }
   SCIPdebugMessage("alpha(%s) = %g\n", SCIPvarGetName(origvar), alpha_i);

   /* ******************************************* *
    * if f > 0, add pair to record                *
    * ******************************************* */
   if( SCIPisGT(scip, alpha_i - SCIPfloor(scip, alpha_i), 0) )
   {
      found = TRUE;
      SCIPdebugMessage("fractional alpha(%s) = %g\n", SCIPvarGetName(origvar), alpha_i);

      /* ******************************************* *
       * compute nu_F                                *
       * ******************************************* */
      nu_F = 0;
      for( l = 0; l < Fsize; ++l )
      {
         if( (isense == GCG_COMPSENSE_GE && SCIPisGE(scip, getGeneratorEntry(scip, F[l], origvar), ivalue)) //F[l]->generator[i], ivalue)) ||
            || (isense == GCG_COMPSENSE_LT && SCIPisLT(scip, getGeneratorEntry(scip, F[l], origvar), ivalue)) )//F[l]->generator[i], ivalue)) )
         {
               nu_F += SCIPgetSolVal(GCGrelaxGetMasterprob(scip), NULL, F[l]);//F[l]->mastervarValue;
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
         copyS[*Ssize].component = origvar;//i;
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

   /** @todo mb: from here that's unclear to me, I suppose he explores only one path */
   /* add something to the end of S */
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
   (*S)[*Ssize-1].component = origvar;//i;
   (*S)[*Ssize-1].sense = GCG_COMPSENSE_GE;
   (*S)[*Ssize-1].bound = median;

   for( k=0; k<Fsize; ++k )
   {
      if( SCIPisGE(scip, getGeneratorEntry(scip, F[k], origvar), median) )//F[k]->generator[i], median) )
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
            ++Cupper;
         else
            ++Clower;
      }
   }

   if( SCIPisLE(scip, alpha_i, 0) && Fupper != 0 )
      Flower = INT_MAX; //Fupper
   if( SCIPisEQ(scip, alpha_i, muF) && Flower != 0 )
      Fupper = INT_MAX; //Flower

   /** @todo mb: one can easily merge these two branches as they are similar */
   //choose smallest partition
   if( ((Fupper <= Flower && Fupper > 0 ) || Flower <= 0) && Fupper != INT_MAX )
   {
      SCIPdebugMessage("chose upper bound Fupper = %d, Cupper = %d\n", Fupper, Cupper);

      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
      for( j = 0, k = 0; k < Fsize; ++k )
      {
         if( SCIPisGE(scip, getGeneratorEntry(scip, F[k], origvar), median) )//F[k]->generator[i], median) )
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
      SCIPdebugMessage("chose lower bound Flower = %d Clower = %d\n", Flower, Clower);
      (*S)[*Ssize-1].sense = GCG_COMPSENSE_LT;
      SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
      j = 0;
      for( k=0; k<Fsize; ++k )
      {
         if( SCIPisLT(scip, getGeneratorEntry(scip, F[k], origvar), median) )//F[k]->generator[i], median) )
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

   Explore( scip, CopyC, Csize, newsequencesizes, p+1, copyF, Fsize, IndexSet, IndexSetSize, S, Ssize, record );

   SCIPfreeMemoryArray(scip, &copyF);
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
 * decides whether Separate or Explore should be done
 */
static
SCIP_RETCODE ChooseSeparateMethod(
   SCIP*                scip,               /**< */
   //GCG_STRIP**          F,                  /**< */
   SCIP_VAR**           F,
   int                  Fsize,              /**< */
   GCG_COMPSEQUENCE**   S,                  /**< */
   int*                 Ssize,              /**< */
   GCG_COMPSEQUENCE**   C,                  /**< */
   int                  Csize,              /**< */
   int*                 CompSizes           /**< */
   )
{
   int i;
   SCIP_VAR** IndexSet;
   int IndexSetSize;
   GCG_RECORD* record;
   int exploreSsize;
   GCG_COMPSEQUENCE* exploreS;

   assert(Fsize > 0);
   assert(F != NULL);
   IndexSetSize = 0;
   exploreSsize = 0;
   exploreS = NULL;
   record = NULL;
   IndexSet = NULL;

   SCIPdebugMessage("Calling Separate\n");

   SCIP_CALL( SCIPallocBuffer(scip, &record) );
   record->recordsize = 0;
   record->record = NULL;
   record->sequencesizes = NULL;

   //calculate IndexSet
   //IndexSetSize = getIndexSetSize(scip, F); //F[0]->generatorsize;
   // set with origvars ??
   //SCIP_CALL( SCIPallocMemoryArray(scip, &IndexSet, IndexSetSize) );
   //for( i=0; i<IndexSetSize; ++i )
      //IndexSet[i] = i;
   SCIP_CALL( InitIndexSet(scip, F, Fsize, &IndexSet, &IndexSetSize) );
   assert(IndexSetSize > 0);
   assert(IndexSet != NULL);

   #ifdef SCIP_DEBUG
   for( i=0; i<Csize; ++i )
   {
      SCIPdebugMessage("Compsizes[%d] = %d\n", i, CompSizes[i]);
   }
   #endif

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

   ChoseS( scip, &record, S, Ssize );
   assert(*S!=NULL);

   #ifdef SCIP_DEBUG
   for( i=0; i<*Ssize; ++i )
   {
      SCIPdebugMessage("S[%d].component = %s\n", i, SCIPvarGetName((*S)[i].component) );
      SCIPdebugMessage("S[%d].sense = %d\n", i, (*S)[i].sense);
      SCIPdebugMessage("S[%d].bound = %.6g\n", i, (*S)[i].bound);
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

   if( (*branchdata)->S != NULL && (*branchdata)->Ssize > 0 )
   {
      SCIPfreeMemoryArray(scip, &((*branchdata)->S));
      (*branchdata)->S = NULL;
   }

   if( (*branchdata)->consS != NULL && (*branchdata)->consSsize > 0 )
   {
      SCIPfreeMemoryArray(scip, &((*branchdata)->consS));
      (*branchdata)->consS = NULL;
   }

   if( (*branchdata)->nchildNodes > 0 )
   {
      if( (*branchdata)->childlhs != NULL )
      {
         SCIPfreeMemoryArray(scip, &((*branchdata)->childlhs));
         (*branchdata)->childlhs = NULL;
      }
//      if((*branchdata)->childbranchdatas != NULL) // this should not happen
//      {
//         SCIPfreeMemoryArray(scip, &((*branchdata)->childbranchdatas));
//         (*branchdata)->childbranchdatas = NULL;
//      }
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

   SCIPdebugMessage("nchildren = %d\n", nchildren);

   for(i=0; i<nchildren; ++i)
   {
      SCIP_CONS* childcons;
      GCG_BRANCHDATA* branchdata;
      SCIP_Bool same;
      int j;

      SCIPdebugMessage("childnr = %d\n", i);

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
   SCIPdebugMessage("child not pruned \n");
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
   int i;

   i = 0;
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

      //cons = GCGconsMasterbranchGetParentcons(cons);  //skip the first ancestor node
      parentdata = GCGconsMasterbranchGetBranchdata(cons);
      if(parentdata == NULL)
      {
         //root node: check children for pruning
         SCIPdebugMessage("parentdata == NULL, check children\n");
         return checkchildconsS(scip, lhs, childS, childSsize, cons, childBlocknr);
      }
      if(strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetbranchrule(cons)), "generic") != 0)
         return checkchildconsS(scip, lhs, childS, childSsize, cons, childBlocknr);

      if( parentdata->blocknr == childBlocknr && parentdata->Ssize >= childSsize && parentdata->S != NULL )
      {
         for( i=0; i<childSsize; ++i )
         {
            if( parentdata->Ssize < childSsize )
               continue;

            if( !SCIPisEQ(scip, parentdata->childlhs[i], lhs) )
               continue;

            if( parentdata->S[i].component == childS[i].component && SCIPisEQ(scip, parentdata->S[i].bound, childS[i].bound) )
            {
               if( parentdata->S[i].sense != childS[i].sense && i < childSsize-1 )
                  break; //subset = FALSE;

               if( i == childSsize-1 )
               {
                  if( childSsize == parentdata->Ssize )
                  {
                     SCIPdebugMessage("child pruned 1\n");
                     return TRUE;
                  }
                  if( parentdata->S[i].sense != childS[i].sense )  //child is induced by parantdata->S
                  {
                     SCIPdebugMessage("child pruned 2\n");
                     return TRUE;
                  }
               }
            }
            else
               break; //subset = FALSE; no subset or other lhs
         }
      }
      cons = GCGconsMasterbranchGetParentcons(cons);
   }

   SCIPdebugMessage("child not pruned\n");
   return FALSE;
}

/** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_RETCODE createChildNodesGeneric(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,        /**< branching rule */
   GCG_COMPSEQUENCE*     S,                 /**< Component Bound Sequence defining the nodes */
   int                   Ssize,             /**< size of S*/
   int                   blocknr,           /**< number of the block */
   //GCG_STRIP**           F,                 /**< strips with mu>0 */
   //int                   Fsize,             /**< number of fractional strips in current column class*/
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

   assert(scip != NULL);
   assert(Ssize > 0);
   assert(S != NULL);
   //assert(F != NULL);
   //assert(Fsize > 0);

   SCIPdebugMessage("Vanderbeck branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, GCGrelaxGetNIdenticalBlocks(scip, blocknr));

   lhs = 0;
   p = 0;
   k = 0;
   i = 0;
   L = 0;
   pL = GCGrelaxGetNIdenticalBlocks(scip, blocknr);
   mu = 0;
   parentdata = NULL;


   /** @todo mb: use initialization */
   if( masterbranchcons != NULL )
   {
      parentdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);
      if( parentdata == NULL )
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
   }
   else
      parentdata = NULL;
   /** @todo mb: END extract **/

   // get variable data of the master problem
   masterscip = GCGrelaxGetMasterprob(scip);
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

      // allocate branchdata for same child and store information
      SCIP_CALL( SCIPallocMemory(scip, &branchchilddata) );
      /** @todo use initialisation ! */
      branchchilddata->consblocknr = blocknr;
      branchchilddata->mastercons = NULL;
      branchchilddata->S = NULL;
      branchchilddata->consS = NULL;
      branchchilddata->C = NULL;
      branchchilddata->sequencesizes = NULL;
      branchchilddata->childlhs = NULL;
//      branchchilddata->childbranchdatas = NULL;
      //branchchilddata->childnr = p;
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
         for( k=0;k<=p;++k )
         {
            ncopymastervars = nmastervars2;
            for( i=0; i<ncopymastervars; ++i )
            {
               SCIP_VAR* swap;
               SCIP_Real generator_i;
//               SCIP_Real* generator;
//               int generatorsize;
//               SCIP_Bool* compisinteger;

  //             generator = NULL;
  //             compisinteger = NULL;

               if( GCGvarGetBlock(mastervars2[i]) == blocknr )
               {
                  //getGenerators(scip, &generator, &generatorsize, &compisinteger, blocknr, mastervars, nmastervars, mastervars2[i]);
                  generator_i = getGeneratorEntry(scip, mastervars2[i], S[k].component); //generator[ (int) SCIPceil(scip, S[k].component-0.5)];

                  if( (S[k].sense == GCG_COMPSENSE_GE && SCIPisGE(scip, generator_i, S[k].bound) && k == p) ||
                      (S[k].sense == GCG_COMPSENSE_LT && SCIPisLT(scip, generator_i, S[k].bound) && k == p) )
                  {
                     mu += SCIPgetSolVal(masterscip, NULL, mastervars2[i]);
                  }
                  else if( ncopymastervars > 0 )
                  {
                     swap = mastervars2[i];
                     mastervars2[i] = mastervars2[ncopymastervars-1];
                     mastervars2[ncopymastervars-1] = swap;
                     --ncopymastervars;
                     --i;
                  }
               }
               else if( ncopymastervars > 0 )
               {
                  swap = mastervars2[i];
                  mastervars2[i] = mastervars2[ncopymastervars-1];
                  mastervars2[ncopymastervars-1] = swap;
                  --ncopymastervars;
                  --i;
               }
            }
            nmastervars2 = ncopymastervars;
         }
         if( p == Ssize-1 )
            L = SCIPceil(scip, mu);
         else
         {
            L = mu;
         }
         lhs = pL-L+1;
         SCIPdebugMessage("pL-L+1 = %g \n", pL-L+1);
      }
      pL = L;

      branchchilddata->lhs = lhs;
      branchchilddata->childlhs = NULL;
      branchchilddata->C = NULL;
      branchchilddata->sequencesizes = NULL;
      SCIPdebugMessage("L = %g \n", L);
      SCIPdebugMessage("lhs set to %g \n", lhs);

      /* define names for origbranch constraints */
      (void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "node(%d,%d, %f) last comp=%s", p+1, blocknr, lhs, SCIPvarGetName(branchchilddata->consS[branchchilddata->consSsize-1].component));

      if( masterbranchcons == NULL || !pruneChildNodeByDominanceGeneric(scip, lhs, branchchilddata->consS, branchchilddata->consSsize, masterbranchcons, blocknr) )
      {
         //++(*nmasternodes);
         if( masterbranchcons != NULL )
         {
            ++(parentdata->nchildNodes);
            if( parentdata->nchildNodes == 1 )
            {
               SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childlhs), parentdata->nchildNodes) );
            }
            else
            {
               SCIP_CALL( SCIPreallocMemoryArray(scip, &(parentdata->childlhs), parentdata->nchildNodes) );
            }
            parentdata->childlhs[parentdata->nchildNodes - 1] = lhs;

            for( i=0; i<branchchilddata->consSsize; ++i )
            {
               SCIPdebugMessage("setting consS[%d].component = %s\n", i, SCIPvarGetName(branchchilddata->consS[i].component) );
               SCIPdebugMessage("setting consS[%d].sense = %d\n", i, branchchilddata->consS[i].sense);
               SCIPdebugMessage("setting consS[%d].bound = %g\n", i, branchchilddata->consS[i].bound);

            }
            SCIP_CALL( SCIPcreateChild(masterscip, &child, 0.0, SCIPgetLocalTransEstimate(masterscip)) );
            SCIP_CALL( GCGcreateConsMasterbranch(masterscip, &childcons, child, GCGconsMasterbranchGetActiveCons(masterscip)) );
            SCIP_CALL( SCIPaddConsNode(masterscip, child, childcons, NULL) );

            assert(branchchilddata != NULL);
            SCIP_CALL( GCGconsMasterbranchSetOrigConsData(masterscip, childcons, childname, branchrule, branchchilddata,
               NULL, 0, FALSE, FALSE, FALSE, NULL, 0, NULL, 0) );

            // release constraints
            SCIP_CALL( SCIPreleaseCons(masterscip, &childcons) );
            SCIPdebugMessage("branchchilddata->consSsize = %d\n", branchchilddata->consSsize);
         }
      }
      else
         SCIPfreeMemory(scip, &branchchilddata);
   }

   SCIPfreeMemoryArray(scip, &mastervars2);
   SCIPfreeMemoryArray(scip, &copymastervars);

   if( parentdata->nchildNodes <= 0 )
   {
      SCIPdebugMessage("node cut off, since all childnodes have been pruned\n");

      *result = SCIP_CUTOFF;
   }

   return SCIP_OKAY;
}

/** prepares informations for using the generic branching scheme */
/** @return SCIP_RETCODE */
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
   int Ssize;
   int Csize;
   int Fsize;
   int* sequencesizes;
   int blocknr;
   int i;
   int c;
//   int k;
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
   branchdata = NULL;
   S = NULL;
   C = NULL;
   F = NULL;
   sequencesizes = NULL;

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

   //k = 0;
   masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
   SCIPdebugMessage("branching in block %d \n", blocknr);

   //calculate F and the strips
   for( i=0; i<nbranchcands; ++i )
   {
      mastervar = branchcands[i];
      assert(GCGvarIsMaster(mastervar));

      if( blocknr == GCGvarGetBlock(mastervar) )
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

            //SCIP_CALL( SCIPallocBuffer(origscip, &(F[Fsize]) ) );
            //F[Fsize]->blocknr = blocknr;
            //F[Fsize]->mastervar = mastervar;
            //F[Fsize]->mastervarValue= mastervarValue;
            F[Fsize] = mastervar;

            //getGenerators(origscip, &(F[Fsize]->generator), &(F[Fsize]->generatorsize), &(F[Fsize]->compisinteger), blocknr, mastervars, nmastervars, mastervar);
            ++Fsize;
            //++k;
         }
      }
   }

   //old data to regard?
   if( masterbranchcons != NULL && GCGconsMasterbranchGetBranchdata(masterbranchcons) != NULL )
   {
      //calculate C
      SCIPdebugMessage("calculating C\n");
      Csize = 1;
      SCIP_CALL( SCIPallocMemoryArray(origscip, &C, Csize) );
      SCIP_CALL( SCIPallocMemoryArray(origscip, &sequencesizes, Csize) );

      parentcons = masterbranchcons;

      if( parentcons != NULL && GCGconsMasterbranchGetBranchdata(parentcons)->consS != NULL && GCGconsMasterbranchGetBranchdata(parentcons)->consblocknr == blocknr )
      {
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
         assert(branchdata != NULL);
         assert(branchdata->consSsize > 0);
         C[0] = NULL;
         SCIP_CALL( SCIPallocMemoryArray(origscip, &(C[0]), branchdata->consSsize) );
         for( i=0; i<branchdata->consSsize; ++i )
         {
            C[0][i] = branchdata->consS[i];
         }
         sequencesizes[0] = branchdata->consSsize;
         SCIPdebugMessage("sequencesizes[0]= %d\n", sequencesizes[0]);

         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
      }
      else
         Csize = 0;
      while( parentcons != NULL && Csize >0 )
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
            SCIPdebugMessage("sequencesizes[%d]= %.d, consSsize = %d\n", Csize-1, sequencesizes[Csize-1], branchdata->consSsize);
         }
         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
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
         SCIP_CALL( ChooseSeparateMethod(origscip, F, Fsize, &S, &Ssize, C, Csize, sequencesizes) );
      }
      else
      {
         SCIPdebugMessage("C == NULL\n");
         //SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
         SCIP_CALL( ChooseSeparateMethod( origscip, F, Fsize, &S, &Ssize, NULL, 0, NULL ) );
      }
      SCIPfreeMemoryArray(origscip, &sequencesizes);
      for( i=0; i<Csize; ++i )
         SCIPfreeMemoryArray(origscip, &(C[i]));
      SCIPfreeMemoryArray(origscip, &C);
   }
   else
   {
      SCIPdebugMessage("root node\n");
      //SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
      SCIP_CALL( ChooseSeparateMethod( origscip, F, Fsize, &S, &Ssize, NULL, 0, NULL ) );
   }
   assert(S!=NULL);

   if( feasible )
   {
      SCIPdebugMessage("Vanderbeck generic branching rule could not find variables to branch on!\n");
      return -1;
   }

   /* create the |S|+1 child nodes in the branch-and-bound tree */
   SCIP_CALL( createChildNodesGeneric(origscip, branchrule, S, Ssize, blocknr, masterbranchcons, result) );// createChildNodesGeneric(origscip, branchrule, S, Ssize, blocknr, F, Fsize, masterbranchcons, result) );

   SCIPdebugMessage("free F\n");
  // for( i=0; i<Fsize; ++i )
  // {
   //   SCIPfreeMemoryArray(origscip, &(F[i]->compisinteger));
    //  SCIPfreeMemoryArray(origscip, &(F[i]->generator));
    //  SCIPfreeBuffer(origscip, &(F[i]));
   //}
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
   SCIPdebugMessage("pricer copy called.\n");
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
//   int oldnmastervars;
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

   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocMemoryArray(origscip, &copymastervars, nmastervars) );

   for( i=0; i<nmastervars; ++i )
      copymastervars[i] = mastervars[i];

  // oldnmastervars = nmastervars;

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
      SCIPdebugMessage("S[%d].component = %s\n", i, SCIPvarGetName(branchdata->consS[i].component) );
      SCIPdebugMessage("S[%d].sense = %d\n", i, branchdata->consS[i].sense);
      SCIPdebugMessage("S[%d].bound = %.6g\n", i, branchdata->consS[i].bound);
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
            SCIP_Real* generator;
            SCIP_Bool* compisinteger;
            int generatorsize;
            SCIP_Real generator_i;

            generatorsize = 0;
            generator = NULL;
            compisinteger = NULL;
            if( GCGvarGetBlock(copymastervars[i]) == branchdata->consblocknr )
            {
               //getGenerators(origscip, &generator, &generatorsize, &compisinteger, branchdata->consblocknr, mastervars, oldnmastervars, copymastervars[i]);

               //assert((int) SCIPceil(origscip, branchdata->consS[p].component-0.5) < generatorsize);
               generator_i = getGeneratorEntry(origscip, copymastervars[i], branchdata->consS[p].component);//generator[(int) SCIPceil(origscip, branchdata->consS[p].component-0.5)];

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
                     copymastervars[i] = copymastervars[nnewmastervars-1];
                     --i;
                     --nnewmastervars;
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
                     copymastervars[i] = copymastervars[nnewmastervars-1];
                     --i;
                     --nnewmastervars;
                  }
               }
            }
            if( generatorsize > 0 )
            {
               SCIPfreeMemoryArray(origscip, &generator);
               generator = NULL;
               SCIPfreeMemoryArray(origscip, &compisinteger);
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
   SCIP_CALL( SCIPaddCons(scip, branchdata->mastercons) );

   SCIPdebugMessage("%d vars added with lhs= %g\n", nvarsadded, branchdata->lhs);

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

/** initialize branchdata at the node */
static
SCIP_RETCODE initNodeBranchdata(
   GCG_BRANCHDATA*      nodebranchdata,     /**< branching data to set */
   GCG_BRANCHDATA*      branchdata          /**< branchdata to copy from */
   )
{
   int j;

   assert( nodebranchdata != NULL);
   assert( branchdata != NULL );

   nodebranchdata->lhs = branchdata->lhs;
   nodebranchdata->nchildNodes = 0;
   nodebranchdata->consblocknr = branchdata->consblocknr;
   nodebranchdata->mastercons = branchdata->mastercons;
   nodebranchdata->S = NULL;
   nodebranchdata->consS = NULL;
//   nodebranchdata->childbranchdatas = NULL;
   nodebranchdata->childlhs = NULL;
   nodebranchdata->C = NULL;
   nodebranchdata->sequencesizes = NULL;
   //nodebranchdata->childnr = branchdata->childnr;
   nodebranchdata->consSsize = branchdata->consSsize;

   SCIPdebugMessage("consSsize = %d, lhs = %g\n", nodebranchdata->consSsize, nodebranchdata->lhs);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(nodebranchdata->consS), nodebranchdata->consSsize) );

   for ( j = 0; j < nodebranchdata->consSsize; ++j )
   {
      nodebranchdata->consS[j].component = branchdata->consS[j].component;
      nodebranchdata->consS[j].sense = branchdata->consS[j].sense;
      nodebranchdata->consS[j].bound = branchdata->consS[j].bound;
      SCIPdebugMessage("consS[%d].component = %s\n", j, SCIPvarGetName(nodebranchdata->consS[j].component) );
      SCIPdebugMessage("consS[%d].sense = %d\n", j, nodebranchdata->consS[j].sense);
      SCIPdebugMessage("consS[%d].bound = %.6g\n", j, nodebranchdata->consS[j].bound);
   }

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpGeneric)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_Bool feasible;
   SCIP_Bool discretization;
   SCIP_CONS* masterbranchcons;
   GCG_BRANCHDATA* branchdata;
//   int i;

   feasible = TRUE;
   branchdata = NULL;

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

   masterbranchcons = GCGconsMasterbranchGetActiveCons(scip);

   if( masterbranchcons != NULL )
      branchdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);

   if( branchdata!=NULL )
   {
      SCIPdebugMessage("branchdata->nchildNodes = %d\n", branchdata->nchildNodes);

//      /**todo clear branchcilddata of parent node if its not already cleared */
//      if( branchdata->childbranchdatas != NULL )
//      {
//         for( i = 0; i < branchdata->nchildNodes; ++i )
//            SCIPfreeMemory(origscip, &(branchdata->childbranchdatas[i]));

//         if( branchdata->nchildNodes > 0 )
//            SCIPfreeMemoryArray(origscip, &(branchdata->childbranchdatas));
//         branchdata->childbranchdatas = NULL;
//      }
   }
   else
   {
      SCIPdebugMessage("Branchdata is NULL!\n");
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
//      SCIPfreeMemoryArray(origscip, &(branchdata->childbranchdatas));
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
   SCIP_CALL( SCIPincludeEventHdlrGenericbranchvaradd(scip) );  // GCGrelaxGetMasterprob(scip)) );

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
   (*branchdata)->Ssize = 0;
   (*branchdata)->S = NULL;
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
