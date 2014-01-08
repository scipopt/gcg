/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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

/**@file   decomp.c
 * @ingroup DECOMP
 * @brief  generic methods for working with different decomposition structures
 * @author Martin Bergner
 *
 * Various methods to work with the decomp structure
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "decomp.h"
#include "gcg.h"
#include "cons_decomp.h"
#include "scip/scip.h"
#include "struct_decomp.h"
#include "scip_misc.h"
#include "relax_gcg.h"

#include <assert.h>

typedef struct {
   SCIP_Real mean;
   SCIP_Real median;
   SCIP_Real max;
   SCIP_Real min;
} DEC_STATISTIC;

#define ELEM_SWAP(a,b) { register SCIP_Real t=(a);(a)=(b);(b)=t; }

static
SCIP_Real quick_select_median(SCIP_Real arr[], unsigned int n)
{
   unsigned int low, high;
   unsigned int median;

   low = 0;
   high = n - 1;
   median = (low + high) / 2;

   for( ;; )
   {
      unsigned int middle, ll, hh;
      if( high <= low ) /* One element only */
         return arr[median];

      if( high == low + 1 ) /* Two elements only */
      {
         if( arr[low] > arr[high] )
            ELEM_SWAP(arr[low], arr[high]);
         return arr[median];
      }

      /* Find median of low, middle and high items; swap into position low */
      middle = (low + high) / 2;
      if( arr[middle] > arr[high] )
         ELEM_SWAP(arr[middle], arr[high]);
      if( arr[low] > arr[high] )
         ELEM_SWAP(arr[low], arr[high]);
      if( arr[middle] > arr[low] )
         ELEM_SWAP(arr[middle], arr[low]);
      /* Swap low item (now in position middle) into position (low+1) */
      ELEM_SWAP(arr[middle], arr[low + 1]);
      /* Nibble from each end towards middle, swapping items when stuck */
      ll = low + 1;
      hh = high;
      for( ;; )
      {
         do
            ll++;
         while( arr[low] > arr[ll] );
         do
            hh--;
         while( arr[hh] > arr[low] );
         if( hh < ll )
            break;
         ELEM_SWAP(arr[ll], arr[hh]);
      }
      /* Swap middle item (in position low) back into correct position */
      ELEM_SWAP(arr[low], arr[hh]);

      /* Re-set active partition */
      if( hh <= median )
         low = ll;
      if( hh >= median )
         high = hh - 1;
   }

   return arr[median];
}


/** fill out subscipvars arrays from the information from vartoblock */
static
SCIP_RETCODE fillOutVarsFromVartoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            haslinking          /**< returns whether there are linking variables */
   )
{
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;

   SCIP_VAR** linkingvars;
   int nlinkingvars;
   int i;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);
   assert(nblocks >= 0);
   assert(vars != NULL);
   assert(nvars > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   nlinkingvars = 0;

   *haslinking = FALSE;

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) ); /*lint !e866*/
      nsubscipvars[i] = 0;
   }

   /* handle variables */
   for( i = 0; i < nvars; ++i )
   {
      int block;
      SCIP_VAR* var;

      var = vars[i];
      assert(var != NULL);
      if( !SCIPhashmapExists(vartoblock, var) )
         block = nblocks+1;
      else
      {
         block = (int)(size_t)SCIPhashmapGetImage(vartoblock, var); /*lint !e507*/
      }

      assert(block > 0 && block <= nblocks+2);

      /* if variable belongs to a block */
      if( block <= nblocks )
      {
         SCIPdebugMessage("var %s in block %d.\n", SCIPvarGetName(var), block-1);
         subscipvars[block-1][nsubscipvars[block-1]] = var;
         ++(nsubscipvars[block-1]);
      }
      else /* variable is linking or master*/
      {
         assert(block == nblocks+1 || block == nblocks+2 );

         if( block == nblocks+2 )
            SCIPdebugMessage("var %s is linking.\n", SCIPvarGetName(var));
         else
            SCIPdebugMessage("var %s is in master only.\n", SCIPvarGetName(var));

         linkingvars[nlinkingvars] = var;
         ++nlinkingvars;
      }
   }

   if( nlinkingvars > 0 )
   {
      SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars) );
      *haslinking = TRUE;
   }

   for( i = 0; i < nblocks; ++i )
   {
      if( nsubscipvars[i] == 0 )
      {
         SCIPfreeBufferArray(scip, &subscipvars[i]);
         subscipvars[i] = NULL;
      }
   }
   if( nblocks > 0 )
   {
      SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
   }
   SCIPfreeBufferArray(scip, &nsubscipvars);

   DECdecompSetVartoblock(decdecomp, vartoblock);

   for( i = 0; i < nblocks; ++i )
   {
     SCIPfreeBufferArrayNull(scip, &subscipvars[i]);
   }

   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &linkingvars);

   return SCIP_OKAY;
}


/** fill out subscipcons arrays from the information from constoblock */
static
SCIP_RETCODE fillOutConsFromConstoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            haslinking          /**< returns whether there are linking constraints */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CONS*** subscipconss;
   int* nsubscipconss;

   SCIP_CONS** linkingconss;
   int nlinkingconss;
   int i;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(constoblock != NULL);
   assert(nblocks >= 0);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(haslinking != NULL);

   DECdecompSetConstoblock(decdecomp, constoblock);

   SCIP_CALL( SCIPallocMemoryArray(scip, &linkingconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );

   *haslinking = FALSE;
   retcode = SCIP_OKAY;
   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &subscipconss[i], nconss) ); /*lint !e866*/
      nsubscipconss[i] = 0;
   }

   nlinkingconss = 0;

   /* handle constraints */
   for( i = 0; i < nconss; ++i )
   {
      int block;
      SCIP_CONS* cons;

      cons = conss[i];
      assert(cons != NULL);
      if( !SCIPhashmapExists(decdecomp->constoblock, cons) )
      {
         block = nblocks+1;
         SCIP_CALL( SCIPhashmapInsert(decdecomp->constoblock, cons, (void*) (size_t) block) );
      }
      else
      {
         block = (int)(size_t)SCIPhashmapGetImage(decdecomp->constoblock, cons); /*lint !e507*/
      }

      assert(block > 0 && block <= nblocks+1);

      /* if constraint belongs to a block */
      if( block <= nblocks )
      {
         SCIPdebugMessage("cons %s in block %d.\n", SCIPconsGetName(cons), block-1);
         subscipconss[block-1][nsubscipconss[block-1]] = cons;
         ++(nsubscipconss[block-1]);
      }
      else /* constraint is linking */
      {
         SCIPdebugMessage("cons %s is linking.\n", SCIPconsGetName(cons));
         assert(block == nblocks+1);
         linkingconss[nlinkingconss] = cons;
         ++nlinkingconss;
      }
   }

   if( nlinkingconss > 0 )
   {
      retcode = DECdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss);
      *haslinking = TRUE;
   }
   if( nblocks > 0 )
   {
      retcode = DECdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss);
   }
   SCIPfreeMemoryArray(scip, &linkingconss);
   SCIPfreeBufferArray(scip, &nsubscipconss);

   for( i = 0; i < nblocks; ++i )
   {
     SCIPfreeMemoryArray(scip, &subscipconss[i]);
   }

   SCIPfreeBufferArray(scip, &subscipconss);

   return retcode;
}


const char *DECgetStrType(
   DEC_DECTYPE type
   )
{
   const char * names[] = { "unknown", "arrowhead", "staircase", "diagonal", "bordered" };
   return names[type];
}

/** initializes the decdecomp structure to absolutely nothing */
SCIP_RETCODE DECdecompCreate(
   SCIP*                 scip,               /**< Pointer to the SCIP instance */
   DEC_DECOMP**          decomp              /**< Pointer to the decdecomp instance */
   )
{
   assert(scip != NULL);
   assert(decomp != NULL);
   SCIP_CALL( SCIPallocMemory(scip, decomp) );

   (*decomp)->type = DEC_DECTYPE_UNKNOWN;
   (*decomp)->constoblock = NULL;
   (*decomp)->vartoblock = NULL;
   (*decomp)->subscipvars = NULL;
   (*decomp)->subscipconss = NULL;
   (*decomp)->nsubscipconss = NULL;
   (*decomp)->nsubscipvars = NULL;
   (*decomp)->linkingconss = NULL;
   (*decomp)->nlinkingconss = 0;
   (*decomp)->linkingvars = NULL;
   (*decomp)->nlinkingvars = 0;
   (*decomp)->stairlinkingvars = NULL;
   (*decomp)->nstairlinkingvars = NULL;
   (*decomp)->nblocks = 0;
   (*decomp)->presolved = (SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING);
   (*decomp)->consindex = NULL;
   (*decomp)->varindex = NULL;
   (*decomp)->detector = NULL;

   return SCIP_OKAY;
}

/** frees the decdecomp structure */
SCIP_RETCODE DECdecompFree(
   SCIP*                 scip,               /**< pointer to the SCIP instance */
   DEC_DECOMP**          decdecomp           /**< pointer to the decdecomp instance */
   )
{
   DEC_DECOMP* decomp;
   int i;
   int j;

   assert( scip!= NULL );
   assert( decdecomp != NULL);
   decomp = *decdecomp;

   assert(decomp != NULL);

   for( i = 0; i < decomp->nblocks; ++i )
   {
      if( decomp->nsubscipvars != NULL )
      {
         for( j = 0; j < decomp->nsubscipvars[i]; ++j )
         {
            if( decomp->subscipvars[i][j] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(decomp->subscipvars[i][j])) );
            }
         }

         SCIPfreeMemoryArrayNull(scip, &(decomp->subscipvars[i]));
      }
      if( decomp->nsubscipconss != NULL )
      {
         for( j = 0; j < decomp->nsubscipconss[i]; ++j )
         {
            if( decomp->subscipconss[i][j] != NULL )
            {
               SCIP_CALL( SCIPreleaseCons(scip, &(decomp->subscipconss[i][j])) );
            }
         }
         SCIPfreeMemoryArrayNull(scip, &decomp->subscipconss[i]);
      }
   }

   if( decomp->linkingvars != NULL )
   {
      for( i = 0; i < decomp->nlinkingvars; ++i )
      {
         if( decomp->linkingvars[i] != NULL )
         {
            if( decomp->linkingvars[i] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(decomp->linkingvars[i])) );
            }
         }
      }
   }

   if( decomp->stairlinkingvars != NULL )
      for( i = 0; i < decomp->nblocks-1; ++i )
      {
         for( j = 0; j < decomp->nstairlinkingvars[i]; ++j )
         {
            if( decomp->stairlinkingvars[i][j] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(decomp->stairlinkingvars[i][j])) );
            }
         }
         SCIPfreeMemoryArrayNull(scip, &decomp->stairlinkingvars[i]);
      }

   /* free hashmaps if they are not NULL */
   if( decomp->constoblock != NULL )
      SCIPhashmapFree(&decomp->constoblock);
   if( decomp->vartoblock != NULL )
      SCIPhashmapFree(&decomp->vartoblock);
   if( decomp->varindex != NULL )
      SCIPhashmapFree(&decomp->varindex);
   if( decomp->consindex != NULL )
      SCIPhashmapFree(&decomp->consindex);

   for( i = 0; i < decomp->nlinkingconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(decomp->linkingconss[i])) );
   }
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->stairlinkingvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nstairlinkingvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingconss);
   SCIPfreeMemoryNull(scip, decdecomp);

   return SCIP_OKAY;
}

/** sets the type of the decomposition */
SCIP_RETCODE DECdecompSetType(
   DEC_DECOMP*           decdecomp,          /**< pointer to the decdecomp instance */
   DEC_DECTYPE           type               /**< type of the decomposition */
   )
{
   SCIP_Bool valid;
   valid = TRUE;
   assert(decdecomp != NULL);
   switch( type )
   {
   case DEC_DECTYPE_DIAGONAL:
      valid = decdecomp->nlinkingconss == 0 && decdecomp->linkingconss == NULL;
      valid = valid && decdecomp->nlinkingvars == 0 && decdecomp->linkingvars == NULL;
      break;
   case DEC_DECTYPE_ARROWHEAD:
      valid = TRUE;
      break;
   case DEC_DECTYPE_UNKNOWN:
      valid = FALSE;
      break;
   case DEC_DECTYPE_BORDERED:
      valid = decdecomp->nlinkingvars == 0 && decdecomp->linkingvars == NULL;
      break;
   case DEC_DECTYPE_STAIRCASE:
      valid = decdecomp->nlinkingconss == 0 && decdecomp->linkingconss == NULL;
      break;
   default:
      valid = FALSE;
      break;
   }

   if( !valid )
   {
      SCIPerrorMessage("The decomposition is not of the given type!\n");
      return SCIP_INVALIDDATA;
   }


   decdecomp->type = type;

   return SCIP_OKAY;
}

/** gets the type of the decomposition */
DEC_DECTYPE DECdecompGetType(
   DEC_DECOMP*           decdecomp           /**< Pointer to the decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->type;
}


/** sets the presolved flag for decomposition */
void DECdecompSetPresolved(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_Bool             presolved           /**< presolved flag for decomposition */
   )
{
   assert(decdecomp != NULL);

   decdecomp->presolved = presolved;
}

/** gets the presolved flag for decomposition */
SCIP_Bool DECdecompGetPresolved(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->presolved;
}

/** sets the number of blocks for decomposition */
void DECdecompSetNBlocks(
   DEC_DECOMP*           decdecomp,          /**< Pointer to the decdecomp instance */
   int                   nblocks             /**< number of blocks for decomposition */
   )
{
   assert(decdecomp != NULL);
   assert(nblocks >= 0);

   decdecomp->nblocks = nblocks;
}

/** gets the number of blocks for decomposition */
int DECdecompGetNBlocks(
   DEC_DECOMP*           decdecomp           /**< Pointer to the decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->nblocks;
}

/** copies the input subscipvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetSubscipvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_VAR***           subscipvars,        /**< Subscipvars array  */
   int*                  nsubscipvars        /**< number of subscipvars per block */
   )
{
   SCIP_Bool valid;
   int i;
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipvars != NULL);
   assert(nsubscipvars != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->subscipvars == NULL);
   assert(decdecomp->nsubscipvars == NULL);

   valid = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipvars, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipvars, decdecomp->nblocks) );

   assert(decdecomp->subscipvars != NULL);
   assert(decdecomp->nsubscipvars != NULL);

   BMSclearMemoryArray(decdecomp->subscipvars, decdecomp->nblocks);
   BMSclearMemoryArray(decdecomp->nsubscipvars, decdecomp->nblocks);

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      assert((subscipvars[b] == NULL) == (nsubscipvars[b] == 0));
      decdecomp->nsubscipvars[b] = nsubscipvars[b];

      if( nsubscipvars[b] < 0 )
      {
         SCIPerrorMessage("Number of variables per subproblem must be nonnegative.\n");
         valid = FALSE;
      }
      else if( nsubscipvars[b] > 0 )
      {
         assert(subscipvars[b] != NULL);
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &(decdecomp->subscipvars[b]), subscipvars[b], nsubscipvars[b]) ); /*lint !e866*/

         for( i = 0; i < nsubscipvars[b]; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, decdecomp->subscipvars[b][i]) );
         }
      }
   }

   if( !valid )
   {
      return SCIP_INVALIDDATA;
   }


   return SCIP_OKAY;
}

/** returns the subscipvars array of the given decdecomp structure */
SCIP_VAR***  DECdecompGetSubscipvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->subscipvars;
}

/** returns the nsubscipvars array of the given decdecomp structure */
int*  DECdecompGetNSubscipvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->nsubscipvars;
}

/** copies the input subscipconss array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetSubscipconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_CONS***          subscipconss,       /**< Subscipconss array  */
   int*                  nsubscipconss       /**< number of subscipconss per block */
   )
{
   SCIP_Bool valid;
   int i;
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipconss != NULL);
   assert(nsubscipconss != NULL);

   assert(decdecomp->nblocks > 0);
   assert(decdecomp->subscipconss == NULL);
   assert(decdecomp->nsubscipconss == NULL);

   valid = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipconss, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipconss, decdecomp->nblocks) );

   assert(decdecomp->subscipconss != NULL);
   assert(decdecomp->nsubscipconss != NULL);

   BMSclearMemoryArray(decdecomp->subscipconss, decdecomp->nblocks);
   BMSclearMemoryArray(decdecomp->nsubscipconss, decdecomp->nblocks);

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      if( nsubscipconss[b] <= 0 || subscipconss[b] == NULL )
      {
         SCIPerrorMessage("Block %d is empty and thus invalid. Each block needs at least one constraint.\n", b);
         valid = FALSE;
      }

      decdecomp->nsubscipconss[b] = nsubscipconss[b];

      if( nsubscipconss[b] > 0 )
      {
         assert(subscipconss[b] != NULL);
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->subscipconss[b], subscipconss[b], nsubscipconss[b]) ); /*lint !e866*/
         for( i = 0; i < nsubscipconss[b]; ++i )
         {
            SCIP_CALL( SCIPcaptureCons(scip, decdecomp->subscipconss[b][i]) );
         }
      }
   }

   if( !valid )
      return SCIP_INVALIDDATA;

   return SCIP_OKAY;
}

/** returns the subscipconss array of the given decdecomp structure */
SCIP_CONS***  DECdecompGetSubscipconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->subscipconss;
}

/** returns the nsubscipconss array of the given decdecomp structure */
int*  DECdecompGetNSubscipconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->nsubscipconss;
}

/** copies the input linkingconss array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetLinkingconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_CONS**           linkingconss,       /**< Linkingconss array  */
   int                   nlinkingconss       /**< number of linkingconss per block */
   )
{
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingconss != NULL);
   assert(nlinkingconss >= 0);

   assert(decdecomp->linkingconss == NULL);
   assert(decdecomp->nlinkingconss == 0);

   decdecomp->nlinkingconss = nlinkingconss;

   if( nlinkingconss > 0 )
   {
      int i;
      assert(linkingconss != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingconss, linkingconss, nlinkingconss) );

      for( i = 0; i < nlinkingconss; ++i )
      {
         SCIP_CALL( SCIPcaptureCons(scip, decdecomp->linkingconss[i]) );
      }
   }

   if( (linkingconss == NULL) !=  (nlinkingconss == 0) )
   {
      SCIPerrorMessage("Number of linking constraints and linking constraint array are inconsistent.\n");
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

/** returns the linkingconss array of the given decdecomp structure */
SCIP_CONS**  DECdecompGetLinkingconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->linkingconss;
}

/** returns the nlinkingconss array of the given decdecomp structure */
int  DECdecompGetNLinkingconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nlinkingconss >= 0);

   return decdecomp->nlinkingconss;
}

/** copies the input linkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetLinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_VAR**            linkingvars,        /**< Linkingvars array  */
   int                   nlinkingvars        /**< number of linkingvars per block */
   )
{
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingvars != NULL || nlinkingvars == 0);

   assert(decdecomp->linkingvars == NULL);
   assert(decdecomp->nlinkingvars == 0);

   decdecomp->nlinkingvars = nlinkingvars;

   if( nlinkingvars > 0 )
   {
      int i;
      assert(linkingvars != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingvars, linkingvars, nlinkingvars) );

      for( i = 0; i < nlinkingvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, decdecomp->linkingvars[i]) );
      }
   }

   if( (linkingvars == NULL) != (nlinkingvars == 0) )
   {
      SCIPerrorMessage("Number of linking variables and linking variable array are inconsistent.\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** returns the linkingvars array of the given decdecomp structure */
SCIP_VAR**  DECdecompGetLinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->linkingvars;
}

/** returns the nlinkingvars array of the given decdecomp structure */
int  DECdecompGetNLinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nlinkingvars >= 0);

   return decdecomp->nlinkingvars;
}

/** copies the input stairlinkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetStairlinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_VAR***           stairlinkingvars,   /**< Linkingvars array  */
   int*                  nstairlinkingvars   /**< number of linkingvars per block */
   )
{
   SCIP_Bool valid;
   int b;
   int i;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(stairlinkingvars != NULL);
   assert(nstairlinkingvars != NULL);

   assert(decdecomp->nblocks > 0);

   assert(decdecomp->stairlinkingvars == NULL);
   assert(decdecomp->nstairlinkingvars == NULL);

   valid = TRUE; /**@todo A valid check needs to be implemented */ /*lint !e774 */

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->stairlinkingvars, decdecomp->nblocks-1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nstairlinkingvars, decdecomp->nblocks-1) );

   assert(decdecomp->stairlinkingvars != NULL);
   assert(decdecomp->nstairlinkingvars != NULL);

   BMSclearMemoryArray(decdecomp->stairlinkingvars, decdecomp->nblocks-1);
   BMSclearMemoryArray(decdecomp->nstairlinkingvars, decdecomp->nblocks-1);

   for( b = 0; b < decdecomp->nblocks-1; ++b )
   {
      assert(nstairlinkingvars[b] > 0 || stairlinkingvars[b] == NULL);
      decdecomp->nstairlinkingvars[b] = nstairlinkingvars[b];
      if( stairlinkingvars[b] != NULL )
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &(decdecomp->stairlinkingvars[b]), stairlinkingvars[b], nstairlinkingvars[b]) ); /*lint !e866 */
      }
      else
      {
         decdecomp->stairlinkingvars[b] = NULL;
      }
   }

   for( b = 0; b < decdecomp->nblocks-1; ++b )
   {
      for( i = 0; i < nstairlinkingvars[b]; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, decdecomp->stairlinkingvars[b][i]) );
      }
   }
   if( !valid )
   {
      SCIPerrorMessage("The staircase linking variables are inconsistent.\n");
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

/** returns the stairlinkingvars array of the given decdecomp structure */
SCIP_VAR***  DECdecompGetStairlinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->stairlinkingvars;
}

/** returns the nstairlinkingvars array of the given decdecomp structure */
int*  DECdecompGetNStairlinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nstairlinkingvars != NULL );
   return decdecomp->nstairlinkingvars;
}

/** sets the vartoblock hashmap of the given decdecomp structure */
void  DECdecompSetVartoblock(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         vartoblock          /**< Vartoblock hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);

   decdecomp->vartoblock = vartoblock;
}

/** returns the vartoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetVartoblock(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->vartoblock;
}

/** sets the constoblock hashmap of the given decdecomp structure */
void  DECdecompSetConstoblock(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         constoblock         /**< Constoblock hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(constoblock != NULL);

   decdecomp->constoblock = constoblock;
}

/** returns the constoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetConstoblock(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->constoblock;
}

/** sets the varindex hashmap of the given decdecomp structure */
void  DECdecompSetVarindex(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         varindex            /**< Varindex hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(varindex != NULL);
   decdecomp->varindex = varindex;
}

/** returns the varindex hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetVarindex(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->varindex;
}

/** sets the consindex hashmap of the given decdecomp structure */
void  DECdecompSetConsindex(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         consindex           /**< Consindex hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(consindex != NULL);
   decdecomp->consindex = consindex;
}

/** returns the consindex hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetConsindex(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->consindex;
}

/** completely initializes decdecomp from the values of the hashmaps */
SCIP_RETCODE DECfillOutDecdecompFromHashmaps(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   )
{
   SCIP_HASHMAP* varindex;
   SCIP_HASHMAP* consindex;
   int* nsubscipconss;
   int* nsubscipvars;
   int* nstairlinkingvars;
   SCIP_VAR*** stairlinkingvars;
   SCIP_CONS*** subscipconss;
   SCIP_Bool success;
   int idx;
   int cindex;
   int cumindex;
   SCIP_Bool haslinking;
   int i;
   int b;
   SCIP_VAR** curvars;
   int ncurvars;
   int j;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);
   assert(constoblock != NULL);
   assert(nblocks >= 0);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(conss != NULL);
   assert(nconss > 0);

   DECdecompSetNBlocks(decdecomp, nblocks);

   SCIP_CALL( DECdecompSetType(decdecomp, DEC_DECTYPE_DIAGONAL) );
   SCIP_CALL_QUIET( fillOutConsFromConstoblock(scip, decdecomp, constoblock, nblocks, conss, nconss, &haslinking) );

   if( haslinking )
   {
      SCIPdebugMessage("Decomposition has linking constraints and is bordered.\n");
      SCIP_CALL( DECdecompSetType(decdecomp, DEC_DECTYPE_BORDERED) );
   }

   SCIP_CALL( fillOutVarsFromVartoblock(scip,  decdecomp, vartoblock, nblocks, vars, nvars, &haslinking) );

   if( haslinking )
   {
      SCIPdebugMessage("Decomposition has linking variables and is arrowhead.\n");
      SCIP_CALL( DECdecompSetType(decdecomp, DEC_DECTYPE_ARROWHEAD) );
   }

   if( !staircase )
   {
      SCIP_CALL( DECdecompCheckConsistency(scip, decdecomp) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&consindex, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &stairlinkingvars, nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nstairlinkingvars, nblocks) );

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(stairlinkingvars[i]), nvars) ); /*lint !e866*/
      nstairlinkingvars[i] = 0;
   }

   nsubscipconss = DECdecompGetNSubscipconss(decdecomp);
   subscipconss = DECdecompGetSubscipconss(decdecomp);
   nsubscipvars = DECdecompGetNSubscipvars(decdecomp);

   idx = 0;
   cindex = 0;
   cumindex = 0;

   /* try to deduce staircase map */
   for( b = 0; b < nblocks; ++b )
   {
      SCIPdebugMessage("block %d (%d vars):\n", b, nsubscipvars[b]);
      idx = 0;

      for( i = 0; i < nsubscipconss[b]; ++i )
      {

         int linkindex = 0;
         SCIP_CONS* cons = subscipconss[b][i];

         SCIP_CALL( SCIPhashmapInsert(consindex, cons, (void*)(size_t)(cindex+1)) );
         ++cindex;
         SCIP_CALL( SCIPgetConsNVars(scip, cons, &ncurvars, &success) );
         assert(success);

         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );

         SCIP_CALL( SCIPgetConsVars(scip, cons, curvars, ncurvars, &success) );
         assert(success);

         for( j = 0; j < ncurvars; ++j )
         {
            SCIP_VAR* probvar = SCIPvarGetProbvar(curvars[j]);

            /* if the variable is linking */
            if( (int)(size_t)SCIPhashmapGetImage(vartoblock, probvar) >= nblocks+1 ) /*lint !e507*/
            {
               /* if it has not been already assigned, it links to the next block */
               if( !SCIPhashmapExists(varindex, probvar) )
               {
                  int vindex = cumindex+nsubscipvars[b]+linkindex+1;
                  SCIPdebugMessage("assigning link var <%s> to index <%d>\n", SCIPvarGetName(probvar), vindex);
                  SCIP_CALL( SCIPhashmapInsert(varindex, probvar, (void*)(size_t)(vindex)) );
                  stairlinkingvars[b][nstairlinkingvars[b]] = probvar;
                  ++(nstairlinkingvars[b]);
                  linkindex++;
               }
            }
            else
            {
               if( !SCIPhashmapExists(varindex, probvar) )
               {
                  int vindex = cumindex+idx+1;
                  assert(((int) (size_t) SCIPhashmapGetImage(vartoblock, probvar)) -1 == b);  /*lint !e507*/
                  SCIPdebugMessage("assigning block var <%s> to index <%d>\n", SCIPvarGetName(probvar), vindex);
                  SCIP_CALL( SCIPhashmapInsert(varindex, probvar, (void*)(size_t)(vindex)) );
                  ++idx;
               }
            }
         }
         SCIPfreeBufferArray(scip, &curvars);
      }
      if( b < nblocks-1 )
      {
         cumindex += nsubscipvars[b] + nstairlinkingvars[b];
      }
      idx += cumindex;
   }
   DECdecompSetVarindex(decdecomp, varindex);
   DECdecompSetConsindex(decdecomp, consindex);
   SCIP_CALL( DECdecompSetType(decdecomp, DEC_DECTYPE_STAIRCASE) );

   for( b = 0; b < nblocks; ++b )
   {
      if( nstairlinkingvars[b] == 0 )
      {
         SCIPfreeMemoryArrayNull(scip, &(stairlinkingvars[b]));
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(stairlinkingvars[b]), nstairlinkingvars[b]) ); /*lint !e866*/
      }
   }

   SCIP_CALL( DECdecompSetStairlinkingvars(scip, decdecomp, stairlinkingvars, nstairlinkingvars) );

   for( b = 0; b < nblocks; ++b )
   {
      SCIPfreeMemoryArrayNull(scip, &stairlinkingvars[b]);
   }
   SCIPfreeMemoryArray(scip, &stairlinkingvars);
   SCIPfreeMemoryArray(scip, &nstairlinkingvars);

   SCIP_CALL( DECdecompCheckConsistency(scip, decdecomp) );

   return SCIP_OKAY;
}

/** completely fills out detector structure from only the constraint partition */
SCIP_RETCODE DECfilloutDecdecompFromConstoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   )
{
   SCIP_HASHMAP* vartoblock;
   int i;
   int j;

   SCIP_VAR** curvars;
   int ncurvars;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(constoblock != NULL);
   assert(nblocks >= 0);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(conss != NULL);

   assert(nconss > 0);

   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nvars) );
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      consblock = (int)(size_t)SCIPhashmapGetImage(constoblock, conss[i]);  /*lint !e507*/

      assert(consblock > 0 && consblock <= nblocks+1);
      if( consblock == nblocks+1 )
      {
         SCIPdebugMessage("cons <%s> is linking and need not be handled\n", SCIPconsGetName(conss[i]));
         continue;
      }
      SCIP_CALL( SCIPgetConsNVars(scip, conss[i], &ncurvars, &success) );
      assert(success);

      SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );

      SCIP_CALL( SCIPgetConsVars(scip, conss[i], curvars, ncurvars, &success) );
      assert(success);
      SCIPdebugMessage("cons <%s> (%d vars) is in block %d.\n", SCIPconsGetName(conss[i]), ncurvars, consblock);

      for( j = 0; j < ncurvars; ++j )
      {
         int varblock;
         SCIP_VAR* probvar = SCIPvarGetProbvar(curvars[j]);
         assert( SCIPvarIsActive(probvar) );
         if( SCIPhashmapExists(vartoblock, probvar) )
            varblock = (int) (size_t) SCIPhashmapGetImage(vartoblock, probvar); /*lint !e507*/
         else
            varblock = nblocks+1;
         /** if the constraint is in a block and the variable is not in the same block */
         if( !SCIPhashmapExists(vartoblock, probvar) && consblock <= nblocks )
         {
            SCIPdebugMessage(" var <%s> not been handled before, adding to block %d\n", SCIPvarGetName(probvar), consblock);
            SCIP_CALL( SCIPhashmapSetImage(vartoblock, probvar, (void*) (size_t) consblock) );
         }
         else if( varblock != consblock && consblock <= nblocks )
         {
            SCIPdebugMessage(" var <%s> has been handled before, adding to linking (%d != %d)\n", SCIPvarGetName(probvar), consblock, varblock);
            SCIP_CALL( SCIPhashmapSetImage(vartoblock, probvar, (void*) (size_t) (nblocks+2)) );
         }
         else if( consblock == nblocks+1 )
         {
            SCIPdebugMessage(" var <%s> not handled and current cons linking, var is master.\n", SCIPvarGetName(probvar));
            SCIP_CALL( SCIPhashmapSetImage(vartoblock, probvar, (void*) (size_t) (nblocks+1)) );
         }
         else
         {
            assert(consblock == varblock);
            SCIPdebugMessage(" var <%s> is handled and in same block as cons (%d == %d).\n", SCIPvarGetName(probvar), consblock, varblock);
         }
      }

      SCIPfreeBufferArray(scip, &curvars);
   }

   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPhashmapExists(vartoblock, vars[i]) )
      {
         SCIPdebugMessage(" var <%s> not handled at all and now in master\n", SCIPvarGetName(vars[i]));
         SCIP_CALL( SCIPhashmapSetImage(vartoblock, vars[i], (void*) (size_t) (nblocks+1)) );
      }
   }

   SCIP_CALL_QUIET( DECfillOutDecdecompFromHashmaps(scip, decdecomp, vartoblock, constoblock, nblocks, vars, nvars, conss, nconss, staircase) );

   return SCIP_OKAY;
}

/** sets the detector for the given decdecomp structure */
void DECdecompSetDetector(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(decdecomp != NULL);

   decdecomp->detector = detector;
}

/** gets the detector for the given decdecomp structure */
DEC_DETECTOR* DECdecompGetDetector(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->detector;
}

/** transforms all constraints and variables, updating the arrays */
SCIP_RETCODE DECdecompTransform(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   )
{
   int b;
   int c;
   int v;
   SCIP_HASHMAP* newconstoblock;
   SCIP_HASHMAP* newvartoblock;
   SCIP_VAR* newvar;

   assert(SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED);

   SCIP_CALL( SCIPhashmapCreate(&newconstoblock, SCIPblkmem(scip), SCIPgetNConss(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&newvartoblock, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* transform all constraints and put them into constoblock */
   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      for( c = 0; c < decdecomp->nsubscipconss[b]; ++c )
      {
         SCIP_CONS* newcons;
         SCIPdebugMessage("%d, %d: %s (%s)\n", b, c, SCIPconsGetName(decdecomp->subscipconss[b][c]), SCIPconsIsTransformed(decdecomp->subscipconss[b][c])?"t":"o" );
         assert(decdecomp->subscipconss[b][c] != NULL);
         newcons = SCIPfindCons(scip, SCIPconsGetName(decdecomp->subscipconss[b][c]));
         if( newcons != decdecomp->subscipconss[b][c] )
         {
            SCIP_CALL( SCIPcaptureCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &(decdecomp->subscipconss[b][c])) );
            decdecomp->subscipconss[b][c] = newcons;
         }
         assert(decdecomp->subscipconss[b][c] != NULL);
         assert(!SCIPhashmapExists(newconstoblock, decdecomp->subscipconss[b][c]));
         SCIP_CALL( SCIPhashmapSetImage(newconstoblock, decdecomp->subscipconss[b][c], (void*) (size_t) (b+1)) );
      }
   }
   /* transform all variables and put them into vartoblock */
   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      int idx;
      for( v = 0, idx = 0; v < decdecomp->nsubscipvars[b]; ++v )
      {
         assert(decdecomp->subscipvars[b][v] != NULL);

         SCIPdebugMessage("%d, %d: %s (%p, %s)\n", b, v, SCIPvarGetName(decdecomp->subscipvars[b][v]),
            (void*)decdecomp->subscipvars[b][v], SCIPvarIsTransformed(decdecomp->subscipvars[b][v])?"t":"o" );

         /* make sure that newvar is a transformed variable */
         SCIP_CALL( SCIPgetTransformedVar(scip, decdecomp->subscipvars[b][v], &newvar) );
         SCIP_CALL( SCIPreleaseVar(scip, &(decdecomp->subscipvars[b][v])) );

         assert(newvar != NULL);
         assert(SCIPvarIsTransformed(newvar));

         newvar = SCIPvarGetProbvar(newvar);
         assert(newvar != NULL);

         /* the probvar can also be fixed, in which case we do not need it in the block; furthermore, multiple variables
          * can resolve to the same active problem variable, so we check whether we already handled the variable
          * @todo: why do we ignore fixed variables? They could still be present in constraints?
          */
         if( SCIPvarIsActive(newvar) && !SCIPhashmapExists(newvartoblock, newvar) )
         {
            decdecomp->subscipvars[b][idx] = newvar;
            SCIP_CALL( SCIPcaptureVar(scip, newvar) );
            SCIPdebugMessage("%d, %d: %s (%p, %s)\n", b, v, SCIPvarGetName(decdecomp->subscipvars[b][idx]),
               (void*)decdecomp->subscipvars[b][idx], SCIPvarIsTransformed(decdecomp->subscipvars[b][idx])?"t":"o" );

            assert(decdecomp->subscipvars[b][idx] != NULL);
            assert(!SCIPhashmapExists(newvartoblock, decdecomp->subscipvars[b][idx]));
            SCIP_CALL( SCIPhashmapSetImage(newvartoblock, decdecomp->subscipvars[b][idx], (void*) (size_t) (b+1)) );
            ++idx;
         }
      }
      decdecomp->nsubscipvars[b] = idx;
   }

   /* transform all linking constraints */
   for( c = 0; c < decdecomp->nlinkingconss; ++c )
   {
      SCIP_CONS* newcons;

      SCIPdebugMessage("m, %d: %s (%s)\n", c, SCIPconsGetName(decdecomp->linkingconss[c]), SCIPconsIsTransformed(decdecomp->linkingconss[c])?"t":"o" );
      assert(decdecomp->linkingconss[c] != NULL);
      newcons = SCIPfindCons(scip, SCIPconsGetName(decdecomp->linkingconss[c]));
      if( newcons != decdecomp->linkingconss[c] )
      {
         SCIP_CALL( SCIPcaptureCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &(decdecomp->linkingconss[c])) );
         decdecomp->linkingconss[c] = newcons;
      }
      SCIP_CALL( SCIPhashmapSetImage(newconstoblock, decdecomp->linkingconss[c],(void*) (size_t) (decdecomp->nblocks+1) ) );

      assert(decdecomp->linkingconss[c] != NULL);
   }

   /* transform all linking variables */
   for( v = 0; v < decdecomp->nlinkingvars; ++v )
   {
      SCIPdebugMessage("m, %d: %s (%p, %s)\n", v, SCIPvarGetName(decdecomp->linkingvars[v]),
         (void*)decdecomp->linkingvars[v], SCIPvarIsTransformed(decdecomp->linkingvars[v])?"t":"o");
      assert(decdecomp->linkingvars[v] != NULL);

      if( !SCIPvarIsTransformed(decdecomp->linkingvars[v]) )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, decdecomp->linkingvars[v], &newvar) );
         newvar = SCIPvarGetProbvar(newvar);
      }
      else
         newvar = decdecomp->linkingvars[v];
      assert(newvar != NULL);
      assert(SCIPvarIsTransformed(newvar));
      SCIP_CALL( SCIPreleaseVar(scip, &(decdecomp->linkingvars[v])) );

      decdecomp->linkingvars[v] = newvar;
      SCIP_CALL( SCIPcaptureVar(scip, decdecomp->linkingvars[v]) );
      SCIP_CALL( SCIPhashmapSetImage(newvartoblock, decdecomp->linkingvars[v], (void*) (size_t) (decdecomp->nblocks+1) ) );
      SCIPdebugMessage("m, %d: %s (%p, %s)\n", v, SCIPvarGetName(decdecomp->linkingvars[v]),
         (void*)decdecomp->linkingvars[v], SCIPvarIsTransformed(decdecomp->linkingvars[v])?"t":"o");
      assert(decdecomp->linkingvars[v] != NULL);
   }

   SCIPhashmapFree(&decdecomp->constoblock);
   decdecomp->constoblock = newconstoblock;
   SCIPhashmapFree(&decdecomp->vartoblock);
   decdecomp->vartoblock = newvartoblock;

   SCIP_CALL( DECdecompCheckConsistency(scip, decdecomp) );

   return SCIP_OKAY;
}

/**
 * Adds all those constraints that were added to the problem after the decomposition as created
 */
SCIP_RETCODE DECdecompAddRemainingConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decomposition data structure */
   )
{
   int c;

   assert(scip != NULL);
   assert(decdecomp != NULL);

   for( c = 0; c < SCIPgetNConss(scip); ++c )
   {
      SCIP_CONS * cons = SCIPgetConss(scip)[c];

      if( !GCGisConsGCGCons(cons) )
      {
         if( !SCIPhashmapExists(DECdecompGetConstoblock(decdecomp), cons) )
         {
            int block;
            SCIP_CALL( DECdetermineConsBlock(scip, decdecomp, cons, &block) );
            SCIPdebugMessage("cons <%s> in block %d/%d\n", SCIPconsGetName(cons), block, DECdecompGetNBlocks(decdecomp) );
            if( block == DECdecompGetNBlocks(decdecomp) )
            {
               SCIP_CALL( SCIPreallocMemoryArray(scip, &decdecomp->linkingconss, decdecomp->nlinkingconss+1) );
               decdecomp->linkingconss[decdecomp->nlinkingconss] = cons;
               decdecomp->nlinkingconss += 1;
               SCIP_CALL( SCIPhashmapInsert(decdecomp->constoblock, cons, (void*) (size_t) (DECdecompGetNBlocks(decdecomp)+1)) );
            }
            else
            {
               SCIP_CALL( SCIPreallocMemoryArray(scip, &decdecomp->subscipconss[block], decdecomp->nsubscipconss[block]+1) );
               decdecomp->subscipconss[block][decdecomp->nsubscipconss[block]] = cons;
               decdecomp->nsubscipconss[block] += 1;
               SCIP_CALL( SCIPhashmapInsert(decdecomp->constoblock, cons, (void*) (size_t) (block+1)) );
            }
            SCIP_CALL( SCIPcaptureCons(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}

/** checks the consistency of the data structure
 *
 *  In particular, it checks whether the redundant information in the structure agree and
 *  whether the variables in the structure are both existant in the arrays and in the problem
 */
SCIP_RETCODE DECdecompCheckConsistency(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decomposition data structure */
   )
{
#ifndef NDEBUG
   int c;
   int b;
   int v;


   SCIPdebugMessage("Problem is %stransformed\n", SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED ? "": "not ");

   for( v = 0; v < SCIPgetNVars(scip); ++v )
   {
      assert(SCIPhashmapExists(DECdecompGetVartoblock(decdecomp), SCIPgetVars(scip)[v]));
   }

   for( c = 0; c < SCIPgetNConss(scip); ++c )
   {
      if( !GCGisConsGCGCons(SCIPgetConss(scip)[c]) )
      {
         assert(SCIPhashmapExists(DECdecompGetConstoblock(decdecomp), SCIPgetConss(scip)[c]));
      }
   }

   /* Check whether subscipcons are correct */
   for( b = 0; b < DECdecompGetNBlocks(decdecomp); ++b )
   {
      for( c = 0; c < DECdecompGetNSubscipconss(decdecomp)[b]; ++c )
      {
         SCIP_VAR** curvars;
         int ncurvars;
         SCIP_CONS* cons = DECdecompGetSubscipconss(decdecomp)[b][c];
         SCIPdebugMessage("Cons <%s> in block %d = %d\n", SCIPconsGetName(cons), b, ((int) (size_t) SCIPhashmapGetImage(DECdecompGetConstoblock(decdecomp), cons)) -1);  /*lint !e507*/
         assert(SCIPfindCons(scip, SCIPconsGetName(cons)) != NULL);
         assert(((int) (size_t) SCIPhashmapGetImage(DECdecompGetConstoblock(decdecomp), cons)) -1 == b); /*lint !e507*/
         ncurvars = SCIPgetNVarsXXX(scip, cons);
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, curvars, ncurvars) );

         for( v = 0; v < ncurvars; ++v )
         {
            int varblock;
            SCIP_VAR* var = SCIPvarGetProbvar(curvars[v]);
            varblock = ((int) (size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), var)) -1;  /*lint !e507*/
            SCIPdebugMessage("\tVar <%s> in block %d = %d\n", SCIPvarGetName(var), b, varblock);
            assert(SCIPfindVar(scip, SCIPvarGetName(var)) != NULL);
            assert(SCIPvarIsActive(var));
            assert(varblock == b || varblock == DECdecompGetNBlocks(decdecomp)+1 );
         }
         SCIPfreeMemoryArray(scip, &curvars);
      }

      for( v = 0; v < DECdecompGetNSubscipvars(decdecomp)[b]; ++v )
      {
         int varblock;
         SCIP_VAR* var = DECdecompGetSubscipvars(decdecomp)[b][v];
         varblock = ((int) (size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), var)) -1; /*lint !e507*/
         SCIPdebugMessage("Var <%s> in block %d = %d\n", SCIPvarGetName(var), b, varblock);
         assert(SCIPfindVar(scip, SCIPvarGetName(var)) != NULL);
         assert(SCIPvarIsActive(var));
         assert(varblock == b || varblock == DECdecompGetNBlocks(decdecomp)+1);
      }
   }

   /* check linking constraints and variables */
   for( v = 0; v < DECdecompGetNLinkingvars(decdecomp); ++v )
   {
      int varblock;
      varblock = (int) (size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), DECdecompGetLinkingvars(decdecomp)[v]);
      assert( varblock == DECdecompGetNBlocks(decdecomp) +1 || varblock == DECdecompGetNBlocks(decdecomp)+2); /*lint !e507*/
   }
   for (c = 0; c < DECdecompGetNLinkingconss(decdecomp); ++c)
   {
      assert(((int) (size_t) SCIPhashmapGetImage(DECdecompGetConstoblock(decdecomp), DECdecompGetLinkingconss(decdecomp)[c])) -1 ==  DECdecompGetNBlocks(decdecomp)); /*lint !e507*/
   }

   switch( DECdecompGetType(decdecomp) )
   {
   case DEC_DECTYPE_UNKNOWN:
         assert(FALSE);
      break;
   case DEC_DECTYPE_ARROWHEAD:
      assert(DECdecompGetNLinkingvars(decdecomp) > 0);
      break;
   case DEC_DECTYPE_BORDERED:
      assert(DECdecompGetNLinkingvars(decdecomp) == 0 && DECdecompGetNLinkingconss(decdecomp) > 0);
      break;
   case DEC_DECTYPE_DIAGONAL:
      assert(DECdecompGetNLinkingvars(decdecomp) == 0 && DECdecompGetNLinkingconss(decdecomp) == 0);
      break;
   case DEC_DECTYPE_STAIRCASE:
      assert(DECdecompGetNLinkingvars(decdecomp) > 0 && DECdecompGetNLinkingconss(decdecomp) == 0);
      break;
   default:
         assert(FALSE);
         break;
   }

#endif
   return SCIP_OKAY;
}

/** creates a decomposition with all constraints in the master */
SCIP_RETCODE DECcreateBasicDecomp(
   SCIP*                 scip,                /**< SCIP data structure */
   DEC_DECOMP**          decomp               /**< decomposition structure */
   )
{
   SCIP_HASHMAP* constoblock;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);

   SCIP_CALL( DECdecompCreate(scip, decomp) );
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );

   for( c = 0; c < nconss; ++c )
   {
      if( GCGisConsGCGCons(conss[c]) )
         continue;

      SCIP_CALL( SCIPhashmapInsert(constoblock, conss[c], (void*) (size_t) 1 ) );
   }

   SCIP_CALL( DECfilloutDecdecompFromConstoblock(scip, *decomp, constoblock, 0, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   return SCIP_OKAY;
}

/**
 * processes block representatives
 *
 * @return returns the number of blocks
 */
static
int processBlockRepresentatives(
   int                   maxblock,           /**< maximal number of blocks */
   int*                  blockrepresentative /**< array blockrepresentatives */
   )
{
   int i;
   int tempblock = 1;

   assert(maxblock >= 1);
   assert(blockrepresentative != NULL );
   SCIPdebugPrintf("Blocks: ");

   /* postprocess blockrepresentatives */
   for (i = 1; i < maxblock; ++i)
   {
      /* forward replace the representatives */
      assert(blockrepresentative[i] >= 0);
      assert(blockrepresentative[i] < maxblock);
      if (blockrepresentative[i] != i)
         blockrepresentative[i] = blockrepresentative[blockrepresentative[i]];
      else
      {
         blockrepresentative[i] = tempblock;
         ++tempblock;
      }
      /* It is crucial that this condition holds */
      assert(blockrepresentative[i] <= i);
      SCIPdebugPrintf("%d ", blockrepresentative[i]);
   }
   SCIPdebugPrintf("\n");
   return tempblock-1;
}

/** */
static
SCIP_RETCODE assignConstraintsToRepresentatives(
   SCIP*                 scip,               /**< */
   SCIP_CONS**           conss,              /**< */
   int                   nconss,             /**< */
   SCIP_Bool*            consismaster,       /**< */
   SCIP_HASHMAP*         constoblock,        /**< */
   int*                  vartoblock,         /**< */
   int*                  nextblock,          /**< */
   int*                  blockrepresentative /**< */
   )
{

   int i;
   int j;
   int k;
   SCIP_VAR** curvars;
   int ncurvars;

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   /* go through the all constraints */
   for (i = 0; i < nconss; ++i)
   {
      int consblock;
      SCIP_CONS* cons = conss[i];

      assert(cons != NULL);
      if (GCGisConsGCGCons(cons))
         continue;

      if (consismaster[i])
         continue;

      /* get variables of constraint */
      ncurvars = SCIPgetNVarsXXX(scip, cons);
      curvars = NULL;
      if (ncurvars > 0)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, curvars, ncurvars) );
      }
      assert(ncurvars >= 0);
      assert(ncurvars <= SCIPgetNVars(scip));
      assert(curvars != NULL || ncurvars == 0);

      assert(SCIPhashmapGetImage(constoblock, cons) == NULL);

      /* if there are no variables, put it in the first block, otherwise put it in the next block */
      if (ncurvars == 0)
         consblock = -1;
      else
         consblock = *nextblock;

      /* go through all variables */
      for (j = 0; j < ncurvars; ++j)
      {
         SCIP_VAR* probvar;
         int varindex;
         int varblock;

         assert(curvars != NULL);
         probvar = SCIPvarGetProbvar(curvars[j]);
         assert(probvar != NULL);

         varindex = SCIPvarGetProbindex(probvar);
         assert(varindex >= 0);
         assert(varindex < SCIPgetNVars(scip));

         /** @todo what about deleted variables? */
         /* get block of variable */
         varblock = vartoblock[varindex];

         SCIPdebugMessage("\tVar %s (%d): ", SCIPvarGetName(probvar), varblock);
         /* if variable is assigned to a block, assign constraint to that block */
         if (varblock > -1 && varblock != consblock)
         {
            consblock = MIN(consblock, blockrepresentative[varblock]);
            SCIPdebugPrintf("still in block %d.\n", varblock);
         }
         else if (varblock == -1)
         {
            /* if variable is free, assign it to the new block for this constraint */
            varblock = consblock;
            assert(varblock > 0);
            assert(varblock <= *nextblock);
            vartoblock[varindex] = varblock;
            SCIPdebugPrintf("new in block %d.\n", varblock);
         } else
         {
            assert((varblock > 0) && (consblock == varblock));
            SCIPdebugPrintf("no change.\n");
         }

         SCIPdebugPrintf("VARINDEX: %d (%d)\n", varindex, vartoblock[varindex]);
      }

      /* if the constraint belongs to a new block, mark it as such */
      if (consblock == *nextblock)
      {
         assert(consblock > 0);
         blockrepresentative[consblock] = consblock;
         assert(blockrepresentative[consblock] > 0);
         assert(blockrepresentative[consblock] <= *nextblock);
         ++(*nextblock);
      }

      SCIPdebugMessage("Cons %s will be in block %d (next %d)\n", SCIPconsGetName(cons), consblock, *nextblock);

      for (k = 0; k < ncurvars; ++k)
      {
         int curvarindex;
         SCIP_VAR* curprobvar;
         int oldblock;
         assert(curvars != NULL);
         curprobvar = SCIPvarGetProbvar(curvars[k]);
         curvarindex = SCIPvarGetProbindex(curprobvar);

         oldblock = vartoblock[curvarindex];
         assert((oldblock > 0) && (oldblock <= *nextblock));
         SCIPdebugMessage("\tVar %s ", SCIPvarGetName(curprobvar));
         if (oldblock != consblock)
         {
            SCIPdebugPrintf("reset from %d to block %d.\n", oldblock, consblock);
            vartoblock[curvarindex] = consblock;
            SCIPdebugPrintf("VARINDEX: %d (%d)\n", curvarindex, consblock);

            if ((blockrepresentative[oldblock] != -1) && (blockrepresentative[oldblock] > blockrepresentative[consblock]))
            {
               int oldrepr;
               oldrepr = blockrepresentative[oldblock];
               SCIPdebugMessage("\t\tBlock representative from block %d changed from %d to %d.\n", oldblock, blockrepresentative[oldblock], consblock);
               assert(consblock > 0);
               blockrepresentative[oldblock] = consblock;
               if ((oldrepr != consblock) && (oldrepr != oldblock))
               {
                  blockrepresentative[oldrepr] = consblock;
                  SCIPdebugMessage("\t\tBlock representative from block %d changed from %d to %d.\n", oldrepr, blockrepresentative[oldrepr], consblock);
               }
            }
         }
         else
         {
            SCIPdebugPrintf("will not be changed from %d to %d.\n", oldblock, consblock);
         }
      }

      SCIPfreeBufferArrayNull(scip, &curvars);
      assert(consblock >= 1 || consblock == -1);
      assert(consblock <= *nextblock);

      /* store the constraint block */
      if (consblock != -1)
      {
         SCIPdebugMessage("cons %s in block %d\n", SCIPconsGetName(cons), consblock);
         SCIP_CALL( SCIPhashmapInsert(constoblock, cons, (void*)(size_t)consblock) );
      }
      else
      {
         SCIPdebugMessage("ignoring %s\n", SCIPconsGetName(cons));
      }
   }

   return SCIP_OKAY;
}

/** */
static
SCIP_RETCODE fillConstoblock(
   SCIP_CONS**           conss,
   int                   nconss,
   SCIP_Bool*            consismaster,       /**< */
   int                   nblocks,            /**< */
   SCIP_HASHMAP*         constoblock,        /**< */
   SCIP_HASHMAP*         newconstoblock,     /**< */
   int*                  blockrepresentative /**< */
   )
{
   int i;

   /* convert temporary data to detectordata */
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      SCIP_CONS* cons = conss[i];

      if( GCGisConsGCGCons(cons) )
         continue;

      if( consismaster[i] )
      {
         SCIP_CALL( SCIPhashmapInsert(newconstoblock, cons, (void*) (size_t) (nblocks+1)) );
         continue;
      }

      if( !SCIPhashmapExists(constoblock, cons) )
         continue;

      consblock = (int) (size_t) SCIPhashmapGetImage(constoblock, cons); /*lint !e507*/
      assert(consblock > 0);
      consblock = blockrepresentative[consblock];
      assert(consblock <= nblocks);
      SCIP_CALL( SCIPhashmapInsert(newconstoblock, cons, (void*)(size_t)consblock) );
      SCIPdebugMessage("%d %s\n", consblock, SCIPconsGetName(cons));
   }
   return SCIP_OKAY;
}

/** creates a decomposition with provided constraints in the master
 * The function will put the remaining constraints in one or more pricing problems
 * depending on whether the subproblems decompose with no variables in common.
 */
SCIP_RETCODE DECcreateDecompFromMasterconss(
   SCIP*                 scip,                /**< SCIP data structure */
   DEC_DECOMP**          decomp,              /**< decomposition structure */
   SCIP_CONS**           masterconss,         /**< constraints to be put in the master */
   int                   nmasterconss         /**< number of constraints in the master */
   )
{
   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* newconstoblock;
   SCIP_CONS** conss;
   int nconss;
   SCIP_VAR** vars;
   int nvars;
   int nblocks;
   int* blockrepresentative;
   int nextblock = 1;
   SCIP_Bool* consismaster;
   int i;
   int* vartoblock;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert((masterconss == NULL) == ( nmasterconss == 0));
   assert(SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED);

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   assert( nmasterconss <= nconss );


   if( GCGisConsGCGCons(conss[nconss-1]) )
      --nconss;

   nblocks = nconss-nmasterconss+1;
   assert(nblocks > 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &blockrepresentative, nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consismaster, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &vartoblock, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&newconstoblock, SCIPblkmem(scip), nconss) );

   for( i = 0; i < nmasterconss; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(constoblock, masterconss[i], (void*)(size_t) (nblocks+1)) );
   }

   for( i = 0; i < nconss; ++i )
   {
      assert(!GCGisConsGCGCons(conss[i]));
      consismaster[i] = SCIPhashmapExists(constoblock, conss[i]);
   }

   for( i = 0; i < nvars; ++i )
   {
      vartoblock[i] = -1;
   }

   for( i = 0; i < nblocks; ++i )
   {
      blockrepresentative[i] = -1;
   }

   SCIP_CALL( assignConstraintsToRepresentatives(scip, conss, nconss, consismaster, constoblock, vartoblock, &nextblock, blockrepresentative) );

   /* postprocess blockrepresentatives */
   nblocks = processBlockRepresentatives(nextblock, blockrepresentative);

   /* convert temporary data to detectordata */
   SCIP_CALL( fillConstoblock(conss, nconss, consismaster, nblocks, constoblock, newconstoblock, blockrepresentative) );
   SCIP_CALL( DECdecompCreate(scip, decomp) );
   SCIP_CALL( DECfilloutDecdecompFromConstoblock(scip, *decomp, newconstoblock, nblocks, vars, nvars, conss, nconss, FALSE) );

   SCIPfreeMemoryArray(scip, &blockrepresentative);
   SCIPfreeMemoryArray(scip, &consismaster);
   SCIPfreeMemoryArray(scip, &vartoblock);
   SCIPhashmapFree(&constoblock);

   return SCIP_OKAY;
}

/** increase the corresponding count of the */
static
void incVarsData(
   SCIP_VAR*              var,                /**< variable to consider */
   int*                   nbinvars,           /**< pointer to array of size nproblems to store number of binary subproblem vars */
   int*                   nintvars,           /**< pointer to array of size nproblems to store number of integer subproblem vars */
   int*                   nimplvars,          /**< pointer to array of size nproblems to store number of implied subproblem vars */
   int*                   ncontvars,          /**< pointer to array of size nproblems to store number of continues subproblem vars */
   int                    nproblems,          /**< size of the arrays*/
   int                    i                   /**< index of the array to increase */
)
{
   assert(var != NULL);
   assert(i >= 0);
   assert(i < nproblems);

   if( nbinvars != NULL && (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarIsBinary(var)) )
   {
      ++(nbinvars[i]);
      assert(nbinvars[i] > 0);
   }
   if( nintvars != NULL && (SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER && !SCIPvarIsBinary(var)) )
   {
      ++(nintvars[i]);
      assert(nintvars[i] > 0);
   }
   if( nimplvars != NULL && (SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT) )
   {
      ++(nimplvars[i]);
      assert(nimplvars[i] > 0);
   }
   if( ncontvars != NULL && (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS) )
   {
      ++(ncontvars[i]);
      assert(ncontvars[i] > 0);
   }
}

/* score methods */

/** return the number of variables and binary, integer, implied integer, continuous variables of all subproblems */
void DECgetSubproblemVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition structure */
   int*                  nvars,              /**< pointer to array of size nproblems to store number of subproblem vars or NULL */
   int*                  nbinvars,           /**< pointer to array of size nproblems to store number of binary subproblem vars or NULL */
   int*                  nintvars,           /**< pointer to array of size nproblems to store number of integer subproblem vars or NULL */
   int*                  nimplvars,          /**< pointer to array of size nproblems to store number of implied subproblem vars or NULL */
   int*                  ncontvars,          /**< pointer to array of size nproblems to store number of continuous subproblem vars or NULL */
   int                   nproblems           /**< size of the arrays*/
)
{
   int i;
   int j;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(nproblems > 0);

   assert(DECdecompGetType(decomp) != DEC_DECTYPE_UNKNOWN);
   if( nvars != NULL )
      BMSclearMemoryArray(nvars, nproblems);
   if( nbinvars != NULL )
      BMSclearMemoryArray(nbinvars, nproblems);
   if( nintvars != NULL )
      BMSclearMemoryArray(nintvars, nproblems);
   if( nimplvars != NULL )
      BMSclearMemoryArray(nimplvars, nproblems);
   if( ncontvars != NULL )
      BMSclearMemoryArray(ncontvars, nproblems);

   for( i = 0; i < nproblems; ++i )
   {
      SCIP_VAR*** subscipvars;
      int* nsubscipvars;

      nsubscipvars = DECdecompGetNSubscipvars(decomp);
      subscipvars = DECdecompGetSubscipvars(decomp);
      if( nvars != NULL )
         nvars[i] = nsubscipvars[i];

      for( j = 0; j < nsubscipvars[i]; ++j )
      {
         incVarsData(subscipvars[i][j], nbinvars, nintvars, nimplvars, ncontvars, nproblems, i);
      }
   }
}


/** return the number of variables and binary, integer, implied integer, continuous variables of the master */
void DECgetLinkingVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition structure */
   int*                  nvars,              /**< pointer to store number of linking vars or NULL */
   int*                  nbinvars,           /**< pointer to store number of binary linking vars or NULL */
   int*                  nintvars,           /**< pointer to store number of integer linking vars or NULL */
   int*                  nimplvars,          /**< pointer to store number of implied linking vars or NULL */
   int*                  ncontvars           /**< pointer to store number of continuous linking vars or NULL */
)
{
   int i;
   SCIP_VAR** linkingvars;
   int nlinkingvars;

   assert(scip != NULL);
   assert(decomp != NULL);

   assert(DECdecompGetType(decomp) != DEC_DECTYPE_UNKNOWN);

   nlinkingvars = DECdecompGetNLinkingvars(decomp);
   linkingvars = DECdecompGetLinkingvars(decomp);

   if( nvars != NULL )
      *nvars = nlinkingvars;
   if( nbinvars != NULL )
      *nbinvars = 0;
   if( nintvars != NULL )
      *nintvars = 0;
   if( nimplvars != NULL )
      *nimplvars = 0;
   if( ncontvars != NULL )
      *ncontvars = 0;


   for( i = 0; i < nlinkingvars; ++i )
   {
      incVarsData(linkingvars[i], nbinvars, nintvars, nimplvars, ncontvars, 1, 0);
   }
}

/**
 * returns the number of nonzeros of each column of the constraint matrix both in the subproblem and in the master
 * @note For linking variables, the number of nonzeros in the subproblems corresponds to the number on nonzeros
 * in the border
 *
 * @note The arrays have to be allocated by the caller
 *
 * @pre This function assumes that constraints are partitioned in the decomp structure, no constraint is present in more than one block
 *
 */
SCIP_RETCODE DECgetDensityData(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition structure */
   SCIP_VAR**            vars,               /**< pointer to array store variables belonging to density */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< pointer to array to store constraints belonging to the density */
   int                   nconss,             /**< number of constraints */
   int*                  varsubproblemdensity, /**< pointer to array to store the nonzeros for the subproblems */
   int*                  varmasterdensity,   /**< pointer to array to store the nonzeros for the master */
   int*                  conssubproblemdensity, /**< pointer to array to store the nonzeros for the subproblems */
   int*                  consmasterdensity   /**< pointer to array to store the nonzeros for the master */
)
{
   int nlinkingconss;
   SCIP_HASHMAP* vartoblock;
   SCIP_CONS** curconss;

   int ncurvars;
   SCIP_VAR** curvars;
   SCIP_Bool success;

   int i;
   int j;
   int v;
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(nvars == SCIPgetNVars(scip));
   assert(conss != NULL);
   assert(nconss == SCIPgetNConss(scip));
   assert(varsubproblemdensity != NULL);
   assert(varmasterdensity != NULL);
   assert(conssubproblemdensity != NULL);
   assert(consmasterdensity != NULL);

   /* make sure the passed data is initialised to 0 */
   BMSclearMemoryArray(vars, nvars);
   BMSclearMemoryArray(conss, nconss);
   BMSclearMemoryArray(varsubproblemdensity, nvars);
   BMSclearMemoryArray(varmasterdensity, nvars);
   BMSclearMemoryArray(conssubproblemdensity, nconss);
   BMSclearMemoryArray(consmasterdensity, nconss);

   BMScopyMemoryArray(vars, SCIPgetVars(scip), nvars);

   vartoblock = DECdecompGetVartoblock(decomp);
   c = 0;
   for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
   {
      curconss = DECdecompGetSubscipconss(decomp)[i];
      assert(curconss != NULL);

      for( j = 0; j < DECdecompGetNSubscipconss(decomp)[i]; ++j )
      {
         assert(c < nconss); /* This assertion and the logic forbids constraints in more than one block */
         conss[c] = curconss[j];

         SCIP_CALL( SCIPgetConsNVars(scip, curconss[j], &ncurvars, &success) );
         assert(success);
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetConsVars(scip, curconss[j], curvars, ncurvars, &success) );
         assert(success);

         for( v = 0; v < ncurvars; ++v )
         {
            SCIP_VAR* var;
            int block;
            int probindex;
            var = curvars[v];
            var = SCIPvarGetProbvar(var);
            probindex = SCIPvarGetProbindex(var);
            assert(probindex >= 0);
            assert(probindex < nvars);
            varsubproblemdensity[probindex] += 1;
            assert(varsubproblemdensity[probindex] > 0);
            block = (int) (size_t) SCIPhashmapGetImage(vartoblock, var); /*lint !e507*/
            assert(block > 0);
            if( block <= DECdecompGetNBlocks(decomp) )
            {
               conssubproblemdensity[c] +=1;
            }
            else
            {
               consmasterdensity[c] += 1;
            }
         }

         SCIPfreeMemoryArray(scip, &curvars);
         c++;
      }
   }

   nlinkingconss = DECdecompGetNLinkingconss(decomp);
   curconss = DECdecompGetLinkingconss(decomp);

   for( j = 0; j < nlinkingconss; ++j )
   {
      assert(c < nconss); /* This assertion and the logic forbids constraints in more than one block */
      SCIP_CALL( SCIPgetConsNVars(scip, curconss[j], &ncurvars, &success) );
      assert(success);
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
      SCIP_CALL( SCIPgetConsVars(scip, curconss[j], curvars, ncurvars, &success) );
      assert(success);

      conss[c] = curconss[j];

      for( v = 0; v < ncurvars; ++v )
      {
         SCIP_VAR* var;
         int probindex;
         var = curvars[v];
         var = SCIPvarGetProbvar(var);
         probindex = SCIPvarGetProbindex(var);
         assert(probindex >= 0);
         assert(probindex < nvars);
         varmasterdensity[probindex] += 1;
         assert(varmasterdensity[probindex] > 0);
         SCIPdebugMessage("Var <%s> appears in cons <%s>, total count: %d\n", SCIPvarGetName(var), SCIPconsGetName(curconss[j]), varmasterdensity[probindex]);
      }

      consmasterdensity[c] = ncurvars;
      c++;

      SCIPfreeMemoryArray(scip, &curvars);
   }

   return SCIP_OKAY;
}

/** helper function to increase correct lock */
static
void increaseLock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lhs,                /**< left side of constraint */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             rhs,                /**< right side of constraint */
   int*                  downlock,           /**< pointer to store downlock */
   int*                  uplock              /**< pointer to store uplock */
   )
{
   assert(scip != NULL);
   assert(downlock != NULL);
   assert(uplock != NULL);

   if( !SCIPisInfinity(scip, -lhs) )
   {
      if( SCIPisPositive(scip, coef) )
         ++(*downlock);
      if( SCIPisNegative(scip, coef) )
         ++(*uplock);

   }
   if( !SCIPisInfinity(scip, rhs) )
   {
      if( SCIPisPositive(scip, coef) )
         ++(*uplock);
      if( SCIPisNegative(scip, coef) )
         ++(*downlock);
   }

}

/**
 *  calculates the number of up and down locks of variables for a given decomposition in both the original problem and the pricingproblems
 *
 *  @note All arrays need to be allocated by the caller
 *
 *  @warning This function needs a lot of memory (nvars*nblocks+1) array entries
 */
SCIP_RETCODE DECgetVarLockData(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition structure */
   SCIP_VAR**            vars,               /**< pointer to array store variables belonging to density */
   int                   nvars,              /**< number of variables */
   int                   nsubproblems,       /**< number of sub problems */
   int**                 subsciplocksdown,   /**< pointer to two dimensional array to store the down locks for the subproblems */
   int**                 subsciplocksup,     /**< pointer to two dimensional array to store the down locks for the subproblems */
   int*                  masterlocksdown,    /**< pointer to array to store the down locks for the master */
   int*                  masterlocksup       /**< pointer to array to store the down locks for the master */
   )
{
   int nlinkingconss;
   SCIP_CONS** curconss;
   SCIP_VAR** curvars;
   SCIP_Real* curvals;
   int ncurvars;
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_Bool success;

   int i;
   int j;
   int v;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(nvars == SCIPgetNVars(scip));
   assert(subsciplocksdown != NULL);
   assert(subsciplocksup != NULL);
   assert(masterlocksdown != NULL);
   assert(masterlocksup != NULL);

   /* make sure the passed data is initialised to 0 */
   BMSclearMemoryArray(vars, nvars);
   BMScopyMemoryArray(vars, SCIPgetVars(scip), nvars);
   BMSclearMemoryArray(masterlocksdown, nvars);
   BMSclearMemoryArray(masterlocksup, nvars);
   for( i = 0; i < nsubproblems; ++i )
   {
      BMSclearMemoryArray(subsciplocksdown[i], nvars);
      BMSclearMemoryArray(subsciplocksup[i], nvars);
   }

   for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
   {
      curconss = DECdecompGetSubscipconss(decomp)[i];
      assert(curconss != NULL);

      for( j = 0; j < DECdecompGetNSubscipconss(decomp)[i]; ++j )
      {

         SCIP_CALL( SCIPgetConsNVars(scip, curconss[j], &ncurvars, &success) );
         assert(success);
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvals, ncurvars) );

         SCIP_CALL( SCIPgetValsXXX(scip, curconss[j], curvals, ncurvars) );
         SCIP_CALL( SCIPgetConsVars(scip, curconss[j], curvars, ncurvars, &success) );
         assert(success);

         rhs = SCIPgetRhsXXX(scip, curconss[j]);
         lhs = SCIPgetLhsXXX(scip, curconss[j]);

         for( v = 0; v < ncurvars; ++v )
         {
            SCIP_VAR* var;
            int probindex;
            var = curvars[v];
            var = SCIPvarGetProbvar(var);
            probindex = SCIPvarGetProbindex(var);
            assert(probindex >= 0);
            assert(probindex < nvars);

            assert((size_t)SCIPhashmapGetImage(DECdecompGetVartoblock(decomp), var) > 0);
            increaseLock(scip, lhs, curvals[v], rhs, &(subsciplocksdown[i][probindex]), &(subsciplocksup[i][probindex]));
         }

         SCIPfreeMemoryArray(scip, &curvars);
         SCIPfreeMemoryArray(scip, &curvals);
      }
   }

   nlinkingconss = DECdecompGetNLinkingvars(decomp);
   curconss = DECdecompGetLinkingconss(decomp);
   for( j = 0; j < nlinkingconss; ++j )
   {
      SCIP_CALL( SCIPgetConsNVars(scip, curconss[j], &ncurvars, &success) );
      assert(success);
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvals, ncurvars) );

      SCIP_CALL( SCIPgetValsXXX(scip, curconss[j], curvals, ncurvars) );
      SCIP_CALL( SCIPgetConsVars(scip, curconss[j], curvars, ncurvars, &success) );
      assert(success);

      rhs = SCIPgetRhsXXX(scip, curconss[j]);
      lhs = SCIPgetLhsXXX(scip, curconss[j]);

      for( v = 0; v < ncurvars; ++v )
      {
         SCIP_VAR* var;
         int probindex;
         var = curvars[v];
         var = SCIPvarGetProbvar(var);
         probindex = SCIPvarGetProbindex(var);
         assert(probindex >= 0);
         assert(probindex < nvars);
         increaseLock(scip, lhs, curvals[v], rhs, &(masterlocksdown[probindex]), &(masterlocksup[probindex]));
      }

      SCIPfreeMemoryArray(scip, &curvars);
      SCIPfreeMemoryArray(scip, &curvals);
   }

   return SCIP_OKAY;
}


/** computes the score of the given decomposition based on the border, the average density score and the ratio of
 * linking variables
 */
SCIP_RETCODE DECevaluateDecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition data structure */
   DEC_SCORES*           score               /**< returns the score of the decomposition */
   )
{
   int matrixarea;
   int borderarea;
   int nvars;
   int nconss;
   int i;
   int j;
   int k;
   /*   int blockarea; */
   SCIP_Real varratio;
   int* nzblocks;
   int nblocks;
   int* nlinkvarsblocks;
   int* nvarsblocks;
   SCIP_Real* blockdensities;
   int* blocksizes;
   SCIP_Real density;

   assert(scip != NULL);
   assert(score != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   nblocks = DECdecompGetNBlocks(decdecomp);

   SCIP_CALL( SCIPallocBufferArray(scip, &nzblocks, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlinkvarsblocks, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockdensities, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocksizes, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvarsblocks, nblocks) );
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nvars*nconss;

   /* calculate slave sizes, nonzeros and linkingvars */
   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CONS** curconss;
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled;

      SCIP_CALL( SCIPallocBufferArray(scip, &ishandled, nvars) );
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;
      for( j = 0; j < nvars; ++j )
      {
         ishandled[j] = FALSE;
      }
      curconss = DECdecompGetSubscipconss(decdecomp)[i];
      ncurconss = DECdecompGetNSubscipconss(decdecomp)[i];

      for( j = 0; j < ncurconss; ++j )
      {
         SCIP_VAR** curvars;
         SCIP_VAR* var;
         int ncurvars;
         ncurvars = SCIPgetNVarsXXX(scip, curconss[j]);
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, curconss[j], curvars, ncurvars) );

         for( k = 0; k < ncurvars; ++k )
         {
            int block;
            if( !SCIPisVarRelevant(curvars[k]) )
               continue;

            var = SCIPvarGetProbvar(curvars[k]);
            assert(var != NULL);
            if( !SCIPisVarRelevant(var) )
               continue;

            assert(SCIPvarIsActive(var));
            assert(!SCIPvarIsDeleted(var));
            ++(nzblocks[i]);
            if( !SCIPhashmapExists(DECdecompGetVartoblock(decdecomp), var) )
            {
               block = (int)(size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), curvars[k]); /*lint !e507*/
            }
            else
            {
               assert(SCIPhashmapExists(DECdecompGetVartoblock(decdecomp), var));
               block = (int)(size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), var); /*lint !e507*/
            }

            if( block == nblocks+1 && ishandled[SCIPvarGetProbindex(var)] == FALSE )
            {
               ++(nlinkvarsblocks[i]);
            }
            ishandled[SCIPvarGetProbindex(var)] = TRUE;
         }

         SCIPfreeBufferArray(scip, &curvars);
      }

      for( j = 0; j < nvars; ++j )
      {
         if( ishandled[j] )
         {
            ++nvarsblock;
         }
      }

      blocksizes[i] = nvarsblock*ncurconss;
      nvarsblocks[i] = nvarsblock;
      if( blocksizes[i] > 0 )
      {
         blockdensities[i] = 1.0*nzblocks[i]/blocksizes[i];
      }
      else
      {
         blockdensities[i] = 0.0;
      }

      assert(blockdensities[i] >= 0 && blockdensities[i] <= 1.0);
      SCIPfreeBufferArray(scip, &ishandled);
   }

   borderarea = DECdecompGetNLinkingconss(decdecomp)*nvars+DECdecompGetNLinkingvars(decdecomp)*(nconss-DECdecompGetNLinkingconss(decdecomp));

   density = 1E20;
   varratio = 1.0;
   for( i = 0; i < nblocks; ++i )
   {
      density = MIN(density, blockdensities[i]);

      if( DECdecompGetNLinkingvars(decdecomp) > 0 )
      {
         varratio *= 1.0*nlinkvarsblocks[i]/DECdecompGetNLinkingvars(decdecomp);
      }
      else
      {
         varratio = 0;
      }
   }

   score->linkingscore = (0.5+0.5*varratio);
   score->borderscore = (1.0*(borderarea)/matrixarea);
   score->densityscore = (1-density);

   switch( DECdecompGetType(decdecomp) )
   {
   case DEC_DECTYPE_ARROWHEAD:
      score->totalscore = score->borderscore*score->linkingscore*score->densityscore;
      break;
   case DEC_DECTYPE_BORDERED:
      score->totalscore = score->borderscore*score->linkingscore*score->densityscore;
      break;
   case DEC_DECTYPE_DIAGONAL:
      score->totalscore = 0.0;
      break;
   case DEC_DECTYPE_STAIRCASE:
      SCIPwarningMessage(scip, "Decomposition type is %s, cannot compute score\n", DECgetStrType(DECdecompGetType(decdecomp)));
      score->totalscore = 0.1;
      break;
   case DEC_DECTYPE_UNKNOWN:
      SCIPerrorMessage("Decomposition type is %s, cannot compute score\n", DECgetStrType(DECdecompGetType(decdecomp)));
      assert(FALSE);
      break;
   default:
      SCIPerrorMessage("No rule for this decomposition type, cannot compute score\n");
      assert(FALSE);
      break;
   }

   SCIPfreeBufferArray(scip, &nzblocks);
   SCIPfreeBufferArray(scip, &nlinkvarsblocks);
   SCIPfreeBufferArray(scip, &blockdensities);
   SCIPfreeBufferArray(scip, &blocksizes);
   SCIPfreeBufferArray(scip, &nvarsblocks);

   return SCIP_OKAY;
}

static
SCIP_RETCODE computeVarDensities(
      SCIP*              scip,               /**< SCIP data structure */
      DEC_DECOMP*        decomp,             /**< decomposition data structure */
      int*               varprobdensity,     /**< density information */
      int*               varmasterdensity,   /**< density information */
      SCIP_VAR**         vars,               /**< */
      int                nvars,              /**< */
      DEC_STATISTIC*     blockvardensities,  /**< */
      DEC_STATISTIC*     mastervardensity,   /**< */
      int                nblocks             /**< */
   )
{
   int v;
   int b;
   SCIP_Real** vardistribution;
   int* nvardistribution;
   SCIP_Real* mastervardistribution;

   SCIP_Real max = 0;
   SCIP_Real min = 1.0;
   SCIP_Real median = 0;
   SCIP_Real mean = 0;

   assert(scip != NULL);
   assert(decomp != NULL);

   assert(vars != NULL);
   assert(nvars > 0);
   assert(blockvardensities != NULL);
   assert(mastervardensity != NULL);
   assert(nblocks >= 0);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vardistribution, nblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nvardistribution, nblocks) );

   BMSclearMemoryArray(vardistribution, nblocks);
   BMSclearMemoryArray(nvardistribution, nblocks);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &mastervardistribution, nvars) );
   BMSclearMemoryArray(mastervardistribution, nvars);

   for( b = 0; b < nblocks; ++b )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vardistribution[b], DECdecompGetNSubscipvars(decomp)[b]) );
      BMSclearMemoryArray(vardistribution[b], DECdecompGetNSubscipvars(decomp)[b]);
   }

   for( v = 0; v < nvars; ++v )
   {
      int block = GCGvarGetBlock(vars[v]);
      SCIPdebugMessage("Var <%s>:", SCIPvarGetName(vars[v]));


      mastervardistribution[v] = 1.0*varmasterdensity[v]/DECdecompGetNLinkingconss(decomp);
      SCIPdebugPrintf("master %d ", varmasterdensity[v]);

      if( block >= 0 )
      {
         vardistribution[block][nvardistribution[block]] = 1.0*varprobdensity[v]/DECdecompGetNSubscipconss(decomp)[block];
         SCIPdebugPrintf("block %d %.3f\n", block, vardistribution[block][nvardistribution[block]]);
         ++(nvardistribution[block]);
      }
      else
      {
         SCIPdebugPrintf("\n");
      }
   }

   for( b = 0; b < nblocks; ++b )
   {
      int ncurvars = DECdecompGetNSubscipvars(decomp)[b];

      max = 0;
      min = 1.0;
      median = 0;
      mean = 0;



      SCIPdebugMessage("block %d:", b);
      for( v = 0; v < ncurvars; ++v )
      {

         SCIPdebugPrintf(" <%s> %.3f", SCIPvarGetName(DECdecompGetSubscipvars(decomp)[b][v]), vardistribution[b][v]);
         max = MAX(max, vardistribution[b][v]);
         min = MIN(min, vardistribution[b][v]);
         mean += 1.0*vardistribution[b][v]/ncurvars;

      }
      if( ncurvars > 0 )
         median = quick_select_median(vardistribution[b], ncurvars);

      SCIPdebugPrintf("\nmin: %.3f, max: %.3f, median: %.3f, mean: %.3f\n", min, max, median, mean);

      blockvardensities[b].max = max;
      blockvardensities[b].min = min;
      blockvardensities[b].median = median;
      blockvardensities[b].mean = mean;
   }
   max = 0;
   min = 1.0;
   mean = 0;

   SCIPdebugMessage("master:");

   for( v = 0; v < nvars; ++v )
   {

      SCIPdebugPrintf(" <%s> %.3f", SCIPvarGetName(vars[v]), mastervardistribution[v]);
      max = MAX(max, mastervardistribution[v]);
      min = MIN(min, mastervardistribution[v]);
      mean += 1.0*mastervardistribution[v]/nvars;

   }
   median = quick_select_median(mastervardistribution, nvars);
   SCIPdebugPrintf("\nmin: %.3f, max: %.3f, median: %.3f, mean: %.3f\n", min, max, median, mean);


   mastervardensity->max = max;
   mastervardensity->min = min;
   mastervardensity->median = median;
   mastervardensity->mean = mean;

   for( b = 0; b < nblocks; ++b )
   {
      SCIPfreeBlockMemoryArray(scip, &vardistribution[b], DECdecompGetNSubscipvars(decomp)[b]);
   }

   SCIPfreeBlockMemoryArray(scip, &mastervardistribution, nvars);

   SCIPfreeBlockMemoryArray(scip, &nvardistribution, nblocks);
   SCIPfreeBlockMemoryArray(scip, &vardistribution, nblocks);

   return SCIP_OKAY;
}


/** display statistics about the decomposition */
SCIP_RETCODE GCGprintDecompStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   )
{
   DEC_DECOMP* decomp;
   DEC_SCORES scores;
   SCIP_VAR** vars;
   SCIP_CONS** conss;

   int nvars;
   int nconss;

   int* nallvars;
   int* nbinvars;
   int* nintvars;
   int* nimplvars;
   int* ncontvars;

   int nblocks;
   int nblocksrelevant;
   int nlinkvars;
   int nlinkbinvar;
   int nlinkintvars;
   int nlinkimplvars;
   int nlinkcontvars;
   int b;

   int* varprobdensity;
   int* varmasterdensity;
   int* consprobsensity;
   int* consmasterdensity;

   DEC_STATISTIC* blockvardensities;
   DEC_STATISTIC* blockconsdensities;
   DEC_STATISTIC mastervardensity;

   assert(scip != NULL);

   decomp = DECgetBestDecomp(scip);
   assert(decomp != NULL);
   nblocks = DECdecompGetNBlocks(decomp);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) < SCIP_STAGE_PRESOLVED )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "No Dantzig-Wolfe reformulation applied. The problem was most likely already solved by the LP in the original problem.\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nallvars, nblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nbinvars, nblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nintvars, nblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nimplvars, nblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ncontvars, nblocks) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blockvardensities, nblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blockconsdensities, nblocks) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varprobdensity, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varmasterdensity, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conss, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consprobsensity, nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consmasterdensity, nconss) );

   SCIP_CALL( DECevaluateDecomposition(scip, decomp, &scores) );

   DECgetSubproblemVarsData(scip, decomp, nallvars, nbinvars, nintvars, nimplvars, ncontvars, nblocks);
   DECgetLinkingVarsData(scip, decomp, &nlinkvars, &nlinkbinvar, &nlinkintvars, &nlinkimplvars, &nlinkcontvars);
   SCIP_CALL( DECgetDensityData(scip, decomp, vars, nvars, conss, nconss, varprobdensity, varmasterdensity, consprobsensity, consmasterdensity) );

   SCIP_CALL( computeVarDensities(scip, decomp, varprobdensity, varmasterdensity, vars, nvars, blockvardensities, &mastervardensity, nblocks) );

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Decomp statistics  :\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  type             : %10s\n", DECgetStrType(DECdecompGetType(decomp)));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  detector         : %10s\n", decomp->detector == NULL? "provided": DECdetectorGetName(decomp->detector));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  blocks           : %10d\n", DECdecompGetNBlocks(decomp));

   nblocksrelevant = nblocks;
   for( b = 0; b < nblocks; ++b )
   {
      if( GCGrelaxGetNIdenticalBlocks(scip, b) == 0 )
         nblocksrelevant -= 1;
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  aggr. blocks     : %10d\n", nblocksrelevant);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Master statistics  :      nvars   nbinvars   nintvars  nimplvars  ncontvars     nconss  min(dens)  max(dens) medi(dens) mean(dens)\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  master           : %10d %10d %10d %10d %10d %10d %10.3f %10.3f %10.3f %10.3f\n", nlinkvars,
         nlinkbinvar, nlinkintvars, nlinkimplvars, nlinkcontvars, DECdecompGetNLinkingconss(decomp),
         mastervardensity.min, mastervardensity.max, mastervardensity.median, mastervardensity.mean);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Pricing statistics :      nvars   nbinvars   nintvars  nimplvars  ncontvars     nconss  min(dens)  max(dens) medi(dens) mean(dens)  identical\n");
   for( b = 0; b < nblocks; ++b )
   {
      if( GCGrelaxIsPricingprobRelevant(scip, b) )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " %10lld        : %10d %10d %10d %10d %10d %10d %10.3f %10.3f %10.3f %10.3f %10d\n", b+1, nallvars[b], nbinvars[b], nintvars[b], nimplvars[b], ncontvars[b],
               DECdecompGetNSubscipconss(decomp)[b], blockvardensities[b].min, blockvardensities[b].max, blockvardensities[b].median, blockvardensities[b].mean, GCGrelaxGetNIdenticalBlocks(scip, b));
      }
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Decomp Scores      :\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  border area      : %10.3f\n", scores.borderscore);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  avg. density     : %10.3f\n", scores.densityscore);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  linking score    : %10.3f\n", scores.linkingscore);

   SCIPfreeBlockMemoryArray(scip, &vars, nvars);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);

   SCIPfreeBlockMemoryArray(scip, &varprobdensity, SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &varmasterdensity, SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &consprobsensity, SCIPgetNConss(scip));
   SCIPfreeBlockMemoryArray(scip, &consmasterdensity, SCIPgetNConss(scip));

   SCIPfreeBlockMemoryArray(scip, &blockvardensities, nblocks);
   SCIPfreeBlockMemoryArray(scip, &blockconsdensities, nblocks);

   SCIPfreeBlockMemoryArray(scip, &nallvars, nblocks);
   SCIPfreeBlockMemoryArray(scip, &nbinvars, nblocks);
   SCIPfreeBlockMemoryArray(scip, &nintvars, nblocks);
   SCIPfreeBlockMemoryArray(scip, &nimplvars, nblocks);
   SCIPfreeBlockMemoryArray(scip, &ncontvars, nblocks);

   return SCIP_OKAY;
}

/** returns whether both structures lead to the same decomposition */
SCIP_Bool DECdecompositionsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp1,            /**< first decomp data structure */
   DEC_DECOMP*           decomp2             /**< second decomp data structure */
)
{
   SCIP_HASHMAP* constoblock1;
   SCIP_HASHMAP* constoblock2;

   SCIP_CONS** conss;
   int nconss;

   SCIP_VAR** vars;
   int i;

   assert(scip != NULL);
   assert(decomp1 != NULL);
   assert(decomp2 != NULL);

   if( DECdecompGetNBlocks(decomp1) != DECdecompGetNBlocks(decomp2) )
   {
      return FALSE;
   }

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   vars = SCIPgetVars(scip);

   constoblock1 = DECdecompGetConstoblock(decomp1);
   constoblock2 = DECdecompGetConstoblock(decomp2);

   for( i = 0; i < nconss; ++i )
   {
      if( SCIPhashmapGetImage(constoblock1, conss[i]) != SCIPhashmapGetImage(constoblock2, conss[i]) )
         return FALSE;
   }

   for( i = 0; i < nconss; ++i )
   {
      if( SCIPhashmapGetImage(constoblock1, vars[i]) != SCIPhashmapGetImage(constoblock2, vars[2]) )
         return FALSE;
   }

   return TRUE;
}

/** filters similar decompositions from a given list and moves them to the end
 * @return the number of unique decompositions
 */
int DECfilterSimilarDecompositions(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP**          decs,               /**< array of decompositions */
   int                   ndecs               /**< number of decompositions */
)
{
   int i;
   int j;
   int nunique;
   assert(scip != NULL);
   assert(decs != NULL);
   assert(ndecs > 0);

   nunique = ndecs;
   for( i = 0; i < nunique; ++i )
   {
      for( j = i+1; j < nunique; ++j )
      {
         DEC_DECOMP* tmp;
         if( DECdecompositionsAreEqual(scip, decs[i], decs[j]) )
         {
            tmp = decs[nunique-1];
            decs[nunique-1] = decs[j];
            decs[j] = tmp;
            --nunique;
            --j;
         }
      }
   }
   return nunique;
}

/** returns the number of the block that the constraint is with respect to the decomposition */
SCIP_RETCODE DECdetermineConsBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition */
   SCIP_CONS*            cons,               /**< constraint to check */
   int                   *block              /**< block of the constraint (or nblocks for master) */
)
{
   SCIP_VAR** curvars = NULL;
   int ncurvars = 0;
   SCIP_Bool success = FALSE;
   int i;
   int nblocks ;

   int nmastervars = 0;
   int npricingvars = 0;

   SCIP_HASHMAP* vartoblock;
   assert(scip != NULL);
   assert(decomp != NULL);
   assert(cons != NULL);
   assert(block != NULL);

   *block = -2;

   SCIP_CALL( SCIPgetConsNVars(scip, cons, &ncurvars, &success) );
   assert(success);

   if( ncurvars == 0 )
      return SCIP_OKAY;

   vartoblock= DECdecompGetVartoblock(decomp);
   assert(vartoblock != NULL);

   nblocks = DECdecompGetNBlocks(decomp);

   SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
   SCIP_CALL( SCIPgetConsVars(scip, cons, curvars, ncurvars, &success) );
   assert(success);

   for( i = 0; i < ncurvars && *block != nblocks; ++i )
   {
      int varblock;

      assert(SCIPhashmapExists(vartoblock, SCIPvarGetProbvar(curvars[i])));
      varblock = ((int) (size_t) SCIPhashmapGetImage(vartoblock, SCIPvarGetProbvar(curvars[i])))-1;

      /* if variable is linking skip*/
      if( varblock == nblocks+1 )
      {
         continue;
      }
      else if( varblock == nblocks )
      {
         ++nmastervars;
         continue;
      }
      else if( *block != varblock )
      {
         ++npricingvars;
         if( *block < 0 )
            *block = varblock;
         else
         {
            assert(*block != nblocks);
            *block = nblocks;
            break;
         }
      }
   }

   SCIPfreeBufferArrayNull(scip, &curvars);

   if( ncurvars > 0 && *block == -2 )
      *block = nblocks;

   if( npricingvars == 0 && nmastervars > 0)
      *block = -1;

   return SCIP_OKAY;
}

/** move a master constraint to pricing problem */
SCIP_RETCODE DECdecompMoveLinkingConsToPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition */
   int                   consindex,          /**< index of constraint to move */
   int                   block               /**< block of the pricing problem where to move */
   )
{
   SCIP_CONS* linkcons;
   SCIP_VAR** curvars = NULL;
   int ncurvars = 0;
   SCIP_Bool success = FALSE;
   int v;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(consindex >= 0 && consindex < decomp->nlinkingconss);
   assert(block >= 0 && block < decomp->nblocks);

   linkcons = decomp->linkingconss[consindex];

   SCIP_CALL( SCIPgetConsNVars(scip, linkcons, &ncurvars, &success) );
   assert(success);
   SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
   SCIP_CALL( SCIPgetConsVars(scip, linkcons, curvars, ncurvars, &success) );
   assert(success);

   decomp->linkingconss[consindex] =  decomp->linkingconss[decomp->nlinkingconss-1];
   decomp->nlinkingconss -= 1;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &decomp->subscipconss[block], decomp->nsubscipconss[block]+1) );
   decomp->subscipconss[block][decomp->nsubscipconss[block]] = linkcons;
   decomp->nsubscipconss[block] += 1;
   SCIP_CALL( SCIPhashmapSetImage(decomp->constoblock, linkcons, (void*) (size_t) (block+1)) );

   for( v = 0; v < ncurvars; ++v )
   {
      SCIP_VAR* probvar = SCIPvarGetProbvar(curvars[v]);
      assert(SCIPhashmapExists(decomp->vartoblock, probvar));
      /* if variable is in master only, move to subproblem */
      if( (int) (size_t) SCIPhashmapGetImage(decomp->vartoblock, probvar) == decomp->nblocks+1 )
      {
         SCIP_CALL( SCIPhashmapSetImage(decomp->vartoblock, probvar, (void*) (size_t) (block+1)) );
         SCIP_CALL( SCIPreallocMemoryArray(scip, &decomp->subscipvars[block], decomp->nsubscipvars[block] + 1) );
         decomp->subscipvars[block][decomp->nsubscipvars[block]] = probvar;
         decomp->nsubscipvars[block] += 1;
         SCIP_CALL( DECdecompRemoveLinkingVar(scip, decomp, probvar, &success) );
         assert(success);
      }
   }

   SCIPfreeBufferArrayNull(scip, &curvars);
   return SCIP_OKAY;
}


/** tries to assign masterconss to pricing problem */
SCIP_RETCODE DECtryAssignMasterconssToExistingPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition */
   int*                  transferred         /**< number of master constraints reassigned */
   )
{
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(transferred != NULL);

   *transferred = 0;

   for( c = 0; c < decomp->nlinkingconss; ++c )
   {
      int block;
      SCIP_CALL( DECdetermineConsBlock(scip, decomp, decomp->linkingconss[c], &block) );

      if( block == DECdecompGetNBlocks(decomp) || block < 0)
      {
         continue;
      }

      SCIP_CALL( DECdecompMoveLinkingConsToPricing(scip, decomp, c, block) );
      --c;
      *transferred += 1;
   }

   if( *transferred > 0 )
   {
      if( decomp->nlinkingconss > 0 )
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &decomp->linkingconss, decomp->nlinkingconss) );
      }
      else
      {
         SCIPfreeMemoryArrayNull(scip, &decomp->linkingconss);
      }
   }

   return SCIP_OKAY;
}

/** Removes a variable from the linking variable array */
SCIP_RETCODE DECdecompRemoveLinkingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition */
   SCIP_VAR*             var,                /**< variable to remove */
   SCIP_Bool*            success             /**< indicates whether the variable was successfully removed */
   )
{
   int v;
   assert(scip != NULL);
   assert(decomp != NULL);
   assert(var != NULL);
   assert(success != NULL);

   *success = FALSE;

   for( v = 0; v < decomp->nlinkingvars; ++v )
   {
      if( decomp->linkingvars[v] == var )
      {
         decomp->linkingvars[v] = decomp->linkingvars[decomp->nlinkingvars-1];
         decomp->nlinkingvars -= 1;
         *success = TRUE;
      }
   }

   if( *success )
   {
      if( decomp->nlinkingvars == 0)
      {
         SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);
         if( DECdecompGetNLinkingconss(decomp) == 0)
         {
            SCIP_CALL( DECdecompSetType(decomp, DEC_DECTYPE_DIAGONAL) );
         }
         else
         {
            SCIP_CALL( DECdecompSetType(decomp, DEC_DECTYPE_BORDERED) );
         }
      }
      else
      {
         SCIPreallocMemoryArray(scip, &decomp->linkingvars, decomp->nlinkingvars);
      }
   }
   return SCIP_OKAY;
}

/** tries to assign masterconss to new and existing pricing problems */
SCIP_RETCODE DECtryAssignMasterconssToNewPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition */
   DEC_DECOMP**          newdecomp,          /**< new decomposition, if successful */
   int*                  transferred         /**< number of master constraints reassigned */
   )
{
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(newdecomp != NULL);
   assert(transferred != NULL);

   *newdecomp = NULL;
   *transferred = 0;

   for( c = 0; c < decomp->nlinkingconss; ++c )
   {
      int block;
      int i;
      int nconss;
      SCIP_HASHMAP* constoblock;
      SCIP_CALL( DECdetermineConsBlock(scip, decomp, decomp->linkingconss[c], &block) );

      if( block >= 0)
      {
         continue;
      }
      SCIPdebugMessage("Cons <%s> in new pricing problem\n", SCIPconsGetName(decomp->linkingconss[c]));
      nconss = SCIPgetNConss(scip);
      SCIP_CALL( DECdecompCreate(scip, newdecomp) );
      SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), SCIPgetNConss(scip)) );

      for( i = 0; i < nconss; ++i )
      {
         int consblock;
         SCIP_CONS* cons = SCIPgetConss(scip)[i];
         assert(SCIPhashmapExists(decomp->constoblock, cons));
         consblock = (int) (size_t) SCIPhashmapGetImage(decomp->constoblock, cons);
         SCIPdebugMessage("Cons <%s> %d -> %d\n", SCIPconsGetName(cons), consblock, consblock+1);

         SCIP_CALL( SCIPhashmapSetImage(constoblock, cons, (void*) (size_t) (consblock+1)) );
      }
      SCIP_CALL( SCIPhashmapSetImage(constoblock, decomp->linkingconss[c], (void*) (size_t) (1)) );
      SCIPdebugMessage("Cons <%s>    -> %d\n", SCIPconsGetName(decomp->linkingconss[c]), 1);

      SCIP_CALL( DECfilloutDecdecompFromConstoblock(scip, *newdecomp, constoblock, decomp->nblocks+1, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );
      *transferred += 1;
      break;
   }

   return SCIP_OKAY;
}

/** polish the decomposition and try to greedily assign master constraints to pricing problem where usefule */
SCIP_RETCODE DECcreatePolishedDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decomp,             /**< decomposition */
   DEC_DECOMP**          newdecomp           /**< new decomposition, if successful */
   )
{
   int transferred = 0;
   DEC_DECOMP* origdecomp = decomp;
   DEC_DECOMP* tempdecomp = NULL;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(newdecomp != NULL);
   *newdecomp = decomp;

   do
   {
      SCIP_CALL( DECtryAssignMasterconssToExistingPricing(scip, *newdecomp, &transferred) );
      SCIPdebugMessage("%d conss transferred to existing pricing\n", transferred);
      SCIP_CALL( DECtryAssignMasterconssToNewPricing(scip, *newdecomp, &tempdecomp, &transferred) );
      SCIPdebugMessage("%d conss transferred to new pricing\n", transferred);
      if( transferred > 0 )
      {
         if( *newdecomp != origdecomp )
         {
            SCIP_CALL( DECdecompFree(scip, newdecomp) );
         }
         *newdecomp = tempdecomp;
      }
   } while( transferred > 0 );

   if( *newdecomp == origdecomp )
   {
      *newdecomp = NULL;
   }

   return SCIP_OKAY;
}
