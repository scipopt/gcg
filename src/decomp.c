/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2012 Operations Research, RWTH Aachen University       */
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
#include "pub_decomp.h"
#include "scip/scip.h"
#include "struct_decomp.h"

#include <assert.h>

/** converts the DEC_DECTYPE enum to a string */
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
   (*decomp)->nblocks = 0;
   (*decomp)->consindex = NULL;
   (*decomp)->varindex = NULL;

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
      for( j = 0; j < decomp->nsubscipvars[i]; ++j )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(decomp->subscipvars[i][j])) );
      }
      SCIPfreeMemoryArrayNull(scip, &(decomp->subscipvars[i]));

      for( j = 0; j < decomp->nsubscipconss[i]; ++j )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(decomp->subscipconss[i][j])) );
      }
      SCIPfreeMemoryArray(scip, &decomp->subscipconss[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipconss);

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
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingconss);

   for( i = 0; i < decomp->nlinkingvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(decomp->linkingvars[i])) );
   }
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);

   SCIPfreeMemory(scip, decdecomp);

   return SCIP_OKAY;
}

/** sets the type of the decomposition */
void DECdecompSetType(
   DEC_DECOMP*           decdecomp,          /**< pointer to the decdecomp instance */
   DEC_DECTYPE           type,               /**< type of the decomposition */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   assert(decdecomp != NULL);
   switch( type )
   {
   case DEC_DECTYPE_DIAGONAL:
      *valid = decdecomp->nlinkingconss == 0 && decdecomp->linkingconss == NULL;
      *valid = *valid && decdecomp->nlinkingvars == 0 && decdecomp->linkingvars == NULL;
      break;
   case DEC_DECTYPE_ARROWHEAD:
      *valid = TRUE;
      break;
   case DEC_DECTYPE_UNKNOWN:
      *valid = FALSE;
      break;
   case DEC_DECTYPE_BORDERED:
      *valid = decdecomp->nlinkingvars == 0 && decdecomp->linkingvars == NULL;
      break;
   case DEC_DECTYPE_STAIRCASE:
      *valid = decdecomp->nlinkingconss == 0 && decdecomp->linkingconss == NULL;
      break;
   default:
      *valid = FALSE;
      break;
   }

   decdecomp->type = type;
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
   int*                  nsubscipvars,       /**< number of subscipvars per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipvars != NULL);
   assert(nsubscipvars != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->subscipvars == NULL);
   assert(decdecomp->nsubscipvars == NULL);

   *valid = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipvars, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipvars, decdecomp->nblocks) );

   assert(decdecomp->subscipvars != NULL);
   assert(decdecomp->nsubscipvars != NULL);

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      assert((subscipvars[b] == NULL) == (nsubscipvars[b] == 0));
      decdecomp->nsubscipvars[b] = nsubscipvars[b];

      if( nsubscipvars[b] < 0 )
         *valid = FALSE;
      else if( nsubscipvars[b] > 0 )
      {
         assert(subscipvars[b] != NULL);
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &(decdecomp->subscipvars[b]), subscipvars[b], nsubscipvars[b]) ); /*lint !e866*/

         for( i = 0; i < nsubscipvars[b]; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, decdecomp->subscipvars[b][i]) );
         }
      }
      else if( nsubscipvars[b] == 0 )
      {
         decdecomp->subscipvars[b] = NULL;
      }
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
   int*                  nsubscipconss,      /**< number of subscipconss per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipconss != NULL);
   assert(nsubscipconss != NULL);
   assert(valid != NULL);

   assert(decdecomp->nblocks > 0);
   assert(decdecomp->subscipconss == NULL);
   assert(decdecomp->nsubscipconss == NULL);

   *valid = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipconss, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipconss, decdecomp->nblocks) );

   assert(decdecomp->subscipconss != NULL);
   assert(decdecomp->nsubscipconss != NULL);

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      if( nsubscipconss[b] <= 0 || subscipconss[b] == NULL )
         *valid = FALSE;

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
   int                   nlinkingconss,      /**< number of linkingconss per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingconss != NULL);
   assert(nlinkingconss >= 0);

   assert(decdecomp->linkingconss == NULL);
   assert(decdecomp->nlinkingconss == 0);

   decdecomp->nlinkingconss = nlinkingconss;

   if( nlinkingconss > 0 )
   {
      assert(linkingconss != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingconss, linkingconss, nlinkingconss) );

      for( i = 0; i < nlinkingconss; ++i )
      {
         SCIP_CALL( SCIPcaptureCons(scip, decdecomp->linkingconss[i]) );
      }
   }

   *valid = linkingconss != NULL || nlinkingconss == 0;

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
   int                   nlinkingvars,       /**< number of linkingvars per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingvars != NULL || nlinkingvars == 0);

   assert(decdecomp->linkingvars == NULL);
   assert(decdecomp->nlinkingvars == 0);

   decdecomp->nlinkingvars = nlinkingvars;

   if( nlinkingvars > 0 )
   {
      assert(linkingvars != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingvars, linkingvars, nlinkingvars) );

      for( i = 0; i < nlinkingvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, decdecomp->linkingvars[i]) );
      }
   }

   *valid = linkingvars != NULL || nlinkingvars == 0;

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

/** sets the vartoblock hashmap of the given decdecomp structure */
void  DECdecompSetVartoblock(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         vartoblock,         /**< Vartoblock hashmap */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);

   *valid = TRUE;

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
   SCIP_HASHMAP*         constoblock,        /**< Constoblock hashmap */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   assert(decdecomp != NULL);
   assert(constoblock != NULL);

   *valid = TRUE;

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
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;
   int i;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);
   assert(constoblock != NULL);
   assert(nblocks > 0);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(valid != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   nlinkingconss = 0;
   nlinkingvars = 0;

   DECdecompSetConstoblock(decdecomp, constoblock, valid);
   assert(*valid);
   DECdecompSetVartoblock(decdecomp, vartoblock, valid);
   assert(*valid);

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], nconss) ); /*lint !e866*/
      nsubscipconss[i] = 0;
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

      assert(block > 0 && block <= nblocks+1);

      /* if variable belongs to a block */
      if( block <= nblocks )
      {
         SCIPdebugMessage("var %s in block %d.\n", SCIPvarGetName(var), block-1);
         subscipvars[block-1][nsubscipvars[block-1]] = var;
         ++(nsubscipvars[block-1]);
      }
      else /* variable is linking */
      {
         SCIPdebugMessage("var %s is linking.\n", SCIPvarGetName(var));
         assert(block == nblocks+1);
         linkingvars[nlinkingvars] = var;
         ++nlinkingvars;
      }
   }

   /* handle constraints */
   for( i = 0; i < nconss; ++i )
   {
      int block;
      SCIP_CONS* cons;

      cons = conss[i];
      assert(cons != NULL);
      if( !SCIPhashmapExists(constoblock, cons) )
         block = nblocks+1;
      else
      {
         block = (int)(size_t)SCIPhashmapGetImage(constoblock, cons); /*lint !e507*/
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
      SCIP_CALL( DECdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss, valid) );
      assert(*valid);
      DECdecompSetType(decdecomp, DEC_DECTYPE_BORDERED, valid);
      assert(*valid);
   }
   if( nlinkingvars > 0 )
   {
      SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars, valid) );
      assert(*valid);
      DECdecompSetType(decdecomp, DEC_DECTYPE_ARROWHEAD, valid);
      assert(*valid);
   }

   DECdecompSetNBlocks(decdecomp, nblocks);

   SCIP_CALL( DECdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss, valid) );
   assert(*valid);
   SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars, valid) );
   assert(*valid);

   SCIPfreeBufferArray(scip, &linkingconss);
   SCIPfreeBufferArray(scip, &linkingvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);

   for( i = 0; i < nblocks; ++i )
   {
     SCIPfreeBufferArray(scip, &subscipconss[i]);
     SCIPfreeBufferArray(scip, &subscipvars[i]);
   }

   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &subscipvars);

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
         SCIP_CALL( SCIPhashmapInsert(newconstoblock, decdecomp->subscipconss[b][c], (void*) (size_t) b) );
      }
   }
   /* transform all variables and put them into vartoblock */
   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      for( v = 0; v < decdecomp->nsubscipvars[b]; ++v )
      {
         SCIPdebugMessage("%d, %d: %s (%p, %s)\n", b, v, SCIPvarGetName(decdecomp->subscipvars[b][v]), decdecomp->subscipvars[b][v], SCIPvarIsTransformed(decdecomp->subscipvars[b][v])?"t":"o" );
         assert(decdecomp->subscipvars[b][v] != NULL);
         if( !SCIPvarIsTransformed(decdecomp->subscipvars[b][v]) )
         {
            SCIP_CALL( SCIPgetTransformedVar(scip, decdecomp->subscipvars[b][v], &newvar) );
            SCIP_CALL( SCIPcaptureVar(scip, newvar) );
            SCIP_CALL( SCIPreleaseVar(scip, &(decdecomp->subscipvars[b][v])) );
         }
         else
            newvar = decdecomp->subscipvars[b][v];
         assert(newvar != NULL);
         assert(SCIPvarIsTransformed(newvar));

         decdecomp->subscipvars[b][v] = newvar;
         SCIPdebugMessage("%d, %d: %s (%p, %s)\n", b, v, SCIPvarGetName(decdecomp->subscipvars[b][v]), decdecomp->subscipvars[b][v], SCIPvarIsTransformed(decdecomp->subscipvars[b][v])?"t":"o" );
         assert(decdecomp->subscipvars[b][v] != NULL);
         assert(!SCIPhashmapExists(newvartoblock, decdecomp->subscipvars[b][v]));
         SCIP_CALL( SCIPhashmapInsert(newvartoblock, decdecomp->subscipvars[b][v], (void*) (size_t) b) );
      }
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
      assert(decdecomp->linkingconss[c] != NULL);
   }

   /* transform all linking variables */
   for( v = 0; v < decdecomp->nlinkingvars; ++v )
   {
      SCIPdebugMessage("m, %d: %s (%p, %s)\n", v, SCIPvarGetName(decdecomp->linkingvars[v]), decdecomp->linkingvars[v], SCIPvarIsTransformed(decdecomp->linkingvars[v])?"t":"o" );
      assert(decdecomp->linkingvars[v] != NULL);

      if( !SCIPvarIsTransformed(decdecomp->linkingvars[v]) )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, decdecomp->linkingvars[v], &newvar) );
      }
      else
         newvar = decdecomp->linkingvars[v];
      assert(newvar != NULL);
      assert(SCIPvarIsTransformed(newvar));

      decdecomp->linkingvars[v] = newvar;
      SCIPdebugMessage("m, %d: %s (%p, %s)\n", v, SCIPvarGetName(decdecomp->linkingvars[v]), decdecomp->linkingvars[v], SCIPvarIsTransformed(decdecomp->linkingvars[v])?"t":"o" );
      assert(decdecomp->linkingvars[v] != NULL);
   }

   SCIPhashmapFree(&decdecomp->constoblock);
   decdecomp->constoblock = newconstoblock;
   SCIPhashmapFree(&decdecomp->vartoblock);
   decdecomp->vartoblock = newvartoblock;
   return SCIP_OKAY;
}
