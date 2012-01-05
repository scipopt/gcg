/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
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
/** initializes the decdecomp structure to absolutely nothing */
SCIP_RETCODE DECdecdecompCreate(
   SCIP* scip,           /**< Pointer to the SCIP instance */
   DECDECOMP** decomp    /**< Pointer to the decdecomp instance */
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
void DECdecdecompFree(
   SCIP* scip,           /**< Pointer to the SCIP instance */
   DECDECOMP** decdecomp /**< Pointer to the decdecomp instance */
   )
{
   DECDECOMP* decomp;
   int i;

   assert( scip!= NULL );
   assert( decdecomp != NULL);
   decomp = *decdecomp;

   assert(decomp != NULL);

   for( i = 0; i < decomp->nblocks; ++i )
   {
      SCIPfreeMemoryArray(scip, &decomp->subscipvars[i]);
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

   SCIPfreeMemoryArrayNull(scip, &decomp->linkingconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);

   SCIPfreeMemory(scip, decdecomp);
}

/** sets the type of the decomposition */
void DECdecdecompSetType(
   DECDECOMP* decdecomp, /**< Pointer to the decdecomp instance */
   DEC_DECTYPE type      /**< type of the decomposition */
   )
{
   assert(decdecomp != NULL);
   decdecomp->type = type;
}

/** gets the type of the decomposition */
DEC_DECTYPE DECdecdecompGetType(
   DECDECOMP* decdecomp  /**< Pointer to the decdecomp instance */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->type;
}


/** sets the number of blocks for decomposition */
void DECdecdecompSetNBlocks(
   DECDECOMP* decdecomp, /**< Pointer to the decdecomp instance */
   int nblocks           /**< number of blocks for decomposition */
   )
{
   assert(decdecomp != NULL);
   assert(nblocks >= 0);
   decdecomp->nblocks = nblocks;
}

/** gets the number of blocks for decomposition */
int DECdecdecompGetNBlocks(
   DECDECOMP* decdecomp  /**< Pointer to the decdecomp instance */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->nblocks;
}

/** Copies the input subscipvars array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetSubscipvars(
   SCIP* scip,                 /**< SCIP data structure */
   DECDECOMP* decdecomp,       /**< DECDECOMP data structure */
   SCIP_VAR*** subscipvars,    /**< Subscipvars array  */
   int* nsubscipvars           /**< number of subscipvars per block */
   )
{
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipvars != NULL);
   assert(nsubscipvars != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->subscipvars == NULL);
   assert(decdecomp->nsubscipvars == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipvars, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipvars, decdecomp->nblocks) );

   assert(decdecomp->subscipvars != NULL);
   assert(decdecomp->nsubscipvars != NULL);

   for( b = 0; b < decdecomp->nblocks; ++b)
   {
      assert(nsubscipvars[b] > 0);
      decdecomp->nsubscipvars[b] = nsubscipvars[b];

      assert(subscipvars[b] != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->subscipvars[b], subscipvars[b], nsubscipvars[b]) );
   }

   return SCIP_OKAY;
}

/** Returns the subscipvars array of the given decdecomp structure */
SCIP_VAR***  DECdecdecompGetSubscipvars(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->subscipvars;
}

/** Returns the nsubscipvars array of the given decdecomp structure */
int*  DECdecdecompGetNSubscipvars(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->nsubscipvars;
}

/** Copies the input subscipconss array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetSubscipconss(
   SCIP* scip,                 /**< SCIP data structure */
   DECDECOMP* decdecomp,       /**< DECDECOMP data structure */
   SCIP_CONS*** subscipconss,    /**< Subscipconss array  */
   int* nsubscipconss           /**< number of subscipconss per block */
   )
{
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipconss != NULL);
   assert(nsubscipconss != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->subscipconss == NULL);
   assert(decdecomp->nsubscipconss == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipconss, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipconss, decdecomp->nblocks) );

   assert(decdecomp->subscipconss != NULL);
   assert(decdecomp->nsubscipconss != NULL);

   for( b = 0; b < decdecomp->nblocks; ++b)
   {
      assert(nsubscipconss[b] > 0);
      decdecomp->nsubscipconss[b] = nsubscipconss[b];

      assert(subscipconss[b] != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->subscipconss[b], subscipconss[b], nsubscipconss[b]) );
   }

   return SCIP_OKAY;
}

/** Returns the subscipconss array of the given decdecomp structure */
SCIP_CONS***  DECdecdecompGetSubscipconss(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->subscipconss;
}

/** Returns the nsubscipconss array of the given decdecomp structure */
int*  DECdecdecompGetNSubscipconss(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->nsubscipconss;
}

/** Copies the input linkingconss array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetLinkingconss(
   SCIP* scip,                 /**< SCIP data structure */
   DECDECOMP* decdecomp,       /**< DECDECOMP data structure */
   SCIP_CONS** linkingconss,    /**< Linkingconss array  */
   int nlinkingconss          /**< number of linkingconss per block */
   )
{
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingconss != NULL);
   assert(nlinkingconss > 0);

   assert(decdecomp->linkingconss == NULL);
   assert(decdecomp->nlinkingconss == 0);

   decdecomp->nlinkingconss = nlinkingconss;

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingconss, linkingconss, nlinkingconss) );

   return SCIP_OKAY;
}

/** Returns the linkingconss array of the given decdecomp structure */
SCIP_CONS**  DECdecdecompGetLinkingconss(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->linkingconss;
}

/** Returns the nlinkingconss array of the given decdecomp structure */
int  DECdecdecompGetNLinkingconss(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nlinkingconss >= 0);
   return decdecomp->nlinkingconss;
}

/** Copies the input linkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetLinkingvars(
   SCIP* scip,                 /**< SCIP data structure */
   DECDECOMP* decdecomp,       /**< DECDECOMP data structure */
   SCIP_VAR** linkingvars,    /**< Linkingvars array  */
   int nlinkingvars          /**< number of linkingvars per block */
   )
{
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingvars != NULL);
   assert(nlinkingvars > 0);

   assert(decdecomp->linkingvars == NULL);
   assert(decdecomp->nlinkingvars == 0);

   decdecomp->nlinkingvars = nlinkingvars;

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingvars, linkingvars, nlinkingvars) );

   return SCIP_OKAY;
}

/** Returns the linkingvars array of the given decdecomp structure */
SCIP_VAR**  DECdecdecompGetLinkingvars(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->linkingvars;
}

/** Returns the nlinkingvars array of the given decdecomp structure */
int  DECdecdecompGetNLinkingvars(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nlinkingvars >= 0);
   return decdecomp->nlinkingvars;
}

/** Sets the vartoblock hashmap of the given decdecomp structure */
void  DECdecdecompSetVartoblock(
   DECDECOMP*    decdecomp,      /**< DECDECOMP data structure */
   SCIP_HASHMAP* vartoblock      /**< Vartoblock hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);
   decdecomp->vartoblock = vartoblock;
}

/** Returns the vartoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecdecompGetVartoblock(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->vartoblock;
}

/** Sets the constoblock hashmap of the given decdecomp structure */
void  DECdecdecompSetConstoblock(
   DECDECOMP*    decdecomp,      /**< DECDECOMP data structure */
   SCIP_HASHMAP* constoblock      /**< Constoblock hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(constoblock != NULL);
   decdecomp->constoblock = constoblock;
}

/** Returns the constoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecdecompGetConstoblock(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->constoblock;
}