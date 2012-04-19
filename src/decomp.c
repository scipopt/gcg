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

/* #define SCIP_DEBUG */

#include "decomp.h"
#include "pub_decomp.h"
#include "scip/scip.h"
#include "struct_decomp.h"

#include <assert.h>

/** Converts the DEC_DECTYPE enum to a string */
const char *DECgetStrType(
   DEC_DECTYPE type
   )
{
   const char * names[] = {"unknown","arrowhead", "staircase", "diagonal", "bordered"};
   return names[type];
}

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
   (*decomp)->stairlinkingvars = NULL;
   (*decomp)->nstairlinkingvars = NULL;
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
   for( i = 0; i < decomp->nblocks-1; ++i )
   {
      SCIPfreeMemoryArray(scip, &decomp->stairlinkingvars[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);
   SCIPfreeMemoryArray(scip, &decomp->stairlinkingvars);
   SCIPfreeMemoryArray(scip, &decomp->nstairlinkingvars);

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

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      //assert(nsubscipvars[b] > 0);
      decdecomp->nsubscipvars[b] = nsubscipvars[b];

      assert(subscipvars[b] != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->subscipvars[b], subscipvars[b], nsubscipvars[b]) ); /*lint !e866*/
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

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      assert(nsubscipconss[b] > 0);
      decdecomp->nsubscipconss[b] = nsubscipconss[b];

      assert(subscipconss[b] != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->subscipconss[b], subscipconss[b], nsubscipconss[b]) ); /*lint !e866*/
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
   int nlinkingvars            /**< number of linkingvars per block */
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

/** Copies the input stairlinkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetStairlinkingvars(
   SCIP* scip,                 /**< SCIP data structure */
   DECDECOMP* decdecomp,       /**< DECDECOMP data structure */
   SCIP_VAR*** stairlinkingvars,    /**< Linkingvars array  */
   int* nstairlinkingvars           /**< number of linkingvars per block */
   )
{
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(stairlinkingvars != NULL);
   assert(nstairlinkingvars != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->stairlinkingvars == NULL);
   assert(decdecomp->nstairlinkingvars == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->stairlinkingvars, decdecomp->nblocks-1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nstairlinkingvars, decdecomp->nblocks-1) );

   assert(decdecomp->stairlinkingvars != NULL);
   assert(decdecomp->nstairlinkingvars != NULL);

   for( b = 0; b < decdecomp->nblocks-1; ++b)
   {
      assert(nstairlinkingvars[b] > 0);
      decdecomp->nstairlinkingvars[b] = nstairlinkingvars[b];

      assert(stairlinkingvars[b] != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->stairlinkingvars[b], stairlinkingvars[b], nstairlinkingvars[b]) );
   }

   return SCIP_OKAY;
}

/** Returns the stairlinkingvars array of the given decdecomp structure */
SCIP_VAR***  DECdecdecompGetStairlinkingvars(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->stairlinkingvars;
}

/** Returns the nstairlinkingvars array of the given decdecomp structure */
int*  DECdecdecompGetNStairlinkingvars(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nstairlinkingvars != NULL );
   return decdecomp->nstairlinkingvars;
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

/** Sets the varindex hashmap of the given decdecomp structure */
void  DECdecdecompSetVarindex(
   DECDECOMP*    decdecomp,      /**< DECDECOMP data structure */
   SCIP_HASHMAP* varindex        /**< Varindex hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(varindex != NULL);
   decdecomp->varindex = varindex;
}

/** Returns the varindexk hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecdecompGetVarindex(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->varindex;
}

/** Sets the consindex hashmap of the given decdecomp structure */
void  DECdecdecompSetConsindex(
   DECDECOMP*    decdecomp,      /**< DECDECOMP data structure */
   SCIP_HASHMAP* consindex       /**< Consindex hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(consindex != NULL);
   decdecomp->consindex = consindex;
}

/** Returns the constoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecdecompGetConsindex(
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->consindex;
}

/** completely initializes decdecomp from the values of the hashmaps */
SCIP_RETCODE DECfillOutDecdecompFromHashmaps(
   SCIP*                 scip,               /**< SCIP data structure */
   DECDECOMP*            decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss              /**< number of constraints */
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

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   nlinkingconss = 0;
   nlinkingvars = 0;

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
      if( !SCIPhashmapExists(decdecomp->constoblock, cons) )
         block = nblocks+1;
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
      SCIP_CALL( DECdecdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss) );
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
   }
   if( nlinkingvars > 0 )
   {
      SCIP_CALL( DECdecdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars) );
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_ARROWHEAD);
   }
   SCIP_CALL( DECdecdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss) );
   SCIP_CALL( DECdecdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
   DECdecdecompSetVartoblock(decdecomp, vartoblock);
   DECdecdecompSetConstoblock(decdecomp, constoblock);

   SCIPfreeBufferArray(scip, &linkingconss);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);

   for( i = 0; i < nblocks; ++i )
   {
     SCIPfreeBufferArray(scip, &subscipconss[i]);
     SCIPfreeBufferArray(scip, &subscipvars[i]);
   }

   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &linkingvars);

   return SCIP_OKAY;
}

/** sets the detector for the given decdecomp structure */
void DECdecdecompSetDetector(
   DECDECOMP*    decdecomp,   /**< decdecomp instance */
   DEC_DETECTOR* detector     /**< detector data structure */
   )
{
   assert(decdecomp != NULL);
   assert(detector != NULL);

   decdecomp->detector = detector;
}

/** gets the detector for the given decdecomp structure */
DEC_DETECTOR* DECdecdecompGetDetector(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->detector;
}

/** transforms all constraints and variables, updating the arrays */
SCIP_RETCODE DECdecdecompTransform(
   SCIP*      scip,           /**< SCIP data structure */
   DECDECOMP* decdecomp       /**< decdecomp instance */
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
         SCIPdebugMessage("%d, %d: %s (%s)\n", b, c, SCIPconsGetName(decdecomp->subscipconss[b][c]), SCIPconsIsTransformed(decdecomp->subscipconss[b][c])?"t":"o" );
         assert(decdecomp->subscipconss[b][c] != NULL);
         decdecomp->subscipconss[b][c] = SCIPfindCons(scip, SCIPconsGetName(decdecomp->subscipconss[b][c]));
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
      SCIPdebugMessage("m, %d: %s (%s)\n", c, SCIPconsGetName(decdecomp->linkingconss[c]), SCIPconsIsTransformed(decdecomp->linkingconss[c])?"t":"o" );
      assert(decdecomp->linkingconss[c] != NULL);
      decdecomp->linkingconss[c] = SCIPfindCons(scip, SCIPconsGetName(decdecomp->linkingconss[c]));
      assert(decdecomp->linkingconss[c] != NULL);
   }

   /** @todo include Friederikes code */

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
