/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_decomp.h
 * @ingroup DECOMP
 * @ingroup PUBLICMETHODS
 * @brief  public methods for working with decomposition structures
 * @author Martin Bergner
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef __SCIP_PUB_DECOMP_H__
#define __SCIP_PUB_DECOMP_H__

#include "type_decomp.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

/** converts the DEC_DECTYPE enum to a string */
const char *DECgetStrType(
   DEC_DECTYPE type           /**< decomposition type */
   );

/** initializes the decdecomp structure to absolutely nothing */
SCIP_RETCODE DECdecdecompCreate(
   SCIP*       scip,          /**< SCIP data structure */
   DECDECOMP** decdecomp      /**< decdecomp instance */
   );

/** frees the decdecomp structure */
void DECdecdecompFree(
   SCIP*       scip,          /**< SCIP data structure */
   DECDECOMP** decdecomp      /**< decdecomp instance */
   );

/** sets the type of the decomposition */
void DECdecdecompSetType(
   DECDECOMP*  decdecomp,     /**< decdecomp instance */
   DEC_DECTYPE type           /**< type of the decomposition */
   );

/** gets the type of the decomposition */
DEC_DECTYPE DECdecdecompGetType(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** sets the number of blocks for decomposition */
void DECdecdecompSetNBlocks(
   DECDECOMP* decdecomp,      /**< decdecomp instance */
   int        nblocks         /**< number of blocks for decomposition */
   );

/** gets the number of blocks for decomposition */
int DECdecdecompGetNBlocks(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** copies the input subscipvars array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetSubscipvars(
   SCIP*       scip,          /**< SCIP data structure */
   DECDECOMP*  decdecomp,     /**< decdecomp instance */
   SCIP_VAR*** subscipvars,   /**< subscipvars array  */
   int*        nsubscipvars   /**< number of subscipvars per block */
   );

/** returns the subscipvars array of the given decdecomp structure */
SCIP_VAR*** DECdecdecompGetSubscipvars(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** returns the nsubscipvars array of the given decdecomp structure */
int* DECdecdecompGetNSubscipvars(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** copies the input subscipconss array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetSubscipconss(
   SCIP*        scip,         /**< SCIP data structure */
   DECDECOMP*   decdecomp,    /**< decdecomp instance */
   SCIP_CONS*** subscipconss, /**< subscipconss array  */
   int*         nsubscipconss /**< number of subscipconss per block */
   );

/** returns the subscipconss array of the given decdecomp structure */
SCIP_CONS*** DECdecdecompGetSubscipconss(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** returns the nsubscipconss array of the given decdecomp structure */
int*  DECdecdecompGetNSubscipconss(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** copies the input linkingconss array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetLinkingconss(
   SCIP*       scip,          /**< SCIP data structure */
   DECDECOMP*  decdecomp,     /**< decdecomp instance */
   SCIP_CONS** linkingconss,  /**< linkingconss array  */
   int         nlinkingconss  /**< number of linkingconss per block */
   );

/** returns the linkingconss array of the given decdecomp structure */
SCIP_CONS**  DECdecdecompGetLinkingconss(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** returns the nlinkingconss array of the given decdecomp structure */
int  DECdecdecompGetNLinkingconss(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** copies the input linkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecdecompSetLinkingvars(
   SCIP*      scip,           /**< SCIP data structure */
   DECDECOMP* decdecomp,      /**< decdecomp instance */
   SCIP_VAR** linkingvars,    /**< linkingvars array  */
   int        nlinkingvars    /**< number of linkingvars per block */
   );

/** returns the linkingvars array of the given decdecomp structure */
SCIP_VAR** DECdecdecompGetLinkingvars(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** returns the nlinkingvars array of the given decdecomp structure */
int DECdecdecompGetNLinkingvars(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** sets the vartoblock hashmap of the given decdecomp structure */
void DECdecdecompSetVartoblock(
   DECDECOMP*    decdecomp,   /**< decdecomp instance */
   SCIP_HASHMAP* vartoblock   /**< Vartoblock hashmap */
   );

/** returns the vartoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP* DECdecdecompGetVartoblock(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** sets the constoblock hashmap of the given decdecomp structure */
void DECdecdecompSetConstoblock(
   DECDECOMP*    decdecomp,   /**< decdecomp instance */
   SCIP_HASHMAP* constoblock  /**< Constoblock hashmap */
   );

/** returns the constoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP* DECdecdecompGetConstoblock(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );


/** completely initializes decdecomp from the values of the hashmaps */
SCIP_RETCODE DECfillOutDecdecompFromHashmaps(
   SCIP*         scip,        /**< SCIP data structure */
   DECDECOMP*    decdecomp,   /**< decdecomp instance */
   SCIP_HASHMAP* vartoblock,  /**< variable to block hashmap */
   SCIP_HASHMAP* constoblock, /**< constraint to block hashmap */
   int           nblocks,     /**< number of blocks */
   SCIP_VAR**    vars,        /**< variable array */
   int           nvars,       /**< number of variables */
   SCIP_CONS**   conss,       /**< constraint array */
   int           nconss       /**< number of constraints */
   );

/** sets the detector for the given decdecomp structure */
void DECdecdecompSetDetector(
   DECDECOMP*    decdecomp,   /**< decdecomp instance */
   DEC_DETECTOR* detector     /**< detector data structure */
   );

/** gets the detector for the given decdecomp structure */
DEC_DETECTOR* DECdecdecompGetDetector(
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

/** transforms all constraints and variables, updating the arrays */
SCIP_RETCODE DECdecdecompTransform(
   SCIP*      scip,           /**< SCIP data structure */
   DECDECOMP* decdecomp       /**< decdecomp instance */
   );

#ifdef __cplusplus
}
#endif
#endif
