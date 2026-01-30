/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   pricingcb.c
 * @ingroup OTHER_CFILES
 * @brief  methods for pricing callback
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "gcg/gcg.h"

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/pub_message.h"

#include "gcg/pricingcb.h"
#include "gcg/struct_pricingcb.h"

/* default parameter settings for the pricing callbacks */
#define GCG_DEFAULT_ENABLED        FALSE
#define GCG_DEFAULT_EXCLUSIVE      FALSE

/** compares two pricing callbacks w.r.t their priority */
SCIP_DECL_SORTPTRCOMP(GCGpricingcbComp)
{  /*lint --e{715}*/
   GCG_PRICINGCB* pricingcb1 = (GCG_PRICINGCB*)elem1;
   GCG_PRICINGCB* pricingcb2 = (GCG_PRICINGCB*)elem2;

   assert(pricingcb1 != NULL);
   assert(pricingcb2 != NULL);

   return pricingcb2->priority - pricingcb1->priority; /* prefer higher priorities */
}

/** comparison method for sorting pricing callbacks w.r.t their name */
SCIP_DECL_SORTPTRCOMP(GCGpricingcbCompName)
{
   return strcmp(GCGpricingcbGetName((GCG_PRICINGCB*)elem1), GCGpricingcbGetName((GCG_PRICINGCB*)elem2));
}

/** internal method for creating a pricing callback structure */
static
SCIP_RETCODE doPricingcbCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_PRICINGCB**       pricingcb,          /**< pointer to the pricing callback data structure */
   const char*           name,               /**< name of the pricing callback */
   const char*           desc,               /**< description of the pricing callback */
   int                   priority,           /**< priority of the the pricing callback */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree)), /**< destructor of the pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit)), /**< initialize the pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit)), /**< deinitialize the pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol)),/**< solving process initialization method of the pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol)),/**< solving process deinitialization method of the pricing callback */
   GCG_DECL_PRICINGCBPREPRICING((*pricingcbprepricing)),/**< pre-pricing method of the pricing callback */
   GCG_DECL_PRICINGCBPOSTPRICING((*pricingcbpostpricing)),/**< post-pricing method of the pricing callback */
   GCG_PRICINGCBDATA*  pricingcbdata         /**< pricing callback data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(pricingcb != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(pricingcb) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*pricingcb)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*pricingcb)->desc, desc, strlen(desc)+1) );
   (*pricingcb)->priority = priority;
   (*pricingcb)->pricingcbfree = pricingcbfree;
   (*pricingcb)->pricingcbinit = pricingcbinit;
   (*pricingcb)->pricingcbexit = pricingcbexit;
   (*pricingcb)->pricingcbinitsol = pricingcbinitsol;
   (*pricingcb)->pricingcbexitsol = pricingcbexitsol;
   (*pricingcb)->pricingcbprepricing = pricingcbprepricing;
   (*pricingcb)->pricingcbpostpricing = pricingcbpostpricing;
   (*pricingcb)->pricingcbdata = pricingcbdata;
   SCIP_CALL( SCIPcreateClock(scip, &(*pricingcb)->setuptime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*pricingcb)->pricingcbclock) );
   (*pricingcb)->nprepricingcalls = 0;
   (*pricingcb)->npostpricingcalls = 0;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricingcb/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of the pricing callback <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, paramname, paramdesc,
                  &(*pricingcb)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricingcb/%s/enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname,
        "are the methods of this pricing callback enabled?", &(*pricingcb)->enabled, FALSE,
        GCG_DEFAULT_ENABLED, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricingcb/%s/exclusive", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname,
        "are the methods of this pricing callback executed exclusively (only takes effect if highest priority callback)?",
        &(*pricingcb)->exclusive, FALSE, GCG_DEFAULT_EXCLUSIVE, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a pricing callback */
SCIP_RETCODE GCGpricingcbCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB**       pricingcb,          /**< pointer to the pricing callback data structure */
   const char*           name,               /**< name of the pricing callback */
   const char*           desc,               /**< description of the pricing callback */
   int                   priority,           /**< priority of the the pricing callback */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree)), /**< destructor of the pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit)), /**< initialize the pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit)), /**< deinitialize the pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol)),/**< solving process initialization method of the pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol)),/**< solving process deinitialization method of the pricing callback */
   GCG_DECL_PRICINGCBPREPRICING((*pricingcbprepricing)),/**< pre-pricing method of the pricing callback */
   GCG_DECL_PRICINGCBPOSTPRICING((*pricingcbpostpricing)),/**< post-pricing method of the pricing callback */
   GCG_PRICINGCBDATA*  pricingcbdata         /**< pricing callback data */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   scip = GCGgetMasterprob(gcg);
   SCIP_CALL_FINALLY( doPricingcbCreate(scip, pricingcb, name, desc, priority,
         pricingcbfree, pricingcbinit, pricingcbexit, pricingcbinitsol, pricingcbexitsol,
         pricingcbprepricing, pricingcbpostpricing, pricingcbdata), (void)GCGpricingcbFree(gcg, pricingcb) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of the pricing callback */
SCIP_RETCODE GCGpricingcbFree(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB**       pricingcb           /**< pointer to the pricing callback data structure */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);
   assert(*pricingcb != NULL);

   scip = GCGgetMasterprob(gcg);

   /* call destructor of the pricing callback */
   if( (*pricingcb)->pricingcbfree != NULL )
   {
      SCIP_CALL( (*pricingcb)->pricingcbfree(gcg, *pricingcb) );
   }

   SCIPfreeClock(scip, &(*pricingcb)->pricingcbclock);
   SCIPfreeClock(scip, &(*pricingcb)->setuptime);
   BMSfreeMemoryArray(&(*pricingcb)->name);
   BMSfreeMemoryArray(&(*pricingcb)->desc);
   BMSfreeMemory(pricingcb);

   return SCIP_OKAY;
}

/** initializes the pricing callback */
SCIP_RETCODE GCGpricingcbInit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   SCIP* scip;
   SCIP_Bool misc_resetstat;

   assert(gcg != NULL);
   assert(pricingcb != NULL);

   scip = GCGgetMasterprob(gcg);
   SCIP_CALL( SCIPgetBoolParam(scip, "misc/resetstat", &misc_resetstat) );

   if( misc_resetstat )
   {
      SCIPresetClock(scip, pricingcb->setuptime);
      SCIPresetClock(scip, pricingcb->pricingcbclock);

      pricingcb->nprepricingcalls = 0;
      pricingcb->npostpricingcalls = 0;
   }

   if( pricingcb->pricingcbinit != NULL )
   {
      /* start timing */
      SCIPstartClock(scip, pricingcb->setuptime);

      SCIP_CALL( pricingcb->pricingcbinit(gcg, pricingcb) );

      /* stop timing */
      SCIPstopClock(scip, pricingcb->setuptime);
   }

   return SCIP_OKAY;
}

/** calls exit method of the pricing callback */
SCIP_RETCODE GCGpricingcbExit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);

   scip = GCGgetMasterprob(gcg);

   if( pricingcb->pricingcbexit != NULL )
   {
      /* start timing */
      SCIPstartClock(scip, pricingcb->setuptime);

      SCIP_CALL( pricingcb->pricingcbexit(gcg, pricingcb) );

      /* stop timing */
      SCIPstopClock(scip, pricingcb->setuptime);
   }

   return SCIP_OKAY;
}

/** informs pricing callback that the branch and bound process is being started */
SCIP_RETCODE GCGpricingcbInitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);

   scip = GCGgetMasterprob(gcg);

   /* call solving process initialization method of the pricing callback */
   if( pricingcb->pricingcbinitsol != NULL )
   {
      /* start timing */
      SCIPstartClock(scip, pricingcb->setuptime);

      SCIP_CALL( pricingcb->pricingcbinitsol(gcg, pricingcb) );

      /* stop timing */
      SCIPstopClock(scip, pricingcb->setuptime);
   }

   return SCIP_OKAY;
}

/** informs pricing callback that the branch and bound process data is being freed */
SCIP_RETCODE GCGpricingcbExitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);

   scip = GCGgetMasterprob(gcg);

   /* call solving process deinitialization method of pricing callback */
   if( pricingcb->pricingcbexitsol != NULL )
   {
      /* start timing */
      SCIPstartClock(scip, pricingcb->setuptime);

      SCIP_CALL( pricingcb->pricingcbexitsol(gcg, pricingcb) );

      /* stop timing */
      SCIPstopClock(scip, pricingcb->setuptime);
   }

   return SCIP_OKAY;
}

/** calls pre-pricing method of the pricing callback */
SCIP_RETCODE GCGpricingcbPrepricing(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_PRICER*          pricer,             /**< the pricer */
   GCG_PRICETYPE         type,               /**< the type of pricing, either redcost or farkas */
   SCIP_Bool*            abort,              /**< should the pricing be aborted? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);
   assert(pricingcb->pricingcbprepricing != NULL);
   assert(pricer != NULL);
   assert(abort != NULL);
   assert(result != NULL);

   /* the pricing callback should only be called if it is enabled */
   assert(pricingcb->enabled);

   scip = GCGgetMasterprob(gcg);
   (*abort) = FALSE;
   (*result) = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "executing the pre-pricing method of pricing callback <%s>\n", pricingcb->name);

   /* start timing */
   SCIPstartClock(scip, pricingcb->pricingcbclock);

   SCIP_CALL( pricingcb->pricingcbprepricing(gcg, pricingcb, pricer, type, abort, result) );

   /* stop timing */
   SCIPstopClock(scip, pricingcb->pricingcbclock);

   /* evaluate result */
   if( (*result) != SCIP_DIDNOTRUN
      && (*result) != SCIP_SUCCESS )
   {
      SCIPerrorMessage("pre-pricing method of pricing callback <%s> returned invalid result <%d>\n",
         pricingcb->name, (*result));
      return SCIP_INVALIDRESULT;
   }

   pricingcb->nprepricingcalls++;

   return SCIP_OKAY;
}

/** calls post-pricing method of the pricing callback */
SCIP_RETCODE GCGpricingcbPostpricing(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_PRICER*          pricer,             /**< the pricer */
   GCG_PRICETYPE         type,               /**< the type of pricing, either redcost or farkas */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   assert(pricingcb != NULL);
   assert(pricingcb->pricingcbpostpricing != NULL);
   assert(pricer != NULL);
   assert(result != NULL);

   scip = GCGgetMasterprob(gcg);

   /* the pricing callback should only be called if it is enabled */
   assert(pricingcb->enabled);

   (*result) = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "executing the post-pricing method of pricing callback <%s>\n", pricingcb->name);

   /* start timing */
   SCIPstartClock(scip, pricingcb->pricingcbclock);

   SCIP_CALL( pricingcb->pricingcbpostpricing(gcg, pricingcb, pricer, type, result) );

   /* stop timing */
   SCIPstopClock(scip, pricingcb->pricingcbclock);

   /* evaluate result */
   if( (*result) != SCIP_DIDNOTRUN
      && (*result) != SCIP_SUCCESS )
   {
      SCIPerrorMessage("post-pricing method of pricing callback <%s> returned invalid result <%d>\n",
         pricingcb->name, (*result));
      return SCIP_INVALIDRESULT;
   }

   pricingcb->npostpricingcalls++;

   return SCIP_OKAY;
}

/** gets user data of the pricing callback */
GCG_PRICINGCBDATA* GCGpricingcbGetData(
   GCG_PRICINGCB*      pricingcb          /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->pricingcbdata;
}

/** sets user data of the pricing callback; user has to free old data in advance! */
void GCGpricingcbSetData(
   GCG_PRICINGCB*       pricingcb,          /**< pricing callback */
   GCG_PRICINGCBDATA*   pricingcbdata       /**< new pricing callback user data */
   )
{
   assert(pricingcb != NULL);

   pricingcb->pricingcbdata = pricingcbdata;
}

/* new callback setter methods */

/** sets priority of the pricing callback */
void GCGpricingcbSetPriority(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   int                   priority            /**< new priority of the pricing callback */
   )
{
   assert(pricingcb != NULL);
   pricingcb->priority = priority;
}

/** sets destructor callback of the pricing callback */
void GCGpricingcbSetFree(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree))  /**< destructor of the pricing callback */
   )
{
   assert(pricingcb != NULL);

   pricingcb->pricingcbfree = pricingcbfree;
}

/** sets initialization callback of the pricing callback */
void GCGpricingcbSetInit(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit))  /**< initialize the pricing callback */
   )
{
   assert(pricingcb != NULL);

   pricingcb->pricingcbinit = pricingcbinit;
}

/** sets deinitialization callback of the pricing callback */
void GCGpricingcbSetExit(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit))  /**< deinitialize the pricing callback */
   )
{
   assert(pricingcb != NULL);

   pricingcb->pricingcbexit = pricingcbexit;
}

/** sets solving process initialization callback of the pricing callback */
void GCGpricingcbSetInitsol(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol))/**< solving process initialization callback of the pricing callback */
   )
{
   assert(pricingcb != NULL);

   pricingcb->pricingcbinitsol = pricingcbinitsol;
}

/** sets solving process deinitialization callback of pricing callback */
void GCGpricingcbSetExitsol(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol))/**< solving process deinitialization callback of the pricing callback */
   )
{
   assert(pricingcb != NULL);

   pricingcb->pricingcbexitsol = pricingcbexitsol;
}

/** gets name of the pricing callback */
const char* GCGpricingcbGetName(
   GCG_PRICINGCB*       pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->name;
}

/** gets description of the pricing callback */
const char* GCGpricingcbGetDesc(
   GCG_PRICINGCB*       pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->desc;
}

/** gets priority of the pricing callback */
int GCGpricingcbGetPriority(
   GCG_PRICINGCB*       pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->priority;
}

/** gets the number of times the pre-pricing method of the pricing callback plugin was called */
SCIP_Longint GCGpricingcbGetNPrepricingCalls(
   GCG_PRICINGCB*       pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->nprepricingcalls;
}

/** gets the number of times the post-pricing method of the pricing callback plugin was called */
SCIP_Longint GCGpricingcbGetNPostpricingCalls(
   GCG_PRICINGCB*       pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->npostpricingcalls;
}

/** gets time in seconds used by this pricing callback for setting up */
SCIP_Real GCGpricingcbGetSetupTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return SCIPgetClockTime(GCGgetMasterprob(gcg), pricingcb->setuptime);
}

/** gets time in seconds used in this pricing callback */
SCIP_Real GCGpricingcbGetTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return SCIPgetClockTime(GCGgetMasterprob(gcg), pricingcb->pricingcbclock);
}

/** sets the enabled flag of the pricing callback method */
void GCGpricingcbSetEnabled(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_Bool             enabled             /**< flag to indicate whether the pricing callback is enabled */
   )
{
   assert(pricingcb != NULL);

   pricingcb->enabled = enabled;
}

/** sets the exclusive flag of the pricing callback plugin method */
void GCGpricingcbSetExclusive(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_Bool             exclusive           /**< flag to indicate whether the pricing callback plugin is executed exclusively */
   )
{
   assert(pricingcb != NULL);

   pricingcb->exclusive = exclusive;
}

/** returns whether the pricing callback is enabled */
SCIP_Bool GCGpricingcbIsEnabled(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->enabled;
}

/** returns whether the methods of this pricing callback should be executed exclusively */
SCIP_Bool GCGpricingcbIsExclusive(
   GCG_PRICINGCB*       pricingcb           /**< pricing callback */
   )
{
   assert(pricingcb != NULL);

   return pricingcb->exclusive;
}
