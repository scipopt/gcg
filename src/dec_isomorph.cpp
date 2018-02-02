/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   dec_isomorph.cpp
 * @ingroup DETECTORS
 * @brief  detector for pricing problems that can be aggregated (uses bliss)
 * @author Martin Bergner
 * @author Daniel Peters
 *
 * This detector finds subproblems that can be aggregated thus reducing the symmetry of the problem using color preserving
 * automorphisms and bliss.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include "dec_isomorph.h"
#include "pub_decomp.h"
#include "cons_decomp.h"
#include "scip_misc.h"

#include "graph.hh"
#include "pub_gcgvar.h"
#include <cstring>
#include <cassert>
#include <algorithm>

#include "pub_bliss.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "isomorph"  /**< name of detector */
#define DEC_DESC                 "Detector for pricing problems suitable for aggregation" /**< description of detector*/
#define DEC_PRIORITY             100         /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'I'         /**< display character of detector */

#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#define DEC_SKIP                 TRUE        /**< should the detector be skipped if others found decompositions */

#define DEFAULT_MAXDECOMPS       1           /**< default maximum number of decompositions */

/*
 * Data structures
 */

/** detector data */
struct DEC_DetectorData
{
   SCIP_RESULT result;                       /**< result pointer to indicate success or failure */
   int maxdecomps;                           /**< maximum number of decompositions */

};

typedef struct struct_hook AUT_HOOK;

/** saves information of the permutation */
struct struct_hook
{
   SCIP_Bool aut;                            /**< true if there is an automorphism */
   unsigned int n;                           /**< number of permutations */
   SCIP* scip;                               /**< scip to search for automorphisms */
   int* conssperm;                           /**< permutations of conss*/

   /** constructor for the hook struct*/
   struct_hook(SCIP_Bool aut,  /**< true if there is an automorphism */
      unsigned int       n,                  /**< number of permutations */
      SCIP*              scip                /**< scip to search for automorphisms */
   );

   ~struct_hook();
   /** getter for the bool aut */
   SCIP_Bool getBool();

   /** setter for the bool aut */
   void setBool(SCIP_Bool aut);

   /** getter for the SCIP */
   SCIP* getScip();
};

SCIP* struct_hook::getScip()
{
   return this->scip;
}


/** constructor of the hook struct */
struct_hook::struct_hook(
   SCIP_Bool             aut_,               /**< true if there is an automorphism */
   unsigned int          n_,                 /**< number of permutations */
   SCIP*                 scip_               /**< array of scips to search for automorphisms */
   ) : conssperm(NULL)
{
   aut = aut_;
   n = n_;
   scip = scip_;
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &conssperm, SCIPgetNConss(scip)) ); /*lint !e666*/

}
struct_hook::~struct_hook()
{   /*lint -esym(1540,struct_hook::conssperm) */
   SCIPfreeMemoryArrayNull(scip, &conssperm);
   scip = NULL;
}
/** hook function to save the permutation of the graph */
static
void fhook(
   void*                 user_param,         /**< data structure to save hashmaps with permutation */
   unsigned int          N,                  /**< number of permutations */
   const unsigned int*   aut                 /**< array of permutations */
   )
{ /*lint -e715*/
   int i;
   int nconss;
   SCIP_CONS** conss;
   AUT_HOOK* hook = (AUT_HOOK*) user_param;
   int auti;
   int ind;

   nconss = SCIPgetNConss(hook->getScip());
   assert(nconss == SCIPgetNConss(hook->getScip()));
   conss = SCIPgetConss(hook->getScip());

   for( i = 0; i < nconss; i++ )
   {
      assert(aut[i] < INT_MAX);
      if( (size_t) i != aut[i])
      {
         auti = (int) aut[i];

         SCIPdebugMessage("%d <%s> <-> %d <%s>\n", i, SCIPconsGetName(conss[i]), auti, SCIPconsGetName(conss[auti]));

         ind = MIN(i, auti);

         if( hook->conssperm[i] != -1)
            ind = MIN(ind, hook->conssperm[i]);
         if( hook->conssperm[auti] != -1 )
            ind = MIN(ind, hook->conssperm[auti]);

         hook->conssperm[i] = ind;
         hook->conssperm[auti] = ind;
         hook->setBool(TRUE);
      }
   }
}

static
SCIP_RETCODE allocMemory(
    SCIP*                scip,               /**< SCIP data structure */
    AUT_COLOR*           colorinfo,          /**< struct to save intermediate information */
    int                  nconss,             /**< number of constraints */
    int                  nvars               /**< number of variables */
    )
{
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarraycoefs, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarrayvars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarrayconss, nconss) );
   colorinfo->alloccoefsarray = nvars;
   return SCIP_OKAY;
}

/** destructor for colorinfo */
static
void freeMemory(
   SCIP*                 scip,               /**< SCIP data structure */
   AUT_COLOR*            colorinfo           /**< struct to save intermediate information */
)
{
   int i;

   for( i = 0; i < colorinfo->lenvarsarray; i++ )
   {
      AUT_VAR* svar = (AUT_VAR*) colorinfo->ptrarrayvars[i];
      delete svar;
   }
   for( i = 0; i < colorinfo->lenconssarray; i++ )
   {
      AUT_CONS* scons = (AUT_CONS*) colorinfo->ptrarrayconss[i];
      delete scons;
   }
   for( i = 0; i < colorinfo->lencoefsarray; i++ )
   {
      AUT_COEF* scoef = (AUT_COEF*) colorinfo->ptrarraycoefs[i];
      delete scoef;
   }

   SCIPfreeMemoryArray(scip, &colorinfo->ptrarraycoefs);
   SCIPfreeMemoryArray(scip, &colorinfo->ptrarrayconss);
   SCIPfreeMemoryArray(scip, &colorinfo->ptrarrayvars);
}

/** set up a help structure for graph creation */
static
SCIP_RETCODE setupArrays(
   SCIP*                 scip,               /**< SCIP to compare */
   AUT_COLOR*            colorinfo,          /**< data structure to save intermediate data */
   SCIP_RESULT*          result              /**< result pointer to indicate success or failure */
   )
{ /*lint -esym(593,scoef) */
   int i;
   int j;
   int nconss;
   int nvars;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   AUT_COEF* scoef;
   AUT_CONS* scons;
   SCIP_Bool added;

   //allocate max n of coefarray, varsarray, and boundsarray in scip
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( allocMemory(scip, colorinfo, nconss, nvars) );

   conss = SCIPgetConss(scip);
   vars = SCIPgetVars(scip);

   //save the properties of variables in a struct array and in a sorted pointer array
   for( i = 0; i < nvars; i++ )
   {
      AUT_VAR* svar = new AUT_VAR(scip, vars[i]);
      //add to pointer array iff it doesn't exist
      SCIP_CALL( colorinfo->insert(svar, &added) );
      SCIPdebugMessage("%s color %d %d\n", SCIPvarGetName(vars[i]), colorinfo->get(*svar), colorinfo->color);
      //otherwise free allocated memory
      if( !added )
         delete svar;
   }

   //save the properties of constraints in a struct array and in a sorted pointer array
   for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
   {
      SCIP_Real* curvals = NULL;
      int ncurvars = GCGconsGetNVars(scip, conss[i]);
      if( ncurvars == 0 )
         continue;
      scons = new AUT_CONS(scip, conss[i]);
      //add to pointer array iff it doesn't exist
      //SCIPdebugMessage("nconss %d %d\n", nconss, *result);
      SCIP_CALL( colorinfo->insert(scons, &added) );
      SCIPdebugMessage("%s color %d %d\n", SCIPconsGetName(conss[i]), colorinfo->get(*scons), colorinfo->color);
      //otherwise free allocated memory
      if( !added )
         delete scons;

      SCIP_CALL( SCIPallocBufferArray(scip, &curvals, ncurvars) );
      SCIP_CALL( GCGconsGetVals(scip, conss[i], curvals, ncurvars) );
      //save the properties of variables of the constraints in a struct array and in a sorted pointer array
      for( j = 0; j < ncurvars; j++ )
      {
         added = FALSE;
         scoef = new AUT_COEF(scip, curvals[j]);
         //test, whether the coefficient is not zero
         if( !SCIPisZero(scip, scoef->getVal()) )
         {
            //add to pointer array iff it doesn't exist
            SCIP_CALL( colorinfo->insert(scoef, &added) );
            SCIPdebugMessage("%f color %d %d\n", scoef->getVal(), colorinfo->get(*scoef), colorinfo->color);
         }
         //otherwise free allocated memory
         if( !added )
            delete scoef;

      }
      SCIPfreeBufferArray(scip, &curvals);
   }
   return SCIP_OKAY;
}

/** create a graph out of an array of scips */
static
SCIP_RETCODE createGraph(
   SCIP*                 scip,               /**< SCIP to compare */
   AUT_COLOR             colorinfo,          /**< data structure to save intermediate data */
   bliss::Graph*         graph,              /**< graph needed for discovering isomorphism */
   SCIP_RESULT*          result              /**< result pointer to indicate success or failure */
   )
{
   int i;
   int j;
   int z;
   int nvars;
   int nconss;
   int ncurvars;
   int curvar;
   int color;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_VAR** curvars = NULL;
   SCIP_Real* curvals = NULL;
   unsigned int nnodes;
   nnodes = 0;
   //building the graph out of the arrays
   bliss::Graph* h = graph;
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   conss = SCIPgetConss(scip);
   vars = SCIPgetVars(scip);
   z = 0;

   //add a node for every constraint
   for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
   {
      ncurvars = GCGconsGetNVars(scip, conss[i]);

      AUT_CONS scons(scip, conss[i]);
      color = colorinfo.get(scons);

      if( color == -1 )
      {
         *result = SCIP_DIDNOTFIND;
         break;
      }

      assert(color >= 0);
      (void)h->add_vertex((unsigned int) color);
      nnodes++;
   }
   //add a node for every variable
   for( i = 0; i < nvars && *result == SCIP_SUCCESS; i++ )
   {
      AUT_VAR svar(scip, vars[i]);
      color = colorinfo.get(svar);

      if( color == -1 )
      {
         *result = SCIP_DIDNOTFIND;
         break;
      }
      (void) h->add_vertex((unsigned int) (colorinfo.getLenCons() + color));
      nnodes++;
   }
   //connecting the nodes with an additional node in the middle
   //it is necessary, since only nodes have colors
   for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
   {
      AUT_CONS scons(scip, conss[i]);
      ncurvars = GCGconsGetNVars(scip, conss[i]);
      if( ncurvars == 0 )
         continue;
      SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
      SCIP_CALL( GCGconsGetVars(scip, conss[i], curvars, ncurvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &curvals, ncurvars) );
      SCIP_CALL( GCGconsGetVals(scip, conss[i], curvals, ncurvars) );

      for( j = 0; j < ncurvars; j++ )
      {
         AUT_COEF scoef(scip, curvals[j]);
         AUT_VAR svar(scip, curvars[j]);

         color = colorinfo.get(scoef);
         if( color == -1 )
         {
            *result = SCIP_DIDNOTFIND;

            break;
         }
         curvar = SCIPvarGetProbindex(curvars[j]);
         (void) h->add_vertex((unsigned int) (colorinfo.getLenCons() + colorinfo.getLenVar() + color)); /*lint !e864 */
         nnodes++;
         h->add_edge((unsigned int)i, (unsigned int) (nconss + nvars + z));
         h->add_edge((unsigned int) (nconss + nvars + z), (unsigned int) (nconss + curvar));
         SCIPdebugMessage(
               "nz: c <%s> (id: %d, colour: %d) -> nz (id: %d) (value: %f, colour: %d) -> var <%s> (id: %d, colour: %d) \n",
               SCIPconsGetName(conss[i]), i, colorinfo.get(scons),
               nconss + nvars + z, scoef.getVal(),
               color + colorinfo.getLenCons() + colorinfo.getLenVar(), /*lint !e864 */
               SCIPvarGetName(curvars[j]), nconss + curvar,
               colorinfo.get(svar) + colorinfo.getLenCons());  /*lint !e864 */
         z++;

      }

      SCIPfreeBufferArray(scip, &curvals);
      SCIPfreeBufferArray(scip, &curvars);

   }
   SCIPdebugMessage("Iteration 1: nnodes = %ud, Cons = %d, Vars = %d\n", nnodes, colorinfo.getLenCons(), colorinfo.getLenVar()); /*lint !e864 */
   assert(*result == SCIP_SUCCESS && nnodes == h->get_nof_vertices());

   //free all allocated memory
   freeMemory(scip, &colorinfo);
   return SCIP_OKAY;
}

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeIsomorph)
{ /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detector initialization method (called after problem was transformed) */
static
DEC_DECL_INITDETECTOR(detectorInitIsomorph)
{ /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** renumbers the permutations from 0 to n-1 and returns the number of permutations
 * @return the number of permutations
 */
int renumberPermutations(
   int*                  permutation,        /**< the permutation */
   int                   permsize            /**< size of the permutation */
)
{
   // renumbering from 0 to number of permutations
   int nperms = -1;

   for( int i = 0; i < permsize; i++ )
   {
      SCIPdebugMessage("%d: %d -> ", i, permutation[i]);
      if( permutation[i] == -1 )
      {
         SCIPdebugPrintf("%d\n", permutation[i]);
         continue;
      }

      if( permutation[i] > nperms && permutation[permutation[i]] > nperms )
      {
         nperms++;
         permutation[i] = nperms;
      }
      else
      {
         permutation[i] = permutation[permutation[i]];
      }
      SCIPdebugPrintf("%d\n", permutation[i]);
   }

   return nperms+1;
}

/** collapses the permutation, if possible */
void collapsePermutation(
   int*                  permutation,        /**< the permutation */
   int                   permsize            /**< size of the permutation */
)
{
   int tmp = 0;
   // assign to a permutation circle only one number
   for( int i = 0; i < permsize; i++ )
   {
      if( permutation[i] != -1 && permutation[i] != i )
      {
         tmp = permutation[i];
         permutation[i] = permutation[tmp];
      }
      SCIPdebugMessage("%d %d\n",i, permutation[i]);

   }
}

/** reorder such that the best permutation is represented by 0, the second best by 1, etc. */
SCIP_RETCODE reorderPermutations(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  permutation,        /**< the permutation */
   int                   permsize,           /**< size of the permutation */
   int                   nperms              /**< number of permutations */
)
{
   int i;
   int* count;
   int* order;
   int* invorder;

   assert(scip != NULL);
   assert(permutation != NULL);
   assert(permsize > 0);
   assert(nperms > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &count, nperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &order, nperms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &invorder, nperms) );
   BMSclearMemoryArray(count, nperms);
   BMSclearMemoryArray(order, nperms);
   BMSclearMemoryArray(invorder, nperms);

   /* initialize order array that will give a mapping from new to old representatives */
   for( i = 0; i < nperms; ++i )
   {
      order[i] = i;
   }

   /* count sizes of orbits */
   for( i = 0; i < permsize; ++i )
   {
      if( permutation[i] >= 0 )
      {
         count[permutation[i]] += 1;

         SCIPdebugMessage("permutation[i] = %d; count %d\n", permutation[i], count[permutation[i]]);
      }
   }

   /* sort count and order array */
   SCIPsortDownIntInt(count, order, nperms);

#ifdef SCIP_DEBUG

   for( i = 0; i < nperms; ++i )
   {
      SCIPdebugMessage("count[%d] = %d, order[%d] = %d\n", i, count[i], i, order[i]);
   }
#endif

   /* compute invorder array that gives a mapping from old to new representatives */
   for( i = 0; i < nperms; ++i )
   {
      invorder[order[i]] = i;
      SCIPdebugMessage("invorder[%d] = %d\n", order[i], invorder[order[i]]);
   }

   SCIPdebugMessage("Best permutation with orbit of size %d, best %d\n", count[0], order[0]);

   /* update representatives of constraints */
   for( i = 0; i < permsize; ++i )
   {
      if( permutation[i] >= 0 )
         permutation[i] = invorder[permutation[i]];
   }

   SCIPfreeBufferArray(scip, &count);
   SCIPfreeBufferArray(scip, &order);
   SCIPfreeBufferArray(scip, &invorder);

   return SCIP_OKAY;
}

/** detector structure detection method, tries to detect a structure in the problem */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectIsomorph)
{ /*lint -esym(429,ptrhook)*/
   bliss::Graph graph;
   bliss::Stats bstats;
   AUT_HOOK *ptrhook;
   AUT_COLOR *colorinfo;

   *ndecdecomps = 0;
   *decdecomps = NULL;

   int nconss = SCIPgetNConss(scip);
   int i;
   int unique;

   colorinfo = new AUT_COLOR();
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting aggregatable structure: ");
   SCIP_CALL( setupArrays(scip, colorinfo, &detectordata->result) );
   SCIP_CALL( createGraph(scip, *colorinfo, &graph, &detectordata->result) );

   ptrhook = new AUT_HOOK(FALSE, graph.get_nof_vertices(), scip);
   for( i = 0; i < nconss; i++ )
   {
      ptrhook->conssperm[i] = -1;
   }

   graph.find_automorphisms(bstats, fhook, ptrhook);

   if( !ptrhook->getBool() )
      detectordata->result = SCIP_DIDNOTFIND;

   if( detectordata->result == SCIP_SUCCESS )
   {
      int nperms;
      DEC_DECOMP* newdecomp;
      int nmasterconss;
      SCIP_CONS** masterconss = NULL;
      int p;

      // assign to a permutation circle only one number
      collapsePermutation(ptrhook->conssperm, nconss);
      // renumbering from 0 to number of permutations
      nperms = renumberPermutations(ptrhook->conssperm, nconss);

      // filter decomposition with largest orbit
      SCIP_CALL( reorderPermutations(scip, ptrhook->conssperm, nconss, nperms) );

      SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, MIN(detectordata->maxdecomps, nperms)) ); /*lint !e506*/

      int pos = 0;
      for( p = 0; p < nperms && pos < detectordata->maxdecomps; ++p )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &masterconss, nconss) );

         nmasterconss = 0;
         for( i = 0; i < nconss; i++ )
         {
            if( p != ptrhook->conssperm[i] )
            {
               masterconss[nmasterconss] = SCIPgetConss(scip)[i];
               SCIPdebugMessage("%s\n", SCIPconsGetName(masterconss[nmasterconss]));
               nmasterconss++;
            }
         }
         SCIPdebugMessage("%d\n", nmasterconss);

         if( nmasterconss < SCIPgetNConss(scip) )
         {
            SCIP_CALL( DECcreateDecompFromMasterconss(scip, &((*decdecomps)[pos]), masterconss, nmasterconss) );

            SCIPfreeMemoryArray(scip, &masterconss);
         }
         else
         {
            SCIPfreeMemoryArray(scip, &masterconss);

            continue;
         }


         SCIP_CALL( DECcreatePolishedDecomp(scip, (*decdecomps)[pos], &newdecomp) );
         if( newdecomp != NULL )
         {
            SCIP_CALL( DECdecompFree(scip, &((*decdecomps)[pos])) );
            (*decdecomps)[pos] = newdecomp;
         }

         ++pos;
      }
      *ndecdecomps = pos;

      if( *ndecdecomps > 0 )
      {
         unique = DECfilterSimilarDecompositions(scip, *decdecomps, *ndecdecomps);
      }
      else
      {
         unique = *ndecdecomps;
      }

      for( p = unique; p < *ndecdecomps; ++p )
      {
         SCIP_CALL( DECdecompFree(scip, &((*decdecomps)[p])) );
         (*decdecomps)[p] = NULL;
      }

      *ndecdecomps = unique;

      if( *ndecdecomps > 0 )
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );
      }

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "found %d decompositions.\n", *ndecdecomps);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "not found.\n");
   }

   if( *ndecdecomps == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, decdecomps);
   }

   delete ptrhook;
   delete colorinfo;
   *result = detectordata->result;

   return SCIP_OKAY;
}

/*
 * detector specific interface methods
 */

/** creates the handler for isomorph subproblems and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorIsomorphism(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP,
      detectordata, detectorDetectIsomorph, detectorFreeIsomorph, detectorInitIsomorph, NULL) );

   /* add isomorph constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/isomorph/maxdecomps",
      "Maximum number of decompositions", &detectordata->maxdecomps, FALSE,
      DEFAULT_MAXDECOMPS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
