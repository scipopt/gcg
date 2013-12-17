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

/**@file   dec_isomorph.c
 * @ingroup DETECTORS
 * @brief  detector of problems that can be aggregated
 * @author Martin Bergner
 * @author Daniel Peters
 *
 * @todo try to merge decompositions with the same number of blocks if possible
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include "dec_isomorph.h"
#include "pub_decomp.h"
#include "cons_decomp.h"
#include "scip_misc.h"

#include "graph.hh"
#include "bliss_automorph.h"
#include "pub_gcgvar.h"
#include <cstring>
#include <cassert>
#include <algorithm>
#include "pub_bliss.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "isomorph"  /**< name of detector */
#define DEC_DESC                 "Detector for pricng problems suitable for aggregation" /**< description of detector*/
#define DEC_PRIORITY             700         /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'I'         /**< display character of detector */

#define DEC_ENABLED              TRUE        /**< should the detection be enabled */
#define DEC_SKIP                 TRUE        /**< should the detector be skipped if others found decompositions */

/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
   SCIP_RESULT result;                       /**< result pointer to indicate success or failure */
   int numofsol;                             /**< number of solutions */
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

   /** getter for the number of nodes */
   int getNNodes();

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
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &conssperm, SCIPgetNConss(scip)) );

}
struct_hook::~struct_hook()
{
   SCIPfreeMemoryArrayNull(scip_, &conssperm);
}
/** hook function to save the permutation of the graph */
static
void hook(
   void*                 user_param,         /**< data structure to save hashmaps with permutation */
   unsigned int          N,                  /**< number of permutations */
   const unsigned int*   aut                 /**< array of permutations */
   )
{
   int i;
   int nconss;
   SCIP_CONS** conss;
   AUT_HOOK* hook = (AUT_HOOK*) user_param;

   nconss = SCIPgetNConss(hook->getScip());
   assert(nconss == SCIPgetNConss(hook->getScip()));
   conss = SCIPgetConss(hook->getScip());

   for( i = 0; i < nconss; i++ )
   {
      assert(aut[i] < INT_MAX);
      if( (size_t) i != aut[i])
      {
         SCIPdebugMessage("%d <%s> <-> %d <%s>\n",i, SCIPconsGetName(conss[i]), aut[i], SCIPconsGetName(conss[aut[i]]));
         int index = MIN(i, aut[i]);
         if( hook->conssperm[i] != -1)
            index = MIN(index, hook->conssperm[i]);
         if( hook->conssperm[aut[i]] != -1 )
            index = MIN(index, hook->conssperm[aut[i]]);

         hook->conssperm[i] = index;
         hook->conssperm[aut[i]] = index;
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

static
SCIP_RETCODE freeMemory(SCIP* scip, /**< SCIP data structure */
AUT_COLOR* colorinfo /**< struct to save intermediate information */
)
{
   int i;
   AUT_COEF* scoef;
   AUT_VAR* svar;
   AUT_CONS* scons;

   for( i = 0; i < colorinfo->lenvarsarray; i++ )
   {
      svar = (AUT_VAR*) colorinfo->ptrarrayvars[i];
      delete svar;
   }
   for( i = 0; i < colorinfo->lenconssarray; i++ )
   {
      scons = (AUT_CONS*) colorinfo->ptrarrayconss[i];
      delete scons;
   }
   for( i = 0; i < colorinfo->lencoefsarray; i++ )
   {
      scoef = (AUT_COEF*) colorinfo->ptrarraycoefs[i];
      delete scoef;
   }

   SCIPfreeMemoryArray(scip, &colorinfo->ptrarraycoefs);
   SCIPfreeMemoryArray(scip, &colorinfo->ptrarrayconss);
   SCIPfreeMemoryArray(scip, &colorinfo->ptrarrayvars);
   return SCIP_OKAY;
}

/** set up a help structure for graph creation */
static
SCIP_RETCODE setuparrays(
   SCIP*                 scip,               /**< SCIP to compare */
   AUT_COLOR*            colorinfo,          /**< data structure to save intermediate data */
   SCIP_RESULT*          result              /**< result pointer to indicate success or failure */
   )
{
   int i;
   int j;
   int ncurvars;
   int nconss;
   int nvars;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   AUT_COEF* scoef;
   AUT_VAR* svar;
   AUT_CONS* scons;
   SCIP_Bool added;

   //allocate max n of coefarray, varsarray, and boundsarray in scip
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   allocMemory(scip, colorinfo, nconss, nvars);

   conss = SCIPgetConss(scip);
   vars = SCIPgetVars(scip);

   //save the properties of variables in a struct array and in a sorted pointer array
   for( i = 0; i < nvars; i++ )
   {
      svar = new AUT_VAR(scip, vars[i]);
      //add to pointer array iff it doesn't exist
      colorinfo->insert(svar, &added);
//      SCIPdebugMessage("%s color %d %d\n", SCIPvarGetName(vars[i]), colorinfo->get(*svar), colorinfo->color);
      //otherwise free allocated memory
      if( !added )
         delete svar;
   }

   //save the properties of constraints in a struct array and in a sorted pointer array
   for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
   {
      SCIP_Real* curvals;
      ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
      if( ncurvars == 0 )
         continue;
      scons = new AUT_CONS(scip, conss[i]);
      //add to pointer array iff it doesn't exist
      //SCIPdebugMessage("nconss %d %d\n", nconss, *result);
      colorinfo->insert(scons, &added);
  //    SCIPdebugMessage("%s color %d %d\n", SCIPconsGetName(conss[i]), colorinfo->get(*scons), colorinfo->color);
      //otherwise free allocated memory
      if( !added )
         delete scons;

      SCIP_CALL( SCIPallocMemoryArray(scip, &curvals, ncurvars) );
      SCIPgetValsXXX(scip, conss[i], curvals, ncurvars);
      //save the properties of variables of the constraints in a struct array and in a sorted pointer array
      for( j = 0; j < ncurvars; j++ )
      {
         scoef = new AUT_COEF(scip, curvals[j]);
         //test, whether the coefficient is not zero
         if( !SCIPisEQ(scip, scoef->getVal(), 0) )
         {
            //add to pointer array iff it doesn't exist
            colorinfo->insert(scoef, &added);
//            SCIPdebugMessage("%f color %d %d\n", scoef->getVal(), colorinfo->get(*scoef), colorinfo->color);
         }
         //otherwise free allocated memory
         if( !added )
            delete scoef;
      }
      SCIPfreeMemoryArray(scip, &curvals);
   }
   return SCIP_OKAY;
}

/** create a graph out of an array of scips */
static SCIP_RETCODE createGraph(
   SCIP*                 scip,               /**< SCIP to compare */
   AUT_COLOR             colorinfo,          /**< result pointer to indicate success or failure */
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
   SCIP_VAR** curvars;
   SCIP_Real* curvals;
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
      ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
      if( ncurvars == 0 )
         continue;

      AUT_CONS scons(scip, conss[i]);
      color = colorinfo.get(scons);

      if( color == -1 )
      {
         *result = SCIP_DIDNOTFIND;
         break;
      }

      h->add_vertex(color);
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
      h->add_vertex(colorinfo.getLenCons() + color);
      nnodes++;
   }
   //connecting the nodes with an additional node in the middle
   //it is necessary, since only nodes have colors
   for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
   {
      AUT_CONS scons(scip, conss[i]);
      ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
      if( ncurvars == 0 )
         continue;
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
      SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars);
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvals, ncurvars) );
      SCIPgetValsXXX(scip, conss[i], curvals, ncurvars);

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
         h->add_vertex(colorinfo.getLenCons() + colorinfo.getLenVar() + color);
         nnodes++;
         h->add_edge(i, nconss + nvars + z);
         h->add_edge(nconss + nvars + z, nconss + curvar);
         SCIPdebugMessage(
               "nz: c <%s> (id: %d, colour: %d) -> nz (id: %d) (value: %f, colour: %d) -> var <%s> (id: %d, colour: %d) \n",
               SCIPconsGetName(conss[i]), i, colorinfo.get(scons),
               nconss + nvars + z, scoef.getVal(),
               color + colorinfo.getLenCons() + colorinfo.getLenVar(),
               SCIPvarGetName(curvars[j]), nconss + curvar,
               colorinfo.get(svar) + colorinfo.getLenCons());
         z++;

      }

      SCIPfreeMemoryArray(scip, &curvals);
      SCIPfreeMemoryArray(scip, &curvars);
   }
   SCIPdebugMessage("Iteration 1: nnodes = %d, Cons = %d, Vars = %d\n", nnodes, colorinfo.getLenCons(), colorinfo.getLenVar());
   assert(*result == SCIP_SUCCESS && nnodes == h->get_nof_vertices());

   //free all allocated memory
   freeMemory(scip, &colorinfo);
   return SCIP_OKAY;
}

/** destructor of detector to free detector data (called when SCIP is exiting) */
static DEC_DECL_EXITDETECTOR(exitIsomorphism)
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

/** detection initialization function of detector (called before solving is about to begin) */
static DEC_DECL_INITDETECTOR(initIsomorphism)
{ /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->result = SCIP_SUCCESS;
   detectordata->numofsol = 10000;

   return SCIP_OKAY;
}

/** renumbers the permutations from 0 to n-1 and returns the number of permutations
 * @return the number of permutations
 */
int renumberPermutations(
   int*                  permutation,          /**< the permutation */
   int                   permsize              /**< size of the permutation */
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
   int*                  permutation,          /**< the permutation */
   int                   permsize              /**< size of the permutation */
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

/** detection function of detector */
static DEC_DECL_DETECTSTRUCTURE(detectIsomorphism)
{
   bliss::Graph graph;
   bliss::Stats bstats;
   AUT_HOOK *ptrhook;
   AUT_COLOR *colorinfo;
   *ndecdecomps = 0;
   *decdecomps = NULL;
   int nconss = SCIPgetNConss(scip);
   SCIP_CONS** masterconss;
   int i;
   int nmasterconss;
   colorinfo = new AUT_COLOR();

   SCIP_CALL( setuparrays(scip, colorinfo, &detectordata->result) );
   SCIP_CALL( createGraph(scip, *colorinfo, &graph, &detectordata->result) );

   ptrhook = new AUT_HOOK(FALSE, graph.get_nof_vertices(), scip);
   for( i = 0; i < nconss; i++ )
   {
      ptrhook->conssperm[i] = -1;
   }

   graph.find_automorphisms(bstats, hook, ptrhook);

   if( !ptrhook->getBool() )
      detectordata->result = SCIP_DIDNOTFIND;

   if( detectordata->result == SCIP_SUCCESS )
   {
      // assign to a permutation circle only one number
      collapsePermutation(ptrhook->conssperm, nconss);
      // renumbering from 0 to number of permutations
      renumberPermutations(ptrhook->conssperm, nconss);
      // create decomposition for every permutation
      SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, detectordata->numofsol) ); /*lint !e506*/
      SCIP_CALL( SCIPallocMemoryArray(scip, &masterconss, nconss) );

      nmasterconss = 0;
      for( i = 0; i < nconss; i++ )
      {
         if( -1 == ptrhook->conssperm[i] )
         {
            masterconss[nmasterconss] = SCIPgetConss(scip)[i];
            SCIPdebugMessage("%s\n", SCIPconsGetName(masterconss[nmasterconss]));
            nmasterconss++;
         }
      }
      SCIPdebugMessage("%d\n", nmasterconss);


      SCIP_CALL( DECcreateDecompFromMasterconss(scip, &((*decdecomps)[0]), masterconss, nmasterconss) );
      *ndecdecomps = 1;
      SCIPfreeMemoryArray(scip, &masterconss);

      for( i = 0; i < *ndecdecomps; i++ )
      {
         assert((*decdecomps)[i] != NULL);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found with %d blocks.\n", DECdecompGetNBlocks((*decdecomps)[i]));
      }
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
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

/** creates the handler for connected constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionIsomorphism(
   SCIP* scip /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create connected constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, DEC_SKIP, detectordata, detectIsomorphism, initIsomorphism, exitIsomorphism) );

   /* add connected constraint handler parameters */
   return SCIP_OKAY;
}
