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


typedef struct struct_cons AUT_CONS;
typedef struct struct_var AUT_VAR;
typedef struct struct_coef AUT_COEF;
typedef struct struct_hook AUT_HOOK;
typedef struct struct_colorinformation AUT_COLOR;

/** saves a constraint with its corresponding scip */
struct struct_cons
{
   SCIP* scip;                               /**< SCIP data structure */
   SCIP_CONS* cons;                          /**< pointer to SCIP constraint */

   /** constructor for the constraint struct */
   struct_cons(SCIP* scip, SCIP_CONS* scons);

   /** getter for the SCIP constraint */
   SCIP_CONS* getCons();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves a variable with its corresponding scip */
struct struct_var
{
   SCIP* scip;                               /**< SCIP data structure */
   SCIP_VAR* var;                            /**< pointer to SCIP variable */

   /** constructor for the variable struct */
   struct_var(SCIP* scip, SCIP_VAR* svar);

   /** getter for the SCIP variable */
   SCIP_VAR* getVar();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves a coefficient with its corresponding scip */
struct struct_coef
{
   SCIP* scip;                               /**< SCIP data structure */
   SCIP_Real val;                            /**< SCIP Real value */

   /** constructor for the coefficient struct */
   struct_coef(SCIP* scip, SCIP_Real val);

   /** getter for the SCIP Real value */
   SCIP_Real getVal();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves helping information for creating the graph */
struct struct_colorinformation
{
   int color;                                /**< color of the nodes of the graph */
   int lenconssarray;                        /**< size of ptrarrayconss */
   int lenvarsarray;                         /**< size of ptrarrayvars */
   int lencoefsarray;                        /**< size of ptrarraycoefs */
   int alloccoefsarray;                      /**< allocated size of ptrarraycoefs */
   void** ptrarraycoefs;                     /**< array of pointers to coefficient */
   void** ptrarrayvars;                      /**< array of pointers to variables */
   void** ptrarrayconss;                     /**< array of pointers to constraints */

   /** constructor for the  colorinformation struct */
   struct_colorinformation(int color,        /**< color of the nodes of the graph */
   int lenvars,                              /**< size of ptrarrayvars */
   int lenconss,                             /**< size of ptrarrayconss */
   int lencoefs                              /**< size of ptrarraycoefs */
   );

   /** insert a variable to its pointer array */
   void insert(AUT_VAR* svar, SCIP_Bool* added);

   /** insert a constraint to its pointer array */
   void insert(AUT_CONS* scons, SCIP_Bool* added);

   /** insert a coefficient to its pointer array */
   void insert(AUT_COEF* scoef, SCIP_Bool* added);

   /** getter for the length of the variable array */
   int getLenVar();

   /** getter for the length of the constraint array */
   int getLenCons();

   /** getter for the variable struct */
   int get(AUT_VAR svar);

   /** getter for the constraint struct */
   int get(AUT_CONS scons);

   /** getter for the coefficient struct */
   int get(AUT_COEF scoef);
};

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

SCIP_CONS* struct_cons::getCons()
{
   return this->cons;
}

SCIP* struct_cons::getScip()
{
   return this->scip;
}

SCIP_VAR* struct_var::getVar()
{
   return this->var;
}

SCIP* struct_var::getScip()
{
   return this->scip;
}

SCIP* struct_coef::getScip()
{
   return this->scip;
}

SCIP_Real struct_coef::getVal()
{
   return this->val;
}

SCIP* struct_hook::getScip()
{
   return this->scip;
}

static
SCIP_DECL_SORTPTRCOMP(sortptrcons);
static
SCIP_DECL_SORTPTRCOMP(sortptrvar);
static
SCIP_DECL_SORTPTRCOMP(sortptrval);

/** inserts a variable to the pointer array of colorinformation */
void struct_colorinformation::insert(
   AUT_VAR*                svar,              /**< variable which is to add */
   SCIP_Bool*              added             /**< true if a variable was added */
   )
{
   int pos;
   if( !SCIPsortedvecFindPtr(this->ptrarrayvars, sortptrvar, svar, this->lenvarsarray, &pos) )
   {
      SCIPsortedvecInsertPtr(this->ptrarrayvars, sortptrvar, svar, &this->lenvarsarray, NULL);
      *added = TRUE;
      this->color++;
   }
   else
      *added = FALSE;
}

/** inserts a constraint to the pointer array of colorinformation */
void struct_colorinformation::insert(
   AUT_CONS*             scons,              /**< constraint which is to add */
   SCIP_Bool*            added               /**< true if a constraint was added */
   )
{
   int pos;
   if( !SCIPsortedvecFindPtr(this->ptrarrayconss, sortptrcons, scons,
         this->lenconssarray, &pos) )
   {
      SCIPsortedvecInsertPtr(this->ptrarrayconss, sortptrcons, scons,
            &this->lenconssarray, NULL);
      *added = TRUE;
      this->color++;
   }
   else
      *added = FALSE;
}

/** inserts a coefficient to the pointer array of colorinformation */
void struct_colorinformation::insert(
   AUT_COEF*             scoef,              /**< coefficient which is to add */
   SCIP_Bool*            added               /**< true if a coefficient was added */
   )
{
   int pos;
   int size;
   if( !SCIPsortedvecFindPtr(this->ptrarraycoefs, sortptrval, scoef, this->lencoefsarray, &pos) )
   {
      size = SCIPcalcMemGrowSize(scoef->getScip(), this->alloccoefsarray+1);
      if( this->lencoefsarray % this->alloccoefsarray == 0 )
      {
         SCIP_CALL_ABORT( SCIPreallocMemoryArray(scip, &this->ptrarraycoefs, size) );
         this->alloccoefsarray = size;
      }

      SCIPsortedvecInsertPtr(this->ptrarraycoefs, sortptrval, scoef, &this->lencoefsarray, NULL);
      *added = TRUE;
      this->color++;
   }
   else
      *added = FALSE;

}

int struct_colorinformation::get(
   AUT_VAR               svar                /**< variable whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   found = SCIPsortedvecFindPtr(this->ptrarrayvars, sortptrvar, &svar, this->lenvarsarray, &pos);
   return found ? pos : -1;
}

int struct_colorinformation::get(
   AUT_CONS              scons               /**< constraint whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   found = SCIPsortedvecFindPtr(this->ptrarrayconss, sortptrcons, &scons, this->lenconssarray, &pos);
   return found ? pos : -1;
}

int struct_colorinformation::get(
   AUT_COEF              scoef               /**< coefficient whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   found = SCIPsortedvecFindPtr(this->ptrarraycoefs, sortptrval, &scoef, this->lencoefsarray, &pos);
   return found ? pos : -1;
}

/** constructor of the variable struct */
struct_var::struct_var(
   SCIP*                 scip_,              /**< SCIP data structure */
   SCIP_VAR*             svar                /**< SCIP variable */
   )
{
   scip = scip_;
   var = svar;
}

/** constructor of the constraint struct */
struct_cons::struct_cons(
   SCIP*                 scip_,              /**< SCIP data structure */
   SCIP_CONS*            scons               /**< SCIP constraint */
   )
{
   scip = scip_;
   cons = scons;
}

/** constructor of the coefficient struct */
struct_coef::struct_coef(
   SCIP*                 scip_,              /**< SCIP data structure */
   SCIP_Real             val_                /**< SCIP value */
   )
{
   scip = scip_;
   val = val_;
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
/** compare two values of two scips */
static
int comp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< value 1 to compare */
   SCIP_Real             val2                /**< value 2 to compare */
   )
{
   if( SCIPisLT(scip, val1, val2) )
      return -1;
   if( SCIPisGT(scip, val1, val2) )
      return 1;
   else
      return 0;
}

/** compare two constraints of two scips */
static
int comp(
   SCIP*                 scip,               /**< SCIP data structure */
   AUT_CONS*             cons1,              /**< constraint 1 to compare */
   AUT_CONS*             cons2               /**< constraint 2 to compare */
   )
{
   if( comp(scip, SCIPgetRhsXXX(scip, cons1->getCons()), SCIPgetRhsXXX(scip, cons2->getCons())) != 0 )
      return comp(scip, SCIPgetRhsXXX(scip, cons1->getCons()), SCIPgetRhsXXX(scip, cons2->getCons()));
   assert(SCIPisEQ(scip, SCIPgetRhsXXX(scip, cons1->getCons()), SCIPgetRhsXXX(scip, cons2->getCons())));

   if( comp(scip, SCIPgetLhsXXX(scip, cons1->getCons()), SCIPgetLhsXXX(scip, cons2->getCons())) != 0 )
      return comp(scip, SCIPgetLhsXXX(scip, cons1->getCons()), SCIPgetLhsXXX(scip, cons2->getCons()));
   assert(SCIPisEQ(scip, SCIPgetLhsXXX(scip, cons1->getCons()), SCIPgetLhsXXX(scip, cons2->getCons())));

   return strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons1->getCons())), SCIPconshdlrGetName(SCIPconsGetHdlr(cons2->getCons())));
}

/** compare two variables of two scips */
static
int comp(
   SCIP*                 scip,               /**< SCIP data structure */
   AUT_VAR*              var1,               /**< variable 1 to compare */
   AUT_VAR*              var2                /**< variable 2 to compare */
   )
{
   if( comp(scip, SCIPvarGetUbGlobal(var1->getVar()), SCIPvarGetUbGlobal(var2->getVar())) != 0 )
      return comp(scip, SCIPvarGetUbGlobal(var1->getVar()), SCIPvarGetUbGlobal(var2->getVar()));
   assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var1->getVar()), SCIPvarGetUbGlobal(var2->getVar())));

   if( comp(scip, SCIPvarGetLbGlobal(var1->getVar()), SCIPvarGetLbGlobal(var2->getVar())) != 0 )
      return comp(scip, SCIPvarGetLbGlobal(var1->getVar()), SCIPvarGetLbGlobal(var2->getVar()));
   assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var1->getVar()), SCIPvarGetLbGlobal(var2->getVar())));

   if( comp(scip, SCIPvarGetObj((var1->getVar())), SCIPvarGetObj(var2->getVar())) != 0 )
      return comp(scip, SCIPvarGetObj(var1->getVar()), SCIPvarGetObj(var2->getVar()));
   assert(SCIPisEQ(scip, SCIPvarGetObj(var1->getVar()), SCIPvarGetObj(var2->getVar())));

   if( SCIPvarGetType(var1->getVar()) < SCIPvarGetType(var2->getVar()) )
      return -1;
   if( SCIPvarGetType(var1->getVar()) > SCIPvarGetType(var2->getVar()) )
      return 1;
   return 0;
}

static
SCIP_DECL_SORTPTRCOMP(sortptrcons)
{
   AUT_CONS* aut1 = (AUT_CONS*) elem1;
   AUT_CONS* aut2 = (AUT_CONS*) elem2;
   return comp(aut1->getScip(), aut1, aut2);
}

static
SCIP_DECL_SORTPTRCOMP(sortptrvar)
{
   AUT_VAR* aut1 = (AUT_VAR*) elem1;
   AUT_VAR* aut2 = (AUT_VAR*) elem2;
   return comp(aut1->getScip(), aut1, aut2);
}

static
SCIP_DECL_SORTPTRCOMP(sortptrval)
{
   AUT_COEF* aut1 = (AUT_COEF*) elem1;
   AUT_COEF* aut2 = (AUT_COEF*) elem2;
   return comp(aut1->getScip(), aut1->getVal(), aut2->getVal());
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
         SCIPdebugMessage("cons <%s> <-> cons <%s>\n", SCIPconsGetName(conss[i]), SCIPconsGetName(conss[aut[i]]));
         SCIPdebugMessage("i <%d> <-> aut[i] <%d> \n", i, aut[i]);
         int index = MIN(i, aut[i]);
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
   AUT_COEF *scoef;
   AUT_VAR *svar;
   AUT_CONS *scons;
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

      scons = new AUT_CONS(scip, conss[i]);
      color = colorinfo.get(*scons);
      delete scons;
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
      svar = new AUT_VAR(scip, vars[i]);
      color = colorinfo.get(*svar);
      delete svar;

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
      scons = new AUT_CONS(scip, conss[i]);
      ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
      if( ncurvars == 0 )
         continue;
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
      SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars);
      SCIP_CALL( SCIPallocMemoryArray(scip, &curvals, ncurvars) );
      SCIPgetValsXXX(scip, conss[i], curvals, ncurvars);
      for( j = 0; j < ncurvars; j++ )
      {
         scoef = new AUT_COEF(scip, curvals[j]);
         color = colorinfo.get(*scoef);
         if( color == -1 )
         {
            *result = SCIP_DIDNOTFIND;
            delete scoef;
            break;
         }
         curvar = SCIPvarGetProbindex(curvars[j]);
         h->add_vertex(colorinfo.getLenCons() + colorinfo.getLenVar() + color);
         nnodes++;
         h->add_edge(i, nconss + nvars + z);
         h->add_edge(nconss + nvars + z, nconss + curvar);
         SCIPdebugMessage(
               "nz: c <%s> (id: %d, colour: %d) -> nz (id: %d) (value: %f, colour: %d) -> var <%s> (id: %d, colour: %d) \n",
               SCIPconsGetName(conss[i]), i, colorinfo.get(*scons),
               nconss + nvars + z, scoef->getVal(),
               color + colorinfo.getLenCons() + colorinfo.getLenVar(),
               SCIPvarGetName(curvars[j]), nconss + curvar,
               colorinfo.get(*svar) + colorinfo.getLenCons());
         z++;
         delete scoef;
      }
      delete scons;

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
   detectordata->numofsol = 100;

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
   int n = 0;
   int tmp = n;
   int i;
   for( i = 0; i < permsize; i++ )
   {
      if( permutation[i] != -1 && permutation[i] != n && permutation[i] != tmp )
      {
         n++;
         tmp = permutation[i];
         permutation[i] = n;
      }
      if( permutation[i] != -1 && permutation[i] != n && permutation[i] == tmp )
      {
         permutation[i] = n;
      }
      SCIPdebugMessage("%d\n", permutation[i]);
   }
   n++;
   return n;
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
      if( permutation[i] != -1 && permutation[i] < i )
      {
         tmp = permutation[i];
         permutation[i] = permutation[tmp];
      }
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
   colorinfo = new AUT_COLOR(0, 0, 0, 0);

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
