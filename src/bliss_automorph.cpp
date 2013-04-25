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

/**@file    bliss_automorph.cpp
 * @brief   automorphism recognition of SCIPs
 * @author  Daniel Peters
 * @author  Martin Bergner
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph.hh"
#include "bliss_automorph.h"
#include "scip_misc.h"
#include "scip/scip.h"
#include "pub_gcgvar.h"
#include "relax_gcg.h"
#include "scip/cons_linear.h"

#include <cstring>


typedef struct Aut_Cons AUT_CONS;
typedef struct Aut_Var AUT_VAR;
typedef struct Aut_Coef AUT_COEF;
typedef struct struct_hook AUT_HOOK;
typedef struct struct_colorinformation AUT_COLOR;

/** saves a constraint with its corresponding scip */
struct Aut_Cons
{
   SCIP* scip;                              /**< SCIP data structure */
   SCIP_CONS* cons;                         /**< pointer to SCIP constraint */

   /** constructor for the constraint struct */
   Aut_Cons( SCIP* scip, SCIP_CONS* scons );

   /** getter for the SCIP constraint */
   SCIP_CONS* getCons();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves a variable with its corresponding scip */
struct Aut_Var
{
   SCIP* scip;                              /**< SCIP data structure */
   SCIP_VAR* var;                           /**< pointer to SCIP variable */

   /** constructor for the variable struct */
   Aut_Var( SCIP* scip, SCIP_VAR* svar );

   /** getter for the SCIP variable */
   SCIP_VAR* getVar();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves a coefficient with its corresponding scip */
struct Aut_Coef
{
   SCIP* scip;                              /**< SCIP data structure */
   SCIP_Real val;                           /**< SCIP Real value */

   /** constructor for the coefficient struct */
   Aut_Coef( SCIP* scip, SCIP_Real val );

   /** getter for the SCIP Real value */
   SCIP_Real getVal();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves helping information for creating the graph */
struct struct_colorinformation
{
   int                  color;              /**< color of the nodes of the graph */
   int                  lenconssarray;      /**< size of ptrarrayconss */
   int                  lenvarsarray;       /**< size of ptrarrayvars */
   int                  lencoefsarray;      /**< size of ptrarraycoefs */
   void**               ptrarraycoefs;      /**< array of pointers to coefficient */
   void**               ptrarrayvars;       /**< array of pointers to variables */
   void**               ptrarrayconss;      /**< array of pointers to constraints */

   /** default constructor */
   struct_colorinformation();

   /** constructor for the  colorinformation struct */
   struct_colorinformation(
      int               color,              /**< color of the nodes of the graph */
      int               lenvars,            /**< size of ptrarrayvars */
      int               lenconss,           /**< size of ptrarrayconss */
      int               lencoefs            /**< size of ptrarraycoefs */
      );

   /** insert a variable to its pointer array */
   void insert( AUT_VAR* svar, SCIP_Bool* added);

   /** insert a constraint to its pointer array */
   void insert( AUT_CONS* scons, SCIP_Bool* added);

   /** insert a coefficient to its pointer array */
   void insert( AUT_COEF* scoef, SCIP_Bool* added);

   /** getter for the length of the variable array */
   int getLenVar();

   /** getter for the length of the constraint array */
   int getLenCons();

   /** getter for the variable struct */
   int get( AUT_VAR svar );

   /** getter for the constraint struct */
   int get( AUT_CONS scons );

   /** getter for the coefficient struct */
   int get( AUT_COEF scoef );
};

/** saves information of the permutation */
struct struct_hook
{
   SCIP_Bool aut;                           /**< true if there is an automorphism */
   unsigned int n;                          /**< number of permutations */
   SCIP_HASHMAP* varmap;                    /**< hashmap for permutated variables */
   SCIP_HASHMAP* consmap;                   /**< hashmap for permutated constraints */
   SCIP** scips;                            /**< array of scips to search for automorphisms */

   /** constructor for the hook struct*/
   struct_hook(
      SCIP_HASHMAP* varmap,                 /**< hashmap for permutated variables */
      SCIP_HASHMAP* consmap,                /**< hashmap for permutated constraints */
      SCIP_Bool aut,                        /**< true if there is an automorphism */
      unsigned int n,                       /**< number of permutations */
      SCIP** scips                          /**< array of scips to search for automorphisms */
      );

   /** getter for the bool aut */
   SCIP_Bool getBool();

   /** setter for the bool aut */
   void setBool(SCIP_Bool aut);

   /** getter for the number of nodes */
   int getNNodes();

   /** getter for the variables hashmap */
   SCIP_HASHMAP* getVarHash();

   /** getter for the constraints hashmap */
   SCIP_HASHMAP* getConsHash();

   /** getter for the SCIPs */
   SCIP** getScips();
};

SCIP_CONS* Aut_Cons::getCons()
{
   return this->cons;
}

SCIP* Aut_Cons::getScip()
{
   return this->scip;
}

SCIP_VAR* Aut_Var::getVar()
{
   return this->var;
}

SCIP* Aut_Var::getScip()
{
   return this->scip;
}

SCIP* Aut_Coef::getScip()
{
   return this->scip;
}

SCIP_Real Aut_Coef::getVal()
{
   return this->val;
}

void struct_hook::setBool( SCIP_Bool aut_ )
{
   aut = aut_;
}

SCIP_Bool struct_hook::getBool()
{
   return this->aut;
}

int struct_hook::getNNodes()
{
   return this->n;
}

SCIP_HASHMAP* struct_hook::getVarHash()
{
   return this->varmap;
}

SCIP_HASHMAP* struct_hook::getConsHash()
{
   return this->consmap;
}

SCIP** struct_hook::getScips()
{
   return this->scips;
}

static
SCIP_DECL_SORTPTRCOMP(sortptrcons);
static
SCIP_DECL_SORTPTRCOMP(sortptrvar);
static
SCIP_DECL_SORTPTRCOMP(sortptrval);

/** inserts a variable to the pointer array of colorinformation */
void struct_colorinformation::insert(
   AUT_VAR* svar,                           /**< variable which is to add */
   SCIP_Bool* added                         /**< true if a variable was added */
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
   AUT_CONS* scons,                         /**< constraint which is to add */
   SCIP_Bool* added                         /**< true if a constraint was added */
   )
{
   int pos;
   if( !SCIPsortedvecFindPtr(this->ptrarrayconss, sortptrcons, scons, this->lenconssarray, &pos) )
   {
      SCIPsortedvecInsertPtr(this->ptrarrayconss, sortptrcons, scons, &this->lenconssarray, NULL);
      *added = TRUE;
      this->color++;
   }
   else
      *added = FALSE;
}

/** inserts a coefficient to the pointer array of colorinformation */
void struct_colorinformation::insert(
   AUT_COEF* scoef,                         /**< coefficient which is to add */
   SCIP_Bool* added                         /**< true if a coefficient was added */
   )
{
   int pos;
   if( !SCIPsortedvecFindPtr(this->ptrarraycoefs, sortptrval, scoef, this->lencoefsarray, &pos) )
   {
      SCIPsortedvecInsertPtr(this->ptrarraycoefs, sortptrval, scoef, &this->lencoefsarray, NULL);
      *added = TRUE;
      this->color++;
   }
   else
      *added = FALSE;
}

int struct_colorinformation::getLenVar()
{
   return lenvarsarray;
}

int struct_colorinformation::getLenCons()
{
   return lenconssarray;
}

int struct_colorinformation::get(
   AUT_VAR svar                             /**< variable whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   found = SCIPsortedvecFindPtr(this->ptrarrayvars, sortptrvar, &svar, this->lenvarsarray, &pos);
   return found ? pos : -1;
}

int struct_colorinformation::get(
   AUT_CONS scons                           /**< constraint whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   found = SCIPsortedvecFindPtr(this->ptrarrayconss, sortptrcons, &scons, this->lenconssarray, &pos);
   return found ? pos : -1;
}

int struct_colorinformation::get(
   AUT_COEF scoef                           /**< coefficient whose pointer you want */
   )
{
   int pos;
   SCIP_Bool found;
   found = SCIPsortedvecFindPtr(this->ptrarraycoefs, sortptrval, &scoef, this->lencoefsarray, &pos);
   return found ? pos : -1;
}

/** constructor of the variable struct */
Aut_Var::Aut_Var(
   SCIP* scip_,                             /**< SCIP data structure */
   SCIP_VAR* svar                           /**< SCIP variable */
   )
{
   scip = scip_;
   var = svar;
}

/** constructor of the constraint struct */
Aut_Cons::Aut_Cons(
   SCIP* scip_,                             /**< SCIP data structure */
   SCIP_CONS* scons                         /**< SCIP constraint */
   )
{
   scip = scip_;
   cons = scons;
}

/** constructor of the coefficient struct */
Aut_Coef::Aut_Coef(
   SCIP* scip_,                             /**< SCIP data structure */
   SCIP_Real val_                           /**< SCIP value */
   )
{
   scip = scip_;
   val = val_;
}

/** default constructor */
struct_colorinformation::struct_colorinformation()
 : color(0), lenconssarray(0), lenvarsarray(0), lencoefsarray(0),
ptrarraycoefs(NULL), ptrarrayvars(NULL), ptrarrayconss(NULL)
{

}

/** constructor of the color struct */
struct_colorinformation::struct_colorinformation(
   int color_,                              /**< color of the nodes of the graph */
   int lenvars,                             /**< length of ptrvarsarray */
   int lenconss,                            /**< length of ptrconsarray */
   int lencoefs                             /**< length of ptrcoefsarray */
   )
: color(color_), lenconssarray(lenconss), lenvarsarray(lenvars), lencoefsarray(lencoefs),
ptrarraycoefs(NULL),  ptrarrayvars(NULL), ptrarrayconss(NULL)
{

}

/** constructor of the hook struct */
struct_hook::struct_hook(
   SCIP_HASHMAP* varmap_,                   /**< hashmap of permutated variables */
   SCIP_HASHMAP* consmap_,                  /**< hahsmap of permutated constraints */
   SCIP_Bool aut_,                          /**< true if there is an automorphism */
   unsigned int n_,                         /**< number of permutations */
   SCIP** scips_                            /**< array of scips to search for automorphisms */
   )
{
   aut = aut_;
   n = n_;
   consmap = consmap_;
   varmap = varmap_;
   scips = scips_;
}

/** compare two values of two scips */
static
int comp(
   SCIP* scip,                              /**< SCIP data structure */
   SCIP_Real val1,                          /**< value 1 to compare */
   SCIP_Real val2                           /**< value 2 to compare */
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
   SCIP* scip,                              /**< SCIP data structure */
   AUT_CONS* cons1,                         /**< constraint 1 to compare */
   AUT_CONS* cons2                          /**< constraint 2 to compare */
)
{
   if(comp(scip, SCIPgetRhsXXX(cons1->getScip(), cons1->getCons()), SCIPgetRhsXXX(cons2->getScip(), cons2->getCons())) != 0)
      return comp(scip, SCIPgetRhsXXX(cons1->getScip(), cons1->getCons()), SCIPgetRhsXXX(cons2->getScip(), cons2->getCons()));
   assert(SCIPisEQ(scip, SCIPgetRhsXXX(cons1->getScip(), cons1->getCons()), SCIPgetRhsXXX(cons2->getScip(), cons2->getCons())));

   if(comp(scip, SCIPgetLhsXXX(cons1->getScip(), cons1->getCons()), SCIPgetLhsXXX(cons2->getScip(), cons2->getCons())) != 0)
      return comp(scip, SCIPgetLhsXXX(cons1->getScip(), cons1->getCons()), SCIPgetLhsXXX(cons2->getScip(), cons2->getCons()));
   assert(SCIPisEQ(scip, SCIPgetLhsXXX(cons1->getScip(), cons1->getCons()), SCIPgetLhsXXX(cons2->getScip(), cons2->getCons())));

   return strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons1->getCons())), SCIPconshdlrGetName(SCIPconsGetHdlr(cons2->getCons())));
}

/** compare two variables of two scips */
static
int comp(
   SCIP* scip,                              /**< SCIP data structure */
   AUT_VAR* var1,                           /**< variable 1 to compare */
   AUT_VAR* var2                            /**< variable 2 to compare */
)
{
   if(comp(scip, SCIPvarGetUbGlobal(var1->getVar()), SCIPvarGetUbGlobal(var2->getVar())) != 0)
      return comp(scip, SCIPvarGetUbGlobal(var1->getVar()), SCIPvarGetUbGlobal(var2->getVar()));
   assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var1->getVar()), SCIPvarGetUbGlobal(var2->getVar())));

   if(comp(scip, SCIPvarGetLbGlobal(var1->getVar()), SCIPvarGetLbGlobal(var2->getVar())) != 0)
      return comp(scip, SCIPvarGetLbGlobal(var1->getVar()), SCIPvarGetLbGlobal(var2->getVar()));
   assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var1->getVar()), SCIPvarGetLbGlobal(var2->getVar())));

   if(comp(scip, SCIPvarGetObj(GCGpricingVarGetOriginalVar(var1->getVar())), SCIPvarGetObj(GCGpricingVarGetOriginalVar(var2->getVar()))) != 0)
      return comp(scip, SCIPvarGetObj(GCGpricingVarGetOriginalVar(var1->getVar())), SCIPvarGetObj(GCGpricingVarGetOriginalVar(var2->getVar())));
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
   AUT_CONS* aut1 = (AUT_CONS*)elem1;
   AUT_CONS* aut2 = (AUT_CONS*)elem2;
   return comp( aut1->getScip(), aut1, aut2);
}

static
SCIP_DECL_SORTPTRCOMP(sortptrvar)
{
   AUT_VAR* aut1 = (AUT_VAR*)elem1;
   AUT_VAR* aut2 = (AUT_VAR*)elem2;
   return comp( aut1->getScip(), aut1, aut2);
}

static
SCIP_DECL_SORTPTRCOMP(sortptrval)
{
   AUT_COEF* aut1 = (AUT_COEF*)elem1;
   AUT_COEF* aut2 = (AUT_COEF*)elem2;
   return comp( aut1->getScip(), aut1->getVal(), aut2->getVal());
}

/** hook function to save the permutation of the graph */
static
void hook(
   void* user_param,                 /**< data structure to save hashmaps with permutation */
   unsigned int N,                          /**< number of permutations */
   const unsigned int *aut                  /**< array of permutations */
   )
{

   int i;
   int j;
   int n;
   int nvars;
   int nconss;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_CONS** conss1;
   SCIP_CONS** conss2;
   AUT_HOOK* hook = (AUT_HOOK*) user_param;

   j = 0;
   n = hook->getNNodes();

   if(hook->getBool())
      return;

   SCIPdebugMessage("Looking for a permutation from [0,%d] bijective to [%d:%d]\n", n/2-1, n/2, n-1);
   for( i = 0; i < n / 2; i++ )
   {
      assert(aut[i] < INT_MAX);
      if( (int) (aut[i]) >= n / 2 )
      {
         SCIPdebugMessage("%d -> %d\n", i, aut[i]);
         j++;
      }
      else
         break;
   }
   if( j == n / 2 )
   {
      hook->setBool(TRUE);
   }

   SCIPdebugMessage("Permutation %s found.\n", hook->getBool() ? "":"not");
   SCIPdebugMessage("j = %d\n", j);

   if( !hook->getBool() )
      return;

   nvars = SCIPgetNVars(hook->getScips()[0]);
   assert(nvars == SCIPgetNVars(hook->getScips()[1]));

   vars1 = SCIPgetVars(hook->getScips()[0]);
   vars2 = SCIPgetVars(hook->getScips()[1]);
   nconss = SCIPgetNConss(hook->getScips()[0]);
   assert(nconss == SCIPgetNConss(hook->getScips()[1]));

   conss1 = SCIPgetConss(hook->getScips()[0]);
   conss2 = SCIPgetConss(hook->getScips()[1]);

   for( i = 0; i < nvars+nconss; i++ )
   {
      /* Assuming the following layout:
       *  0 ... nconss-1 = vertex ids for constraints
       *  nconss ... nconss+nvars-1 = vertex ids for variables
       *  nconss+nvars ... n-1 = nonzero entries (not relevant)
       */
      if( i < nconss )
      {
         int consindex = i;
         int consindex2 = aut[i]-n/2;
         assert( consindex2 >= 0);
         assert( consindex2 < nconss);
         SCIP_CONS* cons1 = conss1[consindex];
         SCIP_CONS* cons2 = conss2[consindex2];
         SCIP_CALL_ABORT( SCIPhashmapInsert(hook->getConsHash(), cons1, cons2) );
         SCIPdebugMessage("cons <%s> <-> cons <%s>\n", SCIPconsGetName(cons1), SCIPconsGetName(cons2));
      }
      else if( i < nvars+nconss )
      {
         int varindex = i-nconss;
         int varindex2 = aut[i]-nconss-n/2;
         assert( varindex2 >= 0);
         assert( varindex2 < nvars);
         SCIP_VAR* var1 = vars1[varindex];
         SCIP_VAR* var2 = vars2[varindex2];
         SCIP_CALL_ABORT( SCIPhashmapInsert(hook->getVarHash(), var1, var2) );
         SCIPdebugMessage("var <%s> <-> var <%s>\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
      }
   }
}

/** tests if two scips have the same number of variables */
static
SCIP_RETCODE testScipVars(
   SCIP* scip1,                             /**< first SCIP data structure */
   SCIP* scip2,                             /**< second SCIP data structure */
   SCIP_RESULT* result                      /**< result pointer to indicate success or failure */
   )
{
   if(SCIPgetNVars(scip1) != SCIPgetNVars(scip2))
   {
      *result = SCIP_DIDNOTFIND;
   }
   return SCIP_OKAY;
}

/** tests if two scips have the same number of constraints */
static
SCIP_RETCODE testScipCons(
   SCIP* scip1,                             /**< first SCIP data structure */
   SCIP* scip2,                             /**< second SCIP data structure */
   SCIP_RESULT* result                      /**< result pointer to indicate success or failure */
   )
{
   if(SCIPgetNConss(scip1) != SCIPgetNConss(scip2))
   {
      *result = SCIP_DIDNOTFIND;
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE allocMemory(
   SCIP* scip,                              /**< SCIP data structure */
   AUT_COLOR* colorinfo,                    /**< struct to save intermediate information */
   int nconss,                              /**< number of constraints */
   int nvars                                /**< number of variables */
   )
{
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarraycoefs, (nconss * nvars)));
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarrayvars, nvars));
   SCIP_CALL( SCIPallocMemoryArray(scip, &colorinfo->ptrarrayconss, nconss));
   return SCIP_OKAY;
}

static SCIP_RETCODE reallocMemory(
   SCIP* scip,                              /**< SCIP data structure */
   AUT_COLOR* colorinfo,                    /**< struct to save intermediate information */
   int nconss,                              /**< number of constraints */
   int nvars                                /**< number of variables */
   )
{
   SCIP_CALL( SCIPreallocMemoryArray(scip, &colorinfo->ptrarraycoefs, colorinfo->lencoefsarray + (nconss * nvars)));
   SCIP_CALL( SCIPreallocMemoryArray(scip, &colorinfo->ptrarrayvars, colorinfo->lenvarsarray + nvars));
   SCIP_CALL( SCIPreallocMemoryArray(scip, &colorinfo->ptrarrayconss, colorinfo->lenconssarray + nconss));
   return SCIP_OKAY;
}

static
SCIP_RETCODE freeMemory(
   SCIP* scip,                              /**< SCIP data structure */
   AUT_COLOR* colorinfo                     /**< struct to save intermediate information */
   )
{
   int i;
   AUT_COEF* scoef;
   AUT_VAR* svar;
   AUT_CONS* scons;

   for(i = 0; i < colorinfo->lenvarsarray; i++ ){
      svar =  (AUT_VAR*) colorinfo->ptrarrayvars[i];
      delete svar;
   }
   for(i = 0; i < colorinfo->lenconssarray; i++ ){
      scons = (AUT_CONS*) colorinfo->ptrarrayconss[i];
      delete scons;
   }
   for(i = 0; i < colorinfo->lencoefsarray; i++ ){
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
   SCIP*                origscip,           /**< SCIP data structure */
   SCIP**               scips,              /**< SCIPs to compare */
   int                  nscips,             /**< number of SCIPs */
   AUT_COLOR*           colorinfo,          /**< data structure to save intermediate data */
   SCIP_RESULT*         result              /**< result pointer to indicate success or failure */
   )
{
   int i;
   int j;
   int s;
   int ncurvars;
   int nconss;
   int nvars;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   AUT_COEF* scoef;
   AUT_VAR* svar;
   AUT_CONS* scons;
   SCIP_Bool added;

   added = FALSE;

   //allocate max n of coefarray, varsarray, and boundsarray in origscip
   nconss = SCIPgetNConss(scips[0]);
   nvars = SCIPgetNVars(scips[0]);
   allocMemory(origscip, colorinfo, nconss, nvars);

   for( s = 0; s < nscips && *result == SCIP_SUCCESS; ++s )
   {
      SCIP *scip = scips[s];
      conss = SCIPgetConss(scip);
      vars = SCIPgetVars(scip);
      SCIPdebugMessage("Handling SCIP %i (%d x %d)\n", s, nconss, nvars);
      //save the properties of variables in a struct array and in a sorted pointer array
      for( i = 0; i < nvars; i++ )
      {
         svar = new AUT_VAR(scip, vars[i]);
         //add to pointer array iff it doesn't exist
         colorinfo->insert( svar, &added );
         if( s > 0 && added)
         {
           *result = SCIP_DIDNOTFIND;
            break;
         }
         //otherwise free allocated memory
         if( !added )
            delete svar;
      }
      //save the properties of constraints in a struct array and in a sorted pointer array
      for( i = 0; i < nconss && *result == SCIP_SUCCESS; i++ )
      {
         SCIP_Real* curvals;
         ncurvars = SCIPgetNVarsXXX(scip, conss[i]);    //SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(origscip), SCIPgetNVars(scip1)+1) );
         if( ncurvars == 0 )
            continue;
         scons = new AUT_CONS(scip, conss[i]);
         //add to pointer array iff it doesn't exist
         colorinfo->insert( scons, &added );
         if( s > 0 && added)
         {
           *result = SCIP_DIDNOTFIND;
           break;
         }
         //otherwise free allocated memory
         if( !added )
            delete scons;

         SCIP_CALL( SCIPallocMemoryArray(origscip, &curvals, ncurvars));
         SCIPgetValsXXX(scip, conss[i], curvals, ncurvars);
         //save the properties of variables of the constraints in a struct array and in a sorted pointer array
         for( j = 0; j < ncurvars; j++ )
         {
            scoef = new AUT_COEF(scip, curvals[j] );
            //test, whether the coefficient is not zero
            if( !SCIPisEQ(scip, scoef->getVal(), 0) )
            {
               //add to pointer array iff it doesn't exist
               colorinfo->insert( scoef, &added );
               if( s > 0 && added)
               {
                  *result = SCIP_DIDNOTFIND;
                  break;
               }
            }
            //otherwise free allocated memory
            if( !added )
               delete scoef;
         }
         SCIPfreeMemoryArray(origscip, &curvals);
      }
      //size of the next instance, in order to allocate memory
      if( s < nscips - 1 )
      {
         nconss = SCIPgetNConss(scips[s + 1]);
         nvars = SCIPgetNVars(scips[s + 1]);
      }
      //set zero, if no next instance exists
      else
      {
         nconss = 0;
         nvars = 0;
      }
      //reallocate memory with size of ptrarray[bounds, vars, coefs] + max of scip[i+1]
      SCIP_CALL( reallocMemory(origscip, colorinfo, nconss, nvars) );
   }

   /* add color information for master constraints */
   SCIP_CONS** origmasterconss = GCGrelaxGetLinearOrigMasterConss(origscip);
   int nmasterconss = GCGrelaxGetNMasterConss(origscip);

   SCIP_CALL( reallocMemory(origscip, colorinfo, nmasterconss, SCIPgetNVars(origscip)) );

   for( i = 0; i < nmasterconss && *result == SCIP_SUCCESS; ++i )
   {
      SCIP_CONS* mastercons = origmasterconss[i];
      SCIP_Real* curvals = SCIPgetValsLinear(origscip, mastercons);
      ncurvars = SCIPgetNVarsLinear(origscip, mastercons);

      /* add right color for master constraint */
      scons = new AUT_CONS(origscip, mastercons);
      colorinfo->insert( scons, &added );

      /* if it hasn't been added, it is already present */
      if(!added)
         delete scons;

      for( j = 0; j < ncurvars; ++j )
      {
         scoef = new AUT_COEF(origscip, curvals[j] );

         added = FALSE;

         if( !SCIPisEQ(origscip, scoef->getVal(), 0) )
         {
            colorinfo->insert( scoef, &added );
         }

         if( !added )
            delete scoef;
      }
   }

   return SCIP_OKAY;
}
/** create a graph out of an array of scips */
static
SCIP_RETCODE createGraph(
   SCIP*                origscip,           /**< SCIP data structure */
   SCIP**               scips,              /**< SCIPs to compare */
   int                  nscips,             /**< number of SCIPs */
   int*                 pricingindices,     /**< indices of the given pricing problems */
   AUT_COLOR            colorinfo,          /**< result pointer to indicate success or failure */
   bliss::Graph*        graph,              /**< graph needed for discovering isomorphism */
   int*                 pricingnodes,             /**< number of pricing nodes without master  */
   SCIP_RESULT*         result              /**< result pointer to indicate success or failure */
   )
{
   int i;
   int j;
   int z;
   int s;
   int nvars;
   int nconss;
   int ncurvars;
   int curvar;
   int* nnodesoffset;
   int color;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_VAR** curvars;
   SCIP_Real* curvals;

   unsigned int nnodes;
   nnodes = 0;
   //building the graph out of the arrays
   bliss::Graph* h = graph;

   int* pricingnonzeros;
   int* mastercoefindex;
   SCIP_CALL( SCIPallocMemoryArray(origscip, &pricingnonzeros, nscips) );
   SCIP_CALL( SCIPallocMemoryArray(origscip, &nnodesoffset, nscips) );
   SCIP_CALL( SCIPallocMemoryArray(origscip, &mastercoefindex, nscips) );
   BMSclearMemoryArray(pricingnonzeros, nscips);
   BMSclearMemoryArray(nnodesoffset, nscips);
   BMSclearMemoryArray(mastercoefindex, nscips);

   SCIP_CONS** origmasterconss = GCGrelaxGetLinearOrigMasterConss(origscip);
   int nmasterconss = GCGrelaxGetNMasterConss(origscip);

   for( s = 0; s < nscips && *result == SCIP_SUCCESS; ++s)
   {
      SCIPdebugMessage("Pricing problem %d\n", pricingindices[s]);
      SCIP *scip = scips[s];
      nnodesoffset[s] = nnodes;
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

         color = colorinfo.get( AUT_CONS(scip, conss[i]) );

         if(color == -1) {
            *result = SCIP_DIDNOTFIND;
            break;
         }

         SCIPdebugMessage("cons <%s> color %d\n", SCIPconsGetName(conss[i]), color);
         h->add_vertex(color);
         nnodes++;
      }
      //add a node for every variable
      for( i = 0; i < nvars && *result == SCIP_SUCCESS; i++ )
      {
         color = colorinfo.get( AUT_VAR(scip, vars[i]) );

         if(color == -1) {
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
         int conscolor = colorinfo.get(AUT_CONS(scip, conss[i]));
         ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
         if( ncurvars == 0 )
            continue;
         SCIP_CALL( SCIPallocMemoryArray(origscip, &curvars, ncurvars));
         SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars);
         SCIP_CALL( SCIPallocMemoryArray(origscip, &curvals, ncurvars));
         SCIPgetValsXXX(scip, conss[i], curvals, ncurvars);
         for( j = 0; j < ncurvars; j++ )
         {
            int varcolor = colorinfo.get( AUT_VAR(scip, curvars[j] ) ) + colorinfo.getLenCons();
            color = colorinfo.get( AUT_COEF(scip, curvals[j] ) );
            if( color == -1 )
            {
               *result = SCIP_DIDNOTFIND;
               break;
            }
            color += colorinfo.getLenCons() + colorinfo.getLenVar();
            curvar = SCIPvarGetProbindex(curvars[j]);
            h->add_vertex(color);
            nnodes++;
            h->add_edge(nnodesoffset[s] + i, nnodesoffset[s] + nconss + nvars + z);
            h->add_edge(nnodesoffset[s] + nconss + nvars + z, nnodesoffset[s]+nconss + curvar);
            SCIPdebugMessage("nz: c <%s> (id: %d, color: %d) -> nz (id: %d) (value: %f, color: %d) -> var <%s> (id: %d, color: %d) \n",
                              SCIPconsGetName(conss[i]),
                              nnodesoffset[s] + i,
                              conscolor,
                              nnodesoffset[s] + nconss + nvars + z,
                              curvals[j],
                              color+colorinfo.getLenCons() + colorinfo.getLenVar(),
                              SCIPvarGetName(curvars[j]),
                              nnodesoffset[s] + nconss + curvar,
                              varcolor);
            z++;
         }
         SCIPfreeMemoryArray(origscip, &curvals);
         SCIPfreeMemoryArray(origscip, &curvars);
      }
      pricingnonzeros[s] = z;

      /* add coefficient nodes for nonzeros in the master */
      for( i = 0; i < nmasterconss && *result == SCIP_SUCCESS; ++i )
      {
         SCIP_CONS* mastercons = origmasterconss[i];
         curvars = SCIPgetVarsLinear(origscip, mastercons);
         curvals = SCIPgetValsLinear(origscip, mastercons);
         ncurvars = SCIPgetNVarsLinear(origscip, mastercons);
         for( j = 0; j < ncurvars; ++j )
         {
            if( GCGvarIsLinking(curvars[j]) )
            {
               SCIPdebugMessage("Var <%s> is linking, abort detection.\n", SCIPvarGetName(curvars[j]));
               *result = SCIP_DIDNOTFIND;
               return SCIP_OKAY;
            }
            int block = GCGvarGetBlock(curvars[j]);

            /* ignore if the variable belongs to a different block */
            if( block != pricingindices[s] )
            {
               SCIPdebugMessage("Var <%s> belongs to a different block (%d)\n", SCIPvarGetName(curvars[j]), block);
               continue;
            }


            color = colorinfo.get(AUT_COEF(origscip, curvals[j]));
            assert(color != -1);
            color += colorinfo.getLenCons() + colorinfo.getLenVar();

            /* add coefficent node for current coeff */
            h->add_vertex(color);
            SCIPdebugMessage("master nz for var <%s> (id: %d) (value: %f, color: %d)\n", SCIPvarGetName(curvars[j]), nnodes, curvals[i], color);
            nnodes++;
         }
      }
      SCIPdebugMessage("Iteration %d: nnodes = %d\n", s, nnodes);
      assert(*result == SCIP_SUCCESS && nnodes == h->get_nof_vertices());
   }
   /* connect the created graphs with nodes for the master problem */

   SCIPdebugMessage( "handling %d masterconss\n", nmasterconss);
   *pricingnodes = nnodes;

   for( i = 0; i < nmasterconss && *result == SCIP_SUCCESS; ++i )
   {
      SCIP_CONS* mastercons = origmasterconss[i];
      curvars = SCIPgetVarsLinear(origscip, mastercons);
      curvals = SCIPgetValsLinear(origscip, mastercons);
      ncurvars = SCIPgetNVarsLinear(origscip, mastercons);

      SCIPdebugMessage("Handling cons <%s>\n", SCIPconsGetName(mastercons));

      /* create node for masterconss and get right color */
      int conscolor = colorinfo.get(AUT_CONS(origscip, mastercons));
      assert(conscolor != -1);
      h->add_vertex(conscolor);
      int masterconsnode = nnodes;
      nnodes++;

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP* pricingscip = NULL;
         if( GCGvarIsLinking(curvars[j]) )
         {
            SCIPdebugMessage("Var <%s> is linking, abort detection.\n", SCIPvarGetName(curvars[j]));
            *result = SCIP_DIDNOTFIND;
            return SCIP_OKAY;
         }
         int block = GCGvarGetBlock(curvars[j]);
         int index = -1;
         SCIPdebugMessage("Var <%s> is in block %d\n", SCIPvarGetName(curvars[j]), block);
         for( s = 0; s < nscips; ++s )
         {
            if( block == pricingindices[s])
            {
               index = s;
               pricingscip = scips[s];
               break;
            }
         }

         /* ignore if the variable belongs to a different block */
         if( pricingscip == NULL )
         {
            SCIPdebugMessage("Var <%s> belongs to a different block (%d)\n", SCIPvarGetName(curvars[j]), block);
            continue;
         }

         color = colorinfo.get(AUT_COEF(origscip, curvals[j]));
         assert(color != -1);
         color += colorinfo.getLenCons() + colorinfo.getLenVar();
         SCIP_VAR* pricingvar = GCGoriginalVarGetPricingVar(curvars[j]);

         /* get coefficient node for current coefficient */
         int coefnodeindex = nnodesoffset[index] + SCIPgetNVars(pricingscip) + SCIPgetNConss(pricingscip) + pricingnonzeros[index] + mastercoefindex[index];
         ++(mastercoefindex[index]);

         int varcolor = colorinfo.get(AUT_VAR(pricingscip, pricingvar));
         assert(varcolor != -1);
         varcolor += colorinfo.getLenCons();

         assert( (uint) masterconsnode < h->get_nof_vertices());
         assert( (uint) coefnodeindex < h->get_nof_vertices());
         /* master constraint and coefficient */
         h->add_edge(masterconsnode, coefnodeindex);
         SCIPdebugMessage("ma: c <%s> (id: %d, color: %d) -> nz (id: %d) (value: <%.6f> , color: %d) -> pricingvar <%s> (id: %d, color: %d)\n",
            SCIPconsGetName(mastercons),
            masterconsnode, conscolor, coefnodeindex, curvals[j], color, SCIPvarGetName(pricingvar),
            nnodesoffset[index] + SCIPgetNConss(pricingscip) + SCIPvarGetProbindex(pricingvar), varcolor);

         /* get node index for pricing variable and connect masterconss, coeff and pricingvar nodes */
         h->add_edge(coefnodeindex, nnodesoffset[index] + SCIPgetNConss(pricingscip) + SCIPvarGetProbindex(pricingvar));
      }
   }


   //free all allocated memory
   SCIP_CALL( freeMemory(origscip, &colorinfo) );
   SCIPfreeMemoryArray(origscip, &mastercoefindex);
   SCIPfreeMemoryArray(origscip, &nnodesoffset);
   SCIPfreeMemoryArray(origscip, &pricingnonzeros);

   return SCIP_OKAY;
}

/** compare two graphs w.r.t. automorphism */
extern
SCIP_RETCODE cmpGraphPair(
   SCIP*                 origscip,            /**< SCIP data structure */
   SCIP*                 scip1,               /**< first SCIP data structure to compare */
   SCIP*                 scip2,               /**< second SCIP data structure to compare */
   int                   prob1,               /**< index of first pricing prob */
   int                   prob2,               /**< index of second pricing prob */
   SCIP_RESULT*          result,              /**< result pointer to indicate success or failure */
   SCIP_HASHMAP*         varmap,              /**< hashmap to save permutation of variables */
   SCIP_HASHMAP*         consmap              /**< hashmap to save permutation of constraints */
   )
{
   bliss::Graph graph;
   bliss::Stats bstats;
   AUT_HOOK *ptrhook;
   AUT_COLOR colorinfo;
   int nscips;
   SCIP* scips[2];
   int pricingindices[2];
   int pricingnodes = 0;
   scips[0] = scip1;
   scips[1] = scip2;
   pricingindices[0] = prob1;
   pricingindices[1] = prob2;
   nscips = 2;
   *result = SCIP_SUCCESS;

   SCIP_CALL( testScipVars(scips[0], scips[1], result) );
   SCIP_CALL( testScipCons(scips[0], scips[1], result) );

   SCIP_CALL( setuparrays(origscip, scips, nscips, &colorinfo, result) );
   SCIP_CALL( createGraph(origscip, scips, nscips, pricingindices, colorinfo, &graph,  &pricingnodes, result) );

   ptrhook = new AUT_HOOK(varmap, consmap, FALSE, pricingnodes, scips);
   graph.find_automorphisms(bstats, hook, ptrhook);

   varmap = ptrhook->getVarHash();
   consmap = ptrhook->getConsHash();

   if( !ptrhook->getBool() )
      *result = SCIP_DIDNOTFIND;

   delete ptrhook;
   return SCIP_OKAY;
}

/** returns bliss version */
extern
const char* GCGgetBlissVersion(void)
{
   return bliss::version;
}
