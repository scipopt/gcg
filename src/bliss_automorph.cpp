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

/**@file    bliss_automorph.cpp
 * @brief   automorphism recognition of SCIPs
 * @author  Daniel Peters
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "graph.hh"
#include "bliss_automorph.h"
#include "scip_misc.h"
#include "scip/scip.h"
#include <cstring>


typedef struct struct_cons AUT_CONS;
typedef struct struct_var AUT_VAR;
typedef struct struct_coef AUT_COEF;
typedef struct struct_hook AUT_HOOK;
typedef struct struct_colorinformation AUT_COLOR;

/** saves a constraint with its corresponding scip */
struct struct_cons
{
   SCIP* scip;                              /**< SCIP data structure */
   SCIP_CONS* cons;                         /**< pointer to SCIP constraint */

   /** constructor for the constraint struct */
   struct_cons( SCIP* scip, SCIP_CONS* scons );

   /** getter for the SCIP constraint */
   SCIP_CONS* getCons();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves a variable with its corresponding scip */
struct struct_var
{
   SCIP* scip;                              /**< SCIP data structure */
   SCIP_VAR* var;                           /**< pointer to SCIP variable */

   /** constructor for the variable struct */
   struct_var( SCIP* scip, SCIP_VAR* svar );

   /** getter for the SCIP variable */
   SCIP_VAR* getVar();

   /** getter for the SCIP itself */
   SCIP* getScip();
};

/** saves a coefficient with its corresponding scip */
struct struct_coef
{
   SCIP* scip;                              /**< SCIP data structure */
   SCIP_Real val;                           /**< SCIP Real value */

   /** constructor for the coefficient struct */
   struct_coef( SCIP* scip, SCIP_Real val );

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
struct_var::struct_var(
   SCIP* scip_,                             /**< SCIP data structure */
   SCIP_VAR* svar                           /**< SCIP variable */
   )
{
   scip = scip_;
   var = svar;
}

/** constructor of the constraint struct */
struct_cons::struct_cons(
   SCIP* scip_,                             /**< SCIP data structure */
   SCIP_CONS* scons                         /**< SCIP constraint */
   )
{
   scip = scip_;
   cons = scons;
}

/** constructor of the coefficient struct */
struct_coef::struct_coef(
   SCIP* scip_,                             /**< SCIP data structure */
   SCIP_Real val_                           /**< SCIP value */
   )
{
   scip = scip_;
   val = val_;
}

/** constructor of the color struct */
struct_colorinformation::struct_colorinformation(
   int color_,                              /**< color of the nodes of the graph */
   int lenvars,                             /**< length of ptrvarsarray */
   int lenconss,                            /**< length of ptrconsarray */
   int lencoefs                             /**< length of ptrcoefsarray */
   )
{
   color = color_;
   lenvarsarray = lenvars;
   lenconssarray = lenconss;
   lencoefsarray = lencoefs;
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
   Scip* scip,                              /**< SCIP data structure */
   AUT_CONS* cons1,                         /**< constraint 1 to compare */
   AUT_CONS* cons2                          /**< constraint 2 to compare */
)
{
   if(comp(scip, SCIPgetRhsXXX(scip, cons1->getCons()), SCIPgetRhsXXX(scip, cons2->getCons())) != 0)
      return comp(scip, SCIPgetRhsXXX(scip, cons1->getCons()), SCIPgetRhsXXX(scip, cons2->getCons()));
   assert(SCIPisEQ(scip, SCIPgetRhsXXX(scip, cons1->getCons()), SCIPgetRhsXXX(scip, cons2->getCons())));

   if(comp(scip, SCIPgetLhsXXX(scip, cons1->getCons()), SCIPgetLhsXXX(scip, cons2->getCons())) != 0)
      return comp(scip, SCIPgetLhsXXX(scip, cons1->getCons()), SCIPgetLhsXXX(scip, cons2->getCons()));
   assert(SCIPisEQ(scip, SCIPgetLhsXXX(scip, cons1->getCons()), SCIPgetLhsXXX(scip, cons2->getCons())));

   return strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons1->getCons())), SCIPconshdlrGetName(SCIPconsGetHdlr(cons2->getCons())));
}

/** compare two variables of two scips */
static
int comp(
   Scip* scip,                              /**< SCIP data structure */
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

   if(comp(scip, SCIPvarGetObj(var1->getVar()), SCIPvarGetObj(var2->getVar())) != 0)
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
void hook(void* user_param,                 /**< data structure to save hashmaps with permutation */
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

   for( i = 0; i < n / 2; i++ )
   {
      assert(aut[i] < INT_MAX);
      if( (int) (aut[i]) >= n / 2 )
         j++;
      else
         break;
   }
   if( j == n / 2 )
      hook->setBool(TRUE);

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

   for(i = 0; i < nvars+nconss; i++)
   {
      if( i < nvars )
      {
         assert( aut[i] < INT_MAX);
         assert( (int) (aut[i]-n/2) < nvars);
         SCIP_CALL_ABORT( SCIPhashmapInsert(hook->getVarHash(), vars1[i], vars2[aut[i]-n/2]) );
      }
      else if (i < nvars+nconss)
      {
         assert( i-nvars >= 0);
         assert( aut[i] < INT_MAX);
         assert( (int) (aut[i]-nvars-n/2) < nconss);
         SCIP_CALL_ABORT( SCIPhashmapInsert(hook->getConsHash(), conss1[i-nvars], conss2[aut[i]-nvars-n/2]) );
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
      SCIPdebugMessage("Handling SCIP %i\n", s);
      conss = SCIPgetConss(scip);
      vars = SCIPgetVars(scip);
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
      reallocMemory(origscip, colorinfo, nconss, nvars);
   }
   return SCIP_OKAY;
}
/** create a graph out of an array of scips */
static
SCIP_RETCODE createGraph(
   SCIP*                origscip,           /**< SCIP data structure */
   SCIP**               scips,              /**< SCIPs to compare */
   int                  nscips,             /**< number of SCIPs */
   AUT_COLOR            colorinfo,          /**< result pointer to indicate success or failure */
   bliss::Graph*        graph,              /**< graph needed for discovering isomorphism */
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
   int nnodesoffset[nscips];
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
   for( s = 0; s < nscips && *result == SCIP_SUCCESS; ++s)
   {
      SCIP *scip = scips[s];
      bliss::Graph* h = graph;
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

         scons = new AUT_CONS(scip, conss[i]);
         color = colorinfo.get( *scons );
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
         svar = new AUT_VAR(scip, vars[i]);
         color = colorinfo.get( *svar );
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
         scons = new AUT_CONS(scip, conss[i]);
         ncurvars = SCIPgetNVarsXXX(scip, conss[i]);
         if( ncurvars == 0 )
            continue;
         SCIP_CALL( SCIPallocMemoryArray(origscip, &curvars, ncurvars));
         SCIPgetVarsXXX(scip, conss[i], curvars, ncurvars);
         SCIP_CALL( SCIPallocMemoryArray(origscip, &curvals, ncurvars));
         SCIPgetValsXXX(scip, conss[i], curvals, ncurvars);
         for( j = 0; j < ncurvars; j++ )
         {
            scoef = new AUT_COEF(scip, curvals[j] );
            color = colorinfo.get( *scoef );
            if( color == -1 )
            {
               *result = SCIP_DIDNOTFIND;
               break;
            }
            curvar = SCIPvarGetProbindex(curvars[j]);
            h->add_vertex(colorinfo.getLenCons() + colorinfo.getLenVar() + color);
            nnodes++;
            h->add_edge(nnodesoffset[s] + i, nnodesoffset[s] + nconss + nvars + z);
            h->add_edge(nnodesoffset[s] + nconss + nvars + z, nnodesoffset[s]+nconss + curvar);
            SCIPdebugMessage("nz: c <%s> (id: %d, colour: %d) -> nz (id: %d) (value: %f, colour: %d) -> var <%s> (id: %d, colour: %d) \n",
                              SCIPconsGetName(conss[i]),
                              nnodesoffset[s] + i,
                              colorinfo.get( *scons ),
                              nnodesoffset[s] + nconss + nvars + z,
                              scoef->getVal(),
                              color+colorinfo.getLenCons() + colorinfo.getLenVar(),
                              SCIPvarGetName(curvars[j]),
                              nnodesoffset[s] + nconss + curvar,
                              colorinfo.get( *svar)+colorinfo.getLenCons());
            z++;
         }
         SCIPfreeMemoryArray(origscip, &curvals);
         SCIPfreeMemoryArray(origscip, &curvars);
      }
      SCIPdebugMessage("Iteration %d: nnodes = %d\n", s, nnodes);
      assert(*result == SCIP_SUCCESS && nnodes == h->get_nof_vertices());
   }
   //free all allocated memory
   freeMemory(origscip, &colorinfo);
   return SCIP_OKAY;
}

/** compare two graphs w.r.t. automorphism */
extern
SCIP_RETCODE cmpGraphPair(
   SCIP*                 origscip,            /**< SCIP data structure */
   SCIP*                 scip1,               /**< first SCIP data structure to compare */
   SCIP*                 scip2,               /**< second SCIP data structure to compare */
   SCIP_RESULT*          result,              /**< result pointer to indicate success or failure */
   SCIP_HASHMAP* 	 varmap,              /**< hashmap to save permutation of variables */
   SCIP_HASHMAP* 	 consmap              /**< hashmap to save permutation of constraints */
   )
{
   bliss::Graph graph;
   bliss::Stats bstats;
   AUT_HOOK *ptrhook;
   AUT_COLOR *colorinfo;
   int nscips;
   SCIP* scips[2];

   colorinfo = new AUT_COLOR(0, 0, 0, 0 );
   scips[0] = scip1;
   scips[1] = scip2;
   nscips = 2;
   *result = SCIP_SUCCESS;

   SCIP_CALL( testScipVars(scips[0], scips[1], result) );
   SCIP_CALL( testScipCons(scips[0], scips[1], result) );

   SCIP_CALL( setuparrays(origscip, scips, nscips, colorinfo, result) );
   SCIP_CALL( createGraph(origscip, scips, nscips, *colorinfo, &graph, result) );

   ptrhook = new AUT_HOOK(varmap, consmap, FALSE, graph.get_nof_vertices(), scips);
   graph.find_automorphisms(bstats, hook, ptrhook);

   varmap = ptrhook->getVarHash();
   consmap = ptrhook->getConsHash();

   return SCIP_OKAY;
}
