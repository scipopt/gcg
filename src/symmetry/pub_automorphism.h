/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file    pub_automorphism.h
 * @ingroup PUBLICCOREAPI
 * @brief   helper functions for automorphism detection
 *
 * @author  Martin Bergner
 * @author  Daniel Peters
 * @author  Jonas Witt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PUB_AUTOMORPHISM_H_
#define PUB_AUTOMORPHISM_H_

#include "gcg/def.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct struct_graph AUT_GRAPH;
typedef struct struct_graph_data AUT_GRAPH_DATA;
typedef struct struct_cons AUT_CONS;
typedef struct struct_var AUT_VAR;
typedef struct struct_coef AUT_COEF;
typedef struct struct_colorinformation AUT_COLOR;
/**
* @ingroup AUTOMORPHISM
* @{
  */

#ifdef WITH_BLISS
/** returns bliss version */
GCG_EXPORT
void GCGgetBlissName(char* buffer, int len);
#endif

#ifdef WITH_NAUTY
/** returns nauty version */
GCG_EXPORT
void GCGgetNautyName(char* buffer, int len);
#endif

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
struct struct_graph
{
   /** initializes the graph */
   SCIP_RETCODE init(
      SCIP* scip,                                                  /**< SCIP struct */
      int nvertices                                                /**< number of vertices */
      );

   /** destroys the graph */
   SCIP_RETCODE destroy();

   /** sets the color of a vertex */
   void setColor(
      int vertex,                                                  /**< vertex */
      int color                                                    /**< color */
      );

   /** adds an edge */
   void addEdge(
      int v1,                                                      /**< first vertex */
      int v2                                                       /**< second vertex */
      );

   /** returns the number of vertices */
   unsigned int getNVertices();

   /** found automorphisms */
   SCIP_RETCODE findAutomorphisms(
      void* userdata,                                              /**< pointer to user data that will be passed to fhook */
      void (*fhook)(void*, unsigned int, const unsigned int*),     /**< pointer to function that will be called for each found generator */
      unsigned int searchnodelimit,                                /**< search node limit (requires patched bliss version) */
      unsigned int generatorlimit                                  /**< generator limit (requires patched bliss version or version >=0.76) */
      );

   /** signals that the search should be terminated */
   void terminateSearch();

   AUT_GRAPH_DATA* graphdata;
};

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
   int                  alloccoefsarray;    /**< allocated size of ptrarraycoefs */

   void**               ptrarraycoefs;      /**< array of pointers to coefficient */
   void**               ptrarrayvars;       /**< array of pointers to variables */
   void**               ptrarrayconss;      /**< array of pointers to constraints */

   SCIP_Bool            onlysign;           /**< use sign of values instead of values? (should be FALSE if we check whether pricin problems can be aggregated) */

   /** default constructor */
   struct_colorinformation();

   /** insert a variable to its pointer array */
   SCIP_RETCODE insert( AUT_VAR* svar, SCIP_Bool* added);

   /** insert a constraint to its pointer array */
   SCIP_RETCODE insert( AUT_CONS* scons, SCIP_Bool* added);

   /** insert a coefficient to its pointer array */
   SCIP_RETCODE insert( AUT_COEF* scoef, SCIP_Bool* added);

   /** getter for the length of the variable array */
   int getLenVar();

   /** getter for the length of the constraint array */
   int getLenCons();

   /** getter for the variable struct */
   int get( AUT_VAR svar);

   /** getter for the constraint struct */
   int get( AUT_CONS scons);

   /** getter for the coefficient struct */
   int get( AUT_COEF scoef);

   /** set onlysign bool */
   SCIP_RETCODE setOnlySign(SCIP_Bool onlysign_);

   /** get onlysign bool */
   SCIP_Bool getOnlySign();
};


#endif
/** @} */
#endif /* PUB_AUTOMORPHISM_H_ */
