/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   struct_event.h
 * @brief  datastructures for managing events
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_VARDATA_H__
#define __SCIP_STRUCT_VARDATA_H__


/** type of the variable */
enum GCG_Vartype
{
   GCG_VARTYPE_ORIGINAL   = 0,       /**< variable belongs to original problem */
   GCG_VARTYPE_PRICING    = 1,       /**< variable belongs to a pricing problem */
   GCG_VARTYPE_MASTER     = 2        /**< variable belongs to the master problem */
};
typedef enum GCG_Vartype GCG_VARTYPE;

/** data for original variables */
struct GCG_OrigVarData
{
   SCIP_VAR*             pricingvar;             /**< corresponding variable in the pricing program */
   SCIP_Real*            coefs;                  /**< coefficiants in the linking constraints of the original program */
   int                   ncoefs;                 /**< number of coefficiants */
   SCIP_VAR**            mastervars;             /**< variables in the master problem that contain the variable */
   SCIP_Real*            mastervals;             /**< value of this variable in the master problem variables */
   int                   nmastervars;            /**< number of corresponding master variables */
   int                   maxmastervars;          /**< length of arrays mastervars and vals */
};
typedef struct GCG_OrigVarData GCG_ORIGVARDATA;

/** data for pricing variables */
struct GCG_PricingVarData
{
   SCIP_VAR**            origvars;                /**< corresponding variables in the original program */
   int                   norigvars;               /**< number of corresponding variables in the original program */
};
typedef struct GCG_PricingVarData GCG_PRICINGVARDATA;


/** @todo don't copy the pointers to the original vars for each master variable, store them in a central place */
/** data for master variables */
struct GCG_MasterVarData
{
   int                    norigvars;             /**< number of variables in the original program corresponding to  the current variable */
   SCIP_VAR**             origvars;              /**< variables in the original program corresponding to  the current variable */
   SCIP_Real*             origvals;                  /**< this variable represents vals[i] times the variable origvars[i] in the
                                                  *   original program */
};
typedef struct GCG_MasterVarData GCG_MASTERVARDATA;

/** variable data structure */
struct SCIP_VarData
{
   union
   {
      GCG_ORIGVARDATA    origvardata;        /**< data for original variables */
      GCG_PRICINGVARDATA pricingvardata;     /**< data for pricing variables */
      GCG_MASTERVARDATA  mastervardata;      /**< data for variable of the master problem */
   } data;
   GCG_VARTYPE        vartype;               /**< type of variable */
   int                blocknr;               /**< number of the block and pricing problem, the variable belongs to */
};


#endif
