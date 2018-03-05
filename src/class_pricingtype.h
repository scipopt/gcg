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

/**@file   class_pricingtype.h
 * @brief  abstraction for SCIP pricing types
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_CLASS_PRICINGTYPE_H__
#define GCG_CLASS_PRICINGTYPE_H__

#include "objscip/objscip.h"
#include "pricer_gcg.h"

class PricingType
{
protected:
   SCIP_CLOCK* clock;

   int calls;
   int maxvarsround;
   int maxvarsroundroot;
   int maxsuccessfulmips;
   int maxrounds;
   double mipsrel;
   double mipsrelroot;
   GCG_PRICETYPE type;
   SCIP *scip_;

public:
   /** constructor */
   PricingType();

   PricingType(SCIP *p_scip);

   /** destructor */
   virtual ~PricingType();

   /** get dual value of a constraint */
   virtual double consGetDual(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONS*         cons                /**< constraint to get dual for */
   ) const =0;

   /** get dual value of a row */
   virtual double rowGetDual(
      SCIP_ROW*          row                 /**< row to get dual value for */
   ) const =0;

   /** get objective value of variable */
   virtual double varGetObj(
      SCIP_VAR*          var                 /**< variable to get objective value for */
   ) const =0;

   /** adds parameters to the SCIP data structure */
   virtual SCIP_RETCODE addParameters() =0;

   /** returns whether the optimal pricing can be aborted */
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull /**< number of non-Null pricing problems*/
   ) const = 0;

   /** returns whether the heuristic pricing can be aborted */
   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull /**< number of non-Null pricing problems*/
   ) const = 0 ;

   /** starts the clock */
   virtual SCIP_RETCODE startClock();

   /** stops the clock */
   virtual SCIP_RETCODE stopClock();

   /** returns the time of the clock */
   virtual double getClockTime() const;

   /** returns the maximal number of rounds */
   int getMaxrounds() const
   {
      return maxrounds;
   }

   /** returns the maximal number of successful mip solutions */
   int getMaxsuccessfulmips() const
   {
      return maxsuccessfulmips;
   }

   /** returns the maximal number of variables at the root node */
   int getMaxvarsroundroot() const
   {
      return maxvarsroundroot;
   }

   /** returns the relative ratio of MIPs to be solved */
   double getMipsrel() const
   {
      return mipsrel;
   }

   /** returns the relative ratio of MIPs to be solved at the root node */
   double getMipsrelroot() const
   {
      return mipsrelroot;
   }

   /** returns the type of this pricing type */
   GCG_PRICETYPE getType() const
   {
      return type;
   }

   /** returns the number of calls so far */
   int getCalls() const
   {
      return calls;
   }

   /** returns the maximal number of vars per pricing round */
   int getMaxvarsround() const
   {
      return maxvarsround;
   }

   /** increases the number of calls */
   virtual void incCalls()
   {
      calls++;
   }

   /** resets the number of calls and the clock for a restart */
   SCIP_RETCODE resetCalls() {
      calls = 0;
      SCIP_CALL( SCIPresetClock(scip_, clock) );
      return SCIP_OKAY;
   }

};

class ReducedCostPricing : public PricingType
{
public:
   /** constructor */
   ReducedCostPricing();

   ReducedCostPricing(
          SCIP* p_scip
          );

    /** destructor */
    virtual ~ReducedCostPricing() {}

    virtual SCIP_RETCODE addParameters();
    virtual SCIP_Real consGetDual(SCIP* scip, SCIP_CONS* cons) const;
    virtual SCIP_Real rowGetDual(SCIP_ROW* row) const;
    virtual SCIP_Real varGetObj(SCIP_VAR *var) const ;
    virtual SCIP_Bool canOptimalPricingBeAborted(
          int                  nfoundvars,         /**< number of variables found so far */
          int                  solvedmips,         /**< number of MIPS solved so far */
          int                  successfulmips,     /**< number of sucessful mips solved so far */
          SCIP_Real            successfulmipsrel,  /**< number of sucessful mips solved so far */
          int                  npricingprobsnotnull /**< number of non-Null pricing problems*/
      ) const ;

    virtual SCIP_Bool canHeuristicPricingBeAborted(
        int                  nfoundvars,         /**< number of variables found so far */
        int                  solvedmips,         /**< number of MIPS solved so far */
        int                  successfulmips,     /**< number of sucessful mips solved so far */
        SCIP_Real            successfulmipsrel,  /**< number of sucessful mips solved so far */
        int                  npricingprobsnotnull /**< number of non-Null pricing problems*/
    ) const;

};

class FarkasPricing : public PricingType
{
public:
   /** constructor */
   FarkasPricing();
   FarkasPricing(
          SCIP* p_scip
          );
   /** destructor */
   virtual ~FarkasPricing() {}
   virtual SCIP_RETCODE addParameters();
   virtual SCIP_Real consGetDual(SCIP *scip, SCIP_CONS *cons) const;
   virtual SCIP_Real rowGetDual(SCIP_ROW* row) const;
   virtual SCIP_Real varGetObj(SCIP_VAR *var) const;
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull /**< number of non-Null pricing problems*/
   ) const;
   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull /**< number of non-Null pricing problems*/
   ) const;
};

#endif /* CLASS_PRICINGTYPE_H_ */
