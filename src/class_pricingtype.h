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

/**@file   class_pricingtype.h
 * @brief  GCG variable pricer
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "objscip/objscip.h"
#include "pricer_gcg.h"
#include "class_instanciated.h"

#ifndef GCG_CLASS_PRICINGTYPE_H_
#define GCG_CLASS_PRICINGTYPE_H_


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
   PricingType(SCIP *p_scip);
   virtual ~PricingType();
   virtual double consGetDual(SCIP *scip, SCIP_CONS *cons) =0;
   virtual double rowGetDual(SCIP_ROW *row) =0;
   virtual double varGetObj(SCIP_VAR *var) =0;
   virtual SCIP_RETCODE addParameters() =0;
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull
   ) = 0;

   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull
   ) = 0;
    virtual SCIP_RETCODE startClock();
    virtual SCIP_RETCODE stopClock();
    virtual double getClockTime();

    int getMaxrounds() const
    {
        return maxrounds;
    }

    int getMaxsuccessfulmips() const
    {
        return maxsuccessfulmips;
    }

    int getMaxvarsroundroot() const
    {
        return maxvarsroundroot;
    }

    double getMipsrel() const
    {
        return mipsrel;
    }

    double getMipsrelroot() const
    {
        return mipsrelroot;
    }

    GCG_PRICETYPE getType() const
    {
        return type;
    }

    int getCalls() const
    {
        return calls;
    }

    int getMaxvarsround() const
    {
        return maxvarsround;
    }

    inline virtual void incCalls()
    {
        calls++;
    }

    SCIP_RETCODE resetCalls() {
       calls = 0;
       SCIP_CALL( SCIPresetClock(scip_, clock) );
       return SCIP_OKAY;
    };

};

class ReducedCostPricing : public PricingType, public Instanciated<ReducedCostPricing>
{
public:
   ReducedCostPricing(
          SCIP* p_scip
          );
    virtual ~ReducedCostPricing() {};
    virtual SCIP_RETCODE addParameters();
    virtual SCIP_Real consGetDual(SCIP* scip, SCIP_CONS* cons);
    virtual SCIP_Real rowGetDual(SCIP_ROW* row);
    virtual SCIP_Real varGetObj(SCIP_VAR *var);
    virtual SCIP_Bool canOptimalPricingBeAborted(
          int                  nfoundvars,         /**< number of variables found so far */
          int                  solvedmips,         /**< number of MIPS solved so far */
          int                  successfulmips,     /**< number of sucessful mips solved so far */
          SCIP_Real            successfulmipsrel,     /**< number of sucessful mips solved so far */
          int                  npricingprobsnotnull
      );

    virtual SCIP_Bool canHeuristicPricingBeAborted(
        int                  nfoundvars,         /**< number of variables found so far */
        int                  solvedmips,         /**< number of MIPS solved so far */
        int                  successfulmips,     /**< number of sucessful mips solved so far */
        SCIP_Real            successfulmipsrel,     /**< number of sucessful mips solved so far */
        int                  npricingprobsnotnull
    );

};

class FarkasPricing : public PricingType, public Instanciated<FarkasPricing>
{
public:
   FarkasPricing(
          SCIP* p_scip
          );
   virtual ~FarkasPricing() {};
   virtual SCIP_RETCODE addParameters();
   virtual SCIP_Real consGetDual(SCIP *scip, SCIP_CONS *cons);
   virtual SCIP_Real rowGetDual(SCIP_ROW* row);
   virtual SCIP_Real varGetObj(SCIP_VAR *var);
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull
   );
   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int               nfoundvars,         /**< number of variables found so far */
      int               solvedmips,         /**< number of MIPS solved so far */
      int               successfulmips,     /**< number of sucessful mips solved so far */
      SCIP_Real         successfulmipsrel,  /**< number of sucessful mips solved so far */
      int               npricingprobsnotnull
   );
};

#endif /* CLASS_PRICINGTYPE_H_ */
