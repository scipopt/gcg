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

/**@file   class_pricingtype.cpp
 * @brief  GCG variable pricer
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "objscip/objscip.h"
#include "pricer_gcg.h"
#ifndef GCG_CLASS_PRICINGTYPE_H_
#define GCG_CLASS_PRICINGTYPE_H_

class PricingType
{
protected:
   SCIP_CLOCK* clock;
   int calls;
   int maxvarsround;
    GCG_PRICETYPE type;
public:
    PricingType();
    virtual ~PricingType() {};
    virtual double consGetDual(SCIP *scip, SCIP_CONS *cons) = 0;
    virtual double rowGetDual(SCIP_ROW* row) = 0;
    virtual double varGetObj(SCIP_VAR *var) = 0;
    virtual SCIP_RETCODE startClock(SCIP *scip);
    virtual SCIP_RETCODE stopClock(SCIP *scip);
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

    void setMaxvarsround(int maxvarsround)
    {
        this->maxvarsround = maxvarsround;
    }

    inline virtual void incCalls()
    {
        calls++;
    }

};

class ReducedCostPricing : public PricingType
{
public:
    ReducedCostPricing();
    virtual ~ReducedCostPricing() {};
    virtual double consGetDual(SCIP* scip, SCIP_CONS* cons);
    virtual double rowGetDual(SCIP_ROW* row);
    virtual double varGetObj(SCIP_VAR *var);

};

class FarkasPricing : public PricingType
{
public:
    FarkasPricing();
    virtual ~FarkasPricing() {};
    virtual double consGetDual(SCIP *scip, SCIP_CONS *cons);
    virtual double rowGetDual(SCIP_ROW* row);
    virtual double varGetObj(SCIP_VAR *var);

};

#endif /* CLASS_PRICINGTYPE_H_ */
