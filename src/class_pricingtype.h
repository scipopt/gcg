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
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_CLASS_PRICINGTYPE_H__
#define GCG_CLASS_PRICINGTYPE_H__

#include "objscip/objscip.h"
#include "pricer_gcg.h"

class PricingType
{
protected:
  SCIP*                 scip_;
  GCG_PRICETYPE         type;
  SCIP_CLOCK*           clock;

  int                   calls;
  int                   maxrounds;
  int                   maxcolsroundroot;
  int                   maxcolsround;
  int                   maxcolsprobroot;
  int                   maxcolsprob;
  int                   maxsuccessfulprobs;
  SCIP_Real             relmaxprobsroot;
  SCIP_Real             relmaxprobs;

public:
   /** constructor */
   PricingType();

   PricingType(
      SCIP                  *p_scip
      );

   /** destructor */
   virtual ~PricingType();

   /** get dual value of a constraint */
   virtual SCIP_Real consGetDual(
      SCIP*                 scip,               /**< SCIP data structure */
      SCIP_CONS*            cons                /**< constraint to get dual for */
      ) const = 0;

   /** get dual value of a row */
   virtual SCIP_Real rowGetDual(
      SCIP_ROW*             row                 /**< row to get dual value for */
      ) const = 0;

   /** get objective value of variable */
   virtual SCIP_Real varGetObj(
      SCIP_VAR*             var                 /**< variable to get objective value for */
      ) const = 0;

   /** adds parameters to the SCIP data structure */
   virtual SCIP_RETCODE addParameters() = 0;

   /** returns whether optimal pricing can be aborted */
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
      int                   nsolvedprobs,       /**< number of pricing problems solved so far */
      int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
      SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
      int                   npricingprobsnotnull /**< number of relevant pricing problems */
      ) const = 0;

   /** returns whether heuristic pricing can be aborted */
   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
      int                   nsolvedprobs,       /**< number of pricing problems solved so far */
      int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
      SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
      int                   npricingprobsnotnull /**< number of relevant pricing problems */
      ) const = 0;

   /** starts the clock */
   virtual SCIP_RETCODE startClock();

   /** stops the clock */
   virtual SCIP_RETCODE stopClock();

   /** returns the time of the clock */
   virtual SCIP_Real getClockTime() const;

   /** returns the maximal number of rounds */
   int getMaxrounds() const
   {
      return maxrounds;
   }

   /** returns the maximal number of columns per pricing round at root node */
   int getMaxcolsroundroot() const
   {
      return maxcolsroundroot;
   }

   /** returns the maximal number of columns per pricing round */
   int getMaxcolsround() const
   {
      return maxcolsround;
   }

   /** returns the maximal number of columns per problem to be generated during pricing at root node */
   int getMaxcolsprobroot() const
   {
      return maxcolsprobroot;
   }

   /** returns the maximal number of columns per problem to be generated during pricing */
   int getMaxcolsprob() const
   {
      return maxcolsprob;
   }

   /** returns the maximal number of successfully solved pricing problems */
   int getMaxsuccessfulprobs() const
   {
      return maxsuccessfulprobs;
   }

   /** returns the maximal percentage of pricing problems that are solved at root node if variables have already been found */
   SCIP_Real getRelmaxprobsroot() const
   {
      return relmaxprobsroot;
   }

   /** returns the maximal percentage of pricing problems that are solved if variables have already been found */
   SCIP_Real getRelmaxprobs() const
   {
      return relmaxprobs;
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

   /** increases the number of calls */
   virtual void incCalls()
   {
      calls++;
   }

   /** resets the number of calls and the clock for a restart */
   SCIP_RETCODE resetCalls()
   {
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
      SCIP*                 p_scip
      );

    /** destructor */
    virtual ~ReducedCostPricing() {}

    virtual SCIP_RETCODE addParameters();

    virtual SCIP_Real consGetDual(
      SCIP*                 scip, 
      SCIP_CONS*            cons
      ) const;

    virtual SCIP_Real rowGetDual(
      SCIP_ROW*             row
      ) const;

    virtual SCIP_Real varGetObj(
      SCIP_VAR*             var
      ) const ;

   /** returns whether optimal pricing can be aborted */
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
      int                   nsolvedprobs,       /**< number of pricing problems solved so far */
      int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
      SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
      int                   npricingprobsnotnull /**< number of relevant pricing problems */
      ) const;

   /** returns whether heuristic pricing can be aborted */
   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
      int                   nsolvedprobs,       /**< number of pricing problems solved so far */
      int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
      SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
      int                   npricingprobsnotnull /**< number of relevant pricing problems */
      ) const;
};

class FarkasPricing : public PricingType
{
public:
   /** constructor */
   FarkasPricing();

   FarkasPricing(
      SCIP*                 p_scip
      );

   /** destructor */
   virtual ~FarkasPricing() {}

   virtual SCIP_RETCODE addParameters();

   virtual SCIP_Real consGetDual(
      SCIP*                 scip, 
      SCIP_CONS*            cons
      ) const;

   virtual SCIP_Real rowGetDual(
      SCIP_ROW*             row
      ) const;

   virtual SCIP_Real varGetObj(
      SCIP_VAR*             var
      ) const;

   /** returns whether optimal pricing can be aborted */
   virtual SCIP_Bool canOptimalPricingBeAborted(
      int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
      int                   nsolvedprobs,       /**< number of pricing problems solved so far */
      int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
      SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
      int                   npricingprobsnotnull /**< number of relevant pricing problems */
      ) const;

   /** returns whether heuristic pricing can be aborted */
   virtual SCIP_Bool canHeuristicPricingBeAborted(
      int                   nfoundcols,         /**< number of negative reduced cost columns found so far */
      int                   nsolvedprobs,       /**< number of pricing problems solved so far */
      int                   nsuccessfulprobs,   /**< number of pricing problems solved successfully so far */
      SCIP_Real             relmaxsuccessfulprobs, /**< maximal percentage of pricing problems that need to be solved successfully */
      int                   npricingprobsnotnull /**< number of relevant pricing problems */
      ) const;
};

#endif /* CLASS_PRICINGTYPE_H_ */
