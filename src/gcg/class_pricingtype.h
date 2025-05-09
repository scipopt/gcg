/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
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
#include "gcg/pricer_gcg.h"
#include "gcg/type_extendedmasterconsdata.h"

/**
 * @ingroup PRICING_PRIV
 * @{
 */

class PricingType
{
protected:
   GCG*                  gcg;                   /**< GCG data structure */
   SCIP*                 masterprob;            /**< SCIP instance (master problem) */
   GCG_PRICETYPE         type;                  /**< type of pricing */
   SCIP_CLOCK*           clock;                 /**< CPU clock */

   int                   calls;                 /**< number of times this type of pricing was called */
   int                   maxrounds;             /**< maximal number of pricing rounds */
   int                   maxcolsroundroot;      /**< maximal number of columns per pricing round at root node */
   int                   maxcolsround;          /**< maximal number of columns per pricing round */
   int                   maxcolsprobroot;       /**< maximal number of columns per problem to be generated at root node */
   int                   maxcolsprob;           /**< maximal number of columns per problem to be generated */
   int                   maxsuccessfulprobs;    /**< maximal number of successfully solved pricing problems until pricing
                                                  *  loop is aborted */
   SCIP_Real             relmaxprobsroot;       /**< maximal percentage of pricing problems that are solved at root node if
                                                  *  variables have already been found */
   SCIP_Real             relmaxprobs;           /**< maximal percentage of pricing problems that are solved if variables
                                                  *  have already been found */
   SCIP_Real             relmaxsuccessfulprobs; /**< maximal percentage of pricing problems that need to be solved successfully */

public:
   /** constructor */
   PricingType() = delete;

   PricingType(
      GCG*                   gcgstruct
      );

   /** destructor */
   virtual ~PricingType();

   /** get dual value of a constraint */
   virtual SCIP_Real consGetDual(
      SCIP_CONS*            cons                /**< constraint to get dual for */
      ) const = 0;

   /** get dual value of a row */
   virtual SCIP_Real rowGetDual(
      SCIP_ROW*             row                 /**< row to get dual value for */
      ) const = 0;

   /** get dual value of an extended master cons */
   virtual SCIP_Real extendedmasterconsGetDual(
      GCG_EXTENDEDMASTERCONSDATA*    extendedmasterconsdata       /**< extended master cons data */
      ) const = 0;

   /** get objective value of variable */
   virtual SCIP_Real varGetObj(
      SCIP_VAR*             var                 /**< variable to get objective value for */
      ) const = 0;

   /** adds parameters to the SCIP data structure */
   virtual SCIP_RETCODE addParameters() = 0;

   /** starts the clock */
   virtual SCIP_RETCODE startClock();

   /** stops the clock */
   virtual SCIP_RETCODE stopClock();

   /** returns the time of the clock */
   virtual SCIP_Real getClockTime() const;

   /** returns the maximal number of rounds */
   virtual int getMaxrounds() const
   {
      return maxrounds;
   }

   /** returns the maximal number of columns per pricing round */
   virtual int getMaxcolsround() const = 0;

   /** returns the maximal number of columns per problem to be generated during pricing */
   virtual int getMaxcolsprob() const = 0;

   /** returns the maximal number of successfully solved pricing problems */
   int getMaxsuccessfulprobs() const
   {
      return maxsuccessfulprobs;
   }

   /** returns the maximal percentage of pricing problems that are solved if variables have already been found */
   virtual SCIP_Real getRelmaxprobs() const = 0;

   /** returns the maximal percentage of pricing problems that need to be solved successfully */
   SCIP_Real getRelmaxsuccessfulprobs() const
   {
      return relmaxsuccessfulprobs;
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
      SCIP_CALL( SCIPresetClock(masterprob, clock) );
      return SCIP_OKAY;
   }

};

class ReducedCostPricing : public PricingType
{
public:
   /** constructor */
   ReducedCostPricing() = delete;

   ReducedCostPricing(
      GCG*                  gcgstruct
      );

   /** destructor */
   virtual ~ReducedCostPricing() {}

   virtual SCIP_RETCODE addParameters();

   virtual SCIP_Real consGetDual(
     SCIP_CONS*            cons
     ) const;

   virtual SCIP_Real rowGetDual(
     SCIP_ROW*             row
     ) const;

   virtual SCIP_Real extendedmasterconsGetDual(
     GCG_EXTENDEDMASTERCONSDATA*    extendedmasterconsdata
     ) const;

   virtual SCIP_Real varGetObj(
     SCIP_VAR*             var
     ) const ;

   /** returns the maximal number of columns per pricing round */
   virtual int getMaxcolsround() const;

   /** returns the maximal number of columns per problem to be generated during pricing */
   virtual int getMaxcolsprob() const;

   /** returns the maximal percentage of pricing problems that are solved if variables have already been found */
   virtual SCIP_Real getRelmaxprobs() const;
};

class FarkasPricing : public PricingType
{
public:
   /** constructor */
   FarkasPricing() = delete;

   FarkasPricing(
      GCG*                  gcgstruct
      );

   /** destructor */
   virtual ~FarkasPricing() {}

   virtual SCIP_RETCODE addParameters();

   virtual SCIP_Real consGetDual(
      SCIP_CONS*            cons
      ) const;

   virtual SCIP_Real rowGetDual(
      SCIP_ROW*             row
      ) const;

   virtual SCIP_Real extendedmasterconsGetDual(
      GCG_EXTENDEDMASTERCONSDATA*    extendedmasterconsdata
      ) const;

   virtual SCIP_Real varGetObj(
      SCIP_VAR*             var
      ) const;

   /** returns the maximal number of columns per pricing round */
   virtual int getMaxcolsround() const;

   /** returns the maximal number of columns per problem to be generated during pricing */
   virtual int getMaxcolsprob() const;

   /** returns the maximal percentage of pricing problems that are solved if variables have already been found */
   virtual SCIP_Real getRelmaxprobs() const;
};

/** @} */
#endif /* CLASS_PRICINGTYPE_H_ */
