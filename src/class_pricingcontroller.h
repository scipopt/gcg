/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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

/**@file   class_pricingcontroller.h
 * @brief  pricing controller managing the pricing strategy
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_PRICINGCONTROLLER_H_
#define CLASS_PRICINGCONTROLLER_H_

#include "class_colpool.h"
#include "class_pricingtype.h"
#include "type_gcgpqueue.h"
#include "type_pricingjob.h"
#include "objscip/objscip.h"

namespace gcg {

class Pricingcontroller
{

private:
   SCIP*                 scip_;              /**< SCIP instance (master problem) */
   GCG_PRICINGJOB**      pricingjobs;        /**< pricing jobs, one per pricing problem */
   int                   npricingprobs;      /**< number of pricing problems */

   /* parameters */
   SCIP_Bool             useheurpricing;     /**< should heuristic pricing be used? */
   int                   sorting;            /**< how should pricing problems be sorted */
   SCIP_Real             successfulmipsrel;  /**< factor of MIPs to be solved successfully until pricing is aborted */
   int                   eagerfreq;          /**< frequency at which all pricingproblems should be solved */

   /* strategy */
   GCG_PQUEUE*           pqueue;             /**< priority queue containing the pricing jobs */
   SCIP_Real*            score;              /**< scores of the pricing problems */
   PricingType*          pricingtype_;       /**< current pricing type */

   /* statistics */
   int                   eagerage;           /**< iterations since last eager iteration */


public:
   /** constructor */
   Pricingcontroller(SCIP* scip);

   /** destructor */
   virtual ~Pricingcontroller();

   SCIP_RETCODE addParameters();

   SCIP_RETCODE initSol();

   SCIP_RETCODE exitSol();

   /** pricing initialization, called right at the beginning of pricing */
   void initPricing(
      PricingType*          pricingtype         /**< type of pricing */
   );

   /** pricing deinitialization, called when pricing is finished */
   void exitPricing();

   /** setup the priority queue (done once per stabilization round): add all pricing jobs to be performed */
   SCIP_RETCODE setupPriorityQueue(
      SCIP_Real*            dualsolconv,        /**< dual solution values / Farkas coefficients of convexity constraints */
      int                   maxcols             /**< maximum number of columns to be generated */
      );

   /** get the next pricing job to be performed */
   GCG_PRICINGJOB* getNextPricingjob();

   /** update statistics of a pricing job, and possibly add it again to the queue with different settings */
   void updatePricingjob(
      GCG_PRICINGJOB*       pricingjob,         /**< pricing job */
      SCIP_STATUS           status,             /**< status after solving the pricing problem */
      SCIP_Real             lowerbound,         /**< lower bound returned by the pricing problem */
      GCG_COL**             cols,               /**< columns found by the last solving of the pricing problem */
      int                   ncols               /**< number of columns found */
      );

   /** decide whether a pricing job must be treated again */
   void evaluatePricingjob(
      GCG_PRICINGJOB*       pricingjob         /**< pricing job */
      );

   /** return whether the reduced cost is valid */
   SCIP_Bool redcostIsValid();

   /** reset the lower bound of a pricing job */
   void resetPricingjobLowerbound(
      GCG_PRICINGJOB*       pricingjob          /**< pricing job */
      );

   /** for all pricing jobs, move their columns to the column pool */
   SCIP_RETCODE moveColsToColpool(
      Colpool*           colpool             /**< column pool */
      );

   /** get best columns found by the pricing jobs */
   void getBestCols(
      GCG_COL**             cols                /**< column array to be filled */
      );

   /** free all columns of the pricing jobs */
   void freeCols();

   /** returns whether pricing can be aborted */
   SCIP_Bool abortPricing(
      PricingType*          pricetype,          /**< type of pricing */
      int                   nfoundvars,         /**< number of variables found so far */
      int                   solvedmips,         /**< number of MIPS solved so far */
      int                   successfulmips,     /**< number of successful mips solved so far */
      SCIP_Bool             optimal             /**< optimal or heuristic pricing */
      ) const;

   void resetEagerage();

   void increaseEagerage();


private:
   /** comparison operator for pricing jobs w.r.t. their solution priority */
   static
   SCIP_DECL_SORTPTRCOMP(comparePricingjobs);
};

} /* namespace gcg */
#endif /* CLASS_PRICINGCONTROLLER_H_ */
