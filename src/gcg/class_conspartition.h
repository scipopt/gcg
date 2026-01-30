/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   class_conspartition.h
 * @brief  class representing a partition of a set of constraints
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_CONSPARTITION_H__
#define GCG_CLASS_CONSPARTITION_H__


#include "gcg/class_indexpartition.h"

namespace gcg
{

/** @ingroup PARTITIONS
 * @{
 */

enum ConsClassDecompInfo
{
   BOTH = 0,                     /**< assign class to master or pricing problem */
   ONLY_MASTER = 1,              /**< assign class only to master problem */
   ONLY_PRICING = 2              /**< assign class only to pricing problem */
};
typedef enum ConsClassDecompInfo CONS_DECOMPINFO;


class ConsPartition : public IndexPartition
{

public:

   /** constructor */
   GCG_EXPORT
   ConsPartition(
      GCG*                 gcgstruct,           /**< GCG data structure */
      const char*          name,                /**< name of partition (will be copied) */
      int                  nClasses,            /**< initial number of classes */
      int                  nConss               /**< number of constraints to be classified */
   );

   /** copy constructor */
   GCG_EXPORT
   ConsPartition(
      const ConsPartition* toCopy              /**< ConsPartition to be copied */
   );


   /** destructor */
   GCG_EXPORT
   ~ConsPartition();


   /** creates a new class, returns index of the class */
   GCG_EXPORT
   int addClass(
      const char* name,                /**< name of the class (will be copied) */
      const char* desc,                /**< description of the class (will be copied) */
      CONS_DECOMPINFO decompInfo            /**< decomposition code of the class */
   );

   /** assigns a constraint to a class */
   GCG_EXPORT
   void assignConsToClass(
      int consindex,                   /**< index of the constraint */
      int classindex                   /**< index of the class */
   );

   /** returns a vector containing all possible subsets of the chosen classindices */
   GCG_EXPORT
   std::vector<std::vector<int>> getAllSubsets(
      bool both,                       /**< true, if BOTH classes should be considered */
      bool only_master,                /**< true, if ONLY_MASTER classes should be considered */
      bool only_pricing                /**< true, if ONLY_PRICING classes should be considered */
   );

   /** returns the decomposition info of a class */
   GCG_EXPORT
   CONS_DECOMPINFO getClassDecompInfo(
      int classindex                   /**< index of class */
   );

   /** returns the name of the class a constraint is assigned to */
   GCG_EXPORT
   const char* getClassNameOfCons(
      int consindex                    /**< index of constraint */
   );


   /** returns the index of the class a constraint is assigned to */
   GCG_EXPORT
   int getClassOfCons(
      int consindex                    /**< index of constraint */
   );

   /** returns vector containing the assigned class of each constraint */
   GCG_EXPORT
   const int* getConssToClasses(
   );

   /** returns the number of constraints */
   GCG_EXPORT
   int getNConss(
   );

   /** returns a vector with the numbers of constraints that are assigned to the classes */
   GCG_EXPORT
   std::vector<int> getNConssOfClasses(
   );

   /** returns whether a constraint is already assigned to a class */
   GCG_EXPORT
   bool isConsClassified(
      int consindex                    /**< index of constraint */
   );


   /** returns partition with reduced number of classes
    *  if the current number of classes is greater than an upper bound
    *  and lower than 2*(upper bound) (returns NULL otherwise) */
   GCG_EXPORT
   ConsPartition* reduceClasses(
      int maxNumberOfClasses           /**< upper bound */
   );

   /** sets the decomposition code of a class */
   GCG_EXPORT
   void setClassDecompInfo(
      int classindex,                  /**< index of class */
      CONS_DECOMPINFO decompInfo            /**< decomposition code of class */
   );

};

/** @} */
} /* namespace gcg */
#endif /* GCG_CLASS_CONSPARTITION_H__ */
