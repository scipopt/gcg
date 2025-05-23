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

/**@file   class_varpartition.h
 * @brief  class representing a partition of a set of variables
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_VARPARTITION_H__
#define GCG_CLASS_VARPARTITION_H__


#include "gcg/class_indexpartition.h"

namespace gcg
{

/** @ingroup PARTITIONS
 * @{
 */

enum VarClassDecompInfo
{
   ALL = 0,
   LINKING = 1,
   MASTER = 2,
   BLOCK = 3
};
typedef enum VarClassDecompInfo VAR_DECOMPINFO;


class VarPartition : public IndexPartition
{

public:

   /** constructor */
   GCG_EXPORT
   VarPartition(
      GCG*                 gcgstruct,           /**< GCG data structure */
      const char*          name,                /**< name of partition (will be copied) */
      int                  nClasses,            /**< initial number of classes */
      int                  nVars                /**< number of variables to be classified */
   );

   /** copy constructor */
   GCG_EXPORT
   VarPartition(
      const VarPartition* toCopy              /**< VarPartition to be copied */
   );


   /** destructor */
   GCG_EXPORT
   ~VarPartition();


   /** creates a new class, returns index of the class */
   GCG_EXPORT
   int addClass(
      const char* name,                /**< name of the class (will be copied) */
      const char* desc,                /**< description of the class (will be copied) */
      VAR_DECOMPINFO decompInfo        /**< decomposition code of the class */
   );

   /** assigns a variable to a class */
   GCG_EXPORT
   void assignVarToClass(
      int varindex,                    /**< index of the variable */
      int classindex                   /**< index of the class */
   );

   /** returns a vector containing all possible subsets of the chosen classindices */
   GCG_EXPORT
   std::vector<std::vector<int>> getAllSubsets(
      bool all,                        /**< true, if ALL classes should be considered */
      bool linking,                    /**< true, if LINKING classes should be considered */
      bool master,                     /**< true, if MASTER classes should be considered */
      bool block                       /**< true, if BLOCK classes should be considered */
   );

   /** returns the decomposition info of a class */
   GCG_EXPORT
   VAR_DECOMPINFO getClassDecompInfo(
      int classindex                   /**< index of class */
   );

   /** returns the name of the class a variable is assigned to */
   GCG_EXPORT
   const char* getClassNameOfVar(
      int varindex                    /**< index of variable */
   );


   /** returns the index of the class a variable is assigned to */
   GCG_EXPORT
   int getClassOfVar(
      int varindex                    /**< index of variable */
   );

   /** returns vector containing the assigned class of each variable */
   GCG_EXPORT
   const int* getVarsToClasses(
   );

   /** returns the number of variables */
   GCG_EXPORT
   int getNVars(
   );

   /** returns a vector with the numbers of variables that are assigned to the classes */
   GCG_EXPORT
   std::vector<int> getNVarsOfClasses(
   );

   /** returns whether a variable is already assigned to a class */
   GCG_EXPORT
   bool isVarClassified(
      int varindex                    /**< index of variable */
   );


   /** returns partition with reduced number of classes
    *  if the current number of classes is greater than an upper bound
    *  and lower than 2*(upper bound) (returns NULL otherwise) */
   GCG_EXPORT
   VarPartition* reduceClasses(
      int maxNumberOfClasses           /**< upper bound */
   );

   /** sets the decomposition code of a class */
   GCG_EXPORT
   void setClassDecompInfo(
      int classindex,                  /**< index of class */
      VAR_DECOMPINFO decompInfo        /**< decomposition code of class */
   );

};

/**@} */
} /* namespace gcg */
#endif /* GCG_CLASS_VARPARTITION_H__ */
