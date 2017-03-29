/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   class_consclassifier.h
 * @brief  class for classifying constraints
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_CONSCLASSIFIER_H__
#define GCG_CLASS_CONSCLASSIFIER_H__

#include "objscip/objscip.h"
#include <vector>
#include <string>

namespace gcg
{

enum ClassDecompInfo
{
   BOTH = 0,                     /**< assign class to master or pricing problem */
   ONLY_MASTER = 1,              /**< assign class only to master problem */
   ONLY_PRICING = 2              /**< assign class only to pricing problem */
};
typedef enum ClassDecompInfo DECOMPINFO;


class ConsClassifier
{

private:
   SCIP*                      scip;                   /**< scip data structure */
   std::string                name;                   /**< name of the classifier */
   int                        nClasses;               /**< number of classes the classifier provides */
   int                        nConss;                 /**< number of constraints */
   std::vector<int>           consToClasses;          /**< constraint i is assigned to class consToClasses[i] (-1 if not assigned)*/
   std::vector<std::string>   classNames;             /**< name of class k is classNames[k] */
   std::vector<std::string>   classDescriptions;      /**< information text describing class k is classDescriptions[k] */
   std::vector<DECOMPINFO>    classDecompInfo;        /**< encodes whether a class k should be assigned to master or pricing problem */

public:

   /** constructor */
   ConsClassifier(
      SCIP*                scip,                /**< scip data structure */
      const char*          name,                /**< name of classifier (will be copied) */
      int                  nClasses,            /**< initial number of classes */
      int                  nConss               /**< number of constraints to be classified */
   );

   /** copy constructor */
   ConsClassifier(
      const ConsClassifier* toCopy              /**< ConsClassifier to be copied */
   );


   /** destructor */
   ~ConsClassifier();


   /** creates a new class, returns index of the class */
   int addClass(
      const char* name,                /**< name of the class (will be copied) */
      const char* desc,                /**< description of the class (will be copied) */
      DECOMPINFO decompInfo            /**< decomposition code of the class */
   );

   /** assigns a constraint to a class */
   void assignConsToClass(
      int consindex,                   /**< index of the constraint */
      int classindex                   /**< index of the class */
   );


   /** returns the decomposition info */
   const DECOMPINFO* getClassDecompInfo(
   );

   /** returns the decomposition info of a class */
   DECOMPINFO getClassDecompInfoOfClass(
      int classindex                   /**< index of class */
   );

   /** returns the information text of a class */
   const char* getClassDescription(
      int classindex                   /**< index of class */
   );

   /** returns the name of a class */
   const char* getClassName(
      int classindex                   /**< index of class */
   );

   /** returns the name of the class a constraint is assigned to */
   const char* getClassNameOfCons(
      int consindex                    /**< index of constraint */
   );


   /** returns the index of the class a constraint is assigned to */
   int getClassOfCons(
      int consindex                    /**< index of constraint */
   );

   /** returns vector containing the assigned class of each constraint */
   const int* getConsToClasses(
   );


   /** returns the name of the classifier */
   const char* getName(
   );


   /** returns the number of classes the classifier provides */
   int getNClasses(
   );

   /** returns the number of constraints */
   int getNConss(
   );


   /** returns whether a constraint is already assigned to a class */
   bool isConsClassified(
      int consindex                    /**< index of constraint */
   );


   /** returns classifier with reduced number of classes
    *  if the current number of classes is greater than an upper bound
    *  and lower than 2*(upper bound) (returns NULL otherwise) */
   ConsClassifier* reduceClasses(
      int maxNumberOfClasses           /**< upper bound */
   );

   /** removes all classes which do not have any assigned constraint (classindices may change)
    *  returns number of removed classes */
   int removeEmptyClasses(
   );

   /** sets the decomposition code of a class */
   void setClassDecompInfo(
      int classindex,                  /**< index of class */
      DECOMPINFO decompInfo            /**< decomposition code of class */
   );

   /** sets the information text of a class */
   void setClassDescription(
      int classindex,                  /**< index of class */
      const char* desc                 /**< description of class (will be copied) */
   );

   /** sets the name of a class */
   void setClassName(
      int classindex,                  /**< index of class */
      const char* name                 /**< name of class (will be copied) */
   );

};


} /* namespace gcg */
#endif /* GCG_CLASS_CONSCLASSIFIER_H__ */
