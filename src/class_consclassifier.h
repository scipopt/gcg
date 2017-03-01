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

#ifndef GCG_CLASS_CONSCLASSIFIER_H_
#define GCG_CLASS_CONSCLASSIFIER_H_

#include <vector>
#include <string>

namespace gcg
{

class ConsClassifier
{

private:
   SCIP*                      scip;                   /**< scip data structure */
   int                        nClasses;               /**< number of classes the classifier provides */
   int                        nConss;                 /**< number of constraints */
   std::vector<int>           consToClasses;          /**< constraint i is assigned to class consToClasses[i] */
   std::string                name;                   /**< name of the classifier */
   std::vector<std::string>   classNames;             /**< name of class k is classNames[k] */
   std::vector<int>           classDecompInfo;        /**< encodes whether a class k should be assigned to master or pricing problem */

public:

   /** constructor */
   ConsClassifier(
      SCIP* scip,                /** scip data structure */
      std::string name,          /** name of classifier */
      int nConss                 /** number of constraints to be classified */
   );

   /** destructor */
   ~ConsClassifier();


   /** creates a new class, returns number of class */
   int addClass(
      std::string name          /**< name of the class */
      // ??? int classDecompInfo        /**< decomposition code of the class */
   );

   /** assigns a constraint to a class */
   void assignConsToClass(
      int consindex,            /**< index of the constraint */
      int classindex            /**< index of the class */
   );


   /** returns the number of classes the classifier provides */
   int getNClasses(
   );


   /** returns vector containing the assigned class of each constraint */
   const int* getConsToClasses(
   );

   /** returns the index of the class a constraint is assigned to */
   int getClassOfCons(
      int consindex              /**< index of constraint */
   );


   /** returns the name of the classifier */
   const std::string getName(
   );


   const std::string* getClassNames(
   );

   /** returns the name of the class a constraint is assigned to */
   const std::string getClassNameOfCons(
      int consindex              /**< index of constraint */
   );


   /** returns the decomposition code */
   const int* getClassDecompInfo(
   );

   /** returns the decomposition code of a class */
   int getClassDecompInfoOfClass(
      int classindex             /**< index of class */
   );


   /** set the decomposition code of a class */
   void setClassDecompInfo(
      int classindex,            /**< index of class */
      int decompInfo             /**< decomposition code of class */
   );

};


} /* namespace gcg */
#endif /* GCG_CLASS_CONSCLASSIFIER_H_ */
