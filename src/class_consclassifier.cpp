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

#include "class_consclassifier.h"

#include <assert.h>


namespace gcg {


/** constructor */
ConsClassifier::ConsClassifier(
   SCIP*          _scip,
   const char*    givenName,
   int            givenNClasses,
   int            givenNCons
) :
   scip(_scip),  name(std::string(givenName)), nClasses(givenNClasses), nConss(givenNCons), consToClasses(givenNCons, -1),
   classNames(givenNClasses, ""), classDescriptions(givenNClasses, ""), classDecompInfo(givenNClasses, BOTH)
{

}

/** destructor */
ConsClassifier::~ConsClassifier()
{
}

/** creates a new class, returns index of the class */
int ConsClassifier::addClass( const char* givenName, const char* givenDesc, DECOMPINFO givenDecompInfo )
{
   assert((int) classNames.size() == nClasses);
   assert((int) classDescriptions.size() == nClasses);
   assert((int) classDecompInfo.size() == nClasses);
   assert(givenDecompInfo == BOTH || givenDecompInfo == ONLY_MASTER || givenDecompInfo == ONLY_PRICING );

   std::string givenClassName = std::string(givenName);

   classNames.push_back( givenClassName );
   classDescriptions.push_back( givenDesc );
   classDecompInfo.push_back( givenDecompInfo );

   ++nClasses;

   return nClasses - 1;
}

/** assigns a constraint to a class */
void ConsClassifier::assignConsToClass( int givenConsindex, int givenClassindex )
{
   assert(0 <= givenConsindex && givenConsindex < nConss);
   assert(-1 <= givenClassindex && givenClassindex < nClasses);

   consToClasses[givenConsindex] = givenClassindex;
}


/** returns the decomposition code */
const DECOMPINFO* ConsClassifier::getClassDecompInfo()
{
   if ( nClasses > 0)
      return &classDecompInfo[0];
   else
      return NULL;
}

/** returns the decomposition code of a class */
DECOMPINFO ConsClassifier::getClassDecompInfoOfClass( int givenClassindex )
{
   assert(0 <= givenClassindex && givenClassindex < nClasses);

   return classDecompInfo[givenClassindex];
}

/** returns the information text of a class */
const char* ConsClassifier::getClassDescription( int givenClassindex )
{
   assert(0 <= givenClassindex && givenClassindex < nClasses);

   return classDescriptions[givenClassindex].c_str();
}

/** returns the name of a class */
const char* ConsClassifier::getClassName( int givenClassindex )
{
   assert(0 <= givenClassindex && givenClassindex < nClasses);

   return classNames[givenClassindex].c_str();
}

/** returns the name of the class a constraint is assigned to */
const char* ConsClassifier::getClassNameOfCons( int givenConsindex )
{
   assert(0 <= givenConsindex && givenConsindex < nConss);
   assert(0 <= consToClasses[givenConsindex] && consToClasses[givenConsindex] < nClasses);

   return classNames[consToClasses[givenConsindex]].c_str();
}


/** returns the index of the class a constraint is assigned to */
int ConsClassifier::getClassOfCons( int givenConsindex )
{
   assert(0 <= givenConsindex && givenConsindex < nConss);

   return consToClasses[givenConsindex];
}

/** returns vector containing the assigned class of each constraint */
const int* ConsClassifier::getConsToClasses()
{
   if ( nConss > 0 )
      return &consToClasses[0];
   else
      return NULL;
}


/** returns the name of the classifier */
const char* ConsClassifier::getName()
{
   return name.c_str();
}


/** returns the number of classes the classifier provides */
int ConsClassifier::getNClasses()
{
   return nClasses;
}

/** returns the number of constraints */
int ConsClassifier::getNConss()
{
   return nConss;
}


/** returns whether a constraint is already assigned to a class */
bool ConsClassifier::isConsClassified( int givenConsindex )
{
   assert(0 <= givenConsindex && givenConsindex < nConss);

   return consToClasses[givenConsindex] != -1;
}


/** sets the decomposition code of a class */
void ConsClassifier::setClassDecompInfo( int givenClassindex, DECOMPINFO givenDecompInfo )
{
   assert(0 <= givenClassindex && givenClassindex < nClasses);
   assert(givenDecompInfo == BOTH || givenDecompInfo == ONLY_MASTER || givenDecompInfo == ONLY_PRICING );

   classDecompInfo[givenClassindex] = givenDecompInfo;
}

/** sets the information text of a class */
void ConsClassifier::setClassDescription( int givenClassindex, const char* givenDesc )
{
   assert(0 <= givenClassindex && givenClassindex < nClasses);

   classDescriptions[givenClassindex] = std::string(givenDesc);
}

/** sets the name of a class */
void ConsClassifier::setClassName( int givenClassindex, const char* givenName )
{
   assert(0 <= givenClassindex && givenClassindex < nClasses);

   classNames[givenClassindex] = std::string(givenName);
}

} /* namespace gcg */
