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

ConsClassifier::ConsClassifier(
   SCIP*          _scip,
   std::string    givenName,
   int            givenNCons
) :
   scip(_scip), name(givenName), nClasses(0), nConss(givenNCons), consToClasses(givenNCons, -1), classNames(0), classDecompInfo(0)
{
}

ConsClassifier::~ConsClassifier()
{
}


int ConsClassifier::addClass( std::string givenName )
{
   assert((int) classNames.size() == nClasses);

   classNames.push_back( givenName );
   classDecompInfo.push_back( -1 );

   ++nClasses;

   return nClasses - 1;
}

void ConsClassifier::assignConsToClass( int givenConsindex, int givenClassindex )
{
   assert(0 <= givenConsindex && givenConsindex <= nConss);
   assert(0 <= givenClassindex && givenClassindex <= nClasses);

   consToClasses[givenConsindex] = givenClassindex;
}


int ConsClassifier::getNClasses()
{
   return nClasses;
}


const int* ConsClassifier::getConsToClasses()
{
   return &consToClasses[0];
}

int ConsClassifier::getClassOfCons( int givenConsindex )
{
   assert(0 <= givenConsindex && givenConsindex <= nConss);

   return consToClasses[givenConsindex];
}


const std::string ConsClassifier::getName()
{
   return name;
}


const std::string* ConsClassifier::getClassNames()
{
   return &classNames[0];
}

const std::string ConsClassifier::getClassNameOfCons( int givenConsindex )
{
   assert(0 <= givenConsindex && givenConsindex <= nConss);

   return classNames[givenConsindex];
}


const int* ConsClassifier::getClassDecompInfo()
{
   return &classDecompInfo[0];
}

int ConsClassifier::getClassDecompInfoOfClass( int givenClassindex )
{
   assert(0 <= givenClassindex && givenClassindex <= nClasses);

   return classDecompInfo[givenClassindex];
}


void ConsClassifier::setClassDecompInfo( int givenClassindex, int givenDecompInfo )
{
   assert(0 <= givenClassindex && givenClassindex <= nClasses);

   classDecompInfo[givenClassindex] = givenDecompInfo;
}

} /* namespace gcg */
