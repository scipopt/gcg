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

/**@file   class_consclassifier.cpp
 * @brief  class for classifying constraints
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_consclassifier.h"

#include <assert.h>
#include <sstream>
#include <algorithm>


namespace gcg {


/** constructor */
ConsClassifier::ConsClassifier(
   SCIP*          _scip,
   const char*    givenName,
   int            givenNClasses,
   int            givenNCons
) :
   IndexClassifier(_scip, givenName, givenNClasses, givenNCons)
{

}

/** copy constructor */
ConsClassifier::ConsClassifier(
   const ConsClassifier* toCopy
) : IndexClassifier( toCopy )
{
}

/** destructor */
ConsClassifier::~ConsClassifier()
{
}

/** creates a new class, returns index of the class */
int ConsClassifier::addClass( const char* givenName, const char* givenDesc, CONS_DECOMPINFO givenDecompInfo )
{
   int classindex = IndexClassifier::addClass( givenName, givenDesc );
   setClassDecompInfo( classindex, givenDecompInfo );

   return classindex;
}

/** assigns a constraint to a class */
void ConsClassifier::assignConsToClass( int givenConsindex, int givenClassindex )
{
   IndexClassifier::assignIndexToClass( givenConsindex, givenClassindex );
}

/** returns a vector containing all possible subsets of the chosen classindices */
std::vector<std::vector<int>> ConsClassifier::getAllSubsets( bool both, bool only_master, bool only_pricing )
{
   std::vector<int> classindices;
   for ( int i = 0; i < getNClasses(); ++i )
   {
      if ( ( both && getClassDecompInfo( i ) == BOTH ) || ( only_master && getClassDecompInfo( i ) == ONLY_MASTER )
            || ( only_pricing && getClassDecompInfo( i ) == ONLY_PRICING ) )
         classindices.push_back( i );
   }
   return IndexClassifier::getAllSubsets( classindices );
}

/** returns the decomposition code of a class */
CONS_DECOMPINFO ConsClassifier::getClassDecompInfo( int givenClassindex )
{
   int decompInfo = IndexClassifier::getClassDecompInfo( givenClassindex );
   CONS_DECOMPINFO interp;

   assert( 0 <= decompInfo && decompInfo <= 2);

   switch ( decompInfo )
   {
   case 0:
      interp = BOTH;
      break;
   case 1:
      interp = ONLY_MASTER;
      break;
   case 2:
      interp = ONLY_PRICING;
      break;
   default:
      interp = BOTH;
      break;
   }

   return interp;
}

/** returns the name of the class a constraint is assigned to */
const char* ConsClassifier::getClassNameOfCons( int givenConsindex )
{
   return IndexClassifier::getClassNameOfIndex( givenConsindex );
}

/** returns the index of the class a constraint is assigned to */
int ConsClassifier::getClassOfCons( int givenConsindex )
{
   return IndexClassifier::getClassOfIndex( givenConsindex );
}

/** returns vector containing the assigned class of each constraint */
const int* ConsClassifier::getConssToClasses()
{
   std::vector<int> conssToClasses = IndexClassifier::getIndicesToClasses();
   if ( conssToClasses.size() > 0 )
      return &conssToClasses[0];
   else
      return NULL;
}

/** returns the number of constraints */
int ConsClassifier::getNConss()
{
   return IndexClassifier::getNIndices();
}

/** returns a vector with the numbers of constraints that are assigned to the classes */
std::vector<int> ConsClassifier::getNConssOfClasses()
{
   return IndexClassifier::getNIndicesOfClasses();
}


/** returns whether a constraint is already assigned to a class */
bool ConsClassifier::isConsClassified( int givenConsindex )
{
   return IndexClassifier::isIndexClassified( givenConsindex );
}

/** returns classifier with reduced number of classes */
ConsClassifier* ConsClassifier::reduceClasses( int givenMaxNumber )
{
   std::vector<int> classindexmapping = IndexClassifier::reduceClasses( givenMaxNumber );
   ConsClassifier* newClassifier;
   std::stringstream newName;
   std::stringstream newClassdesc;

   if ( classindexmapping.empty() )
      return NULL;

   /** create new ConsClassifier */
   newName << getName() << "-red-to-" << givenMaxNumber;
   newClassifier = new ConsClassifier( scip, newName.str().c_str(), givenMaxNumber, getNConss() );

   /** reassign conss */
   for( int i = 0; i < newClassifier->getNConss(); ++i)
   {
      if ( getClassOfCons(i) != -1 )
      {
         newClassifier->assignConsToClass( i, classindexmapping[getClassOfCons(i)] );
      }
   }

   /** set new class names and descriptions (enlarged class has index 0) */
   newClassifier->setClassName( 0, "merged" );
   newClassifier->setClassDecompInfo( 0, BOTH );

   for ( int i = 0; i < getNClasses(); ++i )
   {
     if ( classindexmapping[i] == 0 )
     {
        newClassdesc << getClassDescription( i ) << " - ";
     }
     else
     {
        newClassifier->setClassName( classindexmapping[i], getClassName(i) );
        newClassifier->setClassDescription( classindexmapping[i], getClassDescription(i) );
        newClassifier->setClassDecompInfo( classindexmapping[i], getClassDecompInfo(i) );
     }
   }

   newClassifier->setClassDescription( 0, newClassdesc.str().c_str() );


   return newClassifier;
}

/** sets the decomposition code of a class */
void ConsClassifier::setClassDecompInfo( int givenClassindex, CONS_DECOMPINFO givenDecompInfo )
{
   assert(givenDecompInfo == BOTH || givenDecompInfo == ONLY_MASTER || givenDecompInfo == ONLY_PRICING );

   IndexClassifier::setClassDecompInfo( givenClassindex, (int) givenDecompInfo );
}

} /* namespace gcg */
