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

/**@file   class_varclassifier.cpp
 * @brief  class for classifying variables
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_varclassifier.h"

#include <assert.h>
#include <sstream>
#include <algorithm>


namespace gcg {


/** constructor */
VarClassifier::VarClassifier(
   SCIP*          _scip,
   const char*    givenName,
   int            givenNClasses,
   int            givenNVars
) :
   IndexClassifier(_scip, givenName, givenNClasses, givenNVars)
{

}

/** copy constructor */
VarClassifier::VarClassifier(
   const VarClassifier* toCopy
) : IndexClassifier( toCopy )
{
}

/** destructor */
VarClassifier::~VarClassifier()
{
}

/** creates a new class, returns index of the class */
int VarClassifier::addClass( const char* givenName, const char* givenDesc, VAR_DECOMPINFO givenDecompInfo )
{
   int classindex = IndexClassifier::addClass( givenName, givenDesc );
   setClassDecompInfo( classindex, givenDecompInfo );

   return classindex;
}

/** assigns a variable to a class */
void VarClassifier::assignVarToClass( int givenVarindex, int givenClassindex )
{
   IndexClassifier::assignIndexToClass( givenVarindex, givenClassindex );
}

/** returns a vector containing all possible subsets of the chosen classindices */
std::vector<std::vector<int>> VarClassifier::getAllSubsets( bool all, bool linking, bool master, bool block )
{
   std::vector<int> classindices;
   for ( int i = 0; i < getNClasses(); ++i )
   {
      if ( ( all && getClassDecompInfo( i ) == ALL ) ||
           ( linking && getClassDecompInfo( i ) == LINKING ) ||
           ( master && getClassDecompInfo( i ) == MASTER ) ||
           ( block && getClassDecompInfo( i ) == BLOCK )
          )
         classindices.push_back( i );
   }
   return IndexClassifier::getAllSubsets( classindices );
}

/** returns the decomposition code of a class */
VAR_DECOMPINFO VarClassifier::getClassDecompInfo( int givenClassindex )
{
   int decompInfo = IndexClassifier::getClassDecompInfo( givenClassindex );
   VAR_DECOMPINFO interp;

   assert( 0 <= decompInfo && decompInfo <= 3);

   switch ( decompInfo )
   {
   case 0:
      interp = ALL;
      break;
   case 1:
      interp = LINKING;
      break;
   case 2:
      interp = MASTER;
      break;
   case 3:
      interp = BLOCK;
      break;
   default:
      interp = ALL;
      break;
   }

   return interp;
}

/** returns the name of the class a variable is assigned to */
const char* VarClassifier::getClassNameOfVar( int givenVarindex )
{
   return IndexClassifier::getClassNameOfIndex( givenVarindex );
}

/** returns the index of the class a variable is assigned to */
int VarClassifier::getClassOfVar( int givenVarindex )
{
   return IndexClassifier::getClassOfIndex( givenVarindex );
}

/** returns vector containing the assigned class of each variable */
const int* VarClassifier::getVarsToClasses()
{
   std::vector<int> varsToClasses = IndexClassifier::getIndicesToClasses();
   if ( varsToClasses.size() > 0 )
      return &varsToClasses[0];
   else
      return NULL;
}

/** returns the number of variables */
int VarClassifier::getNVars()
{
   return IndexClassifier::getNIndices();
}

/** returns a vector with the numbers of variables that are assigned to the classes */
std::vector<int> VarClassifier::getNVarsOfClasses()
{
   return IndexClassifier::getNIndicesOfClasses();
}


/** returns whether a variable is already assigned to a class */
bool VarClassifier::isVarClassified( int givenVarindex )
{
   return IndexClassifier::isIndexClassified( givenVarindex );
}

/** returns classifier with reduced number of classes */
VarClassifier* VarClassifier::reduceClasses( int givenMaxNumber )
{
   std::vector<int> classindexmapping = IndexClassifier::reduceClasses( givenMaxNumber );
   VarClassifier* newClassifier;
   std::stringstream newName;
   std::stringstream newClassdesc;

   if ( classindexmapping.empty() )
      return NULL;

   /** create new VarClassifier */
   newName << getName() << "-red-to-" << givenMaxNumber;
   newClassifier = new VarClassifier( scip, newName.str().c_str(), givenMaxNumber, getNVars() );

   /** reassign vars */
   for( int i = 0; i < newClassifier->getNVars(); ++i)
   {
      if ( getClassOfVar(i) != -1 )
      {
         newClassifier->assignVarToClass( i, classindexmapping[getClassOfVar(i)] );
      }
   }

   /** set new class names and descriptions (enlarged class has index 0) */
   newClassifier->setClassName( 0, "merged" );
   newClassifier->setClassDecompInfo( 0, ALL );

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
void VarClassifier::setClassDecompInfo( int givenClassindex, VAR_DECOMPINFO givenDecompInfo )
{
   assert( givenDecompInfo == ALL || givenDecompInfo == LINKING || givenDecompInfo == MASTER || givenDecompInfo == BLOCK );

   IndexClassifier::setClassDecompInfo( givenClassindex, (int) givenDecompInfo );
}

} /* namespace gcg */
