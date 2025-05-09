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

/**@file   class_conspartition.cpp
 * @brief  class representing a partition of a set of constraints
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_conspartition.h"

#include <cassert>
#include <sstream>
#include <algorithm>


namespace gcg {


/** constructor */
ConsPartition::ConsPartition(
   GCG*           gcgstruct,
   const char*    givenName,
   int            givenNClasses,
   int            givenNCons
) :
   IndexPartition(gcgstruct, givenName, givenNClasses, givenNCons)
{

}

/** copy constructor */
ConsPartition::ConsPartition(
   const ConsPartition* toCopy
) : IndexPartition(toCopy )
{
}

/** destructor */
ConsPartition::~ConsPartition()
{
}

/** creates a new class, returns index of the class */
int ConsPartition::addClass(const char* givenName, const char* givenDesc, CONS_DECOMPINFO givenDecompInfo)
{
   int classindex = IndexPartition::addClass(givenName, givenDesc);
   setClassDecompInfo(classindex, givenDecompInfo);

   return classindex;
}

/** assigns a constraint to a class */
void ConsPartition::assignConsToClass(int givenConsindex, int givenClassindex)
{
   IndexPartition::assignIndexToClass(givenConsindex, givenClassindex );
}

/** returns a vector containing all possible subsets of the chosen classindices */
std::vector<std::vector<int>> ConsPartition::getAllSubsets(bool both, bool only_master, bool only_pricing )
{
   std::vector<int> classindices;
   for( int i = 0; i < getNClasses(); ++i )
   {
      if( ( both && getClassDecompInfo( i ) == BOTH ) || ( only_master && getClassDecompInfo( i ) == ONLY_MASTER )
            || ( only_pricing && getClassDecompInfo( i ) == ONLY_PRICING ) )
         classindices.push_back( i );
   }
   return IndexPartition::getAllSubsets(classindices );
}

/** returns the decomposition code of a class */
CONS_DECOMPINFO ConsPartition::getClassDecompInfo(int givenClassindex)
{
   int decompInfo = IndexPartition::getClassDecompInfo(givenClassindex);
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
const char* ConsPartition::getClassNameOfCons(int givenConsindex)
{
   return IndexPartition::getClassNameOfIndex(givenConsindex);
}

/** returns the index of the class a constraint is assigned to */
int ConsPartition::getClassOfCons(int givenConsindex)
{
   return IndexPartition::getClassOfIndex(givenConsindex);
}

/** returns vector containing the assigned class of each constraint */
const int* ConsPartition::getConssToClasses()
{
   std::vector<int>& conssToClasses = IndexPartition::getIndicesToClasses();
   if( !conssToClasses.empty() )
      return &conssToClasses[0];
   else
      return NULL;
}

/** returns the number of constraints */
int ConsPartition::getNConss()
{
   return IndexPartition::getNIndices();
}

/** returns a vector with the numbers of constraints that are assigned to the classes */
std::vector<int> ConsPartition::getNConssOfClasses()
{
   return IndexPartition::getNIndicesOfClasses();
}


/** returns whether a constraint is already assigned to a class */
bool ConsPartition::isConsClassified(int givenConsindex)
{
   return IndexPartition::isIndexClassified(givenConsindex);
}

/** returns partition with reduced number of classes */
ConsPartition* ConsPartition::reduceClasses(int givenMaxNumber)
{
   std::vector<int> classindexmapping = IndexPartition::reduceClasses(givenMaxNumber);
   ConsPartition* newPartition;
   std::stringstream newName;
   std::stringstream newClassdesc;

   if( classindexmapping.empty() )
      return NULL;

   /* create new ConsPartition */
   newName << getName() << "-red-to-" << givenMaxNumber;
   newPartition = new ConsPartition(gcg, newName.str().c_str(), givenMaxNumber, getNConss());

   /* reassign conss */
   for( int i = 0; i < newPartition->getNConss(); ++i)
   {
      if( getClassOfCons(i) != -1 )
      {
         newPartition->assignConsToClass(i, classindexmapping[getClassOfCons(i)]);
      }
   }

   /* set new class names and descriptions (enlarged class has index 0) */
   newPartition->setClassName(0, "merged");
   newPartition->setClassDecompInfo(0, BOTH);

   for( int i = 0; i < getNClasses(); ++i )
   {
     if( classindexmapping[i] == 0 )
     {
        newClassdesc << getClassDescription( i ) << " - ";
     }
     else
     {
        newPartition->setClassName(classindexmapping[i], getClassName(i));
        newPartition->setClassDescription(classindexmapping[i], getClassDescription(i));
        newPartition->setClassDecompInfo(classindexmapping[i], getClassDecompInfo(i));
     }
   }

   newPartition->setClassDescription(0, newClassdesc.str().c_str());

   return newPartition;
}

/** sets the decomposition code of a class */
void ConsPartition::setClassDecompInfo(int givenClassindex, CONS_DECOMPINFO givenDecompInfo)
{
   assert(givenDecompInfo == BOTH || givenDecompInfo == ONLY_MASTER || givenDecompInfo == ONLY_PRICING );

   IndexPartition::setClassDecompInfo(givenClassindex, (int) givenDecompInfo );
}

} /* namespace gcg */
