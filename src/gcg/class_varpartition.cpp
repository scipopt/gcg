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

/**@file   class_varpartition.cpp
 * @brief  class representing a partition of a set of variables
 * @author Julius Hense
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/class_varpartition.h"

#include <cassert>
#include <sstream>
#include <algorithm>


namespace gcg {


/** constructor */
VarPartition::VarPartition(
   GCG*           gcgstruct,
   const char*    givenName,
   int            givenNClasses,
   int            givenNVars
) :
   IndexPartition(gcgstruct, givenName, givenNClasses, givenNVars)
{

}

/** copy constructor */
VarPartition::VarPartition(
   const VarPartition* toCopy
) : IndexPartition(toCopy)
{
}

/** destructor */
VarPartition::~VarPartition()
{
}

/** creates a new class, returns index of the class */
int VarPartition::addClass(const char* givenName, const char* givenDesc, VAR_DECOMPINFO givenDecompInfo)
{
   int classindex = IndexPartition::addClass(givenName, givenDesc);
   setClassDecompInfo(classindex, givenDecompInfo);

   return classindex;
}

/** assigns a variable to a class */
void VarPartition::assignVarToClass(int givenVarindex, int givenClassindex)
{
   IndexPartition::assignIndexToClass(givenVarindex, givenClassindex);
}

/** returns a vector containing all possible subsets of the chosen classindices */
std::vector<std::vector<int>> VarPartition::getAllSubsets(bool all, bool linking, bool master, bool block)
{
   std::vector<int> classindices;
   for ( int i = 0; i < getNClasses(); ++i )
   {
      if( ( all && getClassDecompInfo(i) == ALL ) ||
           ( linking && getClassDecompInfo(i) == LINKING ) ||
           ( master && getClassDecompInfo(i) == MASTER ) ||
           ( block && getClassDecompInfo(i) == BLOCK )
          )
         classindices.push_back(i);
   }
   return IndexPartition::getAllSubsets(classindices);
}

/** returns the decomposition code of a class */
VAR_DECOMPINFO VarPartition::getClassDecompInfo(int givenClassindex)
{
   int decompInfo = IndexPartition::getClassDecompInfo(givenClassindex);
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
const char* VarPartition::getClassNameOfVar(int givenVarindex)
{
   return IndexPartition::getClassNameOfIndex(givenVarindex );
}

/** returns the index of the class a variable is assigned to */
int VarPartition::getClassOfVar(int givenVarindex)
{
   return IndexPartition::getClassOfIndex(givenVarindex);
}

/** returns vector containing the assigned class of each variable */
const int* VarPartition::getVarsToClasses()
{
   std::vector<int>& varsToClasses = IndexPartition::getIndicesToClasses();
   if( !varsToClasses.empty() )
      return &varsToClasses[0];
   else
      return NULL;
}

/** returns the number of variables */
int VarPartition::getNVars()
{
   return IndexPartition::getNIndices();
}

/** returns a vector with the numbers of variables that are assigned to the classes */
std::vector<int> VarPartition::getNVarsOfClasses()
{
   return IndexPartition::getNIndicesOfClasses();
}


/** returns whether a variable is already assigned to a class */
bool VarPartition::isVarClassified(int givenVarindex)
{
   return IndexPartition::isIndexClassified(givenVarindex);
}

/** returns partition with reduced number of classes */
VarPartition* VarPartition::reduceClasses(int givenMaxNumber)
{
   std::vector<int> classindexmapping = IndexPartition::reduceClasses(givenMaxNumber);
   VarPartition* newPartition;
   std::stringstream newName;
   std::stringstream newClassdesc;

   if( classindexmapping.empty() )
      return NULL;

   /* create new VarPartition */
   newName << getName() << "-red-to-" << givenMaxNumber;
   newPartition = new VarPartition(gcg, newName.str().c_str(), givenMaxNumber, getNVars());

   /* reassign vars */
   for( int i = 0; i < newPartition->getNVars(); ++i)
   {
      if( getClassOfVar(i) != -1 )
      {
         newPartition->assignVarToClass(i, classindexmapping[getClassOfVar(i)]);
      }
   }

   /* set new class names and descriptions (enlarged class has index 0) */
   newPartition->setClassName(0, "merged");
   newPartition->setClassDecompInfo(0, ALL);

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
void VarPartition::setClassDecompInfo(int givenClassindex, VAR_DECOMPINFO givenDecompInfo)
{
   assert( givenDecompInfo == ALL || givenDecompInfo == LINKING || givenDecompInfo == MASTER || givenDecompInfo == BLOCK );

   IndexPartition::setClassDecompInfo(givenClassindex, (int) givenDecompInfo);
}

} /* namespace gcg */
