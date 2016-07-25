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

/**@file   class_seeed.cpp
 * @brief  class with functions for seeed (aka incomplete decomposition )
 * @author Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "class_seeed.h"
#include "gcg.h"
#include "class_seeedpool.h"

#include <iostream>
#include <exception>
#include <algorithm>

#define SCIP_CALL_EXC(x)   do                                                                                  \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) !=  SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                                                          \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

namespace gcg {

/** constructor(s) */
 Seeed::Seeed(
     SCIP*             scip,
	  int               givenId,      		   	/**< id that is given to this seeed */
	  int               givenNDetectors,         /**< number of detectors */
	  int				givenNConss,				/**number of constraints */
	  int 				givenNVars				/**number of variables */
    ): scip(scip), id(givenId), nBlocks(0),nVars(givenNVars), nConss(givenNConss), propagatedByDetector(std::vector<bool>(givenNDetectors, false)), openVarsAndConssCalculated(false), completedGreedily(false){
	 }

 Seeed::Seeed(const Seeed *seeedToCopy, Seeedpool* seeedpool)
 {
    scip = (seeedToCopy->scip);
    id = seeedpool->getNewIdForSeeed();
    nBlocks = seeedToCopy->nBlocks;
    nVars = seeedToCopy->nVars;
    nConss = seeedToCopy->nConss;
    masterConss = seeedToCopy->masterConss;
    masterVars = seeedToCopy->masterVars;
    conssForBlocks = seeedToCopy->conssForBlocks;
    varsForBlocks = seeedToCopy->varsForBlocks;
    linkingVars = seeedToCopy->linkingVars;
    stairlinkingVars = seeedToCopy->stairlinkingVars;
    openVars = seeedToCopy->openVars;
    openConss = seeedToCopy->openConss;
    propagatedByDetector = seeedToCopy->propagatedByDetector;
    detectorChain = seeedToCopy->detectorChain;
    openVarsAndConssCalculated = seeedToCopy->openVarsAndConssCalculated;
    completedGreedily = seeedToCopy->completedGreedily;

//    masterConss = std::vector<int>(0);
//    for( size_t i = 0; i < seeedToCopy->masterConss.size(); ++i)
//       masterConss.push_back(seeedToCopy->masterConss[i]);
//
//    masterVars = std::vector<int>(0);
//    for( size_t i = 0; i < seeedToCopy->masterVars.size(); ++i)
//       masterVars.push_back(seeedToCopy->masterVars[i]);
//
//    linkingVars = std::vector<int>(0);
//    for( size_t i = 0; i < seeedToCopy->linkingVars.size(); ++i)
//       linkingVars.push_back(seeedToCopy->linkingVars[i]);
//
//    openVars = std::vector<int>(0);
//    for( size_t i = 0; i < seeedToCopy->openVars.size(); ++i)
//       openVars.push_back(seeedToCopy->openVars[i]);
//
//    openConss = std::vector<int>(0);
//    for( size_t i = 0; i < seeedToCopy->openConss.size(); ++i)
//       openConss.push_back(seeedToCopy->openConss[i]);
//
//    detectorChain = std::vector<int>(0);
//    for( size_t i = 0; i < seeedToCopy->detectorChain.size(); ++i)
//       detectorChain.push_back(seeedToCopy->detectorChain[i]);
//
//    propagatedByDetector = std::vector<bool>(0);
//    for( size_t i = 0; i < seeedToCopy->propagatedByDetector.size(); ++i)
//       propagatedByDetector.push_back(seeedToCopy->propagatedByDetector[i]);
//
//    stairlinkingVars = std::vector<std::vector<int>>(0);
//    std::vector<int> vec = std::vector<int>(0);
//    for( int b = 0; b < nBlocks; ++b)
//       stairlinkingVars.push_back(vec);
//    for( size_t i; i < seeedToCopy->stairlinkingVars.size(); ++i )
//    {
//       for( size_t j; j < seeedToCopy->stairlinkingVars[i].size(); ++j)
//          stairlinkingVars[i].push_back(seeedToCopy->stairlinkingVars[i][j]);
//    }
//
//    conssForBlocks = std::vector<std::vector<int>>(0);
//    //std::vector<int> vec = std::vector<int>(0);
//    for( int b = 0; b < nBlocks; ++b)
//       conssForBlocks.push_back(vec);
//    for( size_t i; i < seeedToCopy->conssForBlocks.size(); ++i )
//    {
//       for( size_t j; j < seeedToCopy->conssForBlocks[i].size(); ++j)
//          conssForBlocks[i].push_back(seeedToCopy->conssForBlocks[i][j]);
//    }
//
//    varsForBlocks = std::vector<std::vector<int>>(0);
//    //std::vector<int> vec = std::vector<int>(0);
//    for( int b = 0; b < nBlocks; ++b)
//       varsForBlocks.push_back(vec);
//    for( size_t i; i < seeedToCopy->varsForBlocks.size(); ++i )
//    {
//       for( size_t j; j < seeedToCopy->varsForBlocks[i].size(); ++j)
//          varsForBlocks[i].push_back(seeedToCopy->varsForBlocks[i][j]);
//    }
 }

 Seeed::~Seeed(){}

/** check the consistency of this seeed */
  bool Seeed::checkConsistency(
  ){

	 // @TODO: check mastervars, stairlinking

	  /**check variables (every variable is assigned at most once) */

	  std::vector<bool> openVarsBool(nVars, true) ;
	  std::vector<int>  openVarsVec(0);
	  std::vector<int>::const_iterator varIter = linkingVars.begin();
	  std::vector<int>::const_iterator varIterEnd = linkingVars.end();
	  for(; varIter != varIterEnd; ++varIter)
	  {
		  if(!openVarsBool[*varIter])
		  {
			  std::cout << "Warning! Variable with index " << *varIter << "is already assigned " << std::endl;
			  return false;
		  }
		  openVarsBool[*varIter] = false;
	  }
	  for(int b =0; b < nBlocks; ++b)
	  {
		  varIter = varsForBlocks[b].begin();
		  varIterEnd = varsForBlocks[b].end();
		  for(; varIter != varIterEnd; ++varIter)
		  {
			  if(!openVarsBool[*varIter])
			  {
				  std::cout << "Warning! Variable with index " << *varIter << "is already assigned " << std::endl;
				  return false;
			  }
			  openVarsBool[*varIter] = false;
		  }
	  }

	  /** check constraints (every constraint is assigned at most once */
	  std::vector<bool> openConssBool(nConss, true) ;
	  std::vector<int>  openConssVec(0);
	  std::vector<int>::const_iterator consIter = masterConss.begin();
	  std::vector<int>::const_iterator consIterEnd = masterConss.end();
	  for(; consIter != consIterEnd; ++consIter)
	  {
		  if(!openConssBool[*consIter])
		  {
			  std::cout << "Warning! Constraint with index " << *consIter << "is already assigned " << std::endl;
			  return false;
		  }
		  openConssBool[*consIter] = false;
	  }

	  for(int b =0; b < nBlocks; ++b)
	  {
		  consIter = conssForBlocks[b].begin();
		  consIterEnd = conssForBlocks[b].end();
		  for(; consIter != consIterEnd; ++consIter)
		  {
			  if(!openConssBool[*consIter])
			  {
				  std::cout << "Warning! Constraint with index " << *consIter << "is already assigned " << std::endl;
				  return false;
			  }
			  openConssBool[*consIter] = false;
		  }
	  }

	  return true;
  }


  /** set-methods */

  /** set number of blocks, atm only increasing number of blocks  */
  SCIP_RETCODE Seeed::setNBlocks(int newNBlocks
  ){
	  assert(newNBlocks >= nBlocks);

	  /** increase number of blocks in conssForBlocks and varsForBlocks */
	  for(int b = nBlocks; b < newNBlocks; ++b )
	  {
		  conssForBlocks.push_back(std::vector<int>(0) );
		  varsForBlocks.push_back(std::vector<int>(0) );
	  }

	  nBlocks = newNBlocks;

	  return SCIP_OKAY;
  }

  /** add a constraint to the master constraints */
  SCIP_RETCODE Seeed::setConsToMaster(
		   int consToMaster
  ){
	  masterConss.push_back(consToMaster);

	  return SCIP_OKAY;
  }

  /** add a variable to the master variables (every constraint consisting it is in master ) */
  SCIP_RETCODE Seeed::setVarToMaster(
		  int varToMaster )
  {
  	  masterVars.push_back(varToMaster);

  	  return SCIP_OKAY;
  }



  /** add a constraint to a block */
  SCIP_RETCODE Seeed::setConsToBlock(
		   int consToBlock,
		   int block
  ){
	  assert(conssForBlocks.size() > block);

	  conssForBlocks[block].push_back(consToBlock);

	  return SCIP_OKAY;
  }

  /** add a variable to a block */
  SCIP_RETCODE Seeed::setVarToBlock(
		   int varToBlock,
		   int block
  )
  {
	  assert(varsForBlocks.size() > block);

	  varsForBlocks[block].push_back(varToBlock);
	  return SCIP_OKAY;
  }

  /** add a variable to the linking variables */
  SCIP_RETCODE Seeed::setVarToLinking(
		   int varToLinking
  ){
	  linkingVars.push_back(varToLinking);
	  return SCIP_OKAY;
  }

  /** add a variable to the stair linking variables */
  SCIP_RETCODE Seeed::setVarToStairlinking(
		   int varToStairLinking, int block1, int block2
  ){
	  stairlinkingVars[block1].push_back(varToStairLinking);
	  stairlinkingVars[block2].push_back(varToStairLinking);

	  return SCIP_OKAY;
  }

  /** finds linking-variables that are actually master-variables. I.e. the variable is adjacent to only master-constraints. */
  SCIP_RETCODE Seeed::findVarsLinkingToMaster(Seeedpool* seeedpool)
  {
    int i;
    int j;
    int* varcons;
    bool isMasterVar;
    int* lvars = getLinkingvars();
    std::vector<int> foundMasterVarIndices;
    
    // sort Master constraints for binary search
    sortMasterconss();

    for (i = 0; i < getNLinkingvars(); ++i)
    {
      isMasterVar = true;
      varcons = seeedpool->getConssForVar(lvars[i]);
      for(j = 0; j < seeedpool->getNConssForVar(lvars[i]); ++j)
      {
        if (!std::binary_search(masterConss.begin(), masterConss.end(), varcons[j])
        {
          isMasterVar = false;
          break;
        }
      }

      if (isMasterVar)
      {
        foundMasterVarIndices.push_back(i);
      }
    }

    for (std::vector<int>::iterator it = foundMasterVarIndices.rbegin(); it != foundMasterVarIndices.rend(); ++it)
    {
      masterVars.push_back(lvars[*it]);
      linkingVars.erase(*it);
    }

    return SCIP_OKAY;
  }
  
  /** finds linking-variables that are actually stairlinking-ariables. I.e. the variable is adjacent to constraints in exactly two block (not adjacent to open cons and master-cons). */
  SCIP_RETCODE Seeed::findVarsLinkingToStairlinking(Seeedpool* seedpool)
  {
    int i;
    int j;
    int k;
   
    int consblock;
    int block1 = -1;
    int block2 = -1;

    int* varcons;
    int* lvars = getLinkingvars();

    std::vector<int> foundMasterVarIndices;

    for(i = 0; i < getNBlocks(); ++i)
    {
      sortConssForBlock(i);
    }

    for(i = 0; i < getNLinkingvars(); ++i)
    {
      block1 = -1; block2 = -1;
      varcons = seeedpool->getConssForVar(lvars[i]); 
      for (j = 0; j < seeedpool->getNConssForVar(lvars[i]); ++j)
      {
        consblock = -1;
        for (k = 0; k < nBlocks; ++k)
        {
          if (std::binary_search(conssForBlocks[k].begin(), conssForBlocks[k].end(), varcons[j]))
          {
            consblock = k;
            break;
          }
        }

        if (consblock == -1)
        {
          block1 = -1;
          block2 = -1;
          break;
        }
        else if (block1 == consblock || block2 == consblock)
        {
          continue;
        }
        else if (block1 == -1)
        {
          block1 = consblock;
          continue;
        }
        else if (block2 == -1)
        {
          block2 = consblock;
          continue;
        }
        else
        {
          block1 = -1;
          block2 = -1;
          break;
        }
      }

      if (block1 != -1 && block2 != -1)
      {
        setVarToStairlinking(lvars[i], block1, block2);
        foundMasterVarIndices.push_back(i);
      }
    }

    for (std::vector<int>::iterator it = foundMasterVarIndices.rbegin(); it != foundMasterVarIndices.rend(); ++it)
    {
      linkingVars.erase(*it);
    }
  }

  /** sets seeed to be propagated by detector with detectorID  */
  SCIP_RETCODE Seeed::setDetectorPropagated(
		   int detectorID
  ){
	  assert(propagatedByDetector.size() > detectorID );
	  propagatedByDetector[detectorID]  = true;

	  return SCIP_OKAY;
  }


  /** get-methods */

  /** returns vector containing master conss */
  const int* Seeed::getMasterconss(
  ){
	  return &masterConss[0];
  }

  /** returns vector containing master conss */
  int Seeed::getNMasterconss(
    ){
  	  return (int) masterConss.size();
    }

  /** sorts master conss */
  void Seeed::sortMasterconss(
    ){
     std::sort(masterConss.begin(), masterConss.end());
     return;
    }

  /** returns vector containing master vars (every constraint containing a master var is in master)*/
  const int* Seeed::getMastervars(
  ){
	  return &masterVars[0];
  }

  /** returns vector containing master vars (every constraint containing a master var is in master)*/
  int Seeed::getNMastervars(
  ){
	  return(int) masterVars.size();
  }

  /** sorts master vars */
  void Seeed::sortMastervars(
  ){
     std::sort(masterVars.begin(), masterVars.end());
     return;
  }

  /** returns number of blocks */
  int Seeed::getNBlocks(
  ){
     return nBlocks;
  }

  /** returns vector containing master conss */
  const int* Seeed::getConssForBlock(
		   int block
  ){
	  return &conssForBlocks[block][0];
  }

  /** returns vector containing master conss */
  int Seeed::getNConssForBlock(
		  int block
  ){
	  return (int)conssForBlocks[block].size();
  }

  /** sorts cons of a certain block */
  void Seeed::sortConssForBlock(
        int block
  ){
     std::sort(conssForBlocks[block].begin(), conssForBlocks[block].end());
     return;
  }

  /** returns vector containing vars of a certain block */
  const int* Seeed::getVarsForBlock(
		   int block
  ){
	  return &varsForBlocks[block][0];
  }

  /** returns vector containing vars of a certain block */
   int Seeed::getNVarsForBlock(
  		   int block
    ){
  	  return (int)varsForBlocks[block].size();
    }

   /** sorts cons of a certain block */
   void Seeed::sortVarsForBlock(
         int block
   ){
      std::sort(varsForBlocks[block].begin(), varsForBlocks[block].end());
      return;
   }

  /** returns vector containing linking vars */
  const int* Seeed::getLinkingvars(
  ){
	  return &linkingVars[0];
  }

  /** returns size of vector containing linking vars */
    int Seeed::getNLinkingvars(
    ){
  	  return (int)linkingVars.size();
    }

    /** sorts linking vars */
    void Seeed::sortLinkingvars(
    ){
       std::sort(linkingVars.begin(), linkingVars.end());
       return;
    }

  /** returns vector containing stairlinking vars */
  const int* Seeed::getStairlinkingvars(
       int block
  ){
	  return &stairlinkingVars[block][0];
  }

  /** returns vector containing stairlinking vars */
  int Seeed::getNStairlinkingvars(
     int block
    ){
  	  return (int) stairlinkingVars[block].size();
    }

  /** sorts stairlinking vars */
  void Seeed::sortStairlinkingvars(
     int block
  ){
     std::sort(stairlinkingVars[block].begin(), stairlinkingVars[block].end());
     return;
  }

  /** returns vector containing variables not assigned yet*/
  const int* Seeed::getOpenvars(
    ){
  	  if(!openVarsAndConssCalculated)
  	  {
  		  calcOpenconss();
  		  calcOpenvars();

  		  openVarsAndConssCalculated = true;
  	  }

  	  return &openVars[0];
    }

  /** returns vector containing constraints not assigned yet */
  const int* Seeed::getOpenconss(
  ){
	  if(!openVarsAndConssCalculated)
	  {
		  calcOpenconss();
		  calcOpenvars();

		  openVarsAndConssCalculated = true;
	  }

	  return &openConss[0];
  }

  /** returns size of vector containing variables not assigned yet */
  int Seeed::getNOpenconss()
  {
	  if(!openVarsAndConssCalculated)
	  {
		  calcOpenconss();
		  calcOpenvars();

		  openVarsAndConssCalculated = true;
	  }
	  return (int) openConss.size();
}

      /** returns size of vector containing constraints not assigned yet */
  int Seeed::getNOpenvars()
   {
 	  if(!openVarsAndConssCalculated)
 	  {
 		  calcOpenconss();
 		  calcOpenvars();

 		  openVarsAndConssCalculated = true;
 	  }
 	  return (int) openVars.size();
 }


  /** constructs vector containing variables not assigned yet */
  void Seeed::calcOpenvars(
  ){

	  openVars = std::vector<int>(0);
	  std::vector<bool> openVarsBool(nVars, true) ;
	  std::vector<int>::const_iterator varIter = linkingVars.begin();
	  std::vector<int>::const_iterator varIterEnd = linkingVars.end();
	  for(; varIter != varIterEnd; ++varIter)
		  openVarsBool[*varIter] = false;
	  for(int b =0; b < nBlocks; ++b)
	  {
		  varIter = varsForBlocks[b].begin();
		  varIterEnd = varsForBlocks[b].end();
		  for(; varIter != varIterEnd; ++varIter)
		  		  openVarsBool[*varIter] = false;
	  }

	  for (int i = 0; i < nVars; ++i)
	  {
		  if(openVarsBool[i])
			  openVars.push_back(i);
	  }


	  return;

  }

  /** calculates vector containing constraints not assigned yet */
  void  Seeed::calcOpenconss(
  ){
	  std::vector<bool> openConssBool(nConss, true) ;
	  openConss = std::vector<int>(0);
	  std::vector<int>::const_iterator consIter = masterConss.begin();
	  std::vector<int>::const_iterator consIterEnd = masterConss.end();
	  for(; consIter != consIterEnd; ++consIter)
	  		  openConssBool[*consIter] = false;
	  for(int b =0; b < nBlocks; ++b)
	  {
		  consIter = conssForBlocks[b].begin();
		  consIterEnd = conssForBlocks[b].end();
		  for(; consIter != consIterEnd; ++consIter)
			  openConssBool[*consIter] = false;
	  }

	  for (int i = 0; i < nConss; ++i)
	  {
		  if(openConssBool[i])
			  openConss.push_back(i);
	  }


	  return;
  }

  /** returns whether this seeed was propagated by certain detector */
  bool Seeed::isPropagatedBy(
		   int detectorID
  ){
	  assert(propagatedByDetector.size() > detectorID);
	  return propagatedByDetector[detectorID];
  }

  /** assigns the open cons and open vars */
  SCIP_RETCODE Seeed::completeGreedily(Seeedpool* seeedpool){

     std::vector<int> assignedOpenvars; /** stores the assigned open vars */

     /** tools to update openVars */
     std::vector<int> oldOpenvars;
     bool found;
     bool notassigned;

	  if(!openVarsAndConssCalculated)
	  {
		  calcOpenconss();
		  calcOpenvars();

		  openVarsAndConssCalculated = true;
	  }

	  if(nBlocks == 0)
	  {
	     nBlocks = 1;
	     std::vector<int> vec = std::vector<int>(0);
	     conssForBlocks.push_back(vec);
	     varsForBlocks.push_back(vec);
	     stairlinkingVars.push_back(vec);
	     if(nBlocks == 0 && openConss.size()!=0)
	     {
	        setConsToBlock(openConss[0], 0);
	        openConss.erase(openConss.begin());
	     }
	     else if(nBlocks == 0 && openConss.size() == 0 && masterConss.size() !=0)
	     {
	        setConsToBlock(masterConss[0], 0);
	        masterConss.erase(masterConss.begin());
	     }
	     else if(openConss.size() == 0 && masterConss.size() ==0)
	     {
	        assert(false);
	     }

	  }



	  /** check if the openVars can found in a constraint yet */
	  for( size_t i = 0; i < openVars.size(); ++i )
	  {
	     bool checkVar = true; /** if the var has to be checked any more */
	     bool varInBlock; /** if the var can be found in the block */
	     std::vector<int> varInBlocks; /** stores, in which block the variable can be found */

	     /** test if the variable can be found in blocks */
	     for( int b = 0; b < nBlocks; ++b )
	     {
	        varInBlock = false;
	        for( size_t k = 0; k < conssForBlocks[b].size() && !varInBlock; ++k)
	        {
	           for( int l = 0; l < seeedpool->getNVarsForCons(conssForBlocks[b][k]); ++l )
	           {
	              if( openVars[i] == *(seeedpool->getVarsForCons(conssForBlocks[b][k])) )
	              {
	                 varInBlocks.push_back(b);
	                 varInBlock = true;
	                 break;
	              }
	           }
	        }
	     }
	     if( varInBlocks.size() == 1) /** if the variable can be found in one block set the variable to a variable of the block*/
	     {
	        setVarToBlock(openVars[i], varInBlocks[0]);
	        assignedOpenvars.push_back(openVars[i]);
	        continue; /** the variable openVars[i] does'nt need to be checked any more */
	     }
	     else if( varInBlocks.size() == 2) /** if the variable can be found in two blocks check if it is a linking var or a stairlinking var*/
	     {
	        if ( varInBlocks[0] + 1 == varInBlocks[1] )
	        {
	           setVarToStairlinking(openVars[i], varInBlocks[0], varInBlocks[1]);
	           setVarToBlock(openVars[i], varInBlocks[0]);
	           setVarToBlock(openVars[i], varInBlocks[1]);
	           assignedOpenvars.push_back(openVars[i]);
	           continue; /** the variable openVars[i] does'nt need to be checked any more */
	        }
	        else
	        {
	           setVarToLinking(openVars[i]);
	           setVarToBlock(openVars[i], varInBlocks[0]);
	           setVarToBlock(openVars[i], varInBlocks[1]);
	           assignedOpenvars.push_back(openVars[i]);
	           continue; /** the variable openVars[i] does'nt need to be checked any more */
	        }
	     }
	     else if( varInBlocks.size() > 2 ) /** if the variable can be found in more than two blocks it is a linking var */
	     {
	        setVarToLinking(openVars[i]);
	        for( size_t k = 0; k < varInBlocks.size(); k++ )
	        {
	           setVarToBlock(openVars[i], varInBlocks[k]);
	        }
	        assignedOpenvars.push_back(openVars[i]);
	        continue; /** the variable openVars[i] does'nt need to be checked any more */
	     }

	     /** if the variable can be found in an open constraint it is still an open var */
	     for( size_t j = 0; j < openConss.size(); ++j )
	     {
	        for( int k = 0; k < seeedpool->getNVarsForCons(j); ++k)
	        {
	           if( openVars[i] == seeedpool->getVarsForCons(j)[k])
	           {
	              checkVar = false;
	              break;
	           }
	        }
	        if ( !checkVar )
	        {
	           break;
	        }
	     }

	     /** test if the variable can be found in a master constraint yet */
	     for( size_t j = 0; j < masterConss.size() && checkVar; ++j )
	     {
	        for ( int k = 0; k < seeedpool->getNVarsForCons(masterConss[j]); ++k )
	        {
	           if( openVars[i] == seeedpool->getVarsForCons(masterConss[j])[k] )
	           {
	              setVarToMaster(openVars[i]);
	              assignedOpenvars.push_back(openVars[i]);
	              checkVar = false; /** the variable openVars[i] does'nt need to be checked any more */
	              break;
	           }
	        }
	     }

	  }
	  /** delete the assigned open vars */
	  oldOpenvars = openVars;
	  openVars.clear();
	  for ( size_t i = 0; i < oldOpenvars.size(); ++i)
	  {
	     found = false;
	     for ( size_t j = 0; j < assignedOpenvars.size(); ++j )
	     {
	        if( oldOpenvars[i] == assignedOpenvars[j] )
	        {
	           found = true;
	           break;
	        }
	     }
	     if( !found ) /** var is still open var */
	     {
	        openVars.push_back(oldOpenvars[i]);
	     }
	  }
	  oldOpenvars.clear();
	  assignedOpenvars.clear();


	  /** assign open conss greedily */
	  for( size_t i = 0; i < openConss.size(); ++i)
	  {
	     std::vector<int> vecOpenvarsOfBlock; /** stores the open vars of the blocks */
	     bool consGotBlockcons = false; /** if the constraint can be assigned to a block */

	     /** check if the constraint can be assigned to a block */
	     for ( int j = 0; j < nBlocks; ++j )
	     {
	        /** check if all vars of the constraint are a block var of the current block, an open var, a linkingvar or a mastervar*/
	        consGotBlockcons = true;
	        for( int k = 0; k < seeedpool->getNVarsForCons(i) ; ++k )
	        {
	           if ( isVarBlockvarOfBlock(seeedpool->getVarsForCons(openConss[i])[k], j) || isVarMastervar(seeedpool->getVarsForCons(openConss[i])[k]) ||
	              isVarOpenvar(seeedpool->getVarsForCons(openConss[i])[k]) )
	           {
	              if ( isVarOpenvar(seeedpool->getVarsForCons(openConss[i])[k]) )
	              {
	                 vecOpenvarsOfBlock.push_back(seeedpool->getVarsForCons(openConss[i])[k]); /**!!!*/
	              }
	           }
	           else
	           {
	              vecOpenvarsOfBlock.clear(); /** the open vars do'nt get vars of the block */
	              consGotBlockcons = false; /** the constraint can't be constraint of the block, check the next block */
	              break;
	           }
	        }
	        if ( consGotBlockcons ) /** the constraint can be assigned to the current block */
	        {
	           setConsToBlock(openConss[i], j);
	           for( size_t k = 0; k < vecOpenvarsOfBlock.size(); ++k) /** the openvars in the constraint get block vars */
	           {
	              setVarToBlock(vecOpenvarsOfBlock[k], j);
	              assignedOpenvars.push_back(vecOpenvarsOfBlock[k]);
	           }
	           vecOpenvarsOfBlock.clear();

	           /** delete the assigned open vars */
	           oldOpenvars = openVars;
	           openVars.clear();
	           for ( size_t l = 0; l < oldOpenvars.size(); ++l)
	           {
	              found = false;
	              for ( size_t m = 0; m < assignedOpenvars.size(); ++m )
	              {
	                 if( oldOpenvars[l] == assignedOpenvars[m] )
	                 {
	                    found = true;
	                    break;
	                 }
	              }
	              if( !found ) /** var is still open var */
	              {
	                 openVars.push_back(oldOpenvars[l]);
	              }
	           }
	           oldOpenvars.clear();
	           assignedOpenvars.clear();

	           break;
	        }
	     }


	     if( !consGotBlockcons ) /** the constraint can not be assigned to a block, set it to master cons */
	     {
	        setConsToMaster(openConss[i]);
	     }
	  }

	  for(size_t i = 0; i < openVars.size(); i++)
	  {
	     notassigned = true;
	     //setVarToMaster(openVars[i]);
	     for(size_t j = 0; j < masterConss.size() && notassigned; ++j)
	     {
	        for(int k = 0; k < seeedpool->getNVarsForCons(masterConss[j]) && notassigned ; ++k)
	        {
	           if(openVars[i] == seeedpool->getVarsForCons(masterConss[j])[k])
	           {
	              setVarToMaster(openVars[i]);
	              assignedOpenvars.push_back(openVars[i]);
	              notassigned = false;
	           }
	        }
	     }
	  }
	  /** delete the assigned open vars */
	  oldOpenvars = openVars;
	  openVars.clear();
	  for ( size_t i = 0; i < oldOpenvars.size(); ++i)
	  {
	     found = false;
	     for ( size_t j = 0; j < assignedOpenvars.size(); ++j )
	     {
	        if( oldOpenvars[i] == assignedOpenvars[j] )
	        {
	           found = true;
	           break;
	        }
	     }
	     if( !found ) /** var is still open var */
	     {
	        openVars.push_back(oldOpenvars[i]);
	     }
	  }
	  oldOpenvars.clear();
	  assignedOpenvars.clear();

	  /** check if the open cons are all assigned */
	  if( ! checkAllConsAssigned() )
	  {
	     std::cout << "ERROR: Something went wrong, there are still open cons, although all should have been assigned ";
	     assert(false);
	  }

	  /** check if the open vars are all assigned */
	  if( ! openVars.empty() )
	  {
	     std::cout << "ERROR: Something went wrong, there are still open vars, although all should have been assigned ";
	     assert(false);
	  }

	  completedGreedily = true;
	  return SCIP_OKAY;

  }

  /** returns whether the var is a linking var */
  bool Seeed::isVarLinkingvar(int var){
     for( size_t i = 0;  i < linkingVars.size(); ++i)
     {
        if( var == linkingVars[i])
        {
           return true;
        }
     }
     return false;
  }

  /** return whether the var is a var of the block */
  bool Seeed::isVarBlockvarOfBlock(int var, int block){
     for( size_t i = 0;  i < varsForBlocks[block].size(); ++i)
     {
        if( var == varsForBlocks[block][i])
        {
           return true;
        }
     }
     return false;
  }

  /** returns whether the var is a master var */
  bool Seeed::isVarMastervar(int var){
     for( size_t i = 0;  i < masterVars.size(); ++i)
     {
        if( var == masterVars[i])
        {
           return true;
        }
     }
     return false;
  }

  /** returns whether the var is an open var */
  bool Seeed::isVarOpenvar(int var){
     for( size_t i = 0;  i < openVars.size(); ++i)
     {
        if( var == openVars[i])
        {
           return true;
        }
     }
     return false;
  }

  /** returns indes of the openvar in the vector */
  int Seeed::getIndexOfOpenvar(int var){
     for( size_t i = 0;  i < openVars.size(); ++i)
     {
        if( var == openVars[i])
        {
           return i;
        }
     }
     return -1;
  }

 /** returns whether all cons are assigned and deletes the vector open cons if all are assigned */
bool Seeed::checkAllConsAssigned(){
   for( size_t i = 0; i < openConss.size(); ++i)
   {
      bool consfound = false;
      for( size_t k = 0; k < masterConss.size(); ++k)
      {
         if( openConss[i] == masterConss[k] )
         {
            consfound = true;
            break;
         }
      }
      for( int b = 0; b < nBlocks && !consfound; ++b )
      {
         for( size_t k = 0; k < conssForBlocks[b].size(); ++k)
         {
            if( openConss[i] == conssForBlocks[b][k])
            {
               consfound = true;
               break;
            }
         }
      }
      if(!consfound)
      {
         return false;
      }
   }
   openConss.clear();
   return true;
}

int* Seeed::getDetectorchain()
{
   return &detectorChain[0];
}

int Seeed::getNDetectors()
{
   return (int) detectorChain.size();
}

int Seeed::getID()
{
   return id;
}

int Seeed::getNConss()
{
  return nConss;
}

int Seeed::getNVars()
{
   return nVars;
}

bool Seeed::isSeeedCompletedGreedily()
{
   return completedGreedily;
}

} /* namespace gcg */
