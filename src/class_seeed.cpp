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

#include <iostream>
#include <exception>

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
	  int               givenId,      		   	/**< id that is given to this seeed */
	  int               givenNDetectors,         /**< number of detectors */
	  int				givenNConss,				/**number of constraints */
	  int 				givenNVars				/**number of variables */
    ): id(givenId), nBlocks(0),nVars(givenNVars), nConss(givenNConss), propagatedByDetector(std::vector<bool>(givenNDetectors, false)), openVarsAndConssCalculated(false){

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

	  conssForBlocks[consToBlock].push_back(consToBlock);

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
		   int varToStairLinking
  ){
	  stairlinkingVars.push_back(varToStairLinking);

	  return SCIP_OKAY;
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


  /** returns vector containing master vars (every constraint containing a master var is in master )*/
  const int* Seeed::getMastervars(
  ){
	  return &masterVars[0];
  }

  /** returns vector containing master vars (every constraint containing a master var is in master )*/
  int Seeed::getNMastervars(
  ){
	  return(int) masterVars.size();
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


  /** returns vector containing stairlinking vars */
  const int* Seeed::getStairlinkingvars(
  ){
	  return &stairlinkingVars[0];
  }

  /** returns vector containing stairlinking vars */
  int Seeed::getNStairlinkingvars(
    ){
  	  return (int) stairlinkingVars.size();
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
		  if(!openVarsBool[i])
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
		  if(!openConssBool[i])
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

  SCIP_RETCODE Seeed::completeGreedily(){

	  if(!openVarsAndConssCalculated)
	  {
		  calcOpenconss();
		  calcOpenvars();

		  openVarsAndConssCalculated = true;
	  }



	  return SCIP_OKAY;

  }



} /* namespace gcg */
