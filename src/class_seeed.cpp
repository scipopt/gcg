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

/**@file   class_colpool.cpp
 * @brief  class with functions for colpool
 * @author Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "class_seeed.h"
#include "gcg.h"


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
	  int               id,      		   	/**< id that is given to this seeed */
	  int               nDetectors,         /**< number of detectors */
	  int				nConss,				/**number of constraints */
	  int 				nVars				/**number of variables */
    ): id(id), nBlocks(0), propagatedByDetector(std::vector<bool>(nDetectors, false)){

	 for(int i = 0; i < nConss; ++i)
	 {
		 openConss.push_back(i);
	 }

	 for(int j = 0; j < nVars; ++j)
	 {
		 openVars.push_back(j);
	 }


 }

 Seeed::~Seeed(){}

/** check the consistency of this seeed */
  bool Seeed::checkConsistency(
  ){
	  return false;
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
	  masterconss.push_back(consToMaster);

	  return SCIP_OKAY;
  }

  /** add a variable to the master variables (every constraint consisting it is in master ) */
  SCIP_RETCODE Seeed::setVarToMaster(
		   int varToMaster
  );

  /** add a constraint to a block */
  SCIP_RETCODE Seeed::setConsToBlock(
		   int consToBlock,
		   int block
  );

  /** add a variable to a block */
  SCIP_RETCODE Seeed::setVarToBlock(
		   int varToBlock,
		   int block
  );

  /** add a variable to the linking variables */
  SCIP_RETCODE Seeed::setVarToLinking(
		   int varToLinking
  );

  /** add a variable to the stair linking variables */
  SCIP_RETCODE Seeed::setVarToStairlinking(
		   int varToStairLinking
  );


  /** get-methods */

  /** returns vector containing master conss */
  std::vector<int> const & Seeed::getMasterconss(
  );


  /** returns vector containing master vars (every constraint containing a master var is in master )*/
  std::vector<int> const & Seeed::getMastervars(
  );

  /** returns vector containing master conss */
  std::vector<int> const & Seeed::getConssForBlock(
		   int block
  );

  /** returns vector containing vars of a certain block */
  std::vector<int> const & Seeed::getVarsForBlock(
		   int block
  );

  /** returns vector containing linking vars */
  std::vector<int> const & Seeed::getLinkingvars(
  );

  /** returns vector containing stairlinking vars */
  std::vector<int> const & Seeed::getStairlinkingvars(
  );

  /** returns vector containing variables not assigned yet */
  std::vector<int> const & Seeed::getOpenvars(
  );

  /** returns vector containing constraints not assigned yet */
  std::vector<int> const & Seeed::getOpenconss(
  );

  /** returns whether this seeed was propagated by certain detector */
  bool Seeed::isPropagatedBy(
		   int detectorID
  );



} /* namespace gcg */
