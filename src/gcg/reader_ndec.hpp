/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file   reader_ndec.hpp
 * @brief  ndec file reader for (nested) structure information
 * @author Erik Muehmer
 * @ingroup FILEREADERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_READER_NDEC_HPP__
#define GCG_READER_NDEC_HPP__

#include <vector>
#include <unordered_map>
#include <string>
#include "class_partialdecomp.h"

namespace gcg
{

struct DecompositionData;

struct BlockData
{
   BlockData(int probnr) : decomposition(NULL), probnr(probnr), symmetricalblock(probnr) {}
   BlockData(const BlockData&) = delete;
   BlockData(BlockData&&) noexcept;
   ~BlockData();
   BlockData& operator=(const BlockData&) = delete;

   std::vector<std::string> constraints;
   DecompositionData* decomposition;
   int symmetricalblock;
   int probnr;
};

struct DecompositionData
{
   DecompositionData() = default;
   ~DecompositionData() = default;
   BLOCK_STRUCTURE* createBlockStructure(SCIP* scip, DETPROBDATA* detprobdata);

   std::vector<std::string> masterconstraints;
   std::vector<BlockData> blocks;
   std::unordered_map<std::string, std::string> symmetrydata;
};

struct NestedDecompositionData
{
   NestedDecompositionData() : version(0), presolved(false), rootdecomposition(NULL) {}
   NestedDecompositionData(const NestedDecompositionData&) = delete;
   ~NestedDecompositionData();
   NestedDecompositionData& operator=(const NestedDecompositionData&) = delete;

   int version;
   std::string name;
   bool presolved;
   std::string description;
   DecompositionData* rootdecomposition;
};

}

#endif //GCG_READER_NDEC_HPP__
