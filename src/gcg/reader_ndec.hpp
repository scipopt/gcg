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

#include <jansson.h>
#include <vector>
#include <unordered_map>
#include <string>
#include "class_partialdecomp.h"

namespace gcg
{

struct DecompositionData;
struct NDecFileHandler;

struct BlockData
{
   BlockData() : decomposition(NULL), symmetricalblock(-1) {}
   BlockData(const BlockData&) = delete;
   BlockData(BlockData&&) noexcept;
   ~BlockData();
   BlockData& operator=(const BlockData&) = delete;

   std::vector<std::string> constraints;
   DecompositionData* decomposition;
   int symmetricalblock;
};

struct DecompositionData
{
   DecompositionData() = default;
   ~DecompositionData() = default;
   BLOCK_STRUCTURE* createBlockStructure(DETPROBDATA* detprobdata);

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

class AbstractElementParser
{
public:
   explicit AbstractElementParser(SCIP* scip, NDecFileHandler& filehandler)
      : scip_(scip), filehandler_(filehandler), error_(false) {}

   virtual ~AbstractElementParser() = default;

   virtual void handleKeyValuePair(const char* name, json_t* value) = 0;

   virtual void handleValue(json_t* value) = 0;

   bool error() const
   {
      return error_;
   }

protected:
   NDecFileHandler& filehandler_;
   SCIP* scip_;
   bool error_;
};

struct NDecFileHandler
{
public:
   NDecFileHandler(SCIP* scip, const char* filename);
   NDecFileHandler(SCIP* scip, FILE* wfile);
   ~NDecFileHandler();

   SCIP_RETCODE initialize();

   bool parseElement(AbstractElementParser& elementparser, json_t* element);

   bool readNDec(NestedDecompositionData& data);

   bool writeNDec(PARTIALDECOMP* decomp);

private:
   bool serializeBlock(json_t* json, PARTIALDECOMP* decomp, int block);
   bool serializeBlockStructure(json_t* json, PARTIALDECOMP* decomp, BLOCK_STRUCTURE* blockstructure);
   bool serializeBlockStructureBlock(json_t* json, PARTIALDECOMP* decomp, BLOCK_STRUCTURE* blockstructure, int block);
   bool serializeDecomposition(json_t* json, PARTIALDECOMP* decomp);
   bool setObjectValue(const char* key, json_t* value, json_t* object = NULL, bool decref = true);
   bool appendArrayValue(json_t* value, json_t* array, bool decref = true);

   static size_t jsonLoadCallback(void* buffer, size_t buflen, void* data);
   static int jsonDumpCallback(const char* buffer, size_t buflen, void* data);

   SCIP_FILE* rfile_;
   FILE* wfile_;
   json_t* json_;
   json_error_t error_;
   SCIP* scip_;
};

class AbstractNestedDecompositionElementParser : public AbstractElementParser
{
public:
   AbstractNestedDecompositionElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data)
      : AbstractElementParser(scip, filehandler), data_(data) {}

   ~AbstractNestedDecompositionElementParser() override = default;

protected:
   DecompositionData* parseDecomposition(json_t* value);

   NestedDecompositionData& data_;
};

class BlockElementParser : public AbstractNestedDecompositionElementParser
{
public:
   BlockElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data,
      BlockData& blockdata) : AbstractNestedDecompositionElementParser(scip, filehandler, data), blockdata_(blockdata),
      parsingconstraints(false) {}

   ~BlockElementParser() override = default;

   void handleKeyValuePair(const char* name, json_t* value) override;

   void handleValue(json_t* value) override;

private:
   BlockData& blockdata_;
   bool parsingconstraints;
};

class DecompositionElementParser : public AbstractNestedDecompositionElementParser
{
public:
   DecompositionElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data,
      DecompositionData& decdata) : AbstractNestedDecompositionElementParser(scip, filehandler, data),
      decdata_(decdata), parsingmasterconstraints(false), parsingblocks(false), parsingsymmetry(false) {}

   ~DecompositionElementParser() override = default;

   void handleKeyValuePair(const char* name, json_t* value) override;

   void handleValue(json_t* value) override;

private:
   DecompositionData& decdata_;
   bool parsingmasterconstraints;
   bool parsingblocks;
   bool parsingsymmetry;
};

class RootElementParser : public AbstractNestedDecompositionElementParser
{
public:
   RootElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data)
      : AbstractNestedDecompositionElementParser(scip, filehandler, data) {}

   ~RootElementParser() override = default;

   void handleKeyValuePair(const char* name, json_t* value) override;

   void handleValue(json_t* value) override;
};

}

#endif //GCG_READER_NDEC_HPP__
