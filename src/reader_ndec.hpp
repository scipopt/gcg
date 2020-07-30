/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       */
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

#include "yaml.h"
#include <vector>
#include <unordered_map>
#include <string>

struct DecompositionData;
struct NDecFileHandler;

struct BlockData
{
   BlockData() : decomposition(NULL), symmetricalblock(-1) {}

   std::vector<std::string> constraints;
   DecompositionData* decomposition;
   int symmetricalblock;
};

struct DecompositionData
{
   DecompositionData() = default;
   ~DecompositionData();

   std::vector<std::string> masterconstraints;
   std::vector<BlockData*> blocks;
};

struct NestedDecompositionData
{
   NestedDecompositionData() : version(0), presolved(false), rootdecomposition(NULL) {}
   ~NestedDecompositionData();

   int version;
   std::string name;
   bool presolved;
   std::string comment;
   std::vector<DecompositionData*> decompositions;
   std::unordered_map<std::string, DecompositionData*> anchors;
   DecompositionData* rootdecomposition;
   std::unordered_map<std::string, std::string> symmetrydata;
};

class AbstractElementParser
{
public:
   explicit AbstractElementParser(SCIP* scip, NDecFileHandler& filehandler)
      : scip_(scip), filehandler_(filehandler), error_(false) {}

   virtual ~AbstractElementParser() = default;

   virtual bool handleMappingStart(const char* name, const char* anchor) = 0;

   virtual void handleMappingEnd() = 0;

   virtual bool handleSequenceStart(const char* name, const char* anchor) = 0;

   virtual void handleSequenceEnd() = 0;

   virtual void handleKeyValuePair(const char* name, const char* value, const char* anchor) = 0;

   virtual void handleKeyAliasPair(const char* name, const char* anchor) = 0;

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

   ~NDecFileHandler();

   void initialize();

   bool parseElement(AbstractElementParser& elementparser);

   bool readNDec(AbstractElementParser& rootparser);

private:
   static int yamlReadHandler(void* data, unsigned char* buffer, size_t size, size_t *size_read);

   SCIP_FILE* file_;
   yaml_parser_t parser_;
   SCIP* scip_;
};

class DummyElementParser : public AbstractElementParser
{
public:
   DummyElementParser(SCIP* scip, NDecFileHandler& filehandler)
      : AbstractElementParser(scip, filehandler) {}

   ~DummyElementParser() override = default;

   bool handleMappingStart(const char* name, const char* anchor) override
   {
      return false;
   }

   void handleMappingEnd() override {}

   bool handleSequenceStart(const char* name, const char* anchor) override
   {
      return false;
   }

   void handleSequenceEnd() override {}

   void handleKeyValuePair(const char* name, const char* value, const char* anchor) override {}

   void handleKeyAliasPair(const char* name, const char* anchor) override {}
};

class AbstractNestedDecompositionElementParser : public AbstractElementParser
{
public:
   AbstractNestedDecompositionElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data)
      : AbstractElementParser(scip, filehandler), data_(data) {}

   ~AbstractNestedDecompositionElementParser() override = default;

protected:
   void parseDecomposition(const char* anchor);

   DecompositionData* getDecompositionData(const char* anchor);

   void skipElement();

   NestedDecompositionData& data_;
};

class BlockElementParser : public AbstractNestedDecompositionElementParser
{
public:
   BlockElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data,
      BlockData& blockdata) : AbstractNestedDecompositionElementParser(scip, filehandler, data), blockdata_(blockdata),
      parsingconstraints(false) {}

   ~BlockElementParser() override = default;

   bool handleMappingStart(const char* name, const char* anchor) override;

   void handleMappingEnd() override {}

   bool handleSequenceStart(const char* name, const char* anchor) override;

   void handleSequenceEnd() override;

   void handleKeyValuePair(const char* name, const char* value, const char* anchor) override;

   void handleKeyAliasPair(const char* name, const char* anchor) override;

private:
   BlockData& blockdata_;
   bool parsingconstraints;
};

class DecompositionElementParser : public AbstractNestedDecompositionElementParser
{
public:
   DecompositionElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data,
      DecompositionData& decdata) : AbstractNestedDecompositionElementParser(scip, filehandler, data),
      decdata_(decdata), parsingmasterconstraints(false), parsingblocks(false) {}

   ~DecompositionElementParser() override = default;

   bool handleMappingStart(const char* name, const char* anchor) override;

   void handleMappingEnd() override {}

   bool handleSequenceStart(const char* name, const char* anchor) override;

   void handleSequenceEnd() override;

   void handleKeyValuePair(const char* name, const char* value, const char* anchor) override;

   void handleKeyAliasPair(const char* name, const char* anchor) override;

private:
   DecompositionData& decdata_;
   bool parsingmasterconstraints;
   bool parsingblocks;
};

class RootElementParser : public AbstractNestedDecompositionElementParser
{
public:
   RootElementParser(SCIP* scip, NDecFileHandler& filehandler, NestedDecompositionData& data)
      : AbstractNestedDecompositionElementParser(scip, filehandler, data), parsingdecomps(false),
      parsingsymmetry(false) {}

   ~RootElementParser() override = default;

   bool handleMappingStart(const char* name, const char* anchor) override;

   void handleMappingEnd() override;

   bool handleSequenceStart(const char* name, const char* anchor) override;

   void handleSequenceEnd() override;

   void handleKeyValuePair(const char* name, const char* value, const char* anchor) override;

   void handleKeyAliasPair(const char* name, const char* anchor) override;

private:
   bool parsingdecomps;
   bool parsingsymmetry;
};

#endif //GCG_READER_NDEC_HPP__
