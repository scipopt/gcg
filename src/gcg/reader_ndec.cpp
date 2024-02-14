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

/**@file   reader_ndec.cpp
 * @brief  ndec file reader for (nested) structure information
 * @author Erik Muehmer
 * @ingroup FILEREADERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// #define SCIP_DEBUG

#include "reader_ndec.h"
#include "reader_ndec.hpp"
#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.hpp"
#include "scip/scip_mem.h"

#include <cassert>
#include <string>

#define READER_NAME             "ndecreader"
#define READER_DESC             "file reader for blocks in ndec format"
#define READER_EXTENSION        "ndec"

#define NDEC_VERSION             1

using namespace gcg;

static constexpr bool checkVersion(int version)
{
   return version <= NDEC_VERSION;
}

static constexpr bool checkJson(int returnvalue)
{
   return returnvalue == 0;
}

/** data for dec reader */
struct SCIP_ReaderData
{
};

/* reads ndec file */
static
SCIP_RETCODE readNDec(
   SCIP*                 scip,
   const char*           filename,
   SCIP_RESULT*          result
   )
{
   NestedDecompositionData data;
   NDecFileHandler filehandler(scip, filename);
   SCIP_CALL( filehandler.initialize() );
   if( filehandler.readNDec(data) )
   {
      if( data.rootdecomposition )
      {
         int nblocks = (int)data.rootdecomposition->blocks.size();

         if( data.presolved && SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
         {
            SCIPinfoMessage(scip, NULL,
               "Reading presolved decomposition but problem is not presolved yet. Calling SCIPpresolve()\n");
            SCIPpresolve(scip);
         }

         PARTIALDECOMP* partialdec = new PARTIALDECOMP(scip, !data.presolved);
         DETPROBDATA* detprobdata = partialdec->getDetprobdata();
         for( auto& cons : data.rootdecomposition->masterconstraints )
         {
            if( !partialdec->fixConsToMasterByName(cons.c_str()) )
               SCIPwarningMessage(scip, "Could not set constraint %s as master constraint.\n", cons.c_str());
         }
         partialdec->setNBlocks(nblocks);
         for( int block = 0; block < nblocks; ++block )
         {
            BlockData& blockdata = data.rootdecomposition->blocks[block];
            for( auto& cons : blockdata.constraints )
            {
               if( !partialdec->fixConsToBlockByName(cons.c_str(), block) )
                  SCIPwarningMessage(scip, "Could not set constraint %s as block constraint.\n", cons.c_str());
            }
            if( blockdata.decomposition )
            {
               BLOCK_STRUCTURE* nestedstructure = blockdata.decomposition->createBlockStructure(detprobdata);
               partialdec->setBlockStructure(block, nestedstructure);
            }
         }
         GCGconshdlrDecompAddPreexisitingPartialDec(scip, partialdec);

         bool success = partialdec->setSymmetryInformation(
            [&data] (int b)
            {
               assert(b < (int)data.rootdecomposition->blocks.size());
               return data.rootdecomposition->blocks[b].symmetricalblock;
            },
            [&data, detprobdata, partialdec] (int b, int vi)
            {
               SCIP_VAR* var = detprobdata->getVar(partialdec->getVarsForBlock(b)[vi]);
               assert(var != NULL);
               assert(data.rootdecomposition->symmetrydata.find(SCIPvarGetName(var)) != data.rootdecomposition->symmetrydata.end());
               int ri = detprobdata->getIndexForVar(data.rootdecomposition->symmetrydata[SCIPvarGetName(var)].c_str());
               assert(partialdec->getVarProbindexForBlock(ri, data.rootdecomposition->blocks[b].symmetricalblock) >= 0);
               return partialdec->getVarProbindexForBlock(ri, data.rootdecomposition->blocks[b].symmetricalblock);
            }
         );
         if( !success )
         {
            SCIPwarningMessage(scip, "Could not set symmetry information.\n");
         }
      }
      else
         SCIPwarningMessage(scip, "No root decomposition is specified.\n");
      *result = SCIP_SUCCESS;
   }
   else
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/* write a ndec file for a given decomposition */
static
SCIP_RETCODE writePartialdec(
   SCIP*                 scip,
   FILE*                 file,
   gcg::PARTIALDECOMP*   partialdec,
   SCIP_RESULT*          result
   )
{
   NDecFileHandler filehandler(scip, file);
   SCIP_CALL( filehandler.initialize() );
   if( filehandler.writeNDec(partialdec) )
   {
      *result = SCIP_SUCCESS;
   }
   else
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_WRITEERROR;
   }
   return SCIP_OKAY;
}

BlockData::BlockData(BlockData&& block) noexcept
{
   symmetricalblock = block.symmetricalblock;
   constraints = std::move(block.constraints);
   decomposition = block.decomposition;
   block.decomposition = NULL;
}

BlockData::~BlockData()
{
   delete decomposition;
}

NestedDecompositionData::~NestedDecompositionData()
{
   delete rootdecomposition;
}

BLOCK_STRUCTURE* DecompositionData::createBlockStructure(
   DETPROBDATA* detprobdata
   )
{
   BLOCK_STRUCTURE* blockstructure = new BLOCK_STRUCTURE();
   int idx;
   for( auto& cons : masterconstraints )
   {
      idx = detprobdata->getIndexForCons(cons.c_str());
      if( idx >= 0 )
         blockstructure->masterconss.push_back(idx);

   }
   for( auto& blockdata : blocks )
   {
      blockstructure->blockconss.emplace_back();
      for( auto& cons : blockdata.constraints )
      {
         idx = detprobdata->getIndexForCons(cons.c_str());
         if( idx >= 0 )
            blockstructure->blockconss.back().push_back(idx);
      }
      if( blockdata.decomposition )
         blockstructure->blockstructures.push_back(blockdata.decomposition->createBlockStructure(detprobdata));
      else
         blockstructure->blockstructures.emplace_back();
   }
   // @todo: set nested symmetry information
   return blockstructure;
}

NDecFileHandler::NDecFileHandler(
   SCIP* scip,
   const char* filename
   ) : scip_(scip), json_(NULL), error_(), wfile_(NULL)
{
   rfile_ = SCIPfopen(filename, "r");
}

NDecFileHandler::NDecFileHandler(
   SCIP* scip,
   FILE* file
   ) : scip_(scip), json_(NULL), error_(), wfile_(file), rfile_(NULL)
{
}

NDecFileHandler::~NDecFileHandler()
{
   if( json_ )
      json_decref(json_);
   if( rfile_ )
      SCIPfclose(rfile_);
}

SCIP_RETCODE NDecFileHandler::initialize()
{
   if( rfile_ )
   {
      json_ = json_load_callback(&jsonLoadCallback, this, 0, &error_);
   }
   else
   {
      json_ = json_object();
   }
   return SCIP_OKAY;
}

bool NDecFileHandler::parseElement(
   AbstractElementParser& elementparser,
   json_t* element
   )
{
   bool error = false;
   json_t* value;

   if( json_is_object(element) )
   {
      const char* key;
      json_object_foreach(element, key, value)
      {
         elementparser.handleKeyValuePair(key, value);
      }
   }
   else if (json_is_array(element))
   {
      size_t index;
      json_array_foreach(element, index, value)
      {
         elementparser.handleValue(value);
      }
   }
   else
   {
      SCIPwarningMessage(scip_, "Unexpected JSON type: %d\n", json_typeof(element));
      error = true;
   }

   error |= elementparser.error();
   return !error;
}

bool NDecFileHandler::readNDec(
   NestedDecompositionData& data
   )
{
   bool error = false;

   if( !rfile_ )
   {
      SCIPwarningMessage(scip_, "JSON parser is not initialized.");
      error = true;
   }
   else if( !json_ )
   {
      SCIPwarningMessage(scip_, "Could not parse JSON, line %d: %s\n", error_.line, error_.text);
      error = true;
   }
   else if( !json_is_object(json_) )
   {
      SCIPwarningMessage(scip_, "Decomposition is invalid (root has to be an object).\n");
      error = true;
   }
   else
   {
      RootElementParser rootparser(scip_, *this, data);
      error = !parseElement(rootparser, json_);
   }

   return !error;
}

bool NDecFileHandler::writeNDec(
   gcg::PARTIALDECOMP* decomp
   )
{
   bool success = true;
   success &= setObjectValue("version", json_integer(NDEC_VERSION));
   success &= setObjectValue("problem_name", json_string(SCIPgetProbName(scip_)));
   success &= setObjectValue("decomp_id", json_integer(decomp->getID()));
   success &= setObjectValue("presolved", json_boolean(!decomp->isAssignedToOrigProb()));
   success &= setObjectValue("num_blocks", json_integer(decomp->getNBlocks()));

   json_t* jsondecomp = json_object();
   success &= serializeDecomposition(jsondecomp, decomp);
   success &= setObjectValue("decomposition", jsondecomp);

   if( success )
   {
      success = checkJson(json_dump_callback(json_, &jsonDumpCallback, this, JSON_INDENT(2)));
   }
   return success;
}

bool NDecFileHandler::serializeBlock(
   json_t* json,
   gcg::PARTIALDECOMP* decomp,
   int block
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();
   json_t* jsonconstraints = json_array();

   for( int i: decomp->getConssForBlock(block) )
   {
      auto* cons = detprobdata->getCons(i);
      success &= appendArrayValue(json_string(SCIPconsGetName(cons)), jsonconstraints);
   }
   success &= setObjectValue("constraints", jsonconstraints, json);

   if( decomp->aggInfoCalculated() )
   {
      success &= setObjectValue("symmetrical_block", json_integer(decomp->getReprBlockForEqClass(decomp->getEqClassForBlock(block))), json);
   }

   if( decomp->isNested() )
   {
      json_t* jsonblockstructure = json_object();
      serializeBlockStructure(jsonblockstructure, decomp, decomp->getBlockStructure(block));
      success &= setObjectValue("decomposition", jsonblockstructure, json);
   }

   return success;
}

bool NDecFileHandler::serializeBlockStructure(
   json_t* json,
   gcg::PARTIALDECOMP* decomp,
   gcg::BLOCK_STRUCTURE* blockstructure
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();

   json_t* jsonmasterconstraints = json_array();
   for( int i: blockstructure->masterconss )
   {
      auto* cons = detprobdata->getCons(i);
      success &= appendArrayValue(json_string(SCIPconsGetName(cons)), jsonmasterconstraints);
   }
   success &= setObjectValue("master_constraints", jsonmasterconstraints, json);

   json_t* jsonblocks = json_array();
   for( int b = 0; b < (int)blockstructure->blockconss.size(); ++b )
   {
      json_t* jsonblock = json_object();
      success &= serializeBlockStructureBlock(jsonblock, decomp, blockstructure, b);
      success &= appendArrayValue(jsonblock, jsonblocks);
   }
   success &= setObjectValue("blocks", jsonblocks, json);

   // @todo: add "symmetry_mapping"

   return success;
}

bool NDecFileHandler::serializeBlockStructureBlock(
   json_t* json,
   PARTIALDECOMP* decomp,
   BLOCK_STRUCTURE* blockstructure,
   int block
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();
   json_t* jsonconstraints = json_array();

   for( int i: blockstructure->blockconss[block] )
   {
      auto* cons = detprobdata->getCons(i);
      success &= appendArrayValue(json_string(SCIPconsGetName(cons)), jsonconstraints);
   }
   success &= setObjectValue("constraints", jsonconstraints, json);

   // @todo: add "symmetrical_block"

   if( blockstructure->blockstructures[block] )
   {
      json_t* jsonblockstructure = json_object();
      serializeBlockStructure(jsonblockstructure, decomp, blockstructure->blockstructures[block]);
      success &= setObjectValue("decomposition", jsonblockstructure, json);
   }

   return success;
}

bool NDecFileHandler::serializeDecomposition(
   json_t* json,
   gcg::PARTIALDECOMP* decomp
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();

   json_t* jsonmasterconstraints = json_array();
   for( int i: decomp->getMasterconss() )
   {
      auto* cons = detprobdata->getCons(i);
      success &= appendArrayValue(json_string(SCIPconsGetName(cons)), jsonmasterconstraints);
   }
   success &= setObjectValue("master_constraints", jsonmasterconstraints, json);

   if( !decomp->aggInfoCalculated() )
   {
      decomp->calcAggregationInformation(true);
   }

   json_t* jsonblocks = json_array();
   for( int b = 0; b < decomp->getNBlocks(); ++b )
   {
      json_t* jsonblock = json_object();
      success &= serializeBlock(jsonblock, decomp, b);
      success &= appendArrayValue(jsonblock, jsonblocks);
   }
   success &= setObjectValue("blocks", jsonblocks, json);

   if( decomp->aggInfoCalculated() )
   {
      json_t* jsonsymmetry = json_object();
      for( int ec = 0; ec < decomp->getNEquivalenceClasses(); ++ec )
      {
         int repblock = decomp->getReprBlockForEqClass(ec);
         auto& eqclassblocks = decomp->getBlocksForEqClass(ec);
         for( int i = 0; i < (int) eqclassblocks.size(); ++i )
         {
            int b = eqclassblocks[i];
            if( b == repblock )
               continue;
            for( int vi = 0; vi < (int)decomp->getRepVarmap(ec, i).size(); ++vi )
            {
               int rvi = decomp->getRepVarmap(ec, i)[vi];
               auto* var = detprobdata->getVar(decomp->getVarsForBlock(b)[vi]);
               auto* repvar = detprobdata->getVar(decomp->getVarsForBlock(repblock)[rvi]);
               success &= setObjectValue(SCIPvarGetName(var), json_string(SCIPvarGetName(repvar)), jsonsymmetry);
            }
         }
      }
      success &= setObjectValue("symmetry_mapping", jsonsymmetry, json);
   }

   return success;
}

bool NDecFileHandler::setObjectValue(
   const char* key,
   json_t* value,
   json_t* object,
   bool decref
   )
{
   bool success;

   if( !object )
      object = json_;

   if( decref )
   {
      success = checkJson(json_object_set_new(object, key, value));
   }
   else
   {
      success = checkJson(json_object_set(object, key, value));
   }

   if( !success )
   {
      SCIPwarningMessage(scip_, "Could not set value with key '%s'\n", key);
   }
   return success;
}

bool NDecFileHandler::appendArrayValue(
   json_t* value,
   json_t* array,
   bool decref
   )
{
   bool success;

   if( decref )
   {
      success = checkJson(json_array_append_new(array, value));
   }
   else
   {
      success = checkJson(json_array_append(array, value));
   }

   if( !success )
   {
      SCIPwarningMessage(scip_, "Could not append value.\n");
   }
   return success;
}

size_t NDecFileHandler::jsonLoadCallback(
   void* buffer,
   size_t buflen,
   void* data
   )
{
   auto* filehandler = (NDecFileHandler*) data;
   size_t size_read = SCIPfread(buffer, 1, buflen, filehandler->rfile_);
   return (size_read == 0 && !SCIPfeof(filehandler->rfile_)) ? (size_t)-1 : size_read;
}

int NDecFileHandler::jsonDumpCallback(
   const char* buffer,
   size_t buflen,
   void* data
   )
{
   auto* filehandler = (NDecFileHandler*) data;
   assert(buflen <= INT_MAX);
   SCIPinfoMessage(filehandler->scip_, filehandler->wfile_, "%.*s", (int)buflen, buffer);
   // size_t size_written = SCIPfwrite(buffer, 1, buflen, filehandler->wfile_);
   // return (size_written == 0) ? 1 : 0;
   return 0;
}

DecompositionData* AbstractNestedDecompositionElementParser::parseDecomposition(json_t* value)
{
   DecompositionData* decompdata = new DecompositionData();
   DecompositionElementParser decompositionparser(scip_, filehandler_, data_, *decompdata);
   if( !filehandler_.parseElement(decompositionparser, value) )
      error_ = true;
   return decompdata;
}

void RootElementParser::handleKeyValuePair(
   const char* name,
   json_t* value
   )
{
   if( strcmp(name, "version") == 0 )
   {
      if( json_is_integer(value) )
      {
         data_.version = (int) json_integer_value(value);
         if( !checkVersion(data_.version))
         {
            SCIPwarningMessage(scip_, "Invalid version.\n");
            error_ = true;
         }
      }
      else
      {
         SCIPwarningMessage(scip_, "Version must be an integer.");
         error_ = true;
      }
   }
   else if( strcmp(name, "name") == 0 )
   {
      if( json_is_string(value) )
      {
         data_.name = std::string(json_string_value(value));
      }
      else
      {
         SCIPwarningMessage(scip_, "Decomposition name must be a string.");
         error_ = true;
      }
   }
   else if( strcmp(name, "description") == 0 )
   {
      if( json_is_string(value) )
      {
         data_.description = std::string(json_string_value(value));
      }
   }
   else if( strcmp(name, "presolved") == 0 )
   {
      if( json_is_string(value) )
      {
         data_.presolved = strcmp(json_string_value(value), "true") == 0 ||
                           strcmp(json_string_value(value), "t") == 0 ||
                           strcmp(json_string_value(value), "yes") == 0 ||
                           strcmp(json_string_value(value), "y") == 0 ||
                           strcmp(json_string_value(value), "1") == 0;
      }
      else if ( json_is_boolean(value) )
      {
         data_.presolved = json_boolean_value(value);
      }
      else
      {
         SCIPwarningMessage(scip_, "Could not parse value of 'presolved'.");
         error_ = true;
      }
   }
   else if( strcmp(name, "decomposition") == 0 )
   {
      if( json_is_object(value) )
      {
         data_.rootdecomposition = parseDecomposition(value);
      }
      else
      {
         SCIPwarningMessage(scip_, "Decomposition must be an object.\n");
         error_ = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown element '%s'.\n", name);
   }
}

void RootElementParser::handleValue(
   json_t* value
   )
{
}

void DecompositionElementParser::handleKeyValuePair(
   const char* name,
   json_t* value
   )
{
   if( parsingsymmetry )
   {
      if( json_is_string(value) )
      {
         decdata_.symmetrydata.emplace(std::string(name), std::string(json_string_value(value)));
      }
      else
      {
         SCIPwarningMessage(scip_, "Symmetry information must consist of strings.");
         error_ = true;
      }
   }
   else
   {
      if( strcmp(name, "master_constraints") == 0 )
      {
         if( json_is_array(value))
         {
            parsingmasterconstraints = true;
            if( !filehandler_.parseElement(*this, value))
               error_ = true;
            parsingmasterconstraints = false;
         }
         else
         {
            SCIPwarningMessage(scip_, "Constraints must be given as an array of strings.\n");
            error_ = true;
         }
      }
      else if( strcmp(name, "blocks") == 0 )
      {
         if( json_is_array(value))
         {
            parsingblocks = true;
            if( !filehandler_.parseElement(*this, value))
               error_ = true;
            parsingblocks = false;
         }
         else
         {
            SCIPwarningMessage(scip_, "Blocks must be given as an array of objects.\n");
            error_ = true;
         }
      }
      else if( strcmp(name, "symmetry_mapping") == 0 )
      {
         if( json_is_object(value) )
         {
            parsingsymmetry = true;
            if( !filehandler_.parseElement(*this, value))
               error_ = true;
            parsingsymmetry = false;
         }
         else
         {
            SCIPwarningMessage(scip_, "Symmetry information must be a mapping of strings.\n");
            error_ = true;
         }
      }
      else
      {
         SCIPdebugMessage("Skipping unknown element '%s'\n", name);
      }
   }
}

void DecompositionElementParser::handleValue(
   json_t* value
   )
{
   if( parsingblocks )
   {
      if( json_is_object(value) )
      {
         decdata_.blocks.emplace_back();
         BlockElementParser blockparser(scip_, filehandler_, data_, decdata_.blocks.back());
         if( !filehandler_.parseElement(blockparser, value))
            error_ = true;
      }
      else
      {
         SCIPwarningMessage(scip_, "Block must be an object.\n");
         error_ = true;
      }
   }
   else if( parsingmasterconstraints )
   {
      if( json_is_string(value) )
      {
         decdata_.masterconstraints.emplace_back(json_string_value(value));
      }
      else
      {
         SCIPwarningMessage(scip_, "Constraints must be given as an array of strings.\n");
         error_ = true;
      }
   }
}

void BlockElementParser::handleKeyValuePair(
   const char* name,
   json_t* value
   )
{
   if( strcmp(name, "symmetrical_block") == 0 )
   {
      if( json_is_integer(value) )
      {
         blockdata_.symmetricalblock = (int)json_integer_value(value);
      }
      else
      {
         SCIPwarningMessage(scip_, "Could not parse block number.\n");
         error_ = true;
      }
   }
   else if( strcmp(name, "decomposition") == 0 )
   {
      if( json_is_object(value) )
      {
         blockdata_.decomposition = parseDecomposition(value);
      }
      else
      {
         SCIPwarningMessage(scip_, "Decomposition must be an object.\n");
         error_ = true;
      }
   }
   else if( strcmp(name, "constraints") == 0 )
   {
      if( json_is_array(value) )
      {
         parsingconstraints = true;
         filehandler_.parseElement(*this, value);
         parsingconstraints = false;
      }
      else
      {
         SCIPwarningMessage(scip_, "Constraints must be given as an array of strings.\n");
         error_ = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown element '%s'\n", name);
   }
}

void BlockElementParser::handleValue(
   json_t* value
   )
{
   if( parsingconstraints )
   {
      if( json_is_string(value) )
      {
         blockdata_.constraints.emplace_back(json_string_value(value));
      }
      else
      {
         SCIPwarningMessage(scip_, "Constraints must be given as an array of strings.\n");
         error_ = true;
      }
   }
}

/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeNDec)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadNDec)
{  /*lint --e{715}*/

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL,
         "Please read in a problem before reading in the corresponding structure file!\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( readNDec(scip, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteNDec)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);

   gcg::PARTIALDECOMP* partialdec = GCGgetPartialdecToWrite(scip, transformed);

   if(partialdec == NULL) {
      SCIPwarningMessage(scip, "There is no writable partialdec!\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( writePartialdec(scip, file, partialdec, result) );

   return SCIP_OKAY;
}

/* includes the ndec file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderNDec(
   SCIP*                 scip
   )
{
   SCIP_READERDATA* readerdata = NULL;

   /* create dec reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include dec reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL, readerFreeNDec,
      readerReadNDec, readerWriteNDec, readerdata));

   return SCIP_OKAY;
}
