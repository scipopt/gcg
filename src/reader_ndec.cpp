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

#define SCIP_DEBUG

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

static bool checkVersion(int version)
{
   return version == NDEC_VERSION;
}

/** data for dec reader */
struct SCIP_ReaderData
{
};

/* reads ndec file */
SCIP_RETCODE readNDec(
   SCIP*                 scip,
   const char*           filename,
   SCIP_RESULT*          result
   )
{
   NestedDecompositionData data;
   NDecFileHandler filehandler(scip, filename);
   RootElementParser rootparser(scip, filehandler, data);
   filehandler.initialize();
   if( filehandler.readNDec(rootparser) )
   {
      if( data.rootdecomposition )
      {
         int nblocks = data.rootdecomposition->blocks.size();

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
         // todo: set symmetry information
         GCGconshdlrDecompAddPreexisitingPartialDec(scip, partialdec);
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
SCIP_RETCODE writePartialdec(
   SCIP*                 scip,
   FILE*                 file,
   gcg::PARTIALDECOMP*   partialdec,
   SCIP_RESULT*          result
   )
{
   return SCIP_NOTIMPLEMENTED;
}

NestedDecompositionData::~NestedDecompositionData()
{
   for( DecompositionData* data : decompositions )
      delete data;
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
   }
   return blockstructure;
}

NDecFileHandler::NDecFileHandler(
   SCIP* scip,
   const char* filename
   ) : scip_(scip), parser_()
{
   file_ = SCIPfopen(filename, "r");
}

NDecFileHandler::~NDecFileHandler()
{
   yaml_parser_delete(&parser_);
   SCIPfclose(file_);
}

void NDecFileHandler::initialize()
{
   yaml_parser_initialize(&parser_);
   yaml_parser_set_input(&parser_, &yamlReadHandler, this);
}

bool NDecFileHandler::parseElement(
   AbstractElementParser& elementparser
   )
{
   yaml_event_t* keyevent = NULL;
   bool error = false;
   int depth = 0;
   do
   {
      auto* event = new yaml_event_t();
      bool handled = false;
      if( !yaml_parser_parse(&parser_, event) )
         break;
      SCIPdebugMessage("YAML parser state: %d\n", parser_.state);
      SCIPdebugMessage("YAML start marker: %lu, %lu, %lu\n",
         event->start_mark.index, event->start_mark.line, event->start_mark.column);
      SCIPdebugMessage("YAML end marker: %lu, %lu, %lu\n",
         event->end_mark.index, event->end_mark.line, event->end_mark.column);
      switch( event->type )
      {
         case YAML_ALIAS_EVENT:
            SCIPdebugMessage("YAML event type: YAML_ALIAS_EVENT, anchor: %s\n", event->data.alias.anchor);
            if( keyevent )
            {
               elementparser.handleKeyAliasPair((char*) keyevent->data.scalar.value,
                  (char*) event->data.alias.anchor);
               yaml_event_delete(keyevent);
               delete keyevent;
               keyevent = NULL;
            }
            else
               elementparser.handleKeyAliasPair(NULL, (char*) event->data.alias.anchor);
            break;
         case YAML_SCALAR_EVENT:
            SCIPdebugMessage("YAML event type: YAML_SCALAR_EVENT, value: %s, anchor: %s\n",
               event->data.scalar.value, event->data.scalar.anchor);
            switch( parser_.state )
            {
               case YAML_PARSE_FLOW_MAPPING_KEY_STATE:
               case YAML_PARSE_BLOCK_MAPPING_KEY_STATE:
               case YAML_PARSE_FLOW_SEQUENCE_ENTRY_MAPPING_END_STATE:
                  if( keyevent )
                  {
                     elementparser.handleKeyValuePair((char*) keyevent->data.scalar.value,
                        (char*) event->data.scalar.value, (char*) event->data.scalar.anchor);
                     yaml_event_delete(keyevent);
                     delete keyevent;
                     keyevent = NULL;
                  }
                  else
                     elementparser.handleKeyValuePair(NULL, (char*) event->data.scalar.value,
                        (char*) event->data.scalar.anchor);
                  break;
               case YAML_PARSE_FLOW_MAPPING_VALUE_STATE:
               case YAML_PARSE_BLOCK_MAPPING_VALUE_STATE:
               case YAML_PARSE_FLOW_SEQUENCE_ENTRY_MAPPING_VALUE_STATE:
                  assert(keyevent == NULL);
                  keyevent = event;
                  event = NULL;
                  break;
               case YAML_PARSE_BLOCK_SEQUENCE_ENTRY_STATE:
               case YAML_PARSE_FLOW_SEQUENCE_ENTRY_STATE:
               case YAML_PARSE_INDENTLESS_SEQUENCE_ENTRY_STATE:
                  assert(keyevent == NULL);
                  elementparser.handleKeyValuePair(NULL, (char*) event->data.scalar.value,
                     (char*) event->data.scalar.anchor);
                  break;
               default:
                  SCIPwarningMessage(scip_, "State of parser is unexpected: %d\n", parser_.state);
                  error = true;
                  break;
            }
            break;
         case YAML_SEQUENCE_START_EVENT:
            SCIPdebugMessage("YAML event type: YAML_SEQUENCE_START_EVENT, anchor: %s\n",
               event->data.sequence_start.anchor);
            if( keyevent )
            {
               handled = elementparser.handleSequenceStart((char *) keyevent->data.scalar.value,
                  (char *) event->data.sequence_start.anchor);
               yaml_event_delete(keyevent);
               delete keyevent;
               keyevent = NULL;
            }
            else
               handled = elementparser.handleSequenceStart(NULL, (char *) event->data.sequence_start.anchor);
            if( !handled )
               depth++;
            break;
         case YAML_SEQUENCE_END_EVENT:
            SCIPdebugMessage("YAML event type: YAML_SEQUENCE_END_EVENT\n");
            depth--;
            elementparser.handleSequenceEnd();
            break;
         case YAML_MAPPING_START_EVENT:
            SCIPdebugMessage("YAML event type: YAML_MAPPING_START_EVENT, anchor: %s\n",
               event->data.mapping_start.anchor);
            if( keyevent )
            {
               handled = elementparser.handleMappingStart((char *) keyevent->data.scalar.value,
                  (char *) event->data.mapping_start.anchor);
               yaml_event_delete(keyevent);
               delete keyevent;
               keyevent = NULL;
            }
            else
               handled = elementparser.handleMappingStart(NULL, (char *) event->data.mapping_start.anchor);
            if( !handled )
               depth++;
            break;
         case YAML_MAPPING_END_EVENT:
            SCIPdebugMessage("YAML event type: YAML_MAPPING_END_EVENT\n");
            depth--;
            elementparser.handleMappingEnd();
            break;
         default:
            SCIPwarningMessage(scip_, "Received unexpected YAML event, type: %d\n", event->type);
            error = true;
            break;
      }
      if( event )
      {
         yaml_event_delete(event);
         delete event;
      }
      error |= elementparser.error();
   }
   while( depth >= 0 && !error && parser_.state != YAML_PARSE_END_STATE );
   return !error;
}

bool NDecFileHandler::readNDec(
   AbstractElementParser& rootparser
   )
{
   yaml_event_t event;
   bool error = false;

   do
   {
      if( !yaml_parser_parse(&parser_, &event) )
         break;
      SCIPdebugMessage("YAML parser state: %d\n", parser_.state);
      SCIPdebugMessage("YAML start marker: %lu, %lu, %lu\n",
         event.start_mark.index, event.start_mark.line, event.start_mark.column);
      SCIPdebugMessage("YAML end marker: %lu, %lu, %lu\n",
         event.end_mark.index, event.end_mark.line, event.end_mark.column);
      switch( event.type )
      {
         case YAML_STREAM_START_EVENT:
            SCIPdebugMessage("YAML event type: YAML_STREAM_START_EVENT\n");
            break;
         case YAML_STREAM_END_EVENT:
            SCIPdebugMessage("YAML event type: YAML_STREAM_END_EVENT\n");
            break;
         case YAML_DOCUMENT_START_EVENT:
            SCIPdebugMessage("YAML event type: YAML_DOCUMENT_START_EVENT\n");
            break;
         case YAML_DOCUMENT_END_EVENT:
            SCIPdebugMessage("YAML event type: YAML_DOCUMENT_END_EVENT\n");
            break;
         case YAML_MAPPING_START_EVENT:
            SCIPdebugMessage("YAML event type: YAML_MAPPING_START_EVENT, anchor: %s\n",
               event.data.mapping_start.anchor);
            error |= !parseElement(rootparser);
            error |= rootparser.error();
            break;
         default:
            SCIPwarningMessage(scip_, "Received unexpected YAML event, type: %d\n", event.type);
            error = true;
            break;
      }
      yaml_event_delete(&event);
   }
   while( !error && parser_.state != YAML_PARSE_END_STATE );

   if( parser_.error != YAML_NO_ERROR )
   {
      error = true;
      SCIPwarningMessage(scip_, "YAML error occurred:\n  problem: %s\n  context: %s\n",
         parser_.problem, parser_.context);
   }
   return !error;
}

int NDecFileHandler::yamlReadHandler(
   void *data,
   unsigned char *buffer,
   size_t size,
   size_t *size_read
   )
{
   auto* filehandler = (NDecFileHandler*) data;
   *size_read = SCIPfread(buffer, 1, size, filehandler->file_);
   return *size_read < 0 ? 0 : 1;
}

void AbstractNestedDecompositionElementParser::parseDecomposition(const char* anchor)
{
   DecompositionData* decompdata = new DecompositionData();
   data_.decompositions.push_back(decompdata);
   if( anchor )
      data_.anchors.emplace(std::string(anchor), decompdata);
   DecompositionElementParser decompositionparser(scip_, filehandler_, data_, *decompdata);
   if( !filehandler_.parseElement(decompositionparser) )
      error_ = true;
}

void AbstractNestedDecompositionElementParser::skipElement()
{
   DummyElementParser dummyparser(scip_, filehandler_);
   if( !filehandler_.parseElement(dummyparser) )
      error_ = true;
}

DecompositionData* AbstractNestedDecompositionElementParser::getDecompositionData(const char* anchor)
{
   DecompositionData* decdata = NULL;
   auto itr = data_.anchors.find(std::string(anchor));
   if( itr != data_.anchors.end() )
      decdata = itr->second;
   else
      SCIPwarningMessage(scip_, "Unknown decomposition anchor: %s\n", anchor);
   return decdata;
}

bool RootElementParser::handleMappingStart(
   const char* name,
   const char* anchor
   )
{
   bool skip = false;
   bool processed = false;
   if( parsingdecomps )
   {
      parseDecomposition(anchor);
      processed = true;
   }
   else if( name )
   {
      if( strcmp(name, "symmetry") == 0 )
      {
         parsingsymmetry = true;
      }
      else if( strcmp(name, "rootdecomposition") == 0 )
      {
         int idx = data_.decompositions.size();
         parseDecomposition(anchor);
         assert(idx < data_.decompositions.size());
         data_.rootdecomposition = data_.decompositions[idx];
         processed = true;
      }
      else
      {
         SCIPdebugMessage("Skipping unknown mapping element '%s'\n", name);
         skip = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown mapping element\n");
      skip = true;
   }

   if( skip )
   {
      skipElement();
      processed = true;
   }
   return processed;
}

void RootElementParser::handleMappingEnd()
{
   if( parsingsymmetry )
      parsingsymmetry = false;
}

bool RootElementParser::handleSequenceStart(
   const char* name,
   const char* anchor
   )
{
   bool skip = false;
   bool processed = false;
   if( name )
   {
      if( strcmp(name, "decompositions") == 0 )
      {
         parsingdecomps = true;
      }
      else
      {
         SCIPdebugMessage("Skipping unknown sequence element '%s'\n", name);
         skip = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown sequence element\n");
      skip = true;
   }

   if( skip )
   {
      skipElement();
      processed = true;
   }
   return processed;
}

void RootElementParser::handleSequenceEnd()
{
   if( parsingdecomps )
      parsingdecomps = false;
}

void RootElementParser::handleKeyValuePair(
   const char* name,
   const char* value,
   const char* anchor
   )
{
   if( parsingsymmetry )
   {
      assert(name);
      data_.symmetrydata.emplace(std::string(name), std::string(value));
   }
   else if( name )
   {
      if( strcmp(name, "version") == 0 )
      {
         try
         {
            data_.version = std::stoi(value);
            if( !checkVersion(data_.version))
            {
               SCIPwarningMessage(scip_, "Invalid version.\n");
               error_ = true;
               return;
            }
         }
         catch( const std::exception &e )
         {
            SCIPwarningMessage(scip_, "Could not parse version: %s\n", value);
            error_ = true;
            return;
         }
      }
      else if( strcmp(name, "name") == 0 )
      {
         data_.name = std::string(value);
      }
      else if( strcmp(name, "comment") == 0 )
      {
         data_.comment = std::string(value);
      }
      else if( strcmp(name, "presolved") == 0 )
      {
         data_.presolved = strcmp(value, "true") == 0 ||
                           strcmp(value, "t") == 0 ||
                           strcmp(value, "yes") == 0 ||
                           strcmp(value, "y") == 0 ||
                           strcmp(value, "1") == 0;
      }
   }
}

void RootElementParser::handleKeyAliasPair(
   const char* name,
   const char* anchor
   )
{
   assert(anchor);
   if( name )
   {
      if( strcmp(name, "rootdecomposition") == 0 )
      {
         data_.rootdecomposition = getDecompositionData(anchor);
      }
      else
         SCIPwarningMessage(scip_, "Only decomposition anchors are allowed.\n");
   }
}

bool DecompositionElementParser::handleMappingStart(
   const char* name,
   const char* anchor
   )
{
   bool skip = false;
   bool processed = false;
   if( parsingblocks )
   {
      decdata_.blocks.emplace_back();
      BlockElementParser blockparser(scip_, filehandler_, data_, decdata_.blocks.back());
      if( !filehandler_.parseElement(blockparser) )
         error_ = true;
      processed = true;
   }
   else if( name )
   {
      SCIPdebugMessage("Skipping unknown mapping element '%s'\n", name);
      skip = true;
   }
   else
   {
      SCIPdebugMessage("Skipping unknown mapping element\n");
      skip = true;
   }

   if( skip )
   {
      skipElement();
      processed = true;
   }
   return processed;
}

bool DecompositionElementParser::handleSequenceStart(
   const char* name,
   const char* anchor
   )
{
   bool skip = false;
   bool processed = false;
   if( name )
   {
      if( strcmp(name, "masterconstraints") == 0 )
      {
         parsingmasterconstraints = true;
      }
      else if( strcmp(name, "blocks") == 0 )
      {
         parsingblocks = true;
      }
      else
      {
         SCIPdebugMessage("Skipping unknown sequence element '%s'\n", name);
         skip = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown sequence element\n");
      skip = true;
   }

   if( skip )
   {
      skipElement();
      processed = true;
   }
   return processed;
}

void DecompositionElementParser::handleSequenceEnd()
{
   if( parsingmasterconstraints )
      parsingmasterconstraints = false;
   else if( parsingblocks )
      parsingblocks = false;
}

void DecompositionElementParser::handleKeyValuePair(
   const char* name,
   const char* value,
   const char* anchor
   )
{
   if( parsingmasterconstraints )
   {
      assert(value);
      decdata_.masterconstraints.emplace_back(std::string(value));
   }
}

void DecompositionElementParser::handleKeyAliasPair(
   const char* name,
   const char* anchor
   )
{
   SCIPwarningMessage(scip_, "Only decomposition anchors are allowed.\n");
}

bool BlockElementParser::handleMappingStart(
   const char* name,
   const char* anchor
   )
{
   bool skip = false;
   bool processed = false;
   if( name )
   {
      if( strcmp(name, "decomposition") == 0 )
      {
         int idx = data_.decompositions.size();
         parseDecomposition(anchor);
         assert(idx < data_.decompositions.size());
         blockdata_.decomposition = data_.decompositions[idx];
         processed = true;
      }
      else
      {
         SCIPdebugMessage("Skipping unknown mapping element '%s'\n", name);
         skip = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown mapping element\n");
      skip = true;
   }

   if( skip )
   {
      skipElement();
      processed = true;
   }
   return processed;
}

bool BlockElementParser::handleSequenceStart(
   const char* name,
   const char* anchor
   )
{
   bool skip = false;
   bool processed = false;
   if( name )
   {
      if( strcmp(name, "constraints") == 0 )
      {
         parsingconstraints = true;
      }
      else
      {
         SCIPdebugMessage("Skipping unknown sequence element '%s'\n", name);
         skip = true;
      }
   }
   else
   {
      SCIPdebugMessage("Skipping unknown sequence element\n");
      skip = true;
   }

   if( skip )
   {
      skipElement();
      processed = true;
   }
   return processed;
}

void BlockElementParser::handleSequenceEnd()
{
   if( parsingconstraints )
      parsingconstraints = false;
}

void BlockElementParser::handleKeyValuePair(
   const char* name,
   const char* value,
   const char* anchor
   )
{
   if( parsingconstraints )
   {
      assert(value);
      blockdata_.constraints.emplace_back(std::string(value));
   }
   else if( name )
   {
      if( strcmp(name, "symmetrical_block") == 0 )
      {
         try
         {
            blockdata_.symmetricalblock = std::stoi(value);
         }
         catch( const std::exception &e )
         {
            SCIPwarningMessage(scip_, "Could not parse block number: %s\n", value);
            error_ = true;
         }
      }
   }
}

void BlockElementParser::handleKeyAliasPair(
   const char* name,
   const char* anchor
   )
{
   assert(anchor);
   if( name && strcmp(name, "decomposition") == 0 )
   {
      blockdata_.decomposition = getDecompositionData(anchor);
   }
   else
      SCIPwarningMessage(scip_, "Only decomposition anchors are allowed.\n");
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

   gcg::PARTIALDECOMP* partialdec = DECgetPartialdecToWrite(scip, transformed);

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
   SCIP_READERDATA* readerdata;

   /* create dec reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include dec reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL, readerFreeNDec,
      readerReadNDec, readerWriteNDec, readerdata));

   return SCIP_OKAY;
}
