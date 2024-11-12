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

/**@file   reader_jdec.cpp
 * @brief  jdec file reader for (JSON formatted) structure information
 * @author Erik Muehmer
 * @ingroup FILEREADERS
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// #define SCIP_DEBUG

#include "reader_jdec.h"
#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "cons_decomp.hpp"
#include "scip/scip_mem.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <jansson.h>

#define READER_NAME             "jdecreader"
#define READER_DESC             "jdec (JSON formatted structure information) file reader"
#define READER_EXTENSION        "jdec"

#define JDEC_VERSION             1

using namespace gcg;

/** checks version */
static constexpr bool checkVersion(int version)
{
   return version > 0 && version <= JDEC_VERSION;
}

/** checks return value of json lib calls */
static constexpr bool checkJson(int returnvalue)
{
   return returnvalue == 0;
}

/** data for dec reader */
struct SCIP_ReaderData
{
};

struct JDecDecompositionData;

/** struct to store block information (read from file) */
struct JDecBlockData
{
   /** constructor */
   JDecBlockData(
      int number                             /**< number of block */
      )
      : decomposition(NULL),
        symmetricalblock(number),
        blocknumber(number) {}

   /** deleted copy constructor */
   JDecBlockData(const JDecBlockData&) = delete;

   /** move constructor */
   JDecBlockData(
      JDecBlockData&& block                  /**< other block data */
      ) noexcept;

   /** destructor */
   ~JDecBlockData();

   /** deleted assignment operator */
   JDecBlockData& operator=(const JDecBlockData&) = delete;

   /** move assignment */
   JDecBlockData& operator=(
      JDecBlockData&& block                  /**< other block data */
      );

   /** less than operator to order block by their indices */
   bool operator<(
      const JDecBlockData& block             /**< other block to compare against */
      ) const;

   std::vector<std::string> constraints;     /**< names of constraints */
   JDecDecompositionData* decomposition;     /**< pointer to decomposition object */
   int symmetricalblock;                     /**< number of representative/symmetrical block */
   int blocknumber;                          /**< number of block/subproblem */
};

/** struct to store (nested) decomposition data (read from file) */
struct JDecDecompositionData
{
   /** constructor */
   JDecDecompositionData() : presolved(false), masterconstraints(), blocks(), symmetryvardata() {}

   /** destructor */
   ~JDecDecompositionData() = default;

   /** creates a block structure object */
   BLOCK_STRUCTURE* createBlockStructure(
      SCIP* scip,                            /**< scip object */
      DETPROBDATA* detprobdata               /**< detprobdata used to create the block structure object */
      );

   bool presolved;                                                /** is this a decomposition of a presolved model? */
   std::vector<std::string> masterconstraints;                    /**< vector containing names of master constraints */
   std::vector<JDecBlockData> blocks;                             /**< vector containing block data of each block */
   std::unordered_map<std::string, std::string> symmetryvardata;  /**< symmetry mapping for variables: name of variable -> name of its representative variable */
};

/** struct that stores the metadata of the read decomposition and a pointer to the actual decomposition */
struct JDecData
{
   /** constructor */
   JDecData() : version(0), rootdecomposition(NULL) {}
   
   /** deleted copy constructor */
   JDecData(const JDecData&) = delete;

   /** destructor */
   ~JDecData();

   /** deleted assignment operator */
   JDecData& operator=(const JDecData&) = delete;

   int version;                              /** version of jdec file */
   std::string name;                         /** decomposition's name */
   std::string description;                  /** decomposition's description */
   JDecDecompositionData* rootdecomposition; /** actual decomposition */
};

class AbstractJDecElementParser;

/** writes and reads jdec files */
class JDecFileHandler
{
public:
   /** constructor creating an object ready to read a jdec file */
   JDecFileHandler(
      SCIP* scip,                            /**< scip data structure */
      const char* filename                   /**< path of the jdec file to be read */
      );

   /** constructor creating an object ready to write a jdec file */
   JDecFileHandler(
      SCIP* scip,                            /**< scip data structure */
      FILE* wfile                            /**< file pointer used to write the jdec file */
      );

   /** destructor */
   ~JDecFileHandler();

   /** parses a json element using an element parser, returns a bool indicating success/failure */
   bool parseElement(
      AbstractJDecElementParser& elementparser, /**< the element parser */
      json_t* element                        /**< pointer to the json element */
      );

   /** 
    * reads a jdec file an stores the information in a data object, returns a bool indicating success/failure
    * 
    * @note must be constructed as reader
    */
   bool readJDec(
      JDecData& data                         /**< JDecData object used to store the read information */
      );

   /**
    * writes a partialldec to a jdec file, returns a bool indicating success/failure, returns a bool indicating success/failure
    * 
    * @note must be constructed as writer
    */
   bool writeJDec(
      PARTIALDECOMP* decomp                  /**< partialdec that will be written to the jdec file */
      );

private:
   /** initialize the handler */
   void initialize();

   /** serializes a specific block of a partialdec and puts it into an existing json object, returns a bool indicating success/failure */
   bool serializeBlock(
      json_t* json,                          /**< pointer to json struct that is already initialized as json object */
      PARTIALDECOMP* decomp,                 /**< pointer to partialdec */
      int block                              /**< number of block that should be serialized */
      );

   /** serializes a block structure and puts it into an existing json object, returns a bool indicating success/failure */
   bool serializeBlockStructure(
      json_t* json,                          /**< pointer to json struct that is already initialized as json object */
      PARTIALDECOMP* decomp,                 /**< pointer to partialdec the block structure belongs to */
      BLOCK_STRUCTURE* blockstructure        /**< pointer to block structure */
      );

   /** serializes a block of a block structure and puts it into an existing json object, returns a bool indicating success/failure */
   bool serializeBlockStructureBlock(
      json_t* json,                          /**< pointer to json struct that is already initialized as json object */
      PARTIALDECOMP* decomp,                 /**< pointer to partialdec the block structure belongs to */
      BLOCK_STRUCTURE* blockstructure,       /**< pointer to block structure */
      int block                              /**< number of block that should be serialized */
      );

   /** serializes a partialdec and puts it into an existing json object, returns a bool indicating success/failure */
   bool serializeDecomposition(
      json_t* json,                          /**< pointer to json struct that is already initialized as json object */
      PARTIALDECOMP* decomp                  /**< pointer to partialdec that should be serialized */
      );

   /** sets a value of the root or a provided json object for a specific key (key-value pair), returns a bool indicating success/failure */
   bool setObjectValue(
      const char* key,                       /**< the key used to store the value */
      json_t* value,                         /**< pointer to a json struct that will be set as value */
      json_t* object = NULL,                 /**< pointer to json object that will be modified, if NULL the root json object of the handler will be modified */
      bool decref = true                     /**< should the reference counter of the value json struct be decreased? (see JANSSON API documentation) */
      );

   /** appends a value to an json array, returns a bool indicating success/failure */
   bool appendArrayValue(
      json_t* value,                         /**< pointer to a json struct that will be set as value */
      json_t* array,                         /**< pointer to the json array */
      bool decref = true                     /**< should the reference counter of the value json struct be decreased? (see JANSSON API documentation) */
      );

   /** (static) callback provided to JANSSON used to read data: when called it writes at most buflen bytes to buffer and returns the number of bytes written */
   static size_t jsonLoadCallback(
      void* buffer,                          /**< pointer to buffer */
      size_t buflen,                         /**< length of buffer */
      void* data                             /**< pointer to user data, will point to the handler itself */
      );

   /** (static) callback provided to JANSSON used to write data: when called it writes the output contained in buffer to the jdec file and returns 0 on success or -1 else */
   static int jsonDumpCallback(
      const char* buffer,                    /**< pointer to buffer */
      size_t buflen,                         /**< length of buffer */
      void* data                             /**< pointer to user data, will point to the handler itself */
      );

   SCIP_FILE* rfile_;                        /**< SCIP file pointer to read from */
   FILE* wfile_;                             /**< file pointer to write to */
   json_t* json_;                            /**< root json object */
   json_error_t error_;                      /**< will contain error information of JANSSON if decoding fails */
   SCIP* scip_;                              /**< SCIP data structure */
};

/** abstract element parser class used to process data read by the jdec file handler */
class AbstractJDecElementParser
{
public:
   /** constructor */
   explicit AbstractJDecElementParser(
      SCIP* scip,                            /**< scip data structure */
      JDecFileHandler& filehandler           /**< jdec file handler that uses this parser*/
      ) : filehandler_(filehandler), scip_(scip), error_(false) {}

   /** destructor */
   virtual ~AbstractJDecElementParser() = default;

   /** abstract method that will be called for each key-value pair read by the handler */
   virtual void handleKeyValuePair(
      const char* name,                      /**< name/key of the value read */
      json_t* value                          /**< pointer to the value */
      ) = 0;

   /** abstract method that will be called for each value read by the handler */
   virtual void handleValue(
      json_t* value                          /**< pointer to the value */
      ) = 0;

   /** returns true if an error occured during parsing or false otherwise */
   bool error() const
   {
      return error_;
   }

protected:
   JDecFileHandler& filehandler_;            /**< the file handler that uses this parser */
   SCIP* scip_;                              /**< scip data structure */
   bool error_;                              /**< should be set to true if an error ocurred */
};

/** abstract decomposition element parser */
class AbstractJDecDecompositionElementParser : public AbstractJDecElementParser
{
public:
   /** constructor */
   AbstractJDecDecompositionElementParser(
      SCIP* scip,                            /**< scip data structure */
      JDecFileHandler& filehandler           /**< jdec file handler using this parser */
      ) : AbstractJDecElementParser(scip, filehandler) {}

   /** destructor */
   ~AbstractJDecDecompositionElementParser() override = default;

protected:
   /** parses/deserializes the decomposition contained by a json object and returns the created decomposition data*/
   JDecDecompositionData* parseDecomposition(
      json_t* value                          /**< the json object */
      );
};

/** block element parser, parses/deserializes blocks of decompositions */
class JDecBlockElementParser : public AbstractJDecDecompositionElementParser
{
public:
   /** constructor */
   JDecBlockElementParser(
      SCIP* scip,                            /**< scip data structure */
      JDecFileHandler& filehandler,          /**< file handler using this parser */
      JDecBlockData& blockdata               /**< block data structure the data is stored in */
      )
      : AbstractJDecDecompositionElementParser(scip, filehandler),
        blockdata_(blockdata),
        parsingconstraints(false)
      {}

   /** destructor */
   ~JDecBlockElementParser() override = default;

   /** implements the corresponding abstract function */
   void handleKeyValuePair(
      const char* name,
      json_t* value
      ) override;

   /** implements the corresponding abstract function */
   void handleValue(
      json_t* value
      ) override;

private:
   JDecBlockData& blockdata_;                /**< block data struct used to store the parsed data */
   bool parsingconstraints;                  /**< are we currently parsing the constraint array? */
};

/** actual decomposition element parser, parses/deserializes a decomposition */
class JDecDecompositionElementParser : public AbstractJDecDecompositionElementParser
{
public:
   /** constructor */
   JDecDecompositionElementParser(
      SCIP* scip,                            /**< scip data structure */
      JDecFileHandler& filehandler,          /**< file handler using this parser */
      JDecDecompositionData& decdata         /**< decomposition data structure the data is stored in */
      )
      : AbstractJDecDecompositionElementParser(scip, filehandler),
        decdata_(decdata),
        parsingmasterconstraints(false),
        parsingblocks(false),
        parsingsymmetry(false)
      {}

   /** destructor */
   ~JDecDecompositionElementParser() override = default;

   /** implements the corresponding abstract function */
   void handleKeyValuePair(
      const char* name,
      json_t* value
      ) override;

   /** implements the corresponding abstract function */
   void handleValue(
      json_t* value
      ) override;

private:
   JDecDecompositionData& decdata_;          /**< decomposition data struct used to store the parsed data */
   bool parsingmasterconstraints;            /**< are we currently parsing the master constraints? */
   bool parsingblocks;                       /**< are we currently parsing the blocks? */
   bool parsingsymmetry;                     /**< are we currently parsing the symmetry mapping of the variables? */
};

/** root element parser, parses/deserializes the root object of a jdec file */
class JDecRootElementParser : public AbstractJDecDecompositionElementParser
{
public:
   /** constructor */
   JDecRootElementParser(
      SCIP* scip,                            /**< scip data structure */
      JDecFileHandler& filehandler,          /**< file handler using this parser */
      JDecData& data                         /**< jdec data structure the data is stored in */
      ) : AbstractJDecDecompositionElementParser(scip, filehandler), data_(data) {}

   /** destructor */
   ~JDecRootElementParser() override = default;

   /** implements the corresponding abstract function */
   void handleKeyValuePair(
      const char* name,
      json_t* value
      ) override;

   /** implements the corresponding abstract function */
   void handleValue(
      json_t* value
      ) override;

private:
   JDecData& data_;                          /**< the jdec data structure used for storing the parsed data */
};

JDecBlockData::JDecBlockData(
   JDecBlockData&& block
   ) noexcept
{
   blocknumber = block.blocknumber;
   symmetricalblock = block.symmetricalblock;
   constraints = std::move(block.constraints);
   decomposition = block.decomposition;
   block.decomposition = NULL;
}

JDecBlockData::~JDecBlockData()
{
   delete decomposition;
}

JDecBlockData& JDecBlockData::operator=(
   JDecBlockData&& block
   )
{
   blocknumber = block.blocknumber;
   symmetricalblock = block.symmetricalblock;
   constraints = std::move(block.constraints);
   if( decomposition != NULL )
      delete decomposition;
   decomposition = block.decomposition;
   block.decomposition = NULL;
   return *this;
}

bool JDecBlockData::operator<(
   const JDecBlockData& block
   ) const
{
   return blocknumber < block.blocknumber;
}

JDecData::~JDecData()
{
   delete rootdecomposition;
}

BLOCK_STRUCTURE* JDecDecompositionData::createBlockStructure(
   SCIP* scip,
   DETPROBDATA* detprobdata
   )
{
   BLOCK_STRUCTURE* blockstructure = new BLOCK_STRUCTURE();
   int idx;

   if( presolved )
   {
      SCIPwarningMessage(scip, "Decomposition of blocks must not belong to a presolved model, ignoring.");
   }

   // master
   for( auto& cons : masterconstraints )
   {
      idx = detprobdata->getIndexForCons(cons.c_str());
      if( idx >= 0 )
         blockstructure->masterconss.push_back(idx);
   }

   // blocks
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
         blockstructure->blockstructures.push_back(blockdata.decomposition->createBlockStructure(scip, detprobdata));
      else
         blockstructure->blockstructures.emplace_back();
   }

   // symmetry
   if( !symmetryvardata.empty() )
   {
      bool success = true;
      for( auto& blockdata : blocks )
      {
         if( blockdata.symmetricalblock >= 0 && blockdata.symmetricalblock < (int)blocks.size() )
         {
            blockstructure->symmetricalblocks.push_back(blockdata.symmetricalblock);
         }
         else
         {
            SCIPwarningMessage(scip, "Got invalid block number: %d.\n", blockdata.symmetricalblock);
            success = false;
            break;
         }
      }

      if( success )
      {
         int idx2;
         blockstructure->symmetryvardata = std::move(blockstructure->symmetryvardata);
         for( auto& it : symmetryvardata )
         {
            idx = detprobdata->getIndexForVar(it.first.c_str());
            idx2 = detprobdata->getIndexForVar(it.second.c_str());
            if( idx >= 0 && idx2 >= 0 )
            {
               blockstructure->symmetryvardata.emplace(idx, idx2);
            }
            else
            {
               SCIPwarningMessage(scip, "Got invalid variable mapping: <%s> -> <%s>.\n", it.first.c_str(), it.second.c_str());
               success = false;
               blockstructure->symmetryvardata.clear();
               break;
            }
         }
      }

      if( !success )
      {
         SCIPwarningMessage(scip, "Could not set nested symmetry information.\n");
      }
   }

   return blockstructure;
}

JDecFileHandler::JDecFileHandler(
   SCIP* scip,
   const char* filename
   )
   : wfile_(NULL),
     json_(NULL),
     error_(),
     scip_(scip)
{
   rfile_ = SCIPfopen(filename, "r");
   initialize();
}

JDecFileHandler::JDecFileHandler(
   SCIP* scip,
   FILE* file
   )
   : rfile_(NULL),
     wfile_(file),
     json_(NULL),
     error_(),
     scip_(scip)
{
   initialize();
}

JDecFileHandler::~JDecFileHandler()
{
   if( json_ )
      json_decref(json_);
   if( rfile_ )
      SCIPfclose(rfile_);
}

void JDecFileHandler::initialize()
{
   if( rfile_ )
   {
      json_ = json_load_callback(&jsonLoadCallback, this, 0, &error_);
   }
   else
   {
      json_ = json_object();
   }
}

bool JDecFileHandler::parseElement(
   AbstractJDecElementParser& elementparser,
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

bool JDecFileHandler::readJDec(
   JDecData& data
   )
{
   bool error = false;

   if( rfile_ == NULL )
   {
      SCIPwarningMessage(scip_, "JSON parser is not initialized.");
      error = true;
   }
   else if( json_ == NULL )
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
      JDecRootElementParser rootparser(scip_, *this, data);
      error = !parseElement(rootparser, json_);
   }

   if( !error && !checkVersion(data.version) )
   {
      error = true;
      SCIPwarningMessage(scip_, "Invalid version.\n");
   }

   return !error;
}

bool JDecFileHandler::writeJDec(
   gcg::PARTIALDECOMP* decomp
   )
{
   bool success = (wfile_ != NULL && json_ != NULL);
   if( success )
   {
      success &= setObjectValue("version", json_integer(JDEC_VERSION));
      success &= setObjectValue("problem_name", json_string(SCIPgetProbName(scip_)));
      success &= setObjectValue("decomposition_id", json_integer(decomp->getID()));

      json_t* jsondecomp = json_object();
      success &= serializeDecomposition(jsondecomp, decomp);
      success &= setObjectValue("decomposition", jsondecomp);
   }
   else
   {
      SCIPwarningMessage(scip_, "JSON parser is not initialized.");
   }

   if( success )
   {
      success = checkJson(json_dump_callback(json_, &jsonDumpCallback, this, JSON_INDENT(2)));
   }
   return success;
}

bool JDecFileHandler::serializeBlock(
   json_t* json,
   gcg::PARTIALDECOMP* decomp,
   int block
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();
   json_t* jsonconstraints = json_array();

   success &= setObjectValue("index", json_integer(block), json);

   for( int i: decomp->getConssForBlock(block) )
   {
      auto* cons = detprobdata->getCons(i);
      success &= appendArrayValue(json_string(SCIPconsGetName(cons)), jsonconstraints);
   }
   success &= setObjectValue("constraints", jsonconstraints, json);

   if( decomp->aggInfoCalculated() )
   {
      success &= setObjectValue("symmetry_representative_block", json_integer(decomp->getReprBlockForEqClass(decomp->getEqClassForBlock(block))), json);
   }

   if( decomp->isNested() && decomp->getBlockStructure(block) != NULL )
   {
      json_t* jsonblockstructure = json_object();
      serializeBlockStructure(jsonblockstructure, decomp, decomp->getBlockStructure(block));
      success &= setObjectValue("decomposition", jsonblockstructure, json);
   }

   return success;
}

bool JDecFileHandler::serializeBlockStructure(
   json_t* json,
   gcg::PARTIALDECOMP* decomp,
   gcg::BLOCK_STRUCTURE* blockstructure
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();

   assert(blockstructure != NULL);

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

   return success;
}

bool JDecFileHandler::serializeBlockStructureBlock(
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

   // @todo: add "symmetry_representative_block"

   if( blockstructure->blockstructures[block] )
   {
      json_t* jsonblockstructure = json_object();
      serializeBlockStructure(jsonblockstructure, decomp, blockstructure->blockstructures[block]);
      success &= setObjectValue("decomposition", jsonblockstructure, json);
   }

   return success;
}

bool JDecFileHandler::serializeDecomposition(
   json_t* json,
   gcg::PARTIALDECOMP* decomp
   )
{
   bool success = true;
   auto* detprobdata = decomp->getDetprobdata();

   success &= setObjectValue("presolved", json_boolean(!decomp->isAssignedToOrigProb()), json);
   success &= setObjectValue("n_blocks", json_integer(decomp->getNBlocks()), json);

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

   if( decomp->aggInfoCalculated() && decomp->getNEquivalenceClasses() < decomp->getNBlocks() )
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
      success &= setObjectValue("symmetry_var_mapping", jsonsymmetry, json);
   }

   return success;
}

bool JDecFileHandler::setObjectValue(
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

bool JDecFileHandler::appendArrayValue(
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

size_t JDecFileHandler::jsonLoadCallback(
   void* buffer,
   size_t buflen,
   void* data
   )
{
   auto* filehandler = (JDecFileHandler*) data;
   size_t size_read = SCIPfread(buffer, 1, buflen, filehandler->rfile_);
   return (size_read == 0 && !SCIPfeof(filehandler->rfile_)) ? (size_t)-1 : size_read;
}

int JDecFileHandler::jsonDumpCallback(
   const char* buffer,
   size_t buflen,
   void* data
   )
{
   auto* filehandler = (JDecFileHandler*) data;
   assert(buflen <= INT_MAX);
   // not sure if this is the best function to call but SCIP's readers use it too
   SCIPinfoMessage(filehandler->scip_, filehandler->wfile_, "%.*s", (int)buflen, buffer);
   // size_t size_written = SCIPfwrite(buffer, 1, buflen, filehandler->wfile_);
   // return (size_written == 0) ? 1 : 0;
   return 0;
}

JDecDecompositionData* AbstractJDecDecompositionElementParser::parseDecomposition(
   json_t* value
   )
{
   JDecDecompositionData* decompdata = new JDecDecompositionData();
   JDecDecompositionElementParser decompositionparser(scip_, filehandler_, *decompdata);
   if( !filehandler_.parseElement(decompositionparser, value) )
      error_ = true;
   return decompdata;
}

void JDecRootElementParser::handleKeyValuePair(
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

void JDecRootElementParser::handleValue(
   json_t* value
   )
{
}

void JDecDecompositionElementParser::handleKeyValuePair(
   const char* name,
   json_t* value
   )
{
   if( parsingsymmetry )
   {
      if( json_is_string(value) )
      {
         decdata_.symmetryvardata.emplace(std::string(name), std::string(json_string_value(value)));
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
            if( !filehandler_.parseElement(*this, value) )
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
            if( !filehandler_.parseElement(*this, value) )
               error_ = true;
            // sort blocks as users can assign indices which may not be sorted
            std::sort(decdata_.blocks.begin(), decdata_.blocks.end());
            parsingblocks = false;
         }
         else
         {
            SCIPwarningMessage(scip_, "Blocks must be given as an array of objects.\n");
            error_ = true;
         }
      }
      else if( strcmp(name, "symmetry_var_mapping") == 0 )
      {
         if( json_is_object(value) )
         {
            parsingsymmetry = true;
            if( !filehandler_.parseElement(*this, value) )
               error_ = true;
            parsingsymmetry = false;
         }
         else
         {
            SCIPwarningMessage(scip_, "Symmetry information must be a mapping of strings.\n");
            error_ = true;
         }
      }
      else if( strcmp(name, "presolved") == 0 )
      {
         if( json_is_string(value) )
         {
            decdata_.presolved = strcmp(json_string_value(value), "true") == 0 ||
                              strcmp(json_string_value(value), "t") == 0 ||
                              strcmp(json_string_value(value), "yes") == 0 ||
                              strcmp(json_string_value(value), "y") == 0 ||
                              strcmp(json_string_value(value), "1") == 0;
         }
         else if ( json_is_boolean(value) )
         {
            decdata_.presolved = json_boolean_value(value);
         }
         else
         {
            SCIPwarningMessage(scip_, "Could not parse value of 'presolved'.");
            error_ = true;
         }
      }
      else
      {
         SCIPdebugMessage("Skipping unknown element '%s'\n", name);
      }
   }
}

void JDecDecompositionElementParser::handleValue(
   json_t* value
   )
{
   if( parsingblocks )
   {
      if( json_is_object(value) )
      {
         decdata_.blocks.emplace_back(decdata_.blocks.size());
         JDecBlockElementParser blockparser(scip_, filehandler_, decdata_.blocks.back());
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

void JDecBlockElementParser::handleKeyValuePair(
   const char* name,
   json_t* value
   )
{
   if( strcmp(name, "index") == 0 )
   {
      if( json_is_integer(value) )
      {
         blockdata_.blocknumber = (int)json_integer_value(value);
         assert(blockdata_.blocknumber >= 0);
      }
      else
      {
         SCIPwarningMessage(scip_, "Could not parse block index.\n");
         error_ = true;
      }
   }
   else if( strcmp(name, "symmetry_representative_block") == 0 )
   {
      if( json_is_integer(value) )
      {
         blockdata_.symmetricalblock = (int)json_integer_value(value);
      }
      else
      {
         SCIPwarningMessage(scip_, "Could not parse block number of representative block (symmetry).\n");
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

void JDecBlockElementParser::handleValue(
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

/* reads jdec file */
static
SCIP_RETCODE readJDec(
   SCIP*                 scip,               /**< scip data structure */
   const char*           filename,           /**< path of file */
   SCIP_RESULT*          result              /**< pointer to scip result */
   )
{
   JDecData data;
   JDecFileHandler filehandler(scip, filename);

   if( filehandler.readJDec(data) )
   {
      if( data.rootdecomposition )
      {
         int nblocks = (int)data.rootdecomposition->blocks.size();

         if( data.rootdecomposition->presolved && SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
         {
            SCIPinfoMessage(scip, NULL,
               "Reading presolved decomposition but problem is not presolved yet. Calling SCIPpresolve()\n");
            SCIPpresolve(scip);
         }

         PARTIALDECOMP* partialdec = new PARTIALDECOMP(scip, !data.rootdecomposition->presolved);
         DETPROBDATA* detprobdata = partialdec->getDetprobdata();
         for( auto& cons : data.rootdecomposition->masterconstraints )
         {
            if( !partialdec->fixConsToMasterByName(cons.c_str()) )
               SCIPwarningMessage(scip, "Could not set constraint %s as master constraint.\n", cons.c_str());
         }
         partialdec->setNBlocks(nblocks);
         for( int block = 0; block < nblocks; ++block )
         {
            JDecBlockData& blockdata = data.rootdecomposition->blocks[block];
            assert(block == blockdata.blocknumber);
            for( auto& cons : blockdata.constraints )
            {
               if( !partialdec->fixConsToBlockByName(cons.c_str(), block) )
                  SCIPwarningMessage(scip, "Could not set constraint %s as block constraint.\n", cons.c_str());
            }
            if( blockdata.decomposition )
            {
               BLOCK_STRUCTURE* nestedstructure = blockdata.decomposition->createBlockStructure(scip, detprobdata);
               partialdec->setBlockStructure(block, nestedstructure);
            }
            else
            {
               partialdec->setBlockStructure(block, NULL);
            }
         }
         GCGconshdlrDecompAddPreexisitingPartialDec(scip, partialdec);

         if( !data.rootdecomposition->symmetryvardata.empty() )
         {
            bool success = true;
            auto& symmetryvardata = data.rootdecomposition->symmetryvardata;

            // check symmetry data
            for( int b = 0; b < partialdec->getNBlocks() && success; ++b )
            {
               int symmetricalblock = data.rootdecomposition->blocks[b].symmetricalblock;
               for( int vi = 0; vi < partialdec->getNVarsForBlock(b) && success; ++vi )
               {
                  SCIP_VAR* var = detprobdata->getVar(partialdec->getVarsForBlock(b)[vi]);
                  assert(var != NULL);
                  const auto& it = symmetryvardata.find(SCIPvarGetName(var));
                  if( symmetricalblock == b )
                  {
                     success = (it == symmetryvardata.end() || it->first == it->second);
                     assert(success);
                  }
                  else if( it != symmetryvardata.end() )
                  {
                     SCIP_VAR* reprvar = SCIPfindVar(scip, it->second.c_str());
                     success = (reprvar != NULL &&
                        partialdec->getVarProbindexForBlock(detprobdata->getIndexForVar(reprvar), symmetricalblock) >= 0);
                     assert(success);
                  }
                  else
                  {
                     success = false;
                     assert(success);
                  }
               }
            }

            // if successful, set symmetry data
            if( success )
            {
               success = partialdec->setSymmetryInformation(
                  [&] (int b)
                  {
                     assert(b < (int)data.rootdecomposition->blocks.size());
                     return data.rootdecomposition->blocks[b].symmetricalblock;
                  },
                  [&] (int b, int vi)
                  {
                     SCIP_VAR* var = detprobdata->getVar(partialdec->getVarsForBlock(b)[vi]);
                     assert(var != NULL);
                     assert(symmetryvardata.find(SCIPvarGetName(var)) != symmetryvardata.end());
                     int ri = detprobdata->getIndexForVar(symmetryvardata[SCIPvarGetName(var)].c_str());
                     assert(partialdec->getVarProbindexForBlock(ri, data.rootdecomposition->blocks[b].symmetricalblock) >= 0);
                     return partialdec->getVarProbindexForBlock(ri, data.rootdecomposition->blocks[b].symmetricalblock);
                  }
               );
            }
            if( !success )
            {
               SCIPwarningMessage(scip, "Could not set symmetry information.\n");
            }
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

/* writes a jdec file for a given decomposition */
static
SCIP_RETCODE writePartialdec(
   SCIP*                 scip,               /**< scip data structure */
   FILE*                 file,               /**< file pointer */
   gcg::PARTIALDECOMP*   partialdec,         /**< partialdec to be written */
   SCIP_RESULT*          result              /**< pointer to scip result */
   )
{
   JDecFileHandler filehandler(scip, file);

   if( filehandler.writeJDec(partialdec) )
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

/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeJDec)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadJDec)
{  /*lint --e{715}*/

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL,
         "Please read in a problem before reading in the corresponding structure file!\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( readJDec(scip, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteJDec)
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

extern
SCIP_RETCODE SCIPincludeReaderJDec(
   SCIP*                 scip
   )
{
   SCIP_READERDATA* readerdata = NULL;

   /* create dec reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include jdec reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL, readerFreeJDec,
      readerReadJDec, readerWriteJDec, readerdata));

   return SCIP_OKAY;
}
