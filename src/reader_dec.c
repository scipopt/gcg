/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file   reader_dec.c
 * @brief  DEC file reader
 * @ingroup FILEREADERS
 * @author Lukas Kirchhart
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h>
#endif
#include <ctype.h>
#include <struct_decomp.h>
#include <type_decomp.h>

#include "reader_dec.h"
#include "relax_gcg.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip_misc.h"

#include "pub_gcgvar.h"
#include "cons_decomp.h"
#include "pub_decomp.h"

#define READER_NAME             "decreader"
#define READER_DESC             "file reader for blocks corresponding to a mip format"
#define READER_EXTENSION        "dec"


/*
 * Data structures
 */
#define DEC_MAX_LINELEN       65536
#define DEC_MAX_PUSHEDTOKENS  2
#define DEC_PRINTLEN          100

/** Section in DEC File */
enum DecSection
{
   DEC_START, DEC_NBLOCKS, DEC_BLOCK, DEC_MASTERCONSS, DEC_END
};
typedef enum DecSection DECSECTION;

enum DecExpType
{
   DEC_EXP_NONE, DEC_EXP_UNSIGNED, DEC_EXP_SIGNED
};
typedef enum DecExpType DECEXPTYPE;

/** DEC reading data */
struct DecInput
{
   SCIP_FILE* file;
   char linebuf[DEC_MAX_LINELEN];
   char* token;
   char* tokenbuf;
   char* pushedtokens[DEC_MAX_PUSHEDTOKENS];
   int npushedtokens;
   int linenumber;
   int linepos;
   int nblocks;
   int blocknr;               /**< number of the currentblock between 0 and Nblocks-1*/
   DECSECTION section;
   SCIP_Bool haserror;
};
typedef struct DecInput DECINPUT;

/** data for gp reader */
struct SCIP_ReaderData
{
   DECDECOMP* decdecomp;
   SCIP_HASHMAP *vartoindex;
   int* varstoblock;          /**< index=var id // value= -1 or blockID or -2 for multipleblocks*/
   int* nblockvars;           /**< n variable per block that are no linkingvars*/
   SCIP_CONS*** blockcons;    /**< [decnr][consid]  */
   int* nblockcons;
   SCIP_HASHMAP* constoblock;
   int nlinkingcons;
   int nlinkingvars;
};

static const int NOVALUE = - 1;
static const int LINKINGVALUE = - 2;
static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=";
static const char commentchars[] = "\\";




/*
 * Local methods (for reading)
 */

/** issues an error message and marks the DEC data to have errors */
static
void syntaxError(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput,  /**< DEC reading data */
   const char* msg      /**< error message */
   )
{
   char formatstr[256];

   assert(decinput != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error in line %d: %s ('%s')\n",
           decinput->linenumber, msg, decinput->token);
   if( decinput->linebuf[strlen(decinput->linebuf) - 1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", decinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", decinput->linebuf);
   }
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", decinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, formatstr, "^");
   decinput->section = DEC_END;
   decinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   assert(decinput != NULL);

   return decinput->haserror;
}

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char c   /**< input character */
   )
{
   return (c == '\0') || (strchr(delimchars, c) != NULL);
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char c /**< input character */
   )
{
   return (strchr(tokenchars, c) != NULL);
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char c,              /**< input character */
   char nextc,          /**< next input character */
   SCIP_Bool firstchar, /**< is the given character the first char of the token? */
   SCIP_Bool* hasdot,   /**< pointer to update the dot flag */
   DECEXPTYPE* exptype  /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit(c) )
      return TRUE;
   else if( (*exptype == DEC_EXP_NONE) && ! (*hasdot) && (c == '.') )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( ! firstchar && (*exptype == DEC_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = DEC_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit(nextc) )
      {
         *exptype = DEC_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == DEC_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = DEC_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   int i;

   assert(decinput != NULL);

   /* clear the line */
   BMSclearMemoryArray(decinput->linebuf, DEC_MAX_LINELEN);

   /* read next line */
   decinput->linepos = 0;
   decinput->linebuf[DEC_MAX_LINELEN - 2] = '\0';
   if( SCIPfgets(decinput->linebuf, sizeof (decinput->linebuf), decinput->file) == NULL )
      return FALSE;
   decinput->linenumber ++;
   if( decinput->linebuf[DEC_MAX_LINELEN - 2] != '\0' )
   {
      SCIPerrorMessage("Error: line %d exceeds %d characters\n", decinput->linenumber, DEC_MAX_LINELEN - 2);
      decinput->haserror = TRUE;
      return FALSE;
   }
   decinput->linebuf[DEC_MAX_LINELEN - 1] = '\0';
   decinput->linebuf[DEC_MAX_LINELEN - 2] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++ i )
   {
      char* commentstart;

      commentstart = strchr(decinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart + 1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
      }
   }

   return TRUE;
}

/** swaps the addresses of two pointers */
static
void swapPointers(
   char** pointer1,     /**< first pointer */
   char** pointer2      /**< second pointer */
   )
{
   char* tmp;

   tmp = * pointer1;
   *pointer1 = * pointer2;
   *pointer2 = tmp;
}

/** reads the next token from the input file into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   SCIP_Bool hasdot;
   DECEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(decinput != NULL);
   assert(decinput->linepos < DEC_MAX_LINELEN);

   /* check the token stack */
   if( decinput->npushedtokens > 0 )
   {
      swapPointers(&decinput->token, &decinput->pushedtokens[decinput->npushedtokens - 1]);
      decinput->npushedtokens --;
      SCIPdebugMessage("(line %d) read token again: '%s'\n", decinput->linenumber, decinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = decinput->linebuf;
   while( isDelimChar(buf[decinput->linepos]) )
   {
      if( buf[decinput->linepos] == '\0' )
      {
         if( ! getNextLine(decinput) )
         {
            decinput->section = DEC_END;
            SCIPdebugMessage("(line %d) end of file\n", decinput->linenumber);
            return FALSE;
         }
         assert(decinput->linepos == 0);
      }
      else
         decinput->linepos ++;
   }
   assert(decinput->linepos < DEC_MAX_LINELEN);
   assert(! isDelimChar(buf[decinput->linepos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = DEC_EXP_NONE;
   if( isValueChar(buf[decinput->linepos], buf[decinput->linepos + 1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < DEC_MAX_LINELEN);
         assert(! isDelimChar(buf[decinput->linepos]));
         decinput->token[tokenlen] = buf[decinput->linepos];
         tokenlen ++;
         decinput->linepos ++;
      }
      while( isValueChar(buf[decinput->linepos], buf[decinput->linepos + 1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < DEC_MAX_LINELEN);
         decinput->token[tokenlen] = buf[decinput->linepos];
         tokenlen ++;
         decinput->linepos ++;
         if( tokenlen == 1 && isTokenChar(decinput->token[0]) )
            break;
      }
      while( ! isDelimChar(buf[decinput->linepos]) && ! isTokenChar(buf[decinput->linepos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1
              && (decinput->token[tokenlen - 1] == '<' || decinput->token[tokenlen - 1] == '>' || decinput->token[tokenlen - 1] == '=')
              && buf[decinput->linepos] == '=' )
      {
         decinput->linepos ++;
      }
      else if( decinput->token[tokenlen - 1] == '=' && (buf[decinput->linepos] == '<' || buf[decinput->linepos] == '>') )
      {
         decinput->token[tokenlen - 1] = buf[decinput->linepos];
         decinput->linepos ++;
      }
   }
   assert(tokenlen < DEC_MAX_LINELEN);
   decinput->token[tokenlen] = '\0';

   SCIPdebugMessage("(line %d) read token: '%s'\n", decinput->linenumber, decinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   assert(decinput != NULL);
   assert(decinput->npushedtokens < DEC_MAX_PUSHEDTOKENS);

   swapPointers(&decinput->pushedtokens[decinput->npushedtokens], &decinput->token);
   decinput->npushedtokens ++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   assert(decinput != NULL);

   swapPointers(&decinput->token, &decinput->tokenbuf);
}

/** returns whether the current token is a value */
static
SCIP_Bool isInt(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput,  /**< DEC reading data */
   int* value           /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(decinput != NULL);
   assert(value != NULL);

   if( strcasecmp(decinput->token, "INFINITY") == 0 || strcasecmp(decinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      long val;
      char* endptr;

      val = strtol(decinput->token, &endptr, 0);
      if( endptr != decinput->token && * endptr == '\0' )
      {
         *value = val;
         return TRUE;
      }
   }

   return FALSE;
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   SCIP_Bool iscolon;

   assert(decinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(decinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(decinput) )
   {
      iscolon = (strcmp(decinput->token, ":") == 0);
      pushToken(decinput);
   }

   /* reinstall the previous token by swapping back the token buffer */
   swapTokenBuffer(decinput);

   /* check for ':' */
   if( iscolon )
      return FALSE;

   if( strcasecmp(decinput->token, "NBLOCKS") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: NBLOCKS\n", decinput->linenumber);
      decinput->section = DEC_NBLOCKS;
      return TRUE;
   }

   if( strcasecmp(decinput->token, "BLOCK") == 0 )
   {
      int blocknr;

      decinput->section = DEC_BLOCK;

      if( getNextToken(decinput) )
      {
         /* read block number */
         if( isInt(scip, decinput, &blocknr) )
         {
            assert(blocknr >= 0);
            assert(blocknr <= decinput->nblocks);

            decinput->blocknr = blocknr - 1;
         }
         else
            syntaxError(scip, decinput, "no block number after block keyword!\n");
      }
      else
         syntaxError(scip, decinput, "no block number after block keyword!\n");

      SCIPdebugMessage("new section: BLOCK %d\n", decinput->blocknr);

      return TRUE;

   }

   if( strcasecmp(decinput->token, "MASTERCONSS") == 0 )
   {
      decinput->section = DEC_MASTERCONSS;

      SCIPdebugMessage("new section: MASTERCONSS\n");

      return TRUE;
   }

   if( strcasecmp(decinput->token, "END") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: END\n", decinput->linenumber);
      decinput->section = DEC_END;
      return TRUE;
   }

   return FALSE;
}

/** reads the header of the file */
static
SCIP_RETCODE readStart(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   assert(decinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( ! getNextToken(decinput) )
         return SCIP_OKAY;
   }
   while( ! isNewSection(scip, decinput) );

   return SCIP_OKAY;
}

/** reads the nblocks section */
static
SCIP_RETCODE readNBlocks(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   int nblocks;

   assert(decinput != NULL);

   while( getNextToken(decinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, decinput) )
         return SCIP_OKAY;

      /* read number of blocks */
      if( isInt(scip, decinput, &nblocks) )
      {
         /* assert(nblocks > 0); */


         if( decinput->nblocks == NOVALUE )
         {
            decinput->nblocks = nblocks;
            /*            GCGrelaxSetNPricingprobs(scip, nblocks);*/
         }
         else
            syntaxError(scip, decinput, "2 integer values in nblocks section");
         SCIPdebugMessage("Number of blocks = %d\n", decinput->nblocks);
      }
   }

   return SCIP_OKAY;
}

/** reads the blocks section */
static
SCIP_RETCODE readBlock(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   int varvalue;
   int varvalueindex;
   int consindex;
   int i;
   int blockid;
   int nvars;
   SCIP_CONS* cons;
   SCIP_VAR** vars;

   assert(decinput != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   while( getNextToken(decinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, decinput) )
         break;

      /* the token must be the name of an existing cons */
      cons = SCIPfindCons(scip, decinput->token);
      if( cons == NULL )
      {
         syntaxError(scip, decinput, "unknown variable in block section");
         break;
      }

      /* get all vars for  the specific constraint */
      vars = SCIPgetVarsXXX(scip, cons);
      nvars = SCIPgetNVarsXXX(scip, cons);
      for( i = 0; i < nvars; i ++ )
      {
         /*
          * saving for dedecomp
          */
         readerdata->nblockvars[decinput->blocknr] = (readerdata->nblockvars[decinput->blocknr]) + 1;

         /* mark each var if it is in none, one or more blocks */
         varvalueindex = SCIPvarGetProbindex(vars[i]);
         varvalue = readerdata->varstoblock[varvalueindex ];
         if( varvalue == NOVALUE )
         {
            readerdata->varstoblock[varvalueindex ] = decinput->blocknr;
         }
         else if( varvalue != NOVALUE || varvalue != LINKINGVALUE )
         {
            readerdata->varstoblock[varvalueindex ] = LINKINGVALUE;
            /* decrease the value again if it is a linking var */
            readerdata->nblockvars[decinput->blocknr] = readerdata->nblockvars[decinput->blocknr] - 1;
            readerdata->nlinkingvars = readerdata->nlinkingvars + 1;
         }
      }
      /*
       * saving block <-> constraint
       */
      /** @todo check if linking constraints are no in the subscipcons */
      blockid = decinput->blocknr;
      consindex = readerdata->nblockcons[blockid];
      readerdata->blockcons[decinput->blocknr][consindex] = cons;
      readerdata->nblockcons[decinput->blocknr] = readerdata->nblockcons[decinput->blocknr] + 1;


      if( SCIPhashmapGetImage(readerdata->constoblock, cons) != (void*) (size_t) LINKINGVALUE )
      {
         SCIP_CALL( SCIPhashmapSetImage(readerdata->constoblock, cons, (void*) (size_t) LINKINGVALUE) );
         readerdata->nlinkingcons --;
      }
      else
      {
         SCIP_CALL( SCIPhashmapSetImage(readerdata->constoblock, cons, (void*) ((size_t) decinput->blocknr)) );
      }
      SCIPfreeMemoryArray(scip, &vars);
   }

   return SCIP_OKAY;
}

/** reads the masterconss section */
static
SCIP_RETCODE readMasterconss(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput   /**< DEC reading data */
   )
{
   assert(decinput != NULL);

   while( getNextToken(decinput) )
   {
      SCIP_CONS* cons;

      /* check if we reached a new section */
      if( isNewSection(scip, decinput) )
         break;

      cons = NULL;

      /* the token must be the name of an existing constraint */
      cons = SCIPfindCons(scip, decinput->token);
      if( cons == NULL )
      {
         syntaxError(scip, decinput, "unknown constraint in masterconss section");
         break;
      }
      else
      {
         /* set the block number of the variable to the number of the current block */
         /*SCIP_CALL( GCGrelaxMarkConsMaster(scip, cons) );*/
      }
   }

   return SCIP_OKAY;
}

/** fills the whole Decomp struct after the dec file has been read */
static
SCIP_RETCODE fillDecompStruct(
   SCIP* scip,                   /**< SCIP data structure */
   DECINPUT* decinput,           /**< DEC reading data */
   SCIP_READERDATA* readerdata   /**< reader data*/
   )
{

   DECDECOMP* decomp;
   SCIP_VAR** allvars;
   SCIP_CONS** allcons;
   int i;
   int j;
   int n;
   int n2;
   int value;
   int ind;
   int nconss;

   assert(scip != NULL);
   assert(decinput != NULL);
   assert(readerdata != NULL);
   assert(readerdata->decdecomp != NULL);
   decomp = readerdata->decdecomp;

   /*
    * alloc
    */

   decomp->nblocks = decinput->nblocks;
   decomp->type = DEC_DECTYPE_ARROWHEAD;

   nconss = SCIPgetNConss(scip);
   /* nvars */
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->nsubscipvars, decinput->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->nsubscipconss, decinput->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->subscipvars, decinput->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->subscipconss, decinput->nblocks) );

   for( i = 0; i < decinput->nblocks; ++i )
   {
      decomp->nsubscipvars[i] = 0;
      decomp->nsubscipconss[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->subscipvars[i], readerdata->nblockvars[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->subscipconss[i], readerdata->nblockcons[i]) );
   }

   /* linking vars */
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->linkingvars, readerdata->nlinkingvars) );
   decomp->nlinkingvars = 0;

   /* linking conss */
   SCIP_CALL( SCIPallocMemoryArray(scip, &decomp->linkingconss, readerdata->nlinkingcons) );
   decomp->nlinkingconss = 0;

   /* hashmaps */
   SCIP_CALL( SCIPhashmapCreate(&decomp->constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&decomp->vartoblock, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* init */
   allvars = SCIPgetVars(scip);
   allcons = SCIPgetConss(scip);

   /*
    * insert values to decomp struct
    */
   /* vars to blocks */
   n = SCIPgetNVars(scip);
   for( i = 0; i < n; i ++ )
   {
      value = readerdata->varstoblock[i];
      if( value == NOVALUE )
      {
         /** @todo What should be done  in this case? */
      }
      else if( value == LINKINGVALUE )
      {
         ind = decomp->nlinkingvars;
         decomp->linkingvars[ind] = allvars[i];
         decomp->nlinkingvars = decomp->nlinkingvars + 1;
         /* hashmap */
         SCIP_CALL( SCIPhashmapInsert(decomp->vartoblock, allvars[i], (void*) (size_t) LINKINGVALUE) );
      }
      else
      { /* value = block id=index of block */
         assert(value >= 0);
         assert(value <= decinput->nblocks);
         /** @todo check if "i" is the id of the variable */
         ind = decomp->nsubscipvars[value];

         assert(ind >= 0);
         assert(ind <= readerdata->nblockvars[value]);

         decomp->subscipvars[value][ind] = allvars[i];
         decomp->nsubscipvars[value] ++;

         /* hashmap */
         SCIP_CALL( SCIPhashmapInsert(decomp->vartoblock, allvars[i], (void*) (size_t) value) );
      }
   }

   ind = 0;
   for( i = 0; i < nconss; i ++ )
   {
      if( SCIPhashmapGetImage(readerdata->constoblock, allcons[i]) == (void*) (size_t) LINKINGVALUE )
      {
         decomp->linkingconss[ind] = allcons[i];
         decomp->nlinkingconss ++;
         ind ++;
      }
   }

   for( i = 0; i < nconss; i ++ )
   {
      SCIP_CALL( SCIPhashmapInsert(decomp->constoblock, allcons[i], (void*) (size_t) LINKINGVALUE) );
   }
   n = decinput->nblocks;
   for( i = 0; i < n; i ++ )
   {
      n2 = readerdata->nblockcons[i];
      for( j = 0; j < n2; j ++ )
      {
         ind = decomp->nsubscipconss[i];
         decomp->subscipconss[i][ind] = readerdata->blockcons[i][j];
         decomp->nsubscipconss[i] ++;
         /* hashmap */
         /* @todo besser machen? */
         SCIP_CALL( SCIPhashmapRemove(decomp->constoblock, readerdata->blockcons[i][j]) );
         SCIP_CALL( SCIPhashmapInsert(decomp->constoblock, readerdata->blockcons[i][j], (void*) (size_t) i) );
      }
   }

   return SCIP_OKAY;
}

/** reads an DEC file */
static
SCIP_RETCODE
readDECFile(
   SCIP* scip,          /**< SCIP data structure */
   DECINPUT* decinput,  /**< DEC reading data */
   const char* filename /**< name of the input file */
   )
{
   int i;
   int nconss;
   int nblocksread;
   int nvars;
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   SCIP_CONS** conss;
   nblocksread = FALSE;
   assert(decinput != NULL);

   /*   SCIP_CALL( GCGcreateOrigVarsData(scip) );*/

   /* open file */
   decinput->file = SCIPfopen(filename, "r");
   if( decinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /** @todo check */
   assert(scip != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   readerdata->nlinkingcons = SCIPgetNConss(scip);
   readerdata->nlinkingvars = 0;
   nvars = SCIPgetNVars(scip);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   /* alloc: var -> block mapping */

   SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->varstoblock, nvars) );

   for( i = 0; i < nvars; i ++ )
   {
      readerdata->varstoblock[i] = NOVALUE;
   }

  SCIP_CALL( SCIPhashmapCreate(&readerdata->constoblock, SCIPblkmem(scip), nconss) );

   for( i = 0; i < SCIPgetNConss(scip); i ++ )
   {
      SCIP_CALL( SCIPhashmapInsert(readerdata->constoblock, conss[i], (void*) (size_t) LINKINGVALUE) );
   }

   /* parse the file */
   decinput->section = DEC_START;
   while( decinput->section != DEC_END && ! hasError(decinput) )
   {
      switch( decinput->section )
      {
         case DEC_START:
            SCIP_CALL( readStart(scip, decinput) );
            break;

         case DEC_NBLOCKS:
            SCIP_CALL( readNBlocks(scip, decinput) );
            break;

         case DEC_BLOCK:
            if( nblocksread == FALSE )
            {
               /* alloc n vars per block */
               SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->nblockvars, decinput->nblocks) );
               SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->nblockcons, decinput->nblocks) );
               SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->blockcons, decinput->nblocks) );
               for( i = 0; i < decinput->nblocks; i ++ )
               {
                  readerdata->nblockvars[i] = 0;
                  readerdata->nblockcons[i] = 0;
                  SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->blockcons[i], nconss) );
               }
               nblocksread = TRUE;
            }
            SCIP_CALL( readBlock(scip, decinput) );
            break;

         case DEC_MASTERCONSS:
            SCIP_CALL( readMasterconss(scip, decinput) );
            break;

         case DEC_END: /* this is already handled in the while() loop */
         default:
            SCIPerrorMessage("invalid DEC file section <%d>\n", decinput->section);
            return SCIP_INVALIDDATA;
      }
   }
   /* fill decomp */
   SCIP_CALL( fillDecompStruct(scip, decinput, readerdata) );

   /* add decomp to cons_decomp */
   SCIP_CALL( SCIPconshdlrDecompAddDecdecomp(scip, readerdata->decdecomp) );


   for( i = 0; i < decinput->nblocks; i ++ )
   {
      SCIPfreeBufferArray(scip, &readerdata->blockcons[i]);
   }
   SCIPfreeBufferArray(scip, &readerdata->blockcons);
   SCIPfreeBufferArray(scip, &readerdata->nblockcons);
   SCIPfreeBufferArray(scip, &readerdata->varstoblock);
   SCIPfreeBufferArray(scip, &readerdata->nblockvars);
   SCIPhashmapFree(&readerdata->constoblock);

   /* close file */
   SCIPfclose(decinput->file);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeDec)
{
   SCIP_READERDATA* readerdata;
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   if( readerdata->decdecomp->type == DEC_DECTYPE_UNKNOWN)
      DECdecdecompFree(scip, &readerdata->decdecomp);
   SCIPfreeMemory(scip, &readerdata);
   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadDec)
{
   SCIP_CALL( SCIPreadDec(scip, reader, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteDec)
{   /*lint --e{715}*/
   SCIP_READERDATA* readerdata;
   assert(scip != NULL);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPwriteDecomp(scip, file, DECgetBestDecomp(scip), TRUE) );
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the dec file reader in SCIP */
SCIP_RETCODE
SCIPincludeReaderDec(
        SCIP * scip /**< SCIP data structure */
        )
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   /* create dec reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );
   SCIP_CALL( DECdecdecompCreate(scip, &readerdata->decdecomp) );

   /* include lp reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeDec, readerReadDec, readerWriteDec, readerdata));

   return SCIP_OKAY;
}

/* reads problem from file */
SCIP_RETCODE SCIPreadDec(
   SCIP* scip,             /**< SCIP data structure */
   SCIP_READER* reader,    /**< the file reader itself */
   const char* filename,   /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT * result    /**< pointer to store the result of the file reading call */
   )
{
   DECINPUT decinput;
   int i;

   /* initialize DEC input data */
   decinput.file = NULL;
   decinput.linebuf[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &decinput.token, DEC_MAX_LINELEN) );
   decinput.token[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &decinput.tokenbuf, DEC_MAX_LINELEN) );
   decinput.tokenbuf[0] = '\0';
   for( i = 0; i < DEC_MAX_PUSHEDTOKENS; ++ i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &decinput.pushedtokens[i], DEC_MAX_LINELEN) );
   }

   decinput.npushedtokens = 0;
   decinput.linenumber = 0;
   decinput.linepos = 0;
   decinput.section = DEC_START;
   decinput.nblocks = NOVALUE;
   decinput.blocknr = - 2;
   decinput.haserror = FALSE;

   /* read the file */
   SCIP_CALL( readDECFile(scip, &decinput, filename) );

   /* free dynamically allocated memory */
   SCIPfreeMemoryArray(scip, &decinput.token);
   SCIPfreeMemoryArray(scip, &decinput.tokenbuf);
   for( i = 0; i < DEC_MAX_PUSHEDTOKENS; ++ i )
   {
      SCIPfreeMemoryArray(scip, &decinput.pushedtokens[i]);
   }

   /* evaluate the result */
   if( decinput.haserror )
      return SCIP_READERROR;
   else
   {
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}

/** write the data optionally using the decomposition data */
static
SCIP_RETCODE writeData(
   SCIP* scip,          /**< SCIP data structure */
   FILE* file,          /**< File pointer to write to */
   DECDECOMP* decdecomp /**< Decomposition pointer */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(file != NULL);
   assert(decdecomp != NULL);

   assert(decdecomp->type == DEC_DECTYPE_ARROWHEAD
           || decdecomp->type == DEC_DECTYPE_BORDERED
           || decdecomp->type == DEC_DECTYPE_DIAGONAL
           || decdecomp->type == DEC_DECTYPE_UNKNOWN
           || decdecomp->type == DEC_DECTYPE_STAIRCASE);
   SCIPdebugMessage("DECDECOMP Type: %d\n",decdecomp->type);

   /* if we don't have staicase, but something else, go through the blocks and create the indices */
   /* conss */
   assert(decdecomp->nsubscipconss != NULL);
   assert(decdecomp->subscipconss != NULL);

   /* linking cons */
   assert(decdecomp->nlinkingconss >= 0 && decdecomp->nlinkingconss < SCIPgetNConss(scip));
   assert(decdecomp->linkingconss != NULL || decdecomp->nlinkingconss == 0 );

   SCIPinfoMessage(scip, file, "NBLOCKS\n");
   SCIPinfoMessage(scip, file, "%d\n", decdecomp->nblocks);

   for( i = 0; i < decdecomp->nblocks; i ++ )
   {
      SCIPinfoMessage(scip, file, "BLOCK %d\n", i + 1);
      for( j = 0; j < decdecomp->nsubscipconss[i]; j ++ )
      {
         SCIPinfoMessage(scip, file, "%s\n", SCIPconsGetName(decdecomp->subscipconss[i][j]));
      }
   }

   if( decdecomp->nlinkingconss > 0 )
   {
      SCIPinfoMessage(scip, file, "MASTERCONSS\n");
      for( i = 0; i < decdecomp->nlinkingconss; i ++ )
      {
         SCIPinfoMessage(scip, file, "%s\n", SCIPconsGetName(decdecomp->linkingconss[i]));
      }
   }

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPwriteDecomp(
   SCIP* scip,                   /**< SCIP data structure */
   FILE* file,                   /**< File pointer to write to */
   DECDECOMP* decdecomp,         /**< Decomposition pointer */
   SCIP_Bool writeDecomposition  /**< whether to write decomposed problem */
   )
{
   char outname[SCIP_MAXSTRLEN];
   assert(scip != NULL);
   assert(file != NULL);

   if( writeDecomposition )
   {
      if( decdecomp == NULL )
      {
         SCIPwarningMessage("Cannot write decomposed problem if decomposition structure empty!");
         writeDecomposition = FALSE;
         /* return SCIP_INVALIDDATA; */
      }
   }

   /* print header */
   if( decdecomp == NULL )
   {
      SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   }
   else
   {
      SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d", SCIPgetProbName(scip), decdecomp->nblocks);
   }

   if( writeDecomposition )
   {

      /* write data */
      SCIP_CALL( writeData(scip, file, decdecomp) );

   }


   return SCIP_OKAY;
}

