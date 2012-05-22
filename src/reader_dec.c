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
 * @author Lukas Kirchhart
 * @author Martin Bergner
 * @author Gerald Gamrath
 * This reader reads in a dec-file that defines the structur to be used for the decomposition.
 * The structure is defined constraint-wise, i.e., the number of blocks and the constraints belonging
 * to each block are  defined.  If needed, constraints can also be  forced into the master, even if
 * they could be transferred to one block.
 *
 * The keywords are:
 * - NBlocks: to be followed by a line giving the number of blocks
 * - Block i with 1 <= i <= nblocks: to be followed by the names of the constraints belonging to block i,
                  one per line.
 * - Masterconss: to be followed by names of constraints, one per line, that should go into the master,
 *                even if they only contain variables of one block and could thus be added to this block.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

/* @todo really needed? #include <stdlib.h> */
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strcasecmp() */
#endif
#include <ctype.h>

#include "reader_dec.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"

#include "cons_decomp.h"
#include "pub_decomp.h"

#define READER_NAME             "decreader"
#define READER_DESC             "file reader for blocks in dec format"
#define READER_EXTENSION        "dec"


/*
 * Data structures
 */
#define DEC_MAX_LINELEN       65536
#define DEC_MAX_PUSHEDTOKENS  2

/** section in DEC File */
enum DecSection
{
   DEC_START, DEC_NBLOCKS, DEC_BLOCK, DEC_MASTERCONSS, DEC_END
};
typedef enum DecSection DECSECTION;

/** exponent indicator of the a value */
enum DecExpType
{
   DEC_EXP_NONE, DEC_EXP_UNSIGNED, DEC_EXP_SIGNED
};
typedef enum DecExpType DECEXPTYPE;

/** DEC reading data */
struct DecInput
{
   SCIP_FILE* file;                          /**< file to read */
   char linebuf[DEC_MAX_LINELEN];            /**< line buffer */
   char* token;                              /**< current token */
   char* tokenbuf;                           /**< token buffer */
   char* pushedtokens[DEC_MAX_PUSHEDTOKENS]; /**< token stack */
   int npushedtokens;                        /**< size of token buffer */
   int linenumber;                           /**< current line number */
   int linepos;                              /**< current line position (column) */
   int nblocks;                              /**< number of blocks */
   int blocknr;                              /**< number of the currentblock between 0 and Nblocks-1*/
   DECSECTION section;                       /**< current section */
   SCIP_Bool haserror;                       /**< flag to indicate an error occurence */
};
typedef struct DecInput DECINPUT;

/** data for dec reader */
struct SCIP_ReaderData
{
   DECDECOMP* decdecomp;      /**< decomposition data structure*/
   int* varstoblock;          /**< index=var id // value= -1 or blockID or -2 for multipleblocks*/
   int* nblockvars;           /**< n variable per block that are not linkingvars*/
   SCIP_CONS*** blockconss;   /**< array of assignments from constraints to their blocks [blocknr][consid]  */
   int* nblockconss;          /**< number of block-constraints for blockID*/
   SCIP_HASHMAP* constoblock; /**< hashmap key=constaint value=block*/
   int nlinkingconss;         /**< number of linking constraints*/
   int nlinkingvars;          /**< number of linking vars*/
};
static const int NOVALUE = -1;
static const int LINKINGVALUE = -2;
static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=";
static const char commentchars[] = "\\";




/*
 * Local methods (for reading)
 */

/** issues an error message and marks the DEC data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput,           /**< DEC reading data */
   const char*           msg                 /**< error message */
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
   DECINPUT*             decinput            /**< DEC reading data */
   )
{
   assert(decinput != NULL);

   return decinput->haserror;
}

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char                  c                   /**< input character */
   )
{
   return (c == '\0') || (strchr(delimchars, c) != NULL);
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char                  c                   /**< input character */
   )
{
   return (strchr(tokenchars, c) != NULL);
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char                 c,                   /**< input character */
   char                 nextc,               /**< next input character */
   SCIP_Bool            firstchar,           /**< is the given character the first char of the token? */
   SCIP_Bool*           hasdot,              /**< pointer to update the dot flag */
   DECEXPTYPE*          exptype              /**< pointer to update the exponent type */
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
   DECINPUT*             decinput            /**< DEC reading data */
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
   char**                pointer1,           /**< first pointer */
   char**                pointer2            /**< second pointer */
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
   DECINPUT*             decinput            /**< DEC reading data */
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
   if( isValueChar(buf[decinput->linepos], buf[decinput->linepos + 1], TRUE, &hasdot, &exptype) ) /*lint !e679*/
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < DEC_MAX_LINELEN);
         assert(! isDelimChar(buf[decinput->linepos]));
         decinput->token[tokenlen] = buf[decinput->linepos];
         ++tokenlen;
         ++(decinput->linepos);
         assert(decinput->linepos < DEC_MAX_LINELEN-1);
      }
      while( isValueChar(buf[decinput->linepos], buf[decinput->linepos + 1], FALSE, &hasdot, &exptype) ); /*lint !e679*/
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
   DECINPUT*             decinput            /**< DEC reading data */
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
   DECINPUT*             decinput            /**< DEC reading data */
   )
{
   assert(decinput != NULL);

   swapPointers(&decinput->token, &decinput->tokenbuf);
}

/** returns whether the current token is a value */
static
SCIP_Bool isInt(
   SCIP*                 scip,                /**< SCIP data structure */
   DECINPUT*             decinput,            /**< DEC reading data */
   int*                  value                /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   long val;
   char* endptr;

   assert(decinput != NULL);
   assert(value != NULL);
   assert(!(strcasecmp(decinput->token, "INFINITY") == 0) && !(strcasecmp(decinput->token, "INF") == 0));

   val = strtol(decinput->token, &endptr, 0);
   if( endptr != decinput->token && * endptr == '\0' )
   {
      if( val < INT_MIN || val > INT_MAX ) /*lint !e685*/
         return FALSE;

      *value = (int) val;
      return TRUE;
   }

   return FALSE;
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput            /**< DEC reading data */
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
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput            /**< DEC reading data */
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
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput            /**< DEC reading data */
   )
{
   int nblocks;

   assert(scip != NULL);
   assert(decinput != NULL);

   while( getNextToken(decinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, decinput) )
      {
         if( decinput->nblocks == NOVALUE )
            syntaxError(scip, decinput, "no integer value in nblocks section");
         else
            return SCIP_OKAY;
      }

      /* read number of blocks */
      if( isInt(scip, decinput, &nblocks) )
      {
         if( decinput->nblocks == NOVALUE )
            decinput->nblocks = nblocks;
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
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput,           /**< DEC reading data */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   int oldblock;
   int varidx;
   int consindex;
   int i;
   int blockid;
   int nvars;
   SCIP_CONS* cons;
   SCIP_VAR** vars;

   assert(decinput != NULL);
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
         syntaxError(scip, decinput, "unknown constraint in block section");
         break;
      }

      /* get all vars for the specific constraint */
      nvars = SCIPgetNVarsXXX(scip, cons);
      vars = NULL;
      if( nvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, vars, nvars) );
      }

      blockid = decinput->blocknr;

      for( i = 0; i < nvars; i ++ )
      {
         assert(vars != NULL); /* for flexelint */

         /* store for each var whether it is in none, one or more blocks */
         varidx = SCIPvarGetProbindex(vars[i]);
         oldblock = readerdata->varstoblock[varidx];

         /* variable was assigned to no block before, just assign it to the new block */
         if( oldblock == NOVALUE )
         {
            SCIPdebugMessage("\tVar %s temporary in block %d.\n", SCIPvarGetName(vars[i]), blockid);
            readerdata->varstoblock[varidx] = blockid;
            ++(readerdata->nblockvars[blockid]);
         }
         /* variable was assigned to another (non-linking) block before, so it becomes a linking variable, now */
         else if( (oldblock != LINKINGVALUE) && oldblock != blockid )
         {
            SCIPdebugMessage("\tVar %s is linking (old %d != %d new).\n", SCIPvarGetName(vars[i]), oldblock, blockid);
            assert(oldblock != blockid);

            readerdata->varstoblock[varidx] = LINKINGVALUE;

            /* decrease the value again if it is a linking var */
            --(readerdata->nblockvars[oldblock]);
            ++(readerdata->nlinkingvars);
         }
      }

      /*
       * saving block <-> constraint
       */

      /** @todo check if linking constraints are not in the subscipcons */
      consindex = readerdata->nblockconss[blockid];
      readerdata->blockconss[blockid][consindex] = cons;
      ++(readerdata->nblockconss[blockid]);

      assert(SCIPhashmapGetImage(readerdata->constoblock, cons) == (void*)(size_t) LINKINGVALUE);

      SCIPdebugMessage("cons %s is in block %d\n", SCIPconsGetName(cons), blockid);
      SCIP_CALL( SCIPhashmapSetImage(readerdata->constoblock, cons, (void*) ((size_t) blockid)) );
      --(readerdata->nlinkingconss);

      SCIPfreeBufferArrayNull(scip, &vars);
   }

   return SCIP_OKAY;
}

/** reads the masterconss section */
static
SCIP_RETCODE readMasterconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput,           /**< DEC reading data */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   assert(scip != NULL);
   assert(decinput != NULL);
   assert(readerdata != NULL);

   while( getNextToken(decinput) )
   {
      SCIP_CONS* cons;

      /* check if we reached a new section */
      if( isNewSection(scip, decinput) )
         break;

      /* the token must be the name of an existing constraint */
      cons = SCIPfindCons(scip, decinput->token);
      if( cons == NULL )
      {
         syntaxError(scip, decinput, "unknown constraint in masterconss section");
         break;
      }
      else
      {
         assert(SCIPhashmapGetImage(readerdata->constoblock, cons) == (void*)(size_t) LINKINGVALUE);

         SCIPdebugMessage("cons %s is linking constraint\n", decinput->token);
      }
   }

   return SCIP_OKAY;
}

/** fills the whole Decomp struct after the dec file has been read */
static
SCIP_RETCODE fillDecompStruct(
   SCIP*                 scip,               /**< SCIP data structure */
   DECINPUT*             decinput,           /**< DEC reading data */
   SCIP_READERDATA*      readerdata          /**< reader data*/
   )
{
   DECDECOMP* decomp;
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;
   SCIP_VAR** allvars;
   SCIP_CONS** allcons;
   SCIP_CONS*** subscipconss;
   SCIP_CONS** linkingconss;
   SCIP_VAR*** subscipvars;
   SCIP_VAR** linkingvars;
   int* nsubscipconss;
   int* nsubscipvars;
   int nlinkingconss;
   int nlinkingvars;
   int i;
   int j;
   int nvars;
   int nblockconss;
   int blocknr;
   int ind;
   int nconss;
   int nblocks;

   assert(scip != NULL);
   assert(decinput != NULL);
   assert(readerdata != NULL);
   assert(readerdata->decdecomp != NULL);
   decomp = readerdata->decdecomp;

   printf("decomp = %p\n", decomp);

   allvars = SCIPgetVars(scip);
   allcons = SCIPgetConss(scip);
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   nblocks = decinput->nblocks;

   DECdecdecompSetNBlocks(decomp, nblocks);
   DECdecdecompSetType(decomp, DEC_DECTYPE_ARROWHEAD);

   /* get memory for subscip variables and constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );

   for( i = 0; i < nblocks; ++i )
   {
      nsubscipvars[i] = 0;
      nsubscipconss[i] = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], readerdata->nblockvars[i]) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], readerdata->nblockconss[i]) ); /*lint !e866*/
   }

   /* get memory for linking variables and constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, readerdata->nlinkingvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, readerdata->nlinkingconss) );
   nlinkingvars = 0;
   nlinkingconss = 0;

   /* assign variables to blocks or as linking variables according to the varstoblock structure */
   for( i = 0; i < nvars; i ++ )
   {
      SCIPdebugMessage("var %s ", SCIPvarGetName(allvars[i]));
      blocknr = readerdata->varstoblock[i];

      if( blocknr == NOVALUE )
      {
         SCIPdebugMessage("is unknown\n" );
         /** @todo What should be done in this case? gg: copy directly into master */
      }
      else if( blocknr == LINKINGVALUE )
      {
         /* add variable to array of linking variables */
         ind = nlinkingvars;
         linkingvars[ind] = allvars[i];
         ++nlinkingvars;

         SCIPdebugMessage("is linking\n" );
      }
      else
      {
         assert(blocknr >= 0);
         assert(blocknr <= nblocks);
         assert(SCIPvarGetProbindex(allvars[i]) == i);

         /* get current number of variables in the block */
         ind = nsubscipvars[blocknr];
         assert(ind >= 0);
         assert(ind <= readerdata->nblockvars[blocknr]);

         /* add variable to array of variables in the block */
         subscipvars[blocknr][ind] = allvars[i];
         nsubscipvars[blocknr] ++;

         SCIPdebugMessage("is in block %d\n", blocknr);
      }
   }

   /* set subscip and linking variables in decomposition structure */
   SCIP_CALL( DECdecdecompSetSubscipvars(scip, decomp, subscipvars, nsubscipvars) );
   SCIP_CALL( DECdecdecompSetLinkingvars(scip, decomp, linkingvars, nlinkingvars) );

   /* copy linking constraints and set them in decomposition data */
   for( i = 0; i < nconss; i ++ )
   {
      if( SCIPhashmapGetImage(readerdata->constoblock, allcons[i]) == (void*) (size_t) LINKINGVALUE )
      {
         linkingconss[nlinkingconss] = allcons[i];
         ++nlinkingconss;
         SCIPdebugMessage("cons %s is linking\n", SCIPconsGetName(allcons[i]));
      }
   }
   SCIP_CALL( DECdecdecompSetLinkingconss(scip, decomp, linkingconss, nlinkingconss) );

   /* hashmaps */
   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nvars) );

   for( i = 0; i < nconss; i ++ )
   {
      SCIP_CALL( SCIPhashmapInsert(constoblock, allcons[i], (void*) (size_t) LINKINGVALUE) );
   }
   for( i = 0; i < nblocks; i ++ )
   {
      nblockconss = readerdata->nblockconss[i];
      for( j = 0; j < nblockconss; j ++ )
      {
         ind = nsubscipconss[i];
         subscipconss[i][ind] = readerdata->blockconss[i][j];
         ++nsubscipconss[i];

         /* hashmap */
         SCIPdebugMessage("cons %s is in block %d\n", SCIPconsGetName(readerdata->blockconss[i][j]), i);
         SCIP_CALL( SCIPhashmapSetImage(constoblock, readerdata->blockconss[i][j], (void*) (size_t) i) );
      }
   }
   SCIP_CALL( DECdecdecompSetSubscipconss(scip, decomp, subscipconss, nsubscipconss) );
   DECdecdecompSetConstoblock(decomp, constoblock);
   DECdecdecompSetVartoblock(decomp, vartoblock);

   SCIPfreeBufferArray(scip, &linkingconss);
   SCIPfreeBufferArray(scip, &linkingvars);
   for( i = nblocks - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &subscipconss[i]);
      SCIPfreeBufferArray(scip, &subscipvars[i]);
   }
   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);

   return SCIP_OKAY;
}

/** reads a DEC file */
static
SCIP_RETCODE readDECFile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< Reader data structure */
   DECINPUT*             decinput,           /**< DEC reading data */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_CONS** conss;
   int nconss;
   int nblocksread;
   int nvars;
   int i;

   assert(decinput != NULL);
   assert(scip != NULL);
   assert(reader != NULL);

   nblocksread = FALSE;

   /* open file */
   decinput->file = SCIPfopen(filename, "r");
   if( decinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   readerdata->nlinkingconss = SCIPgetNConss(scip);
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

   /* cons -> block mapping */
   SCIP_CALL( SCIPhashmapCreate(&readerdata->constoblock, SCIPblkmem(scip), nconss) );
   for( i = 0; i < nconss; i ++ )
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
               SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->nblockconss, decinput->nblocks) );
               SCIP_CALL( SCIPallocBufferArray(scip, &readerdata->blockconss, decinput->nblocks) );
               for( i = 0; i < decinput->nblocks; i ++ )
               {
                  readerdata->nblockvars[i] = 0;
                  readerdata->nblockconss[i] = 0;
                  SCIP_CALL( SCIPallocBufferArray(scip, &(readerdata->blockconss[i]), nconss) ); /*lint !e866*/
               }
               nblocksread = TRUE;
            }
            SCIP_CALL( readBlock(scip, decinput, readerdata) );
            break;

         case DEC_MASTERCONSS:
            SCIP_CALL( readMasterconss(scip, decinput, readerdata) );
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

   for( i = decinput->nblocks - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &readerdata->blockconss[i]);
   }
   SCIPfreeBufferArray(scip, &readerdata->blockconss);
   SCIPfreeBufferArray(scip, &readerdata->nblockconss);
   SCIPfreeBufferArray(scip, &readerdata->nblockvars);
   SCIPfreeBufferArray(scip, &readerdata->varstoblock);
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

   /* free decomp structure and readerdata */
   if( DECdecdecompGetType(readerdata->decdecomp) == DEC_DECTYPE_UNKNOWN )
      DECdecdecompFree(scip, &readerdata->decdecomp);
   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadDec)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPreadDec(scip, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteDec)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);

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
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create dec reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );
   SCIP_CALL( DECdecdecompCreate(scip, &readerdata->decdecomp) );

   /* include dec reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeDec, readerReadDec, readerWriteDec, readerdata));

   return SCIP_OKAY;
}

/* reads problem from file */
SCIP_RETCODE SCIPreadDec(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   SCIP_READER* reader;
   DECINPUT decinput;
   int i;

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   /* initialize DEC input data */
   decinput.file = NULL;
   decinput.linebuf[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &decinput.token, DEC_MAX_LINELEN) ); /*lint !e506*/
   decinput.token[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &decinput.tokenbuf, DEC_MAX_LINELEN) ); /*lint !e506*/
   decinput.tokenbuf[0] = '\0';
   for( i = 0; i < DEC_MAX_PUSHEDTOKENS; ++ i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &decinput.pushedtokens[i], DEC_MAX_LINELEN) ); /*lint !e506 !e866*/
   }

   decinput.npushedtokens = 0;
   decinput.linenumber = 0;
   decinput.linepos = 0;
   decinput.section = DEC_START;
   decinput.nblocks = NOVALUE;
   decinput.blocknr = - 2;
   decinput.haserror = FALSE;

   /* read the file */
   SCIP_CALL( readDECFile(scip, reader, &decinput, filename) );

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
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DECDECOMP*            decdecomp           /**< Decomposition pointer */
   )
{
   SCIP_CONS*** subscipconss;
   SCIP_CONS** linkingconss;
   int* nsubscipconss;
   int nlinkingconss;
   int nblocks;
   int i;
   int j;

   assert(scip != NULL);
   assert(file != NULL);
   assert(decdecomp != NULL);

   assert(DECdecdecompGetType(decdecomp) == DEC_DECTYPE_ARROWHEAD
           || DECdecdecompGetType(decdecomp) == DEC_DECTYPE_BORDERED
           || DECdecdecompGetType(decdecomp) == DEC_DECTYPE_DIAGONAL
           || DECdecdecompGetType(decdecomp) == DEC_DECTYPE_UNKNOWN
           || DECdecdecompGetType(decdecomp) == DEC_DECTYPE_STAIRCASE);
   SCIPdebugMessage("DECDECOMP Type: %s\n", DECgetStrType(DECdecdecompGetType(decdecomp)));

   /* if we don't have staicase, but something else, go through the blocks and create the indices */
   /* subscip conss */
   subscipconss = DECdecdecompGetSubscipconss(decdecomp);
   nsubscipconss = DECdecdecompGetNSubscipconss(decdecomp);
   assert(subscipconss != NULL);
   assert(nsubscipconss != NULL);

   /* linking cons */
   linkingconss = DECdecdecompGetLinkingconss(decdecomp);
   nlinkingconss = DECdecdecompGetNLinkingconss(decdecomp);
   assert(nlinkingconss >= 0 && nlinkingconss < SCIPgetNConss(scip));
   assert(linkingconss != NULL || nlinkingconss == 0 );

   nblocks = DECdecdecompGetNBlocks(decdecomp);

   SCIPinfoMessage(scip, file, "NBLOCKS\n");
   SCIPinfoMessage(scip, file, "%d\n", nblocks);

   for( i = 0; i < nblocks; i ++ )
   {
      SCIPinfoMessage(scip, file, "BLOCK %d\n", i + 1);
      for( j = 0; j < nsubscipconss[i]; j ++ )
      {
         SCIPinfoMessage(scip, file, "%s\n", SCIPconsGetName(subscipconss[i][j]));
      }
   }

   if( nlinkingconss > 0 )
   {
      assert(linkingconss != NULL); /* for flexelint */
      SCIPinfoMessage(scip, file, "MASTERCONSS\n");
      for( i = 0; i < nlinkingconss; i ++ )
      {
         SCIPinfoMessage(scip, file, "%s\n", SCIPconsGetName(linkingconss[i]));
      }
   }

   return SCIP_OKAY;
}

/** write a DEC file for a given decomposition */
SCIP_RETCODE SCIPwriteDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DECDECOMP*            decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition  /**< whether to write decomposed problem */
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
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   }
   else
   {
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d", SCIPgetProbName(scip), DECdecdecompGetNBlocks(decdecomp));
   }

   if( writeDecomposition )
   {
      /* write data */
      SCIP_CALL( writeData(scip, file, decdecomp) );
   }

   return SCIP_OKAY;
}

