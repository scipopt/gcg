/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   reader_lpb.c
 * @ingroup FILEReaders 
 * @brief  LPB file reader
 * @author Gerald Gamrath
 *
 * @todo Test for uniqueness of variable names (after cutting down).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h>
#endif
#include <ctype.h>

#include "scip/reader_lpb.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/pub_misc.h"

#define READER_NAME             "lpbreader"
#define READER_DESC             "file reader for MIPs in ILOG's LPB file format"
#define READER_EXTENSION        "lpb"


/*
 * Data structures
 */
#define LPB_MAX_LINELEN       65536
#define LPB_MAX_PUSHEDTOKENS  2
#define LPB_INIT_COEFSSIZE    8192
#define LPB_MAX_PRINTLEN      561       /**< the maximum length of any line is 560 + '\\0' = 561*/
#define LPB_MAX_NAMELEN       256       /**< the maximum length for any name is 255 + '\\0' = 256 */
#define LPB_PRINTLEN          100

/** Section in LPB File */
enum LpbSection
{
   LPB_START, LPB_OBJECTIVE, LPB_CONSTRAINTS, LPB_BOUNDS, LPB_GENERALS, LPB_BINARIES, LPB_SEMICONTINUOUS, LPB_SOS, LPB_END
};
typedef enum LpbSection LPBSECTION;

enum LpbExpType
{
   LPB_EXP_NONE, LPB_EXP_UNSIGNED, LPB_EXP_SIGNED
};
typedef enum LpbExpType LPBEXPTYPE;

enum LpbSense
{
   LPB_SENSE_NOTHING, LPB_SENSE_LE, LPB_SENSE_GE, LPB_SENSE_EQ
};
typedef enum LpbSense LPBSENSE;

/** LPB reading data */
struct LpbInput
{
   SCIP_FILE*           file;
   char                 linebuf[LPB_MAX_LINELEN];
   char                 probname[LPB_MAX_LINELEN];
   char                 objname[LPB_MAX_LINELEN];
   char*                token;
   char*                tokenbuf;
   char*                pushedtokens[LPB_MAX_PUSHEDTOKENS];
   int                  npushedtokens;
   int                  linenumber;
   int                  linepos;
   LPBSECTION           section;
   SCIP_OBJSENSE        objsense;
   SCIP_Bool            inlazyconstraints;
   SCIP_Bool            inusercuts;
   SCIP_Bool            haserror;
};
typedef struct LpbInput LPBINPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=";
static const char commentchars[] = "\\";




/*
 * Local methods (for reading)
 */

/** issues an error message and marks the LPB data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput,            /**< LPB reading data */
   const char*           msg                 /**< error message */
   )
{
   char formatstr[256];

   assert(lpbinput != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error in line %d: %s ('%s')\n",
      lpbinput->linenumber, msg, lpbinput->token);
   if( lpbinput->linebuf[strlen(lpbinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", lpbinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", lpbinput->linebuf);
   }
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", lpbinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, formatstr, "^");
   lpbinput->section  = LPB_END;
   lpbinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   return lpbinput->haserror;
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
   char                  c,                  /**< input character */
   char                  nextc,              /**< next input character */
   SCIP_Bool             firstchar,          /**< is the given character the first char of the token? */
   SCIP_Bool*            hasdot,             /**< pointer to update the dot flag */
   LPBEXPTYPE*            exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit(c) )
      return TRUE;
   else if( (*exptype == LPB_EXP_NONE) && !(*hasdot) && (c == '.') )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == LPB_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = LPB_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit(nextc) )
      {
         *exptype = LPB_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == LPB_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = LPB_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   int i;

   assert(lpbinput != NULL);

   /* clear the line */
   BMSclearMemoryArray(lpbinput->linebuf, LPB_MAX_LINELEN);

   /* read next line */
   lpbinput->linepos = 0;
   lpbinput->linebuf[LPB_MAX_LINELEN-2] = '\0';
   if( SCIPfgets(lpbinput->linebuf, sizeof(lpbinput->linebuf), lpbinput->file) == NULL )
      return FALSE;
   lpbinput->linenumber++;
   if( lpbinput->linebuf[LPB_MAX_LINELEN-2] != '\0' )
   {
      SCIPerrorMessage("Error: line %d exceeds %d characters\n", lpbinput->linenumber, LPB_MAX_LINELEN-2);
      lpbinput->haserror = TRUE;
      return FALSE;
   }
   lpbinput->linebuf[LPB_MAX_LINELEN-1] = '\0';
   lpbinput->linebuf[LPB_MAX_LINELEN-2] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(lpbinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
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

   tmp = *pointer1;
   *pointer1 = *pointer2;
   *pointer2 = tmp;
}

/** reads the next token from the input file into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   SCIP_Bool hasdot;
   LPBEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(lpbinput != NULL);
   assert(lpbinput->linepos < LPB_MAX_LINELEN);

   /* check the token stack */
   if( lpbinput->npushedtokens > 0 )
   {
      swapPointers(&lpbinput->token, &lpbinput->pushedtokens[lpbinput->npushedtokens-1]);
      lpbinput->npushedtokens--;
      SCIPdebugMessage("(line %d) read token again: '%s'\n", lpbinput->linenumber, lpbinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = lpbinput->linebuf;
   while( isDelimChar(buf[lpbinput->linepos]) )
   {
      if( buf[lpbinput->linepos] == '\0' )
      {
         if( !getNextLine(lpbinput) )
         {
            lpbinput->section = LPB_END;
            SCIPdebugMessage("(line %d) end of file\n", lpbinput->linenumber);
            return FALSE;
         }
         assert(lpbinput->linepos == 0);
      }
      else
         lpbinput->linepos++;
   }
   assert(lpbinput->linepos < LPB_MAX_LINELEN);
   assert(!isDelimChar(buf[lpbinput->linepos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = LPB_EXP_NONE;
   if( isValueChar(buf[lpbinput->linepos], buf[lpbinput->linepos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < LPB_MAX_LINELEN);
         assert(!isDelimChar(buf[lpbinput->linepos]));
         lpbinput->token[tokenlen] = buf[lpbinput->linepos];
         tokenlen++;
         lpbinput->linepos++;
      }
      while( isValueChar(buf[lpbinput->linepos], buf[lpbinput->linepos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < LPB_MAX_LINELEN);
         lpbinput->token[tokenlen] = buf[lpbinput->linepos];
         tokenlen++;
         lpbinput->linepos++;
         if( tokenlen == 1 && isTokenChar(lpbinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[lpbinput->linepos]) && !isTokenChar(buf[lpbinput->linepos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1
         && (lpbinput->token[tokenlen-1] == '<' || lpbinput->token[tokenlen-1] == '>' || lpbinput->token[tokenlen-1] == '=')
         && buf[lpbinput->linepos] == '=' )
      {
         lpbinput->linepos++;
      }
      else if( lpbinput->token[tokenlen-1] == '=' && (buf[lpbinput->linepos] == '<' || buf[lpbinput->linepos] == '>') )
      {
         lpbinput->token[tokenlen-1] = buf[lpbinput->linepos];
         lpbinput->linepos++;
      }
   }
   assert(tokenlen < LPB_MAX_LINELEN);
   lpbinput->token[tokenlen] = '\0';

   SCIPdebugMessage("(line %d) read token: '%s'\n", lpbinput->linenumber, lpbinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);
   assert(lpbinput->npushedtokens < LPB_MAX_PUSHEDTOKENS);

   swapPointers(&lpbinput->pushedtokens[lpbinput->npushedtokens], &lpbinput->token);
   lpbinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);
   assert(lpbinput->npushedtokens < LPB_MAX_PUSHEDTOKENS);

   swapPointers(&lpbinput->pushedtokens[lpbinput->npushedtokens], &lpbinput->tokenbuf);
   lpbinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   swapPointers(&lpbinput->token, &lpbinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   SCIP_Bool iscolon;

   assert(lpbinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(lpbinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(lpbinput) )
   {
      iscolon = (strcmp(lpbinput->token, ":") == 0);
      pushToken(lpbinput);
   }

   /* reinstall the previous token by swapping back the token buffer */
   swapTokenBuffer(lpbinput);

   /* check for ':' */
   if( iscolon )
      return FALSE;

   if( strcasecmp(lpbinput->token, "MINIMIZE") == 0
      || strcasecmp(lpbinput->token, "MINIMUM") == 0
      || strcasecmp(lpbinput->token, "MIN") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: OBJECTIVE\n", lpbinput->linenumber);
      lpbinput->section = LPB_OBJECTIVE;
      lpbinput->objsense = SCIP_OBJSENSE_MINIMIZE;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "MAXIMIZE") == 0
      || strcasecmp(lpbinput->token, "MAXIMUM") == 0
      || strcasecmp(lpbinput->token, "MAX") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: OBJECTIVE\n", lpbinput->linenumber);
      lpbinput->section = LPB_OBJECTIVE;
      lpbinput->objsense = SCIP_OBJSENSE_MAXIMIZE;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "SUBJECT") == 0 )
   {
      /* check if the next token is 'TO' */
      swapTokenBuffer(lpbinput);
      if( getNextToken(lpbinput) )
      {
         if( strcasecmp(lpbinput->token, "TO") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", lpbinput->linenumber);
            lpbinput->section = LPB_CONSTRAINTS;
            lpbinput->inlazyconstraints = FALSE;
            lpbinput->inusercuts = FALSE;
            return TRUE;
         }
         else
            pushToken(lpbinput);
      }
      swapTokenBuffer(lpbinput);
   }

   if( strcasecmp(lpbinput->token, "SUCH") == 0 )
   {
      /* check if the next token is 'THAT' */
      swapTokenBuffer(lpbinput);
      if( getNextToken(lpbinput) )
      {
         if( strcasecmp(lpbinput->token, "THAT") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", lpbinput->linenumber);
            lpbinput->section = LPB_CONSTRAINTS;
            lpbinput->inlazyconstraints = FALSE;
            lpbinput->inusercuts = FALSE;
            return TRUE;
         }
         else
            pushToken(lpbinput);
      }
      swapTokenBuffer(lpbinput);
   }

   if( strcasecmp(lpbinput->token, "st") == 0
      || strcasecmp(lpbinput->token, "S.T.") == 0
      || strcasecmp(lpbinput->token, "ST.") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: CONSTRAINTS\n", lpbinput->linenumber);
      lpbinput->section = LPB_CONSTRAINTS;
      lpbinput->inlazyconstraints = FALSE;
      lpbinput->inusercuts = FALSE;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "LAZY") == 0 )
   {
      /* check if the next token is 'CONSTRAINTS' */
      swapTokenBuffer(lpbinput);
      if( getNextToken(lpbinput) )
      {
         if( strcasecmp(lpbinput->token, "CONSTRAINTS") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS (lazy)\n", lpbinput->linenumber);
            lpbinput->section = LPB_CONSTRAINTS;
            lpbinput->inlazyconstraints = TRUE;
            lpbinput->inusercuts = FALSE;
            return TRUE;
         }
         else
            pushToken(lpbinput);
      }
      swapTokenBuffer(lpbinput);
   }

   if( strcasecmp(lpbinput->token, "USER") == 0 )
   {
      /* check if the next token is 'CUTS' */
      swapTokenBuffer(lpbinput);
      if( getNextToken(lpbinput) )
      {
         if( strcasecmp(lpbinput->token, "CUTS") == 0 )
         {
            SCIPdebugMessage("(line %d) new section: CONSTRAINTS (user cuts)\n", lpbinput->linenumber);
            lpbinput->section = LPB_CONSTRAINTS;
            lpbinput->inlazyconstraints = FALSE;
            lpbinput->inusercuts = TRUE;
            return TRUE;
         }
         else
            pushToken(lpbinput);
      }
      swapTokenBuffer(lpbinput);
   }

   if( strcasecmp(lpbinput->token, "BOUNDS") == 0
      || strcasecmp(lpbinput->token, "BOUND") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: BOUNDS\n", lpbinput->linenumber);
      lpbinput->section = LPB_BOUNDS;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "GENERAL") == 0
      || strcasecmp(lpbinput->token, "GENERALS") == 0
      || strcasecmp(lpbinput->token, "GEN") == 0
      || strcasecmp(lpbinput->token, "INTEGER") == 0
      || strcasecmp(lpbinput->token, "INTEGERS") == 0
      || strcasecmp(lpbinput->token, "INT") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: GENERALS\n", lpbinput->linenumber);
      lpbinput->section = LPB_GENERALS;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "BINARY") == 0
      || strcasecmp(lpbinput->token, "BINARIES") == 0
      || strcasecmp(lpbinput->token, "BIN") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: BINARIES\n", lpbinput->linenumber);
      lpbinput->section = LPB_BINARIES;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "SEMI-CONTINUOUS") == 0
      || strcasecmp(lpbinput->token, "SEMIS") == 0
      || strcasecmp(lpbinput->token, "SEMI") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: SEMICONTINUOUS\n", lpbinput->linenumber);
      lpbinput->section = LPB_SEMICONTINUOUS;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "SOS") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: SOS\n", lpbinput->linenumber);
      lpbinput->section = LPB_SOS;
      return TRUE;
   }

   if( strcasecmp(lpbinput->token, "END") == 0 )
   {
      SCIPdebugMessage("(line %d) new section: END\n", lpbinput->linenumber);
      lpbinput->section = LPB_END;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   LPBINPUT*              lpbinput,            /**< LPB reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(lpbinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( lpbinput->token[1] == '\0' )
   {
      if( *lpbinput->token == '+' )
         return TRUE;
      else if( *lpbinput->token == '-' )
      {
         *sign *= -1;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isValue(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput,            /**< LPB reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(lpbinput != NULL);
   assert(value != NULL);

   if( strcasecmp(lpbinput->token, "INFINITY") == 0 || strcasecmp(lpbinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(lpbinput->token, &endptr);
      if( endptr != lpbinput->token && *endptr == '\0' )
      {
         *value = val;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is an equation sense */
static
SCIP_Bool isSense(
   LPBINPUT*              lpbinput,            /**< LPB reading data */
   LPBSENSE*              sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(lpbinput != NULL);

   if( strcmp(lpbinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = LPB_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(lpbinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = LPB_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(lpbinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = LPB_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns the variable with the given name, or creates a new variable if it does not exist */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 name,               /**< name of the variable */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Bool*            created             /**< pointer to store whether a new variable was created, or NULL */
   )
{
   assert(name != NULL);
   assert(var != NULL);

   *var = SCIPfindVar(scip, name);
   if( *var == NULL )
   {
      SCIP_VAR* newvar;
      SCIP_Bool dynamiccols;
      SCIP_Bool initial;
      SCIP_Bool removable;

      SCIP_CALL( SCIPgetBoolpbaram(scip, "reading/lpbreader/dynamiccols", &dynamiccols) );
      initial = !dynamiccols;
      removable = dynamiccols;

      /* create new variable of the given name */
      SCIPdebugMessage("creating new variable: <%s>\n", name);
      SCIP_CALL( GCGcreateVar(scip, &newvar, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
            initial, removable, NULL, NULL, NULL, NULL) );
      SCIP_CALL( GCGaddOriginalVar(scip, newvar) );
      *var = newvar;

      if( created != NULL )
         *created = TRUE;
   }
   else if( created != NULL )
      *created = FALSE;

   return SCIP_OKAY;
}

/** reads the header of the file */
static
SCIP_RETCODE readStart(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( !getNextToken(lpbinput) )
         return SCIP_OKAY;
   }
   while( !isNewSection(lpbinput) );

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput,            /**< LPB reading data */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   LPB_MAX_LINELEN */
   SCIP_VAR***           vars,               /**< pointer to store the array with variables (must be freed by caller) */
   SCIP_Real**           coefs,              /**< pointer to store the array with coefficients (must be freed by caller) */
   int*                  ncoefs,             /**< pointer to store the number of coefficients */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;
   int coefssize;

   assert(lpbinput != NULL);
   assert(name != NULL);
   assert(vars != NULL);
   assert(coefs != NULL);
   assert(ncoefs != NULL);
   assert(newsection != NULL);

   *vars = NULL;
   *coefs = NULL;
   *name = '\0';
   *ncoefs = 0;
   *newsection = FALSE;

   /* read the first token, which may be the name of the line */
   if( getNextToken(lpbinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(lpbinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* remember the token in the token buffer */
      swapTokenBuffer(lpbinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(lpbinput) )
      {
         if( strcmp(lpbinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the line name */
            strncpy(name, lpbinput->tokenbuf, LPB_MAX_LINELEN);
            name[LPB_MAX_LINELEN - 1] = '\0';
            SCIPdebugMessage("(line %d) read constraint name: '%s'\n", lpbinput->linenumber, name);
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            pushToken(lpbinput);
            pushBufferToken(lpbinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(lpbinput);
      }
   }

   /* initialize buffers for storing the coefficients */
   coefssize = LPB_INIT_COEFSSIZE;
   SCIP_CALL( SCIPallocMemoryArray(scip, vars, coefssize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, coefs, coefssize) );

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   *ncoefs = 0;
   while( getNextToken(lpbinput) )
   {
      SCIP_VAR* var;

      /* check if we reached a new section */
      if( isNewSection(lpbinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* check if we reached an equation sense */
      if( isSense(lpbinput, NULL) )
      {
         /* put the sense back onto the token stack */
         pushToken(lpbinput);
         break;
      }

      /* check if we read a sign */
      if( isSign(lpbinput, &coefsign) )
      {
         SCIPdebugMessage("(line %d) read coefficient sign: %+d\n", lpbinput->linenumber, coefsign);
         havesign = TRUE;
         continue;
      }

      /* all but the first coefficient need a sign */
      if( *ncoefs > 0 && !havesign )
      {
         syntaxError(scip, lpbinput, "expected sign ('+' or '-') or sense ('<' or '>')");
         return SCIP_OKAY;
      }

      /* check if we read a value */
      if( isValue(scip, lpbinput, &coef) )
      {
         SCIPdebugMessage("(line %d) read coefficient value: %g with sign %+d\n", lpbinput->linenumber, coef, coefsign);
         if( havevalue )
         {
            syntaxError(scip, lpbinput, "two consecutive values");
            return SCIP_OKAY;
         }
         havevalue = TRUE;
         continue;
      }

      /* the token is a variable name: get the corresponding variable (or create a new one) */
      SCIP_CALL( getVariable(scip, lpbinput->token, &var, NULL) );

      /* insert the coefficient */
      SCIPdebugMessage("(line %d) read coefficient: %+g<%s>\n", lpbinput->linenumber, coefsign * coef, SCIPvarGetName(var));
      if( !SCIPisZero(scip, coef) )
      {
         /* resize the vars and coefs array if needed */
         if( *ncoefs >= coefssize )
         {
            coefssize *= 2;
            coefssize = MAX(coefssize, (*ncoefs)+1);
            SCIP_CALL( SCIPreallocMemoryArray(scip, vars, coefssize) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, coefs, coefssize) );
         }
         assert(*ncoefs < coefssize);

         /* add coefficient */
         (*vars)[*ncoefs] = var;
         (*coefs)[*ncoefs] = coefsign * coef;
         (*ncoefs)++;
      }

      /* reset the flags and coefficient value for the next coefficient */
      coefsign = +1;
      coef = 1.0;
      havesign = FALSE;
      havevalue = FALSE;
   }

   return SCIP_OKAY;
}

/** reads the objective section */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*             lpbinput            /**< LPB reading data */
   )
{
   char name[LPB_MAX_LINELEN];
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int ncoefs;
   SCIP_Bool newsection;

   assert(lpbinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, lpbinput, name, &vars, &coefs, &ncoefs, &newsection) );
   if( !hasError(lpbinput) )
   {
      int i;

      /* set the objective values */
      for( i = 0; i < ncoefs; ++i )
      {
         SCIP_CALL( SCIPchgVarObj(GCGprobGetOrigprob(scip), vars[i], coefs[i]) );
      }
   }

   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &vars);
   SCIPfreeMemoryArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** reads the constraints section */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   char name[LPB_MAX_LINELEN];
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Bool newsection;
   LPBSENSE sense;
   SCIP_Real sidevalue;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamicrows;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   int ncoefs;
   int sidesign;

   assert(lpbinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, lpbinput, name, &vars, &coefs, &ncoefs, &newsection) );
   if( hasError(lpbinput) )
      goto TERMINATE;
   if( newsection )
   {
      if( ncoefs > 0 )
         syntaxError(scip, lpbinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if( !getNextToken(lpbinput) || !isSense(lpbinput, &sense) )
   {
      syntaxError(scip, lpbinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if( !getNextToken(lpbinput) )
   {
      syntaxError(scip, lpbinput, "missing right hand side");
      goto TERMINATE;
   }
   if( isSign(lpbinput, &sidesign) )
   {
      if( !getNextToken(lpbinput) )
      {
         syntaxError(scip, lpbinput, "missing value of right hand side");
         goto TERMINATE;
      }
   }
   if( !isValue(scip, lpbinput, &sidevalue) )
   {
      syntaxError(scip, lpbinput, "expected value as right hand side");
      goto TERMINATE;
   }
   sidevalue *= sidesign;

   /* assign the left and right hand side, depending on the constraint sense */
   switch( sense )
   {
   case LPB_SENSE_GE:
      lhs = sidevalue;
      rhs = SCIPinfinity(scip);
      break;
   case LPB_SENSE_LE:
      lhs = -SCIPinfinity(scip);
      rhs = sidevalue;
      break;
   case LPB_SENSE_EQ:
      lhs = sidevalue;
      rhs = sidevalue;
      break;
   case LPB_SENSE_NOTHING:
   default:
      SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
      return SCIP_INVALIDDATA;
   }

   /* create and add the linear constraint */
   SCIP_CALL( SCIPgetBoolpbaram(scip, "reading/lpbreader/dynamicconss", &dynamicconss) );
   SCIP_CALL( SCIPgetBoolpbaram(scip, "reading/lpbreader/dynamicrows", &dynamicrows) );
   initial = !dynamicrows && !lpbinput->inlazyconstraints && !lpbinput->inusercuts;
   separate = TRUE;
   enforce = !lpbinput->inusercuts;
   check = !lpbinput->inusercuts;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = dynamicconss;
   removable = dynamicrows || lpbinput->inusercuts;
   SCIP_CALL( GCGcreateConsLinear(scip, name, ncoefs, vars, coefs, lhs, rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE, -1) );

 TERMINATE:
   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &vars);
   SCIPfreeMemoryArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** reads the bounds section */
static
SCIP_RETCODE readBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   while( getNextToken(lpbinput) )
   {
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Real lb;
      SCIP_Real ub;
      int sign;
      SCIP_Bool hassign;
      LPBSENSE leftsense;

      /* check if we reached a new section */
      if( isNewSection(lpbinput) )
         return SCIP_OKAY;

      /* default bounds are [0,+inf] */
      lb = 0.0;
      ub = SCIPinfinity(scip);
      leftsense = LPB_SENSE_NOTHING;

      /* check if the first token is a sign */
      sign = +1;
      hassign = isSign(lpbinput, &sign);
      if( hassign && !getNextToken(lpbinput) )
      {
         syntaxError(scip, lpbinput, "expected value");
         return SCIP_OKAY;
      }

      /* the first token must be either a value or a variable name */
      if( isValue(scip, lpbinput, &value) )
      {
         /* first token is a value: the second token must be a sense */
         if( !getNextToken(lpbinput) || !isSense(lpbinput, &leftsense) )
         {
            syntaxError(scip, lpbinput, "expected bound sense '<=', '=', or '>='");
            return SCIP_OKAY;
         }

         /* update the bound corresponding to the sense */
         switch( leftsense )
         {
         case LPB_SENSE_GE:
            ub = sign * value;
            break;
         case LPB_SENSE_LE:
            lb = sign * value;
            break;
         case LPB_SENSE_EQ:
            lb = sign * value;
            ub = sign * value;
            break;
         case LPB_SENSE_NOTHING:
         default:
            SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
            return SCIP_INVALIDDATA;
         }
      }
      else if( hassign )
      {
         syntaxError(scip, lpbinput, "expected value");
         return SCIP_OKAY;
      }
      else
         pushToken(lpbinput);

      /* the next token must be a variable name */
      if( !getNextToken(lpbinput) )
      {
         syntaxError(scip, lpbinput, "expected variable name");
         return SCIP_OKAY;
      }
      SCIP_CALL( getVariable(scip, lpbinput->token, &var, NULL) );

      /* the next token might be another sense, or the word "free" */
      if( getNextToken(lpbinput) )
      {
         LPBSENSE rightsense;

         if( isSense(lpbinput, &rightsense) )
         {
            /* check, if the senses fit */
            if( leftsense == LPB_SENSE_NOTHING
               || (leftsense == LPB_SENSE_LE && rightsense == LPB_SENSE_LE)
               || (leftsense == LPB_SENSE_GE && rightsense == LPB_SENSE_GE) )
            {
               if( !getNextToken(lpbinput) )
               {
                  syntaxError(scip, lpbinput, "expected value or sign");
                  return SCIP_OKAY;
               }

               /* check if the next token is a sign */
               sign = +1;
               hassign = isSign(lpbinput, &sign);
               if( hassign && !getNextToken(lpbinput) )
               {
                  syntaxError(scip, lpbinput, "expected value");
                  return SCIP_OKAY;
               }

               /* the next token must be a value */
               if( !isValue(scip, lpbinput, &value) )
               {
                  syntaxError(scip, lpbinput, "expected value");
                  return SCIP_OKAY;
               }

               /* update the bound corresponding to the sense */
               switch( rightsense )
               {
               case LPB_SENSE_GE:
                  lb = sign * value;
                  break;
               case LPB_SENSE_LE:
                  ub = sign * value;
                  break;
               case LPB_SENSE_EQ:
                  lb = sign * value;
                  ub = sign * value;
                  break;
               case LPB_SENSE_NOTHING:
               default:
                  SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
                  return SCIP_INVALIDDATA;
               }
            }
            else
            {
               syntaxError(scip, lpbinput, "the two bound senses do not fit");
               return SCIP_OKAY;
            }
         }
         else if( strcasecmp(lpbinput->token, "FREE") == 0 )
         {
            if( leftsense != LPB_SENSE_NOTHING )
            {
               syntaxError(scip, lpbinput, "variable with bound is marked as 'free'");
               return SCIP_OKAY;
            }
            lb = -SCIPinfinity(scip);
            ub = SCIPinfinity(scip);
         }
         else
         {
            /* the token was no sense: push it back to the token stack */
            pushToken(lpbinput);
         }
      }

      /* change the bounds of the variable if bounds have been given (do not destroy earlier specification of bounds) */
      if ( lb != 0.0 )
	 SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      if ( ub != SCIPinfinity(scip) )
	 SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      SCIPdebugMessage("(line %d) new bounds: <%s>[%g,%g]\n", lpbinput->linenumber, SCIPvarGetName(var),
	 SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

/** reads the generals section */
static
SCIP_RETCODE readGenerals(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   while( getNextToken(lpbinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;

      /* check if we reached a new section */
      if( isNewSection(lpbinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpbinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpbinput, "unknown variable in generals section");
         return SCIP_OKAY;
      }

      /* mark the variable to be integral */
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER) );
   }

   return SCIP_OKAY;
}

/** reads the binaries section */
static
SCIP_RETCODE readBinaries(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   while( getNextToken(lpbinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;

      /* check if we reached a new section */
      if( isNewSection(lpbinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, lpbinput->token, &var, &created) );
      if( created )
      {
         syntaxError(scip, lpbinput, "unknown variable in binaries section");
         return SCIP_OKAY;
      }

      /* mark the variable to be binary and change its bounds appropriately */
      if( SCIPvarGetLbGlobal(var) < 0.0 )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );
      }
      if( SCIPvarGetUbGlobal(var) > 1.0 )
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 1.0) );
      }
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY) );
   }

   return SCIP_OKAY;
}

/** reads the semicontinuous section */
static
SCIP_RETCODE readSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput             /**< LPB reading data */
   )
{
   assert(lpbinput != NULL);

   while( getNextToken(lpbinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(lpbinput) )
         return SCIP_OKAY;

      /* semi-continuous variables are not yet supported by SCIP */
      syntaxError(scip, lpbinput, "semi-continuous variables not yet supported by SCIP");
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** reads an LPB file */
static
SCIP_RETCODE readLPBFile(
   SCIP*                 scip,               /**< SCIP data structure */
   LPBINPUT*              lpbinput,            /**< LPB reading data */
   const char*           filename            /**< name of the input file */
   )
{
   assert(lpbinput != NULL);

   /* open file */
   lpbinput->file = SCIPfopen(filename, "r");
   if( lpbinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProbGCG(scip, filename, npricingprobs) );

   /* parse the file */
   lpbinput->section = LPB_START;
   while( lpbinput->section != LPB_END && !hasError(lpbinput) )
   {
      switch( lpbinput->section )
      {
      case LPB_START:
         SCIP_CALL( readStart(scip, lpbinput) );
         break;

      case LPB_OBJECTIVE:
         SCIP_CALL( readObjective(scip, lpbinput) );
         break;

      case LPB_CONSTRAINTS:
         SCIP_CALL( readConstraints(scip, lpbinput) );
         break;

      case LPB_BOUNDS:
         SCIP_CALL( readBounds(scip, lpbinput) );
         break;

      case LPB_GENERALS:
         SCIP_CALL( readGenerals(scip, lpbinput) );
         break;

      case LPB_BINARIES:
         SCIP_CALL( readBinaries(scip, lpbinput) );
         break;

      case LPB_SEMICONTINUOUS:
         SCIP_CALL( readSemicontinuous(scip, lpbinput) );
         break;

      case LPB_END: /* this is already handled in the while() loop */
      default:
         SCIPerrorMessage("invalid LPB file section <%d>\n", lpbinput->section);
         return SCIP_INVALIDDATA;
      }
   }

   /* close file */
   SCIPfclose(lpbinput->file);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeLpb NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadLpb)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadLpb(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteLpb NULL

/*
 * reader specific interface methods
 */

/** includes the lpb file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderLpb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create lpb reader data */
   readerdata = NULL;

   /* include lp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeLpb, readerReadLpb, readerWriteLpb, readerdata) );

   /* add lpb reader parameters */
   SCIP_CALL( SCIPaddBoolparam(scip,
         "reading/lpbreader/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/lpbreader/dynamiccols", "should columns be added and removed dynamically to the LPB?",
         NULL, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/lpbreader/dynamicrows", "should rows be added and removed dynamically to the LP?",
         NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}


/* reads problem from file */
SCIP_RETCODE SCIPreadLpb(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_READER*       reader,             /**< the file reader itself */
   const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   LPBINPUT lpbinput;
   int i;

   /* initialize LPB input data */
   lpbinput.file = NULL;
   lpbinput.linebuf[0] = '\0';
   lpbinput.probname[0] = '\0';
   lpbinput.objname[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpbinput.token, LPB_MAX_LINELEN) );
   lpbinput.token[0] = '\0';
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpbinput.tokenbuf, LPB_MAX_LINELEN) );
   lpbinput.tokenbuf[0] = '\0';
   for( i = 0; i < LPB_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &lpbinput.pushedtokens[i], LPB_MAX_LINELEN) );
   }

   lpbinput.npushedtokens = 0;
   lpbinput.linenumber = 0;
   lpbinput.linepos = 0;
   lpbinput.section = LPB_START;
   lpbinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   lpbinput.inlazyconstraints = FALSE;
   lpbinput.inusercuts = FALSE;
   lpbinput.haserror = FALSE;

   /* read the file */
   SCIP_CALL( readLPBFile(scip, &lpbinput, filename) );

   /* free dynamically allocated memory */
   SCIPfreeMemoryArray(scip, &lpbinput.token);
   SCIPfreeMemoryArray(scip, &lpbinput.tokenbuf);
   for( i = 0; i < LPB_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeMemoryArray(scip, &lpbinput.pushedtokens[i]);
   }

   /* evaluate the result */
   if( lpbinput.haserror )
      return SCIP_PARSEERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, lpbinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}
