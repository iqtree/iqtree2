#include <stdio.h>
#include "lexer.h"

static const char * rawtext;
static int          rawtext_size;
static int          pos = 0;

int lex_table[SIZE_ASCII] = {
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN,     SYMBOL_TAB,      SYMBOL_CR, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN,      SYMBOL_LF, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*      */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/*  !"# */   SYMBOL_SPACE, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* $%&' */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* ()*+ */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* ,-./ */ SYMBOL_UNKNOWN,    SYMBOL_CHAR, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* 0123 */   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,
/* 4567 */   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,   SYMBOL_DIGIT,
/* 89:; */   SYMBOL_DIGIT,   SYMBOL_DIGIT, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,
/* <=>? */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,    SYMBOL_CHAR,
/* @ABC */ SYMBOL_UNKNOWN,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* DEFG */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* HIJK */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* LMNO */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* PQRS */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* TUVW */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* XYZ[ */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR, SYMBOL_UNKNOWN,
/* \]^_ */ SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN,    SYMBOL_CHAR,
/* `abc */ SYMBOL_UNKNOWN,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* defg */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* hijk */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* lmno */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* pqrs */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* tuvw */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR,
/* xyz{ */    SYMBOL_CHAR,    SYMBOL_CHAR,    SYMBOL_CHAR, SYMBOL_UNKNOWN,
/* |}~  */    SYMBOL_CHAR, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN, SYMBOL_UNKNOWN
 };

extern int debug ;

int 
get_next_byte (void)
{
  if (pos == rawtext_size) 
   {
     ++pos;
     return (EOS);
   }

  return (rawtext[pos++]);
}

int
get_next_symbol ()
{
  int ch, sym;

  ch = get_next_byte ();

  if (ch == EOS) return (SYMBOL_EOF);
  if (ch >= SIZE_ASCII) return SYMBOL_UNKNOWN;

  sym = lex_table[ch];

  if (sym == SYMBOL_LF)
   {
     if (get_next_byte() == '\n')
      {
        sym = SYMBOL_LFCR;
      }
     else
      {
       --pos;
      }
   }

  return sym;
  
}

struct lex_token
get_token (int * input)
{
  struct lex_token  token;
  int               start_pos;

  token.lexeme = rawtext + pos - 1;
  start_pos    = pos;

  switch (*input)
   {
     case SYMBOL_SPACE:
     case SYMBOL_TAB:
       do
        {
          *input = get_next_symbol();
        } while (*input == SYMBOL_SPACE || *input == SYMBOL_TAB);
       token.len   = pos - start_pos;
       token.class = LEX_WHITESPACE; 
       if (*input == SYMBOL_LFCR) --token.len;
       break;
       
     case SYMBOL_DIGIT:
       do
        {
          *input = get_next_symbol();   
        } while (*input == SYMBOL_DIGIT);
       
       if (*input != SYMBOL_CHAR)
        {
          token.len   = pos - start_pos;
          token.class = LEX_NUMBER;
        }
       else
        {
          do {
            *input = get_next_symbol();
          } while (*input == SYMBOL_CHAR || *input == SYMBOL_DIGIT);
          token.len   = pos - start_pos;
          token.class = LEX_STRING;
        }
       if (*input == SYMBOL_LFCR) --token.len;
       break;

     case SYMBOL_CHAR:
 //    if (debug == 1) printf ("HERE!!!!! %d \n", rawtext[pos - 1]);
       do
        {
          *input = get_next_symbol();
 //         if (debug == 1) printf ("TEST: %d\n", *input);
        } while (*input == SYMBOL_CHAR || *input == SYMBOL_DIGIT);
  //      if (debug == 1) printf ("pos startpod %d %d\n", pos, start_pos);
       token.len   = pos - start_pos;
       token.class = LEX_STRING;
       if (*input == SYMBOL_LFCR) --token.len;
       break;
       
     case SYMBOL_EOF:
       token.class = LEX_EOF;
       break;

     case SYMBOL_CR:
     case SYMBOL_LF:
     case SYMBOL_LFCR:
       do
        {
          *input = get_next_symbol();
        } while (*input == SYMBOL_CR || *input == SYMBOL_LFCR || *input == SYMBOL_LF);
       token.class = LEX_NEWLINE;
       break;
     case SYMBOL_UNKNOWN:
       token.class = LEX_UNKNOWN;
       break;
   }

  return (token);
}

void
init_lexan (const char * text, int n)
{
  rawtext      = text;
  rawtext_size = n;
  pos = 0;
}
