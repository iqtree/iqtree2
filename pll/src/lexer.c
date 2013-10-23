#include <stdio.h>
#include "lexer.h"

static const char * rawtext;
static long rawtext_size;
static long pos = 0;

int lex_table[PLL_ASCII_SIZE] = {
/*      */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN,     PLL_SYM_TAB,      PLL_SYM_CR,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN,      PLL_SYM_LF, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*      */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/*  !"# */   PLL_SYM_SPACE, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/* $%&' */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/* ()*+ */  PLL_SYM_OPAREN,  PLL_SYM_CPAREN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN,
/* ,-./ */   PLL_SYM_COMMA,    PLL_SYM_DASH,     PLL_SYM_DOT,     PLL_SYM_SLASH,
/* 0123 */   PLL_SYM_DIGIT,   PLL_SYM_DIGIT,   PLL_SYM_DIGIT,     PLL_SYM_DIGIT,
/* 4567 */   PLL_SYM_DIGIT,   PLL_SYM_DIGIT,   PLL_SYM_DIGIT,     PLL_SYM_DIGIT,
/* 89:; */   PLL_SYM_DIGIT,   PLL_SYM_DIGIT,   PLL_SYM_COLON, PLL_SYM_SEMICOLON,
/* <=>? */ PLL_SYM_UNKNOWN,   PLL_SYM_EQUAL, PLL_SYM_UNKNOWN,      PLL_SYM_CHAR,
/* @ABC */ PLL_SYM_UNKNOWN,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* DEFG */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* HIJK */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* LMNO */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* PQRS */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* TUVW */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* XYZ[ */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,   PLL_SYM_UNKNOWN,
/* \]^_ */ PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,      PLL_SYM_CHAR,
/* `abc */ PLL_SYM_UNKNOWN,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* defg */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* hijk */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* lmno */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* pqrs */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* tuvw */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,      PLL_SYM_CHAR,
/* xyz{ */    PLL_SYM_CHAR,    PLL_SYM_CHAR,    PLL_SYM_CHAR,   PLL_SYM_UNKNOWN,
/* |}~  */    PLL_SYM_CHAR, PLL_SYM_UNKNOWN, PLL_SYM_UNKNOWN,   PLL_SYM_UNKNOWN
 };

int 
get_next_byte (void)
{
  if (pos == rawtext_size) 
   {
     ++pos;
     return (PLL_EOS);
   }

  return (rawtext[pos++]);
}

int
get_next_symbol (void)
{
  int ch, sym;

  ch = get_next_byte ();

  if (ch == PLL_EOS) return (PLL_SYM_EOF);
  if (ch >= PLL_ASCII_SIZE) return (PLL_SYM_UNKNOWN);

  sym = lex_table[ch];

  if (sym == PLL_SYM_LF)
   {
     if (get_next_byte() == '\n')
      {
        sym = PLL_SYM_LFCR;
      }
     else
      {
        --pos;
      }
   }

  return sym;
}

pllLexToken
get_token (int * input)
{
  pllLexToken       token;
  int               start_pos;
  int               floating = 0;

  token.lexeme = rawtext + pos - 1;
  start_pos    = pos;

  switch (*input)
   {
     case PLL_SYM_SLASH:
       token.tokenType = PLL_TOKEN_SLASH;
       *input = get_next_symbol();
       break;

     case PLL_SYM_DASH:
       token.tokenType = PLL_TOKEN_DASH;
       *input = get_next_symbol();
       break;

     case PLL_SYM_EQUAL:
       token.tokenType = PLL_TOKEN_EQUAL;
       *input = get_next_symbol();
       break;

     case PLL_SYM_SEMICOLON:
       token.tokenType = PLL_TOKEN_SEMICOLON;
       *input = get_next_symbol();
       break;

     case PLL_SYM_COMMA:
       token.tokenType = PLL_TOKEN_COMMA;
       *input = get_next_symbol();
       break;

     case PLL_SYM_COLON:
       token.tokenType = PLL_TOKEN_COLON;
       *input = get_next_symbol();
       break;

     case PLL_SYM_OPAREN:
       token.tokenType = PLL_TOKEN_OPAREN;
       *input = get_next_symbol();
       break;

     case PLL_SYM_CPAREN:
       token.tokenType = PLL_TOKEN_CPAREN;
       *input = get_next_symbol();
       break;

     case PLL_SYM_SPACE:
     case PLL_SYM_TAB:
       do
        {
          *input = get_next_symbol();
        } while (*input == PLL_SYM_SPACE || *input == PLL_SYM_TAB);
       token.len   = pos - start_pos;
       token.tokenType = PLL_TOKEN_WHITESPACE; 
       if (*input == PLL_SYM_LFCR) --token.len;
       break;
       
     case PLL_SYM_DIGIT:
       do
        {
          *input = get_next_symbol();   
        } while (*input == PLL_SYM_DIGIT);

       if (*input == PLL_SYM_DOT)
        {
          floating = 1;
          do
           {
             *input = get_next_symbol ();
           } while (*input == PLL_SYM_DIGIT);
        }

       if (*input != PLL_SYM_CHAR)
        {
          token.len   = pos - start_pos;
          if (!floating)
            token.tokenType = PLL_TOKEN_NUMBER;
          else
            token.tokenType = PLL_TOKEN_FLOAT;
        }
       else
        {
          do {
            *input = get_next_symbol();
          } while (*input == PLL_SYM_CHAR || *input == PLL_SYM_DIGIT || *input == PLL_SYM_DOT);
          token.len   = pos - start_pos;
          token.tokenType = PLL_TOKEN_STRING;
        }

       if (*input == PLL_SYM_LFCR) --token.len;
       break;

     case PLL_SYM_CHAR:
       do
        {
          *input = get_next_symbol();
        } 
       while (*input == PLL_SYM_CHAR  || 
              *input == PLL_SYM_DIGIT || 
              *input == PLL_SYM_DASH  ||
              *input == PLL_SYM_DOT);
       token.len   = pos - start_pos;
       token.tokenType = PLL_TOKEN_STRING;
       if (*input == PLL_SYM_LFCR) --token.len;
       break;
       
     case PLL_SYM_EOF:
       token.tokenType = PLL_TOKEN_EOF;
       break;

     case PLL_SYM_CR:
     case PLL_SYM_LF:
     case PLL_SYM_LFCR:
       do
        {
          *input = get_next_symbol();
        } while (*input == PLL_SYM_CR || *input == PLL_SYM_LFCR || *input == PLL_SYM_LF);
       token.tokenType = PLL_TOKEN_NEWLINE;
       break;
     case PLL_SYM_UNKNOWN:
       token.tokenType = PLL_TOKEN_UNKNOWN;
       break;
     default:
       break;
   }

  return (token);
}

void
lex_table_amend_phylip (void)
{
  lex_table['-'] = lex_table['.'] = PLL_SYM_CHAR; 
}

void
lex_table_amend_fasta (void)
{
  lex_table['-'] = lex_table['.'] = lex_table['>'] = PLL_SYM_CHAR; 
}

void
lex_table_restore (void)
{
  lex_table['-'] = PLL_SYM_DASH;
  lex_table['.'] = PLL_SYM_DOT; 
  lex_table['>'] = PLL_SYM_UNKNOWN;
}

void
init_lexan (const char * text, long n)
{
  rawtext      = text;
  rawtext_size = n;
  pos          = 0;
}
