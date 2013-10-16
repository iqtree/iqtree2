#ifndef LEXER_H
#define LEXER_H

#define  SIZE_ASCII              128
#define  EOS                     0x00000200

#define  SYMBOL_CR               0x00000001
#define  SYMBOL_LF               0x00000002
#define  SYMBOL_LFCR             0x00000004
#define  SYMBOL_DIGIT            0x00000008
#define  SYMBOL_CHAR             0x00000010
#define  SYMBOL_SPACE            0x00000020
#define  SYMBOL_TAB              0x00000040
#define  SYMBOL_EOF              0x00000080
#define  SYMBOL_UNKNOWN          0x00000100

#define  LEX_NUMBER              0x00000001
#define  LEX_STRING              0x00000002
#define  LEX_EOF                 0x00000004
#define  LEX_WHITESPACE          0x00000008
#define  LEX_NEWLINE             0x00000010
#define  LEX_UNKNOWN             0x00000020

struct lex_token
 {
   int          class;
   const char * lexeme;
   int          len;
 };

int get_next_byte (void);
int get_next_symbol (void);
struct lex_token get_token (int * input);
void init_lexan (const char * text, int n);

#endif
