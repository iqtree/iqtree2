#ifndef __pll_LEXER__
#define __pll_LEXER__

#define  PLL_ASCII_SIZE                128
#define  PLL_EOS                       0x00000200

#define  PLL_SYM_CR                    1 << 0
#define  PLL_SYM_LF                    1 << 1
#define  PLL_SYM_LFCR                  1 << 2
#define  PLL_SYM_DIGIT                 1 << 3
#define  PLL_SYM_CHAR                  1 << 4
#define  PLL_SYM_SPACE                 1 << 5
#define  PLL_SYM_TAB                   1 << 6
#define  PLL_SYM_EOF                   1 << 7
#define  PLL_SYM_UNKNOWN               1 << 8
#define  PLL_SYM_DOT                   1 << 9
#define  PLL_SYM_COLON                 1 << 10
#define  PLL_SYM_OPAREN                1 << 11
#define  PLL_SYM_CPAREN                1 << 12
#define  PLL_SYM_COMMA                 1 << 13
#define  PLL_SYM_SEMICOLON             1 << 14
#define  PLL_SYM_EQUAL                 1 << 15
#define  PLL_SYM_DASH                  1 << 16
#define  PLL_SYM_SLASH                 1 << 17
#define  PLL_SYM_PLUS                  1 << 18

#define  PLL_TOKEN_NUMBER              1 << 0
#define  PLL_TOKEN_STRING              1 << 1
#define  PLL_TOKEN_EOF                 1 << 2
#define  PLL_TOKEN_WHITESPACE          1 << 3
#define  PLL_TOKEN_NEWLINE             1 << 4
#define  PLL_TOKEN_UNKNOWN             1 << 5
#define  PLL_TOKEN_COLON               1 << 6
#define  PLL_TOKEN_OPAREN              1 << 7
#define  PLL_TOKEN_CPAREN              1 << 8
#define  PLL_TOKEN_FLOAT               1 << 9
#define  PLL_TOKEN_COMMA               1 << 10
#define  PLL_TOKEN_SEMICOLON           1 << 11
#define  PLL_TOKEN_EQUAL               1 << 12
#define  PLL_TOKEN_DASH                1 << 13
#define  PLL_TOKEN_SLASH               1 << 14

#define CONSUME(x)         while (token.tokenType & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

typedef struct
 {
   int 	        tokenType;
   const char * lexeme;
   int          len;
 } pllLexToken;

int get_next_byte (void);
int get_next_symbol (void);
pllLexToken get_token (int * input);
void init_lexan (const char * text, long n);
void lex_table_amend_phylip (void);
void lex_table_amend_fasta (void);
void lex_table_restore (void);
#endif
