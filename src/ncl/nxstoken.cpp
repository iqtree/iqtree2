//	Copyright (C) 1999-2003 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library) version 2.0.
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Sets atEOF and atEOL to false, comment and token to the empty string, filecol and fileline to 1, filepos to 0, 
|	labileFlags to 0 and saved and special to the null character. Initializes the istream reference data 
|	member in to the supplied istream `i'.
*/
NxsToken::NxsToken(
  istream &i)	/* the istream object to which the token is to be associated */
  : in(i)
	{
	atEOF		= false;
	atEOL		= false;
	comment.clear();
	filecol		= 1L;
	fileline	= 1L;
	filepos		= 0L;
	labileFlags	= 0;
	saved		= '\0';
	special		= '\0';
	
	whitespace[0]  = ' ';
	whitespace[1]  = '\t';
	whitespace[2]  = '\n';
	whitespace[3]  = '\0';

	punctuation[0]	= '(';
	punctuation[1]	= ')';
	punctuation[2]	= '[';
	punctuation[3]	= ']';
	punctuation[4]	= '{';
	punctuation[5]	= '}';
	punctuation[6]	= '/';
	punctuation[7]	= '\\';
	punctuation[8]	= ',';
	punctuation[9]	= ';';
	punctuation[10]	= ':';
	punctuation[11]	= '=';
	punctuation[12]	= '*';
	punctuation[13]	= '\'';
	punctuation[14]	= '"';
	punctuation[15]	= '`';
	punctuation[16]	= '+';
	punctuation[17]	= '-';
	punctuation[18]	= '<';
	punctuation[19]	= '>';
	punctuation[20]	= '\0';
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing needs to be done; all objects take care of deleting themselves.
*/
NxsToken::~NxsToken()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads rest of comment (starting '[' already input) and acts accordingly. If comment is an output comment, and if 
|	an output stream has been attached, writes the output comment to the output stream. Otherwise, output comments are 
|	simply ignored like regular comments. If the labileFlag bit saveCommandComments is in effect, the comment (without 
|	the square brackets) will be stored in token. 
*/
void NxsToken::GetComment()
	{
	// Set comment level to 1 initially.  Every ']' encountered reduces
	// level by one, so that we know we can stop when level becomes 0.
	//
	int level = 1;

	// Get first character
	//
	char ch = GetNextChar();
	if (atEOF)
		{
		errormsg = "Unexpected end of file inside comment";
		throw NxsException( errormsg, GetFilePosition(), GetFileLine(), GetFileColumn());
		}

	// See if first character is the output comment symbol ('!')
	// or command comment symbol (&)
	//
	int printing = 0;
	int command = 0;
	if (ch == '!')
		printing = 1;
	else if (ch == '&' && labileFlags & saveCommandComments)
		{
		command = 1;
		AppendToToken(ch);
		}
	else if (ch == ']')
		return;

	// Now read the rest of the comment
	//
	for(;;)
		{
		ch = GetNextChar();
		if (atEOF)
			break;

		if (ch == ']')
			level--;
		else if (ch == '[')
			level++;

		if (level == 0)
			break;

		if (printing)
			AppendToComment(ch);
		else if (command)
			AppendToToken(ch);
		}

	if (printing)
		{
		// Allow output comment to be printed or displayed in most appropriate
		// manner for target operating system
		//
		OutputComment(comment);

		// Now that we are done with it, free the memory used to store the comment
		//
		//comment;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads rest of a token surrounded with curly brackets (the starting '{' has already been input) up to and including
|	the matching '}' character. All nested curly-bracketed phrases will be included.
*/
void NxsToken::GetCurlyBracketedToken()
	{
	// Set level to 1 initially.  Every '}' encountered reduces
	// level by one, so that we know we can stop when level becomes 0.
	//
	int level = 1;

	char ch;
	for(;;)
		{
		ch = GetNextChar();
		if (atEOF)
			break;

		if (ch == '}')
			level--;
		else if (ch == '{')
			level++;

		AppendToToken(ch);

		if (level == 0)
			break;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Gets remainder of a double-quoted NEXUS word (the first double quote character was read in already by GetNextToken).
|	This function reads characters until the next double quote is encountered. Tandem double quotes within a 
|	double-quoted NEXUS word are not allowed and will be treated as the end of the first word and the beginning of the 
|	next double-quoted NEXUS word. Tandem single quotes inside a double-quoted NEXUS word are saved as two separate 
|	single quote characters; to embed a single quote inside a double-quoted NEXUS word, simply use the single quote by 
|	itself (not paired with another tandem single quote).
*/
void NxsToken::GetDoubleQuotedToken()
	{
	char ch;

	for(;;)
		{
		ch = GetNextChar();
		if (atEOF)
			break;

		if (ch == '\"')
			break;
		else
			AppendToToken(ch);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Gets remainder of a quoted NEXUS word (the first single quote character was read in already by GetNextToken). This
|	function reads characters until the next single quote is encountered. An exception occurs if two single quotes occur
|	one after the other, in which case the function continues to gather characters until an isolated single quote is
|	found. The tandem quotes are stored as a single quote character in the token NxsString.
*/
void NxsToken::GetQuoted()
	{
	char ch;

	for(;;)
		{
		ch = GetNextChar();
		if (atEOF)
			break;

		if (ch == '\'' && saved == '\'')
			{
			// Paired single quotes, save as one single quote
			//
			AppendToToken(ch);
			saved = '\0';
			}
		else if (ch == '\'' && saved == '\0')
			{
			// Save the single quote to see if it is followed by another
			//
			saved = '\'';
			}
		else if (saved == '\'')
			{
			// Previously read character was single quote but this is something else, save current character so that it will
			// be the first character in the next token read
			//
			saved = ch;
			break;
			}
		else
			AppendToToken(ch);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads rest of parenthetical token (starting '(' already input) up to and including the matching ')' character.  All
|	nested parenthetical phrases will be included.
*/
void NxsToken::GetParentheticalToken()
	{
	// Set level to 1 initially.  Every ')' encountered reduces
	// level by one, so that we know we can stop when level becomes 0.
	//
	int level = 1;

	char ch;
	for(;;)
		{
		ch = GetNextChar();
		if (atEOF)
			break;

		if (ch == ')')
			level--;
		else if (ch == '(')
			level++;

		AppendToToken(ch);

		if (level == 0)
			break;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if token begins with the capitalized portion of `s' and, if token is longer than `s', the remaining 
|	characters match those in the lower-case portion of `s'. The comparison is case insensitive. This function should be
|	used instead of the Begins function if you wish to allow for abbreviations of commands and also want to ensure that 
|	user does not type in a word that does not correspond to any command.
*/
bool NxsToken::Abbreviation(
  NxsString s)	/* the comparison string */
	{
	int k;
	int slen = s.size();
	int tlen = token.size();
	char tokenChar, otherChar;

	// The variable mlen refers to the "mandatory" portion
	// that is the upper-case portion of s
	//
	int mlen;
	for (mlen = 0; mlen < slen; mlen++)
		{
		if (!isupper(s[mlen]))
			break;
		}

	// User must have typed at least mlen characters in
	// for there to even be a chance at a match
	//
	if (tlen < mlen)
		return false;

	// If user typed in more characters than are contained in s,
	// then there must be a mismatch
	//
	if (tlen > slen)
		return false;

	// Check the mandatory portion for mismatches
	//
	for (k = 0; k < mlen; k++)
		{
		tokenChar = (char)toupper( token[k]);
		otherChar = s[k];
		if (tokenChar != otherChar)
			return false;
		}

	// Check the auxiliary portion for mismatches (if necessary)
	//
	for (k = mlen; k < tlen; k++)
		{
		tokenChar = (char)toupper( token[k]);
		otherChar = (char)toupper( s[k]);
		if (tokenChar != otherChar)
			return false;
		}

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if token NxsString begins with the NxsString `s'. This function should be used instead of the Equals 
|	function if you wish to allow for abbreviations of commands.
*/
bool NxsToken::Begins(
  NxsString s,			/* the comparison string */
  bool respect_case)	/* determines whether comparison is case sensitive */
	{
	unsigned k;
	char tokenChar, otherChar;

	unsigned slen = s.size();
	if (slen > token.size())
		return false;

	for (k = 0; k < slen; k++)
		{
		if (respect_case)
			{
			tokenChar = token[k];
			otherChar = s[k];
			}
		else
			{
			tokenChar = (char)toupper( token[k]);
			otherChar = (char)toupper( s[k]);
			}

		if (tokenChar != otherChar)
			return false;
		}

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if token NxsString exactly equals `s'. If abbreviations are to be allowed, either Begins or 
|	Abbreviation should be used instead of Equals.
*/
bool NxsToken::Equals(
  NxsString s,			/* the string for comparison to the string currently stored in this token */
  bool respect_case)	/* if true, comparison will be case-sensitive */
	{
	unsigned k;
	char tokenChar, otherChar;

	unsigned slen = s.size();
	if (slen != token.size())
		return false;

	for (k = 0; k < token.size(); k++)
		{
		if (respect_case)
			{
			tokenChar = token[k];
			otherChar = s[k];
			}
		else
			{
			tokenChar = (char)toupper( token[k]);
			otherChar = (char)toupper( s[k]);
			}
		if (tokenChar != otherChar)
			return false;
		}

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads characters from in until a complete token has been read and stored in token. GetNextToken performs a number 
|	of useful operations in the process of retrieving tokens:
|~
|	o any underscore characters encountered are stored as blank spaces (unless the labile flag bit preserveUnderscores
|	  is set)
|	o if the first character of the next token is an isolated single quote, then the entire quoted NxsString is saved 
|	  as the next token
|	o paired single quotes are automatically converted to single quotes before being stored
|	o comments are handled automatically (normal comments are treated as whitespace and output comments are passed to 
|	  the function OutputComment which does nothing in the NxsToken class but can be overridden in a derived class to 
|	  handle these in an appropriate fashion)
|	o leading whitespace (including comments) is automatically skipped
|	o if the end of the file is reached on reading this token, the atEOF flag is set and may be queried using the AtEOF 
|	  member function
|	o punctuation characters are always returned as individual tokens (see the Maddison, Swofford, and Maddison paper 
|	  for the definition of punctuation characters) unless the flag ignorePunctuation is set in labileFlags,
|	  in which case the normal punctuation symbols are treated just like any other darkspace character.
|~
|	The behavior of GetNextToken may be altered by using labile flags. For example, the labile flag saveCommandComments 
|	can be set using the member function SetLabileFlagBit. This will cause comments of the form [&X] to be saved as 
|	tokens (without the square brackets), but only for the aquisition of the next token. Labile flags are cleared after 
|	each application.
*/
void NxsToken::GetNextToken()
	{
	ResetToken();

	char ch = ' ';
	if (saved == '\0' || IsWhitespace(saved))
		{
		// Skip leading whitespace
		//
		while( IsWhitespace(ch) && !atEOF)
			ch = GetNextChar();
		saved = ch;
		}

	for(;;)
		{
		// Break now if singleCharacterToken mode on and token length > 0.
		//
		if (labileFlags & singleCharacterToken && token.size() > 0)
			break;

		// Get next character either from saved or from input stream.
		//
		if (saved != '\0')
			{
			ch = saved;
			saved = '\0';
			}
		else
			ch = GetNextChar();

		// Break now if we've hit EOF.
		//
		if (atEOF)
			break;

		if (ch == '\n' && labileFlags & newlineIsToken)
			{
			if (token.size() > 0)
				{
				// Newline came after token, save newline until next time when it will be 
				// reported as a separate token.
				//
				atEOL = 0;
				saved = ch;
				}
			else
				{
				atEOL = 1;
				AppendToToken(ch);
				}
			break;
			}

		else if (IsWhitespace(ch))
			{
			// Break only if we've begun adding to token (remember, if we hit a comment before a token,
			// there might be further white space between the comment and the next token).
			//
			if (token.size() > 0)
				break;
			}

		else if (ch == '_')
			{
			// If underscores are discovered in unquoted tokens, they should be 
			// automatically converted to spaces.
			//
			if (!(labileFlags & preserveUnderscores))
				ch = ' ';
			AppendToToken(ch);
			}

		else if (ch == '[')
			{
			// Get rest of comment and deal with it, but notice that we only break if the comment ends a token,
			// not if it starts one (comment counts as whitespace). In the case of command comments 
			// (if saveCommandComment) GetComment will add to the token NxsString, causing us to break because
			// token.size() will be greater than 0.
			//
			GetComment();
			if (token.size() > 0)
			break;
			}

		else if (ch == '(' && labileFlags & parentheticalToken)
			{
			AppendToToken(ch);

			// Get rest of parenthetical token.
			//
			GetParentheticalToken();
			break;
			}

		else if (ch == '{' && labileFlags & curlyBracketedToken)
			{
			AppendToToken(ch);

			// Get rest of curly-bracketed token.
			//
			GetCurlyBracketedToken();
			break;
			}

		else if (ch == '\"' && labileFlags & doubleQuotedToken)
			{
			// Get rest of double-quoted token.
			//
			GetDoubleQuotedToken();
			break;
			}

		else if (ch == '\'')
			{
			if (token.size() > 0)
				{
				// We've encountered a single quote after a token has
				// already begun to be read; should be another tandem
				// single quote character immediately following.
				//
				ch = GetNextChar();
				if (ch == '\'')
					AppendToToken(ch);
				else
					{
					errormsg = "Expecting second single quote character";
					throw NxsException( errormsg, GetFilePosition(), GetFileLine(), GetFileColumn());
					}
				}
			else
				{
				// Get rest of quoted NEXUS word and break, since
				// we will have eaten one token after calling GetQuoted.
				//
				GetQuoted();
				}
			break;
			}

		else if (IsPunctuation(ch))
			{
			if (token.size() > 0)
				{
				// If we've already begun reading the token, encountering
				// a punctuation character means we should stop, saving
				// the punctuation character for the next token.
				//
				saved = ch;
				break;
				}
			else
				{
				// If we haven't already begun reading the token, encountering
				// a punctuation character means we should stop and return
				// the punctuation character as this token (i.e., the token
				// is just the single punctuation character.
				//
				AppendToToken(ch);
				break;
				}
			}

		else
			{
			AppendToToken(ch);
			}

		}

	labileFlags = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Strips whitespace from currently-stored token. Removes leading, trailing, and embedded whitespace characters.
*/
void NxsToken::StripWhitespace()
	{
	NxsString s;
	for (unsigned j = 0; j < token.size(); j++)
		{
		if (IsWhitespace( token[j]))
			continue;
		s += token[j];
		}
	token = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts all alphabetical characters in token to upper case.
*/
void NxsToken::ToUpper()
	{
	for (unsigned i = 0; i < token.size(); i++)
		token[i] = (char)toupper(token[i]);
	}

