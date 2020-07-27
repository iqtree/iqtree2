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

#ifndef NCL_NXSTOKEN_H
#define NCL_NXSTOKEN_H

#include <iostream>	      //for std::ostream
#include "nxsstring.h"    //for NxsString and NxsStringVector
#include "nxsexception.h" //for NxsException
#include "nxsdefs.h"      //for file_pos
/*----------------------------------------------------------------------------------------------------------------------
|	NxsToken objects are used by NxsReader to extract words (tokens) from a NEXUS data file. NxsToken objects know to
|	correctly skip NEXUS comments and understand NEXUS punctuation, making reading a NEXUS file as simple as repeatedly
|	calling the GetNextToken() function and then interpreting the token returned. If the token object is not attached 
|	to an input stream, calls to GetNextToken() will have no effect. If the token object is not attached to an output
|	stream, output comments will be discarded (i.e., not output anywhere) and calls to Write or Writeln will be 
|	ineffective. If input and output streams have been attached to the token object, however, tokens are read one at a
|	time from the input stream, and comments are correctly read and either written to the output stream (if an output
|	comment) or ignored (if not an output comment). Sequences of characters surrounded by single quotes are read in as
|	single tokens. A pair of adjacent single quotes are stored as a single quote, and underscore characters are stored
|	as blanks.
*/
class NxsToken
	{
	public:

		enum NxsTokenFlags	/* For use with the variable labileFlags */
			{
			saveCommandComments		= 0x0001,	/* if set, command comments of the form [&X] are not ignored but are instead saved as regular tokens (without the square brackets, however) */
			parentheticalToken		= 0x0002,	/* if set, and if next character encountered is a left parenthesis, token will include everything up to the matching right parenthesis */
			curlyBracketedToken		= 0x0004,	/* if set, and if next character encountered is a left curly bracket, token will include everything up to the matching right curly bracket */
			doubleQuotedToken		= 0x0008,	/* if set, grabs entire phrase surrounded by double quotes */
			singleCharacterToken	= 0x0010,	/* if set, next non-whitespace character returned as token */
			newlineIsToken			= 0x0020,	/* if set, newline character treated as a token and atEOL set if newline encountered */
			tildeIsPunctuation		= 0x0040,	/* if set, tilde character treated as punctuation and returned as a separate token */
			useSpecialPunctuation	= 0x0080,	/* if set, character specified by the data member special is treated as punctuation and returned as a separate token */
			hyphenNotPunctuation	= 0x0100,	/* if set, the hyphen character is not treated as punctutation (it is normally returned as a separate token) */
			preserveUnderscores		= 0x0200,	/* if set, underscore characters inside tokens are not converted to blank spaces (normally, all underscores are automatically converted to blanks) */
			ignorePunctuation		= 0x0400	/* if set, the normal punctuation symbols are treated the same as any other darkspace characters */
			};

		NxsString		errormsg;

						NxsToken(std::istream &i);
		virtual			~NxsToken();

		bool			AtEOF();
		bool			AtEOL();
		bool			Abbreviation(NxsString s);
		bool			Begins(NxsString s, bool respect_case = false);
		void			BlanksToUnderscores();
		bool			Equals(NxsString s, bool respect_case = false);
		long			GetFileColumn() const;
		file_pos		GetFilePosition() const;
		long			GetFileLine() const;
		void			GetNextToken();
		NxsString		GetToken(bool respect_case = true);
		const char		*GetTokenAsCStr(bool respect_case = true);
		const NxsString	&GetTokenReference();
		int				GetTokenLength() const;
		bool			IsPlusMinusToken();
		bool			IsPunctuationToken();
		bool			IsWhitespaceToken();
		void			ReplaceToken(const NxsString s);
		void			ResetToken();
		void			SetSpecialPunctuationCharacter(char c);
		void			SetLabileFlagBit(int bit);
		bool			StoppedOn(char ch);
		void			StripWhitespace();
		void			ToUpper();
		void			Write(std::ostream &out);
		void			Writeln(std::ostream &out);

		virtual void	OutputComment(const NxsString &msg);
		void GetNextContiguousToken(char stop_char); // Added by BQM
	protected:

		void			AppendToComment(char ch);
		void			AppendToToken(char ch);
		char			GetNextChar();
		void			GetComment();
		void			GetCurlyBracketedToken();
		void			GetDoubleQuotedToken();
		void			GetQuoted();
		void			GetParentheticalToken();
		bool			IsPunctuation(char ch);
		bool			IsWhitespace(char ch);

	private:

		std::istream	&in;				/* reference to input stream from which tokens will be read */
		file_pos		filepos;			/* current file position (for Metrowerks compiler, type is streampos rather than long) */
		long			fileline;			/* current file line */
		long			filecol;			/* current column in current line (refers to column immediately following token just read) */
		NxsString		token;				/* the character buffer used to store the current token */
		NxsString		comment;			/* temporary buffer used to store output comments while they are being built */
		char			saved;				/* either '\0' or is last character read from input stream */
		bool			atEOF;				/* true if end of file has been encountered */
		bool			atEOL;				/* true if newline encountered while newlineIsToken labile flag set */
		char			special;			/* ad hoc punctuation character; default value is '\0' */
		int				labileFlags;		/* storage for flags in the NxsTokenFlags enum */
		char			punctuation[21];	/* stores the 20 NEXUS punctuation characters */
		char			whitespace[4];		/* stores the 3 whitespace characters: blank space, tab and newline */
	};

typedef NxsToken NexusToken;

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the token for functions that only need read only access - faster than GetToken.
*/
inline const NxsString &NxsToken::GetTokenReference()
	{
	return token;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	This function is called whenever an output comment (i.e., a comment beginning with an exclamation point) is found 
|	in the data file. This version of OutputComment does nothing; override this virtual function to display the output 
|	comment in the most appropriate way for the platform you are supporting.
*/
inline void NxsToken::OutputComment(
  const NxsString &msg)	/* the contents of the printable comment discovered in the NEXUS data file */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(msg)
#	endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `ch' to end of comment NxsString.
*/
inline void NxsToken::AppendToComment(
  char ch)	/* character to be appended to comment */
	{
	comment += ch;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds `ch' to end of current token.
*/
inline void NxsToken::AppendToToken(
  char ch)	/* character to be appended to token */
	{
	// First three lines proved necessary to keep Borland's implementation of STL from crashing
	// under some circumstances (may no longer be necessary)
	//
	char s[2];
	s[0] = ch;
	s[1] = '\0';

	token += s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads next character from in and does all of the following before returning it to the calling function:
|~
|	o if character read is either a carriage return or line feed, the variable line is incremented by one and the
|	  variable col is reset to zero
|	o if character read is a carriage return, and a peek at the next character to be read reveals that it is a line
|	  feed, then the next (line feed) character is read
|	o if either a carriage return or line feed is read, the character returned to the calling function is '\n' if 
|	  character read is neither a carriage return nor a line feed, col is incremented by one and the character is
|	  returned as is to the calling function
|	o in all cases, the variable filepos is updated using a call to the tellg function of istream.
|~
*/
inline char NxsToken::GetNextChar()
	{
	int ch = in.get();
	int failed = in.bad();
	if (failed)
		{
		errormsg = "Unknown error reading data file (check to make sure file exists)";
		throw NxsException(errormsg, *this);
		}

	if (ch == 13 || ch == 10)
		{
		fileline++;
		filecol = 1L;

		if (ch == 13 && (int)in.peek() == 10) 
			ch = in.get();

		atEOL = 1;
		}
	else if (ch == EOF)
		atEOF = 1;
	else
		{
		filecol++;
		atEOL = 0;
		}

#	if defined(__DECCXX)
		filepos = 0L;
#	else
    // BQM this cause crash compiling with clang under Windows!
//		filepos = in.tellg();
    filepos += 1;
#	endif

	if (atEOF)
		return '\0';
	else if (atEOL)
		return '\n';
	else
		return (char)ch;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if character supplied is considered a punctuation character. The following twenty characters are 
|	considered punctuation characters:
|>
|	()[]{}/\,;:=*'"`+-<>
|>
|	Exceptions:
|~
|	o The tilde character ('~') is also considered punctuation if the tildeIsPunctuation labile flag is set
|	o The special punctuation character (specified using the SetSpecialPunctuationCharacter) is also considered 
|	  punctuation if the useSpecialPunctuation labile flag is set
|	o The hyphen (i.e., minus sign) character ('-') is not considered punctuation if the hyphenNotPunctuation 
|	  labile flag is set
|~
|	Use the SetLabileFlagBit method to set one or more NxsLabileFlags flags in `labileFlags'
*/
inline bool NxsToken::IsPunctuation(
  char ch)	/* the character in question */
	{
	// PAUP 4.0b10 
	//  o allows ]`<> inside taxon names
	//  o allows `<> inside taxset names
	//
	bool is_punctuation = false;
	if (strchr(punctuation, ch))
		is_punctuation = true;
	if (labileFlags & tildeIsPunctuation  && ch == '~')
		is_punctuation = true;
	if (labileFlags & useSpecialPunctuation  && ch == special)
		is_punctuation = true;
	if (labileFlags & hyphenNotPunctuation  && ch == '-')
		is_punctuation = false;

	return is_punctuation;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if character supplied is considered a whitespace character. Note: treats '\n' as darkspace if labile
|	flag newlineIsToken is in effect.
*/
inline bool NxsToken::IsWhitespace(
  char ch)	/* the character in question */
	{
	bool ws = false;

	// If ch is found in the whitespace array, it's whitespace
	//
	if (strchr(whitespace, ch))
		ws = true;

	// Unless of course ch is the newline character and we're currently
	// treating newlines as darkspace!
	//
	if (labileFlags & newlineIsToken && ch == '\n')
		ws = false;

	return ws;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if last call to GetNextToken encountered the end-of-file character (or for some reason the 
|	input stream is now out of commission).
*/
inline bool NxsToken::AtEOF()
	{
	return atEOF;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if last call to GetNextToken encountered the newline character while the newlineIsToken 
|	labile flag was in effect.
*/
inline bool NxsToken::AtEOL()
	{
	return atEOL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts all blanks in token to underscore characters. Normally, underscores found in the tokens read from a NEXUS
|	file are converted to blanks automatically as they are read; this function reverts the blanks back to underscores. 
*/
inline void NxsToken::BlanksToUnderscores()
	{
	token.BlanksToUnderscores();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns value stored in `filecol', which keeps track of the current column in the data file (i.e., number of 
|	characters since the last new line was encountered).
*/
inline long  NxsToken::GetFileColumn() const
	{
	return filecol;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value stored in filepos, which keeps track of the current position in the data file (i.e., number of 
|	characters since the beginning of the file).  Note: for Metrowerks compiler, you must use the offset() method of 
|	the streampos class to use the value returned.
*/
inline file_pos  NxsToken::GetFilePosition() const
	{
	return filepos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value stored in `fileline', which keeps track of the current line in the data file (i.e., number of new 
|	lines encountered thus far).
*/
inline long  NxsToken::GetFileLine() const
	{
	return fileline;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the data member `token'. Specifying false for`respect_case' parameter causes all characters in `token'
|	to be converted to upper case before `token' is returned. Specifying true results in GetToken returning exactly 
|	what it read from the file.
*/
inline NxsString NxsToken::GetToken(
  bool respect_case)	/* determines whether token is converted to upper case before being returned */
	{
	if (!respect_case)
		ToUpper();

	return token;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the data member `token' as a C-style string. Specifying false for`respect_case' parameter causes all 
|	characters in `token' to be converted to upper case before the `token' C-string is returned. Specifying true 
|	results in GetTokenAsCStr returning exactly what it read from the file.
*/
inline const char *NxsToken::GetTokenAsCStr(
  bool respect_case)	/* determines whether token is converted to upper case before being returned */
	{
	if (!respect_case)
		ToUpper();

	return token.c_str();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns token.size().
*/
inline int NxsToken::GetTokenLength() const
	{
	return token.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if current token is a single character and this character is either '+' or '-'.
*/
inline bool NxsToken::IsPlusMinusToken()
	{
	if (token.size() == 1 && ( token[0] == '+' || token[0] == '-') )
		return true;
	else
		return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if current token is a single character and this character is a punctuation character (as defined in 
|	IsPunctuation function).
*/
inline bool NxsToken::IsPunctuationToken()
	{
	if (token.size() == 1 && IsPunctuation( token[0]))
		return true;
	else
		return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if current token is a single character and this character is a whitespace character (as defined in 
|	IsWhitespace function).
*/
inline bool NxsToken::IsWhitespaceToken()
	{
	if (token.size() == 1 && IsWhitespace( token[0]))
		return true;
	else
		return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces current token NxsString with s.
*/
inline void NxsToken::ReplaceToken(
  const NxsString s)	/* NxsString to replace current token NxsString */
	{
	token = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets token to the empty NxsString ("").
*/
inline void NxsToken::ResetToken()
	{
	token.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the special punctuation character to `c'. If the labile bit useSpecialPunctuation is set, this character will 
|	be added to the standard list of punctuation symbols, and will be returned as a separate token like the other 
|	punctuation characters.
*/
inline void NxsToken::SetSpecialPunctuationCharacter(
  char c)	/* the character to which `special' is set */
	{
	special = c;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the bit specified in the variable `labileFlags'. The available bits are specified in the NxsTokenFlags enum.
|	All bits in `labileFlags' are cleared after each token is read.
*/
inline void NxsToken::SetLabileFlagBit(
  int bit)	/* the bit (see NxsTokenFlags enum) to set in `labileFlags' */
	{
	labileFlags |= bit;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Checks character stored in the variable saved to see if it matches supplied character `ch'. Good for checking such 
|	things as whether token stopped reading characters because it encountered a newline (and labileFlags bit 
|	newlineIsToken was set):
|>
|	StoppedOn('\n');
|>
|	or whether token stopped reading characters because of a punctuation character such as a comma:
|>
|	StoppedOn(',');
|>
*/
inline bool NxsToken::StoppedOn(
  char ch)	/* the character to compare with saved character */
	{
	if (saved == ch)
		return true;
	else
		return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Simply outputs the current NxsString stored in `token' to the output stream `out'. Does not send a newline to the 
|	output stream afterwards.
*/
inline void NxsToken::Write(
  std::ostream &out)	/* the output stream to which to write token NxsString */
	{
	out << token;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Simply outputs the current NxsString stored in `token' to the output stream `out'. Sends a newline to the output 
|	stream afterwards.
*/
inline void NxsToken::Writeln(
  std::ostream &out)	/* the output stream to which to write `token' */
	{
	out << token << std::endl;
	}

/**
 * Added by BQM: return the contiguous string (including white space) as token
 * until hitting stop_char
 * @param stop_char a character to stop reading in
 */
inline void NxsToken::GetNextContiguousToken(char stop_char) {
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
	for (;;) {

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
		if (ch == stop_char) {
			saved = ch;
			break;
		}
		AppendToToken(ch);
	}
	// Skip ending whitespace
	if (token.empty()) return;
	NxsString::iterator last = token.end();
	while (last != token.begin() && IsWhitespace(*(last-1))) {
		last--;
	}
	if (last != token.end()) token.erase(last, token.end());
}

#endif
