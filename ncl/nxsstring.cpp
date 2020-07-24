//	Copyright (C) 1999-2003 Paul O. Lewis and Mark T. Holder
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

/*--------------------------------------------------------------------------------------------------------------------------
|	Capitalizes every character in the stored string.
*/
NxsString &NxsString::ToUpper()
	{
	for (NxsString::iterator sIt = begin(); sIt != end(); sIt++)
		*sIt = (char) toupper(*sIt);
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Appends a string representation of the supplied double to the stored string and returns a reference to itself.
*/
NxsString &NxsString::operator+=(
  const double d)	/* the double value to append */
	{
	char tmp[81];

	// Create a C-string representing the supplied double value. 
	// The # causes a decimal point to always be output.
	//
	sprintf(tmp, "%#3.6f", d);
	unsigned tmplen = (unsigned)strlen(tmp);

	// If the C-string has a lot of trailing zeros, lop them off
	//
	for (;;)
		{
		if (tmplen < 3 || tmp[tmplen-1] != '0' || tmp[tmplen-2] == '.')
			break;
		tmp[tmplen-1] = '\0';
		tmplen--;
		}

	append(tmp);
	return *this;
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	Adds `n' copies of the character `c' to the end of the stored string and returns a reference to itself.
*/
NxsString &NxsString::AddTail(
  char c,		/* the character to use in the appended tail */
  unsigned n)	/* the number of times `c' is to be appended */
	{
	char s[2];
	s[0] = c;
	s[1] = '\0';

	for (unsigned i = 0; i < n; i++)
		append(s);

	return *this;
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	Replaces the stored string with a copy of itself surrounded by single quotes (single quotes inside the string are 
|	converted to the '' pair of characters that signify a single quote). Returns a reference to itself.
*/
NxsString &NxsString::AddQuotes() 
	{
	NxsString withQuotes;
	int len = length();
	withQuotes.reserve(len + 4);
	withQuotes += '\'';
	for (NxsString::const_iterator sIt = begin(); sIt != end(); sIt++)
		{
		withQuotes += *sIt;
		if (*sIt == '\'')
			withQuotes += '\'';
		}
	withQuotes += '\'';
	*this = withQuotes;

	return *this;
	}
 
/*--------------------------------------------------------------------------------------------------------------------------
|	Appends a printf-style formatted string onto the end of this NxsString and returns the number of characters added to the 
|	string. For example, the following code would result in the string s being set to "ts-tv rate ratio = 4.56789":
|>
|	double kappa = 4.56789;
|	NxsString s;
|	s.PrintF("ts-tv rate ratio = %.5f", kappa);
|>
*/
int NxsString::PrintF(
  const char *formatStr,	/* the printf-style format string */
  ...)						/* other arguments referred to by the format string */
  	{
  	const int kInitialBufferSize = 256;
  	char buf[kInitialBufferSize];

	// Create a pointer to the list of optional arguments
	//
  	va_list argList;

	// Set arg_ptr to the first optional argument in argList. The
	// second argument (formatStr) is the last non-optional argument.
	// 
  	va_start(argList, formatStr);

	// If vsnprintf returns -1, means kInitialBufferSize was not large enough.
	// In this case, only kInitialBufferSize bytes are written.
	//
  	int nAdded = vsnprintf(buf, kInitialBufferSize, formatStr, argList);

	// Reset the argument list pointer
	//
  	va_end(argList);

	// Currently, if formatted string is too long to fit into the supplied buf,
	// just adding a terminating '\0' and returning the truncated string
	// Need to think of a better solution
	//
	if (nAdded  < 0 || nAdded >= kInitialBufferSize)
		buf[kInitialBufferSize - 1] = '\0';

	*this << buf;

#if 0
	// This part not being used anymore because there seems to be some differences
	// between compilers in what is returned from the vsnprintf function. VC returns
	// -1 if string is too long, Metrowerks returns the number of bytes that it would
	// have used had there been enough space! 
	//

  	if (nAdded >= kInitialBufferSize)
  		{
  		char *tempbuf = new char[nAdded + 2];

  		va_list argList;
  		va_start(argList, formatStr);

  		unsigned newNAdded = vsnprintf(tempbuf, nAdded + 1, formatStr, argList);

  		va_end(argList);

  		assert(nAdded == newNAdded);

  		*this << tempbuf;
  		delete [] tempbuf;
  		}
  	else
  		*this << buf;

#endif

  	return nAdded;
  	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the string is a abbreviation (or complete copy) of the argument `s'.
*/
bool NxsString::IsStdAbbreviation(
  const NxsString &s,	/* the string for which the stored string is potentially an abbreviation */
  bool respectCase)		/* if true, comparison will be case-sensitive */
  const
	{
	if (empty())
		return false;

	// s is the unabbreviated comparison string
	//
	const unsigned slen = static_cast<unsigned long>(s.size());

	// t is the stored string
	//
	const unsigned tlen = static_cast<unsigned long>(size());

	// t cannot be an abbreviation of s if it is longer than s
	//
	if (tlen > slen)
		return false;

	// Examine each character in t and return false (meaning "not an abbreviation")
	// if at any point the corresponding character in s is different
	//
	for (unsigned k = 0; k < tlen; k++)
		{
		if (respectCase)
			{
			if ((*this)[k] != s[k])
				return false;
			}
		else if (toupper((*this)[k]) != toupper(s[k]))
			return false;
		}

	return true;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the stored string is a case-insensitive abbreviation (or complete copy) of `s' and the stored string 
| 	has all of the characters that are in the initial capitalized portion of `s'. For example if `s' is "KAPpa" then 
|	"kappa", "kapp", or "kap" (with any capitalization pattern) will return true and all other strings will return false. 
|	Always returns false if the stored string has length of zero.
*/
bool NxsString::IsCapAbbreviation(
  const NxsString &s)	/* the string for which the stored string is potentially an abbreviation */
  const
	{
	if (empty())
		return false;

	// s is the unabbreviated comparison string
	//
	const unsigned slen = static_cast<unsigned>(s.size());

	// t is the stored string
	//
	const unsigned tlen = static_cast<unsigned>(size());

	// If the stored string is longer than s then it cannot be an abbreviation of s
	//
	if (tlen > slen)
		return false;
	
	unsigned k = 0;
	for (; k < slen; k++) 
		{
		if (isupper(s[k]))	
			{
			// If still in the uppercase portion of s and we've run out of characters
			// in t, then t is not a valid abbrevation of s
			//
			if (k >= tlen)
				return false;

			// If kth character in t is not equal to kth character in s, then
			// t is not an abbrevation of s
			//
			char tokenChar = (char)toupper((*this)[k]);
			if (tokenChar != s[k])
				return false;
			}
		else if (!isalpha(s[k]))
			{
			// Get here if we are no longer in the upper case portion of s and 
			// s[k] is not an alphabetic character. This section is necessary because
			// we are dealing with a section of s that is not alphabetical and thus
			// we cannot tell whether this should be part of the abbrevation or not
			// (i.e. we cannot tell if it is capitalized or not). In this case, we
			// pretend that we are still in the upper case portion of s and return
			// false if we have run out of characters in t (meaning that the abbreviation
			// was too short) or we find a mismatch.
			//
			if (k >= tlen)
				return false;

			if ((*this)[k] != s[k])
				return false;
			}
		else
			{
			// Get here if we are no longer in the upper case portion of s and
			// s[k] is an alphabetic character. Just break because we have determined
			// that t is in fact a valid abbreviation of s.
			//
			break;
			}
		}

	// Check the lower case portion of s and any corresponding characters in t for mismatches
	// Even though the abbreviation is valid up to this point, it will become invalid if
	// any mismatches are found beyond the upper case portion of s
	//
	for (; k < tlen; k++)
		{
  		const char tokenChar = (char)toupper((*this)[k]);
  		const char otherChar = (char)toupper(s[k]);
		if (tokenChar != otherChar)
			return false;
		}

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Right-justifies `x' in a field `w' characters wide, using blank spaces to fill in unused portions on the left-hand 
|	side of the field. Specify true for `clear_first' to first empty the string. Assumes `w' is large enough to 
|	accommodate the string representation of `x'.
*/
NxsString &NxsString::RightJustifyLong(
  long x,			/* long value to right justify */
  unsigned int w,	/* width of field */
  bool clear_first)	/* if true, initialize string first to empty string */
	{
	bool x_negative = (x < 0L ? true : false);
	unsigned long xabs = (x_negative ? (-x) : x);
	unsigned num_spaces = w;

	// If w = 10 and x = 123, we need 7 blank spaces before x
	// log10(123) is 2.09, indicating that x is at least 10^2 = 100 but not
	// 10^3 = 1000, thus x requires at least 3 characters to display
	// 
	unsigned x_width = (x == 0 ? 1 :1 + (int)log10((double)xabs));
	if (x_negative)
		x_width++;	// for the minus sign

	assert(x_width <= num_spaces);
	num_spaces -= x_width;

	if (clear_first)
		erase();

	for (unsigned k = 0; k < num_spaces; k++)
		*this += ' ';

	if (x_negative)
		*this += '-';

	*this += xabs;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Right-justifies `x' in a field `w' characters wide with precision `p', using blank spaces to fill in unused 
|	portions on the left-hand side of the field. Specify true for `clear_first' to first empty the string. Assumes that
|	the specified width is enough to accommodate the string representation of `x'.
*/
NxsString &NxsString::RightJustifyDbl(
  double x,				/* double value to right justify */
  unsigned w,			/* width of field */
  unsigned p,			/* precision to use when displaying `x' */
  bool clear_first)		/* if true, initialize stored string first to the empty string */
	{
	if (clear_first)
		erase();

	char fmtstr[81];
	sprintf(fmtstr, "%%.%df", p);
	NxsString tmp;
	tmp.PrintF(fmtstr, x);

	unsigned num_spaces = w - tmp.length();
	assert(num_spaces >= 0);

	for (unsigned k = 0; k < num_spaces; k++)
		*this += ' ';

	*this += tmp;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Right-justifies `s' in a field `w' characters wide, using blank spaces to fill in unused portions on the left-hand
|	side of the field. Specify true for `clear_first' to first empty the string. Assumes that the specified width is 
|	enough to accommodate `s'.
*/
NxsString &NxsString::RightJustifyString(
  const NxsString &s,	/* string to right justify */
  unsigned w,			/* width of field */
  bool clear_first)		/* if true, initialize string first to the empty string */
	{
	if (clear_first)
		erase();

	unsigned num_spaces = w - s.length();
	assert(num_spaces >= 0);

	for (unsigned k = 0; k < num_spaces; k++)
		*this += ' ';

	*this += s;
	return *this;
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	Returns true if the string needs to be surrounded by single-quotes to make it a single nexus token.
*/
bool NxsString::QuotesNeeded() const
	{
	bool quotes_needed = false;

	for (NxsString::const_iterator sIt = begin(); sIt != end(); sIt++)
		{
		char c = (*sIt);

		if (!isgraph(c))
			{
			// The standard C function isgraph returns zero if c is either a space or is not a printable character.
			//
			quotes_needed = true;
			}
		else if (strchr("(){}\"-]/\\,;:=*`+<>", c) != NULL)
			{
			// Get here if c is any NEXUS punctuation mark except left square bracket ([) or apostrophe (').
			// Left square bracket characters and apostrophes never get returned as punctuation by NxsToken,
			// so we should never encounter them here. 
			//
			
			if (length() > 1)
				quotes_needed = true;
			}
		else if (c == '\'' || c == '[')
			{
			// Get here if c is either an apostrophe or left square bracket. Quotes are needed if one of these
			// characters is all there is to this string
			//
			//@POL Mark, I'm confused.
			//
			quotes_needed = true;
			}

		if (quotes_needed)
			break;
		}

	return quotes_needed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts any blank spaces found in the stored string to the underscore character.
*/
NxsString &NxsString::BlanksToUnderscores()
	{
	unsigned len = length();
	for (unsigned k = 0; k < len; k++)
		{
		char &ch = at(k);
		if (ch == ' ')
			ch = '_';
		}
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts any underscore characters found in the stored string to blank spaces.
*/
NxsString &NxsString::UnderscoresToBlanks()
	{
	unsigned len = length();
	for (unsigned k = 0; k < len; k++)
		{
		char &ch = at(k);
		if (ch == '_')
			ch = ' ';
		}
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Shortens stored string to `n' - 3 characters, making the last three characters "...". If string is already less than 
|	`n' characters in length, this function has no effect. This is useful when it is desirable to show some of the
|	contents of a string, even when the string will not fit in its entirety into the space available for displaying it.
|	Assumes that `n' is at least 4.
*/
NxsString &NxsString::ShortenTo(
  unsigned n)	/* maximum number of characters available for displaying the string */
	{
	assert(n > 3);
	if (length() <= static_cast<unsigned>(n))
		return *this;
	
	NxsString s;
	for (NxsString::iterator sIt = begin(); sIt != end(); sIt++)
		{
		s += (*sIt);
		if (s.length() >= n - 3)
			break;
		}
	s += "...";

	*this = s;
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts every character in the stored string to its lower case equivalent.
*/
NxsString &NxsString::ToLower()
	{
	for (NxsString::iterator sIt = begin(); sIt != end(); sIt++)
		{
		char c = (char)tolower(*sIt);
		*sIt = c;
		}
	return *this;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the stored string can be interpreted as a double value, and returns false otherwise.
*/
bool NxsString::IsADouble() const
	{
	const char	*str			= c_str();
	unsigned	i				= 0;
	bool		hadDecimalPt	= false;
	bool		hadExp			= false;
	bool		hadDigit		= false;
	bool		hadDigitInExp	= false;

	//	First char can be -
	//
	if (str[i]=='-')
		i++;
		
	while (str[i])
		{
		if (isdigit(str[i]))
			{
			//	Digits are always OK
			//
			if (hadExp)
				hadDigitInExp = true;
			else
				hadDigit = true;
			}
		else if (str[i] == '.')
			{
			//	One decimal point is allowed and it must be before the exponent
			//
			if (hadExp || hadDecimalPt)	
				return false;
			hadDecimalPt = true;
			}
		else if (str[i] == 'e' || str[i] == 'E')
			{
			//	One e is allowed, but it must be after at least one digit
			//
			if (hadExp || !hadDigit)
				return false;
			hadExp = true;
			}
		else if (str[i] == '-')
			{
			//	Another - is allowed if it is preceded by e
			//
			if (!hadExp || (str[i-1] != 'e' && str[i-1] != 'E')	)
				return false;
			}
		else	
			return false;
		i++;
		}
		
	if (hadExp)
		{
		if (hadDigitInExp)
			return true;
		return false;
		}

	if (hadDigit)
		return true;
	return false;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if stored string can be interpreted as a long integer.
*/
bool NxsString::IsALong() const
	{
	const char *str	= c_str();
	unsigned i		= 0;

	//	First char can be -
	//
	if (str[i]=='-')
		i++;

	if (!isdigit(str[i]))
		return false;

	while (str[i])
		{
		if (!isdigit(str[i]))
			return false;
		i++;
		}

	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the stored string is a non-case-sensitive copy of the argument `s'. Note: will return true if both the
|	stored string and `s' are empty strings.
*/
bool NxsString::EqualsCaseInsensitive(
  const NxsString &s)	/* the comparison string */
  const
	{
	unsigned k;
	unsigned slen = s.size();
	unsigned tlen = size();
	if (slen != tlen)
		return false;

	for (k = 0; k < tlen; k++)
		{
  		if ((char)toupper((*this)[k]) != (char)toupper(s[k]))
			return false;
		}

	return true;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Creates a string representation of the hexadecimal version of the long integer `p'. For example, if `p' equals 123,
|	and if 2 was specified for `nFours', the resulting string would be "7B". If 4 was specified for `nFours', then the
|	resulting string would be "007B". 
*/
NxsString NxsString::ToHex(
  long p,			/* the value to display in hexadecimal */
  unsigned nFours)	/* the number of hexadecimal digits to display */
	{
	NxsString s;
	char decod[] = "0123456789ABCDEF";
	for (int i = nFours - 1; i >= 0 ; i--)
		{
		unsigned long k = (p >> (4*i));
		unsigned long masked = (k & 0x000f);
		s += decod[masked];
		}
	return s;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Checks to see if the stored string begins with upper case letters and, if so, returns all of the contiguous capitalized
|	prefix. If the stored string begins with lower case letters, an empty string is returned.
*/
NxsString NxsString::UpperCasePrefix() const
	{
	NxsString x;
	unsigned i = 0;
	while (i < size() && isupper((*this)[i]))
		x += (*this)[i++];
	return x;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Converts the stored string to an unsigned int using the standard C function strtol, throwing NxsX_NotANumber if the 
|	conversion fails. Returns UINT_MAX if the number is too large to fit in an unsigned (or was a negative number).
*/
unsigned NxsString::ConvertToUnsigned() const
	{
	long l = ConvertToLong();
	if (l < 0 || l >UINT_MAX) 
		return UINT_MAX;
	return static_cast<unsigned> (l);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Converts the stored string to an int using the standard C function strtol, throwing NxsX_NotANumber if the conversion 
|	fails. Returns INT_MAX if the number is too large to fit in an int or -INT_MAX if it is too small.
*/
int NxsString::ConvertToInt() const
	{
	long l = ConvertToLong();
	if (l == LONG_MAX || l > INT_MAX)
		return INT_MAX;
	if (l == -LONG_MAX || l <-INT_MAX)
		return -INT_MAX;
	return static_cast<int> (l);
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Converts the stored string to a long using the standard C function strtol, throwing NxsX_NotANumber if the conversion 
|	fails.
*/
long NxsString::ConvertToLong() const
	{
	if (length() == 0 || !(isdigit(at(0)) || at(0) == '-'))
		throw NxsX_NotANumber();
	const char *b = c_str();
	char *endP;
	long l = strtol(b, &endP, 10);
	if (l == 0 && endP == b)
		throw NxsX_NotANumber();
	return l;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Converts the stored string to a double using the standard C function strtod, throwing NxsX_NotANumber if the conversion
|	fails. Returns DBL_MAX or -DBL_MAX if the number is out of bounds.
*/
double NxsString::ConvertToDouble() const
	{
	if (length() == 0)
		throw NxsX_NotANumber();

	char ch = at(0);
	if (isdigit(ch) || ch == '-' || ch == '.'|| toupper(ch) == 'E')
		{
		const char *b = c_str();
		char *endP;
		double d = strtod(b, &endP);
		if (d == 0.0 && endP == b)
			throw NxsX_NotANumber();
		if (d == HUGE_VAL)
			return DBL_MAX;
		if (d == -HUGE_VAL)
			return -DBL_MAX;
		return d;
		}
	throw NxsX_NotANumber();
#if defined (DEMANDS_UNREACHABLE_RETURN)
	return DBL_MAX;
#endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Transforms the vector of NxsString objects by making them all lower case and then capitalizing the first portion of 
|	them so that the capitalized portion is enough to uniquely specify each. Returns true if the strings are long enough 
|	to uniquely specify each. Horrendously bad algorithm, but shouldn't be called often.
*/
bool SetToShortestAbbreviation(
  NxsStringVector	&strVec,		/* vector of NxsString objects */
  bool 				allowTooShort)	/* */
	{
	NxsStringVector upperCasePortion;
	unsigned i;
	for (i = 0; i < strVec.size(); i++)
		{
		// Change the next string to lower case
		//
		strVec[i].ToLower();

		unsigned prefLen = 0;
		NxsString pref;
		
		if (prefLen >= strVec[i].size())
			return false;
		pref += (char) toupper(strVec[i][prefLen++]);
		bool moreChars = true;

		// Keep adding letters from the current string until pref is unique.
		// Then add this pref to upperCasePortion (vector of previous prefs)
		//
		for (;moreChars;)
			{
			size_t prevInd = 0;
			for (; prevInd < upperCasePortion.size(); prevInd++)
				{
				if (pref == upperCasePortion[prevInd])
					{
					// 	Conflict  - both abbreviations need to grow
					//
					if (prefLen >= strVec[i].size())
						{
						if (allowTooShort)
							{
							if (prefLen < strVec[prevInd].size())
								upperCasePortion[prevInd] += (char) toupper(strVec[prevInd][prefLen]);
							moreChars = false;
							break;
							}
						else
							return false;
						}
					pref += (char) toupper(strVec[i][prefLen]);
					if (prefLen >= strVec[prevInd].size())
						{
						if (allowTooShort)
							{
							prevInd = 0;
							prefLen++;
							break;
							}
						else
							return false;
						}
					upperCasePortion[prevInd] += (char) toupper(strVec[prevInd][prefLen++]);
					prevInd = 0;
					break;
					}
				else
					{
					unsigned j;
					for (j = 0; j < prefLen; j++)
						{
						if (pref[j] != upperCasePortion[prevInd][j])
							break;
						}
					if (j == prefLen)
						{
						//	pref agrees with the first part of another abbreviation, lengthen it.
						//
						if (prefLen >= strVec[i].size())
							{
							if (allowTooShort)
								{
								moreChars = false;
								break;
								}
							else
								return false;
							}
						pref += (char) toupper(strVec[i][prefLen++]);
						break;
						}
					}
				}
			if (prevInd == upperCasePortion.size() || !moreChars)
				{
				// Made it all the way through with no problems, add this 
				// prefix as command i's upper case portion
				//
				upperCasePortion.push_back(pref);
				break;
				}
			}
		}

	for (i = 0; i < strVec.size(); i++)
		{
		for (size_t j = 0; j < upperCasePortion[i].size(); j++)
			strVec[i][j] = upperCasePortion[i][j];
		}

	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns a vector of NxsString objects that match the entire `testStr'.
*/
NxsStringVector GetVecOfPossibleAbbrevMatches(
  const NxsString		&testStr,		/* string to match */
  const NxsStringVector	&possMatches)	/* vector of possible matches */
	{
	NxsStringVector matches;
	for (size_t i = 0; i < possMatches.size(); i++)
		{
		if (testStr.Abbreviates(possMatches[i]))
			matches.push_back(possMatches[i]);
		}
	return matches;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Written to make it easy to initialize a vector of strings. Similar to the perl split function. Converts a string like
|	this -- "A|bro|ken strin|g" -- to a vector of strings with four elements:  "A", "bro", "ken string", and "g".
*/
NxsStringVector BreakPipeSeparatedList(
  const NxsString &strList)	/* the string submitted for splitting */
  	{
	NxsString::const_iterator p = strList.begin();
	NxsString ss;
	NxsStringVector retVec;
	for (;;)
		{
		bool done = (p == strList.end());
		if (done || (*p == '|')) 
			{
			retVec.push_back(ss);
			ss.clear();
			if (done)
				break;
			p++;
			}
		ss += *p;
		p++;
		}
	return retVec;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns true if the Equals comparison function is true for this or any element in the vector `s'.
|   (James B. 23-Jul-2020, moved this here, from nxsstring.h, because the references to NxsStringVector 
|    were a problem when compiling in Visual Studio).
*/
bool NxsString::IsInVector(
	const NxsStringVector& s, /* the vector of NxsString objects to be searched */
	NxsString::CmpEnum mode)	/* the argument passed to the Equals function, which is called for every element in the vector `s' */
	const
{
	for (NxsStringVector::const_iterator sIt = s.begin(); sIt != s.end(); sIt++)
	{
		if (Equals(*sIt, mode))
			return true;
	}
	return false;
}

