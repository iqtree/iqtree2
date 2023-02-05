/***************************************************************************
 *   Copyright (C) 2017-2020 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   James Barbetti <james_barbetti@yahoo.com>                             *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
//
//  safeGetLine (the string version): 
//  created by BUI Quang Minh <minh.bui@univie.ac.at> 18-Mar-2017.
//  This file:   Created by James Barbetti on 03-Aug-2020.
//  safeGetTrimmedLine: James Barbetti 03-Aug-2020.
//

#ifndef safe_io_h
#define safe_io_h
#include <string>  //for std::string
#include <sstream> //for std::stringstream

/**
 * @brief  Read a line from a stream into a string buffer
 * @tparam S  - the type of the input stream (S must have
 *              an S::entry type exposing a rdbuf() member
 *              function returning a std::sreambuf pointer)
 *              (in practice, S will likely be a 
 *               std:ifstream).
 * @param  is - reference to the input stream
 * @param  t  - reference to a std::istringstream to read
 * @return S& a reference to the input file stream
 *              (will be positioned after the linefeed at
 *               the end of the line; or at end of file,
 *               if the file ended before the line did).
 * @note   The content of t is cleared before reading
 *         behins.  This function does not append the 
 *         std::istringstream.
 * @note   Trailing "\n" or "\r\n" will be skipped over, 
 *         and will NOT be read into t.
 * @note   The characters in the stream are read one-by-one 
 *         using a std::streambuf. That is faster than 
 *         reading them one-by-one using the std::istream
 *         (or whatever).
 *         Code that uses streambuf this way must be 
 *         guarded by a sentry object. The sentry object 
 *         performs various tasks, such as thread 
 *         synchronization and updating the stream state.
 */
template <class S=std::istringstream> 
S& safeGetLine(S& is, std::string& t)
{
    t.clear();

    typename S::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n') {
                sb->sbumpc();
            }
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty()) {
                is.setstate(std::ios::eofbit);
            }
            return is;
        default:
            t += (char)c;
        }
    }
}

/**
 * @brief  Get a line from a stream, and trim off leading
 *         and trailing white space.
 * @tparam S - the type of the input stream (see safeGetLine
 *             for a summary of what is required of S)
 * @param  is   - reference to the input file stream
 * @param  line - reference to the string variable into
 *                which the trimmed line is to be read.
 * @return S& a reference to the input stream
 *         (will be positioned after the linefeed at
 *          the end of the line; or at end of file,
 *          if the file ended before the line did). 
 * @note   Implemented in terms of safeGetLine().
 * @note   Leading ' ' and '\t' characters are trimmed.
 *         Trailing ' ', '\t', '\r', and '\n' characters
 *         are trimmed.
 */
template <class S=std::istringstream> 
S& safeGetTrimmedLine(S& is, std::string& line) {
    safeGetLine(is, line);
    size_t start = 0;
    size_t stop  = line.length();
    if (0<stop) {
        while (line[start]==' ' || line[start]=='\t') {
            ++start;
            if (start == stop) {
                break;
            }
        }
    }
    if (start<stop) {
        while (line[stop-1]==' ' || line[stop-1]=='\t'
               || line[stop-1]=='\r' || line[stop-1]=='\n') {
            --stop;
            if (start == stop) {
                break;
            }
        }
    }
    if (start==stop) {
        line.clear();
    } else {
        line = line.substr(start, stop-start);
    }
    return is;
}

/**
 * @brief  Read a line from an input stream and load it up
 *         into a std::stringstream (passed in by reference).
 * @tparam S - the type of the input stream.
 * @param  is - reference to the input file stream
 * @param  lineStream - reference to the output string stream
 * @return S& - reference to the input file stream
 *              (will be positioned after the linefeed at
 *               the end of the line; or at end of file,
 *               if the file ended before the line did).
 */
template <class S=std::istringstream> 
S& safeGetTrimmedLineAsStream
    (S& is, std::stringstream& lineStream) {
    std::string lineString;
    safeGetTrimmedLine<S>(is, lineString);
    lineStream.str(lineString);
    return is;
}

/**
 * @brief  Returns true if a stream's read pointer
 *         is positioned at the end of the stream
 * @tparam S  - the type of the input stream.   
 *         S must support exceptions(), operator>>,
 *         eof() and ungetc() all working just as
 *         std::ifstream's same-named member functions do.
 * @param  in - a reference to the input stream
 * @return true  - if read pointer positioned at end
 * @return false - otherwise
 */
template <class S> bool isAtEndOfFile(S& in) {
    char ch;
    in.exceptions(std::ios::goodbit);
    (in) >> ch;
    if (in.eof()) {
        return true;
    }
    in.unget();
    in.exceptions(std::ios::failbit | std::ios::badbit);
    return false;
}

#endif /* safe_io_h */
