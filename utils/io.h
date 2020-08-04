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

#ifndef io_h
#define io_h
#include <string>

template <class S> S& safeGetLine(S& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    typename S::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

template <class S> S& safeGetTrimmedLine(S& is, std::string& line) {
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

template <class S> S& safeGetTrimmedLineAsStream(S& is, std::stringstream& lineStream) {
    std::string lineString;
    safeGetTrimmedLine(is, lineString);
    lineStream.str(lineString);
    return is;
}

#endif /* io_h */
