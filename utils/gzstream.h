// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.h
// Revision      : $Revision: 1.5 $
// Revision_date : $Date: 2002/04/26 23:30:15 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#ifndef GZSTREAM_H
#define GZSTREAM_H 1

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
//#include "zlib-1.2.7/zlib.h"
#include <zlib.h>
#include "progress.h"

#define GZ_NO_COMPRESSION (1L << 11)

#ifdef GZSTREAM_NAMESPACE
namespace GZSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See below for user classes.
// ----------------------------------------------------------------------------

class gzstreambuf : public std::streambuf {
	friend class gzstreambase;
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for igzstream at the end.

    gzFile           file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    char             opened;             // open/close state of stream
    int              mode;               // I/O mode

    size_t           compressed_length;
    size_t           compressed_position; //only tracked for read (input) streams
    
    int flush_buffer();
    
    progress_display* progress;
    
public:
    gzstreambuf() : opened(0), compressed_length(0), compressed_position(0), progress(nullptr) {
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position      
        // ASSERT: both input & output capabilities will not be used together
    }
    int is_open() const { return opened; }
    gzstreambuf* open( const char* name, int open_mode, int compression_level=9);
    gzstreambuf* close();
    ~gzstreambuf() { close(); }
    
    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();

    size_t getCompressedLength()   const;
    size_t getCompressedPosition() const;
    void setProgress(progress_display* p);
};

class gzstreambase : virtual public std::ios {
protected:
    gzstreambuf buf;
public:
    gzstreambase() { init(&buf); }
    gzstreambase( const char* name, int open_mode, int compression_level = 9);
    ~gzstreambase();
    void open( const char* name, int open_mode, int compression_level=9);
    void close();
	z_off_t get_raw_bytes(); // BQM: return number of uncompressed bytes

    const gzstreambuf* rdbuf() const { return &buf; }
    gzstreambuf* rdbuf() { return &buf; }
    size_t getCompressedLength() const  {
        return buf.getCompressedLength(); }
    size_t getCompressedPosition() const {
        return buf.getCompressedPosition(); }
};

// ----------------------------------------------------------------------------
// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz* 
// function interface of the zlib. Files are compatible with gzip compression.
// ----------------------------------------------------------------------------

class igzstream : public gzstreambase, public std::istream {
public:
    igzstream() : gzstreambase(), std::istream( &buf) {}
    igzstream( const char* name, int open_mode = std::ios::in)
        : gzstreambase( name, open_mode), std::istream( &buf) {}
    const gzstreambuf* rdbuf() const { return gzstreambase::rdbuf(); }
    gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::in) {
        gzstreambase::open( name, open_mode);
    }
    bool is_open() const {
        const gzstreambuf* buf = rdbuf();
        return buf!=nullptr && buf->is_open();
    }
};

class pigzstream: public igzstream {
protected:
    std::string format_name;
    progress_display* progress;
public:
    typedef igzstream super;
    pigzstream(const char* format);
    void open( const char* name, int open_mode = std::ios::in);
    void close();
    void done();
    
    const gzstreambuf* rdbuf() const;
    gzstreambuf* rdbuf();
    
    void hideProgress();
    void showProgress();
    ~pigzstream();
};

class ogzstream : public gzstreambase, public std::ostream {
public:
    ogzstream() : gzstreambase(), std::ostream( &buf) {
    }
    ogzstream( const char* name, int mode = std::ios::out,
               int compression_level = 9)
        : gzstreambase( name, mode, compression_level)
        , std::ostream( &buf) {
    }
    gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::out,
               int compression_level = 9) {
        gzstreambase::open( name, open_mode, compression_level);
    }
};

#ifdef GZSTREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif

#endif // GZSTREAM_H
// ============================================================================
// EOF //

