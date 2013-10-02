//
// Copyright 2012 Jeff Verkoeyen
// Originally ported from https://github.com/ingenuitas/python-tesseract/blob/master/fmemopen.c
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#if defined __APPLE__ || defined __MACH__


/*--------------------------------------------------------------*/
/* portable version for fmemopen for MAC OSX */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

struct fmem {
  size_t pos;
  size_t size;
  char *buffer;
};
typedef struct fmem fmem_t;

static int readfn(void *handler, char *buf, int size) {
  fmem_t *mem = handler;
  size_t available = mem->size - mem->pos;
  
  if (size > available) {
    size = available;
  }
  memcpy(buf, mem->buffer + mem->pos, sizeof(char) * size);
  mem->pos += size;
  
  return size;
}

static int writefn(void *handler, const char *buf, int size) {
  fmem_t *mem = handler;
  size_t available = mem->size - mem->pos;

  if (size > available) {
    size = available;
  }
  memcpy(mem->buffer + mem->pos, buf, sizeof(char) * size);
  mem->pos += size;

  return size;
}

static fpos_t seekfn(void *handler, fpos_t offset, int whence) {
  size_t pos;
  fmem_t *mem = handler;

  switch (whence) {
    case SEEK_SET: pos = offset; break;
    case SEEK_CUR: pos = mem->pos + offset; break;
    case SEEK_END: pos = mem->size + offset; break;
    default: return -1;
  }

  if (pos > mem->size) {
    return -1;
  }

  mem->pos = pos;
  return (fpos_t)pos;
}

static int closefn(void *handler) {
  free(handler);
  return 0;
}

FILE *fmemopen(void *buf, size_t size, const char *mode) {
  // This data is released on fclose.
  fmem_t* mem = (fmem_t *) malloc(sizeof(fmem_t));

  // Zero-out the structure.
  memset(mem, 0, sizeof(fmem_t));

  mem->size = size;
  mem->buffer = buf;

  // funopen's man page: https://developer.apple.com/library/mac/#documentation/Darwin/Reference/ManPages/man3/funopen.3.html
  return funopen(mem, readfn, writefn, seekfn, closefn);
}

#elif defined WIN32 || defined _WIN32 || defined __WIN32__

/*--------------------------------------------------------------*/
/* SLOW portable version for fmemopen for WIN32 (using temp file) */
/*--------------------------------------------------------------*/
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *fmemopen(void *buf, size_t size, const char *mode) {
    char temppath[MAX_PATH - 13];
    if (0 == GetTempPath(sizeof(temppath), temppath)) {
                puts("Can't get temp path");
        return NULL;
        }
    char filename[MAX_PATH + 1];
    if (0 == GetTempFileName(temppath, "IQT", 0, filename)) {
        puts("Can't get file name");
        return NULL;
    }
        printf("file::%s\n",filename);
    /* FILE *f = fopen(filename, "wb");
      if (NULL == f)
        return NULL;
        */
    FILE *f;
    errno_t err;

        if( (err  = fopen_s( &f, filename, "wb" )) !=0 )
      printf( "The file '%s' was not opened\n", filename );

    fwrite(buf, size, 1, f);
    fclose(f);



    /* return fopen(filename, mode); */
        FILE *f2;
        if( (err  = fopen_s( &f2, filename, mode )) !=0 )
                return f2;


}


#endif
