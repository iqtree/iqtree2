/**
    Computing with GPUs template
    Copyright (C) 2012 Nikolaos Alachiotis, Fernando Izquierdo, and 
    Solon P. Pissis.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifdef _USE_GPU
#include "CL/opencl.h"
#endif


unsigned int add_two_vectors_cpu ( unsigned int n, int * a, int * b, int * c );

#ifdef _USE_GPU

unsigned int add_two_vectors_gpu ( unsigned int n, int * a, int * b, int * c );
unsigned int kernel_launch ( cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, unsigned int n, int * a, int * b, int * c );
cl_platform_id get_gpu_id(int * error);
cl_device_id get_dev_id(cl_platform_id  gpu_id, int * error);
cl_context create_context(cl_device_id dev_id, int * error);
cl_command_queue create_cmd_queue (cl_device_id dev_id, cl_context context, int * error);
cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id dev_id, cl_context context, int * error);
cl_mem malloc_device (cl_context context, size_t size, int * error);
void init_device_mem_int (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, int * mem, size_t size, int * error);
void init_device_mem_uint (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, unsigned int * mem, size_t size, int * error);
void init_device_mem_float (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, float * mem, size_t size, int * error);
void set_kernel_arguments (cl_kernel kernel, cl_command_queue cmd_queue, cl_mem cl_mem0, cl_mem cl_mem1, cl_mem cl_mem2 );
void read_device_mem_int (cl_command_queue cmd_queue, size_t size, int * mem, cl_mem dev_mem, int * error);
void fill_aVec ( unsigned int n, int * a, int * aVec );
void fill_bVec ( unsigned int n, int * b, int * bVec );

#endif
