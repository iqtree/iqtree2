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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

unsigned int add_two_vectors_cpu ( unsigned int n, int * a, int * b, int * c )
{
	unsigned int i;
	for ( i = 0; i < n; i ++ )
		c[i] = a[i] + b[i]; 
	return ( 1 );
}

#ifdef _USE_GPU

unsigned int add_two_vectors_gpu ( unsigned int n, int * a, int * b, int * c )
{

	int err = -1;

	/* get the GPU id */
	cl_platform_id gpu_id = get_gpu_id( &err );
	if( err ) 
	{	
      		return ( 0 );
	}
	
	/* get the device id */
	cl_device_id dev_id = get_dev_id( gpu_id, &err );
	if( err )
	{	
      		return ( 0 );
	}

	/* create the context using dev_id */
	cl_context context = create_context( dev_id, &err );
	if( err )
	{	
      		return ( 0 );
	}

	/* create a list to hold the commands to be executed by GPU */
	cl_command_queue cmd_queue = create_cmd_queue ( dev_id, context, &err );
	if( err )
	{	
      		return ( 0 );
	}

	/* create a kernel */
	cl_kernel kernel;

	/* load the kernel ``kernel.cl'' with name ``vectors_kernel''*/
	kernel = load_kernel ( "kernel.cl", "vectors_kernel", dev_id, context, &err );
	if( err )
	{	
      		return ( 0 );
	}

	/* run the kernel */
	if( ! ( kernel_launch ( kernel, context, cmd_queue, n, a, b, c ) ) )
		return (0);			

        clReleaseContext ( context );
	clReleaseCommandQueue ( cmd_queue );
        clReleaseKernel ( kernel );

	return ( 1 );
 }


unsigned int kernel_launch ( cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, unsigned int n, int * a, int * b, int * c )
{
        int error = 1;

	int * aVec = calloc( n , sizeof ( int ) );
	int * bVec = calloc( n , sizeof ( int ) );
        int * cVec = calloc( n , sizeof ( int ) ); 

	cl_int err;	

	if( aVec == NULL   || bVec == NULL      || cVec == NULL )
	{	
      		return ( 0 );
	}

        /* Here it is not needed */
	fill_aVec ( n, a, aVec );
	fill_bVec ( n, b, bVec );

	cl_mem aVec_device = malloc_device (context, n * sizeof( int ), &error);
	if(error)
	{
      		return ( 0 );	
	}
	
	init_device_mem_int (context, cmd_queue, aVec_device, aVec, n, &error);
	if(error)
	{
      		return ( 0 );	
	}
	
	cl_mem bVec_device = malloc_device (context, n * sizeof( int ), &error);
	if(error)
	{
      		return ( 0 );	
	}

	init_device_mem_int (context, cmd_queue, bVec_device, bVec, n, &error);
	if(error)
	{
      		return ( 0 );	
	}

	cl_mem cVec_device = malloc_device (context, n * sizeof( int ), &error);
	if(error)
	{
      		return ( 0 );	
	}

	init_device_mem_int (context, cmd_queue, cVec_device, cVec, n, &error);
	if(error)
	{
      		return ( 0 );	
	}

	err = clFinish( cmd_queue );
	if( err != CL_SUCCESS )
	{
      		return ( 0 );	
	}

	set_kernel_arguments ( kernel, cmd_queue, aVec_device, bVec_device, cVec_device );

	err = clFinish( cmd_queue );
	if( err != CL_SUCCESS )
	{
      		return ( 0 );	
	}	

	size_t WorkSizeGlobal[] = {n};
	size_t WorkSizeLocal[]  = {1};

	err = clEnqueueNDRangeKernel( cmd_queue, kernel, 1, NULL, WorkSizeGlobal, WorkSizeLocal, 0, NULL, NULL);
	if( err != CL_SUCCESS )
	{
      		return ( 0 );	
	}

	err=clFinish( cmd_queue );
	if( err != CL_SUCCESS )
	{
      		return ( 0 );	
	}	

	read_device_mem_int (cmd_queue, n, c, cVec_device, &error);
	if( error )
	{
      		return ( 0 );	
	}

	/*Here c should contain the result */
	free (aVec);
	free (bVec);
	free (cVec);

	clReleaseMemObject(aVec_device);
	clReleaseMemObject(bVec_device);
	clReleaseMemObject(cVec_device);

	return ( 1 );
}

cl_platform_id get_gpu_id(int * error)
{
	cl_platform_id gpu_id = NULL;

	cl_uint platforms_total=0;
	
	cl_int err;

	err = clGetPlatformIDs (0, NULL, &platforms_total);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}
	if(platforms_total<=0)
	{
		*error = 1;
		return NULL;
	}
	cl_platform_id * gpu_id_vec = malloc(sizeof(cl_platform_id) * platforms_total);

	err = clGetPlatformIDs (platforms_total, gpu_id_vec, NULL);

   // printf( "platforms: %d\n", platforms_total );
    int i;
    int use_platform = -1;
    
    
    // choose the nvidia platform
    for( i = 0; i < platforms_total; ++i ) {
        char str[256];
        
        clGetPlatformInfo( gpu_id_vec[i], CL_PLATFORM_VENDOR, sizeof(str), str, NULL );
        
     //   printf( "vendor: %s\n", str );
        
        if( strstr( str, "NVIDIA" ) != NULL ) {
            use_platform = i;
            break;
        }
    }
    
	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return NULL;
	}
	
	if( use_platform == -1 ) {
        *error = 1;
        return NULL;
    }
	gpu_id = gpu_id_vec[use_platform];

	free( gpu_id_vec );
	*error = 0;

	return gpu_id;
}

cl_device_id get_dev_id(cl_platform_id  gpu_id, int * error)
{
	cl_device_id dev_id = NULL;

	cl_uint devices_total;

	cl_int err;
	
	err = clGetDeviceIDs(gpu_id, CL_DEVICE_TYPE_GPU, 0, NULL, &devices_total);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	if(devices_total<=0)
	{
		*error = 1;
		return NULL;
	}

	cl_device_id * dev_id_vec = malloc(sizeof(cl_device_id) * devices_total);

	err = clGetDeviceIDs(gpu_id, CL_DEVICE_TYPE_ALL, devices_total, dev_id_vec, NULL);

    
    
    
	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	dev_id = dev_id_vec[0];

	free( dev_id_vec );
	*error = 0;
	return dev_id;
}

cl_context create_context(cl_device_id dev_id, int * error)
{


	cl_context context;

	cl_int err;
    
	context = clCreateContext (0,1, &dev_id, NULL,NULL, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}
	
	*error = 0;
	return context;
}

cl_command_queue create_cmd_queue (cl_device_id dev_id, cl_context context, int * error)
{
	cl_int err;

	cl_command_queue cmd_queue;

	cmd_queue = clCreateCommandQueue(context, dev_id, 0, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	return cmd_queue;
}

cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id dev_id, cl_context context, int * error)
{
	cl_kernel kernel;

	FILE * fp = fopen(name, "r");

	if(fp==NULL)
	{
		*error = 1;
		return NULL;
	}

	char * source;
	size_t size;

	fseek(fp, 0, SEEK_END);
	size = ftell(fp);
	fseek(fp, 0, SEEK_SET);
		
	source = malloc((size+1)*sizeof(char));
	
	size = fread(source, 1, size, fp);

	fclose(fp);
	
	source[size] = '\0';

	cl_int err;

	cl_program program = clCreateProgramWithSource (context, 1, (const char **) &source, &size, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

	cl_build_status status;

	clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);

	if(status!=CL_BUILD_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	kernel = clCreateKernel(program, kernel_name , &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	*error = 0;
        free ( source );
	clReleaseProgram(program);
	return kernel;
}

cl_mem malloc_device (cl_context context, size_t size, int * error)
{
	cl_mem mem = NULL;

	cl_int err;

	mem = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &err);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return NULL;
	}

	*error=0;
	return mem;	
}

void init_device_mem_int (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, int * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(int), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

void init_device_mem_uint (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, unsigned int * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(unsigned int), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

void init_device_mem_float (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, float * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(float), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

void set_kernel_arguments (cl_kernel kernel, cl_command_queue cmd_queue, cl_mem cl_mem0, cl_mem cl_mem1, cl_mem cl_mem2 )
{
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &cl_mem0);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &cl_mem1);
	clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &cl_mem2);
	clFinish(cmd_queue);
}

void read_device_mem_int (cl_command_queue cmd_queue, size_t size, int * mem, cl_mem dev_mem, int * error)
{
	cl_int err;

	err = clEnqueueReadBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof( int ), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}
	
	*error = 0;
	return;
}

void fill_aVec ( unsigned int n, int * a, int * aVec )
{
	unsigned int i;
	for( i = 0; i < n; i++ )
		aVec[i] = a[i];
}

void fill_bVec ( unsigned int n, int * b, int * bVec )
{
	unsigned int i;
	for( i = 0; i < n; i++ )
		bVec[i] = b[i];
}

#endif
