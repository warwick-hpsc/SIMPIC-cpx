#include <iostream>
#include <cuda.h>
#include<cuda_runtime.h>                                                                                                                        
#include<cusparse_v2.h> 

/*cuSPARSE status and handle definition*/
cusparseStatus_t status;
cusparseHandle_t handle=0;

/*CUDA error wrapper*/
static void CUDA_ERROR( cudaError_t err)
{
  if(err!= cudaSuccess){
    printf("CUDA ERROR: %s, exiting\n", cudaGetErrorString(err));
    exit(-1);
    }
}

/*Right and left hand side voltages*/
Scalar rhsV(Scalar t)
{
  return 0.;
}
Scalar lhsV(Scalar t)
{
  return lhsvoltage;
}


/*################################################################/*
 * Fields initializers *
 * ###############################################################*/

/*Array allocation in CPU and GPU*/
void allocate_arrays(int dealloc=0)
{
  if(dealloc)
    {
      /*Delete CPU arrays*/
      delete [] narray;
      delete [] phiarray;
      delete [] Earray;
      
      /*Free GPU memory*/
      cudaFree(Earray_gpu);
      cudaFree(narray_gpu);
      cudaFree(nc_gpu);
      cudaFree(d_a);
      cudaFree(d_b);
      cudaFree(d_c);
      return;
    }
  
  /*Allocate GPU matrix arrays*/ 
  cudaMalloc((void**)&Earray_gpu, ng*sizeof(double));
  cudaMalloc((void**)&narray_gpu, ng*sizeof(double));
  cudaMalloc((void**)&nc_gpu, sizeof(int));
  cudaMalloc((void **)&d_a, ng*sizeof(double));
  cudaMalloc((void **)&d_b, ng*sizeof(double));
  cudaMalloc((void **)&d_c, ng*sizeof(double));
  
  /*Create CPU arrays*/  
  narray = new double[ng];
  phiarray = new double[ng];
  Earray = new double[ng];
}

/*Sets the matrix coefficients in each thread*/
__global__ void setcoeffs(Scalar scale, double *a,double *b,double *c, int *nc_gpu)
{
  /*the id of the thread*/
  int tid=blockIdx.x*blockDim.x+threadIdx.x;           

  /*Each thread writes three values depending on its tid.
   a is lower diagonal, b diagonal, c upper diagonal*/
  if(tid == 0 )
  {
      a[tid] = 0.0;
      b[tid] = 1.0;
      c[tid] = 0.0;
	}
  else if (tid == *nc_gpu)
  {
      a[tid] = 0.0;
      b[tid] = 1.0;
      c[tid] = 0.0;
  } 
  else if( (tid>0)&&(tid<*nc_gpu))
  {
      a[tid] = scale;
      b[tid] = -2.*scale;
      c[tid] = scale;
  }
}

/*Initialization of all variables and arrays*/
void init_fields()
{
  /*Array allocation*/
  allocate_arrays();
  
  /*Setting initial values*/
  cudaMemcpy(nc_gpu, &nc, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemset(Earray_gpu, 0, ng*sizeof(double));
  
  /*Definition of GPU threads and blocks for the correct
   behaviour of the field solver. Nthreads_fields may vary 
   with hardware*/
  Nblocks_fields = ng/Nthreads_max +1;
  /*Set matrix coefficients*/
  setcoeffs<<<Nblocks_fields,Nthreads_max>>>(-epsilon/(q*dx*dx), d_a, d_b, d_c, nc_gpu);

  /*Create a cuSPARSE handle. Needed to call sparse functions*/
  status=cusparseCreate(&handle);
}

/*################################################################/*
 * Laplace and gradient functions *
 * ###############################################################*/

/*Sums the Laplace equation for rhs and lhs voltages*/
__global__ void sumLaplace_GPU(double *pphh, Scalar dx, Scalar rv, Scalar lv, Scalar xl, Scalar L, int ng)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;           
  Scalar frac, xlocal;

  if( (tid >= 0) && ( tid < ng)){
    xlocal = xl + tid*dx;
    frac = xlocal/L;
    pphh[tid] += (frac*rv + (1. - frac)*lv) ;
  }
}

/*Computes the gradient in GPU*/
__global__ void gradient_GPU(double *grad, double *arr, int n, Scalar scale)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;           

  /*Intermediate points calculation*/
  if( (tid>0) && (tid < n))
  {
    grad[tid] = scale*(arr[tid+1] - arr[tid-1]);
  }
  
  /*Fix up end points*/
  if(tid==0)
  {
    grad[0] = 2*scale*(arr[1] - arr[0]);
  }
  if(tid==n)
  {
    grad[n] = 2*scale*(arr[n] - arr[n-1]);
  }
}

/*################################################################/*
 * Main function *
 * ###############################################################*/

/*Updates phi and E values in dt*/
void advancefields(Scalar ddt)
{
  starttime(FIELDS);

  /*Correct density array for computation*/
  starttime(MEMORY_GPU);
  cudaMemset(&(narray_gpu[0]), 0, sizeof(double)); 
  cudaMemset(&(narray_gpu[nc]), 0, sizeof(double)); 
  endtime(MEMORY_GPU);

  starttime(PERFORMANCE);
  /*Trimatrix solver
   *Function from cuSPARSE for trimatrix solution*
   *Overwrites n_array with phi values*/
  status=cusparseDgtsv(handle,nc,1,d_a,d_b,d_c,narray_gpu,nc);
  if (status != CUSPARSE_STATUS_SUCCESS)
  {
     std::cout << status << std::endl;
  }  
  endtime(PERFORMANCE);

  /*Laplace correction and gradient calculation with values
   * that are already in GPU*/
  sumLaplace_GPU<<<Nblocks_fields,Nthreads_max>>>(narray_gpu, dx, rhsV(t), lhsV(t), xl, L, ng);
  gradient_GPU<<<Nblocks_fields,Nthreads_max>>>( Earray_gpu,  narray_gpu,  nc,  -0.5/dx);
  
  /*Copy data back to CPU*/
  if(diag_flag){
    starttime(MEMORY_GPU);
    cudaMemcpy(phiarray, narray_gpu, ng*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(Earray, Earray_gpu,ng*sizeof(double),cudaMemcpyDeviceToHost);
    endtime(MEMORY_GPU);
  }
  
  endtime(FIELDS);
}
