#include <starpu.h>
#include <cuda.h>
//#include<cusparse_v2.h> 

#define LAPLACE

TriMatrix *A;
void advancefields(void *buffers[], void *_args);
void advancefields_GPU(void *buffers[], void*_args);
void solveTridiagonal(void *buffers[], void *_cl_args);
void setcoeffs(Scalar scale);
 //double* d_a;
 //double* d_b;
 //double* d_c;
 /*cuSPARSE status and handle definition*/
 //cusparseStatus_t status;

struct cl_args_f {
  Scalar qm;
  Scalar dt;
  Scalar dx;
  Scalar xl;
  Scalar xr;
  int nl, nr;
  int nlp, nrp; // counters for sending particle buffers
  int nc;
  TriMatrix* A;
  double* d_a;
  double* d_b;
  double* d_c;

} cl_args_fields;

struct cl_args_s {
  int nc;
  TriMatrix* A;
  //double* d_a;
  //double* d_b;
  //double* d_c;

} cl_args_matrix_solver;


/*###########################################
 * #######   Initialisers ##############
 *########################################### */

void allocate_arrays(int dealloc=0)
{
  if(dealloc)
    {
      starpu_data_unregister(narray_handle);
      starpu_data_unregister(Earray_handle);
      starpu_data_unregister(phiarray_handle);
      return;
    }

  starpu_malloc((void**)&narray, ng*sizeof(double));	
  starpu_malloc((void**)&phiarray, ng*sizeof(double));	
  starpu_malloc((void**)&Earray,ng*sizeof(double));	
  starpu_malloc((void**)&nback,ng*sizeof(double));	


  starpu_vector_data_register(&narray_handle, STARPU_MAIN_RAM, (uintptr_t)narray, ng, sizeof(narray[0]));
  starpu_vector_data_register(&phiarray_handle, STARPU_MAIN_RAM, (uintptr_t)phiarray, ng, sizeof(phiarray[0]));
  starpu_vector_data_register(&Earray_handle, STARPU_MAIN_RAM, (uintptr_t)Earray, ng, sizeof(Earray[0]));

  task->handles[4]=narray_handle;
  task->handles[5]=Earray_handle;

  static struct starpu_codelet fields_cl =
 {
   //.where = STARPU_CPU | STARPU_CUDA,
   //.cpu_funcs = { advancefields },
   .where = STARPU_CUDA,
   .cuda_funcs = { advancefields_GPU },
   //.cpu_funcs_name = { "advancefields" },
   .nbuffers = 3,
   .modes = {STARPU_RW, STARPU_RW, STARPU_RW}
  };

  static struct starpu_codelet solver_cl=
  {
    .where = STARPU_CPU,
    .cpu_funcs = {solveTridiagonal},
    .nbuffers = 2,
    .modes = {STARPU_RW, STARPU_RW}
  };

  task_fields = starpu_task_create();
  task_fields->cl = &fields_cl;
  task_fields->handles[0] = narray_handle;
  task_fields->handles[1] = phiarray_handle;
  task_fields->handles[2] = Earray_handle;
  task_fields->cl_arg = &cl_args_fields;
  task_fields->cl_arg_size = sizeof(cl_args_f);
  task_fields->synchronous = 1;
  task_fields->destroy = 0;
  
  task_solver = starpu_task_create();
  task_solver->cl = &solver_cl;
  task_solver->handles[0] = narray_handle;
  task_solver->handles[1] = phiarray_handle;
  task_solver->cl_arg = &cl_args_matrix_solver;
  task_solver->cl_arg_size = sizeof(cl_args_s);
  task_solver->synchronous = 1;
  task_solver->destroy = 0;
  //starpu_cusparse_init() ;
}

void allocate_fieldbuffers(int dealloc=0)
{
  if(dealloc)
    {
      delete [] fsbuffer;
      delete [] frbuffer;
      return;
    }

  fsbuffer = new Scalar[2*nproc];
  frbuffer = new Scalar[2*nproc];
}


void init_fields()
{
  allocate_arrays();
  allocate_fieldbuffers();

  memset(nback, 0, ng*sizeof(Scalar));

  setcoeffs(-epsilon/(q*cl_args.dx*cl_args.dx));


  cl_args_fields.A = A;
  cl_args_matrix_solver.A = A;
  cl_args_matrix_solver.nc = nc;

  cl_args_fields.d_a = A->tri_a;
  cl_args_fields.d_b = A->tri_b;
  cl_args_fields.d_c = A->tri_c;

}

/*###########################################
 * #######  General functions ##############
 *########################################### */


Scalar rhsV(Scalar t)
{
  return 0.;
}

Scalar lhsV(Scalar t)
{
  return lhsvoltage;
}


void setcoeffs(Scalar scale)
{
  A = new TriMatrix(ng);

  // lhs BC
  if(rank == 0)
    {
      A->a(0) = 0.0;
      A->b(0) = 1.0;
      A->c(0) = 0.0;
    }
  else
    {
      //      A->a(0) = ;
      A->b(0) = -scale*(1.0 + cl_args.dx/cl_args.xl);
      A->c(0) = scale;

    }

  // rhs BC
  if(rank == last)
    {
      A->a(nc) = 0.0;
      A->b(nc) = 1.0;
      A->c(nc) = 0.0;
    }
  else
    {
      A->a(nc) = scale;
      A->b(nc) = -scale*(1.0 + cl_args.dx/(L - cl_args.xr));
    }

  for(int i=1; i < nc; i++)
    {
      A->a(i) = scale;
      A->b(i) = -2.*scale;
      A->c(i) = scale;
    }

}


/*###########################################
 * #######   CPU functions ##############
 *########################################### */

void sumLaplace(Scalar *pphh)
{
  register int i;
  Scalar rv, lv, frac,  xlocal;

  rv = rhsV(t);
  lv = lhsV(t);

  for(i=0, xlocal = cl_args.xl; i < ng; i++, xlocal +=cl_args.dx)
    {
      frac = xlocal/L;
      pphh[i] += frac*rv + (1. - frac)*lv;
    }
}


void communicate_fields(void)
{
  int i, j, k;
  Scalar lcoord, rcoord;
  Scalar xval, phival;

  Scalar ttemp[2];

#ifdef DEBUG
  //  fprintf(stderr, "proc %d entering communicate_fields\n", rank);
#endif

  if(nproc > 1)
    {
      ttemp[0] = phiarray[0];
      ttemp[1] = phiarray[nc];

      MPI_Allgather(ttemp, 2, MPI_SCALAR, frbuffer, 2, MPI_SCALAR, MPI_COMM_WORLD);


      for(i=0, j=0; i < nproc; i++, j+=2)
	{
#ifdef DEBUG
	  //	  fprintf(stderr, "rank %d %d: %g %g;\n",
	  //		  rank, i, frbuffer[j], frbuffer[j+1]);
#endif
	}

      // compute contribution from lhs
      for(i=0, j=0; i < rank; i++, j+=2)
	{
	  phival = frbuffer[j+1];
	  lcoord = L*(i+1)/Scalar(nproc);
	  lcoord = L - lcoord;
	  for(k=0, xval=L-cl_args.xl; k < ng; k++, xval -= cl_args.dx)
	    {
	      phiarray[k] += phival*(xval/lcoord);
	    }
	}

      // compute contribution from rhs
      for(i=rank+1; i < nproc; i++)
	{
	  j = 2*i;
	  phival = frbuffer[j];
	  rcoord = L*i/Scalar(nproc);
	  for(k=0, xval=cl_args.xl; k < ng; k++, xval += cl_args.dx)
	    {
	      phiarray[k] += phival*(xval/rcoord);
	    }
	}

    }

#ifdef DEBUG
  //      fprintf(stderr, "proc %d done with communicate_fields\n", rank);
#endif
}

void advancefields(void *buffers[], void *_cl_args)
{
  Scalar nlold, nrold;

  Scalar *narray = (Scalar *)STARPU_VECTOR_GET_PTR(buffers[0]);
  Scalar *phiarray = (Scalar *)STARPU_VECTOR_GET_PTR(buffers[1]);
  Scalar *Earray = (Scalar *)STARPU_VECTOR_GET_PTR(buffers[2]);

  struct cl_args_t *cl_args = (struct cl_args_t *)_cl_args;


  // modify density array
  if(rank == 0)
    {
      nlold = narray[0];
      narray[0] = 0.;
    }
  else
    {
      narray[0] *= 2;
    }
  if(rank == last)
    {
      nrold = narray[nc];
      narray[nc] = 0.;
    }
  else
    {
      narray[nc] *= 2;
    }

  A->solve(narray, phiarray);

  // restore density array
  if(rank == 0)
    {
      narray[0] = 2*nlold;
    }
  if(rank == last)
    {
      narray[nc] = 2*nrold;
    }

  starttime(MPITIME);
  communicate_fields();
  endtime(MPITIME);

  sumLaplace(phiarray);
  gradient(Earray, phiarray, ng, -0.5/cl_args->dx);

  // "fix up" endpoints
  Earray[0] = -(phiarray[1] - phiarray[0])/cl_args->dx;
  Earray[nc] = -(phiarray[nc] - phiarray[nc-1])/cl_args->dx;


}

/*###########################################
 * #######   GPU functions ##############
 *########################################### */

/*CPU tridiagonal solver*/
void solveTridiagonal(void *buffers[], void *_cl_args)
{
  struct cl_args_s *cl_args = (struct cl_args_s*)_cl_args;


  Scalar nlold, nrold;

  double *narray = (double *)STARPU_VECTOR_GET_PTR(buffers[0]);
  double *phiarray = (double *)STARPU_VECTOR_GET_PTR(buffers[1]);

  A = cl_args->A;
  nlold = narray[0];
  nrold = narray[cl_args->nc];
  narray[0] = 0.;
  narray[cl_args->nc] = 0.;
  A->solve(narray, phiarray);
  narray[0] = 2.*nlold;
  narray[cl_args->nc] = 2.*nrold;

  std::cout<<"Tridiagonal solver"<<std::endl;
}


/*Sums the Laplace equation for rhs and lhs voltages*/
__global__ void sumLaplace_GPU(double *pphh, Scalar dx, Scalar rv, Scalar lv, Scalar xl, Scalar L, int ng)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;           
  Scalar frac, xlocal;
  /*
  if(tid ==0){
  for(int i = 0; i<ng; i++){
   printf("%g , ", pphh[i]);
    }}
    */
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

__global__ void calculateE(double *pphh, double *grad, Scalar scale, Scalar dx, Scalar rv, Scalar lv, Scalar xl, Scalar L, int ng)
{
  int tid=blockIdx.x*blockDim.x+threadIdx.x;           
  Scalar frac, xlocal;
  
  if( (tid >= 0) && ( tid < ng)){
    xlocal = xl + tid*dx;
    frac = xlocal/L;
    pphh[tid] += (frac*rv + (1. - frac)*lv) ;
    if (tid==0){
      atomicAdd_double(&(grad[tid]), -2*scale*pphh[tid]);
      atomicAdd_double(&(grad[tid+1]), -scale*pphh[tid]);
    }
    else if (tid==ng -1){
      atomicAdd_double(&(grad[tid]), 2*scale*pphh[tid]);
      atomicAdd_double(&(grad[tid-1]), scale*pphh[tid]);
    }
    else if (tid == 1){
      atomicAdd_double(&(grad[tid+1]), -scale*pphh[tid]);
      atomicAdd_double(&(grad[tid-1]), 2*scale*pphh[tid]);
      }
    else if (tid == ng-2){
      atomicAdd_double(&(grad[tid+1]), -2*scale*pphh[tid]);
      atomicAdd_double(&(grad[tid-1]), scale*pphh[tid]);
    }
      else{
      atomicAdd_double(&(grad[tid+1]), -scale*pphh[tid]);
      atomicAdd_double(&(grad[tid-1]), scale*pphh[tid]);
    }
  }
}


/*Updates phi and E values in dt*/
void advancefields_GPU(void *buffers[], void *_cl_args)
{

  struct cl_args_f *cl_args = (struct cl_args_f*)_cl_args;
  
  Scalar nlold, nrold;

  std::cout<<"runs GPU fields"<<std::endl;

  double *narray = (double *)STARPU_VECTOR_GET_PTR(buffers[0]);
  double *phiarray = (double *)STARPU_VECTOR_GET_PTR(buffers[1]);
  double *Earray = (double *)STARPU_VECTOR_GET_PTR(buffers[2]);

  /*Correct density array for computation*/

  
  A = cl_args->A;

  double narray_cpu[ng];
  double phiarray_cpu[ng];
  cudaMemcpy(phiarray_cpu, phiarray, ng*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(narray_cpu, narray, ng*sizeof(double),cudaMemcpyDeviceToHost);
  nlold = narray_cpu[0];
  nrold = narray_cpu[cl_args->nc];
  narray_cpu[0] = 0.;
  narray_cpu[cl_args->nc] = 0.;
  A->solve(narray_cpu, phiarray_cpu);
  narray_cpu[0] = 2.*nlold;
  narray_cpu[cl_args->nc] = 2.*nrold;
  cudaMemcpy(phiarray, phiarray_cpu, ng*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(narray, narray_cpu, ng*sizeof(double),cudaMemcpyHostToDevice);


  /*Need to figure out how to call cursparse from starPU
   * Raises [starpu][execute_job_on_cuda][assert failure] Unless when using the STARPU_CUDA_ASYNC flag, CUDA codelets have to wait for termination of their kernels on the starpu_cuda_get_local_stream() stream
   * a.out: ../../src/drivers/cuda/driver_cuda.c:663: execute_job_on_cuda: Assertion `cudaStreamQuery(starpu_cuda_get_local_stream()) == cudaSuccess' failed.
  */
  /*
  status = cusparseDgtsv(starpu_cusparse_get_local_handle(), cl_args->nc ,1, cl_args->d_a , cl_args->d_b, cl_args->d_c, narray, cl_args->nc);

  std::cout<<"runs GPU after cusparse"<<std::endl;
  cudaStreamSynchronize(starpu_cuda_get_local_stream());
  if (status != CUSPARSE_STATUS_SUCCESS)
  {
     std::cout << status << std::endl;
  }  
  /*

  /*Laplace correction and gradient calculation with values
   * that are already in GPU*/

  int Nthreads_max = prop.maxThreadsPerBlock;
  int Nblocks_fields = ng/Nthreads_max +1;
  
  sumLaplace_GPU<<<Nblocks_fields,Nthreads_max, 0, starpu_cuda_get_local_stream()>>>(phiarray, cl_args->dx, rhsV(t), lhsV(t), cl_args->xl, cl_args->xr - cl_args->xl, cl_args->nc +1);
  cudaStreamSynchronize(starpu_cuda_get_local_stream());

  gradient_GPU<<<Nblocks_fields,Nthreads_max, 0, starpu_cuda_get_local_stream()>>>( Earray,  phiarray,  cl_args->nc,  -0.5/cl_args->dx);
  //cudaMemset(Earray, 0, ng*sizeof(Scalar));
  //calculateE<<<Nblocks_fields,Nthreads_max, 0, starpu_cuda_get_local_stream()>>>(phiarray, Earray, -0.5/cl_args->dx, cl_args->dx, rhsV(t), lhsV(t), cl_args->xl, cl_args->xr - cl_args->xl, cl_args->nc +1);
  cudaStreamSynchronize(starpu_cuda_get_local_stream());

}
