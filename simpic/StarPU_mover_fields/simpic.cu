#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <iostream>

#include "global.h"
#include <mpi.h>
#include <starpu_mpi.h>
#include "timing.cpp"
#include "diagnostics.cpp"
#include "trimatrix.cpp"
#include "fields.cu"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

void moveparticlesGPU(void *buffers[], void *_cl_args);
void moveparticles(void *buffers[], void * _cl_args);

/* CUDA error wraper */
static void CUDA_ERROR( cudaError_t err) 
{
    if (err != cudaSuccess) {
        printf("CUDA ERROR: %s, exiting\n", cudaGetErrorString(err));
        exit(-1);
    }
}

/*###########################################
 * #######  General functions ##############
 *########################################### */
inline void lhs(int j)
{
  lhsbuf[cl_args.nlp] = pdata[j];
  lhsbuf[cl_args.nlp+1] = pdata[j+1];
  cl_args.nlp += 2;

#ifdef DEBUG
  fprintf(stderr, "rank %d nlp => %d\n", rank, cl_args.nlp);
#endif

}

inline void rhs(int j)
{
  rhsbuf[cl_args.nrp] = pdata[j];
  rhsbuf[cl_args.nrp+1] = pdata[j+1];
  cl_args.nrp += 2;

#ifdef DEBUG
  fprintf(stderr, "rank %d nrp => %d\n", rank, cl_args.nrp);
#endif
}



/*###########################################
 * #######  Initialisers ##############
 *########################################### */

void allocate_particles(void)
{
  
  // nproc + 1 -- safety factor for this test
  // (should be more dynamic)
  starpu_malloc((void**)&pdata, (nproc+1)*2*npart*sizeof(Scalar));	
  starpu_malloc((void**)&lhsbuf, (nproc+1)*2*npart*sizeof(Scalar));	
  starpu_malloc((void**)&rhsbuf,(nproc+1)*2*npart*sizeof(Scalar));	
  starpu_malloc((void**)&check_alive,(nproc+1)*npart*sizeof(bool));	
  for(int i=0, j=0; i < npart; i++, j+=2)	
    {	
      pdata[j] = cl_args.xl + ((i + 1)*(cl_args.xr-cl_args.xl))/Scalar(npart + 1);	
      pdata[j+1] = 0.;	
      check_alive[i] = true;	
    }	
 	
  starpu_vector_data_register(&pdata_handle, STARPU_MAIN_RAM,	
			      (uintptr_t)pdata, (nproc+1)*2*npart,	
			      sizeof(pdata[0]));	
  starpu_vector_data_register(&lhsbuf_handle, STARPU_MAIN_RAM,	
			      (uintptr_t)lhsbuf, (nproc+1)*2*npart,	
			      sizeof(lhsbuf[0]));	
  starpu_vector_data_register(&rhsbuf_handle, STARPU_MAIN_RAM,	
			      (uintptr_t)rhsbuf, (nproc+1)*2*npart,	
			      sizeof(rhsbuf[0]));	
  starpu_vector_data_register(&alive_handle, STARPU_MAIN_RAM, (uintptr_t)check_alive,(nproc+1)*npart,sizeof(check_alive[0]));	
  	
  
  static struct starpu_codelet cl = 
    {
      .where = STARPU_CUDA,
      .cuda_funcs = {moveparticlesGPU},
      .nbuffers = 6,
      .modes = {STARPU_RW, STARPU_RW, STARPU_RW, STARPU_RW, STARPU_RW, STARPU_R}
    };
  task = starpu_task_create();
  task->cl = &cl;
  task->handles[0] = pdata_handle;
  task->handles[1] = alive_handle;
  task->handles[2] = lhsbuf_handle;
  task->handles[3] = rhsbuf_handle;
  task->cl_arg = &cl_args;
  task->cl_arg_size = sizeof(cl_args_t);
  task->synchronous = 1;
  task->destroy = 0;

}	

void parsecmdline(int argc, char **argv)
{
  int i;
  int nnppcc;   // number of particles per cell
  int nnccpppp; // number of cells per processor
  int nnnttt;   // number of time steps

  double ddttff, lhsdefault;

  #ifdef DEBUG
  //  fprintf(stderr, "Program name: %s\n", argv[0]);
  #endif

  ddttff = lhsdefault = 0.;
  nnppcc = nnccpppp = nnnttt = 0;

  for(i=1; i < argc; i++)
    {
      #ifdef DEBUG
      //      fprintf(stderr, "arg %d : %s\n", i, argv[i]);
      #endif

      if(strcmp(argv[i], "-ppc") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &nnppcc);
	}
      else if(strcmp(argv[i], "-ncpp") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &nnccpppp);
	}
      else if(strcmp(argv[i], "-nt") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &nnnttt);
	}
      else if(strcmp(argv[i], "-dtfactor") == 0)
	{
	  i++;
	  sscanf(argv[i], "%lf", &ddttff);
	}
      else if(strcmp(argv[i], "-lhsv") == 0)
	{
	  i++;
	  sscanf(argv[i], "%lf", &lhsdefault);
	}
      else // default case
	{
	  fprintf(stderr, "\nError:\n");
	  fprintf(stderr, "Unrecognized argument \"%s\"\n", argv[i]);
	  Quit();
	}
    }

  #ifdef DEBUG
  //  fprintf(stderr, "\ncmd line args parsed: ");
  //  fprintf(stderr, "%d particles/cell, %d cells/processor, %d timesteps\n",
  //	  nnppcc, nnccpppp, nnnttt);
  #endif

  if((nnppcc < 1) || (nnccpppp < 1) || (nnnttt < 1))
    {
      fprintf(stderr, "\nError, input arguments must be entered!\n");
      Quit();
    }

  if(ddttff <= 0.)
    {
      fprintf(stderr, "\nError, dtfactor MUST BE positive!\n");
      Quit();
    }

  // set simulation variables from input data
  ntimesteps = nnnttt;
  nc = nnccpppp;
  ng = nc + 1;
  npart = nc*nnppcc;
  dtfactor = ddttff;
  lhsvoltage = lhsdefault;
}

void init(void)
{

  // remove old .dat files
  system("rm -f *.dat");

  // get MPI data
  MPI_Get_version(&ver, &subver);
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  wtres = MPI_Wtick();

  // print processor string name
  sprintf(procname, "%03d", rank);

  // reset timer variables
  for(int j=0; j < NTIMES; j++)
    {
      tttimes[j] = 0.;
    }
  // label timer counters -- do this here because
  // I can't find any documentation on how to declare
  // multidimensional variables at compile time
  // aka C is stupid
  sprintf(names[PARTICLES], "Particles");
  sprintf(names[FIELDS], "Fields");
  sprintf(names[MPITIME], "MPI");
  sprintf(names[DIAGNOSTICS], "Diagnostics");
  sprintf(names[INIT], "Initialization");

  nparttot = 0.;

  density = 1.E13;
  epsilon = 8.85E-12;
  area = 1.;
  L = 1.;
  q = 1.602E-19;
  m = 9.1E-31;


  cl_args.xl = rank/Scalar(nproc);
  cl_args.xr = (rank + 1.)/Scalar(nproc);
  cl_args.xl *= L;
  cl_args.xr *= L;
  cl_args.nc = nc*nproc;
  Ll = cl_args.xr - cl_args.xl;
  cl_args.dx = Ll / Scalar(nc);

  cl_args.nl = rank - 1;
  cl_args.nr = rank + 1;
  
  cl_args_fields.xl =  cl_args.xl;
  cl_args_fields.xr =  cl_args.xr;
  cl_args_fields.nc =  cl_args.nc;
  cl_args_fields.dx =  cl_args.dx;
  cl_args_fields.qm =  cl_args.qm;

  
  last = nproc - 1;
  if(rank == 0)
    {
      cl_args.nl = last;
    }
  if(rank == last)
    {
      cl_args.nr = 0;
    }

  np2c = density*area*Ll/Scalar(npart);

  allocate_particles();
  init_fields();

  // calculate time step from plasma frequency
  wp = density*q*q/(epsilon*m);
  wp = sqrt(wp);

  cl_args.dt = dtfactor/wp;

  qscale = np2c/(area*cl_args.dx);

  t = ntimesteps*cl_args.dt;
  cl_args.qm = q*cl_args.dt/m;

  diagnosticsflag = true;
}

/*###########################################
 * #######  CPU particle mover functions ##############
 *########################################### */

inline Scalar gridx(Scalar x)
{
  return ((x - cl_args.xl)/cl_args.dx);
}

inline Scalar getfrac(Scalar x)
{
  static int n;

  n = int(x);
  return x-n;
}

inline void weightparticle(Scalar x, Scalar weight, Scalar *array)
{
  static int n;
  static Scalar frac, xx;

  xx = gridx(x);
  n = int(xx);
  frac = xx-n;

  array[n] += weight*(1. - frac);
  array[n+1] += weight*frac;
}

inline Scalar particleweight(Scalar x, Scalar *array)
{
  int n;
  static Scalar frac, xx;

  xx = gridx(x);
  n = int(xx);
  frac = xx-n;

  return frac*array[n+1] + (1.0 - frac)*array[n];
}


void moveparticles(void *buffers[], void * _cl_args)	
{	
  int i, j;	
  Scalar E;	
  Scalar *pdata = (Scalar *) STARPU_VECTOR_GET_PTR(buffers[0]);	
  Scalar *lhsbuf = (Scalar *) STARPU_VECTOR_GET_PTR(buffers[2]);	
  Scalar *rhsbuf = (Scalar *) STARPU_VECTOR_GET_PTR(buffers[3]);	
  bool *alive = (bool *) STARPU_VECTOR_GET_PTR(buffers[1]);	
  Scalar *narray = (Scalar *) STARPU_VECTOR_GET_PTR(buffers[4]);	
  Scalar *Earray = (Scalar *) STARPU_VECTOR_GET_PTR(buffers[5]);	
  struct cl_args_t *cl_args = (struct cl_args_t *) _cl_args;	
  int n_partition = STARPU_VECTOR_GET_NX(buffers[1]);	
  //return;	
  starttime(PARTICLES);	
  for(i=0 ; i < n_partition; i++)	
    {	
      j = 2*i;	
      if(alive[i]==true)	
     { 	
      E = particleweight(pdata[j], Earray);	
      pdata[j+1] += cl_args->qm*E;	
      pdata[j] += pdata[j+1]*cl_args->dt;	
      #ifdef DEBUG	
      //      fprintf(stderr, "particle %d (%d): E=%g, x=%g, v=%g, qm=%g\n", i, j,	
      //	      E, pdata[j], pdata[j+1], qm);	
      #endif	
      if((pdata[j] > cl_args->xl) && (pdata[j] < cl_args->xr))	
	{	
	  weightparticle(pdata[j], qscale, narray);	
	  	
	}	
      else	
	{	
	  #ifdef DEBUG	
	  //  fprintf(stderr, "\nParticle %d[of %d] lost, E=%g\n", i, npart, E);	
	  #endif	
	  alive[i]=false;	
	}	
      }		
    }	
  endtime(PARTICLES);	
}	

void communicateparticles(void)
{
  MPI_Request lstat, rstat;
  MPI_Request plstat, prstat;
  MPI_Request llstat, rrstat;
  MPI_Request pllstat, prrstat;
  MPI_Status stat;

#ifdef DEBUG
  fprintf(stderr,
	  "Rank %d communicateparticles(), t=%g, nlp=%d, nrp=%d, npart=%d\n",
	  rank, tt, cl_args.nlp, cl_args.nrp, npart);

#endif

  // recv/send from/to left
  if(rank != 0)
    {

#ifdef DEBUG
      fprintf(stderr, "left: Processor %d sending %d particles to proc %d\n",
	      rank, cl_args.nlp, cl_args.nl);
#endif

      MPI_Irecv(&nlr, 1, MPI_INT, cl_args.nl, 0, MPI_COMM_WORLD, &lstat);
      MPI_Isend(&cl_args.nlp, 1, MPI_INT, cl_args.nl, 0, MPI_COMM_WORLD, &llstat);
    }

  // recv/send from/to right
  if(rank != last)
    {

#ifdef DEBUG
      fprintf(stderr, "right: Processor %d sending %d particles to proc %d\n",
	      rank, cl_args.nrp, cl_args.nr);
#endif

      MPI_Irecv(&nrr, 1, MPI_INT, cl_args.nr, 0, MPI_COMM_WORLD, &rstat);
      MPI_Isend(&cl_args.nrp, 1, MPI_INT, cl_args.nr, 0, MPI_COMM_WORLD, &rrstat);
    }

#ifdef DEBUG
  fprintf(stderr, "proc %d waiting for # of particles, t=%g\n", rank, tt);
#endif

  if(rank != 0)
    {
#ifdef DEBUG
      fprintf(stderr, "proc %d: rank != 0, t=%g\n", rank, tt);
#endif

      MPI_Wait(&lstat, &stat);
#ifdef DEBUG
      fprintf(stderr, "left: Processor %d receiving %d particles from proc %d\n",
	      rank, nlr, cl_args.nl);
#endif
    }
  if(rank != last)
    {
#ifdef DEBUG
      fprintf(stderr, "proc %d: rank != last, t=%g\n", rank, tt);
#endif

      MPI_Wait(&rstat, &stat);
#ifdef DEBUG
      fprintf(stderr, "right: Processor %d receving %d particles from proc %d\n",
	      rank, nrr, cl_args.nr);
#endif
    }

#ifdef DEBUG
  fprintf(stderr, "Rank %d ready for particle communication, t=%g:\n", rank, tt);
#endif


  // here, one would reallocate, if necessary
  // for this simple test, we use unrealistically
  // large buffers to avoid this possibility

  if(rank != 0)
    {
      // receive from left
      MPI_Irecv(&(pdata[2*npart]), nlr, MPI_SCALAR, cl_args.nl, 0, MPI_COMM_WORLD, &plstat);

      npart += (nlr/2);

      // send to left
      MPI_Isend(lhsbuf, cl_args.nlp, MPI_SCALAR, cl_args.nl, 0, MPI_COMM_WORLD, &pllstat);

    }
  if(rank != last)
    {
      MPI_Irecv(&(pdata[2*npart]), nrr, MPI_SCALAR, cl_args.nr, 0, MPI_COMM_WORLD, &prstat);
      npart += (nrr/2);

      // send to right
      MPI_Isend(rhsbuf, cl_args.nrp, MPI_SCALAR, cl_args.nr, 0, MPI_COMM_WORLD, &prrstat);

      MPI_Wait(&prstat, &stat);
    }
  //  if(rank != 0) { MPI_Wait(&plstat, &stat); }

  #ifdef DEBUG
  fprintf(stderr, "rank %d waiting for MPI_Barrier(), t=%g\n", rank, tt);
  #endif

  // wait for all processors to catch up
  MPI_Barrier(MPI_COMM_WORLD);

  #ifdef DEBUG
  fprintf(stderr, "rank %d done with communicateparticles(), t=%g.\n", rank, tt);
  #endif

}

/*###########################################
 * #######  GPU particle mover functions ##############
 *########################################### */

__device__ double gridxCUDA(double x, double xl, double dx)
{
  return ((x - xl)/dx);
}

__device__ double particleweightCUDA(double x, double *array, double xl, double dx)
{
  int n;
  double frac, xx;

  xx = gridxCUDA(x, xl, dx);
  n = int(xx);
  frac = xx-n;

  return frac*array[n+1] + (1.0 - frac)*array[n];
}

__device__ void weightparticleCUDA(double x, double weight, double *array, double xl, double dx)
{
  int n;
  //double frac, xx;
  double frac, xx;

  xx = gridxCUDA(x, xl, dx);
  n = int(xx);
  frac = xx-n;

  atomicAdd_double(&(array[n]),weight*(1.0 - frac));
  atomicAdd_double(&(array[n+1]),weight*frac);
}

__global__ void moveparticlesCUDA(double *pdata, bool *alive, double *Earray, 	
                                      double *narray, double xl, double xr, double dx, double dt, double  qm, double qscale,	
                                      int npart)
{
  double E;



  int j = blockIdx.x*blockDim.x+threadIdx.x;

if(j < npart && alive[j]==true)	
  {	
    int index=2*j;	
    E = particleweightCUDA(pdata[index], Earray, xl, dx);	
      pdata[index+1] = pdata[index+1] + qm*E;	
      pdata[index] = pdata[index] + pdata[index+1]*dt;	
      if((pdata[index] > xl) && (pdata[index] < xr))	
      {	
        weightparticleCUDA(pdata[index], qscale, narray, xl, dx);	
	//printf("%g %g , ", pdata[index], pdata[index+1]);
      }	
      else	
      {	
        alive[j]=false;	
      }	
  }	
}

/* GPU (CUDA) version of  moveparticles */
void moveparticlesGPU(void *buffers[], void *_cl_args)
{

  struct cl_args_t *cl_args = (struct cl_args_t*)_cl_args;
  /* definition of host (CPU) arrays */
  double *pdata = (double *)STARPU_VECTOR_GET_PTR(buffers[0]);
  bool *alive = (bool*)STARPU_VECTOR_GET_PTR(buffers[1]);	
  double *lhsbuf = (double *)STARPU_VECTOR_GET_PTR(buffers[2]);
  double *rhsbuf = (double *)STARPU_VECTOR_GET_PTR(buffers[3]); 
  double *narray = (double *)STARPU_VECTOR_GET_PTR(buffers[4]);	
  double *Earray = (double *)STARPU_VECTOR_GET_PTR(buffers[5]);	


  starttime(PARTICLES);
  cudaMemset(narray, 0, ng*sizeof(Scalar));

  /* Executing kernel */
  int threadsPerBlock = 256;
  //int blocksPerGrid =(2*npart + threadsPerBlock - 1) / threadsPerBlock;
  int blocksPerGrid =(npart + threadsPerBlock - 1) / threadsPerBlock;	
  //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
  moveparticlesCUDA<<<blocksPerGrid, threadsPerBlock, 0, starpu_cuda_get_local_stream()>>>(pdata, alive, Earray, narray,	
                                                            cl_args->xl, cl_args->xr, cl_args->dx, cl_args->dt, cl_args->qm, qscale, npart);
  cudaStreamSynchronize(starpu_cuda_get_local_stream());
  endtime(PARTICLES);


}


/*###########################################
 * #######  Printing functions ##############
 *########################################### */

void print_diagnostics(char *fmod, int close=0)
{
  static int init = 0;
  char ffname[STR_SIZE];
  char temp[STR_SIZE];

  if(!init)
    {
      #ifdef DEBUG
      //      fprintf(stderr, "\nInitializing diagnostics...\n");
      #endif

      sprintf(diags[NDENSITY], "density");
      thediags[NDENSITY] = narray;
      sprintf(diags[NEFIELD], "E");
      thediags[NEFIELD] = Earray;
      sprintf(diags[NVOLTAGE], "phi");
      thediags[NVOLTAGE] = phiarray;
      sprintf(diags[NVXPHASE], "vxx");
      thediags[NVXPHASE] = pdata;

      for(int i=0; i < NDIAGNOSTICS; i++)
	{
	  sprintf(temp, "%s%s", diags[i], procname);
	  strcpy(diags[i], temp);
	}

      init = 1;
    }

  if(close)
    {
      return;
    }

  for(int i=0; i < NDIAGNOSTICS; i++)
    {
      sprintf(ffname, "%s%s%s", diags[i], fmod, fext);
      if(i < NVXPHASE)
	{
	  printtarray(i, ng, cl_args.dx);
	}
      else
	{
	  printphase(ffname);
	}
    }
  printt();
}


void printdata(void)
{
  printf("rank %d of %d processors\n", rank, nproc);
  fprintf(stdout, "MPI version: %d.%d, ", ver, subver);
  fprintf(stdout, "MPI_Wtime precision: %g\n", wtres);
  printf("t=%.3g, dt=%.3g, ntimesteps=%d, ", t, cl_args.dt, ntimesteps);
  printf("density=%.3g, wp=%.3g, np2c=%.3g\n", density, wp, np2c);
  printf("%d particles, ng=%d, nc=%d, ", npart, ng, nc);
  printf("L=%g, Ll=%g, dx=%g, area=%g\n", L, Ll, cl_args.dx, area);

  fflush(stdout);
}


/*###########################################
 * #######  End functions ##############
 *########################################### */


// free fields arrays
void free_fields(void)
{
  allocate_arrays(1);
}




void free_particles(void)
{
  starpu_data_unregister(pdata_handle);
  starpu_data_unregister(alive_handle);
  starpu_data_unregister(lhsbuf_handle);
  starpu_data_unregister(rhsbuf_handle);

  starpu_free(pdata);
  starpu_free(check_alive);
  starpu_free(lhsbuf);
  starpu_free(rhsbuf);
   
}

void Quit(void)
{
  free_particles();
  free_fields();

  allocate_fieldbuffers(1);

  if(diagnosticsflag)
    {
      printt(1);
      printtarray(0, 0, 0., 1);
    }

  #ifdef DEBUG
  //  fprintf(stderr, "proc %d ready to finalize MPI\n", rank);
  #endif


  //exit(0);
  //MPI_Abort(MPI_COMM_WORLD, 0);

   MPI_Finalize();
   exit(0);
}




/*###########################################
 * #######  Main functions ##############
 *########################################### */


void mainloop(Scalar tmax, Scalar dt)
{
  tt = 0.;
  while(tt < tmax)
    {
      #ifdef DEBUG
      fprintf(stderr, "rank %d mainloop: t=%g, npart=%d{\n",
	      rank, tt, npart);
      #endif
      
      // reset needed variables
      cl_args.nlp = cl_args.nrp = 0;

      nparttot += npart;

      /*Particle mover task*/
      starpu_task_submit(task);

 
      #ifdef DEBUG
      //      fprintf(stderr, "after move particles, npart=%d\n", npart);
      #endif

      if(nproc > 1)
	{
	  starttime(MPITIME);
	  communicateparticles();
	  endtime(MPITIME);
	}

      /*Advance fields mover tasks*/
            starttime(FIELDS);

      //starpu_task_submit(task_solver);

      starpu_task_submit(task_fields);
      starpu_task_wait_for_all();
            endtime(FIELDS);
    //starpu_cusparse_shutdown();

      #ifdef DEBUG
      //      fprintf(stderr, "after fields, npart=%d\n", npart);
      #endif

      if(diagnosticsflag)
	{
	  starpu_data_acquire(narray_handle, STARPU_R);
	  starpu_data_acquire(Earray_handle, STARPU_R);
	  starpu_data_acquire(phiarray_handle, STARPU_R);

	  char msg[] = "";
	  starttime(DIAGNOSTICS);
	  print_diagnostics(msg);
	  endtime(DIAGNOSTICS);
	  starpu_data_release(narray_handle);
	  starpu_data_release(Earray_handle);
	  starpu_data_release(phiarray_handle);

	}

      #ifdef DEBUG
      //      fprintf(stderr, "after diagnostics, npart=%d\n", npart);
      #endif

      tt += dt;

      MPI_Barrier(MPI_COMM_WORLD);

      #ifdef DEBUG
      fprintf(stderr, "}[%d]\n", rank);
      #endif
    }
}



int main(int argc, char **argv)
{
  /*get info on our GPU, defaulting to first one*/
  CUDA_ERROR(cudaGetDeviceProperties(&prop,0));
  printf("Found GPU '%s' with %g Gb of global memory, max %d threads per block, and %d multiprocessors\n", 
      prop.name, prop.totalGlobalMem/(1024.0*1024.0),
      prop.maxThreadsPerBlock,prop.multiProcessorCount);
      

  /*init CUDA*/
  CUDA_ERROR(cudaSetDevice(0));

  starpu_init(NULL);

  MPI_Init(&argc, &argv);

  tstart = MPI_Wtime();

  starttime(INIT);

  parsecmdline(argc, argv);

  init();

  #ifdef DEBUG
  if(rank == 0)
    {
      printdata();
    }
  #endif

  endtime(INIT);

  mainloop(t, cl_args.dt);

  #ifdef DEBUG
  //  fprintf(stderr, "Proc %d done with mainloop\n", rank);
  #endif

  MPI_Barrier(MPI_COMM_WORLD);

  tend = MPI_Wtime();

  if(rank == 0)
    {
      printtimes(stdout);
    }

  #ifdef DEBUG
  //  fprintf(stderr, "Proc %d ready to Quit()\n", rank);
  #endif

  Quit();

#ifdef DEBUG
  //  fprintf(stderr, "Proc %d DONE!!!\n", rank);
#endif
  starpu_shutdown();
  return 0;
}
