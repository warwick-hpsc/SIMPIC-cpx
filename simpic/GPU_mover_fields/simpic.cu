#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>


#include "global.h"
#include <mpi.h>
#include "timing.cpp"
#include "diagnostics.cpp"
#include "fields.cu"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>


/*################################################################/*
 * Initializers *
 * ###############################################################*/

/*Allocates memory for particles*/
void allocate_particles(void)
{

  pdata = new double[2*npart*sizeof(double)];

  /*GPU data allocation*/
  CUDA_ERROR(cudaMalloc((void**)&pdata_gpu,2*npart*sizeof(double)));
  CUDA_ERROR(cudaMalloc((void**)&check_alive,npart*sizeof(bool))); 
 
 
  /*Number of blocks used for GPU particle mover*/
  Nblocks_particles = int(npart/Nthreads_max) + 1;
 
  /*Initial distribution*/
  for(int i=0, j=0; i < npart; i++, j+=2)
    {
      pdata[j] = xl + ((i + 1)*(xr-xl))/Scalar(npart + 1);
      pdata[j+1] = 0.;
    }

  /*Copy and set data to GPU*/
  CUDA_ERROR(cudaMemcpy(pdata_gpu, pdata, 2*npart*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_ERROR(cudaMemset(check_alive,1, npart*sizeof(bool)));
}

/*Gets input data from command line*/
void parsecmdline(int argc, char **argv)
{
  int i;
  int nnppcc;   // number of particles per cell
  int nnccpppp; // number of cells per processor 
  int nnnttt;   // number of time steps

  double ddttff, lhsdefault;

  ddttff = lhsdefault = 0.;
  nnppcc = nnccpppp = nnnttt = 0;
  dg_steps= 1;
  
  for(i=1; i < argc; i++)
    {
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
      
      else if(strcmp(argv[i], "-diagsteps") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &dg_steps);
	}
      else // default case
	{
	  fprintf(stderr, "\nError:\n");
	  fprintf(stderr, "Unrecognized argument \"%s\"\n", argv[i]);
	  Quit();
	}
    }


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

  /*Set simulation variables from input data*/
  ntimesteps = nnnttt;
  nc = nnccpppp;
  ng = nc + 1;
  npart = nc*nnppcc;
  dtfactor = ddttff;
  lhsvoltage = lhsdefault;
}

/*Init function*/
void init(void)
{

  /*Remove old .dat files*/
  system("rm -f *.dat");

  /*Get MPI data*/
  MPI_Get_version(&ver, &subver);  
  wtres = MPI_Wtick();

  /*Reset timer variables*/
  for(int j=0; j < NTIMES; j++)
    {
      tttimes[j] = 0.;
    }

  sprintf(names[PARTICLES], "Particles");
  sprintf(names[FIELDS], "Fields");
  sprintf(names[DIAGNOSTICS], "Diagnostics");
  sprintf(names[INIT], "Initialization"); 
  sprintf(names[MEMORY_GPU], "GPU memory");
  sprintf(names[PERFORMANCE], "Tridiagonal calculation");

  nparttot = 0.;

  /*Physical constants for simulation*/
  density = 1.E13;
  epsilon = 8.85E-12;
  area = 1.;
  L = 1.;
  q = 1.602E-19;
  m = 9.1E-31;

  /*space variables*/
  xl = 0; 	//Left boundary position
  xr = L;	//Right boundary position
  Ll = xr - xl;	//Length of the spatial domain
  dx = Ll / Scalar(nc);	//Distance between space grid points

  np2c = density*area*Ll/Scalar(npart);

  /*Initialise particles and fields*/
  allocate_particles();
  init_fields();

  /*Calculate time step from plasma frequency*/
  wp = density*q*q/(epsilon*m);
  wp = sqrt(wp);

  /*Simulation constants*/
  dt = dtfactor/wp;
  qscale = np2c/(area*dx);
  t = ntimesteps*dt;
  qm = q*dt/m;

}

/*Free CPU and GPU data*/
void free_particles(void)
{
  delete [] pdata;

  CUDA_ERROR(cudaFree(pdata_gpu));
  CUDA_ERROR(cudaFree(check_alive));
}


/*Free fields arrays*/
void free_fields(void)
{
  allocate_arrays(1);
}

/*Free memory at the end of the program*/
void Quit(void)
{
  free_particles();
  free_fields();

   MPI_Finalize();
   exit(0);
}

/*################################################################/*
 * Printers *
 * ###############################################################*/

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
	  printtarray(i, ng, dx);
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
  fprintf(stdout, "MPI version: %d.%d, ", ver, subver);
  fprintf(stdout, "MPI_Wtime precision: %g\n", wtres);
  printf("t=%.3g, dt=%.3g, ntimesteps=%d, ", t, dt, ntimesteps);
  printf("density=%.3g, wp=%.3g, np2c=%.3g\n", density, wp, np2c);
  printf("%d particles, ng=%d, nc=%d, ", npart, ng, nc);
  printf("L=%g,  dx=%g, area=%g\n", L,  dx, area);

  fflush(stdout);
}

/*################################################################/*
 * Particle mover functions *
 * ###############################################################*/
 

/*Extrapolates values from particles to grid points.
 * Used to calculate density at the grid points*/
__device__ void weightparticleCUDA(double x, double weight, double *array, double xl, double dx)
{
  int n;
  float frac, xx;

  xx = (x-xl)/dx;
  n = int(xx);
  frac = xx-n;
  
  atomicAdd_double(&(array[n]),weight*(1.0 - frac));
  atomicAdd_double(&(array[n+1]),weight*frac);
}

/*Interpolates values for particle between two grid points.
 * Used to compute E field in particle position*/
__host__ __device__ double particleweight(double x, double *array, double xl, double dx, int nc)
{
  int n;
  double frac, xx;

  xx = (x-xl)/dx;
  n = int(xx);
  frac = xx-n;

  return frac*array[n+1] + (1.0 - frac)*array[n];
}

/*Computes new position and velocity for the particle. 
 * Calculates contribution to charge density for this particle */
__global__ void moveparticlesCUDA(double *pdata, double *Earray, bool *alive, double *narray, double xl, double xr, double dx, double dt, double qm, double qscale, int npart, int nc)
{
  double E;
  int j = blockIdx.x*blockDim.x+threadIdx.x;
  

  if(j<npart && alive[j]==true)
 {
   int index = 2*j;
   E = particleweight(pdata[index],Earray,xl,dx, nc);
   pdata[index+1] += qm*E;
   pdata[index] += pdata[index+1]*dt;

   if(pdata[index]>xl && pdata[index]<xr)
   {
     weightparticleCUDA(pdata[index],qscale,narray,xl,dx);
   }
   else
   {
     alive[j]=false;
   } 
 }
} 

/*Main particle mover function.
 * Moves all particles at dt */
void moveparticlesGPU(double dt)
{
  starttime(PARTICLES);

  starttime(MEMORY_GPU);
  CUDA_ERROR(cudaMemset(narray_gpu, 0, ng*sizeof(double)));
  endtime(MEMORY_GPU);
  
  moveparticlesCUDA<<<Nblocks_particles, Nthreads_max>>>(pdata_gpu,Earray_gpu, check_alive, narray_gpu,xl,xr,dx,dt,qm,qscale,npart, nc);
  
  if(diag_flag){
    starttime(MEMORY_GPU);
    CUDA_ERROR(cudaMemcpy(narray,narray_gpu,ng*sizeof(double),cudaMemcpyDeviceToHost));
    CUDA_ERROR(cudaMemcpy(pdata, pdata_gpu, 2*npart*sizeof(double),cudaMemcpyDeviceToHost));
    endtime(MEMORY_GPU);
    /*fix density end points because of half volume*/
    narray[0] *= 2.;
    narray[nc] *= 2.;
  }

  endtime(PARTICLES);
}


/*################################################################/*
 * Main functions *
 * ###############################################################*/

void mainloop(Scalar tmax, Scalar dt)
{
  tt = 0.;
  int nsteps = 0;
  while(tt < tmax)
    {
      #ifdef DEBUG
      fprintf(stderr, "Completed %g %: t=%g, npart=%d \n", 
	       float(nsteps)/(float(ntimesteps))*100.,tt, npart);
      #endif
      
      nsteps+=1;
      nparttot += npart;
      diag_flag = ((nsteps % dg_steps == 0)|| (nsteps == ntimesteps-1));
      /*Recalculate particle positions, velocity and density*/
      moveparticlesGPU(dt);    

      /*With new density recalculate E field*/
      advancefields(dt);

      /*Save diagnostics each dg_steps*/
      if( diag_flag )
      {
	  starttime(DIAGNOSTICS);
	  char msg[] = "";
	  print_diagnostics(msg);
	  endtime(DIAGNOSTICS);
      }

      tt += dt;
    }
}

int main(int argc, char **argv)
{
  /*Get CUDA info*/
  CUDA_ERROR(cudaGetDeviceProperties(&prop,0));
  printf("Found GPU '%s' with %g Gb of global memory, max %d threads per block, and %d multiprocessors\n", prop.name, prop.totalGlobalMem/(1024.0*1024.0), prop.maxThreadsPerBlock, prop.multiProcessorCount);
  Nthreads_max = prop.maxThreadsPerBlock;  //Max threads to be used
  CUDA_ERROR(cudaSetDevice(0));		//Set the CUDA device

  //Init MPI
  MPI_Init(&argc, &argv);
  tstart = MPI_Wtime();

  starttime(INIT);
  parsecmdline(argc, argv); 	//Get input data

  init();	//Initialise

  #ifdef DEBUG
    printdata();
    std::cout<<"Max GPU threads per block "<<Nthreads_max;
    std::cout<<", using "<< Nblocks_particles<<" blocks for particles and ";
    std::cout<<Nblocks_fields<< " blocks for fields \n";
    std::cout<<"At least 1 thread is per particle is needed and 1 thread per grid point \n";
  #endif
  endtime(INIT);
  
  mainloop(t, dt);

  tend = MPI_Wtime();

  printtimes(stdout);

  Quit();

  return 0;
}
