// global variables
#define Scalar double
#define MPI_SCALAR MPI_DOUBLE

#define STR_SIZE 501


#define DEBUG

#define NG 100
#define NP 10

// definitions for diagnostics
#define NDENSITY 0
#define NEFIELD  1
#define NVOLTAGE 2
#define NVXPHASE 3
#define NDIAGNOSTICS 4


//Particle mover arrays
double *pdata;
double *pdata_gpu;
bool *check_alive;
int npart;

// field arrays
double *Earray;
double *phiarray;
double *narray;
double *nback;
double *Earray_gpu;
double *narray_gpu;
double *d_a,*d_b,*d_c,*d_r;
int *nc_gpu;

//GPU values
int Nthreads_max;
int Nblocks_particles;
int Nblocks_fields;

cudaDeviceProp prop;

//Physical and simulation values
Scalar area;
Scalar density;
Scalar np2c;
Scalar q, m, qm;
Scalar qscale;
Scalar epsilon;
Scalar wp;
Scalar El, Er;
Scalar dx, dt, t;
double xl, xr, L, Ll;
int ntimesteps;
Scalar lhsvoltage;
int ng, nc, ngglobal;


int ver, subver;

Scalar dtfactor;
int dg_steps;
bool diag_flag;



// variables for timing
#define PARTICLES 0
#define FIELDS 1
#define DIAGNOSTICS 2
#define INIT 3
#define MEMORY_GPU 4
#define PERFORMANCE 5
#define NTIMES 6
Scalar tttimes[NTIMES];
char names[NTIMES][STR_SIZE];
Scalar tstart, tend;  // global times
Scalar wtres;
Scalar nparttot;

Scalar tt;  // universal time



char diags[NDIAGNOSTICS][STR_SIZE];
Scalar * thediags[NDIAGNOSTICS];

char fext[] = ".dat";

char procname[STR_SIZE];

// global functions
void Quit(void);


//atomicAdd for doubles
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd_double(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif



