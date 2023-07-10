#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* no of grid points */
const unsigned int scale = 2;
const unsigned int nx = 32 * scale;
const unsigned int ny = nx;

/* no of directions */
const unsigned int ndir = 9;

/* memory size for grid points */
const size_t mem_size_ndir = sizeof(double) * nx * ny * ndir;
const size_t mem_size_scalar = sizeof(double) * nx * ny;

/* weights */
const double w0 = 4.0 / 9.0;  // zero weight
const double ws = 1.0 / 9.0;  // adjacent weight
const double wd = 1.0 / 36.0; // diagonal weight
const double wi[] = {w0, ws, ws, ws, ws, wd, wd, wd, wd};

/* velocities */
const int dirx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int diry[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

/* kinematic viscosity */
const double nu = 1.0 / 6.0;

/* taylor green parameters */
const double u_max = 0.04 / scale;
const double rho0 = 1.0;

/* no of time steps */
const unsigned int nsteps = 200 * scale * scale;
const unsigned int nsave = 50 * scale * scale;

/* indexing 2d scalar array to 1d array */
inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return nx * y + x;
}

/* indexing 2d distribution function for each velocity direction to 1d array */
inline size_t field_index(unsigned int x, unsigned int y, unsigned int d)
{
    return nx * (y + ny * d) + x;
}

/* function declarations */
void taylor_green(unsigned int, unsigned int, unsigned int, double *, double *, double *);
void initialize_flow(unsigned int, double *, double *, double *);
void init_equilibrium(double *, double *, double *, double *);
void compute_rho_u(double *, double *, double *, double *);
void stream(double*,double*);
void collide(double *, double *, double *, double *);
void write_vtk(double *, double *, double *, unsigned int);

int main()
{
    /* allocate memory for distribution function, density and velocity field */
    double *f1 = (double *)malloc(mem_size_ndir);
    double *f2 = (double *)malloc(mem_size_ndir);
    double *rho = (double *)malloc(mem_size_scalar);
    double *ux = (double *)malloc(mem_size_scalar);
    double *uy = (double *)malloc(mem_size_scalar);

    /* compute total memory allocated for all the fields */
    double bytesPerMiB = 1024.0 * 1024.0;
    size_t total_mem_bytes = 2 * mem_size_ndir + 3 * mem_size_scalar;

    if (f1 == NULL || f2 == NULL || rho == NULL || ux == NULL || uy == NULL)
    {
        fprintf(stderr, "Error: unable to allocate required memory (%.1f MiB).\n", total_mem_bytes / bytesPerMiB);
        exit(-1);
    }

    /* compute Taylor-Green flow at t=0 to initialise rho, ux, uy fields. */
    initialize_flow(0, rho, ux, uy);

    /* compute equilibrium distribution from rho, ux, uy */
    init_equilibrium(f1, rho, ux, uy);

    for (unsigned int n = 0; n < nsteps; ++n)
    {
        /* stream from f1 and storing to f2 */
        stream(f1, f2);

        /* calculate post-streaming density and velocity */
        compute_rho_u(f2, rho, ux, uy);

        /* perform collision on f2 */
        collide(f2, rho, ux, uy);

        if ((n + 1) % nsave == 0)
        {
            write_vtk(rho, ux, uy, n);
        }

        // swap pointers
        double *temp = f1;
        f1 = f2;
        f2 = temp;
    }

    /* deallocate memory */
    free(f1);
    free(f2);
    free(rho);
    free(ux);
    free(uy);
}

/* taylor green analytical solution of incompressible NS equation */
void taylor_green(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v)
{
    double kx = 2.0 * M_PI / nx;
    double ky = 2.0 * M_PI / ny;
    double td = 1.0 / (nu * (kx * kx + ky * ky));

    double X = x + 0.5;
    double Y = y + 0.5;
    double ux = -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y) * exp(-1.0 * t / td);
    double uy = u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y) * exp(-1.0 * t / td);
    double P = -0.25 * rho0 * u_max * u_max * ((ky / kx) * cos(2.0 * kx * X) + (kx / ky) * cos(2.0 * ky * Y)) * exp(-2.0 * t / td);
    double rho = rho0 + 3.0 * P;

    *r = rho;
    *u = ux;
    *v = uy;
}

/* set velocity and pressure at each grid point */
void initialize_flow(unsigned int t, double *r, double *u, double *v)
{
    for (unsigned int x = 0; x < nx; ++x)
    {
        for (unsigned int y = 0; y < ny; ++y)
        {
            size_t index = scalar_index(x, y);
            taylor_green(t, x, y, &r[index], &u[index], &v[index]);
        }
    }
}

/* get equilibrium distribution function from density and velocities */
void init_equilibrium(double *f, double *r, double *u, double *v)
{
    for (unsigned int y = 0; y < ny; ++y)
    {
        for (unsigned int x = 0; x < nx; ++x)
        {
            double rho = r[scalar_index(x, y)];
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];

            for (unsigned int i = 0; i < ndir; ++i)
            {
                double cidotu = dirx[i] * ux + diry[i] * uy;

                f[field_index(x, y, i)] = wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));
            }
        }
    }
}

/* write density and velocity in vtk format */
void write_vtk(double *r, double *u, double *v, unsigned int n)
{
    char filename[50];
    sprintf(filename, "output_%04d.vtk", n);

    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        fprintf(stderr, "Failed to open file for writing: %s\n", filename);
        exit(-1);
    }

    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "Scalar Data Example\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET STRUCTURED_GRID\n");
    fprintf(f, "DIMENSIONS %d %d 1\n", nx, ny);
    fprintf(f, "POINTS %d double\n", nx * ny);

    /* write grid points */
    for (unsigned int x = 0; x < nx; ++x)
        for (unsigned int y = 0; y < ny; ++y)
            fprintf(f, "%d %d 0\n", x, y);

    fprintf(f, "POINT_DATA %d\n", nx * ny);

    /* write density */
    fprintf(f, "SCALARS density double\n");
    fprintf(f, "LOOKUP_TABLE default\n");

    for (unsigned int x = 0; x < nx; ++x)
    {
        for (unsigned int y = 0; y < ny; ++y)
        {
            size_t index = scalar_index(x, y);
            fprintf(f, "%.15e\n", r[index]);
        }
    }

    /* write velocity */
    fprintf(f, "VECTORS velocity double\n");

    for (unsigned int x = 0; x < nx; ++x)
    {
        for (unsigned int y = 0; y < ny; ++y)
        {
            size_t index = scalar_index(x, y);
            fprintf(f, "%.15e %.15e 0.0\n", u[index], v[index]);
        }
    }

    fclose(f);
}

/* compute flow properties */
void compute_rho_u(double *f, double *r, double *u, double *v)
{
    for (unsigned int x = 0; x < nx; ++x)
    {
        for (unsigned int y = 0; y < ny; ++y)
        {
            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;

            for (unsigned int i = 0; i < ndir; ++i)
            {
                rho += f[field_index(x, y, i)];
                ux += dirx[i] * f[field_index(x, y, i)];
                uy += diry[i] * f[field_index(x, y, i)];
            }

            r[scalar_index(x, y)] = rho;
            u[scalar_index(x, y)] = ux / rho;
            v[scalar_index(x, y)] = uy / rho;
        }
    }
}

/* collision operation */
void collide(double *f, double *r, double *u, double *v)
{
    // useful constants:
    const double tauinv = 2.0 / (6.0 * nu + 1.0); // 1/tau
    const double omtauinv = 1.0 - tauinv;         // 1 - 1/tau

    for (unsigned int x = 0; x < nx; ++x)
    {
        for (unsigned int y = 0; y < ny; ++y)
        {
            double rho = r[scalar_index(x, y)];
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];

            for (unsigned int i = 0; i < ndir; ++i)
            {
                // calculate dot product
                double cidotu = dirx[i] * ux + diry[i] * uy;

                // calculate equilibrium
                double feq = wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));

                // relax to equilibrium
                f[field_index(x, y, i)] = omtauinv * f[field_index(x, y, i)] + tauinv * feq;
            }
        }
    }
}

/* streaming operation */
void stream(double *f_src, double *f_dst)
{
    for (unsigned int x = 0; x < nx; ++x)
    {
        for (unsigned int y = 0; y < ny; ++y)
        {
            for (unsigned int i = 0; i < ndir; ++i)
            {
                // enforce periodicity
                unsigned int xmd = (nx + x - dirx[i]) % nx;
                unsigned int ymd = (ny + y - diry[i]) % ny;

                f_dst[field_index(x, y, i)] = f_src[field_index(xmd, ymd, i)];
            }
        }
    }
}
