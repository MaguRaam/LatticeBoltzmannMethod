#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <execution>
#include <chrono>
#include <cmath>

// sqrt(velocity_X*velocity_X + velocity_Y*velocity_Y + velocity_Z*velocity_Z)

using namespace std;
using namespace std::chrono;

double Re = 100;    // Reynolds number
double ulb = 0.1;   // inlet Velocity in lattice units
std::vector<double> uinlet{0, 0, ulb}; //inlet velocity vector in lattice units
double rholb = 1.0; // Density in lattice units
int Nx = 100;        // Number of nodes in x-direction or length in lattice units
int Ny = 100;        // Number of nodes in y-direction or length in lattice units
int Nz = 500;       // Number of nodes in z-direction or length in lattice units
double max_t = 5000; // Total time in lattice units
int vtk_freq = 4500;  // Frequency in LU for output of terminal message and profiles (use 0 for no messages)

enum class CellType : uint8_t
{
    bounce_back,
    inlet,
    bulk,
    outlet
};

inline auto d3q19_constants()
{
    // The discrete velocities of the d3q19 mesh.
    vector<array<int, 3>> c_vect = {
        {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {-1, -1, 0}, {-1, 1, 0}, {-1, 0, -1}, {-1, 0, 1}, {0, -1, -1}, {0, -1, 1}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {0, 1, 1}, {0, 1, -1}};

    // The opposite of a given direction.
    vector<int> opp_vect =
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8};

    // The lattice weights.
    vector<double> t_vect =
        {
            1. / 18., 1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
            1. / 3.,
            1. / 18., 1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};
    return make_tuple(c_vect, opp_vect, t_vect);
}

// Stores the number of elements in a rectangular-shaped simulation.
struct Dim
{
    Dim(int nx_, int ny_, int nz_)
        : nx(nx_), ny(ny_), nz(nz_),
          nelem(static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz)),
          npop(19 * nelem)
    {
    }
    int nx, ny, nz;
    size_t nelem, npop;
};

// Compute lattice-unit variables and discrete space and time step for a given lattice
// velocity (acoustic scaling).
auto lbParameters(double ulb, int lref, double Re)
{
    double nu = ulb * static_cast<double>(lref) / Re;
    double omega = 1. / (3. * nu + 0.5);
    return make_tuple(nu, omega);
}

// Print the simulation parameters to the terminal.
void printParameters(double Re, double nu, double omega, double ulb, double max_t)
{
    cout << "Lid-driven 3D cavity, " << endl;
    cout << "Re = " << Re << endl;
    cout << "Kinematic viscosity = " << nu << endl;
    cout << "relaxation time = " << 1.0/omega << endl;
    cout << "ulb = " << ulb << endl;
    cout << "max_t = " << max_t << endl;
    cout << "Ma = " << fabs(ulb/sqrt(1.0/3.0)) << std::endl;
    
}

// Instances of this class are function objects: the function-call operator executes a
// collision-streaming cycle.
struct LBM
{
    using CellData = double;
    static size_t sizeOfLattice(size_t nelem) { return 2 * 19 * nelem; }

    CellData *lattice;
    CellType *flag;
    int *parity;
    std::array<int, 3> *c;
    int *opp;
    double *t;
    double omega;
    Dim dim;

    // Convert linear index to Cartesian indices.
    auto i_to_xyz(int i)
    {
        int iX = i / (dim.ny * dim.nz);
        int remainder = i % (dim.ny * dim.nz);
        int iY = remainder / dim.nz;
        int iZ = remainder % dim.nz;
        return std::make_tuple(iX, iY, iZ);
    };

    // Convert Cartesian indices to linear index.
    size_t xyz_to_i(int x, int y, int z)
    {
        return z + dim.nz * (y + dim.ny * x);
    };

    // The "f" operator always refers to the population with parity 0, in which
    // the momentum-exchange term is stored.
    double &f(int i, int k)
    {
        return lattice[k * dim.nelem + i];
    }

    // Get the pre-collision populations of the current time step (their location depends on parity).
    double &fin(int i, int k)
    {
        return lattice[*parity * dim.npop + k * dim.nelem + i];
    }

    // Get the post-collision, streamed populations at the end of the current step.
    double &fout(int i, int k)
    {
        return lattice[(1 - *parity) * dim.npop + k * dim.nelem + i];
    }

    // Initialize the lattice with zero velocity and density 1.
    auto iniLattice(double &f0)
    {
        auto i = &f0 - lattice;
        for (int k = 0; k < 19; ++k)
        {
            fin(i, k) = t[k];
        }
    };

    // Compute the macroscopic variables density and velocity on a cell.
    auto macro(double &f0)
    {
        auto i = &f0 - lattice;
        double X_M1 = fin(i, 0) + fin(i, 3) + fin(i, 4) + fin(i, 5) + fin(i, 6);
        double X_P1 = fin(i, 10) + fin(i, 13) + fin(i, 14) + fin(i, 15) + fin(i, 16);
        double X_0 = fin(i, 9) + fin(i, 1) + fin(i, 2) + fin(i, 7) + fin(i, 8) + fin(i, 11) + fin(i, 12) + fin(i, 17) + fin(i, 18);

        double Y_M1 = fin(i, 1) + fin(i, 3) + fin(i, 7) + fin(i, 8) + fin(i, 14);
        double Y_P1 = fin(i, 4) + fin(i, 11) + fin(i, 13) + fin(i, 17) + fin(i, 18);

        double Z_M1 = fin(i, 2) + fin(i, 5) + fin(i, 7) + fin(i, 16) + fin(i, 18);
        double Z_P1 = fin(i, 6) + fin(i, 8) + fin(i, 12) + fin(i, 15) + fin(i, 17);

        double rho = X_M1 + X_P1 + X_0;
        std::array<double, 3> u{(X_P1 - X_M1) / rho, (Y_P1 - Y_M1) / rho, (Z_P1 - Z_M1) / rho};
        return std::make_pair(rho, u);
    }

    // Execute the streaming step on a cell.
    void stream(int i, int k, int iX, int iY, int iZ, double pop_out, double rho, std::array<double, 3> const &u, double usqr)
    {
        int XX = iX + c[k][0];
        int YY = iY + c[k][1];
        int ZZ = iZ + c[k][2];
        size_t nb = xyz_to_i(XX, YY, ZZ);

        double ck_u = c[k][0] * u[0] + c[k][1] * u[1] + c[k][2] * u[2];

        switch (flag[nb])
        {
        case CellType::bounce_back:
            fout(i, opp[k]) = pop_out;
            break;

        case CellType::inlet:
            fout(i, opp[k]) = pop_out - 6. * t[k] * rho * (c[k][0] * uinlet[0] + c[k][1] * uinlet[1] + c[k][2] * uinlet[2] );
            break;

        case CellType::outlet:
            fout(i, opp[k]) =  -pop_out + 2.0*t[k]*rholb*(1.0 + 4.5*ck_u*ck_u - usqr);
            break;

        case CellType::bulk:
            fout(nb, k) = pop_out;
            break;

        }
    };

    // Execute second-order BGK collision on a cell.
    auto collideBgk(int i, int k, double rho, std::array<double, 3> const &u, double usqr)
    {
        double ck_u = c[k][0] * u[0] + c[k][1] * u[1] + c[k][2] * u[2];
        double eq = rho * t[k] * (1. + 3. * ck_u + 4.5 * ck_u * ck_u - usqr);
        double eqopp = eq - 6. * rho * t[k] * ck_u;
        double pop_out = (1. - omega) * fin(i, k) + omega * eq;
        double pop_out_opp = (1. - omega) * fin(i, opp[k]) + omega * eqopp;
        return std::make_pair(pop_out, pop_out_opp);
    }

    // Execute an iteration (collision and streaming) on a cell.
    void operator()(double &f0)
    {
        int i = &f0 - lattice;
        if (flag[i] == CellType::bulk)
        {
            auto [rho, u] = macro(f0);
            auto [iX, iY, iZ] = i_to_xyz(i);
                   
            double usqr = 1.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

            for (int k = 0; k < 9; ++k)
            {
                auto [pop_out, pop_out_opp] = collideBgk(i, k, rho, u, usqr);

                stream(i, k, iX, iY, iZ, pop_out, rho, u, usqr);
                stream(i, opp[k], iX, iY, iZ, pop_out_opp, rho, u, usqr);
            }

            // Treat the case of the "rest-population" (c[k] = {0, 0, 0,}).
            for (int k : {9})
            {
                double eq = rho * t[k] * (1. - usqr);
                fout(i, k) = (1. - omega) * fin(i, k) + omega * eq;
            }
        }
    }
};

// Compute the average kinetic energy in the domain.
double computeEnergy(LBM &lbm)
{
    Dim const &dim = lbm.dim;
    double energy = 0.;
    for (int iX = 0; iX < dim.nx; ++iX)
    {
        for (int iY = 0; iY < dim.ny; ++iY)
        {
            for (int iZ = 0; iZ < dim.nz; ++iZ)
            {
                size_t i = lbm.xyz_to_i(iX, iY, iZ);
                if (lbm.flag[i] != CellType::bounce_back)
                {
                    auto [rho, v] = lbm.macro(lbm.lattice[i]);
                    energy += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                }
            }
        }
    }
    energy *= 0.5;
    return energy;
}

// Initializes the flag matrix and the populations to set up a cavity flow. In particular, the
// momentum-exchange term is assigned to the wall cells.
void setupFlags(LBM &lbm)
{
    Dim const &dim = lbm.dim;
    for (size_t i = 0; i < dim.nelem; ++i)
    {
        auto [iX, iY, iZ] = lbm.i_to_xyz(i);
        
        if (iX == 0 || iX == dim.nx - 1 || iY == 0 || iY == dim.ny - 1)
            lbm.flag[i] = CellType::bounce_back;
        
        else if (iZ == 0)
            lbm.flag[i] = CellType::inlet;

        else if (iZ == dim.nz - 1)
            lbm.flag[i] = CellType::outlet;

        else
            lbm.flag[i] = CellType::bulk;
    }
}

// Save the macroscopic variables in a VTK file (for post-processing for
// example with Paraview).
template <class LBM>
void saveVtkFields(LBM &lbm, int time_iter)
{
    using namespace std;
    Dim const &dim = lbm.dim;
    double dx = 1.0;

    stringstream fNameStream;
    fNameStream << "sol_" << setfill('0') << setw(7) << time_iter << ".vtk";
    string fName = fNameStream.str();
    ofstream fStream;
    fStream.open(fName.c_str());
    fStream << "# vtk DataFile Version 2.0" << endl;
    fStream << "iteration " << time_iter << endl;
    fStream << "ASCII" << endl;
    fStream << endl;
    fStream << "DATASET STRUCTURED_POINTS" << endl;
    fStream << "DIMENSIONS " << dim.nx - 2 << " " << dim.ny - 2 << " " << dim.nz - 2 << endl;
    fStream << "ORIGIN 0 0 0" << endl;
    fStream << "SPACING " << dx << " " << dx << " " << dx << endl;
    fStream << endl;
    fStream << "POINT_DATA " << (dim.nx - 2) * (dim.ny - 2) * (dim.nz - 2) << endl;
    // Density
    fStream << "SCALARS Density float 1" << endl;
    fStream << "LOOKUP_TABLE default" << endl;
    for (int iZ = 1; iZ < dim.nz - 1; ++iZ)
    {
        for (int iY = 1; iY < dim.ny - 1; ++iY)
        {
            for (int iX = 1; iX < dim.nx - 1; ++iX)
            {
                size_t i = lbm.xyz_to_i(iX, iY, iZ);
                auto [rho, v] = lbm.macro(lbm.lattice[i]);
                fStream << rho << " ";
            }
            fStream << endl;
        }
        fStream << endl;
    }
    fStream << endl;
    // velocity_X
    fStream << "SCALARS velocity_X float 1" << endl;
    fStream << "LOOKUP_TABLE default" << endl;
    for (int iZ = 1; iZ < dim.nz - 1; ++iZ)
    {
        for (int iY = 1; iY < dim.ny - 1; ++iY)
        {
            for (int iX = 1; iX < dim.nx - 1; ++iX)
            {
                size_t i = lbm.xyz_to_i(iX, iY, iZ);
                auto [rho, v] = lbm.macro(lbm.lattice[i]);
                fStream << v[0] << " ";
            }

            fStream << endl;
        }
        fStream << endl;
    }
    fStream << endl;
    // velocity_Y
    fStream << "SCALARS velocity_Y float 1" << endl;
    fStream << "LOOKUP_TABLE default" << endl;
    for (int iZ = 1; iZ < dim.nz - 1; ++iZ)
    {
        for (int iY = 1; iY < dim.ny - 1; ++iY)
        {
            for (int iX = 1; iX < dim.nx - 1; ++iX)
            {
                size_t i = lbm.xyz_to_i(iX, iY, iZ);
                auto [rho, v] = lbm.macro(lbm.lattice[i]);
                fStream << v[1] << " ";
            }
            fStream << endl;
        }
        fStream << endl;
    }
    fStream << endl;
    // velocity_Z
    fStream << "SCALARS velocity_Z float 1" << endl;
    fStream << "LOOKUP_TABLE default" << endl;
    for (int iZ = 1; iZ < dim.nz - 1; ++iZ)
    {
        for (int iY = 1; iY < dim.ny - 1; ++iY)
        {
            for (int iX = 1; iX < dim.nx - 1; ++iX)
            {
                size_t i = lbm.xyz_to_i(iX, iY, iZ);
                auto [rho, v] = lbm.macro(lbm.lattice[i]);
                fStream << v[2] << " ";
            }
            fStream << endl;
        }
        fStream << endl;
    }
}


template <class LBM>
void saveVtkFlags(LBM &lbm)
{
    using namespace std;
    Dim const &dim = lbm.dim;
    double dx = 1.0;

    stringstream fNameStream;
    fNameStream << "flags" << ".vtk";
    string fName = fNameStream.str();
    ofstream fStream;
    fStream.open(fName.c_str());
    fStream << "# vtk DataFile Version 2.0" << endl;
    fStream << "iteration " << 0 << endl;
    fStream << "ASCII" << endl;
    fStream << endl;
    fStream << "DATASET STRUCTURED_POINTS" << endl;
    fStream << "DIMENSIONS " << dim.nx << " " << dim.ny << " " << dim.nz << endl;
    fStream << "ORIGIN 0 0 0" << endl;
    fStream << "SPACING " << dx << " " << dx << " " << dx << endl;
    fStream << endl;
    fStream << "POINT_DATA " << dim.nx * dim.ny * dim.nz << endl;
    // Density
    fStream << "SCALARS Flags int 1" << endl;
    fStream << "LOOKUP_TABLE default" << endl;
    for (int iZ = 0; iZ < dim.nz; ++iZ)
    {
        for (int iY = 0; iY < dim.ny; ++iY)
        {
            for (int iX = 0; iX < dim.nx; ++iX)
            {
                size_t i = lbm.xyz_to_i(iX, iY, iZ);
                
                fStream << static_cast<int>(lbm.flag[i]) << " ";
            }
            fStream << endl;
        }
        fStream << endl;
    }
    fStream << endl;
}


int main()
{
    // CellData is either a double (structure-of-array) or an array<double, 19> (array-of-structure).
    using CellData = typename LBM::CellData;

    Dim dim{Nx, Ny, Nz};

    auto [nu, omega] = lbParameters(ulb, Nx - 2, Re);
    printParameters(Re, nu, omega, ulb, max_t);

    vector<CellData> lattice_vect(LBM::sizeOfLattice(dim.nelem));
    CellData *lattice = &lattice_vect[0];

    // The "vector" is used as a convenient way to allocate the flag array on the heap.
    vector<CellType> flag_vect(dim.nelem);
    CellType *flag = &flag_vect[0];

    // The "vector" is used as a convenient way to allocate the integer value "parity"
    // on the heap. That's a necessity when the code is used on GPUs, because the
    // parity (which is 0/1, depending on whether the current time iteration is even/odd)
    // is modified on the host, and accessed on the device. The parity flag is used
    // to decide in which of the two population arrays the pre-collision LB variables are located.
    vector<int> parity_vect{0};
    int *parity = &parity_vect[0];

    // The lattice constants (discrete velocities, opposite indices, weghts) are mostly
    // used in the not-unrolled versions of the code. They must be allocated on the
    // heap so they can be shared with the device in case of a GPU execution.
    auto [c, opp, t] = d3q19_constants();

    // Instantiate the function object for the for_each call.
    LBM lbm{lattice, flag, parity, &c[0], &opp[0], &t[0], omega, dim};

    // Initialize the populations.
    for_each(lattice, lattice + dim.nelem, [&lbm](CellData &f0)
             { lbm.iniLattice(f0); });

    // Set up the geometry of the cavity.
    setupFlags(lbm);
    //saveVtkFlags(lbm);

    for (int time_iter = 0; time_iter < max_t; ++time_iter)
    {

        std::cout << "iter = " << time_iter << std::endl;

        // write flow data in vtk
        if (time_iter % vtk_freq == 0)
        {
          std::cout << "write vtk \n";
          saveVtkFields(lbm, time_iter);
	}

        // With the double-population scheme, collision and streaming are fused and executed in the following loop.
        for_each(execution::par_unseq, lattice, lattice + dim.nelem, lbm);

        // After a collision-streaming cycle, swap the parity for the next iteration.
        *parity = 1 - *parity;
    }
    
    
    
}
