// simulate Poiseullie flow in 3D

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

using namespace std;

// set simulation parameters
double Re = 10;
int N = 32;
double ulb = 0.02; // Velocity in lattice units

double max_t = 0.01; // Total time in dim.less units
int vtk_freq = 1;     // Data writing frequency



enum class CellType : uint8_t
{
    bounce_back,
    inlet,
    outlet,
    bulk
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
    double dx = 1. / static_cast<double>(lref);
    double dt = dx * ulb;
    return make_tuple(nu, omega, dx, dt);
}

// Instances of this class are function objects: the function-call operator executes a
// collision-streaming cycle.
struct LBM
{
    static size_t sizeOfLattice(size_t nelem) { return 2 * 19 * nelem; }

    double *lattice;
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

    // Compute unkown distribution function for inlet and outlet
    // Reference: "Implementation of on-site velocity boundary conditions for D3Q19 lattice Boltzmann simulations"
    double onsite(int i, int k, double rho, std::array<double, 3> const &u, const std::array<double, 3> &n)
    {
        double rhoby6 = rho / 6.0;
        double rhoby3 = rho / 3.0;

        // lambda to compute dot product of 3d vectors
        auto dot = [](const auto &x, const auto &y)
        { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; };

        // compute dot products
        double ck_u = dot(c[k], u);
        double ck_n = dot(c[k], n);

        // compute tangent vector
        std::array<double, 3> tk{c[k][0] - ck_n * n[0], c[k][1] - ck_n * n[1], c[k][2] - ck_n * n[2]};

        double tk_u = dot(tk, u);

        double sum = 0.0;

        for (int j = 0; j < 19; ++j)
        {
            double tk_cj = dot(tk, c[j]);
            double cj_n = dot(c[j], n);

            sum += fin(i, j) * tk_cj * (1.0 - fabs(cj_n));
        }

        return -rhoby6 * ck_u - rhoby3 * tk_u + 0.5 * sum;
    }

    // Execute the streaming step on a cell.
    void stream(int i, int k, int iX, int iY, int iZ, double pop_out, double rho, std::array<double, 3> const &u)
    {
        int XX = iX + c[k][0];
        int YY = iY + c[k][1];
        int ZZ = iZ + c[k][2];
        size_t nb = xyz_to_i(XX, YY, ZZ);
        size_t nb_opp = xyz_to_i(iX - c[k][0], iY - c[k][1], iZ - c[k][2]);

        switch (flag[nb])
        {
        case CellType::bounce_back:
            fout(i, opp[k]) = pop_out;
            break;

        case CellType::outlet:
            fout(i, opp[k]) = pop_out + onsite(i, k, 1.0, {0, 0, u[2]}, {0, 0, 1});
            break;

        case CellType::inlet:
            fout(i, opp[k]) = pop_out + onsite(i, k, rho, {0, 0, 1.0}, {0, 0, -1});
            break;

        case CellType::bulk:
            fout(nb, k) = pop_out;
            break;
        }
    };

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

                // TODO Need to remove rho after testing
                stream(i, k, iX, iY, iZ, pop_out, rho, u);
                stream(i, opp[k], iX, iY, iZ, pop_out_opp, rho, u);
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

// Initializes the flag matrix to set up the Poiseuille flow
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
            lbm.flag[i] = CellType::bounce_back;
        else
            lbm.flag[i] = CellType::bulk;
    }
}

// Save the macroscopic variables in a VTK file (for post-processing for
// example with Paraview).
template <class LBM>
void saveVtkFields(LBM &lbm, int time_iter, double dx = 0.)
{
    using namespace std;
    Dim const &dim = lbm.dim;
    if (dx == 0.)
    {
        dx = 1. / dim.nx;
    }
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

int main()
{
    // domain
    Dim dim{N, N, N};

    // compute LB parameters
    auto [nu, omega, dx, dt] = lbParameters(ulb, N - 2, Re);

    // allocate memory for distribution function
    vector<double> lattice_vect(LBM::sizeOfLattice(dim.nelem));
    double *lattice = &lattice_vect[0];

    // flag the cells
    vector<CellType> flag_vect(dim.nelem);
    CellType *flag = &flag_vect[0];

    // parity to decide in which of the two population arrays the pre-collision LB variables are located.
    vector<int> parity_vect{0};
    int *parity = &parity_vect[0];

    // discrete velocities, opposite indices, weights
    auto [c, opp, t] = d3q19_constants();

    // Instantiate the function object for the for_each call.
    LBM lbm{lattice, flag, parity, &c[0], &opp[0], &t[0], omega, dim};

    // Initialize the populations.
    for_each(lattice, lattice + dim.nelem, [&lbm](double &f0)
             { lbm.iniLattice(f0); });

    // Set up the geometry of the cavity.
    setupFlags(lbm);

    // Maximum number of time iterations
    int max_time_iter = static_cast<int>(max_t / dt);

    // time iteration
    for (int time_iter = 0; time_iter < max_time_iter; ++time_iter)
    {
        std::cout << "iter = " << time_iter << "\n";

        // write flow data in vtk
        if (time_iter % vtk_freq == 0)
            saveVtkFields(lbm, time_iter);

        // With the double-population scheme, collision and streaming are fused and executed in the following loop.
        for_each(execution::par_unseq, lattice, lattice + dim.nelem, lbm);

        // After a collision-streaming cycle, swap the parity for the next iteration.
        *parity = 1 - *parity;
    }
}
