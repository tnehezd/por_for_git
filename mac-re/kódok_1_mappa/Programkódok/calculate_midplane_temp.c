#include <stdio.h>
#include <math.h>

// Physical and numerical constants
#define RHO 1.501780435e-3           // Density [g/cm^3]
#define SIGMA0 0.0001                // Surface density constant [units?]
#define NGRID 100                    // Number of grid points for radial discretization
#define RMIN 0.1                    // Minimum radius of the disk [AU]
#define RMAX 15.0                   // Maximum radius of the disk [AU]
#define DGRID ((RMAX - RMIN) / (NGRID - 1)) // Radial step size (AU)

#define ALPHA_VISC 0.01              // Alpha parameter for disk viscosity (dimensionless)
#define ASP_RATIO 0.05               // Aspect ratio (disk thickness / radius)
#define G 1.0                       // Gravitational constant (normalized units)
#define SB_CONST 2.533732811E14      // Stefan-Boltzmann constant (W AU^-2 K^-4)
#define T_BG 10.0                   // Background temperature (K)
#define M_STAR 1.0                  // Mass of central star (normalized units)

// FILE pointer for output file
FILE *fil;

/**
 * Calculates the radial-dependent kinematic viscosity ν(r) with a "dead zone" profile.
 *
 * The dead zone models regions of reduced turbulence within the disk where the effective
 * viscosity is lowered. This is implemented by modulating alpha viscosity parameter
 * with smooth transitions (using hyperbolic tangent functions).
 *
 * @param r Radial distance from the star (in AU)
 * @return Effective kinematic viscosity at radius r
 */
double viscosity(double r) {
    // Dead zone parameters defining inner and outer boundaries (in AU)
    const double r_dze_i = 0.0;
    const double r_dze_o = 10.0;
    const double Dr_dze_i = 0.0;
    const double Dr_dze_o = 0.5;

    // Reduction factor inside the dead zone
    const double a_mod = 0.01;

    // Smooth transition of alpha parameter across dead zone edges
    double alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * 
                     (tanh((r - r_dze_i) / (Dr_dze_i + 1e-10)) + tanh((r_dze_o - r) / Dr_dze_o));

    // Compute viscosity based on alpha, aspect ratio squared, gravitational parameter, and radius
    return ALPHA_VISC * alpha_r * ASP_RATIO * ASP_RATIO * G * sqrt(r);
}

/**
 * Returns a constant coefficient related to thermal radiation calculations.
 * Possibly used for cooling or radiative transfer terms in the disk.
 */
double Coeff_m() {
    return 2 * SB_CONST;
}

/**
 * Keplerian angular frequency Ω_K(r) assuming central mass M_STAR and gravitational constant G.
 *
 * @param r Radius [AU]
 * @return Keplerian angular velocity at radius r
 */
double kep_freq(double r) {
    return sqrt((G * M_STAR) / pow(r, 3));
}

/**
 * Coefficient function corresponding to thermal diffusion term's second derivative component.
 *
 * @param r Radius [AU]
 * @return Coefficient for the radial second derivative in thermal conduction equation
 */
double Coeff_1(double r) {
    return 3.0 * viscosity(r);
}

/**
 * Coefficient function corresponding to the first derivative (advection-like) term in the
 * thermal diffusion equation.
 *
 * @param r Radius [AU]
 * @return Coefficient for the radial first derivative in thermal conduction equation
 */
double Coeff_2(double r) {
    return 9.0 * viscosity(r) / (2.0 * r);
}

int main() {
    double temp, temp_old;
    double a, b, k0;
    double r;
    double t_midplane[NGRID + 2]; // Temperature at the disk midplane; includes ghost zones for boundaries
    double t = 0.0, dt = 0.1, tmax = 10.0;
    int i;

    // Open output file to write final temperature profile
    fil = fopen("temperature.dat", "w");
    if (!fil) {
        fprintf(stderr, "Error: Unable to open output file.\n");
        return 1;
    }

    // Initialize boundary conditions (ghost cells)
    // Fixed temperature at boundaries, assuming ambient/background disk temperature
    t_midplane[0] = T_BG;       
    t_midplane[NGRID + 1] = T_BG;

    // Initial guess and iterative solver for temperature at each radius based on local opacity laws
    for (i = 1; i <= NGRID; i++) {
        r = RMIN + (i - 1) * DGRID;

        temp = 1000.0;  // Starting temperature guess [K]
        do {
            temp_old = temp;

            // Opacity model: piecewise function defining opacity coefficient (k0), density exponent (a), and temperature exponent (b)
            // These values are likely fitted to dust/gas opacity tables for protoplanetary disks or stellar disks.
            if (temp <= 170.0) {
                k0 = 2.250419625E3;
                a = 0.0;
                b = 2.0;
            } else if (temp > 170.0 && temp <= 210.0) {
                k0 = 2.250419637E23;
                a = 0.0;
                b = -7.0;
            } else if (temp > 210.0 && temp <= pow(RHO, 1.0 / 15.0)) {
                k0 = 5.626049092E29;
                a = 0.0;
                b = 1.0;
            } else if (temp > pow(RHO, 1.0 / 15.0) && temp <= 3000.0) {
                k0 = 2.250419637E41;
                a = 2.0 / 3.0;
                b = -9.0;
            } else if (temp > 3000.0 && temp <= pow(RHO, 1.0 / 21.0)) {
                k0 = 2.250419637e-1;
                a = 2.0 / 3.0;
                b = 3.0;
            } else if (temp > pow(RHO, 1.0 / 21.0) && temp <= pow(RHO, 4.0 / 75.0)) {
                k0 = 1.125209818e-29;
                a = 1.0 / 3.0;
                b = 10.0;
            } else if (temp > pow(RHO, 4.0 / 75.0)) {
                k0 = 1.687814744E28;
                a = 1.0;
                b = -5.0 / 2.0;
            }

            // Calculate opacity using the piecewise power law
            double kappa = k0 * pow(RHO, a) * pow(temp_old, b);

            // Optical depth estimation (surface density SIGMA0 assumed)
            double optical_depth_disk = kappa * SIGMA0 / 2.0;

            // Effective optical depth accounting for radiative transfer approximation
            double optical_depth_effective = (3.0 * optical_depth_disk / 8.0) + (sqrt(3.0) / 4.0) + (1.0 / (4.0 * optical_depth_disk));

            // Temperature update from radiative equilibrium balance (simplified)
            // Includes viscous heating (through viscosity and Keplerian frequency) and background temperature floor
            temp = pow(
                (optical_depth_effective * (9.0 / 4.0 * SIGMA0 * viscosity(r) * pow(kep_freq(r), 2)) + Coeff_m() * pow(T_BG, 4))
                / Coeff_m(),
                0.25
            );

        } while (fabs(temp_old - temp) > 1.0e-6);  // iterate until convergence below tolerance

        t_midplane[i] = temp;

        printf("Radius = %.7lf AU, Midplane temperature = %.7lf K\n", r, t_midplane[i]);
    }

    // Time evolution of temperature via explicit finite difference solution of heat conduction
    // This loop updates temperature until final simulation time tmax
    while (t < tmax) {
        double t_midplane_new[NGRID + 2];

        // Copy current temperature profile before update
        for (i = 0; i <= NGRID + 1; i++) {
            t_midplane_new[i] = t_midplane[i];
        }

        // Update temperature at each radial grid point using discretized diffusion equation
        for (i = 1; i < NGRID; i++) {
            r = RMIN + (i - 1) * DGRID;

            // Radial second derivative term multiplied by coefficient Coeff_1 (diffusion term)
            double second_derivative = (t_midplane[i + 1] - 2.0 * t_midplane[i] + t_midplane[i - 1]) / (DGRID * DGRID);

            // Radial first derivative term multiplied by coefficient Coeff_2 (advection-like term)
            double first_derivative = (t_midplane[i + 1] - t_midplane[i - 1]) / (2.0 * DGRID);

            // Explicit Euler timestep for temperature evolution
            t_midplane_new[i] = t_midplane[i] + dt * (Coeff_1(r) * second_derivative + Coeff_2(r) * first_derivative);
        }

        // Update temperature profile for next timestep
        for (i = 1; i < NGRID; i++) {
            t_midplane[i] = t_midplane_new[i];
        }

        t += dt;
    }

    // Write final temperature profile to file
    for (i = 1; i <= NGRID; i++) {
        fprintf(fil, "%.7lf %.7lf\n", RMIN + (i - 1) * DGRID, t_midplane[i]);
    }

    fclose(fil);

    return 0;
}
