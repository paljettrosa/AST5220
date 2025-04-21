#ifndef _POWERSPECTRUM_HEADER
#define _POWERSPECTRUM_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <functional>
#include <utility> 
#include <fstream> 
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class PowerSpectrum {
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
    Perturbations *pert        = nullptr;

    // Parameters defining the primordial power-spectrum
    double A_s        = 2.1e-9;
    double n_s        = 0.965;
    double kpivot_Mpc = 0.05;
    
    // The ells's we will compute quantities for
    Vector ells{ 
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};
   
    //=====================================================================
    // [1] Create bessel function splines needed for the LOS integration
    //=====================================================================

    // Splines of bessel-functions for each value of ell in the array above
    std::vector<Spline> j_ell_spline;
    
    // Generate splines of bessel-functions for each ell
    void generate_bessel_function_splines(
        const double z_max, 
        const int npts_z); 
    
    //=====================================================================
    // [2] Do the line of sight integration and spline the result
    //=====================================================================
    
    // Do LOS integration for all ells and all k's in the given k_array
    void line_of_sight_integration(
        Vector &x_array,
        Vector &k_array);
  
    // Do the line of sight integration for a single quantity
    Vector2D line_of_sight_integration_single(
        Vector &x_array,
        Vector &k_array, 
        std::function<double(double,double)> &source_function);
    
    // Splines of the reusult of the LOS integration
    std::vector<Spline> ThetaT_ell_of_k_spline;
    std::vector<Spline> ThetaE_ell_of_k_spline;
    std::vector<Spline> Nu_ell_of_k_spline;
    std::vector<Spline> ThetaL_ell_of_k_spline;
    
    //=====================================================================
    // [3] Integrate to get power-spectrum
    //=====================================================================

    Vector solve_for_C_ell(
        Vector & logk_array,
        std::vector<Spline> & f_ell, 
        std::vector<Spline> & g_ell);

    // Splines with the power-spectra
    Spline C_ell_TT_spline{"C_ell_TT_spline"};
    Spline C_ell_TE_spline{"C_ell_TE_spline"};
    Spline C_ell_EE_spline{"C_ell_EE_spline"};
    Spline C_ell_Nu_spline{"C_ell_Nu_spline"};
    Spline C_ell_lens_spline{"C_ell_lens_spline"};

  public:

    // Constructors
    PowerSpectrum() = delete;
    PowerSpectrum(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec, 
        Perturbations *pert,
        double A_s,
        double n_s,
        double kpivot_Mpc);
    
    // Do all the solving: bessel functions, LOS integration and then compute C_ells
    void solve(
        const double x_start,
        const double x_end,
        const int npts_x,
        const double k_start,
        const double k_end,
        const int npts_k);  // TODO: reasonable?

    // The dimensionless primordial power-spectrum Delta = 2pi^2/k^3 P(k) TODO which one is this? make private?
    double primordial_power_spectrum(const double k_Mpc) const;

    // Get P(k, x) for a given x in units of (Mpc)^3
    double get_matter_power_spectrum(const double x, const double k_Mpc) const;

    // Get the quantities we have computed TODO: int ell instead?
    double get_C_ell_TT(const double ell) const;
    double get_C_ell_TE(const double ell) const;
    double get_C_ell_EE(const double ell) const;
    double get_C_ell_Nu(const double ell) const;
    double get_C_ell_lens(const double ell) const;

    //TODO: get_k_eq / print_k_eq?

    // Print some useful info about the class
    void info() const;

    // Output matter power spectrum with conventional normalization
    void output_P_k(
        const double k_Mpc_min, 
        const double k_Mpc_max, 
        const std::string filename) const;

    // Output C_ells with conventional normalizations
    void output_C_ells(std::string filename) const;

    //TODO: correlation functions, CMB map, etc.
};

#endif
