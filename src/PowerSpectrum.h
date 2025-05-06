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
    double A_s;
    double n_s;
    double kpivot_Mpc;
    
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
        const double k_min,
        const double k_max,
        const bool only_TT);
  
    // Do the line of sight integration for a single quantity
    Vector2D line_of_sight_integration_single(
        Vector &x_array,
        Vector &k_array, 
        std::function<double(double,double)> &source_function);
    
    // Splines of the reusult of the LOS integration
    std::vector<Spline> ThetaT_ell_of_k_spline;
    std::vector<Spline> ThetaE_ell_of_k_spline;
    std::vector<Spline> Nu_ell_of_k_spline;
    std::vector<Spline> Psi_ell_of_k_spline;
    
    //=====================================================================
    // [3] Integrate to get power-spectrum
    //=====================================================================

    Vector solve_for_C_ell(
        std::vector<Spline> & f_ell, 
        std::vector<Spline> & g_ell,
        const double k_min,
        const double k_max);

    // Splines with the power-spectra
    Spline C_ell_TT_spline{"C_ell_TT_spline"};
    Spline C_ell_TE_spline{"C_ell_TE_spline"};
    Spline C_ell_EE_spline{"C_ell_EE_spline"};
    Spline C_ell_nu_spline{"C_ell_nu_spline"};
    Spline C_ell_Psi_spline{"C_ell_Psi_spline"};


    //=====================================================================
    // [4] Compute the angular correlation functions
    //=====================================================================
    void compute_angular_correlation(Vector &theta_array);
    void compute_lensed_angular_correlation(Vector &theta_array);
    
    // Spline with the angular correlation functions
    Spline C_of_theta_spline{"C_of_theta_spline"};
    Spline C_of_theta_lensed_spline{"C_of_theta_lensed_spline"};


    //=====================================================================
    // [5] Compute the lensed CMB spectrum
    //=====================================================================
    void solve_lensed_spectrum(const int npts);
    
    // Spline with the lensed TT power-spectrum
    Spline C_ell_Theta_spline{"C_ell_Theta_spline"};


    //=====================================================================
    // [6] Compute the correlation function
    //=====================================================================
    void solve_correlation_function(
        Vector &r_array,
        const double k_min,
        const double k_max);
    
    // Spline with the correlation function
    Spline xi_spline{"xi_spline"};

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
        const double k_min,
        const double k_max,
        const bool only_TT = false,
        const bool angular_correlation = false,
        const bool lensed_TT = false,
        const bool correlation_function = false,
        const double r_min = Constants.Mpc,
        const double r_max = 500.0*Constants.Mpc,
        const int npts_r = 1000);  // TODO: reasonable?

    // The dimensionless primordial power-spectrum Delta = 2pi^2/k^3 P(k) TODO which one is this? make private?
    double primordial_power_spectrum(const double k_Mpc) const;

    // Get P(k, x) for a given x in units of (Mpc)^3 (or dimensionless)
    double get_matter_power_spectrum(
        const double x, 
        const double k_Mpc, 
        const bool dimensionful = true,
        const bool total = true,
        const int component = 0) const;

    // Get the quantities we have computed 
    double get_ThetaT_ell(const double k, const int ell_idx) const;
    double get_ThetaE_ell(const double k, const int ell_idx) const;
    double get_Nu_ell(const double k, const int ell_idx) const;
    double get_Psi_ell(const double k, const int ell_idx) const;
    double get_C_of_theta(const double theta) const;
    double get_C_of_theta_lensed(const double theta) const;
    double get_C_ell_TT(const double ell) const;
    double get_C_ell_TE(const double ell) const;
    double get_C_ell_EE(const double ell) const;
    double get_C_ell_nu(const double ell) const;
    double get_C_ell_Psi(const double ell) const;
    double get_C_ell_Theta(const double ell) const;
    double get_xi(const double r) const;

    // Print some useful info about the class
    void info() const;

    // Print the equality scale
    void print_equality_scale() const;    

    // Output the transfer functions
    void output_transfer_functions(
        const double k_min, 
        const double k_max,
        const int ell,
        const std::string filename) const; 
    
    // Output angular correlation functions
    void output_C_of_theta(
        const std::string filename,
        const bool lensed_C = false) const; 

    // Output C_ells with conventional normalizations
    void output_C_ells(
        std::string filename,
        const bool only_TT = false,
        const bool lensed_TT = false) const;

    // Output the matter power spectrum with conventional normalization
    void output_P_k(
        const double k_Mpc_min, 
        const double k_Mpc_max,
        const std::string filename,
        const bool components = false) const; 
    
    // Output correlation function with conventional normalization
    void output_xi(
        const double r_Mpc_min, 
        const double r_Mpc_max,
        const std::string filename) const; 

    //TODO: correlation functions, CMB map, etc.
};

#endif
