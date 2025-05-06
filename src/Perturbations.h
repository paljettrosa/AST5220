#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    // The cosmology and recombination history we use
    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;

    // Number of multipoles
    int n_ell_Theta;
    int n_ell_Theta_P;
    int n_ell_Nu;

    // Compute source function for lensing potential?
    bool lensing;

    // Include polarization and/or neutrinos?
    bool polarization = false;
    bool neutrinos    = false;
    
    // Number of scalars 
    int n_scalars = 5;

    // Indices
    int idx_Phi         = 0;
    int idx_delta_CDM   = 1; 
    int idx_delta_b     = 2;
    int idx_v_CDM       = 3;
    int idx_v_b         = 4;
    int idx_start_Theta = n_scalars;

    // Set by input parameters
    int n_ell_tot;
    int idx_start_Theta_P;
    int idx_start_Nu;
    
    // Neutrino fraction
    double f_nu;

    // Integrate perturbations and spline the result
    void integrate_perturbations(
        const double x_start, 
        const double x_end,
        const int npts_x,
        const double k_min,
        const double k_max,
        const int npts_k);
    
    // Compute source functions and spline the result
    void compute_source_functions(
        const double x_start, 
        const double x_end,
        const int npts_x,
        const double k_min,
        const double k_max,
        const int npts_k,
        const bool SW,
        const bool ISW,
        const bool Doppler,
        const bool pol);
    
    // Set the initial conditions in the very beginning
    Vector set_ic(
        const double x, 
        const double k) const;

    // Set the initial conditions after tight coupling ends
    Vector set_ic_after_tight_coupling(
        const Vector &y_tc,
        const double x, 
        const double k) const;
    
    // Find the index corresponding to when tight coupling ends
    int get_tight_coupling_index(
        const Vector &x_array, 
        const double k) const;
    
    // Right hand side of the ODE in the tight coupling regime
    int rhs_tight_coupling_ode(
        double x, 
        double k, 
        const double *y,
        double *dydx);
    
    // Right hand side of the ODE in the full regime
    int rhs_full_ode(
        double x, 
        double k, 
        const double *y, 
        double *dydx);

    // Splines of scalar perturbations quantities
    Spline2D Phi_spline{"Phi_spline"};
    Spline2D Psi_spline{"Psi_spline"};
    Spline2D Pi_spline{"Pi_spline"};
    Spline2D delta_CDM_spline{"delta_CDM_spline"};
    Spline2D delta_b_spline{"delta_b_spline"};
    Spline2D v_CDM_spline{"v_CDM_spline"};
    Spline2D v_b_spline{"v_b_spline"};

    // Splines of mulipole quantities
    std::vector<Spline2D> Theta_spline;
    std::vector<Spline2D> Theta_P_spline;
    std::vector<Spline2D> Nu_spline;
   
    // Splines of source functions
    Spline2D source_T_spline{"source_T_spline"};      // Temperature
    Spline2D source_E_spline{"source_E_spline"};      // Polarization
    Spline2D source_nu_spline{"source_nu_spline"};    // Neutrinos
    Spline2D source_Psi_spline{"source_Psi_spline"};  // Lensing

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec,
        int n_ell_Theta,
        int n_ell_Theta_P,
        int n_ell_Nu,
        bool lensing); 

    // Do all the solving
    void solve(
        const double x_start, 
        const double x_end,
        const int npts_x,
        const double k_min,
        const double k_max,
        const int npts_k,
        const bool SW = true,
        const bool ISW = true,
        const bool Doppler = true,
        const bool pol = true);

    // Get the quantities we have integrated
    double get_Phi(const double x, const double k) const;
    double get_dPhidx(const double x, const double k) const;
    double get_Psi(const double x, const double k) const;
    double get_dPsidx(const double x, const double k) const;
    double get_ddPsiddx(const double x, const double k) const;
    double get_Pi(const double x, const double k) const;
    double get_dPidx(const double x, const double k) const;
    double get_ddPiddx(const double x, const double k) const;
    double get_delta_CDM(const double x, const double k) const;
    double get_delta_b(const double x, const double k) const;
    double get_v_CDM(const double x, const double k) const;
    double get_v_b(const double x, const double k) const;
    double get_dv_bdx(const double x, const double k) const;
    double get_Theta(const double x, const double k, const int ell) const;
    double get_Theta_P(const double x, const double k, const int ell) const;
    double get_Nu(const double x, const double k, const int ell) const;

    // Get the source functions
    double get_source_T(const double x, const double k) const;
    double get_source_E(const double x, const double k) const;
    double get_source_nu(const double x, const double k) const;
    double get_source_Psi(const double x, const double k) const;

    // Get the class parameters
    bool get_polarization_bool() const;
    bool get_neutrinos_bool() const;
    bool get_lensing_bool() const;
    double get_neutrino_fraction() const;

    // Print some useful info about the class
    void info() const;

    // Compute and print the end of the tight-coupling regime for a given mode
    void print_tight_coupling_time(const double k) const;

    // Print the time when a mode enters the horizon
    void print_horizon_entry_time(const double k) const;

    // Output info to file
    void output(
        const double x_min, 
        const double x_max, 
        const double k, 
        const std::string filename,
        const bool source = false) const;
};

#endif
