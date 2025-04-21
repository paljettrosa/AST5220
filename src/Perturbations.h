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

    // Include polarization and/or neutrinos?
    bool polarization;
    bool neutrinos;

    // Compute source function for lensing potential?
    bool lensing;
    
    // Number of scalars and multipoles
    int n_scalars        = Constants.n_scalars;
    int n_ell_Theta      = Constants.n_ell_Theta;
    int n_ell_Thetap     = Constants.n_ell_Thetap;
    int n_ell_Nu         = Constants.n_ell_Nu;

    // Indices
    int idx_Phi          = Constants.idx_Phi;
    int idx_deltaCDM     = Constants.idx_deltaCDM; 
    int idx_deltab       = Constants.idx_deltab;
    int idx_vCDM         = Constants.idx_vCDM;
    int idx_vb           = Constants.idx_vb;
    int idx_start_Theta  = Constants.idx_start_Theta;

    // Set by booleans
    int n_ell_tot;
    int idx_start_Thetap;
    int idx_start_Nu;
    
    // Neutrino fraction
    double f_nu;

    // Integrate perturbations and spline the result
    void integrate_perturbations(
        const double x_start, 
        const double x_end,
        const int npts_x,
        const double k_start,
        const double k_end,
        const int npts_k);
    
    // Compute source functions and spline the result
    void compute_source_functions(
        const double x_start, 
        const double x_end,
        const int npts_x,
        const double k_start,
        const double k_end,
        const int npts_k);
    
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
    std::vector<Spline2D> Theta_p_spline;
    std::vector<Spline2D> Nu_spline;
   
    // Splines of source functions
    Spline2D ST_spline{"ST_spline"};    // Temperature
    Spline2D SE_spline{"SE_spline"};    // Polarization
    Spline2D SN_spline{"SN_spline"};    // Neutrinos
    Spline2D SL_spline{"SL_spline"};    // Lensing

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec,
        bool polarization,
        bool neutrinos,
        bool lensing); 

    // Do all the solving
    void solve(
        const double x_start, 
        const double x_end,
        const int npts_x,
        const double k_start,
        const double k_end,
        const int npts_k);

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
    double get_Theta_p(const double x, const double k, const int ell) const;
    double get_Nu(const double x, const double k, const int ell) const;

    // Get the source functions
    double get_Source_T(const double x, const double k) const;
    double get_Source_E(const double x, const double k) const;
    double get_Source_N(const double x, const double k) const;
    double get_Source_L(const double x, const double k) const;

    // Get the class parameters
    bool get_polarization_bool() const;
    bool get_neutrinos_bool() const;
    bool get_lensing_bool() const;
    double get_neutrino_fraction() const;

    // Print some useful info about the class
    void info() const;

    // Output info to file
    void output(
        const double x_min, 
        const double x_max, 
        const double k, 
        const std::string filename) const;
};

#endif
