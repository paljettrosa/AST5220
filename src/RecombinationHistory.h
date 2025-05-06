#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;

    // Reionization parameters
    double z_reion;
    double Delta_z_reion;
    double z_Hereion;
    double Delta_z_Hereion;

    // Derived reionization parameters
    double y_reion;
    double Delta_y_reion;
    double f_He;

    // Include Helium and/or reionization?
    bool Helium       = false;
    bool reionization = false;

    // x at recombination (determined when solving for Xe)
    double x_recombination        = 0.0;

    // Xe at recombination
    const double Xe_recombination = 0.1;


    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_Saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_Peebles_ode(
        double x, 
        const double *y, 
        double *dydx, 
        const bool baryon_temp,
        const double x_tol);
    
    // Solve for Xe and ne
    void solve_number_density_electrons(
        const double x_start, 
        const double x_end, 
        const int npts,
        const double Xe_Saha_limit, 
        const bool baryon_temp,
        const double x_tol,
        const bool only_Saha);
  

    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================
    void solve_optical_depth_tau(
        const double x_start, 
        const double x_end, 
        const int npts, 
        const bool baryon_tau, 
        const bool only_Saha);


    //===============================================================
    // [3] Compute sound horizon s
    //===============================================================
    void solve_sound_horizon_s(
        const double x_start, 
        const double x_end, 
        const int npts);


    // Splines contained in this class
    Spline Xe_of_x_spline{"Xe"};
    Spline log_ne_of_x_spline{"ne"};
    Spline Xe_noreion_of_x_spline{"Xe_noreion"};
    Spline log_ne_noreion_of_x_spline{"ne_noreion"};
    Spline tau_of_x_spline{"tau"};
    Spline tau_b_of_x_spline{"tau_b"}; 
    Spline x_of_tau_spline{"x_of_tau"};
    Spline x_of_tau_b_spline{"x_of_tau_b"};
    Spline g_tilde_of_x_spline{"g"};  
    Spline g_tilde_b_of_x_spline{"g_b"}; 
    Spline s_of_x_spline{"s"};
    Spline y_of_x_spline{"y"}; 

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp,
        double z_reion,
        double Delta_z_reion,
        double z_Hereion,
        double Delta_z_Hereion);

    // Do all the solving
    void solve(
        const double x_start, 
        const double x_end, 
        const int npts, 
        const double Xe_Saha_limit = 0.9999,
        const bool tau_g = true, 
        const bool sound_horizon = false, 
        const bool baryon_temp = false, 
        const double x_tol = -7.0,
        const bool baryon_tau = false, 
        const bool only_Saha = false);

    // Get functions
    double Xe_of_x(const double x, const bool no_reionization = false) const;
    double ne_of_x(const double x, const bool no_reionization = false) const;
    double nb_of_x(const double x) const;
    double tau_of_x(const double x, const bool baryon_tau = false) const;
    double dtaudx_of_x(const double x, const bool baryon_tau = false) const;
    double ddtauddx_of_x(const double x, const bool baryon_tau = false) const;
    double x_of_tau(const double tau, const bool baryon_tau = false) const;
    double g_tilde_of_x(const double x, const bool baryon_tau = false) const;
    double dgdx_tilde_of_x(const double x, const bool baryon_tau = false) const;
    double ddgddx_tilde_of_x(const double x, const bool baryon_tau = false) const;
    double s_of_x(const double x) const;
    double T_b_of_x(const double x, const bool baryon_temp = false) const;
    double get_Yp() const;
    double get_z_reion() const;
    double get_Delta_z_reion() const;
    double get_z_Hereion() const;
    double get_Delta_z_Hereion() const;
    double get_x_recombination() const;

    // Print some useful info about the class
    void info() const;

    // Compute and print the freeze-out abundance of free electrons
    void print_freeze_out_abundance() const;

    // Print the optical debth(s) at reionization
    void print_tau_reionization(const bool baryon_tau = false) const;

    // Print the times and horizon sizes for decoupling (also for baryons if drag = true) and recombination
    void print_decoupling_and_recombination(const bool drag = false, const bool Saha = false) const;

    // Output some data to file
    void output(
        const double x_min, 
        const double x_max, 
        const std::string filename, 
        const bool tau_g = true,
        const bool sound_horizon = false, 
        const bool baryon_temp = false, 
        const bool baryon_tau = false) const;
};

#endif
