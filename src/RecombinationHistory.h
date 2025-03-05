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
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;


    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe and ne
    void solve_number_density_electrons(const double x_start, const double x_end, const int npts, bool timing = true);
  

    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // Solve for tau and g_tilde
    void solve_optical_depth_tau(const double x_start, const double x_end, const int npts, bool timing = true);


    //===============================================================
    // [3] Compute sound horizon s
    //===============================================================

    // Solve for s
    void solve_sound_horizon_s(const double x_start, const double x_end, const int npts, bool timing = true);


    // Splines contained in this class
    Spline Xe_of_x_spline{"Xe"};
    Spline log_ne_of_x_spline{"ne"};
    Spline tau_of_x_spline{"tau"}; 
    Spline g_tilde_of_x_spline{"g"};  
    Spline s_of_x_spline{"s"};

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp);

    // Do all the solving
    void solve(const double x_start, const double x_end, const int npts, bool Xe_ne = true, bool tau_g = true, bool s = false, bool timing = true);
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const double x_min, const double x_max, const std::string filename, bool s = false) const;

    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double Xe_of_x(double x) const;
    double ne_of_x(double x) const;
    double s_of_x(double x) const;
    double nb_of_x(double x) const;
    double get_Yp() const;
};

#endif
