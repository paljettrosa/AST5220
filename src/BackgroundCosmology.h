#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double Omega_b0;                // Baryon density today
    double Omega_CDM0;              // CDM density today
    double Omega_k0;                // Curvature density today
    double N_eff;                   // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double T_CMB0;                  // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double H_0;                     // The Hubble parameter today H0 = 100h km/s/Mpc
    double Omega_gamma0;            // Photon density today (follows from TCMB)
    double Omega_nu0;               // Neutrino density today (follows from TCMB and Neff)
    double Omega_Lambda0;           // Dark energy density today  = 1 - Omega_m0 - Omega_r0 - Omega_k0

    // Derived expressions for x
    double x_rm;                    // Radiation-matter equality
    double x_acc;                   // Onset of acceleration
    double x_mLambda;               // Matter-dark energy equality

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"t"};
 
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double Omega_b0, 
        double OmegaCDM0, 
        double Omega_k0,
        double N_eff, 
        double T_CMB0
        );

    // Do all the solving
    void solve(
        const double x_start, 
        const double x_end, 
        const int npts, 
        bool eta = true, 
        bool t = false, 
        bool timing = true);

    // Get functions
    double eta_of_x(double x) const;
    double detadx_of_x(double x) const;
    double t_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double get_Omega_b(double x = 0.0) const; 
    double get_Omega_gamma(double x = 0.0) const;
    double get_Omega_nu(double x = 0.0) const;
    double get_Omega_CDM(double x = 0.0) const; 
    double get_Omega_Lambda(double x = 0.0) const; 
    double get_Omega_k(double x = 0.0) const; 
    double get_H_0() const;
    double get_h() const;
    double get_N_eff() const;
    double get_T_CMB(double x = 0.0) const;
    double get_r_of_x(double x) const;

    // Distance measures
    double get_comoving_distance_of_x(double x) const;
    double get_luminosity_distance_of_x(double x) const;
    double get_angular_diameter_distance_of_x(double x) const;

    // Print some useful info about the class
    void info() const;

    // For printing the cosmic and conformal times at important values of x
    void print_times() const;

    // Output some results to file
    void output(
        const double x_min, 
        const double x_max, 
        const std::string filename, 
        bool t = false, 
        bool detadx = false, 
        bool distances = false, 
        bool TCMB = false) const;
};

#endif
