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
    double OmegaB;                  // Baryon density today
    double OmegaCDM;                // CDM density today
    double OmegaLambda;             // Dark energy density today
    double Neff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB;                    // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double OmegaR;                  // Photon density today (follows from TCMB)
    double OmegaNu;                 // Neutrino density today (follows from TCMB and Neff)
    double OmegaK;                  // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc

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
        double OmegaB, 
        double OmegaCDM, 
        double OmegaK,
        double Neff, 
        double TCMB
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve(const double x_start, const double x_end, const int npts, bool eta = true, bool t = false, bool timing = true);

    // For printing the cosmic and conformal times at important values of x
    void print_times() const;

    // Output some results to file
    void output(const double x_min, const double x_max, const std::string filename, bool t = false, bool detadx = false, bool distances = false, bool TCMB = false) const;

    // Get functions that we must implement
    double eta_of_x(double x) const;
    double detadx_of_x(double x) const;
    double t_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double get_OmegaB(double x = 0.0) const; 
    double get_OmegaM(double x = 0.0) const; 
    double get_OmegaR(double x = 0.0) const;
    double get_OmegaRtot(double x = 0.0) const; 
    double get_OmegaNu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const; 
    double get_OmegaLambda(double x = 0.0) const; 
    double get_OmegaK(double x = 0.0) const; 
    double get_OmegaMnu(double x = 0.0) const; 
    double get_H0() const;
    double get_h() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;
    double get_r_of_x(double x) const;

    // Distance measures
    double get_comoving_distance_of_x(double x) const;
    double get_luminosity_distance_of_x(double x) const;
    double get_angular_diameter_distance_of_x(double x) const;

};

#endif
