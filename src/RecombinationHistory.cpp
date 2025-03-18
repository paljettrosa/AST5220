#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp,
    double z_reion,
    double Delta_z_reion,
    double z_Hereion,
    double Delta_z_Hereion) :
  cosmo(cosmo),
  Yp(Yp),
  z_reion(z_reion),
  Delta_z_reion(Delta_z_reion),
  z_Hereion(z_Hereion),
  Delta_z_Hereion(Delta_z_Hereion)
{

  // Derived parameters
  y_reion       = pow(1.0 + z_reion, 3.0/2.0);
  Delta_y_reion = 3.0/2.0 * sqrt(1.0 + z_reion) * Delta_z_reion;
  f_He          = Yp / (4.0*(1.0 - Yp));

  // Booleans
  if (Yp      != 0.0)       Helium = true;
  if (z_reion != 0.0) reionization = true;
}


//====================================================
// Main solver method
//====================================================
void RecombinationHistory::solve(
    const double x_start, 
    const double x_end, 
    const int npts,
    bool tau_g, 
    bool sound_horizon, 
    bool baryon_temp,
    double x_tol, 
    bool baryon_tau, 
    bool only_Saha,
    bool timing)
{
  // Compute and spline Xe, ne
  solve_number_density_electrons(x_start, x_end, npts, baryon_temp, x_tol, only_Saha, timing);
  
  if (tau_g) {
    // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx (also for baryons if baryon_tau is true)
    solve_optical_depth_tau(x_start, x_end, npts, baryon_tau, only_Saha, timing);
  }

  if (sound_horizon) {
    // Compute and spline s
    solve_sound_horizon_s(x_start, x_end, npts, timing);
  }
}


//=============================================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//=============================================================================
void RecombinationHistory::solve_number_density_electrons(
    const double x_start, 
    const double x_end, 
    const int npts,
    bool baryon_temp,
    double x_tol,
    bool only_Saha,
    bool timing)
{
  if (timing) Utils::StartTiming("Xe and ne");
  
  // Set up x-array and make arrays to store X_e(x) and n_e(x) (and y(x) in case we want this)
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector Xe_array (npts);
  Vector ne_array (npts);
  Vector y_array  (npts);

  // Set up arrays for no reionization if we include this
  Vector Xe_noreion_array  (npts);
  Vector ne_noreion_array  (npts);

  // Arrays for solving Peebles equation
  Vector Peebles_x_array;
  Vector Xe_ini;

  //======================================
  // Calculate recombination history
  //======================================
  bool Saha_regime           = true;
  bool switched_regimes      = false;
  bool reached_recombination = false;

  for (int i = 0; i < npts; i++) {

    // Get X_e from solving the Saha equation
    auto Xe_ne_data = electron_fraction_from_Saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Has recombination happened yet?
    if (Xe_current < Xe_recombination && !reached_recombination) {
      x_recombination       = x_array[i];
      reached_recombination = true;
    }

    // Are we still in the Saha regime?
    if (Xe_current < Xe_Saha_limit && !only_Saha) {
      Saha_regime = false;
      if (!switched_regimes && timing) {
        std::cout << "Switched to Peebles regime at x = " << x_array[i] << "\n";
        switched_regimes = true;
      }
    }

    if (Saha_regime) {
      // Store the result we got from the Saha equation
      Xe_array[i] = Xe_current;
      ne_array[i] = log(ne_current);

      // Store array without reionization (the same as with in Saha regime)
      if (reionization) {
        Xe_noreion_array[i] = Xe_current;
        ne_noreion_array[i] = log(ne_current);
      }

      // Store baryon temperature (still tightly coupled to photons)
      if (baryon_temp) y_array[i] = 0.0;
    } 
    
    else {
      // Compute X_e from current time til today by solving the Peebles equation

      ODESolver Peebles_Xe_ode;
      // Adjust accuracy if solving for baryon temperature
      if (baryon_temp) Peebles_Xe_ode.set_accuracy(1e-3, 1e-12, 1e-12);

      // The Peebles ODE equation
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_Peebles_ode(x, Xe, dXedx, baryon_temp, x_tol);
      };
      
      // Set up IC, solve the ODE and fetch the result 
      if (baryon_temp) Xe_ini = {Xe_array[i-1], 0.0};                  // Second IC corresponds to y = Tb/T - 1
      else             Xe_ini = {Xe_array[i-1]};
      Peebles_x_array         = Vector(x_array.begin() + i-1, x_array.end());

      Peebles_Xe_ode.solve(dXedx, Peebles_x_array, Xe_ini);

      auto Peebles_Xe_array = Peebles_Xe_ode.get_data_by_component(0);

      // Add the results to the arrays (and add reionization)
      double z;
      double y;
      double Xe_addition;
      double nH;

      for (int j = i; j < npts; j++) {
        Xe_array[j]           = Peebles_Xe_array[j-i+1];
        if (reionization) {
          Xe_noreion_array[j] = Peebles_Xe_array[j-i+1];
          z                   = 1.0/exp(x_array[j]) - 1.0;
          y                   = pow(1.0 + z, 3.0/2.0);
          Xe_addition         = (1.0 + f_He)/2.0 * (1.0 + tanh((y_reion - y)/Delta_y_reion));
          if (Helium) 
            Xe_addition      += f_He/2.0 * (1.0 + tanh((z_Hereion - z)/Delta_z_Hereion));
          Xe_array[j]        += Xe_addition;
        }

        // Has recombination happened yet?
        if (Xe_array[j] < Xe_recombination && !reached_recombination) {
          x_recombination       = x_array[j];
          reached_recombination = true;
        }

        nH                    = (1.0 - Yp)*nb_of_x(x_array[j]);
        ne_array[j]           = log(nH*Xe_array[j]);
        if (reionization) 
          ne_noreion_array[j] = log(nH*Xe_noreion_array[j]);
      }

      // Add the baryon temperature results if we solved for it
      if (baryon_temp) {
        auto Peebles_y_array = Peebles_Xe_ode.get_data_by_component(1);
        for (int j = i; j < npts; j++) {
          y_array[j] = Peebles_y_array[j-i+1];
        }
      }
      
      break; 
    }

    // Slow and unnecessary to solve for the entire range of x when using only Saha,
    // as we only need it to compute decoupling and recombination estimates
    if (only_Saha && reached_recombination) {
      if (abs(Xe_array[i] - Xe_array[i-1]) < 1e-7) {
        for (int j = i+1; j < npts; j++) {
          Xe_array[j]  = Xe_array[j-1] - abs(Xe_array[i] - Xe_array[i-1]);
          ne_array[j]  = ne_array[j-1] - abs(ne_array[i] - ne_array[i-1]); // Realistic decrease, for plotting purposes
        }
        break;
      }
    }
  }

  // Spline the result
  Xe_of_x_spline.create(x_array, Xe_array, "Xe_of_x");
  log_ne_of_x_spline.create(x_array, ne_array, "ne_of_x");
  if (reionization) {
    Xe_noreion_of_x_spline.create(x_array, Xe_noreion_array, "Xe_noreion_of_x");
    log_ne_noreion_of_x_spline.create(x_array, ne_noreion_array, "ne_noreion_of_x");
  }
  if (baryon_temp && !only_Saha) y_of_x_spline.create(x_array, y_array, "y_of_x");

  if (timing) Utils::EndTiming("Xe and ne");
}


//====================================================
// Solve the Saha equation(s) to get X_e and n_e
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_Saha_equation(double x) const{
  // Physical constants
  const double k_b       = Constants.k_b;
  const double m_e       = Constants.m_e;
  const double hbar      = Constants.hbar;
  const double epsilon_0 = Constants.epsilon_0;
  const double chi_0     = Constants.chi_0;
  const double chi_1     = Constants.chi_1;

  // Fetch cosmological parameters
  const double Tb     = cosmo->get_TCMB(x);

  // Compute baryon, Hydrogen and Helium number densities
  const double nb     = nb_of_x(x);
  const double nH     = (1.0 - Yp)*nb;
  const double nHe    = Yp*nb/4.0;

  // Compute common factor
  const double factor = pow(m_e*k_b*Tb / (2.0*pow(hbar, 2)*M_PI), 3.0/2.0);

  // Electron fraction and number density
  double Xe;
  double ne;

  if (Helium) {
    // Compute right-hand sides of Saha equations
    const double rhs_H   = factor * exp(-epsilon_0/(k_b*Tb));
    const double rhs_He  = 2.0 * factor * exp(-chi_0/(k_b*Tb));
    const double rhs_He2 = 4.0 * factor * exp(-chi_1/(k_b*Tb));

    // Iteratively solve for f_e
    double fe     = 1.0;
    double fe_old = 0.0;
    double tol    = 1e-10;

    double x_H;
    double x_He;
    double x_He2;
    while (abs(fe - fe_old) >= tol) {
      x_H    = 1.0 / (1.0 + fe*nb/rhs_H);
      x_He   = 1.0 / (1.0 + fe*nb/rhs_He + rhs_He2/(fe*nb));
      x_He2  = rhs_He2*x_He / (fe*nb);
      fe_old = fe;
      fe     = (2.0*x_He2 + x_He)*Yp/4.0 + x_H*(1.0 - Yp);
    }
    
    // Compute electron fraction and number density
    Xe = fe/(1.0 - Yp);
    ne = fe*nb;
  }

  else {
    // Compute coefficients (a = 1, c = -b)
    const double b   = 1.0/nb * factor * exp(-epsilon_0/(k_b*Tb));

    // Set tolerance for Taylor expanding the numerator
    const double tol = 1e-3;

    // Compute Xe and ne from the Saha equation
    Xe;
    if (4.0/b < tol) Xe = 1.0;
    else             Xe = b/2.0 * (-1.0 + sqrt(1.0 + 4.0/b));
    ne = Xe*nH;
  }

  return std::pair<double,double>(Xe, ne);
}


//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_Peebles_ode(
    double x, 
    const double *Xe, 
    double *dXedx, 
    bool baryon_temp,
    double x_tol)
{
  // Current value of Xe
  const double Xe_current   = Xe[0]; 

  // Physical constants in SI units
  const double k_b          = Constants.k_b;
  const double m_e          = Constants.m_e;
  const double c            = Constants.c;
  const double hbar         = Constants.hbar;
  const double sigma_T      = Constants.sigma_T;
  const double lambda_2s1s  = Constants.lambda_2s1s;
  const double epsilon_0    = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double H            = cosmo->H_of_x(x);
  const double TCMB         = cosmo->get_TCMB(x);

  double Tb                 = TCMB;
  // Use precise baryon temperature if we solve for it
  // The ODE for y is very unstable in the beginning of the Peebles regime
  if (baryon_temp && x > x_tol) Tb = (Xe[1] + 1.0)*TCMB;

  // Compute relevant quantities
  const double nH           = (1.0 - Yp)*nb_of_x(x);
  const double n_1s         = (1.0 - Xe_current)*nH;
  const double lambda_alpha = H * pow(3.0*epsilon_0*hbar/c, 3) / (pow(8.0*M_PI, 2)*n_1s);
  const double phi_2        = 0.448 * log(epsilon_0/(k_b*Tb));
  const double alpha        = m_e*c/hbar * sqrt(3*sigma_T / (8*M_PI));
  const double alpha2       = 64.0*M_PI/sqrt(27.0*M_PI) * pow(alpha*hbar/m_e, 2)/c * sqrt(epsilon_0/(k_b*Tb)) * phi_2;
  const double beta         = alpha2 * pow(m_e*k_b*Tb / (2.0*pow(hbar, 2)*M_PI), 3.0/2.0) * exp(-epsilon_0/(k_b*Tb));
  const double beta2        = alpha2 * pow(m_e*k_b*Tb / (2.0*pow(hbar, 2)*M_PI), 3.0/2.0) * exp(-epsilon_0/(4.0*k_b*Tb));
  const double C_r          = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2);

  // Update value for dXedx
  dXedx[0]                  = C_r/H * (beta*(1.0 - Xe_current) - nH*alpha2*pow(Xe_current, 2));

  // Update value for dydx if we solve for the baryon temperature
  if (baryon_temp) {
    const double y_current = Xe[1]; 
    const double mu        = Constants.m_H;
    const double OmegaR    = cosmo->get_OmegaR(0.0);
    const double OmegaB    = cosmo->get_OmegaB(0.0);

    double ne              = nH*Xe_current;
    double R               = 4.0*OmegaR / (3.0*OmegaB*exp(x));
    double dtaudx          = -c*ne*sigma_T/H;

    dXedx[1] = - 1.0 + y_current * (2.0*mu/m_e * R*dtaudx - 1.0);
  }

  return GSL_SUCCESS;
}


//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================
void RecombinationHistory::solve_optical_depth_tau(
    const double x_start, 
    const double x_end, 
    const int npts, 
    bool baryon_tau, 
    bool only_Saha,
    bool timing)
{
  if (timing) Utils::StartTiming("tau and g_tilde");

  // Physical constants in SI units
  const double c       = Constants.c;
  const double sigma_T = Constants.sigma_T;

  // Set up x-arrays to integrate over
  Vector x_array = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tau_ode;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Compute relevant quantities
    double H  = cosmo->H_of_x(x);
    double ne = ne_of_x(x);
    double ne_noreion;
    if (reionization) ne_noreion = ne_of_x(x, true);

    // Find correct index for tau_b
    int tau_b_idx;
    if (reionization) tau_b_idx = 2;
    else              tau_b_idx = 1;

    // Set the derivative for photon and baryon optical depth
    dtaudx[0]   = -c*ne*sigma_T/H;
    if (reionization) 
      dtaudx[1] = -c*ne_noreion*sigma_T/H;

    if (baryon_tau) { 
      // Fetch cosmological parameters
      const double OmegaR = cosmo->get_OmegaR(0.0);
      const double OmegaB = cosmo->get_OmegaB(0.0);

      double R  = 4.0*OmegaR / (3.0*OmegaB*exp(x));
      dtaudx[tau_b_idx] = R*dtaudx[0];
      if (reionization) 
        dtaudx[tau_b_idx+1] = R*dtaudx[1];
    }

    return GSL_SUCCESS;
  };

  // Set up IC, solve the ODE and fetch the result 
  Vector tau_ini;
  if (reionization && baryon_tau)
    tau_ini = {0.0, 0.0, 0.0, 0.0};  
  else if (reionization || baryon_tau)           
    tau_ini = {0.0, 0.0};
  else                
    tau_ini = {0.0};

  tau_ode.solve(dtaudx, x_array, tau_ini);
  auto tau_array = tau_ode.get_data_by_component(0);

  // Create spline of x(tau_noreion) for estimating decoupling time
  if (reionization) {
    auto tau_noreion_array = tau_ode.get_data_by_component(1);
    x_of_tau_spline.create(tau_noreion_array, x_array, "x_of_tau");
  }
  else
    x_of_tau_spline.create(tau_array, x_array, "x_of_tau");

  // Repeat process for x(tau_b_noreion) if we want to estime x_drag
  if (baryon_tau) {
    if (reionization) {
      auto tau_b_noreion_array = tau_ode.get_data_by_component(3);
      x_of_tau_b_spline.create(tau_b_noreion_array, x_array, "x_of_tau_b");
    }
    else {
      auto tau_b_array = tau_ode.get_data_by_component(1);
      x_of_tau_b_spline.create(tau_b_array, x_array, "x_of_tau_b");
    }
  }

  // Reverse order of arrays to make splines tau(x)
  std::reverse(x_array.begin(), x_array.end());
  std::reverse(tau_array.begin(), tau_array.end());

  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");

  // Repeat process for the baryon optical depth if needed
  if (baryon_tau) {
    int tau_b_idx;
    if (reionization) tau_b_idx = 2;
    else              tau_b_idx = 1;
    auto tau_b_array = tau_ode.get_data_by_component(tau_b_idx);

    std::reverse(tau_b_array.begin(), tau_b_array.end());
    tau_b_of_x_spline.create(x_array, tau_b_array, "tau_b_of_x");
  }

  // Compute visibility functions and spline everything
  if (!only_Saha) { // No need to do this when only using Saha
    Vector g_tilde_array (npts);
    Vector g_tilde_b_array (npts);
    for (int i = 0; i < npts; i++) {
      g_tilde_array[i] = -dtaudx_of_x(x_array[i]) * exp(-tau_of_x(x_array[i]));
      if (baryon_tau)
        g_tilde_b_array[i] = -dtaudx_of_x(x_array[i], true) * exp(-tau_of_x(x_array[i], true));
    }
    g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde_of_x");
    if (baryon_tau)
      g_tilde_b_of_x_spline.create(x_array, g_tilde_b_array, "g_tilde_b_of_x");


  }

  if (timing) Utils::EndTiming("tau and g_tilde");
}


//===================================================
// Solve for the sound horizon and spline the result
//===================================================
void RecombinationHistory::solve_sound_horizon_s(
    const double x_start, 
    const double x_end, 
    const int npts, 
    bool timing)
{
  if (timing) Utils::StartTiming("s");

  // Physical constants in SI units
  const double c      = Constants.c;

  // Fetch cosmological parameters
  const double OmegaR = cosmo->get_OmegaR(0.0);
  const double OmegaB = cosmo->get_OmegaB(0.0);

  // Set up x-array to integrate over
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // Declare relevant quantities
  double Hp;
  double R;
  double c_s;

  // The ODE system ds/dxx
  ODESolver s_ode;
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    // Compute relevant quantities
    Hp  = cosmo->Hp_of_x(x);
    R   = 4.0*OmegaR / (3.0*OmegaB*exp(x));
    c_s = c * sqrt(R / (3*(1 + R)));

    // Set the derivative for the sound horizon
    dsdx[0] = c_s/Hp;

    return GSL_SUCCESS;
  };

  // Set up and solve the ODE and make s splines
  Hp  = cosmo->Hp_of_x(x_start);
  R   = 4.0*OmegaR / (3.0*OmegaB*exp(x_start));
  c_s = c * sqrt(R / (3*(1 + R)));
  Vector s_ini {c_s/Hp};

  s_ode.solve(dsdx, x_array, s_ini);
  auto s_array = s_ode.get_data_by_component(0);

  // Make spline
  s_of_x_spline.create(x_array, s_array, "s_of_x");

  if (timing) Utils::EndTiming("s");
}


//====================================================
// Get methods
//====================================================
double RecombinationHistory::Xe_of_x(double x, bool no_reionization) const{
  if (no_reionization) return Xe_noreion_of_x_spline(x);
  else                 return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x, bool no_reionization) const{
  if (no_reionization) return exp(log_ne_noreion_of_x_spline(x));
  else                 return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::nb_of_x(double x) const{
  // Physical constants in SI units
  const double G      = Constants.G;
  const double m_H    = Constants.m_H;

  // Fetch cosmological parameters
  const double H0     = cosmo->get_H0();
  const double OmegaB = cosmo->get_OmegaB(0.0);

  // Compute baryon number density
  double rho_c0       = 3.0*pow(H0, 2) / (8.0*M_PI*G);
  double nb           = OmegaB*rho_c0 / (m_H*exp(3*x));
  return nb;
}

double RecombinationHistory::tau_of_x(double x, bool baryon_tau) const{
  if (baryon_tau) return tau_b_of_x_spline(x);
  else            return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x, bool baryon_tau) const{
  if (baryon_tau) return tau_b_of_x_spline.deriv_x(x);
  else            return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x, bool baryon_tau) const{
  if (baryon_tau) return tau_b_of_x_spline.deriv_xx(x);
  else            return tau_of_x_spline.deriv_xx(x); 
}

double RecombinationHistory::x_of_tau(double tau, bool baryon_tau) const{
  if (baryon_tau) return x_of_tau_b_spline(tau);
  else            return x_of_tau_spline(tau);
}

double RecombinationHistory::g_tilde_of_x(double x, bool baryon_tau) const{
  if (baryon_tau) return g_tilde_b_of_x_spline(x);
  else            return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x, bool baryon_tau) const{
  if (baryon_tau) return g_tilde_b_of_x_spline.deriv_x(x);
  else            return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x, bool baryon_tau) const{
  if (baryon_tau) return g_tilde_b_of_x_spline.deriv_xx(x);
  else            return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
}

double RecombinationHistory::Tb_of_x(double x, bool baryon_temp) const{
  double TCMB = cosmo->get_TCMB(x);
  if (baryon_temp) return (y_of_x_spline(x) + 1.0)*TCMB;
  else             return TCMB;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

double RecombinationHistory::get_z_reion() const{
  return z_reion;
}

double RecombinationHistory::get_Delta_z_reion() const{
  return Delta_z_reion;
}

double RecombinationHistory::get_z_Hereion() const{
  return z_Hereion;
}

double RecombinationHistory::get_Delta_z_Hereion() const{
  return Delta_z_Hereion;
}

double RecombinationHistory::get_x_recombination() const{
  if (x_recombination == 0.0) std::cout << "\nYou need to solve the system first!\n";
  return x_recombination;             
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:               " << Yp               << "\n";
  std::cout << "z_reion:          " << z_reion          << "\n";
  std::cout << "Delta_z_reion:    " << Delta_z_reion    << "\n";
  std::cout << "z_Hereion:        " << z_Hereion        << "\n";
  std::cout << "Delta_z_Hereion:  " << Delta_z_Hereion  << "\n";
  std::cout << std::endl;
} 

//====================================================
// Print out freeze-out abundance of free electrons
//====================================================
void RecombinationHistory::print_freeze_out_abundance() const{
  RecombinationHistory rec(cosmo, Yp, 0.0, 0.0, 0.0, 0.0);
  rec.solve(-13.0, 0.0, 1000, false, false, false, -7.0, false, false, false);

  std::cout << "\n";
  std::cout << "Freeze-out abundance of free electrons:\n";
  std::cout << "X_e = " << rec.Xe_of_x(0.0) << "\n";
  std::cout << "n_e = " << rec.ne_of_x(0.0) << "\n";
}

//====================================================
// Print out optical(s) debth at reionization
//====================================================
void RecombinationHistory::print_tau_reionization(bool baryon_tau) const{
  double x_reion = log(1.0 / (z_reion + 1.0));

  std::cout << "\n";
  std::cout << "Optical debth at reionization:\n";
  std::cout << "tau   = "                 << tau_of_x(x_reion)       << "\n";
  if (baryon_tau) std::cout << "tau_b = " << tau_of_x(x_reion, true) << "\n";
}

//====================================================
// Print out times and horizon sizes at decoupling 
// and recombination
//====================================================
void RecombinationHistory::print_decoupling_and_recombination(bool drag, bool Saha) const{
  double x_decoup;
  double x_drag;
  double x_recomb;
  if (Saha) {
      std::cout << "\n";
      std::cout << "==================================\n";
      std::cout << "Using only the Saha approximation:\n";
      std::cout << "==================================";
      RecombinationHistory rec(cosmo, Yp, 0.0, 0.0, 0.0, 0.0);
      if (drag) {
        rec.solve(-13.0, 0.0, 1000, true, false, false, -7.0, true, true, false);
        x_drag         = rec.x_of_tau(1.0, true);
      }
      else
        rec.solve(-13.0, 0.0, 1000, true, false, false, -7.0, false, true, false);
      x_decoup         = rec.x_of_tau(1.0);
      x_recomb         = rec.get_x_recombination();
  }
  else {  
      std::cout << "\n";
      std::cout << "==================================\n";
      std::cout << "Using both Saha and Peebles:\n";
      std::cout << "==================================";
      x_decoup         = x_of_tau(1.0);
      x_recomb         = x_recombination;
      if (drag) x_drag = x_of_tau(1.0, true);
  }

  // Photon decoupling
  double z_decoup      = 1.0/exp(x_decoup) - 1.0;
  double t_decoup      = cosmo->t_of_x(x_decoup);
  double eta_decoup    = cosmo->eta_of_x(x_decoup);
  double r_s           = s_of_x(x_decoup);

  std::cout << "\n";
  std::cout << "Photon decoupling:\n";
  std::cout << "x:        " << x_decoup                             << "\n";
  std::cout << "z:      "   << z_decoup                             << "\n";
  std::cout << "t:       "  << t_decoup/Constants.kyr               << " kyr\n";
  std::cout << "eta/c:   "  << eta_decoup/Constants.c/Constants.Myr << " Myr\n";
  std::cout << "eta:     "  << eta_decoup/Constants.Mpc             << " Mpc\n";
  std::cout << "s:       "  << r_s/Constants.Mpc                    << " Mpc\n";

  // Baryon decoupling
  if (drag) {
    double z_drag      = 1.0/exp(x_drag) - 1.0;
    double t_drag      = cosmo->t_of_x(x_drag);
    double eta_drag    = cosmo->eta_of_x(x_drag);
    double s_drag      = s_of_x(x_drag);

    std::cout << "\n";
    std::cout << "Baryon decoupling:\n";
    std::cout << "x:        " << x_drag                             << "\n";
    std::cout << "z:      "   << z_drag                             << "\n";
    std::cout << "t:       "  << t_drag/Constants.kyr               << " kyr\n";
    std::cout << "eta/c:   "  << eta_drag/Constants.c/Constants.Myr << " Myr\n";
    std::cout << "eta:     "  << eta_drag/Constants.Mpc             << " Mpc\n";
    std::cout << "s:       "  << s_drag/Constants.Mpc               << " Mpc\n";

    std::cout << "\n";
    std::cout << "Changes during the drag epoch:\n";
    std::cout << "dx:       " << x_drag-x_decoup                                 << "\n";
    std::cout << "dz:     "   << z_drag-z_decoup                                 << "\n";
    std::cout << "dt:      "  << (t_drag-t_decoup)/Constants.kyr                 << " kyr\n";
    std::cout << "deta/c:  "  << (eta_drag-eta_decoup)/Constants.c/Constants.Myr << " Myr\n";
    std::cout << "deta:     " << (eta_drag-eta_decoup)/Constants.Mpc             << " Mpc\n";
    std::cout << "ds:       " << (s_drag-r_s)/Constants.Mpc                      << " Mpc\n";
  }

  // Recombination
  double z_recomb      = 1.0/exp(x_recomb) - 1.0;
  double t_recomb      = cosmo->t_of_x(x_recomb);
  double eta_recomb    = cosmo->eta_of_x(x_recomb);
  double s_recomb      = s_of_x(x_recomb);

  std::cout << "\n";
  std::cout << "Recombination:\n";
  std::cout << "x:        " << x_recomb                             << "\n";
  std::cout << "z:      "   << z_recomb                             << "\n";
  std::cout << "t:       "  << t_recomb/Constants.kyr               << " kyr\n";
  std::cout << "eta/c:   "  << eta_recomb/Constants.c/Constants.Myr << " Myr\n";
  std::cout << "eta:     "  << eta_recomb/Constants.Mpc             << " Mpc\n";
  std::cout << "s:       "  << s_recomb/Constants.Mpc               << " Mpc\n";
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(
    const double x_min, 
    const double x_max, 
    const std::string filename, 
    bool tau_g, 
    bool sound_horizon, 
    bool baryon_temp, 
    bool baryon_tau) const
{
  std::ofstream fp(filename.c_str());
  const int npts = static_cast<int>(x_max - x_min)*10000 + 1;
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  auto print_data = [&] (const double x) {
    fp << x                          << " ";
    fp << Xe_of_x(x)                 << " ";
    fp << ne_of_x(x)                 << " ";
    if (tau_g) {
      fp << tau_of_x(x)                << " ";
      fp << dtaudx_of_x(x)             << " ";
      fp << ddtauddx_of_x(x)           << " ";
      fp << g_tilde_of_x(x)      << " ";
      fp << dgdx_tilde_of_x(x)   << " ";
      fp << ddgddx_tilde_of_x(x) << " ";
    }
    if (sound_horizon) {
      fp << s_of_x(x)                << " ";
    }
    if (baryon_temp) {
      fp << Tb_of_x(x)               << " ";
      fp << Tb_of_x(x, true)         << " ";
    }
    if (baryon_tau) {
      fp << tau_of_x(x, true)        << " ";
      fp << dtaudx_of_x(x, true)     << " ";
      fp << ddtauddx_of_x(x, true)   << " ";
      fp << g_tilde_of_x(x, true)        << " ";
      fp << dgdx_tilde_of_x(x, true)     << " ";
      fp << ddgddx_tilde_of_x(x, true)   << " ";
    }
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

