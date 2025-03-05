#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(const double x_start, const double x_end, const int npts, bool Xe_ne, bool tau_g, bool s, bool timing){
  
  if (Xe_ne) {
    // Compute and spline Xe, ne
    solve_number_density_electrons(x_start, x_end, npts, timing);
  }
  
  if (tau_g) {
    // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx
    solve_optical_depth_tau(x_start, x_end, npts, timing);
  }

  if (s) {
    // Compute and spline s
    solve_sound_horizon_s(x_start, x_end, npts, timing);
  }
}

//=============================================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//=============================================================================

void RecombinationHistory::solve_number_density_electrons(const double x_start, const double x_end, const int npts, bool timing){
  if (timing) Utils::StartTiming("Xe and ne");
  
  //=============================================================================
  // Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector peebles_x_array; // TODO: integrate point by point instead?
  Vector Xe_array (npts);
  Vector ne_array (npts);

  // Calculate recombination history
  bool saha_regime = true;

  for (int i = 0; i < npts; i++) {

    //==============================================================
    // Get X_e from solving the Saha equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if (Xe_current < Xe_saha_limit)
      saha_regime = false;

    if (saha_regime) {
      //=============================================================================
      // Store the result we got from the Saha equation
      //=============================================================================
      Xe_array[i] = Xe_current;
      ne_array[i] = log(ne_current);
    } 
    
    else {
      //=============================================================================
      // Compute X_e from current time til today by solving the Peebles equation
      //=============================================================================

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      Vector Xe_ini = {Xe_array[i-1]};

      peebles_x_array = Vector(x_array.begin() + i-1, x_array.end());
      peebles_Xe_ode.solve(dXedx, peebles_x_array, Xe_ini);
      auto peebles_Xe_array = peebles_Xe_ode.get_data_by_component(0);
      // TODO integrate point by point instead?

      double nH;
      for (int j = i; j < npts; j++) {
        Xe_array[j] = peebles_Xe_array[j-i+1];
        nH = (1 + Yp)*nb_of_x(peebles_x_array[j-i+1]);
        ne_array[j] = log(peebles_Xe_array[j-i+1]*nH);
      }

      break;
    }
  }

  //=============================================================================
  // Spline the result.
  //=============================================================================
  Xe_of_x_spline.create(x_array, Xe_array, "Xe_of_x");
  log_ne_of_x_spline.create(x_array, ne_array, "ne_of_x");

  if (timing) Utils::EndTiming("Xe and ne");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  // Physical constants
  const double k_b       = Constants.k_b;
  const double m_e       = Constants.m_e;
  const double hbar      = Constants.hbar;
  const double epsilon_0 = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double Tb  = cosmo->get_TCMB(x);

  // Compute baryon and Hydrogen number densities
  const double nb  = nb_of_x(x);
  const double nH  = (1 - Yp)*nb;

  // Compute coefficients (a = 1, c = -b)
  const double b   = 1.0/nb * pow(m_e*k_b*Tb / (2.0*pow(hbar, 2)*M_PI), 3.0/2.0) * exp(-epsilon_0/(k_b*Tb));

  // Tolerance for Taylor expanding the numerator
  const double tol = 1e-3;

  // Electron fraction and number density
  double Xe;
  double n_e;
  //=============================================================================
  // Compute Xe and ne from the Saha equation
  //=============================================================================
  if (4.0/b < tol) { 
    Xe = 1.0;
  }
  else {
    Xe = b/2.0 * (-1.0 + sqrt(1.0 + 4.0/b));
  }
  n_e = Xe*nH;

  return std::pair<double,double>(Xe, n_e);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){
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
  const double Tb           = cosmo->get_TCMB(x);

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

  //=============================================================================
  // Update value for dXedx
  //=============================================================================
  dXedx[0] = C_r/H * (beta*(1 - Xe_current) - nH*alpha2*pow(Xe_current, 2));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_optical_depth_tau(const double x_start, const double x_end, const int npts, bool timing){
  if (timing) Utils::StartTiming("tau and g_tilde");

  // Physical constants in SI units
  const double c       = Constants.c;
  const double sigma_T = Constants.sigma_T;

  // Declare relevant quantities
  double H;
  double ne;

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionization
  // TODO: split into regions?
  Vector x_array = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tau_ode;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Compute relevant quantities
    H  = cosmo->H_of_x(x);
    ne = ne_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -c*ne*sigma_T/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE and make tau splines
  //=============================================================================
  Vector tau_ini {0.0};

  tau_ode.solve(dtaudx, x_array, tau_ini);
  auto tau_array = tau_ode.get_data_by_component(0);

  // Reverse order of arrays to make splines
  std::reverse(x_array.begin(), x_array.end());
  std::reverse(tau_array.begin(), tau_array.end());
  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");

  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================
  Vector g_tilde_array (npts);
  for (int i = 0; i < npts; i++) {
    g_tilde_array[i] = -dtaudx_of_x(x_array[i]) * exp(-tau_of_x(x_array[i]));
  }
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g_tilde_of_x");

  if (timing) Utils::EndTiming("tau and g_tilde");
}

//=============================================================================
// Solve for the sound horizon and spline the result
//=============================================================================

void RecombinationHistory::solve_sound_horizon_s(const double x_start, const double x_end, const int npts, bool timing){
  if (timing) Utils::StartTiming("s");

  // Physical constants in SI units
  const double c      = Constants.c;

  // Fetch cosmological parameters
  const double OmegaR = cosmo->get_OmegaR(0.0);
  const double OmegaB = cosmo->get_OmegaB(0.0);

  // Declare relevant quantities
  double Hp;
  double R;
  double c_s;

  // Set up x-array to integrate over
  Vector x_array = Utils::linspace(x_start, x_end, npts);

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

  //=============================================================================
  // Set up and solve the ODE and make s splines
  //=============================================================================
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

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x); 
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x); 
}

double RecombinationHistory::Xe_of_x(double x) const{
  // return exp(log_Xe_of_x_spline(x));
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
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

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const double x_min, const double x_max, const std::string filename, bool s) const{
  std::ofstream fp(filename.c_str());
  const int npts = static_cast<int>(x_max - x_min)*100 + 1;
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    if (s) {
      fp << s_of_x(x)          << " ";
    }
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

