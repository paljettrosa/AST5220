#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================
PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_Mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_Mpc(kpivot_Mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(
    const double x_start,
    const double x_end,
    const int npts_x,
    const double k_start,
    const double k_end,
    const int npts_k)
{
  // Make k-arrays
  // Vector k_array     = Utils::linspace(k_start, k_end, npts_k);
  // Vector logk_array  = log(k_array); //TODO should this be dimensionless?
  // Set up x-array and the k-array
  Vector x_array    = Utils::linspace(x_start, x_end, npts_x);
  Vector logk_array = Utils::linspace(log10(k_start), log10(k_end), npts_k); //TODO should this be dimensionless?
  Vector k_array(npts_k);
  for (int ik = 0; ik < npts_k; ik++)
    k_array[ik]     = pow(10.0, logk_array[ik]);

  // Compute z_max from x_start TODO correct?
  const double eta_0 = cosmo->eta_of_x(0.0);
  const double eta   = cosmo->eta_of_x(x_start);
  const double z_max = k_end*(eta_0 - eta);
  if (z_max > 40000.0) {
    std::cout << "Current value of z_max is " << z_max << ", but it should not exceed 40000! Use a smaller value of |x_start|." << std::endl;
    exit(1);
  }

  // Implement generate_bessel_function_splines
  generate_bessel_function_splines(z_max, npts_x);

  // TODO: maybe pass x_array?
  // Line of sight integration to get Theta_ell(k)
  line_of_sight_integration(x_array, k_array);

  //=========================================================================
  // TODO: Integration to get C_ell by solving dC_ell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_C_ell
  //=========================================================================
  auto C_ell_TT = solve_for_C_ell(logk_array, ThetaT_ell_of_k_spline, ThetaT_ell_of_k_spline);
  C_ell_TT_spline.create(ells, C_ell_TT, "C_ell_TT_of_ell");
  
  if (pert->get_polarization_bool()) {
    auto C_ell_TE = solve_for_C_ell(logk_array, ThetaT_ell_of_k_spline, ThetaE_ell_of_k_spline);
    C_ell_TE_spline.create(ells, C_ell_TE, "C_ell_TE_of_ell");

    auto C_ell_EE = solve_for_C_ell(logk_array, ThetaE_ell_of_k_spline, ThetaE_ell_of_k_spline);
    C_ell_EE_spline.create(ells, C_ell_EE, "C_ell_EE_of_ell");
  }

  if (pert->get_neutrinos_bool()) {
    auto C_ell_Nu = solve_for_C_ell(logk_array, Nu_ell_of_k_spline, Nu_ell_of_k_spline);
    C_ell_Nu_spline.create(ells, C_ell_Nu, "C_ell_Nu_of_ell");
  }

  if (pert->get_lensing_bool()) {
    auto C_ell_lens = solve_for_C_ell(logk_array, ThetaL_ell_of_k_spline, ThetaL_ell_of_k_spline);
    C_ell_lens_spline.create(ells, C_ell_lens, "C_ell_lens_of_ell");
  }
  
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================
void PowerSpectrum::generate_bessel_function_splines(const double z_max, const int npts_z){
  Utils::StartTiming("besselspline");
  
  // Create z-array for splining
  Vector z_array = Utils::linspace(0.0, z_max, npts_z);

  // Declare array for j_ell(z)
  Vector j_ell_array(npts_z);

  // Make storage for the splines
  j_ell_spline = std::vector<Spline>(ells.size());

  // TODO: correct?
  for (size_t i = 0; i < ells.size(); i++) {
    const int ell = ells[i];

    for (int iz = 0; iz < npts_z; iz++) 
      j_ell_array[iz] = Utils::j_ell(ell, z_array[iz]);
    
    // Make the j_ell_splines[i] spline
    j_ell_spline[i].create(z_array, j_ell_array, "j_" + std::to_string(ell) + "_spline");
  }

  Utils::EndTiming("besselspline");
}

//=========================================================================
// Do the line of sight integration for a single source function
//=========================================================================
Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector &x_array,
    Vector &k_array, 
    std::function<double(double,double)> &source_function)
{
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  for (size_t i = 0; i < ells.size(); i++) {
    const int ell = ells[i];

    for (size_t ik = 0; ik < k_array.size(); ik++) {
      // TODO: correct? evaluate integral directly?
      // Fetch current value of k
      const double k = k_array[ik];

      // Fetch cosmological parameters
      const double eta_0 = cosmo->eta_of_x(0.0);

      // Declare relevant quantities
      double S_tilde;
      double eta;
      double j_ell;

      // The line of sight ODE system
      ODESolver LOS_ode;
      ODEFunction df_elldx = [&](double x, const double *f_ell, double *df_elldx){
        // Compute relevant quantities
        S_tilde = source_function(x, k);
        eta     = cosmo->eta_of_x(x);
        j_ell   = j_ell_spline[i](k*(eta_0 - eta));

        // Set the derivative
        df_elldx[0] = S_tilde*j_ell;

        return GSL_SUCCESS;
      };

      // Set up initial conditions and integrate
      Vector f_ell_ini = {0.0};
      LOS_ode.solve(df_elldx, x_array, f_ell_ini);

      // Store results for splining
      result[i][ik] = LOS_ode.get_data_by_component(0)[x_array.size()-1]; //TODO: correct?
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(
    Vector &x_array,
    Vector &k_array)
{
  // Fetch booleans
  const bool polarization = pert->get_polarization_bool();
  const bool neutrinos    = pert->get_neutrinos_bool();
  const bool lensing      = pert->get_lensing_bool();

  // TODO: is this correct?
  const int n_ells = ells.size();
  
  // Make storage for the splines we are to create
  ThetaT_ell_of_k_spline   = std::vector<Spline>(n_ells);
  if (polarization)
    ThetaE_ell_of_k_spline = std::vector<Spline>(n_ells);
  if (neutrinos)
    Nu_ell_of_k_spline     = std::vector<Spline>(n_ells);
  if (lensing)
    ThetaL_ell_of_k_spline = std::vector<Spline>(n_ells);

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x, k);
  };

  // Do the line of sight integration
  Vector2D ThetaT_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_T);

  // Solve for ThetaE_ell(k)
  Vector2D ThetaE_ell_of_k;
  if (polarization) {
    std::function<double(double,double)> source_function_E = [&](double x, double k){
    return pert->get_Source_E(x, k);
    };
    ThetaE_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_E);
  }

  // Solve for Nu_ell(k)
  Vector2D Nu_ell_of_k;
  if (neutrinos) {
    std::function<double(double,double)> source_function_N = [&](double x, double k){
    return pert->get_Source_N(x, k);
    };
    Nu_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_N);
  }

  // Solve for ThetaL_ell(k)
  Vector2D ThetaL_ell_of_k;
  if (lensing) {
    std::function<double(double,double)> source_function_L = [&](double x, double k){
    return pert->get_Source_L(x, k);
    };
    ThetaL_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_L);
  }

  // Spline the results
  for (int i = 0; i < n_ells; i++) {
    ThetaT_ell_of_k_spline[i].create(k_array, ThetaT_ell_of_k[i], "ThetaT_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
    if (polarization)
      ThetaE_ell_of_k_spline[i].create(k_array, ThetaE_ell_of_k[i], "ThetaE_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
    if (neutrinos)
      Nu_ell_of_k_spline[i].create(k_array, Nu_ell_of_k[i], "Nu_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
    if (lensing)
      ThetaL_ell_of_k_spline[i].create(k_array, ThetaL_ell_of_k[i], "ThetaL_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
  }
}

//====================================================
// Compute the power spectrum
//====================================================
Vector PowerSpectrum::solve_for_C_ell(
    Vector &logk_array,
    std::vector<Spline> &f_ell_spline,
    std::vector<Spline> &g_ell_spline)
{
  Utils::StartTiming("powerspectrum");

  const int n_ells = ells.size();
  const int npts_k = logk_array.size();

  // Make storage for the results
  Vector C_ell(n_ells);

  for (int i = 0; i < n_ells; i++) {
    const int ell = ells[i];

    // Declare relevant quantities
    double k;
    double Delta;

    // The line of sight ODE system
    ODESolver powspec_ode;
    ODEFunction dC_elldlogk = [&](double logk, const double *C_ell, double *dC_elldlogk){
      // Compute k
      k = pow(10.0, logk);

      // Compute the dimensionless primordial power spectrum
      Delta = primordial_power_spectrum(k*Constants.Mpc);

      // Set the derivative
      dC_elldlogk[0] = 4.0*M_PI*Delta*f_ell_spline[i](k)*g_ell_spline[i](k);

      return GSL_SUCCESS;
    };

    // Set up initial conditions and integrate
    Vector C_ell_ini = {0.0}; //TODO: correct initial conditions?
    powspec_ode.solve(dC_elldlogk, logk_array, C_ell_ini);

    // Store last value for splining for splining
    C_ell[i] = powspec_ode.get_data_by_component(0)[npts_k-1]; //TODO correct? use last value since that is the value of the integral over all k's?
  }

  Utils::EndTiming("powerspectrum");

  return C_ell;
}

//====================================================
// The dimensionless primordial power-spectrum
//====================================================
double PowerSpectrum::primordial_power_spectrum(const double k_Mpc) const{
  return A_s * pow(k_Mpc/kpivot_Mpc, n_s - 1.0);
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_matter_power_spectrum(
    const double x, 
    const double k_Mpc) const
{
  // Compute the dimensionful primordial power spectrum
  double P_prim   = 2.0*M_PI*M_PI/pow(k_Mpc/Constants.Mpc, 3.0) * primordial_power_spectrum(k_Mpc);

  // Fetch quantities
  double Omega_m0 = cosmo->get_Omega_b() + cosmo->get_Omega_CDM();
  double H_0      = cosmo->get_H_0();
  double Phi      = pert->get_Phi(x, k_Mpc/Constants.Mpc);

  // Compute the gauge invariant density perturbation
  double Delta_m  = pow(Constants.c*(k_Mpc/Constants.Mpc)/H_0, 2.0) * 2.0*Phi/(3.0*Omega_m0) * exp(x);

  // Compute the matter power spectrum
  return pow(fabs(Delta_m), 2.0) * P_prim;
}

double PowerSpectrum::get_C_ell_TT(const double ell) const{
  return C_ell_TT_spline(ell);
}
double PowerSpectrum::get_C_ell_TE(const double ell) const{
  return C_ell_TE_spline(ell);
}
double PowerSpectrum::get_C_ell_EE(const double ell) const{
  return C_ell_EE_spline(ell);
}
double PowerSpectrum::get_C_ell_Nu(const double ell) const{
  return C_ell_Nu_spline(ell);
}
double PowerSpectrum::get_C_ell_lens(const double ell) const{
  return C_ell_lens_spline(ell);
}

//====================================================
// Print some useful info about the class
//====================================================

void PowerSpectrum::info() const{ //TODO: change?
  std::cout << "\n";
  std::cout << "Info about power spectrum class:\n";
  std::cout << "n_ells:         " << ells.size() << "\n";
  std::cout << "A_s:            " << A_s         << "\n";
  std::cout << "n_s:            " << n_s         << "\n";
  std::cout << "kpivot (1/Mpc): " << kpivot_Mpc  << "\n";
}

//====================================================
// Output matter power spectrum to file
//====================================================
void PowerSpectrum::output_P_k(
        const double k_Mpc_min, 
        const double k_Mpc_max, 
        const std::string filename) const
{
  // Output k in units of h/Mpc and P(k) in units of (Mpc/h)^3
  std::ofstream fp(filename.c_str());
  const double h     = cosmo->get_h();
  const int npts     = static_cast<int>(k_Mpc_max - k_Mpc_min)*10000 + 1; //TODO: divide by h here to get the correct amount of points?
  Vector k_Mpc_array = Utils::linspace(k_Mpc_min, k_Mpc_max, npts);

  auto print_data = [&] (const double k_Mpc) {
    fp << k_Mpc/h                                             << " ";
    fp << get_matter_power_spectrum(0.0, k_Mpc) * pow(h, 3.0) << " ";
    fp << "\n";
  };
  std::for_each(k_Mpc_array.begin(), k_Mpc_array.end(), print_data);
}

//====================================================
// Output the C_ells to file TODO: all ell-dependent quantities
//====================================================
void PowerSpectrum::output_C_ells(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ell_max = int(ells[ells.size()-1]);
  auto ell_array    = Utils::linspace(2, ell_max, ell_max-1);

  auto print_data = [&] (const double ell) {
    double normfactor  = (ell*(ell+1)) / (2.0*M_PI) 
                         * pow(1e6*cosmo->get_T_CMB(), 2.0);
    double normfactorN = (ell*(ell+1)) / (2.0*M_PI) 
                         * pow(1e6*cosmo->get_T_CMB() *  pow(4.0/11.0, 1.0/3.0), 2.0); //TODO: correct?
    double normfactorL = (ell*(ell+1))*(ell*(ell+1)) / (2.0*M_PI); //TODO: correct? what about cmb temp?
    fp << ell                                    << " ";
    fp << C_ell_TT_spline(ell) * normfactor      << " ";
    if (pert->get_polarization_bool()) {
      fp << C_ell_EE_spline(ell) * normfactor    << " ";
      fp << C_ell_TE_spline(ell) * normfactor    << " ";
    }
    if (pert->get_neutrinos_bool())
      fp << C_ell_Nu_spline(ell) * normfactorN   << " ";
    if (pert->get_lensing_bool())
      fp << C_ell_lens_spline(ell) * normfactorL << " ";
    fp << "\n";
  };
  std::for_each(ell_array.begin(), ell_array.end(), print_data);
}