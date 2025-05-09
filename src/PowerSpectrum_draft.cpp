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
    const double k_min,
    const double k_max,
    const bool only_TT,
    const bool angular_correlation,
    const bool lensed_TT,
    const bool correlation_function,
    const double r_min,
    const double r_max,
    const int npts_r)
{
  // Set up the x-array
  Vector x_array = Utils::linspace(x_start, x_end, npts_x);

  // Compute z_max from x_start TODO correct?
  const double eta_0 = cosmo->eta_of_x(0.0);
  const double eta   = cosmo->eta_of_x(x_start);
  const double z_max = k_max*(eta_0 - eta);
  if (z_max > 40000.0) {
    std::cout << "Current value of z_max is " << z_max << ", but it should not exceed 40000! Use a smaller value of |x_start|." << std::endl;
    exit(1);
  }

  // Implement generate_bessel_function_splines
  generate_bessel_function_splines(z_max, npts_x);

  // Line of sight integration to get Theta_ell(k)
  line_of_sight_integration(x_array, k_min, k_max, only_TT);

  // Integration to get C_ell 
  auto C_ell_TT = solve_for_C_ell(ThetaT_ell_of_k_spline, ThetaT_ell_of_k_spline, k_min, k_max);
  C_ell_TT_spline.create(ells, C_ell_TT, "C_ell_TT_of_ell");

  if (!only_TT) {
    // TODO comment
    if (pert->get_polarization_bool()) {
      auto C_ell_TE = solve_for_C_ell(ThetaT_ell_of_k_spline, ThetaE_ell_of_k_spline, k_min, k_max);
      C_ell_TE_spline.create(ells, C_ell_TE, "C_ell_TE_of_ell");

      auto C_ell_EE = solve_for_C_ell(ThetaE_ell_of_k_spline, ThetaE_ell_of_k_spline, k_min, k_max);
      C_ell_EE_spline.create(ells, C_ell_EE, "C_ell_EE_of_ell");
    }

    if (pert->get_neutrinos_bool()) {
      auto C_ell_nu = solve_for_C_ell(Nu_ell_of_k_spline, Nu_ell_of_k_spline, k_min, k_max);
      C_ell_nu_spline.create(ells, C_ell_nu, "C_ell_nu_of_ell");
    }

    // if (pert->get_lensing_bool()) { TODO add back
    //   auto C_ell_Psi = solve_for_C_ell(Psi_ell_of_k_spline, Psi_ell_of_k_spline, k_min, k_max);
    //   C_ell_Psi_spline.create(ells, C_ell_Psi, "C_ell_Psi_of_ell");
    // }

    if (angular_correlation) {
      // Set up the theta-array
      Vector theta_array = Utils::linspace(0.0, M_PI, 100000); //TODO: 

      // Compute the angular correlation function C(theta)
      compute_angular_correlation(theta_array);

      if (pert->get_lensing_bool() && lensed_TT) {
        // Compute the lensed angular correlation function C^Theta(theta)
        compute_lensed_angular_correlation(theta_array);

        // Integrate to get the lensed TT-spectrum C^Theta_ell
        solve_lensed_spectrum(100000);
      }
    }
  }

  if (correlation_function) {
    // Set up the r-array
    Vector r_array = Utils::linspace(r_min, r_max, npts_r);

    // Integrate to get the correlation function xi(r)
    solve_correlation_function(r_array);
  }
}

//=========================================================================
// Generate splines of j_ell(z) needed for LOS integration
//=========================================================================
void PowerSpectrum::generate_bessel_function_splines(const double z_max, const int npts_z){
  Utils::StartTiming("besselspline");
  
  // Create z-array for splining
  Vector z_array = Utils::linspace(0.0, z_max, npts_z);

  // Declare array for j_ell(z)
  Vector j_ell_array(npts_z);

  // Make storage for the splines
  j_ell_spline = std::vector<Spline>(ells.size());

  // #pragma omp parallel for schedule(dynamic, 1) //TODO: 
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

  int n_ells = ells.size();
  int npts_k = k_array.size();
  int npts_x = x_array.size();

  // Make storage for the results
  Vector2D result = Vector2D(n_ells, Vector(npts_k));

  // Fetch cosmological parameters
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Fetch time step
  double dx = x_array[1] - x_array[0];

  #pragma omp parallel for schedule(dynamic, 1) //TODO: experiment with this
  for (int i = 0; i < n_ells; i++) {
    // Declare relevant quantities
    double S_tilde;
    double eta;
    double j_ell;
    double integral;

    for (int ik = 0; ik < npts_k; ik++) {
      // Fetch current value of k
      const double k = k_array[ik];

      // Set integral to zero
      integral = 0.0;

      // // Compute the integral
      for (int ix = 0; ix < npts_x; ix++) {
        S_tilde = source_function(x_array[ix], k);
        eta     = cosmo->eta_of_x(x_array[ix]);
        j_ell   = j_ell_spline[i](k*(eta_0 - eta));

        if (ix == 0 || ix == npts_x-1)
          integral += S_tilde*j_ell/2.0;
        else
          integral += S_tilde*j_ell;
      }

      // Store the result
      result[i][ik] = integral*dx;
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
    const double k_min,
    const double k_max,
    const bool only_TT)
{
  // Fetch booleans
  const bool polarization = pert->get_polarization_bool();
  const bool neutrinos    = pert->get_neutrinos_bool();
  const bool lensing      = pert->get_lensing_bool();

  // Fetch cosmological parameters
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Set the k-array
  double dk      = M_PI / (3.0*eta_0);
  int npts_k     = static_cast<int>((k_max - k_min)/dk);
  Vector k_array = Utils::linspace(k_min, k_max, npts_k);
  
  // Make storage for the splines we are to create
  const int n_ells = ells.size();
  ThetaT_ell_of_k_spline   = std::vector<Spline>(n_ells);
  if (!only_TT) {
    if (polarization)
      ThetaE_ell_of_k_spline = std::vector<Spline>(n_ells);
    if (neutrinos)
      Nu_ell_of_k_spline     = std::vector<Spline>(n_ells);
    if (lensing)
      Psi_ell_of_k_spline    = std::vector<Spline>(n_ells);
  }

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_source_T(x, k);
  };

  // Do the line of sight integration
  Vector2D ThetaT_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_T);

  Vector2D ThetaE_ell_of_k;
  Vector2D Nu_ell_of_k;
  Vector2D Psi_ell_of_k;
  if (!only_TT) {
    // TODO comment
    // Solve for ThetaE_ell(k)
    if (polarization) {
      std::function<double(double,double)> source_function_E = [&](double x, double k){
      return pert->get_source_E(x, k);
      };
      ThetaE_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_E);
    }

    // Solve for Nu_ell(k)
    if (neutrinos) {
      std::function<double(double,double)> source_function_nu = [&](double x, double k){
      return pert->get_source_nu(x, k);
      };
      Nu_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_nu);

      // Add the term (Nu_0 + Psi)delta(eta) TODO: correct?
      double Nu_0;
      double Psi;
      double arg;
      for (int ik = 0; ik < npts_k; ik++) {
        Nu_0 = pert->get_Nu(x_array[0], k_array[ik], 0);
        Psi  = pert->get_Psi(x_array[0], k_array[ik]);
        arg  = k_array[ik] * (eta_0 - cosmo->eta_of_x(x_array[0]));
        for (int i = 0; i < n_ells; i++)
          Nu_ell_of_k[i][ik] += (Nu_0 + Psi) * j_ell_spline[i](arg);
      }
    }

    // // Solve for Psi_ell(k) TODO add back
    // if (lensing) {
    //   std::function<double(double,double)> source_function_Psi = [&](double x, double k){
    //   return pert->get_source_Psi(x, k);
    //   };
    //   Psi_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_Psi);
    // }
  }

  // Spline the results
  for (int i = 0; i < n_ells; i++) {
    ThetaT_ell_of_k_spline[i].create(k_array, ThetaT_ell_of_k[i], "ThetaT_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
    if (!only_TT) {
      // TODO comment
      if (polarization)
        ThetaE_ell_of_k_spline[i].create(k_array, ThetaE_ell_of_k[i], "ThetaE_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
      if (neutrinos)
        Nu_ell_of_k_spline[i].create(k_array, Nu_ell_of_k[i], "Nu_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
      // if (lensing) TODO add back
      //   Psi_ell_of_k_spline[i].create(k_array, Psi_ell_of_k[i], "Psi_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
    }
  }
}

//====================================================
// Solve for the power spectrum
//====================================================
Vector PowerSpectrum::solve_for_C_ell(
    std::vector<Spline> &f_ell_spline,
    std::vector<Spline> &g_ell_spline,
    const double k_min,
    const double k_max)
{
  Utils::StartTiming("powerspectrum");

  // Make storage for the results
  int n_ells = ells.size();
  Vector C_ell(n_ells);

  // Fetch cosmological parameters
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Set the k-array
  double dk  = M_PI / (64.0*eta_0); //TODO: maybe change back to 32
  int npts_k = static_cast<int>((k_max - k_min)/dk);
  Vector logk_array = Utils::linspace(log(k_min), log(k_max), npts_k);

  // Get dlogk
  double dlogk = logk_array[1] - logk_array[0];

  // #pragma omp parallel for schedule(dynamic, 1) //TODO
  for (int i = 0; i < n_ells; i++) {
    // Fetch current ell
    const int ell = ells[i];

    // Declare relevant quantities
    double k;
    double Delta;

    // Set integral to zero
    double integral = 0.0;

    // Compute the integral
    for (int ik = 0; ik < npts_k; ik++) {
      // Fetch k
      k = exp(logk_array[ik]);

      // Compute the dimensionless primordial power spectrum
      Delta = primordial_power_spectrum(k*Constants.Mpc);

      // Set the derivative
      if (ik == 0 or ik == npts_k-1)
        integral += 2.0*M_PI*Delta*f_ell_spline[i](k)*g_ell_spline[i](k); //TODO why wrong with fabs?
      else
        integral += 4.0*M_PI*Delta*f_ell_spline[i](k)*g_ell_spline[i](k); //TODO why wrong with fabs?
    }

    // Store the result
    C_ell[i] = integral*dlogk;
  }

  Utils::EndTiming("powerspectrum");

  return C_ell;
}


//====================================================
// Solve for the angular correlation function
//====================================================
void PowerSpectrum::compute_angular_correlation(Vector &theta_array)
{
  Utils::StartTiming("angularcorrelation");

  // // Fetch number of ells
  // const int n_ells = ells.size();

  // Make array of ells
  const int ell_max = int(ells[ells.size()-1]);
  Vector ell_array  = Utils::linspace(2, ell_max, ell_max-1);

  // Make storage for the result
  const int npts_theta = theta_array.size();
  Vector C_array(npts_theta);

  // // Declare relevant quantities TODO
  // double costheta;
  // double sum;
  // int ell;
  // double C_ell;
  // double P_ell;

  #pragma omp parallel for schedule(dynamic, 1) //TODO
  for (int it = 0; it < npts_theta; it++) {
    // Compute cos(theta)
    double costheta = cos(theta_array[it]);

    // Set the sum to zero
    double sum = 0.0;

    // Declare relevant quantities
    double C_ell;
    double P_ell;

    // Compute the sum
    for (int ell = 2; ell < ell_max+1; ell++) {
      // Fetch the current value of the power spectrum
      C_ell = get_C_ell_TT(ell);

      // Fetch the current value of the Legendre polynomial
      P_ell = gsl_sf_legendre_Pl(ell, costheta);

      // Add to the sum
      sum += (2.0*ell + 1.0) * C_ell * P_ell;
    }

    // Store the result
    C_array[it] = sum / (4.0*M_PI);
  }

  // Spline the result
  C_of_theta_spline.create(theta_array, C_array, "C_of_theta_spline");

  Utils::EndTiming("angularcorrelation");
}


//====================================================
// Solve for the lensed angular correlation function
//====================================================
void PowerSpectrum::compute_lensed_angular_correlation(Vector &theta_array)
{
  Utils::StartTiming("lensedangularcorrelation");

  // // Fetch number of ells and largest ell TODO
  // const int n_ells  = ells.size();
  // const int ell_max = ells[n_ells-1];

  // Make array of ells
  const int ell_max = int(ells[ells.size()-1]);
  Vector ell_array  = Utils::linspace(2, ell_max, ell_max-1);

  // Make storage for the results
  const int npts_theta = theta_array.size();
  Vector C_array(npts_theta);

  // Declare relevant quantities TODO
  double theta;
  double costheta;   
  double sum;
  double C_gl_of_0;
  double C_gl;
  double C_gl2;
  double sigma_squared;
  // double C_ell_TT;
  // double C_ell_Psi;
  // double P_ell;
  // double d_pp; // d^ell_1,1    
  // double d_pm; // d^ell_1,-1  
  // double d_mp; // d^ell_-1,1 

  // Compute C_gl and C_gl,2, and C(theta) from those
  for (int it = 0; it < npts_theta; it++) {
    // Fetch current theta and compute cos(theta)
    theta      = theta_array[it];
    costheta   = cos(theta);

    // Set the sums to zero
    C_gl  = 0.0;
    C_gl2 = 0.0;

    // Compute the sum
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : C_gl, C_gl2)//TODO
    for (int ell = 2; ell < ell_max+1; ell++) {
      // Fetch the current value of the lensing power spectrum
      // double C_ell_Psi = get_C_ell_Psi(ell);
      double C_ell_Psi = get_C_ell_Psi(ell, true); //TODO change back

      // Fetch the current value of the Legendre polynomial
      double P_ell = gsl_sf_legendre_Pl(ell, costheta);

      // Fetch the reduced Wigner functions
      double d_pp = get_reduced_wigner_d(ell, 1, 1, theta);
      double d_mp = get_reduced_wigner_d(ell, -1, 1, theta);
      // std::cout << "d_pp(" << ell << ") = " << d_pp << std::endl;
      // std::cout << "d_mp(" << ell << ") = " << d_mp << std::endl;

      // Add to the sums
      C_gl  += (2.0*ell + 1.0)*ell*(ell + 1.0) * C_ell_Psi * d_pp;
      C_gl2 += (2.0*ell + 1.0)*ell*(ell + 1.0) * C_ell_Psi * d_mp;
    }

    // Divide by final factor
    C_gl  /= (4.0*M_PI);
    C_gl2 /= (4.0*M_PI);

    // Store C_gl(0)
    if (it == 0)
      C_gl_of_0 = C_gl;
    
    // Compute sigma^2
    sigma_squared = C_gl_of_0 - C_gl;
    
    // Compute C(theta)
    sum = 0.0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : sum)//TODO
    for (int ell = 2; ell < ell_max+1; ell++) {
      // Fetch the current value of the unlensed power spectrum
      double C_ell_TT = get_C_ell_TT(ell);

      // Fetch the current value of the Legendre polynomial
      double P_ell = gsl_sf_legendre_Pl(ell, costheta);

      // Fetch the reduced Wigner function
      double d_pm = get_reduced_wigner_d(ell, 1, -1, theta);
      // std::cout << "d_pm(" << ell << ") = " << d_pm << std::endl;

      // Add to the sum
      sum += (2.0*ell + 1.0)*C_ell_TT*exp(-ell*(ell + 1.0)*sigma_squared/2.0) *
             (P_ell + ell*(ell + 1.0)/2.0 * C_gl2 * d_pm);
    }

    // Store the result
    C_array[it] = sum / (4.0*M_PI);
    // std::cout << "C_array(" << theta_array[it] << ") = " << C_array[it] << std::endl;
  }

  // Spline the result
  C_of_theta_lensed_spline.create(theta_array, C_array, "C_of_theta_lensed_spline");

  Utils::EndTiming("lensedangularcorrelation");
}


//====================================================
// Solve for the lensed CMB power spectrum
//====================================================
void PowerSpectrum::solve_lensed_spectrum(const int npts)
{
  Utils::StartTiming("lensedspectrum");

  // TODO
  // // Make costheta_array
  // Vector costheta_array = Utils::linspace(-1.0, 1.0, npts);

  // // Get dcos(theta)
  // double dcostheta = costheta_array[1] - costheta_array[0];

  // Make theta_array
  Vector theta_array = Utils::linspace(0.0, M_PI, npts); //TODO

  // Get dtheta
  double dtheta = theta_array[1] - theta_array[0];

  // Make storage for the result //TODO
  // const int n_ells = ells.size();
  // Vector C_ell_lensed_array(n_ells);
  const int ell_max = int(ells[ells.size()-1]);
  Vector ell_array  = Utils::linspace(2, ell_max, ell_max-1);
  Vector C_ell_lensed_array(ell_max-1);

  // // Declare relevant quantities TODO
  // double theta;
  // double C_of_theta;
  // double P_ell;
  // double integral;

  #pragma omp parallel for schedule(dynamic, 1) //TODO
  for (int ell = 2; ell < ell_max+1; ell++) {
    // Set the integral to zero
    double integral = 0.0;

    // Compute the integral
    for (int it = 0; it < npts; it++) { 
      // Compute current theta
      // theta = acos(costheta_array[it]);
      // std::cout << "theta = " << theta*180/M_PI << std::endl;

      double theta = theta_array[it];

      // Fetch the current values of C(theta) and P_ell(theta)
      double C_of_theta = get_C_of_theta_lensed(theta);
      // if (theta > 3.0) //TODO: remove
      //   C_of_theta = get_C_of_theta(theta);
      // std::cout << "C(theta) = " << C_of_theta*1e10 << std::endl;
      // P_ell = gsl_sf_legendre_Pl(ell, costheta_array[it]);
      double P_ell = gsl_sf_legendre_Pl(ell, cos(theta));
      // std::cout << "P_ell = " <<  P_ell << std::endl;

      // Add to the integral
      integral += C_of_theta * P_ell * sin(theta);
    }

    // Store the result
    // C_ell_lensed_array[i] = 2.0*M_PI * integral * dcostheta;
    // C_ell_lensed_array[i] = 2.0*M_PI * integral * dtheta; //TODO
    C_ell_lensed_array[ell-2] = 2.0*M_PI * integral * dtheta;
  }

  // Spline the result
  C_ell_lensed_spline.create(ell_array, C_ell_lensed_array, "C_ell_lensed_spline");

  Utils::EndTiming("lensedspectrum");
}


//====================================================
// Solve for the correlation function
//====================================================
void PowerSpectrum::solve_correlation_function(
  Vector &r_array)
{
  Utils::StartTiming("correlationfunction");

  // Get k_min and set a large k_max
  double k_min = pert->get_k_min_max().first;
  double k_max = 10.0/Constants.Mpc;

  // Make storage for the results
  int npts_r = r_array.size();
  Vector xi_array(npts_r);

  // TODO remove?
  // Compute the distortion parameter
  double Omega_m0 = cosmo->get_Omega_CDM() + cosmo->get_Omega_b();
  double b_gal    = 1.5; //TODO: between 1.5 and 2.0
  double beta     = pow(Omega_m0, 0.55) / b_gal;

  // #pragma omp parallel for schedule(dynamic, 1) //TODO
  for (int ir = 0; ir < npts_r; ir++) {
    // Fetch current r
    const double r = r_array[ir];

    // Fetch cosmological parameters TODO: remove?
    const double eta_0 = cosmo->eta_of_x(0.0);

    // Set the k-array
    double dk  = M_PI / (64.0*r); //TODO: change to use eta_0?
    int npts_k = static_cast<int>((k_max - k_min)/dk);
    Vector logk_array = Utils::linspace(log(k_min), log(k_max), npts_k);

    // Get dlogk
    double dlogk = logk_array[1] - logk_array[0];

    // Declare relevant quantities
    double k;
    double P_k;
    double P_k_0; //TODO

    // Set integral to zero
    double integral = 0.0;

    // Compute the integral
    for (int ik = 0; ik < npts_k; ik++) {
      // Compute current k
      k = exp(logk_array[ik]);

      // Compute the dimensionless matter power spectrum in real space
      P_k = get_matter_power_spectrum(0.0, k*Constants.Mpc, false);

      // Set the derivative
      if (ik == 0 || ik == npts_k-1)
        integral += sin(k*r)/(k*r) * P_k/2.0;
      else
        integral += sin(k*r)/(k*r) * P_k;
      
      // // Compute the monopole of the redshift-space power spectrum
      // P_k_0 = (1.0 + 2.0*beta/3.0 + beta*beta/5.0) * P_k;

      // // Fetch the bessel function  
      // double j_0 = gsl_sf_bessel_j0(k * r);

      // // Set the derivative
      // if (ik == 0 || ik == npts_k-1)
      //   integral += pow(k, 3.0) * P_k_0 * j_0/2.0;
      // else
      //   integral += pow(k, 3.0) * P_k_0 * j_0;
    }

    // Store the result
    // xi_array[ir] = integral*dlogk / (2.0*M_PI*M_PI); //TODO change back
    xi_array[ir] = integral*dlogk;
  }

  // Spline the result
  xi_spline.create(r_array, xi_array, "xi_spline");

  Utils::EndTiming("correlationfunction");
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
    const double k_Mpc,
    const bool dimensionful,
    const bool total,
    const int component) const
{
  // Compute the dimensionful primordial power spectrum in Mpc^3
  double Delta    = primordial_power_spectrum(k_Mpc);
  double P_prim   = 2.0*M_PI*M_PI * pow(k_Mpc, -3.0) * Delta;

  if (total) {
    // Fetch quantities
    double Omega_m0 = cosmo->get_Omega_b() + cosmo->get_Omega_CDM();
    double H_0      = cosmo->get_H_0();
    double k_max = pert->get_k_min_max().second;
    double Phi;
    if (k_Mpc/Constants.Mpc > k_max) {
      // Use second order Taylor expansion to approximate Phi(x, k)
      Phi = pert->get_Phi(x, k_max);
      double dPhidk =   pert->get_dPhidk(x, k_max);
      double ddPhiddk = pert->get_ddPhiddk(x, k_max);
      Phi += dPhidk*(k_Mpc/Constants.Mpc - k_max) + ddPhiddk*pow(k_Mpc/Constants.Mpc - k_max, 2.0)/2.0;
      if (Phi < 0.0)
        Phi = 0.0;
    } 
    else
      Phi = pert->get_Phi(x, k_Mpc/Constants.Mpc);

    // Compute the gauge invariant density perturbation
    double Delta_m  = pow(Constants.c*(k_Mpc/Constants.Mpc)/H_0, 2.0) * 2.0*Phi/(3.0*Omega_m0) * exp(x);

    // Compute the total matter power spectrum
    if (dimensionful)
      return pow(Delta_m, 2.0) * P_prim;
    // Return k^3P(k)/2pi^2 for computing xi
    else
      return pow(Delta_m, 2.0) * Delta;
  }
  else 
  { //TODO: rewrite this to fit k_Mpc > k_max, or maybe just remove
    if (component > 3) {
      std::cout << "The value of component must either be 0 (CDM), 1 (baryons), 2 (photons) or 3 (neutrinos). Returning 0.0" << std::endl;
      return 0.0;
    }

    // Declare quantities
    double w_i;
    double delta_i;
    double v_i;
    
    if (component == 0) {
      // CDM
      w_i     = 0.0;
      delta_i = pert->get_delta_CDM(x, k_Mpc/Constants.Mpc);
      v_i     = pert->get_v_CDM(x, k_Mpc/Constants.Mpc);
    }
    else if (component == 1) {
      // Baryons
      w_i     = 0.0;
      delta_i = pert->get_delta_b(x, k_Mpc/Constants.Mpc);
      v_i     = pert->get_v_b(x, k_Mpc/Constants.Mpc);
    }
    else if (component == 2) {
      // Photons
      w_i     = 1.0/3.0;
      delta_i = 4.0*pert->get_Theta(x, k_Mpc/Constants.Mpc, 0);
      v_i     = -3.0*pert->get_Theta(x, k_Mpc/Constants.Mpc, 1);
    }
    else if (component == 3) {
      // Neutrinos
      w_i     = 1.0/3.0;
      delta_i = 4.0*pert->get_Nu(x, k_Mpc/Constants.Mpc, 0);
      v_i     = -3.0*pert->get_Nu(x, k_Mpc/Constants.Mpc, 1);
    }

    // Fetch cosmological quantities
    double Hp = cosmo->Hp_of_x(x);

    // Compute the gauge invariant density perturbation
    double Delta_i = delta_i - 3.0*(1.0 + w_i)*Hp / (Constants.c*k_Mpc/Constants.Mpc) * v_i;

    // Compute the contribution to the matter power spectrum
    return pow(Delta_i, 2.0) * P_prim;
  }
}

// TODO: ask Hans
double PowerSpectrum::get_reduced_wigner_d(int ell, int m, int mp, double theta) {
    if (abs(m) > 1 || abs(mp) > 1) {
        throw std::invalid_argument("Only implemented for |m|,|m'| ≤ 1");
    }

    // const double costheta = cos(theta/2.0);
    // const double sintheta = sin(theta/2.0);
    // double sum = 0.0;

    // // Loop over summation index i with bounds depending on l, m, mp
    // for (int i = std::max(0, m - mp); i <= std::min(ell + m, ell - mp); ++i) { // TODO: maybe change to gsl function
    //   double factor = exp((std::lgamma(ell + m + 1) + std::lgamma(ell - m + 1) + 
    //                        std::lgamma(ell + mp + 1) + std::lgamma(ell - mp + 1)) / 2.0 -
    //                        std::lgamma(i + 1) - std::lgamma(ell + m - i + 1) -
    //                        std::lgamma(ell - mp - i + 1) - std::lgamma(i + mp - m + 1));
    //   sum += pow(-1, i) * factor * pow(costheta, 2*ell + m - mp - 2*i) * pow(sintheta, 2*i + mp - m);
    // }

    // return sum;


    const int ell_max = 50;

    const double sintheta_min = 1e-10;
    const double delta = 1e-11;

    if (ell <= ell_max) {
      // Use the WignerDMatrix operator
      Quaternions::Quaternion R = Quaternions::Quaternion(cos(theta/2.0), 0.0, sin(theta/2.0), 0.0); // y-axis rotation
      SphericalFunctions::WignerDMatrix D(R);
      return std::real(D(ell, mp, m)); //TODO: mp before m?

      // double sum = 0.0;

      // // Loop over summation index i with bounds depending on l, m, mp
      // for (int i = std::max(0, m - mp); i <= std::min(ell + m, ell - mp); ++i) { // TODO: maybe change to gsl function
      //   double factor = exp((std::lgamma(ell + m + 1) + std::lgamma(ell - m + 1) + 
      //                       std::lgamma(ell + mp + 1) + std::lgamma(ell - mp + 1)) / 2.0 -
      //                       std::lgamma(i + 1) - std::lgamma(ell + m - i + 1) -
      //                       std::lgamma(ell - mp - i + 1) - std::lgamma(i + mp - m + 1));
      //   sum += pow(-1, i) * factor * pow(cos(theta/2.0), 2*ell + m - mp - 2*i) * pow(sin(theta/2.0), 2*i + mp - m);
      // }

      // return sum;
    }
    else {
      double sintheta = sin(theta);
      if (sintheta < 1e-14) {
        // Approximate as 0 or pi
        sintheta = theta;
        Quaternions::Quaternion R = (theta < M_PI_2) ? Quaternions::Quaternion(1.0, 0.0, 0.0, 0.0) : Quaternions::Quaternion(-1.0, 0.0, 0.0, 0.0);
        SphericalFunctions::WignerDMatrix D(R);
        return std::real(D(ell, mp, m)); //TODO: mp before m?;
      }
      else {
        // Use the Debye approximation
        double phi     = M_PI_4 * (2.0*m + 2.0*mp + 1.0);        
        double d_Debye = cos((ell + 1.0/2.0)*theta - phi) / sqrt(2.0*M_PI*(ell + 1.0/2.0)*sintheta);

        if (sintheta > sintheta_min)
          return d_Debye;
        else {
          // Blend with small-angle approximation
          double blend   = 1.0 / (1.0 + exp((sintheta_min - abs(sintheta)) / delta)); //TODO where does this come from?
          double d_small = (m == mp) ? (1.0 - 0.5*ell*(ell + 1.0) * theta*theta) : 0.0; //TODO where does this come from?
          return (1.0 - blend)*d_small + blend*d_Debye;
        }
      }
    // else {
    //   // Use the recurrence relation from Varshalovich et al.
    //   double d_ell_1 = std::real(D(ell_max, m, mp));   // d^{ell-1}_{m,m'} TODO: mp before m?
    //   double d_ell_2 = std::real(D(ell_max-1, m, mp)); // d^{ell-2}_{m,m'} TODO: mp before m?

    //   double d_ell;
    //   double A;
    //   double B;
    //   double C;
    //   for (int l = ell_max+1; l < ell+1; l++) {
    //     // Compute coefficients
    //     A = - sqrt(((l - 1.0)*(l - 1.0) - m*m)*((l - 1.0)*(l - 1.0) - mp*mp)) / ((l - 1.0)*(2.0*l - 1.0));  // coefficient for d^{ell-2}_{m,m'}
    //     B = cos(theta) - m*mp / ((l - 1.0)*l);                                                              // coefficient for d^{ell-1}_{m,m'}
    //     C = sqrt((l*l - m*m)*(l*l - mp*mp)) / (l*(2.0*l - 1.0));                                            // coefficient for d^{ell}_{m,m'}

    //     // Avoid dividing by zero
    //     // if (C < 1e-5)
    //     //   d_ell = 0.0;
    //     // else
    //     //   d_ell = (B * d_ell_1 + A * d_ell_2) / C;
    //     // if (C < 1e-2)
    //     //   C = 1e-2;
    //     d_ell = (B * d_ell_1 + A * d_ell_2) / C;

    //     // // Update previous Wigner functions
    //     // d_ell_2 = d_ell_1;
    //     // d_ell_1 = d_ell;

    //     if (std::isnan(d_ell)) {
    //       // Approximate with previous value
    //       d_ell   = d_ell_1;
    //       d_ell_2 = d_ell_1;
    //     }
    //     else {
    //       // Update previous Wigner functions
    //       d_ell_2 = d_ell_1;
    //       d_ell_1 = d_ell;
    //     }
    //   }
    //   return d_ell;
    // }
    }
}

// TODO remove this
void PowerSpectrum::read_Planck_data(const std::string filename) {
    // Open the file
    std::ifstream fp(filename.c_str());

    // String variable to store the read data
    std::string line;

    // Skip the first line
    std::getline(fp, line); 

    // Make storage for data
    Vector ell_array;
    Vector C_ell_Psi_array;

    double ell, C_ell_TT, C_ell_TE, C_ell_EE, C_ell_BB, C_ell_Psi;
    while (std::getline(fp, line)) {
        // Parse input TODO: correct?
        std::istringstream iss(line);
        iss >> ell >> C_ell_TT >> C_ell_TE >> C_ell_EE >> C_ell_BB >> C_ell_Psi;

        // Add to arrays
        double normfactor = (ell * (ell+1.0)) * (ell * (ell+1.0)) / (2.0*M_PI);
        ell_array.push_back(ell);
        C_ell_Psi_array.push_back(C_ell_Psi / normfactor);
    }

    // Close file
    fp.close();

    // Make spline
    C_ell_Psi_Planck_spline.create(ell_array, C_ell_Psi_array);
}

double PowerSpectrum::get_ThetaT_ell(const double k, const int ell_idx) const{
  return ThetaT_ell_of_k_spline[ell_idx](k);
}
double PowerSpectrum::get_ThetaE_ell(const double k, const int ell_idx) const{
  return ThetaE_ell_of_k_spline[ell_idx](k);
}
double PowerSpectrum::get_Nu_ell(const double k, const int ell_idx) const{
  return Nu_ell_of_k_spline[ell_idx](k);
}
double PowerSpectrum::get_Psi_ell(const double k, const int ell_idx) const{
  return Psi_ell_of_k_spline[ell_idx](k);
}
double PowerSpectrum::get_C_of_theta(const double theta) const{
  return C_of_theta_spline(theta);
}
double PowerSpectrum::get_C_of_theta_lensed(const double theta) const{
  return C_of_theta_lensed_spline(theta);
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
double PowerSpectrum::get_C_ell_nu(const double ell) const{
  return C_ell_nu_spline(ell);
}
double PowerSpectrum::get_C_ell_Psi(const double ell, const bool Planck) const{ //TODO change back to normal
  if (Planck)
    return C_ell_Psi_Planck_spline(ell);
  else
    return C_ell_Psi_spline(ell);
}
double PowerSpectrum::get_C_ell_lensed(const double ell) const{
  return C_ell_lensed_spline(ell);
}
double PowerSpectrum::get_xi(const double r) const{
  return xi_spline(r);
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
// Print the equality scale in Mpc^-1
//====================================================
void PowerSpectrum::print_equality_scale() const{
  // Fetch cosmological parameters
  double Omega_b0     = cosmo->get_Omega_b();
  double Omega_CDM0   = cosmo->get_Omega_CDM();
  double Omega_gamma0 = cosmo->get_Omega_gamma();
  double Omega_nu0    = cosmo->get_Omega_nu();

  // Compute the scale factor at radiation-matter equality
  double a_eq = (Omega_gamma0 + Omega_nu0) / (Omega_b0 + Omega_CDM0);

  // Compute the Hubble parameter
  double H_eq = cosmo->H_of_x(log(a_eq));

  // Compute and print the equality scale
  double k_eq = a_eq * H_eq / Constants.c;
  std::cout << "\n";
  std::cout << "Equality scale in 1/Mpc:\n";
  std::cout << "k_eq = " << k_eq*Constants.Mpc << "\n";
}   


//====================================================
// Output transfer functions to file
//====================================================
void PowerSpectrum::output_transfer_functions(
        const double k_min, 
        const double k_max, 
        const int ell,
        const std::string filename) const
{
  std::ofstream fp(filename.c_str());
  const double h = cosmo->get_h();
  const int npts = static_cast<int>(log(k_max) - log(k_min))*10000 + 1; 
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), npts));

  // Output the product k*eta_0 for plotting the transfer function
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Fetch the correct index for ell
  int ell_idx = ells.size();
  for (size_t i; i < ells.size(); i++) {
    if (ell == ells[i]) {
      ell_idx = i;
      break;
    }
  }
  if (ell_idx == ells.size()) {
    std::cout << "No transfer function computed for ell = " << ell << ", defaulting to ell = 2." << std::endl;
    ell_idx = 0;
  }

  auto print_data = [&] (const double k) {
    fp << k                            << " ";
    fp << k*eta_0                      << " ";
    fp << get_ThetaT_ell(k, ell_idx)   << " ";
    if (pert->get_polarization_bool())
      fp << get_ThetaE_ell(k, ell_idx) << " ";
    if (pert->get_neutrinos_bool())
      fp << get_Nu_ell(k, ell_idx)     << " ";
    // if (pert->get_lensing_bool()) TODO add back
    //   fp << get_Psi_ell(k, ell_idx)    << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

//====================================================
// Output the angular correlation functions to file
//====================================================
void PowerSpectrum::output_C_of_theta(
        const std::string filename,
        const bool lensed_C) const
{
  std::ofstream fp(filename.c_str());
  Vector theta_array = Utils::linspace(0.0, M_PI, 1001);

  auto print_data = [&] (const double theta) {
    fp << theta                            << " ";
    fp << get_C_of_theta(theta)            << " ";
    if (lensed_C)
      fp << get_C_of_theta_lensed(theta)   << " ";
    fp << "\n";
  };
  std::for_each(theta_array.begin(), theta_array.end(), print_data);
}


//====================================================
// Output the C_ells to file
//====================================================
void PowerSpectrum::output_C_ells(
        std::string filename,
        const bool only_TT,
        const bool lensed_TT) const
{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ell_max = int(ells[ells.size()-1]);
  Vector ell_array  = Utils::linspace(2, ell_max, ell_max-1);

  auto print_data = [&] (const double ell) {
    double normfactor      = (ell * (ell+1.0)) / (2.0*M_PI) 
                             * pow(1e6*cosmo->get_T_CMB(), 2.0);
    double normfactor_TE   = sqrt((ell+2.0) * (ell+1.0) * ell * (ell-1.0)) * normfactor;
    double normfactor_EE   = (ell+2.0) * (ell+1.0) * ell * (ell-1.0)
                             * 1e5 * pow(1e6*cosmo->get_T_CMB(), 2.0);
    double normfactor_nu   = pow(1e6*cosmo->get_T_CMB() * pow(4.0/11.0, 1.0/3.0), 2.0);
    double normfactor_lens = 1e7 * (ell * (ell+1.0)) * (ell * (ell+1.0)) / (2.0*M_PI);
    fp << ell                                         << " ";
    fp << C_ell_TT_spline(ell) * normfactor           << " ";
    if (!only_TT) {
      if (pert->get_polarization_bool()) {
        fp << C_ell_TE_spline(ell) * normfactor_TE    << " ";
        fp << C_ell_EE_spline(ell) * normfactor_EE    << " ";
      }
      if (pert->get_neutrinos_bool())
        fp << C_ell_nu_spline(ell) * normfactor_nu    << " ";
      if (pert->get_lensing_bool()) {
        // fp << C_ell_Psi_spline(ell) * normfactor_lens << " "; //TODO add back
        fp << C_ell_Psi_Planck_spline(ell) * normfactor_lens << " ";
        if (lensed_TT)
          fp << C_ell_lensed_spline(ell) * normfactor << " "; //TODO: correct norm-factor?
      }
    }
    fp << "\n";
  };
  std::for_each(ell_array.begin(), ell_array.end(), print_data);
}


//====================================================
// Output matter power spectrum to file
//====================================================
void PowerSpectrum::output_P_k(
        const double k_Mpc_min, 
        const double k_Mpc_max, 
        const std::string filename,
        const bool components) const
{
  // Output k in units of h/Mpc and P(k) in units of (Mpc/h)^3
  std::ofstream fp(filename.c_str());
  const double h     = cosmo->get_h();
  const int npts     = static_cast<int>(log(k_Mpc_max) - log(k_Mpc_min))*10000 + 1; 
  Vector k_Mpc_array = exp(Utils::linspace(log(k_Mpc_min), log(k_Mpc_max), npts));

  auto print_data = [&] (const double k_Mpc) {
    fp << k_Mpc/h                                             << " ";
    fp << get_matter_power_spectrum(0.0, k_Mpc) * pow(h, 3.0) << " ";
    if (components) 
      for (int i = 0; i < 4; i++)
        fp << get_matter_power_spectrum(0.0, k_Mpc, true, false, i) * pow(h, 3.0) << " ";
    fp << "\n";
  };
  std::for_each(k_Mpc_array.begin(), k_Mpc_array.end(), print_data);
}


//====================================================
// Output correlation function to file
//====================================================
void PowerSpectrum::output_xi(
        const double r_Mpc_min, 
        const double r_Mpc_max, 
        const std::string filename) const
{
  // Output r in units of Mpc/h and xi(r) in units of (Mpc)^2 TODO: or (Mpc/h)^2?
  std::ofstream fp(filename.c_str());
  const double h     = cosmo->get_h();
  const int npts     = static_cast<int>(log(r_Mpc_max) - log(r_Mpc_min))*10000 + 1; 
  Vector r_Mpc_array = exp(Utils::linspace(log(r_Mpc_min), log(r_Mpc_max), npts)); //TODO: maybe not log?

  auto print_data = [&] (const double r_Mpc) {
    fp << r_Mpc*h                                     << " "; 
    fp << get_xi(r_Mpc*Constants.Mpc)*pow(r_Mpc, 2.0) << " "; //TODO: other normalization? multiply by h? ask Hans
    fp << "\n";
  };
  std::for_each(r_Mpc_array.begin(), r_Mpc_array.end(), print_data);
}