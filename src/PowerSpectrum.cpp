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
    double kpivot_Mpc,
    int ell_max) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_Mpc(kpivot_Mpc),
  ell_max(ell_max)
{
  // Fill the first ells with tigther spacing
  for (size_t i = 0; i < first_ells.size(); i++) {
    if (first_ells[i] > ell_max) {
      ells.push_back(ell_max);
      break;
    }
    else
      ells.push_back(first_ells[i]);
  }

  // Fill the remaining ells 
  if (ell_max > 300) {
    for (int ell = 350; ell < ell_max; ell += 50)
      ells.push_back(ell);
    ells.push_back(ell_max);
  }

  // Fetch number of ells
  n_ells = ells.size();
}

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

  // Compute z_max from x_start and k_max
  const double eta_0 = cosmo->eta_of_x(0.0);
  const double eta   = cosmo->eta_of_x(x_start);
  const double z_max = k_max*(eta_0 - eta);
  if (z_max > 40000.0) {
    std::cout << "Current value of z_max is " << z_max << ", but it should not exceed 40000! Use a smaller value of |x_start| and/or k_max." << std::endl;
    exit(1);
  }

  // Generate splines for the Bessel functions
  generate_bessel_function_splines(z_max);

  // Line of sight integration to get the transfer functions Theta_ell(k)
  line_of_sight_integration(x_array, k_min, k_max, only_TT);

  // Integration to get C_ell 
  auto C_ell_TT = solve_for_C_ell(ThetaT_ell_of_k_spline, ThetaT_ell_of_k_spline, k_min, k_max, "TT");
  C_ell_TT_spline.create(ells, C_ell_TT, "C_ell_TT_of_ell");

  if (!only_TT) {
    // TODO comment
    if (pert->get_polarization_bool()) {
      auto C_ell_TE = solve_for_C_ell(ThetaT_ell_of_k_spline, ThetaE_ell_of_k_spline, k_min, k_max, "TE");
      C_ell_TE_spline.create(ells, C_ell_TE, "C_ell_TE_of_ell");

      auto C_ell_EE = solve_for_C_ell(ThetaE_ell_of_k_spline, ThetaE_ell_of_k_spline, k_min, k_max, "EE");
      C_ell_EE_spline.create(ells, C_ell_EE, "C_ell_EE_of_ell");
    }

    if (pert->get_neutrinos_bool()) {
      auto C_ell_nu = solve_for_C_ell(Nu_ell_of_k_spline, Nu_ell_of_k_spline, k_min, k_max, "nu");
      C_ell_nu_spline.create(ells, C_ell_nu, "C_ell_nu_of_ell");
    }

    if (pert->get_lensing_bool()) {
      auto C_ell_Psi = solve_for_C_ell(Psi_ell_of_k_spline, Psi_ell_of_k_spline, k_min, k_max, "Psi");
      C_ell_Psi_spline.create(ells, C_ell_Psi, "C_ell_Psi_of_ell");
    }

    if (angular_correlation) {
      // Set up the theta-array 
      Vector theta_array = Utils::linspace(0.0, M_PI, 10000); 

      // Compute the angular correlation function C(theta)
      compute_angular_correlation(theta_array);

      if (pert->get_lensing_bool() && lensed_TT) {
        // Generate splines for the reduced Wigner functions
        generate_reduced_Wigner_d_splines(theta_array);

        // Compute the lensed angular correlation function C^lensed(theta)
        compute_lensed_angular_correlation(theta_array);

        // Integrate to get the lensed TT-spectrum C^lensed_ell
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
void PowerSpectrum::generate_bessel_function_splines(const double z_max){
  Utils::StartTiming("besselspline");
  
  // Set the z-array
  double dz      = M_PI / 16.0;
  int npts_z     = static_cast<int>(z_max/dz);
  Vector z_array = Utils::linspace(0.0, z_max, npts_z);

  // Make storage for the splines
  j_ell_spline = std::vector<Spline>(n_ells);

  #pragma omp parallel for schedule(dynamic, 1) 
  for (int i = 0; i < n_ells; i++) {
    // Fetch current ell
    const int ell = ells[i];

    // Declare and fill array for j_ell(z)
    Vector j_ell_array(npts_z);
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
    std::function<double(double,double)> &source_function,
    std::string source)
{
  std::string label = "lineofsight(" + source + ")";
  Utils::StartTiming(label);

  int npts_k = k_array.size();
  int npts_x = x_array.size();

  // Make storage for the results
  Vector2D result = Vector2D(n_ells, Vector(npts_k));

  // Fetch cosmological parameters
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Fetch time step
  double dx = x_array[1] - x_array[0];

  #pragma omp parallel for schedule(dynamic, 1) 
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

      // Compute the integral
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

  Utils::EndTiming(label);
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
  ThetaT_ell_of_k_spline     = std::vector<Spline>(n_ells);
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
  Vector2D ThetaT_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_T, "temperature");

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
      ThetaE_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_E, "polarization");
    }

    // Solve for Nu_ell(k)
    if (neutrinos) {
      std::function<double(double,double)> source_function_nu = [&](double x, double k){
      return pert->get_source_nu(x, k);
      };
      Nu_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_nu, "neutrinos");

      // Add the term (Nu_0 + Psi)delta(eta) 
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

    // Solve for Psi_ell(k) 
    if (lensing) {
      std::function<double(double,double)> source_function_Psi = [&](double x, double k){
      return pert->get_source_Psi(x, k);
      };
      Psi_ell_of_k = line_of_sight_integration_single(x_array, k_array, source_function_Psi, "lensing");
    }
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
      if (lensing) 
        Psi_ell_of_k_spline[i].create(k_array, Psi_ell_of_k[i], "Psi_" + std::to_string(static_cast<int>(ells[i])) + "_of_k_spline");
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
    const double k_max,
    std::string spectrum)
{
  std::string label = "powerspectrum(" + spectrum + ")";
  Utils::StartTiming(label);

  // Make storage for the results
  Vector C_ell(n_ells);

  // Fetch cosmological parameters
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Set the k-array
  double dk  = M_PI / (64.0*eta_0); 
  int npts_k = static_cast<int>((k_max - k_min)/dk);
  Vector logk_array = Utils::linspace(log(k_min), log(k_max), npts_k);

  // Get dlogk
  double dlogk = logk_array[1] - logk_array[0];

  for (int i = 0; i < n_ells; i++) {
    // Declare relevant quantities
    double k;
    double P_prim;

    // Set integral to zero
    double integral = 0.0;

    // Compute the integral
    for (int ik = 0; ik < npts_k; ik++) {
      // Fetch k
      k = exp(logk_array[ik]);

      // Fetch the dimensionless primordial power spectrum
      P_prim = primordial_power_spectrum(k);

      // Set the derivative
      if (ik == 0 or ik == npts_k-1)
        integral += 2.0*M_PI*P_prim*f_ell_spline[i](k)*g_ell_spline[i](k);
      else
        integral += 4.0*M_PI*P_prim*f_ell_spline[i](k)*g_ell_spline[i](k);
    }

    // Store the result
    C_ell[i] = integral*dlogk;
  }

  Utils::EndTiming(label);
  return C_ell;
}


//====================================================
// Solve for the angular correlation function
//====================================================
void PowerSpectrum::compute_angular_correlation(Vector &theta_array)
{
  Utils::StartTiming("angularcorrelation");

  // Make array of ells 
  Vector ell_array = Utils::linspace(2, ell_max, ell_max-1);

  // Make storage for the result
  const int npts_theta = theta_array.size();
  Vector C_array(npts_theta);

  #pragma omp parallel for schedule(dynamic, 1)
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
// Generate splines for the reduced Wigner functions
//====================================================
void PowerSpectrum::generate_reduced_Wigner_d_splines(Vector &theta_array) {
  // Lambda function to make the rest less repetitive
  auto low_ell_Wigner_d = [=] (int ell, int m, int n, double theta) -> double {
    // Precompute log(cos(theta/2)) and log(sin(theta/2))
    const double log_cos_half = log(std::max(std::cos(theta / 2.0), 1e-15));
    const double log_sin_half = log(std::max(std::sin(theta / 2.0), 1e-15));

    // Precompute constant term
    const double log_term = (std::lgamma(ell + m + 1) + std::lgamma(ell - m + 1) +
                             std::lgamma(ell + n + 1) + std::lgamma(ell - n + 1)) / 2.0;
    
    // Declare relevant quantities
    double d_contrib;
    double d_pos = 0.0;
    double d_neg = 0.0;

    // Loop over summation index i with bounds depending on l, m, n
    for (int i = std::max(0, m - n); i <= std::min(ell + m, ell - n); ++i) { 
      d_contrib = exp(log_term -
                      std::lgamma(i + 1) - 
                      std::lgamma(ell + m - i + 1) -
                      std::lgamma(ell - n - i + 1) - 
                      std::lgamma(i + n - m + 1) +
                      (2*ell + m - n - 2*i) * log_cos_half + 
                      (2*i + n - m) * log_sin_half);

      // Sum positive and negative contributions separately to avoid cancellation error
      if (i % 2 == 0)
        d_pos += d_contrib;
      else
        d_neg += d_contrib;
    }
    return d_pos - d_neg;    
  };

  Utils::StartTiming("reducedWigner");

  // Set up arrays for computing the splines
  Vector ell_array = Utils::linspace(2, ell_max, ell_max-1);
  const int npts_theta = theta_array.size();
  Vector2D d_array(3);
  for (int i = 0; i < 3; i++)
    d_array[i] = Vector(npts_theta * (ell_max-1));

  // Maximum value of ell before we switch from log to recursion
  const int ell_max_log = 50;

  #pragma omp parallel for schedule(dynamic, 1) 
  for (int it = 0; it < npts_theta; it++) {
    // Fetch current theta and compute cosine
    const double theta    = theta_array[it];
    const double costheta = cos(theta);

    // Declare private thread index 
    int idx;

    // Declare the recursion coefficients
    double A;
    double B;
    double C;
    for (int ell = 2; ell < ell_max+1; ell++) {
      // Compute the index
      idx = ell-2 + (ell_max-1)*it;

      if (ell <= ell_max_log) {
        // Log-norm
        d_array[0][idx] = low_ell_Wigner_d(theta, ell, 1, 1);
        d_array[1][idx] = low_ell_Wigner_d(theta, ell, 1, -1);
        d_array[2][idx] = low_ell_Wigner_d(theta, ell, -1, 1);
      }
      else {
        // Compute recursion coefficients
        A = (ell*ell - 1.0) / (ell*(2.0*ell - 1.0));
        B = 1.0 / (ell*(ell - 1.0));
        C = ((ell - 1.0)*(ell - 1.0) - 1.0) / ((ell - 1.0)*(2.0*ell - 1.0));

        // Fill the values (d_-11 and d_1-1 have m*n = -1, hence + in front of B)
        d_array[0][idx] = ((costheta - B)*d_array[0][idx-1] - C*d_array[0][idx-2]) / A;
        d_array[1][idx] = ((costheta + B)*d_array[1][idx-1] - C*d_array[1][idx-2]) / A;
        d_array[2][idx] = ((costheta + B)*d_array[2][idx-1] - C*d_array[2][idx-2]) / A;
      }
    }
  }

  // Spline the results
  d_ell_pp_spline.create(ell_array, theta_array, d_array[0], "d_ell_pp_spline");
  d_ell_pm_spline.create(ell_array, theta_array, d_array[1], "d_ell_pm_spline");
  d_ell_mp_spline.create(ell_array, theta_array, d_array[2], "d_ell_mp_spline");

  Utils::EndTiming("reducedWigner");
}


//====================================================
// Solve for the lensed angular correlation function
//====================================================
void PowerSpectrum::compute_lensed_angular_correlation(Vector &theta_array)
{
  Utils::StartTiming("lensedangularcorrelation");

  // Make array of ells
  Vector ell_array = Utils::linspace(2, ell_max, ell_max-1);

  // Make storage for the results
  const int npts_theta = theta_array.size();
  Vector C_array(npts_theta);

  // Compute C_gl(0) beforehand
  double C_gl_of_0 = get_C_gl(0.0).first;

  // Compute C_gl and C_gl,2, and C(theta) from those
  #pragma omp parallel for schedule(dynamic, 1)
  for (int it = 0; it < npts_theta; it++) {
    // Fetch current theta and compute cos(theta)
    double theta    = theta_array[it];
    double costheta = cos(theta);

    // Compute C_gl and C_gl2
    auto C_gl_pair = get_C_gl(theta);
    double C_gl    = C_gl_pair.first;
    double C_gl2   = C_gl_pair.second;
    
    // Compute sigma^2 
    double sigma_squared = C_gl_of_0 - C_gl;
    
    // Compute C(theta)
    double C_ell_TT;
    double P_ell;
    double d_pm;
    double sum = 0.0;
    for (int ell = 2; ell < ell_max+1; ell++) {
      // Fetch the current value of the unlensed power spectrum
      C_ell_TT = get_C_ell_TT(ell);

      // Fetch the current value of the Legendre polynomial
      P_ell = gsl_sf_legendre_Pl(ell, costheta);

      // Fetch the reduced Wigner function
      d_pm = get_reduced_Wigner_d(ell, 1, -1, theta);

      // Add to the sum
      sum += (2.0*ell + 1.0)*C_ell_TT*exp(-ell*(ell + 1.0)*sigma_squared/2.0) *
             (P_ell + ell*(ell + 1.0)/2.0 * C_gl2 * d_pm);
    }

    // Store the result
    C_array[it] = sum / (4.0*M_PI);
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

  // Make array with thetas
  Vector theta_array = Utils::linspace(0.0, M_PI, npts);

  // Get dtheta
  double dtheta = theta_array[1] - theta_array[0];

  // Make storage for the result 
  Vector C_ell_lensed_array(n_ells);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < n_ells; i++) { 
    int ell = ells[i];

    // Set the integral to zero
    double integral = 0.0;

    // Compute the integral
    for (int it = 0; it < npts; it++) { 
      // Fetch current theta
      double theta = theta_array[it]; //TODO

      // Fetch the current values of C(theta) and P_ell(cos(theta))
      double C_of_theta = get_C_of_theta_lensed(theta);
      double P_ell = gsl_sf_legendre_Pl(ell, cos(theta)); 

      // Set the derivative
      if (it == 0 || it == npts-1)
        integral += C_of_theta * P_ell * sin(theta) / 2.0;
      else
        integral += C_of_theta * P_ell * sin(theta) ;
    }

    // Store the result 
    C_ell_lensed_array[i] = 2.0*M_PI * integral * dtheta;
  }

  // Spline the result
  C_ell_lensed_spline.create(ells, C_ell_lensed_array, "C_ell_lensed_spline");

  Utils::EndTiming("lensedspectrum");
}


//====================================================
// Solve for the correlation function
//====================================================
void PowerSpectrum::solve_correlation_function(Vector &r_array) {
  Utils::StartTiming("correlationfunction");

  // Cover a wide range of k's to avoid ringing
  double k_min = 1e-10/Constants.Mpc;
  double k_max = 1e2/Constants.Mpc;

  // Make storage for the results
  int npts_r = r_array.size();
  Vector xi_array(npts_r);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int ir = 0; ir < npts_r; ir++) {
    // Fetch current r
    const double r = r_array[ir];

    // Set the k-array
    double dk  = M_PI / (64.0*r); 
    int npts_k = static_cast<int>((k_max - k_min)/dk);
    Vector logk_array = Utils::linspace(log(k_min), log(k_max), npts_k);

    // Get dlogk
    double dlogk = logk_array[1] - logk_array[0];

    // Declare relevant quantities
    double k;
    double P_k;

    // Set integral to zero
    double integral = 0.0;

    // Compute the integral
    for (int ik = 0; ik < npts_k; ik++) {
      // Compute current k
      k = exp(logk_array[ik]);

      // Compute the dimensionless matter power spectrum in real space
      P_k = get_matter_power_spectrum(k, 0.0, false);

      // Set the derivative
      if (ik == 0 || ik == npts_k-1)
        integral += sin(k*r)/(k*r) * P_k/2.0;
      else
        integral += sin(k*r)/(k*r) * P_k;
    }

    // Store the result
    xi_array[ir] = integral*dlogk;
  }

  // Spline the result
  xi_spline.create(r_array, xi_array, "xi_spline");

  Utils::EndTiming("correlationfunction");
}


//====================================================
// The dimensionless primordial power-spectrum
//====================================================
double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow(k*Constants.Mpc/kpivot_Mpc, n_s - 1.0);
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_matter_power_spectrum(
    double k,
    const double x,
    const bool dimensionful) const
{
  // Fetch minimum and maximum k-values
  auto k_pair  = pert->get_k_min_max();
  double k_min = k_pair.first;
  double k_max = k_pair.second;

  // Use min/max values if outside k-range
  const double k_orig = k;
  if (k < k_min)
    k = k_min;
  else if (k > k_max)
    k = k_max;

    // Fetch quantities
  double Omega_m0 = cosmo->get_Omega_b() + cosmo->get_Omega_CDM();
  double H_0 = cosmo->get_H_0();
  double Phi = pert->get_Phi(x, k);

  // Compute the gauge invariant density perturbation
  double Delta_m = pow(Constants.c*k/H_0, 2.0) * 2.0*Phi/(3.0*Omega_m0) * exp(x);

  // Fetch the dimensionless primordial power spectrum
  double P_prim = primordial_power_spectrum(k);

  // Compute the total matter power spectrum
  double P_k = pow(Delta_m, 2.0) * P_prim;
  if (dimensionful)
    P_k *= 2.0*M_PI*M_PI * pow(k, -3.0);
  
  // Approximate extrapolations beyond k-range
  if (k_orig < k_min)
    P_k *= pow(k_orig / k_min, n_s); 
  else if (k_orig > k_max) {
    double P_k_prev = get_matter_power_spectrum(0.99*k_max, x, dimensionful);
    double log_slope = log(P_k / P_k_prev) / log(1.0 / 0.99);
    P_k *= pow(k_orig / k_max, log_slope);
  }
  
  return P_k;
}

std::pair<double,double> PowerSpectrum::get_C_gl(double theta) {  
  // Compute cos(theta)
  double costheta = cos(theta);

  // Set the sums to zero
  double C_gl  = 0.0;
  double C_gl2 = 0.0;

  // Compute the sum
  for (int ell = 2; ell < ell_max+1; ell++) { 
    // Fetch the current value of the lensing potential power spectrum
    double C_ell_Psi = get_C_ell_Psi(ell); 

    // Fetch the current value of the Legendre polynomial
    double P_ell = gsl_sf_legendre_Pl(ell, costheta);

    // Fetch the reduced Wigner functions
    double d_pp = get_reduced_Wigner_d(ell, 1, 1, theta);
    double d_mp = get_reduced_Wigner_d(ell, -1, 1, theta);

    double sgn_pp = 1.0;
    if (d_pp < 0.0)
      sgn_pp = -1.0;
    double sgn_mp = 1.0;
    if (d_mp < 0.0)
      sgn_mp = -1.0;  

    // Add to the sums
    C_gl  += (2.0*ell + 1.0)*ell*(ell + 1.0) * C_ell_Psi * d_pp;
    C_gl2 += (2.0*ell + 1.0)*ell*(ell + 1.0) * C_ell_Psi * d_mp; 
  }

  // Divide by final factor
  C_gl  /= (4.0*M_PI);
  C_gl2 /= (4.0*M_PI);

  return std::pair(C_gl, C_gl2);
}

double PowerSpectrum::get_reduced_Wigner_d(const double ell, const int m, const int n, const double theta) const{
  if (m == n)
    return d_ell_pp_spline(ell, theta);
  else
    if (m == 1)
      return d_ell_pm_spline(ell, theta);
    else
      return d_ell_mp_spline(ell, theta);
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
double PowerSpectrum::get_C_ell_Psi(const double ell) const{
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
void PowerSpectrum::info() const{
  std::cout << "\n";
  std::cout << "Info about power spectrum class:\n";
  std::cout << "A_s:            " << A_s         << "\n";
  std::cout << "n_s:            " << n_s         << "\n";
  std::cout << "kpivot (1/Mpc): " << kpivot_Mpc  << "\n";
  std::cout << "ell_max:        " << ell_max     << "\n";
  std::cout << "n_ells:         " << ells.size() << "\n";
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
// Output methods
//====================================================
void PowerSpectrum::output_transfer_functions(
        const double k_min, 
        const double k_max, 
        const int ell,
        const std::string filename,
        const bool only_TT) const
{
  std::ofstream fp(filename.c_str());
  const double h = cosmo->get_h();
  const int npts = static_cast<int>(log(k_max) - log(k_min))*10000 + 1; 
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), npts));

  // Output the product k*eta_0 for plotting the transfer function
  const double eta_0 = cosmo->eta_of_x(0.0);

  // Fetch the correct index for ell
  int ell_idx = ells.size();
  for (size_t i = 0; i < ells.size(); i++) {
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
    if (!only_TT) {
      if (pert->get_polarization_bool())
        fp << get_ThetaE_ell(k, ell_idx) << " ";
      if (pert->get_neutrinos_bool())
        fp << get_Nu_ell(k, ell_idx)     << " ";
      if (pert->get_lensing_bool())
        fp << get_Psi_ell(k, ell_idx)    << " ";
    }
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

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


void PowerSpectrum::output_C_ells(
        std::string filename,
        const bool only_TT,
        const bool lensed_TT) const
{
  std::ofstream fp(filename.c_str());
  Vector ell_array = Utils::linspace(2, ell_max, ell_max-1);

  auto print_data = [&] (const double ell) {
    double normfactor     = (ell * (ell+1.0)) / (2.0*M_PI) 
                            * pow(1e6*cosmo->get_T_CMB(), 2.0);
    double normfactor_TE  = sqrt((ell+2.0) * (ell+1.0) * ell * (ell-1.0)) * normfactor;
    double normfactor_EE  = (ell+2.0) * (ell+1.0) * ell * (ell-1.0)
                            * 1e5 * pow(1e6*cosmo->get_T_CMB(), 2.0);
    double normfactor_nu  = (ell * (ell+1.0)) / (2.0*M_PI)
                            * pow(1e6*cosmo->get_T_CMB() * pow(4.0/11.0, 1.0/3.0), 2.0);
    double normfactor_Psi = 1e7 * (ell * (ell+1.0)) * (ell * (ell+1.0)) / (2.0*M_PI);
    fp << ell                                         << " ";
    fp << C_ell_TT_spline(ell) * normfactor           << " ";
    if (!only_TT) {
      // TODO comment
      if (pert->get_polarization_bool()) {
        fp << C_ell_TE_spline(ell) * normfactor_TE    << " ";
        fp << C_ell_EE_spline(ell) * normfactor_EE    << " ";
      }
      if (pert->get_neutrinos_bool())
        fp << C_ell_nu_spline(ell) * normfactor_nu    << " ";
      if (pert->get_lensing_bool()) {
        fp << C_ell_Psi_spline(ell) * normfactor_Psi  << " "; 
        if (lensed_TT)
          fp << C_ell_lensed_spline(ell) * normfactor << " "; 
      }
    }
    fp << C_ell_lensed_spline(ell) * normfactor << " ";
    fp << "\n";
  };
  std::for_each(ell_array.begin(), ell_array.end(), print_data);
}


void PowerSpectrum::output_P_k(
        const double k_min_Mpc, 
        const double k_max_Mpc,  
        const std::string filename,
        const double x) const
{
  // Output k in units of h/Mpc and P(k) in units of (Mpc/h)^3
  std::ofstream fp(filename.c_str());
  const double h = cosmo->get_h();
  const int npts = static_cast<int>(log(k_max_Mpc) - log(k_min_Mpc))*10000 + 1; 
  Vector k_Mpc_array = exp(Utils::linspace(log(k_min_Mpc), log(k_max_Mpc), npts));

  auto print_data = [&] (const double k_Mpc) {
    fp << k_Mpc/h                                                                       << " ";
    fp << get_matter_power_spectrum(k_Mpc/Constants.Mpc, x) * pow(h/Constants.Mpc, 3.0) << " ";
    fp << "\n";
  };
  std::for_each(k_Mpc_array.begin(), k_Mpc_array.end(), print_data);
}

void PowerSpectrum::output_xi(
        const double r_min_Mpc, 
        const double r_max_Mpc, 
        const std::string filename) const
{
  // Output r in units of Mpc/h and xi(r) in units of (Mpc/h)^2
  std::ofstream fp(filename.c_str());
  const double h = cosmo->get_h();
  const int npts = static_cast<int>(r_max_Mpc - r_min_Mpc)*10000 + 1; 
  Vector r_Mpc_array = Utils::linspace(r_min_Mpc, r_max_Mpc, npts); 

  auto print_data = [&] (const double r_Mpc) {
    fp << r_Mpc*h                                       << " "; 
    fp << get_xi(r_Mpc*Constants.Mpc)*pow(r_Mpc*h, 2.0) << " "; 
    fp << "\n";
  };
  std::for_each(r_Mpc_array.begin(), r_Mpc_array.end(), print_data);
}