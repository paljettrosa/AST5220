#include "Perturbations.h"

//====================================================
// Constructors
//====================================================
Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec,
    int n_ell_Theta,
    int n_ell_Theta_P,
    int n_ell_Nu,
    double k_min,
    double k_max,
    bool lensing) : // TODO: maybe just have lensing in solve-func?
  cosmo(cosmo), 
  rec(rec),
  n_ell_Theta(n_ell_Theta),
  n_ell_Theta_P(n_ell_Theta_P),
  n_ell_Nu(n_ell_Nu),
  k_min(k_min),
  k_max(k_max),
  lensing(lensing)
{
  // Compute total number of quantities
  n_ell_tot         = n_scalars + n_ell_Theta + n_ell_Theta_P + n_ell_Nu;

  // Booleans
  if (n_ell_Theta_P != 0) polarization = true;
  if (n_ell_Nu      != 0)    neutrinos = true;

  // Compute indices
  idx_start_Theta_P = idx_start_Theta   + n_ell_Theta   - (1 - polarization);
  idx_start_Nu      = idx_start_Theta_P + n_ell_Theta_P + (neutrinos - polarization);

  // Neutrino fraction
  double Omega_gamma0 = cosmo->get_Omega_gamma();
  double Omega_nu0    = cosmo->get_Omega_nu();
  f_nu = Omega_nu0 / (Omega_gamma0 + Omega_nu0) * neutrinos;
}

//====================================================
// Do all the solving
//====================================================
void Perturbations::solve(
    const double x_start, 
    const double x_end,
    const int npts_x,
    const int npts_k,
    const bool SW,
    const bool ISW,
    const bool Doppler,
    const bool pol)
{
  // Integrate all the perturbation equation and spline the result
  integrate_perturbations(x_start, x_end, npts_x, npts_k);

  // Compute source functions and spline the result
  compute_source_functions(x_start, x_end, npts_x, npts_k, SW, ISW, Doppler, pol);
}

//====================================================
// Integrate all perturbations and spline the results
//====================================================
void Perturbations::integrate_perturbations(
    const double x_start, 
    const double x_end,
    const int npts_x,
    const int npts_k)
{
  Utils::StartTiming("integrateperturbation");

  // Physical constants
  double c = Constants.c;

  // Set up arrays for computing the splines
  Vector2D y_array(n_ell_tot+2);
  for (int i = 0; i < n_ell_tot+2; i++)
    y_array[i] = Vector(npts_x * npts_k);

  // Set up x-array and the k-array
  Vector x_array = Utils::linspace(x_start, x_end, npts_x);
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), npts_k)); 

  // Loop over all wavenumbers
  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < npts_k; ik++){

    // Progress bar...
    if( (10*ik) / npts_k != (10*ik+10) / npts_k ) {
      std::cout << (100*ik+100)/npts_k << "% " << std::flush;
      if(ik == npts_k-1) std::cout << std::endl;
    }

    // Fetch current value of k
    const double k = k_array[ik];

    // Set up x-arrays for the two regimes
    int idx_tight_end    = get_tight_coupling_index(x_array, k);
    Vector x_array_tight = Vector(x_array.begin(), x_array.begin() + idx_tight_end+1);
    Vector x_array_full  = Vector(x_array.begin() + idx_tight_end, x_array.end());

    // The tight coupling ODE system
    ODESolver tight_ode;
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // The full ODE system
    ODESolver full_ode;
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Set up initial conditions in the tight coupling regime
    auto y_ini_tight = set_ic(x_start, k);

    // Integrate from x_start -> x_end_tight
    tight_ode.solve(dydx_tight_coupling, x_array_tight, y_ini_tight);

    // Store results for splining and set up initial conditions in the full regime
    Vector y_tc(n_ell_tot);
    for (int i = 0; i < n_ell_tot; i++) {
      auto y_tight_array = tight_ode.get_data_by_component(i);
      for (int ix = 0; ix < idx_tight_end; ix++)
        y_array[i][ix + npts_x*ik] = y_tight_array[ix];
      y_tc[i] = y_tight_array[idx_tight_end];
    }
    auto y_ini_full = set_ic_after_tight_coupling(y_tc, x_array_full[0], k);

    // Compute Theta_ell for ell>1 and Theta_P from initial conditions
    double Hp;
    double dtaudx;
    for (int ix = 0; ix < idx_tight_end; ix++) {
        Hp     = cosmo->Hp_of_x(x_array[ix]);
        dtaudx = rec->dtaudx_of_x(x_array[ix]);

        if (polarization)
          y_array[idx_start_Theta+2][ix + npts_x*ik]       = - 8.0*c*k / (15.0*Hp*dtaudx) * y_array[idx_start_Theta+1][ix + npts_x*ik];
        else
          y_array[idx_start_Theta+2][ix + npts_x*ik]       = - 20.0*c*k / (45.0*Hp*dtaudx) * y_array[idx_start_Theta+1][ix + npts_x*ik];
        for (int ell = 3; ell < n_ell_Theta; ell++) 
          y_array[idx_start_Theta+ell][ix + npts_x*ik]     = - ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * y_array[idx_start_Theta+ell-1][ix + npts_x*ik];

        if (polarization) {
          y_array[idx_start_Theta_P][ix + npts_x*ik]       = 5.0/4.0 * y_array[idx_start_Theta+2][ix + npts_x*ik];
          y_array[idx_start_Theta_P+1][ix + npts_x*ik]     = - c*k / (4.0*Hp*dtaudx) * y_array[idx_start_Theta+2][ix + npts_x*ik];
          y_array[idx_start_Theta_P+2][ix + npts_x*ik]     = 1.0/4.0 * y_array[idx_start_Theta+2][ix + npts_x*ik];
          for (int ell = 3; ell < n_ell_Theta_P; ell++)
            y_array[idx_start_Theta_P+ell][ix + npts_x*ik] = - ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * y_array[idx_start_Theta_P+ell-1][ix + npts_x*ik];
        }
    }

    // Integrate from x_end_tight -> x_end
    full_ode.solve(dydx_full, x_array_full, y_ini_full);

    // Store results for splining
    for (int i = 0; i < n_ell_tot; i++) {
      auto y_full_array  = full_ode.get_data_by_component(i);
      for (int ix = idx_tight_end; ix < npts_x; ix++)
          y_array[i][ix + npts_x*ik] = y_full_array[ix-idx_tight_end];
    }

    // Compute Psi and Pi
    double Omega_gamma0 = cosmo->get_Omega_gamma();
    double Omega_nu0    = cosmo->get_Omega_nu();
    double H_0          = cosmo->get_H_0();
    for (int ix = 0; ix < npts_x; ix++) {
      // Compute Psi
      if (neutrinos)
        y_array[n_ell_tot][ix + npts_x*ik] = - y_array[idx_Phi][ix + npts_x*ik]
                                             - 12.0*pow(H_0/(c*k), 2.0) * exp(-2.0*x_array[ix])
                                             * (Omega_gamma0 * y_array[idx_start_Theta+2][ix + npts_x*ik] 
                                             + Omega_nu0 * y_array[idx_start_Nu+2][ix + npts_x*ik]);
      else
        y_array[n_ell_tot][ix + npts_x*ik] = - y_array[idx_Phi][ix + npts_x*ik]
                                             - 12.0*pow(H_0/(c*k), 2.0) * exp(-2.0*x_array[ix])
                                             * (Omega_gamma0 * y_array[idx_start_Theta+2][ix + npts_x*ik]);

      // Compute Pi
      y_array[n_ell_tot+1][ix + npts_x*ik] = y_array[idx_start_Theta+2][ix + npts_x*ik];
      if (polarization)
        y_array[n_ell_tot+1][ix + npts_x*ik] += y_array[idx_start_Theta_P][ix + npts_x*ik]
                                                + y_array[idx_start_Theta_P+2][ix + npts_x*ik];
    }

  }
  Utils::EndTiming("integrateperturbation");

  //====================================================
  // Make all splines needed
  //====================================================
  Phi_spline.create(x_array, k_array, y_array[idx_Phi], "Phi_spline");
  Psi_spline.create(x_array, k_array, y_array[n_ell_tot], "Psi_spline");
  Pi_spline.create(x_array, k_array, y_array[n_ell_tot+1], "Pi_spline");
  delta_CDM_spline.create(x_array, k_array, y_array[idx_delta_CDM], "delta_CDM_spline");
  delta_b_spline.create(x_array, k_array, y_array[idx_delta_b], "delta_b_spline");
  v_CDM_spline.create(x_array, k_array, y_array[idx_v_CDM], "v_CDM_spline");
  v_b_spline.create(x_array, k_array, y_array[idx_v_b], "v_b_spline");

  Theta_spline = std::vector<Spline2D>(n_ell_Theta);
  for (int ell = 0; ell < n_ell_Theta; ell++)
    Theta_spline[ell].create(x_array, k_array, y_array[idx_start_Theta+ell], "Theta_" + std::to_string(ell) + "_spline");

  if (polarization)
    Theta_P_spline = std::vector<Spline2D>(n_ell_Theta_P);
    for (int ell = 0; ell < n_ell_Theta_P; ell++)
      Theta_P_spline[ell].create(x_array, k_array, y_array[idx_start_Theta_P+ell], "Theta_P_" + std::to_string(ell) + "_spline");
   
  if (neutrinos) {
    Nu_spline = std::vector<Spline2D>(n_ell_Nu);
    for (int ell = 0; ell < n_ell_Nu; ell++)
      Nu_spline[ell].create(x_array, k_array, y_array[idx_start_Nu+ell], "Nu_" + std::to_string(ell) + "_spline");
  }
}

//====================================================
// Compute the source function(s)
//====================================================
void Perturbations::compute_source_functions(
    const double x_start, 
    const double x_end,
    const int npts_x,
    const int npts_k,
    const bool SW,
    const bool ISW,
    const bool Doppler,
    const bool pol)
{
  Utils::StartTiming("source");

  // Physical constants
  double c = Constants.c;

  // Set up x-array and the k-array
  Vector x_array = Utils::linspace(x_start, x_end, npts_x);
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), npts_k));

  // Make storage for the source functions
  Vector source_T_array(npts_x * npts_k);    // Temperature
  Vector source_E_array(npts_x * npts_k);    // Polarization
  Vector source_nu_array(npts_x * npts_k);   // Neutrinos
  Vector source_Psi_array(npts_x * npts_k);  // Lensing

  // Fetch time-independent quantities
  double Omega_k0 = cosmo->get_Omega_k();
  double x_rec    = rec->get_x_recombination();
  double chi_rec  = cosmo->get_comoving_distance_of_x(x_rec);

  // Compute source functions
  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < npts_k; ik++){
    // Fetch current value of k
    const double k = k_array[ik];

    for(int ix = 0; ix < npts_x; ix++){
      // Fetch current value of x
      const double x = x_array[ix];

      // Cosmological variables
      double Hp           = cosmo->Hp_of_x(x);
      double dHpdx        = cosmo->dHpdx_of_x(x);
      double ddHpddx      = cosmo->ddHpddx_of_x(x);
      double chi          = cosmo->get_comoving_distance_of_x(x);

      // Recombination variables
      double tau          = rec->tau_of_x(x);
      double g_tilde      = rec->g_tilde_of_x(x);
      double dgdx_tilde   = rec->dgdx_tilde_of_x(x);
      double ddgddx_tilde = rec->ddgddx_tilde_of_x(x);

      // Perturbation variables
      double Phi          = get_Phi(x, k);
      double dPhidx       = get_dPhidx(x, k);
      double Psi          = get_Psi(x, k);
      double dPsidx       = get_dPsidx(x, k);
      double ddPsiddx     = get_ddPsiddx(x, k);
      double Pi           = get_Pi(x, k);
      double dPidx        = get_dPidx(x, k);
      double ddPiddx      = get_ddPiddx(x, k);
      double v_b          = get_v_b(x, k);
      double dv_bdx       = get_dv_bdx(x, k);
      double Theta_0      = get_Theta(x, k, 0);

      // Temperature source
      source_T_array[ix + npts_x*ik] = 0.0;
      if (SW)
        source_T_array[ix + npts_x*ik] += g_tilde * (Theta_0 + Psi + Pi/4.0);
      if (ISW)
        source_T_array[ix + npts_x*ik] += exp(-tau) * (dPsidx - dPhidx);
      if (Doppler)
        source_T_array[ix + npts_x*ik] -= 1.0/(c*k) * (dHpdx*g_tilde*v_b + Hp*dgdx_tilde*v_b + Hp*g_tilde*dv_bdx);
      if (pol)
        source_T_array[ix + npts_x*ik] += 3.0/4.0 * pow(c*k, -2.0) *
                                          (dHpdx * (dHpdx*g_tilde*Pi + Hp*dgdx_tilde*Pi + Hp*g_tilde*dPidx) +
                                           Hp * (ddHpddx*g_tilde*Pi + Hp*ddgddx_tilde*Pi + Hp*g_tilde*ddPiddx +
                                                 2.0 * (dHpdx*dgdx_tilde*Pi + dHpdx*g_tilde*dPidx + Hp*dgdx_tilde*dPidx)));

      // Polarization source
      if(polarization){
        if (x > -3.0) //TODO: Hans used -0.01
          source_E_array[ix + npts_x*ik] = source_E_array[ix-1 + npts_x*ik]; 
        else
          source_E_array[ix + npts_x*ik] = 3.0/4.0 * g_tilde*Pi / pow(k*chi, 2.0);
      }

      // Neutrino source (the delta-term is added after integrating to get N_ell)
      if (neutrinos)
        source_nu_array[ix + npts_x*ik] = dPsidx - dPhidx;

      // Lensing potential source
      if (lensing) {
        double source_Psi; 
        double W = 0.0;
        if (x >= x_rec) {
          if (Omega_k0 == 0)
            W = (chi -  chi_rec) / chi_rec;
          else if (Omega_k0 < 0)
            W = sin(chi - chi_rec) / sin(chi_rec);
          else
            W = sinh(chi - chi_rec) / sinh(chi_rec);
        }
        source_Psi = - 2.0*c*Psi / (Hp*chi) * W;
        if (std::isfinite(source_Psi)) 
          source_Psi_array[ix + npts_x*ik] = source_Psi;
        else
          source_Psi_array[ix + npts_x*ik] = source_Psi_array[ix-1 + npts_x*ik]; 
      }
    }
  }

  // Spline the source functions
  source_T_spline.create(x_array, k_array, source_T_array, "source_T_x_k");
  if (polarization)
    source_E_spline.create(x_array, k_array, source_E_array, "source_E_x_k");
  if (neutrinos)
    source_nu_spline.create(x_array, k_array, source_nu_array, "source_nu_x_k");
  if (lensing)
    source_Psi_spline.create(x_array, k_array, source_Psi_array, "source_Psi_x_k");

  Utils::EndTiming("source");
}

//====================================================
// Set the initial conditions in the very beginning
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{
  // The vector we are going to fill
  Vector y_ini(n_ell_tot);

  // References to the quantities we are going to set
  double &Phi               =  y_ini[idx_Phi];
  double &delta_CDM         =  y_ini[idx_delta_CDM];
  double &delta_b           =  y_ini[idx_delta_b];
  double &v_CDM             =  y_ini[idx_v_CDM];
  double &v_b               =  y_ini[idx_v_b];
  double *Theta             = &y_ini[idx_start_Theta];
  double *Theta_P           = &y_ini[idx_start_Theta_P];
  double *Nu                = &y_ini[idx_start_Nu];

  // Physical constants
  double c = Constants.c;

  // Cosmological quantities
  double Omega_nu0 = cosmo->get_Omega_nu();
  double H_0       = cosmo->get_H_0();
  double Hp        = cosmo->Hp_of_x(x);

  // Scalar quantities (Gravitational potental, baryons and CDM)
  double Psi     = - 1.0 / (3.0/2.0 + 2.0*f_nu/5.0);
  Phi            = - (1.0 + 2.0*f_nu/5.0) * Psi;
  delta_CDM      = - 3.0/2.0 * Psi;
  delta_b        = delta_CDM;
  v_CDM          = - c*k / (2.0*Hp) * Psi;
  v_b            = v_CDM; 

  // Photon temperature perturbations (ell >= 2 not solved for in tight-coupling)
  *Theta         = - 1.0/2.0 * Psi;
  *(Theta+1)     = c*k / (6.0*Hp) * Psi;
  for (int ell = 2; ell < n_ell_Theta; ell++) 
    *(Theta+ell) = 0.0;

  // Photon polarization perturbations (not solved for in tight-coupling)
  if (polarization) {
    for (int ell = 0; ell < n_ell_Theta_P; ell++)
      *(Theta_P+ell) = 0.0;
  }

  // Neutrino perturbations
  if (neutrinos) {
    *Nu     = *Theta;
    *(Nu+1) = *(Theta+1);
    *(Nu+2) = - pow(c*k/H_0, 2.0) * exp(2.0*x) * (Phi + Psi) / (12.0*Omega_nu0);
    for (int ell = 3; ell < n_ell_Nu; ell++)
      *(Nu+ell) = c*k / ((2.0*ell + 1.0)*Hp) * *(Nu+ell-1);
  }

  return y_ini;
}

//====================================================
// Set the initial conditions after tight coupling
//====================================================
Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const
{
  // The vector we are going to fill
  Vector y_ini(n_ell_tot);

  // References to the tight coupling quantities
  const double &Phi_tc       =  y_tc[idx_Phi];
  const double &delta_CDM_tc =  y_tc[idx_delta_CDM];
  const double &delta_b_tc   =  y_tc[idx_delta_b];
  const double &v_CDM_tc     =  y_tc[idx_v_CDM];
  const double &v_b_tc       =  y_tc[idx_v_b];
  const double *Theta_tc     = &y_tc[idx_start_Theta];
  const double *Nu_tc        = &y_tc[idx_start_Nu];

  // References to the quantities we are going to set
  double &Phi                =  y_ini[idx_Phi];
  double &delta_CDM          =  y_ini[idx_delta_CDM];
  double &delta_b            =  y_ini[idx_delta_b];
  double &v_CDM              =  y_ini[idx_v_CDM];
  double &v_b                =  y_ini[idx_v_b];
  double *Theta              = &y_ini[idx_start_Theta];
  double *Theta_P            = &y_ini[idx_start_Theta_P];
  double *Nu                 = &y_ini[idx_start_Nu];

  // Physical constants
  double c = Constants.c;

  // Cosmological and recombination quantities
  double Hp      = cosmo->Hp_of_x(x);
  double dtaudx  = rec->dtaudx_of_x(x);

  // Scalar quantities (Gravitational potental, baryons and CDM)
  Phi            = Phi_tc;
  delta_CDM      = delta_CDM_tc;
  delta_b        = delta_b_tc;
  v_CDM          = v_CDM_tc;
  v_b            = v_b_tc; 

  // Photon temperature perturbations 
  *Theta         = *Theta_tc;
  *(Theta+1)     = *(Theta_tc+1);
  if (polarization)
    *(Theta+2)   = - 8.0*c*k / (15.0*Hp*dtaudx) * *(Theta+1);
  else
    *(Theta+2)   = - 20.0*c*k / (45.0*Hp*dtaudx) * *(Theta+1);
  for (int ell = 3; ell < n_ell_Theta; ell++) 
    *(Theta+ell) = - ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * *(Theta+ell-1);

  // Photon polarization perturbations
  if (polarization) {
    *Theta_P     = 5.0/4.0 * *(Theta+2);
    *(Theta_P+1) = - c*k / (4.0*Hp*dtaudx) * *(Theta+2);
    *(Theta_P+2) = 1.0/4.0 * *(Theta+2);
    for (int ell = 3; ell < n_ell_Theta_P; ell++)
      *(Theta_P+ell) = - ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * *(Theta_P+ell-1);
  }

  // Neutrino perturbations
  if (neutrinos) {
    for (int ell = 0; ell < n_ell_Nu; ell++)
      *(Nu+ell)  = *(Nu_tc+ell);
  }

  return y_ini;
}

//====================================================
// Compute the index when tight coupling ends
//====================================================
int Perturbations::get_tight_coupling_index(
    const Vector &x_array, 
    const double k) const
{
  double c      = Constants.c;
  double Hp     = cosmo->Hp_of_x(x_array[0]);
  double dtaudx = rec->dtaudx_of_x(x_array[0]);

  int i = 1;
  while (fabs(dtaudx) > 10.0 * std::min(1.0, c*k/Hp) && x_array[i] < -8.3) {
    Hp     = cosmo->Hp_of_x(x_array[i]);
    dtaudx = rec->dtaudx_of_x(x_array[i]);
    i++;
  }

  return i;
}


//====================================================
// The right hand side of the tight coupling ODE
//====================================================
int Perturbations::rhs_tight_coupling_ode(
    double x, 
    double k, 
    const double *y, 
    double *dydx)
{ 
  // The different quantities in the y array
  const double &Phi         =  y[idx_Phi];
  const double &delta_CDM   =  y[idx_delta_CDM];
  const double &delta_b     =  y[idx_delta_b];
  const double &v_CDM       =  y[idx_v_CDM];
  const double &v_b         =  y[idx_v_b];
  const double *Theta       = &y[idx_start_Theta];
  const double *Theta_P     = &y[idx_start_Theta_P];
  const double *Nu          = &y[idx_start_Nu];

  // References to the quantities we are going to set in the dydx array
  double &dPhidx            =  dydx[idx_Phi];
  double &ddelta_CDMdx      =  dydx[idx_delta_CDM];
  double &ddelta_bdx        =  dydx[idx_delta_b];
  double &dv_CDMdx          =  dydx[idx_v_CDM];
  double &dv_bdx            =  dydx[idx_v_b];
  double *dThetadx          = &dydx[idx_start_Theta];
  double *dTheta_Pdx        = &dydx[idx_start_Theta_P];
  double *dNudx             = &dydx[idx_start_Nu];

  // Physical constants
  double c = Constants.c;

  // Cosmological parameters and variables
  double Omega_CDM0   = cosmo->get_Omega_CDM();
  double Omega_b0     = cosmo->get_Omega_b();
  double Omega_gamma0 = cosmo->get_Omega_gamma();
  double Omega_nu0    = cosmo->get_Omega_nu();
  double R            = 4.0*Omega_gamma0/(3.0*Omega_b0*exp(x));
  double H_0          = cosmo->get_H_0();
  double Hp           = cosmo->Hp_of_x(x);
  double dHpdx        = cosmo->dHpdx_of_x(x);
  double eta          = cosmo->eta_of_x(x);

  // Recombination variables
  double dtaudx       = rec->dtaudx_of_x(x);
  double ddtauddx     = rec->ddtauddx_of_x(x);

  // Compute Theta_2 from initial conditions
  double Theta_2;
  if (polarization)
    Theta_2     = - 8.0*c*k / (15.0*Hp*dtaudx) * *(Theta+1);
  else
    Theta_2     = - 20.0*c*k / (45.0*Hp*dtaudx) * *(Theta+1);

  // Compute Psi from other quantities
  double Psi    = - Phi;
  if (neutrinos)
    Psi        -= 12.0*pow(H_0/(c*k), 2.0) * exp(-2.0*x) * (Omega_gamma0 * Theta_2 + Omega_nu0 * *(Nu+2));
  else
    Psi        -= 12.0*pow(H_0/(c*k), 2.0) * exp(-2.0*x) * Omega_gamma0 * Theta_2;

  // Tight-coupling quantities
  dPhidx        = Psi - pow(c*k/Hp, 2.0)/3.0 * Phi;
  if (neutrinos)
    dPhidx     += pow(H_0/Hp, 2.0)/2.0 * ((Omega_CDM0*delta_CDM + Omega_b0*delta_b)*exp(-x) + 4.0*(Omega_gamma0*(*Theta) + Omega_nu0*(*Nu))*exp(-2.0*x));
  else
    dPhidx     += pow(H_0/Hp, 2.0)/2.0 * ((Omega_CDM0*delta_CDM + Omega_b0*delta_b)*exp(-x) + 4.0*Omega_gamma0*(*Theta)*exp(-2.0*x));
  *dThetadx     = - c*k/Hp * *(Theta+1) - dPhidx;
  double q      = (- ((1.0 - R)*dtaudx + (1.0 + R)*ddtauddx) * (3.0*(*(Theta+1)) + v_b)
                   - c*k/Hp * Psi
                   + (1.0 - dHpdx/Hp) * c*k/Hp * (- *Theta + 2.0*Theta_2)
                   - c*k/Hp * *dThetadx) /
                  ((1.0 + R)*dtaudx + dHpdx/Hp - 1.0);
  dv_bdx        = 1.0/(1.0 + R) *
                  (- v_b - c*k/Hp * Psi 
                   + R * (q + c*k/Hp * (- *Theta + 2.0*Theta_2 - Psi)));
  *(dThetadx+1) = (q - dv_bdx)/3.0;

  // Remaining scalar quantities
  ddelta_CDMdx  = c*k/Hp * v_CDM - 3.0*dPhidx;
  ddelta_bdx    = c*k/Hp * v_b - 3.0*dPhidx;
  dv_CDMdx      = - v_CDM - c*k/Hp * Psi;

  // Remaining photon temperature perturbations (set by IC's)
  for (int ell = 2; ell < n_ell_Theta; ell++)
    *(dThetadx+ell) = 0.0;

  // Photon polarization perturbations (set by IC's)
  if (polarization) {
    for (int ell = 0; ell < n_ell_Theta_P; ell++) 
      *(dTheta_Pdx+ell) = 0.0;
  }

  // Neutrino perturbations
  if (neutrinos) {
    *dNudx     = - c*k/Hp * *(Nu+1) - dPhidx;
    *(dNudx+1) = c*k/(3.0*Hp) * (*Nu - 2.0*(*(Nu+2)) + Psi);
    for (int ell = 2; ell < n_ell_Nu-1; ell++)
      *(dNudx+ell) = c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Nu+ell-1)) - (ell+1.0)*(*(Nu+ell+1)));
    *(dNudx+n_ell_Nu-1) = c*k/Hp * (*(Nu+n_ell_Nu-2))
                          - c*(n_ell_Nu)/(Hp*eta) * (*(Nu+n_ell_Nu-1));
  }

  return GSL_SUCCESS;
}


//====================================================
// The right hand side of the full ODE
//====================================================
int Perturbations::rhs_full_ode(
    double x, 
    double k, 
    const double *y, 
    double *dydx) 
{
  // The different quantities in the y array
  const double &Phi         =  y[idx_Phi];
  const double &delta_CDM   =  y[idx_delta_CDM];
  const double &delta_b     =  y[idx_delta_b];
  const double &v_CDM       =  y[idx_v_CDM];
  const double &v_b         =  y[idx_v_b];
  const double *Theta       = &y[idx_start_Theta];
  const double *Theta_P     = &y[idx_start_Theta_P];
  const double *Nu          = &y[idx_start_Nu];

  // References to the quantities we are going to set in the dydx array
  double &dPhidx            =  dydx[idx_Phi];
  double &ddelta_CDMdx      =  dydx[idx_delta_CDM];
  double &ddelta_bdx        =  dydx[idx_delta_b];
  double &dv_CDMdx          =  dydx[idx_v_CDM];
  double &dv_bdx            =  dydx[idx_v_b];
  double *dThetadx          = &dydx[idx_start_Theta];
  double *dTheta_Pdx        = &dydx[idx_start_Theta_P];
  double *dNudx             = &dydx[idx_start_Nu];

  // Physical constants
  double c = Constants.c;

  // Cosmological parameters and variables
  double Omega_CDM0   = cosmo->get_Omega_CDM();
  double Omega_b0     = cosmo->get_Omega_b();
  double Omega_gamma0 = cosmo->get_Omega_gamma();
  double Omega_nu0    = cosmo->get_Omega_nu();
  double R            = 4.0*Omega_gamma0/(3.0*Omega_b0*exp(x));
  double H_0          = cosmo->get_H_0();
  double Hp           = cosmo->Hp_of_x(x);
  double dHpdx        = cosmo->dHpdx_of_x(x);
  double eta          = cosmo->eta_of_x(x);

  // Recombination variables
  double dtaudx       = rec->dtaudx_of_x(x);

  // Compute Psi and Pi from other quantities
  double Psi = - Phi;
  double Pi  = *(Theta+2);
  if (neutrinos)
    Psi     -= 12.0*pow(H_0/(c*k), 2.0) * exp(-2.0*x) * (Omega_gamma0 * *(Theta+2) + Omega_nu0 * *(Nu+2));
  else
    Psi     -= 12.0*pow(H_0/(c*k), 2.0) * exp(-2.0*x) * Omega_gamma0 * *(Theta+2);
  if (polarization)
    Pi      += *Theta_P + *(Theta_P+2); 

  // Scalar quantities (Gravitational potental, baryons and CDM)
  dPhidx        = Psi - pow(c*k/Hp, 2.0)/3.0 * Phi;
  if (neutrinos)
    dPhidx     += pow(H_0/Hp, 2.0)/2.0 * ((Omega_CDM0*delta_CDM + Omega_b0*delta_b)*exp(-x) + 4.0*(Omega_gamma0*(*Theta) + Omega_nu0*(*Nu))*exp(-2.0*x));
  else
    dPhidx     += pow(H_0/Hp, 2.0)/2.0 * ((Omega_CDM0*delta_CDM + Omega_b0*delta_b)*exp(-x) + 4.0*Omega_gamma0*(*Theta)*exp(-2.0*x));
  ddelta_CDMdx  = c*k/Hp * v_CDM - 3.0*dPhidx;
  ddelta_bdx    = c*k/Hp * v_b - 3.0*dPhidx;
  dv_CDMdx      = - v_CDM - c*k/Hp * Psi;
  dv_bdx        = - v_b - c*k/Hp * Psi + dtaudx*R*(3.0*(*(Theta+1)) + v_b);

  // Photon temperature perturbations
  *dThetadx     = - c*k/Hp * *(Theta+1) - dPhidx;
  *(dThetadx+1) = c*k/(3.0*Hp) * (*Theta - 2.0*(*(Theta+2)) + Psi) + dtaudx*(*(Theta+1) + v_b/3.0);
  for (int ell = 2; ell < n_ell_Theta-1; ell++) {
    *(dThetadx+ell)    = c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Theta+ell-1)) - (ell+1.0)*(*(Theta+ell+1))) + dtaudx * *(Theta+ell);
    if (ell == 2) 
      *(dThetadx+ell) -= dtaudx*Pi/10.0;
  }
  *(dThetadx+n_ell_Theta-1) = c*k/Hp * (*(Theta+n_ell_Theta-2))
                              - c*n_ell_Theta/(Hp*eta) * (*(Theta+n_ell_Theta-1))
                              + dtaudx * *(Theta+n_ell_Theta-1);
  
  // Photon polarization perturbations 
  if (polarization) {
    *dTheta_Pdx = - c*k/Hp * *(Theta_P+1) + dtaudx*(*Theta_P - Pi/2.0);
    for (int ell = 1; ell < n_ell_Theta_P-1; ell++) {
      *(dTheta_Pdx+ell)    = c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Theta_P+ell-1)) - (ell+1.0)*(*(Theta_P+ell+1))) + dtaudx * *(Theta_P+ell);
      if (ell == 2) 
        *(dTheta_Pdx+ell) -= dtaudx*Pi/10.0;
    }
    *(dTheta_Pdx+n_ell_Theta_P-1) = c*k/Hp * (*(Theta_P+n_ell_Theta_P-2))
                                    - c*n_ell_Theta_P/(Hp*eta) * (*(Theta_P+n_ell_Theta_P-1))
                                    + dtaudx * *(Theta_P+n_ell_Theta_P-1);
  }

  // Neutrino perturbations
  if (neutrinos) {
    *dNudx     = - c*k/Hp * *(Nu+1) - dPhidx;
    *(dNudx+1) = c*k/(3.0*Hp) * (*Nu - 2.0*(*(Nu+2)) + Psi);
    for (int ell = 2; ell < n_ell_Nu-1; ell++)
      *(dNudx+ell) = c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Nu+ell-1)) - (ell+1.0)*(*(Nu+ell+1)));
    *(dNudx+n_ell_Nu-1) = c*k/Hp * (*(Nu+n_ell_Nu-2))
                          - c*n_ell_Nu/(Hp*eta) * (*(Nu+n_ell_Nu-1));
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x, k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x, k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x, k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x, k);
}
double Perturbations::get_ddPsiddx(const double x, const double k) const{
  return Psi_spline.deriv_xx(x, k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x, k);
}
double Perturbations::get_dPidx(const double x, const double k) const{
  return Pi_spline.deriv_x(x, k);
}
double Perturbations::get_ddPiddx(const double x, const double k) const{
  return Pi_spline.deriv_xx(x, k);
}
double Perturbations::get_delta_CDM(const double x, const double k) const{
  return delta_CDM_spline(x, k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x, k);
}
double Perturbations::get_v_CDM(const double x, const double k) const{
  return v_CDM_spline(x, k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x, k);
}
double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x, k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x, k);
}
double Perturbations::get_Theta_P(const double x, const double k, const int ell) const{
  return Theta_P_spline[ell](x, k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x, k);
}

double Perturbations::get_source_T(const double x, const double k) const{
  return source_T_spline(x, k);
}
double Perturbations::get_source_E(const double x, const double k) const{
  return source_E_spline(x, k);
}
double Perturbations::get_source_nu(const double x, const double k) const{
  return source_nu_spline(x, k);
}
double Perturbations::get_source_Psi(const double x, const double k) const{
  return source_Psi_spline(x, k);
}

std::pair<double,double> Perturbations::get_k_min_max() const{
  return std::pair<double,double>(k_min, k_max);
}

bool Perturbations::get_polarization_bool() const{
  return polarization;
}
bool Perturbations::get_neutrinos_bool() const{
  return neutrinos;
}
bool Perturbations::get_lensing_bool() const{
  return lensing;
}
double Perturbations::get_neutrino_fraction() const{
  return f_nu;
}

//====================================================
// Print some useful info about the class
//====================================================
void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  if(polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";
  if(lensing)
    std::cout << "We compute the lensing potential source\n";
  else
    std::cout << "We do not compute the lensing potential source\n";
  std::cout << "We solve from k_min = " << k_min*Constants.Mpc << "/Mpc to k_max = " << k_max*Constants.Mpc << "/Mpc\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "idx_Phi:            " << idx_Phi              << "\n";
  std::cout << "idx_delta_CDM:      " << idx_delta_CDM        << "\n";
  std::cout << "idx_delta_b:        " << idx_delta_b          << "\n";
  std::cout << "idx_v_CDM:          " << idx_v_CDM            << "\n";
  std::cout << "idx_v_b:            " << idx_v_b              << "\n";
  std::cout << "idx_start_Theta:    " << idx_start_Theta      << "\n";
  std::cout << "n_ell_Theta:        " << n_ell_Theta          << "\n";
  if(polarization){
    std::cout << "idx_start_Theta_P:  " << idx_start_Theta_P  << "\n";
    std::cout << "n_ell_Theta_P:      " << n_ell_Theta_P      << "\n";
  }
  if(neutrinos){
    std::cout << "idx_start_Nu:       " << idx_start_Nu       << "\n";
    std::cout << "n_ell_Nu:           " << n_ell_Nu           << "\n";
  }
  std::cout << "n_ell_tot:          " << n_ell_tot            << "\n";
}


//====================================================
// Print the tight coupling time for a given mode
//====================================================
void Perturbations::print_tight_coupling_time(const double k) const{
  // Make an x-array for finding the time
  Vector x_array = Utils::linspace(-15.0, -8.0, 1000);

  // Find the correct index
  int idx_tight_end = get_tight_coupling_index(x_array, k);

  // Print result
  std::cout << "\n";
  std::cout << "End of tight-coupling for k = " << k*Constants.Mpc << "/Mpc:\n";
  std::cout << "x_tc = " << x_array[idx_tight_end-1]      << "\n";
  std::cout << "a_tc = " << exp(x_array[idx_tight_end-1]) << "\n";
}   


//====================================================
// Print the time of horizon entry for a given mode
//====================================================
void Perturbations::print_horizon_entry_time(const double k) const{
  // Make an x-array for finding the time
  Vector x_array = Utils::linspace(-20.0, 0.0, 10000);

  // Find x_enter
  double Hp;
  double x_enter;
  for (int i = 0; i < 10000; i++) {
    x_enter = x_array[i];
    Hp = cosmo->Hp_of_x(x_enter);
    if (Constants.c*k/Hp >= 1.0)  // TODO: maybe change to k*eta instead
      break;
  }

  // Print result
  std::cout << "\n";
  std::cout << "Time of horizon entry for k = " << k*Constants.Mpc << "/Mpc:\n";
  std::cout << "x_enter = " << x_enter      << "\n";
  std::cout << "a_enter = " << exp(x_enter) << "\n";
}   


//====================================================
// Output some results to file for a given value of k
//====================================================
void Perturbations::output(
    const double x_min, 
    const double x_max, 
    const double k, 
    const std::string filename,
    const bool source) const
{
  std::ofstream fp(filename.c_str());
  const int npts = static_cast<int>(x_max - x_min)*10000 + 1;
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  auto print_data = [&] (const double x) {
    fp << x                        << " ";
    fp << get_Phi(x, k)            << " ";
    fp << get_Psi(x, k)            << " ";
    fp << get_Pi(x, k)             << " ";
    fp << get_delta_CDM(x, k)      << " ";
    fp << get_delta_b(x, k)        << " ";
    fp << get_v_CDM(x, k)          << " ";
    fp << get_v_b(x, k)            << " ";
    fp << get_Theta(x, k, 0)       << " ";
    fp << get_Theta(x, k, 1)       << " ";
    fp << get_Theta(x, k, 2)       << " ";
    if (polarization) {
      fp << get_Theta_P(x, k, 0)   << " ";
      fp << get_Theta_P(x, k, 1)   << " ";
      fp << get_Theta_P(x, k, 2)   << " ";
    }
    if (neutrinos) {
      fp << get_Nu(x, k, 0)        << " ";
      fp << get_Nu(x, k, 1)        << " ";
      fp << get_Nu(x, k, 2)        << " ";
    }
    if (source) {
      double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
      fp << get_source_T(x, k)     << " ";
      fp << get_source_T(x, k) * Utils::j_ell(5,   arg)     << " ";
      fp << get_source_T(x, k) * Utils::j_ell(50,  arg)     << " ";
      fp << get_source_T(x, k) * Utils::j_ell(500, arg)     << " ";
      if (polarization) {
        fp << get_source_E(x, k)   << " ";
        fp << get_source_E(x, k) * Utils::j_ell(5,   arg)   << " ";
        fp << get_source_E(x, k) * Utils::j_ell(50,  arg)   << " ";
        fp << get_source_E(x, k) * Utils::j_ell(500, arg)   << " ";
      }
      if (neutrinos) {
        fp << get_source_nu(x, k)  << " ";
        fp << get_source_nu(x, k) * Utils::j_ell(5,   arg)  << " ";
        fp << get_source_nu(x, k) * Utils::j_ell(50,  arg)  << " ";
        fp << get_source_nu(x, k) * Utils::j_ell(500, arg)  << " ";
      }
      if (lensing) {
        fp << get_source_Psi(x, k) << " ";
        fp << get_source_Psi(x, k) * Utils::j_ell(5,   arg) << " ";
        fp << get_source_Psi(x, k) * Utils::j_ell(50,  arg) << " ";
        fp << get_source_Psi(x, k) * Utils::j_ell(500, arg) << " ";
      }
    }
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

