#include "Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec,
    bool polarization,
    bool neutrinos,
    bool lensing) : 
  cosmo(cosmo), 
  rec(rec),
  polarization(polarization),
  neutrinos(neutrinos),
  lensing(lensing)
{
  // Adjust number of polarization and/or neutrino multipoles
  n_ell_Thetap     *= polarization;
  n_ell_Nu         *= neutrinos;
  n_ell_tot         = n_scalars + n_ell_Theta + n_ell_Thetap + n_ell_Nu;

  // Compute indices
  idx_start_Thetap  = idx_start_Theta  + n_ell_Theta  - (1 - polarization);
  idx_start_Nu      = idx_start_Thetap + n_ell_Thetap + (neutrinos - polarization);

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
    const double k_start,
    const double k_end,
    const int npts_k)
{
  // Integrate all the perturbation equation and spline the result
  integrate_perturbations(x_start, x_end, npts_x, k_start, k_end, npts_k);

  // Compute source functions and spline the result
  compute_source_functions(x_start, x_end, npts_x, k_start, k_end, npts_k);
}

//====================================================
// Integrate all perturbations and spline the results
//====================================================
void Perturbations::integrate_perturbations(
    const double x_start, 
    const double x_end,
    const int npts_x,
    const double k_start,
    const double k_end,
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
  Vector k_array = Utils::linspace(log10(k_start), log10(k_end), npts_k); // These values are corrected in the loop

  // Loop over all wavenumbers
  for(int ik = 0; ik < npts_k; ik++){

    // Progress bar...
    if( (10*ik) / npts_k != (10*ik+10) / npts_k ) {
      std::cout << (100*ik+100)/npts_k << "% " << std::flush;
      if(ik == npts_k-1) std::cout << std::endl;
    }

    // Compute current value of k
    const double k = pow(10.0, k_array[ik]);

    // Change k so that the k_array has logarithmic spacing
    k_array[ik] = k;

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
    // TODO: maybe bad way to do this? compute the dependent multipoles or keep them at zero?
    Vector y_tc(n_ell_tot);
    for (int i = 0; i < n_ell_tot; i++) {
      auto y_tight_array = tight_ode.get_data_by_component(i);
      for (int ix = 0; ix < idx_tight_end; ix++) 
        y_array[i][ix + npts_x*ik] = y_tight_array[ix];
      y_tc[i] = y_tight_array[idx_tight_end];
    }
    auto y_ini_full = set_ic_after_tight_coupling(y_tc, x_array_full[0], k);

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
        y_array[n_ell_tot+1][ix + npts_x*ik] += y_array[idx_start_Thetap][ix + npts_x*ik]
                                                + y_array[idx_start_Thetap+2][ix + npts_x*ik];
    }

  }
  Utils::EndTiming("integrateperturbation");

  //====================================================
  // Make all splines needed
  //====================================================
  Phi_spline.create(x_array, k_array, y_array[idx_Phi], "Phi_spline");
  Psi_spline.create(x_array, k_array, y_array[n_ell_tot], "Psi_spline");
  Pi_spline.create(x_array, k_array, y_array[n_ell_tot+1], "Pi_spline");
  delta_CDM_spline.create(x_array, k_array, y_array[idx_deltaCDM], "delta_CDM_spline");
  delta_b_spline.create(x_array, k_array, y_array[idx_deltab], "delta_b_spline");
  v_CDM_spline.create(x_array, k_array, y_array[idx_vCDM], "v_CDM_spline");
  v_b_spline.create(x_array, k_array, y_array[idx_vb], "v_b_spline");

  Theta_spline = std::vector<Spline2D>(n_ell_Theta);
  for (int ell = 0; ell < n_ell_Theta; ell++)
    Theta_spline[ell].create(x_array, k_array, y_array[idx_start_Theta+ell], "Theta_" + std::to_string(ell) + "_spline");

  if (polarization)
    Theta_p_spline = std::vector<Spline2D>(n_ell_Thetap);
    for (int ell = 0; ell < n_ell_Thetap; ell++)
      Theta_p_spline[ell].create(x_array, k_array, y_array[idx_start_Thetap+ell], "Theta_p_" + std::to_string(ell) + "_spline");
   
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
    const double k_start,
    const double k_end,
    const int npts_k)
{
  Utils::StartTiming("source");

  // Physical constants
  double c = Constants.c;

  // Set up x-array and the k-array
  Vector x_array = Utils::linspace(x_start, x_end, npts_x);
  Vector k_array = Utils::linspace(log10(k_start), log10(k_end), npts_k); // These values are corrected in the loop

  // Make storage for the source functions
  Vector ST_array(npts_x * npts_k);   // Temperature
  Vector SE_array(npts_x * npts_k);   // Polarization
  Vector SN_array(npts_x * npts_k);   // Neutrinos
  Vector SL_array(npts_x * npts_k);   // Lensing

  // Fetch time-independent quantities
  double Omega_k0 = cosmo->get_Omega_k();
  double x_rec    = rec->get_x_recombination();
  double chi_rec  = cosmo->get_comoving_distance_of_x(x_rec);

  // Compute source functions
  for(int ik = 0; ik < npts_k; ik++){
    // Compute current value of k
    const double k = pow(10.0, k_array[ik]);

    // Change k so that the k_array has logarithmic spacing
    k_array[ik] = k;

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

      // Temperature source TODO: correct to use chain rule? rewrite in a different way?
      ST_array[ix + npts_x*ik] = g_tilde * (Theta_0 + Psi + Pi/4.0)
                                 + exp(-tau) * (dPsidx - dPhidx)
                                 - 1.0/(c*k) * (dHpdx*g_tilde*v_b + Hp*dgdx_tilde*v_b + Hp*g_tilde*dv_bdx)
                                 + 3.0/4.0 * pow(c*k, -2.0) *
                                 (dHpdx * (dHpdx*g_tilde*Pi + Hp*dgdx_tilde*Pi + Hp*g_tilde*dPidx)
                                  + Hp * (ddHpddx*g_tilde*Pi + dHpdx*dgdx_tilde*Pi + dHpdx*g_tilde*dPidx 
                                          + dHpdx*dgdx_tilde*Pi + Hp*ddgddx_tilde*Pi + Hp*dgdx_tilde*dPidx 
                                          + dHpdx*g_tilde*dPidx + Hp*dgdx_tilde*dPidx + Hp*g_tilde*ddPiddx));

      // Polarization source TODO: correct expression? correct to multiply by c?
      if(polarization){
        SE_array[ix + npts_x*ik] = 3.0/4.0 * g_tilde*Pi * pow(c*k*chi, -2.0);
      }

      // Neutrino source
      if (neutrinos) {
        // TODO: correct to multiply by chi? what is delta(eta)?
        double Nu_0 = get_Nu(x, k, 0);
        SN_array[ix + npts_x*ik] = (Nu_0 + Psi)*chi + dPsidx - dPhidx;
      }

      // Lensing source TODO: correct definition of S_k? okay to include here?
      if (lensing) {
        double W = 0.0;
        if (x >= x_rec) {
          if (Omega_k0 == 0)
            W = (chi -  chi_rec) / chi;
          else if (Omega_k0 <= 0) //TODO: double check
            W = sin(chi - chi_rec) / sin(chi);
          else
            W = sinh(chi - chi_rec) / sinh(chi);
        }
        SL_array[ix + npts_x*ik] = - 2.0*c*Psi / (Hp*chi) * W;
      }
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  if (polarization) {
    SE_spline.create(x_array, k_array, SE_array, "Source_Pol_x_k");
  }
  if (neutrinos) {
    SN_spline.create(x_array, k_array, SN_array, "Source_Neu_x_k");
  }
  if (lensing) {
    SL_spline.create(x_array, k_array, SL_array, "Source_Lens_x_k");
  }

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
  double &delta_CDM         =  y_ini[idx_deltaCDM];
  double &delta_b           =  y_ini[idx_deltab];
  double &v_CDM             =  y_ini[idx_vCDM];
  double &v_b               =  y_ini[idx_vb];
  double *Theta             = &y_ini[idx_start_Theta];
  double *Theta_p           = &y_ini[idx_start_Thetap];
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
    for (int ell = 0; ell < n_ell_Thetap; ell++)
      *(Theta_p+ell) = 0.0;
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
  const double &delta_CDM_tc =  y_tc[idx_deltaCDM];
  const double &delta_b_tc   =  y_tc[idx_deltab];
  const double &v_CDM_tc     =  y_tc[idx_vCDM];
  const double &v_b_tc       =  y_tc[idx_vb];
  const double *Theta_tc     = &y_tc[idx_start_Theta];
  const double *Nu_tc        = &y_tc[idx_start_Nu];

  // References to the quantities we are going to set
  double &Phi                =  y_ini[idx_Phi];
  double &delta_CDM          =  y_ini[idx_deltaCDM];
  double &delta_b            =  y_ini[idx_deltab];
  double &v_CDM              =  y_ini[idx_vCDM];
  double &v_b                =  y_ini[idx_vb];
  double *Theta              = &y_ini[idx_start_Theta];
  double *Theta_p            = &y_ini[idx_start_Thetap];
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
    *Theta_p     = 5.0/4.0 * *(Theta+2);
    *(Theta_p+1) = - c*k / (4.0*Hp*dtaudx) * *(Theta+2);
    *(Theta_p+2) = 1.0/4.0 * *(Theta+2);
    for (int ell = 3; ell < n_ell_Thetap; ell++)
      *(Theta_p+ell) = - ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * *(Theta_p+ell-1);
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
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================
int Perturbations::rhs_tight_coupling_ode(
    double x, 
    double k, 
    const double *y, 
    double *dydx)
{ 
  // The different quantities in the y array
  const double &Phi         =  y[idx_Phi];
  const double &delta_CDM   =  y[idx_deltaCDM];
  const double &delta_b     =  y[idx_deltab];
  const double &v_CDM       =  y[idx_vCDM];
  const double &v_b         =  y[idx_vb];
  const double *Theta       = &y[idx_start_Theta];
  const double *Theta_p     = &y[idx_start_Thetap];
  const double *Nu          = &y[idx_start_Nu];

  // References to the quantities we are going to set in the dydx array
  double &dPhidx            =  dydx[idx_Phi];
  double &ddelta_CDMdx      =  dydx[idx_deltaCDM];
  double &ddelta_bdx        =  dydx[idx_deltab];
  double &dv_CDMdx          =  dydx[idx_vCDM];
  double &dv_bdx            =  dydx[idx_vb];
  double *dThetadx          = &dydx[idx_start_Theta];
  double *dTheta_pdx        = &dydx[idx_start_Thetap];
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
    for (int ell = 0; ell < n_ell_Thetap; ell++) 
      *(dTheta_pdx+ell) = 0.0;
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
  const double &delta_CDM   =  y[idx_deltaCDM];
  const double &delta_b     =  y[idx_deltab];
  const double &v_CDM       =  y[idx_vCDM];
  const double &v_b         =  y[idx_vb];
  const double *Theta       = &y[idx_start_Theta];
  const double *Theta_p     = &y[idx_start_Thetap];
  const double *Nu          = &y[idx_start_Nu];

  // References to the quantities we are going to set in the dydx array
  double &dPhidx            =  dydx[idx_Phi];
  double &ddelta_CDMdx      =  dydx[idx_deltaCDM];
  double &ddelta_bdx        =  dydx[idx_deltab];
  double &dv_CDMdx          =  dydx[idx_vCDM];
  double &dv_bdx            =  dydx[idx_vb];
  double *dThetadx          = &dydx[idx_start_Theta];
  double *dTheta_pdx        = &dydx[idx_start_Thetap];
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
    Pi      += *Theta_p + *(Theta_p+2); 

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
    *dTheta_pdx = - c*k/Hp * *(Theta_p+1) + dtaudx*(*Theta_p - Pi/2.0);
    for (int ell = 1; ell < n_ell_Thetap-1; ell++) {
      *(dTheta_pdx+ell)    = c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Theta_p+ell-1)) - (ell+1.0)*(*(Theta_p+ell+1))) + dtaudx * *(Theta_p+ell);
      if (ell == 2) 
        *(dTheta_pdx+ell) -= dtaudx*Pi/10.0;
    }
    *(dThetadx+n_ell_Thetap-1) = c*k/Hp * (*(Theta+n_ell_Thetap-2))
                                - c*n_ell_Thetap/(Hp*eta) * (*(Theta_p+n_ell_Thetap-1))
                                + dtaudx * *(Theta_p+n_ell_Thetap-1);
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
  return Phi_spline(x,k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_ddPsiddx(const double x, const double k) const{
  return Psi_spline.deriv_xx(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_dPidx(const double x, const double k) const{
  return Pi_spline.deriv_x(x,k);
}
double Perturbations::get_ddPiddx(const double x, const double k) const{
  return Pi_spline.deriv_xx(x,k);
}
double Perturbations::get_delta_CDM(const double x, const double k) const{
  return delta_CDM_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_CDM(const double x, const double k) const{
  return v_CDM_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Source_N(const double x, const double k) const{
  return SN_spline(x,k);
}
double Perturbations::get_Source_L(const double x, const double k) const{
  return SL_spline(x,k);
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
  // TODO: maybe change this back
  // std::cout << "x_start:       " << x_start                << "\n";
  // std::cout << "x_end:         " << x_end                  << "\n";
  // std::cout << "n_x:     " << n_x              << "\n";
  // std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  // std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  // std::cout << "n_k:     " << n_k              << "\n";
  if(polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "idx_Phi:            " << idx_Phi              << "\n";
  std::cout << "idx_deltaCDM:       " << idx_deltaCDM         << "\n";
  std::cout << "idx_deltab:         " << idx_deltab           << "\n";
  std::cout << "idx_v_CDM:          " << idx_vCDM             << "\n";
  std::cout << "idx_v_b:            " << idx_vb               << "\n";
  std::cout << "idx_start_Theta:    " << idx_start_Theta      << "\n";
  std::cout << "n_ell_Theta:        " << n_ell_Theta          << "\n";
  if(polarization){
    std::cout << "idx_start_Thetap:   " << idx_start_Thetap   << "\n";
    std::cout << "n_ell_Thetap:       " << n_ell_Thetap       << "\n";
  }
  if(neutrinos){
    std::cout << "idx_start_Nu:       " << idx_start_Nu       << "\n";
    std::cout << "n_ell_Nu:           " << n_ell_Nu           << "\n";
  }
  std::cout << "n_ell_tot:          " << n_ell_tot            << "\n";
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(
    const double x_min, 
    const double x_max, 
    const double k, 
    const std::string filename) const
{
  std::ofstream fp(filename.c_str());
  const int npts = static_cast<int>(x_max - x_min)*10000 + 1;
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                   << " ";
    fp << get_Phi(x, k)       << " ";
    fp << get_Psi(x, k)       << " ";
    fp << get_Pi(x, k)        << " ";
    fp << get_delta_CDM(x, k) << " ";
    fp << get_delta_b(x, k)   << " ";
    fp << get_v_CDM(x, k)     << " ";
    fp << get_v_b(x, k)       << " ";
    fp << get_Theta(x, k, 0)  << " ";
    fp << get_Theta(x, k, 1)  << " ";
    fp << get_Theta(x, k, 2)  << " ";
    // TODO: change this?
    fp << get_Source_T(x, k)  << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

