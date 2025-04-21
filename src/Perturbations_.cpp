#include "Perturbations_.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{
  if (Constants.neutrinos) {
    double OmegaR  = cosmo->get_OmegaR();
    double OmegaNu = cosmo->get_OmegaNu();
    f_nu = OmegaNu / (OmegaR + OmegaNu);
  }
  else
    f_nu = 0.0;
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
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
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

  // Set up arrays for computing the splines TODO: correct?
  Vector2D y_array(Constants.n_ell_tot_full);
  for (int i = 0; i < Constants.n_ell_tot_full; i++)
    y_array[i] = Vector(npts_x * npts_k);

  // Set up x-array and the k-array
  Vector x_array = Utils::linspace(x_start, x_end, npts_x);
  Vector k_array = Utils::linspace(log10(k_start), log10(k_end), npts_k);

  // Loop over all wavenumbers
  for(int ik = 0; ik < npts_k; ik++){

    // Progress bar...
    if( (10*ik) / npts_k != (10*ik+10) / npts_k ) {
      std::cout << (100*ik+100)/npts_k << "% " << std::flush;
      if(ik == npts_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = pow(10.0, k_array[ik]);

    // Set up x-arrays for the two regimes
    Vector tight_x_array;
    Vector full_x_array;

    // Find value to integrate to and divide x_array accordingly
    double x_end_tight = get_tight_coupling_time(x_start, k);
    int idx_end_tight;
    for (int i = 0; i < npts_x; i++) {
      if (x_array[i] >= x_end_tight) {
        // TODO: are all points included here?
        // TODO: maybe do this in function
        idx_end_tight = i-1;
        tight_x_array = Vector(x_array.begin(), x_array.begin() + i-1);
        full_x_array  = Vector(x_array.begin() + i, x_array.end());
        break;
      }
    }

    // The tight coupling ODE system
    ODESolver tight_ode;
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // Integrate from x_start -> x_end_tight
    tight_ode.solve(dydx_tight_coupling, tight_x_array, y_tight_coupling_ini);
    Vector y_tight_coupling(Constants.n_ell_tot_tc);
    for (int i = 0; i < Constants.n_ell_tot_tc; i++)
      y_tight_coupling[i] = tight_ode.get_data_by_component(i)[idx_end_tight];
    // TODO: fill arrays with all values? 

    // The full ODE system
    ODESolver full_ode;
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Set up initial conditions 
    // auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k); //TODO maybe change to this
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, full_x_array[0], k);

    // Integrate from x_end_tight -> x_end
    full_ode.solve(dydx_full, full_x_array, y_full_ini);

    // Store results for splining
    for (int i = 0; i < Constants.n_ell_tot_full; i++) {
      auto y_tight_array = tight_ode.get_data_by_component(i);
      auto y_full_array  = full_ode.get_data_by_component(i);
      for (int ix = 0; ix < npts_x; ix++) {
        if (ix <= idx_end_tight)
          y_array[i][ix + npts_x*ik] = y_tight_array[ix];
        else
          y_array[i][ix + npts_x*ik] = y_full_array[ix];
      }
    }
    // TODO: correct way to do this?? 

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // Set the initial conditions in the tight coupling regime
  //=============================================================================
  double Hp  = cosmo->Hp_of_x(x);
  double Psi = - 1.0 / (3.0/2.0 + 2.0*f_nu/5.0);

  // Scalar quantities (Gravitational potental, baryons and CDM)
  Phi        = - (1.0 + 2.0*f_nu/5.0) * Psi;
  delta_cdm  = - 3.0/2.0 * Psi;
  delta_b    = delta_cdm;
  v_cdm      = - Constants.c*k / (2.0*Hp) * Psi;
  v_b        = v_cdm; 

  // Photon temperature perturbations 
  *Theta     = - 1.0/2.0 * Psi;
  *(Theta+1) = Constants.c*k / (6.0*Hp) * Psi;

  // Neutrino perturbations
  if (neutrinos) {
    *Nu      = *Theta;
    *(Nu+1)  = *(Theta+1);
    *(Nu+2)  = - pow(Constants.c*k/cosmo->get_H0(), 2.0) * exp(2.0*x) * (Phi + Psi) / (12.0*cosmo->get_OmegaNu());
    for (int ell = 3; ell < n_ell_tot_tc; ell++)
      *(Nu+ell) = Constants.c*k / ((2.0*ell + 1.0)*Hp) * *(Nu+ell-1);
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================
Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &Phi             =  y[Constants.ind_Phi];
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double *Theta           = &y[Constants.ind_start_theta];
  double *Theta_p         = &y[Constants.ind_start_thetap];
  double *Nu              = &y[Constants.ind_start_nu];

  //=============================================================================
  // Fill in the initial conditions for the full equation system
  //=============================================================================
  double Hp      = cosmo->Hp_of_x(x);
  double dtaudx  = rec->dtaudx_of_x(x);

  // Scalar quantities (Gravitational potental, baryons and CDM)
  Phi            = Phi_tc;
  delta_cdm      = delta_cdm_tc;
  delta_b        = delta_b_tc;
  v_cdm          = v_cdm_tc;
  v_b            = v_b_tc; 

  // Photon temperature perturbations 
  *Theta         = *Theta_tc;
  *(Theta+1)     = *(Theta_tc+1);
  if (polarization)
    *(Theta+2)   = - 8.0*Constants.c*k / (15.0*Hp*dtaudx) * *(Theta+1);
  else
    *(Theta+2)   = - 20.0*Constants.c*k / (45.0*Hp*dtaudx) * *(Theta+1);
  for (int ell = 3; ell < n_ell_theta; ell++) 
    *(Theta+ell) = - Constants.c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * *(Theta+ell-1);

  // Photon polarization perturbations
  if (polarization) {
    *Theta_p     = 5.0/4.0 * *(Theta+2);
    *(Theta_p+1) = - Constants.c*k / (4.0*Hp*dtaudx) * *(Theta+2);
    *(Theta_p+2) = 1.0/4.0 * *(Theta+2);
    for (int ell = 3; ell < n_ell_thetap; ell++)
      *(Theta_p+ell) = - Constants.c*k / ((2.0*ell + 1.0)*Hp*dtaudx) * *(Theta_p+ell-1);
  }

  // Neutrino perturbations
  if (neutrinos) {
    for (int ell = 0; ell < n_ell_neutrinos_tc; ell++)
      *(Nu+ell) = *(Nu_tc+ell);
    if (n_ell_neutrinos_tc < n_ell_neutrinos) {
      for (int ell = n_ell_neutrinos_tc; ell < n_ell_neutrinos; ell++)
        *(Nu+ell) = Constants.c*k / ((2.0*ell + 1.0)*Hp) * *(Nu+ell-1);
    }
  }

  return y;
}

//====================================================
// The time when tight coupling ends
//====================================================
double Perturbations::get_tight_coupling_time(const double x_start, const double k) const{
  //TODO do differently?
  double x_end_tight = -8.3;
  Vector x_vals = Utils::linspace(x_start, x_end_tight, 1000);

  for (int i = 0; i < 1000; i++) {
    double Hp     = cosmo->Hp_of_x(x_vals[i]);
    double dtaudx = rec->dtaudx_of_x(x_vals[i]);
    if (fabs(dtaudx) <= 10.0 * std::min(1.0, Constants.c*k/Hp)) {
      x_end_tight = x_vals[i];
      break;
    }
  }

  return x_end_tight;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      // const int index = ix + n_x * ik;
      const int index = ix + 1000 * ik; //TODO: maybe change this back

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperature source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  // Cosmological parameters and variables
  double OmegaCDM = cosmo->get_OmegaCDM();
  double OmegaB   = cosmo->get_OmegaB();
  double OmegaR   = cosmo->get_OmegaR();
  double OmegaNu  = cosmo->get_OmegaNu();
  double R        = 4.0*OmegaR/(3.0*OmegaB*exp(x));
  double H0       = cosmo->get_H0();
  double Hp       = cosmo->Hp_of_x(x);
  double dHpdx    = cosmo->dHpdx_of_x(x);

  // Recombination variables
  double dtaudx   = rec->dtaudx_of_x(x);
  double ddtauddx = rec->ddtauddx_of_x(x);

  // Compute Psi from other quantities
  double Psi      = - Phi;
  if (neutrinos)
    Psi -= 12.0*pow(H0/(Constants.c*k), 2.0) * exp(-2.0*x) * (OmegaR * *(Theta+2) + OmegaNu * *(Nu+2));
  else
    Psi -= 12.0*pow(H0/(Constants.c*k), 2.0) * exp(-2.0*x) * OmegaR * *(Theta+2);

  // Tight-coupling quantities
  dPhidx        = Psi - pow(Constants.c*k/Hp, 3.0) * Phi/3.0;
  if (neutrinos)
    dPhidx     += pow(H0/Hp, 2.0)/2.0 * ((OmegaCDM*delta_cdm + OmegaB*delta_b)*exp(-x) + 4.0*(OmegaR*(*Theta) + OmegaNu*(*Nu))*exp(-2.0*x));
  else
    dPhidx     += pow(H0/Hp, 2.0)/2.0 * ((OmegaCDM*delta_cdm + OmegaB*delta_b)*exp(-x) + 4.0*OmegaR*(*Theta)*exp(-2.0*x));
  *dThetadx     = - Constants.c*k/Hp * *(Theta+1) - dPhidx;
  double q      = (- ((1.0 - R)*dtaudx + (1.0 + R)*ddtauddx) * (3.0*(*(Theta+1)) + v_b)
                   - Constants.c*k/Hp * Psi
                   + (1.0 - dHpdx/Hp) * Constants.c*k/Hp * (- *Theta + 2.0*(*(Theta+2)))
                   - Constants.c*k/Hp * *dThetadx) /
                  ((1.0 + R)*dtaudx + dHpdx/Hp - 1.0);
  dv_bdx        = 1.0/(1.0 + R) *
                  (- v_b - Constants.c*k/Hp * Psi 
                   + R * (q + Constants.c*k/Hp * (- *Theta + 2.0*(*(Theta+2)) - Psi)));
  *(dThetadx+1) = (q - dv_bdx)/3.0;

  // Remaining scalar quantities
  ddelta_cdmdx  = Constants.c*k/Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx    = Constants.c*k/Hp * v_b - 3.0*dPhidx;
  dv_cdmdx      = - v_cdm - Constants.c*k/Hp * Psi;

  // Neutrino perturbations
  if (neutrinos) {
    *dNudx     = - Constants.c*k/Hp * *(Nu+1) - dPhidx;
    *(dNudx+1) = Constants.c*k/(3.0*Hp) * (*Nu - 2.0*(*(Nu+2)) + Psi);
    for (int ell = 2; ell < n_ell_neutrinos_tc-1; ell++)
      *(dNudx+ell) = Constants.c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Nu+ell-1)) - (ell+1.0)*(*(Nu+ell+1)));
    *(dNudx+n_ell_neutrinos_tc-1) = Constants.c*k/Hp * (*(Nu+n_ell_neutrinos_tc-2))
                                    - Constants.c*(n_ell_neutrinos_tc)/(Hp*cosmo->eta_of_x(x)) * (*(Nu+n_ell_neutrinos_tc-1));
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================
int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &Phi             =  y[Constants.ind_Phi];
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  double OmegaCDM = cosmo->get_OmegaCDM();
  double OmegaB   = cosmo->get_OmegaB();
  double OmegaR   = cosmo->get_OmegaR();
  double OmegaNu  = cosmo->get_OmegaNu();
  double R        = 4.0*OmegaR/(3.0*OmegaB*exp(x));
  double H0       = cosmo->get_H0();
  double Hp       = cosmo->Hp_of_x(x);
  double dHpdx    = cosmo->dHpdx_of_x(x);

  // Recombination variables
  double dtaudx   = rec->dtaudx_of_x(x);

  // Compute Psi and Pi from other quantities
  double Psi      = - Phi;
  double Pi       = *(Theta+2);
  if (neutrinos)
    Psi -= 12.0*pow(H0/(Constants.c*k), 2.0) * exp(-2.0*x) * (OmegaR * *(Theta+2) + OmegaNu * *(Nu+2));
  else
    Psi -= 12.0*pow(H0/(Constants.c*k), 2.0) * exp(-2.0*x) * OmegaR * *(Theta+2);
  if (polarization)
    Pi  += *Theta_p + *(Theta_p+2); 

  // Scalar quantities (Gravitational potental, baryons and CDM)
  dPhidx       = Psi - pow(Constants.c*k/Hp, 3.0) * Phi/3.0;
  if (neutrinos)
    dPhidx    += pow(H0/Hp, 2.0)/2.0 * ((OmegaCDM*delta_cdm + OmegaB*delta_b)*exp(-x) + 4.0*(OmegaR*(*Theta) + OmegaNu*(*Nu))*exp(-2.0*x));
  else
    dPhidx    += pow(H0/Hp, 2.0)/2.0 * ((OmegaCDM*delta_cdm + OmegaB*delta_b)*exp(-x) + 4.0*OmegaR*(*Theta)*exp(-2.0*x));
  ddelta_cdmdx = Constants.c*k/Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx   = Constants.c*k/Hp * v_b - 3.0*dPhidx;
  dv_cdmdx     = - v_cdm - Constants.c*k/Hp * Psi;
  dv_bdx       = - v_b - Constants.c*k/Hp * Psi + dtaudx*R*(3.0*(*(Theta+1)) + v_b);

  // Photon temperature perturbations
  *dThetadx       = - Constants.c*k/Hp * *(Theta+1) - dPhidx;
  *(dThetadx+1)   = Constants.c*k/(3.0*Hp) * (*Theta - 2.0*(*(Theta+2)) + Psi) + dtaudx*(*(Theta+1) + v_b/3.0);
  for (int ell = 2; ell < n_ell_theta-1; ell++) {
    *(dThetadx+ell)    = Constants.c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Theta+ell-1)) - (ell+1.0)*(*(Theta+ell+1))) + dtaudx * *(Theta+ell);
    if (ell == 2) 
      *(dThetadx+ell) -= dtaudx*Pi/10.0;
  }
  *(dThetadx+n_ell_theta-1) = Constants.c*k/Hp * (*(Theta+n_ell_theta-2))
                              - Constants.c*(n_ell_theta)/(Hp*cosmo->eta_of_x(x)) * (*(Theta+n_ell_theta-1))
                              + dtaudx * *(Theta+n_ell_theta-1);

  // Photon polarization perturbations 
  if (polarization) {
    *dTheta_pdx = - Constants.c*k/Hp * *(Theta_p+1) + dtaudx*(*Theta_p - Pi/2.0);
    for (int ell = 1; ell < n_ell_thetap-1; ell++) {
      *(dTheta_pdx+ell)    = Constants.c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Theta_p+ell-1)) - (ell+1.0)*(*(Theta_p+ell+1))) + dtaudx * *(Theta_p+ell);
      if (ell == 2) 
        *(dTheta_pdx+ell) -= dtaudx*Pi/10.0;
    }
    *(dThetadx+n_ell_thetap-1) = Constants.c*k/Hp * (*(Theta+n_ell_thetap-2))
                                - Constants.c*(n_ell_thetap)/(Hp*cosmo->eta_of_x(x)) * (*(Theta_p+n_ell_thetap-1))
                                + dtaudx * *(Theta_p+n_ell_thetap-1);
  }

  // Neutrino perturbations
  if (neutrinos) {
    *dNudx     = - Constants.c*k/Hp * *(Nu+1) - dPhidx;
    *(dNudx+1) = Constants.c*k/(3.0*Hp) * (*Nu - 2.0*(*(Nu+2)) + Psi);
    for (int ell = 2; ell < n_ell_neutrinos-1; ell++)
      *(dNudx+ell) = Constants.c*k / ((2.0*ell + 1.0)*Hp) * (ell*(*(Nu+ell-1)) - (ell+1.0)*(*(Nu+ell+1)));
    *(dNudx+n_ell_neutrinos-1) = Constants.c*k/Hp * (*(Nu+n_ell_neutrinos-2))
                                 - Constants.c*(n_ell_neutrinos)/(Hp*cosmo->eta_of_x(x)) * (*(Nu+n_ell_neutrinos-1));
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
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
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double x_min, const double x_max, const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  // const int npts = 5000;
  // auto x_array = Utils::linspace(x_start, x_end, npts);
  const int npts = static_cast<int>(x_max - x_min)*100 + 1; // TODO: maybe change
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x, k, 0) << " ";
    fp << get_Theta(x, k, 1) << " ";
    fp << get_Theta(x, k, 2) << " ";
    fp << get_Phi(x, k)      << " ";
    fp << get_Psi(x, k)      << " ";
    fp << get_Pi(x, k)       << " ";
    fp << get_Source_T(x, k) << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

