#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // Derived parameters
  //=============================================================================
  H0 = Constants.H0_over_h*h;
  OmegaR = 2.*pow(M_PI, 2)/30. * 
           pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 5)) * 
           8.*M_PI*Constants.G/(3.*pow(H0, 2));
  OmegaNu = Neff * 7./8. * pow(4./11., 4./3.) * OmegaR;
  OmegaLambda = 1. - (OmegaB + OmegaCDM) - (OmegaR + OmegaNu) - OmegaK;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
// Correct to compute t as well? Move to own solve?
void BackgroundCosmology::solve(){
  Utils::StartTiming("eta");
  Utils::StartTiming("t");
    
  Vector x_array;

  // The ODE for deta/dx
  // TODO: Correct to change name?
  ODEFunction detadx_func = [&](double x, const double *eta, double *detadx){
    //=============================================================================
    // TODO: Correct? How does this function work?
    //=============================================================================
    double Hp = Hp_of_x(x);
    detadx[0] = Constants.c/Hp;

    return GSL_SUCCESS;
  };

  ODEFunction dtdx_func = [&](double x, const double *t, double *dtdx){
    double H = H_of_x(x);
    dtdx[0] = 1/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  
  //=============================================================================
  // TODO: Reasonable amount of points? Correct way to do this?
  //=============================================================================
  const int npts = 100'000;
  x_array = Utils::linspace(x_start, x_end, npts);

  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  Vector eta_array(npts);
  Vector t_array(npts);
  double Hp = Hp_of_x(x_start);
  double H = H_of_x(x_start);
  eta_array[0] = Constants.c/Hp;
  t_array[0] = 1/(2*H);

  double detadx;
  double dtdx;
  const double dx = (x_end - x_start) / npts;

  // Evolve forward in time
  for (int i = 1; i < npts; i++) {
    detadx_func(x_array[i-1], &eta_array[i-1], &detadx);
    eta_array[i] = eta_array[i-1] + detadx*dx;

    dtdx_func(x_array[i-1], &t_array[i-1], &dtdx);
    t_array[i] = t_array[i-1] + dtdx*dx;
  }

  // TODO: Necessary to give name since it is defined in header?
  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");
  t_of_x_spline.create(x_array, t_array, "t_of_x");
  Utils::EndTiming("eta");
  Utils::EndTiming("t");
}


//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  double H = H0 * 
             sqrt((OmegaB + OmegaCDM)*exp(-3.*x) + 
             (OmegaR + OmegaNu)*exp(-4.*x) + 
             OmegaK*exp(-2.*x) + 
             OmegaLambda);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  double Hp = H0 * 
              sqrt((OmegaB + OmegaCDM)*exp(-x) + 
                   (OmegaR + OmegaNu)*exp(-2.*x) + 
                   OmegaK + 
                   OmegaLambda*exp(2.*x));
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Correct? Okay to call on other functions?
  //=============================================================================
  double Hp = Hp_of_x(x);
  double dHpdx = - pow(H0, 2)/(2.*Hp) * 
                 ((OmegaB + OmegaCDM)*exp(-x) + 
                  2.*(OmegaR + OmegaNu)*exp(-2.*x) - 
                  2.*OmegaLambda*exp(2.*x));

  // // TODO: ALTERNATIVE
  // double dHpdx = H0/2 * 
  //                (-(OmegaB + OmegaCDM)*exp(-x) - 
  //                 2.*(OmegaR + OmegaNu)*exp(-2.*x) - 
  //                 2.*OmegaLambda*exp(2.*x)) /
  //                sqrt((OmegaB + OmegaCDM)*exp(-x) + 
  //                     (OmegaR + OmegaNu)*exp(-2.*x) + 
  //                     OmegaK + 
  //                     OmegaLambda*exp(2.*x));

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Correct? Okay to call on other functions?
  //=============================================================================
  double Hp = Hp_of_x(x);
  double dHpdx = dHpdx_of_x(x);
  double ddHpddx = pow(H0, 2)/(2.*Hp) *
                   ((OmegaB + OmegaCDM)*exp(-x) + 
                    4.*(OmegaR + OmegaNu)*exp(-2.*x) + 
                    4.*OmegaLambda*exp(2.*x) -
                    dHpdx/(2*Hp));

  // // TODO: ALTERNATIVE
  // double denominator = (OmegaB + OmegaCDM)*exp(-x) + 
  //                      (OmegaR + OmegaNu)*exp(-2.*x) + 
  //                      OmegaK + 
  //                      OmegaLambda*exp(2.*x);
  // double ddHpddx = H0/2 * 
  //                (((OmegaB + OmegaCDM)*exp(-x) + 
  //                  4.*(OmegaR + OmegaNu)*exp(-2.*x) + 
  //                  4.*OmegaLambda*exp(2.*x)) /
  //                 sqrt(denominator) +
  //                 1/2 *
  //                 ((OmegaB + OmegaCDM)*exp(-x) + 
  //                  2.*(OmegaR + OmegaNu)*exp(-2.*x) - 
  //                  2.*OmegaLambda*exp(2.*x)) / 
  //                 pow(denominator, 3./2.));

  return ddHpddx;
}

//=============================================================================
// TODO: Correct to call on H_of_x on these?
//=============================================================================
double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  double H = H_of_x(x);
  double OmegaB_of_x = OmegaB / (exp(3.*x) * pow(H/H0, 2));
  return OmegaB_of_x;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  double H = H_of_x(x);
  double OmegaR_of_x = OmegaR / (exp(4.*x) * pow(H/H0, 2));
  return OmegaR_of_x;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  double H = H_of_x(x);
  double OmegaNu_of_x = OmegaNu / (exp(4.*x) * pow(H/H0, 2));
  return OmegaNu_of_x;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  double H = H_of_x(x);
  double OmegaCDM_of_x = OmegaCDM / (exp(3.*x) * pow(H/H0, 2));
  return OmegaCDM_of_x;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  double H = H_of_x(x);
  double OmegaLambda_of_x = OmegaLambda / pow(H/H0, 2);
  return OmegaLambda_of_x;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;
  double H = H_of_x(x);
  double OmegaK_of_x = OmegaK / (exp(2.*x) * pow(H/H0, 2));
  return OmegaK_of_x;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  double eta0 = eta_of_x(x_start);
  double eta = eta_of_x(x);
  double chi = eta0 - eta;
  return chi;
}

double BackgroundCosmology::get_r_of_x(double x) const{
  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  double chi = get_comoving_distance_of_x(x);
  double r;
  if (OmegaK == 0.) {
    r = chi;
  }
  else if (OmegaK < 0) {
    r = chi *
        sin(sqrt(abs(OmegaK))*H0*chi/Constants.c) /
        sqrt(abs(OmegaK))*H0*chi/Constants.c;
  }
  else {
    r = chi *
        sinh(sqrt(abs(OmegaK))*H0*chi/Constants.c) /
        sqrt(abs(OmegaK))*H0*chi/Constants.c;
  }
  return r;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  double a = exp(x);
  double r = get_r_of_x(x);
  double d_L = r/a;
  return d_L;
}

double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Correct?
  //=============================================================================
  double a = exp(x);
  double r = get_r_of_x(x);
  double d_A = a*r;
  return d_A;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

