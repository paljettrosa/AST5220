#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double Omega_b0, 
    double Omega_CDM0, 
    double Omega_k0,
    double N_eff, 
    double T_CMB0) :
  h(h),
  Omega_b0(Omega_b0),
  Omega_CDM0(Omega_CDM0),
  Omega_k0(Omega_k0),
  N_eff(N_eff), 
  T_CMB0(T_CMB0)
{

  // Derived parameters
  H_0           = Constants.H_0_over_h*h;
  Omega_gamma0  = 2.0*pow(M_PI, 2)/30.0 * 
                pow(Constants.k_b*T_CMB0, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 5)) * 
                8.0*M_PI*Constants.G/(3.0*pow(H_0, 2));
  Omega_nu0     = N_eff * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0) * Omega_gamma0;
  Omega_Lambda0 = 1.0 - (Omega_b0 + Omega_CDM0) - (Omega_gamma0 + Omega_nu0) - Omega_k0;

  // Derived analytical expressions for x
  x_rm      = log((Omega_gamma0 + Omega_nu0) / (Omega_b0 + Omega_CDM0));
  x_acc     = (1.0/3.0) * log((Omega_b0 + Omega_CDM0) / (2.0*Omega_Lambda0));
  x_mLambda = (1.0/3.0) * log((Omega_b0 + Omega_CDM0) / Omega_Lambda0);
}

//====================================================
// Solve the background
//====================================================
void BackgroundCosmology::solve(
    const double x_start, 
    const double x_end, 
    const int npts, 
    bool eta, 
    bool t, 
    bool timing)
{
  ODESolver ode;
  
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  if (eta) {
    if (timing) Utils::StartTiming("eta");

    // The ODE for deta/dx
    ODEFunction detadx = [&](double x, const double *eta, double *detadx){
      double Hp = Hp_of_x(x);
      detadx[0] = Constants.c/Hp;

      return GSL_SUCCESS;
    };
    
    double Hp = Hp_of_x(x_start);
    Vector eta_ini {Constants.c/Hp};

    ode.solve(detadx, x_array, eta_ini);
    auto eta_array = ode.get_data_by_component(0);

    eta_of_x_spline.create(x_array, eta_array, "eta_of_x");
    if (timing) Utils::EndTiming("eta");
  };

  if (t) {
    if (timing) Utils::StartTiming("t");

    // The ODE for dt/dx
    ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
      double H = H_of_x(x);
      dtdx[0]  = 1.0/H;

      return GSL_SUCCESS;
    };
    
    double H = H_of_x(x_start);
    Vector t_ini {1.0/(2.0*H)};

    ode.solve(dtdx, x_array, t_ini);
    auto t_array = ode.get_data_by_component(0);

    t_of_x_spline.create(x_array, t_array, "t_of_x");
    if (timing) Utils::EndTiming("t");
  }
}


//====================================================
// Get methods
//====================================================
double BackgroundCosmology::H_of_x(double x) const{
  double H = H_0 * 
             sqrt((Omega_b0 + Omega_CDM0)*exp(-3.0*x) + 
             (Omega_gamma0 + Omega_nu0)*exp(-4.0*x) + 
             Omega_k0*exp(-2.0*x) + 
             Omega_Lambda0);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double Hp = H_0 * 
              sqrt((Omega_b0 + Omega_CDM0)*exp(-x) + 
                   (Omega_gamma0 + Omega_nu0)*exp(-2.0*x) + 
                   Omega_k0 + 
                   Omega_Lambda0*exp(2.0*x));
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  double Hp    = Hp_of_x(x);
  double dHpdx = - pow(H_0, 2)/(2.0*Hp) * 
                 ((Omega_b0 + Omega_CDM0)*exp(-x) + 
                  2.0*(Omega_gamma0 + Omega_nu0)*exp(-2.0*x) - 
                  2.0*Omega_Lambda0*exp(2.0*x));
  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double Hp      = Hp_of_x(x);
  double dHpdx   = dHpdx_of_x(x);
  double ddHpddx = pow(H_0, 2)/Hp *
                   ((Omega_b0 + Omega_CDM0)*exp(-x)/2.0 + 
                    2.0*(Omega_gamma0 + Omega_nu0)*exp(-2.0*x) + 
                    2.0*Omega_Lambda0*exp(2.0*x) -
                    pow(dHpdx, 2)/pow(H_0, 2));
  return ddHpddx;
}

double BackgroundCosmology::get_Omega_b(double x) const{ 
  if(x == 0.0) return Omega_b0;
  double H       = H_of_x(x);
  double Omega_b = Omega_b0 / (exp(3.0*x) * pow(H/H_0, 2));
  return Omega_b;
}

double BackgroundCosmology::get_Omega_gamma(double x) const{ 
  if(x == 0.0) return Omega_gamma0;
  double H           = H_of_x(x);
  double Omega_gamma = Omega_gamma0 / (exp(4.0*x) * pow(H/H_0, 2));
  return Omega_gamma;
}

double BackgroundCosmology::get_Omega_nu(double x) const{ 
  if(x == 0.0) return Omega_nu0;
  double H        = H_of_x(x);
  double Omega_nu = Omega_nu0 / (exp(4.0*x) * pow(H/H_0, 2));
  return Omega_nu;
}

double BackgroundCosmology::get_Omega_CDM(double x) const{ 
  if(x == 0.0) return Omega_CDM0;
  double H         = H_of_x(x);
  double Omega_CDM = Omega_CDM0 / (exp(3.0*x) * pow(H/H_0, 2));
  return Omega_CDM;
}

double BackgroundCosmology::get_Omega_Lambda(double x) const{ 
  if(x == 0.0) return Omega_Lambda0;
  double H            = H_of_x(x);
  double Omega_Lambda = Omega_Lambda0 / pow(H/H_0, 2);
  return Omega_Lambda;
}

double BackgroundCosmology::get_Omega_k(double x) const{ 
  if(x == 0.0) return Omega_k0;
  double H       = H_of_x(x);
  double Omega_k = Omega_k0 / (exp(2.0*x) * pow(H/H_0, 2));
  return Omega_k;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  double eta_0 = eta_of_x(0.0);
  double eta   = eta_of_x(x);
  double chi   = eta_0 - eta;
  return chi;
}

double BackgroundCosmology::get_r_of_x(double x) const{
  double chi = get_comoving_distance_of_x(x);
  double r;
  if (Omega_k0 == 0.0) {
    r = chi;
  }
  else if (Omega_k0 < 0.0) {
    r = chi *
        sin(sqrt(abs(Omega_k0))*H_0*chi/Constants.c) /
        (sqrt(abs(Omega_k0))*H_0*chi/Constants.c);
  }
  else {
    r = chi *
        sinh(sqrt(abs(Omega_k0))*H_0*chi/Constants.c) /
        (sqrt(abs(Omega_k0))*H_0*chi/Constants.c);
  }
  return r;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double a   = exp(x);
  double r   = get_r_of_x(x);
  double d_L = r/a;
  return d_L;
}

double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
  double a   = exp(x);
  double r   = get_r_of_x(x);
  double d_A = a*r;
  return d_A;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::detadx_of_x(double x) const{
  return eta_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H_0() const{ 
  return H_0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_N_eff() const{ 
  return N_eff; 
}

double BackgroundCosmology::get_T_CMB(double x) const{ 
  if(x == 0.0) return T_CMB0;
  return T_CMB0 * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "h:             " << h             << "\n";
  std::cout << "Omega_b0:      " << Omega_b0      << "\n";
  std::cout << "Omega_CDM0:    " << Omega_CDM0    << "\n";
  std::cout << "Omega_Lambda0: " << Omega_Lambda0 << "\n";
  std::cout << "Omega_k0:      " << Omega_k0      << "\n";
  std::cout << "Omega_nu0:     " << Omega_nu0     << "\n";
  std::cout << "Omega_gamma0:  " << Omega_gamma0  << "\n";
  std::cout << "N_eff:         " << N_eff         << "\n";
  std::cout << "T_CMB0:        " << T_CMB0        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Print out important values of t and eta
//====================================================
void BackgroundCosmology::print_times() const{
  std::cout << "\n";
  std::cout << "Radiation-matter equality:\n";
  std::cout << "t:      " << t_of_x(x_rm)/Constants.kyr << " kyr\n";
  std::cout << "eta/c:  " << eta_of_x(x_rm)/Constants.c/Constants.Myr << " Myr\n";
  std::cout << "eta:    " << eta_of_x(x_rm)/Constants.Mpc << " Mpc\n";

  std::cout << "\n";
  std::cout << "Onset of acceleration:\n";
  std::cout << "t:      " << t_of_x(x_acc)/Constants.Gyr << " Gyr\n";
  std::cout << "eta/c:  " << eta_of_x(x_acc)/Constants.c/Constants.Gyr << " Gyr\n";
  std::cout << "eta:    " << eta_of_x(x_acc)/Constants.Gpc << " Gpc\n";

  std::cout << "\n";
  std::cout << "Matter-dark energy equality:\n";
  std::cout << "t:      " << t_of_x(x_mLambda)/Constants.Gyr << " Gyr\n";
  std::cout << "eta/c:  " << eta_of_x(x_mLambda)/Constants.c/Constants.Gyr << " Gyr\n";
  std::cout << "eta:    " << eta_of_x(x_mLambda)/Constants.Gpc << " Gpc\n";

  std::cout << "\n";
  std::cout << "Today:\n";
  std::cout << "t:      " << t_of_x(0.0)/Constants.Gyr << " Gyr\n";
  std::cout << "eta/c:  " << eta_of_x(0.0)/Constants.c/Constants.Gyr << " Gyr\n";
  std::cout << "eta:    " << eta_of_x(0.0)/Constants.Gpc << " Gpc\n\n";
}

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(
    const double x_min, 
    const double x_max, 
    const std::string filename, 
    bool t, 
    bool detadx, 
    bool distances, 
    bool T_CMB) const
{
  const int npts = static_cast<int>(x_max - x_min)*100 + 1; 
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                   << " ";
    fp << eta_of_x(x)         << " ";
    fp << Hp_of_x(x)          << " ";
    fp << dHpdx_of_x(x)       << " ";
    fp << ddHpddx_of_x(x)     << " ";
    fp << get_Omega_b(x)      << " ";
    fp << get_Omega_CDM(x)    << " ";
    fp << get_Omega_Lambda(x) << " ";
    fp << get_Omega_gamma(x)  << " ";
    fp << get_Omega_nu(x)     << " ";
    fp << get_Omega_k(x)      << " ";
    if (t) {
      fp << t_of_x(x)         << " ";
    }
    if (detadx) {
      fp << detadx_of_x(x)    << " ";
    }
    if (distances) {
      fp << get_comoving_distance_of_x(x)         << " ";
      fp << get_luminosity_distance_of_x(x)       << " ";
      fp << get_angular_diameter_distance_of_x(x) << " ";
    }
    if (T_CMB) {
      fp << get_T_CMB(x)      << " ";
    }
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

