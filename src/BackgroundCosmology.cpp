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
  OmegaR = 2.0*pow(M_PI, 2)/30.0 * 
           pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar, 3)*pow(Constants.c, 5)) * 
           8.0*M_PI*Constants.G/(3.0*pow(H0, 2));
  OmegaNu = Neff * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0) * OmegaR;
  OmegaLambda = 1.0 - (OmegaB + OmegaCDM) - (OmegaR + OmegaNu) - OmegaK;

  //=============================================================================
  // Derived analytical expressions for x
  //=============================================================================

  x_rm = log((OmegaR + OmegaNu) / (OmegaB + OmegaCDM));
  x_acc = (1.0/3.0) * log((OmegaB + OmegaCDM) / (2.0*OmegaLambda));
  x_mLambda = (1.0/3.0) * log((OmegaB + OmegaCDM) / OmegaLambda);
}

//====================================================
// Solve the background
//====================================================

void BackgroundCosmology::solve(const double x_start, const double x_end, const int npts, bool eta, bool t, bool timing){
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
      dtdx[0] = 1.0/H;

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
  double H = H0 * 
             sqrt((OmegaB + OmegaCDM)*exp(-3.0*x) + 
             (OmegaR + OmegaNu)*exp(-4.0*x) + 
             OmegaK*exp(-2.0*x) + 
             OmegaLambda);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double Hp = H0 * 
              sqrt((OmegaB + OmegaCDM)*exp(-x) + 
                   (OmegaR + OmegaNu)*exp(-2.0*x) + 
                   OmegaK + 
                   OmegaLambda*exp(2.0*x));
  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  double Hp = Hp_of_x(x);
  double dHpdx = - pow(H0, 2)/(2.0*Hp) * 
                 ((OmegaB + OmegaCDM)*exp(-x) + 
                  2.0*(OmegaR + OmegaNu)*exp(-2.0*x) - 
                  2.0*OmegaLambda*exp(2.0*x));
  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double Hp = Hp_of_x(x);
  double dHpdx = dHpdx_of_x(x);
  double ddHpddx = pow(H0, 2)/Hp *
                   ((OmegaB + OmegaCDM)*exp(-x)/2.0 + 
                    2.0*(OmegaR + OmegaNu)*exp(-2.0*x) + 
                    2.0*OmegaLambda*exp(2.0*x) -
                    pow(dHpdx, 2)/pow(H0, 2));
  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  double H = H_of_x(x);
  double OmegaB_of_x = OmegaB / (exp(3.0*x) * pow(H/H0, 2));
  return OmegaB_of_x;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  double H = H_of_x(x);
  double OmegaR_of_x = OmegaR / (exp(4.0*x) * pow(H/H0, 2));
  return OmegaR_of_x;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;
  double H = H_of_x(x);
  double OmegaNu_of_x = OmegaNu / (exp(4.0*x) * pow(H/H0, 2));
  return OmegaNu_of_x;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  double H = H_of_x(x);
  double OmegaCDM_of_x = OmegaCDM / (exp(3.0*x) * pow(H/H0, 2));
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
  double OmegaK_of_x = OmegaK / (exp(2.0*x) * pow(H/H0, 2));
  return OmegaK_of_x;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  double eta0 = eta_of_x(0.0);
  double eta = eta_of_x(x);
  double chi = eta0 - eta;
  return chi;
}

double BackgroundCosmology::get_r_of_x(double x) const{
  double chi = get_comoving_distance_of_x(x);
  double r;
  if (OmegaK == 0.0) {
    r = chi;
  }
  else if (OmegaK < 0.0) {
    r = chi *
        sin(sqrt(abs(OmegaK))*H0*chi/Constants.c) /
        (sqrt(abs(OmegaK))*H0*chi/Constants.c);
  }
  else {
    r = chi *
        sinh(sqrt(abs(OmegaK))*H0*chi/Constants.c) /
        (sqrt(abs(OmegaK))*H0*chi/Constants.c);
  }
  return r;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  double a = exp(x);
  double r = get_r_of_x(x);
  double d_L = r/a;
  return d_L;
}

double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
  double a = exp(x);
  double r = get_r_of_x(x);
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
void BackgroundCosmology::output(const double x_min, const double x_max, const std::string filename, bool t, bool detadx, bool distances, bool TCMB) const{
  const int npts = static_cast<int>(x_max - x_min)*100 + 1; 
  Vector x_array = Utils::linspace(x_min, x_max, npts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    if (t) {
      fp << t_of_x(x)        << " ";
    }
    if (detadx) {
      fp << detadx_of_x(x)   << " ";
    }
    if (distances) {
      fp << get_comoving_distance_of_x(x)         << " ";
      fp << get_luminosity_distance_of_x(x)       << " ";
      fp << get_angular_diameter_distance_of_x(x) << " ";
    }
    if (TCMB) {
      fp << get_TCMB(x)      << " ";
    }
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

