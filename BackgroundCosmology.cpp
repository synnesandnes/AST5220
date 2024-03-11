#include "BackgroundCosmology.h"

putchar('\n');

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h,               // H0/(100km/s/Mpc)
    double OmegaB,          // Baryon density (today)
    double OmegaCDM,        // CDM density (today)
    double OmegaK,          // Dark energy density (today)
    double Neff,            // Effective number of relativistic species
    double TCMB) :          // Temperature of CMB (today) [K]
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  TCMB(TCMB),
  Neff(Neff)
{
    
    H0 = Constants.H0_over_h*h; // Hubble parameter today [1/s]
    
    // Constants
    
    const double kb    = Constants.k_b;
    const double hbar  = Constants.hbar;
    const double c      = Constants.c;
    const double G      = Constants.G;
    
    
    OmegaR = pow(M_PI, 2) / 15 * pow((kb * TCMB),4) / pow(hbar, 3) / pow(c,5) * 8 * M_PI * G /3 / H0 / H0;
    
    OmegaNu = Neff * (7./8.) * pow(4./11., 4./3.) * OmegaR;
    
    const double K = 0;
    
    OmegaK = - K * pow(Constants.c, 2) / H0 / H0;  // Not said to compute?
    
    OmegaLambda = 1 - (OmegaR + OmegaNu + OmegaB + OmegaCDM + OmegaK);
    
    // Sum of Omegas today
    //OmegaSum = OmegaB + OmegaCDM + OmegaLambda + OmegaR + OmegaK + OmegaNu;
    
    
}
    // Solving the background

    void BackgroundCosmology::solve(){
        Utils::StartTiming("Eta");
        
        double npts = 1000; //  number of points for the splines
        
        // Range of x
        
        Vector x_array = Utils::linspace(x_start, x_end, npts);
        
        
        // ODE for d eta / dx
        ODEFunction detadx = [&](double x, const double *eta, double *detadx){
            
            detadx[0] = Constants.c / Hp_of_x(x); // 1/ H_of_x
            
            return GSL_SUCCESS;
        };
        
        
        
        // Setting the initial condition and solving the ODE
        //Vector eta_initial{Constants.c / Hp_of_x(x)};
        Vector eta_initial{0.0};  //
        ODESolver ode;
        ode.solve(detadx, x_array, eta_initial);
        
        // Getting the results
        auto eta_array = ode.get_data_by_component(0);
        
        // Creating the spline
        eta_of_x_spline.create(x_array, eta_array, "eta_of_x");
        Utils::EndTiming("Eta");
        
        
        Utils::StartTiming("t");
        
        // ODE for dt / dx
        
        ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
            
            dtdx[0] = 1 / H_of_x(x);
            
            return GSL_SUCCESS;
        };
        
        
        // Setting the initial condition and solving the ODE
        
        Vector t_initial{1./(2.*H_of_x(x_start))};
        
        ode.solve(dtdx, x_array, t_initial);
        
        auto t_array = ode.get_data_by_component(0);
        
        t_of_x_spline.create(x_array, t_array, "t");
        
        Utils::EndTiming("t");
    }
        
      
// Get methods
        
// Using the Friedmann eq, we return the Hubble param. as a function of x
        
double BackgroundCosmology::H_of_x(double x) const{
    double a = exp(x);
    double res = H0 * sqrt( (OmegaB + OmegaCDM) / pow(a, 3) +  (OmegaR + OmegaNu) / pow(a, 4) + OmegaK/pow(a, 2) + OmegaLambda);
    
    return res;
    
}
        
// Returns the derivative of the Hubble parameter, dH(x)/dx.
// Also calculates the double derivative of Hubble parameter, d^2H(x)/dx^2
        
Doublepair BackgroundCosmology::derivatives_Hdx(double x) const{
    double a = exp(x);
    
    double root = (OmegaB + OmegaCDM) / pow(a, 3) +  (OmegaR + OmegaNu) / pow(a, 4) + OmegaK / pow(a, 2) + OmegaLambda);
    
    double drootdx = -3 * (OmegaB+OmegaCDM) / pow(a, 3) - 4 * (OmegaR+OmegaNu) / pow(a, 4) - 2 * OmegaK / pow(a, 2);
    
    double ddrootdx = 9 * (OmegaB + OmegaCDM) / pow(a, 3) + 16 * (OmegaR + OmegaNu) / pow(a, 4) + 4 * OmegaK / pow(a, 2);
    
    res_d = drootdx
    res_dd = ddrootdx
    return Doublepair(res_d, res_dd)
}

// Returning H prime (Hp = aH, and x = ln a) as a func of x and H_of_x
double BackgroundCosmology::Hp_of_x(double x) const{
    double a = exp(x);
    double res =  a * H_of_x(x);
    
    return res;
    
}
    
// Returning the derivative of H prime, dHp(x)/dx
double BackgroundCosmology::dHpdx_of_x(double x) const{
    double a = exp(x);
    double root = (OmegaB + OmegaCDM) / pow(a, 3) +  (OmegaR + OmegaNu) / pow(a, 4) + OmegaK / pow(a, 2) + OmegaLambda);
    Doublepair der = derivatives_Hdx(x);
    double res = (a * H0 * der.first / (2 * pow(root, 1./2.))) + Hp_of_x(x) ;
    return res;
    
}


// Returning double derivative of Hubble prime, d^2 Hp(x)/dx^2
     
double BackgroundCosmology::ddHpddx_of_x(double x) const{
    Doublepair der = derivatives_Hdx(x);
    double a = exp(x);
    double dHdx = der.first;
    double H = H_of_x(x);
    double ddHdx = der.second;

    return 2 * dHpdx_of_x(x) - Hp_of_x(x) + (a * H0 / 2)*(ddHdx / pow(root, 1./2.) - dHdx / (2 * pow(root, 3./2.)));
}
 
        
        
double BackgroundCosmology::get_OmegaB(double x) const{
  if(x == 0.0) return OmegaB;
    double a = exp(x);
    double Omega = pow(H0, 2) * OmegaB / (pow(a, 3) * pow(Hp_of_x(x), 2));

  return Omega;
}
        
        
double BackgroundCosmology::get_OmegaR(double x) const{
  if(x == 0.0) return OmegaR;
    double a = exp(x);
    double Omega = pow(H0, 2) * OmegaR /(pow(a, 4) * pow(Hp_of_x(x), 2));

  return Omega;
}
        
// Not sure if either: Will always return zero as we are computing without neutrinos

// Or if I should use formula for Omega_Nu as defined in report
double BackgroundCosmology::get_OmegaNu(double x) const{
    double a = exp(x);
    double Omega = pow(H0, 2) * OmegaNu /(pow(a, 4) * pow(Hp_of_x(x), 2));

    return Omega;
}
        
double BackgroundCosmology::get_OmegaCDM(double x) const{
  if(x == 0.0) return OmegaCDM;
    double a = exp(x);
    double Omega = pow(H0, 2) * OmegaCDM / (pow(a, 3) * pow(Hp_of_x(x), 2));;

  return Omega;
}
    
        
double BackgroundCosmology::get_OmegaLambda(double x) const{
  if(x == 0.0) return OmegaLambda;
    double Omega = pow(H0, 2) * OmegaLambda / pow(Hp_of_x(x), 2);

  return Omega;
}



// Not sure if: Will always return zero as we are computing without curvature


double BackgroundCosmology::get_OmegaK(double x) const{

  if(x == 0.0) return OmegaK;
    double a = exp(x);
    double Omega = OmegaK *pow(H0, 2)/(pow(a, 2) * pow(Hp_of_x(x), 2);
    
  return Omega;
}





double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
double res = get_comoving_distance_of_x(x) / exp(x);

return res;
}


double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
    double res = get_comoving_distance_of_x(x) * exp(x);
    return res;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
    double res = eta_of_x(0) - eta_of_x(x);
return res;
}



        
double BackgroundCosmology::get_rho_crit(double x) const{
  if(x == 0.0) return 3*H0*H0/8/M_PI/Constants.G;

  double H_temp = H_of_x(x);
  double res = 3*H_temp*H_temp/8/M_PI/Constants.G;

  return res;
}
        
        
 // Splining

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
// Printing out info about the class
//====================================================
void BackgroundCosmology::info() const{
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB       << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM     << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda  << "\n";
  std::cout << "OmegaK:      " << OmegaK       << "\n";
  std::cout << "OmegaNu:     " << OmegaNu      << "\n";
  std::cout << "OmegaR:      " << OmegaR       << "\n";
  std::cout << "Neff:        " << Neff         << "\n";
  std::cout << "h:           " << h            << "\n";
  std::cout << "TCMB:        " << TCMB         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some data to file
//====================================================
        
void BackgroundCosmology::output(const std::string filename) const{
  // Create x_array to write to file using splines. We chose a smaller interval than the one used to solve the equations (wanted to avoid boundary problems)
  const double x_min = -10;  // -10
  const double x_max = 0.0;  // 0.0
  const int    n_pts = 1000; // 100
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";

    fp << H_of_x(x)          << " ";

    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";

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

