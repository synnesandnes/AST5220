#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "utils.h"

using namespace std;

using Vector = std::vector<double>;

using Doublepair = std::pair<double,double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // H0/(100km/s/Mpc)
    double OmegaB;                  // Baryon density (today)
    double OmegaCDM;                // CDM density (today)
    double OmegaLambda;             // Dark energy density (today)
    double Neff;                    // Effective number of relativistic species
    double TCMB;                    // Temperature of CMB (today) [K]
   
    // Derived parameters
    double OmegaR;                  // Photon density (today)
    double OmegaNu;                 // Neutrino density (today)
    double OmegaK;                  // Curvature density
    double H0;                      // The Hubble parameter (today)

    // Start and end of x-int
    double x_start = Constants.x_start;
    double x_end   = Constants.x_end;

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"t"};
 
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double OmegaB, 
        double OmegaCDM, 
        double OmegaK,
        double Neff, 
        double TCMB
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve();

    // Output some results to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double t_of_x(double x) const;
    double eta_of_x(double x) const;
    double H_of_x(double x) const;
    Doublepair BackgroundCosmology::derivatives_Hdx(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double get_OmegaB(double x = 0.0) const; 
    double get_OmegaM(double x = 0.0) const; 
    double get_OmegaR(double x = 0.0) const;
    double get_OmegaRtot(double x = 0.0) const; 
    double get_OmegaNu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const; 
    double get_OmegaLambda(double x = 0.0) const; 
    double get_OmegaK(double x = 0.0) const; 
    double get_OmegaMnu(double x = 0.0) const; 
    double get_H0() const;
    double get_h() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;

    // Distance measures
    double get_luminosity_distance_of_x(double x) const;
    double get_comoving_distance_of_x(double x) const;
    double get_rho_crit(double x) const;
    double get_angular_diameter_distance_of_x(double x) const;

};

#endif
