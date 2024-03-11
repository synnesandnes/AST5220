
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
    OmegaSum = OmegaB + OmegaCDM + OmegaLambda + OmegaR + OmegaK + OmegaNu;
    
}
