# PRF2020
This repository contains data and code from Mandel et al. (2020).

Data, analysis, and interpretation by T.L. Mandel and colleagues at the University of California, Merced. Code written by T.L. Mandel (current contact: tracy.mandel@unh.edu)

For full information, see:
* Mandel, T.L., Zhou, D.Z., Waldrop, L., Theillard, M., Kleckner, D., & Khatri, S. (2020). Retention of rising droplets in density stratification. Physical Review Fluids 12(5): 124803. https://link.aps.org/doi/10.1103/PhysRevFluids.5.124803

# CODE #

The code posted here includes:
* demo_droplet_code:  A MATLAB script that shows an example of calling the dropTimescales function with experimental data.
* dropTimescales.m:   A MATLAB function that computes the entrainment and retention times associated with an object crossing through a relatively sharp density    interface.
* intersections.m:    This function is called within dropTimescales.m. The intersections.m code was developed by Douglas M. Schwartz; the original license is included in this repository. His code can be obtained directly from the MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections

  
# DATA #

The data uploaded here includes:
* PRF2020_droplet_data.mat: A .mat file containing the kinematic tracking data as well as postprocessed data used in this paper.

Variables contained in this .mat file include:

drops:  A structure array containing data for each of the 179 experimental cases. The fields are as follows:
 * rho:        Ambient fluid density [g/cm^3]
 * z_centered: Vertical positions at which these density values are defined [cm]
 * h_95:       The transition region thickness, defined as the distance accounting for 95% of the density difference [cm]
 * h_lo:       Vertical position of bottom of transition region [cm]
 * h_hi:       Vertical position of top of transition region [cm]
 * N:          Brunt-Vaisala or buoyancy frequency computed at the center of the transition region [s^-1]
 * d:          Droplet diameter [cm]
 * rho_drop:   Droplet density [g/cm^3]
 * rhoU:       Upper-layer fluid density [g/cm^3]
 * rhoL:       Lower-layer fluid density [g/cm^3]
 * t:          Time vector for kinematic tracking data [s]
 * y:          Vertical position of droplet over time [cm]
 * v:          Vertical velocity of droplet over time [cm/s]
 * vRsq:       R^2 of best-fit lines used to compute velocity (see full description in Mandel et al. 2020)
 * u_upper:    Upper-layer terminal velocity [cm/s]
 * u_lower:    Lower-layer terminal velocity [cm/s]
 * u_min:      Minimum velocity reached (for velocity values with R^2 > 0.9) [cm/s]
 * z_min:      Vertical position of velocity minimum [cm]
 * t_min:      Time between drop's entry to transition region and velocity minimum (based on Srdic-Mitrovic et al. 1999) [sec]
 * y_ll:       Approx. upper threshold of constant terminal velocity in the lower layer [cm]
 * y_ul:       Approx. lower threshold of constant terminal velocity in the upper layer [cm]
 * t_e:        Entrainment time [sec]
 * t_r:        Retention time [sec]
 * d_r:        Retention distance [cm]

Beyond the "drops" structure array, vectors containing the values of various nondimensional numbers for each experimental case are also included:
* Arl:    Archimedes number based on lower-layer properties
* Aru:    Archimedes number based on upper-layer properties
* drhoU:  Normalized ifference between upper- and lower-layer ambient fluid densitites, (rhoU-rhoL)/rhoU
* dU:     Difference between upper- and lower-layer terminal velocities [cm/s]
* dU_norm:  Normalized difference, (U_u - U_l)/U_u
* Fr:     Froude number (based on upper-layer properties)
* hond:   Ratio of transition region thickness (h) to droplet diameter (d)
* Rel:    Reynolds number based on lower-layer properties
* Reu:    Reynolds number based on upper-layer properties
* tau_e:  Nondimensionalized entrainment time (t_e x N)
* tau_r:  Nondimensionalized retention time (t_r x N)

The date these data were processed is saved in the variable "date_processed".
