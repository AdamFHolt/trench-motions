import scipy, math, os, numpy as np
import scipy.special
from scipy.optimize import fsolve
import scipy.integrate

def compute_plate_buoyancy(age,DT,k,alpha,rho0):

	delz=0.5;    # km
	zmin=0;
	zmax=400.0;
	T0=273.;     # K
	T1=T0+DT;    # K

	z=np.arange(0,zmax,delz)*1e3
	T_erf= T1 - DT * scipy.special.erfc(z/(2*np.sqrt(k*age*1e6*365*24*60*60)));
	B_erf= (T1 - T_erf) * rho0 * alpha;
	B_erf_int=scipy.integrate.simps(y=B_erf, x=z, even='avg')

	return B_erf_int

def compute_plate_isotherm(age,DT,k,Tiso):

	T0=273.;     # K
	T1=T0+DT;    # K

	z_iso=scipy.special.erfcinv((T1-Tiso)/DT)*(2*np.sqrt(k*age*1e6*365*24*60*60));
	z_iso=z_iso/1e3;
	
	return z_iso

def compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n,pre,v_trans):

	vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
	g = 9.81;

	if formulation == 1:
		vsp = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h))  ) * (1.0/vel_converter)
	elif formulation == 2:
		vsp = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g)  - \
			(2.0 * vc * slabL * (visc_asthen/h))  ) * (1.0/vel_converter)
	elif formulation == 3:
		vsp = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((1./6.) * (H**2/Rmin) * yield_stress) - \
			(2.0 * vc * slabL * (visc_asthen/h))  ) * (1.0/vel_converter)
	elif formulation == 4:
		vsp = (h/(visc_asthen*Lsp)) * (pre * (slabL * oceanic_buoy * g) - \
			(2.0 * vc * slabL * (visc_asthen/h))  ) * (1.0/vel_converter)
	elif formulation == 5 or formulation == 6:
        power_law_visc_slab = visc_asthen * ((vc)/(2 * h * trans_strain_rate))**((1-n)/n)
        composite_visc_slab = ((1.0/power_law_visc_slab) + (1.0/visc_asthen))**(-1.0)
        if formulation == 5:
			vsp_initial_guess = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - (2.0 * vc * slabL * (visc_asthen/h))) 
			force_term = (h/(visc_asthen*Lsp)) *  ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - (2.0 * vc * slabL * (composite_visc_slab/h)))
		else:
			vsp_initial_guess = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((1./6.) * (H**2/Rmin) * yield_stress) - (2.0 * vc * slabL * (visc_asthen/h))) 
			force_term = (h/(visc_asthen*Lsp)) *  ((slabL * oceanic_buoy * g) - ((1./6.) * (H**2/Rmin) * yield_stress) - (2.0 * vc * slabL * (composite_visc_slab/h)))
		vsp_min = -100000.0*vc; vsp_max = 100000.0*vc;
		max_num_its = 3000; misfit = 1.0;
		dvsp = 100000.0*vc;
		for i in range(0,max_num_its):
			print "Power law iteration %.0f" % (i+1)
			vsps = np.linspace(vsp_min,vsp_max,50);
			for k in range(0,len(vsps)):
				trial_vsp = vsps[k];
				val = force_term + ((abs(trial_vsp)/(2.0*h*trans_strain_rate))**((n-1)/n))*force_term - trial_vsp;
				if abs(val) < misfit:
					misfit = abs(val);
					actual_vsp = trial_vsp;
			dmisfit = abs(val) - misfit;
			if abs(dmisfit) < abs(0.0001*vc):
				break;
			else:
				dvsp = 0.01 * dvsp;
				vsp_min = actual_vsp - dvsp
				vsp_max = actual_vsp + dvsp
		vsp = actual_vsp * (1.0/vel_converter)
	elif formulation == 7:
		v_trans = v_trans * vel_converter
		visc_asthen_slab = visc_asthen * ((vc/v_trans)**((1-n)/n))
		vsp = (((h/(visc_asthen*Lsp))**n) * (1.0/(v_trans**(n-1))) * (pre * (slabL * oceanic_buoy * g) - (2.0 * vc * slabL * (visc_asthen_slab/h)) )**n) * (1.0/vel_converter)
		

	return vsp


def compute_vsp_withDP(vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,DP_ref,visc_asthen_ref,w_ref,vt_ref,w,slabD,yield_stress):

	vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
	g = 9.81;

	# for derivation, see page ~81 in triangle notepad

	if formulation == 1: # viscous bending, with DP
		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) ) * (1.0/vel_converter)
	elif formulation == 2: # no bending, with DP
		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabL * oceanic_buoy * g) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) ) * (1.0/vel_converter)
	elif formulation == 3: # plastic bending, with DP
		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabL * oceanic_buoy * g) - ((1./6.) * (H**2/Rmin) * yield_stress) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) ) * (1.0/vel_converter)

	return vsp







