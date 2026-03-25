import scipy, numpy as np
import scipy.special
import scipy.integrate


def integrate_simpson(y, x):
	if hasattr(scipy.integrate, 'simps'):
		return scipy.integrate.simps(y=y, x=x, even='avg')
	return scipy.integrate.simpson(y=y, x=x)


def compute_plate_buoyancy(age,DT,k,alpha,rho0):

	delz=0.5;    # km
	zmin=0;
	zmax=1000.0;
	T0=273.;     # K
	T1=T0+DT;    # K

	z=np.arange(0,zmax,delz)*1e3
	T_erf= T1 - DT * scipy.special.erfc(z/(2*np.sqrt(k*age*1e6*365*24*60*60)));
	B_erf= (T1 - T_erf) * rho0 * alpha;
	B_erf_int=integrate_simpson(y=B_erf, x=z)

	return B_erf_int

def compute_plate_isotherm(age,DT,k,Tiso):

	T0=273.;     # K
	T1=T0+DT;    # K

	z_iso=scipy.special.erfcinv((T1-Tiso)/DT)*(2*np.sqrt(k*age*1e6*365*24*60*60));
	z_iso=z_iso/1e3;
	
	return z_iso


def compute_vsp_withDP(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,slabL_buoy,dip,oceanic_buoy,DP_ref,visc_asthen_ref,w_ref,vt_ref,\
	w,slabD,yield_stress,pre,ridge_push):

	vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s
	g = 9.81;

	#------------- NEWTONIAN FORMULATIONS -------------------
	if formulation == 1: # viscous bending, with DP

		prefactor =  visc_asthen * ( (2.0*Lsp/h) + (slabL/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) + \
			vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	elif formulation == 2: # plastic bending, with DP

		prefactor =  visc_asthen * ( (2.0*Lsp/h) + (slabL/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push  - ((1./6.) * (H**2/Rmin) * yield_stress) + \
			vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	else:
		raise ValueError("unsupported formulation: {}".format(formulation))
		
	vsp = vsp *(1.0/vel_converter)

	return vsp
