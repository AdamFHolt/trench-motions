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
	w,slabD,yield_stress,pre,external_force_factor,PSP_force_transmitted,ridge_push,Lop):

	vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
	g = 9.81;

	# for derivations: see page 81 in triangle notepad for newtonian

	#------------- NEWTONIAN FORMULATIONS -------------------
	if formulation == 1: # viscous bending, with DP

		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	elif formulation == 2: # plastic bending, with DP

		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push  - ((1./6.) * (H**2/Rmin) * yield_stress) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	elif formulation == 3: # just slab pull (variable pre-factor), with DP

		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * (pre * (slabD * oceanic_buoy * g)  + ridge_push - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	elif formulation == 7: # viscous bending, with h proportional to Lsp

		Lsp_ref = 5000.e3; 
		hsp_ref = 250.e3;
		h0 = 100.e3;
		hsp_eff = h0 +  Lsp*((hsp_ref - h0)/Lsp_ref)
		prefactor =  visc_asthen * ( (Lsp/hsp_eff) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/hsp_eff)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	elif formulation == 8: # viscous bending, with h proportional to Vsp

		vsp_ref = 5.0 * vel_converter # 5 cm/yr 
		hsp_ref = 200.e3; 
		h0 = 150.e3;

		max_num_its = 5000; tolerance=0.001 * vel_converter
		vsp_init = vsp_ref

		for i in range(0,max_num_its):

			hSP_eff = h0 + vsp_init*((hsp_ref - h0)/vsp_ref)

			prefactor =  visc_asthen * ( (Lsp/hSP_eff) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
			vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
				(2.0 * vc * slabL * (visc_asthen/hSP_eff)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))
			
			if abs(vsp_init - vsp) < tolerance:
				break
			else:
				vsp_init = vsp_init + 0.5*(vsp-vsp_init)
				if vsp < 0:
					break

	elif formulation == 9: # viscous bending, with OP drag

		prefactor =  visc_asthen * ( ((Lsp+Lop)/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + (vc * Lop * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)))

	else:
		raise ValueError("unsupported formulation: {}".format(formulation))
		
	vsp = vsp *(1.0/vel_converter)

	return vsp
