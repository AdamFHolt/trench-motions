import scipy, math, os, numpy as np
import scipy.special
from scipy.optimize import fsolve
import scipy.integrate

def compute_plate_buoyancy(age,DT,k,alpha,rho0):

	delz=0.5;    # km
	zmin=0;
	zmax=1000.0;
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

def compute_vsp(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,dip,oceanic_buoy,yield_stress,trans_strain_rate,n,pre,external_force_factor,PSP_force_transmitted):

	vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
	g = 9.81;

	if formulation == 1:
		vsp = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + (external_force_factor * PSP_force_transmitted) )
	elif formulation == 2:
		vsp = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((1./6.) * (H**2/Rmin) * yield_stress) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + (external_force_factor * PSP_force_transmitted) )
	elif formulation == 3:
		vsp = (h/(visc_asthen*Lsp)) * (pre * (slabL * oceanic_buoy * g) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + (external_force_factor * PSP_force_transmitted) )
	elif formulation >= 4:
		power_law_visc_slab = visc_asthen * ((vc)/(2 * h * trans_strain_rate))**((1-n)/n)
		composite_visc_slab = ((1.0/power_law_visc_slab) + (1.0/visc_asthen))**(-1.0)

		if formulation == 4:
			force_term        = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - (2.0 * vc * slabL * (composite_visc_slab/h)) + (external_force_factor * PSP_force_transmitted) )
		elif formulation == 5:
			force_term        = (h/(visc_asthen*Lsp)) * ((slabL * oceanic_buoy * g) - ((1./6.) * (H**2/Rmin) * yield_stress) - (2.0 * vc * slabL * (composite_visc_slab/h)) + (external_force_factor * PSP_force_transmitted))
		elif formulation == 6:
			force_term        = (h/(visc_asthen*Lsp)) * (pre*(slabL * oceanic_buoy * g) - (2.0 * vc * slabL * (composite_visc_slab/h)) + (external_force_factor * PSP_force_transmitted))

		vsp_min = -20000.0*vc; vsp_max = 20000.0*vc;
		max_num_its = 30; misfit = 1.0;
		dvsp = 20000.0*vc;
		for i in range(0,max_num_its):
			trial_vsps = np.linspace(vsp_min,vsp_max,50);
			for k in range(0,len(trial_vsps)):
				trial_vsp = trial_vsps[k];
				val = force_term + ((abs(trial_vsp)/(2.0*h*trans_strain_rate))**((n-1)/n))*force_term - trial_vsp;
				if abs(val) < misfit:
					misfit = abs(val);
					vsp = trial_vsp;
			if abs(misfit) < abs(0.0005*vc):
				break;
			else:
				dvsp = 0.1 * dvsp;
				vsp_min = vsp - dvsp
				vsp_max = vsp + dvsp

	vsp = vsp * (1.0/vel_converter)

	return vsp


def compute_vsp_withDP(formulation,vc,h,visc_asthen,visc_lith,H,Lsp,Rmin,slabL,slabL_buoy,dip,oceanic_buoy,DP_ref,visc_asthen_ref,w_ref,vt_ref,\
	w,slabD,yield_stress,n,pre,trans_strain_rate,composite,external_force_factor,PSP_force_transmitted,ridge_push):

	vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 
	g = 9.81;

	# for derivation, see page 81 in triangle notepad for newtonian, page 27 for power-law
	if formulation == 1: # viscous bending, with DP
		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
	elif formulation == 2: # plastic bending, with DP
		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push  - ((1./6.) * (H**2/Rmin) * yield_stress) - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
	elif formulation == 3: # just slab pull (variable pre-factor), with DP
		prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * (pre * (slabD * oceanic_buoy * g)  + ridge_push - \
			(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
	elif formulation == 7: # viscous bending, with h proportional to Lsp
		Lsp_ref = 5000.e3; 
		hsp_ref = 250.e3;
		h0 = 100.e3;
		hsp_eff = h0 +  Lsp*((hsp_ref - h0)/Lsp_ref)
		prefactor =  visc_asthen * ( (Lsp/hsp_eff) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
		vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
			(2.0 * vc * slabL * (visc_asthen/h0)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
	elif formulation == 8: # viscous bending, with h proportional to Vsp

		vsp_ref = 5.0 * vel_converter # 5 cm/yr 
		max_num_its = 10000; tolerance=0.001 * vel_converter
		vsp_init = vsp_ref

		for i in range(0,max_num_its):
			hSP_eff = 200.e3 + (vsp_init/vsp_ref)
			prefactor =  visc_asthen * ( (Lsp/hSP_eff) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
			vsp = (1/prefactor) * ((slabD * oceanic_buoy * g) + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
				(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
			
			if abs(vsp_init - vsp) < tolerance:
				# print "success after %.0f iterations..." % i
				# print "%.6f, %.6f" % (vsp_init/vel_converter,vsp/vel_converter)
				break
			else:
				vsp_init = vsp_init + 0.5*(vsp-vsp_init)
				if i == max_num_its-1:
					print "%.6f, %.6f %.10f" % (vsp_init/vel_converter,vsp/vel_converter,(vsp_init - vsp)/vel_converter)
					print "COMPLETE DISASTER -- no convergence!"
					exit()


	else: # power-law, viscous bending / plastic bending

		power_law_visc_slab = visc_asthen * ((vc)/(2 * h * trans_strain_rate))**((1-n)/n)
		if composite == 1:
			composite_visc_slab = ((1.0/power_law_visc_slab) + (1.0/visc_asthen))**(-1.0)
		else:
			composite_visc_slab = power_law_visc_slab

		if formulation == 4: 	# viscous bending
			main_term        = ((slabD * oceanic_buoy * g)  + ridge_push - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - (2.0 * vc * slabL * (composite_visc_slab/h)) \
				+ vc*((slabD*DP_ref*composite_visc_slab*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
		elif formulation == 5: 	# plastic bending
			main_term        = ((slabD * oceanic_buoy * g)  + ridge_push - ((1./6.) * (H**2/Rmin) * yield_stress)      - (2.0 * vc * slabL * (composite_visc_slab/h)) \
				+ vc*((slabD*DP_ref*composite_visc_slab*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))
		elif formulation == 6: 	# constant slab pull prefactor
			main_term        = ((pre *slabD * oceanic_buoy * g)  + ridge_push - (2.0 * vc * slabL * (composite_visc_slab/h)) \
				+ vc*((slabD*DP_ref*composite_visc_slab*w)/(visc_asthen_ref*w_ref*vt_ref)) + (external_force_factor * PSP_force_transmitted))

		vsp_min = -20000.0*vc; vsp_max = 20000.0*vc;
		max_num_its = 30; misfit = 1.0;
		dvsp = 20000.0*vc;
		for i in range(0,max_num_its):

			trial_vsps = np.linspace(vsp_min,vsp_max,50);
			for k in range(0,len(trial_vsps)):
				trial_vsp = trial_vsps[k];
				power_law_visc_sp = visc_asthen * (abs(trial_vsp)/(2 * h * trans_strain_rate))**((1-n)/n)
				if composite == 1:
					composite_visc_sp = ((1.0/power_law_visc_sp) + (1.0/visc_asthen))**(-1.0)
				else:
					composite_visc_sp = power_law_visc_sp
				prefactor = (composite_visc_sp*(Lsp/h) +  composite_visc_slab*((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
				val = (1.0/prefactor)*main_term - trial_vsp;
				if abs(val) < misfit:
					misfit = abs(val);
					vsp = trial_vsp;
					composite_visc_sp_pref = composite_visc_sp

			if abs(misfit) < abs(0.0005*vc):
				break;
			else:
				dvsp = 0.1 * dvsp;
				vsp_min = vsp - dvsp
				vsp_max = vsp + dvsp

		
	vsp = vsp *(1.0/vel_converter)

	return vsp
