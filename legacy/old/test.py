
n = 3.5
visc_asthen = 1.e20
vc = 1.3001e-9
h = 200.3e3
trans_strain_rate = 5e-15
slabL = 1106000.
oceanic_buoy = 5880829.40855227
g = 9.81
H = 52003.77279504
Rmin = 318000.
visc_lith = 1e+20 
vc = 1.30010147e-09
slabD = 1200000.
DP_ref = 23500000.0 
w = 4750000.
visc_asthen_ref = 3e+20 
w_ref = 4448000.
vt_ref = 1.58548959919e-09
Lsp=4000000.
vel_converter = 0.01/(365. * 24. * 60. * 60.) ; # cm/yr to m/s 

power_law_visc_slab = visc_asthen * ((vc)/(2 * h * trans_strain_rate))**((1-n)/n)
composite_visc_slab = ((1.0/power_law_visc_slab) + (1.0/visc_asthen))**(-1.0)
print composite_visc_slab
main_term        = ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - (2.0 * vc * slabL * (composite_visc_slab/h)) + vc*((slabD*DP_ref*composite_visc_slab*w)/(visc_asthen_ref*w_ref*vt_ref)) )

v_sp=2.981751e-07
power_law_visc_sp = visc_asthen * (abs(v_sp)/(2 * h * trans_strain_rate))**((1-n)/n)
composite_visc_sp = ((1.0/power_law_visc_sp) + (1.0/visc_asthen))**(-1.0)
print composite_visc_sp
prefactor = composite_visc_sp * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )

print v_sp * (1/vel_converter), vc * (1/vel_converter)
print (1/prefactor) * main_term * (1/vel_converter)


####

prefactor =  visc_asthen * ( (Lsp/h) + ((slabD * DP_ref * w)/(visc_asthen_ref * w_ref * vt_ref)) )
vsp = (1/prefactor) * ((slabL * oceanic_buoy * g) - ((2./3.) * (H**3/Rmin**3) * visc_lith * vc) - \
		(2.0 * vc * slabL * (visc_asthen/h)) + vc*((slabD*DP_ref*visc_asthen*w)/(visc_asthen_ref*w_ref*vt_ref)) )