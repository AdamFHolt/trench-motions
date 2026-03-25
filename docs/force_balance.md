# Force Balance: Trench Motion Formulations

## Notation

| Symbol | Code variable | Description | Units |
|--------|--------------|-------------|-------|
| $v_{sp}$ | `vsp` | Subducting plate velocity / slab pull rate (output of `compute_vsp_withDP`) | cm/yr |
| $v_c$ | `vc` | Convergence velocity (subducting − overriding) | m/s |
| $v_t = v_c - v_{sp}$ | — | Trench velocity (positive = retreat). The code computes `vt_estimate = vc/vel_converter − vsp = vt` (both in cm/yr); a positive `vt_estimate` means retreat. | cm/yr |
| $\eta_A$ | `visc_asthen` | Asthenosphere viscosity | Pa·s |
| $\eta_L$ | `visc_lith` | Lithosphere viscosity (F1 only) | Pa·s |
| $h$ | `h` | Asthenosphere channel thickness | m |
| $H$ | `H` | Lithosphere thickness (thermal isotherm depth) | m |
| $R$ | `Rmin` | Minimum radius of curvature at trench | m |
| $D$ | `slabD` | Slab depth | m |
| $L_p$ | `Lsp` | Plate drag length: length of subducting plate, ridge-to-trench (col 26, `LSP, simplified`) | m |
| $L_s$ | `slabL` | Slab drag length: measured slab length (Lallemand col 8, `L`) | m |
| $B$ | `oceanic_buoy` | Integrated oceanic buoyancy ($\int \rho_0 \alpha \Delta T  \mathrm{erfc}\left(\frac{z}{2\sqrt{\kappa t}}\right) dz$) | kg/m² |
| $F_R$ | `ridge_push` | Ridge push force per unit trench length | N/m |
| $g$ | — | Gravitational acceleration (9.81 m/s²) | m/s² |
| $\sigma_Y$ | `yield_stress` | Lithospheric yield stress (F2 only) | Pa |
| $\Delta P_0$ | `DP_ref` | Reference dynamic pressure (from mantle flow computation) | Pa |
| $w$ | `w` | Along-trench width of subduction zone | m |
| $C_{DP}$ | — | Dimensionless DP scaling coefficient (see below) | — |

### Thermal quantities

$H$ and $B$ are computed from the oceanic plate age via the error-function cooling model:

$$H = \mathrm{erfcinv}\left(\frac{T_1 - T_{iso}}{\Delta T}\right) \cdot 2\sqrt{\kappa  t}$$

$$B = \int_0^\infty \rho_0  \alpha  (T_1 - T(z))  dz \quad \text{where} \quad T(z) = T_1 - \Delta T \cdot \mathrm{erfc}\left(\frac{z}{2\sqrt{\kappa t}}\right)$$

The ridge push force integrates the horizontal pressure gradient from the elevated ridge:

$$F_R = g \rho_0 \alpha \Delta T \left(1 + \frac{2\rho_0 \alpha \Delta T}{\pi(\rho_0 - \rho_w)}\right) \kappa t$$

---

## Physical Overview

The force balance equates **driving forces** (slab pull, ridge push) against **resisting forces** (bending, plate drag, slab drag, mantle back-pressure) per unit trench length. All formulations share the same structure — they differ only in how the bending term is treated.

The unknowns are $v_{sp}$ (and equivalently $v_t = v_c - v_{sp}$). The convergence velocity $v_c$ is an observable input; all other terms are either observed (geometry, age) or swept parameters (viscosities).

### Velocity conventions

- **Plate drag** (Couette flow of asthenosphere below the flat subducting plate) scales with the **absolute plate velocity $v_{sp}$** over the ridge-to-trench length $L_p$.
- **Slab drag** (asthenosphere sheared alongside the inclined slab) also scales with $v_{sp}$, over the measured slab length $L_s$ (from Lallemand 2005).
- **DP** scales with the **trench velocity $v_t = v_c - v_{sp}$** and acts as a **driving force** on $v_{sp}$ when the trench retreats ($v_t > 0$).
- **Bending** scales with the **convergence velocity $v_c$** (viscous) or is velocity-independent (plastic).

### Dynamic pressure parameterisation

Instead of solving the full Stokes problem, dynamic pressure scales linearly from a reference solution computed analytically for standard conditions $(\eta_{A,\text{ref}}, w_\text{ref}, v_{t,\text{ref}})$:

$$\Delta P = \Delta P_0 \cdot \frac{\eta_A}{\eta_{A,\text{ref}}} \cdot \frac{w}{w_\text{ref}} \cdot \frac{v_t}{v_{t,\text{ref}}}$$

The force per unit trench length acting over the slab face (depth $D$) is $F_{DP} = \Delta P \cdot D$. Defining the dimensionless coefficient:

$$C_{DP} = \frac{D \cdot \Delta P_0 \cdot w}{\eta_{A,\text{ref}} \cdot w_\text{ref} \cdot v_{t,\text{ref}}}$$

the dynamic pressure force $F_{DP} = \eta_A  C_{DP}  v_t$ acts in the direction of plate motion, driving $v_{sp}$ when the trench retreats ($v_t > 0$). Reference values: $\eta_{A,\text{ref}} = 3\times10^{20}$ Pa·s, $w_\text{ref} = 4448$ km, $v_{t,\text{ref}} = 5$ cm/yr. The canonical $\Delta P_0 = 23.5$ MPa is the analytically-derived maximum for a free-slip mantle base using Royden and Holt (2020, doi:10.1029/2019GC008771)

---

## Formulation 1 — Viscous Bending

### Force balance

$$\boxed{D B g + F_R + \underbrace{\eta_A C_{DP}  (v_c - v_{sp})}_{\text{dynamic pressure}} = \underbrace{\frac{2}{3}\frac{H^3}{R^3}\eta_L  v_c}_{\text{viscous bending}} + \underbrace{2\eta_A \frac{L_p}{h} v_{sp}}_{\text{plate drag}} + \underbrace{\eta_A \frac{L_s}{h} v_{sp}}_{\text{slab drag}}}$$

- **Plate drag** and **slab drag** both scale with the absolute plate velocity $v_{sp}$.
- **Dynamic pressure** scales with the trench velocity $v_t = v_c - v_{sp}$ and drives $v_{sp}$ when the trench retreats ($v_t > 0$).
- **Bending** scales with convergence $v_c$.

### Solved for $v_{sp}$

Substituting $v_t = v_c - v_{sp}$ into $F_{DP} = \eta_A C_{DP} v_t$, moving $F_{DP}$ to the right and collecting $v_{sp}$ terms:

$$\boxed{v_{sp} = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L  v_c + \eta_A C_{DP} v_c}{\eta_A\left(\dfrac{2L_p + L_s}{h} + C_{DP}\right)}}$$

The denominator $\eta_A((2L_p + L_s)/h + C_{DP})$ is the total resistance to slab motion. The $+\eta_A C_{DP} v_c$ in the numerator is the $v_c$ component of the DP driving term.

### Trench velocity

$$v_t = v_c - v_{sp}$$

---

## Formulation 2 — Plastic Bending

### Force balance

$$\boxed{D B g + F_R + \underbrace{\eta_A C_{DP}  (v_c - v_{sp})}_{\text{dynamic pressure}} = \underbrace{\frac{1}{6}\frac{H^2}{R}\sigma_Y}_{\text{plastic bending}} + \underbrace{2\eta_A \frac{L_p}{h} v_{sp}}_{\text{plate drag}} + \underbrace{\eta_A \frac{L_s}{h} v_{sp}}_{\text{slab drag}}}$$

The only difference from F1 is the bending term. Once the lithosphere yields, the bending moment is set by the yield stress $\sigma_Y$ and geometry alone — it is **independent of velocity** and **independent of $\eta_L$**. This has a qualitatively different scaling: in F1 bending resistance grows with $v_c$, in F2 it is fixed.

### Solved for $v_{sp}$

$$\boxed{v_{sp} = \frac{D B g + F_R - \dfrac{1}{6}\dfrac{H^2}{R}\sigma_Y + \eta_A C_{DP} v_c}{\eta_A\left(\dfrac{2L_p + L_s}{h} + C_{DP}\right)}}$$

$$v_t = v_c - v_{sp}$$

---

## Comparison Summary

| | F1 Viscous | F2 Plastic |
|--|--|--|
| Bending term | $\frac{2}{3}\frac{H^3}{R^3}\eta_L v_c$ | $\frac{1}{6}\frac{H^2}{R}\sigma_Y$ |
| Bending depends on $v_c$? | yes | **no** |
| Channel thickness $h$ | fixed, 200 km | fixed, 200 km |
| Plate drag length | $L_p$ (ridge–trench) | same |
| Slab drag length | $L_s$ (Lallemand `L`) | same |
| Denominator | $\eta_A((2L_p+L_s)/h + C_{DP})$ | same as F1 |
