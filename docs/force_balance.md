# Force Balance: Trench Motion Formulations

## Notation

| Symbol | Code variable | Description | Units |
|--------|--------------|-------------|-------|
| $v_{sp}$ | `vsp` | Subducting plate velocity / slab pull rate (output of `compute_vsp_withDP`) | m/s |
| $v_c$ | `vc` | Convergence velocity (subducting − overriding) | m/s |
| $v_t = v_{sp} - v_c$ | — | Physical trench retreat velocity (positive = retreat). **Note:** the code computes `vt_estimate = vc − vsp = −vt` for comparison with observations, so a negative `vt_estimate` means retreat. | m/s |
| $\eta_A$ | `visc_asthen` | Asthenosphere viscosity | Pa·s |
| $\eta_L$ | `visc_lith` | Lithosphere viscosity | Pa·s |
| $h$ | `h` | Asthenosphere channel thickness | m |
| $H$ | `H` | Lithosphere thickness (thermal isotherm depth) | m |
| $R$ | `Rmin` | Minimum radius of curvature at trench | m |
| $D$ | `slabD` | Slab depth | m |
| $L$ | `slabL` | Horizontal plate length (plate drag) | m |
| $\mathcal{L}_{sp}$ | `Lsp` | Slab length in asthenosphere (slab drag) | m |
| $\mathcal{L}_{op}$ | `Lop` | Overriding plate length (F5 only) | m |
| $B$ | `oceanic_buoy` | Integrated oceanic buoyancy ($\int \rho_0 \alpha \Delta T \, \mathrm{erfc}\!\left(\frac{z}{2\sqrt{\kappa t}}\right) dz$) | kg/m² |
| $F_R$ | `ridge_push` | Ridge push force per unit trench length | N/m |
| $g$ | — | Gravitational acceleration (9.81 m/s²) | m/s² |
| $\sigma_Y$ | `yield_stress` | Lithospheric yield stress (F2 only) | Pa |
| $\Delta P_0$ | `DP_ref` | Reference dynamic pressure (from mantle flow computation) | Pa |
| $w$ | `w` | Along-trench width of subduction zone | m |
| $C_{DP}$ | — | Dimensionless DP scaling coefficient (see below) | — |

### Thermal quantities

$H$ and $B$ are computed from the oceanic plate age via the error-function cooling model:

$$H = \mathrm{erfcinv}\!\left(\frac{T_1 - T_{iso}}{\Delta T}\right) \cdot 2\sqrt{\kappa \, t}$$

$$B = \int_0^\infty \rho_0 \, \alpha \, (T_1 - T(z)) \, dz \quad \text{where} \quad T(z) = T_1 - \Delta T \cdot \mathrm{erfc}\!\left(\frac{z}{2\sqrt{\kappa t}}\right)$$

The ridge push force integrates the horizontal pressure gradient from the elevated ridge:

$$F_R = g \, \rho_0 \, \alpha \, \Delta T \left(1 + \frac{2\rho_0 \alpha \Delta T}{\pi(\rho_0 - \rho_w)}\right) \kappa \, t$$

---

## Physical Overview

The force balance equates **driving forces** (slab pull, ridge push) against **resisting forces** (bending, plate drag, asthenosphere drag, mantle back-pressure) per unit trench length. All formulations share the same structure — they differ only in how the bending term and asthenosphere channel thickness $h$ are treated.

The unknowns are $v_{sp}$ (and equivalently $v_t = v_{sp} - v_c$). The convergence velocity $v_c$ is an observable input; all other terms are either observed (geometry, age) or swept parameters (viscosities).

### Dynamic pressure parameterisation

Instead of solving the full Stokes problem, dynamic pressure scales linearly from a reference solution computed analytically for standard conditions $(\eta_{A,\text{ref}},\, w_\text{ref},\, v_{t,\text{ref}})$:

$$\Delta P = \Delta P_0 \cdot \frac{\eta_A}{\eta_{A,\text{ref}}} \cdot \frac{w}{w_\text{ref}} \cdot \frac{v_t}{v_{t,\text{ref}}}$$

The force per unit trench length acting over the slab face (depth $D$) is $F_{DP} = \Delta P \cdot D$. Defining the dimensionless coefficient:

$$C_{DP} = \frac{D \cdot \Delta P_0 \cdot w}{\eta_{A,\text{ref}} \cdot w_\text{ref} \cdot v_{t,\text{ref}}}$$

the DP force becomes $\eta_A \, C_{DP} \, v_t$. Reference values: $\eta_{A,\text{ref}} = 3\times10^{20}$ Pa·s, $w_\text{ref} = 4448$ km, $v_{t,\text{ref}} = 5$ cm/yr. The canonical $\Delta P_0 = 23.5$ MPa is the analytically-derived maximum for a free-slip mantle base.

---

## Formulation 1 — Viscous Bending

### Force balance

$$\boxed{D B g + F_R \;=\; \underbrace{\frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c}_{\text{viscous bending}} + \underbrace{2\eta_A \frac{L}{h} v_c}_{\text{plate drag}} + \underbrace{\eta_A \frac{\mathcal{L}_{sp}}{h} v_{sp}}_{\text{slab drag}} + \underbrace{\eta_A C_{DP} \, v_t}_{\text{DP back-pressure}}}$$

- **Slab drag** scales with the absolute plate velocity $v_{sp}$ (Couette flow of asthenosphere sheared by the moving slab).
- **DP back-pressure** scales with the trench velocity $v_t$ (mantle displaced by trench retreat).
- **Bending** and **plate drag** scale with convergence $v_c$ only.

### Solved for $v_t$

Substituting $v_{sp} = v_t + v_c$ and collecting terms:

$$\boxed{v_t = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L \, v_c - \eta_A \dfrac{2L + \mathcal{L}_{sp}}{h} v_c}{\eta_A\!\left(\dfrac{\mathcal{L}_{sp}}{h} + C_{DP}\right)}}$$

The denominator is the total resistance to trench motion (slab drag + DP). The numerator is the net driving force minus the velocity-independent resistances.

---

## Formulation 2 — Plastic Bending

### Force balance

$$\boxed{D B g + F_R \;=\; \underbrace{\frac{1}{6}\frac{H^2}{R}\sigma_Y}_{\text{plastic bending}} + \underbrace{2\eta_A \frac{L}{h} v_c}_{\text{plate drag}} + \underbrace{\eta_A \frac{\mathcal{L}_{sp}}{h} v_{sp}}_{\text{slab drag}} + \underbrace{\eta_A C_{DP} \, v_t}_{\text{DP back-pressure}}}$$

The only difference from F1 is the bending term. Once the lithosphere yields, the bending moment is set by the yield stress $\sigma_Y$ and geometry alone — it is **independent of velocity** and **independent of $\eta_L$**. This has a qualitatively different scaling: in F1 bending resistance grows with $v_c$, in F2 it is fixed.

### Solved for $v_t$

$$\boxed{v_t = \frac{D B g + F_R - \dfrac{1}{6}\dfrac{H^2}{R}\sigma_Y - \eta_A \dfrac{2L + \mathcal{L}_{sp}}{h} v_c}{\eta_A\!\left(\dfrac{\mathcal{L}_{sp}}{h} + C_{DP}\right)}}$$

---

## Formulation 3 — Viscous Bending, $h \propto \mathcal{L}_{sp}$

Channel thickness scales linearly with slab length, reflecting longer slabs maintaining a thicker shear zone:

$$h_\text{eff} = h_0 + \mathcal{L}_{sp} \cdot \frac{h_\text{ref} - h_0}{\mathcal{L}_{sp,\text{ref}}} \qquad h_0 = 100\text{ km},\; h_\text{ref} = 250\text{ km},\; \mathcal{L}_{sp,\text{ref}} = 5000\text{ km}$$

### Force balance

Identical structure to F1 with $h \to h_\text{eff}$:

$$\boxed{D B g + F_R \;=\; \frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c + 2\eta_A \frac{L}{h_\text{eff}} v_c + \eta_A \frac{\mathcal{L}_{sp}}{h_\text{eff}} v_{sp} + \eta_A C_{DP} \, v_t}$$

Note: $C_{DP}$ is independent of $h$ and is unchanged.

### Solved for $v_t$

$$\boxed{v_t = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L \, v_c - \eta_A \dfrac{2L + \mathcal{L}_{sp}}{h_\text{eff}} v_c}{\eta_A\!\left(\dfrac{\mathcal{L}_{sp}}{h_\text{eff}} + C_{DP}\right)}}$$

---

## Formulation 4 — Viscous Bending, $h \propto v_{sp}$ (iterative)

Channel thickness scales linearly with the slab velocity — faster slabs maintain a thicker shear zone:

$$h_\text{eff}(v_{sp}) = h_0 + v_{sp} \cdot \frac{h_\text{ref} - h_0}{v_{sp,\text{ref}}} \qquad h_0 = 150\text{ km},\; h_\text{ref} = 200\text{ km},\; v_{sp,\text{ref}} = 5\text{ cm/yr}$$

### Force balance

$$\boxed{D B g + F_R \;=\; \frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c + 2\eta_A \frac{L}{h_\text{eff}(v_{sp})} v_c + \eta_A \frac{\mathcal{L}_{sp}}{h_\text{eff}(v_{sp})} v_{sp} + \eta_A C_{DP} \, v_t}$$

### Solved for $v_t$ (implicit)

Because $h_\text{eff}$ depends on $v_{sp}$, the equation is nonlinear and has no closed-form solution. Written implicitly:

$$\boxed{v_t = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L \, v_c - \eta_A \dfrac{2L + \mathcal{L}_{sp}}{h_\text{eff}(v_t + v_c)} v_c}{\eta_A\!\left(\dfrac{\mathcal{L}_{sp}}{h_\text{eff}(v_t + v_c)} + C_{DP}\right)}}$$

Solved by fixed-point iteration, initialising at $v_{sp} = v_{sp,\text{ref}}$ and relaxing with step 0.5 until $|v_{sp}^{(n+1)} - v_{sp}^{(n)}| < 10^{-3}$ cm/yr. If $v_{sp}$ goes negative, iteration stops early.

---

## Formulation 5 — Viscous Bending + Overriding Plate Drag

The asthenosphere under the overriding plate (length $\mathcal{L}_{op}$) is also sheared by the trench motion, adding a drag term coupled to $v_{sp}$.

### Force balance

$$\boxed{D B g + F_R \;=\; \frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c + \underbrace{\left(2\eta_A \frac{L}{h} - \eta_A\frac{\mathcal{L}_{op}}{h}\right)}_{\text{net plate drag}} v_c + \underbrace{\eta_A \frac{\mathcal{L}_{sp} + \mathcal{L}_{op}}{h} v_{sp}}_{\text{slab + OP drag}} + \eta_A C_{DP} \, v_t}$$

The overriding plate drag adds $\eta_A(\mathcal{L}_{op}/h) v_{sp}$ to the resistance and partially offsets the $v_c$-dependent plate drag by $\eta_A(\mathcal{L}_{op}/h) v_c$.

### Solved for $v_t$

$$\boxed{v_t = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L \, v_c - \eta_A \dfrac{2L + \mathcal{L}_{sp}}{h} v_c}{\eta_A\!\left(\dfrac{\mathcal{L}_{sp} + \mathcal{L}_{op}}{h} + C_{DP}\right)}}$$

Crucially, the numerator is **identical to F1**. The overriding plate drag only increases the denominator. A wider overriding plate reduces $v_t$ by increasing the total resistance to trench motion, without changing the net driving force or convergence-dependent resistances.

---

## Comparison Summary

| | F1 Viscous | F2 Plastic | F3 $h(\mathcal{L}_{sp})$ | F4 $h(v_{sp})$ | F5 OP drag |
|--|--|--|--|--|--|
| Bending term | $\frac{2}{3}\frac{H^3}{R^3}\eta_L v_c$ | $\frac{1}{6}\frac{H^2}{R}\sigma_Y$ | same as F1 | same as F1 | same as F1 |
| Bending depends on $v_c$? | yes | **no** | yes | yes | yes |
| Channel thickness $h$ | fixed, 200 km | fixed, 200 km | $h(\mathcal{L}_{sp})$ | $h(v_{sp})$ (iterative) | fixed, 200 km |
| Denominator | $\eta_A(\mathcal{L}_{sp}/h + C_{DP})$ | same as F1 | $\eta_A(\mathcal{L}_{sp}/h_\text{eff} + C_{DP})$ | same structure | $\eta_A((\mathcal{L}_{sp}+\mathcal{L}_{op})/h + C_{DP})$ |
| Numerator vs F1 | — | replace bending | replace $h$ | replace $h$ | **identical** |
| Closed form? | yes | yes | yes | **no** | yes |
