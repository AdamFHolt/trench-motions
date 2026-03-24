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
| $L_p$ | `slabL` | Plate drag length: ridge-to-trench distance (Lallemand col 21, `f`) | m |
| $L_s$ | `Lsp` | Slab drag length: along-dip slab length $= D/\sin\delta$ | m |
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

The force balance equates **driving forces** (slab pull, ridge push) against **resisting forces** (bending, plate drag, slab drag, mantle back-pressure) per unit trench length. All formulations share the same structure — they differ only in how the bending term and asthenosphere channel thickness $h$ are treated.

The unknowns are $v_{sp}$ (and equivalently $v_t = v_{sp} - v_c$). The convergence velocity $v_c$ is an observable input; all other terms are either observed (geometry, age) or swept parameters (viscosities).

### Velocity conventions

- **Plate drag** (Couette flow of asthenosphere below the flat subducting plate) scales with the **absolute plate velocity $v_{sp}$** over the ridge-to-trench length $L_p$.
- **Slab drag** (asthenosphere sheared alongside the inclined slab) also scales with $v_{sp}$, over the along-dip slab length $L_s = D/\sin\delta$.
- **DP back-pressure** scales with the **trench retreat velocity $v_t = v_{sp} - v_c$**.
- **Bending** scales with the **convergence velocity $v_c$** (viscous) or is velocity-independent (plastic).

### Dynamic pressure parameterisation

Instead of solving the full Stokes problem, dynamic pressure scales linearly from a reference solution computed analytically for standard conditions $(\eta_{A,\text{ref}},\, w_\text{ref},\, v_{t,\text{ref}})$:

$$\Delta P = \Delta P_0 \cdot \frac{\eta_A}{\eta_{A,\text{ref}}} \cdot \frac{w}{w_\text{ref}} \cdot \frac{v_t}{v_{t,\text{ref}}}$$

The force per unit trench length acting over the slab face (depth $D$) is $F_{DP} = \Delta P \cdot D$. Defining the dimensionless coefficient:

$$C_{DP} = \frac{D \cdot \Delta P_0 \cdot w}{\eta_{A,\text{ref}} \cdot w_\text{ref} \cdot v_{t,\text{ref}}}$$

the DP force becomes $\eta_A \, C_{DP} \, v_t$. Reference values: $\eta_{A,\text{ref}} = 3\times10^{20}$ Pa·s, $w_\text{ref} = 4448$ km, $v_{t,\text{ref}} = 5$ cm/yr. The canonical $\Delta P_0 = 23.5$ MPa is the analytically-derived maximum for a free-slip mantle base.

---

## Formulation 1 — Viscous Bending

### Force balance

$$\boxed{D B g + F_R \;=\; \underbrace{\frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c}_{\text{viscous bending}} + \underbrace{2\eta_A \frac{L_p}{h} v_{sp}}_{\text{plate drag}} + \underbrace{\eta_A \frac{L_s}{h} v_{sp}}_{\text{slab drag}} + \underbrace{\eta_A C_{DP} \, v_t}_{\text{DP back-pressure}}}$$

- **Plate drag** and **slab drag** both scale with the absolute plate velocity $v_{sp}$.
- **DP back-pressure** scales with the trench velocity $v_t$ (mantle displaced by trench retreat).
- **Bending** scales with convergence $v_c$.

### Solved for $v_{sp}$

Substituting $F_{DP} = \eta_A C_{DP}(v_{sp} - v_c)$ and collecting $v_{sp}$ terms:

$$\boxed{v_{sp} = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L \, v_c + \eta_A C_{DP} v_c}{\eta_A\!\left(\dfrac{2L_p + L_s}{h} + C_{DP}\right)}}$$

The denominator $\eta_A((2L + \mathcal{L}_{sp})/h + C_{DP})$ is the total resistance to slab motion. The $+\eta_A C_{DP} v_c$ in the numerator is an algebraic artifact of expanding $F_{DP} = \eta_A C_{DP}(v_{sp} - v_c)$; it is not a physical driving force.

### Trench velocity

$$v_t = v_{sp} - v_c \qquad \text{(positive = trench retreats / slab outpaces convergence)}$$

---

## Formulation 2 — Plastic Bending

### Force balance

$$\boxed{D B g + F_R \;=\; \underbrace{\frac{1}{6}\frac{H^2}{R}\sigma_Y}_{\text{plastic bending}} + \underbrace{2\eta_A \frac{L_p}{h} v_{sp}}_{\text{plate drag}} + \underbrace{\eta_A \frac{L_s}{h} v_{sp}}_{\text{slab drag}} + \underbrace{\eta_A C_{DP} \, v_t}_{\text{DP back-pressure}}}$$

The only difference from F1 is the bending term. Once the lithosphere yields, the bending moment is set by the yield stress $\sigma_Y$ and geometry alone — it is **independent of velocity** and **independent of $\eta_L$**. This has a qualitatively different scaling: in F1 bending resistance grows with $v_c$, in F2 it is fixed.

### Solved for $v_{sp}$

$$\boxed{v_{sp} = \frac{D B g + F_R - \dfrac{1}{6}\dfrac{H^2}{R}\sigma_Y + \eta_A C_{DP} v_c}{\eta_A\!\left(\dfrac{2L_p + L_s}{h} + C_{DP}\right)}}$$

### Trench velocity

$$v_t = v_{sp} - v_c \qquad \text{(positive = trench retreats / slab outpaces convergence)}$$

---

## Formulation 3 — Viscous Bending, $h \propto \mathcal{L}_{sp}$

Channel thickness scales linearly with slab length, reflecting longer slabs maintaining a thicker shear zone:

$$h_\text{eff} = h_0 + \mathcal{L}_{sp} \cdot \frac{h_\text{ref} - h_0}{\mathcal{L}_{sp,\text{ref}}} \qquad h_0 = 100\text{ km},\; h_\text{ref} = 250\text{ km},\; \mathcal{L}_{sp,\text{ref}} = 5000\text{ km}$$

### Force balance

Identical structure to F1 with $h \to h_\text{eff}$:

$$\boxed{D B g + F_R \;=\; \frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c + 2\eta_A \frac{L}{h_\text{eff}} v_{sp} + \eta_A \frac{L_s}{h_\text{eff}} v_{sp} + \eta_A C_{DP} \, v_t}$$

Note: $C_{DP}$ is independent of $h$ and is unchanged.

### Solved for $v_{sp}$

$$\boxed{v_{sp} = \frac{D B g + F_R - \dfrac{2}{3}\dfrac{H^3}{R^3}\eta_L \, v_c + \eta_A C_{DP} v_c}{\eta_A\!\left(\dfrac{2L_p + L_s}{h_\text{eff}} + C_{DP}\right)}}$$

### Trench velocity

$$v_t = v_{sp} - v_c \qquad \text{(positive = trench retreats / slab outpaces convergence)}$$

---

## Formulation 4 — Viscous Bending, $h \propto v_{sp}$ (closed-form quadratic)

Channel thickness scales linearly with the slab velocity — faster slabs maintain a thicker shear zone:

$$h_\text{eff}(v_{sp}) = h_0 + v_{sp} \cdot \frac{h_\text{ref} - h_0}{v_{sp,\text{ref}}} \qquad h_0 = 150\text{ km},\; h_\text{ref} = 200\text{ km},\; v_{sp,\text{ref}} = 5\text{ cm/yr}$$

### Force balance

$$\boxed{D B g + F_R \;=\; \frac{2}{3}\frac{H^3}{R^3}\eta_L \, v_c + 2\eta_A \frac{L}{h_\text{eff}(v_{sp})} v_{sp} + \eta_A \frac{L_s}{h_\text{eff}(v_{sp})} v_{sp} + \eta_A C_{DP} \, v_t}$$

### Solved for $v_{sp}$ (quadratic)

Let $P = D B g + F_R - \tfrac{2}{3}(H^3/R^3)\eta_L v_c + \eta_A C_{DP} v_c$ collect all terms that do not depend on $v_{sp}$ (the $+\eta_A C_{DP} v_c$ is the algebraic artifact from expanding $v_t = v_{sp} - v_c$).

Substituting $h_\text{eff} = h_0 + m \, v_{sp}$ (where $m = (h_\text{ref}-h_0)/v_{sp,\text{ref}}$) and multiplying through by $h_\text{eff}$ yields $a \, v_{sp}^2 + b \, v_{sp} + c = 0$, with:

$$a = \eta_A C_{DP} m, \qquad b = \eta_A(2L_p + L_s + C_{DP} h_0) - P m, \qquad c = -P h_0$$

Taking the positive root:

$$\boxed{v_{sp} = \frac{-b + \sqrt{b^2 - 4ac}}{2a}}$$

### Trench velocity

$$v_t = v_{sp} - v_c \qquad \text{(positive = trench retreats / slab outpaces convergence)}$$

---

## Comparison Summary

| | F1 Viscous | F2 Plastic | F3 $h(\mathcal{L}_{sp})$ | F4 $h(v_{sp})$ |
|--|--|--|--|--|
| Bending term | $\frac{2}{3}\frac{H^3}{R^3}\eta_L v_c$ | $\frac{1}{6}\frac{H^2}{R}\sigma_Y$ | same as F1 | same as F1 |
| Bending depends on $v_c$? | yes | **no** | yes | yes |
| Channel thickness $h$ | fixed, 200 km | fixed, 200 km | $h(\mathcal{L}_{sp})$ | $h(v_{sp})$ (quadratic) |
| Plate drag length | $L_p$ (ridge–trench) | same | same | same |
| Slab drag length | $L_s = D/\sin\delta$ | same | same | same |
| Denominator | $\eta_A((2L_p+L_s)/h + C_{DP})$ | same as F1 | $\eta_A((2L_p+L_s)/h_\text{eff} + C_{DP})$ | same structure |
| Numerator vs F1 | — | replace bending | replace $h$ | replace $h$ |
| Closed form? | yes | yes | yes | yes (quadratic) |
