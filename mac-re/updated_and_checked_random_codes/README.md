### updated_and_checked_random_codes

A collection of random codes made by me in the past... I have a bunch of them, so no I'm checking all of them to find out, which of them are useful or not... So, let's clear my repo out....

## T_calc_midplane
May 15 2025

It contains to .c codes: 
calc_temp_origin.c : the orginal midplane temperature calculator -- fuzzy and buggy, I have to say
calculate_temperature_midplane.c : some syntax errors are fixed, more comments are made, a bit fancier, but still not working fine:


# Protoplanetary Disk Temperature Model Notes

This document summarizes important points and potential improvements regarding the opacity-based temperature calculation model for a protoplanetary disk.

---

## 1. Opacity (κ) and Regimes

- The current opacity (`κ`) formula is piecewise-defined with temperature thresholds (e.g. 170 K, 210 K).
- **Check the source and physical justification** of the constants (`k0`, `a`, `b`) for each regime. They might come from simplified or outdated literature.
- Ensure **continuity at regime boundaries** to avoid unphysical jumps in opacity.
- Recommended references for opacity data:
  - Bell & Lin (1994)
  - Semenov et al. (2003)
- Consider plotting κ(T) to visually inspect for irregularities or discontinuities.

---

## 2. Temperature Calculation and Energy Balance

- Temperature is calculated based on an energy balance involving viscous heating and background radiation:
  
  \[
  T = \left( \frac{\text{opteff} \cdot \left(\frac{9}{4} \sigma_0 \nu \Omega_K^2 + 2 \sigma_\text{SB} T_\text{bg}^4 \right)}{2 \sigma_\text{SB}} \right)^{1/4}
  \]

- `opteff` (effective optical depth) is derived from the optical depth `optd`, which depends on surface density and opacity.
- The current formula for `opteff` is an approximate interpolation between optically thin and thick limits.
- Large jumps in κ or surface density (Σ) may cause instability in temperature.
- For more robust modeling, use separate formulae for optically thin and thick regimes with conditional branching.

---

## 3. Time Evolution Scheme

- The temperature update uses an explicit finite difference scheme for diffusion and advection terms:

  \[
  T_i^{n+1} = T_i^n + \Delta t \left( A \frac{T_{i+1} - 2T_i + T_{i-1}}{\Delta r^2} + B \frac{T_{i+1} - T_{i-1}}{2 \Delta r} \right)
  \]

- Explicit schemes require very small time steps (`dt`) for numerical stability.
- Consider switching to an implicit or semi-implicit scheme (e.g., Crank–Nicolson) for better stability and larger time steps.
- Make sure that the radius array `r` is updated correctly inside the time evolution loop; otherwise, calculations use stale values.

---

## 4. Recommendations for Further Development

- Export κ(T) and optical depth profiles to files for post-processing and visualization.
- Plot equilibrium and time-evolving temperature profiles to monitor the model’s behavior.
- Experiment with varying initial surface density (`Σ_0`) to see its effect on temperature distribution.
- Review and validate all physical constants and formulae with current literature to ensure realistic modeling.

---

## 5. Useful References

- Bell, K. R., & Lin, D. N. C. (1994). “Using opacity regimes to model protoplanetary disks.”
- Semenov, D., et al. (2003). “Rosseland and Planck mean opacities for protoplanetary disks.”

---

*This note is for personal reference and future improvements of the protoplanetary disk temperature model.*
