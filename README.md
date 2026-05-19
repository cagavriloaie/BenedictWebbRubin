// =============================================================================
// BWRS — Gas flow calculation through throttling devices
// Revision 3.0  |  05.2026
// Eng. Agavriloaie Constantin  (original R 01.2004)
// ELCOST Impex
// ============================================================================= 
//
// DESCRIPTION
//   Calculates density, dynamic viscosity, and flow rates for a gas mixture
//   of up to 35 components, based on:
//     • Benedict-Webb-Rubin-Starling (BWRS) equation of state for density;
//     • Chapman-Enskog model with high-density correction for viscosity;
//     • ISO 5167-2/3/4:2003 (Reader-Harris/Gallagher equation)
//       for flow rate through throttling devices.
//
// INPUTS
//   1. Mixture composition — molar fractions for each of the
//      35 components (methane, ethane, propane, ... acetylene).
//   2. Throttling device type (1–9):
//        [1] Orifice plate — corner taps
//        [2] Orifice plate — flange taps
//        [3] Orifice plate — D and D/2 taps
//        [4] ISA 1932 nozzle
//        [5] Long-radius nozzle
//        [6] Classical Venturi tube — rough-cast convergent
//        [7] Classical Venturi tube — machined convergent
//        [8] Classical Venturi tube — rough-welded sheet-metal convergent
//        [9] Venturi nozzle
//   3. Pipe inner diameter D [mm] and orifice diameter d [mm]
//      (at 20 °C reference temperature).
//   4. Temperature T [°C], absolute pressure p [kPa], and differential
//      pressure Δp [kPa] at the measurement point.
//
// OUTPUTS (for each T / p / Δp set)
//   • Mixture density ρ(T, p)                  [kg/m³]
//   • Dynamic viscosity η(T, p)                [μPa·s]
//   • Mass flow rate Qm                        [kg/s]
//   • Volumetric flow at selectable reference conditions
//       (0 °C / 101.325 kPa → Nm³/h; 15 °C / 101.325 kPa → Sm³/h; etc.)
//   • Volumetric flow at T and p               [m³/h]
//   • Mean gas velocity in pipe                [m/s]
//   • Permanent pressure loss                  [kPa]
//   • Throttling ratio β and Reynolds number Re
//
// VALIDATION
//   Checks ISO 5167 applicability limits for pipe diameter,
//   orifice diameter, throttling ratio β, and Reynolds number Re —
//   displays an error message when a limit is exceeded.
//
// DENSITY ALGORITHM
//   Bisection on the BWRS equation (Starling 1973, 11 parameters) until 
//   convergence |p_calc − p| < 5×10⁻⁴ atm.  Equation of state:
//     p = ρRT + (B₀RT − A₀ − C₀/T² + D₀/T³ − E₀/T⁴)ρ²
//             + (bRT − a − d/T)ρ³ + α(a + d/T)ρ⁶
//             + (c/T²)ρ³(1 + γρ²)exp(−γρ²)
//   Mixture constants are computed once from composition using quadratic
//   mixing rules (A₀–E₀, γ) and cubic mixing rules (a–d, α).
//
// FLOW ALGORITHM
//   Reynolds iteration until relative convergence |ΔQm/Qm| < 10⁻⁶.
//   Discharge coefficient C and expansibility factor ε are recalculated
//   at each iteration as a function of Re and β.
//
// REFERENCES
//   • ISO 5167-2:2003 — Orifice plates (Reader-Harris/Gallagher equation)
//   • ISO 5167-3:2003 — Nozzles and Venturi nozzles
//   • ISO 5167-4:2003 — Venturi tubes
//   • Starling K.E. (1973) — Fluid Thermodynamic Properties for Light
//     Petroleum Systems, Gulf Publishing Co., Houston
//   • Nishiumi H., Saito S. (1975) — J. Chem. Eng. Japan 8(5), 356–360
//   • Benedict, Webb, Rubin (1940) — J. Chem. Phys. 8, 334
// =============================================================================
