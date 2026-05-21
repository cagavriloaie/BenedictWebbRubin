// =============================================================================
// BWRS — Gas flow calculation through throttling devices
// Revision 3.0  |  05.2026
// Author  Constantin Agavriloaie  ·  CJ-RO  (original R 01.2004)
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
//   5. Reference conditions — one of five presets or a custom temperature:
//        Romania / EU  (DIN 1343)  :  0 °C and 15 °C / 101.325 kPa
//        ISO 13443 / UK / Italy    : 15 °C / 101.325 kPa
//        USA — AGA-3  (60 °F)      : 15.56 °C / 101.325 kPa
//        Russia — GOST 30319-1     : 20 °C / 101.325 kPa
//        Custom                    : user-entered T [°C]
//
// OUTPUTS (for each T / p / Δp set)
//   • Mixture density ρ(T, p)                  [kg/m³]
//   • Dynamic viscosity η(T, p)                [μPa·s]
//   • Compressibility factor Z(T, p)           [—]
//   • Mass flow rate Qm                        [kg/s  and  kg/h]
//   • Volumetric flow at reference conditions  [Nm³/h, Sm³/h, m³/h]
//       (0 °C / 101.325 kPa → Nm³/h; 15 °C / 101.325 kPa → Sm³/h; etc.)
//   • Volumetric flow at T and p               [m³/h]
//   • Mean gas velocity in pipe                [m/s]
//   • Permanent pressure loss                  [kPa]
//   • Throttling ratio β = d(t)/D(t)           [—]
//   • Reynolds number Re                       [—]
//   • Discharge coefficient C                  [—]
//   • Expansibility factor ε                   [—]
//   • Isentropic exponent κ                    [—]
//   • Pipe and orifice diameters at T          [mm]
//
// VALIDATION
//   Checks ISO 5167 applicability limits for pipe diameter,
//   orifice diameter, throttling ratio β, and Reynolds number Re —
//   displays an error message when a limit is exceeded.
//   Thermal expansion corrections applied to D and d before use:
//     pipe (carbon steel)        α_D = 12.2×10⁻⁶ /°C  (ISO 5167)
//     orifice (stainless steel)  α_d = 16.5×10⁻⁶ /°C
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
// VISCOSITY ALGORITHM
//   Chapman-Enskog dilute-gas model with high-density correction
//   (Lucas / Chung-Lee-Starling):
//     η = η₀(T, x) + kA/ξ · [exp(k₁·ρ_r) − exp(−k₂·ρ_r)]^k₃
//   where η₀ is the low-pressure mixture viscosity from kinetic theory,
//   ρ_r = ρ/ρ_c the reduced density, and ξ the viscosity reducing
//   parameter derived from Tc, M, Pc of the mixture.
//
// ISENTROPIC EXPONENT
//   Computed from the ideal-gas molar heat capacity at 20 °C:
//     κ = Cp° / (Cp° − R),   Cp° = Σ xᵢ·Cp°ᵢ,   R = 8.314 J/(mol·K)
//   Used in the ISO 5167 expansibility factor ε.
//
// FLOW ALGORITHM
//   Reynolds iteration until relative convergence |ΔQm/Qm| < 10⁻⁶.
//   Discharge coefficient C and expansibility factor ε are recalculated
//   at each iteration as a function of Re and β.
//
// UNCERTAINTY (ISO 5167-1)
//   Optional measurement uncertainty budget after each calculation.
//   Relative standard uncertainties [%] entered by user (or defaults):
//     u(C)  — discharge coefficient, ISO 5167 default per device type
//     u(ε)  — expansibility factor
//     u(d)  — orifice diameter
//     u(D)  — pipe diameter
//     u(Δp) — differential pressure transmitter
//     u(ρ)  — gas density (BWRS equation of state)
//   Combined relative standard uncertainty (law of propagation):
//     u²(Qm)/Qm² = u²(C) + u²(ε) + [2/(1−β⁴)]²·u²(d)
//                         + [2β⁴/(1−β⁴)]²·u²(D)
//                         + ¼·u²(Δp) + ¼·u²(ρ)
//   Expanded uncertainty: U(Qm) = 2·u(Qm)  at k=2, 95% confidence.
//   Absolute ± values reported for Qm [kg/s, kg/h] and Qv [ref. unit].
//
// PERSISTENCE
//   Composition saved to / loaded from  BWRSflow_comp.dat  (35 molar fractions).
//   Device configuration saved to / loaded from  BWRSflow_conf.dat
//   (device type integer, D [mm], d [mm]).
//
// SELF-TEST
//   Optional validation suite at startup: 8 gas compositions (CH4, C2H6,
//   CO2, H2, N2, Ar, standard NG, rich NG, synthetic air, H2S) tested for
//   density, viscosity, Z-factor, molar mass, and Qm across all 9 ISO 5167
//   device types.  Results checked against NIST WebBook, AGA-8, ISO 5167.
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

#include "BWRSflow.h"
#include "../Src/Ui/console.h"


// Navigation flow inside main():
//
//   ┌──────────────────────────────────────────────────────────────────────────┐
//   │  OUTER LOOP  — picks a new device type each iteration                    │
//   │                                                                          │
//   │  Load BWRSflow_conf.dat → reuse? ──yes──► DIAMETER LOOP (skip prompt)         │
//   │         │ no                                                             │
//   │         ▼                                                                │
//   │  Select device type (1–9)                                                │
//   │         │                                                                │
//   │  ┌──────┴───────────────────────────────────────────────────────────┐    │
//   │  │  DIAMETER LOOP  — same device type; re-entered when nav = 2      │    │
//   │  │                                                                  │    │
//   │  │  selectDiameters()  ←── also re-entered on nav = 2               │    │
//   │  │    PipeD ESC ──────────────────────────────────► outer loop      │    │
//   │  │    OrificeD ESC ──► back to PipeD                                │    │
//   │  │         │                                                        │    │
//   │  │  ┌──────┴───────────────────────────────────────────────────┐    │    │
//   │  │  │  INNER LOOP  — same device + diameters; vary T/P/ΔP      │    │    │
//   │  │  │                                                          │    │    │
//   │  │  │  readConditions()                                        │    │    │
//   │  │  │    ESC on T ──────────────────────────────► outer loop   │    │    │
//   │  │  │    ESC on P ──► back to T                                │    │    │
//   │  │  │    ESC on ΔP ─► back to P                                │    │    │
//   │  │  │         │                                                │    │    │
//   │  │  │  calcDensity / calcViscosity / calcMassFlow              │    │    │
//   │  │  │    diverge ──────────────────────────────────► outer loop│    │    │
//   │  │  │         │                                                │    │    │
//   │  │  │  PrintFlowResults / PrintUncertainty                     │    │    │
//   │  │  │         │                                                │    │    │
//   │  │  │  nav:  1 ──► top of inner loop (new T/P/ΔP)              │    │    │
//   │  │  │        2 ──► diameter loop (new D/d, same device type)   │    │    │
//   │  │  │        3 ──► outer loop  (new device type)               │    │    │
//   │  │  └──────────────────────────────────────────────────────────┘    │    │
//   │  └──────────────────────────────────────────────────────────────────┘    │
//   └──────────────────────────────────────────────────────────────────────────┘
//
int main() {
  InitConsole();

  // ── 1. Gas composition ──────────────────────────────────────────────────────
  double x[ARRAY_SIZE] = {};
  PrintBanner();
  ReadComposition(x);  // reads BWRSflow_comp.dat or prompts the user

  // ── 2. BWRS mixture constants (computed once for the whole session) ─────────
  using Clock = std::chrono::high_resolution_clock;
  using Us    = std::chrono::duration<double, std::micro>;

  auto t_mix0 = Clock::now();
  BwrConst      bwr       = calcBwrMixture(x);        // 11-parameter mixing rules
  double        kappa_mix = calcKappa(x);              // isentropic exponent for ISO 5167 ε
  MixtureThermo thermo    = calcMixtureThermo(x, bwr); // Tc, Pc, ξ  used by viscosity model
  double t_mix_us = Us(Clock::now() - t_mix0).count();

  // ── 3. Reference conditions (Romania/EU, ISO, AGA-3, GOST, or custom) ───────
  CountryRef ref;
  double ror_ref[2] = {};  // [0] density at chosen ref T,  [1] density at 15 °C
  SelectRefConditions(bwr, &ref, ror_ref);

  // ── Wizard: pipe D then orifice d, with ESC back-navigation ─────────────────
  // Returns false when the user presses ESC at the pipe-D prompt,
  // signalling that the caller should return to device-type selection.
  auto selectDiameters = [](double* d_int, double* d_orif) -> bool {
    enum class Step { PipeD, OrificeD };
    Step step = Step::PipeD;
    *d_int = *d_orif = 0.0;
    std::printf("\n");

    for (;;) {
      if (step == Step::PipeD) {
        // ── Pipe inner diameter ─────────────────────────────────────────────
        std::printf("  %s>%s %sPipe D    (20" U_DEG"C) [mm]%s"
                    "  [" U_LARR " reselect device]%s : %s",
                    COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_YELLOW,
                    COLOR_BOLD_WHITE, COLOR_RESET);
        if (!ReadDouble(d_int)) {
          std::printf("%s", COLOR_RESET);
          *d_int = 0.0;
          return false;  // ESC at D → caller returns to device-type selection
        }
        std::printf("%s", COLOR_RESET);
        if (*d_int > 0.0) {
          step = Step::OrificeD;
        } else {
          std::printf("%s  Invalid value " U_MDASH " must be positive.%s\n",
                      COLOR_BOLD_RED, COLOR_RESET);
        }

      } else {
        // ── Orifice diameter ────────────────────────────────────────────────
        std::printf("  %s>%s %sOrifice D (20" U_DEG"C) [mm]%s"
                    "  [" U_LARR " reenter pipe D]%s : %s",
                    COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_YELLOW,
                    COLOR_BOLD_WHITE, COLOR_RESET);
        if (!ReadDouble(d_orif)) {
          std::printf("%s", COLOR_RESET);
          *d_orif = *d_int = 0.0;
          step = Step::PipeD;  // ESC at d → back to pipe D
          continue;
        }
        std::printf("%s", COLOR_RESET);

        if (*d_orif <= 0.0) {
          std::printf("%s  Invalid value " U_MDASH " must be positive.%s\n",
                      COLOR_BOLD_RED, COLOR_RESET);
        } else if (*d_orif >= *d_int) {
          std::printf("%s  Orifice D must be less than pipe D (%g mm).%s\n",
                      COLOR_BOLD_RED, *d_int, COLOR_RESET);
        } else {
          double beta = *d_orif / *d_int;
          if (beta < 0.10 || beta > 0.80) {
            std::printf("%s  " U_BETA " = %.4f " U_MDASH
                        " outside ISO 5167 range [0.10, 0.80].%s\n",
                        COLOR_BOLD_RED, beta, COLOR_RESET);
          } else {
            return true;  // valid β — both diameters accepted
          }
        }
      }
    }
  };

  // ── Wizard: temperature, pressure, differential pressure with ESC back-nav ──
  // Returns false when the user presses ESC at the T prompt,
  // signalling that the caller should exit the inner loop and pick a new device.
  auto readConditions = [](double* T, double* P, double* dP) -> bool {
    enum class Step { Temperature, Pressure, DeltaP };
    Step step = Step::Temperature;
    *T = *P = *dP = 0.0;

    std::printf("\n\n%s  " U_HLINE U_HLINE " Measurement conditions%s\n", COLOR_YELLOW, COLOR_RESET);

    for (;;) {
      if (step == Step::Temperature) {
        // ── Gas temperature ───────────────────────────────────────────────
        std::printf("  %s>%s %sTemperature         [" U_DEG"C]%s"
                    "  [" U_LARR " change device]%s : %s",
                    COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_YELLOW,
                    COLOR_BOLD_WHITE, COLOR_RESET);
        if (!ReadDouble(T)) {
          std::printf("%s", COLOR_RESET);
          return false;  // ESC at T → caller exits inner loop, picks new device
        }
        std::printf("%s", COLOR_RESET);
        if (*T <= -273.15) {
          std::printf("%s  Temperature below absolute zero (-273.15" U_DEG"C).%s\n",
                      COLOR_BOLD_RED, COLOR_RESET);
        } else {
          step = Step::Pressure;
        }

      } else if (step == Step::Pressure) {
        // ── Absolute pressure ─────────────────────────────────────────────
        std::printf("  %s>%s %sPressure            [kPa]%s"
                    "  [" U_LARR " reenter T]%s : %s",
                    COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_YELLOW,
                    COLOR_BOLD_WHITE, COLOR_RESET);
        if (!ReadDouble(P)) {
          std::printf("%s", COLOR_RESET);
          step = Step::Temperature;  // ESC at P → back to T
        } else {
          std::printf("%s", COLOR_RESET);
          if (*P <= 0.0) {
            std::printf("%s  Pressure must be positive.%s\n",
                        COLOR_BOLD_RED, COLOR_RESET);
          } else {
            step = Step::DeltaP;
          }
        }

      } else {
        // ── Differential pressure ─────────────────────────────────────────
        std::printf("  %s>%s %sDifferential pressure [kPa]%s"
                    "  [" U_LARR " reenter P]%s : %s",
                    COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_YELLOW,
                    COLOR_BOLD_WHITE, COLOR_RESET);
        if (!ReadDouble(dP)) {
          std::printf("%s", COLOR_RESET);
          step = Step::Pressure;  // ESC at ΔP → back to P
        } else {
          std::printf("%s", COLOR_RESET);
          if (*dP <= 0.0) {
            std::printf("%s  Differential pressure must be positive.%s\n",
                        COLOR_BOLD_RED, COLOR_RESET);
          } else if (*dP >= *P) {
            std::printf("%s  Differential must be less than p = %g kPa.%s\n",
                        COLOR_BOLD_RED, *P, COLOR_RESET);
          } else {
            return true;  // T, P, ΔP all valid
          }
        }
      }
    }
  };

  // ══════════════════════════════════════════════════════════════════════════
  // OUTER LOOP — picks a new throttling device type each iteration.
  // ══════════════════════════════════════════════════════════════════════════
  for (;;) {
    int    tip_raw = 0;
    double d_int   = 0.0;
    double d_orif  = 0.0;
    bool   diameters_ready = false;  // true when D/d come from a saved config

    // ── 4a. Offer to reuse the saved device configuration ─────────────────────
    {
      int sv_tip; double sv_d_int, sv_d_orif;
      if (LoadConfig(&sv_tip, &sv_d_int, &sv_d_orif)
          && sv_tip >= TIP_MIN && sv_tip <= TIP_MAX) {
        std::printf(
            "\n%s  " U_HLINE U_HLINE " Saved configuration%s\n"
            "%s  Device      : %s%s\n"
            "%s  Pipe D      : %s%g mm%s\n"
            "%s  Orifice D   : %s%g mm%s\n",
            COLOR_YELLOW, COLOR_RESET,
            COLOR_BOLD_WHITE, tipName(static_cast<DeviceType>(sv_tip)), COLOR_RESET,
            COLOR_BOLD_WHITE, COLOR_BOLD_GREEN, sv_d_int, COLOR_RESET,
            COLOR_BOLD_WHITE, COLOR_BOLD_GREEN, sv_d_orif, COLOR_RESET);
        std::printf("  %s>%s %sReuse this configuration? [y/n] : %s",
                    COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET);
        if (AskYesNo()) {
          tip_raw        = sv_tip;
          d_int          = sv_d_int;
          d_orif         = sv_d_orif;
          diameters_ready = true;  // skip manual diameter entry below
        }
      }
    }

    // ── 4b. Select device type (skipped when reusing saved config) ────────────
    if (!diameters_ready) {
      std::printf("\n%s  " U_HLINE U_HLINE " Throttling device%s\n", COLOR_YELLOW, COLOR_RESET);
      std::printf("%s%s", COLOR_BOLD_WHITE, TIP_DISP);
      std::printf("%s", COLOR_RESET);
      tip_raw = ReadChoice(TIP_MIN, TIP_MAX);
    }

    DeviceType tip = static_cast<DeviceType>(tip_raw);

    // ══════════════════════════════════════════════════════════════════════════
    // DIAMETER LOOP — same device type; re-entered when the user picks nav = 2
    // ("change D / d") after seeing results.
    // ══════════════════════════════════════════════════════════════════════════
    bool need_new_diameters = !diameters_ready;  // first pass: skip if from config

    for (;;) {  // diameter loop

      // ── 4c. Select pipe D and orifice d ─────────────────────────────────────
      if (need_new_diameters) {
        if (!selectDiameters(&d_int, &d_orif)) break;  // ESC at D → outer loop

        std::printf("%s\n  Save configuration? [y/n] %s>%s ",
                    COLOR_BOLD_WHITE, COLOR_CYAN, COLOR_RESET);
        if (AskYesNo()) SaveConfig(tip_raw, d_int, d_orif);
      }
      need_new_diameters = true;  // subsequent passes always prompt (nav = 2)

      // ════════════════════════════════════════════════════════════════════════
      // INNER LOOP — same device and diameters; user varies T, P, ΔP.
      // ESC on T breaks this loop and lets the diameter loop exit to the outer
      // loop so the user can choose a new device type.
      // ════════════════════════════════════════════════════════════════════════
      bool change_diameters = false;

      for (;;) {  // inner loop
        double temperature   = 0.0;
        double pressure      = 0.0;
        double pressure_diff = 0.0;

        // ── 5. Read T, P, ΔP ────────────────────────────────────────────────
        if (!readConditions(&temperature, &pressure, &pressure_diff))
          break;  // ESC at T → exit inner loop → exit diameter loop → outer loop

        std::printf("  %c", 7);  // terminal bell — signals start of computation

        // ── 6a. Gas density via BWRS bisection ──────────────────────────────
        int  iter_rho = 0;
        auto t_rho0   = Clock::now();
        double ro     = calcDensity(temperature, pressure / KPA_PER_ATM, bwr, &iter_rho);
        double t_rho_us = Us(Clock::now() - t_rho0).count();

        if (!std::isfinite(ro) || ro <= 0.0) {
          std::printf(
              "%s\n  BWRS density did not converge (T=%.1f" U_DEG"C, P=%.1f kPa)"
              " " U_MDASH " check that conditions are in the gas-phase region.%s\n",
              COLOR_BOLD_RED, temperature, pressure, COLOR_RESET);
          break;  // abort this T/P set → outer loop (new device)
        }

        // ── 6b. Dynamic viscosity (Chapman-Enskog + high-density correction) ─
        auto t_eta0 = Clock::now();
        double eta  = calcViscosity(temperature, ro, x, thermo.m_rocCrit, thermo.m_csi);
        double t_eta_us = Us(Clock::now() - t_eta0).count();

        // ── 6c. Mass flow rate (ISO 5167 Reynolds iteration) ─────────────────
        int        iter_qm = 0;
        FlowResult flow;
        auto t_qm0  = Clock::now();
        double qm   = calcMassFlow(pressure_diff, pressure, temperature,
                                   tip, d_int, d_orif, ro, eta,
                                   &flow, kappa_mix, &iter_qm);
        double t_qm_us = Us(Clock::now() - t_qm0).count();

        if (qm == 0.0) break;  // ISO 5167 applicability limit violated → outer loop

        // ── 7. Display results ───────────────────────────────────────────────
        PrintFlowResults(tip, d_int, d_orif, temperature, pressure, pressure_diff,
                         ro, eta, qm, flow, bwr, ref, ror_ref, x);
        PrintCpuProfile(t_mix_us, t_rho_us, t_eta_us, t_qm_us, iter_rho, iter_qm);
        PrintUncertainty(tip, flow, qm, ref, ror_ref);

        // ── 8. Navigation ────────────────────────────────────────────────────
        //   1 → new T/P/ΔP, same device and diameters   (continue inner loop)
        //   2 → new D and d, same device type            (restart diameter loop)
        //   3 → new device type                          (exit both loops)
        std::printf("\n  %s1%s new conditions   %s2%s change D / d"
                    "   %s3%s change device   %sESC%s quit\n\n",
                    COLOR_BOLD_GREEN, COLOR_RESET, COLOR_BOLD_GREEN, COLOR_RESET,
                    COLOR_BOLD_GREEN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET);
        std::printf("  %s>%s ", COLOR_CYAN, COLOR_RESET);

        int nav = ReadChoice(1, 3);
        if (nav == 2) { change_diameters = true; break; }  // restart diameter loop
        if (nav == 3) break;  // exit inner loop → exit diameter loop → outer loop
        // nav == 1: continue inner loop (new T/P/ΔP, same device and diameters)
      }  // end inner loop

      if (!change_diameters) break;  // exit diameter loop → outer loop (new device)
    }  // end diameter loop
  }  // end outer loop
}
