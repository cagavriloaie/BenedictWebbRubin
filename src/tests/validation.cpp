#include "validation.h"
#include "../Common/constants.h"
#include "../Common/types.h"
#include "../Data/components.h"
#include "../Formulas/bwrs.h"
#include "../Formulas/viscosity.h"
#include "../Standards/iso5167.h"
#include "../Ui/console.h"

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <io.h>
#else
#  include <unistd.h>
#endif

#include <cmath>
#include <cstdio>
#include <cstring>

void runValidationTest() {
  std::printf("\n%s  " U_HLINE U_HLINE " Implementation validation test%s\n", COLOR_YELLOW, COLOR_RESET);

  int npass = 0, ntotal = 0;

  // chk: compares a value against an interval and prints PASS/FAIL
  // src = reference source label, printed at the end of each line
  auto chk = [&](const char* label, double v, double lo, double hi, const char* src) {
    ntotal++;
    bool ok = (v >= lo && v <= hi);
    if (ok) npass++;
    int vlen = (int)std::strlen(label) - Utf8ExtraBytes(label);
    std::printf("    %s%s", COLOR_BOLD_WHITE, label);
    for (int k = vlen; k < 54; k++) std::putchar(' ');
    char ivl[32];
    std::snprintf(ivl, sizeof(ivl), "[%.4f, %.4f]", lo, hi);
    std::printf("%s%-10.5f%s  %-24s  %s%-4s%s  %s[%s]%s\n",
                COLOR_BOLD_GREEN, v, COLOR_RESET,
                ivl,
                ok ? COLOR_BOLD_GREEN : COLOR_BOLD_RED, ok ? "PASS" : "FAIL", COLOR_RESET,
                COLOR_BOLD_WHITE, src, COLOR_RESET);
  };

  // test_comp: runs the test suite (density x2, viscosity, mass flow rate)
  //   for composition xc with the given parameters and reference intervals.
  //   d_flow_mm = 0 → skip flow test.
  auto test_comp = [&](const char* comp_name,
                       const double* xc,
                       double rho20_lo, double rho20_hi,
                       double rho0_lo,  double rho0_hi,
                       double eta_lo,   double eta_hi,
                       double D_mm,     double d_flow_mm,
                       double p_kpa,    double dp_kpa,
                       bool leading_nl = true) {
    std::printf(leading_nl ? "\n  %s" U_RAQUO " %s%s\n" : "  %s" U_RAQUO " %s%s\n",
                COLOR_YELLOW, comp_name, COLOR_RESET);

    BwrConst bt = calcBwrMixture(xc);
    MixtureThermo thermo_c = calcMixtureThermo(xc, bt);

    double rho1 = calcDensity(20.0, 1.0, bt);
    double rho2 = calcDensity( 0.0, 1.0, bt);
    chk("Density " U_RHO " 20" U_DEG"C, 101.325 kPa [kg/m" U_SUP3 "]", rho1, rho20_lo, rho20_hi, "BWRS 1973");
    chk("Density " U_RHO "  0" U_DEG"C, 101.325 kPa [kg/m" U_SUP3 "]", rho2, rho0_lo,  rho0_hi,  "BWRS 1973");
    chk("Viscosity " U_ETA " 20" U_DEG"C, 101.325 kPa [" U_MICRO "Pa" U_CDOT "s]",
        calcViscosity(20.0, rho1, xc, thermo_c.m_rocCrit, thermo_c.m_csi) * PA_TO_MICRO_PA, eta_lo, eta_hi, "CE/NIST");
    double T_K1 = 20.0 + KELVIN_OFFSET;
    double Z1   = 1.0 * bt.m_molarMass / (rho1 * GAS_CONSTANT_R * T_K1);
    chk("Compressibility factor Z 20" U_DEG"C, 101.325 kPa [-]", Z1, 0.975, 1.003, "BWRS/ideal");
    double mx_ref = 0.0;
    for (int ii = 1; ii <= NUM_COMPONENTS; ii++) mx_ref += xc[ii] * MOLAR_MASS_TABLE[ii];
    chk("Molar mass M [g/mol]", bt.m_molarMass,
        mx_ref * (1.0 - 1e-8), mx_ref * (1.0 + 1e-8), "mixing rule");

    if (d_flow_mm > 0.0) {
      double rho_f = calcDensity(20.0, p_kpa / KPA_PER_ATM, bt);
      double eta_f = calcViscosity(20.0, rho_f, xc, thermo_c.m_rocCrit, thermo_c.m_csi);
      double A_o = PI / 4.0 * (d_flow_mm * 1e-3) * (d_flow_mm * 1e-3);
      static const struct { int tip; const char* abv; double alpha_nom; double eps_nom; }
        tipuri[] = {
          {1, "Dia-U",   0.609, 0.993},
          {2, "Dia-F",   0.608, 0.993},
          {3, "Dia-D/2", 0.608, 0.993},
          {4, "Aj-ISA",  0.997, 0.997},
          {5, "Aj-RL",   1.005, 0.997},
          {6, "Ven-B",   0.997, 0.997},
          {7, "Ven-P",   1.008, 0.997},
          {8, "Ven-T",   0.998, 0.997},
          {9, "Aj-V",    0.996, 0.997},
        };
      double qm_diau = 0.0;
      for (const auto& t : tipuri) {
        FlowResult fr;
        double qm_f = calcMassFlow(dp_kpa, p_kpa, 20.0,
                                   static_cast<DeviceType>(t.tip),
                                   D_mm, d_flow_mm, rho_f, eta_f, &fr);
        if (t.tip == 1) qm_diau = qm_f;
        if (qm_f > 0.0) {
          double qm_est = t.alpha_nom * t.eps_nom * A_o
                        * std::sqrt(2000.0 * dp_kpa * rho_f);
          char lbl[80];
          std::snprintf(lbl, sizeof(lbl),
                        "Mass flow Qm D=%g,d=%g,p=%g," U_DELTA "p=%g %s [kg/s]",
                        D_mm, d_flow_mm, p_kpa, dp_kpa, t.abv);
          chk(lbl, qm_f, qm_est * 0.85, qm_est * 1.15, "ISO 5167");
        } else {
          ntotal++;
          std::printf("    %sQm test %s: ISO 5167 ERROR%s\n", COLOR_BOLD_RED, t.abv, COLOR_RESET);
        }
      }
      if (qm_diau > 0.0 && rho2 > 0.0) {
        double qv_n     = qm_diau / rho2 * SECONDS_PER_HOUR;
        double qm_est_u = tipuri[0].alpha_nom * tipuri[0].eps_nom * A_o
                        * std::sqrt(2000.0 * dp_kpa * rho_f);
        double qv_est   = qm_est_u / rho2 * SECONDS_PER_HOUR;
        chk("Volumetric flow Qv Dia-U, 0" U_DEG"C ref [Nm" U_SUP3 "/h]",
            qv_n, qv_est * 0.85, qv_est * 1.15, "ISO 5167/BWRS");
      }
    }
  };

  // sec_box: prints a description box before a test section.
  auto sec_box = [](const char* sec, const char* desc, int n) {
    auto hl = [](int cnt) { for (int k = 0; k < cnt; ++k) std::printf(U_HLINE); };
    std::printf("\n  " U_TL); hl(22); std::printf(U_MT); hl(74);
    std::printf(U_MT); hl(7); std::printf(U_TR "\n");
    std::printf("  " U_VLINE " %-20s " U_VLINE " %-72s " U_VLINE " %-5s " U_VLINE "\n",
                "Section", "Checks", "Tests");
    std::printf("  " U_ML); hl(22); std::printf(U_CROSS); hl(74);
    std::printf(U_CROSS); hl(7); std::printf(U_MR "\n");

    char nt[8]; std::snprintf(nt, sizeof(nt), "+%d", n);
    std::printf("  " U_VLINE " %-20s " U_VLINE " %-72s " U_VLINE " %-5s " U_VLINE "\n",
                sec, desc, nt);
    std::printf("  " U_BL); hl(22); std::printf(U_MB); hl(74);
    std::printf(U_MB); hl(7); std::printf(U_BR "\n");
  };

  // ── Test compositions ──────────────────────────────────────────────────────
  // Component indices: CH4=1, C2H6=2, C3H8=3, i-C4=4, n-C4=5, H2=24,
  //   CO=25, H2S=26, He=27, Ar=28, N2=29, O2=30, CO2=31, Ammonia=34
  //
  // ── Density references (BWRS vs ideal gas at 1 atm) ─────────────────────
  // Ideal gas 20°C: rho = M / (R × 293.15) = M / 24.0554 [kg/m³]
  // Ideal gas  0°C: rho = M / (R × 273.15) = M / 22.4136 [kg/m³]
  //
  //  Composition    M [g/mol] rho_id_20  rho_id_0   B₂(20°C)[L/mol]  rho_BWRS_20
  //  CH4             16.043    0.6672     0.7157     −0.0447           0.668
  //  C2H6            30.070    1.2500     1.3415     −0.195            1.260
  //  CO2             44.011    1.8296     1.9635     −0.147            1.841
  //  H2               2.016    0.0838     0.0899     +0.022            0.0837
  //  N2              28.016    1.1647     1.2498     −0.0044           1.162
  //  Ar              39.944    1.6604     1.7821     −0.018            1.661
  //
  // B₂ source: computed from BWRS Starling 1973 parameters via:
  //   B₂(T) = B₀ − A₀/(RT) − C₀/(RT³) + D₀/(RT⁴) − E₀/(RT⁵)
  //
  // ── Viscosity references (Chapman-Enskog + density correction) ───────────
  //  Composition  η_calc [μPa·s]   η_NIST 20°C [μPa·s]   CE error
  //  CH4           11.46            10.99                   +4.3 %
  //  C2H6           9.68             9.36                   +3.4 %
  //  CO2           15.49            14.89                   +4.0 %
  //  H2             8.79             8.91                   −1.3 %
  //  N2            18.47            17.54                   +5.3 %
  //  Ar            22.41            22.74                   −1.4 %
  //
  // NIST source: https://webbook.nist.gov (Transport Properties, 100 kPa, 20°C)
  //
  // ── Mass flow references ─────────────────────────────────────────────────
  // All tests: D=200 mm, d=80 mm, β=0.4, T=20°C, p=500 kPa, Δp=5 kPa
  // β=0.4 chosen to satisfy all ISO 5167 constraints simultaneously.
  // Qm_est = alpha_nom × eps_nom × A_o × √(2000 × Δp_kPa × ρ)
  //   A_o = π/4 × (0.08)² = 5.027e−3 m²
  // Accepted range: ±15% of Qm_est (covers Re and BWRS variation)

  double xCH4[ARRAY_SIZE] = {}; xCH4[1] = 1.0;   // pure methane

  double xC2H6[ARRAY_SIZE] = {}; xC2H6[2] = 1.0;  // pure ethane

  double xCO2[ARRAY_SIZE] = {}; xCO2[31] = 1.0;   // pure CO2

  double xH2[ARRAY_SIZE] = {}; xH2[24] = 1.0;     // pure hydrogen

  double xN2[ARRAY_SIZE] = {}; xN2[29] = 1.0;     // pure nitrogen

  double xAr[ARRAY_SIZE] = {}; xAr[28] = 1.0;     // pure argon

  double xGN[ARRAY_SIZE] = {};                     // standard natural gas
  xGN[1] = 0.900; xGN[2] = 0.060; xGN[3] = 0.020;
  xGN[29] = 0.015; xGN[31] = 0.005;

  double xGNB[ARRAY_SIZE] = {};                    // rich natural gas
  xGNB[1] = 0.850; xGNB[2] = 0.100; xGNB[3] = 0.030;
  xGNB[29] = 0.010; xGNB[31] = 0.010;

  sec_box("Test compositions",
          "8 compositions: rho(20/0C), eta, Z, molar mass, Qm 9 ISO 5167 types", 120);
  std::printf("%s  " U_HLINE U_HLINE " Legend%s\n", COLOR_YELLOW, COLOR_RESET);

  std::printf("\n%s  Physical quantities:%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("    " U_RHO "   [kg/m" U_SUP3 "]   density\n");
  std::printf("    " U_ETA "   [" U_MICRO "Pa" U_CDOT "s]   dynamic viscosity\n");
  std::printf("    Z   [-]       compressibility factor\n");
  std::printf("    M   [g/mol]   molar mass\n");
  std::printf("    Qm  [kg/s]    mass flow rate\n");
  std::printf("    Qv  [Nm" U_SUP3 "/h]   volumetric flow (normal conditions: 0" U_DEG"C, 101.325 kPa)\n");

  std::printf("\n%s  Flow device parameters:%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("    D   [mm]      pipe diameter\n");
  std::printf("    d   [mm]      orifice / throat diameter\n");
  std::printf("    " U_BETA "   [-]       diameter ratio  " U_BETA " = d/D\n");
  std::printf("    p   [kPa]     absolute pressure\n");
  std::printf("    " U_DELTA "p  [kPa]     differential pressure\n");

  std::printf("\n%s  ISO 5167 devices:%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("    Dia-U   = Orifice plate, corner taps\n");
  std::printf("    Dia-F   = Orifice plate, flange taps\n");
  std::printf("    Dia-D/2 = Orifice plate, D and D/2 taps\n");
  std::printf("    Aj-ISA  = ISA 1932 nozzle\n");
  std::printf("    Aj-RL   = Long-radius nozzle (ASME long-radius)\n");
  std::printf("    Ven-B   = Venturi tube " U_MDASH " rough-cast convergent      (C " U_APPROX " 0.984)\n");
  std::printf("    Ven-P   = Venturi tube " U_MDASH " machined convergent        (C " U_APPROX " 0.995)\n");
  std::printf("    Ven-T   = Venturi tube " U_MDASH " welded sheet-metal conv.   (C " U_APPROX " 0.985)\n");
  std::printf("    Aj-V    = Venturi nozzle                             (C " U_APPROX " 0.983)\n");

  std::printf("\n%s  Reference sources [right column of each test]:%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("    BWRS 1973     = Benedict-Webb-Rubin-Starling (Starling, K.E., 1973)\n");
  std::printf("    CE/NIST       = Chapman-Enskog + density correction / NIST WebBook\n");
  std::printf("    BWRS/ideal    = BWRS vs ideal gas at 101.325 kPa\n");
  std::printf("    BWRS/NIST     = BWRS compared with NIST WebBook tabulated data\n");
  std::printf("    NIST/BWRS     = NIST reference, calculated with BWRS\n");
  std::printf("    BWRS/AGA-8    = BWRS compared with AGA-8 data\n");
  std::printf("    BWRS/estim.   = BWRS compared with engineering estimate\n");
  std::printf("    BWRS/Z        = " U_RHO " monotonicity with Z-factor correction\n");
  std::printf("    ISO 5167      = ISO 5167-2/3/4:2003 (Reader-Harris/Gallagher)\n");
  std::printf("    ISO 5167/BWRS = Qv = Qm(ISO 5167) / " U_RHO "(BWRS, 0" U_DEG"C, 101.325 kPa)\n");
  std::printf("    mixing rule   = linear rule  M = " U_SIGMA " xi" U_CDOT "Mi\n");
  std::printf("    ideal gas     = monotonicity test  " U_RHO "(2p)/" U_RHO "(p) " U_APPROX " 2\n");
  std::printf("    ideal gas" U_TIMES "Z   = as above, with Z-factor correction\n");
  std::printf("    single-phase  = Z " U_IN " [0.75, 1.02] check (gas phase)\n");
  std::printf("    quasi-ideal   = " U_RHO "/" U_RHO "_ideal " U_IN " [0.990, 1.025] at p " U_LEQ " 500 kPa\n");
  std::printf("    monotone T    = " U_RHO "(T1) > " U_RHO "(T2) at T1 < T2  (density increases with cooling)\n");
  std::printf("    monotone P    = " U_RHO "(p1) < " U_RHO "(p2) at p1 < p2  (density increases with pressure)\n");
  std::printf("    M + Z         = density ratio " U_APPROX " (Ma/Mb) " U_TIMES " (Zb/Za)\n");
  std::printf("    CE/kinetic    = Chapman-Enskog from kinetic theory (monotone " U_ETA " vs. T)\n");
  std::printf("\n");

  //           Composition   rho_20[kg/m3]  rho_0[kg/m3]   eta[µPa·s]   D    d   p    dp
  test_comp("Methane (CH4)",  xCH4, 0.660,0.672, 0.712,0.724, 10.5,12.5, 200,80, 500,5, false);
  test_comp("Ethane (C2H6)", xC2H6,1.220,1.280, 1.310,1.375,  8.5,10.8, 200,80, 500,5);
  test_comp("CO2",           xCO2, 1.800,1.860, 1.930,2.000, 14.0,16.8, 200,80, 500,5);
  test_comp("Hydrogen (H2)", xH2,  0.082,0.085, 0.088,0.091,  8.2, 9.5, 200,80, 500,5);
  test_comp("Nitrogen (N2)", xN2,  1.155,1.175, 1.240,1.260, 16.5,20.5, 200,80, 500,5);
  test_comp("Argon (Ar)",    xAr,  1.650,1.672, 1.772,1.794, 21.0,23.5, 200,80, 500,5);
  test_comp("Std NG",        xGN,  0.730,0.748, 0.784,0.802, 10.5,12.5, 200,80, 500,5);
  test_comp("Rich NG",       xGNB, 0.765,0.790, 0.820,0.848, 10.2,12.5, 200,80, 500,5);

  sec_box("Pure CH4 standalone",
          "Pure CH4: rho monotonicity vs P (x2); Z factor at 20 atm", 3);
  std::printf("  %s" U_RAQUO " Standalone tests pure CH4%s\n", COLOR_YELLOW, COLOR_RESET);
  {
    BwrConst bt_a = calcBwrMixture(xCH4);
    double T_a    = 20.0 + KELVIN_OFFSET;
    double rho_05 = calcDensity(20.0,  0.5, bt_a);   // 50.66 kPa
    double rho_10 = calcDensity(20.0,  1.0, bt_a);   // 101.325 kPa
    double rho_20 = calcDensity(20.0, 20.0, bt_a);   // ≈ 2.026 MPa
    chk("Monotone " U_RHO "(1atm)/" U_RHO "(0.5atm) CH4 [-]",
        rho_10 / rho_05, 1.90, 2.10, "ideal gas");
    chk("Monotone " U_RHO "(20atm)/" U_RHO "(1atm) CH4 [-]",
        rho_20 / rho_10, 17.0, 23.0, "ideal gas" U_TIMES "Z");
    double Z_hi = 20.0 * bt_a.m_molarMass / (rho_20 * GAS_CONSTANT_R * T_a);
    chk("Z factor CH4 20" U_DEG"C, 20 atm [-]", Z_hi, 0.88, 0.98, "BWRS/NIST");
  }

  sec_box("High-pressure rho",
          "CH4/N2/CO2/GN at 500-5066 kPa: rho BWRS, Z, monotonicity, consistency", 8);
  // References:
  //   CH4  500 kPa, 293.15 K: Z≈0.9986 (NIST)       → ρ≈3.30  kg/m³ (≈ideal)
  //   CH4 5066 kPa, 293.15 K: Z≈0.911  (BWRS/AGA-8) → ρ≈36.6  kg/m³
  //   N2  5066 kPa, 293.15 K: Z≈0.991  (NIST)       → ρ≈58.8  kg/m³
  //   CO2  500 kPa, 293.15 K: Z≈0.974  (NIST)       → ρ≈9.27  kg/m³ (non-ideal)
  //   Std NG 5066 kPa, 20°C : BWRS estimate          → ρ≈38-40 kg/m³
  std::printf("  %s" U_RAQUO " Density at high pressure%s\n", COLOR_YELLOW, COLOR_RESET);
  {
    BwrConst bCH4 = calcBwrMixture(xCH4);
    BwrConst bN2  = calcBwrMixture(xN2);
    BwrConst bCO2 = calcBwrMixture(xCO2);
    BwrConst bGN  = calcBwrMixture(xGN);

    double T_K   = 20.0 + KELVIN_OFFSET;
    double p500  = 500.0  / KPA_PER_ATM;   //  4.935 atm
    double p5066 = 5066.0 / KPA_PER_ATM;   // 50.00  atm

    double rCH4_500  = calcDensity(20.0, p500,  bCH4);
    double rCH4_5066 = calcDensity(20.0, p5066, bCH4);
    chk("Density CH4 20" U_DEG"C,  500 kPa [kg/m" U_SUP3 "]",
        rCH4_500,  3.24, 3.35, "NIST/BWRS");
    chk("Density CH4 20" U_DEG"C, 5066 kPa [kg/m" U_SUP3 "]",
        rCH4_5066, 34.5, 39.0, "BWRS/AGA-8");

    // Z at 50 atm: B₂(CH4,293K)≈−44 cm³/mol → Z_virial≈0.909; BWRS→0.911
    double Z50 = p5066 * bCH4.m_molarMass
               / (rCH4_5066 * GAS_CONSTANT_R * T_K);
    chk("Z factor CH4 20" U_DEG"C, 50 atm [-]",
        Z50, 0.88, 0.95, "BWRS/AGA-8");

    double rN2_5066 = calcDensity(20.0, p5066, bN2);
    chk("Density N2  20" U_DEG"C, 5066 kPa [kg/m" U_SUP3 "]",
        rN2_5066, 56.0, 62.0, "NIST/BWRS");

    double rCO2_500 = calcDensity(20.0, p500, bCO2);
    chk("Density CO2 20" U_DEG"C,  500 kPa [kg/m" U_SUP3 "]",
        rCO2_500, 8.90, 9.70, "NIST/BWRS");

    double rGN_5066 = calcDensity(20.0, p5066, bGN);
    chk("Density std NG 20" U_DEG"C, 5066 kPa [kg/m" U_SUP3 "]",
        rGN_5066, 35.0, 42.0, "BWRS/estim.");

    chk("Monotone " U_RHO " CH4: " U_RHO "(5066)/" U_RHO "(500) [-]",
        rCH4_5066 / rCH4_500, 9.5, 11.5, "BWRS/Z");

    double rN2_500 = calcDensity(20.0, p500, bN2);
    chk("Consistency: " U_RHO "(CO2)/" U_RHO "(N2) at 500 kPa [-]",
        rCO2_500 / rN2_500, 1.50, 1.65, "M + Z");
  }

  sec_box("Rich NG 4Tx4P",
          "Rich NG 4Tx4P: Z, rho/rho_ideal, monotonicity vs T and vs P", 30);
  // Temperatures: -10°C (winter), 10°C, 30°C, 60°C (summer/compressor)
  // Pressures:    101 kPa (ref.), 500 kPa, 2000 kPa, 7000 kPa (pipeline)
  // T_min = -10°C = 263.15 K > Tc_pseudo(GNB) ≈ 208 K → single-phase gas.
  std::printf("  %s" U_RAQUO " Rich NG: T " U_TIMES " P matrix  (4 " U_TIMES " 4 conditions)%s\n",
              COLOR_YELLOW, COLOR_RESET);
  {
    BwrConst bGNB = calcBwrMixture(xGNB);

    const double T_c[]   = { -10.0,   10.0,   30.0,   60.0 };
    const double P_kpa[] = { 101.325, 500.0, 2000.0, 7000.0 };
    const char*  T_lbl[] = { "-10",  " 10",  " 30",  " 60"  };
    const char*  P_lbl[] = { "  101", "  500", " 2000", " 7000" };
    constexpr int nT = 4, nP = 4;

    double rho[nT][nP] = {};

    for (int ti = 0; ti < nT; ti++) {
      for (int pi = 0; pi < nP; pi++) {
        double T_K      = T_c[ti] + KELVIN_OFFSET;
        double p_atm    = P_kpa[pi] / KPA_PER_ATM;
        rho[ti][pi]     = calcDensity(T_c[ti], p_atm, bGNB);
        double Z        = p_atm * bGNB.m_molarMass
                        / (rho[ti][pi] * GAS_CONSTANT_R * T_K);
        char lbl[64];
        std::snprintf(lbl, sizeof(lbl),
                      "Z GNB T=%s" U_DEG"C P=%s kPa",
                      T_lbl[ti], P_lbl[pi]);
        chk(lbl, Z, 0.75, 1.02, "single-phase");
      }
    }

    // At p ≤ 500 kPa: ρ/ρ_ideal ∈ [0.990, 1.025] (quasi-ideal)
    // Upper limit 1.025 covers T=-10°C where B_mix≈-85 cm³/mol → Z≈0.981
    for (int ti = 0; ti < nT; ti++) {
      for (int pi = 0; pi < 2; pi++) {
        double T_K    = T_c[ti] + KELVIN_OFFSET;
        double rho_id = P_kpa[pi] * 1e3 * bGNB.m_molarMass * 1e-3
                      / (8.314 * T_K);
        char lbl[64];
        std::snprintf(lbl, sizeof(lbl),
                      U_RHO "/" U_RHO "_id GNB T=%s" U_DEG"C P=%s kPa",
                      T_lbl[ti], P_lbl[pi]);
        chk(lbl, rho[ti][pi] / rho_id, 0.990, 1.025, "quasi-ideal");
      }
    }

    // Monotonicity vs T at P = 2000 kPa: ρ(T_i)/ρ(T_{i+1}) ≈ 1.07–1.10
    chk("GNB: " U_RHO "(-10" U_DEG"C)/" U_RHO "( 10" U_DEG"C) at 2000 kPa",
        rho[0][2] / rho[1][2], 1.03, 1.15, "monotone T");
    chk("GNB: " U_RHO "( 10" U_DEG"C)/" U_RHO "( 30" U_DEG"C) at 2000 kPa",
        rho[1][2] / rho[2][2], 1.03, 1.15, "monotone T");
    chk("GNB: " U_RHO "( 30" U_DEG"C)/" U_RHO "( 60" U_DEG"C) at 2000 kPa",
        rho[2][2] / rho[3][2], 1.03, 1.15, "monotone T");

    // Monotonicity vs P at T = 10°C
    chk("GNB: " U_RHO "(500)/" U_RHO "(101) at 10" U_DEG"C",
        rho[1][1] / rho[1][0], 4.70, 5.10, "monotone P");
    chk("GNB: " U_RHO "(2000)/" U_RHO "(500) at 10" U_DEG"C",
        rho[1][2] / rho[1][1], 3.70, 4.30, "monotone P");
    chk("GNB: " U_RHO "(7000)/" U_RHO "(2000) at 10" U_DEG"C",
        rho[1][3] / rho[1][2], 3.50, 4.80, "monotone P");
  }

  sec_box("Negative tests",
          "d<12.5 Dia-U; beta<0.23 Dia-U; D>760 Dia-F; D>500 Aj-ISA; d=50 Aj-V", 9);
  std::printf("  %s" U_RAQUO " ISO 5167 negative tests%s\n", COLOR_YELLOW, COLOR_RESET);
  {
    BwrConst bt_n = calcBwrMixture(xCH4);
    double rho_n  = calcDensity(20.0, 500.0 / KPA_PER_ATM, bt_n);
    double eta_n  = 11.5e-6;  // CH4 viscosity [Pa*s]
    FlowResult fr_n;
    // Suppress printError messages during negative tests
#ifdef _WIN32
    int saved_fd = _dup(1);
    FILE* nul_f  = std::fopen("NUL", "w");
    if (nul_f) { _dup2(_fileno(nul_f), 1); std::fclose(nul_f); }
#else
    int saved_fd = dup(1);
    FILE* nul_f  = std::fopen("/dev/null", "w");
    if (nul_f) { dup2(fileno(nul_f), 1); std::fclose(nul_f); }
#endif
    double qn1 = calcMassFlow(5.0,500.0,20.0, DeviceType::ORIFICE_CORNER,
                              200.0,170.0, rho_n,eta_n, &fr_n);  // β=0.85>0.80
    double qn2 = calcMassFlow(5.0,500.0,20.0, DeviceType::ORIFICE_CORNER,
                               30.0, 15.0, rho_n,eta_n, &fr_n);  // D=30<50mm
    double qn3 = calcMassFlow(5.0,500.0,20.0, DeviceType::VENTURI_MACHINED,
                              200.0, 60.0, rho_n,eta_n, &fr_n);  // β=0.30<0.40
    double qn4 = calcMassFlow(5.0,500.0,20.0, DeviceType::VENTURI_WELDED_SHEET,
                              200.0,150.0, rho_n,eta_n, &fr_n);  // β=0.75>0.70
    double qn5 = calcMassFlow(5.0,500.0,20.0, DeviceType::ORIFICE_CORNER,
                              200.0, 10.0, rho_n,eta_n, &fr_n);  // d=10<12.5mm
    double qn6 = calcMassFlow(5.0,500.0,20.0, DeviceType::ORIFICE_CORNER,
                              200.0, 44.0, rho_n,eta_n, &fr_n);  // β=0.22<0.23
    double qn7 = calcMassFlow(5.0,500.0,20.0, DeviceType::ORIFICE_FLANGE,
                              800.0,400.0, rho_n,eta_n, &fr_n);  // D=800>760mm
    double qn8 = calcMassFlow(5.0,500.0,20.0, DeviceType::NOZZLE_ISA,
                              550.0,220.0, rho_n,eta_n, &fr_n);  // D=550>500mm
    double qn9 = calcMassFlow(5.0,500.0,20.0, DeviceType::VENTURI_NOZZLE,
                              200.0, 50.0, rho_n,eta_n, &fr_n);  // d=50<=50mm
    std::fflush(stdout);
#ifdef _WIN32
    _dup2(saved_fd, 1);
    _close(saved_fd);
#else
    dup2(saved_fd, 1);
    close(saved_fd);
#endif
    const char* ns = "ISO 5167";
    chk(U_BETA "=0.85 Dia-U (max 0.80) " U_RARR " Qm=0",    qn1, -0.001, 0.001, ns);
    chk("D=30mm Dia-U (min 50mm) " U_RARR " Qm=0",            qn2, -0.001, 0.001, ns);
    chk(U_BETA "=0.30 Ven-P (min 0.40) " U_RARR " Qm=0",    qn3, -0.001, 0.001, ns);
    chk(U_BETA "=0.75 Ven-T (max 0.70) " U_RARR " Qm=0",    qn4, -0.001, 0.001, ns);
    chk("d=10mm Dia-U (min 12.5mm) " U_RARR " Qm=0",          qn5, -0.001, 0.001, ns);
    chk(U_BETA "=0.22 Dia-U (min 0.23) " U_RARR " Qm=0",    qn6, -0.001, 0.001, ns);
    chk("D=800mm Dia-F (max 760mm) " U_RARR " Qm=0",          qn7, -0.001, 0.001, ns);
    chk("D=550mm Aj-ISA (max 500mm) " U_RARR " Qm=0",         qn8, -0.001, 0.001, ns);
    chk("d=50mm Aj-V (limit d>50mm) " U_RARR " Qm=0",      qn9, -0.001, 0.001, ns);
  }

  sec_box("Monotone eta vs T",
          "eta(CH4) increases with T: ratios -10C/20C/60C; sqrt(T) kinetic", 3);
  // Gases: η ∝ √T at dilute limit (kinetic theory). CE includes factor (1+k·ln(T/cs))·√T.
  // NIST CH4: η(-10°C)≈10.5 μPa·s; η(20°C)≈11.0; η(60°C)≈11.9 μPa·s.
  std::printf("  %s" U_RAQUO " Viscosity monotonicity vs T (CH4)%s\n",
              COLOR_YELLOW, COLOR_RESET);
  {
    BwrConst bt_mt = calcBwrMixture(xCH4);
    MixtureThermo thermo_mt = calcMixtureThermo(xCH4, bt_mt);
    double eN10 = calcViscosity(-10.0, calcDensity(-10.0, 1.0, bt_mt), xCH4, thermo_mt.m_rocCrit, thermo_mt.m_csi) * PA_TO_MICRO_PA;
    double eP20 = calcViscosity( 20.0, calcDensity( 20.0, 1.0, bt_mt), xCH4, thermo_mt.m_rocCrit, thermo_mt.m_csi) * PA_TO_MICRO_PA;
    double eP60 = calcViscosity( 60.0, calcDensity( 60.0, 1.0, bt_mt), xCH4, thermo_mt.m_rocCrit, thermo_mt.m_csi) * PA_TO_MICRO_PA;
    chk("Monotone " U_ETA " CH4: " U_ETA "(60" U_DEG"C)/" U_ETA "(-10" U_DEG"C) [-]",
        eP60 / eN10, 1.05, 1.25, "CE/kinetic");
    chk("Monotone " U_ETA " CH4: " U_ETA "(20" U_DEG"C)/" U_ETA "(-10" U_DEG"C) [-]",
        eP20 / eN10, 1.01, 1.10, "CE/kinetic");
    chk("Monotone " U_ETA " CH4: " U_ETA "(60" U_DEG"C)/" U_ETA "(20" U_DEG"C) [-]",
        eP60 / eP20, 1.01, 1.12, "CE/kinetic");
  }

  sec_box("Gas eta ordering",
          "eta(N2)/eta(CH4)~1.61; eta(CO2)/eta(CH4)~1.35 at 20C, 1 atm", 2);
  // CE: η(CH4)≈11.5; η(CO2)≈15.5; η(N2)≈18.5 μPa·s at 20°C, 1 atm.
  std::printf("  %s" U_RAQUO " Viscosity ordering between gases at 20" U_DEG"C%s\n",
              COLOR_YELLOW, COLOR_RESET);
  {
    auto eta20 = [&](const double* xc) -> double {
      BwrConst bt = calcBwrMixture(xc);
      MixtureThermo thermo = calcMixtureThermo(xc, bt);
      double ro20 = calcDensity(20.0, 1.0, bt);
      return calcViscosity(20.0, ro20, xc, thermo.m_rocCrit, thermo.m_csi) * PA_TO_MICRO_PA;
    };
    double eCH4 = eta20(xCH4);
    double eN2  = eta20(xN2);
    double eCO2 = eta20(xCO2);
    chk("Ordering: " U_ETA "(N2)/" U_ETA "(CH4) at 20" U_DEG"C [-]",
        eN2  / eCH4, 1.40, 1.80, "CE/NIST");
    chk("Ordering: " U_ETA "(CO2)/" U_ETA "(CH4) at 20" U_DEG"C [-]",
        eCO2 / eCH4, 1.20, 1.55, "CE/NIST");
  }

  sec_box("Synthetic air",
          "79% N2+21% O2: rho(20C)~1.199 kg/m3; eta~18.2 muPa*s; Qm for 9 types", 15);
  // M_air = 0.79×28.016 + 0.21×32.000 = 28.853 g/mol
  // ρ_ideal(20°C)=28.853/24.055=1.199 kg/m³; ρ_ideal(0°C)=28.853/22.414=1.287
  // η_NIST air at 20°C ≈ 18.2 μPa·s; Z ≈ 1.000 at 1 atm (quasi-ideal)
  {
    double xAer[ARRAY_SIZE] = {};
    xAer[29] = 0.79;  // N2
    xAer[30] = 0.21;  // O2
    test_comp("Synthetic air (79% N2 + 21% O2)", xAer,
              1.182, 1.218,
              1.270, 1.305,
              15.5,  22.0,
              200, 80, 500, 5, false);
  }

  sec_box("H2S density",
          "Pure H2S (Tc=373.6K, Tr=0.785): rho at 1 atm and 500 kPa, N2 consistency", 3);
  // H2S: Tc=373.6 K, Pc=88.9 atm; at 20°C → Tr=0.785 → strongly non-ideal.
  // ρ_ideal(1 atm)=34.082/24.055=1.417 kg/m³; Z≈0.993 → ρ≈1.427
  // ρ_ideal(500 kPa)=4.935×1.417=6.99 kg/m³; Z≈0.965 → ρ≈7.24
  std::printf("  %s" U_RAQUO " H2S: non-ideal density and consistency%s\n",
              COLOR_YELLOW, COLOR_RESET);
  {
    double xH2S[ARRAY_SIZE] = {};  xH2S[26] = 1.0;
    BwrConst bH2S   = calcBwrMixture(xH2S);
    BwrConst bN2_hs = calcBwrMixture(xN2);
    double rH2S_1   = calcDensity(20.0, 1.0,                  bH2S);
    double rH2S_500 = calcDensity(20.0, 500.0 / KPA_PER_ATM,   bH2S);
    double rN2_500  = calcDensity(20.0, 500.0 / KPA_PER_ATM,   bN2_hs);
    chk("Density H2S 20" U_DEG"C, 101.325 kPa [kg/m" U_SUP3 "]",
        rH2S_1,   1.38, 1.48, "BWRS/ideal");
    chk("Density H2S 20" U_DEG"C,  500 kPa   [kg/m" U_SUP3 "]",
        rH2S_500, 6.80, 7.80, "BWRS/NIST");
    chk("Consistency: " U_RHO "(H2S)/" U_RHO "(N2) at 500 kPa [-]",
        rH2S_500 / rN2_500, 1.20, 1.45, "M + Z");
  }

  // ── Summary ───────────────────────────────────────────────────────────────
  std::printf("\n  %sOverall result: %d/%d tests passed.%s\n",
              (npass == ntotal) ? COLOR_BOLD_GREEN : COLOR_BOLD_RED,
              npass, ntotal, COLOR_RESET);
  if (npass < ntotal)
    std::printf("  %sWARNING: Some tests failed " U_MDASH " check the implementation!%s\n",
                COLOR_BOLD_RED, COLOR_RESET);
  std::printf("\n  %sPress any key to continue...%s", COLOR_BOLD_WHITE, COLOR_RESET);
  WaitKey();
  std::printf("\n");
}
