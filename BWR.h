// BWR.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace {


constexpr double kPi = 3.141592653589;

enum class DeviceType {
  kOrificeCorner   = 1,
  kOrificeFlange  = 2,
  kOrificeDD2     = 3,
  kNozzleIsa        = 4,
  kNozzleLongRadius  = 5,
  kVenturiRoughCast      = 6,
  kVenturiMachined = 7,
  kVenturiWeldedSheet     = 8,
  kVenturiNozzle    = 9,
};

enum class ErrorCode {
  kPipeDiameter   = 1,
  kOrificeDiameter = 2,
  kThrottleRatio  = 3,
  kReynolds           = 4,
  kNumeric            = 5,
};

constexpr int kNumComponents = 35;
constexpr int kArraySize     = kNumComponents + 1;

constexpr double kKelvinOffset      = 273.15;
constexpr double kGasConstantR      = 0.082057366;  // CODATA 2018
constexpr double kRefTempCelsius    = 20.0;
constexpr double kKpaPerAtm         = 101.325;
constexpr double kReynoldsIsoRef    = 1.0e6;  // ISO 5167 Re normalization
constexpr double kInitialReynolds   = 1.0e6;  // starting guess for Re iteration
constexpr double kReynoldsTolerance = 1.0e-6;  // relative convergence |ΔQm/Qm|
constexpr double kDensityTolerance  = 5.0e-4;
constexpr double kRoMin             = 1.0e-4;
constexpr double kRoMax             = 25.0;
constexpr int    kMaxReynoldsIter   = 100;
constexpr int    kMaxDensityIter    = 200;
constexpr double kSumTolerance      = 1.0e-3;  // tolerance for Σx = 1 check

// Thermal expansion coefficients (linear, per °C, at 20 °C reference)
constexpr double kThermalExpPipe    = 12.2e-6;  // carbon steel pipe (D), ISO 5167
constexpr double kThermalExpOrifice = 16.5e-6;  // stainless steel orifice plate (d)

// ISO 5167 flanged-tap offset [m] = 25.4 mm = 1 inch
constexpr double kFlangeTapM = 0.0254;

// Factor in Qm formula: 2 × (Pa/kPa) so Δp[kPa]·ρ enters as SI
constexpr double k2kPaFactor = 2000.0;

// Isentropic exponent (default air/typical gas mixture; ISO 5167 ε formula)
constexpr double kKappa = 1.31;

// BWRS quadratic mixing-rule coefficient (Starling 1973)
constexpr double kBwrsMixCoef = 8.0;

// Chapman-Enskog polar-correction parameter (viscosity dilute-gas term)
constexpr double kChapEnskog = 0.323;

// High-density viscosity correction (Lucas / Chung-Lee-Starling)
constexpr double kViscHighA    = 10.8e-8;  // prefactor
constexpr double kViscHighExp1 = 1.439;    // positive-exponential argument coefficient
constexpr double kViscHighExp2 = 1.111;    // negative-exponential argument coefficient
constexpr double kViscHighPow  = 1.358;    // outer power

// ANSI color codes
constexpr const char* kReset      = "\033[0;37m";
constexpr const char* kBoldYellow = "\033[1;33m";
constexpr const char* kBoldWhite  = "\033[0;37m";
constexpr const char* kBoldGreen  = "\033[0;32m";
constexpr const char* kBoldRed    = "\033[0;31m";
constexpr const char* kYellow     = "\033[0;33m";
constexpr const char* kCyan       = "\033[0;36m";

// Header block colors (no bold)
constexpr const char* kHdrYellow  = "\033[0;33m";
constexpr const char* kHdrGreen   = "\033[0;32m";
constexpr const char* kHdrCyan    = "\033[0;36m";

constexpr int kTipMin =
    static_cast<int>(DeviceType::kOrificeCorner);
constexpr int kTipMax =
    static_cast<int>(DeviceType::kVenturiNozzle);

const char kTipDisp[] =
    "\n  Throttling device:\n\n"
    "  1.  Orifice plate \xe2\x80\x94 corner taps\n"
    "  2.  Orifice plate \xe2\x80\x94 flange taps\n"
    "  3.  Orifice plate \xe2\x80\x94 D and D/2 taps\n"
    "  4.  ISA 1932 nozzle\n"
    "  5.  Long-radius nozzle\n"
    "  6.  Classical Venturi tube \xe2\x80\x94 rough-cast convergent\n"
    "  7.  Classical Venturi tube \xe2\x80\x94 machined convergent\n"
    "  8.  Classical Venturi tube \xe2\x80\x94 rough-welded sheet-metal convergent\n"
    "  9.  Venturi nozzle\n\n"
    "  Select (1\xe2\x80\x93" "9) >";

// Returns display name for a throttling device type; "" for unknown.
constexpr const char* TipName(DeviceType tip) {
  switch (tip) {
    case DeviceType::kOrificeCorner:   return "Orifice plate \xe2\x80\x94 corner taps";
    case DeviceType::kOrificeFlange:  return "Orifice plate \xe2\x80\x94 flange taps";
    case DeviceType::kOrificeDD2:     return "Orifice plate \xe2\x80\x94 D and D/2 taps";
    case DeviceType::kNozzleIsa:        return "ISA 1932 nozzle";
    case DeviceType::kNozzleLongRadius:  return "Long-radius nozzle";
    case DeviceType::kVenturiRoughCast:      return "Classical Venturi tube \xe2\x80\x94 rough-cast convergent";
    case DeviceType::kVenturiMachined: return "Classical Venturi tube \xe2\x80\x94 machined convergent";
    case DeviceType::kVenturiWeldedSheet:     return "Classical Venturi tube \xe2\x80\x94 rough-welded sheet-metal convergent";
    case DeviceType::kVenturiNozzle:    return "Venturi nozzle";
    default:                               return "";
  }
}

// BWRS mixture constants, computed once from composition.
// Starling (1973) equation of state:
//   p = ρRT + (B₀RT−A₀−C₀/T²+D₀/T³−E₀/T⁴)ρ²
//           + (bRT−a−d/T)ρ³ + α(a+d/T)ρ⁶ + (c/T²)ρ³(1+γρ²)exp(−γρ²)
struct BwrConst {
  double a0         = 0.0;  // A₀ — 2nd-order attractive forces
  double b0         = 0.0;  // B₀ — 2nd-order volume exclusion
  double c0         = 0.0;  // C₀ — thermal correction 1/T²
  double d0         = 0.0;  // D₀ — thermal correction 1/T³ (Starling)
  double e0         = 0.0;  // E₀ — thermal correction 1/T⁴ (Starling)
  double a          = 0.0;  // a  — 3rd-order attractive forces
  double b          = 0.0;  // b  — 3rd-order volume exclusion
  double c          = 0.0;  // c  — 3rd-order thermal correction
  double dv         = 0.0;  // d  — 1/T correction at third virial (Starling)
  double alpha      = 0.0;  // α  — amplitude of ρ⁶ term
  double gamma      = 0.0;  // γ  — Gaussian parameter in exp(−γρ²)
  double molar_mass = 0.0;  // M  — mixture molar mass [g/mol]
};

// Hydraulic results returned by CalcMassFlow()
struct FlowResult {
  double velocity      = 0.0;  // mean gas velocity in pipe [m/s]
  double pressure_loss = 0.0;  // permanent pressure loss through throttling [kPa]
  double beta          = 0.0;  // throttling ratio d/D at measurement temperature [-]
  double reynolds      = 0.0;  // Reynolds number after iteration [-]
  double coef_c        = 0.0;  // discharge coefficient C [-]
  double epsilon       = 0.0;  // expansibility factor ε [-]
  double d_working_mm  = 0.0;  // orifice diameter at working temperature [mm]
  double D_working_mm  = 0.0;  // pipe diameter at working temperature [mm]
  double kappa         = 1.31; // isentropic exponent κ used for ε [-]
};

// Volumetric reference conditions (one preset = 1 or 2 T/101.325 kPa pairs)
struct CountryRef {
  const char* name;      // country/standard name
  int         n;         // number of conditions (1 or 2)
  double      t[2];      // reference temperatures [°C]
  const char* label[2];  // unit labels (e.g. "Nm³/h", "Sm³/h")
};

// Returns true if tip is an orifice plate (corner / flange / D-D/2 taps).
constexpr bool IsDiaphragm(DeviceType tip) {
  return static_cast<int>(tip) <=
         static_cast<int>(DeviceType::kOrificeDD2);
}

// Discharge coefficient C per ISO 5167-2/3/4 (Reader-Harris/Gallagher for orifice plates).
double DischargeCoefficient(DeviceType tip, double d_m, double beta, double re) {
  double coef = 0.0;
  switch (tip) {
    case DeviceType::kOrificeCorner:
    case DeviceType::kOrificeFlange:
    case DeviceType::kOrificeDD2: {
      // Tap positions: L1 upstream, L2p downstream (normalised by D)
      double L1, L2p;
      if (tip == DeviceType::kOrificeCorner) {
        L1 = 0.0;    L2p = 0.0;            // corner taps
      } else if (tip == DeviceType::kOrificeFlange) {
        L1 = kFlangeTapM / d_m;  L2p = L1;   // flange taps: 25.4 mm / D
      } else {
        L1 = 1.0;    L2p = 0.50;           // D and D/2 taps (ISO 5167-2)
      }
      double A  = std::pow(19000.0 * beta / re, 0.8);
      double M2 = 2.0 * L2p / (1.0 - beta);
      coef = 0.5961
           + 0.0261 * beta * beta
           - 0.216  * std::pow(beta, 8.0)
           + 0.000521 * std::pow(1.0e6 * beta / re, 0.7)
           + (0.0188 + 0.0063 * A)
             * std::pow(beta, 3.5) * std::pow(1.0e6 / re, 0.3)
           + (0.043 + 0.080 * std::exp(-10.0 * L1)
                     - 0.123 * std::exp( -7.0 * L1))
             * (1.0 - 0.11 * A) * std::pow(beta, 4.0)
             / (1.0 - std::pow(beta, 4.0))
           - 0.031 * (M2 - 0.8 * std::pow(M2, 1.1))
             * std::pow(beta, 1.3);
      if (d_m < 0.07112)  // correction for D < 71.12 mm
        coef += 0.011 * (0.75 - beta) * (2.8 - d_m / kFlangeTapM);
      break;
    }
    case DeviceType::kNozzleIsa:
      coef = 0.99 - 0.2262 * std::pow(beta, 4.1)
           + (0.000215 - 0.001125 * beta + 0.00249 * std::pow(beta, 4.7))
           * std::pow(kReynoldsIsoRef / re, 1.15);
      break;
    case DeviceType::kNozzleLongRadius:
      coef = 0.9965
           - 0.00653 * std::sqrt(beta)
           * std::sqrt(kReynoldsIsoRef / re);
      break;
    case DeviceType::kVenturiRoughCast:
      coef = 0.984;
      break;
    case DeviceType::kVenturiMachined:
      coef = 0.995;
      break;
    case DeviceType::kVenturiWeldedSheet:
      coef = 0.985;
      break;
    case DeviceType::kVenturiNozzle:
      coef = 0.9858 - 0.196 * std::pow(beta, 4.5);
      break;
    default:
      break;
  }
  return coef;
}

// Velocity coefficient α = C / √(1 − β⁴) used in the mass-flow equation.
double VelocityCoefficient(DeviceType tip, double d_m, double beta, double re) {
  return DischargeCoefficient(tip, d_m, beta, re)
       / std::sqrt(1 - std::pow(beta, 4));
}

// Prints the ISO 5167 error message for the given error code.
void PrintError(ErrorCode code, double red) {
  std::printf("%s", kBoldRed);
  switch (code) {
    case ErrorCode::kPipeDiameter:
      std::printf("\n  Pipe inner diameter is out of range\n");
      break;
    case ErrorCode::kOrificeDiameter:
      std::printf("\n  Throttling device orifice diameter is out of range\n");
      break;
    case ErrorCode::kThrottleRatio:
      std::printf("\n  Throttling ratio \xce\xb2 is out of range\n");
      break;
    case ErrorCode::kReynolds:
      std::printf("\n  Reynolds number (%g) is out of range\n", red);
      break;
    case ErrorCode::kNumeric:
      std::printf("\n  Numerical error — check that inputs are in a physically valid range\n");
      break;
    default:
      break;
  }
  std::printf("%s", kReset);
}

// Calculates mass flow rate Qm [kg/s] via ISO 5167 Reynolds iteration.
// Returns 0.0 on range or convergence error; fills *out with hydraulic results.
double CalcMassFlow(double dp, double p, double t,
            DeviceType tip, double d_int, double d_orif,
            double ro, double eta, FlowResult* out, double kappa = kKappa, int* iters = nullptr) {
  if (!out || !std::isfinite(dp) || !std::isfinite(p) || !std::isfinite(t)
      || !std::isfinite(d_int) || !std::isfinite(d_orif)
      || !std::isfinite(ro) || !std::isfinite(eta)
      || dp <= 0.0 || p <= 0.0 || dp >= p || ro <= 0.0 || eta <= 0.0) {
    PrintError(ErrorCode::kNumeric, 0);
    return 0;
  }

  double d_i = d_int  * (1 + kThermalExpPipe    * (t - kRefTempCelsius));
  double d_o = d_orif * (1 + kThermalExpOrifice * (t - kRefTempCelsius));

  if ((d_i < 50)
      || (d_i > 1000 && tip == DeviceType::kOrificeCorner)
      || (d_i > 760  && (tip == DeviceType::kOrificeFlange
                         || tip == DeviceType::kOrificeDD2))
      || (d_i > 500  && tip == DeviceType::kNozzleIsa)
      || (d_i > 630  && tip == DeviceType::kNozzleLongRadius)) {
    PrintError(ErrorCode::kPipeDiameter, 0);
    return 0;
  }

  if ((d_i > 800  && tip == DeviceType::kVenturiRoughCast)
      || (d_i < 100 && tip == DeviceType::kVenturiRoughCast)
      || (d_i > 250 && tip == DeviceType::kVenturiMachined)
      || (d_i < 50  && tip == DeviceType::kVenturiMachined)
      || (d_i > 1200 && tip == DeviceType::kVenturiWeldedSheet)
      || (d_i < 200  && tip == DeviceType::kVenturiWeldedSheet)
      || (d_i > 500  && tip == DeviceType::kVenturiNozzle)
      || (d_i < 65   && tip == DeviceType::kVenturiNozzle)) {
    PrintError(ErrorCode::kPipeDiameter, 0);
    return 0;
  }

  if ((IsDiaphragm(tip) && d_o < 12.5)
      || (tip == DeviceType::kVenturiNozzle && d_o <= 50)) {
    PrintError(ErrorCode::kOrificeDiameter, 0);
    return 0;
  }

  double beta = d_o / d_i;

  if (((beta < 0.23  || beta > 0.8)  && tip == DeviceType::kOrificeCorner)
      || ((beta < 0.20 || beta > 0.75) && (tip == DeviceType::kOrificeFlange
                                             || tip == DeviceType::kOrificeDD2))
      || ((beta < 0.3  || beta > 0.8) && tip == DeviceType::kNozzleIsa)
      || ((beta < 0.2  || beta > 0.8) && tip == DeviceType::kNozzleLongRadius)) {
    PrintError(ErrorCode::kThrottleRatio, 0);
    return 0;
  }
  if (((beta < 0.3   || beta > 0.75)  && tip == DeviceType::kVenturiRoughCast)
      || ((beta < 0.4   || beta > 0.75) && tip == DeviceType::kVenturiMachined)
      || ((beta < 0.4   || beta > 0.7)  && tip == DeviceType::kVenturiWeldedSheet)
      || ((beta < 0.316 || beta > 0.775) && tip == DeviceType::kVenturiNozzle)) {
    PrintError(ErrorCode::kThrottleRatio, 0);
    return 0;
  }

  double eps;
  if (IsDiaphragm(tip)) {
    eps = 1 - (0.41 + 0.35 * std::pow(beta, 4)) * dp / p / kappa;
  } else {
    double y     = 1 - dp / p;
    double exp2k = 2.0 / kappa;              // exponent y^(2/κ)
    double expk1 = (kappa - 1.0) / kappa;   // exponent y^((κ-1)/κ)
    double b4    = std::pow(beta, 4);
    double yk    = std::pow(y, exp2k);
    eps = std::sqrt(kappa * yk / (kappa - 1.0)
        * (1.0 - b4) / (1.0 - b4 * yk)
        * (1.0 - std::pow(y, expk1)) / (1.0 - y));
  }
  if (!std::isfinite(eps) || eps <= 0.0) {
    PrintError(ErrorCode::kNumeric, 0);
    return 0;
  }

  d_i /= 1000;
  d_o /= 1000;

  double red  = kInitialReynolds;
  double qn   = 0;
  double q0   = 0;
  double alfa = 0;
  int nre = 0;
  bool converged = false;
  for (; nre < kMaxReynoldsIter; ) {
    nre++;
    q0   = qn;
    alfa = VelocityCoefficient(tip, d_i, beta, red);
    qn   = alfa * eps * kPi / 4 * std::pow(d_o, 2) * std::sqrt(k2kPaFactor * dp * ro);
    red  = 4 * qn / (d_i * kPi * eta);  // Re = 4*Qm / (pi*D*mu)
    if (!std::isfinite(alfa) || !std::isfinite(qn) || !std::isfinite(red)) {
      PrintError(ErrorCode::kNumeric, 0);
      return 0;
    }
    if (std::fabs((qn - q0) / std::max(qn, 1e-10)) <= kReynoldsTolerance) {
      converged = true;
      break;
    }
  }
  if (iters) *iters = nre;
  if (!converged) {
    PrintError(ErrorCode::kNumeric, 0);
    return 0;
  }

  bool reynolds_valid = false;
  switch (tip) {
    case DeviceType::kOrificeCorner:
      if ((5000  <= red) && (red <= 1e8) && (0.23 <= beta) && (beta < 0.45))
        reynolds_valid = true;
      if ((10000 <= red) && (red <= 1e8) && (0.45 <= beta) && (beta < 0.77))
        reynolds_valid = true;
      if ((20000 <= red) && (red <= 1e8) && (0.77 <= beta) && (beta <= 0.80))
        reynolds_valid = true;
      break;
    case DeviceType::kOrificeFlange:
    case DeviceType::kOrificeDD2:
      if ((1.26e6 * beta * beta * d_i <= red) && (red <= 1e8))
        reynolds_valid = true;
      break;
    case DeviceType::kNozzleIsa:
      if ((70000 <= red) && (red <= 1e7) && (0.30 <= beta) && (beta < 0.44))
        reynolds_valid = true;
      if ((20000 <= red) && (red <= 1e7) && (0.44 <= beta) && (beta <= 0.80))
        reynolds_valid = true;
      break;
    case DeviceType::kNozzleLongRadius:
      if ((10000  <= red) && (red <= 1e7))
        reynolds_valid = true;
      break;
    case DeviceType::kVenturiRoughCast:
      if ((200000 <= red) && (red <= 2000000))
        reynolds_valid = true;
      break;
    case DeviceType::kVenturiMachined:
      if ((200000 <= red) && (red <= 1000000))
        reynolds_valid = true;
      break;
    case DeviceType::kVenturiWeldedSheet:
      if ((200000 <= red) && (red <= 2000000))
        reynolds_valid = true;
      break;
    case DeviceType::kVenturiNozzle:
      if ((150000 <= red) && (red <= 2000000))
        reynolds_valid = true;
      break;
    default:
      break;
  }

  if (!reynolds_valid) {
    PrintError(ErrorCode::kReynolds, red);
    return 0;
  }

  out->velocity      = 4 * qn / kPi / d_i / d_i / ro;
  // Exact ISO Annex A for diaphragms/nozzles; empirical fraction for classical Venturi
  // (Venturi recovers 80-95% of Δp; orifice formula would give ~60% — far too high)
  double loss_frac;
  if (tip == DeviceType::kVenturiRoughCast) {
    loss_frac = 0.15;
  } else if (tip == DeviceType::kVenturiMachined) {
    loss_frac = 0.08;
  } else if (tip == DeviceType::kVenturiWeldedSheet) {
    loss_frac = 0.15;
  } else {
    double C  = out->coef_c;
    double b2 = beta * beta;
    double sq = std::sqrt(1.0 - b2 * b2 * (1.0 - C * C));
    loss_frac = (sq - C * b2) / (sq + C * b2);
  }
  out->pressure_loss = loss_frac * dp;
  out->beta        = beta;
  out->reynolds    = red;
  out->coef_c      = DischargeCoefficient(tip, d_i, beta, red);
  out->epsilon     = eps;
  out->d_working_mm  = d_o * 1000.0;
  out->D_working_mm  = d_i * 1000.0;
  out->kappa       = kappa;
  return qn;
}

// Solves the BWRS equation for mixture density ρ [kg/m³] via bisection on the
// gas-phase branch. Adaptive upper bracket avoids the van der Waals loop
// (gases below Tc, e.g. H2S, CO2, NH3). Converges when |p_calc − p| < kDensityTolerance.
double CalcDensity(double t, double p, const BwrConst& bwr, int* iters = nullptr) {
  double T   = t + kKelvinOffset;
  double R   = kGasConstantR;
  if (!std::isfinite(T) || !std::isfinite(p) || !std::isfinite(bwr.molar_mass)
      || T <= 0.0 || p <= 0.0 || bwr.molar_mass <= 0.0) {
    if (iters) *iters = 0;
    return 0.0;
  }
  auto pressure_at = [&](double ro) {
    return R * T * ro
         + (bwr.b0 * R * T - bwr.a0 - bwr.c0 / T / T
            + bwr.d0 / (T * T * T) - bwr.e0 / (T * T * T * T)) * ro * ro
         + (bwr.b  * R * T - bwr.a - bwr.dv / T) * std::pow(ro, 3)
         + bwr.alpha * (bwr.a + bwr.dv / T) * std::pow(ro, 6)
         + (bwr.c / T / T) * std::pow(ro, 3)
         * (1 + bwr.gamma * ro * ro)
         * std::exp(-bwr.gamma * ro * ro);
  };

  // Adaptive bracket: start from ideal-gas estimate and double ro2 until
  // p(ro2) >= p_target, stopping if pressure starts decreasing (van der Waals
  // spinodal) to stay on the gas-phase branch.
  double ro1 = kRoMin;
  double p1  = pressure_at(ro1);
  double ro2 = std::max(p / (R * T), ro1 * 2.0);  // ideal-gas starting point
  double p2  = pressure_at(ro2);
  double p_lo = p1;
  for (int k = 0; k < 60; k++) {
    if (!std::isfinite(p2) || p2 < p_lo) break;  // spinodal or non-finite
    if (p2 >= p) break;                           // upper bracket found
    p_lo = p2;
    ro1  = ro2;  p1 = p2;
    ro2  = std::min(ro2 * 2.0, kRoMax);
    p2   = pressure_at(ro2);
  }

  if (!std::isfinite(p1) || !std::isfinite(p2) || p < p1 || p > p2) {
    if (iters) *iters = 0;
    return 0.0;
  }
  double ro   = 0.0;
  double pcal = 0.0;
  int n = 0;
  for (; n < kMaxDensityIter; ) {
    n++;
    ro   = (ro1 + ro2) / 2;
    pcal = pressure_at(ro);
    if (!std::isfinite(pcal)) {
      if (iters) *iters = n;
      return 0.0;
    }
    if (std::fabs(p - pcal) < kDensityTolerance) {
      if (iters) *iters = n;
      return bwr.molar_mass * ro;
    }
    if (pcal > p) ro2 = ro;
    else          ro1 = ro;
  }
  if (iters) *iters = n;
  return 0.0;
}

}  // namespace
