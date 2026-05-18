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

enum class TipDispozitiv {
  kDiafragmaUnghi   = 1,
  kDiafragmaFlansa  = 2,
  kDiafragmaDD2     = 3,
  kAjutajIsa        = 4,
  kAjutajRazaLunga  = 5,
  kVenturiBrut      = 6,
  kVenturiPrelucrat = 7,
  kVenturiTabla     = 8,
  kAjutajVenturi    = 9,
};

enum class ErrorCode {
  kDiametruInterior   = 1,
  kOrificuStrangulare = 2,
  kRaportStrangulare  = 3,
  kReynolds           = 4,
};

constexpr int kNumComponents = 35;
constexpr int kArraySize     = kNumComponents + 1;

constexpr double kKelvinOffset      = 273.15;
constexpr double kGasConstantR      = 0.082055;
constexpr double kRefTempCelsius    = 20.0;
constexpr double kKpaPerAtm         = 101.325;
constexpr double kReynoldsIsoRef    = 1.0e6;  // ISO 5167 Re normalization
constexpr double kInitialReynolds   = 1.0e6;  // starting guess for Re iteration
constexpr double kReynoldsTolerance = 1.0e-4;
constexpr double kDensityTolerance  = 5.0e-4;
constexpr double kRoMin             = 1.0e-4;
constexpr double kRoMax             = 5.0;
constexpr double kSumTolerance      = 1.0e-3;  // tolerance for Σx = 1 check

// Thermal expansion coefficients (linear, per °C, at 20 °C reference)
constexpr double kThermalExpPipe    = 12.2e-6;  // carbon steel pipe (D), ISO 5167
constexpr double kThermalExpOrifice = 16.5e-6;  // stainless steel orifice plate (d)

// ISO 5167 flanged-tap offset [m] = 25.4 mm = 1 inch
constexpr double kFlangeTapM = 0.0254;

// Factor in Qm formula: 2 × (Pa/kPa) so Δp[kPa]·ρ enters as SI
constexpr double k2kPaFactor = 2000.0;

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
constexpr const char* kBoldYellow = "\033[33m";
constexpr const char* kBoldWhite  = "\033[37m";
constexpr const char* kBoldGreen  = "\033[32m";
constexpr const char* kBoldRed    = "\033[31m";
constexpr const char* kYellow     = "\033[33m";
constexpr const char* kCyan       = "\033[36m";

// Bold — exclusiv pentru blocul antet
constexpr const char* kHdrYellow  = "\033[1;33m";
constexpr const char* kHdrGreen   = "\033[1;32m";
constexpr const char* kHdrCyan    = "\033[1;36m";

constexpr int kTipMin =
    static_cast<int>(TipDispozitiv::kDiafragmaUnghi);
constexpr int kTipMax =
    static_cast<int>(TipDispozitiv::kAjutajVenturi);

const char kTipDisp[] =
    "\n  Dispozitiv de strangulare:\n\n"
    "  1.  Diafragmă cu prize în unghi\n"
    "  2.  Diafragmă cu prize la flanșă\n"
    "  3.  Diafragmă cu prize la D și D/2\n"
    "  4.  Ajutaj ISA 1932\n"
    "  5.  Ajutaj cu rază lungă\n"
    "  6.  Tub Venturi clasic — convergent brut turnat\n"
    "  7.  Tub Venturi clasic — convergent prelucrat\n"
    "  8.  Tub Venturi clasic — convergent brut din tablă sudată\n"
    "  9.  Ajutaj Venturi\n\n"
    "  Selectați (1–9) >";

// Funcție: TipName
// Intrări: tip — codul dispozitivului de strangulare (TipDispozitiv)
// Ieșire:  șir cu denumirea dispozitivului; "" pentru valoare necunoscută
// Scop:    conversie enum → text pentru afișare în consolă
constexpr const char* TipName(TipDispozitiv tip) {
  switch (tip) {
    case TipDispozitiv::kDiafragmaUnghi:   return "Diafragmă cu prize în unghi";
    case TipDispozitiv::kDiafragmaFlansa:  return "Diafragmă cu prize la flanșă";
    case TipDispozitiv::kDiafragmaDD2:     return "Diafragmă cu prize la D și D/2";
    case TipDispozitiv::kAjutajIsa:        return "Ajutaj ISA 1932";
    case TipDispozitiv::kAjutajRazaLunga:  return "Ajutaj cu rază lungă";
    case TipDispozitiv::kVenturiBrut:      return "Tub Venturi clasic — convergent brut turnat";
    case TipDispozitiv::kVenturiPrelucrat: return "Tub Venturi clasic — convergent prelucrat";
    case TipDispozitiv::kVenturiTabla:     return "Tub Venturi clasic — convergent brut din tablă sudată";
    case TipDispozitiv::kAjutajVenturi:    return "Ajutaj Venturi";
    default:                               return "";
  }
}

// Constantele BWRS ale amestecului, calculate o singură dată din compoziție.
// Ecuația Starling (1973):
//   p = ρRT + (B₀RT−A₀−C₀/T²+D₀/T³−E₀/T⁴)ρ²
//           + (bRT−a−d/T)ρ³ + α(a+d/T)ρ⁶ + (c/T²)ρ³(1+γρ²)exp(−γρ²)
struct BwrConst {
  double a0         = 0.0;  // A₀ — forțe atractive de ordin 2
  double b0         = 0.0;  // B₀ — excludere de volum de ordin 2
  double c0         = 0.0;  // C₀ — corecție termică 1/T²
  double d0         = 0.0;  // D₀ — corecție termică 1/T³ (Starling)
  double e0         = 0.0;  // E₀ — corecție termică 1/T⁴ (Starling)
  double a          = 0.0;  // a  — forțe atractive de ordin 3
  double b          = 0.0;  // b  — excludere de volum de ordin 3
  double c          = 0.0;  // c  — corecție termică de ordin 3
  double dv         = 0.0;  // d  — corecție 1/T la virial terț (Starling)
  double alpha      = 0.0;  // α  — amplitudine termen ρ⁶
  double gamma      = 0.0;  // γ  — parametru Gaussian în exp(−γρ²)
  double molar_mass = 0.0;  // M  — masa molară a amestecului [g/mol]
};

// Rezultatele hidraulice returnate de CalcMassFlow()
struct FlowResult {
  double viteza      = 0.0;  // viteza medie a gazului în conductă [m/s]
  double pierderea   = 0.0;  // pierderea de presiune prin strangulare [kPa]
  double beta        = 0.0;  // raportul de strangulare d/D la temperatura de măsurare [-]
  double reynolds    = 0.0;  // numărul Reynolds după iterație [-]
  double coef_c      = 0.0;  // coeficientul de debit C [-]
  double epsilon     = 0.0;  // factorul de expansibilitate ε [-]
  double d_lucru_mm  = 0.0;  // diametrul orificiului la temperatura de lucru [mm]
  double D_lucru_mm  = 0.0;  // diametrul conductei la temperatura de lucru [mm]
};

// Condiții de referință volumetrică (un preset = 1 sau 2 perechi T/101.325 kPa)
struct CountryRef {
  const char* tara;      // denumire țară/standard
  int         n;         // număr de condiții (1 sau 2)
  double      t[2];      // temperaturi de referință [°C]
  const char* label[2];  // etichete unitate (ex. "Nm³/h", "Sm³/h")
};

// Funcție: IsDiaphragm
// Intrări: tip — codul dispozitivului de strangulare
// Ieșire:  true dacă tipul este diafragmă (unghi / flanșă / D-D/2), false altfel
// Scop:    selecția formulei corecte pentru factorul de expansibilitate ε
constexpr bool IsDiaphragm(TipDispozitiv tip) {
  return static_cast<int>(tip) <=
         static_cast<int>(TipDispozitiv::kDiafragmaDD2);
}

// Funcție: DischargeCoefficient
// Intrări: tip — tipul dispozitivului; d_m — diametrul interior la locul de măsurare [m];
//          beta — β = d/D [-]; re — numărul Reynolds curent din iterație [-]
// Ieșire:  coeficientul de debit C [-]
// Scop:    calculul C conform ISO 5167-2/3/4 (ecuația Reader-Harris/Gallagher pt. diafragme)
double DischargeCoefficient(TipDispozitiv tip, double d_m, double beta, double re) {
  double coef = 0.0;
  switch (tip) {
    case TipDispozitiv::kDiafragmaUnghi:
    case TipDispozitiv::kDiafragmaFlansa:
    case TipDispozitiv::kDiafragmaDD2: {
      // Poziția prizelor de presiune: L1 amonte, L2p aval (adimensionalizate cu D)
      double L1, L2p;
      if (tip == TipDispozitiv::kDiafragmaUnghi) {
        L1 = 0.0;    L2p = 0.0;            // prize în unghi
      } else if (tip == TipDispozitiv::kDiafragmaFlansa) {
        L1 = kFlangeTapM / d_m;  L2p = L1;   // prize la flanșă: 25,4 mm / D
      } else {
        L1 = 1.0;    L2p = 0.47;           // prize la D și D/2
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
      if (d_m < 0.07112)  // corecție pentru D < 71,12 mm
        coef += 0.011 * (0.75 - beta) * (2.8 - d_m / kFlangeTapM);
      break;
    }
    case TipDispozitiv::kAjutajIsa:
      coef = 0.99 - 0.2262 * std::pow(beta, 4.1)
           + (0.000215 - 0.001125 * beta + 0.00249 * std::pow(beta, 4.7))
           * std::pow(kReynoldsIsoRef / re, 1.15);
      break;
    case TipDispozitiv::kAjutajRazaLunga:
      coef = 0.9965
           - 0.00653 * std::sqrt(beta)
           * std::sqrt(kReynoldsIsoRef / re);
      break;
    case TipDispozitiv::kVenturiBrut:
      coef = 0.984;
      break;
    case TipDispozitiv::kVenturiPrelucrat:
      coef = 0.995;
      break;
    case TipDispozitiv::kVenturiTabla:
      coef = 0.985;
      break;
    case TipDispozitiv::kAjutajVenturi:
      coef = 0.9858 - 0.196 * std::pow(beta, 4.5);
      break;
    default:
      break;
  }
  return coef;
}

// Funcție: VelocityCoefficient
// Intrări: tip, d_m, beta, re — aceleași ca DischargeCoefficient
// Ieșire:  coeficientul de viteză α = C / √(1 − β⁴) [-]
// Scop:    calculul factorului de amplitudine din ecuația debitului masic
double VelocityCoefficient(TipDispozitiv tip, double d_m, double beta, double re) {
  return DischargeCoefficient(tip, d_m, beta, re)
       / std::sqrt(1 - std::pow(beta, 4));
}

// Funcție: PrintError
// Intrări: code — codul de eroare (ErrorCode); red — valoarea Re (relevantă doar pt. kReynolds)
// Ieșire:  —
// Scop:    afișează pe consolă mesajul de eroare ISO 5167 corespunzător codului
void PrintError(ErrorCode code, double red) {
  std::printf("%s", kBoldRed);
  switch (code) {
    case ErrorCode::kDiametruInterior:
      std::printf(
          "\n  Diametrul interior al conductei este necorespunzător\n");
      break;
    case ErrorCode::kOrificuStrangulare:
      std::printf(
          "\n  Orificiul dispozitivului de strangulare este necorespunzător\n");
      break;
    case ErrorCode::kRaportStrangulare:
      std::printf("\n  Raportul de strangulare este necorespunzător\n");
      break;
    case ErrorCode::kReynolds:
      std::printf("\n  Număr Reynolds (%g) necorespunzător\n", red);
      break;
    default:
      break;
  }
  std::printf("%s", kReset);
}

// Funcție: CalcMassFlow
// Intrări: dp — Δp [kPa]; p — presiune absolută [kPa]; t — temperatură [°C];
//          tip — tipul dispozitivului; d_int — D la 20 °C [mm]; d_orif — d la 20 °C [mm];
//          ro — densitate (t,p) [kg/m³]; eta — vâscozitate dinamică (t,p) [Pa·s]; out — structură ieșire hidraulică
// Ieșire:  debitul masic Qm [kg/s]; câmpurile FlowResult prin *out; 0.0 la eroare
// Scop:    validare domeniu ISO 5167, calcul Qm prin iterație pe Re
double CalcMassFlow(double dp, double p, double t,
            TipDispozitiv tip, double d_int, double d_orif,
            double ro, double eta, FlowResult* out, int* iters = nullptr) {
  double d_i = d_int  * (1 + kThermalExpPipe    * (t - kRefTempCelsius));
  double d_o = d_orif * (1 + kThermalExpOrifice * (t - kRefTempCelsius));

  if ((d_i < 50)
      || (d_i > 1000 && tip == TipDispozitiv::kDiafragmaUnghi)
      || (d_i > 760  && (tip == TipDispozitiv::kDiafragmaFlansa
                         || tip == TipDispozitiv::kDiafragmaDD2))
      || (d_i > 500  && tip == TipDispozitiv::kAjutajIsa)
      || (d_i > 630  && tip == TipDispozitiv::kAjutajRazaLunga)) {
    PrintError(ErrorCode::kDiametruInterior, 0);
    return 0;
  }

  if ((d_i > 800  && tip == TipDispozitiv::kVenturiBrut)
      || (d_i < 100 && tip == TipDispozitiv::kVenturiBrut)
      || (d_i > 250 && tip == TipDispozitiv::kVenturiPrelucrat)
      || (d_i < 50  && tip == TipDispozitiv::kVenturiPrelucrat)
      || (d_i > 1200 && tip == TipDispozitiv::kVenturiTabla)
      || (d_i < 200  && tip == TipDispozitiv::kVenturiTabla)
      || (d_i > 500  && tip == TipDispozitiv::kAjutajVenturi)
      || (d_i < 65   && tip == TipDispozitiv::kAjutajVenturi)) {
    PrintError(ErrorCode::kDiametruInterior, 0);
    return 0;
  }

  if ((IsDiaphragm(tip) && d_o < 12.5)
      || (tip == TipDispozitiv::kAjutajVenturi && d_o <= 50)) {
    PrintError(ErrorCode::kOrificuStrangulare, 0);
    return 0;
  }

  double beta = d_o / d_i;

  if (((beta < 0.23  || beta > 0.8)  && tip == TipDispozitiv::kDiafragmaUnghi)
      || ((beta < 0.20 || beta > 0.75) && (tip == TipDispozitiv::kDiafragmaFlansa
                                             || tip == TipDispozitiv::kDiafragmaDD2))
      || ((beta < 0.3  || beta > 0.8) && tip == TipDispozitiv::kAjutajIsa)
      || ((beta < 0.2  || beta > 0.8) && tip == TipDispozitiv::kAjutajRazaLunga)) {
    PrintError(ErrorCode::kRaportStrangulare, 0);
    return 0;
  }
  if (((beta < 0.3   || beta > 0.75)  && tip == TipDispozitiv::kVenturiBrut)
      || ((beta < 0.4   || beta > 0.75) && tip == TipDispozitiv::kVenturiPrelucrat)
      || ((beta < 0.4   || beta > 0.7)  && tip == TipDispozitiv::kVenturiTabla)
      || ((beta < 0.316 || beta > 0.775) && tip == TipDispozitiv::kAjutajVenturi)) {
    PrintError(ErrorCode::kRaportStrangulare, 0);
    return 0;
  }

  double eps;
  if (IsDiaphragm(tip)) {
    eps = 1 - (0.41 + 0.35 * std::pow(beta, 4)) * dp / p / 1.31;
  } else {
    double y = 1 - dp / p;
    eps = std::sqrt(1.31 * std::pow(y, 1.52671) / 0.31
        * (1 - std::pow(beta, 4))
        / (1 - std::pow(beta, 4) * std::pow(y, 1.52671))
        * (1 - std::pow(y, 0.236641)) / (1 - y));
  }

  d_i /= 1000;
  d_o /= 1000;

  double red  = kInitialReynolds;
  double qn   = 0;
  double q0   = 0;
  double alfa = 0;
  int nre = 0;
  do {
    nre++;
    q0   = qn;
    alfa = VelocityCoefficient(tip, d_i, beta, red);
    qn   = alfa * eps * kPi / 4 * std::pow(d_o, 2) * std::sqrt(k2kPaFactor * dp * ro);
    red  = 4 * qn / (d_i * kPi * eta);  // Re = 4*Qm / (pi*D*mu)
  } while (std::fabs(qn - q0) > kReynoldsTolerance);
  if (iters) *iters = nre;

  bool reynolds_valid = false;
  switch (tip) {
    case TipDispozitiv::kDiafragmaUnghi:
      if ((5000  <= red) && (red <= 1e8) && (0.23 <= beta) && (beta < 0.45))
        reynolds_valid = true;
      if ((10000 <= red) && (red <= 1e8) && (0.45 <= beta) && (beta < 0.77))
        reynolds_valid = true;
      if ((20000 <= red) && (red <= 1e8) && (0.77 <= beta) && (beta <= 0.80))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kDiafragmaFlansa:
    case TipDispozitiv::kDiafragmaDD2:
      if ((1.26e6 * beta * beta * d_i <= red) && (red <= 1e8))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kAjutajIsa:
      if ((70000 <= red) && (red <= 1e7) && (0.30 <= beta) && (beta < 0.44))
        reynolds_valid = true;
      if ((20000 <= red) && (red <= 1e7) && (0.44 <= beta) && (beta <= 0.80))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kAjutajRazaLunga:
      if ((10000  <= red) && (red <= 1e7))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kVenturiBrut:
      if ((200000 <= red) && (red <= 2000000))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kVenturiPrelucrat:
      if ((200000 <= red) && (red <= 1000000))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kVenturiTabla:
      if ((200000 <= red) && (red <= 2000000))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kAjutajVenturi:
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

  out->viteza      = 4 * qn / kPi / d_i / d_i / ro;
  out->pierderea   = (1 - alfa * beta * beta) / (1 + alfa * beta * beta) * dp;
  out->beta        = beta;
  out->reynolds    = red;
  out->coef_c      = DischargeCoefficient(tip, d_i, beta, red);
  out->epsilon     = eps;
  out->d_lucru_mm  = d_o * 1000.0;
  out->D_lucru_mm  = d_i * 1000.0;
  return qn;
}

// Funcție: CalcDensity
// Intrări: t — temperatura [°C]; p — presiunea [atm]; bwr — constantele BWRS ale amestecului
// Ieșire:  densitatea amestecului ρ [kg/m³]
// Scop:    rezolvarea ecuației BWRS prin bisecție în [kRoMin, kRoMax] până la |p_calc − p| < kDensityTolerance
double CalcDensity(double t, double p, const BwrConst& bwr, int* iters = nullptr) {
  double T   = t + kKelvinOffset;
  double R   = kGasConstantR;
  double ro1 = kRoMin;
  double ro2 = kRoMax;
  double ro  = 0.0;
  double pcal = 0.0;
  int n = 0;
  do {
    n++;
    ro   = (ro1 + ro2) / 2;
    pcal = R * T * ro
         + (bwr.b0 * R * T - bwr.a0 - bwr.c0 / T / T
            + bwr.d0 / (T * T * T) - bwr.e0 / (T * T * T * T)) * ro * ro
         + (bwr.b  * R * T - bwr.a - bwr.dv / T) * std::pow(ro, 3)
         + bwr.alpha * (bwr.a + bwr.dv / T) * std::pow(ro, 6)
         + (bwr.c / T / T) * std::pow(ro, 3)
         * (1 + bwr.gamma * ro * ro)
         * std::exp(-bwr.gamma * ro * ro);
    if (pcal > p) ro2 = ro;
    else          ro1 = ro;
  } while (std::fabs(p - pcal) >= kDensityTolerance);
  if (iters) *iters = n;
  return bwr.molar_mass * ro;
}

}  // namespace
