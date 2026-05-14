// BWR.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
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

constexpr int kNumComponents = 36;
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

// ANSI color codes
constexpr const char* kReset      = "\033[0m";
constexpr const char* kBoldYellow = "\033[1;33m";
constexpr const char* kBoldCyan   = "\033[1;36m";
constexpr const char* kBoldGreen  = "\033[1;32m";
constexpr const char* kBoldRed    = "\033[1;31m";
constexpr const char* kCyan       = "\033[36m";
constexpr const char* kYellow     = "\033[33m";
constexpr const char* kMagenta    = "\033[35m";

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
    "  Selectați (1–9)  >";

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
  }
  return "";
}

// Constantele BWR ale amestecului, calculate o singură dată din compoziție.
// Ecuația: p = ρRT + (B₀RT−A₀−C₀/T²)ρ² + (bRT−a)ρ³ + aαρ⁶ + (c/T²)ρ³(1+γρ²)exp(−γρ²)
struct BwrConst {
  double a0         = 0.0;  // A₀ — forțe atractive de ordin 2  (termen −A₀ρ²)
  double b0         = 0.0;  // B₀ — excludere de volum de ordin 2 (termen +B₀RTρ²)
  double c0         = 0.0;  // C₀ — corecție termică de ordin 2  (termen −C₀ρ²/T²)
  double a          = 0.0;  // a  — forțe atractive de ordin 3  (termen −aρ³ și +aαρ⁶)
  double b          = 0.0;  // b  — excludere de volum de ordin 3 (termen +bRTρ³)
  double c          = 0.0;  // c  — corecție termică de ordin 3  (termen +cρ³/T²·exp)
  double alpha      = 0.0;  // α  — amplitudine termen de densitate ρ⁶
  double gamma      = 0.0;  // γ  — parametru Gaussian în exp(−γρ²)
  double molar_mass = 0.0;  // M  — masa molară a amestecului [g/mol]
};

// Rezultatele hidraulice returnate de CalcMassFlow()
struct FlowResult {
  double viteza    = 0.0;  // viteza medie a gazului în conductă [m/s]
  double pierderea = 0.0;  // pierderea de presiune prin strangulare [kPa]
};

// Condiții de referință volumetrică (un preset = 1 sau 2 perechi T/101.325 kPa)
struct CountryRef {
  const char* tara;      // denumire țară/standard
  int         n;         // număr de condiții (1 sau 2)
  double       t[2];      // temperaturi de referință [°C]
  const char* label[2];  // etichete unitate (ex. "Nm³/h", "Sm³/h")
};

// Returnează true dacă tipul este o diafragmă (prize unghi, flanșă sau D-D/2).
// Folosit pentru a alege formula factorului de expansibilitate ε.
constexpr bool IsDiaphragm(TipDispozitiv tip) {
  return static_cast<int>(tip) <=
         static_cast<int>(TipDispozitiv::kDiafragmaDD2);
}

// Calculează coeficientul de debit C [-] al dispozitivului de strangulare.
// Pentru diafragme: ecuația Reader-Harris/Gallagher (ISO 5167-2:2003, §8.3.2).
// d_m  — diametrul interior al conductei la locul de măsurare [m]
// beta — raportul de strangulare β = d/D [-]
// re   — numărul Reynolds în conductă (valoarea curentă din iterație) [-]
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
        L1 = 0.0254 / d_m;  L2p = L1;       // prize la flanșă: 25,4 mm / D
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
        coef += 0.011 * (0.75 - beta) * (2.8 - d_m / 0.0254);
      break;
    }
    case TipDispozitiv::kAjutajIsa:
      coef = 0.99 - 0.2262 * std::pow(beta, 4.1)
           + (0.000215 - 0.001125 * beta + 0.00249 * std::pow(beta, 4.7))
           * std::pow(kReynoldsIsoRef / re, 1.15);
      break;
    case TipDispozitiv::kAjutajRazaLunga:
      coef = 0.9965
           - 0.00653 * std::pow(beta, 0.5)
           * std::pow(kReynoldsIsoRef / re, 0.5);
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
  }
  return coef;
}

// Calculează coeficientul de viteză α = C / √(1 − β⁴) [-].
// Înglobează atât coeficientul de debit C cât și factorul geometric 1/√(1−β⁴)
// pentru a obține direct factorul de amplitudine din ecuația debitului masic.
double VelocityCoefficient(TipDispozitiv tip, double d_m, double beta, double re) {
  return std::pow(1 - std::pow(beta, 4), -0.5)
       * DischargeCoefficient(tip, d_m, beta, re);
}

// Afișează pe consolă mesajul de eroare corespunzător codului de eroare.
// Parametrul red este folosit doar pentru ErrorCode::kReynolds, pentru a
// indica valoarea numerică a lui Re care a depășit domeniul ISO 5167.
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
  }
  std::printf("%s", kReset);
}

// Calculează debitul masic Qm [kg/s] prin dispozitivul de strangulare.
// Înainte de calcul validează domeniile ISO 5167 pentru D, d, β și Re.
// Algoritmul iterează corecția cu Re până la |Qm_k − Qm_{k-1}| < kReynoldsTolerance.
// Returnează 0.0 și afișează eroarea dacă vreun parametru depășește domeniul.
// dp     — presiunea diferențială Δp [kPa]
// p      — presiunea absolută [kPa]
// t      — temperatura fluidului [°C]
// tip    — tipul dispozitivului de strangulare
// d_int  — diametrul interior al conductei la temperatura de referință 20 °C [mm]
// d_orif — diametrul orificiului la temperatura de referință 20 °C [mm]
// ro     — densitatea gazului la condiții (t, p) [kg/m³]
// eta    — viscozitatea dinamică la condiții (t, p) [Pa·s]
// out    — ieșire: viteza medie [m/s] și pierderea de presiune [kPa]
double CalcMassFlow(double dp, double p, double t,
            TipDispozitiv tip, double d_int, double d_orif,
            double ro, double eta, FlowResult* out) {
  double d_i = d_int  * (1 + 0.0000122 * (t - kRefTempCelsius));
  double d_o = d_orif * (1 + 0.0000165 * (t - kRefTempCelsius));

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

  double eps, y;
  if (IsDiaphragm(tip)) {
    eps = 1 - (0.41 + 0.35 * std::pow(beta, 4)) * dp / p / 1.31;
  } else {
    y   = 1 - dp / p;
    eps = std::sqrt(1.31 * std::pow(y, 1.52671) / 0.31
        * (1 - std::pow(beta, 4))
        / (1 - std::pow(beta, 4) * std::pow(y, 1.52671))
        * (1 - std::pow(y, 0.236641)) / (1 - y));
  }

  d_i /= 1000;
  d_o /= 1000;

  double red = kInitialReynolds, qn = 0, q0 = 0, alfa = 0;
  do {
    q0   = qn;
    alfa = VelocityCoefficient(tip, d_i, beta, red);
    qn   = alfa * eps * kPi / 4 * std::pow(d_o, 2) * std::sqrt(2000 * dp * ro);
    red  = 4 * qn / (d_i * kPi * eta);  // Re = 4*Qm / (pi*D*mu)
  } while (std::fabs(qn - q0) > kReynoldsTolerance);

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
  }

  if (!reynolds_valid) {
    PrintError(ErrorCode::kReynolds, red);
    return 0;
  }

  out->viteza    = 4 * qn / kPi / d_i / d_i / ro;
  out->pierderea = (1 - alfa * beta * beta) / (1 + alfa * beta * beta) * dp;
  return qn;
}

// Calculează densitatea amestecului ρ [kg/m³] la temperatura t [°C]
// și presiunea p [atm] prin rezolvarea ecuației de stare BWR cu metoda bisecției.
// Caută ρ în intervalul [kRoMin, kRoMax] [mol/L] până când
// |p_BWR(ρ) − p| < kDensityTolerance [atm].
double CalcDensity(double t, double p, const BwrConst& bwr) {
  double T   = t + kKelvinOffset;
  double R   = kGasConstantR;
  double ro1 = kRoMin;
  double ro2 = kRoMax;
  double ro  = 0.0;
  double pcal = 0.0;
  do {
    ro   = (ro1 + ro2) / 2;
    pcal = R * T * ro
         + (bwr.b0 * R * T - bwr.a0 - bwr.c0 / T / T) * ro * ro
         + (bwr.b  * R * T - bwr.a) * std::pow(ro, 3)
         + bwr.a * bwr.alpha * std::pow(ro, 6)
         + (bwr.c / T / T) * std::pow(ro, 3)
         * (1 + bwr.gamma * ro * ro)
         * std::exp(-bwr.gamma * ro * ro);
    if (pcal > p) ro2 = ro;
    else          ro1 = ro;
  } while (std::fabs(p - pcal) >= kDensityTolerance);
  return bwr.molar_mass * ro;
}

}  // namespace
