// BWR.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

namespace {


constexpr float kPi = 3.141592653589f;

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

constexpr float kKelvinOffset      = 273.15f;
constexpr float kGasConstantR      = 0.082055f;
constexpr float kRefTempCelsius    = 20.0f;
constexpr float kKpaPerAtm         = 101.325f;
constexpr float kReynoldsIsoRef    = 1.0e6f;  // ISO 5167 Re normalization
constexpr float kInitialReynolds   = 1.0e6f;  // starting guess for Re iteration
constexpr float kReynoldsTolerance = 1.0e-4f;
constexpr float kDensityTolerance  = 5.0e-4f;
constexpr float kRoMin             = 1.0e-4f;
constexpr float kRoMax             = 5.0f;
constexpr float kSumTolerance      = 1.0e-3f;  // tolerance for Σx = 1 check

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
  float a0         = 0.0f;  // A₀ — forțe atractive de ordin 2  (termen −A₀ρ²)
  float b0         = 0.0f;  // B₀ — excludere de volum de ordin 2 (termen +B₀RTρ²)
  float c0         = 0.0f;  // C₀ — corecție termică de ordin 2  (termen −C₀ρ²/T²)
  float a          = 0.0f;  // a  — forțe atractive de ordin 3  (termen −aρ³ și +aαρ⁶)
  float b          = 0.0f;  // b  — excludere de volum de ordin 3 (termen +bRTρ³)
  float c          = 0.0f;  // c  — corecție termică de ordin 3  (termen +cρ³/T²·exp)
  float alpha      = 0.0f;  // α  — amplitudine termen de densitate ρ⁶
  float gamma      = 0.0f;  // γ  — parametru Gaussian în exp(−γρ²)
  float molar_mass = 0.0f;  // M  — masa molară a amestecului [g/mol]
};

// Rezultatele hidraulice returnate de CalcMassFlow()
struct FlowResult {
  float viteza    = 0.0f;  // viteza medie a gazului în conductă [m/s]
  float pierderea = 0.0f;  // pierderea de presiune prin strangulare [kPa]
};

// Returnează true dacă tipul este o diafragmă (prize unghi, flanșă sau D-D/2).
// Folosit pentru a alege formula factorului de expansibilitate ε.
constexpr bool IsDiaphragm(TipDispozitiv tip) {
  return static_cast<int>(tip) <=
         static_cast<int>(TipDispozitiv::kDiafragmaDD2);
}

// Calculează coeficientul de debit C [-] al dispozitivului de strangulare
// conform formulelor ISO 5167 / STAS 7347-90.
// d_m  — diametrul interior al conductei la locul de măsurare [m]
// beta — raportul de strangulare β = d/D [-]
// re   — numărul Reynolds în conductă (valoarea curentă din iterație) [-]
float DischargeCoefficient(TipDispozitiv tip, float d_m, float beta, float re) {
  float coef = 0.0f;
  switch (tip) {
    case TipDispozitiv::kDiafragmaUnghi:
      coef = 0.5959f + 0.0312f * std::pow(beta, 2.1)
           - 0.184f * std::pow(beta, 8)
           + 0.0029f * std::pow(beta, 2.5)
           * std::pow(kReynoldsIsoRef / re, 0.75);
      break;
    case TipDispozitiv::kDiafragmaFlansa:
      if (1000 * d_m > 58.62f) {
        coef = 0.5959f + 0.0312f * std::pow(beta, 2.1)
             - 0.184f * std::pow(beta, 8)
             + 0.002286f / d_m * std::pow(beta, 4)
             / (1 - std::pow(beta, 4));
        coef = coef - 0.00085598f / d_m * std::pow(beta, 3)
             + 0.0029f * std::pow(beta, 2.5)
             * std::pow(kReynoldsIsoRef / re, 0.75);
      } else {
        coef = 0.5959f + 0.0312f * std::pow(beta, 2.1)
             - 0.184f * std::pow(beta, 8)
             + 0.039f * std::pow(beta, 4) / (1 - std::pow(beta, 4));
        coef = coef - 0.039f * std::pow(beta, 3)
             + 0.0029f * std::pow(beta, 2.5)
             * std::pow(kReynoldsIsoRef / re, 0.75);
      }
      break;
    case TipDispozitiv::kDiafragmaDD2:
      coef = 0.5959f + 0.0312f * std::pow(beta, 2.1)
           - 0.184f * std::pow(beta, 8)
           + 0.039f * std::pow(beta, 4) / (1 - std::pow(beta, 4))
           - 0.015839f * std::pow(beta, 3);
      coef = coef + 0.0029f * std::pow(beta, 2.5)
           * std::pow(kReynoldsIsoRef / re, 0.75);
      break;
    case TipDispozitiv::kAjutajIsa:
      coef = 0.99f - 0.2262f * std::pow(beta, 4.1)
           + (0.000215f - 0.001125f * beta + 0.00249f * std::pow(beta, 4.7))
           * std::pow(kReynoldsIsoRef / re, 1.15);
      break;
    case TipDispozitiv::kAjutajRazaLunga:
      coef = 0.9965f
           - 0.00653f * std::pow(beta, 0.5)
           * std::pow(kReynoldsIsoRef / re, 0.5);
      break;
    case TipDispozitiv::kVenturiBrut:
      coef = 0.984f;
      break;
    case TipDispozitiv::kVenturiPrelucrat:
      coef = 0.995f;
      break;
    case TipDispozitiv::kVenturiTabla:
      coef = 0.985f;
      break;
    case TipDispozitiv::kAjutajVenturi:
      coef = 0.9858f - 0.196f * std::pow(beta, 4.5);
      break;
  }
  return coef;
}

// Calculează coeficientul de viteză α = C / √(1 − β⁴) [-].
// Înglobează atât coeficientul de debit C cât și factorul geometric 1/√(1−β⁴)
// pentru a obține direct factorul de amplitudine din ecuația debitului masic.
float VelocityCoefficient(TipDispozitiv tip, float d_m, float beta, float re) {
  return std::pow(1 - std::pow(beta, 4), -0.5)
       * DischargeCoefficient(tip, d_m, beta, re);
}

// Afișează pe consolă mesajul de eroare corespunzător codului de eroare.
// Parametrul red este folosit doar pentru ErrorCode::kReynolds, pentru a
// indica valoarea numerică a lui Re care a depășit domeniul STAS 7347-90.
void PrintError(ErrorCode code, float red) {
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
// Înainte de calcul validează domeniile STAS 7347-90 pentru D, d, β și Re.
// Algoritmul iterează corecția cu Re până la |Qm_k − Qm_{k-1}| < kReynoldsTolerance.
// Returnează 0.0f și afișează eroarea dacă vreun parametru depășește domeniul.
// dp     — presiunea diferențială Δp [kPa]
// p      — presiunea absolută [kPa]
// t      — temperatura fluidului [°C]
// tip    — tipul dispozitivului de strangulare
// d_int  — diametrul interior al conductei la temperatura de referință 20 °C [mm]
// d_orif — diametrul orificiului la temperatura de referință 20 °C [mm]
// ro     — densitatea gazului la condiții (t, p) [kg/m³]
// eta    — viscozitatea dinamică la condiții (t, p) [Pa·s]
// out    — ieșire: viteza medie [m/s] și pierderea de presiune [kPa]
float CalcMassFlow(float dp, float p, float t,
            TipDispozitiv tip, float d_int, float d_orif,
            float ro, float eta, FlowResult* out) {
  float d_i = d_int  * (1 + 0.0000122f * (t - kRefTempCelsius));
  float d_o = d_orif * (1 + 0.0000165f * (t - kRefTempCelsius));

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
  if ((IsDiaphragm(tip) && d_o < 12.5f)
      || (tip == TipDispozitiv::kAjutajVenturi && d_o <= 50)) {
    PrintError(ErrorCode::kOrificuStrangulare, 0);
    return 0;
  }

  float beta = d_o / d_i;

  if (((beta < 0.23f  || beta > 0.8f)  && tip == TipDispozitiv::kDiafragmaUnghi)
      || ((beta < 0.20f || beta > 0.75f) && (tip == TipDispozitiv::kDiafragmaFlansa
                                             || tip == TipDispozitiv::kDiafragmaDD2))
      || ((beta < 0.3f  || beta > 0.8f) && tip == TipDispozitiv::kAjutajIsa)
      || ((beta < 0.2f  || beta > 0.8f) && tip == TipDispozitiv::kAjutajRazaLunga)) {
    PrintError(ErrorCode::kRaportStrangulare, 0);
    return 0;
  }
  if (((beta < 0.3f   || beta > 0.75f)  && tip == TipDispozitiv::kVenturiBrut)
      || ((beta < 0.4f   || beta > 0.75f) && tip == TipDispozitiv::kVenturiPrelucrat)
      || ((beta < 0.4f   || beta > 0.7f)  && tip == TipDispozitiv::kVenturiTabla)
      || ((beta < 0.316f || beta > 0.775f) && tip == TipDispozitiv::kAjutajVenturi)) {
    PrintError(ErrorCode::kRaportStrangulare, 0);
    return 0;
  }

  float eps, y;
  if (IsDiaphragm(tip)) {
    eps = 1 - (0.41f + 0.35f * std::pow(beta, 4)) * dp / p / 1.31f;
  } else {
    y   = 1 - dp / p;
    eps = std::sqrt(1.31f * std::pow(y, 1.52671f) / 0.31f
        * (1 - std::pow(beta, 4))
        / (1 - std::pow(beta, 4) * std::pow(y, 1.52671f))
        * (1 - std::pow(y, 0.236641f)) / (1 - y));
  }

  d_i /= 1000;
  d_o /= 1000;

  float red = kInitialReynolds, qn = 0, q0 = 0, alfa = 0;
  do {
    q0   = qn;
    alfa = VelocityCoefficient(tip, d_i, beta, red);
    qn   = alfa * eps * kPi / 4 * std::pow(d_o, 2) * std::sqrt(2000 * dp * ro);
    red  = 4 * qn / (d_i * kPi * eta);  // Re = 4*Qm / (pi*D*mu)
  } while (std::fabs(qn - q0) > kReynoldsTolerance);

  bool reynolds_valid = false;
  switch (tip) {
    case TipDispozitiv::kDiafragmaUnghi:
      if ((5000  <= red) && (red <= 1e8) && (0.23f <= beta) && (beta < 0.45f))
        reynolds_valid = true;
      if ((10000 <= red) && (red <= 1e8) && (0.45f <= beta) && (beta < 0.77f))
        reynolds_valid = true;
      if ((20000 <= red) && (red <= 1e8) && (0.77f <= beta) && (beta <= 0.80f))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kDiafragmaFlansa:
    case TipDispozitiv::kDiafragmaDD2:
      if ((1.26e6f * beta * beta * d_i <= red) && (red <= 1e8))
        reynolds_valid = true;
      break;
    case TipDispozitiv::kAjutajIsa:
      if ((70000 <= red) && (red <= 1e7) && (0.30f <= beta) && (beta < 0.44f))
        reynolds_valid = true;
      if ((20000 <= red) && (red <= 1e7) && (0.44f <= beta) && (beta <= 0.80f))
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
float CalcDensity(float t, float p, const BwrConst& bwr) {
  float T   = t + kKelvinOffset;
  float R   = kGasConstantR;
  float ro1 = kRoMin;
  float ro2 = kRoMax;
  float ro  = 0.0f;
  float pcal = 0.0f;
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
