#include "bwrs.h"
#include "../common/constants.h"
#include "../common/types.h"
#include "../data/components.h"
#include <algorithm>
#include <cmath>

// Computes BWRS mixture constants from molar fractions x[1..NUM_COMPONENTS]
// using Starling (1973) quadratic (A₀…E₀, γ) and cubic (a, b, c, d, α) mixing rules.
BwrConst calcBwrMixture(const double* x) {
  BwrConst bt;

  // Quadratic (2nd-order) mixing rules for A₀, B₀, C₀, D₀, E₀, γ
  for (int i = 1; i <= NUM_COMPONENTS; i++) {
    for (int j = 1; j <= NUM_COMPONENTS; j++) {
      if (x[i] <= 0.0 || x[j] <= 0.0) continue;
      double kk = 1 - BWRS_MIX_COEF * std::sqrt(BWRS_V[i] * BWRS_V[j])
               / std::pow(std::cbrt(BWRS_V[i]) + std::cbrt(BWRS_V[j]), 3);
      bt.m_a0    += x[i] * x[j] * std::sqrt(BWRS_A0[i] * BWRS_A0[j]) * (1 - kk);
      bt.m_b0    += x[i] * x[j] * std::sqrt(BWRS_B0[i] * BWRS_B0[j]);
      bt.m_c0    += x[i] * x[j] * std::sqrt(BWRS_C0[i] * BWRS_C0[j]) * std::pow(1 - kk, 3);
      bt.m_d0    += x[i] * x[j] * std::sqrt(BWRS_D0[i] * BWRS_D0[j]) * std::pow(1 - kk, 4);
      bt.m_e0    += x[i] * x[j] * std::sqrt(BWRS_E0[i] * BWRS_E0[j]) * std::pow(1 - kk, 5);
      bt.m_gamma += x[i] * x[j] * std::sqrt(BWRS_GAMA[i] * BWRS_GAMA[j]);
    }
  }

  // Cubic (3rd-order) mixing rules for a, b, c, d, α
  for (int i = 1; i <= NUM_COMPONENTS; i++) {
    for (int j = 1; j <= NUM_COMPONENTS; j++) {
      for (int l = 1; l <= NUM_COMPONENTS; l++) {
        if (x[i] <= 0.0 || x[j] <= 0.0 || x[l] <= 0.0) continue;
        double k1 = 1 - BWRS_MIX_COEF * std::sqrt(BWRS_V[i] * BWRS_V[j])
                 / std::pow(std::cbrt(BWRS_V[i]) + std::cbrt(BWRS_V[j]), 3);
        double k2 = 1 - BWRS_MIX_COEF * std::sqrt(BWRS_V[i] * BWRS_V[l])
                 / std::pow(std::cbrt(BWRS_V[i]) + std::cbrt(BWRS_V[l]), 3);
        double k3 = 1 - BWRS_MIX_COEF * std::sqrt(BWRS_V[j] * BWRS_V[l])
                 / std::pow(std::cbrt(BWRS_V[j]) + std::cbrt(BWRS_V[l]), 3);
        bt.m_a     += x[i] * x[j] * x[l]
                    * std::cbrt(BWRS_A[i]*BWRS_A[j]*BWRS_A[l]) * (1-k1)*(1-k2)*(1-k3);
        bt.m_b     += x[i] * x[j] * x[l] * std::cbrt(BWRS_B[i]*BWRS_B[j]*BWRS_B[l]);
        bt.m_c     += x[i] * x[j] * x[l]
                    * std::cbrt(BWRS_C[i]*BWRS_C[j]*BWRS_C[l]) * (1-k1)*(1-k2)*(1-k3);
        bt.m_dv    += x[i] * x[j] * x[l]
                    * std::cbrt(BWRS_DV[i]*BWRS_DV[j]*BWRS_DV[l]) * (1-k1)*(1-k2)*(1-k3);
        bt.m_alpha += x[i] * x[j] * x[l]
                    * std::cbrt(BWRS_ALFA[i]*BWRS_ALFA[j]*BWRS_ALFA[l]);
      }
    }
  }

  for (int i = 1; i <= NUM_COMPONENTS; i++)
    bt.m_molarMass += x[i] * MOLAR_MASS_TABLE[i];

  return bt;
}

// Solves the BWRS equation for gas-phase molar density [kg/m³] via bisection on pressure.
// t [°C], p [atm]; returns 0 on invalid input or convergence failure.
double calcDensity(double t, double p, const BwrConst& bwr, int* iters) {
  double tKelvin = t + KELVIN_OFFSET;
  if (!std::isfinite(tKelvin) || !std::isfinite(p) || !std::isfinite(bwr.m_molarMass)
      || tKelvin <= 0.0 || p <= 0.0 || bwr.m_molarMass <= 0.0) {
    if (iters) *iters = 0;
    return 0.0;
  }

  auto pressureAt = [&](double ro) {
    return GAS_CONSTANT_R * tKelvin * ro
         + (bwr.m_b0 * GAS_CONSTANT_R * tKelvin - bwr.m_a0 - bwr.m_c0 / tKelvin / tKelvin
            + bwr.m_d0 / (tKelvin * tKelvin * tKelvin)
            - bwr.m_e0 / (tKelvin * tKelvin * tKelvin * tKelvin)) * ro * ro
         + (bwr.m_b  * GAS_CONSTANT_R * tKelvin - bwr.m_a - bwr.m_dv / tKelvin)
           * std::pow(ro, 3)
         + bwr.m_alpha * (bwr.m_a + bwr.m_dv / tKelvin) * std::pow(ro, 6)
         + (bwr.m_c / tKelvin / tKelvin) * std::pow(ro, 3)
         * (1 + bwr.m_gamma * ro * ro)
         * std::exp(-bwr.m_gamma * ro * ro);
  };

  // Adaptive bracket: start from ideal-gas estimate, double ro2 until p(ro2) >= p.
  // Stop if pressure starts decreasing (van der Waals spinodal) to stay gas-phase.
  double ro1 = RO_MIN;
  double pLow = pressureAt(ro1);
  double ro2  = std::max(p / (GAS_CONSTANT_R * tKelvin), ro1 * 2.0);
  double pHigh = pressureAt(ro2);
  double pPrev = pLow;
  for (int k = 0; k < 60; k++) {
    if (!std::isfinite(pHigh) || pHigh < pPrev) break;
    if (pHigh >= p) break;
    pPrev = pHigh;
    ro1   = ro2;  pLow = pHigh;
    ro2   = std::min(ro2 * 2.0, RO_MAX);
    pHigh = pressureAt(ro2);
  }

  if (!std::isfinite(pLow) || !std::isfinite(pHigh) || p < pLow || p > pHigh) {
    if (iters) *iters = 0;
    return 0.0;
  }

  double ro   = 0.0;
  double pCal = 0.0;
  int nIter = 0;
  for (; nIter < MAX_DENSITY_ITER; ) {
    nIter++;
    ro   = (ro1 + ro2) / 2;
    pCal = pressureAt(ro);
    if (!std::isfinite(pCal)) {
      if (iters) *iters = nIter;
      return 0.0;
    }
    if (std::fabs(p - pCal) < DENSITY_TOLERANCE) {
      if (iters) *iters = nIter;
      return bwr.m_molarMass * ro;
    }
    if (pCal > p) ro2 = ro;
    else          ro1 = ro;
  }
  if (iters) *iters = nIter;
  return 0.0;
}

// Computes pseudo-critical mixture properties (critical density ρc, viscosity parameter ξ)
// using Kay's linear mixing rule on Tc, Pc, Zc.
MixtureThermo calcMixtureThermo(const double* x, const BwrConst& bwr) {
  double tCritMix = 0.0, pCritMix = 0.0, zCritMix = 0.0;
  for (int i = 1; i <= NUM_COMPONENTS; i++) {
    tCritMix += x[i] * TEMP_CRITICAL_TABLE[i];
    pCritMix += x[i] * PRESS_CRITICAL_TABLE[i];
    zCritMix += x[i] * Z_CRITICAL_TABLE[i];
  }
  MixtureThermo thermo;
  thermo.m_rocCrit = pCritMix / (GAS_CONSTANT_R * zCritMix * tCritMix);
  thermo.m_csi     = std::pow(tCritMix, 6)
                   / std::sqrt(bwr.m_molarMass)
                   / std::cbrt(pCritMix * pCritMix);
  return thermo;
}

// Returns the isentropic exponent κ = Cp°/(Cp°−R) for the mixture from ideal-gas heat capacities.
double calcKappa(const double* x) {
  double cpMix = 0.0;
  for (int i = 1; i <= NUM_COMPONENTS; i++)
    cpMix += x[i] * CP0_IDEAL[i];
  return cpMix / (cpMix - GAS_CONSTANT_RSI);
}
