#include "viscosity.h"
#include "../Common/constants.h"
#include "../Data/components.h"
#include <cmath>

// Computes gas mixture dynamic viscosity [mPa·s] using Chapman-Enskog dilute-gas theory
// with a Chung-Lee-Starling high-density correction term.
// tCelsius [°C], ro [kg/m³], rocCrit [kg/m³], csi — viscosity parameter from calcMixtureThermo.
double calcViscosity(double tCelsius, double ro, const double* x,
                     double rocCrit, double csi) {
  double tKelvin  = tCelsius + KELVIN_OFFSET;
  double rocRed   = ro / rocCrit;

  // Σ xᵢ √Mᵢ — denominator of the Chapman-Enskog mixing rule
  double mxSum = 0.0;
  for (int i = 1; i <= NUM_COMPONENTS; i++)
    mxSum += x[i] * std::sqrt(MOLAR_MASS_TABLE[i]);

  // Dilute-gas (zero-density) mixture viscosity via Chapman-Enskog kinetic theory
  double viscosity = 0.0;
  for (int i = 1; i <= NUM_COMPONENTS; i++) {
    viscosity += (1 + CHAP_ENSKOG * std::log(tKelvin / VISC_CS_PARAM[i]))
               / (1 + CHAP_ENSKOG * std::log(KELVIN_OFFSET / VISC_CS_PARAM[i]))
               * std::sqrt(tKelvin / KELVIN_OFFSET)
               * VISC_ET_PARAM[i] * x[i] * std::sqrt(MOLAR_MASS_TABLE[i]);
  }
  viscosity /= mxSum;

  // High-density correction term
  viscosity += VISC_HIGH_A / csi
             * std::pow(
                 std::exp(VISC_HIGH_EXP1 * rocRed) - std::exp(-VISC_HIGH_EXP2 * rocRed),
                 VISC_HIGH_POW);

  return viscosity;
}
