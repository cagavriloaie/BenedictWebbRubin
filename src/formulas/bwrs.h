#pragma once
#include "../common/constants.h"
#include "../common/types.h"

// Pseudo-critical properties needed for the viscosity correlation.
struct MixtureThermo {
  double m_rocCrit;  // pseudo-critical molar density [mol/L]
  double m_csi;      // viscosity reducing parameter ξ
};

// Computes BWRS mixture constants from molar composition x[1..NUM_COMPONENTS].
// Uses Starling (1973) quadratic (2nd-order) and cubic (3rd-order) mixing rules.
BwrConst calcBwrMixture(const double* x);

// Solves the BWRS equation for mixture density ρ [kg/m³] via bisection on the
// gas-phase branch.  p is absolute pressure [atm].
// Returns 0 on convergence failure.
double calcDensity(double t, double p, const BwrConst& bwr, int* iters = nullptr);

// Computes mixture pseudo-critical properties needed for the viscosity model.
MixtureThermo calcMixtureThermo(const double* x, const BwrConst& bwr);

// Computes the mixture isentropic exponent κ = Cp°/(Cp° − R) at 20 °C
// from ideal-gas molar heat capacities.
double calcKappa(const double* x);
