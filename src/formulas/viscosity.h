#pragma once

// Computes dynamic viscosity η [Pa·s] at temperature tCelsius [°C] and
// density ro [kg/m³] for the mixture given by molar fractions x[1..NUM_COMPONENTS].
//
// Uses the Chapman-Enskog dilute-gas model with a high-density correction
// (Lucas / Chung-Lee-Starling):
//   η = η₀(T, x) + kA/ξ × [exp(k₁·ρᵣ) − exp(−k₂·ρᵣ)]^k₃
//
// rocCrit [mol/L]  and  csi [K^6 / (g/mol)^0.5 / atm^(2/3)]  are
// obtained from calcMixtureThermo().
double calcViscosity(double tCelsius, double ro, const double* x,
                     double rocCrit, double csi);
