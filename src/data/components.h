#pragma once
#include "../Common/constants.h"

// BWRS equation-of-state parameters (Starling 1973).
// Arrays indexed 1..NUM_COMPONENTS; index 0 is unused (= 0).
// BWRS_B0, BWRS_B, BWRS_ALFA, BWRS_GAMA, BWRS_V, VISC_ET_PARAM are stored pre-scaled (see components.cpp).

extern const double BWRS_A0   [ARRAY_SIZE];  // A₀
extern const double BWRS_B0   [ARRAY_SIZE];  // B₀  (original ÷ 100)
extern const double BWRS_C0   [ARRAY_SIZE];  // C₀  (original × 1e5)
extern const double BWRS_A    [ARRAY_SIZE];  // a
extern const double BWRS_B    [ARRAY_SIZE];  // b   (original ÷ 100)
extern const double BWRS_C    [ARRAY_SIZE];  // c   (original × 1e5)
extern const double BWRS_ALFA [ARRAY_SIZE];  // α   (original ÷ 1000)
extern const double BWRS_GAMA [ARRAY_SIZE];  // γ   (original ÷ 100)
extern const double BWRS_V    [ARRAY_SIZE];  // V   (original ÷ 1000)
extern const double BWRS_D0   [ARRAY_SIZE];  // D₀  [atm·L²·K³/mol²]
extern const double BWRS_E0   [ARRAY_SIZE];  // E₀  [atm·L²·K⁴/mol²]
extern const double BWRS_DV   [ARRAY_SIZE];  // d   [atm·L³·K/mol³]
extern const double MOLAR_MASS_TABLE [ARRAY_SIZE];  // molar mass M [g/mol]

// Transport-property parameters
extern const double VISC_CS_PARAM     [ARRAY_SIZE];  // Chapman-Enskog σ*ε [K]
extern const double VISC_ET_PARAM     [ARRAY_SIZE];  // low-pressure viscosity factor (original ÷ 1e4)
extern const double TEMP_CRITICAL_TABLE  [ARRAY_SIZE];  // critical temperature [K]
extern const double PRESS_CRITICAL_TABLE [ARRAY_SIZE];  // critical pressure [atm]
extern const double Z_CRITICAL_TABLE     [ARRAY_SIZE];  // critical compressibility factor

// Ideal-gas molar heat capacity at 20 °C [J/(mol·K)]
extern const double CP0_IDEAL[ARRAY_SIZE];

// Gross (superior) calorific value per component [MJ/mol] — ISO 6976:2016
// Inert components (He, Ar, N2, O2, CO2) = 0.
extern const double HHV_TABLE[ARRAY_SIZE];
