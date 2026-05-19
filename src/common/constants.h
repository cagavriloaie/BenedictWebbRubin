#pragma once

constexpr double PI = 3.141592653589;

constexpr int NUM_COMPONENTS = 35;
constexpr int ARRAY_SIZE     = NUM_COMPONENTS + 1;

constexpr double KELVIN_OFFSET       = 273.15;
constexpr double GAS_CONSTANT_R      = 0.082057366;  // CODATA 2018 [atm·L/(mol·K)]
constexpr double GAS_CONSTANT_RSI    = 8.31446;       // R [J/(mol·K)] for κ = Cp/(Cp−R)
constexpr double REF_TEMP_CELSIUS    = 20.0;
constexpr double KPA_PER_ATM         = 101.325;
constexpr double REYNOLDS_ISO_REF    = 1.0e6;
constexpr double INITIAL_REYNOLDS    = 1.0e6;
constexpr double REYNOLDS_TOLERANCE  = 1.0e-6;
constexpr double DENSITY_TOLERANCE   = 5.0e-4;
constexpr double RO_MIN              = 1.0e-4;
constexpr double RO_MAX              = 25.0;
constexpr int    MAX_REYNOLDS_ITER   = 100;
constexpr int    MAX_DENSITY_ITER    = 200;
constexpr double SUM_TOLERANCE       = 1.0e-3;

// Thermal expansion coefficients (linear, per °C, at 20 °C reference)
constexpr double THERMAL_EXP_PIPE    = 12.2e-6;  // carbon steel pipe (D), ISO 5167
constexpr double THERMAL_EXP_ORIFICE = 16.5e-6;  // stainless steel orifice plate (d)

// ISO 5167 flanged-tap offset [m] = 25.4 mm = 1 inch
constexpr double FLANGE_TAP_M = 0.0254;

// Factor in Qm formula: 2 × (Pa/kPa) so Δp[kPa]·ρ enters as SI
constexpr double FACTOR_2KPA = 2000.0;

// Isentropic exponent (default air/typical gas mixture; ISO 5167 ε formula)
constexpr double KAPPA_DEFAULT = 1.31;

// BWRS quadratic mixing-rule coefficient (Starling 1973)
constexpr double BWRS_MIX_COEF = 8.0;

// Chapman-Enskog polar-correction parameter (viscosity dilute-gas term)
constexpr double CHAP_ENSKOG = 0.323;

// High-density viscosity correction (Lucas / Chung-Lee-Starling)
constexpr double VISC_HIGH_A    = 10.8e-8;
constexpr double VISC_HIGH_EXP1 = 1.439;
constexpr double VISC_HIGH_EXP2 = 1.111;
constexpr double VISC_HIGH_POW  = 1.358;

// ANSI color codes
constexpr const char* COLOR_RESET       = "\033[0;37m";
constexpr const char* COLOR_BOLD_YELLOW = "\033[1;33m";
constexpr const char* COLOR_BOLD_WHITE  = "\033[0;37m";
constexpr const char* COLOR_BOLD_GREEN  = "\033[0;32m";
constexpr const char* COLOR_BOLD_RED    = "\033[0;31m";
constexpr const char* COLOR_YELLOW      = "\033[0;33m";
constexpr const char* COLOR_CYAN        = "\033[0;36m";
constexpr const char* COLOR_HDR_YELLOW  = "\033[0;33m";
constexpr const char* COLOR_HDR_GREEN   = "\033[0;32m";
constexpr const char* COLOR_HDR_CYAN    = "\033[0;36m";
