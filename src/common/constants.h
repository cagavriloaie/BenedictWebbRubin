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

constexpr double SECONDS_PER_HOUR   = 3600.0;   // s/h, for Qm [kg/s] → Q [m³/h]
constexpr double PA_TO_MICRO_PA     = 1.0e6;    // μPa·s per Pa·s
constexpr int    COMP_NAME_WIDTH    = 25;        // column width for component names

// UTF-8 string constants — use #define so they splice into string literals
// Greek / physics
#define U_RHO    "\xcf\x81"          // ρ
#define U_ETA    "\xce\xb7"          // η
#define U_BETA   "\xce\xb2"          // β
#define U_EPS    "\xce\xb5"          // ε
#define U_DELTA  "\xCE\x94"          // Δ
#define U_SIGMA  "\xce\xa3"          // Σ
#define U_MICRO  "\xC2\xB5"          // µ
// Math / punctuation
#define U_MDASH  "\xe2\x80\x94"      // —
#define U_NDASH  "\xe2\x80\x93"      // –
#define U_CDOT   "\xC2\xB7"          // ·
#define U_PLUSMN "\xc2\xb1"          // ±
#define U_APPROX "\xe2\x89\x88"      // ≈
#define U_NEQ    "\xe2\x89\xa0"      // ≠
#define U_LEQ    "\xe2\x89\xa4"      // ≤
#define U_IN     "\xe2\x88\x88"      // ∈
#define U_TIMES  "\xc3\x97"          // ×
#define U_DEG    "\xC2\xB0"          // °
#define U_MINUS  "\xe2\x88\x92"      // − (math minus)
#define U_COPY   "\xc2\xa9"          // ©
#define U_RAQUO  "\xc2\xbb"          // »
// Superscripts / subscripts
#define U_SUP2   "\xc2\xb2"          // ²
#define U_SUP3   "\xC2\xB3"          // ³
#define U_SUP4   "\xe2\x81\xb4"      // ⁴
#define U_SUB2   "\xe2\x82\x82"      // ₂
#define U_SUB4   "\xe2\x82\x84"      // ₄
#define U_SUB6   "\xe2\x82\x86"      // ₆
// Arrows
#define U_RARR   "\xE2\x86\x92"      // →
#define U_LARR   "\xe2\x86\x90"      // ←
// Box-drawing (table borders)
#define U_HLINE  "\xe2\x94\x80"      // ─
#define U_DLINE  "\xe2\x95\x90"      // ═
#define U_VLINE  "\xe2\x94\x82"      // │
#define U_TL     "\xe2\x94\x8c"      // ┌
#define U_TR     "\xe2\x94\x90"      // ┐
#define U_BL     "\xe2\x94\x94"      // └
#define U_BR     "\xe2\x94\x98"      // ┘
#define U_ML     "\xe2\x94\x9c"      // ├
#define U_MR     "\xe2\x94\xa4"      // ┤
#define U_MT     "\xe2\x94\xac"      // ┬
#define U_MB     "\xe2\x94\xb4"      // ┴
#define U_CROSS  "\xe2\x94\xbc"      // ┼

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
