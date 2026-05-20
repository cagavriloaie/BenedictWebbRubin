#pragma once
#include "constants.h"

enum class DeviceType {
  ORIFICE_CORNER     = 1,
  ORIFICE_FLANGE     = 2,
  ORIFICE_DD2        = 3,
  NOZZLE_ISA         = 4,
  NOZZLE_LONG_RADIUS = 5,
  VENTURI_ROUGH_CAST = 6,
  VENTURI_MACHINED   = 7,
  VENTURI_WELDED_SHEET = 8,
  VENTURI_NOZZLE     = 9,
};

enum class ErrorCode {
  PIPE_DIAMETER    = 1,
  ORIFICE_DIAMETER = 2,
  THROTTLE_RATIO   = 3,
  REYNOLDS         = 4,
  NUMERIC          = 5,
};

// BWRS mixture constants, computed once from composition.
// Starling (1973) equation of state:
//   p = ρRT + (B₀RT−A₀−C₀/T²+D₀/T³−E₀/T⁴)ρ²
//           + (bRT−a−d/T)ρ³ + α(a+d/T)ρ⁶ + (c/T²)ρ³(1+γρ²)exp(−γρ²)
struct BwrConst {
  double m_a0        = 0.0;
  double m_b0        = 0.0;
  double m_c0        = 0.0;
  double m_d0        = 0.0;
  double m_e0        = 0.0;
  double m_a         = 0.0;
  double m_b         = 0.0;
  double m_c         = 0.0;
  double m_dv        = 0.0;
  double m_alpha     = 0.0;
  double m_gamma     = 0.0;
  double m_molarMass = 0.0;  // M [g/mol]
};

// Hydraulic results returned by calcMassFlow()
struct FlowResult {
  double m_velocity      = 0.0;             // mean gas velocity in pipe [m/s]
  double m_pressureLoss  = 0.0;             // permanent pressure loss [kPa]
  double m_beta          = 0.0;             // throttling ratio d/D at measurement temperature [-]
  double m_reynolds      = 0.0;             // Reynolds number after iteration [-]
  double m_coefC         = 0.0;             // discharge coefficient C [-]
  double m_epsilon       = 0.0;             // expansibility factor ε [-]
  double m_dWorkingMm    = 0.0;             // orifice diameter at working temperature [mm]
  double m_DWorkingMm    = 0.0;             // pipe diameter at working temperature [mm]
  double m_kappa         = KAPPA_DEFAULT;   // isentropic exponent κ used for ε [-]
};

// Volumetric reference conditions (one preset = 1 or 2 T/101.325 kPa pairs)
struct CountryRef {
  const char* m_name;
  int         m_n;
  double      m_t[2];
  const char* m_label[2];
};

constexpr bool isDiaphragm(DeviceType tip) {
  return static_cast<int>(tip) <= static_cast<int>(DeviceType::ORIFICE_DD2);
}

constexpr const char* tipName(DeviceType tip) {
  switch (tip) {
    case DeviceType::ORIFICE_CORNER:      return "Orifice plate " U_MDASH " corner taps";
    case DeviceType::ORIFICE_FLANGE:      return "Orifice plate " U_MDASH " flange taps";
    case DeviceType::ORIFICE_DD2:         return "Orifice plate " U_MDASH " D and D/2 taps";
    case DeviceType::NOZZLE_ISA:          return "ISA 1932 nozzle";
    case DeviceType::NOZZLE_LONG_RADIUS:  return "Long-radius nozzle";
    case DeviceType::VENTURI_ROUGH_CAST:  return "Classical Venturi tube " U_MDASH " rough-cast convergent";
    case DeviceType::VENTURI_MACHINED:    return "Classical Venturi tube " U_MDASH " machined convergent";
    case DeviceType::VENTURI_WELDED_SHEET: return "Classical Venturi tube " U_MDASH " rough-welded sheet-metal convergent";
    case DeviceType::VENTURI_NOZZLE:      return "Venturi nozzle";
    default:                              return "";
  }
}

constexpr int TIP_MIN = static_cast<int>(DeviceType::ORIFICE_CORNER);
constexpr int TIP_MAX = static_cast<int>(DeviceType::VENTURI_NOZZLE);

inline const char TIP_DISP[] =
    "\n  Throttling device:\n\n"
    "  1.  Orifice plate " U_MDASH " corner taps\n"
    "  2.  Orifice plate " U_MDASH " flange taps\n"
    "  3.  Orifice plate " U_MDASH " D and D/2 taps\n"
    "  4.  ISA 1932 nozzle\n"
    "  5.  Long-radius nozzle\n"
    "  6.  Classical Venturi tube " U_MDASH " rough-cast convergent\n"
    "  7.  Classical Venturi tube " U_MDASH " machined convergent\n"
    "  8.  Classical Venturi tube " U_MDASH " rough-welded sheet-metal convergent\n"
    "  9.  Venturi nozzle\n\n"
    "  Select (1" U_NDASH "9) >";
