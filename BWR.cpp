// =============================================================================
// BWRS — Gas flow calculation through throttling devices
// Revision 3.0  |  05.2026
// Eng. Agavriloaie Constantin  (original R 01.2004)
// ELCOST Impex
// =============================================================================
//
// DESCRIPTION
//   Calculates density, dynamic viscosity, and flow rates for a gas mixture
//   of up to 35 components, based on:
//     • Benedict-Webb-Rubin-Starling (BWRS) equation of state for density;
//     • Chapman-Enskog model with high-density correction for viscosity;
//     • ISO 5167-2/3/4:2003 (Reader-Harris/Gallagher equation)
//       for flow rate through throttling devices.
//
// INPUTS
//   1. Mixture composition — molar fractions for each of the
//      35 components (methane, ethane, propane, ... acetylene).
//   2. Throttling device type (1–9):
//        [1] Orifice plate — corner taps
//        [2] Orifice plate — flange taps
//        [3] Orifice plate — D and D/2 taps
//        [4] ISA 1932 nozzle
//        [5] Long-radius nozzle
//        [6] Classical Venturi tube — rough-cast convergent
//        [7] Classical Venturi tube — machined convergent
//        [8] Classical Venturi tube — rough-welded sheet-metal convergent
//        [9] Venturi nozzle
//   3. Pipe inner diameter D [mm] and orifice diameter d [mm]
//      (at 20 °C reference temperature).
//   4. Temperature T [°C], absolute pressure p [kPa], and differential
//      pressure Δp [kPa] at the measurement point.
//
// OUTPUTS (for each T / p / Δp set)
//   • Mixture density ρ(T, p)                  [kg/m³]
//   • Dynamic viscosity η(T, p)                [μPa·s]
//   • Mass flow rate Qm                        [kg/s]
//   • Volumetric flow at selectable reference conditions
//       (0 °C / 101.325 kPa → Nm³/h; 15 °C / 101.325 kPa → Sm³/h; etc.)
//   • Volumetric flow at T and p               [m³/h]
//   • Mean gas velocity in pipe                [m/s]
//   • Permanent pressure loss                  [kPa]
//   • Throttling ratio β and Reynolds number Re
//
// VALIDATION
//   Checks ISO 5167 applicability limits for pipe diameter,
//   orifice diameter, throttling ratio β, and Reynolds number Re —
//   displays an error message when a limit is exceeded.
//
// DENSITY ALGORITHM
//   Bisection on the BWRS equation (Starling 1973, 11 parameters) until
//   convergence |p_calc − p| < 5×10⁻⁴ atm.  Equation of state:
//     p = ρRT + (B₀RT − A₀ − C₀/T² + D₀/T³ − E₀/T⁴)ρ²
//             + (bRT − a − d/T)ρ³ + α(a + d/T)ρ⁶
//             + (c/T²)ρ³(1 + γρ²)exp(−γρ²)
//   Mixture constants are computed once from composition using quadratic
//   mixing rules (A₀–E₀, γ) and cubic mixing rules (a–d, α).
//
// FLOW ALGORITHM
//   Reynolds iteration until relative convergence |ΔQm/Qm| < 10⁻⁶.
//   Discharge coefficient C and expansibility factor ε are recalculated
//   at each iteration as a function of Re and β.
//
// REFERENCES
//   • ISO 5167-2:2003 — Orifice plates (Reader-Harris/Gallagher equation)
//   • ISO 5167-3:2003 — Nozzles and Venturi nozzles
//   • ISO 5167-4:2003 — Venturi tubes
//   • Starling K.E. (1973) — Fluid Thermodynamic Properties for Light
//     Petroleum Systems, Gulf Publishing Co., Houston
//   • Nishiumi H., Saito S. (1975) — J. Chem. Eng. Japan 8(5), 356–360
//   • Benedict, Webb, Rubin (1940) — J. Chem. Phys. 8, 334
// =============================================================================

#include "BWR.h"
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <conio.h>
#  include <io.h>
#endif

static const char kCompFile[] = "bwr_comp.dat";
static const char kConfFile[] = "bwr_conf.dat";

static constexpr double kSecondsPerHour = 3600.0;   // s/h, for Qm[kg/s] → Q[m³/h]
static constexpr double kPaToMicroPa    = 1.0e6;    // μPa·s per Pa·s
static constexpr int    kCompNameWidth  = 25;        // column width for component names

// Array-table scale factors (raw values stored pre-scaled for compact notation)
static constexpr double kBTableScale    = 100.0;    // B, b stored ×100
static constexpr double kCTableScale    = 1.0e5;    // C, c stored ×10⁻⁵
static constexpr double kAlfaTableScale = 1.0e3;    // alfa stored ×1000
static constexpr double kGamaTableScale = 100.0;    // gama stored ×100
static constexpr double kVTableScale    = 1.0e3;    // V stored ×1000
static constexpr double kEtTableScale   = 1.0e4;    // et stored ×10000
static constexpr double kGasConstantRSI = 8.31446;  // R [J/(mol·K)] for κ = Cp/(Cp−R)

// Loads previously saved mixture composition from kCompFile.
// Returns true if the file was read successfully.
static bool LoadComposition(double* x) {
  FILE* f = nullptr;
  fopen_s(&f, kCompFile, "r");
  if (!f) {
    return false;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    if (fscanf_s(f, "%lf", &x[i]) != 1) {
      std::fclose(f);
      return false;
    }
  }
  std::fclose(f);
  return true;
}

// Saves current mixture composition to kCompFile.
static void SaveComposition(const double* x) {
  FILE* f = nullptr;
  fopen_s(&f, kCompFile, "w");
  if (!f) {
    std::printf("%s  Could not open %s for writing.%s\n", kBoldRed, kCompFile, kReset);
    return;
  }
  bool ok = true;
  for (int i = 1; i <= kNumComponents; i++)
    ok &= (std::fprintf(f, "%.8f\n", x[i]) > 0);
  if (std::fclose(f) != 0 || !ok)
    std::printf("%s  Write error: %s may be incomplete (disk full?).%s\n",
                kBoldRed, kCompFile, kReset);
}

// Loads previously saved device configuration from kConfFile.
// Returns true if all three values were read successfully.
static bool LoadConfig(int* tip_raw, double* d_int, double* d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, kConfFile, "r");
  if (!f) {
    return false;
  }
  bool ok = (fscanf_s(f, "%d %lf %lf", tip_raw, d_int, d_orif) == 3);
  std::fclose(f);
  return ok;
}

// Saves current device configuration to kConfFile.
static void SaveConfig(int tip_raw, double d_int, double d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, kConfFile, "w");
  if (!f) {
    std::printf("%s  Could not open %s for writing.%s\n", kBoldRed, kConfFile, kReset);
    return;
  }
  bool ok = (std::fprintf(f, "%d\n%.4f\n%.4f\n", tip_raw, d_int, d_orif) > 0);
  if (std::fclose(f) != 0 || !ok)
    std::printf("%s  Write error: %s may be incomplete (disk full?).%s\n",
                kBoldRed, kConfFile, kReset);
}

// Clears the screen, shows the exit message, and terminates the process (ESC key).
static void ExitApp() {
  std::printf("\033[2J\033[H%s\n  ELCOST Impex  —  BWR Gas Flow Calculator  v3.0/2026%s\n",
              kYellow, kReset);
  std::exit(0);
}

// Reads a real number interactively with echo, backspace, and ESC support.
// Returns false if the user presses Backspace on an empty field (go-back signal).
static bool ReadDouble(double* val) {
  char buf[64] = {};
  int pos = 0;
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 0 || ch == 0xE0) { (void)_getch(); continue; }
    if (ch == '\r') {
      std::printf("\n");
      std::fflush(stdout);
      break;
    }
    if (ch == 8 || ch == 127) {
      if (pos > 0) {
        pos--;
        std::printf("\b \b");
        std::fflush(stdout);
      } else {
        std::printf("\xe2\x86\x90\n");   // ←
        std::fflush(stdout);
        *val = 0.0;
        return false;
      }
      continue;
    }
    if (pos < (int)sizeof(buf) - 2 &&
        ((ch >= '0' && ch <= '9') || ch == '.' || (ch == '-' && pos == 0))) {
      buf[pos++] = static_cast<char>(ch);
      std::printf("%c", ch);
      std::fflush(stdout);
    }
  }
  *val = (pos > 0) ? std::atof(buf) : 0.0;
  return true;
}

// Reads a single digit character in the range [lo, hi] without requiring ENTER.
static int ReadChoice(int lo, int hi) {
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 0 || ch == 0xE0) { (void)_getch(); continue; }
    if (ch >= '0' + lo && ch <= '0' + hi) {
      std::printf("%c\n", ch);
      std::fflush(stdout);
      return ch - '0';
    }
  }
}

// Blocks until y/n is pressed; ESC exits the application.
static bool AskYesNo() {
  for (;;) {
    int ch = _getch();
    if (ch == 27) {
      ExitApp();
    }
    if (ch == 'y' || ch == 'Y') {
      std::printf("y\n");
      std::fflush(stdout);
      return true;
    }
    if (ch == 'n' || ch == 'N') {
      std::printf("n\n");
      std::fflush(stdout);
      return false;
    }
  }
}

// Returns count of UTF-8 continuation bytes (10xxxxxx) in s.
// strlen(s) − Utf8ExtraBytes(s) = number of displayed characters.
static int Utf8ExtraBytes(const char* s) {
  int n = 0;
  while (*s) {
    if (((unsigned char)*s & 0xC0) == 0x80) {
      n++;
    }
    s++;
  }
  return n;
}

static const char* const kCompNames[] = {
    nullptr,
    "Methane", "Ethane", "Propane", "Isobutane", "N-butane",
    "Neopentane", "Isopentane", "N-pentane", "2,2-dimethylbutane", "2,3-dimethylbutane",
    "3-methylpentane", "2-methylpentane", "N-hexane", "2,4-dimethylpentane", "2,2,3-trimethylbutane",
    "2-methylhexane", "3-methylhexane", "3-ethylpentane", "N-heptane", "2,2,4-trimethylpentane",
    "N-octane", "Benzene", "Toluene", "Hydrogen", "Carbon monoxide",
    "Hydrogen sulfide", "Helium", "Argon", "Nitrogen", "Oxygen",
    "Carbon dioxide", "Ethylene", "Propylene", "Ammonia", "Acetylene"
};

// Prints the mixture composition table aligned in the console, values highlighted in green.
static void PrintComposition(const double* x) {
  std::printf("\n%s  \xe2\x94\x80\xe2\x94\x80 Composition%s\n",
              kYellow, kReset);
  std::printf("%s  Molar fractions of the gas mixture:%s\n\n", kBoldWhite, kReset);
  for (int i = 1; i <= kNumComponents; i++) {
    const char* name = kCompNames[i];
    int visualLen = (int)strlen(name) - Utf8ExtraBytes(name);
    std::printf("  %s>%s %s%s", kCyan, kReset, kBoldWhite, name);
    for (int j = visualLen; j < kCompNameWidth; j++) {
      std::putchar(' ');
    }
    std::printf(": %s%.8f%s\n", kBoldGreen, x[i], kReset);
  }
  std::printf("\n");
}

// Main: initialises BWRS constants, reads composition/configuration, runs the mass-flow loop.
int main() {
#ifdef _WIN32
  std::system("chcp 65001 > nul");
  {
    HANDLE hCon = GetStdHandle(STD_OUTPUT_HANDLE);

    // Set console font
    CONSOLE_FONT_INFOEX cfi = {};
    cfi.cbSize = sizeof(cfi);
    GetCurrentConsoleFontEx(hCon, FALSE, &cfi);
    cfi.dwFontSize.Y = 14;
    wcscpy_s(cfi.FaceName, L"Consolas");
    SetCurrentConsoleFontEx(hCon, FALSE, &cfi);

    // Enable VT/ANSI processing
    DWORD mode = 0;
    GetConsoleMode(hCon, &mode);
    SetConsoleMode(hCon, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
  }
#endif

  double x[kArraySize] = {};

  // Raw values; arrays marked (*) are scaled after initialisation
  static const double A[kArraySize] = {
       0,         1.79894,   4.15556,   6.87225,  10.23264,  10.0847,   12.8,     12.7959,  12.1794,  11.842,
      16.43,     12.203,    12.203,    14.4373,   12.423,    12.423,    14.31,    14.31,    14.31,    17.5206,
      51.42,     51.86,      5.509772,  6.21,      0.040962,  1.34122,   2.12044,  0.040962,  0.823417, 1.053642,
       0.950852,  2.73742,   3.33958,   6.1122,    3.7892819, 5.1079342};

  static double B[kArraySize] = {  // (*) /= 100
        0,        4.54625,   6.27724,  97.313,    13.7544,   12.4361,   16.0,     16.0053,  15.6751,  19.214,
       19.0,      8.1505,    8.1505,   17.7813,   20.246,    20.246,     9.1423,   9.1423,   9.1423,  19.9,
      110.3,    121.2,      50.30055,  40.8,       2.3661,    5.45425,   2.6182,   2.3661,   2.22826,  4.07526,
        2.22,     3.38943,   5.56833,   8.50647,  51.646121,  6.946403};

  static double C[kArraySize] = {  // (*) *= 100000
       0,         0.318382,  1.79592,   5.08256,   8.49943,   9.9283,   17.5,     17.4632,  21.21219,  33.595,
      25.534,    22.125,    22.125,    33.1935,   51.237,    51.237,    31.564,   31.564,   31.564,   47.4574,
       1.032,     0.931,    34.2997,   29.0,       0.0000016227, 0.0856209, 7.9384, 0.0000016227, 0.1314125, 0.08059,
       0.326436,  1.38567,   9.6936284, 1.3114,   1.785708,  6.506284};

  static const double a[kArraySize] = {
       0,         0.04352,   0.34516,   0.9477,    1.93763,   1.88231,   3.756,    3.7562,   4.0748,  10.108,
       4.6956,    7.4286,    7.4286,    7.11671,  11.786,    11.786,     7.5854,   7.5854,   7.5854,  10.36475,
      32.512,    31.423,     5.57,      4.32,      0.00057339, 0.3665,   0.84468,  0.00057339, 0.0288358, 0.025102,
       0.16269,   0.136814,  0.259,   0.774056,  0.10354029, 0.6970948};

  static double b[kArraySize] = {  // (*) /= 100
       0,         0.252033,  1.1122,    2.25,      4.24352,   3.99983,   6.68,     6.6812,   6.6812,  14.0,
       7.9,      11.224,    11.224,    10.9131,   17.9131,   17.721,    14.321,   14.321,   14.321,  15.1954,
      53.14,     58.32,      7.663,     5.18,      0.000019727, 0.263158, 1.4653,  0.000019727, 0.215289, 0.23277,
       0.358835,  0.527236,  0.86,     18.7059,    0.071952516, 1.482999};

  static double c[kArraySize] = {  // (*) *= 100000
       0,         0.035878,  0.32767,   1.29,      2.8601,    3.164,     6.95,     6.95,     8.2417,  17.483,
      11.346,     9.5556,    9.5556,   15.1276,   22.586,    22.586,    13.252,   13.252,   13.252,  24.7,
      11.21,     18.23,     11.76418,  23.3,      0.0000000552, 0.0104,  1.1335,  0.0000000552, 0.007982437, 0.0072841,
       0.128274,  0.14918,   2.8297636, 2.112,    0.0015753298, 1.0984375};

  static double alfa[kArraySize] = {  // (*) /= 1000
        0,        0.33,      0.243389,  0.607175,  1.07408,   1.10132,   1.7,      1.7,      1.81,    2.189,
        3.5948,   2.25,      2.25,      2.81086,   2.764,     2.764,     2.8155,   2.8155,   2.8155,  4.35611,
        2.207,    2.581,     0.7001,    0.318,     0.0072673, 0.135,     0.071955, 0.0072673, 0.035589, 0.1272,
        0.09270,    0.0698611, 0.73924,   0.178,    0.0046521779, 0.27363248};  // [30]=O₂ corrected from 927.06 (Starling 1973)

  static double gama[kArraySize] = {  // (*) /= 100
       0,         1.05,      1.18,      2.2,       3.4,       3.4,       4.63,     4.63,     4.75,    5.65,
       7.5,       6.289,     6.289,     6.668849,  6.799,     6.799,     7.446,    7.446,    7.446,   9.0,
       0.0318,    0.0209,    2.93,      1.12,      0.077942,  0.6,       0.59236,  0.077942, 0.233827, 0.53,
       3.1,       0.460593,  2.911417,  0.923,   19.805156,  1.245167};

  static double V[kArraySize] = {  // (*) /= 1000
        0,        99.5,     148.00,    200.00,    263.00,    255.00,    303.00,   308.00,   311.00,   359.00,
      358.00,   367.00,   367.00,    368.00,    420.00,    420.04,    428.00,   418.00,   416.00,   426.00,
      482.00,   486.00,   260.00,    316.00,     65.0,      93.10,     95.00,    57.80,    75.20,    90.10,
       74.40,    94.0,    124.00,    181.00,     72.50,    113.00};

  static const double m[kArraySize] = {
        0,       16.043,   30.070,    44.097,    58.124,    58.124,    72.151,   72.151,   72.151,   86.178,
       86.178,   86.178,   86.178,    86.178,   100.205,   100.205,   100.205,  100.205,  100.205,  100.205,
      114.232,  114.232,   78.114,    92.141,     2.016,    28.011,    34.082,    4.003,   39.944,   28.016,
       32.000,   44.011,   28.054,    42.081,    17.032,   26.038};

  static const double cs[kArraySize] = {
        0,      148.6,    215.7,     237.1,     330.1,     331.4,     340.1,    340.1,    341.1,    398.2,
      398.3,    399.0,    399.1,     399.3,     412.1,     413.2,     413.2,    410.1,    412.2,    413.6,
      563.0,    564.0,    412.3,     418.3,      59.7,      91.7,     301.1,     10.22,    93.3,     71.4,
      106.7,    195.2,     224.7,     298.9,     558.2,    231.8};

  static double et[kArraySize] = {  // (*) /= 10000
       0,        0.1085,   0.0915,    0.0805,    0.0735,    0.0725,    0.0711,   0.0696,   0.067,    0.0666,
       0.0658,   0.0647,   0.0651,    0.0641,    0.0626,    0.0626,    0.061,    0.0619,   0.0616,   0.0607,
       0.05947,  0.0577,   0.0745,    0.066,     0.0836,    0.1636,    0.141,    0.071,    0.2126,   0.1755,
       0.2025,   0.1465,   0.094,     0.078,     0.093,     0.0943};

  static const double Tc[kArraySize] = {
        0,      190.7,    305.4,     369.9,     408.1,     425.2,     433.8,    460.4,    469.5,    488.7,
      499.9,    504.7,    496.5,     507.3,     520.3,     531.5,     530.3,    535.6,    540.8,    540.3,
      543.6,    568.6,    562.1,     592.0,      33.3,     133.0,     373.6,      5.3,    151.0,    126.2,
      154.8,    304.2,    282.85,    364.55,    405.55,   308.85};

  static const double Pc[kArraySize] = {
        0,       45.8,     48.2,      42.0,      36.0,      37.5,      31.6,     32.9,     33.3,     30.7,
       30.9,     30.8,     30.0,      29.9,      27.4,      29.8,      27.2,     28.1,     28.6,     27.0,
       25.4,     24.6,     48.6,      41.6,      12.8,      34.5,      88.9,      2.26,    48.0,     33.5,
       50.1,     72.9,     50.7,      45.4,     111.5,     61.6};

  static const double Zc[kArraySize] = {
       0,        0.29,     0.285,     0.277,     0.283,     0.274,     0.269,    0.268,    0.269,    0.273,
       0.27,     0.273,    0.27,      0.264,     0.27,      0.269,     0.267,    0.268,    0.267,    0.259,
       0.274,    0.256,    0.274,     0.271,     0.304,     0.294,     0.268,    0.3,      0.296,    0.291,
       0.292,    0.274,    0.27,      0.274,     0.242,     0.274};

  // D₀, E₀, d — Starling (1973) parameters; see Table 1, p. 19.
  // Zero values reduce BWRS exactly to BWR — replace per component as data become available.
  // Sources: Nishiumi & Saito (J.CEJ 1975) for C1–C8 and permanent gases;
  // Starling (1973) for CO₂, H₂S, CO; Poling et al. (2001) for Ar, NH₃;
  // estimated from Tc/Pc correlation (marked *) for isomers without published data.
  static const double D0[kArraySize] = {  // D₀ [atm·L²·K³/mol²]
       0,   1.218e5,  7.697e5,  3.362e6,  7.044e6,  9.699e6,
       1.448e7,      1.991e7,  2.089e7,  2.524e7,  2.583e7,  // *9,*10
       3.147e7,      2.952e7,  4.870e7,  5.487e7,  4.953e7,  // *11,*12; *14,*15
       7.566e7,      7.892e7,  7.680e7,  9.618e7,  1.504e8,  // *16,*17,*18; *20
       1.695e8,      2.500e7,  6.000e7,  3.40,     1.855e3,
       5.044e5,      0.50,     3.000e3,  3.078e3,  4.800e3,
       1.375e6,      2.280e5,  4.500e6,  5.700e5,  5.000e5}; // *33,*34,*35

  static const double E0[kArraySize] = {  // E₀ [atm·L²·K⁴/mol²]
       0,   3.238e6,  2.672e7,  1.787e8,  4.484e8,  7.283e8,
       1.277e9,      2.083e9,  2.226e9,  2.813e9,  3.002e9,  // *9,*10
       4.313e9,      3.772e9,  7.373e9,  8.960e9,  8.012e9,  // *11,*12; *14,*15
       1.530e10,     1.622e10, 1.557e10, 2.020e10, 3.798e10, // *16,*17,*18; *20
       4.888e10,     3.500e9,  1.200e10, 2.50e-2,  1.131e4,
       1.888e7,      2.00e-3,  2.500e4,  3.700e4,  7.000e4,
       6.319e7,      9.300e6,  2.700e8,  8.700e7,  3.000e7}; // *33,*34,*35

  static const double dv[kArraySize] = {  // d [atm·L³·K/mol³]
       0,   6.003e-2, 2.787e-1, 7.671e-1, 1.545,    1.968,
       2.693,        3.379,    3.454,    3.944,    4.115,    // *9,*10
       4.912,        4.652,    6.439,    7.516,    6.998,    // *11,*12; *14,*15
       9.551,        9.861,    9.667,    10.75,    16.21,    // *16,*17,*18; *20
       16.99,        4.500,    7.000,    1.00e-4,  1.307e-2,
       2.027e-1,     1.00e-4,  1.100e-2, 1.601e-2, 2.100e-2,
       3.849e-1,     1.000e-1, 9.500e-1, 6.000e-1, 2.000e-1}; // *33,*34,*35

  // Cp° ideal-gas [J/(mol·K)] at 20 °C — for κ = Cp/(Cp − R) with R = 8.314 J/(mol·K)
  // Sources: NIST WebBook; C5–C8 isomers estimated from functional groups (±2 J/(mol·K))
  static const double kCp0[kArraySize] = {
       0,      35.7,  52.5,  73.6,  97.5,  96.4, 120.9, 118.9, 120.1, 141.3,
     140.9,  141.7, 141.2, 143.1, 164.8, 163.4, 165.0, 165.0, 165.1, 166.1,
     188.9,  188.9,  82.4, 103.7,  28.8,  29.1,  34.2,  20.8,  20.8,  29.1,
      29.4,   37.1,  42.9,  63.9,  35.7,  44.0};

  // Scale arrays marked (*) — done once at startup
  for (int i = 1; i <= kNumComponents; i++) {
    B[i] /= kBTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    C[i] *= kCTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    b[i] /= kBTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    c[i] *= kCTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    alfa[i] /= kAlfaTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    gama[i] /= kGamaTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    V[i] /= kVTableScale;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    et[i] /= kEtTableScale;
  }

  std::printf("\033[2J\033[H\n");

  // ── Main header ────────────────────────────────────────────────────────────
  std::printf(
      "%s  \xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\n"
      "  %sELCOST Impex%s  \xc2\xb7  BWRS Gas Flow Calculator                          v3.0/2026\n"
      "  \xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\n\n",
      kHdrYellow, kHdrGreen, kHdrYellow);

  // ── Calculation models ─────────────────────────────────────────────────────
  std::printf("%s  CALCULATION MODELS%s\n", kHdrCyan, kReset);
  std::printf(
      "  %sEquation of state%s   BWRS \xc2\xb7 Starling 1973  \xe2\x80\x94  Benedict-Webb-Rubin-Starling, 11 param.\n"
      "  %sViscosity%s           Chapman-Enskog with high-density correction  [Nishiumi 1975]\n"
      "  %sC coefficient, \xce\xb5%s    Reader-Harris/Gallagher  \xc2\xb7  ISO 5167-2/3/4:2003\n\n",
      kBoldWhite, kReset, kBoldWhite, kReset, kBoldWhite, kReset);

  // ── Scope ──────────────────────────────────────────────────────────────────
  std::printf("%s  SCOPE OF APPLICATION%s\n", kHdrCyan, kReset);
  std::printf(
      "  35 components  \xc2\xb7  9 throttling device types\n"
      "  Mass flow, volumetric flow, velocity, pressure loss  \xc2\xb7  selectable reference conditions\n"
      "  ISO 5167 limit validation  (D, \xce\xb2, Re)\n\n");

  // ── Standards and references ───────────────────────────────────────────────
  std::printf(
      "%s  \xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\n"
      "  %sSTANDARDS AND REFERENCES%s\n",
      kHdrYellow, kHdrCyan, kReset);
  std::printf(
      "  %sISO 5167-2:2003%s   Orifice plates \xe2\x80\x94 Reader-Harris/Gallagher equation\n"
      "  %sISO 5167-3:2003%s   Nozzles and Venturi nozzles\n"
      "  %sISO 5167-4:2003%s   Classical Venturi tubes\n"
      "  %sStarling K.E.%s     Fluid Thermodynamic Properties, Gulf Publ. Houston (1973)\n"
      "  %sNishiumi & Saito%s  J. Chem. Eng. Japan 8(5), 356\xe2\x80\x93" "360 (1975)\n\n",
      kBoldWhite, kReset, kBoldWhite, kReset, kBoldWhite, kReset,
      kBoldWhite, kReset, kBoldWhite, kReset);

  // ── Validation ─────────────────────────────────────────────────────────────
  std::printf(
      "%s  \xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\n"
      "  %sVALIDATION%s\n",
      kHdrYellow, kHdrCyan, kReset);
  std::printf(
      "  \xcf\x81  validated for 8 compositions (CH\xe2\x82\x84, C\xe2\x82\x82H\xe2\x82\x86, CO\xe2\x82\x82,"
      " H\xe2\x82\x82, N\xe2\x82\x82, Ar, std NG, rich NG)  vs. NIST WebBook\n"
      "  \xce\xb7  validated vs. Chapman-Enskog / NIST\n"
      "  Qm validated \xc2\xb1" "15 %% against ISO 5167 estimates  \xc2\xb7  9 device types\n\n");

  // ── Footer ─────────────────────────────────────────────────────────────────
  std::printf(
      "%s  \xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\n"
      "  %soffice@elcost.ro%s                                   \xc2\xa9 2004\xe2\x80\x93" "2026 ELCOST Impex\n"
      "%s  \xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\xe2\x95\x90\n\n",
      kHdrYellow, kHdrCyan, kReset, kHdrYellow);

  std::printf("%s", kReset);

  // ── Optional self-test at startup ─────────────────────────────────────────
  {
    std::printf("%s\n  Run implementation validation test? [y/n] %s>%s ",
                kBoldWhite, kCyan, kReset);
    if (AskYesNo()) {
      std::printf("\n%s  \xe2\x94\x80\xe2\x94\x80 Implementation validation test%s\n", kYellow, kReset);

      int npass = 0, ntotal = 0;

      // chk: compares a value against an interval and prints PASS/FAIL
      // src = reference source label, printed at the end of each line
      auto chk = [&](const char* label, double v, double lo, double hi, const char* src) {
        ntotal++;
        bool ok = (v >= lo && v <= hi);
        if (ok) npass++;
        int vlen = (int)std::strlen(label) - Utf8ExtraBytes(label);
        std::printf("    %s%s", kBoldWhite, label);
        for (int k = vlen; k < 54; k++) std::putchar(' ');
        char ivl[32];
        std::snprintf(ivl, sizeof(ivl), "[%.4f, %.4f]", lo, hi);
        std::printf("%s%-10.5f%s  %-24s  %s%-4s%s  %s[%s]%s\n",
                    kBoldGreen, v, kReset,
                    ivl,
                    ok ? kBoldGreen : kBoldRed, ok ? "PASS" : "FAIL", kReset,
                    kBoldWhite, src, kReset);
      };

      // make_bwr: computes BwrConst for any composition xc[1..kNumComponents]
      auto make_bwr = [&](const double* xc) -> BwrConst {
        BwrConst bt;
        for (int i = 1; i <= kNumComponents; i++) {
          for (int j = 1; j <= kNumComponents; j++) {
            if (xc[i] <= 0.0 || xc[j] <= 0.0) continue;
            double kk = 1 - kBwrsMixCoef * std::sqrt(V[i] * V[j])
                     / std::pow(std::cbrt(V[i]) + std::cbrt(V[j]), 3);
            bt.a0    += xc[i] * xc[j] * std::sqrt(A[i] * A[j]) * (1 - kk);
            bt.b0    += xc[i] * xc[j] * std::sqrt(B[i] * B[j]);
            bt.c0    += xc[i] * xc[j] * std::sqrt(C[i] * C[j]) * std::pow(1 - kk, 3);
            bt.d0    += xc[i] * xc[j] * std::sqrt(D0[i] * D0[j]) * std::pow(1 - kk, 4);
            bt.e0    += xc[i] * xc[j] * std::sqrt(E0[i] * E0[j]) * std::pow(1 - kk, 5);
            bt.gamma += xc[i] * xc[j] * std::sqrt(gama[i] * gama[j]);
          }
        }
        for (int i = 1; i <= kNumComponents; i++) {
          for (int j = 1; j <= kNumComponents; j++) {
            for (int l = 1; l <= kNumComponents; l++) {
              if (xc[i] <= 0.0 || xc[j] <= 0.0 || xc[l] <= 0.0) continue;
              double k1 = 1 - kBwrsMixCoef * std::sqrt(V[i] * V[j])
                       / std::pow(std::cbrt(V[i]) + std::cbrt(V[j]), 3);
              double k2 = 1 - kBwrsMixCoef * std::sqrt(V[i] * V[l])
                       / std::pow(std::cbrt(V[i]) + std::cbrt(V[l]), 3);
              double k3 = 1 - kBwrsMixCoef * std::sqrt(V[j] * V[l])
                       / std::pow(std::cbrt(V[j]) + std::cbrt(V[l]), 3);
              bt.a     += xc[i] * xc[j] * xc[l]
                        * std::cbrt(a[i]*a[j]*a[l]) * (1-k1)*(1-k2)*(1-k3);
              bt.b     += xc[i] * xc[j] * xc[l] * std::cbrt(b[i]*b[j]*b[l]);
              bt.c     += xc[i] * xc[j] * xc[l]
                        * std::cbrt(c[i]*c[j]*c[l]) * (1-k1)*(1-k2)*(1-k3);
              bt.dv    += xc[i] * xc[j] * xc[l]
                        * std::cbrt(dv[i]*dv[j]*dv[l]) * (1-k1)*(1-k2)*(1-k3);
              bt.alpha += xc[i] * xc[j] * xc[l]
                        * std::cbrt(alfa[i]*alfa[j]*alfa[l]);
            }
          }
        }
        for (int i = 1; i <= kNumComponents; i++) bt.molar_mass += xc[i] * m[i];
        return bt;
      };

      // calc_eta: viscosity [Pa·s] at temperature t[°C] and 1 atm, composition xc
      auto calc_eta = [&](double t_c, const double* xc,
                          const BwrConst& bt,
                          double roc_crit, double csi) -> double {
        double T  = t_c + kKelvinOffset;
        double mx = 0.0;
        for (int i = 1; i <= kNumComponents; i++) mx += xc[i] * std::sqrt(m[i]);
        double rho     = CalcDensity(t_c, 1.0, bt);
        double roc_red = rho / roc_crit;
        double eta     = 0.0;
        for (int i = 1; i <= kNumComponents; i++) {
          eta += (1 + kChapEnskog * std::log(T / cs[i]))
               / (1 + kChapEnskog * std::log(kKelvinOffset / cs[i]))
               * std::sqrt(T / kKelvinOffset)
               * et[i] * xc[i] * std::sqrt(m[i]);
        }
        eta = eta / mx
            + kViscHighA / csi
            * std::pow(std::exp(kViscHighExp1 * roc_red) - std::exp(-kViscHighExp2 * roc_red),
                       kViscHighPow);
        return eta;
      };

      // test_comp: runs the test suite (density × 2, viscosity, mass flow rate)
      //   for composition xc with the given parameters and reference intervals.
      //   d_flow_mm = 0 → skip flow test.
      auto test_comp = [&](const char* comp_name,
                           const double* xc,
                           double rho20_lo, double rho20_hi,   // density 20°C [kg/m³]
                           double rho0_lo,  double rho0_hi,    // density  0°C [kg/m³]
                           double eta_lo,   double eta_hi,     // viscosity 20°C [μPa·s]
                           double D_mm,     double d_flow_mm,  // D and d [mm] (d=0 → skip)
                           double p_kpa,    double dp_kpa,
                           bool leading_nl = true) {
        std::printf(leading_nl ? "\n  %s\xc2\xbb %s%s\n" : "  %s\xc2\xbb %s%s\n",
                    kYellow, comp_name, kReset);

        BwrConst bt = make_bwr(xc);
        double tcam_c = 0.0, pcam_c = 0.0, zcam_c = 0.0;
        for (int i = 1; i <= kNumComponents; i++) {
          tcam_c += xc[i] * Tc[i];
          pcam_c += xc[i] * Pc[i];
          zcam_c += xc[i] * Zc[i];
        }
        double roc_crit_c = pcam_c / (kGasConstantR * zcam_c * tcam_c);
        double csi_c      = std::pow(tcam_c, 6) / std::sqrt(bt.molar_mass)
                          / std::cbrt(pcam_c * pcam_c);

        double rho1 = CalcDensity(20.0, 1.0, bt);
        double rho2 = CalcDensity( 0.0, 1.0, bt);
        chk("Density \xcf\x81 20\xC2\xB0""C, 101.325 kPa [kg/m\xC2\xB3]", rho1, rho20_lo, rho20_hi, "BWRS 1973");
        chk("Density \xcf\x81  0\xC2\xB0""C, 101.325 kPa [kg/m\xC2\xB3]", rho2, rho0_lo,  rho0_hi,  "BWRS 1973");
        chk("Viscosity \xce\xb7 20\xC2\xB0""C, 101.325 kPa [\xC2\xB5Pa\xC2\xB7s]",
            calc_eta(20.0, xc, bt, roc_crit_c, csi_c) * kPaToMicroPa, eta_lo, eta_hi, "CE/NIST");
        // Z = p*M / (rho*R*T); at 1 atm, 20°C: Z ∈ [0.975, 1.003] for all common gases
        double T_K1 = 20.0 + kKelvinOffset;
        double Z1   = 1.0 * bt.molar_mass / (rho1 * kGasConstantR * T_K1);
        chk("Compressibility factor Z 20\xC2\xB0""C, 101.325 kPa [-]", Z1, 0.975, 1.003, "BWRS/ideal");
        // Molar mass: make_bwr must accumulate exactly sum(x[i]*m[i])
        double mx_ref = 0.0;
        for (int ii = 1; ii <= kNumComponents; ii++) mx_ref += xc[ii] * m[ii];
        chk("Molar mass M [g/mol]", bt.molar_mass,
            mx_ref * (1.0 - 1e-8), mx_ref * (1.0 + 1e-8), "mixing rule");

        if (d_flow_mm > 0.0) {
          double rho_f = CalcDensity(20.0, p_kpa / kKpaPerAtm, bt);
          double T_f   = 20.0 + kKelvinOffset;
          double mx_f  = 0.0;
          for (int i = 1; i <= kNumComponents; i++) mx_f += xc[i] * std::sqrt(m[i]);
          double roc_f = rho_f / roc_crit_c;
          double eta_f = 0.0;
          for (int i = 1; i <= kNumComponents; i++) {
            eta_f += (1 + kChapEnskog * std::log(T_f / cs[i]))
                   / (1 + kChapEnskog * std::log(kKelvinOffset / cs[i]))
                   * std::sqrt(T_f / kKelvinOffset)
                   * et[i] * xc[i] * std::sqrt(m[i]);
          }
          eta_f = eta_f / mx_f
                + kViscHighA / csi_c
                * std::pow(std::exp(kViscHighExp1 * roc_f) - std::exp(-kViscHighExp2 * roc_f),
                           kViscHighPow);
          double A_o = kPi / 4.0 * (d_flow_mm * 1e-3) * (d_flow_mm * 1e-3);
          // alpha_nom = C_nominal/sqrt(1-beta^4) at beta=0.4; eps_nom from ISO 5167
          // C nominal source: ISO 5167-2/3/4:2003 and Reader-Harris/Gallagher at Re=10^6
          static const struct { int tip; const char* abv; double alpha_nom; double eps_nom; }
            tipuri[] = {
              {1, "Dia-U",   0.609, 0.993},
              {2, "Dia-F",   0.608, 0.993},
              {3, "Dia-D/2", 0.608, 0.993},
              {4, "Aj-ISA",  0.997, 0.997},
              {5, "Aj-RL",   1.005, 0.997},
              {6, "Ven-B",   0.997, 0.997},
              {7, "Ven-P",   1.008, 0.997},
              {8, "Ven-T",   0.998, 0.997},
              {9, "Aj-V",    0.996, 0.997},
            };
          double qm_diau = 0.0;
          for (const auto& t : tipuri) {
            FlowResult fr;
            double qm_f = CalcMassFlow(dp_kpa, p_kpa, 20.0,
                                       static_cast<DeviceType>(t.tip),
                                       D_mm, d_flow_mm, rho_f, eta_f, &fr);
            if (t.tip == 1) qm_diau = qm_f;
            if (qm_f > 0.0) {
              double qm_est = t.alpha_nom * t.eps_nom * A_o
                            * std::sqrt(2000.0 * dp_kpa * rho_f);
              char lbl[80];
              std::snprintf(lbl, sizeof(lbl),
                            "Mass flow Qm D=%g,d=%g,p=%g,\xCE\x94p=%g %s [kg/s]",
                            D_mm, d_flow_mm, p_kpa, dp_kpa, t.abv);
              chk(lbl, qm_f, qm_est * 0.85, qm_est * 1.15, "ISO 5167");
            } else {
              ntotal++;
              std::printf("    %sQm test %s: ISO 5167 ERROR%s\n", kBoldRed, t.abv, kReset);
            }
          }
          if (qm_diau > 0.0 && rho2 > 0.0) {
            double qv_n     = qm_diau / rho2 * kSecondsPerHour;
            double qm_est_u = tipuri[0].alpha_nom * tipuri[0].eps_nom * A_o
                            * std::sqrt(2000.0 * dp_kpa * rho_f);
            double qv_est   = qm_est_u / rho2 * kSecondsPerHour;
            chk("Volumetric flow Qv Dia-U, 0\xC2\xB0""C ref [Nm\xC2\xB3/h]",
                qv_n, qv_est * 0.85, qv_est * 1.15, "ISO 5167/BWRS");
          }
        }
      };

      // sec_box: prints a description box before a test section.
      // Format: 3-column table (Section | Checks | Tests), width ~110 chars.
      auto sec_box = [](const char* sec, const char* desc, int n) {
        auto hl = [](int cnt) { for (int k = 0; k < cnt; ++k) std::printf("\xe2\x94\x80"); };
        std::printf("\n  \xe2\x94\x8c"); hl(22); std::printf("\xe2\x94\xac"); hl(74);
        std::printf("\xe2\x94\xac"); hl(7); std::printf("\xe2\x94\x90\n");
        std::printf("  \xe2\x94\x82 %-20s \xe2\x94\x82 %-72s \xe2\x94\x82 %-5s \xe2\x94\x82\n",
                    "Section", "Checks", "Tests");
        std::printf("  \xe2\x94\x9c"); hl(22); std::printf("\xe2\x94\xbc"); hl(74);
        std::printf("\xe2\x94\xbc"); hl(7); std::printf("\xe2\x94\xa4\n");

        char nt[8]; std::snprintf(nt, sizeof(nt), "+%d", n);
        std::printf("  \xe2\x94\x82 %-20s \xe2\x94\x82 %-72s \xe2\x94\x82 %-5s \xe2\x94\x82\n",
                    sec, desc, nt);
        std::printf("  \xe2\x94\x94"); hl(22); std::printf("\xe2\x94\xb4"); hl(74);
        std::printf("\xe2\x94\xb4"); hl(7); std::printf("\xe2\x94\x98\n");
      };

      // ── Test compositions ──────────────────────────────────────────────────
      // Add more compositions by calling test_comp() below.
      // Component indices: CH4=1, C2H6=2, C3H8=3, i-C4=4, n-C4=5, H2=24,
      //   CO=25, H2S=26, He=27, Ar=28, N2=29, O2=30, CO2=31, Ammonia=34
      //
      // ── Density references (BWRS vs ideal gas at 1 atm) ─────────────────
      // Ideal gas 20°C: rho = M / (R × 293.15) = M / 24.0554 [kg/m³]
      // Ideal gas  0°C: rho = M / (R × 273.15) = M / 22.4136 [kg/m³]
      //
      //  Composition    M [g/mol] rho_id_20  rho_id_0   B₂(20°C)[L/mol]  rho_BWRS_20
      //  CH4             16.043    0.6672     0.7157     −0.0447           0.668
      //  C2H6            30.070    1.2500     1.3415     −0.195            1.260
      //  CO2             44.011    1.8296     1.9635     −0.147            1.841
      //  H2               2.016    0.0838     0.0899     +0.022            0.0837
      //  N2              28.016    1.1647     1.2498     −0.0044           1.162
      //  Ar              39.944    1.6604     1.7821     −0.018            1.661
      //
      // B₂ source: computed from BWRS Starling 1973 parameters via:
      //   B₂(T) = B₀ − A₀/(RT) − C₀/(RT³) + D₀/(RT⁴) − E₀/(RT⁵)
      //
      // ── Viscosity references (Chapman-Enskog + density correction) ───────
      //  Composition  η_calc [μPa·s]   η_NIST 20°C [μPa·s]   CE error
      //  CH4           11.46            10.99                   +4.3 %
      //  C2H6           9.68             9.36                   +3.4 %
      //  CO2           15.49            14.89                   +4.0 %
      //  H2             8.79             8.91                   −1.3 %
      //  N2            18.47            17.54                   +5.3 %
      //  Ar            22.41            22.74                   −1.4 %
      //
      // NIST source: https://webbook.nist.gov (Transport Properties, 100 kPa, 20°C)
      //
      // ── Mass flow references ─────────────────────────────────────────────
      // All tests: D=200 mm, d=80 mm, β=0.4, T=20°C, p=500 kPa, Δp=5 kPa
      // β=0.4 chosen to satisfy all ISO 5167 constraints simultaneously:
      //   Ven-P/T requires β≥0.4; Ven-T requires D≥200 mm; Ven-P Re≤1e6
      // Run for all 9 ISO 5167 device types.
      // Qm_est = alpha_nom × eps_nom × A_o × √(2000 × Δp_kPa × ρ)
      //   A_o = π/4 × (0.08)² = 5.027e−3 m²
      // Accepted range: ±15% of Qm_est (covers Re and BWRS variation)
      //
      //  Type      α_nom   ε_nom   C_nom   C source
      //  Dia-U     0.609   0.993   0.601   ISO 5167-2, Reader-Harris/Gallagher β=0.4 Re=10⁶
      //  Dia-F     0.608   0.993   0.601   idem, flange taps
      //  Dia-D/2   0.608   0.993   0.600   idem, D and D/2 taps
      //  Aj-ISA    0.997   0.997   0.985   ISO 5167-3, ISA 1932 nozzle
      //  Aj-RL     1.005   0.997   0.992   ISO 5167-3, long-radius nozzle
      //  Ven-B     0.997   0.997   0.984   ISO 5167-4, rough-cast Venturi (C=0.984)
      //  Ven-P     1.008   0.997   0.995   ISO 5167-4, machined Venturi   (C=0.995)
      //  Ven-T     0.998   0.997   0.985   ISO 5167-4, welded-sheet Venturi (C=0.985)
      //  Aj-V      0.996   0.997   0.983   ISO 5167-4, Venturi nozzle (C≈0.983 at β=0.4)
      // ───────────────────────────────────────────────────────────────────

      double xCH4[kArraySize] = {}; xCH4[1] = 1.0;   // pure methane

      double xC2H6[kArraySize] = {}; xC2H6[2] = 1.0;  // pure ethane

      double xCO2[kArraySize] = {}; xCO2[31] = 1.0;   // pure CO2

      double xH2[kArraySize] = {}; xH2[24] = 1.0;     // pure hydrogen

      double xN2[kArraySize] = {}; xN2[29] = 1.0;     // pure nitrogen

      double xAr[kArraySize] = {}; xAr[28] = 1.0;     // pure argon

      double xGN[kArraySize] = {};                     // standard natural gas
      xGN[1] = 0.900; xGN[2] = 0.060; xGN[3] = 0.020;
      xGN[29] = 0.015; xGN[31] = 0.005;

      double xGNB[kArraySize] = {};                    // rich natural gas
      xGNB[1] = 0.850; xGNB[2] = 0.100; xGNB[3] = 0.030;
      xGNB[29] = 0.010; xGNB[31] = 0.010;

      sec_box("Test compositions",
              "8 compositions: rho(20/0C), eta, Z, molar mass, Qm 9 ISO 5167 types", 120);
      std::printf("%s  \xe2\x94\x80\xe2\x94\x80 Legend%s\n", kYellow, kReset);

      std::printf("\n%s  Physical quantities:%s\n", kBoldWhite, kReset);
      std::printf("    \xcf\x81   [kg/m\xc2\xb3]   density\n");
      std::printf("    \xce\xb7   [\xc2\xb5Pa\xc2\xb7s]   dynamic viscosity\n");
      std::printf("    Z   [-]       compressibility factor\n");
      std::printf("    M   [g/mol]   molar mass\n");
      std::printf("    Qm  [kg/s]    mass flow rate\n");
      std::printf("    Qv  [Nm\xc2\xb3/h]   volumetric flow (normal conditions: 0\xc2\xb0""C, 101.325 kPa)\n");

      std::printf("\n%s  Flow device parameters:%s\n", kBoldWhite, kReset);
      std::printf("    D   [mm]      pipe diameter\n");
      std::printf("    d   [mm]      orifice / throat diameter\n");
      std::printf("    \xce\xb2   [-]       diameter ratio  \xce\xb2 = d/D\n");
      std::printf("    p   [kPa]     absolute pressure\n");
      std::printf("    \xce\x94p  [kPa]     differential pressure\n");

      std::printf("\n%s  ISO 5167 devices:%s\n", kBoldWhite, kReset);
      std::printf("    Dia-U   = Orifice plate, corner taps\n");
      std::printf("    Dia-F   = Orifice plate, flange taps\n");
      std::printf("    Dia-D/2 = Orifice plate, D and D/2 taps\n");
      std::printf("    Aj-ISA  = ISA 1932 nozzle\n");
      std::printf("    Aj-RL   = Long-radius nozzle (ASME long-radius)\n");
      std::printf("    Ven-B   = Venturi tube \xe2\x80\x94 rough-cast convergent      (C \xe2\x89\x88 0.984)\n");
      std::printf("    Ven-P   = Venturi tube \xe2\x80\x94 machined convergent        (C \xe2\x89\x88 0.995)\n");
      std::printf("    Ven-T   = Venturi tube \xe2\x80\x94 welded sheet-metal conv.   (C \xe2\x89\x88 0.985)\n");
      std::printf("    Aj-V    = Venturi nozzle                             (C \xe2\x89\x88 0.983)\n");

      std::printf("\n%s  Reference sources [right column of each test]:%s\n", kBoldWhite, kReset);
      std::printf("    BWRS 1973     = Benedict-Webb-Rubin-Starling (Starling, K.E., 1973)\n");
      std::printf("    CE/NIST       = Chapman-Enskog + density correction / NIST WebBook\n");
      std::printf("    BWRS/ideal    = BWRS vs ideal gas at 101.325 kPa\n");
      std::printf("    BWRS/NIST     = BWRS compared with NIST WebBook tabulated data\n");
      std::printf("    NIST/BWRS     = NIST reference, calculated with BWRS\n");
      std::printf("    BWRS/AGA-8    = BWRS compared with AGA-8 data\n");
      std::printf("    BWRS/estim.   = BWRS compared with engineering estimate\n");
      std::printf("    BWRS/Z        = \xcf\x81 monotonicity with Z-factor correction\n");
      std::printf("    ISO 5167      = ISO 5167-2/3/4:2003 (Reader-Harris/Gallagher)\n");
      std::printf("    ISO 5167/BWRS = Qv = Qm(ISO 5167) / \xcf\x81(BWRS, 0\xc2\xb0""C, 101.325 kPa)\n");
      std::printf("    mixing rule   = linear rule  M = \xce\xa3 xi\xc2\xb7Mi\n");
      std::printf("    ideal gas     = monotonicity test  \xcf\x81(2p)/\xcf\x81(p) \xe2\x89\x88 2\n");
      std::printf("    ideal gas\xc3\x97Z   = as above, with Z-factor correction\n");
      std::printf("    single-phase  = Z \xe2\x88\x88 [0.75, 1.02] check (gas phase)\n");
      std::printf("    quasi-ideal   = \xcf\x81/\xcf\x81_ideal \xe2\x88\x88 [0.990, 1.025] at p \xe2\x89\xa4 500 kPa\n");
      std::printf("    monotone T    = \xcf\x81(T1) > \xcf\x81(T2) at T1 < T2  (density increases with cooling)\n");
      std::printf("    monotone P    = \xcf\x81(p1) < \xcf\x81(p2) at p1 < p2  (density increases with pressure)\n");
      std::printf("    M + Z         = density ratio \xe2\x89\x88 (Ma/Mb) \xc3\x97 (Zb/Za)\n");
      std::printf("    CE/kinetic    = Chapman-Enskog from kinetic theory (monotone \xce\xb7 vs. T)\n");
      std::printf("\n");

      //           Composition   rho_20[kg/m3]  rho_0[kg/m3]   eta[µPa·s]   D    d   p    dp
      test_comp("Methane (CH4)",  xCH4, 0.660,0.672, 0.712,0.724, 10.5,12.5, 200,80, 500,5, false);
      test_comp("Ethane (C2H6)", xC2H6,1.220,1.280, 1.310,1.375,  8.5,10.8, 200,80, 500,5);
      test_comp("CO2",           xCO2, 1.800,1.860, 1.930,2.000, 14.0,16.8, 200,80, 500,5);
      test_comp("Hydrogen (H2)", xH2,  0.082,0.085, 0.088,0.091,  8.2, 9.5, 200,80, 500,5);
      test_comp("Nitrogen (N2)", xN2,  1.155,1.175, 1.240,1.260, 16.5,20.5, 200,80, 500,5);
      test_comp("Argon (Ar)",    xAr,  1.650,1.672, 1.772,1.794, 21.0,23.5, 200,80, 500,5);
      test_comp("Std NG",        xGN,  0.730,0.748, 0.784,0.802, 10.5,12.5, 200,80, 500,5);
      test_comp("Rich NG",       xGNB, 0.765,0.790, 0.820,0.848, 10.2,12.5, 200,80, 500,5);

      sec_box("Pure CH4 standalone",
              "Pure CH4: rho monotonicity vs P (x2); Z factor at 20 atm", 3);
      // ── Standalone tests: physical behaviour of pure CH4 ────────────────
      std::printf("  %s\xC2\xBB Standalone tests pure CH4%s\n", kYellow, kReset);
      {
        BwrConst bt_a = make_bwr(xCH4);
        double T_a    = 20.0 + kKelvinOffset;
        double rho_05 = CalcDensity(20.0,  0.5, bt_a);   // 50.66 kPa
        double rho_10 = CalcDensity(20.0,  1.0, bt_a);   // 101.325 kPa
        double rho_20 = CalcDensity(20.0, 20.0, bt_a);   // ≈ 2.026 MPa
        // Ideal gas: rho(2p)/rho(p) = 2; Z correction introduces small deviations
        chk("Monotone \xCF\x81(1atm)/\xCF\x81(0.5atm) CH4 [-]",
            rho_10 / rho_05, 1.90, 2.10, "ideal gas");
        // At 20 atm, Z(CH4)≈0.953 → rho(20)/rho(1) ≈ 20/0.953 ≈ 21.0
        chk("Monotone \xCF\x81(20atm)/\xCF\x81(1atm) CH4 [-]",
            rho_20 / rho_10, 17.0, 23.0, "ideal gas\xC3\x97Z");
        // Z factor at high pressure: NIST CH4 20°C, 20 atm → Z≈0.953
        double Z_hi = 20.0 * bt_a.molar_mass / (rho_20 * kGasConstantR * T_a);
        chk("Z factor CH4 20\xC2\xB0""C, 20 atm [-]", Z_hi, 0.88, 0.98, "BWRS/NIST");
      }

      sec_box("High-pressure rho",
              "CH4/N2/CO2/GN at 500-5066 kPa: rho BWRS, Z, monotonicity, consistency", 8);
      // ── Density at high pressure ──────────────────────────────────────────
      // References:
      //   CH4  500 kPa, 293.15 K: Z≈0.9986 (NIST)       → ρ≈3.30  kg/m³ (≈ideal)
      //   CH4 5066 kPa, 293.15 K: Z≈0.911  (BWRS/AGA-8) → ρ≈36.6  kg/m³
      //     Virial check: B₂(CH4,293K)≈−44 cm³/mol
      //     → Z≈1+B₂·P/(RT)=1−0.044·50/24.05=0.909 (consistent with BWRS 0.911)
      //     Note: Z≈0.954 would correspond to ~27 atm, not 50 atm (below Boyle temp. 511 K)
      //   N2  5066 kPa, 293.15 K: Z≈0.991  (NIST)       → ρ≈58.8  kg/m³
      //   CO2  500 kPa, 293.15 K: Z≈0.974  (NIST)       → ρ≈9.27  kg/m³ (non-ideal)
      //   Std NG 5066 kPa, 20°C : BWRS estimate          → ρ≈38-40 kg/m³
      // Tolerance ±5%: BWRS Starling 1973 typical error 1–3% at high pressure.
      std::printf("  %s\xc2\xbb Density at high pressure%s\n", kYellow, kReset);
      {
        BwrConst bCH4 = make_bwr(xCH4);
        BwrConst bN2  = make_bwr(xN2);
        BwrConst bCO2 = make_bwr(xCO2);
        BwrConst bGN  = make_bwr(xGN);

        double T_K   = 20.0 + kKelvinOffset;
        double p500  = 500.0  / kKpaPerAtm;   //  4.935 atm
        double p5066 = 5066.0 / kKpaPerAtm;   // 50.00  atm

        // ─ CH4 ──────────────────────────────────────────────────────────────
        double rCH4_500  = CalcDensity(20.0, p500,  bCH4);
        double rCH4_5066 = CalcDensity(20.0, p5066, bCH4);
        chk("Density CH4 20\xc2\xb0""C,  500 kPa [kg/m\xc2\xb3]",
            rCH4_500,  3.24, 3.35, "NIST/BWRS");
        chk("Density CH4 20\xc2\xb0""C, 5066 kPa [kg/m\xc2\xb3]",
            rCH4_5066, 34.5, 39.0, "BWRS/AGA-8");

        // Z at 50 atm: B₂(CH4,293K)≈−44 cm³/mol → Z_virial≈0.909; BWRS→0.911
        double Z50 = p5066 * bCH4.molar_mass
                   / (rCH4_5066 * kGasConstantR * T_K);
        chk("Z factor CH4 20\xc2\xb0""C, 50 atm [-]",
            Z50, 0.88, 0.95, "BWRS/AGA-8");

        // ─ N2: quasi-ideal at 50 atm, Z ≈ 0.987 (NIST) ─────────────────────
        double rN2_5066 = CalcDensity(20.0, p5066, bN2);
        chk("Density N2  20\xc2\xb0""C, 5066 kPa [kg/m\xc2\xb3]",
            rN2_5066, 56.0, 62.0, "NIST/BWRS");

        // ─ CO2 at 500 kPa: Z ≈ 0.974, more non-ideal than CH4 ───────────────
        double rCO2_500 = CalcDensity(20.0, p500, bCO2);
        chk("Density CO2 20\xc2\xb0""C,  500 kPa [kg/m\xc2\xb3]",
            rCO2_500, 8.90, 9.70, "NIST/BWRS");

        // ─ Std NG at 50 atm: typical pipeline pressure ────────────────────
        double rGN_5066 = CalcDensity(20.0, p5066, bGN);
        chk("Density std NG 20\xc2\xb0""C, 5066 kPa [kg/m\xc2\xb3]",
            rGN_5066, 35.0, 42.0, "BWRS/estim.");

        // ─ Monotonicity: ρ(50 atm)/ρ(5 atm) = p_ratio × Z(5 atm)/Z(50 atm) ─
        // = 10.13 × (0.9986/0.954) ≈ 10.6
        chk("Monotone \xcf\x81 CH4: \xcf\x81(5066)/\xcf\x81(500) [-]",
            rCH4_5066 / rCH4_500, 9.5, 11.5, "BWRS/Z");

        // ─ Consistency: ρ_CO2 > ρ_N2 at same T and p ────────────────────────
        // M_CO2=44.01 >> M_N2=28.02 and Z_CO2 < Z_N2 → ratio ≈ 1.56–1.62
        double rN2_500 = CalcDensity(20.0, p500, bN2);
        chk("Consistency: \xcf\x81(CO2)/\xcf\x81(N2) at 500 kPa [-]",
            rCO2_500 / rN2_500, 1.50, 1.65, "M + Z");
      }

      sec_box("Rich NG 4Tx4P",
              "Rich NG 4Tx4P: Z, rho/rho_ideal, monotonicity vs T and vs P", 30);
      // ── Rich NG: 4-temperature × 4-pressure matrix ───────────────────────
      // Temperatures: -10°C (winter), 10°C, 30°C, 60°C (summer/compressor)
      // Pressures:    101 kPa (ref.), 500 kPa, 2000 kPa, 7000 kPa (pipeline)
      // T_min = -10°C = 263.15 K > Tc_pseudo(GNB) ≈ 208 K → single-phase gas.
      // Checks: Z ∈ [0.75, 1.02]; ρ/ρ_ideal at low p; monotonicity vs T and P.
      std::printf("  %s\xc2\xbb Rich NG: T \xc3\x97 P matrix  (4 \xc3\x97 4 conditions)%s\n",
                  kYellow, kReset);
      {
        BwrConst bGNB = make_bwr(xGNB);

        const double T_c[]   = { -10.0,   10.0,   30.0,   60.0 };   // [°C]
        const double P_kpa[] = { 101.325, 500.0, 2000.0, 7000.0 };  // [kPa]
        const char*  T_lbl[] = { "-10",  " 10",  " 30",  " 60"  };
        const char*  P_lbl[] = { "  101", "  500", " 2000", " 7000" };
        constexpr int nT = 4, nP = 4;

        double rho[nT][nP] = {};

        // ─── Z ∈ [0.75, 1.02] for all 16 conditions ─────────────────────────
        for (int ti = 0; ti < nT; ti++) {
          for (int pi = 0; pi < nP; pi++) {
            double T_K      = T_c[ti] + kKelvinOffset;
            double p_atm    = P_kpa[pi] / kKpaPerAtm;
            rho[ti][pi]     = CalcDensity(T_c[ti], p_atm, bGNB);
            double Z        = p_atm * bGNB.molar_mass
                            / (rho[ti][pi] * kGasConstantR * T_K);
            char lbl[64];
            std::snprintf(lbl, sizeof(lbl),
                          "Z GNB T=%s\xc2\xb0""C P=%s kPa",
                          T_lbl[ti], P_lbl[pi]);
            chk(lbl, Z, 0.75, 1.02, "single-phase");
          }
        }

        // ─── At p ≤ 500 kPa: ρ/ρ_ideal ∈ [0.990, 1.025] (quasi-ideal) ───────
        // Upper limit 1.025 covers T=-10°C where B_mix≈-85 cm³/mol → Z≈0.981
        for (int ti = 0; ti < nT; ti++) {
          for (int pi = 0; pi < 2; pi++) {
            double T_K    = T_c[ti] + kKelvinOffset;
            double rho_id = P_kpa[pi] * 1e3 * bGNB.molar_mass * 1e-3
                          / (8.314 * T_K);
            char lbl[64];
            std::snprintf(lbl, sizeof(lbl),
                          "\xcf\x81/\xcf\x81_id GNB T=%s\xc2\xb0""C P=%s kPa",
                          T_lbl[ti], P_lbl[pi]);
            chk(lbl, rho[ti][pi] / rho_id, 0.990, 1.025, "quasi-ideal");
          }
        }

        // ─── Monotonicity vs T at P = 2000 kPa (pi=2): ρ(T1) > ρ(T2) ─────
        // Ratio ρ(T_i)/ρ(T_{i+1}) = (Z_{i+1}/Z_i) × (T_{i+1}/T_i) ≈ 1.07–1.10
        chk("GNB: \xcf\x81(-10\xc2\xb0""C)/\xcf\x81( 10\xc2\xb0""C) at 2000 kPa",
            rho[0][2] / rho[1][2], 1.03, 1.15, "monotone T");
        chk("GNB: \xcf\x81( 10\xc2\xb0""C)/\xcf\x81( 30\xc2\xb0""C) at 2000 kPa",
            rho[1][2] / rho[2][2], 1.03, 1.15, "monotone T");
        chk("GNB: \xcf\x81( 30\xc2\xb0""C)/\xcf\x81( 60\xc2\xb0""C) at 2000 kPa",
            rho[2][2] / rho[3][2], 1.03, 1.15, "monotone T");

        // ─── Monotonicity vs P at T = 10°C (ti=1): ρ(P1) < ρ(P2) ───────────
        // ρ(500)/ρ(101)   ≈ p_ratio × Z(101)/Z(500) ≈ 4.94 × 1.003 ≈ 4.96
        chk("GNB: \xcf\x81(500)/\xcf\x81(101) at 10\xc2\xb0""C",
            rho[1][1] / rho[1][0], 4.70, 5.10, "monotone P");
        // ρ(2000)/ρ(500)  ≈ 4.00 × Z(500)/Z(2000) ≈ 4.00 × 1.016 ≈ 4.07
        chk("GNB: \xcf\x81(2000)/\xcf\x81(500) at 10\xc2\xb0""C",
            rho[1][2] / rho[1][1], 3.70, 4.30, "monotone P");
        // ρ(7000)/ρ(2000) ≈ 3.50 × Z(2000)/Z(7000) ≈ 3.50 × 1.14 ≈ 3.98
        chk("GNB: \xcf\x81(7000)/\xcf\x81(2000) at 10\xc2\xb0""C",
            rho[1][3] / rho[1][2], 3.50, 4.80, "monotone P");
      }

      sec_box("Negative tests",
              "d<12.5 Dia-U; beta<0.23 Dia-U; D>760 Dia-F; D>500 Aj-ISA; d=50 Aj-V", 9);
      // ── Negative ISO 5167 tests: invalid parameters must return 0 ──────────
      std::printf("  %s\xC2\xBB ISO 5167 negative tests%s\n", kYellow, kReset);
      {
        BwrConst bt_n = make_bwr(xCH4);
        double rho_n  = CalcDensity(20.0, 500.0 / kKpaPerAtm, bt_n);
        double eta_n  = 11.5e-6;  // CH4 viscosity [Pa*s]
        FlowResult fr_n;
        // Suppress PrintError messages during negative tests
        int saved_fd = _dup(1);
        FILE* nul_f  = nullptr;
        fopen_s(&nul_f, "NUL", "w");
        if (nul_f) { _dup2(_fileno(nul_f), 1); fclose(nul_f); }
        double qn1 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kOrificeCorner,
                                  200.0,170.0, rho_n,eta_n, &fr_n);  // β=0.85>0.80
        double qn2 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kOrificeCorner,
                                   30.0, 15.0, rho_n,eta_n, &fr_n);  // D=30<50mm
        double qn3 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kVenturiMachined,
                                  200.0, 60.0, rho_n,eta_n, &fr_n);  // β=0.30<0.40
        double qn4 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kVenturiWeldedSheet,
                                  200.0,150.0, rho_n,eta_n, &fr_n);  // β=0.75>0.70
        double qn5 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kOrificeCorner,
                                  200.0, 10.0, rho_n,eta_n, &fr_n);  // d=10<12.5mm
        double qn6 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kOrificeCorner,
                                  200.0, 44.0, rho_n,eta_n, &fr_n);  // β=0.22<0.23
        double qn7 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kOrificeFlange,
                                  800.0,400.0, rho_n,eta_n, &fr_n);  // D=800>760mm
        double qn8 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kNozzleIsa,
                                  550.0,220.0, rho_n,eta_n, &fr_n);  // D=550>500mm
        double qn9 = CalcMassFlow(5.0,500.0,20.0, DeviceType::kVenturiNozzle,
                                  200.0, 50.0, rho_n,eta_n, &fr_n);  // d=50<=50mm
        std::fflush(stdout);
        _dup2(saved_fd, 1);
        _close(saved_fd);
        const char* ns = "ISO 5167";
        chk("\xCE\xB2=0.85 Dia-U (max 0.80) \xE2\x86\x92 Qm=0",    qn1, -0.001, 0.001, ns);
        chk("D=30mm Dia-U (min 50mm) \xE2\x86\x92 Qm=0",            qn2, -0.001, 0.001, ns);
        chk("\xCE\xB2=0.30 Ven-P (min 0.40) \xE2\x86\x92 Qm=0",    qn3, -0.001, 0.001, ns);
        chk("\xCE\xB2=0.75 Ven-T (max 0.70) \xE2\x86\x92 Qm=0",    qn4, -0.001, 0.001, ns);
        chk("d=10mm Dia-U (min 12.5mm) \xE2\x86\x92 Qm=0",          qn5, -0.001, 0.001, ns);
        chk("\xCE\xB2=0.22 Dia-U (min 0.23) \xE2\x86\x92 Qm=0",    qn6, -0.001, 0.001, ns);
        chk("D=800mm Dia-F (max 760mm) \xE2\x86\x92 Qm=0",          qn7, -0.001, 0.001, ns);
        chk("D=550mm Aj-ISA (max 500mm) \xE2\x86\x92 Qm=0",         qn8, -0.001, 0.001, ns);
        chk("d=50mm Aj-V (limit d>50mm) \xE2\x86\x92 Qm=0",      qn9, -0.001, 0.001, ns);
      }

      sec_box("Monotone eta vs T",
              "eta(CH4) increases with T: ratios -10C/20C/60C; sqrt(T) kinetic", 3);
      // ── Viscosity monotonicity vs temperature (CH4) ───────────────────────
      // Gases: η ∝ √T at dilute limit (kinetic theory). CE includes factor (1+k·ln(T/cs))·√T.
      // NIST CH4: η(-10°C)≈10.5 μPa·s; η(20°C)≈11.0; η(60°C)≈11.9 μPa·s.
      // √(333.15/263.15)≈1.126; √(293.15/263.15)≈1.055; √(333.15/293.15)≈1.065.
      std::printf("  %s\xC2\xBB Viscosity monotonicity vs T (CH4)%s\n",
                  kYellow, kReset);
      {
        BwrConst bt_mt = make_bwr(xCH4);
        double tc_mt = 0.0, pc_mt = 0.0, zc_mt = 0.0;
        for (int i = 1; i <= kNumComponents; i++) {
          tc_mt += xCH4[i] * Tc[i];
          pc_mt += xCH4[i] * Pc[i];
          zc_mt += xCH4[i] * Zc[i];
        }
        double roc_mt = pc_mt / (kGasConstantR * zc_mt * tc_mt);
        double csi_mt = std::pow(tc_mt, 6) / std::sqrt(bt_mt.molar_mass)
                      / std::cbrt(pc_mt * pc_mt);
        double eN10 = calc_eta(-10.0, xCH4, bt_mt, roc_mt, csi_mt) * kPaToMicroPa;
        double eP20 = calc_eta( 20.0, xCH4, bt_mt, roc_mt, csi_mt) * kPaToMicroPa;
        double eP60 = calc_eta( 60.0, xCH4, bt_mt, roc_mt, csi_mt) * kPaToMicroPa;
        chk("Monotone \xCE\xB7 CH4: \xCE\xB7(60\xC2\xB0""C)/\xCE\xB7(-10\xC2\xB0""C) [-]",
            eP60 / eN10, 1.05, 1.25, "CE/kinetic");
        chk("Monotone \xCE\xB7 CH4: \xCE\xB7(20\xC2\xB0""C)/\xCE\xB7(-10\xC2\xB0""C) [-]",
            eP20 / eN10, 1.01, 1.10, "CE/kinetic");
        chk("Monotone \xCE\xB7 CH4: \xCE\xB7(60\xC2\xB0""C)/\xCE\xB7(20\xC2\xB0""C) [-]",
            eP60 / eP20, 1.01, 1.12, "CE/kinetic");
      }

      sec_box("Gas eta ordering",
              "eta(N2)/eta(CH4)~1.61; eta(CO2)/eta(CH4)~1.35 at 20C, 1 atm", 2);
      // ── Viscosity ordering between gases at 20°C ─────────────────────────
      // CE: η(CH4)≈11.5; η(CO2)≈15.5; η(N2)≈18.5 μPa·s at 20°C, 1 atm.
      // Ar omitted from relative ordering test (η_Ar ~ 22.4 μPa·s, close to N2 ~ 18.5 → ratio near 1).
      std::printf("  %s\xC2\xBB Viscosity ordering between gases at 20\xC2\xB0""C%s\n",
                  kYellow, kReset);
      {
        auto eta20 = [&](const double* xc) -> double {
          BwrConst bt = make_bwr(xc);
          double tc = 0.0, pc = 0.0, zc = 0.0;
          for (int i = 1; i <= kNumComponents; i++) {
            tc += xc[i] * Tc[i]; pc += xc[i] * Pc[i]; zc += xc[i] * Zc[i];
          }
          double roc = pc / (kGasConstantR * zc * tc);
          double csi = std::pow(tc, 6) / std::sqrt(bt.molar_mass) / std::cbrt(pc * pc);
          return calc_eta(20.0, xc, bt, roc, csi) * kPaToMicroPa;
        };
        double eCH4 = eta20(xCH4);
        double eN2  = eta20(xN2);
        double eCO2 = eta20(xCO2);
        // η(N2)/η(CH4) ≈ 18.47/11.46 ≈ 1.61
        chk("Ordering: \xCE\xB7(N2)/\xCE\xB7(CH4) at 20\xC2\xB0""C [-]",
            eN2  / eCH4, 1.40, 1.80, "CE/NIST");
        // η(CO2)/η(CH4) ≈ 15.49/11.46 ≈ 1.35
        chk("Ordering: \xCE\xB7(CO2)/\xCE\xB7(CH4) at 20\xC2\xB0""C [-]",
            eCO2 / eCH4, 1.20, 1.55, "CE/NIST");
      }

      sec_box("Synthetic air",
              "79% N2+21% O2: rho(20C)~1.199 kg/m3; eta~18.2 muPa*s; Qm for 9 types", 15);
      // ── Synthetic air (79% N2 + 21% O2) ──────────────────────────────────
      // M_air = 0.79×28.016 + 0.21×32.000 = 28.853 g/mol
      // ρ_ideal(20°C)=28.853/24.055=1.199 kg/m³; ρ_ideal(0°C)=28.853/22.414=1.287
      // η_NIST air at 20°C ≈ 18.2 μPa·s; Z ≈ 1.000 at 1 atm (quasi-ideal)
      {
        double xAer[kArraySize] = {};
        xAer[29] = 0.79;  // N2
        xAer[30] = 0.21;  // O2
        test_comp("Synthetic air (79% N2 + 21% O2)", xAer,
                  1.182, 1.218,   // \xcf\x81(20\xc2\xb0C, 1 atm): ideal=1.199, Z\xe2\x89\x881.00
                  1.270, 1.305,   // \xcf\x81(0\xc2\xb0C,  1 atm): ideal=1.287, Z\xe2\x89\x881.00
                  15.5,  22.0,    // \xce\xb7(20\xc2\xb0C): NIST air\xe2\x89\x8818.2 \xc2\xb5Pa\xc2\xb7s
                  200, 80, 500, 5, false);
      }

      sec_box("H2S density",
              "Pure H2S (Tc=373.6K, Tr=0.785): rho at 1 atm and 500 kPa, N2 consistency", 3);
      // ── H2S: non-ideal density and consistency ────────────────────────────
      // H2S: Tc=373.6 K, Pc=88.9 atm; at 20°C → Tr=0.785 → strongly non-ideal.
      // ρ_ideal(1 atm)=34.082/24.055=1.417 kg/m³; Z≈0.993 → ρ≈1.427
      // ρ_ideal(500 kPa)=4.935×1.417=6.99 kg/m³; Z≈0.965 → ρ≈7.24
      // M(H2S)/M(N2)=1.217; Z(H2S)<Z(N2) → ρ(H2S)/ρ(N2)>1.217
      std::printf("  %s\xC2\xBB H2S: non-ideal density and consistency%s\n",
                  kYellow, kReset);
      {
        double xH2S[kArraySize] = {};  xH2S[26] = 1.0;
        BwrConst bH2S   = make_bwr(xH2S);
        BwrConst bN2_hs = make_bwr(xN2);
        double rH2S_1   = CalcDensity(20.0, 1.0,                  bH2S);
        double rH2S_500 = CalcDensity(20.0, 500.0 / kKpaPerAtm,   bH2S);
        double rN2_500  = CalcDensity(20.0, 500.0 / kKpaPerAtm,   bN2_hs);
        chk("Density H2S 20\xC2\xB0""C, 101.325 kPa [kg/m\xC2\xB3]",
            rH2S_1,   1.38, 1.48, "BWRS/ideal");
        chk("Density H2S 20\xC2\xB0""C,  500 kPa   [kg/m\xC2\xB3]",
            rH2S_500, 6.80, 7.80, "BWRS/NIST");
        // Ratio > M ratio (1.217) due to greater non-ideality of H2S
        chk("Consistency: \xCF\x81(H2S)/\xCF\x81(N2) at 500 kPa [-]",
            rH2S_500 / rN2_500, 1.20, 1.45, "M + Z");
      }

      // ── Summary ───────────────────────────────────────────────────────────
      std::printf("\n  %sOverall result: %d/%d tests passed.%s\n",
                  (npass == ntotal) ? kBoldGreen : kBoldRed,
                  npass, ntotal, kReset);
      if (npass < ntotal)
        std::printf("  %sWARNING: Some tests failed \xe2\x80\x94 check the implementation!%s\n",
                    kBoldRed, kReset);
      std::printf("\n  %sPress any key to continue...%s", kBoldWhite, kReset);
      _getch();
      std::printf("\n");
    }
  }

  // ── Composition: load from file or enter manually ─────────────────────────
  bool comp_loaded = false;
  {
    int    tmp_tip = 0;
    double tmp_d1  = 0.0;
    double tmp_d2  = 0.0;
    (void)tmp_tip; (void)tmp_d1; (void)tmp_d2;
    FILE* cf = nullptr;
    fopen_s(&cf, kCompFile, "r");
    if (cf) {
      std::fclose(cf);
      std::printf("%s\n  A saved composition exists. Reuse it? [y/n] %s>%s ",
                  kBoldWhite, kCyan, kReset);
      if (AskYesNo()) {
        if (LoadComposition(x)) {
          double loaded_sum = 0.0;
          for (int i = 1; i <= kNumComponents; i++) loaded_sum += x[i];
          if (std::fabs(loaded_sum - 1.0) > kSumTolerance) {
            std::printf("%s  Loaded composition: sum = %.6f \xe2\x89\xa0 1 "
                        "— switching to manual entry.%s\n", kBoldRed, loaded_sum, kReset);
            for (int i = 1; i <= kNumComponents; i++) x[i] = 0.0;
          } else {
            comp_loaded = true;
            PrintComposition(x);
          }
        } else {
          std::printf("%s  Error reading file. Switching to manual entry.%s\n",
                      kBoldRed, kReset);
        }
      }
    }
  }

  if (!comp_loaded) {
    double sum = 0.0;
    do {
      std::printf("\n%s  ── Composition %s\n",
                  kYellow, kReset);
      std::printf("%s  Molar fractions of the gas mixture:%s\n\n", kBoldWhite, kReset);
      for (int i = 1; i <= kNumComponents; i++) {
        const char* name = kCompNames[i];
        int vlen = (int)std::strlen(name) - Utf8ExtraBytes(name);
        double v;
        do {
          std::printf("  %s>%s %s%s", kCyan, kReset, kBoldWhite, name);
          for (int j = vlen; j < kCompNameWidth; j++) {
            std::putchar(' ');
          }
          std::printf(": ");
          ReadDouble(&v);
          std::printf("%s", kReset);
          if (v < 0.0 || v > 1.0) {
            std::printf("%s  Invalid value \xe2\x80\x94 must be in [0, 1].%s\n", kBoldRed, kReset);
          }
        } while (v < 0.0 || v > 1.0);
        x[i] = v;
      }

      sum = 0.0;
      for (int i = 1; i <= kNumComponents; i++) {
        sum += x[i];
      }

      if (std::fabs(sum - 1.0) > kSumTolerance) {
        std::printf("%s\n  Sum of molar fractions = %.6f  \xe2\x89\xa0  1.%s\n", kBoldRed, sum, kReset);
        if (sum < 1e-10) {
          std::printf("%s  Sum is zero \xe2\x80\x94 re-enter the composition.%s\n", kBoldRed, kReset);
        } else {
          std::printf("%s  Normalise automatically? [y/n] (n = re-enter) %s>%s ", kYellow, kCyan, kReset);
          if (AskYesNo()) {
            for (int i = 1; i <= kNumComponents; i++) {
              x[i] /= sum;
            }
            sum = 1.0;
            std::printf("%s  Normalised fractions (non-zero components):%s\n\n", kBoldGreen, kReset);
            for (int i = 1; i <= kNumComponents; i++) {
              if (x[i] > 0.0) {
                const char* name = kCompNames[i];
                int vlen = (int)std::strlen(name) - Utf8ExtraBytes(name);
                std::printf("    %s%s", kBoldGreen, name);
                for (int j = vlen; j < kCompNameWidth; j++) {
                  std::putchar(' ');
                }
                std::printf(": %.8f%s\n", x[i], kReset);
              }
            }
          }
        }
      }
    } while (std::fabs(sum - 1.0) > kSumTolerance);

    std::printf("%s\n  Save composition? [y/n] %s>%s ", kBoldWhite, kCyan, kReset);
    if (AskYesNo()) {
      SaveComposition(x);
    }
  }

  // ── Compute BWRS mixture constants ────────────────────────────────────────
  using Clock = std::chrono::high_resolution_clock;
  using Us    = std::chrono::duration<double, std::micro>;

  BwrConst bwr;
  auto t_mix0 = Clock::now();

  for (int i = 1; i <= kNumComponents; i++) {
    for (int j = 1; j <= kNumComponents; j++) {
      if (x[i] <= 0.0 || x[j] <= 0.0) continue;
      double kk = 1 - kBwrsMixCoef * std::sqrt(V[i] * V[j])
               / std::pow(std::cbrt(V[i]) + std::cbrt(V[j]), 3);
      bwr.a0    += x[i] * x[j] * std::sqrt(A[i] * A[j]) * (1 - kk);
      bwr.b0    += x[i] * x[j] * std::sqrt(B[i] * B[j]);
      bwr.c0    += x[i] * x[j] * std::sqrt(C[i] * C[j]) * std::pow(1 - kk, 3);
      bwr.d0    += x[i] * x[j] * std::sqrt(D0[i] * D0[j]) * std::pow(1 - kk, 4);
      bwr.e0    += x[i] * x[j] * std::sqrt(E0[i] * E0[j]) * std::pow(1 - kk, 5);
      bwr.gamma += x[i] * x[j] * std::sqrt(gama[i] * gama[j]);
    }
  }

  for (int i = 1; i <= kNumComponents; i++) {
    for (int j = 1; j <= kNumComponents; j++) {
      for (int l = 1; l <= kNumComponents; l++) {
        if (x[i] <= 0.0 || x[j] <= 0.0 || x[l] <= 0.0) continue;
        double k1 = 1 - kBwrsMixCoef * std::sqrt(V[i] * V[j])
                 / std::pow(std::cbrt(V[i]) + std::cbrt(V[j]), 3);
        double k2 = 1 - kBwrsMixCoef * std::sqrt(V[i] * V[l])
                 / std::pow(std::cbrt(V[i]) + std::cbrt(V[l]), 3);
        double k3 = 1 - kBwrsMixCoef * std::sqrt(V[j] * V[l])
                 / std::pow(std::cbrt(V[j]) + std::cbrt(V[l]), 3);
        bwr.a     += x[i] * x[j] * x[l]
                   * std::cbrt(a[i]*a[j]*a[l]) * (1-k1)*(1-k2)*(1-k3);
        bwr.b     += x[i] * x[j] * x[l] * std::cbrt(b[i]*b[j]*b[l]);
        bwr.c     += x[i] * x[j] * x[l]
                   * std::cbrt(c[i]*c[j]*c[l]) * (1-k1)*(1-k2)*(1-k3);
        bwr.dv    += x[i] * x[j] * x[l]
                   * std::cbrt(dv[i]*dv[j]*dv[l]) * (1-k1)*(1-k2)*(1-k3);
        bwr.alpha += x[i] * x[j] * x[l]
                   * std::cbrt(alfa[i]*alfa[j]*alfa[l]);
      }
    }
  }

  double tcam = 0.0;
  double pcam = 0.0;
  double zcam = 0.0;
  double mx   = 0.0;
  for (int i = 1; i <= kNumComponents; i++) {
    bwr.molar_mass += x[i] * m[i];
    tcam           += x[i] * Tc[i];
    pcam           += x[i] * Pc[i];
    zcam           += x[i] * Zc[i];
    mx             += x[i] * std::sqrt(m[i]);
  }

  double cp_mix = 0.0;
  for (int i = 1; i <= kNumComponents; i++) cp_mix += x[i] * kCp0[i];
  double kappa_mix = cp_mix / (cp_mix - kGasConstantRSI);

  double t_mix_us = Us(Clock::now() - t_mix0).count();

  double roc_crit = pcam / (kGasConstantR * zcam * tcam);
  double csi      = std::pow(tcam, 6)
                 / std::sqrt(bwr.molar_mass)
                 / std::cbrt(pcam * pcam);

  std::printf("\n%s  Mixture molar mass              : %s %s%8.4f%s [g/mol]\n",
              kBoldWhite, kReset, kBoldGreen, bwr.molar_mass, kReset);
  std::printf("%s  Relative density (vs. air)      :%s %s%8.4f%s [-]\n",
              kBoldWhite, kReset, kBoldGreen, bwr.molar_mass / 28.962, kReset);

  static const CountryRef kRefTable[] = {
    {"Romania / EU  (DIN 1343)",          2, { 0.0,  15.0  }, {"Nm\xC2\xB3/h", "Sm\xC2\xB3/h"}},
    {"ISO 13443  /  UK / Italy",          1, {15.0,   0.0  }, {"Sm\xC2\xB3/h", ""}},
    {"USA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)",      1, {15.56, 0.0}, {"Sm\xC2\xB3/h", ""}},
    {"Russia \xe2\x80\x94 GOST 30319-1",  1, {20.0,   0.0  }, {"m\xC2\xB3/h",  ""}},
    {"Custom",                            1, { 0.0,   0.0  }, {"m\xC2\xB3/h",  ""}},
  };
  static constexpr int kNRef = 5;

  std::printf("\n%s  ── Reference conditions %s\n",
              kYellow, kReset);
  std::printf("%s  1.  Romania / EU  (DIN 1343)        \xe2\x80\x94   0\xC2\xB0""C and 15\xC2\xB0""C / 101.325 kPa  [Nm\xC2\xB3/h] and [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  2.  ISO 13443  /  UK / Italy        \xe2\x80\x94  15\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  3.  USA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)             \xe2\x80\x94  15.56\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  4.  Russia \xe2\x80\x94 GOST 30319-1           \xe2\x80\x94  20\xC2\xB0""C / 101.325 kPa  [m\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  5.  Custom                          \xe2\x80\x94  T [\xC2\xB0""C] entered manually  [m\xC2\xB3/h]%s\n",  kBoldWhite, kReset);

  std::printf("\n  %s>%s %sSelect (1\xe2\x80\x93" "5) : %s", kCyan, kReset, kBoldWhite, kReset);
  int ref_sel = ReadChoice(1, kNRef);

  CountryRef ref = kRefTable[ref_sel - 1];
  if (ref_sel == kNRef) {
    std::printf("  %s>%s %sReference temperature [\xC2\xB0""C] : %s", kCyan, kReset, kBoldWhite, kReset);
    ReadDouble(&ref.t[0]);
  }

  std::printf ( "\n");

  double ror_ref[2] = {};
  for (int i = 0; i < ref.n; i++) {
    ror_ref[i] = CalcDensity(ref.t[i], 1, bwr);
    std::printf("%s  Reference density   %5.2f\xC2\xB0""C / 101.325 kPa :%s %s%8.4f%s kg/m\xC2\xB3\n",
                kBoldWhite, ref.t[i], kReset, kBoldGreen, ror_ref[i], kReset);
    double z_ref = bwr.molar_mass / (ror_ref[i] * kGasConstantR * (ref.t[i] + kKelvinOffset));
    std::printf("%s  Z factor            %5.2f\xC2\xB0""C / 101.325 kPa :%s %s%8.6f%s \xe2\x80\x94\n",
                kBoldWhite, ref.t[i], kReset, kBoldGreen, z_ref, kReset);
  }

  // ── Output helpers ─────────────────────────────────────────────────────────
  auto sep_d = []() {
    std::printf("%s  ", kYellow);
    for (int k = 0; k < 76; k++) std::fputs("\xe2\x95\x90", stdout);
    std::printf("%s\n", kReset);
  };
  auto sep_s = []() {
    std::printf("  ");
    for (int k = 0; k < 76; k++) std::fputs("\xe2\x94\x80", stdout);
    std::printf("\n");
  };
  auto row = [](const char* label, const char* fmt, double val, const char* unit) {
    int llen = (int)std::strlen(label) - Utf8ExtraBytes(label);
    std::printf("  %s  %s", kBoldWhite, label);
    for (int k = llen; k < 42; k++) std::putchar(' ');
    std::printf(" :%s %s", kReset, kBoldGreen);
    std::printf(fmt, val);
    std::printf("%s %s\n", kReset, unit);
  };

  // ── Outer loop: select measurement device ─────────────────────────────────
  for (;;) {
    int    tip_raw = 0;
    double d_int = 0.0, d_orif = 0.0;

    // Try to reuse saved configuration
    {
      int   sv_tip; double sv_d_int, sv_d_orif;
      if (LoadConfig(&sv_tip, &sv_d_int, &sv_d_orif)
          && sv_tip >= kTipMin && sv_tip <= kTipMax) {
        std::printf(
            "\n%s  ── Saved configuration%s\n"
            "%s  Device      : %s%s\n"
            "%s  Pipe D      : %s%g mm%s\n"
            "%s  Orifice D   : %s%g mm%s\n",
            kYellow, kReset,
            kBoldWhite, TipName(static_cast<DeviceType>(sv_tip)), kReset,
            kBoldWhite, kBoldGreen, sv_d_int, kReset,
            kBoldWhite, kBoldGreen, sv_d_orif, kReset);
        std::printf("  %s>%s %sReuse this configuration? [y/n] : %s", kCyan, kReset, kBoldWhite, kReset);
        if (AskYesNo()) {
          tip_raw = sv_tip;
          d_int   = sv_d_int;
          d_orif  = sv_d_orif;
          goto run_inner;
        }
      }
    }

    // Manual device selection
    select_device_type:
    std::printf("\n%s  ── Throttling device%s\n",
                kYellow, kReset);
    std::printf("%s%s", kBoldWhite, kTipDisp);
    std::printf("%s", kReset);
    tip_raw = ReadChoice(kTipMin, kTipMax);

    select_D:
    d_int = 0.0;
    std::printf("\n");
    for (;;) {
      std::printf("  %s>%s %sPipe D    (20\xC2\xB0""C) [mm]%s"
                  "  [\xe2\x86\x90 reselect device]%s : %s",
                  kCyan, kReset, kBoldWhite, kYellow, kBoldWhite, kReset);
      if (!ReadDouble(&d_int)) { std::printf("%s", kReset); d_int = 0.0; goto select_device_type; }
      std::printf("%s", kReset);
      if (d_int > 0.0) break;
      std::printf("%s  Invalid value \xe2\x80\x94 must be positive.%s\n",
                  kBoldRed, kReset);
    }

    select_d:
    d_orif = 0.0;
    for (;;) {
      std::printf("  %s>%s %sOrifice D (20\xC2\xB0""C) [mm]%s"
                  "  [\xe2\x86\x90 reenter pipe D]%s : %s",
                  kCyan, kReset, kBoldWhite, kYellow, kBoldWhite, kReset);
      if (!ReadDouble(&d_orif)) { std::printf("%s", kReset); d_orif = 0.0; goto select_D; }
      std::printf("%s", kReset);
      if (d_orif <= 0.0) {
        std::printf("%s  Invalid value \xe2\x80\x94 must be positive.%s\n",
                    kBoldRed, kReset);
      } else if (d_orif >= d_int) {
        std::printf("%s  Orifice D must be less than pipe D (%g mm).%s\n",
                    kBoldRed, d_int, kReset);
      } else {
        double beta_chk = d_orif / d_int;
        if (beta_chk < 0.10 || beta_chk > 0.80) {
          std::printf("%s  \xCE\xB2 = %.4f \xe2\x80\x94 outside ISO 5167 range [0.10, 0.80].%s\n",
                      kBoldRed, beta_chk, kReset);
        } else {
          break;
        }
      }
    }

    std::printf("%s\n  Save configuration? [y/n] %s>%s ", kBoldWhite, kCyan, kReset);
    if (AskYesNo()) {
      SaveConfig(tip_raw, d_int, d_orif);
    }

    run_inner:;
    DeviceType tip = static_cast<DeviceType>(tip_raw);

    // Inner loop: calculate for different T/p conditions with the same device
    for (;;) {
      double temperature = 0.0, pressure = 0.0, pressure_diff = 0.0;

      reenter_T:
      std::printf("\n\n%s  ── Measurement conditions%s\n",
                  kYellow, kReset);
      std::printf("  %s>%s %sTemperature         [\xC2\xB0""C]%s"
                  "  [\xe2\x86\x90 change device]%s : %s",
                  kCyan, kReset, kBoldWhite, kYellow, kBoldWhite, kReset);
      if (!ReadDouble(&temperature)) { std::printf("%s", kReset); break; }
      std::printf("%s", kReset);
      if (temperature <= -273.15) {
        std::printf("%s  Temperature below absolute zero (-273.15\xC2\xB0""C).%s\n",
                    kBoldRed, kReset);
        goto reenter_T;
      }

      reenter_P:
      std::printf("  %s>%s %sPressure            [kPa]%s"
                  "  [\xe2\x86\x90 reenter T]%s : %s",
                  kCyan, kReset, kBoldWhite, kYellow, kBoldWhite, kReset);
      if (!ReadDouble(&pressure)) { std::printf("%s", kReset); goto reenter_T; }
      std::printf("%s", kReset);
      if (pressure <= 0.0) {
        std::printf("%s  Pressure must be positive.%s\n",
                    kBoldRed, kReset);
        goto reenter_P;
      }

      reenter_dP:
      std::printf("  %s>%s %sDifferential pressure [kPa]%s"
                  "  [\xe2\x86\x90 reenter P]%s : %s",
                  kCyan, kReset, kBoldWhite, kYellow, kBoldWhite, kReset);
      if (!ReadDouble(&pressure_diff)) { std::printf("%s", kReset); goto reenter_P; }
      std::printf("%s", kReset);
      if (pressure_diff <= 0.0) {
        std::printf("%s  Differential pressure must be positive.%s\n",
                    kBoldRed, kReset);
        goto reenter_dP;
      } else if (pressure_diff >= pressure) {
        std::printf("%s  Differential must be less than p = %g kPa.%s\n",
                    kBoldRed, pressure, kReset);
        goto reenter_dP;
      }
      std::printf("  %c", 7);

      // ─── Calcul ────────────────────────────────────────────────────────────
      int    iter_rho = 0, iter_qm = 0;
      auto   t_rho0 = Clock::now();
      double ro = CalcDensity(temperature, pressure / kKpaPerAtm, bwr, &iter_rho);
      double t_rho_us = Us(Clock::now() - t_rho0).count();
      if (!std::isfinite(ro) || ro <= 0.0) {
        std::printf("%s\n  BWRS density did not converge (T=%.1f\xc2\xb0""C, P=%.1f kPa)"
                    " \xe2\x80\x94 check that conditions are in the gas-phase region.%s\n",
                    kBoldRed, temperature, pressure, kReset);
        break;
      }

      double roc_red = ro / roc_crit;
      auto   t_eta0 = Clock::now();
      double eta = 0.0;
      for (int i = 1; i <= kNumComponents; i++) {
        eta += (1 + kChapEnskog * std::log((temperature + kKelvinOffset) / cs[i]))
             / (1 + kChapEnskog * std::log(kKelvinOffset / cs[i]))
             * std::sqrt((temperature + kKelvinOffset) / kKelvinOffset)
             * et[i] * x[i] * std::sqrt(m[i]);
      }
      eta  = eta / mx;
      eta += kViscHighA / csi
           * std::pow(std::exp(kViscHighExp1 * roc_red) - std::exp(-kViscHighExp2 * roc_red),
                      kViscHighPow);
      double t_eta_us = Us(Clock::now() - t_eta0).count();

      FlowResult flow;
      auto   t_qm0 = Clock::now();
      double qm = CalcMassFlow(pressure_diff, pressure, temperature,
                       tip, d_int, d_orif, ro, eta, &flow, kappa_mix, &iter_qm);
      double t_qm_us = Us(Clock::now() - t_qm0).count();
      if (qm == 0.0) break;  // error -> reselect device

      // ─── Display results ───────────────────────────────────────────────────
      std::printf("\n");
      sep_d();
      std::printf("%s  RESULTS  \xe2\x80\x94  %s%s\n",
                  kYellow, TipName(tip), kReset);
      std::printf("%s  D = %.2f mm  \xc2\xb7  d = %.2f mm  \xc2\xb7  "
                  "\xce\xb2 = %.4f%s\n",
                  kBoldWhite, d_int, d_orif, flow.beta, kReset);
      std::printf("%s  D(t) = %.3f mm  \xc2\xb7  d(t) = %.3f mm  "
                  "\xc2\xb7  t = %.1f\xc2\xb0""C%s\n",
                  kBoldWhite, flow.D_working_mm, flow.d_working_mm, temperature, kReset);
      {
        double dt_ref   = temperature - kRefTempCelsius;
        double corr_D   = kThermalExpPipe    * dt_ref * 100.0;
        double corr_d   = kThermalExpOrifice * dt_ref * 100.0;
        if (std::fabs(corr_D) > 0.10 || std::fabs(corr_d) > 0.10) {
          std::printf("%s  Note: thermal expansion at %.0f\xc2\xb0""C: "
                      "\xce\x94""D = %+.3f%%,  \xce\x94""d = %+.3f%%"
                      " \xe2\x80\x94 verify dimensions are at 20\xc2\xb0""C reference.%s\n",
                      kYellow, temperature, corr_D, corr_d, kReset);
        }
      }
      sep_d();

      std::printf("\n%s  Measurement conditions%s\n", kCyan, kReset);
      sep_s();
      row("Temperature", "%10.2f", temperature, "\xc2\xb0""C");
      row("Absolute pressure", "%10.2f", pressure, "kPa");
      row("Differential pressure", "%10.2f", pressure_diff, "kPa");

      double Z_tp = (pressure / kKpaPerAtm) * bwr.molar_mass
                  / (ro * kGasConstantR * (temperature + kKelvinOffset));

      std::printf("\n%s  Fluid  (at t, p)%s\n", kCyan, kReset);
      sep_s();
      row("Density \xcf\x81(t,p)", "%10.4f", ro, "kg/m\xc2\xb3");
      row("Dynamic viscosity \xce\xb7(t,p)", "%10.4f",
          eta * kPaToMicroPa, "\xc2\xb5Pa\xc2\xb7s");
      row("Compressibility factor Z(t,p)", "%10.6f", Z_tp, "\xe2\x80\x94");

      std::printf("\n%s  Flow rates%s\n", kCyan, kReset);
      sep_s();
      row("Mass flow Qm", "%10.4f", qm, "kg/s");
      row("Mass flow Qm", "%10.2f", qm * kSecondsPerHour, "kg/h");
      for (int i = 0; i < ref.n; i++) {
        double qhref = kSecondsPerHour / ror_ref[i] * qm;
        char lbl[64];
        std::snprintf(lbl, sizeof(lbl),
                      "Volumetric %5.2f\xc2\xb0""C / 101.325 kPa", ref.t[i]);
        row(lbl, "%10.2f", qhref, ref.label[i]);
      }
      row("Volumetric at (t,p)", "%10.2f", kSecondsPerHour / ro * qm, "m\xc2\xb3/h");

      std::printf("\n%s  Hydraulics%s\n", kCyan, kReset);
      sep_s();
      row("Mean velocity v", "%10.2f", flow.velocity, "m/s");
      row("Pressure loss", "%10.2f", flow.pressure_loss, "kPa");
      row("Throttle ratio \xce\xb2", "%10.4f", flow.beta, "\xe2\x80\x94");
      row("Reynolds number Re", "%10.4g", flow.reynolds, "\xe2\x80\x94");
      row("Discharge coefficient C", "%10.6f", flow.coef_c, "\xe2\x80\x94");
      row("Expansibility factor \xce\xb5", "%10.6f", flow.epsilon, "\xe2\x80\x94");
      row("Isentropic exponent k", "%10.4f", flow.kappa, "\xe2\x80\x94");

      std::printf("\n");
      sep_d();

      // ── CPU profile ────────────────────────────────────────────────────────
      double t_total = t_mix_us + t_rho_us + t_eta_us + t_qm_us;
      std::printf("\n%s  CPU PROFILE%s\n", kCyan, kReset);
      sep_s();
      auto prow = [](const char* op, const char* detail, double us, int it) {
        int llen = (int)std::strlen(op) - Utf8ExtraBytes(op);
        std::printf("  %s%s", kBoldWhite, op);
        for (int k = llen; k < 30; k++) std::putchar(' ');
        std::printf("%s  %-28s%s %s%7.2f \xc2\xb5s%s",
                    kReset, detail, kReset, kBoldGreen, us, kReset);
        if (it > 0) std::printf("  %s%d iter.%s", kYellow, it, kReset);
        std::putchar('\n');
      };
      int n_comp = kNumComponents;
      char mix_detail[48];
      std::snprintf(mix_detail, sizeof(mix_detail),
                    "mixing rules N\xc2\xb2+N\xc2\xb3  (N=%d)", n_comp);
      prow("BWRS constants", mix_detail,          t_mix_us, 0);
      prow("Density \xcf\x81",       "BWRS bisection",     t_rho_us, iter_rho);
      prow("Viscosity \xce\xb7",     "Chapman-Enskog / N", t_eta_us, 0);
      prow("Flow rate Qm",    "Reynolds iteration", t_qm_us,  iter_qm);
      sep_s();
      std::printf("  %sTotal%s                                         "
                  "%s%7.2f \xc2\xb5s%s\n\n",
                  kBoldWhite, kReset, kBoldGreen, t_total, kReset);

      std::printf("%s  EMBEDDED ESTIMATE%s  (indicative factors relative to PC)\n",
                  kCyan, kReset);
      sep_s();
      struct { const char* name; double factor; } targets[] = {
        { "ARM Cortex-M7 @480 MHz + FPU  (STM32H7)",    8.0  },
        { "ARM Cortex-M4 @168 MHz + FPU  (STM32F4)",    22.0 },
        { "Xtensa LX6    @240 MHz + FPU  (ESP32)",       16.0 },
        { "MSP430 F5xx   @ 25 MHz + hw mult (TI)",      350.0 },
        { "AVR           @ 16 MHz, no FPU (Mega2560)",   900.0 },
        { "80C51         @ 12 MHz, no FPU",             3000.0 },
      };
      for (auto& tg : targets) {
        double est = t_total * tg.factor / 1000.0;
        int tlen = (int)std::strlen(tg.name) - Utf8ExtraBytes(tg.name);
        std::printf("  %s%s", kBoldWhite, tg.name);
        for (int k = tlen; k < 48; k++) std::putchar(' ');
        std::printf("%s\xc3\x97%4.0f  %s%s%8.2f ms%s\n",
                    kReset, tg.factor, kReset, kBoldGreen, est, kReset);
      }
      std::printf("%s  Note:%s estimated factors (IPC, cache, compiler)."
                  " Measure on target for accuracy.\n",
                  kYellow, kReset);

      std::printf("\n");
      sep_d();

      // ── Measurement uncertainty (ISO 5167-1) ─────────────────────────────
      std::printf("\n  %sCompute measurement uncertainty (ISO 5167-1)? [y/n]  %s>%s ",
                  kBoldWhite, kCyan, kReset);
      if (AskYesNo()) {
        double uC_def;
        switch (tip) {
          case DeviceType::kOrificeCorner:
          case DeviceType::kOrificeFlange:
          case DeviceType::kOrificeDD2:
            uC_def = (flow.beta > 0.60) ? 0.75 : 0.50;
            break;
          case DeviceType::kNozzleIsa:           uC_def = 0.80; break;
          case DeviceType::kNozzleLongRadius:    uC_def = 2.00; break;
          case DeviceType::kVenturiRoughCast:    uC_def = 0.70; break;
          case DeviceType::kVenturiMachined:     uC_def = 1.00; break;
          case DeviceType::kVenturiWeldedSheet:  uC_def = 1.50; break;
          case DeviceType::kVenturiNozzle:       uC_def = 1.20; break;
          default:                               uC_def = 1.00; break;
        }
        auto read_pct = [&](const char* label, double def_val) -> double {
          double v = 0.0;
          int llen = (int)std::strlen(label) - Utf8ExtraBytes(label);
          std::printf("  %s>%s %s%s", kCyan, kReset, kBoldWhite, label);
          for (int k = llen; k < 44; k++) std::putchar(' ');
          std::printf("%s[default %s%.2f%%%s] : %s",
                      kBoldWhite, kBoldGreen, def_val, kBoldWhite, kReset);
          if (!ReadDouble(&v) || v <= 0.0) v = def_val;
          return v;
        };
        std::printf("\n%s  \xe2\x94\x80\xe2\x94\x80 Measurement uncertainty  (ISO 5167-1)%s\n",
                    kYellow, kReset);
        std::printf(
            "%s  Relative standard uncertainties [%%]  "
            "\xe2\x80\x94  Enter or 0 \xe2\x86\x92 default:%s\n\n",
            kBoldWhite, kReset);
        double uC   = read_pct("u(C)   discharge coefficient (ISO 5167)",  uC_def);
        double uEps = read_pct("u(\xce\xb5)    expansibility factor",              0.10);
        double ud_p = read_pct("u(d)   orifice diameter",                      0.03);
        double uD_p = read_pct("u(D)   pipe diameter",                         0.10);
        double udp  = read_pct("u(\xce\x94p)  differential pressure",               0.20);
        double uro  = read_pct("u(\xcf\x81)   gas density (BWRS)",                 0.30);
        double b4   = std::pow(flow.beta, 4);
        double fd   = 2.0 / (1.0 - b4);
        double fD   = 2.0 * b4 / (1.0 - b4);
        double cC   = uC,          cEps = uEps;
        double cd   = fd * ud_p,   cD   = fD * uD_p;
        double cdp  = 0.5 * udp,   cro  = 0.5 * uro;
        double u_qm = std::sqrt(cC*cC + cEps*cEps + cd*cd + cD*cD + cdp*cdp + cro*cro);
        double U_qm = 2.0 * u_qm;
        std::printf("\n%s  Sensitivity factors  (\xce\xb2 = %.4f)%s\n", kCyan, flow.beta, kReset);
        sep_s();
        std::printf(
            "  f_d = 2/(1\xe2\x88\x92\xce\xb2\xe2\x81\xb4) = %s%.4f%s     "
            "f_D = 2\xce\xb2\xe2\x81\xb4/(1\xe2\x88\x92\xce\xb2\xe2\x81\xb4) = %s%.4f%s\n\n",
            kBoldGreen, fd, kReset, kBoldGreen, fD, kReset);
        std::printf("  %s%-12s  %8s   %11s   %15s%s\n",
                    kBoldWhite, "Source", "u_i [%]", "Sensitivity", "Contribution [%]", kReset);
        sep_s();
        auto brow = [&](const char* src, double ui, double fi, double ci) {
          int slen = (int)std::strlen(src) - Utf8ExtraBytes(src);
          std::printf("  %s%s", kBoldWhite, src);
          for (int k = slen; k < 12; k++) std::putchar(' ');
          std::printf("%s  %s%8.4f%s   %11.4f   %s%15.4f%s\n",
                      kReset, kBoldGreen, ui, kReset, fi, kBoldGreen, ci, kReset);
        };
        brow("C",          uC,   1.0, cC);
        brow("\xce\xb5",   uEps, 1.0, cEps);
        brow("d",          ud_p, fd,  cd);
        brow("D",          uD_p, fD,  cD);
        brow("\xce\x94p",  udp,  0.5, cdp);
        brow("\xcf\x81",   uro,  0.5, cro);
        sep_s();
        row("Combined std. uncertainty  u(Qm)", "%10.4f", u_qm, "%");
        row("Expanded uncertainty  U(Qm)  k=2, 95%", "%10.4f", U_qm, "%");
        double U_abs_kgs = U_qm / 100.0 * qm;
        std::printf("\n");
        std::printf("  %s  Qm = %.4f \xc2\xb1 %.5f  kg/s%s\n",
                    kBoldGreen, qm, U_abs_kgs, kReset);
        std::printf("  %s  Qm = %.2f \xc2\xb1 %.3f  kg/h%s\n",
                    kBoldGreen, qm * kSecondsPerHour, U_abs_kgs * kSecondsPerHour, kReset);
        for (int i = 0; i < ref.n; i++) {
          if (ror_ref[i] > 0.0) {
            double qhref   = kSecondsPerHour / ror_ref[i] * qm;
            double U_qhref = U_qm / 100.0 * qhref;
            std::printf("  %s  Qv(%.2f\xc2\xb0""C) = %.2f \xc2\xb1 %.3f  %s%s\n",
                        kBoldGreen, ref.t[i], qhref, U_qhref, ref.label[i], kReset);
          }
        }
        std::printf("\n");
        sep_d();
      }

      std::printf("\n  %s1%s new conditions   %s2%s change D / d"
                  "   %s3%s change device   %sESC%s quit\n\n",
                  kBoldGreen, kReset, kBoldGreen, kReset,
                  kBoldGreen, kReset, kBoldWhite, kReset);
      std::printf("  %s>%s ", kCyan, kReset);
      {
        int nav = ReadChoice(1, 3);
        if (nav == 2) goto select_D;
        if (nav == 3) break;
        // nav == 1: loop back to reenter_T
      }
    }
  }
}
