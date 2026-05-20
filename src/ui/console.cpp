#include "console.h"
#include "../tests/validation.h"
#include "../formulas/bwrs.h"

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <conio.h>
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

const char* const COMP_NAMES[ARRAY_SIZE] = {
    nullptr,
    "Methane", "Ethane", "Propane", "Isobutane", "N-butane",
    "Neopentane", "Isopentane", "N-pentane", "2,2-dimethylbutane", "2,3-dimethylbutane",
    "3-methylpentane", "2-methylpentane", "N-hexane", "2,4-dimethylpentane", "2,2,3-trimethylbutane",
    "2-methylhexane", "3-methylhexane", "3-ethylpentane", "N-heptane", "2,2,4-trimethylpentane",
    "N-octane", "Benzene", "Toluene", "Hydrogen", "Carbon monoxide",
    "Hydrogen sulfide", "Helium", "Argon", "Nitrogen", "Oxygen",
    "Carbon dioxide", "Ethylene", "Propylene", "Ammonia", "Acetylene"
};

void ExitApp() {
  std::printf("\033[2J\033[H%s\n  ELCOST Impex  —  BWR Gas Flow Calculator  v3.0/2026%s\n",
              COLOR_YELLOW, COLOR_RESET);
  std::exit(0);
}

bool ReadDouble(double* val) {
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
        std::printf(U_LARR "\n");   // ←
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

int ReadChoice(int lo, int hi) {
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

bool AskYesNo() {
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
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

int Utf8ExtraBytes(const char* s) {
  int n = 0;
  while (*s) {
    if (((unsigned char)*s & 0xC0) == 0x80) n++;
    s++;
  }
  return n;
}

void PrintComposition(const double* x) {
  std::printf("\n%s  " U_HLINE U_HLINE " Composition%s\n",
              COLOR_YELLOW, COLOR_RESET);
  std::printf("%s  Molar fractions of the gas mixture:%s\n\n", COLOR_BOLD_WHITE, COLOR_RESET);
  for (int i = 1; i <= NUM_COMPONENTS; i++) {
    const char* name = COMP_NAMES[i];
    int visualLen = (int)std::strlen(name) - Utf8ExtraBytes(name);
    std::printf("  %s>%s %s%s", COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, name);
    for (int j = visualLen; j < COMP_NAME_WIDTH; j++) std::putchar(' ');
    std::printf(": %s%.8f%s\n", COLOR_BOLD_GREEN, x[i], COLOR_RESET);
  }
  std::printf("\n");
}

bool LoadComposition(double* x) {
  FILE* f = nullptr;
  fopen_s(&f, COMP_FILE, "r");
  if (!f) return false;
  for (int i = 1; i <= NUM_COMPONENTS; i++) {
    if (fscanf_s(f, "%lf", &x[i]) != 1) {
      std::fclose(f);
      return false;
    }
  }
  std::fclose(f);
  return true;
}

void SaveComposition(const double* x) {
  FILE* f = nullptr;
  fopen_s(&f, COMP_FILE, "w");
  if (!f) {
    std::printf("%s  Could not open %s for writing.%s\n", COLOR_BOLD_RED, COMP_FILE, COLOR_RESET);
    return;
  }
  bool ok = true;
  for (int i = 1; i <= NUM_COMPONENTS; i++)
    ok &= (std::fprintf(f, "%.8f\n", x[i]) > 0);
  if (std::fclose(f) != 0 || !ok)
    std::printf("%s  Write error: %s may be incomplete (disk full?).%s\n",
                COLOR_BOLD_RED, COMP_FILE, COLOR_RESET);
}

bool LoadConfig(int* tip_raw, double* d_int, double* d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, CONF_FILE, "r");
  if (!f) return false;
  bool ok = (fscanf_s(f, "%d %lf %lf", tip_raw, d_int, d_orif) == 3);
  std::fclose(f);
  return ok;
}

void SaveConfig(int tip_raw, double d_int, double d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, CONF_FILE, "w");
  if (!f) {
    std::printf("%s  Could not open %s for writing.%s\n", COLOR_BOLD_RED, CONF_FILE, COLOR_RESET);
    return;
  }
  bool ok = (std::fprintf(f, "%d\n%.4f\n%.4f\n", tip_raw, d_int, d_orif) > 0);
  if (std::fclose(f) != 0 || !ok)
    std::printf("%s  Write error: %s may be incomplete (disk full?).%s\n",
                COLOR_BOLD_RED, CONF_FILE, COLOR_RESET);
}

static void hdr_thin(int n = 80) {
  std::printf("  ");
  for (int k = 0; k < n; k++) std::fputs(U_HLINE, stdout);
  std::putchar('\n');
}

static void hdr_dbl(int n = 80) {
  std::printf("  ");
  for (int k = 0; k < n; k++) std::fputs(U_DLINE, stdout);
  std::putchar('\n');
}

void PrintBanner() {
  std::printf("\033[2J\033[H\n");

  // ── Main header ────────────────────────────────────────────────────────────
  std::printf("%s", COLOR_HDR_YELLOW);
  hdr_dbl(81);
  std::printf(
      "  %sELCOST Impex%s  " U_CDOT "  BWRS Gas Flow Calculator v3.0/2026\n"
      "  Author  %sConstantin Agavriloaie%s  " U_CDOT "  CJ-RO\n"
      "  Email   office@elcost.ro\n",
      COLOR_HDR_GREEN, COLOR_HDR_YELLOW, COLOR_HDR_YELLOW, COLOR_HDR_YELLOW);
  hdr_dbl(81);
  std::printf("\n");

  // ── Calculation models ─────────────────────────────────────────────────────
  std::printf("%s  CALCULATION MODELS%s\n", COLOR_HDR_CYAN, COLOR_RESET);
  std::printf(
      "  %sEquation of state%s   BWRS " U_CDOT " Starling 1973  " U_MDASH "  Benedict-Webb-Rubin-Starling, 11 param.\n"
      "  %sViscosity%s           Chapman-Enskog with high-density correction  [Nishiumi 1975]\n"
      "  %sC coefficient, " U_EPS "%s    Reader-Harris/Gallagher  " U_CDOT "  ISO 5167-2/3/4:2003\n\n",
      COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET);

  // ── Scope ──────────────────────────────────────────────────────────────────
  std::printf("%s  SCOPE OF APPLICATION%s\n", COLOR_HDR_CYAN, COLOR_RESET);
  std::printf(
      "  35 components  " U_CDOT "  9 throttling device types\n"
      "  Mass flow, volumetric flow, velocity, pressure loss  " U_CDOT "  selectable reference conditions\n"
      "  ISO 5167 limit validation  (D, " U_BETA ", Re)\n\n");

  // ── Standards and references ───────────────────────────────────────────────
  std::printf("%s", COLOR_HDR_YELLOW);
  hdr_thin();
  std::printf("  %sSTANDARDS AND REFERENCES%s\n", COLOR_HDR_CYAN, COLOR_RESET);
  std::printf(
      "  %sISO 5167-2:2003%s   Orifice plates " U_MDASH " Reader-Harris/Gallagher equation\n"
      "  %sISO 5167-3:2003%s   Nozzles and Venturi nozzles\n"
      "  %sISO 5167-4:2003%s   Classical Venturi tubes\n"
      "  %sStarling K.E.%s     Fluid Thermodynamic Properties, Gulf Publ. Houston (1973)\n"
      "  %sNishiumi & Saito%s  J. Chem. Eng. Japan 8(5), 356" U_NDASH "360 (1975)\n\n",
      COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET,
      COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET);

  // ── Validation ─────────────────────────────────────────────────────────────
  std::printf("%s", COLOR_HDR_YELLOW);
  hdr_thin();
  std::printf("  %sVALIDATION%s\n", COLOR_HDR_CYAN, COLOR_RESET);
  std::printf(
      "  " U_RHO "  validated for 8 compositions (CH" U_SUB4 ", C" U_SUB2 "H" U_SUB6 ", CO" U_SUB2 ","
      " H" U_SUB2 ", N" U_SUB2 ", Ar, std NG, rich NG)  vs. NIST WebBook\n"
      "  " U_ETA "  validated vs. Chapman-Enskog / NIST\n"
      "  Qm validated " U_PLUSMN "15 %% against ISO 5167 estimates  " U_CDOT "  9 device types\n\n");

  // ── Footer ─────────────────────────────────────────────────────────────────
  std::printf("%s", COLOR_HDR_YELLOW);
  hdr_thin();
  std::printf(U_COPY " 2004" U_NDASH "2026 ELCOST Impex\n");
  std::printf("%s", COLOR_HDR_YELLOW);
  hdr_dbl();
  std::printf("\n");

  std::printf("%s", COLOR_RESET);

  // ── Optional self-test at startup ─────────────────────────────────────────
  std::printf("%s\n  Run implementation validation test? [y/n] %s>%s ",
              COLOR_BOLD_WHITE, COLOR_CYAN, COLOR_RESET);
  if (AskYesNo()) {
    runValidationTest();
  }
}

// ── Static output helpers (used only within this file) ────────────────────

static void sep_d() {
  std::printf("%s  ", COLOR_YELLOW);
  for (int k = 0; k < 76; k++) std::fputs(U_DLINE, stdout);
  std::printf("%s\n", COLOR_RESET);
}

static void sep_s() {
  std::printf("  ");
  for (int k = 0; k < 76; k++) std::fputs(U_HLINE, stdout);
  std::printf("\n");
}

static void row(const char* label, const char* fmt, double val, const char* unit) {
  int llen = (int)std::strlen(label) - Utf8ExtraBytes(label);
  std::printf("  %s  %s", COLOR_BOLD_WHITE, label);
  for (int k = llen; k < 42; k++) std::putchar(' ');
  std::printf(" :%s %s", COLOR_RESET, COLOR_BOLD_GREEN);
  std::printf(fmt, val);
  if (unit && unit[0])
    std::printf("%s [%s]\n", COLOR_RESET, unit);
  else
    std::printf("%s\n", COLOR_RESET);
}

// ─────────────────────────────────────────────────────────────────────────────

void InitConsole() {
#ifdef _WIN32
  std::system("chcp 65001 > nul");
  {
    HANDLE hCon = GetStdHandle(STD_OUTPUT_HANDLE);

    CONSOLE_FONT_INFOEX cfi = {};
    cfi.cbSize = sizeof(cfi);
    GetCurrentConsoleFontEx(hCon, FALSE, &cfi);
    cfi.dwFontSize.Y = 14;
    wcscpy_s(cfi.FaceName, L"Consolas");
    SetCurrentConsoleFontEx(hCon, FALSE, &cfi);

    DWORD mode = 0;
    GetConsoleMode(hCon, &mode);
    SetConsoleMode(hCon, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
  }
#endif
}

void ReadComposition(double* x) {
  bool comp_loaded = false;
  {
    FILE* cf = nullptr;
    fopen_s(&cf, COMP_FILE, "r");
    if (cf) {
      std::fclose(cf);
      std::printf("%s\n  A saved composition exists. Reuse it? [y/n] %s>%s ",
                  COLOR_BOLD_WHITE, COLOR_CYAN, COLOR_RESET);
      if (AskYesNo()) {
        if (LoadComposition(x)) {
          double loaded_sum = 0.0;
          for (int i = 1; i <= NUM_COMPONENTS; i++) loaded_sum += x[i];
          if (std::fabs(loaded_sum - 1.0) > SUM_TOLERANCE) {
            std::printf("%s  Loaded composition: sum = %.6f " U_NEQ " 1 "
                        "— switching to manual entry.%s\n", COLOR_BOLD_RED, loaded_sum, COLOR_RESET);
            for (int i = 1; i <= NUM_COMPONENTS; i++) x[i] = 0.0;
          } else {
            comp_loaded = true;
            PrintComposition(x);
          }
        } else {
          std::printf("%s  Error reading file. Switching to manual entry.%s\n",
                      COLOR_BOLD_RED, COLOR_RESET);
        }
      }
    }
  }

  if (!comp_loaded) {
    double sum = 0.0;
    do {
      std::printf("\n%s  " U_HLINE U_HLINE " Composition %s\n",
                  COLOR_YELLOW, COLOR_RESET);
      std::printf("%s  Molar fractions of the gas mixture:%s\n\n", COLOR_BOLD_WHITE, COLOR_RESET);
      for (int i = 1; i <= NUM_COMPONENTS; i++) {
        const char* name = COMP_NAMES[i];
        int vlen = (int)std::strlen(name) - Utf8ExtraBytes(name);
        double v;
        do {
          std::printf("  %s>%s %s%s", COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, name);
          for (int j = vlen; j < COMP_NAME_WIDTH; j++) std::putchar(' ');
          std::printf(": ");
          ReadDouble(&v);
          std::printf("%s", COLOR_RESET);
          if (v < 0.0 || v > 1.0) {
            std::printf("%s  Invalid value " U_MDASH " must be in [0, 1].%s\n", COLOR_BOLD_RED, COLOR_RESET);
          }
        } while (v < 0.0 || v > 1.0);
        x[i] = v;
      }

      sum = 0.0;
      for (int i = 1; i <= NUM_COMPONENTS; i++) sum += x[i];

      if (std::fabs(sum - 1.0) > SUM_TOLERANCE) {
        std::printf("%s\n  Sum of molar fractions = %.6f  " U_NEQ "  1.%s\n", COLOR_BOLD_RED, sum, COLOR_RESET);
        if (sum < 1e-10) {
          std::printf("%s  Sum is zero " U_MDASH " re-enter the composition.%s\n", COLOR_BOLD_RED, COLOR_RESET);
        } else {
          std::printf("%s  Normalise automatically? [y/n] (n = re-enter) %s>%s ", COLOR_YELLOW, COLOR_CYAN, COLOR_RESET);
          if (AskYesNo()) {
            for (int i = 1; i <= NUM_COMPONENTS; i++) x[i] /= sum;
            sum = 1.0;
            std::printf("%s  Normalised fractions (non-zero components):%s\n\n", COLOR_BOLD_GREEN, COLOR_RESET);
            for (int i = 1; i <= NUM_COMPONENTS; i++) {
              if (x[i] > 0.0) {
                const char* name = COMP_NAMES[i];
                int vlen = (int)std::strlen(name) - Utf8ExtraBytes(name);
                std::printf("    %s%s", COLOR_BOLD_GREEN, name);
                for (int j = vlen; j < COMP_NAME_WIDTH; j++) std::putchar(' ');
                std::printf(": %.8f%s\n", x[i], COLOR_RESET);
              }
            }
          }
        }
      }
    } while (std::fabs(sum - 1.0) > SUM_TOLERANCE);

    std::printf("%s\n  Save composition? [y/n] %s>%s ", COLOR_BOLD_WHITE, COLOR_CYAN, COLOR_RESET);
    if (AskYesNo()) SaveComposition(x);
  }
}

void SelectRefConditions(const BwrConst& bwr, CountryRef* out_ref, double ror_ref[2]) {
  std::printf("\n%s  Mixture molar mass              : %s %s%8.4f%s [g/mol]\n",
              COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_GREEN, bwr.m_molarMass, COLOR_RESET);
  std::printf("%s  Relative density (vs. air)      :%s %s%8.4f%s [-]\n",
              COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_GREEN, bwr.m_molarMass / 28.962, COLOR_RESET);

  static const CountryRef kRefTable[] = {
    {"Romania / EU  (DIN 1343)",         2, { 0.0,  15.0  }, {"Nm" U_SUP3 "/h", "Sm" U_SUP3 "/h"}},
    {"ISO 13443  /  UK / Italy",         1, {15.0,   0.0  }, {"Sm" U_SUP3 "/h", ""}},
    {"USA " U_MDASH " AGA-3  (60" U_DEG"F)",     1, {15.56, 0.0}, {"Sm" U_SUP3 "/h", ""}},
    {"Russia " U_MDASH " GOST 30319-1", 1, {20.0,   0.0  }, {"m" U_SUP3 "/h",  ""}},
    {"Custom",                           1, { 0.0,   0.0  }, {"m" U_SUP3 "/h",  ""}},
  };
  static constexpr int kNRef = 5;

  std::printf("\n%s  " U_HLINE U_HLINE " Reference conditions %s\n", COLOR_YELLOW, COLOR_RESET);
  std::printf("%s  1.  Romania / EU  (DIN 1343)        " U_MDASH "   0" U_DEG"C and 15" U_DEG"C / 101.325 kPa  [Nm" U_SUP3 "/h] and [Sm" U_SUP3 "/h]%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("%s  2.  ISO 13443  /  UK / Italy        " U_MDASH "  15" U_DEG"C / 101.325 kPa  [Sm" U_SUP3 "/h]%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("%s  3.  USA " U_MDASH " AGA-3  (60" U_DEG"F)             " U_MDASH "  15.56" U_DEG"C / 101.325 kPa  [Sm" U_SUP3 "/h]%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("%s  4.  Russia " U_MDASH " GOST 30319-1           " U_MDASH "  20" U_DEG"C / 101.325 kPa  [m" U_SUP3 "/h]%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("%s  5.  Custom                          " U_MDASH "  T [" U_DEG"C] entered manually  [m" U_SUP3 "/h]%s\n", COLOR_BOLD_WHITE, COLOR_RESET);
  std::printf("\n  %s>%s %sSelect (1" U_NDASH"5) : %s", COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET);
  int ref_sel = ReadChoice(1, kNRef);

  *out_ref = kRefTable[ref_sel - 1];
  if (ref_sel == kNRef) {
    std::printf("  %s>%s %sReference temperature [" U_DEG"C] : %s", COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, COLOR_RESET);
    ReadDouble(&out_ref->m_t[0]);
  }

  std::printf("\n");
  for (int i = 0; i < out_ref->m_n; i++) {
    ror_ref[i] = calcDensity(out_ref->m_t[i], 1, bwr);
    std::printf("%s  Reference density   %5.2f" U_DEG"C / 101.325 kPa :%s %s%8.4f%s [kg/m" U_SUP3 "]\n",
                COLOR_BOLD_WHITE, out_ref->m_t[i], COLOR_RESET, COLOR_BOLD_GREEN, ror_ref[i], COLOR_RESET);
    double z_ref = bwr.m_molarMass / (ror_ref[i] * GAS_CONSTANT_R * (out_ref->m_t[i] + KELVIN_OFFSET));
    std::printf("%s  Z factor            %5.2f" U_DEG"C / 101.325 kPa :%s %s%8.6f%s [" U_MDASH "]\n",
                COLOR_BOLD_WHITE, out_ref->m_t[i], COLOR_RESET, COLOR_BOLD_GREEN, z_ref, COLOR_RESET);
  }
}

void PrintFlowResults(DeviceType tip, double d_int, double d_orif,
                      double temperature, double pressure, double pressure_diff,
                      double ro, double eta, double qm,
                      const FlowResult& flow, const BwrConst& bwr,
                      const CountryRef& ref, const double* ror_ref) {
  std::printf("\n");
  sep_d();
  std::printf("%s  RESULTS  " U_MDASH "  %s%s\n",
              COLOR_YELLOW, tipName(tip), COLOR_RESET);
  std::printf("%s  D = %.2f mm  " U_CDOT "  d = %.2f mm  " U_CDOT "  "
              U_BETA " = %.4f%s\n",
              COLOR_BOLD_WHITE, d_int, d_orif, flow.m_beta, COLOR_RESET);
  std::printf("%s  D(t) = %.3f mm  " U_CDOT "  d(t) = %.3f mm  "
              U_CDOT "  t = %.1f" U_DEG"C%s\n",
              COLOR_BOLD_WHITE, flow.m_DWorkingMm, flow.m_dWorkingMm, temperature, COLOR_RESET);
  {
    double dt_ref = temperature - REF_TEMP_CELSIUS;
    double corr_D = THERMAL_EXP_PIPE    * dt_ref * 100.0;
    double corr_d = THERMAL_EXP_ORIFICE * dt_ref * 100.0;
    if (std::fabs(corr_D) > 0.10 || std::fabs(corr_d) > 0.10) {
      std::printf("%s  Note: thermal expansion at %.0f" U_DEG"C: "
                  U_DELTA"D = %+.3f%%,  " U_DELTA"d = %+.3f%%"
                  " " U_MDASH " verify dimensions are at 20" U_DEG"C reference.%s\n",
                  COLOR_YELLOW, temperature, corr_D, corr_d, COLOR_RESET);
    }
  }
  sep_d();

  std::printf("\n%s  Measurement conditions%s\n", COLOR_CYAN, COLOR_RESET);
  sep_s();
  row("Temperature", "%10.2f", temperature, U_DEG"C");
  row("Absolute pressure", "%10.2f", pressure, "kPa");
  row("Differential pressure", "%10.2f", pressure_diff, "kPa");

  double Z_tp = (pressure / KPA_PER_ATM) * bwr.m_molarMass
              / (ro * GAS_CONSTANT_R * (temperature + KELVIN_OFFSET));

  std::printf("\n%s  Fluid  (at t, p)%s\n", COLOR_CYAN, COLOR_RESET);
  sep_s();
  row("Density " U_RHO "(t,p)", "%10.4f", ro, "kg/m" U_SUP3);
  row("Dynamic viscosity " U_ETA "(t,p)", "%10.4f",
      eta * PA_TO_MICRO_PA, U_MICRO "Pa" U_CDOT "s");
  row("Compressibility factor Z(t,p)", "%10.6f", Z_tp, U_MDASH);

  std::printf("\n%s  Flow rates%s\n", COLOR_CYAN, COLOR_RESET);
  sep_s();
  row("Mass flow Qm", "%10.4f", qm, "kg/s");
  row("Mass flow Qm", "%10.2f", qm * SECONDS_PER_HOUR, "kg/h");
  for (int i = 0; i < ref.m_n; i++) {
    double qhref = SECONDS_PER_HOUR / ror_ref[i] * qm;
    char lbl[64];
    std::snprintf(lbl, sizeof(lbl),
                  "Volumetric %5.2f" U_DEG"C / 101.325 kPa", ref.m_t[i]);
    row(lbl, "%10.2f", qhref, ref.m_label[i]);
  }
  row("Volumetric at (t,p)", "%10.2f", SECONDS_PER_HOUR / ro * qm, "m" U_SUP3 "/h");

  std::printf("\n%s  Hydraulics%s\n", COLOR_CYAN, COLOR_RESET);
  sep_s();
  row("Mean velocity v", "%10.2f", flow.m_velocity, "m/s");
  row("Pressure loss", "%10.2f", flow.m_pressureLoss, "kPa");
  row("Throttle ratio " U_BETA, "%10.4f", flow.m_beta, U_MDASH);
  row("Reynolds number Re", "%10.4g", flow.m_reynolds, U_MDASH);
  row("Discharge coefficient C", "%10.6f", flow.m_coefC, U_MDASH);
  row("Expansibility factor " U_EPS, "%10.6f", flow.m_epsilon, U_MDASH);
  row("Isentropic exponent k", "%10.4f", flow.m_kappa, U_MDASH);

  std::printf("\n");
  sep_d();
}

void PrintCpuProfile(double t_mix_us, double t_rho_us, double t_eta_us,
                     double t_qm_us, int iter_rho, int iter_qm) {
  auto prow = [](const char* op, const char* detail, double us, int it) {
    int llen = (int)std::strlen(op) - Utf8ExtraBytes(op);
    std::printf("  %s%s", COLOR_BOLD_WHITE, op);
    for (int k = llen; k < 30; k++) std::putchar(' ');
    int dlen = (int)std::strlen(detail) - Utf8ExtraBytes(detail);
    std::printf("%s  %s%s", COLOR_RESET, detail, COLOR_RESET);
    for (int k = dlen; k < 28; k++) std::putchar(' ');
    std::printf("%s%7.2f%s [" U_MICRO "s]", COLOR_BOLD_GREEN, us, COLOR_RESET);
    if (it > 0) std::printf("  %s%d iter.%s", COLOR_YELLOW, it, COLOR_RESET);
    std::putchar('\n');
  };
  auto pcol = [](const char* s, int w) {
    std::printf("%s", s);
    int len = (int)std::strlen(s) - Utf8ExtraBytes(s);
    for (int k = len; k < w; k++) std::putchar(' ');
  };

  double t_total = t_mix_us + t_rho_us + t_eta_us + t_qm_us;
  std::printf("\n%s  CPU PROFILE%s\n", COLOR_CYAN, COLOR_RESET);
  sep_s();
  int n_comp = NUM_COMPONENTS;
  char mix_detail[48];
  std::snprintf(mix_detail, sizeof(mix_detail),
                "mixing rules N" U_SUP2 "+N" U_SUP3 "  (N=%d)", n_comp);
  prow("BWRS constants", mix_detail,          t_mix_us, 0);
  prow("Density " U_RHO,       "BWRS bisection",     t_rho_us, iter_rho);
  prow("Viscosity " U_ETA,     "Chapman-Enskog / N", t_eta_us, 0);
  prow("Flow rate Qm",    "Reynolds iteration", t_qm_us,  iter_qm);
  sep_s();
  std::printf("  %sTotal%s                                                       "
              "%s%7.2f%s [" U_MICRO "s]\n\n",
              COLOR_BOLD_WHITE, COLOR_RESET, COLOR_BOLD_GREEN, t_total, COLOR_RESET);

  std::printf("%s  EMBEDDED ESTIMATE%s  (indicative factors relative to PC)\n",
              COLOR_CYAN, COLOR_RESET);
  sep_s();
  struct { const char* proc; const char* freq; const char* fpu;
           const char* chip; double factor; } targets[] = {
    { "ARM Cortex-M7",  "@600 MHz", "+ FPU",     "(i.MX RT1062)",   6.0 },
    { "ARM Cortex-M7",  "@480 MHz", "+ FPU",     "(STM32H7)",       8.0 },
    { "Xtensa LX7",     "@240 MHz", "+ FPU",     "(ESP32-S3)",     12.0 },
    { "Xtensa LX6",     "@240 MHz", "+ FPU",     "(ESP32)",        16.0 },
    { "ARM Cortex-M4",  "@168 MHz", "+ FPU",     "(STM32F4)",      22.0 },
    { "ARM Cortex-M4",  "@120 MHz", "+ FPU",     "(STM32F3)",      32.0 },
    { "ARM Cortex-M33", "@ 64 MHz", "+ FPU",     "(nRF9160)",      55.0 },
    { "ARM Cortex-M3",  "@ 72 MHz", "no FPU",    "(STM32F103)",   200.0 },
    { "RISC-V RV32",    "@160 MHz", "no FPU",    "(ESP32-C3)",    180.0 },
    { "MSP430 F5xx",    "@ 25 MHz", "+ hw mult", "(TI)",          350.0 },
    { "AVR",            "@ 16 MHz", "no FPU",    "(Mega2560)",    900.0 },
    { "80C51",          "@ 12 MHz", "no FPU",    "",             3000.0 },
    { "Zilog Z80",      "@ 4-8 MHz","no FPU",    "(Elster)",     5000.0 },
  };
  std::printf("  %s", COLOR_BOLD_WHITE);
  pcol("Processor", 14); std::printf("  ");
  pcol("Freq.",      9); std::printf("  ");
  pcol("FPU",        9); std::printf("  ");
  pcol("Chip",      13);
  std::printf("  " U_TIMES "Mult  Est.[ms]%s\n", COLOR_RESET);
  sep_s();
  for (auto& tg : targets) {
    double est = t_total * tg.factor / 1000.0;
    std::printf("  %s", COLOR_BOLD_WHITE); pcol(tg.proc, 14);
    std::printf("%s  ", COLOR_RESET);     pcol(tg.freq,  9);
    std::printf("  ");               pcol(tg.fpu,   9);
    std::printf("  ");               pcol(tg.chip, 13);
    std::printf("  " U_TIMES "%4.0f  %s%8.2f%s [ms]\n",
                tg.factor, COLOR_BOLD_GREEN, est, COLOR_RESET);
  }
  sep_s();
  std::printf("%s  Note:%s estimated factors (IPC, cache, compiler)."
              " Measure on target for accuracy.\n",
              COLOR_YELLOW, COLOR_RESET);

  std::printf("\n");
  sep_d();
}

void PrintUncertainty(DeviceType tip, const FlowResult& flow, double qm,
                      const CountryRef& ref, const double* ror_ref) {
  std::printf("\n  %sCompute measurement uncertainty (ISO 5167-1)? [y/n]  %s>%s ",
              COLOR_BOLD_WHITE, COLOR_CYAN, COLOR_RESET);
  if (!AskYesNo()) return;

  double uC_def;
  switch (tip) {
    case DeviceType::ORIFICE_CORNER:
    case DeviceType::ORIFICE_FLANGE:
    case DeviceType::ORIFICE_DD2:
      uC_def = (flow.m_beta > 0.60) ? 0.75 : 0.50;
      break;
    case DeviceType::NOZZLE_ISA:            uC_def = 0.80; break;
    case DeviceType::NOZZLE_LONG_RADIUS:    uC_def = 2.00; break;
    case DeviceType::VENTURI_ROUGH_CAST:    uC_def = 0.70; break;
    case DeviceType::VENTURI_MACHINED:      uC_def = 1.00; break;
    case DeviceType::VENTURI_WELDED_SHEET:  uC_def = 1.50; break;
    case DeviceType::VENTURI_NOZZLE:        uC_def = 1.20; break;
    default:                                uC_def = 1.00; break;
  }

  auto read_pct = [&](const char* label, double def_val) -> double {
    double v = 0.0;
    int llen = (int)std::strlen(label) - Utf8ExtraBytes(label);
    std::printf("  %s>%s %s%s", COLOR_CYAN, COLOR_RESET, COLOR_BOLD_WHITE, label);
    for (int k = llen; k < 44; k++) std::putchar(' ');
    std::printf("%s[default %s%.2f%%%s] : %s",
                COLOR_BOLD_WHITE, COLOR_BOLD_GREEN, def_val, COLOR_BOLD_WHITE, COLOR_RESET);
    if (!ReadDouble(&v) || v <= 0.0) v = def_val;
    return v;
  };

  std::printf("\n%s  " U_HLINE U_HLINE " Measurement uncertainty  (ISO 5167-1)%s\n",
              COLOR_YELLOW, COLOR_RESET);
  std::printf(
      "%s  Relative standard uncertainties [%%]  "
      U_MDASH "  Enter or 0 " U_RARR " default:%s\n\n",
      COLOR_BOLD_WHITE, COLOR_RESET);
  double uC   = read_pct("u(C)   discharge coefficient (ISO 5167)",  uC_def);
  double uEps = read_pct("u(" U_EPS ")    expansibility factor",              0.10);
  double ud_p = read_pct("u(d)   orifice diameter",                      0.03);
  double uD_p = read_pct("u(D)   pipe diameter",                         0.10);
  double udp  = read_pct("u(" U_DELTA "p)  differential pressure",               0.20);
  double uro  = read_pct("u(" U_RHO ")   gas density (BWRS)",                 0.30);

  double b4   = std::pow(flow.m_beta, 4);
  double fd   = 2.0 / (1.0 - b4);
  double fD   = 2.0 * b4 / (1.0 - b4);
  double cC   = uC,        cEps = uEps;
  double cd   = fd * ud_p, cD   = fD * uD_p;
  double cdp  = 0.5 * udp, cro  = 0.5 * uro;
  double u_qm = std::sqrt(cC*cC + cEps*cEps + cd*cd + cD*cD + cdp*cdp + cro*cro);
  double U_qm = 2.0 * u_qm;

  std::printf("\n%s  Sensitivity factors  (" U_BETA " = %.4f)%s\n", COLOR_CYAN, flow.m_beta, COLOR_RESET);
  sep_s();
  std::printf(
      "  f_d = 2/(1" U_MINUS U_BETA U_SUP4 ") = %s%.4f%s     "
      "f_D = 2" U_BETA U_SUP4 "/(1" U_MINUS U_BETA U_SUP4 ") = %s%.4f%s\n\n",
      COLOR_BOLD_GREEN, fd, COLOR_RESET, COLOR_BOLD_GREEN, fD, COLOR_RESET);
  std::printf("  %s%-12s  %8s   %11s   %15s%s\n",
              COLOR_BOLD_WHITE, "Source", "u_i [%]", "Sensitivity", "Contribution [%]", COLOR_RESET);
  sep_s();

  auto brow = [&](const char* src, double ui, double fi, double ci) {
    int slen = (int)std::strlen(src) - Utf8ExtraBytes(src);
    std::printf("  %s%s", COLOR_BOLD_WHITE, src);
    for (int k = slen; k < 12; k++) std::putchar(' ');
    std::printf("%s  %s%8.4f%s   %11.4f   %s%15.4f%s\n",
                COLOR_RESET, COLOR_BOLD_GREEN, ui, COLOR_RESET, fi, COLOR_BOLD_GREEN, ci, COLOR_RESET);
  };
  brow("C",         uC,   1.0, cC);
  brow(U_EPS,  uEps, 1.0, cEps);
  brow("d",         ud_p, fd,  cd);
  brow("D",         uD_p, fD,  cD);
  brow(U_DELTA "p", udp,  0.5, cdp);
  brow(U_RHO,  uro,  0.5, cro);
  sep_s();
  row("Combined std. uncertainty  u(Qm)", "%10.4f", u_qm, "%");
  row("Expanded uncertainty  U(Qm)  k=2, 95%", "%10.4f", U_qm, "%");

  double U_abs_kgs = U_qm / 100.0 * qm;
  std::printf("\n");
  std::printf("  %s  Qm = %.4f " U_PLUSMN " %.5f  [kg/s]%s\n",
              COLOR_BOLD_GREEN, qm, U_abs_kgs, COLOR_RESET);
  std::printf("  %s  Qm = %.2f " U_PLUSMN " %.3f  [kg/h]%s\n",
              COLOR_BOLD_GREEN, qm * SECONDS_PER_HOUR, U_abs_kgs * SECONDS_PER_HOUR, COLOR_RESET);
  for (int i = 0; i < ref.m_n; i++) {
    if (ror_ref[i] > 0.0) {
      double qhref   = SECONDS_PER_HOUR / ror_ref[i] * qm;
      double U_qhref = U_qm / 100.0 * qhref;
      std::printf("  %s  Qv(%.2f" U_DEG"C) = %.2f " U_PLUSMN " %.3f  [%s]%s\n",
                  COLOR_BOLD_GREEN, ref.m_t[i], qhref, U_qhref, ref.m_label[i], COLOR_RESET);
    }
  }
  std::printf("\n");
  sep_d();
}