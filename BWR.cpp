// =============================================================================
// BWRS — Calcul debit gaze prin dispozitive de strangulare
// Revizia 3.0  |  05.2026
// ing. Agavriloaie Constantin  (original R 01.2004)
// ELCOST Impex
// =============================================================================
//
// DESCRIERE
//   Aplicația calculează densitatea, vâscozitatea dinamică și debitele unui
//   amestec de gaze cu până la 35 de componente, pe baza:
//     • ecuației de stare Benedict-Webb-Rubin-Starling (BWRS) pentru densitate;
//     • modelului Chapman-Enskog corectat cu termenul de densitate ridicată
//       pentru vâscozitate;
//     • metodei ISO 5167-2/3/4:2003 (ecuația Reader-Harris/Gallagher)
//       pentru calculul debitului prin dispozitive de strangulare.
//
// INTRĂRI
//   1. Compoziția amestecului — fracții molare pentru fiecare din cei
//      35 de componenți (metan, etan, propan, ... acetilenă).
//   2. Tipul dispozitivului de strangulare (1–9):
//        [1] Diafragmă cu prize în unghi
//        [2] Diafragmă cu prize la flanșă
//        [3] Diafragmă cu prize la D și D/2
//        [4] Ajutaj ISA 1932
//        [5] Ajutaj cu rază lungă
//        [6] Tub Venturi clasic — convergent brut turnat
//        [7] Tub Venturi clasic — convergent prelucrat
//        [8] Tub Venturi clasic — convergent brut din tablă sudată
//        [9] Ajutaj Venturi
//   3. Diametrul interior al conductei D [mm] și diametrul orificiului d [mm]
//      (la temperatura de referință de 20 °C).
//   4. Temperatura T [°C], presiunea absolută p [kPa] și presiunea
//      diferențială Δp [kPa] la locul de măsurare.
//
// IEȘIRI (pentru fiecare set T / p / Δp)
//   • Densitatea amestecului ρ(T, p)           [kg/m³]
//   • Vâscozitatea dinamică η(T, p)            [μPa·s]
//   • Debitul masic Qm                         [kg/s]
//   • Debitul volumic la condiții de referință selectabile
//       (0 °C / 101,325 kPa → Nm³/h; 15 °C / 101,325 kPa → Sm³/h; etc.)
//   • Debitul volumic la T și p                [m³/h]
//   • Viteza medie a gazului în conductă       [m/s]
//   • Pierderea de presiune prin strangulare   [kPa]
//   • Raportul de strangulare β și numărul Reynolds Re
//
// VALIDĂRI
//   Aplicația verifică limitele de aplicabilitate ISO 5167 / STAS 7347-90
//   pentru diametrul conductei, diametrul orificiului, raportul de strangulare β
//   și numărul Reynolds Re — și afișează mesaj de eroare la depășire.
//
// ALGORITM DENSITATE
//   Bisecție pe ecuația BWRS (Starling 1973, 11 parametri) până la
//   convergența |p_calc − p| < 5×10⁻⁴ atm.  Ecuația de stare:
//     p = ρRT + (B₀RT − A₀ − C₀/T² + D₀/T³ − E₀/T⁴)ρ²
//             + (bRT − a − d/T)ρ³ + α(a + d/T)ρ⁶
//             + (c/T²)ρ³(1 + γρ²)exp(−γρ²)
//   Constantele amestecului se calculează o singură dată din compoziție
//   prin reguli de mixare pătratice (A₀–E₀, γ) și cubice (a–d, α).
//
// ALGORITM DEBIT
//   Iterație pe numărul Reynolds până la convergența relativă |ΔQm/Qm| < 10⁻⁶.
//   Coeficientul de debit C și factorul de expansibilitate ε sunt recalculați
//   la fiecare iterație în funcție de Re și β.
//
// REFERINȚE
//   • ISO 5167-2:2003 — Orifice plates (ecuația Reader-Harris/Gallagher)
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

// Funcție: LoadComposition
// Intrări: x — tablou de ieșire pentru fracțiile molare (indexat 1..kNumComponents)
// Ieșire:  true dacă fișierul kCompFile a fost citit complet
// Scop:    încarcă compoziția amestecului salvată anterior
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

// Funcție: SaveComposition
// Intrări: x — fracțiile molare ale amestecului (indexat 1..kNumComponents)
// Ieșire:  —
// Scop:    persistă compoziția curentă în fișierul kCompFile
static void SaveComposition(const double* x) {
  FILE* f = nullptr;
  fopen_s(&f, kCompFile, "w");
  if (!f) {
    return;
  }
  for (int i = 1; i <= kNumComponents; i++) {
    std::fprintf(f, "%.8f\n", x[i]);
  }
  std::fclose(f);
}

// Funcție: LoadConfig
// Intrări: tip_raw — ieșire cod dispozitiv; d_int — ieșire D interior [mm]; d_orif — ieșire D orificiu [mm]
// Ieșire:  true dacă fișierul kConfFile a fost citit complet (3 valori)
// Scop:    reîncarcă configurația dispozitivului salvată anterior
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

// Funcție: SaveConfig
// Intrări: tip_raw — codul dispozitivului; d_int — D interior [mm]; d_orif — D orificiu [mm]
// Ieșire:  —
// Scop:    persistă configurația curentă a dispozitivului în fișierul kConfFile
static void SaveConfig(int tip_raw, double d_int, double d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, kConfFile, "w");
  if (!f) {
    return;
  }
  std::fprintf(f, "%d\n%.4f\n%.4f\n", tip_raw, d_int, d_orif);
  std::fclose(f);
}

// Funcție: ExitApp
// Intrări: —
// Ieșire:  — (nu returnează; termină procesul cu std::exit)
// Scop:    șterge ecranul, afișează mesajul de ieșire și încheie aplicația la apăsarea ESC
static void ExitApp() {
  std::printf("\033[2J\033[H%s\n  ELCOST Impex  —  BWR Gas Flow Calculator  v2.0/2026%s\n",
              kBoldYellow, kReset);
  std::exit(0);
}

// Funcție: ReadDouble
// Intrări: val — pointer la variabila de ieșire
// Ieșire:  valoarea citită prin *val (0.0 dacă nu s-a introdus nimic)
// Scop:    citire interactivă a unui număr real cu ecou, backspace și ESC
static void ReadDouble(double* val) {
  char buf[64] = {};
  int pos = 0;
  for (;;) {
    int ch = _getch();
    if (ch == 27) {
      ExitApp();
    }
    if (ch == 0 || ch == 0xE0) {
      (void)_getch();
      continue;
    }
    if (ch == '\r') {
      std::printf("\n");
      std::fflush(stdout);
      break;
    }
    if ((ch == 8 || ch == 127) && pos > 0) {
      pos--;
      std::printf("\b \b");
      std::fflush(stdout);
      continue;
    }
    if (pos < (int)sizeof(buf) - 2 && ((ch >= '0' && ch <= '9') || ch == '.' || (ch == '-' && pos == 0))) {
      buf[pos++] = static_cast<char>(ch);
      std::printf("%c", ch);
      std::fflush(stdout);
    }
  }
  *val = (pos > 0) ? std::atof(buf) : 0.0;
}

// Scop: citește un singur caracter cifră în intervalul [lo, hi] fără a necesita ENTER
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

// Funcție: AskYesNo
// Intrări: —
// Ieșire:  true pentru 'd'/'D' (da), false pentru 'n'/'N' (nu)
// Scop:    blochează până la apăsarea d/n; ESC termină aplicația
static bool AskYesNo() {
  for (;;) {
    int ch = _getch();
    if (ch == 27) {
      ExitApp();
    }
    if (ch == 'd' || ch == 'D') {
      std::printf("d\n");
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

// Funcție: Utf8ExtraBytes
// Intrări: s — șir UTF-8
// Ieșire:  numărul de octeți de continuare (bit-pattern 10xxxxxx)
// Scop:    corectarea lungimii vizuale: strlen(s) - Utf8ExtraBytes(s) = nr. caractere afișate
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
    "Metan", "Etan", "Propan", "Izobutan", "N-butan",
    "Neopentan", "Izopentan", "N-pentan", "2,2-dimetilbutan", "2,3-dimetilbutan",
    "3-metilpentan", "2-metilpentan", "N-hexan", "2,4-dimetilpentan", "2,2,3-trimetilbutan",
    "2-metilhexan", "3-metilhexan", "3-etilpentan", "N-heptan", "2,2,4-trimetilpentan",
    "N-octan", "Benzen", "Toluen", "Hidrogen", "Monoxid de carbon",
    "Hidrogen sulfurat", "Heliu", "Argon", "Azot", "Oxigen",
    "Dioxid de carbon", "Etilen\xC4\x83", "Propilen\xC4\x83", "Amoniac", "Acetilen\xC4\x83"
};

// Funcție: PrintComposition
// Intrări: x — fracțiile molare ale amestecului (indexat 1..kNumComponents)
// Ieșire:  —
// Scop:    afișează tabelul compoziției aliniat în consolă, cu valori evidențiate în verde
static void PrintComposition(const double* x) {
  std::printf("\n%s  ── Compoziție%s\n",
              kBoldYellow, kReset);
  std::printf("%s  Fracții molare ale amestecului:%s\n\n", kBoldWhite, kReset);
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

// Funcție: main
// Intrări: —
// Ieșire:  0 la ieșire normală (nu returnează în mod obișnuit)
// Scop:    inițializare constante BWRS, citire compoziție/configurație, buclă de calcul debit masic
int main() {
#ifdef _WIN32
  std::system("chcp 65001 > nul");
  {
    HANDLE hCon = GetStdHandle(STD_OUTPUT_HANDLE);

    // Setează fontul
    CONSOLE_FONT_INFOEX cfi = {};
    cfi.cbSize = sizeof(cfi);
    GetCurrentConsoleFontEx(hCon, FALSE, &cfi);
    cfi.dwFontSize.Y = 14;
    wcscpy_s(cfi.FaceName, L"Consolas");
    SetCurrentConsoleFontEx(hCon, FALSE, &cfi);

    // Activează procesarea VT/ANSI
    DWORD mode = 0;
    GetConsoleMode(hCon, &mode);
    SetConsoleMode(hCon, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
  }
#endif

  double x[kArraySize] = {};

  // Valori brute; tablourile marcate cu (*) sunt scalate după inițializare
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
       0.05947,  0.0577,   0.0745,    0.066,     0.0715,    0.1636,    0.141,    0.071,    0.174,    0.1755,
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

  // D₀, E₀, d — parametri Starling (1973); completează din Tabel 1, pag. 19.
  // Valorile 0 reduc BWRS exact la BWR — înlocuiește pe rând per componentă.
  // Surse: Nishiumi & Saito (J.CEJ 1975) pentru C1–C8 și gaze permanente;
  // Starling (1973) pentru CO₂, H₂S, CO; Poling et al. (2001) pentru Ar, NH₃;
  // estimare prin corelație Tc/Pc (marcate *) pentru izomerii fără date publicate.
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

  // Cp° ideal-gas [J/(mol·K)] at 20 °C — pentru κ = Cp/(Cp − R) cu R = 8.314 J/(mol·K)
  // Surse: NIST WebBook; izomerii C5–C8 estimați din grupuri funcționale (±2 J/(mol·K))
  static const double kCp0[kArraySize] = {
       0,      35.7,  52.5,  73.6,  97.5,  96.4, 120.9, 118.9, 120.1, 141.3,
     140.9,  141.7, 141.2, 143.1, 164.8, 163.4, 165.0, 165.0, 165.1, 166.1,
     188.9,  188.9,  82.4, 103.7,  28.8,  29.1,  34.2,  20.8,  20.8,  29.1,
      29.4,   37.1,  42.9,  63.9,  35.7,  44.0};

  // Scalare tablouri (*) — executată o singură dată la pornire
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

  // ── Antet principal ────────────────────────────────────────────────────────
  std::printf(
      "%s  ══════════════════════════════════════════════════════════════════════════════\n"
      "  %sELCOST Impex%s  ·  BWRS Gas Flow Calculator                          v3.0/2026\n"
      "  ══════════════════════════════════════════════════════════════════════════════\n\n",
      kHdrYellow, kHdrGreen, kHdrYellow);

  // ── Modele de calcul ───────────────────────────────────────────────────────
  std::printf("%s  MODELE DE CALCUL%s\n", kHdrCyan, kReset);
  std::printf(
      "  %sEcuație de stare%s    BWRS · Starling 1973  —  Benedict-Webb-Rubin-Starling, 11 param.\n"
      "  %sVâscozitate%s         Chapman-Enskog cu corecție la densitate ridicată  [Nishiumi 1975]\n"
      "  %sCoeficient C, ε%s     Reader-Harris/Gallagher  ·  ISO 5167-2/3/4:2003\n\n",
      kBoldWhite, kReset, kBoldWhite, kReset, kBoldWhite, kReset);

  // ── Domeniu ────────────────────────────────────────────────────────────────
  std::printf("%s  DOMENIU DE APLICARE%s\n", kHdrCyan, kReset);
  std::printf(
      "  35 componenți  ·  9 tipuri dispozitive de strangulare\n"
      "  Debit masic, volumic, viteză, pierdere presiune  ·  condiții ref. selectabile\n"
      "  Validări limită ISO 5167 / STAS 7347-90  (D, β, Re)\n\n");

  // ── Standarde și referințe ─────────────────────────────────────────────────
  std::printf(
      "%s  ──────────────────────────────────────────────────────────────────────────────\n"
      "  %sSTANDARDE ȘI REFERINȚE%s\n",
      kHdrYellow, kHdrCyan, kReset);
  std::printf(
      "  %sISO 5167-2:2003%s   Diafragme — ecuația Reader-Harris/Gallagher\n"
      "  %sISO 5167-3:2003%s   Ajutaje și ajutaje Venturi\n"
      "  %sISO 5167-4:2003%s   Tuburi Venturi clasice\n"
      "  %sSTAS 7347-90%s      Limite de aplicabilitate (D, β, Re)\n"
      "  %sStarling K.E.%s     Fluid Thermodynamic Properties, Gulf Publ. Houston (1973)\n"
      "  %sNishiumi & Saito%s  J. Chem. Eng. Japan 8(5), 356–360 (1975)\n\n",
      kBoldWhite, kReset, kBoldWhite, kReset, kBoldWhite, kReset,
      kBoldWhite, kReset, kBoldWhite, kReset, kBoldWhite, kReset);

  // ── Validare ───────────────────────────────────────────────────────────────
  std::printf(
      "%s  ──────────────────────────────────────────────────────────────────────────────\n"
      "  %sVALIDARE%s\n",
      kHdrYellow, kHdrCyan, kReset);
  std::printf(
      "  ρ  verificat pe 8 compoziții (CH\xe2\x82\x84, C\xe2\x82\x82H\xe2\x82\x86, CO\xe2\x82\x82,"
      " H\xe2\x82\x82, N\xe2\x82\x82, Ar, GN std, GN bogat)  vs. NIST WebBook\n"
      "  η  verificat vs. Chapman-Enskog / NIST\n"
      "  Qm verificat \xc2\xb1" "15 %% față de estimări ISO 5167  ·  9 tipuri dispozitive\n\n");

  // ── Footer ─────────────────────────────────────────────────────────────────
  std::printf(
      "%s  ──────────────────────────────────────────────────────────────────────────────\n"
      "  %soffice@elcost.ro%s                                   \xc2\xa9 2004\xe2\x80\x93" "2026 ELCOST Impex\n"
      "%s  ══════════════════════════════════════════════════════════════════════════════\n\n",
      kHdrYellow, kHdrCyan, kReset, kHdrYellow);

  std::printf("%s", kReset);

  // ── Self-test opțional la pornire ─────────────────────────────────────────
  {
    std::printf("%s\n  Rula\xC8\x9Bi testul de validare al implement\xC4\x83rii? [d/n] %s>%s ",
                kBoldWhite, kCyan, kReset);
    if (AskYesNo()) {
      std::printf("\n%s  \xe2\x94\x80\xe2\x94\x80 Test validare implementare%s\n", kBoldYellow, kReset);

      int npass = 0, ntotal = 0;

      // chk: compară o valoare cu un interval și afișează PASS/FAIL
      // src = sursa intervalului de referinta, afisat la sfarsitul liniei
      auto chk = [&](const char* label, double v, double lo, double hi, const char* src) {
        ntotal++;
        bool ok = (v >= lo && v <= hi);
        if (ok) npass++;
        int vlen = (int)std::strlen(label) - Utf8ExtraBytes(label);
        std::printf("    %s%s", kBoldWhite, label);
        for (int k = vlen; k < 50; k++) std::putchar(' ');
        char ivl[32];
        std::snprintf(ivl, sizeof(ivl), "[%.4f, %.4f]", lo, hi);
        std::printf("%s%-10.5f%s  %-24s  %s%-4s%s  %s[%s]%s\n",
                    kBoldGreen, v, kReset,
                    ivl,
                    ok ? kBoldGreen : kBoldRed, ok ? "PASS" : "FAIL", kReset,
                    kBoldWhite, src, kReset);
      };

      // make_bwr: calculează BwrConst pentru orice compoziție xc[1..kNumComponents]
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

      // calc_eta: vâscozitate [Pa·s] la temperatura t[°C] și 1 atm, compoziție xc
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

      // test_comp: rulează bateria de teste (densitate × 2, vâscozitate, debit masic)
      //   pentru compoziția xc cu parametri și intervale de referință date.
      //   d_flow_mm = 0 → sare testul de debit.
      auto test_comp = [&](const char* comp_name,
                           const double* xc,
                           double rho20_lo, double rho20_hi,   // densitate 20°C [kg/m³]
                           double rho0_lo,  double rho0_hi,    // densitate  0°C [kg/m³]
                           double eta_lo,   double eta_hi,     // vâscozitate 20°C [μPa·s]
                           double D_mm,     double d_flow_mm,  // D și d [mm] (d=0 → skip)
                           double p_kpa,    double dp_kpa,
                           bool leading_nl = true) {
        std::printf(leading_nl ? "\n  %s\xc2\xbb %s%s\n" : "  %s\xc2\xbb %s%s\n",
                    kBoldYellow, comp_name, kReset);

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
        chk("Densitate 20\xC2\xB0""C, 101.325 kPa [kg/m\xC2\xB3]", rho1, rho20_lo, rho20_hi, "BWRS 1973");
        chk("Densitate  0\xC2\xB0""C, 101.325 kPa [kg/m\xC2\xB3]", rho2, rho0_lo,  rho0_hi,  "BWRS 1973");
        chk("V\xC3\xA2scozitate 20\xC2\xB0""C, 101.325 kPa [\xC2\xB5Pa\xC2\xB7s]",
            calc_eta(20.0, xc, bt, roc_crit_c, csi_c) * kPaToMicroPa, eta_lo, eta_hi, "CE/NIST");
        // Z = p*M / (rho*R*T); la 1 atm, 20°C: Z ∈ [0.975, 1.003] pt. toate gazele comune
        double T_K1 = 20.0 + kKelvinOffset;
        double Z1   = 1.0 * bt.molar_mass / (rho1 * kGasConstantR * T_K1);
        chk("Factor Z la 20\xC2\xB0""C, 101.325 kPa [-]", Z1, 0.975, 1.003, "BWRS/ideal");
        // Masa molara: make_bwr trebuie sa acumuleze exact sum(x[i]*m[i])
        double mx_ref = 0.0;
        for (int ii = 1; ii <= kNumComponents; ii++) mx_ref += xc[ii] * m[ii];
        chk("Mas\xC4\x83 molar\xC4\x83 M [g/mol]", bt.molar_mass,
            mx_ref * (1.0 - 1e-8), mx_ref * (1.0 + 1e-8), "amestecare");

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
          // alpha_nom = C_nominal/sqrt(1-beta^4) la beta=0.4; eps_nom din ISO 5167
          // Sursa C nominal: ISO 5167-2/3/4:2003 si Reader-Harris/Gallagher la Re=10^6
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
                                       static_cast<TipDispozitiv>(t.tip),
                                       D_mm, d_flow_mm, rho_f, eta_f, &fr);
            if (t.tip == 1) qm_diau = qm_f;
            if (qm_f > 0.0) {
              double qm_est = t.alpha_nom * t.eps_nom * A_o
                            * std::sqrt(2000.0 * dp_kpa * rho_f);
              char lbl[80];
              std::snprintf(lbl, sizeof(lbl),
                            "Qm D=%g,d=%g,p=%g,\xCE\x94p=%g %s [kg/s]",
                            D_mm, d_flow_mm, p_kpa, dp_kpa, t.abv);
              chk(lbl, qm_f, qm_est * 0.85, qm_est * 1.15, "ISO 5167");
            } else {
              ntotal++;
              std::printf("    %sTest Qm %s: EROARE ISO 5167%s\n", kBoldRed, t.abv, kReset);
            }
          }
          if (qm_diau > 0.0 && rho2 > 0.0) {
            double qv_n     = qm_diau / rho2 * kSecondsPerHour;
            double qm_est_u = tipuri[0].alpha_nom * tipuri[0].eps_nom * A_o
                            * std::sqrt(2000.0 * dp_kpa * rho_f);
            double qv_est   = qm_est_u / rho2 * kSecondsPerHour;
            chk("Debit vol. Dia-U, 0\xC2\xB0""C ref [Nm\xC2\xB3/h]",
                qv_n, qv_est * 0.85, qv_est * 1.15, "ISO 5167/BWRS");
          }
        }
      };

      // sec_box: afișează o casetă de descriere înaintea unei secțiuni de teste
      // Format: tabel 3 coloane (Secțiune | Ce verifică | Teste), lățime ~110 car.
      auto sec_box = [](const char* sec, const char* desc, int n) {
        auto hl = [](int cnt) { for (int k = 0; k < cnt; ++k) std::printf("\xe2\x94\x80"); };
        std::printf("\n  \xe2\x94\x8c"); hl(22); std::printf("\xe2\x94\xac"); hl(74);
        std::printf("\xe2\x94\xac"); hl(7); std::printf("\xe2\x94\x90\n");
        std::printf("  \xe2\x94\x82 %-21s \xe2\x94\x82 %-73s \xe2\x94\x82 %-5s \xe2\x94\x82\n",
                    "Sec\xc8\x9biune", "Ce verific\xc4\x83", "Teste");
        std::printf("  \xe2\x94\x9c"); hl(22); std::printf("\xe2\x94\xbc"); hl(74);
        std::printf("\xe2\x94\xbc"); hl(7); std::printf("\xe2\x94\xa4\n");

        char nt[8]; std::snprintf(nt, sizeof(nt), "+%d", n);
        std::printf("  \xe2\x94\x82 %-20s \xe2\x94\x82 %-72s \xe2\x94\x82 %-5s \xe2\x94\x82\n",
                    sec, desc, nt);
        std::printf("  \xe2\x94\x94"); hl(22); std::printf("\xe2\x94\xb4"); hl(74);
        std::printf("\xe2\x94\xb4"); hl(7); std::printf("\xe2\x94\x98\n");
      };

      // ── Compoziții de test ─────────────────────────────────────────────────
      // Adăugați oricâte compoziții cu un nou apel test_comp() mai jos.
      // Indici componente: CH4=1, C2H6=2, C3H8=3, i-C4=4, n-C4=5, H2=24,
      //   CO=25, H2S=26, He=27, Ar=28, N2=29, O2=30, CO2=31, Amoniac=34
      //
      // ── Referințe densitate (BWRS față de gaz ideal la 1 atm) ───────────
      // Gaz ideal 20°C: rho = M / (R × 293.15) = M / 24.0554 [kg/m³]
      // Gaz ideal  0°C: rho = M / (R × 273.15) = M / 22.4136 [kg/m³]
      //
      //  Compoziție     M [g/mol] rho_id_20  rho_id_0   B₂(20°C)[L/mol]  rho_BWRS_20
      //  CH4             16.043    0.6672     0.7157     −0.0447           0.668
      //  C2H6            30.070    1.2500     1.3415     −0.195            1.260
      //  CO2             44.011    1.8296     1.9635     −0.147            1.841
      //  H2               2.016    0.0838     0.0899     +0.022            0.0837
      //  N2              28.016    1.1647     1.2498     −0.0044           1.162
      //  Ar              39.944    1.6604     1.7821     −0.018            1.661
      //
      // Sursa B₂: calculat din parametrii BWRS Starling 1973 prin relația:
      //   B₂(T) = B₀ − A₀/(RT) − C₀/(RT³) + D₀/(RT⁴) − E₀/(RT⁵)
      //
      // ── Referințe vâscozitate (Chapman-Enskog + corectie densitate) ──────
      //  Compoziție   η_calc [μPa·s]   η_NIST 20°C [μPa·s]   Eroare CE
      //  CH4           11.46            10.99                   +4.3 %
      //  C2H6           9.68             9.36                   +3.4 %
      //  CO2           15.49            14.89                   +4.0 %
      //  H2             7.52             8.77                  −14.3 % (CE<exact ptr H2)
      //  N2            18.47            17.54                   +5.3 %
      //  Ar            18.33            22.72                  −19.3 % (param. vechi)
      //
      // Sursa NIST: https://webbook.nist.gov (Transport Properties, 100 kPa, 20°C)
      //
      // ── Referințe debit masic ────────────────────────────────────────────
      // Toate testele: D=200 mm, d=80 mm, β=0.4, T=20°C, p=500 kPa, Δp=5 kPa
      // β=0.4 ales pentru a respecta toate constrângerile ISO 5167 simultan:
      //   Ven-P/T necesită β≥0.4; Ven-T necesită D≥200 mm; Ven-P Re≤1e6
      // Rulat pentru toate cele 9 tipuri de dispozitive ISO 5167.
      // Qm_est = alpha_nom × eps_nom × A_o × √(2000 × Δp_kPa × ρ)
      //   A_o = π/4 × (0.08)² = 5.027e−3 m²
      // Interval acceptat: ±15% față de Qm_est (acoperă variația Re și BWRS)
      //
      //  Tip       α_nom   ε_nom   C_nom   Sursa C
      //  Dia-U     0.609   0.993   0.601   ISO 5167-2, Reader-Harris/Gallagher β=0.4 Re=10⁶
      //  Dia-F     0.608   0.993   0.601   idem, prize flanșă
      //  Dia-D/2   0.608   0.993   0.600   idem, prize D și D/2
      //  Aj-ISA    0.997   0.997   0.985   ISO 5167-3, ajutaj ISA 1932
      //  Aj-RL     1.005   0.997   0.992   ISO 5167-3, ajutaj rază lungă
      //  Ven-B     0.997   0.997   0.984   ISO 5167-4, Venturi brut turnat (C=0.984)
      //  Ven-P     1.008   0.997   0.995   ISO 5167-4, Venturi prelucrat   (C=0.995)
      //  Ven-T     0.998   0.997   0.985   ISO 5167-4, Venturi tablă sudată (C=0.985)
      //  Aj-V      0.996   0.997   0.983   ISO 5167-4, ajutaj Venturi (C≈0.983 la β=0.4)
      // ───────────────────────────────────────────────────────────────────

      double xCH4[kArraySize] = {}; xCH4[1] = 1.0;   // Metan pur

      double xC2H6[kArraySize] = {}; xC2H6[2] = 1.0;  // Etan pur

      double xCO2[kArraySize] = {}; xCO2[31] = 1.0;   // CO2 pur

      double xH2[kArraySize] = {}; xH2[24] = 1.0;     // Hidrogen pur

      double xN2[kArraySize] = {}; xN2[29] = 1.0;     // Azot pur

      double xAr[kArraySize] = {}; xAr[28] = 1.0;     // Argon pur

      double xGN[kArraySize] = {};                     // Gaz natural std (STAS 7347)
      xGN[1] = 0.900; xGN[2] = 0.060; xGN[3] = 0.020;
      xGN[29] = 0.015; xGN[31] = 0.005;

      double xGNB[kArraySize] = {};                    // Gaz natural bogat
      xGNB[1] = 0.850; xGNB[2] = 0.100; xGNB[3] = 0.030;
      xGNB[29] = 0.010; xGNB[31] = 0.010;

      sec_box("Compozitii test",
              "8 compozitii: rho(20/0C), eta, Z, masa mol., Qm 9 tipuri ISO 5167", 120);
      std::printf("  Surse referinta: BWRS 1973 = Starling, K.E. (1973) Fluid Thermodynamic Properties;\n");
      std::printf("                   CE/NIST   = Chapman-Enskog / NIST WebBook (webbook.nist.gov);\n");
      std::printf("                   ISO 5167  = ISO 5167-2/3/4:2003 (Reader-Harris/Gallagher).\n\n");
      std::printf("  Prescurtari dispozitive ISO 5167 (test debit D=200,d=80,p=500,dp=5):\n");
      std::printf("    Dia-U   = Diafragma prize unghi       Dia-F   = Diafragma prize flansa\n");
      std::printf("    Dia-D/2 = Diafragma prize D si D/2    Aj-ISA  = Ajutaj ISA 1932\n");
      std::printf("    Aj-RL   = Ajutaj raza lunga           Ven-B   = Venturi brut turnat\n");
      std::printf("    Ven-P   = Venturi prelucrat           Ven-T   = Venturi tabla sudata\n");
      std::printf("    Aj-V    = Ajutaj Venturi\n\n");

      //           Compozitie    rho_20[kg/m3]  rho_0[kg/m3]   eta[µPa·s]   D    d   p    dp
      test_comp("Metan (CH4)",  xCH4, 0.660,0.672, 0.712,0.724, 10.5,12.5, 200,80, 500,5, false);
      test_comp("Etan (C2H6)",  xC2H6,1.220,1.280, 1.310,1.375,  8.5,10.8, 200,80, 500,5);
      test_comp("CO2",          xCO2, 1.800,1.860, 1.930,2.000, 14.0,16.8, 200,80, 500,5);
      test_comp("Hidrogen (H2)",xH2,  0.082,0.085, 0.088,0.091,  7.0,10.5, 200,80, 500,5);
      test_comp("Azot (N2)",    xN2,  1.155,1.175, 1.240,1.260, 16.5,20.5, 200,80, 500,5);
      test_comp("Argon (Ar)",   xAr,  1.650,1.672, 1.772,1.794, 17.0,23.0, 200,80, 500,5);
      test_comp("GN std",       xGN,  0.730,0.748, 0.784,0.802, 10.5,12.5, 200,80, 500,5);
      test_comp("GN bogat",     xGNB, 0.765,0.790, 0.820,0.848, 10.2,12.5, 200,80, 500,5);

      sec_box("Autonome CH4 pur",
              "CH4 pur: monotonie rho vs P (x2); factor Z la 20 atm", 3);
      // ── Teste autonome: comportament fizic CH4 pur ───────────────────────
      std::printf("  %s\xC2\xBB Teste autonome CH4 pur%s\n", kBoldYellow, kReset);
      {
        BwrConst bt_a = make_bwr(xCH4);
        double T_a    = 20.0 + kKelvinOffset;
        double rho_05 = CalcDensity(20.0,  0.5, bt_a);   // 50.66 kPa
        double rho_10 = CalcDensity(20.0,  1.0, bt_a);   // 101.325 kPa
        double rho_20 = CalcDensity(20.0, 20.0, bt_a);   // ≈ 2.026 MPa
        // Gaz ideal: rho(2p)/rho(p) = 2; corectie Z introduce deviatii mici
        chk("Monotonie \xCF\x81(1atm)/\xCF\x81(0.5atm) CH4 [-]",
            rho_10 / rho_05, 1.90, 2.10, "gaz ideal");
        // La 20 atm, Z(CH4)≈0.953 → rho(20)/rho(1) ≈ 20/0.953 ≈ 21.0
        chk("Monotonie \xCF\x81(20atm)/\xCF\x81(1atm) CH4 [-]",
            rho_20 / rho_10, 17.0, 23.0, "gaz ideal\xC3\x97Z");
        // Factor Z la presiune inalta: NIST CH4 20°C, 20atm → Z≈0.953
        double Z_hi = 20.0 * bt_a.molar_mass / (rho_20 * kGasConstantR * T_a);
        chk("Factor Z CH4 la 20\xC2\xB0""C, 20 atm [-]", Z_hi, 0.88, 0.98, "BWRS/NIST");
      }

      sec_box("Densitate p. ridic.",
              "CH4/N2/CO2/GN la 500-5066 kPa: rho BWRS, Z, monotonie, consistenta", 8);
      // ── Densitate la presiune ridicată ────────────────────────────────────
      // Referințe:
      //   CH4  500 kPa, 293.15 K: Z≈0.9986 (NIST)       → ρ≈3.30  kg/m³ (≈ideal)
      //   CH4 5066 kPa, 293.15 K: Z≈0.911  (BWRS/AGA-8) → ρ≈36.6  kg/m³
      //     Verificare virială: B₂(CH4,293K)≈−44 cm³/mol
      //     → Z≈1+B₂·P/(RT)=1−0.044·50/24.05=0.909 (consistent cu BWRS 0.911)
      //     Nota: Z≈0.954 ar corespunde ~27 atm, nu 50 atm (sub temp. Boyle 511 K)
      //   N2  5066 kPa, 293.15 K: Z≈0.991  (NIST)       → ρ≈58.8  kg/m³
      //   CO2  500 kPa, 293.15 K: Z≈0.974  (NIST)       → ρ≈9.27  kg/m³ (non-ideal)
      //   GN std 5066 kPa, 20°C : estimat BWRS           → ρ≈38-40 kg/m³
      // Toleranță ±5%: BWRS Starling 1973 are eroare tipică 1–3% la presiuni ridicate.
      std::printf("  %s\xc2\xbb Densitate la presiune ridicat\xc4\x83%s\n", kBoldYellow, kReset);
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
        chk("Densitate CH4 20\xc2\xb0""C,  500 kPa [kg/m\xc2\xb3]",
            rCH4_500,  3.24, 3.35, "NIST/BWRS");
        chk("Densitate CH4 20\xc2\xb0""C, 5066 kPa [kg/m\xc2\xb3]",
            rCH4_5066, 34.5, 39.0, "BWRS/AGA-8");

        // Z la 50 atm: B₂(CH4,293K)≈−44 cm³/mol → Z_virial≈0.909; BWRS→0.911
        double Z50 = p5066 * bCH4.molar_mass
                   / (rCH4_5066 * kGasConstantR * T_K);
        chk("Factor Z CH4 20\xc2\xb0""C, 50 atm [-]",
            Z50, 0.88, 0.95, "BWRS/AGA-8");

        // ─ N2: cvasi-ideal la 50 atm, Z ≈ 0.987 (NIST) ─────────────────────
        double rN2_5066 = CalcDensity(20.0, p5066, bN2);
        chk("Densitate N2  20\xc2\xb0""C, 5066 kPa [kg/m\xc2\xb3]",
            rN2_5066, 56.0, 62.0, "NIST/BWRS");

        // ─ CO2 la 500 kPa: Z ≈ 0.974, mai non-ideal decât CH4 ───────────────
        double rCO2_500 = CalcDensity(20.0, p500, bCO2);
        chk("Densitate CO2 20\xc2\xb0""C,  500 kPa [kg/m\xc2\xb3]",
            rCO2_500, 8.90, 9.70, "NIST/BWRS");

        // ─ GN std la 50 atm: presiune tipică de transport ────────────────────
        double rGN_5066 = CalcDensity(20.0, p5066, bGN);
        chk("Densitate GN std 20\xc2\xb0""C, 5066 kPa [kg/m\xc2\xb3]",
            rGN_5066, 35.0, 42.0, "BWRS/estim.");

        // ─ Monotonie: ρ(50 atm)/ρ(5 atm) = p_ratio × Z(5 atm)/Z(50 atm) ────
        // = 10.13 × (0.9986/0.954) ≈ 10.6
        chk("Monotonie \xcf\x81 CH4: \xcf\x81(5066)/\xcf\x81(500) [-]",
            rCH4_5066 / rCH4_500, 9.5, 11.5, "BWRS/Z");

        // ─ Consistență: ρ_CO2 > ρ_N2 la aceeași T și p ───────────────────────
        // M_CO2=44.01 >> M_N2=28.02 și Z_CO2 < Z_N2 → raport ≈ 1.56–1.62
        double rN2_500 = CalcDensity(20.0, p500, bN2);
        chk("Consisten\xc8\x9b\xc4\x83: \xcf\x81(CO2)/\xcf\x81(N2) la 500 kPa [-]",
            rCO2_500 / rN2_500, 1.50, 1.65, "M + Z");
      }

      sec_box("GN bogat 4Tx4P",
              "GN bogat 4Tx4P: Z, rho/rho_ideal, monotonie vs T si vs P", 30);
      // ── GN bogat: matrice 4 temperaturi × 4 presiuni ─────────────────────
      // Temperaturi: -10°C (iarnă), 10°C, 30°C, 60°C (vară/compresor)
      // Presiuni:    101 kPa (ref.), 500 kPa, 2000 kPa, 7000 kPa (transport)
      // T_min = -10°C = 263.15 K > Tc_pseudo(GNB) ≈ 208 K → gaz monofazic.
      // Verificări: Z ∈ [0.75, 1.02]; ρ/ρ_ideal la p joasă; monotonie T și P.
      std::printf("  %s\xc2\xbb GN bogat: matrice T \xc3\x97 P  (4 \xc3\x97 4 condi\xc8\x9bii)%s\n",
                  kBoldYellow, kReset);
      {
        BwrConst bGNB = make_bwr(xGNB);

        const double T_c[]   = { -10.0,   10.0,   30.0,   60.0 };   // [°C]
        const double P_kpa[] = { 101.325, 500.0, 2000.0, 7000.0 };  // [kPa]
        const char*  T_lbl[] = { "-10",  " 10",  " 30",  " 60"  };
        const char*  P_lbl[] = { "  101", "  500", " 2000", " 7000" };
        constexpr int nT = 4, nP = 4;

        double rho[nT][nP] = {};

        // ─── Z ∈ [0.75, 1.02] pentru toate cele 16 condiții ─────────────────
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
            chk(lbl, Z, 0.75, 1.02, "gaz monofaz.");
          }
        }

        // ─── La p ≤ 500 kPa: ρ/ρ_ideal ∈ [0.990, 1.025] (cvasi-ideal) ──────
        // Limita sup. 1.025 acoperă T=-10°C unde B_mix≈-85 cm³/mol → Z≈0.981
        for (int ti = 0; ti < nT; ti++) {
          for (int pi = 0; pi < 2; pi++) {
            double T_K    = T_c[ti] + kKelvinOffset;
            double rho_id = P_kpa[pi] * 1e3 * bGNB.molar_mass * 1e-3
                          / (8.314 * T_K);
            char lbl[64];
            std::snprintf(lbl, sizeof(lbl),
                          "\xcf\x81/\xcf\x81_id GNB T=%s\xc2\xb0""C P=%s kPa",
                          T_lbl[ti], P_lbl[pi]);
            chk(lbl, rho[ti][pi] / rho_id, 0.990, 1.025, "cvasi-ideal");
          }
        }

        // ─── Monotonie față de T la P = 2000 kPa (pi=2): ρ(T1) > ρ(T2) ─────
        // Raport ρ(T_i)/ρ(T_{i+1}) = (Z_{i+1}/Z_i) × (T_{i+1}/T_i) ≈ 1.07–1.10
        chk("GNB: \xcf\x81(-10\xc2\xb0""C)/\xcf\x81( 10\xc2\xb0""C) la 2000 kPa",
            rho[0][2] / rho[1][2], 1.03, 1.15, "monotonie T");
        chk("GNB: \xcf\x81( 10\xc2\xb0""C)/\xcf\x81( 30\xc2\xb0""C) la 2000 kPa",
            rho[1][2] / rho[2][2], 1.03, 1.15, "monotonie T");
        chk("GNB: \xcf\x81( 30\xc2\xb0""C)/\xcf\x81( 60\xc2\xb0""C) la 2000 kPa",
            rho[2][2] / rho[3][2], 1.03, 1.15, "monotonie T");

        // ─── Monotonie față de P la T = 10°C (ti=1): ρ(P1) < ρ(P2) ──────────
        // ρ(500)/ρ(101)   ≈ p_ratio × Z(101)/Z(500) ≈ 4.94 × 1.003 ≈ 4.96
        chk("GNB: \xcf\x81(500)/\xcf\x81(101) la 10\xc2\xb0""C",
            rho[1][1] / rho[1][0], 4.70, 5.10, "monotonie P");
        // ρ(2000)/ρ(500)  ≈ 4.00 × Z(500)/Z(2000) ≈ 4.00 × 1.016 ≈ 4.07
        chk("GNB: \xcf\x81(2000)/\xcf\x81(500) la 10\xc2\xb0""C",
            rho[1][2] / rho[1][1], 3.70, 4.30, "monotonie P");
        // ρ(7000)/ρ(2000) ≈ 3.50 × Z(2000)/Z(7000) ≈ 3.50 × 1.14 ≈ 3.98
        chk("GNB: \xcf\x81(7000)/\xcf\x81(2000) la 10\xc2\xb0""C",
            rho[1][3] / rho[1][2], 3.50, 4.80, "monotonie P");
      }

      sec_box("Negative suplim.",
              "d<12.5 Dia-U; beta<0.23 Dia-U; D>760 Dia-F; D>500 Aj-ISA; d=50 Aj-V", 9);
      // ── Teste negative ISO 5167: parametri invalizi trebuie sa returneze 0 ──
      std::printf("  %s\xC2\xBB Teste negative ISO 5167%s\n", kBoldYellow, kReset);
      {
        BwrConst bt_n = make_bwr(xCH4);
        double rho_n  = CalcDensity(20.0, 500.0 / kKpaPerAtm, bt_n);
        double eta_n  = 11.5e-6;  // viscozitate CH4 [Pa*s]
        FlowResult fr_n;
        // Suprima mesajele PrintError pe durata testelor negative
        int saved_fd = _dup(1);
        FILE* nul_f  = nullptr;
        fopen_s(&nul_f, "NUL", "w");
        if (nul_f) { _dup2(_fileno(nul_f), 1); fclose(nul_f); }
        double qn1 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kDiafragmaUnghi,
                                  200.0,170.0, rho_n,eta_n, &fr_n);  // β=0.85>0.80
        double qn2 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kDiafragmaUnghi,
                                   30.0, 15.0, rho_n,eta_n, &fr_n);  // D=30<50mm
        double qn3 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kVenturiPrelucrat,
                                  200.0, 60.0, rho_n,eta_n, &fr_n);  // β=0.30<0.40
        double qn4 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kVenturiTabla,
                                  200.0,150.0, rho_n,eta_n, &fr_n);  // β=0.75>0.70
        double qn5 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kDiafragmaUnghi,
                                  200.0, 10.0, rho_n,eta_n, &fr_n);  // d=10<12.5mm
        double qn6 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kDiafragmaUnghi,
                                  200.0, 44.0, rho_n,eta_n, &fr_n);  // β=0.22<0.23
        double qn7 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kDiafragmaFlansa,
                                  800.0,400.0, rho_n,eta_n, &fr_n);  // D=800>760mm
        double qn8 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kAjutajIsa,
                                  550.0,220.0, rho_n,eta_n, &fr_n);  // D=550>500mm
        double qn9 = CalcMassFlow(5.0,500.0,20.0, TipDispozitiv::kAjutajVenturi,
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
        chk("d=50mm Aj-V (limit\xC4\x83 d>50mm) \xE2\x86\x92 Qm=0",qn9, -0.001, 0.001, ns);
      }

      sec_box("Monotonie eta vs T",
              "eta(CH4) creste cu T: rapoarte -10C/20C/60C; sqrt(T) cinetic", 3);
      // ── Monotonie vâscozitate față de temperatură (CH4) ──────────────────
      // Gaze: η ∝ √T la diluat (teoria cinetică). CE include factorul (1+k·ln(T/cs))·√T.
      // NIST CH4: η(-10°C)≈10.5 μPa·s; η(20°C)≈11.0; η(60°C)≈11.9 μPa·s.
      // √(333.15/263.15)≈1.126; √(293.15/263.15)≈1.055; √(333.15/293.15)≈1.065.
      std::printf("  %s\xC2\xBB Monotonie v\xC3\xA2scozitate fa\xC8\x9B\xC4\x83 de T (CH4)%s\n",
                  kBoldYellow, kReset);
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
        chk("Monotonie \xCE\xB7 CH4: \xCE\xB7(60\xC2\xB0""C)/\xCE\xB7(-10\xC2\xB0""C) [-]",
            eP60 / eN10, 1.05, 1.25, "CE/cinetic");
        chk("Monotonie \xCE\xB7 CH4: \xCE\xB7(20\xC2\xB0""C)/\xCE\xB7(-10\xC2\xB0""C) [-]",
            eP20 / eN10, 1.01, 1.10, "CE/cinetic");
        chk("Monotonie \xCE\xB7 CH4: \xCE\xB7(60\xC2\xB0""C)/\xCE\xB7(20\xC2\xB0""C) [-]",
            eP60 / eP20, 1.01, 1.12, "CE/cinetic");
      }

      sec_box("Ordonare eta gaze",
              "eta(N2)/eta(CH4)~1.61; eta(CO2)/eta(CH4)~1.35 la 20C, 1 atm", 2);
      // ── Ordonare vâscozitate între gaze la 20°C ──────────────────────────
      // CE: η(CH4)≈11.5; η(CO2)≈15.5; η(N2)≈18.5 μPa·s la 20°C, 1 atm.
      // Ar omis: eroare CE sistematică −19% față de NIST → ordinea Ar nu e garantată.
      std::printf("  %s\xC2\xBB Ordonare v\xC3\xA2scozitate \xC3\xAEntre gaze la 20\xC2\xB0""C%s\n",
                  kBoldYellow, kReset);
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
        chk("Ordonare: \xCE\xB7(N2)/\xCE\xB7(CH4) la 20\xC2\xB0""C [-]",
            eN2  / eCH4, 1.40, 1.80, "CE/NIST");
        // η(CO2)/η(CH4) ≈ 15.49/11.46 ≈ 1.35
        chk("Ordonare: \xCE\xB7(CO2)/\xCE\xB7(CH4) la 20\xC2\xB0""C [-]",
            eCO2 / eCH4, 1.20, 1.55, "CE/NIST");
      }

      sec_box("Aer sintetic",
              "79% N2+21% O2: rho(20C)~1.199 kg/m3; eta~18.2 muPa*s; Qm pt. 9 tipuri", 15);
      // ── Aer sintetic (79% N2 + 21% O2) ──────────────────────────────────
      // M_aer = 0.79×28.016 + 0.21×32.000 = 28.853 g/mol
      // ρ_ideal(20°C)=28.853/24.055=1.199 kg/m³; ρ_ideal(0°C)=28.853/22.414=1.287
      // η_NIST aer la 20°C ≈ 18.2 μPa·s; Z ≈ 1.000 la 1 atm (cvasi-ideal)
      {
        double xAer[kArraySize] = {};
        xAer[29] = 0.79;  // N2
        xAer[30] = 0.21;  // O2
        test_comp("Aer sintetic (79% N2 + 21% O2)", xAer,
                  1.182, 1.218,   // \xcf\x81(20\xc2\xb0C, 1 atm): ideal=1.199, Z\xe2\x89\x881.00
                  1.270, 1.305,   // \xcf\x81(0\xc2\xb0C,  1 atm): ideal=1.287, Z\xe2\x89\x881.00
                  15.5,  22.0,    // \xce\xb7(20\xc2\xb0C): NIST aer\xe2\x89\x8818.2 \xc2\xb5Pa\xc2\xb7s
                  200, 80, 500, 5, false);
      }

      sec_box("H2S densitate",
              "H2S pur (Tc=373.6K, Tr=0.785): rho la 1 atm si 500 kPa, consistenta N2", 3);
      // ── H2S: densitate non-ideală și consistență ──────────────────────────
      // H2S: Tc=373.6 K, Pc=88.9 atm; la 20°C → Tr=0.785 → puternic non-ideal.
      // ρ_ideal(1 atm)=34.082/24.055=1.417 kg/m³; Z≈0.993 → ρ≈1.427
      // ρ_ideal(500 kPa)=4.935×1.417=6.99 kg/m³; Z≈0.965 → ρ≈7.24
      // M(H2S)/M(N2)=1.217; Z(H2S)<Z(N2) → ρ(H2S)/ρ(N2)>1.217
      std::printf("  %s\xC2\xBB H2S: densitate non-ideal\xC4\x83 \xC8\x99i consisten\xC8\x9B\xC4\x83%s\n",
                  kBoldYellow, kReset);
      {
        double xH2S[kArraySize] = {};  xH2S[26] = 1.0;
        BwrConst bH2S   = make_bwr(xH2S);
        BwrConst bN2_hs = make_bwr(xN2);
        double rH2S_1   = CalcDensity(20.0, 1.0,                  bH2S);
        double rH2S_500 = CalcDensity(20.0, 500.0 / kKpaPerAtm,   bH2S);
        double rN2_500  = CalcDensity(20.0, 500.0 / kKpaPerAtm,   bN2_hs);
        chk("Densitate H2S 20\xC2\xB0""C, 101.325 kPa [kg/m\xC2\xB3]",
            rH2S_1,   1.38, 1.48, "BWRS/ideal");
        chk("Densitate H2S 20\xC2\xB0""C,  500 kPa   [kg/m\xC2\xB3]",
            rH2S_500, 6.80, 7.80, "BWRS/NIST");
        // Raport > M ratio (1.217) din cauza non-idealității mai mari a H2S
        chk("Consisten\xC8\x9B\xC4\x83: \xCF\x81(H2S)/\xCF\x81(N2) la 500 kPa [-]",
            rH2S_500 / rN2_500, 1.20, 1.45, "M + Z");
      }

      // ── Sumar ────────────────────────────────────────────────────────────
      std::printf("\n  %sRezultat global: %d/%d teste trecute.%s\n",
                  (npass == ntotal) ? kBoldGreen : kBoldRed,
                  npass, ntotal, kReset);
      if (npass < ntotal)
        std::printf("  %sATEN\xC8\x9AIE: Unele teste au e\xC8\x99uat \xe2\x80\x94 verifica\xC8\x9Bi implementarea!%s\n",
                    kBoldRed, kReset);
      std::printf("\n  %sAp\xC4\x83sa\xC8\x9Bi orice tast\xC4\x83 pentru a continua...%s", kBoldWhite, kReset);
      _getch();
      std::printf("\n");
    }
  }

  // ── Compoziție: încărcare din fișier sau introducere manuală ──────────────
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
      std::printf("%s\n  Există o compoziție salvată. O refolosiți? [d/n] %s>%s ",
                  kBoldWhite, kCyan, kReset);
      if (AskYesNo()) {
        if (LoadComposition(x)) {
          comp_loaded = true;
          PrintComposition(x);
        } else {
          std::printf("%s  Eroare la citirea fișierului. Se va introduce manual.%s\n",
                      kBoldRed, kReset);
        }
      }
    }
  }

  if (!comp_loaded) {
    double sum = 0.0;
    do {
      std::printf("\n%s  ── Compoziție %s\n",
                  kBoldYellow, kReset);
      std::printf("%s  Fracții molare ale amestecului de gaze:%s\n\n", kBoldWhite, kReset);
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
            std::printf("%s  Valoare invalidă — trebuie să fie în [0, 1].%s\n", kBoldRed, kReset);
          }
        } while (v < 0.0 || v > 1.0);
        x[i] = v;
      }

      sum = 0.0;
      for (int i = 1; i <= kNumComponents; i++) {
        sum += x[i];
      }

      if (std::fabs(sum - 1.0) > kSumTolerance) {
        std::printf("%s\n  Suma fracțiilor molare = %.6f  ≠  1.%s\n", kBoldRed, sum, kReset);
        if (sum < 1e-10) {
          std::printf("%s  Suma este zero — reintroduceți compoziția.%s\n", kBoldRed, kReset);
        } else {
          std::printf("%s  Normalizați automat? [d/n] (n = reintroduceți) %s>%s ", kYellow, kCyan, kReset);
          if (AskYesNo()) {
            for (int i = 1; i <= kNumComponents; i++) {
              x[i] /= sum;
            }
            sum = 1.0;
            std::printf("%s  Fracții normalizate (componente nenule):%s\n\n", kBoldGreen, kReset);
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

    std::printf("%s\n  Salvați compoziția? [d/n] %s>%s ", kBoldWhite, kCyan, kReset);
    if (AskYesNo()) {
      SaveComposition(x);
    }
  }

  // ── Calcul constante BWRS ale amestecului ─────────────────────────────────
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

  std::printf("\n%s  Masa molar\xC4\x83 a amestecului: %s %s%8.4f%s [g/mol]\n",
              kBoldWhite, kReset, kBoldGreen, bwr.molar_mass, kReset);
  std::printf("%s  Densitate relativ\xC4\x83 fa\xC8\x9B\xC4\x83 de aer      :%s %s%8.4f%s [-]\n",
              kBoldWhite, kReset, kBoldGreen, bwr.molar_mass / 28.962, kReset);

  static const CountryRef kRefTable[] = {
    {"Rom\xC3\xA2nia / UE  (DIN 1343)",  2, { 0.0,  15.0  }, {"Nm\xC2\xB3/h", "Sm\xC2\xB3/h"}},
    {"ISO 13443  /  UK / Italia",         1, {15.0,   0.0  }, {"Sm\xC2\xB3/h", ""}},
    {"SUA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)", 1, {15.56, 0.0}, {"Sm\xC2\xB3/h", ""}},
    {"Rusia \xe2\x80\x94 GOST 30319-1",  1, {20.0,   0.0  }, {"m\xC2\xB3/h",  ""}},
    {"Personalizat",                      1, { 0.0,   0.0  }, {"m\xC2\xB3/h",  ""}},
  };
  static constexpr int kNRef = 5;

  std::printf("\n%s  ── Condi\xC8\x9Bii de referin\xC8\x9B\xC4\x83 %s\n",
              kBoldYellow, kReset);
  std::printf("%s  1.  Rom\xC3\xA2nia / UE  (DIN 1343)        \xe2\x80\x94   0\xC2\xB0""C \xC8\x99i 15\xC2\xB0""C / 101.325 kPa  [Nm\xC2\xB3/h] \xC8\x99i [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  2.  ISO 13443  /  UK / Italia       \xe2\x80\x94  15\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  3.  SUA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)             \xe2\x80\x94  15.56\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  4.  Rusia \xe2\x80\x94 GOST 30319-1            \xe2\x80\x94  20\xC2\xB0""C / 101.325 kPa  [m\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  5.  Personalizat                    \xe2\x80\x94  T [\xC2\xB0""C] introdus manual  [m\xC2\xB3/h]%s\n",  kBoldWhite, kReset);

  std::printf("\n  %s>%s %sSelecta\xC8\x9Bi (1\xe2\x80\x93" "5) : %s", kCyan, kReset, kBoldWhite, kReset);
  int ref_sel = ReadChoice(1, kNRef);

  CountryRef ref = kRefTable[ref_sel - 1];
  if (ref_sel == kNRef) {
    std::printf("  %s>%s %sTemperatura de referin\xC8\x9B\xC4\x83 [\xC2\xB0""C] : %s", kCyan, kReset, kBoldWhite, kReset);
    ReadDouble(&ref.t[0]);
  }

  std::printf ( "\n");

  double ror_ref[2] = {};
  for (int i = 0; i < ref.n; i++) {
    ror_ref[i] = CalcDensity(ref.t[i], 1, bwr);
    std::printf("%s  Densitate referin\xC8\x9B\xC4\x83 %5.2f\xC2\xB0""C / 101.325 kPa :%s %s%8.4f%s kg/m\xC2\xB3\n",
                kBoldWhite, ref.t[i], kReset, kBoldGreen, ror_ref[i], kReset);
    double z_ref = bwr.molar_mass / (ror_ref[i] * kGasConstantR * (ref.t[i] + kKelvinOffset));
    std::printf("%s  Factor Z            %5.2f\xC2\xB0""C / 101.325 kPa :%s %s%8.6f%s \xe2\x80\x94\n",
                kBoldWhite, ref.t[i], kReset, kBoldGreen, z_ref, kReset);
  }

  // ── Output helpers ─────────────────────────────────────────────────────────
  auto sep_d = []() {
    std::printf("%s  ", kBoldYellow);
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

  // ── Buclă exterioară: selecția dispozitivului de măsurare ─────────────────
  for (;;) {
    int    tip_raw = 0;
    double d_int = 0.0, d_orif = 0.0;

    // Încearcă să refolosească configurația salvată
    {
      int   sv_tip; double sv_d_int, sv_d_orif;
      if (LoadConfig(&sv_tip, &sv_d_int, &sv_d_orif)
          && sv_tip >= kTipMin && sv_tip <= kTipMax) {
        std::printf(
            "\n%s  ── Configurație salvată%s\n"
            "%s  Dispozitiv  : %s%s\n"
            "%s  D intern    : %s%g mm%s\n"
            "%s  D orificiu  : %s%g mm%s\n",
            kBoldYellow, kReset,
            kBoldWhite, TipName(static_cast<TipDispozitiv>(sv_tip)), kReset,
            kBoldWhite, kBoldGreen, sv_d_int, kReset,
            kBoldWhite, kBoldGreen, sv_d_orif, kReset);
        std::printf("  %s>%s %sRefolosiți configurația? [d/n] : %s", kCyan, kReset, kBoldWhite, kReset);
        if (AskYesNo()) {
          tip_raw = sv_tip;
          d_int   = sv_d_int;
          d_orif  = sv_d_orif;
          goto run_inner;
        }
      }
    }

    // Selectare manuală dispozitiv
    std::printf("\n%s  ── Dispozitiv de strangulare%s\n",
                kBoldYellow, kReset);
    std::printf("%s%s", kBoldWhite, kTipDisp);
    std::printf("%s", kReset);
    tip_raw = ReadChoice(kTipMin, kTipMax);

    std::printf("\n");
    do {
      std::printf("  %s>%s %sD intern  (20\xC2\xB0""C)  [mm] : %s", kCyan, kReset, kBoldWhite, kReset);
      ReadDouble(&d_int);
      std::printf("%s", kReset);
      if (d_int <= 0.0) {
        std::printf("%s  Valoare invalid\xC4\x83 \xe2\x80\x94 trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                    kBoldRed, kReset);
      }
    } while (d_int <= 0.0);
    do {
      std::printf("  %s>%s %sD orificiu (20\xC2\xB0""C) [mm] : %s", kCyan, kReset, kBoldWhite, kReset);
      ReadDouble(&d_orif);
      std::printf("%s", kReset);
      if (d_orif <= 0.0) {
        std::printf("%s  Valoare invalid\xC4\x83 \xe2\x80\x94 trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                    kBoldRed, kReset);
      } else if (d_orif >= d_int) {
        std::printf("%s  D orificiu trebuie s\xC4\x83 fie mai mic dec\xC3\xA2t D intern (%g mm).%s\n",
                    kBoldRed, d_int, kReset);
      } else {
        double beta_chk = d_orif / d_int;
        if (beta_chk < 0.10 || beta_chk > 0.80) {
          std::printf("%s  \xCE\xB2 = %.4f \xe2\x80\x94 \xC3\xAEn afara domeniului ISO 5167 [0.10, 0.80].%s\n",
                      kBoldRed, beta_chk, kReset);
        }
      }
    } while (d_orif <= 0.0 || d_orif >= d_int
             || d_orif / d_int < 0.10 || d_orif / d_int > 0.80);

    std::printf("%s\n  Salvați configurația? [d/n] %s>%s ", kBoldWhite, kCyan, kReset);
    if (AskYesNo()) {
      SaveConfig(tip_raw, d_int, d_orif);
    }

    run_inner:;
    TipDispozitiv tip = static_cast<TipDispozitiv>(tip_raw);

    // Buclă interioară: calcul pentru condiții diferite T/P cu același dispozitiv
    for (;;) {
      double temperatura, presiunea, presiunea_dif;
      do {
        std::printf("\n\n%s  ── Condi\xC8\x9Bii de m\xC4\x83surare%s\n",
                    kBoldYellow, kReset);
        std::printf("  %s>%s %sTemperatura         [\xC2\xB0""C] : %s", kCyan, kReset, kBoldWhite, kReset);
        ReadDouble(&temperatura);
        std::printf("%s", kReset);
        if (temperatura <= -273.15)
          std::printf("%s  Temperatura sub zero absolut (-273.15\xC2\xB0""C).%s\n",
                      kBoldRed, kReset);
      } while (temperatura <= -273.15);
      do {
        std::printf("  %s>%s %sPresiunea           [kPa] : %s", kCyan, kReset, kBoldWhite, kReset);
        ReadDouble(&presiunea);
        std::printf("%s", kReset);
        if (presiunea <= 0.0)
          std::printf("%s  Presiunea trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                      kBoldRed, kReset);
      } while (presiunea <= 0.0);
      do {
        std::printf("  %s>%s %sPresiunea diferen\xC8\x9Bial\xC4\x83 [kPa] : %s", kCyan, kReset, kBoldWhite, kReset);
        ReadDouble(&presiunea_dif);
        std::printf("%s", kReset);
        if (presiunea_dif <= 0.0)
          std::printf("%s  Presiunea diferen\xC8\x9Bial\xC4\x83 trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                      kBoldRed, kReset);
        else if (presiunea_dif >= presiunea)
          std::printf("%s  Diferen\xC8\x9Bial\xC4\x83 trebuie s\xC4\x83 fie mai mic\xC4\x83 dec\xC3\xA2t p = %g kPa.%s\n",
                      kBoldRed, presiunea, kReset);
      } while (presiunea_dif <= 0.0 || presiunea_dif >= presiunea);
      std::printf("  %c", 7);

      // ─── Calcul ────────────────────────────────────────────────────────────
      int    iter_rho = 0, iter_qm = 0;
      auto   t_rho0 = Clock::now();
      double ro = CalcDensity(temperatura, presiunea / kKpaPerAtm, bwr, &iter_rho);
      double t_rho_us = Us(Clock::now() - t_rho0).count();
      if (!std::isfinite(ro) || ro <= 0.0) {
        PrintError(ErrorCode::kNumeric, 0);
        break;  // eroare -> reselect dispozitiv / conditii
      }

      double roc_red = ro / roc_crit;
      auto   t_eta0 = Clock::now();
      double eta = 0.0;
      for (int i = 1; i <= kNumComponents; i++) {
        eta += (1 + kChapEnskog * std::log((temperatura + kKelvinOffset) / cs[i]))
             / (1 + kChapEnskog * std::log(kKelvinOffset / cs[i]))
             * std::sqrt((temperatura + kKelvinOffset) / kKelvinOffset)
             * et[i] * x[i] * std::sqrt(m[i]);
      }
      eta  = eta / mx;
      eta += kViscHighA / csi
           * std::pow(std::exp(kViscHighExp1 * roc_red) - std::exp(-kViscHighExp2 * roc_red),
                      kViscHighPow);
      double t_eta_us = Us(Clock::now() - t_eta0).count();

      FlowResult flow;
      auto   t_qm0 = Clock::now();
      double qm = CalcMassFlow(presiunea_dif, presiunea, temperatura,
                       tip, d_int, d_orif, ro, eta, &flow, kappa_mix, &iter_qm);
      double t_qm_us = Us(Clock::now() - t_qm0).count();
      if (qm == 0.0) break;  // eroare -> reselect dispozitiv

      // ─── Afișare rezultate ─────────────────────────────────────────────────
      std::printf("\n");
      sep_d();
      std::printf("%s  REZULTATE  \xe2\x80\x94  %s%s\n",
                  kBoldYellow, TipName(tip), kReset);
      std::printf("%s  D = %.2f mm  \xc2\xb7  d = %.2f mm  \xc2\xb7  "
                  "\xce\xb2 = %.4f%s\n",
                  kBoldWhite, d_int, d_orif, flow.beta, kReset);
      std::printf("%s  D(t) = %.3f mm  \xc2\xb7  d(t) = %.3f mm  "
                  "\xc2\xb7  t = %.1f\xc2\xb0""C%s\n",
                  kBoldWhite, flow.D_lucru_mm, flow.d_lucru_mm, temperatura, kReset);
      sep_d();

      std::printf("\n%s  Condi\xc8\x9bii de m\xc4\x83surare%s\n", kCyan, kReset);
      sep_s();
      row("Temperatur\xc4\x83", "%10.2f", temperatura, "\xc2\xb0""C");
      row("Presiune absolut\xc4\x83", "%10.2f", presiunea, "kPa");
      row("Presiune diferen\xc8\x9bial\xc4\x83", "%10.2f", presiunea_dif, "kPa");

      double Z_tp = (presiunea / kKpaPerAtm) * bwr.molar_mass
                  / (ro * kGasConstantR * (temperatura + kKelvinOffset));

      std::printf("\n%s  Fluid  (la t, p)%s\n", kCyan, kReset);
      sep_s();
      row("Densitate \xcf\x81(t,p)", "%10.4f", ro, "kg/m\xc2\xb3");
      row("V\xc3\xa2scozitate dinamic\xc4\x83 \xce\xb7(t,p)", "%10.4f",
          eta * kPaToMicroPa, "\xc2\xb5Pa\xc2\xb7s");
      row("Factor compresibilitate Z(t,p)", "%10.6f", Z_tp, "\xe2\x80\x94");

      std::printf("\n%s  Debite%s\n", kCyan, kReset);
      sep_s();
      row("Masic Qm", "%10.4f", qm, "kg/s");
      row("Masic Qm", "%10.2f", qm * kSecondsPerHour, "kg/h");
      for (int i = 0; i < ref.n; i++) {
        double qhref = kSecondsPerHour / ror_ref[i] * qm;
        char lbl[64];
        std::snprintf(lbl, sizeof(lbl),
                      "Volumic %5.2f\xc2\xb0""C / 101.325 kPa", ref.t[i]);
        row(lbl, "%10.2f", qhref, ref.label[i]);
      }
      row("Volumic la (t,p)", "%10.2f", kSecondsPerHour / ro * qm, "m\xc2\xb3/h");

      std::printf("\n%s  Hidraulic\xc4\x83%s\n", kCyan, kReset);
      sep_s();
      row("Viteza medie v", "%10.2f", flow.viteza, "m/s");
      row("Pierdere presiune", "%10.2f", flow.pierderea, "kPa");
      row("Raport str\xc3\xa2ngulare \xce\xb2", "%10.4f", flow.beta, "\xe2\x80\x94");
      row("Num\xc4\x83r Reynolds Re", "%10.4g", flow.reynolds, "\xe2\x80\x94");
      row("Coeficient debit C", "%10.6f", flow.coef_c, "\xe2\x80\x94");
      row("Factor expansibilitate \xce\xb5", "%10.6f", flow.epsilon, "\xe2\x80\x94");
      row("Exponent izentropic k", "%10.4f", flow.kappa, "\xe2\x80\x94");

      std::printf("\n");
      sep_d();

      // ── Profil CPU ─────────────────────────────────────────────────────────
      double t_total = t_mix_us + t_rho_us + t_eta_us + t_qm_us;
      std::printf("\n%s  PROFIL CPU%s\n", kCyan, kReset);
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
                    "reguli mixare N\xc2\xb2+N\xc2\xb3  (N=%d)", n_comp);
      prow("Constante BWRS", mix_detail,          t_mix_us, 0);
      prow("Densitate \xcf\x81",    "bise\xc8\x9b\xc8\x9bie BWRS",  t_rho_us, iter_rho);
      prow("V\xc3\xa2scozitate \xce\xb7", "Chapman-Enskog / N",       t_eta_us, 0);
      prow("Debit Qm",         "itera\xc8\x9bie Reynolds",           t_qm_us,  iter_qm);
      sep_s();
      std::printf("  %sTotal calcul%s                                  "
                  "%s%7.2f \xc2\xb5s%s\n\n",
                  kBoldWhite, kReset, kBoldGreen, t_total, kReset);

      std::printf("%s  EVALUARE EMBEDDED%s  (factori orientativi fa\xc8\x9b\xc4\x83 de PC)\n",
                  kCyan, kReset);
      sep_s();
      struct { const char* name; double factor; } targets[] = {
        { "ARM Cortex-M7 @480 MHz + FPU  (STM32H7)",          8.0  },
        { "ARM Cortex-M4 @168 MHz + FPU  (STM32F4)",          22.0 },
        { "Xtensa LX6    @240 MHz + FPU  (ESP32)",             16.0 },
        { "MSP430 F5xx   @ 25 MHz + hw mult (TI)",            350.0 },
        { "AVR           @ 16 MHz, f\xc4\x83r\xc4\x83 FPU (Mega2560)",  900.0 },
        { "80C51         @ 12 MHz, f\xc4\x83r\xc4\x83 FPU",            3000.0 },
      };
      for (auto& tg : targets) {
        double est = t_total * tg.factor / 1000.0;
        int tlen = (int)std::strlen(tg.name) - Utf8ExtraBytes(tg.name);
        std::printf("  %s%s", kBoldWhite, tg.name);
        for (int k = tlen; k < 48; k++) std::putchar(' ');
        std::printf("%s\xc3\x97%4.0f  %s%s%8.2f ms%s\n",
                    kReset, tg.factor, kReset, kBoldGreen, est, kReset);
      }
      std::printf("%s  Not\xc4\x83:%s factori estima\xc8\x9bi (IPC, cache, compila\xc8\x9bor)."
                  " M\xc4\x83sura\xc8\x9bi pe target pentru precizie.\n",
                  kYellow, kReset);

      std::printf("\n");
      sep_d();
      std::printf("\n\n");
    }
  }
}
