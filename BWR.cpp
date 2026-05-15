// =============================================================================
// BWR — Calcul debit gaze prin dispozitive de strangulare
// Revizia 2.0  |  05.2026
// ing. Agavriloaie Constantin  (original R 01.2004)
// ELCOST Impex
// =============================================================================
//
// DESCRIERE
//   Aplicația calculează densitatea, viscozitatea dinamică și debitele unui
//   amestec de gaze cu până la 36 de componente, pe baza:
//     • ecuației de stare Benedict-Webb-Rubin (BWR) pentru densitate;
//     • modelului Chapman-Enskog corectat cu termenul de densitate ridicată
//       pentru viscozitate;
//     • metodei ISO 5167-2:2003 (ecuația Reader-Harris/Gallagher)
//       pentru calculul debitului prin dispozitive de strangulare.
//
// INTRĂRI
//   1. Compoziția amestecului — fracții molare pentru fiecare din cei
//      36 de componenți (metan, etan, propan, ... acetilenă).
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
//   • Viscozitatea dinamică η(T, p)            [μPa·s]
//   • Debitul masic Qm                         [kg/s]
//   • Debitul volumic în condiții normale
//       (0 °C / 101,325 kPa)                   [Nm³/h]
//   • Debitul volumic în condiții standard
//       (15 °C / 101,325 kPa)                  [Sm³/h]
//   • Debitul volumic la T și p                [m³/h]
//   • Viteza medie a gazului în conductă       [m/s]
//   • Pierderea de presiune prin strangulare   [kPa]
//
// VALIDĂRI
//   Aplicația verifică limitele de aplicabilitate STAS 7347-90 pentru:
//   diametrul conductei, diametrul orificiului, raportul de strangulare β
//   și numărul Reynolds Re — și afișează mesaj de eroare la depășire.
//
// ALGORITM DENSITATE
//   Bisecție pe ecuația BWR până la convergența |p_calc − p| < 5×10⁻⁴ atm.
//   Constantele BWR ale amestecului se calculează o singură dată din
//   compoziție, folosind reguli de mixare pătratice și cubice.
//
// ALGORITM DEBIT
//   Iterație pe numărul Reynolds până la convergența |Qm_k − Qm_{k-1}| < 10⁻⁴.
//   Coeficientul de debit C și factorul de expansibilitate ε sunt recalculați
//   la fiecare iterație în funcție de Re și β.
//
// REFERINȚE
//   • ISO 5167-2:2003 — Orifice plates (ecuația Reader-Harris/Gallagher)
//   • ISO 5167-3:2003 — Nozzles and Venturi nozzles
//   • ISO 5167-4:2003 — Venturi tubes
//   • Benedict, Webb, Rubin (1940) — J. Chem. Phys. 8, 334
// =============================================================================

#include "BWR.h"
#ifdef _WIN32
#  include <conio.h>
#endif

static const char kCompFile[] = "bwr_comp.dat";
static const char kConfFile[] = "bwr_conf.dat";

static bool LoadComposition(double* x) {
  FILE* f = nullptr;
  fopen_s(&f, kCompFile, "r");
  if (!f) return false;
  for (int i = 1; i <= kNumComponents; i++) {
    if (fscanf_s(f, "%lf", &x[i]) != 1) { std::fclose(f); return false; }
  }
  std::fclose(f);
  return true;
}

static void SaveComposition(const double* x) {
  FILE* f = nullptr;
  fopen_s(&f, kCompFile, "w");
  if (!f) return;
  for (int i = 1; i <= kNumComponents; i++) std::fprintf(f, "%.8f\n", x[i]);
  std::fclose(f);
}

static bool LoadConfig(int* tip_raw, double* d_int, double* d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, kConfFile, "r");
  if (!f) return false;
  bool ok = (fscanf_s(f, "%d %lf %lf", tip_raw, d_int, d_orif) == 3);
  std::fclose(f);
  return ok;
}

static void SaveConfig(int tip_raw, double d_int, double d_orif) {
  FILE* f = nullptr;
  fopen_s(&f, kConfFile, "w");
  if (!f) return;
  std::fprintf(f, "%d\n%.4f\n%.4f\n", tip_raw, d_int, d_orif);
  std::fclose(f);
}

static void ExitApp() {
  std::printf("\033[2J\033[H");
  std::printf("%s\n  ELCOST Impex  —  BWR Gas Flow Calculator  v2.0/2026%s\n",
              kBoldYellow, kReset);
  std::exit(0);
}

static void ReadDouble(double* val) {
  char buf[64] = {};
  int pos = 0;
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 0 || ch == 0xE0) { (void)_getch(); continue; }
    if (ch == '\r') { std::printf("\n"); std::fflush(stdout); break; }
    if ((ch == 8 || ch == 127) && pos > 0) {
      pos--;
      std::printf("\b \b"); std::fflush(stdout);
      continue;
    }
    if (pos < 62 && (ch >= '0' && ch <= '9' || ch == '.' || (ch == '-' && pos == 0))) {
      buf[pos++] = static_cast<char>(ch);
      std::printf("%c", ch); std::fflush(stdout);
    }
  }
  *val = (pos > 0) ? std::atof(buf) : 0.0;
}

static void ReadInt(int* val) {
  char buf[32] = {};
  int pos = 0;
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 0 || ch == 0xE0) { (void)_getch(); continue; }
    if (ch == '\r') { std::printf("\n"); std::fflush(stdout); break; }
    if ((ch == 8 || ch == 127) && pos > 0) {
      pos--;
      std::printf("\b \b"); std::fflush(stdout);
      continue;
    }
    if (pos < 30 && ch >= '0' && ch <= '9') {
      buf[pos++] = static_cast<char>(ch);
      std::printf("%c", ch); std::fflush(stdout);
    }
  }
  *val = (pos > 0) ? std::atoi(buf) : 0;
}

static bool AskYesNo() {
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 'd' || ch == 'D') { std::printf("d\n"); std::fflush(stdout); return true; }
    if (ch == 'n' || ch == 'N') { std::printf("n\n"); std::fflush(stdout); return false; }
  }
}

static int Utf8ExtraBytes(const char* s) {
  int n = 0;
  while (*s) { if (((unsigned char)*s & 0xC0) == 0x80) n++; s++; }
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
    "Dioxid de carbon", "Aer", "Etilen\xC4\x83", "Propilen\xC4\x83", "Amoniac", "Acetilen\xC4\x83"
};

static void PrintComposition(const double* x) {
  std::printf("\n%s    Compoziția amestecului (fracții molare):%s\n", kBoldWhite, kReset);
  for (int i = 1; i <= kNumComponents; i++) {
    const char* name = kCompNames[i];
    int visualLen = (int)strlen(name) - Utf8ExtraBytes(name);
    std::printf("%s    %s", kBoldWhite, name);
    for (int j = visualLen; j < 25; j++) std::putchar(' ');
    std::printf(": %s%.8f%s\n", kBoldGreen, x[i], kReset);
  }
  std::printf("\n");
}

int main() {
#ifdef _WIN32
  std::system("chcp 65001 > nul");  // UTF-8 pentru caractere românești
#endif

  double x[kArraySize] = {};

  // Valori brute; tablourile marcate cu (*) sunt scalate după inițializare
  static const double A[kArraySize] = {
       0,         1.79894,   4.15556,   6.87225,  10.23264,  10.0847,   12.8,     12.7959,  12.1794,  11.842,
      16.43,     12.203,    12.203,    14.4373,   12.423,    12.423,    14.31,    14.31,    14.31,    17.5206,
      51.42,     51.86,      5.509772,  6.21,      0.040962,  1.34122,   2.12044,  0.040962,  0.823417, 1.053642,
       0.950852,  2.73742,   1.0292916, 3.33958,   6.1122,    3.7892819, 5.1079342};

  static double B[kArraySize] = {  // (*) /= 100
        0,        4.54625,   6.27724,  97.313,    13.7544,   12.4361,   16.0,     16.0053,  15.6751,  19.214,
       19.0,      8.1505,    8.1505,   17.7813,   20.246,    20.246,     9.1423,   9.1423,   9.1423,  19.9,
      110.3,    121.2,      50.30055,  40.8,       2.3661,    5.45425,   2.6182,   2.3661,   2.22826,  4.07526,
        2.22,     3.38943,   3.6216378, 5.56833,   8.50647,  51.646121,  6.946403};

  static double C[kArraySize] = {  // (*) *= 100000
       0,         0.318382,  1.79592,   5.08256,   8.49943,   9.9283,   17.5,     17.4632,  21.21219,  33.595,
      25.534,    22.125,    22.125,    33.1935,   51.237,    51.237,    31.564,   31.564,   31.564,   47.4574,
       1.032,     0.931,    34.2997,   29.0,       0.0000016227, 0.0856209, 7.9384, 0.0000016227, 0.1314125, 0.08059,
       0.326436,  1.38567,   0.11882507, 9.6936284, 1.3114,   1.785708,  6.506284};

  static const double a[kArraySize] = {
       0,         0.04352,   0.34516,   0.9477,    1.93763,   1.88231,   3.756,    3.7562,   4.0748,  10.108,
       4.6956,    7.4286,    7.4286,    7.11671,  11.786,    11.786,     7.5854,   7.5854,   7.5854,  10.36475,
      32.512,    31.423,     5.57,      4.32,      0.00057339, 0.3665,   0.84468,  0.00057339, 0.0288358, 0.025102,
       0.16269,   0.136814,  0.041402895, 0.259,   0.774056,  0.10354029, 0.6970948};

  static double b[kArraySize] = {  // (*) /= 100
       0,         0.252033,  1.1122,    2.25,      4.24352,   3.99983,   6.68,     6.6812,   6.6812,  14.0,
       7.9,      11.224,    11.224,    10.9131,   17.9131,   17.721,    14.321,   14.321,   14.321,  15.1954,
      53.14,     58.32,      7.663,     5.18,      0.000019727, 0.263158, 1.4653,  0.000019727, 0.215289, 0.23277,
       0.358835,  0.527236,  0.25625,   0.86,     18.7059,    0.071952516, 1.482999};

  static double c[kArraySize] = {  // (*) *= 100000
       0,         0.035878,  0.32767,   1.29,      2.8601,    3.164,     6.95,     6.95,     8.2417,  17.483,
      11.346,     9.5556,    9.5556,   15.1276,   22.586,    22.586,    13.252,   13.252,   13.252,  24.7,
      11.21,     18.23,     11.76418,  23.3,      0.0000000552, 0.0104,  1.1335,  0.0000000552, 0.007982437, 0.0072841,
       0.128274,  0.14918,   0.1729187,  2.8297636, 2.112,    0.0015753298, 1.0984375};

  static double alfa[kArraySize] = {  // (*) /= 1000
        0,        0.33,      0.243389,  0.607175,  1.07408,   1.10132,   1.7,      1.7,      1.81,    2.189,
        3.5948,   2.25,      2.25,      2.81086,   2.764,     2.764,     2.8155,   2.8155,   2.8155,  4.35611,
        2.207,    2.581,     0.7001,    0.318,     0.0072673, 0.135,     0.071955, 0.0072673, 0.035589, 0.1272,
      927.06,     0.0698611, 14.483933,  0.73924,   0.178,    0.0046521779, 0.27363248};

  static double gama[kArraySize] = {  // (*) /= 100
       0,         1.05,      1.18,      2.2,       3.4,       3.4,       4.63,     4.63,     4.75,    5.65,
       7.5,       6.289,     6.289,     6.668849,  6.799,     6.799,     7.446,    7.446,    7.446,   9.0,
       0.0318,    0.0209,    2.93,      1.12,      0.077942,  0.6,       0.59236,  0.077942, 0.233827, 0.53,
       3.1,       0.460593,  0.88722417, 2.911417,  0.923,   19.805156,  1.245167};

  static double V[kArraySize] = {  // (*) /= 1000
        0,        99.5,     148.00,    200.00,    263.00,    255.00,    303.00,   308.00,   311.00,   359.00,
      358.00,   367.00,   367.00,    368.00,    420.00,    420.04,    428.00,   418.00,   416.00,   426.00,
      482.00,   486.00,   260.00,    316.00,     65.0,      93.10,     95.00,    57.80,    75.20,    90.10,
       74.40,    94.0,    90.52,    124.00,    181.00,     72.50,    113.00};

  static const double m[kArraySize] = {
        0,       16.043,   30.070,    44.097,    58.124,    58.124,    72.151,   72.151,   72.151,   86.178,
       86.178,   86.178,   86.178,    86.178,   100.205,   100.205,   100.205,  100.205,  100.205,  100.205,
      114.232,  114.232,   78.114,    92.141,     2.016,    28.011,    34.082,    4.003,   39.944,   28.016,
       32.000,   44.011,   28.788,    28.054,    42.081,    17.032,   26.038};

  static const double cs[kArraySize] = {
        0,      148.6,    215.7,     237.1,     330.1,     331.4,     340.1,    340.1,    341.1,    398.2,
      398.3,    399.0,    399.1,     399.3,     412.1,     413.2,     413.2,    410.1,    412.2,    413.6,
      563.0,    564.0,    412.3,     418.3,      59.7,      91.7,     301.1,     10.22,    93.3,     71.4,
      106.7,    195.2,     78.6,     224.7,     298.9,     558.2,    231.8};

  static double et[kArraySize] = {  // (*) /= 10000
       0,        0.1085,   0.0915,    0.0805,    0.0735,    0.0725,    0.0711,   0.0696,   0.067,    0.0666,
       0.0658,   0.0647,   0.0651,    0.0641,    0.0626,    0.0626,    0.061,    0.0619,   0.0616,   0.0607,
       0.05947,  0.0577,   0.0745,    0.066,     0.0715,    0.1636,    0.141,    0.071,    0.174,    0.1755,
       0.2025,   0.1465,   0.1815,    0.094,     0.078,     0.093,     0.0943};

  static const double Tc[kArraySize] = {
        0,      190.7,    305.4,     369.9,     408.1,     425.2,     433.8,    460.4,    469.5,    488.7,
      499.9,    504.7,    496.5,     507.3,     520.3,     531.5,     530.3,    535.6,    540.8,    540.3,
      543.6,    568.6,    562.1,     592.0,      33.3,     133.0,     373.6,      5.3,    151.0,    126.2,
      154.8,    304.2,    132.5,     282.85,    364.55,    405.55,   308.85};

  static const double Pc[kArraySize] = {
        0,       45.8,     48.2,      42.0,      36.0,      37.5,      31.6,     32.9,     33.3,     30.7,
       30.9,     30.8,     30.0,      29.9,      27.4,      29.8,      27.2,     28.1,     28.6,     27.0,
       25.4,     24.6,     48.6,      41.6,      12.8,      34.5,      88.9,      2.26,    48.0,     33.5,
       50.1,     72.9,     37.17,     50.7,      45.4,     111.5,     61.6};

  static const double Zc[kArraySize] = {
       0,        0.29,     0.285,     0.277,     0.283,     0.274,     0.269,    0.268,    0.269,    0.273,
       0.27,     0.273,    0.27,      0.264,     0.27,      0.269,     0.267,    0.268,    0.267,    0.259,
       0.274,    0.256,    0.274,     0.271,     0.304,     0.294,     0.268,    0.3,      0.296,    0.291,
       0.292,    0.274,    0.291,     0.27,      0.274,     0.242,     0.274};

  // Scalare tablouri (*) — executată o singură dată la pornire
  for (int i = 0; i <= kNumComponents; i++) { B[i]    /= 100.0;    }
  for (int i = 0; i <= kNumComponents; i++) { C[i]    *= 100000.0; }
  for (int i = 0; i <= kNumComponents; i++) { b[i]    /= 100.0;    }
  for (int i = 0; i <= kNumComponents; i++) { c[i]    *= 100000.0; }
  for (int i = 0; i <= kNumComponents; i++) { alfa[i] /= 1000.0;   }
  for (int i = 0; i <= kNumComponents; i++) { gama[i] /= 100.0;    }
  for (int i = 0; i <= kNumComponents; i++) { V[i]    /= 1000.0;   }
  for (int i = 0; i <= kNumComponents; i++) { et[i]   /= 10000.0;  }

  std::printf("\033[2J\033[H");   // clear screen (ANSI)
  std::printf("%s", kBoldYellow);
  std::printf(
      "\n"
      "  ╔════════════════════════════════════════════════════════╗\n"
      "  ║                                                        ║\n"
      "  ║               E L C O S T   I m p e x                  ║\n"
      "  ║                                                        ║\n"
      "  ╠════════════════════════════════════════════════════════╣\n"
      "  ║                                                        ║\n"
      "  ║   BWR Gas Flow Calculator            v 2.0 / 2026      ║\n"
      "  ║                                                        ║\n"
      "  ║   ▸ Ecuație de stare   —  Benedict · Webb · Rubin      ║\n"
      "  ║   ▸ Viscozitate        —  Chapman-Enskog (corectat)    ║\n"
      "  ║   ▸ Debit gaze         —  ISO 5167-2:2003              ║\n"
      "  ║   ▸ Coef. debit C      —  Reader-Harris / Gallagher    ║\n"
      "  ║   ▸ 36 componenți      —  9 dispozitive de strangulare ║\n"
      "  ║   ▸ Ref. volumetrice   —  selectabile per țară         ║\n"
      "  ║                                                        ║\n"
      "  ╠════════════════════════════════════════════════════════╣\n"
      "  ║                    office@elcost.ro                    ║\n"
      "  ╚════════════════════════════════════════════════════════╝\n"
      "\n");
  std::printf("%s", kReset);

  // ── Compoziție: încărcare din fișier sau introducere manuală ──────────────
  bool comp_loaded = false;
  {
    int tmp_tip; double tmp_d1, tmp_d2;
    (void)tmp_tip; (void)tmp_d1; (void)tmp_d2;
    FILE* cf = nullptr;
    fopen_s(&cf, kCompFile, "r");
    if (cf) {
      std::fclose(cf);
      std::printf("%s\n  Există o compoziție salvată. O refolosiți? [d/n]: %s",
                  kBoldWhite, kReset);
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
    double sum;
    do {
      std::printf("\n%s  Compoziția în fracții molare a amestecului de gaze:%s\n",
                  kBoldWhite, kReset);
      for (int i = 1; i <= kNumComponents; i++) {
        const char* name = kCompNames[i];
        int vlen = (int)std::strlen(name) - Utf8ExtraBytes(name);
        double v;
        do {
          std::printf("%s  %s", kBoldWhite, name);
          for (int j = vlen; j < 22; j++) std::putchar(' ');
          std::printf(": ");
          ReadDouble(&v);
          std::printf("%s", kReset);
          if (v < 0.0 || v > 1.0)
            std::printf("%s  Valoare invalidă — trebuie să fie în [0, 1].%s\n", kBoldRed, kReset);
        } while (v < 0.0 || v > 1.0);
        x[i] = v;
      }

      sum = 0.0;
      for (int i = 1; i <= kNumComponents; i++) sum += x[i];

      if (std::fabs(sum - 1.0) > kSumTolerance) {
        std::printf("%s\n  Suma fracțiilor molare = %.6f  ≠  1.%s\n", kBoldRed, sum, kReset);
        if (sum < kSumTolerance) {
          std::printf("%s  Suma este zero — reintroduceți compoziția.%s\n", kBoldRed, kReset);
        } else {
          std::printf("%s  Normalizați automat? [d/n]  (n = reintroduceți): %s", kYellow, kReset);
          if (AskYesNo()) {
            for (int i = 1; i <= kNumComponents; i++) x[i] /= sum;
            sum = 1.0;
            std::printf("%s  Fracții normalizate (componente nenule):%s\n", kBoldGreen, kReset);
            for (int i = 1; i <= kNumComponents; i++) {
              if (x[i] > 0.0) {
                const char* name = kCompNames[i];
                int vlen = (int)std::strlen(name) - Utf8ExtraBytes(name);
                std::printf("%s    %s", kBoldGreen, name);
                for (int j = vlen; j < 22; j++) std::putchar(' ');
                std::printf(": %.8f%s\n", x[i], kReset);
              }
            }
          }
        }
      }
    } while (std::fabs(sum - 1.0) > kSumTolerance);

    std::printf("%s\n  Salvați compoziția? [d/n]: %s", kBoldWhite, kReset);
    if (AskYesNo()) SaveComposition(x);
  }

  // ── Calcul constante BWR ale amestecului ──────────────────────────────────
  BwrConst bwr;
  double tcam = 0.0, pcam = 0.0, zcam = 0.0, mx = 0.0;

  for (int i = 1; i <= kNumComponents; i++) {
    for (int j = 1; j <= kNumComponents; j++) {
      double kk = 1 - 8.0 * std::sqrt(V[i] * V[j])
               / std::pow(std::pow(V[i], 1.0/3) + std::pow(V[j], 1.0/3), 3);
      bwr.a0    += x[i] * x[j] * std::sqrt(A[i] * A[j]) * (1 - kk);
      bwr.b0    += x[i] * x[j] * std::sqrt(B[i] * B[j]);
      bwr.c0    += x[i] * x[j] * std::sqrt(C[i] * C[j]) * std::pow(1 - kk, 3);
      bwr.gamma += x[i] * x[j] * std::sqrt(gama[i] * gama[j]);
    }
  }

  for (int i = 1; i <= kNumComponents; i++) {
    for (int j = 1; j <= kNumComponents; j++) {
      for (int l = 1; l <= kNumComponents; l++) {
        double k1 = 1 - 8.0 * std::sqrt(V[i] * V[j])
                 / std::pow(std::pow(V[i], 1.0/3) + std::pow(V[j], 1.0/3), 3);
        double k2 = 1 - 8.0 * std::sqrt(V[i] * V[l])
                 / std::pow(std::pow(V[i], 1.0/3) + std::pow(V[l], 1.0/3), 3);
        double k3 = 1 - 8.0 * std::sqrt(V[j] * V[l])
                 / std::pow(std::pow(V[j], 1.0/3) + std::pow(V[l], 1.0/3), 3);
        bwr.a     += x[i] * x[j] * x[l]
                   * std::pow(a[i]*a[j]*a[l] * (1-k1)*(1-k2)*(1-k3), 1.0/3);
        bwr.b     += x[i] * x[j] * x[l] * std::pow(b[i]*b[j]*b[l], 1.0/3);
        bwr.c     += x[i] * x[j] * x[l]
                   * std::pow(c[i]*c[j]*c[l], 1.0/3) * (1-k1)*(1-k2)*(1-k3);
        bwr.alpha += x[i] * x[j] * x[l]
                   * std::pow(alfa[i]*alfa[j]*alfa[l], 1.0/3);
      }
    }
  }

  for (int i = 1; i <= kNumComponents; i++) {
    bwr.molar_mass += x[i] * m[i];
    tcam           += x[i] * Tc[i];
    pcam           += x[i] * Pc[i];
    zcam           += x[i] * Zc[i];
    mx             += x[i] * std::sqrt(m[i]);
  }

  double roc_crit = pcam / (kGasConstantR * zcam * tcam);
  double csi      = std::pow(tcam, 6)
                 / std::pow(bwr.molar_mass, 0.5)
                 / std::pow(pcam, 2.0/3.0);

  std::printf("  Masa molar\xC4\x83 a amestecului             : %s%8.4f%s [g/mol]\n",
              kBoldGreen, bwr.molar_mass, kReset);

  static const CountryRef kRefTable[] = {
    {"Rom\xC3\xA2nia / UE  (DIN 1343)",  2, { 0.0,  15.0  }, {"Nm\xC2\xB3/h", "Sm\xC2\xB3/h"}},
    {"ISO 13443  /  UK / Italia",         1, {15.0,   0.0  }, {"Sm\xC2\xB3/h", ""}},
    {"SUA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)", 1, {15.56, 0.0}, {"Sm\xC2\xB3/h", ""}},
    {"Rusia \xe2\x80\x94 GOST 30319-1",  1, {20.0,   0.0  }, {"m\xC2\xB3/h",  ""}},
    {"Personalizat",                      1, { 0.0,   0.0  }, {"m\xC2\xB3/h",  ""}},
  };
  static constexpr int kNRef = 5;

  std::printf("\n%s  Condi\xC8\x9Bii de referin\xC8\x9B\xC4\x83 volumetric\xC4\x83:%s\n\n",
              kBoldWhite, kReset);
  std::printf("%s  1.  Rom\xC3\xA2nia / UE  (DIN 1343)        \xe2\x80\x94   0\xC2\xB0""C \xC8\x99i 15\xC2\xB0""C / 101.325 kPa  [Nm\xC2\xB3/h] \xC8\x99i [Sm\xC2\xB3/h]%s\n", kBoldWhite, kReset);
  std::printf("%s  2.  ISO 13443  /  UK / Italia       \xe2\x80\x94  15\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n",                                                                    kBoldWhite, kReset);
  std::printf("%s  3.  SUA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)             \xe2\x80\x94  15.56\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n",                                           kBoldWhite, kReset);
  std::printf("%s  4.  Rusia \xe2\x80\x94 GOST 30319-1            \xe2\x80\x94  20\xC2\xB0""C / 101.325 kPa  [m\xC2\xB3/h]%s\n",                                                         kBoldWhite, kReset);
  std::printf("%s  5.  Personalizat                    \xe2\x80\x94  T [\xC2\xB0""C] introdus manual  [m\xC2\xB3/h]%s\n",                                                                 kBoldWhite, kReset);

  int ref_sel = 0;
  do {
    std::printf("%s\n  Selecta\xC8\x9Bi (1\xe2\x80\x93" "5)  > ", kBoldWhite);
    ReadInt(&ref_sel);
    std::printf("%s", kReset);
  } while (ref_sel < 1 || ref_sel > kNRef);

  CountryRef ref = kRefTable[ref_sel - 1];
  if (ref_sel == kNRef) {
    std::printf("%s  Temperatura de referin\xC8\x9B\xC4\x83 [\xC2\xB0""C]  : %s", kBoldWhite, kReset);
    ReadDouble(&ref.t[0]);
  }

  double ror_ref[2] = {};
  for (int i = 0; i < ref.n; i++) {
    ror_ref[i] = CalcDensity(ref.t[i], 1, bwr);
    std::printf("%s  Densitatea la %5.2f\xC2\xB0""C / 101.325 kPa     :%s %s%f%s %s[kg/m\xC2\xB3]%s\n",
                kBoldWhite, ref.t[i], kReset, kBoldGreen, ror_ref[i], kReset, kBoldWhite, kReset);
  }

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
            "%s\n  Configurație salvată:\n"
            "    Dispozitiv  : %s\n"
            "    D intern    : %g mm\n"
            "    D orificiu  : %g mm%s\n",
            kBoldWhite,
            TipName(static_cast<TipDispozitiv>(sv_tip)),
            sv_d_int, sv_d_orif, kReset);
        std::printf("%s  Refolosiți configurația? [d/n]: %s", kBoldWhite, kReset);
        if (AskYesNo()) {
          tip_raw = sv_tip;
          d_int   = sv_d_int;
          d_orif  = sv_d_orif;
          goto run_inner;
        }
      }
    }

    // Selectare manuală dispozitiv
    do {
      std::printf("\n%s%s", kBoldWhite, kTipDisp);
      ReadInt(&tip_raw);
      std::printf("%s", kReset);
      if (tip_raw < kTipMin || tip_raw > kTipMax)
        std::printf("%s\n  Dispozitiv de strangulare nedisponibil!%s\n", kBoldRed, kReset);
    } while (tip_raw < kTipMin || tip_raw > kTipMax);

    std::printf("\n");
    do {
      std::printf("%s  D intern (20\xC2\xB0""C)   [mm] : ", kBoldWhite);
      ReadDouble(&d_int);
      std::printf("%s", kReset);
      if (d_int <= 0.0)
        std::printf("%s  Valoare invalid\xC4\x83 \xe2\x80\x94 trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                    kBoldRed, kReset);
    } while (d_int <= 0.0);
    do {
      std::printf("%s  D orificiu (20\xC2\xB0""C) [mm] : ", kBoldWhite);
      ReadDouble(&d_orif);
      std::printf("%s", kReset);
      if (d_orif <= 0.0)
        std::printf("%s  Valoare invalid\xC4\x83 \xe2\x80\x94 trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                    kBoldRed, kReset);
      else if (d_orif >= d_int)
        std::printf("%s  D orificiu trebuie s\xC4\x83 fie mai mic dec\xC3\xA2t D intern (%g mm).%s\n",
                    kBoldRed, d_int, kReset);
      else {
        double beta_chk = d_orif / d_int;
        if (beta_chk < 0.10 || beta_chk > 0.80)
          std::printf("%s  \xCE\xB2 = %.4f \xe2\x80\x94 \xC3\xAEn afara domeniului ISO 5167 [0.10, 0.80].%s\n",
                      kBoldRed, beta_chk, kReset);
      }
    } while (d_orif <= 0.0 || d_orif >= d_int
             || d_orif / d_int < 0.10 || d_orif / d_int > 0.80);

    std::printf("%s\n  Salvați configurația? [d/n]: %s", kBoldWhite, kReset);
    if (AskYesNo()) SaveConfig(tip_raw, d_int, d_orif);

    run_inner:;
    TipDispozitiv tip = static_cast<TipDispozitiv>(tip_raw);

    // Buclă interioară: calcul pentru condiții diferite T/P cu același dispozitiv
    for (;;) {
      double temperatura, presiunea, presiunea_dif;
      do {
        std::printf("\n\n%s  Temperatura [\xC2\xB0""C]             : ", kBoldWhite);
        ReadDouble(&temperatura);
        std::printf("%s", kReset);
        if (temperatura <= -273.15)
          std::printf("%s  Temperatura sub zero absolut (-273.15\xC2\xB0""C).%s\n",
                      kBoldRed, kReset);
      } while (temperatura <= -273.15);
      do {
        std::printf("%s  Presiunea [kPa]              : ", kBoldWhite);
        ReadDouble(&presiunea);
        std::printf("%s", kReset);
        if (presiunea <= 0.0)
          std::printf("%s  Presiunea trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                      kBoldRed, kReset);
      } while (presiunea <= 0.0);
      do {
        std::printf("%s  Presiunea diferen\xC8\x9Bial\xC4\x83 [kPa] : ", kBoldWhite);
        ReadDouble(&presiunea_dif);
        std::printf("%s", kReset);
        if (presiunea_dif <= 0.0)
          std::printf("%s  Presiunea diferen\xC8\x9Bial\xC4\x83 trebuie s\xC4\x83 fie pozitiv\xC4\x83.%s\n",
                      kBoldRed, kReset);
        else if (presiunea_dif >= presiunea)
          std::printf("%s  Diferen\xC8\x9Bial\xC4\x83 trebuie s\xC4\x83 fie mai mic\xC4\x83 dec\xC3\xA2t p = %g kPa.%s\n",
                      kBoldRed, presiunea, kReset);
      } while (presiunea_dif <= 0.0 || presiunea_dif >= presiunea);
      std::printf("  %c\n", 7);

      double ro = CalcDensity(temperatura, presiunea / kKpaPerAtm, bwr);
      std::printf("  Densitatea (t,p)                       : %s%8.4f%s [kg/m\xC2\xB3]\n",
                  kBoldGreen, ro, kReset);

      double roc_red = ro / roc_crit;
      double eta = 0.0;
      for (int i = 1; i <= kNumComponents; i++) {
        eta += (1 + 0.323 * std::log((temperatura + kKelvinOffset) / cs[i]))
             / (1 + 0.323 * std::log(kKelvinOffset / cs[i]))
             * std::sqrt((temperatura + kKelvinOffset) / kKelvinOffset)
             * et[i] * x[i] * std::sqrt(m[i]);
      }
      eta  = eta / mx;
      eta += 10.8e-8 / csi
           * std::pow(std::exp(1.439 * roc_red) - std::exp(-1.111 * roc_red),
                      1.358);
      std::printf("  Viscozitatea dinamic\xC4\x83 (t,p)            : %s%8.4f%s [\xC2\xB5Pa\xC2\xB7s]\n",
                  kBoldGreen, eta * 1000000, kReset);

      FlowResult flow;
      double qm = CalcMassFlow(presiunea_dif, presiunea, temperatura,
                       tip, d_int, d_orif, ro, eta, &flow);
      if (qm == 0.0) break;  // eroare -> reselect dispozitiv

      std::printf("  Debitul masic (t,p)                    : %s%7.4f%s [kg/s]\n",
                  kBoldGreen, qm, kReset);
      for (int i = 0; i < ref.n; i++) {
        double qhref = 3600.0 / ror_ref[i] * qm;
        std::printf("  Debitul vol.  (%5.2f\xC2\xB0""C / 101.325 kPa)  : %s%7.2f%s [%s]\n",
                    ref.t[i], kBoldGreen, qhref, kReset, ref.label[i]);
      }
      std::printf("  Debitul vol.  (t,p)                    : %s%7.2f%s [m³/h]\n",
                  kBoldGreen, 3600 / ro * qm, kReset);
      std::printf("  Viteza medie a gazului                 : %s%.2f%s [m/s]\n",
                  kBoldGreen, flow.viteza, kReset);
      std::printf("  Pierderea de presiune prin strangulare : %s%.2f%s [kPa]\n",
                  kBoldGreen, flow.pierderea, kReset);
      std::printf("  Raportul de str\xC3\xA2ngulare \xCE\xB2              : %s%.4f%s [-]\n",
                  kBoldGreen, flow.beta, kReset);
      std::printf("  Num\xC4\x83rul Reynolds                       : %s%.4g%s [-]\n",
                  kBoldGreen, flow.reynolds, kReset);
    }
  }
}
