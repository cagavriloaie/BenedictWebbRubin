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

static bool LoadComposition(float* x) {
  FILE* f = std::fopen(kCompFile, "r");
  if (!f) return false;
  for (int i = 1; i <= kNumComponents; i++) {
    if (std::fscanf(f, "%f", &x[i]) != 1) { std::fclose(f); return false; }
  }
  std::fclose(f);
  return true;
}

static void SaveComposition(const float* x) {
  FILE* f = std::fopen(kCompFile, "w");
  if (!f) return;
  for (int i = 1; i <= kNumComponents; i++) std::fprintf(f, "%.8f\n", x[i]);
  std::fclose(f);
}

static bool LoadConfig(int* tip_raw, float* d_int, float* d_orif) {
  FILE* f = std::fopen(kConfFile, "r");
  if (!f) return false;
  bool ok = (std::fscanf(f, "%d %f %f", tip_raw, d_int, d_orif) == 3);
  std::fclose(f);
  return ok;
}

static void SaveConfig(int tip_raw, float d_int, float d_orif) {
  FILE* f = std::fopen(kConfFile, "w");
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

static bool ReadFloat(float* val) {
  char buf[64] = {};
  int pos = 0;
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 0 || ch == 0xE0) { _getch(); continue; }
    if (ch == '\r') { std::printf("\n"); std::fflush(stdout); break; }
    if ((ch == 8 || ch == 127) && pos > 0) {
      pos--;
      std::printf("\b \b"); std::fflush(stdout);
      continue;
    }
    if (pos < 62 && (ch >= '0' && ch <= '9' || ch == '.' || ch == '-')) {
      buf[pos++] = static_cast<char>(ch);
      std::printf("%c", ch); std::fflush(stdout);
    }
  }
  *val = (pos > 0) ? static_cast<float>(std::atof(buf)) : 0.0f;
  return true;
}

static bool ReadInt(int* val) {
  char buf[32] = {};
  int pos = 0;
  for (;;) {
    int ch = _getch();
    if (ch == 27) ExitApp();
    if (ch == 0 || ch == 0xE0) { _getch(); continue; }
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
  return true;
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

static void PrintComposition(const float* x) {
  static const char* const kNames[] = {
      nullptr,
      "Metan", "Etan", "Propan", "Izobutan", "N-butan",
      "Neopentan", "Izopentan", "N-pentan", "2,2-dimetilbutan", "2,3-dimetilbutan",
      "3-metilpentan", "2-metilpentan", "N-hexan", "2,4-dimetilpentan", "2,2,3-trimetilbutan",
      "2-metilhexan", "3-metilhexan", "3-etilpentan", "N-heptan", "2,2,3-trimetilpentan",
      "N-octan", "Benzen", "Toluen", "Hidrogen", "Monoxid de carbon",
      "Hidrogen sulfurat", "Heliu", "Argon", "Azot", "Oxigen",
      "Dioxid de carbon", "Aer", "Etilen\xC4\x83", "Propilen\xC4\x83", "Amoniac", "Acetilen\xC4\x83"
  };
  std::printf("\n%s    Compoziția amestecului (fracții molare):%s\n", kBoldCyan, kReset);
  for (int i = 1; i <= kNumComponents; i++) {
    const char* name = kNames[i];
    int visualLen = (int)strlen(name) - Utf8ExtraBytes(name);
    std::printf("%s    %s", kCyan, name);
    for (int j = visualLen; j < 25; j++) std::putchar(' ');
    std::printf(": %s%.8f%s\n", kBoldGreen, x[i], kReset);
  }
  std::printf("\n");
}

int main() {
#ifdef _WIN32
  std::system("chcp 65001 > nul");  // UTF-8 pentru caractere românești
#endif

  float x[kArraySize] = {};

  // Valori brute; tablourile marcate cu (*) sunt scalate după inițializare
  static const float A[kArraySize] = {
       0,          1.79894f,   4.15556f,   6.87225f,  10.23264f,  10.0847f,   12.8f,     12.7959f,  12.1794f,  11.842f,
      16.43f,     12.203f,    12.203f,    14.4373f,   12.423f,    12.423f,    14.31f,    14.31f,    14.31f,    17.5206f,
      51.42f,     51.86f,      5.509772f,  6.21f,      0.040962f,  1.34122f,   2.12044f,  0.040962f,  0.823417f, 1.053642f,
       0.950852f,  2.73742f,   1.0292916f, 3.33958f,   6.1122f,    3.7892819f, 5.1079342f};

  static float B[kArraySize] = {  // (*) /= 100
        0,         4.54625f,   6.27724f,  97.313f,    13.7544f,   12.4361f,   16.0f,     16.0053f,  15.6751f,  19.214f,
       19.0f,      8.1505f,    8.1505f,   17.7813f,   20.246f,    20.246f,     9.1423f,   9.1423f,   9.1423f,  19.9f,
      110.3f,    121.2f,      50.30055f,  40.8f,       2.3661f,    5.45425f,   2.6182f,   2.3661f,   2.22826f,  4.07526f,
        2.22f,     3.38943f,   3.6216378f, 5.56833f,   8.50647f,  51.646121f,  6.946403f};

  static float C[kArraySize] = {  // (*) *= 100000
       0,          0.318382f,  1.79592f,   5.08256f,   8.49943f,   9.9283f,   17.5f,     17.4632f,  21.21219f,  33.595f,
      25.534f,    22.125f,    22.125f,    33.1935f,   51.237f,    51.237f,    31.564f,   31.564f,   31.564f,   47.4574f,
       1.032f,     0.931f,    34.2997f,   29.0f,       0.0000016227f, 0.0856209f, 7.9384f, 0.0000016227f, 0.1314125f, 0.08059f,
       0.326436f,  1.38567f,   0.11882507f, 9.6936284f, 1.3114f,   1.785708f,  6.506284f};

  static const float a[kArraySize] = {
       0,          0.04352f,   0.34516f,   0.9477f,    1.93763f,   1.88231f,   3.756f,    3.7562f,   4.0748f,  10.108f,
       4.6956f,    7.4286f,    7.4286f,    7.11671f,  11.786f,    11.786f,     7.5854f,   7.5854f,   7.5854f,  10.36475f,
      32.512f,    31.423f,     5.57f,      4.32f,      0.00057339f, 0.3665f,   0.84468f,  0.00057339f, 0.0288358f, 0.025102f,
       0.16269f,   0.136814f,  0.041402895f, 0.259f,   0.774056f,  0.10354029f, 0.6970948f};

  static float b[kArraySize] = {  // (*) /= 100
       0,          0.252033f,  1.1122f,    2.25f,      4.24352f,   3.99983f,   6.68f,     6.6812f,   6.6812f,  14.0f,
       7.9f,      11.224f,    11.224f,    10.9131f,   17.9131f,   17.721f,    14.321f,   14.321f,   14.321f,  15.1954f,
      53.14f,     58.32f,      7.663f,     5.18f,      0.000019727f, 0.263158f, 1.4653f,  0.000019727f, 0.215289f, 0.23277f,
       0.358835f,  0.527236f,  0.25625f,   0.86f,     18.7059f,    0.071952516f, 1.482999f};

  static float c[kArraySize] = {  // (*) *= 100000
       0,          0.035878f,  0.32767f,   1.29f,      2.8601f,    3.164f,     6.95f,     6.95f,     8.2417f,  17.483f,
      11.346f,     9.5556f,    9.5556f,   15.1276f,   22.586f,    22.586f,    13.252f,   13.252f,   13.252f,  24.7f,
      11.21f,     18.23f,     11.76418f,  23.3f,      0.0000000552f, 0.0104f,  1.1335f,  0.0000000552f, 0.007982437f, 0.0072841f,
       0.128274f,  0.14918f,   0.1729187f,  2.8297636f, 2.112f,    0.0015753298f, 1.0984375f};

  static float alfa[kArraySize] = {  // (*) /= 1000
        0,         0.33f,      0.243389f,  0.607175f,  1.07408f,   1.10132f,   1.7f,      1.7f,      1.81f,    2.189f,
        3.5948f,   2.25f,      2.25f,      2.81086f,   2.764f,     2.764f,     2.8155f,   2.8155f,   2.8155f,  4.35611f,
        2.207f,    2.581f,     0.7001f,    0.318f,     0.0072673f, 0.135f,     0.071955f, 0.0072673f, 0.035589f, 0.1272f,
      927.06f,     0.0698611f, 14.483933f,  0.73924f,   0.178f,    0.0046521779f, 0.27363248f};

  static float gama[kArraySize] = {  // (*) /= 100
       0,          1.05f,      1.18f,      2.2f,       3.4f,       3.4f,       4.63f,     4.63f,     4.75f,    5.65f,
       7.5f,       6.289f,     6.289f,     6.668849f,  6.799f,     6.799f,     7.446f,    7.446f,    7.446f,   9.0f,
       0.0318f,    0.0209f,    2.93f,      1.12f,      0.077942f,  0.6f,       0.59236f,  0.077942f, 0.233827f, 0.53f,
       3.1f,       0.460593f,  0.88722417f, 2.911417f,  0.923f,   19.805156f,  1.245167f};

  static float V[kArraySize] = {  // (*) /= 1000
        0,        99.5f,     148.00f,    200.00f,    263.00f,    255.00f,    303.00f,   308.00f,   311.00f,   359.00f,
      358.00f,   367.00f,   367.00f,    368.00f,    420.00f,    420.04f,    428.00f,   418.00f,   416.00f,   426.00f,
      482.00f,   486.00f,   260.00f,    316.00f,     65.0f,      93.10f,     95.00f,    57.80f,    75.20f,    90.10f,
       74.40f,    94.0f,    90.52f,    124.00f,    181.00f,     72.50f,    113.00f};

  static const float m[kArraySize] = {
        0,        16.043f,   30.070f,    44.097f,    58.124f,    58.124f,    72.151f,   72.151f,   72.151f,   86.178f,
       86.178f,   86.178f,   86.178f,    86.178f,   100.205f,   100.205f,   100.205f,  100.205f,  100.205f,  100.205f,
      114.232f,  114.232f,   78.114f,    92.141f,     2.016f,    28.011f,    34.082f,    4.003f,   39.944f,   28.016f,
       32.000f,   44.011f,   28.788f,    28.054f,    42.081f,    17.032f,   26.038f};

  static const float cs[kArraySize] = {
        0,       148.6f,    215.7f,     237.1f,     330.1f,     331.4f,     340.1f,    340.1f,    341.1f,    398.2f,
      398.3f,    399.0f,    399.1f,     399.3f,     412.1f,     413.2f,     413.2f,    410.1f,    412.2f,    413.6f,
      563.0f,    564.0f,    412.3f,     418.3f,      59.7f,      91.7f,     301.1f,     10.22f,    93.3f,     71.4f,
      106.7f,    195.2f,     78.6f,     224.7f,     298.9f,     558.2f,    231.8f};

  static float et[kArraySize] = {  // (*) /= 10000
       0,         0.1085f,   0.0915f,    0.0805f,    0.0735f,    0.0725f,    0.0711f,   0.0696f,   0.067f,    0.0666f,
       0.0658f,   0.0647f,   0.0651f,    0.0641f,    0.0626f,    0.0626f,    0.061f,    0.0619f,   0.0616f,   0.0607f,
       0.05947f,  0.0577f,   0.0745f,    0.066f,     0.0715f,    0.1636f,    0.141f,    0.071f,    0.174f,    0.1755f,
       0.2025f,   0.1465f,   0.1815f,    0.094f,     0.078f,     0.093f,     0.0943f};

  static const float Tc[kArraySize] = {
        0,       190.7f,    305.4f,     369.9f,     408.1f,     425.2f,     433.8f,    460.4f,    469.5f,    488.7f,
      499.9f,    504.7f,    496.5f,     507.3f,     520.3f,     531.5f,     530.3f,    535.6f,    540.8f,    540.3f,
      543.6f,    568.6f,    562.1f,     592.0f,      33.3f,     133.0f,     373.6f,      5.3f,    151.0f,    126.2f,
      154.8f,    304.2f,    132.5f,     282.85f,    364.55f,    405.55f,   308.85f};

  static const float Pc[kArraySize] = {
        0,        45.8f,     48.2f,      42.0f,      36.0f,      37.5f,      31.6f,     32.9f,     33.3f,     30.7f,
       30.9f,     30.8f,     30.0f,      29.9f,      27.4f,      29.8f,      27.2f,     28.1f,     28.6f,     27.0f,
       25.4f,     24.6f,     48.6f,      41.6f,      12.8f,      34.5f,      88.9f,      2.26f,    48.0f,     33.5f,
       50.1f,     72.9f,     37.17f,     50.7f,      45.4f,     111.5f,     61.6f};

  static const float Zc[kArraySize] = {
       0,         0.29f,     0.285f,     0.277f,     0.283f,     0.274f,     0.269f,    0.268f,    0.269f,    0.273f,
       0.27f,     0.273f,    0.27f,      0.264f,     0.27f,      0.269f,     0.267f,    0.268f,    0.267f,    0.259f,
       0.274f,    0.256f,    0.274f,     0.271f,     0.304f,     0.294f,     0.268f,    0.3f,      0.296f,    0.291f,
       0.292f,    0.274f,    0.291f,     0.27f,      0.274f,     0.242f,     0.274f};

  // Scalare tablouri (*) — executată o singură dată la pornire
  for (int i = 0; i <= kNumComponents; i++) { B[i]    /= 100.0f;    }
  for (int i = 0; i <= kNumComponents; i++) { C[i]    *= 100000.0f; }
  for (int i = 0; i <= kNumComponents; i++) { b[i]    /= 100.0f;    }
  for (int i = 0; i <= kNumComponents; i++) { c[i]    *= 100000.0f; }
  for (int i = 0; i <= kNumComponents; i++) { alfa[i] /= 1000.0f;   }
  for (int i = 0; i <= kNumComponents; i++) { gama[i] /= 100.0f;    }
  for (int i = 0; i <= kNumComponents; i++) { V[i]    /= 1000.0f;   }
  for (int i = 0; i <= kNumComponents; i++) { et[i]   /= 10000.0f;  }

  std::printf("\033[2J\033[H");   // clear screen (ANSI)
  std::printf("%s", kBoldYellow);
  std::printf(
      "\n"
      "  ╔════════════════════════════════════════════════════════╗\n"
      "  ║                                                        ║\n"
      "  ║           E L C O S T   I m p e x                      ║\n"
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
      "  ║   ▸ Ref. volumetrice   —  selectabile per \xC8\x9B" "ar\xC4\x83         \xe2\x95\x91\n"
      "  ║                                                        ║\n"
      "  ╚════════════════════════════════════════════════════════╝\n"
      "\n");
  std::printf("%s", kReset);

  // ── Compoziție: încărcare din fișier sau introducere manuală ──────────────
  bool comp_loaded = false;
  {
    int tmp_tip; float tmp_d1, tmp_d2;
    (void)tmp_tip; (void)tmp_d1; (void)tmp_d2;
    FILE* cf = std::fopen(kCompFile, "r");
    if (cf) {
      std::fclose(cf);
      std::printf("%s\n  Există o compoziție salvată. O refolosiți? [d/n]: %s",
                  kCyan, kReset);
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
    float sum;
    do {
      std::printf("\n%s Compoziția în fracții molare a amestecului de gaze:%s\n",
                  kBoldCyan, kReset);
      std::printf("%s  Metan                 : %s", kCyan, kReset); ReadFloat(&x[1]);
      std::printf("%s  Etan                  : %s", kCyan, kReset); ReadFloat(&x[2]);
      std::printf("%s  Propan                : %s", kCyan, kReset); ReadFloat(&x[3]);
      std::printf("%s  Izobutan              : %s", kCyan, kReset); ReadFloat(&x[4]);
      std::printf("%s  N-butan               : %s", kCyan, kReset); ReadFloat(&x[5]);
      std::printf("%s  Neopentan             : %s", kCyan, kReset); ReadFloat(&x[6]);
      std::printf("%s  Izopentan             : %s", kCyan, kReset); ReadFloat(&x[7]);
      std::printf("%s  N-pentan              : %s", kCyan, kReset); ReadFloat(&x[8]);
      std::printf("%s  2,2-dimetilbutan      : %s", kCyan, kReset); ReadFloat(&x[9]);
      std::printf("%s  2,3-dimetilbutan      : %s", kCyan, kReset); ReadFloat(&x[10]);
      std::printf("%s  3-metilpentan         : %s", kCyan, kReset); ReadFloat(&x[11]);
      std::printf("%s  2-metilpentan         : %s", kCyan, kReset); ReadFloat(&x[12]);
      std::printf("%s  N-hexan               : %s", kCyan, kReset); ReadFloat(&x[13]);
      std::printf("%s  2,4-dimetilpentan     : %s", kCyan, kReset); ReadFloat(&x[14]);
      std::printf("%s  2,2,3-trimetilbutan   : %s", kCyan, kReset); ReadFloat(&x[15]);
      std::printf("%s  2-metilhexan          : %s", kCyan, kReset); ReadFloat(&x[16]);
      std::printf("%s  3-metilhexan          : %s", kCyan, kReset); ReadFloat(&x[17]);
      std::printf("%s  3-etilpentan          : %s", kCyan, kReset); ReadFloat(&x[18]);
      std::printf("%s  N-heptan              : %s", kCyan, kReset); ReadFloat(&x[19]);
      std::printf("%s  2,2,3-trimetilpentan  : %s", kCyan, kReset); ReadFloat(&x[20]);
      std::printf("%s  N-octan               : %s", kCyan, kReset); ReadFloat(&x[21]);
      std::printf("%s  Benzen                : %s", kCyan, kReset); ReadFloat(&x[22]);
      std::printf("%s  Toluen                : %s", kCyan, kReset); ReadFloat(&x[23]);
      std::printf("%s  Hidrogen              : %s", kCyan, kReset); ReadFloat(&x[24]);
      std::printf("%s  Monoxid de carbon     : %s", kCyan, kReset); ReadFloat(&x[25]);
      std::printf("%s  Hidrogen sulfurat     : %s", kCyan, kReset); ReadFloat(&x[26]);
      std::printf("%s  Heliu                 : %s", kCyan, kReset); ReadFloat(&x[27]);
      std::printf("%s  Argon                 : %s", kCyan, kReset); ReadFloat(&x[28]);
      std::printf("%s  Azot                  : %s", kCyan, kReset); ReadFloat(&x[29]);
      std::printf("%s  Oxigen                : %s", kCyan, kReset); ReadFloat(&x[30]);
      std::printf("%s  Dioxid de carbon      : %s", kCyan, kReset); ReadFloat(&x[31]);
      std::printf("%s  Aer                   : %s", kCyan, kReset); ReadFloat(&x[32]);
      std::printf("%s  Etilenă               : %s", kCyan, kReset); ReadFloat(&x[33]);
      std::printf("%s  Propilenă             : %s", kCyan, kReset); ReadFloat(&x[34]);
      std::printf("%s  Amoniac               : %s", kCyan, kReset); ReadFloat(&x[35]);
      std::printf("%s  Acetilenă             : %s", kCyan, kReset); ReadFloat(&x[36]);

      sum = 0.0f;
      for (int i = 1; i <= kNumComponents; i++) sum += x[i];

      if (std::fabs(sum - 1.0f) > kSumTolerance) {
        std::printf("%s\n  Suma fracțiilor molare = %.6f  ≠  1.%s\n", kBoldRed, sum, kReset);
        std::printf("%s  Normalizați automat? [d/n]  (n = reintroduceți): %s", kYellow, kReset);
        if (AskYesNo()) {
          for (int i = 1; i <= kNumComponents; i++) x[i] /= sum;
          sum = 1.0f;
          std::printf("%s  Fracțiile au fost normalizate.%s\n", kBoldGreen, kReset);
        }
      }
    } while (std::fabs(sum - 1.0f) > kSumTolerance);

    std::printf("%s\n  Salvați compoziția? [d/n]: %s", kCyan, kReset);
    if (AskYesNo()) SaveComposition(x);
  }

  // ── Calcul constante BWR ale amestecului ──────────────────────────────────
  BwrConst bwr;
  float tcam = 0.0f, pcam = 0.0f, zcam = 0.0f, mx = 0.0f;

  for (int i = 1; i <= kNumComponents; i++) {
    for (int j = 1; j <= kNumComponents; j++) {
      float kk = 1 - 8.0f * std::sqrt(V[i] * V[j])
               / std::pow(std::pow(V[i], 1.0f/3) + std::pow(V[j], 1.0f/3), 3);
      bwr.a0    += x[i] * x[j] * std::sqrt(A[i] * A[j]) * (1 - kk);
      bwr.b0    += x[i] * x[j] * std::sqrt(B[i] * B[j]);
      bwr.c0    += x[i] * x[j] * std::sqrt(C[i] * C[j]) * std::pow(1 - kk, 3);
      bwr.gamma += x[i] * x[j] * std::sqrt(gama[i] * gama[j]);
    }
  }

  for (int i = 1; i <= kNumComponents; i++) {
    for (int j = 1; j <= kNumComponents; j++) {
      for (int l = 1; l <= kNumComponents; l++) {
        float k1 = 1 - 8.0f * std::sqrt(V[i] * V[j])
                 / std::pow(std::pow(V[i], 1.0f/3) + std::pow(V[j], 1.0f/3), 3);
        float k2 = 1 - 8.0f * std::sqrt(V[i] * V[l])
                 / std::pow(std::pow(V[i], 1.0f/3) + std::pow(V[l], 1.0f/3), 3);
        float k3 = 1 - 8.0f * std::sqrt(V[j] * V[l])
                 / std::pow(std::pow(V[j], 1.0f/3) + std::pow(V[l], 1.0f/3), 3);
        bwr.a     += x[i] * x[j] * x[l]
                   * std::pow(a[i]*a[j]*a[l] * (1-k1)*(1-k2)*(1-k3), 1.0f/3);
        bwr.b     += x[i] * x[j] * x[l] * std::pow(b[i]*b[j]*b[l], 1.0f/3);
        bwr.c     += x[i] * x[j] * x[l]
                   * std::pow(c[i]*c[j]*c[l], 1.0f/3) * (1-k1)*(1-k2)*(1-k3);
        bwr.alpha += x[i] * x[j] * x[l]
                   * std::pow(alfa[i]*alfa[j]*alfa[l], 1.0f/3);
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

  float roc_crit = pcam / (kGasConstantR * zcam * tcam);
  float csi      = std::pow(tcam, 6)
                 / std::pow(bwr.molar_mass, 0.5f)
                 / std::pow(pcam, 2.0f/3.0f);

  static const CountryRef kRefTable[] = {
    {"Rom\xC3\xA2nia / UE  (DIN 1343)",  2, { 0.0f,  15.0f  }, {"Nm\xC2\xB3/h", "Sm\xC2\xB3/h"}},
    {"ISO 13443  /  UK / Italia",         1, {15.0f,   0.0f  }, {"Sm\xC2\xB3/h", ""}},
    {"SUA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)", 1, {15.56f, 0.0f}, {"Sm\xC2\xB3/h", ""}},
    {"Rusia \xe2\x80\x94 GOST 30319-1",  1, {20.0f,   0.0f  }, {"m\xC2\xB3/h",  ""}},
    {"Personalizat",                      1, { 0.0f,   0.0f  }, {"m\xC2\xB3/h",  ""}},
  };
  static constexpr int kNRef = 5;

  std::printf("\n%s  Condi\xC8\x9Bii de referin\xC8\x9B\xC4\x83 volumetric\xC4\x83:%s\n\n",
              kBoldCyan, kReset);
  std::printf("%s  1.  Rom\xC3\xA2nia / UE  (DIN 1343)    \xe2\x80\x94   0\xC2\xB0""C \xC8\x99i 15\xC2\xB0""C / 101.325 kPa  [Nm\xC2\xB3/h] \xC8\x99i [Sm\xC2\xB3/h]%s\n", kCyan, kReset);
  std::printf("%s  2.  ISO 13443  /  UK / Italia       \xe2\x80\x94  15\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n",       kCyan, kReset);
  std::printf("%s  3.  SUA \xe2\x80\x94 AGA-3  (60\xC2\xB0""F)          \xe2\x80\x94  15.56\xC2\xB0""C / 101.325 kPa  [Sm\xC2\xB3/h]%s\n",  kCyan, kReset);
  std::printf("%s  4.  Rusia \xe2\x80\x94 GOST 30319-1         \xe2\x80\x94  20\xC2\xB0""C / 101.325 kPa  [m\xC2\xB3/h]%s\n",           kCyan, kReset);
  std::printf("%s  5.  Personalizat                    \xe2\x80\x94  T [\xC2\xB0""C] introdus manual  [m\xC2\xB3/h]%s\n",           kCyan, kReset);

  int ref_sel = 0;
  do {
    std::printf("%s\n  Selecta\xC8\x9Bi (1\xe2\x80\x93" "5)  > %s", kMagenta, kReset);
    ReadInt(&ref_sel);
  } while (ref_sel < 1 || ref_sel > kNRef);

  CountryRef ref = kRefTable[ref_sel - 1];
  if (ref_sel == kNRef) {
    std::printf("%s  Temperatura de referin\xC8\x9B\xC4\x83 [\xC2\xB0""C]  : %s", kBoldCyan, kReset);
    ReadFloat(&ref.t[0]);
  }

  float ror_ref[2] = {};
  for (int i = 0; i < ref.n; i++) {
    ror_ref[i] = CalcDensity(ref.t[i], 1, bwr);
    std::printf("%s  Densitatea la %5.2f\xC2\xB0""C / 101.325 kPa     :%s %s%f%s %s[kg/m\xC2\xB3]%s\n",
                kBoldCyan, ref.t[i], kReset, kBoldGreen, ror_ref[i], kReset, kCyan, kReset);
  }

  // ── Buclă exterioară: selecția dispozitivului de măsurare ─────────────────
  for (;;) {
    int   tip_raw;
    float d_int, d_orif;

    // Încearcă să refolosească configurația salvată
    {
      int   sv_tip; float sv_d_int, sv_d_orif;
      if (LoadConfig(&sv_tip, &sv_d_int, &sv_d_orif)
          && sv_tip >= kTipMin && sv_tip <= kTipMax) {
        std::printf(
            "%s\n  Configurație salvată:\n"
            "    Dispozitiv  : %s\n"
            "    D intern    : %g mm\n"
            "    D orificiu  : %g mm%s\n",
            kCyan,
            TipName(static_cast<TipDispozitiv>(sv_tip)),
            sv_d_int, sv_d_orif, kReset);
        std::printf("%s  Refolosiți configurația? [d/n]: %s", kCyan, kReset);
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
      std::printf("\n%s%s%s", kMagenta, kTipDisp, kReset);
      ReadInt(&tip_raw);
      if (tip_raw < kTipMin || tip_raw > kTipMax)
        std::printf("%s\n  Dispozitiv de strangulare nedisponibil!%s\n", kBoldRed, kReset);
    } while (tip_raw < kTipMin || tip_raw > kTipMax);

    std::printf("\n%s  D intern (20°C)   [mm] : %s", kBoldCyan, kReset); ReadFloat(&d_int);
    std::printf(  "%s  D orificiu (20°C) [mm] : %s", kBoldCyan, kReset); ReadFloat(&d_orif);

    std::printf("%s\n  Salvați configurația? [d/n]: %s", kCyan, kReset);
    if (AskYesNo()) SaveConfig(tip_raw, d_int, d_orif);

    run_inner:;
    TipDispozitiv tip = static_cast<TipDispozitiv>(tip_raw);

    // Buclă interioară: calcul pentru condiții diferite T/P cu același dispozitiv
    for (;;) {
      float temperatura, presiunea, presiunea_dif;
      std::printf("\n\n%s  Temperatura [°C]             : %s", kCyan, kReset);
      ReadFloat(&temperatura);
      std::printf("%s  Presiunea [kPa]              : %s", kCyan, kReset);
      ReadFloat(&presiunea);
      if (presiunea <= 0) {
        std::printf("\033[2J\033[H");
        std::printf("%s\n  ELCOST Impex  —  BWR Gas Flow Calculator  v2.0/2026%s\n",
                    kBoldYellow, kReset);
        return 0;
      }
      std::printf("%s  Presiunea diferențială [kPa] : %s", kCyan, kReset);
      ReadFloat(&presiunea_dif);
      std::printf("  %c\n", 7);

      float ro = CalcDensity(temperatura, presiunea / kKpaPerAtm, bwr);
      std::printf("  Densitatea (t,p)                       : %s%f%s [kg/m³]\n",
                  kBoldGreen, ro, kReset);

      float roc_red = ro / roc_crit;
      float eta = 0.0f;
      for (int i = 1; i <= kNumComponents; i++) {
        eta += (1 + 0.323f * std::log((temperatura + kKelvinOffset) / cs[i]))
             / (1 + 0.323f * std::log(kKelvinOffset / cs[i]))
             * std::sqrt((temperatura + kKelvinOffset) / kKelvinOffset)
             * et[i] * x[i] * std::sqrt(m[i]);
      }
      eta  = eta / mx;
      eta += 10.8e-8f / csi
           * std::pow(std::exp(1.439f * roc_red) - std::exp(-1.111f * roc_red),
                      1.358f);
      std::printf("  Viscozitatea dinamică (t,p)            : %s%f%s [μPa·s]\n",
                  kBoldGreen, eta * 1000000, kReset);

      FlowResult flow;
      float qm = CalcMassFlow(presiunea_dif, presiunea, temperatura,
                       tip, d_int, d_orif, ro, eta, &flow);
      if (qm == 0.0f) break;  // eroare -> reselect dispozitiv

      std::printf("  Debitul masic (t,p)                    : %s%7.4f%s [kg/s]\n",
                  kBoldGreen, qm, kReset);
      for (int i = 0; i < ref.n; i++) {
        float qhref = 3600.0f / ror_ref[i] * qm;
        std::printf("  Debitul vol.  (%5.2f\xC2\xB0""C / 101.325 kPa)  : %s%7.2f%s [%s]\n",
                    ref.t[i], kBoldGreen, qhref, kReset, ref.label[i]);
      }
      std::printf("  Debitul vol.  (t,p)                    : %s%7.2f%s [m³/h]\n",
                  kBoldGreen, 3600 / ro * qm, kReset);
      std::printf("  Viteza medie a gazului                 : %s%.2f%s [m/s]\n",
                  kBoldGreen, flow.viteza, kReset);
      std::printf("  Pierderea de presiune prin strangulare : %s%.2f%s [kPa]\n",
                  kBoldGreen, flow.pierderea, kReset);
    }
  }
}
