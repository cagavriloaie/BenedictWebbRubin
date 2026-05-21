#pragma once
#include "../Common/constants.h"
#include "../Common/types.h"

// File names used for persistence
inline constexpr const char* COMP_FILE = "BWRSflow_comp.dat";
inline constexpr const char* CONF_FILE = "BWRSflow_conf.dat";

// Component display names (1-indexed; index 0 is nullptr)
extern const char* const COMP_NAMES[ARRAY_SIZE];

// Clears the screen, prints the startup banner (header, models, references),
// and optionally runs the validation test suite if the user confirms.
void PrintBanner();

// Clears the screen, shows the exit banner, and terminates the process (ESC key).
void ExitApp();

// Reads a real number interactively with echo, backspace, and ESC support.
// Returns false if the user presses Backspace on an empty field (go-back signal).
bool ReadDouble(double* val);

// Reads a single digit in [lo, hi] without requiring ENTER; ESC exits.
int ReadChoice(int lo, int hi);

// Blocks until y or n is pressed; ESC exits the application.
bool AskYesNo();

// Blocks until any key is pressed.
void WaitKey();

// Returns the count of UTF-8 continuation bytes (10xxxxxx) in s.
// strlen(s) - Utf8ExtraBytes(s) == number of displayed characters.
int Utf8ExtraBytes(const char* s);

// Prints the mixture composition table aligned in the console.
void PrintComposition(const double* x);

// Loads mixture composition from COMP_FILE.
// Returns true if all NUM_COMPONENTS values were read successfully.
bool LoadComposition(double* x);

// Saves mixture composition to COMP_FILE.
void SaveComposition(const double* x);

// Loads device configuration from CONF_FILE.
// Returns true if tip, pipe diameter, and orifice diameter were read successfully.
bool LoadConfig(int* tip_raw, double* d_int, double* d_orif);

// Saves device configuration to CONF_FILE.
void SaveConfig(int tip_raw, double d_int, double d_orif);

// Initialises the Windows console (UTF-8, Consolas font, VT/ANSI).
void InitConsole();

// Loads or prompts for the mixture molar fractions; saves on request.
void ReadComposition(double* x);

// Displays molar mass, prompts for reference conditions, computes reference densities.
void SelectRefConditions(const BwrConst& bwr, CountryRef* out_ref, double ror_ref[2]);

// Prints the full results block for one T/p/Δp calculation.
void PrintFlowResults(DeviceType tip, double d_int, double d_orif,
                      double temperature, double pressure, double pressure_diff,
                      double ro, double eta, double qm,
                      const FlowResult& flow, const BwrConst& bwr,
                      const CountryRef& ref, const double* ror_ref,
                      const double* x);

// Prints the CPU profile and embedded-platform estimate table.
void PrintCpuProfile(double t_mix_us, double t_rho_us, double t_eta_us,
                     double t_qm_us, int iter_rho, int iter_qm);

// Prompts for uncertainty confirmation, then prints the ISO 5167-1 budget.
void PrintUncertainty(DeviceType tip, const FlowResult& flow, double qm,
                      const CountryRef& ref, const double* ror_ref);
