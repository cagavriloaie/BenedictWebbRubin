#pragma once
#include "../Common/constants.h"
#include "../Common/types.h"

// Prints the ISO 5167 error message for the given error code.
void printError(ErrorCode code, double reynoldsNum);

// Discharge coefficient C per ISO 5167-2/3/4 (Reader-Harris/Gallagher for orifices).
// pipeDiamM  — pipe inner diameter at measurement temperature [m]
// beta       — throttling ratio d/D [-]
// reynoldsNum — Reynolds number [-]
double dischargeCoefficient(DeviceType tip, double pipeDiamM, double beta, double reynoldsNum);

// Velocity coefficient α = C / √(1 − β⁴) used in the mass-flow equation.
double velocityCoefficient(DeviceType tip, double pipeDiamM, double beta, double reynoldsNum);

// Calculates mass flow rate Qm [kg/s] via ISO 5167 Reynolds iteration.
// dp, p  — differential and absolute pressure [kPa]
// t      — temperature [°C]
// dInt, dOrif — pipe and orifice diameters at 20 °C reference [mm]
// ro     — gas density [kg/m³]
// eta    — dynamic viscosity [Pa·s]
// Returns 0.0 on range or convergence error; fills *out with hydraulic results.
double calcMassFlow(double dp, double p, double t,
                    DeviceType tip, double dInt, double dOrif,
                    double ro, double eta, FlowResult* out,
                    double kappa = KAPPA_DEFAULT, int* iters = nullptr);
