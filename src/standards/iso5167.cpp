#include "iso5167.h"
#include "../Common/constants.h"
#include "../Common/types.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

// Prints a coloured error message for the given ISO 5167 validation failure.
void printError(ErrorCode code, double reynoldsNum) {
  std::printf("%s", COLOR_BOLD_RED);
  switch (code) {
    case ErrorCode::PIPE_DIAMETER:
      std::printf("\n  Pipe inner diameter is out of range\n");
      break;
    case ErrorCode::ORIFICE_DIAMETER:
      std::printf("\n  Throttling device orifice diameter is out of range\n");
      break;
    case ErrorCode::THROTTLE_RATIO:
      std::printf("\n  Throttling ratio " U_BETA " is out of range\n");
      break;
    case ErrorCode::REYNOLDS:
      std::printf("\n  Reynolds number (%g) is out of range\n", reynoldsNum);
      break;
    case ErrorCode::NUMERIC:
      std::printf("\n  Numerical error — check that inputs are in a physically valid range\n");
      break;
    default:
      break;
  }
  std::printf("%s", COLOR_RESET);
}

// Returns the discharge coefficient C for the given throttling device per ISO 5167-2/3/4.
// pipeDiamM [m], beta = d/D, reynoldsNum — pipe Reynolds number.
double dischargeCoefficient(DeviceType tip, double pipeDiamM, double beta, double reynoldsNum) {
  double dischCoef = 0.0;
  switch (tip) {
    case DeviceType::ORIFICE_CORNER:
    case DeviceType::ORIFICE_FLANGE:
    case DeviceType::ORIFICE_DD2: {
      double tapL1, tapL2Prime;
      if (tip == DeviceType::ORIFICE_CORNER) {
        tapL1 = 0.0;  tapL2Prime = 0.0;
      } else if (tip == DeviceType::ORIFICE_FLANGE) {
        tapL1 = FLANGE_TAP_M / pipeDiamM;  tapL2Prime = tapL1;
      } else {
        tapL1 = 1.0;  tapL2Prime = 0.50;
      }
      double rhgA  = std::pow(19000.0 * beta / reynoldsNum, 0.8);
      double m2Tap = 2.0 * tapL2Prime / (1.0 - beta);
      dischCoef = 0.5961
               + 0.0261 * beta * beta
               - 0.216  * std::pow(beta, 8.0)
               + 0.000521 * std::pow(1.0e6 * beta / reynoldsNum, 0.7)
               + (0.0188 + 0.0063 * rhgA)
                 * std::pow(beta, 3.5) * std::pow(1.0e6 / reynoldsNum, 0.3)
               + (0.043 + 0.080 * std::exp(-10.0 * tapL1)
                         - 0.123 * std::exp( -7.0 * tapL1))
                 * (1.0 - 0.11 * rhgA) * std::pow(beta, 4.0)
                 / (1.0 - std::pow(beta, 4.0))
               - 0.031 * (m2Tap - 0.8 * std::pow(m2Tap, 1.1))
                 * std::pow(beta, 1.3);
      if (pipeDiamM < 0.07112)  // correction for D < 71.12 mm
        dischCoef += 0.011 * (0.75 - beta) * (2.8 - pipeDiamM / FLANGE_TAP_M);
      break;
    }
    case DeviceType::NOZZLE_ISA:
      dischCoef = 0.99 - 0.2262 * std::pow(beta, 4.1)
               + (0.000215 - 0.001125 * beta + 0.00249 * std::pow(beta, 4.7))
               * std::pow(REYNOLDS_ISO_REF / reynoldsNum, 1.15);
      break;
    case DeviceType::NOZZLE_LONG_RADIUS:
      dischCoef = 0.9965
               - 0.00653 * std::sqrt(beta)
               * std::sqrt(REYNOLDS_ISO_REF / reynoldsNum);
      break;
    case DeviceType::VENTURI_ROUGH_CAST:
      dischCoef = 0.984;
      break;
    case DeviceType::VENTURI_MACHINED:
      dischCoef = 0.995;
      break;
    case DeviceType::VENTURI_WELDED_SHEET:
      dischCoef = 0.985;
      break;
    case DeviceType::VENTURI_NOZZLE:
      dischCoef = 0.9858 - 0.196 * std::pow(beta, 4.5);
      break;
    default:
      break;
  }
  return dischCoef;
}

// Returns the velocity coefficient α = C / √(1−β⁴) used in the mass-flow formula.
double velocityCoefficient(DeviceType tip, double pipeDiamM, double beta, double reynoldsNum) {
  return dischargeCoefficient(tip, pipeDiamM, beta, reynoldsNum)
       / std::sqrt(1 - std::pow(beta, 4));
}

// Computes mass flow rate [kg/s] per ISO 5167 via Reynolds iteration.
// dp [kPa], p [kPa], t [°C], dInt/dOrif [mm], ro [kg/m³], eta [mPa·s].
// Validates geometry and Reynolds range; writes detailed results to *out.
double calcMassFlow(double dp, double p, double t,
                    DeviceType tip, double dInt, double dOrif,
                    double ro, double eta, FlowResult* out,
                    double kappa, int* iters) {
  if (!out || !std::isfinite(dp) || !std::isfinite(p) || !std::isfinite(t)
      || !std::isfinite(dInt) || !std::isfinite(dOrif)
      || !std::isfinite(ro) || !std::isfinite(eta)
      || dp <= 0.0 || p <= 0.0 || dp >= p || ro <= 0.0 || eta <= 0.0) {
    printError(ErrorCode::NUMERIC, 0);
    return 0;
  }

  double dPipeT  = dInt  * (1 + THERMAL_EXP_PIPE    * (t - REF_TEMP_CELSIUS));
  double dOrifT  = dOrif * (1 + THERMAL_EXP_ORIFICE * (t - REF_TEMP_CELSIUS));

  if ((dPipeT < 50)
      || (dPipeT > 1000 && tip == DeviceType::ORIFICE_CORNER)
      || (dPipeT > 760  && (tip == DeviceType::ORIFICE_FLANGE
                            || tip == DeviceType::ORIFICE_DD2))
      || (dPipeT > 500  && tip == DeviceType::NOZZLE_ISA)
      || (dPipeT > 630  && tip == DeviceType::NOZZLE_LONG_RADIUS)) {
    printError(ErrorCode::PIPE_DIAMETER, 0);
    return 0;
  }

  if ((dPipeT > 800  && tip == DeviceType::VENTURI_ROUGH_CAST)
      || (dPipeT < 100 && tip == DeviceType::VENTURI_ROUGH_CAST)
      || (dPipeT > 250 && tip == DeviceType::VENTURI_MACHINED)
      || (dPipeT < 50  && tip == DeviceType::VENTURI_MACHINED)
      || (dPipeT > 1200 && tip == DeviceType::VENTURI_WELDED_SHEET)
      || (dPipeT < 200  && tip == DeviceType::VENTURI_WELDED_SHEET)
      || (dPipeT > 500  && tip == DeviceType::VENTURI_NOZZLE)
      || (dPipeT < 65   && tip == DeviceType::VENTURI_NOZZLE)) {
    printError(ErrorCode::PIPE_DIAMETER, 0);
    return 0;
  }

  if ((isDiaphragm(tip) && dOrifT < 12.5)
      || (tip == DeviceType::VENTURI_NOZZLE && dOrifT <= 50)) {
    printError(ErrorCode::ORIFICE_DIAMETER, 0);
    return 0;
  }

  double beta = dOrifT / dPipeT;

  if (((beta < 0.23  || beta > 0.8)  && tip == DeviceType::ORIFICE_CORNER)
      || ((beta < 0.20 || beta > 0.75) && (tip == DeviceType::ORIFICE_FLANGE
                                            || tip == DeviceType::ORIFICE_DD2))
      || ((beta < 0.3  || beta > 0.8) && tip == DeviceType::NOZZLE_ISA)
      || ((beta < 0.2  || beta > 0.8) && tip == DeviceType::NOZZLE_LONG_RADIUS)) {
    printError(ErrorCode::THROTTLE_RATIO, 0);
    return 0;
  }
  if (((beta < 0.3   || beta > 0.75)  && tip == DeviceType::VENTURI_ROUGH_CAST)
      || ((beta < 0.4   || beta > 0.75) && tip == DeviceType::VENTURI_MACHINED)
      || ((beta < 0.4   || beta > 0.7)  && tip == DeviceType::VENTURI_WELDED_SHEET)
      || ((beta < 0.316 || beta > 0.775) && tip == DeviceType::VENTURI_NOZZLE)) {
    printError(ErrorCode::THROTTLE_RATIO, 0);
    return 0;
  }

  double expandFact;
  if (isDiaphragm(tip)) {
    expandFact = 1 - (0.41 + 0.35 * std::pow(beta, 4)) * dp / p / kappa;
  } else {
    double pressureRatio = 1 - dp / p;
    double exp2Kappa     = 2.0 / kappa;
    double expKappaM1    = (kappa - 1.0) / kappa;
    double betaPow4      = std::pow(beta, 4);
    double yKappa        = std::pow(pressureRatio, exp2Kappa);
    expandFact = std::sqrt(kappa * yKappa / (kappa - 1.0)
        * (1.0 - betaPow4) / (1.0 - betaPow4 * yKappa)
        * (1.0 - std::pow(pressureRatio, expKappaM1)) / (1.0 - pressureRatio));
  }
  if (!std::isfinite(expandFact) || expandFact <= 0.0) {
    printError(ErrorCode::NUMERIC, 0);
    return 0;
  }

  dPipeT /= 1000;
  dOrifT /= 1000;

  double reynoldsNum = INITIAL_REYNOLDS;
  double qmCurr = 0;
  double qmPrev = 0;
  double velCoef = 0;
  int nReynolds = 0;
  bool converged = false;
  for (; nReynolds < MAX_REYNOLDS_ITER; ) {
    nReynolds++;
    qmPrev  = qmCurr;
    velCoef = velocityCoefficient(tip, dPipeT, beta, reynoldsNum);
    qmCurr  = velCoef * expandFact * PI / 4 * std::pow(dOrifT, 2)
              * std::sqrt(FACTOR_2KPA * dp * ro);
    reynoldsNum = 4 * qmCurr / (dPipeT * PI * eta);
    if (!std::isfinite(velCoef) || !std::isfinite(qmCurr) || !std::isfinite(reynoldsNum)) {
      printError(ErrorCode::NUMERIC, 0);
      return 0;
    }
    if (std::fabs((qmCurr - qmPrev) / std::max(qmCurr, 1e-10)) <= REYNOLDS_TOLERANCE) {
      converged = true;
      break;
    }
  }
  if (iters) *iters = nReynolds;
  if (!converged) {
    printError(ErrorCode::NUMERIC, 0);
    return 0;
  }

  bool reynoldsValid = false;
  switch (tip) {
    case DeviceType::ORIFICE_CORNER:
      if ((5000  <= reynoldsNum) && (reynoldsNum <= 1e8) && (0.23 <= beta) && (beta < 0.45))
        reynoldsValid = true;
      if ((10000 <= reynoldsNum) && (reynoldsNum <= 1e8) && (0.45 <= beta) && (beta < 0.77))
        reynoldsValid = true;
      if ((20000 <= reynoldsNum) && (reynoldsNum <= 1e8) && (0.77 <= beta) && (beta <= 0.80))
        reynoldsValid = true;
      break;
    case DeviceType::ORIFICE_FLANGE:
    case DeviceType::ORIFICE_DD2:
      if ((1.26e6 * beta * beta * dPipeT <= reynoldsNum) && (reynoldsNum <= 1e8))
        reynoldsValid = true;
      break;
    case DeviceType::NOZZLE_ISA:
      if ((70000 <= reynoldsNum) && (reynoldsNum <= 1e7) && (0.30 <= beta) && (beta < 0.44))
        reynoldsValid = true;
      if ((20000 <= reynoldsNum) && (reynoldsNum <= 1e7) && (0.44 <= beta) && (beta <= 0.80))
        reynoldsValid = true;
      break;
    case DeviceType::NOZZLE_LONG_RADIUS:
      if ((10000 <= reynoldsNum) && (reynoldsNum <= 1e7))
        reynoldsValid = true;
      break;
    case DeviceType::VENTURI_ROUGH_CAST:
      if ((200000 <= reynoldsNum) && (reynoldsNum <= 2000000))
        reynoldsValid = true;
      break;
    case DeviceType::VENTURI_MACHINED:
      if ((200000 <= reynoldsNum) && (reynoldsNum <= 1000000))
        reynoldsValid = true;
      break;
    case DeviceType::VENTURI_WELDED_SHEET:
      if ((200000 <= reynoldsNum) && (reynoldsNum <= 2000000))
        reynoldsValid = true;
      break;
    case DeviceType::VENTURI_NOZZLE:
      if ((150000 <= reynoldsNum) && (reynoldsNum <= 2000000))
        reynoldsValid = true;
      break;
    default:
      break;
  }

  if (!reynoldsValid) {
    printError(ErrorCode::REYNOLDS, reynoldsNum);
    return 0;
  }

  out->m_velocity = 4 * qmCurr / PI / dPipeT / dPipeT / ro;

  // Permanent pressure loss fraction
  double lossFraction;
  if (tip == DeviceType::VENTURI_ROUGH_CAST) {
    lossFraction = 0.15;
  } else if (tip == DeviceType::VENTURI_MACHINED) {
    lossFraction = 0.08;
  } else if (tip == DeviceType::VENTURI_WELDED_SHEET) {
    lossFraction = 0.15;
  } else {
    double dischCoefC = out->m_coefC;
    double betaSq     = beta * beta;
    double sqRoot     = std::sqrt(1.0 - betaSq * betaSq * (1.0 - dischCoefC * dischCoefC));
    lossFraction = (sqRoot - dischCoefC * betaSq) / (sqRoot + dischCoefC * betaSq);
  }
  out->m_pressureLoss = lossFraction * dp;
  out->m_beta         = beta;
  out->m_reynolds     = reynoldsNum;
  out->m_coefC        = dischargeCoefficient(tip, dPipeT, beta, reynoldsNum);
  out->m_epsilon      = expandFact;
  out->m_dWorkingMm   = dOrifT * 1000.0;
  out->m_DWorkingMm   = dPipeT * 1000.0;
  out->m_kappa        = kappa;
  return qmCurr;
}
