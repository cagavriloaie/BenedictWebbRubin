# BWRSflow — Gas Flow Calculator

**Version 3.0 | 05.2026**  
Eng. Agavriloaie Constantin *(original R 01.2004)*  
ELCOST Impex · office@elcost.ro  
Source: https://github.com/cagavriloaie/BenedictWebbRubin

---

## Table of Contents

1. [General Description](#1-general-description)
2. [Input Data](#2-input-data)
3. [Output Data](#3-output-data)
4. [Calculation Models](#4-calculation-models)
   - 4.1 [BWRS Equation of State](#41-bwrs-equation-of-state)
   - 4.2 [Mixing Rules](#42-bwrs-mixing-rules)
   - 4.3 [Density Calculation](#43-density-calculation)
   - 4.4 [Compressibility Factor Z](#44-compressibility-factor-z)
   - 4.5 [Dynamic Viscosity](#45-dynamic-viscosity)
   - 4.6 [Mass Flow Rate ISO 5167](#46-mass-flow-rate--iso-5167)
   - 4.7 [Discharge Coefficients](#47-discharge-coefficients-by-device-type)
   - 4.8 [Expansibility Factor](#48-expansibility-factor-ε)
   - 4.9 [Thermal Corrections for Diameters](#49-thermal-corrections-for-diameters)
   - 4.10 [Calorific Quantities](#410-calorific-quantities)
5. [ISO 5167 Validation Limits](#5-iso-5167-validation-limits)
6. [Measurement Uncertainty](#6-measurement-uncertainty)
7. [How to Use the Application](#7-how-to-use-the-application)
8. [Data Persistence](#8-data-persistence)
9. [Available Reference Conditions](#9-available-reference-conditions)
10. [Component Database](#10-component-database)
11. [Standards and References](#11-standards-and-references)
12. [Frequently Asked Questions](#12-frequently-asked-questions)

---

## 1. General Description

**BWRSflow** calculates the density, dynamic viscosity, and flow rate of gases through differential pressure devices (orifice plates, nozzles, Venturi tubes), for gas mixtures of up to 35 components.

**Scope of application:**
- Natural gas, refinery gas, industrial gases, air, hydrogen, CO₂
- Operating conditions: pressures up to several hundred bar, temperatures −50 … +200 °C
- 9 device types in accordance with ISO 5167-2/3/4:2003

**Implemented methods:**
| Property | Model |
|----------|-------|
| Density ρ | BWRS — Benedict–Webb–Rubin–Starling (Starling 1973, 11 parameters) |
| Viscosity η | Chapman–Enskog (dilute gas) + Lucas/Chung–Lee–Starling high-density correction |
| Mass flow rate Qm | ISO 5167-2/3/4:2003 — Reader–Harris/Gallagher equation + Reynolds iteration |
| Calorific value | ISO 6976:2016 — gross higher heating value per component |

---

## 2. Input Data

### 2.1 Mixture Composition
Molar fractions of the **35 components** (sum = 1.000000).  
The application verifies the sum and offers automatic normalisation in case of deviation.  
The composition can be saved/loaded from the file `BWRSflow_comp.dat`.

### 2.2 Differential Pressure Device Type
| No. | Device |
|-----|--------|
| 1 | Orifice plate — corner taps |
| 2 | Orifice plate — flange taps |
| 3 | Orifice plate — D and D/2 taps |
| 4 | ISA 1932 nozzle |
| 5 | Long-radius nozzle |
| 6 | Classical Venturi tube — rough-cast convergent |
| 7 | Classical Venturi tube — machined convergent |
| 8 | Classical Venturi tube — welded sheet-iron convergent |
| 9 | Venturi nozzle |

### 2.3 Device Geometry
- **D** — internal pipe diameter at 20 °C [mm]
- **d** — orifice/throat diameter at 20 °C [mm]
- The configuration can be saved/loaded from `BWRSflow_conf.dat`

### 2.4 Measurement Conditions
- **T** — gas temperature [°C]
- **p** — absolute pressure [kPa]
- **Δp** — differential pressure [kPa]

---

## 3. Output Data

### Measurement Conditions
| Quantity | Symbol | Unit |
|----------|--------|------|
| Density | ρ(T, p) | kg/m³ |
| Dynamic viscosity | η(T, p) | μPa·s |
| Kinematic viscosity | ν = η/ρ | mm²/s |
| Compressibility factor | Z(T, p) | — |
| Speed of sound | c = √(κ·p/ρ) | m/s |
| Mach number | Ma = v/c | — |
| Pressure ratio | Δp/p | % |

### Flow Rates
| Quantity | Unit |
|----------|------|
| Mass flow rate Qm | kg/s and kg/h |
| Volumetric flow at reference conditions | Nm³/h, Sm³/h (country-dependent) |
| Volumetric flow at (T, p) | m³/h |

### Hydraulics
| Quantity | Unit |
|----------|------|
| Mean pipe velocity v | m/s |
| Permanent pressure loss | kPa |
| Diameter ratio β | — |
| Reynolds number Re | — |
| Discharge coefficient C | — |
| Expansibility factor ε | — |
| Isentropic exponent κ | — |

### Calorific (ISO 6976)
| Quantity | Unit |
|----------|------|
| Mixture HHV (per mass) | MJ/kg |
| Mixture HHV (volume at 0 °C) | MJ/m³ |
| Relative density d | — |
| Wobbe index Ws (0 °C) | MJ/m³ |
| Energy flow rate Qe | kW and MJ/h |

---

## 4. Calculation Models

### 4.1 BWRS Equation of State

The Benedict–Webb–Rubin–Starling equation (Starling, 1973):

```
p = ρRT
  + (B₀RT − A₀ − C₀/T² + D₀/T³ − E₀/T⁴) · ρ²
  + (bRT − a − d/T) · ρ³
  + α(a + d/T) · ρ⁶
  + (c/T²) · ρ³ · (1 + γρ²) · exp(−γρ²)
```

where R = 0.082057366 L·atm/(mol·K), T in Kelvin, ρ in mol/L, p in atm.

**Parameters per pure component:**

| Parameter | Meaning |
|-----------|---------|
| A₀, B₀, C₀, D₀, E₀ | Second-virial coefficients and thermal corrections |
| a, b, c, d | Third-virial coefficients |
| α | Amplitude of the ρ⁶ term |
| γ | Gaussian parameter |

### 4.2 BWRS Mixing Rules

**Quadratic rule** (for A₀, B₀, C₀, D₀, E₀, γ):

Binary interaction coefficient:
```
kᵢⱼ = 1 − 8·√(Vᵢ·Vⱼ) / (∛Vᵢ + ∛Vⱼ)³
```

```
A₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(A₀ᵢ·A₀ⱼ) · (1 − kᵢⱼ)
B₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(B₀ᵢ·B₀ⱼ)
C₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(C₀ᵢ·C₀ⱼ) · (1 − kᵢⱼ)³
...
```

**Cubic rule** (for a, b, c, d, α):
```
a_mix = Σᵢ Σⱼ Σₗ xᵢ·xⱼ·xₗ · ∛(aᵢ·aⱼ·aₗ) · (1−kᵢⱼ)·(1−kᵢₗ)·(1−kⱼₗ)
```

### 4.3 Density Calculation

The BWRS equation is solved for ρ by **bisection**:

```
ρ_lo = 10⁻⁴ mol/L,   ρ_hi = 25 mol/L

Repeat:
  ρ_mid  = (ρ_lo + ρ_hi) / 2
  p_calc = BWRS(ρ_mid, T, mixture_params)
  if p_calc > p  →  ρ_hi = ρ_mid
  else           →  ρ_lo = ρ_mid
until |p − p_calc| < 5×10⁻⁴ atm   (max 200 iterations)
```

### 4.4 Compressibility Factor Z

```
Z(T, p) = p · M / (ρ · R_SI · T)

where:
  M     — molar mass of the mixture [kg/mol]
  R_SI  = 8.31446 J/(mol·K)
  ρ     — density [kg/m³]
  p     — pressure [Pa]
  T     — temperature [K]
```

### 4.5 Dynamic Viscosity

**Dilute gas — Chapman–Enskog:**
```
η_dilute[i] = √(T/Tcᵢ) · (1 + 0.323·ln(T/Tcᵢ)) / (1 + 0.323·ln(273.15/Tcᵢ)) · etᵢ · √mᵢ

η_mix_dilute = Σᵢ (η_dilute[i] · xᵢ · √mᵢ) / √M_mix
```

**High-density correction — Lucas/Chung–Lee–Starling:**
```
η = η_mix_dilute + (10.8×10⁻⁸ / ξ) · [exp(1.439·ρ_red) − exp(−1.111·ρ_red)]^1.358

where:
  ρ_red = ρ / ρ_crit         (reduced density)
  ξ     = Tc^6 / (√M · Pc^(2/3))
  ρ_crit = Pc / (R · Zc · Tc)
```

### 4.6 Mass Flow Rate — ISO 5167

**Primary equation:**
```
Qm = α · ε · (π/4) · d² · √(2000 · Δp · ρ)   [kg/s]

where:
  α = C / √(1 − β⁴)    — velocity of approach factor
  ε                     — expansibility factor
  d                     — orifice diameter at temperature T [m]
  β = d/D               — diameter ratio
```

**Reynolds iteration:**
```
Re₀ = 10⁶

At iteration k:
  αₖ  = C(Reₖ) / √(1 − β⁴)
  Qmₖ = αₖ · ε · (π/4) · d² · √(2000·Δp·ρ)
  Reₖ₊₁ = 4·Qmₖ / (π·D·η)

Convergence: |ΔQm/Qm| < 10⁻⁶   (max 100 iterations)
```

**Mean pipe velocity and permanent pressure loss:**
```
v = 4·Qm / (π·D²·ρ)   [m/s]

Orifice plates / nozzles:
  Δp_perm = (√(1−β⁴(1−C²)) − C·β²) / (√(1−β⁴(1−C²)) + C·β²) · Δp

Classical Venturi:  ~8–15% · Δp  (depends on type)
```

### 4.7 Discharge Coefficients by Device Type

**Orifice plate — Reader–Harris/Gallagher (ISO 5167-2):**
```
C = 0.5961 + 0.0261·β² − 0.216·β⁸
  + 0.000521·(10⁶·β/Re)^0.7
  + (0.0188 + 0.0063·A)·β^3.5·(10⁶/Re)^0.3
  + (0.043 + 0.080·e^(−10L₁) − 0.123·e^(−7L₁))·(1−0.11A)·β⁴/(1−β⁴)
  − 0.031·(M₂ − 0.8·M₂^1.1)·β^1.3
  [+ correction for D < 71.12 mm]

A = (19000·β/Re)^0.8
M₂ = 2·L₂' / (1−β)
```

| Tap type | L₁ | L₂' |
|----------|----|-----|
| Corner | 0 | 0 |
| Flange | 25.4/D [mm→—] | 25.4/D |
| D and D/2 | 1.0 | 0.47 |

**ISA 1932 nozzle (ISO 5167-3):**
```
C = 0.99 − 0.2262·β^4.1 + (0.000215 − 0.001125·β + 0.00249·β^4.7)·(10⁶/Re)^1.15
```

**Long-radius nozzle (ISO 5167-3):**
```
C = 0.9965 − 0.00653·√β·√(10⁶/Re)
```

**Venturi tubes — constant C (ISO 5167-4):**
| Type | C |
|------|---|
| Rough-cast | 0.984 |
| Machined | 0.995 |
| Welded sheet-iron | 0.985 |
| Venturi nozzle | 0.9858 − 0.196·β^4.5 |

### 4.8 Expansibility Factor ε

**Orifice plates (ISO 5167-2):**
```
ε = 1 − (0.41 + 0.35·β⁴) · Δp / (κ·p)
```

**Nozzles and Venturi tubes (ISO 5167-3/4):**
```
y = 1 − Δp/p

ε = √[ κ·y^(2/κ) / (κ−1) · (1−β⁴)/(1−β⁴·y^(2/κ)) · (1−y^((κ−1)/κ))/(1−y) ]
```

**Isentropic exponent κ:**
```
κ = Cp_mix / (Cp_mix − R_SI)
Cp_mix = Σᵢ xᵢ · Cp°ᵢ   [J/(mol·K) at 20 °C]
```

### 4.9 Thermal Corrections for Diameters

```
D(T) = D_ref · (1 + 12.2×10⁻⁶ · (T − 20))   [carbon steel pipe]
d(T) = d_ref · (1 + 16.5×10⁻⁶ · (T − 20))   [stainless steel orifice plate]
```

A warning is displayed if the correction exceeds ±0.10%.

### 4.10 Calorific Quantities

**Higher heating value of the mixture (ISO 6976:2016):**
```
HHV_mol = Σᵢ xᵢ · HHV_i   [MJ/mol]

HHV_mass = HHV_mol / (M_mix × 10⁻³)           [MJ/kg]
HHV_vol  = HHV_mol / 0.022414                  [MJ/m³ at 0 °C / 101.325 kPa]
```

**Wobbe index (ISO 13686):**
```
d_rel = M_mix / 28.9625   (relative density with respect to air)
Ws    = HHV_vol / √d_rel  [MJ/m³]
```

**Energy flow rate:**
```
Qe = Qm · HHV_mass · 1000   [kW]
Qe = Qe · 3.6               [MJ/h]
```

---

## 5. ISO 5167 Validation Limits

The application checks automatically and displays an error if the limits are exceeded.

### Pipe Diameter D [mm]
| Device | D_min | D_max |
|--------|-------|-------|
| Orifice (corner / flange) | 50 | 1000 |
| Orifice (D and D/2) | 50 | 760 |
| ISA 1932 nozzle | 50 | 500 |
| Long-radius nozzle | 50 | 630 |
| Rough-cast Venturi | 100 | 800 |
| Machined Venturi | 50 | 250 |
| Welded Venturi | 200 | 1200 |
| Venturi nozzle | 65 | 500 |

### Diameter Ratio β = d/D
| Device | β_min | β_max |
|--------|-------|-------|
| Orifice plate | 0.23 | 0.80 |
| ISA 1932 nozzle | 0.30 | 0.80 |
| Long-radius nozzle | 0.20 | 0.80 |
| Rough-cast / welded Venturi | 0.30–0.40 | 0.70–0.75 |
| Machined Venturi | 0.40 | 0.75 |
| Venturi nozzle | 0.316 | 0.775 |

### Reynolds Number Re
| Device | Re_min |
|--------|--------|
| Orifice β < 0.45 | 5 000 |
| Orifice 0.45 ≤ β | 10 000–20 000 |
| ISA 1932 nozzle | 20 000–70 000 |
| Long-radius nozzle | 10 000 |
| Classical Venturi | 200 000 |
| Venturi nozzle | 150 000 |

---

## 6. Measurement Uncertainty

The uncertainty calculation is optional (the application prompts the user at the end of each result set).

**Law of propagation of uncertainties (ISO 5167-1):**
```
u²(Qm)/Qm² = u²(C) + u²(ε)
            + [2/(1−β⁴)]² · u²(d)
            + [2β⁴/(1−β⁴)]² · u²(D)
            + (1/2)² · u²(Δp)
            + (1/2)² · u²(ρ)
```

**Expanded uncertainty:** U(Qm) = 2·u(Qm) at k = 2, 95 % confidence level.

**Default values** (user-adjustable):

| Source | Default |
|--------|---------|
| u(C) — discharge coefficient | 0.50 % (orifice β ≤ 0.60) |
| u(ε) — expansibility factor | 0.10 % |
| u(d) — orifice diameter | 0.03 % |
| u(D) — pipe diameter | 0.10 % |
| u(Δp) — differential pressure | 0.20 % |
| u(ρ) — BWRS density | 0.30 % |

---

## 7. How to Use the Application

### Starting
```
Windows:  Build\Release\BWRSflow.exe
Linux:    ./Build/BWRSflow
```

On first launch, an optional **validation test** is offered with 8 reference compositions.

### Workflow

```
1. Gas composition
   ├── If BWRSflow_comp.dat exists → reuse? [y/n]
   └── Otherwise → manual entry (35 molar fractions)
       └── Automatic normalisation if sum ≠ 1

2. Reference conditions
   ├── 1. Romania / EU  (0 °C and 15 °C / 101.325 kPa)
   ├── 2. ISO 13443     (0 °C / 101.325 kPa)
   ├── 3. AGA-3 / USA   (60 °F / 14.73 psia)
   ├── 4. GOST (Russia) (20 °C / 101.325 kPa)
   └── 5. Custom        (user-selected temperature)

3. Device type (1–9)
   └── If BWRSflow_conf.dat exists → reuse? [y/n]

4. Diameters D and d [mm] at 20 °C
   └── Automatic check β ∈ [0.10, 0.80]

5. Measurement conditions T, p, Δp
   └── Calculate → Display results

6. Navigation:
   ├── 1 → new T/p/Δp (same D, d, device type)
   ├── 2 → new D, d (same device type)
   └── 3 → new device type
```

### Special Keys
| Key | Action |
|-----|--------|
| ESC | Exit the application (at any time) |
| ← (Backspace on empty field) | Return to the previous field |
| Enter / 0 | Accept default value (in uncertainty prompts) |

---

## 8. Data Persistence

| File | Contents | Created when |
|------|----------|-------------|
| `BWRSflow_comp.dat` | 35 molar fractions | On request, after composition entry |
| `BWRSflow_conf.dat` | Device type, D [mm], d [mm] | Automatically after each calculation |

Files are created in the **current working directory** (from which the application is launched).

---

## 9. Available Reference Conditions

| Option | Temperatures | Pressure | Flow unit |
|--------|-------------|----------|-----------|
| Romania / EU | 0 °C and 15 °C | 101.325 kPa | Nm³/h and Sm³/h |
| ISO 13443 | 0 °C | 101.325 kPa | Nm³/h |
| AGA-3 / USA | 60 °F (15.56 °C) | 101.325 kPa | Sm³/h |
| GOST (Russia) | 20 °C | 101.325 kPa | Sm³/h |
| Custom | User-selected T | 101.325 kPa | m³/h |

---

## 10. Component Database

35 components, indexed 1–35:

| No. | Component | Formula |
|-----|-----------|---------|
| 1 | Methane | CH₄ |
| 2 | Ethane | C₂H₆ |
| 3 | Propane | C₃H₈ |
| 4 | Isobutane | i-C₄H₁₀ |
| 5 | n-Butane | n-C₄H₁₀ |
| 6 | Neopentane | neo-C₅H₁₂ |
| 7 | Isopentane | i-C₅H₁₂ |
| 8 | n-Pentane | n-C₅H₁₂ |
| 9–13 | C₆ isomers | C₆H₁₄ |
| 14–19 | C₇ isomers | C₇H₁₆ |
| 20–21 | C₈ isomers | C₈H₁₈ |
| 22 | Benzene | C₆H₆ |
| 23 | Toluene | C₇H₈ |
| 24 | Hydrogen | H₂ |
| 25 | Carbon monoxide | CO |
| 26 | Hydrogen sulphide | H₂S |
| 27 | Helium | He |
| 28 | Argon | Ar |
| 29 | Nitrogen | N₂ |
| 30 | Oxygen | O₂ |
| 31 | Carbon dioxide | CO₂ |
| 32 | Ethylene | C₂H₄ |
| 33 | Propylene | C₃H₆ |
| 34 | Ammonia | NH₃ |
| 35 | Acetylene | C₂H₂ |

Per component the following are stored: 11 BWRS parameters, Tc, Pc, Zc, M, 2 viscosity parameters, Cp°, HHV.

---

## 11. Standards and References

| Reference | Description |
|-----------|-------------|
| **ISO 5167-2:2003** | Orifice plates — Reader–Harris/Gallagher equation |
| **ISO 5167-3:2003** | Nozzles and Venturi nozzles |
| **ISO 5167-4:2003** | Classical Venturi tubes |
| **ISO 6976:2016** | Higher and lower calorific values of natural gases |
| **ISO 13443:1996** | Standard reference conditions for natural gas |
| **ISO 13686:1998** | Natural gas quality — Wobbe index |
| **Starling K.E. (1973)** | *Fluid Thermodynamic Properties for Light Petroleum Systems*, Gulf Publishing Co., Houston — BWRS parameters for 11 components C1–C8 |
| **Nishiumi H., Saito S. (1975)** | J. Chem. Eng. Japan **8**(5), 356–360 — BWRS extension with D₀ and E₀ parameters |
| **Benedict M., Webb G.B., Rubin L.C. (1940)** | J. Chem. Phys. **8**, 334 — original BWR equation |
| **Chapman S., Cowling T.G.** | *The Mathematical Theory of Non-Uniform Gases* — dilute-gas viscosity |
| **Lucas K. (1980)** | High-density viscosity correction |
| **Chung T.H., Lee L.L., Starling K.E. (1984)** | High-density viscosity correction |

---

## 12. Frequently Asked Questions

**Q: Why does the sum of molar fractions not need to be exactly 1.000000?**  
A: The application accepts a tolerance of ±0.001. If the sum deviates further, automatic normalisation is offered (each fraction is divided by the total sum).

**Q: What happens if the density calculation does not converge?**  
A: If the bisection fails to converge within 200 iterations or the density falls outside the physical range, the application displays an error and returns to device type selection. Possible causes: pressure/temperature outside the BWRS domain, mixture with a predominantly heavy component at high pressure.

**Q: Why is the mass flow rate zero?**  
A: The application returns Qm = 0 when one of the ISO 5167 limits is exceeded (D too small/large, β out of range, Re too low). The specific error message is displayed in the console.

**Q: Can I use the application for liquids?**  
A: No. BWRS is an equation of state for the gas phase. The application does not detect condensation and results in the liquid phase or near the dew point are incorrect.

**Q: What does the thermal correction for diameters mean?**  
A: D and d are entered at 20 °C (the standard reference temperature). The application automatically corrects to the operating temperature T using the thermal expansion coefficients of carbon steel (pipe) and stainless steel (orifice plate). If the correction exceeds ±0.10%, a warning is displayed.

**Q: How is the Wobbe index interpreted?**  
A: The Wobbe index (Ws) characterises the calorific power per unit of volumetric flow at constant pressure. Two gases with the same Ws release the same energy at the same Δp through a burner — the interchangeability criterion (ISO 13686).

**Q: How accurate is the BWRS density compared to AGA-8?**  
A: For typical natural gases (CH₄ > 70 %) at industrial conditions (p < 200 bar, T = −10…60 °C), BWRS reproduces density with errors below 0.3–0.5 %. AGA-8 DC92 is more accurate at high pressures. BWRS is preferred in embedded applications due to its lower complexity.

**Q: Does the application run on Linux?**  
A: Yes. The code is portable C++20. Build and installation instructions: `Linux/INSTALL.txt` and `Linux/install.sh`.

**Q: How do I save the composition for future sessions?**  
A: At the end of composition entry, the application asks whether you want to save it. Answer `y`. At the next launch, if `BWRSflow_comp.dat` exists, the application offers to reuse it.

**Q: Can I enter the same composition and calculate for multiple T/p/Δp sets?**  
A: Yes. At the end of a calculation, option `1` returns to entering a new T/p/Δp set with the same device and diameters. The composition and BWRS constants are computed once per session.

**Q: What are normal conditions (Nm³/h) vs. standard conditions (Sm³/h)?**  
A: Nm³/h = volumetric flow referenced to 0 °C / 101.325 kPa. Sm³/h = volumetric flow referenced to 15 °C / 101.325 kPa (European convention) or 60 °F / 14.73 psia (US/AGA-3 convention).

**Q: Why does the warning "Mach number Ma > 0.1–0.2" appear?**  
A: The ISO 5167 equation assumes subsonic compressible flow. Above Ma > 0.2, compressibility effects become significant and the expansibility factor ε may no longer be sufficiently accurate. It is recommended to review the operating conditions.
