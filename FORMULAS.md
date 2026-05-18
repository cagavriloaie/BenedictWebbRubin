# Formule matematice utilizate în calculatorul BWR

## Cuprins

1. [Constante fizice și conversii](#1-constante-fizice-și-conversii)
2. [Ecuația de stare BWRS](#2-ecuația-de-stare-bwrs)
3. [Reguli de amestecare BWRS](#3-reguli-de-amestecare-bwrs)
4. [Calculul densității (metoda bisecției)](#4-calculul-densității-metoda-bisecției)
5. [Calculul vâscozității dinamice](#5-calculul-vâscozității-dinamice)
6. [Calculul debitului masic – ISO 5167](#6-calculul-debitului-masic--iso-5167)
7. [Coeficienți de debit pentru dispozitive de laminare](#7-coeficienți-de-debit-pentru-dispozitive-de-laminare)
8. [Factorul de expansibilitate](#8-factorul-de-expansibilitate)
9. [Corecții termice pentru diametre](#9-corecții-termice-pentru-diametre)
10. [Limite de validare ISO 5167 / STAS 7347-90](#10-limite-de-validare-iso-5167--stas-7347-90)
11. [Standarde și referințe bibliografice](#11-standarde-și-referințe-bibliografice)

---

## 1. Constante fizice și conversii

| Simbol | Valoare | Unitate | Descriere |
|--------|---------|---------|-----------|
| π | 3.141592653589 | — | Constanta matematică pi |
| R | 0.082057366 | L·atm/(mol·K) | Constanta universală a gazelor (CODATA 2018) |
| T₀ | 273.15 | K | Offset conversie °C → K |
| T_ref | 20.0 | °C | Temperatura de referință |
| k_kPa/atm | 101.325 | kPa/atm | Conversie presiune |
| α_oțel_C | 12.2 × 10⁻⁶ | 1/°C | Coef. de dilatare termică – oțel carbon (conductă) |
| α_inox | 16.5 × 10⁻⁶ | 1/°C | Coef. de dilatare termică – inox (placă de orificiu) |
| k₂kPa | 2000.0 | — | Factor de unități: 2 × (Pa/kPa) pentru SI |

---

## 2. Ecuația de stare BWRS

Ecuația Benedict–Webb–Rubin–Starling (Starling, 1973):

```
p = ρRT
  + (B₀RT − A₀ − C₀/T² + D₀/T³ − E₀/T⁴) · ρ²
  + (bRT − a − d/T) · ρ³
  + α(a + d/T) · ρ⁶
  + (c/T²) · ρ³ · (1 + γρ²) · exp(−γρ²)
```

**Semnificația parametrilor pentru fiecare componentă pură:**

| Parametru | Descriere |
|-----------|-----------|
| A₀ | Al doilea coeficient virial (forțe atractive) |
| B₀ | Al doilea coeficient virial (excludere de volum) |
| C₀ | Corecție termică de ordinul 1/T² |
| D₀ | Corecție termică de ordinul 1/T³ (extensie Starling) |
| E₀ | Corecție termică de ordinul 1/T⁴ (extensie Starling) |
| a | Al treilea coeficient virial (forțe atractive) |
| b | Al treilea coeficient virial (excludere de volum) |
| c | Al treilea coeficient virial (corecție termică) |
| d | Corecție termică de ordinul 1/T în termenul trei |
| α | Amplitudinea termenului ρ⁶ |
| γ | Parametrul gaussian din termenul exponențial |

Baza de date conține 35 de componente gazoase pure (metan, etan, propan, butani, pentani, hexani, heptani, benzen, toluen, H₂, CO, H₂S, He, Ar, N₂, O₂, CO₂, C₂H₄, C₃H₆, NH₃, C₂H₂).

---

## 3. Reguli de amestecare BWRS

### 3.1 Regulă pătratică (parametrii A₀, B₀, C₀, D₀, E₀, γ)

Coeficientul de interacțiune binar:

```
kᵢⱼ = 1 − 8·√(Vᵢ·Vⱼ) / (∛Vᵢ + ∛Vⱼ)³
```

Parametrii de amestec:

```
A₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(A₀ᵢ·A₀ⱼ) · (1 − kᵢⱼ)
B₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(B₀ᵢ·B₀ⱼ)
C₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(C₀ᵢ·C₀ⱼ) · (1 − kᵢⱼ)³
D₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(D₀ᵢ·D₀ⱼ) · (1 − kᵢⱼ)⁴
E₀_mix = Σᵢ Σⱼ xᵢ·xⱼ · √(E₀ᵢ·E₀ⱼ) · (1 − kᵢⱼ)⁵
γ_mix  = Σᵢ Σⱼ xᵢ·xⱼ · √(γᵢ·γⱼ)
```

### 3.2 Regulă cubică (parametrii a, b, c, d, α)

```
a_mix = Σᵢ Σⱼ Σₗ xᵢ·xⱼ·xₗ · ∛( aᵢ·aⱼ·aₗ ) · (1−kᵢⱼ)·(1−kᵢₗ)·(1−kⱼₗ)
b_mix = Σᵢ Σⱼ Σₗ xᵢ·xⱼ·xₗ · ∛( bᵢ·bⱼ·bₗ )
c_mix = Σᵢ Σⱼ Σₗ xᵢ·xⱼ·xₗ · ∛( cᵢ·cⱼ·cₗ ) · (1−kᵢⱼ)·(1−kᵢₗ)·(1−kⱼₗ)
d_mix = Σᵢ Σⱼ Σₗ xᵢ·xⱼ·xₗ · ∛( dᵢ·dⱼ·dₗ ) · (1−kᵢⱼ)·(1−kᵢₗ)·(1−kⱼₗ)
α_mix = Σᵢ Σⱼ Σₗ xᵢ·xⱼ·xₗ · ∛( αᵢ·αⱼ·αₗ )
```

---

## 4. Calculul densității (metoda bisecției)

Se rezolvă ecuația BWRS față de ρ prin metoda bisecției:

```
T = t [°C] + 273.15

ρ_inf = 10⁻⁴ kg/m³   (limita inferioară)
ρ_sup = 5.0  kg/m³   (limita superioară)

Repetă:
  ρ = (ρ_inf + ρ_sup) / 2
  p_calc = BWRS(ρ, T, parametrii_amestec)
  Dacă p_calc > p_dat  →  ρ_sup = ρ
  Altfel               →  ρ_inf = ρ
Până când |p_dat − p_calc| < 5×10⁻⁴ atm
```

---

## 5. Calculul vâscozității dinamice

### 5.1 Vâscozitate gaz diluat – Chapman–Enskog

Pentru fiecare componentă i:

```
η_diluat[i] = √(T/Tcᵢ) · (1 + 0.323·ln(T/Tcᵢ)) / (1 + 0.323·ln(273.15/Tcᵢ)) · etᵢ · √mᵢ
```

Vâscozitatea amestecului în regim diluat:

```
η_diluat_mix = Σᵢ (η_diluat[i] · xᵢ · √mᵢ) / √M_mix

unde M_mix = Σᵢ (xᵢ · mᵢ)   [masa molară a amestecului, g/mol]
```

### 5.2 Corecție pentru densitate ridicată – Lucas / Chung–Lee–Starling

```
η = η_diluat_mix + (10.8×10⁻⁸ / ξ) · [exp(1.439·ρ_red) − exp(−1.111·ρ_red)]^1.358

unde:
  ρ_red = ρ / ρ_crit            (densitate redusă)
  ξ     = Tc^6 / (√M · Pc^(2/3))     (formulare calibrată — nu forma standard Tc^(1/6))
  ρ_crit = Pc / (R · Zc · Tc)   [mol/L]
```

---

## 6. Calculul debitului masic – ISO 5167

### 6.1 Ecuația principală de debit masic

```
Qm = α · ε · (π/4) · d² · √(2000 · Δp · ρ)

unde:
  Qm  – debit masic [kg/s]
  α   – coeficientul de viteză [-]
  ε   – factorul de expansibilitate [-]
  d   – diametrul orificiului la temperatura de măsurare [m]
  Δp  – presiunea diferențială [kPa]
  ρ   – densitatea gazului la condițiile de măsurare [kg/m³]
  2000 = factor de conversie unități SI (√(2·1000))
```

### 6.2 Coeficientul de viteză

```
α = C / √(1 − β⁴)

unde:
  C – coeficientul de debit
  β = d/D – raportul de laminare
```

### 6.3 Numărul Reynolds și iterația convergentă

```
Inițial: Re₀ = 10⁶

La iterația k:
  αₖ  = f(Reₖ)                        (coeficientul de viteză)
  Qmₖ = αₖ · ε · (π/4) · d² · √(2000·Δp·ρ)
  Re_{k+1} = 4·Qmₖ / (π·D·η)

Convergență: |Qmₖ − Qm_{k-1}| / max(Qmₖ, 10⁻¹⁰) < 10⁻⁶  (toleranță relativă)
```

### 6.4 Viteza medie în conductă

```
v = 4·Qm / (π·D²·ρ)   [m/s]
```

### 6.5 Căderea de presiune la dispozitivul de laminare

Diafragme și ajutaje (ISO 5167-2/3 Anexa A):
```
  Δp_laminare = (√(1−β⁴(1−C²)) − C·β²) / (√(1−β⁴(1−C²)) + C·β²) · Δp
```

Tuburi Venturi clasice (fracție empirică din Δp):
```
  brut turnat: ~15%  ·  Δp
  prelucrat:   ~ 8%  ·  Δp
  tablă sudată: ~15% ·  Δp
```

### 6.6 Debit volumic

```
La condițiile de măsurare:    Qv     = 3600 · Qm / ρ       [m³/h]
La condițiile de referință:   Qv_ref = 3600 · Qm / ρ_ref   [m³/h]
```

---

## 7. Coeficienți de debit pentru dispozitive de laminare

### 7.1 Plăci de orificiu – ecuația Reader–Harris/Gallagher (ISO 5167-2)

```
C = 0.5961
  + 0.0261·β²
  − 0.216·β⁸
  + 0.000521·(10⁶·β/Re)^0.7
  + (0.0188 + 0.0063·A)·β^3.5·(10⁶/Re)^0.3
  + (0.043 + 0.080·exp(−10·L₁) − 0.123·exp(−7·L₁))·(1 − 0.11·A)·β⁴/(1 − β⁴)
  − 0.031·(M₂ − 0.8·M₂^1.1)·β^1.3
  + corecție_D_mic     (dacă D < 71.12 mm)

unde:
  A   = (19000·β/Re)^0.8
  M₂  = 2·L₂ₚ / (1 − β)
```

**Poziția prizelor de presiune (L₁, L₂ₚ):**

| Tip priză | L₁ | L₂ₚ |
|-----------|----|-----|
| Unghi (corner taps) | 0 | 0 |
| Flanșă (flange taps) | 0.0254/D | 0.0254/D |
| D și D/2 | 1.0 | 0.50 |

### 7.2 Ajutaje (nozzles) – ISO 5167-3

**Ajutaj ISA 1932:**
```
C = 0.99 − 0.2262·β^4.1
  + (0.000215 − 0.001125·β + 0.00249·β^4.7)·(10⁶/Re)^1.15
```

**Ajutaj cu rază lungă (Long Radius Nozzle):**
```
C = 0.9965 − 0.00653·√β·√(10⁶/Re)
```

### 7.3 Tuburi Venturi – ISO 5167-4

| Tip | C |
|-----|---|
| Venturi clasic (rugos) | 0.984 |
| Venturi clasic (prelucrat mecanic) | 0.995 |
| Venturi clasic (sudat) | 0.985 |
| Venturi-nozzle | 0.9858 − 0.196·β^4.5 |

---

## 8. Factorul de expansibilitate

### 8.1 Plăci de orificiu (ISO 5167-2)

```
ε = 1 − (0.41 + 0.35·β⁴) · Δp / (κ·p)

unde κ = Cp_mix / (Cp_mix − R_SI)   cu  R_SI = 8.314 J/(mol·K)
         Cp_mix = Σᵢ xᵢ · Cp°ᵢ      (capacitate calorică ideală la 20 °C)
```

### 8.2 Ajutaje și tuburi Venturi (ISO 5167-3/4)

```
y = 1 − Δp/p

ε = √[ κ·y^(2/κ) / (κ−1)
      · (1 − β⁴) / (1 − β⁴·y^(2/κ))
      · (1 − y^((κ−1)/κ)) / (1 − y) ]
```

---

## 9. Corecții termice pentru diametre

Diametrele conductei (D) și orificiului (d) se corectează față de temperatura de referință de 20 °C:

```
D_actual = D_ref · (1 + α_exp · (T − 20))
d_actual = d_ref · (1 + α_exp · (T − 20))

unde:
  T        – temperatura de măsurare [°C]
  α_exp    – coeficientul de dilatare termică liniară [1/°C]
             · oțel carbon (conductă): 12.2 × 10⁻⁶ /°C
             · inox (placă orificiu):  16.5 × 10⁻⁶ /°C
```

---

## 10. Limite de validare ISO 5167 / STAS 7347-90

### 10.1 Diametrul conductei D [mm]

| Dispozitiv | D_min | D_max |
|------------|-------|-------|
| Orificiu – unghi / flanșă | 50 | 1000 |
| Orificiu – D și D/2 | 50 | 760 |
| Ajutaj ISA 1932 | 50 | 500 |
| Ajutaj cu rază lungă | 50 | 630 |
| Venturi clasic rugos | 100 | 800 |
| Venturi clasic prelucrat | 50 | 250 |
| Venturi clasic sudat | 200 | 1200 |
| Venturi-nozzle | 65 | 500 |

### 10.2 Raportul de laminare β = d/D

| Dispozitiv | β_min | β_max |
|------------|-------|-------|
| Orificiu (orice priză) | 0.20 | 0.75 |
| Ajutaj ISA 1932 | 0.30 | 0.80 |
| Ajutaj cu rază lungă | 0.20 | 0.80 |
| Venturi clasic rugos | 0.30 | 0.75 |
| Venturi clasic prelucrat | 0.40 | 0.75 |
| Venturi clasic sudat | 0.40 | 0.70 |
| Venturi-nozzle | 0.316 | 0.775 |

### 10.3 Numărul Reynolds Re

| Dispozitiv | Condiție | Re_min | Re_max |
|------------|----------|--------|--------|
| Orificiu unghi (β < 0.45) | — | 5 000 | 10⁸ |
| Orificiu unghi (0.45 ≤ β < 0.77) | — | 10 000 | 10⁸ |
| Orificiu unghi (0.77 ≤ β ≤ 0.80) | — | 20 000 | 10⁸ |
| Orificiu flanșă / D-D/2 | — | 1.26×10⁶·β²·D | 10⁸ |
| Ajutaj ISA (β < 0.44) | — | 70 000 | 10⁷ |
| Ajutaj ISA (0.44 ≤ β ≤ 0.80) | — | 20 000 | 10⁷ |
| Ajutaj cu rază lungă | — | 10 000 | 10⁷ |
| Venturi clasic | — | 200 000 | 2×10⁶ |
| Venturi-nozzle | — | 150 000 | 2×10⁶ |

---

## 11. Standarde și referințe bibliografice

| Referință | Descriere |
|-----------|-----------|
| **ISO 5167-2:2003** | Plăci de orificiu – ecuația Reader–Harris/Gallagher |
| **ISO 5167-3:2003** | Ajutaje și Venturi-nozzle |
| **ISO 5167-4:2003** | Tuburi Venturi clasice |
| **STAS 7347-90** | Standard român echivalent ISO 5167 |
| **Starling K.E. (1973)** | *Fluid Thermodynamic Properties for Light Petroleum Systems* – parametrii BWRS |
| **Nishiumi H., Saito S. (1975)** | Parametrii D₀, E₀ pentru C1–C8 și gaze permanente |
| **Benedict M., Webb G.B., Rubin L.C. (1940)** | Ecuația originală BWR, *J. Chem. Phys.* 8, 334 |
| **Chapman S., Cowling T.G.** | *The Mathematical Theory of Non-Uniform Gases* – vâscozitate gaz diluat |
| **Lucas K. (1980)** | Corecție vâscozitate la densitate ridicată |
| **Chung T.H., Lee L.L., Starling K.E. (1984)** | Corecție vâscozitate la densitate ridicată |
