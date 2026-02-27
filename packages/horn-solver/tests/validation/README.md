# Validation Test Suite

This suite validates the horn-solver FEM pipeline against known analytical solutions and published data. It answers the question: **are the numbers physically correct?**

## Running the Suite

```bash
# All validation tests (inside Docker)
docker run horn-solver:test pytest packages/horn-solver/tests/validation/ -v

# Just validation tests (with marker)
pytest packages/horn-solver/tests/validation/ -v -m validation

# A single case
pytest packages/horn-solver/tests/validation/test_straight_tube.py -v
```

## Test Cases

### V1: Straight Tube with Robin BC (Tier 1 — exact analytical)

**File:** `test_straight_tube.py`

**Geometry:** Cylinder (radius=50mm, length=500mm) — a straight tube with no flare.

**Physics:** With Dirichlet p=1 at the inlet and Robin BC (∂p/∂n = −jkp) at the outlet, the analytical solution is a pure forward-traveling wave p(z) = exp(−jkz). The pressure magnitude at the outlet is |p(L)| = 1.0 at all frequencies.

**Expected:** SPL = 20·log10(1/20µPa) = 93.98 dB, constant across all frequencies.

**Tolerance:** 0.5 dB — only mesh discretisation error should contribute.

**Why it matters:** If this fails, the solver has a fundamental bug in the Helmholtz equation, boundary conditions, or SPL calculation.

---

### V2: Conical Horn — Webster Equation (Tier 1 — analytical, approximate)

**File:** `test_conical_horn_webster.py`

**Geometry:** Conical frustum (throat=25mm, mouth=100mm, length=500mm, 4:1 expansion).

**Physics:** The Webster horn equation for a conical frustum has an analytical solution in spherical wave coordinates. The FEM solves the full 3D Helmholtz equation; Webster is a 1D plane-wave approximation.

**Expected:** FEM SPL within 3 dB of Webster prediction. SPL should increase with frequency in the passband.

**Tolerance:** 3 dB — Webster is an inherently approximate 1D model compared to full 3D FEM.

**Implementation:** The `webster_conical_spl()` function in the test file computes the analytical reference by solving the 2x2 system for forward/backward wave coefficients with matching boundary conditions.

---

### V3: Community Horn Cross-Validation (Tier 2 — published geometry)

**File:** `test_conical_horn_community.py`

**Geometry:** IJERT Type A air horn (throat=5mm, mouth=28.5mm, length=250mm, 32:1 area expansion). Geometry from Choudhari et al., IJERT Vol.3 Issue 2, 2014.

**Reference data:** Webster equation SPL computed at 20 log-spaced frequencies (200-4000 Hz) for the published geometry. Stored as CSV to decouple from V2's Webster implementation.

**Why this geometry:** It has a high expansion ratio and small throat, exercising different solver behaviour than V2. The IJERT paper includes ANSYS simulation and experimental data for validation context.

**Tolerance:** 6 dB — accounts for 3D vs 1D model mismatch and higher-order mode effects at the small throat.

---

### V4: Exponential Horn — Webster Equation (Tier 1 — analytical, approximate)

**File:** `test_exponential_horn_webster.py`

**Geometry:** Exponential horn (throat=5mm, mouth=22.95mm, length=120mm). Geometry from Rasetshwane & Neely (JASA 2012, PMC3316681).

**Physics:** The Webster horn equation for an exponential horn has an exact analytical solution. The pressure is p(x) = exp(-m·x) · [A·exp(-jk_e·x) + B·exp(jk_e·x)] where k_e = sqrt(k² - m²) and m is the flare rate. Below the cutoff frequency fc ≈ 693 Hz, waves are evanescent.

**Expected:** FEM SPL within 3 dB of Webster prediction. High-pass behavior with cutoff around 693 Hz.

**Tolerance:** 3 dB — Webster is an inherently approximate 1D model compared to full 3D FEM.

**Implementation:** The `webster_exponential_spl()` function solves the 2x2 system with Dirichlet/Robin BCs, matching the pattern from V2.

---

### V5: Hyperbolic Horn — Webster Equation (Tier 1 — analytical, approximate)

**File:** `test_hyperbolic_horn_webster.py`

**Geometry:** Hyperbolic (cosh) horn (throat=12.5mm, mouth=100mm, length=300mm, 8:1 radius ratio). Matches V2's expansion ratio but with hyperbolic flare to isolate profile shape effects.

**Physics:** The Webster horn equation for a cosh horn has an exact analytical solution. The pressure is p(x) = [A·exp(-jk_h·x) + B·exp(jk_h·x)] / cosh(β·x) where k_h = sqrt(k² - β²). Below the cutoff frequency fc ≈ 525 Hz, waves are evanescent.

**Expected:** FEM SPL within 3 dB of Webster prediction. High-pass behavior with cutoff around 525 Hz.

**Tolerance:** 3 dB — Webster is an inherently approximate 1D model compared to full 3D FEM.

**Implementation:** The `webster_hyperbolic_spl()` function solves the 2x2 system with tanh(β·L) terms in the Robin condition.

---

### V6: Exponential Horn Community Cross-Validation (Tier 2 — published geometry)

**File:** `test_exponential_horn_community.py`

**Geometry:** IJERT Type B air horn (throat=5mm, mouth=28.5mm, length=165mm), modeled as exponential profile. Geometry from Choudhari et al., IJERT Vol.3 Issue 2, 2014 — same paper as V3 but different horn and profile.

**Reference data:** Webster equation SPL computed at 20 log-spaced frequencies (200-4000 Hz) for the published geometry. Stored as CSV to decouple from V4's Webster implementation.

**Why this geometry:** It uses the same published source as V3 but a different horn (Type B vs Type A) and different profile (exponential vs conical), providing independent cross-validation of the exponential horn path.

**Tolerance:** 6 dB — accounts for 3D vs 1D model mismatch and higher-order mode effects at the small throat.

## Tolerance Rationale

| Tier | Tolerance | Source of error |
|------|-----------|-----------------|
| Tier 1 exact (V1) | 0.5 dB | Mesh discretisation only |
| Tier 1 approximate (V2, V4, V5) | 3 dB | Webster is 1D approximation of 3D |
| Tier 2 cross-validation (V3, V6) | 6 dB | Different physics, driver models, measurement setup |

## Adding a New Validation Case

1. **Choose a reference:** Identify a horn geometry with published or analytically known frequency response data.

2. **Prepare reference data:** Create a CSV file in `reference_data/` with columns `frequency_hz,spl_db`. Add metadata as header comments:
   ```csv
   # Source: [paper/datasheet URL or citation]
   # Geometry: [type], throat_radius=[m], mouth_radius=[m], length=[m]
   # Measurement conditions: [1W/1m, free-field, etc.]
   # Date digitized: [YYYY-MM-DD]
   frequency_hz,spl_db
   200,85.3
   300,87.1
   ...
   ```

3. **Create the test:** Add a new test file or method that:
   - Generates matching geometry (or loads a STEP file)
   - Runs the solver with appropriate mesh size (lambda/6 at max frequency)
   - Loads the reference data
   - Calls `assert_spl_within_tolerance()` with the appropriate tier tolerance

4. **Mark with `@pytest.mark.validation`** so it can be run separately.

## Reference Data Sources

| File | Case | Source | Status |
|------|------|--------|--------|
| `straight_tube_analytical.json` | V1 | Analytical plane wave solution | Complete |
| `conical_horn_webster.json` | V2 | Webster horn equation | Complete |
| `conical_horn_hornresp.csv` | V3 | Webster SPL for IJERT Type A geometry | Complete |
| `exponential_horn_webster.json` | V4 | Webster horn equation (exponential) | Complete |
| `hyperbolic_horn_webster.json` | V5 | Webster horn equation (hyperbolic) | Complete |
| `exponential_horn_ijert.csv` | V6 | Webster SPL for IJERT Type B geometry | Complete |
