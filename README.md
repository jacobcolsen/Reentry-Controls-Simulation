# Reentry Simulation — Engineering Reference

**Version:** V69 (2026-03-03)
**Files:** `LAUNCH_reentry_simulation.html` · `reentry_trajectory_viz_integrated.html`

---

## Overview

A 6-degree-of-freedom (6-DOF) atmospheric reentry trajectory simulator built in browser-native JavaScript/HTML. The physics engine integrates the full equations of motion in spherical coordinates using a 5th-order Runge-Kutta integrator (Cash-Karp RK45). Outputs include full trajectory state histories, aerodynamic heating profiles, and deceleration loads. A companion 3D visualization tool renders the trajectory on a globe with a real-time HUD, analog gauges, and selectable vehicle models.

Apollo 10 historical data (entry conditions + bank angle time history) is bundled for validation.

---

## Table of Contents

1. [Equations of Motion](#1-equations-of-motion)
2. [Atmospheric Model](#2-atmospheric-model)
3. [Gravity Model](#3-gravity-model)
4. [Aerodynamic Heating Model](#4-aerodynamic-heating-model)
5. [Deceleration Model](#5-deceleration-model)
6. [Numerical Integration](#6-numerical-integration)
7. [Vehicle Aerodynamic Parameters](#7-vehicle-aerodynamic-parameters)
8. [How to Use the Physics Simulator](#8-how-to-use-the-physics-simulator)
9. [How to Use the 3D Visualization Tool](#9-how-to-use-the-3d-visualization-tool)

---

## 1. Equations of Motion

The simulation uses a **3D point-mass trajectory model in spherical coordinates** with lift, drag, and gravity. This is sometimes called the "flat-rotating-Earth" approximation's spherical analog — no Coriolis or centrifugal terms are included, but the full spherical geometry is preserved.

### State Vector

```
x = [ r, θ, φ, V, γ, ψ ]
```

| Symbol | Name | Units |
|--------|------|-------|
| r | Radial distance from Earth's center | km |
| θ | Longitude | rad |
| φ | Latitude | rad |
| V | Velocity magnitude (relative to atmosphere) | km/s |
| γ | Flight path angle (positive = climbing) | rad |
| ψ | Heading angle (measured from North, clockwise) | rad |

### Control Input

| Symbol | Name | Units |
|--------|------|-------|
| σ | Bank angle (roll about velocity vector) | rad |

The bank angle σ is the primary guidance input. It rotates the lift vector about the velocity axis, splitting lift between the vertical (out-of-plane pull) and horizontal (cross-range steering) components.

### The Six Differential Equations

**1. Radial rate (altitude change):**
```
dr/dt = V · sin(γ)
```
Straightforward kinematic relationship — the radial velocity is the component of V along the local vertical.

**2. Longitude rate:**
```
dθ/dt = V · cos(γ) · cos(ψ) / (r · cos(φ))
```
The northward component of horizontal velocity `V·cos(γ)·cos(ψ)` maps to a longitude rate via the spherical metric. The `cos(φ)` in the denominator accounts for meridian convergence — longitude degrees are shorter at higher latitudes.

**3. Latitude rate:**
```
dφ/dt = V · cos(γ) · sin(ψ) / r
```
The eastward component of horizontal velocity maps directly to latitude rate via the spherical arc length.

**4. Velocity (speed) rate:**
```
dV/dt = -(1/2) · ρ · C_D · S · V² / m  -  g · sin(γ)
```
The two deceleration terms are:
- **Aerodynamic drag:** `D = (1/2)·ρ·C_D·S·V²` acting opposite to the velocity vector.
- **Gravity component along velocity:** `g·sin(γ)` — positive γ means climbing, so gravity decelerates; negative γ (descending) means gravity accelerates the vehicle.

There is no thrust term — this is a purely ballistic/glide descent.

**5. Flight path angle rate:**
```
dγ/dt = (1/V) · [ (1/2)·ρ·C_L·S·V²/m · cos(σ)  +  (V²/r - g) · cos(γ) ]
```
Three terms contribute to changing γ:
- **Lift (vertical component):** `(1/2)·ρ·C_L·S·V²/m · cos(σ)` — the component of lift force that rotates the velocity vector toward the vertical. At zero bank angle (σ = 0) all lift is vertical. At σ = 90° no lift acts vertically (all cross-range).
- **Centrifugal relief:** `V²/r · cos(γ)` — at orbital and near-orbital speeds the vehicle "outfalls" the Earth's curvature; this acts like reduced gravity on γ.
- **Gravity (perpendicular to velocity):** `-g · cos(γ)` — at shallow angles this is the dominant force pulling γ negative (diving into the atmosphere).

**6. Heading angle rate:**
```
dψ/dt = (1/V) · [ (1/2)·ρ·C_L·S·V²/m · sin(σ) / cos(γ)  -  V²/r · cos(γ) · cos(ψ) · tan(φ) ]
```
Two terms:
- **Lift (lateral/cross-range component):** `(1/2)·ρ·C_L·S·V²/m · sin(σ) / cos(γ)` — the lateral lift drives cross-range steering. This is the primary mechanism for range control and skip-entry management.
- **Spherical Earth heading drift:** `-V²/r · cos(γ) · cos(ψ) · tan(φ)` — a purely geometric term that arises because lines of constant heading diverge on a sphere. It is equivalent to a course correction needed to maintain a great-circle path.

---

## 2. Atmospheric Model

### Model Type
Exponential (isothermal) density profile. This is the standard analytic model used in reentry trajectory analysis when a full NRLMSISE-00 table is not required.

### Equation
```
ρ(h) = ρ₀ · exp(−β · h)
```

Where altitude `h = r − r_E` (km), and:

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Sea-level density | ρ₀ | 1.225 | kg/m³ |
| Scale height coefficient | β | 0.14 | km⁻¹ |
| Earth radius | r_E | 6378.137 | km |

The **scale height** H = 1/β ≈ 7.14 km. This means density falls by a factor of e ≈ 2.718 for every 7.14 km of altitude gain, broadly consistent with the real atmosphere between ~10–80 km.

### Notes
- The model is isothermal (constant temperature assumed). Real atmosphere has temperature structure (troposphere, stratosphere, mesosphere) that causes the scale height to vary with altitude. For reentry heating, the critical altitude band is 40–80 km where most deceleration occurs — the exponential model is a reasonable fit there.
- No winds are modeled; atmosphere is assumed stationary relative to the rotating Earth frame (i.e., the "atmosphere-relative" and "inertial" velocity are treated as equal, which is adequate for the altitudes involved).
- Above the Kármán line (~100 km) the atmosphere remains in the model but densities become negligible (ρ < 10⁻⁶ kg/m³), so aerodynamic forces correctly approach zero.

---

## 3. Gravity Model

### Model Type
Inverse-square point-mass gravity. No J2 or higher-order geopotential terms are included.

### Equation
```
g(r) = g₀ · (r₀/r)²
```

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Reference gravity | g₀ | 0.00981 | km/s² (= 9.81 m/s²) |
| Reference radius | r₀ | 6500 | km |

The reference radius r₀ = 6500 km corresponds to an altitude of approximately 121.9 km above Earth's surface — close to the nominal entry interface altitude used for Apollo. Setting the reference there (rather than sea level) slightly improves accuracy in the upper atmosphere where the simulation begins.

### Notes
- At sea level (r = 6378 km): g ≈ 9.78 m/s² (within 0.3% of standard gravity — acceptable)
- At entry interface (r = 6500 km): g ≈ 9.81 m/s² exactly
- The model does not account for Earth's oblateness (J2 ≈ 1.08×10⁻³). For reentry trajectories of 500–2000 s duration, J2 effects on the trajectory are small (sub-kilometer level).

---

## 4. Aerodynamic Heating Model

### Model Type
Empirical stagnation-point convective heat flux, based on the Detra-Kemp-Riddell (DKR) correlation adapted for entry vehicles.

### Heat Flux Equation
```
q̇ = C · (ρ/ρ_ref)^0.5 · (V/V_ref)^3.15    [W/m²]
```

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Scaling constant | C | 1.9983 × 10⁸ | W/m² |
| Reference density | ρ_ref | 1.225 | kg/m³ |
| Reference velocity | V_ref | 7.905 | km/s |

The reference velocity V_ref ≈ 7.905 km/s is the circular orbital velocity at Earth's surface — a natural normalization for entry vehicle speeds.

### Why These Exponents?

The theoretical basis for heating correlations is the **Fay-Riddell stagnation point heating** relation:

```
q̇_stag ∝ ρ^n · V^m
```

- The **density exponent (0.5)** comes from boundary layer heat transfer theory: heat flux scales with the square root of density because the convective heating rate depends on the density in the shock layer, and the boundary layer thickness scales with 1/√ρ.
- The **velocity exponent (3.15)** reflects that most of the energy of a hypervelocity flow ends up as enthalpy rise across the shock. The classical DKR result gives approximately V³, with the 0.15 above 3.0 coming from dissociation and real-gas effects at the speeds and altitudes typical of entry.

### Stagnation Temperature
```
T_stag = T_ref · [1 + (γ-1)/2 · M²]
```

Computed using isentropic relations at the stagnation point:
- T_ref = 288 K (standard sea-level temperature)
- γ = 1.4 (ratio of specific heats for air — appropriate for lower altitudes; real gas effects at high altitude would give γ → 1.2, but this is a first-order approximation)
- M = V / a_ref where a_ref = √(γ · R · T_ref) ≈ 340 m/s (standard speed of sound at sea level)

This gives the theoretical maximum temperature the stagnation point gas can reach if all kinetic energy is converted to enthalpy. Actual heat shield surface temperatures are lower by an emissivity/radiation balance factor.

### Peak Heating Reference Values (Apollo 10 conditions)
At peak heating (approximately t ≈ 80–100 s post-entry interface):
- V ≈ 9–10 km/s
- h ≈ 55–65 km
- q̇_peak ≈ 300–600 W/cm² (3–6 MW/m²)

---

## 5. Deceleration Model

### Formulation
Deceleration in units of Earth g's is derived directly from the velocity rate equation:

```
n_g = |dV/dt| / g₀
```

Where `dV/dt` is the total scalar deceleration (km/s²) from Equation 4 above, and g₀ = 0.00981 km/s².

This gives the **felt load factor** on the crew/payload — the net deceleration they would experience relative to free fall. At peak deceleration, the dominant term is aerodynamic drag (the gravity term along the velocity vector is small at the shallow FPA angles of a lifting entry).

### Dynamic Pressure
Dynamic pressure is a useful intermediate quantity:

```
q̄ = (1/2) · ρ · V²    [Pa]
```

Peak dynamic pressure (max-q) occurs when the rate of change of q̄ is zero, which happens when the rate of density increase (falling through the atmosphere) balances the rate of velocity decrease. For Apollo-class entries this typically occurs around 40–50 km altitude.

### Aerodynamic Force Relations

From the state equations:

```
Drag:    D = (1/2) · ρ · C_D · S · V²    [N, converted from km units]
Lift:    L = (1/2) · ρ · C_L · S · V²    [N]
L/D    = C_L / C_D
```

The **ballistic coefficient** β_BC is a key derived parameter:

```
β_BC = m / (C_D · S)    [kg/m²]
```

Higher β_BC means the vehicle decelerates at lower altitude (and higher dynamic pressure), which increases peak heating and peak g-load. Lower β_BC (blunt, draggy vehicles) decelerate higher and slower.

---

## 6. Numerical Integration

### Method: Runge-Kutta 4/5 (Cash-Karp RK45)

The simulation uses a **6-stage, 5th-order** Runge-Kutta method. The Cash-Karp coefficients provide a 4th-order solution and an embedded 5th-order error estimate.

```
k₁ = f(t, x)
k₂ = f(t + c₂·h,  x + h·(a₂₁·k₁))
k₃ = f(t + c₃·h,  x + h·(a₃₁·k₁ + a₃₂·k₂))
k₄ = f(t + c₄·h,  x + h·(a₄₁·k₁ + ...))
k₅ = f(t + c₅·h,  x + h·(a₅₁·k₁ + ...))
k₆ = f(t + c₆·h,  x + h·(a₆₁·k₁ + ...))

x_{n+1} = xₙ + h · Σ bᵢ · kᵢ    (5th order)
```

| Parameter | Value |
|-----------|-------|
| Time step dt | 0.1 seconds (fixed) |
| Max simulation time | 2500 seconds |
| Termination: ground impact | h < 0 km |
| Termination: velocity | V ≤ 0 km/s |
| Termination: numerical failure | NaN detection |

The fixed 0.1 s timestep is conservative — at entry velocities of 7–12 km/s the vehicle travels ~0.7–1.2 km per step, resolving all physically significant dynamics. The RK45 method provides good energy conservation and phase accuracy for this class of problem.

---

## 7. Vehicle Aerodynamic Parameters

Vehicles are selected automatically based on the L/D ratio entered:

| Vehicle | L/D Range | Default C_D | Default S (m²) | Default m (kg) |
|---------|-----------|-------------|-----------------|----------------|
| School Bus (ballistic) | 0.00–0.05 | 1.5 | 15 | 5000 |
| Stardust SRC | 0.05–0.175 | 1.5 | 15 | 5000 |
| Mercury | 0.175–0.30 | 1.5 | 15 | 5000 |
| Apollo CM | 0.30–0.525 | 1.5 | 15 | 5498.2 |
| Dream Chaser | 0.525–0.75 | 1.5 | 15 | 5000 |
| X-37B | 0.75–0.925 | 1.5 | 15 | 5000 |
| Shuttle Orbiter | 0.925+ | 1.5 | 15 | 5000 |

**Apollo 10 Historical Mode** (hardcoded preset):

| Parameter | Value |
|-----------|-------|
| Mass | 5498.2 kg |
| Entry velocity | 11.06715 km/s |
| Flight path angle | −6.62° |
| Entry longitude | 174.24°E |
| Entry latitude | 23.51°S |
| Entry altitude | 120.125 km |
| Bank angle source | `bank_angle_time_history_0_to_2000s.csv` |

---

## 8. How to Use the Physics Simulator

**File:** `LAUNCH_reentry_simulation.html` — open directly in a browser (no server needed).

### Input Parameters

| Control | Description | Range | Default |
|---------|-------------|-------|---------|
| Vehicle Mass | Total entry mass | 2000–10000 kg | 5000 kg |
| Entry Velocity | Inertial speed at entry interface | 7–12 km/s | 9.0 km/s |
| Flight Path Angle | FPA at entry interface (negative = descending) | −20° to +5° | −5.0° |
| L/D Ratio | Lift-to-drag ratio | 0–1 | 0.30 |
| Drag Coefficient C_D | Aerodynamic drag coefficient | 0.5–2.5 | 1.50 |
| Reference Area S | Aerodynamic reference area | 5–50 m² | 15 m² |

### Bank Angle Control

Four input methods are available — select one at a time:

**1. Slider (constant bank)**
- Single slider sets a constant bank angle for the entire trajectory.
- Range: 0°–180°.
- Use this for quick sensitivity checks.

**2. Draw Profile**
- Click and drag on the bank angle canvas to draw a time-history profile.
- X-axis = time (0 to max simulation time).
- Y-axis = bank angle (0°–180°).
- The simulation linearly interpolates between drawn points.

**3. Waypoints**
- Enter time-angle pairs manually (time in seconds, angle in degrees).
- Add multiple waypoints to define a piecewise-linear profile.
- The simulation linearly interpolates between consecutive waypoints.

**4. CSV Import**
- Upload a two-column CSV: `time_seconds, bank_angle_degrees`.
- No header required.
- The bundled `bank_angle_time_history_0_to_2000s.csv` is the Apollo 10 historical profile.

### Running a Simulation

1. Set parameters and bank angle profile.
2. Click **Run Simulation**.
3. Simulation runs instantly (no real-time waiting).
4. Results display as time-history plots:
   - Altitude vs. time
   - Velocity vs. time
   - Flight path angle vs. time
   - Deceleration (g-load) vs. time
   - Heat flux vs. time
   - Latitude/Longitude ground track
5. Click **Export Trajectory CSV** to save the full state vector for import into the 3D visualizer.

### Apollo 10 Historical Mode

Click the **"Load Apollo 10"** preset button to automatically populate all fields with historical entry conditions and load the bank angle CSV profile. Run the simulation and compare outputs against historical records as a validation baseline.

---

## 9. How to Use the 3D Visualization Tool

**File:** `reentry_trajectory_viz_integrated.html` — open in a browser. Requires a trajectory CSV generated by the physics simulator (or import a pre-existing one).

### Loading a Trajectory

1. Click **"Import Trajectory"** (or drag-and-drop a CSV onto the window).
2. The tool parses the state vector columns and renders the trajectory arc over the globe.
3. Press **Play** (or `P`) to begin playback.

---

### Playback Controls

| Control | Action |
|---------|--------|
| `P` or `SPACE` | Play / Pause |
| `R` | Reset to t = 0 |
| `←` / `→` arrow keys | Step backward / forward one frame |
| Time warp buttons | Set playback speed: **1×, 5×, 10×, 25×, 50×** |

---

### Camera Modes

Cycle with `C` or click the camera mode button:

| Mode | Description |
|------|-------------|
| **Orbit** | Free rotation — click and drag to rotate globe, scroll to zoom |
| **Chase** | Camera follows the vehicle from behind at fixed offset |
| **Nadir** | Camera points straight down from the vehicle (ground track view) |
| **Inertial** | Fixed inertial frame — Earth rotates beneath the camera |

**Mouse Controls (Orbit mode):**
| Action | Control |
|--------|---------|
| Rotate globe | Left-click + drag |
| Zoom in/out | Scroll wheel |
| Pan | Right-click + drag |

---

### HUD and Display Overlays

| Key | Toggle |
|-----|--------|
| `H` | **HUD** — fighter-jet style flight data overlay (velocity, altitude, FPA, heading, g-load, time) |
| `T` | **Thermal view** — heat flux color mapping on vehicle surface and plasma glow |
| `I` | **Info panel** — text telemetry readout (top-left corner), collapsible |

---

### Vehicle Viewer Pane

Located bottom-left (450×320 px). Shows the selected vehicle in 3D with real-time effects:

| Effect | Trigger |
|--------|---------|
| Plasma glow | Activates when heat flux exceeds threshold |
| Bow shock visualization | Appears at hypersonic speeds (M > ~5) |
| Heat shield coloration | Color maps to stagnation temperature |
| Ablation particle effects | Activates at peak heating |

Click inside the vehicle viewer to rotate the model independently of the globe camera.

---

### Gauges

Two analog dial gauges are displayed on-screen:

| Gauge | Range | Units |
|-------|-------|-------|
| **Thermal** | 0 – peak heat flux | W/cm² |
| **Deceleration** | 0 – 40+ | g's |

Both update in real-time during playback.

---

### 3D Visualization Keyboard Reference (Quick Card)

```
P / SPACE     Play / Pause
R             Reset
C             Cycle camera modes
H             Toggle HUD overlay
T             Toggle thermal view
I             Toggle info panel
← →           Step one frame
Mouse drag    Rotate (Orbit mode)
Scroll        Zoom
```

---

## File Structure

```
MOAB_Reentry_Viewer_V69/
├── LAUNCH_reentry_simulation.html          # Physics engine + UI
└── reentry_trajectory_viz_integrated.html  # 3D visualization
```

---

## References

- Vinh, N. X., Busemann, A., & Culp, R. D. — *Hypersonic and Planetary Entry Flight Mechanics*, University of Michigan Press, 1980. (EOM derivation)
- Detra, R. W., Kemp, N. H., & Riddell, F. R. — "Addendum to Heat Transfer to Satellite Vehicles Re-entering the Atmosphere," *Jet Propulsion*, 1957. (Heating correlation)
- U.S. Standard Atmosphere, 1976 — NOAA/NASA/USAF. (Atmospheric reference)
- Cash, J. R. & Karp, A. H. — "A Variable Order Runge-Kutta Method for Initial Value Problems with Rapidly Varying Right-Hand Sides," *ACM Transactions on Mathematical Software*, 1990. (Integrator)
- Apollo 10 Mission Report — NASA MSC-00126, 1969. (Historical validation data)
