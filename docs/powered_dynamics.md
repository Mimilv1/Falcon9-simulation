# Powered integration with moving CG

The implementation uses open-control-volume linear/angular momentum balances,
specialized to axial CG migration and a balanced axial exhaust centroid. The
[RocketPy v1 explicit linear and angular equations](https://docs.rocketpy.org/en/latest/technical/equations_of_motion_v1.html)
provide the reference-point formulation. The code below uses current-CG variables
and remaps the longitudinal axis to Body +x. No RocketPy implementation code was copied.

## Definitions

- m: current vehicle mass; J: inertia about current CG, in Body components.
- c: CG station from the structural base; ċ and c̈: migration through the structure.
- qᵢ ≥ 0: port outflow; μᵢ = −qᵢ and μ̇ᵢ = −dqᵢ/dt.
- dᵢ: nozzle-exit station minus current CG station.
- Gᵢ: local exit gyration per unit mass about its own centroid, in m².
- S(d) = (d·d)I − ddᵀ, the per-unit-mass parallel-axis tensor.
- F: net nongravitational force including thrust reaction; M: net CG moment.
- g, ω, a: Body gravity, angular velocity and inertial/geometric CG acceleration.

The CG-form balances implemented in `dynamics.accelerations` are:

**m a = F + m g + Σ[2μᵢ(ω × dᵢ − ċ) + μ̇ᵢdᵢ]**

**Jω̇ = M − ω × (Jω) − J̇ω + Σ μᵢ[S(dᵢ) + Gᵢ]ω**

The angular outflow tensor is consistent with the nozzle angular-momentum flux
in the [RocketPy control-volume derivation](https://docs.rocketpy.org/en/latest/technical/equations_of_motion.html).
The explicit reference-point terms can also be inspected in
[Flight.u_dot_generalized](https://github.com/RocketPy-Team/RocketPy/blob/master/rocketpy/simulation/flight.py).
References were consulted on 2026-10-27.

These CG equations are used under the following restricted closure: c/ċ/c̈
are axial, internal relative-flow angular momentum is omitted, exit stations
are fixed in the structure, and the transverse first moment of outflow and its
derivative vanish. Port validation rejects an unbalanced transverse centroid.
Finite exit gyration is supplied by B; uniform disk flux is an explicit A.
This does not resolve tank/pipe transients or a moving nozzle boundary.

Thrust already contains relative exhaust momentum and pressure reaction. The
additional terms are convective/migration and spin outflow terms; adding q times
exhaust speed as a second thrust would count propulsion twice.

## Body-only state derivatives

Let rᴮ be the vector from a fixed initial origin to CG, expressed in CURRENT
Body axes. Let vᴮ be CG velocity in the same axes. Then:

- ṙᴮ = vᴮ − ω × rᴮ.
- v̇ᴮ = aᴮ − ω × vᴮ.
- q̇ = ½ q ⊗ [0,ω]; quaternion q maps current Body to initial Body.
- Each stage fuel derivative is minus its consumption rate.

The initial Body basis is a mathematical reference with no specified geographic
orientation. Initial q=[1,0,0,0]. ENU/ECEF conversion is deferred.

## Geometric derivatives and sensors

For a settled column of area A and density ρ: h=m/(ρA), ḣ=ṁ/(ρA),
ḧ=m̈/(ρA), and liquid center b+h/2. Combined CG derivatives follow
differentiation of the mass-weighted first moment, including mass second
derivatives. These are geometric model derivatives, not independent thrusts.

For a sensor at fixed structural station s, define l=s−c. Its acceleration is:

**a_sensor = a_CG + ω̇ × l + ω × (ω × l) − 2ω × ċ − c̈**

All terms here are current Body components. Ideal accelerometer specific force
is a_sensor−g. The default station in acceleration_truth is coincident with the
current CG at that instant; C should explicitly supply a fixed s for an IMU.

## Smooth transitions and fuel limits

Mass flow and supplied dq/dt must be consistent within each integration segment.
The μ̇d term matters during a smooth ramp. An abrupt jump in mass flow makes
this term impulsive and requires an explicit impulse/state map; simply switching
q while setting dq/dt=0 is not this smooth model. `propagate` splits scheduled
breakpoints and shortens failed fuel-boundary trial steps. No mass is clipped
to a fabricated post-depletion state.

## Verification and scope

Checks include the independent Tsiolkovsky velocity/position solution; equality
with a separately assembled fixed-reference 6×6 coupled balance; geometric CG
finite differences; roll angular momentum leaving with a finite exit; accepted
base/CG state mapping; fixed-sensor kinematics; and limiting-tank stopping.
A separate local 30 s prescribed-load case compares total linear/angular
momentum budgets and halved integration steps. It is not flight-data validation.

Runnable project configurations now include early ascent, attitude verification
and [smooth shutdown/separation](shutdown_separation.md). Generated plots and
outputs remain outside version control. The first B models are research
surrogates; measured transients, separate-stage aerodynamics, fluid settling/slosh,
variable O/F, navigation, guidance/control, geographic outputs and full mission
validation remain separate.
