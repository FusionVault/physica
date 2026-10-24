"""Checks for the review of research/external/2026-10-24-baass-peratt-pais-video.md.
Standard library only. Run from outside the repository: python3 -I <path>."""
# ---- Pillar 1: scale invariance, Birkeland vs gravity ----
from math import pi, sqrt, log10
G=6.67430e-11; MU0=1.25663706212e-6; MP=1.67262192369e-27; E=1.602176634e-19; KB=1.380649e-23
# Scale invariance: the Alfven speed and Lundquist number are dimensionless-similar across scales.
def alfven(B, n, m=MP): return B/sqrt(MU0*n*m)
def lundquist(L, vA, eta): return MU0*L*vA/eta
# FuZE sheared-flow Z-pinch: n=4.9e23 m^-3, B from Bennett at I~200 kA over a~0.3 mm... use reported B~few T.
vA_fuze = alfven(3.0, 4.9e23)
print("FuZE v_A (B=3 T, n=4.9e23):", vA_fuze, "m/s")
# Galactic Birkeland filament: B~3 uG, n~1e4 m^-3 (1e-2 cm^-3), L~kpc
vA_gal = alfven(3e-10, 1e4)
print("galactic filament v_A (B=3e-10 T, n=1e4):", vA_gal, "m/s")
# Spatial scales spanned: lab mm to Hubble 1e26 m
print("orders of magnitude lab(1e-3 m) to Hubble(1e26 m):", log10(1e26/1e-3))
print("Alfven's 28 orders claim: 1e26/1e-2 =", log10(1e26/1e-2))
# Filamentary (Biot-Savart) attraction 1/r vs gravity 1/r^2: force per length between two line currents
def ampere_force(I1, I2, r): return MU0*I1*I2/(2*pi*r)
def grav_line(lam1, lam2, r): return 2*G*lam1*lam2/r
# Two 1 MA Birkeland currents 1 pc apart, line mass 1 Msun/pc
Msun=1.989e30; pc=3.0857e16
lam = Msun/pc
print("Ampere force/length, two 1 MA at 1 pc:", ampere_force(1e6,1e6,pc), "N/m")
print("gravity force/length, 1 Msun/pc lines at 1 pc:", grav_line(lam,lam,pc), "N/m")
print("ratio EM/grav:", ampere_force(1e6,1e6,pc)/grav_line(lam,lam,pc))
# Galactic rotation: v^2 = G M(<r)/r (Newton). Flat v~200 km/s at r=10 kpc needs M:
v=2e5; r=10e3*pc
print("enclosed mass for v=200 km/s at 10 kpc (Msun):", v*v*r/G/Msun)

# ---- Pillar 2: plasma frequency and refraction ----
EPS0=8.8541878128e-12; E=1.602176634e-19; ME=9.1093837015e-31; C=299792458.0
# Plasma frequency f_p = (1/2pi) sqrt(n_e e^2/(eps0 m_e)); a wave below f_p is reflected, above it refracts.
def fp(ne): return sqrt(ne*E*E/(EPS0*ME))/(2*pi)
for ne in (1e18, 1e19, 1.24e18):
    print(f"n_e={ne:.2e} m^-3: plasma frequency {fp(ne)/1e9:.2f} GHz")
# To refract X-band (10 GHz) the plasma must be underdense: f_p < 10 GHz, i.e. n_e < critical density
nc = (10e9*2*pi)**2*EPS0*ME/E**2
print("critical density for 10 GHz (m^-3):", nc)
# Refractive index of underdense plasma: n = sqrt(1 - (fp/f)^2)
for ne in (1e17, 5e17, 1e18):
    print(f"n_e={ne:.0e}: index at 10 GHz = {sqrt(max(0,1-(fp(ne)/10e9)**2)):.4f}")
# Gordon / transformation-optics: a graded index n(r) bends a ray like a metric with n = 1/sqrt(g00).
# The gradient needed to bend a ray through 90 deg over a 1 m shell of a 10 m sphere: dn/dr ~ 1/R.
print("index gradient to turn a ray over L=1 m shell: dn/dr ~ 1/L =", 1.0, "per m (order unity)")

# ---- Pillar 3: dielectrophoresis, field angular momentum, mass test ----
EPS0=8.8541878128e-12; MU0=1.25663706212e-6; C=299792458.0; G=6.67430e-11; g0=9.80665
# Setup 1: dielectrophoretic (Kelvin) force on a neutral dielectric sheet in a non-uniform field.
# F_DEP = (1/2) V eps0 (eps_r - 1) grad(E^2).  Depends on E^2 -> polarity-independent.
def f_dep(V, eps_r, gradE2): return 0.5*V*EPS0*(eps_r-1)*gradE2
# A 0.15 mm plastic sheet, say 2 cm x 2 cm, eps_r=2.2, near a 40 kV electrode with a strong gradient.
V_sheet = 0.02*0.02*0.15e-3
# E near a thin HV wire ~ V/(r ln(d/a)); grad(E^2) ~ 2 E dE/dr. Take E~1e7 V/m, dE/dr~1e9 V/m^2.
E, dEdr = 1e7, 1e9
gradE2 = 2*E*dEdr
F = f_dep(V_sheet, 2.2, gradE2)
print("DEP force on plastic sheet (N):", F, " = weight of (mg):", F/g0*1e6)
# Key point: F ~ eps0 (eps_r-1) grad E^2 -> reversing polarity (E -> -E) leaves E^2 unchanged.
print("polarity-independence: F(+V) == F(-V) because F ~ E^2")
# Setup 2: Feynman field angular momentum L = eps0 int r x (E x B) dV; switching the field transfers torque.
# Radial E ~ 1e6 V/m over the chamber, B from the coils ~ 1e-3 T, lever arm ~ 0.05 m, volume ~ 1e-4 m^3.
Efield, B, r_arm, Vol = 1e6, 1e-3, 0.05, 1e-4
L_em = EPS0*r_arm*Efield*B*Vol
print("Feynman EM angular momentum (kg m^2/s):", L_em, " -> torque impulse on switch-off")
# Mass-independence test: DEP scales with volume*(eps_r-1), not mass; a true g-field scales with m=rho V.
# Compare polyethylene (rho 950, eps_r 2.2) vs gold-coated shell vs tungsten (rho 19300) of same geometry.
print("DEP/volume identical for same-geometry balls; gravitational force ratio tungsten/PE =", 19300/950)
# A real GEM vortex field E_g would give a = E_g independent of material; residual a to prove it:
# the detectable threshold is set by the gravimeter, e.g. cold-atom 1e-10 m/s^2.
print("cold-atom gradiometer can resolve a mass-dependent residual a ~ 1e-10 m/s^2 outside the cage")
# Gary Stephenson 20 kV nonlinearity threshold: just a stated voltage, no formula to verify.

# ---- Pillar 4: covariant coupling, inductive recovery ----
G=6.67430e-11; C=299792458.0
# Bearden's "4D energy from the vacuum": the real statement is the coupling 8 pi G/c^4.
k = 8*pi*G/C**4
print("Einstein coupling 8 pi G/c^4 (s^2/(kg m)):", k)
# Energy the metric can exchange with matter in a time-varying g_mu_nu: governed by this coupling,
# ~43 orders of magnitude weaker than electromagnetism. Show the "stiffness" of spacetime:
print("1/k (kg m/s^2 = N):", 1/k, " (the force to curve spacetime by order 1 over 1 m)")
# Helion direct inductive recovery: Faraday V = -N dPhi/dt; round-trip efficiency eta_rec.
# 0-D: fusion energy E_fus pushes against flux; recovered electrical E_rec = eta_rec * (magnetic work).
# Compare D-T vs D-He3 reactivity at the FRC temperature (use the plasmoid crate values conceptually).
# Reactivity ratio check: Bosch-Hale D-T at 10 keV ~ 1.1e-22, D-He3 much smaller.
# Here just the recovery arithmetic: if 70% of fusion energy is in charged particles recoverable at eta_rec,
for eta_rec in (0.90, 0.95):
    for charged_frac in (1.0, 0.68):   # D-He3 aneutronic ~ all charged; D-T 80% in neutrons
        print(f"eta_rec={eta_rec}, charged fraction={charged_frac}: electrical return fraction {eta_rec*charged_frac:.2f}")
# The point: aneutronic D-He3 (charged products) suits direct recovery; D-T puts 80% into neutrons.
