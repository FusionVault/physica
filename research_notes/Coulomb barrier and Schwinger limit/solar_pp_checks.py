import math
alpha=7.2973525693e-3; kB=8.617333262e-8; mp=938272.08816; mu=mp/2
T=1.55e7; kT=kB*T; T9=T/1e9
E0=kT*(math.pi*alpha/math.sqrt(2))**(2/3)*(mu/kT)**(1/3); D=4*kT*math.sqrt(E0/(3*kT))
S0=4.01e-22; S1=4.49e-24  # keV b, b  (SF-II Table I)
Seff=S0*(1+5*kT/(36*E0))+S1*E0*(1+35*kT/(36*E0))   # keV b
f0=1.049  # computed Salpeter value (assumed composition)
A=0.5
sv=1.301e-14*(1/A)**(1/3)*f0*(Seff/1000)*T9**(-2/3)*math.exp(-3*E0/kT)
# exact form sqrt(2/(mu kT)) dE0/kT f0 Seff exp(-3E0/kT), c in cm/s, barn=1e-24 cm2
c=2.99792458e10
sv2=math.sqrt(2/(mu*kT))*c*(D/kT)*f0*(Seff*1e-24)*math.exp(-3*E0/kT)
print(f"E0={E0:.4f} keV D={D:.4f} keV Seff={Seff:.4e} keV b")
print(f"<sv>_pp (1.301e-14 form) = {sv:.4e} cm3/s ; exact-prefactor form = {sv2:.4e} cm3/s")
n=150*0.34/1.66053906660e-24
tau=1/(n*sv)/3.15576e7
print(f"n_p={n:.3e} cm^-3 ; proton lifetime vs pp = {tau:.3e} yr")
