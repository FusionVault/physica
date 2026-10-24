# Reproduce published HFGW numbers (Li-Baker, JASON, Baker, Pais). SI unless noted.
from math import pi, sqrt
G=6.67430e-11; c=2.99792458e8; mu0=1.25663706e-6; eps0=8.8541878e-12; hbar=1.054571817e-34; h_pl=2*pi*hbar
Z0=mu0*c; kB=1.380649e-23; sigSB=5.670374e-8; Msun=1.98847e30
def p(label, val, pub=None):
    s=f"{label:62s} {val:11.3e}"
    if pub is not None: s+=f"   published {pub:9.3e}   ratio ours/pub {val/pub:7.3f}"
    print(s)
print("== Li-Baker EPJC 56 (2008) / PRD 80 064013 (2009) ==")
P,W0,nu=10.0,0.05,5e9; w=2*pi*nu; k=w/c; Eph=hbar*w
p("GB peak field psi0 = sqrt(4 Z0 P/(pi W0^2)) [V/m]", sqrt(4*Z0*P/(pi*W0**2)), 1.26e3)
A=1e-30
F_li =c**3*w**2*A**2/(8*pi*G); F_std=c**3*w**2*A**2/(32*pi*G)
p("GW flux Li c^3 w^2 A^2/(8 pi G) [W/m^2]", F_li)
p("GW flux standard, one polarisation c^3w^2A^2/(32piG)", F_std)
p("Graviton flux thru pi W0^2 (Li formula) [1/s]", F_li*pi*W0**2/Eph, 3.77e16)
p("PRD09 eq9 U_gw thru 0.01 m^2 [W]", F_li*0.01, 1.6e-7)
p("PRD09 eq10 N_g [1/s]", F_li*0.01/Eph, 4.8e16)
B,L,ds=10.0,10.0,0.01
for kk,lab in ((100.0,'k=100 as paper'),(k,'k=2pi/lambda')):
    U2=(A*B*kk*L)**2*c*ds/mu0
    p(f"PRD09 eq7 U_em^(2)=(A B k L)^2 c ds/mu0 [W] ({lab})", U2, 2.3e-40)
U2=(A*B*100*L)**2*c*ds/mu0
p("PRD09 eq8 N^(2)=U2/hbar w [1/s] (paper uses h nu=3.3e-24 J)", U2/3.3e-24, 7.0e-17)
p("PRD09 eq11 conversion U2/Ugw", U2/(F_li*ds), 1.4e-33)
p("Standard Gertsenshtein P=4piG B^2 L^2/(mu0 c^4), B=10T,L=10m", 4*pi*G*B**2*L**2/(mu0*c**4))
E0=sqrt(2*Z0*P/ds); p("PRD09 eq14 plane-wave E0 [V/m]", E0, 8.7e2)
N0=P/3.3e-24; p("PRD09 eq15 N0 [1/s]", N0, 3.0e24)
E1=A*B*100*c*L; U1=E1*E0*ds/Z0
p("PRD09 eq17 U1=(1/mu0 c)E1 E0 ds [W]", U1, 6.9e-20)
N1=U1/3.3e-24; p("PRD09 eq18 N1 [1/s]", N1, 2.1e4)
p("  cross term 2 sqrt(N0 N2) (JASON eq 3-32 form) [1/s]", 2*sqrt(N0*U2/3.3e-24))
N1_26=N1*1e4; p("PRD09 eq23 dt>N0/N1^2 at h=1e-26 [s]", N0/N1_26**2, 6.8e7)
print("-- Gaussian-beam scheme, B=3T, l=6m, ds=0.01 m^2, 5 GHz, h=1e-30 --")
B,L=3.0,6.0
U2=(A*B*k*L)**2*c*ds/mu0; Ngw=U2/Eph
p("2nd-order (pure inverse-G) photon rate N_GW, Li field formula", Ngw)
Ust=(0.5*A*B*w*L)**2/(2*Z0)*ds; p("  same, standard E=(1/2)hB w L, <S>=E^2/2Z0", Ust/Eph)
p("cross term 2 sqrt(N_BPF(x=3.5cm)=1.24e22 * N_GW) [1/s]", 2*sqrt(1.24e22*Ngw), 3.54e2)
p("Li dt_min = N_other/N1^2, N_other=2.1e9, N1=8.2e2 [s]", 2.1e9/8.2e2**2, 3.0e3)
p("Li dt_min with N_other=2.1e11 [s]", 2.1e11/8.2e2**2, 3.0e5)
p("JASON quantum limit t > 1/N_GW (Li field formula) [s]", 1/Ngw)
p("JASON limit t > 1/N_GW (standard field) [s]", Eph/Ust)
p("  ratio JASON-limit / Li dt_min(3e3 s)", (1/Ngw)/3.0e3)
print("== JASON JSR-08-506 (Gaussian units where stated) ==")
Gc=6.674e-8; cc=2.998e10
p("eq3-14 F=4piG B^2 L^2/c^4, B=1e5 G, L^2=1e3 cm^2", 4*pi*Gc*1e10*1e3/cc**4, 1e-35)
p("  standard F=G B_G^2 L^2/c^4 (= 4piG B^2L^2/(mu0 c^4) SI)", Gc*1e10*1e3/cc**4)
p("  SI check 4piG B^2 L^2/(mu0 c^4), 10 T, 0.1 m^2", 4*pi*G*100*0.1/(mu0*c**4))
p("eq3-16 U = P L/c, P=10 kW, L=30 cm [J]", 1e4*0.3/c, 1e-5)
Umax=1e6*2e3/1e10; p("eq3-18 Umax = Hdot Q/omega, 1 MW, Q=2e3, w=1e10 [J]", Umax, 2e-1)
PGW=4*pi*Gc*1e10*(Umax*1e7)*30/cc**3*1e-7; p("eq3-15/19 P_GW=4piG B^2 U L/c^3, L=30cm [W]", PGW, 2e-27)
p("  gravitons/s at hbar*w, w=1e10 rad/s", PGW/(hbar*1e10), 2e-2)
E=50e6/1e2/299.79; p("eq3-21 u=E^2/8pi for E=50 MV/m [erg/cm^3]", E**2/(8*pi), 1e7)
p("eq3-24 P_GW with U=3e6 J, L=30 cm [W]", 4*pi*Gc*1e10*3e13*30/cc**3*1e-7, 1e-20)
p("eq4-42 E_GW = h B w L/c (Gaussian) -> V/m, h=3e-32,B=1e5 G,wL/c=10", 3e-32*1e5*10*2.998e4, 1e-22)
p("  ratio 2 E0/E_GW with E0=1.17e3 V/m (JASON: 1e25)", 2*1.17e3/1e-22, 1e25)
p("sec4.2 2 K graviton bath sigma T^4 x 0.01 m^2 [W]", sigSB*2**4*0.01, 1e-8)
p("  conversion 1e-38 -> photon power [W]; wait time at 4kT(2K)", 1e-46)
p("  wait time 1/((1e-46 W)/(4 kB 2K)) [s]", 4*kB*2/1e-46, 1e24)
p("sec4.4 Stephenson U=1e5 W x 1e3 s x Q 1e5 [J]", 1e5*1e3*1e5, 1e13)
p("  in kilotons TNT", 1e13/4.184e12, 2.4)
p("  correct stored U = P Q/omega, 1e5 W, Q=1e5, 5 GHz [J]", 1e5*1e5/(2*pi*5e9), 2.0)
print("== Baker ==")
p("Baker constant 32G/(5c^5) [W s^2/(N^2 m^2)... ]", 32*G/(5*c**5), 1.76e-52)
p("rotating-rod 2G/(45c^5)*4 /(32G/5c^5) i.e. LHS/RHS of laser Eq(1) for binary", (8/45)/(128/5))
p("Physics Procedia 2012: P=1.76e-52(0.5*3.1e13/4e-10)^2 [W]", 1.76e-52*(0.5*3.1e13/4e-10)**2, 2.64e-7)
F0,r,wd=2.0,0.25,2*pi*2.5e9
P1std=16*G/(15*c**5)*r**2*F0**2*wd**2; P1b=1.76e-52*(0.5*2/4e-10)**2
p("single FBAR pair, standard quadrupole (16G/15c^5) r^2 F0^2 w^2", P1std)
p("single FBAR pair, Baker 1.76e-52(2r dF/dt)^2", P1b)
N=1.55e13
p("  N^2 coherent (needs all sources within ~lambda), standard", N**2*P1std)
p("  N incoherent, standard", N*P1std)
amp=lambda F,nu: sqrt(8*pi*G*F/(c**3*(2*pi*nu)**2))
p("Baker A=(8piG F/c^3 w^2)^1/2 coefficient 1.28e-18", sqrt(8*pi*G/c**3)/(2*pi), 1.28e-18)
p("FBAR SPESIF09: F=3.53 W/m^2, lambda=6.1 cm -> A", amp(3.53,c/0.061), 4.9e-28)
p("  same flux, standard one-pol A=(32piG F/c^3w^2)^1/2", 2*amp(3.53,c/0.061))
p("PhysProc12: S=3.55e-9 W/m^2 at 5 GHz -> A", amp(3.55e-9,5e9), 1.52e-32)
p("Laser STAIF06: df=(1+R)(h/lambda)(P dt/(hc/l))/dt [N]", (1.95)*(h_pl/800e-9)*(23e12*33.9e-15/(h_pl*c/800e-9))/33.9e-15, 1.5e5)
m1,m2,a,e=1.4398*Msun,1.3886*Msun,1.95e9,0.6171
fe=(1+73/24*e**2+37/96*e**4)/(1-e**2)**3.5
p("PSR1913+16 Peters-Mathews P_GW [W] (Baker quotes 1.01e25)", 32/5*G**4*m1**2*m2**2*(m1+m2)/(c**5*a**5)*fe)
p("  Baker jerk form 1.76e-52(2 r nu df)^2, r=2.05e9", 1.76e-52*(2*2.05e9*7.2e-5*1.16e33)**2, 1.01e25)
print("== Pais ==")
sg,R,om=5e4,2.0,30000*2*pi/60
p("US10144532 S=sigma^2/eps0 R omega, 5e4 C/m^2, 2 m, 30k RPM", sg**2/eps0*R*om, 1e24)
p("  surface field sigma/eps0 [V/m] (vacuum limit 1e9 per US10135366)", sg/eps0)
smax=eps0*1e9; p("  same formula with sigma at 1e9 V/m [W/m^2]", smax**2/eps0*R*om)
p("US10322827 eq2 sigma=Rv=t_op=fG=1, nu=1e7 [W/m^2]", 1/eps0*(1e7)**3, 1e33)
p("  ... with nu taken as angular 2pi*1e7", 1/eps0*(2*pi*1e7)**3)
p("  surface field for sigma=1 C/m^2 [V/m]", 1/eps0)
p("US10322827 eq4 as written (4piG/c^4)B^2L^2 P, B=1e10 T,L=1m,P=1e33 W", 4*pi*G/c**4*1e20*1*1e33, 1e10)
p("  SI-correct 4piG B^2L^2 P/(mu0 c^4)", 4*pi*G*1e20*1e33/(mu0*c**4))
p("  gravitons/s at JASON's implied 1e-25 J/graviton", 1e10/1e-25, 1e35)
p("  B=1e10 T energy density B^2/2mu0 [J/m^3] (patent: >1e25)", 1e20/(2*mu0))
Es=1.32e18; nuv=(eps0*Es**2*8*pi**2*c**3/hbar)**0.25
p("US10144532 eq4 eps0 E_S^2 = hbar nu^4/(8pi^2c^3) -> nu [Hz]", nuv, 1e22)
p("Schwinger field m_e^2 c^3/(e hbar) [V/m]", (9.1093837e-31)**2*c**3/(1.602176634e-19*hbar))
