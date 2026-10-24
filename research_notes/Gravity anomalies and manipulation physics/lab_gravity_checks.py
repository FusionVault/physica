# Recompute published numbers for lab-gravity notes. Run: python3 -I labgrav.py
import math, cmath
G=6.67430e-11; c=2.99792458e8; hbar=1.054571817e-34; e=1.602176634e-19; eps0=8.8541878128e-12
me=9.1093837e-31; Msun=1.98892e30; amu=1.66053907e-27
def pr(k,v): print(f"{k:58s} {v:.4g}" if isinstance(v,float) else f"{k:58s} {v}")

print("== A. LIGO NCal (Ross et al. 2021) ==")
m=1.0558; M=39.680; rq=0.06033; rh=0.10478; d=1.1802; Phi=math.radians(52.24); z=0.010; L=3995.0
zeta=1+(z/d)**2
def amp_ross2(Phi):  # eq.2 amplitude of cos terms combined
    A=5/(6*zeta)-2/3; B=5/(6*zeta)
    return 9*G*M*m*rq**2/(2*d**4)*zeta**-2.5*abs(A*cmath.exp(-1j*Phi)+B*cmath.exp(-3j*Phi))
def amp_ross3(Phi):
    A=7/(8*zeta)-3/4; B=7/(8*zeta)
    return 15*G*M*m*rh**3/(2*d**5)*zeta**-3.5*abs(A*cmath.exp(-2j*Phi)+B*cmath.exp(-4j*Phi))
F2=amp_ross2(Phi); F3=amp_ross3(Phi)
pr("eq2 point-mass F_x(2f) [pN] (paper 19.04)",F2*1e12)
pr("eq3 point-mass F_x(3f) [pN] (paper 8.97)",F3*1e12)
# brute force point-mass sum, mirror at origin, rotor centre at (d cosPhi, d sinPhi, z)
def Fx(theta, radii_angles):
    s=0.0
    for r,a in radii_angles:
        x=d*math.cos(Phi)+r*math.cos(theta+a); y=d*math.sin(Phi)+r*math.sin(theta+a)
        R=math.sqrt(x*x+y*y+z*z); s+=G*M*m*x/R**3
    return s
N=720
def harmonic(k, ra):
    acc=0j
    for i in range(N):
        th=2*math.pi*i/N; acc+=Fx(th,ra)*cmath.exp(-1j*k*th)
    return 2*abs(acc)/N
quad=[(rq,0),(rq,math.pi)]; hexa=[(rh,0),(rh,2*math.pi/3),(rh,4*math.pi/3)]
pr("brute-force 2f [pN]",harmonic(2,quad)*1e12); pr("brute-force 3f [pN]",harmonic(3,hexa)*1e12)
for f1 in (7.80,8.55,9.56):
    for k,F in ((2,19.16e-12),(3,9.06e-12)):
        fs=k*f1; x=F/(M*(2*math.pi*fs)**2)
        pr(f"  rotor {f1} Hz, line {fs:.2f} Hz: x [m], h=x/L",f"{x:.3e}  h={x/L:.3e}")
pr("G from NCal published",6.56e-11)

print("== B. Virgo NCal (Estevez 2018; Aubin 2024) ==")
def R_sector(b,rho,alpha,rmax,rmin): return G*b*rho*math.sin(alpha)*(rmax**4-rmin**4)/(4*8*math.pi**2)
def P_geom(d,Phi,z):
    Z=1+(z/d)**2
    return math.sqrt((15/Z-6)**2*math.cos(Phi)**2+36*math.sin(Phi)**2)*Z**-2.5/d**4
R=R_sector(0.052,2800,math.radians(45),0.190,0.055); P=P_geom(1.32,0.241,0.643)
pr("Estevez2018 x(f)*f^2 [m Hz^2] (paper 3.3e-14)",R*P)
pr("  sector mass each [kg] (paper ~2.06)",2800*0.052*math.radians(45)/2*(0.190**2-0.055**2))
# O4 Far rotor: Al, two 90deg sectors, D=208mm, b=104.4mm, rmin~0 assumed, d=2.1m, phi=34.5deg, z=0 assumed
for name,rho,dd in (("Far Al",2810,2.1),("Near PVC",1400,1.7)):
    R=R_sector(0.1044,rho,math.pi/2,0.104,0.0); P=P_geom(dd,math.radians(34.5),0.0)
    x=R*P/36**2
    pr(f"  O4 {name}: rotor mass [kg], x(36Hz) [m], h=x/3000",f"{rho*0.1044*math.pi*0.104**2/2:.2f} kg  x={x:.3e}  h={x/3000:.3e}")

print("== C. Westphal 2021 ==")
ms=92.1e-6; mt=90.7e-6
pr("static F at 2.5 mm [N] (paper 9e-14)",G*ms*mt/2.5e-3**2)
pr("a=Gms/d^2 at d0=4.2mm [m/s2]",G*ms/4.2e-3**2)
pr("linear modulation 2Gms dm/d0^3, dm=1.6mm [m/s2] (paper ~3e-10)",2*G*ms*1.6e-3/4.2e-3**3)
pr("human 80 kg at 2.5 m, F on mt [N]",G*80*mt/2.5**2)
pr("Planck mass sqrt(hbar c/G) [ug]",math.sqrt(hbar*c/G)*1e9)

print("== D. Rosi 2014 ==")
keff=2*2*math.pi/780.241e-9; T=0.160
pr("phase per accel k_eff T^2 [rad/(m/s2)]",keff*T**2)
pr("0.547870 rad -> differential accel [m/s2]",0.547870/(keff*T**2))

print("== E. Yukawa / Lee 2020 ==")
rho_d=3.8e3*e/1e-6  # J/m^3
lam_d=(hbar*c/rho_d)**0.25; pr("lambda_d=(hbar c/rho_d)^(1/4) [um] (paper ~85)",lam_d*1e6)
pr("lambda_d with factor 4^(1/4)? (4 hbar c/rho)^(1/4) [um]",(4*hbar*c/rho_d)**0.25*1e6)
pr("boson mass hbar/(lambda c) for 38.6 um [meV] (paper 5.1)",hbar*c/38.6e-6/e*1e3)
def yuk_force_factor(r,alpha,lam): return 1+alpha*(1+r/lam)*math.exp(-r/lam)
pr("force factor r=lam, alpha=1",yuk_force_factor(1,1,1))

print("== F. BMV (Bose et al. 2017) ==")
m=1e-14; dd=450e-6; dx=250e-6; tau=2.5
ph=lambda r: G*m*m*tau/(hbar*r)
p0,pRL,pLR=ph(dd),ph(dd-dx),ph(dd+dx)
pr("phi(d), phi_RL(d-dx), phi_LR(d+dx) [rad]",f"{p0:.3f} {pRL:.3f} {pLR:.3f}")
def witness(dLR,dRL):
    # state: (1/2)[|u>(|u>+e^{i dLR}|d>) + |d>(e^{i dRL}|u>+|d>)], basis uu,ud,du,dd
    psi=[0.5,0.5*cmath.exp(1j*dLR),0.5*cmath.exp(1j*dRL),0.5]
    sx=[[0,1],[1,0]]; sy=[[0,-1j],[1j,0]]; sz=[[1,0],[0,-1]]
    def ev(A,B):
        s=0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        s+=psi[2*i+j].conjugate()*A[i][k]*B[j][l]*psi[2*k+l]
        return s.real
    return abs(ev(sx,sz)-ev(sy,sz))
pr("W with paper phases dLR=-0.2,dRL=0.7 (paper ~1.16)",witness(-0.2,0.7))
pr("W with computed phases",witness(pLR-p0,pRL-p0))
pr("dphi_RL, dphi_LR computed",f"{pRL-p0:.3f} {pLR-p0:.3f}")
pr("mutual accel Gm/(d-dx)^2 [m/s2] (paper ~1e-16)",G*m/(dd-dx)**2)

print("== G. Diosi-Penrose radiation (Donadi 2021) ==")
NA=6.02214076e23; massGe=375*5.323; Na=massGe/72.630*NA; Nz=32
for R0 in (0.54e-10,0.05e-10):
    pref=2*G*e**2*Nz**2*Na/(3*math.pi**1.5*eps0*c**3*R0**3)  # dGamma/domega = pref/omega
    rate=pref*math.log(3800/1000)  # photons/s, 1-3.8 MeV, no efficiency
    pr(f"  R0={R0:.2e} m: photons/day 1-3.8 MeV (eff=1)",rate*86400)
pr("Ge atoms in 375 cm3",Na)
def EDelta_sphere(Mm,Rr,dd):  # Penrose self-energy of difference, uniform sphere, d>=2R
    return G*Mm**2*(6/(5*Rr)-1/dd)
for Mm,Rr in ((1e-14,1e-6),(1e-12,(1e-12/(3*19300)*4/(4*math.pi))**(1/3)),(1e-3,(1e-3/19300*3/(4*math.pi))**(1/3))):
    Ed=EDelta_sphere(Mm,Rr,10*Rr); pr(f"  uniform sphere M={Mm:g} kg R={Rr:.2e} m: hbar/E_Delta [s]",hbar/Ed)

print("== H. Tobar 2024 ==")
Gs=lambda M,vs,h: M*vs**2*h**2/(4*math.pi**3*hbar)
pr("Gamma_stim Al 1800 kg, vs 5.4 km/s, h 5e-22 [Hz] (paper ~1)",Gs(1800,5.4e3,5e-22))
pr("Gamma_spon Nb rho 8570, vs 5e3, R .5, L 1 [Hz] (paper ~1e-33)",8*math.pi*G*8570*(5e3)**4*0.25/(1*c**5))
Mc=1.188*Msun; w=2*math.pi*100
Mopt=24*math.pi**2/5*hbar/((2e-22)**2*(13e3)**2)*(G*Mc/(2*c**3))**(5/3)*w**(8/3)
pr("optimal mass Be, 100 Hz, h0 2e-22, GW170817 [kg] (paper ~15)",Mopt)
pr("graviton energy at 100 Hz [J]",hbar*w)
F=c**3*(2*math.pi*150)**2*(1e-21)**2/(32*math.pi*G)
pr("GW flux h=1e-21 @150 Hz [W/m2]",F); pr("  gravitons per m2 in 0.1 s",F*0.1/(hbar*2*math.pi*150))

print("== I. Overstreet AB ==")
mRb=86.909*amu
for r in (0.07,0.25):
    pr(f"  m*GM/r/hbar, M=1.25 kg, r={r} m [rad/s]",mRb*G*1.25/r/hbar)

print("== J. Frame dragging / Forward ==")
pr("G/c^2 [m/kg]",G/c**2)
rhoW=19300; Rs=0.5; Ms=rhoW*4/3*math.pi*Rs**3; vmax=math.sqrt(1e9/rhoW); om=vmax/Rs
J=0.4*Ms*Rs**2*om
pr("W sphere R=0.5 m mass [kg]",Ms); pr("  omega at v=sqrt(1GPa/rho) [rad/s]",om)
pr("  Omega_LT at pole surface 2GJ/(c^2R^3) [rad/s]",2*G*J/(c**2*Rs**3))
pr("  B_g-type ratio Omega_LT/omega",2*G*J/(c**2*Rs**3)/om)
pr("shell interior 4GM w/(3c^2R), M=1e4,R=1,w=300 [rad/s]",4*G*1e4*300/(3*c**2*1))
gpb=39.2e-3/206264.806/(365.25*86400); pr("GP-B frame drag 39.2 mas/yr [rad/s]",gpb)
pr("induced accel 4Gm a_f/(c^2 r): m=1e4 kg,a_f=1e4,r=1 [m/s2]",4*G*1e4*1e4/(c**2*1))
pr("m/r needed for 9.8 m/s2 with a_f=1e6 [kg/m]",9.8/(4*G/c**2*1e6))

print("== K. Rotating bar GW ==")
Mb=4.9e5; Lb=20; rb=1.0; Om=28
I1=Mb*(Lb**2/12+rb**2/4); I2=Mb*rb**2/2
P=32/5*G/c**5*(I1-I2)**2*Om**6
pr("MTW beam P [W] (sister note 2.3e-29)",P)
pr("  far-zone h at r=lambda=pi c/Om",4*G*(I1-I2)*Om**2/(c**4*(math.pi*c/Om)))
pr("  formal h at 1 m (near zone; not radiative)",4*G*(I1-I2)*Om**2/(c**4*1))
# NCal-size bar: 1 kg, 0.2 m, 100 Hz rotation
I=1*(0.2**2/12); Om=2*math.pi*100
pr("1 kg 0.2 m bar at 100 Hz rotation, P [W]",32/5*G/c**5*I**2*Om**6)

print("== L. Chiao / Kiefer-Weber ==")
pr("Z_g=4piG/c [m2/(kg s)] (paper 2.79e-18)",4*math.pi*G/c)
pr("Z_0 [ohm]",math.sqrt(4e-7*math.pi/eps0))
pr("e^2/(4pi eps0 G me^2)",e**2/(4*math.pi*eps0*G*me**2))
pr("(c/vs)^4 with vs=1e-5 c",1e20)
Mbar=1000; vs=5e3; Lbar=1.5; Q=1e6; w0=2*math.pi*1e3
pr("Weinberg eta=64GMvs^4/(15L^2c^5 Gamma_heat), Gamma=w/Q",64*G*Mbar*vs**4/(15*Lbar**2*c**5*(w0/Q)))

print("== M. Yukawa coupling <-> alpha, BMV witness optimum ==")
mN=1.67262192e-27
g=5e-14; pr("alpha = g^2 hbar c/(4 pi G mN^2) for g_N=5e-14",g*g*hbar*c/(4*math.pi*G*mN**2))
pr("lambda for scalar mass 1 eV: hbar c/(1 eV) [m]",hbar*c/e)
# Optimum over local bases of |<sx sz> - <sy sy>| for the BMV state is 1 + C, C=|sin((dLR+dRL)/2)| (checked numerically in witness_opt section)
for dLR,dRL in ((-0.2,0.7),(pLR-p0,pRL-p0)):
    pr(f"  1+|sin((dLR+dRL)/2)| for ({dLR:.3f},{dRL:.3f})",1+abs(math.sin((dLR+dRL)/2)))
