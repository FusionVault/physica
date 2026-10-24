# Recompute Coulomb-barrier / Gamow / Bosch-Hale / Lawson numbers (CODATA 2018 constants).
import math
alpha = 7.2973525693e-3
hbarc = 197.3269804          # MeV fm
e2 = alpha*hbarc             # MeV fm  (= 1.439964...)
keV_J = 1.602176634e-16
u = 931494.10242             # keV
me = 510.99895               # keV
# nuclear masses (keV), CODATA 2018 where available
m = {'p':938272.08816,'d':1875612.9424,'t':2808921.1325,'h':2808391.6074,
     'b11': 11.009305166*u - 5*me + 0.6704}   # 11B atomic -> nuclear (approx. e- binding 0.67 keV)
Z = {'p':1,'d':1,'t':1,'h':2,'b11':5}
A = {'p':1,'d':2,'t':3,'h':3,'b11':11}
pairs = [('d','d'),('d','t'),('d','h'),('p','b11'),('p','p')]
def mu(a,b): return m[a]*m[b]/(m[a]+m[b])
print("e2 [MeV fm] =", e2)
print("\n# reduced mass, B_G, E_G")
for a,b in pairs:
    mr = mu(a,b); BG = math.pi*alpha*Z[a]*Z[b]*math.sqrt(2*mr)
    print(f"{a}-{b}: mr c2 = {mr:.1f} keV  B_G = {BG:.4f} keV^1/2  E_G = {BG**2:.1f} keV")
print("\n# Coulomb barrier V = Z1Z2 e2 / R")
for a,b in pairs[:4]:
    for r0 in (1.2,1.4):
        R = r0*(A[a]**(1/3)+A[b]**(1/3))
        print(f"{a}-{b}: r0={r0} R={R:.3f} fm  V_C={Z[a]*Z[b]*e2/R*1000:.0f} keV")
print("\n# classical turning point r_tp = Z1Z2 e2/E (fm), eta, exp(-2 pi eta)")
for a,b in pairs[:4]:
    mr = mu(a,b); BG = math.pi*alpha*Z[a]*Z[b]*math.sqrt(2*mr)
    for E in (1,10,100):
        eta = Z[a]*Z[b]*alpha*math.sqrt(mr/(2*E))
        print(f"{a}-{b} E_cm={E} keV: r_tp={Z[a]*Z[b]*e2/(E/1000):.1f} fm eta={eta:.4f} 2pi*eta={2*math.pi*eta:.4f} (B_G/sqrtE={BG/math.sqrt(E):.4f}) P={math.exp(-2*math.pi*eta):.4e}")
print("\n# Gamow peak E0=(E_G kT^2/4)^(1/3), width (1/e full) D=4 sqrt(E0 kT/3)")
kB = 8.617333262e-8  # keV/K
cases = [('p','p',15.7e6),('p','p',1.55e7),('d','t',10/kB),('d','t',15/kB),('d','d',10/kB),('d','h',10/kB),('p','b11',15/kB),('p','b11',300/kB)]
for a,b,T in cases:
    kT = kB*T; EG = 2*mu(a,b)*(math.pi*alpha*Z[a]*Z[b])**2
    E0 = (EG*kT**2/4)**(1/3); D = 4*math.sqrt(E0*kT/3)
    T9 = T/1e9; muA = A[a]*A[b]/(A[a]+A[b])
    E0i = 0.1220*(Z[a]**2*Z[b]**2*muA*T9**2)**(1/3)*1000; Di = 0.2368*(Z[a]**2*Z[b]**2*muA*T9**5)**(1/6)*1000
    tau = 3*E0/kT
    print(f"{a}-{b} T={T:.4g} K (kT={kT:.4g} keV): E0={E0:.3f} keV  D={D:.3f} keV | Iliadis approx E0={E0i:.3f} D={Di:.3f} | tau=3E0/kT={tau:.2f} exp(-tau)={math.exp(-tau):.3e}")

# ---------------- Bosch-Hale ----------------
BH = {  # BG, mrc2, C1..C7  (Bosch & Hale 1992 Table VII via UWFDM-1268 / MFE notebook)
 'DT':  (34.3827,1124656,[1.17302e-9,1.51361e-2,7.51886e-2,4.60643e-3,1.35000e-2,-1.06750e-4,1.36600e-5]),
 'DDn': (31.3970,937814,[5.43360e-12,5.85778e-3,7.68222e-3,0,-2.96400e-6,0,0]),
 'DDp': (31.3970,937814,[5.65718e-12,3.41267e-3,1.99167e-3,0,1.05060e-5,0,0]),
 'DHe3':(68.7508,1124572,[5.51036e-10,6.41918e-3,-2.02896e-3,-1.91080e-5,1.35776e-4,0,0]),
}
def bh_sv(r,T):
    BG,mrc2,C = BH[r]; C1,C2,C3,C4,C5,C6,C7 = C
    th = T/(1 - T*(C2+T*(C4+T*C6))/(1+T*(C3+T*(C5+T*C7))))
    xi = (BG**2/(4*th))**(1/3)
    return C1*th*math.sqrt(xi/(mrc2*T**3))*math.exp(-3*xi)   # cm^3/s
uw = {  # UWFDM-1268 Table 6 (Bosch-Hale model values, cm3/s)
 1:(6.86e-21,9.93e-23,1.02e-22,3.05e-26),2:(2.98e-19,3.11e-21,3.15e-21,1.40e-23),5:(1.37e-17,9.13e-20,9.02e-20,6.36e-21),
 10:(1.14e-16,6.02e-19,5.78e-19,2.12e-19),20:(4.33e-16,2.60e-18,2.40e-18,3.48e-18),50:(8.65e-16,1.13e-17,9.84e-18,5.55e-17),100:(8.45e-16,2.68e-17,2.24e-17,1.72e-16)}
print("\n# Bosch-Hale <sigma v> (cm3/s) computed vs UWFDM-1268 Table 6")
for T,ref in uw.items():
    vals = [bh_sv(r,T) for r in ('DT','DDn','DDp','DHe3')]
    print(f"T={T:>3} keV: " + "  ".join(f"{v:.4e}({rv:.2e},{(v/rv-1)*100:+.1f}%)" for v,rv in zip(vals,ref)))
print("MFE notebook D(d,n)3He endpoints (m3/s): T=0.1 ->", bh_sv('DDn',0.1)*1e-6, " T=100 ->", bh_sv('DDn',100)*1e-6)
# peak of D-T reactivity
Ts=[1+0.01*i for i in range(0,9900)]
Tp=max(Ts,key=lambda T:bh_sv('DT',T)); print(f"DT reactivity max (within fit to 100 keV) at T={Tp:.1f} keV: {bh_sv('DT',Tp):.4e} cm3/s")

# Bosch-Hale S-factor cross sections (keV, mb)
SF = {
 'DT':  (34.3827,[6.927e4,7.454e8,2.050e6,5.2002e4,0],[6.38e1,-9.95e-1,6.981e-5,1.728e-4]),
 'DDp': (31.3970,[5.5576e4,2.1054e2,-3.2638e-2,1.4987e-6,1.8181e-10],[0,0,0,0]),
 'DDn': (31.3970,[5.3701e4,3.3027e2,-1.2706e-1,2.9327e-5,-2.5151e-9],[0,0,0,0]),
 'DHe3':(68.7508,[5.7501e6,2.5226e3,4.5566e1,0,0],[-3.1995e-3,-8.5530e-6,5.9014e-8,0]),
}
def S(r,E):
    BG,Aa,Bb = SF[r]
    num = Aa[0]+E*(Aa[1]+E*(Aa[2]+E*(Aa[3]+E*Aa[4]))); den = 1+E*(Bb[0]+E*(Bb[1]+E*(Bb[2]+E*Bb[3])))
    return num/den
def sigma(r,E): return S(r,E)/E*math.exp(-SF[r][0]/math.sqrt(E))   # mb
print("\n# Bosch-Hale cross sections")
for r,E in (('DT',10),('DT',64),('DDn',50),('DDp',50),('DHe3',250)):
    print(f"{r} E_cm={E} keV: S={S(r,E):.4e} keV mb  sigma={sigma(r,E):.4e} mb")
Es=[1+0.1*i for i in range(5490)]
Ep=max(Es,key=lambda E:sigma('DT',E)); print(f"DT sigma peak at E_cm={Ep:.1f} keV (E_lab,d={Ep*5/3:.1f}): {sigma('DT',Ep):.1f} mb")
# Maxwellian average of S-factor fit as independent check of reactivity fit
def maxw(r,T,mrc2):
    # <sv> = sqrt(8/(pi mu)) (kT)^-3/2 int sigma E exp(-E/kT) dE ; c in cm/s
    c=2.99792458e10; N=20000; Emax=min(40*T, 549 if r=='DT' else 4900)
    h=Emax/N; s=0
    for i in range(1,N+1):
        E=i*h; w=4 if i%2 else 2
        if i==N: w=1
        s+=w*sigma(r,E)*1e-27*E*math.exp(-E/T)
    s*=h/3
    return math.sqrt(8/(math.pi*mrc2))*c*T**-1.5*s
print("\n# Maxwellian integral of S-fit vs reactivity fit (cm3/s)")
for r,mr in (('DT',1124656),('DDn',937814),('DDp',937814),('DHe3',1124572)):
    for T in (2,10,50):
        print(f"{r} T={T}: integral={maxw(r,T,mr):.4e}  fit={bh_sv(r,T):.4e}")

# ---------------- Lawson (Wurzel & Hsu 2022 eqs 16, 22) ----------------
CB=5.34e-37; eF=17.6e3*keV_J; fc=3.5/17.6
def ntau(T,Q=math.inf):
    q=0 if Q==math.inf else 1/Q
    den=(fc+q)*bh_sv('DT',T)*1e-6*eF/4 - CB*math.sqrt(T)
    return 3*T*keV_J/den if den>0 else math.inf
def nTtau(T,Q=math.inf): return T*ntau(T,Q)
Ts=[1+0.05*i for i in range(1981)]
for Q in (math.inf,1):
    Tm=min(Ts,key=lambda T:nTtau(T,Q)); Tn=min(Ts,key=lambda T:ntau(T,Q))
    print(f"\nQ_fuel={Q}: min nT tau = {nTtau(Tm,Q):.3e} keV s m^-3 at T={Tm:.2f} keV ; min n tau = {ntau(Tn,Q):.3e} s m^-3 at T={Tn:.2f} keV")
# ideal ignition temperature: fc <sv> eF/4 = CB sqrt(T)
lo,hi=2.0,10.0
for _ in range(60):
    mid=(lo+hi)/2
    if fc*bh_sv('DT',mid)*1e-6*eF/4 - CB*math.sqrt(mid) > 0: hi=mid
    else: lo=mid
print(f"ideal ignition T (alpha heating = bremsstrahlung) = {hi:.3f} keV")

# ---------------- screening ----------------
print("\n# electron screening enhancement f=exp(pi eta Ue/E) and exact-shift version")
for a,b,Ue in (('d','d',0.020),('d','d',0.025),('d','d',0.300),('d','h',0.119),('d','h',0.219)):
    mr=mu(a,b)
    for E in (1,2,5,10,20):
        eta=Z[a]*Z[b]*alpha*math.sqrt(mr/(2*E)); f=math.exp(math.pi*eta*Ue/E)
        g = lambda x: math.exp(-2*math.pi*Z[a]*Z[b]*alpha*math.sqrt(mr/(2*x)))
        fx = (E/(E+Ue))*g(E+Ue)/g(E)*(E+Ue)/E  # sigma ~ S/E * P ; shift E->E+Ue in P only
        fx = g(E+Ue)/g(E)
        print(f"{a}-{b} Ue={Ue*1000:.0f} eV E={E} keV: f={f:.4f}  exact P(E+Ue)/P(E)={fx:.4f}")
# Salpeter weak screening for p+p, solar core (assumed rho=150 g/cc, X=0.34, Y=0.64, Zm=0.02)
X,Y,Zm=0.34,0.64,0.02
zeta=math.sqrt(X*1+Y*1+Zm*3.5 + 0.92*(X+Y/2+Zm/2))   # metals: Z^2/A ~ 3.5 (C,N,O avg) (assumption)
for rho,T6 in ((150,15.7),(150,15.5)):
    print(f"Salpeter f0(pp) rho={rho} T6={T6}: zeta={zeta:.3f}  f0={math.exp(0.188*zeta*math.sqrt(rho)*T6**-1.5):.4f}")

# ---------------- muon-catalysed fusion ----------------
lam0=1/2.1969811e-6
print("\n# muCF: X = 1/(ws + lam0/lamc)")
print("lambda0 =",lam0)
print("Kelly et al: ws=0.47%, tc/tau=0.18% ->", 1/(0.0047+0.0018), " lamc =",lam0/0.0018)
for ws in (0.00912,0.0047,0.004,0.003):
    print(f"sticking-only cap 1/ws ws={ws*100:.3f}% -> {1/ws:.0f}")
for Emu in (4.70e3,5e3,6e3):
    print(f"E_mu={Emu/1e3} GeV: fusions for energy breakeven = {Emu/17.6:.0f}; at eta_el=0.4 -> {Emu/17.6/0.4:.0f}; at 0.33 -> {Emu/17.6/0.33:.0f}")
print("150 fusions x 17.6 MeV =",150*17.6/1e3,"GeV ; / m_mu c2 =",150*17.6/105.6583755)

# ---------------- small-device energy balances ----------------
MeV_J=1.602176634e-13
Edd=(3.268+4.033)/2   # mean Q of the two D-D branches, MeV
print("\n# energy balances (computed)")
print("pyro (Naranjo): beam 100 keV x 4 nA =",100e3*4e-9,"W ; 1000 fusions/s ->",1000*Edd*MeV_J,"W")
print("LCF (Steinetz): beam 2.9 MV x 15 mA =",2.9e6*15e-3,"W ; 1.5e3 n/s (x2 branches) ->",2*1.5e3*Edd*MeV_J,"W; ratio",2*1.5e3*Edd*MeV_J/(2.9e6*15e-3))
print("Fusor (Eindhoven): 43 kV x 30 mA =",43e3*30e-3,"W ; 1e6 n/s (x2) ->",2e6*Edd*MeV_J,"W ; Q=",2e6*Edd*MeV_J/(43e3*30e-3))
print("Hirsch ~1e8 DD fusions/s ->",1e8*Edd*MeV_J,"W ; same rate DT ->",1e8*17.59*MeV_J,"W")
print("NIF gains: 3.15/2.05=",3.15/2.05," 3.88/2.05=",3.88/2.05," 5.0/2.05=",5.0/2.05," 8.6/2.08=",8.6/2.08," 3.5/2.065=",3.5/2.065)
# ---------------- Koonin-Nauenberg ----------------
NA=6.02214076e23; yr=3.15576e7
for lab,rate in (('d+d',10**-63.5),('d+d (3e-64)',3e-64),('p+d',1e-55),('d+t',10**-68.9)):
    print(f"{lab}: rate {rate:.2e}/s per molecule; per mole of molecules {rate*NA:.2e}/s ; mean wait for one fusion in a mole {1/(rate*NA)/yr:.2e} yr")

# ---------------- p-11B (Bosch-Hale-form HT fit; Nevins&Swain 2000 vs Tentori&Belloni 2023 Table 2) ----------------
EGb=2*mu('p','b11')*(math.pi*alpha*5)**2; mrb=mu('p','b11')
NS=[4.4467e-14,-5.9357e-2,2.0165e-1,1.0404e-3,2.7621e-3,-9.1653e-6,9.8305e-7]
TB=[9.7827e-14,-5.1610e-2,1.3240e-1,3.8446e-4,1.2499e-3,-5.6715e-6,1.1615e-6]
def pb(C,T,res=False):
    C1,C2,C3,C4,C5,C6,C7=C
    th=T/(1-T*(C2+T*(C4+T*C6))/(1+T*(C3+T*(C5+T*C7)))); xi=(EGb/(4*th))**(1/3)
    v=C1*th*math.sqrt(xi/(mrb*T**3))*math.exp(-3*xi)
    if res: v+=5.41e-21*T**-1.5*math.exp(-148/T)   # N&S resonance term (m3/s), assumption
    return v
print("\n# p-11B reactivity (m3/s); E_G =",EGb," mrc2=",mrb)
for T in (50,70,100,150,200,300,400,500):
    a=pb(NS,T); ar=pb(NS,T,True); b=pb(TB,T)
    print(f"T={T}: NS={a:.4e} NS+res={ar:.4e} TB={b:.4e}  TB/NS-1={(b/a-1)*100:+.1f}%  TB/(NS+res)-1={(b/ar-1)*100:+.1f}%")
Ts=[50+i for i in range(451)]
Tm=max(Ts,key=lambda T:pb(TB,T)); print("TB max in 50-500 keV at",Tm,pb(TB,Tm))
Tm=max(Ts,key=lambda T:pb(NS,T,True)); print("NS(+res) max in 50-500 keV at",Tm,pb(NS,Tm,True))
