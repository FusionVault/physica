# Gravity-field mathematics for satellite gravity data: an implementable formula reference

Scope: these are formulas a developer can code directly, in Rust or another language, to turn spherical-harmonic (SH) gravity models (static models such as EGM2008 or EIGEN, and monthly GRACE/GRACE-FO Level-2 models) into gravity-field functionals and hydrology products. Each section has three parts. *Cited Findings* are taken from the source named on each line. *Inferences* are my own derivations or synthesis; where I checked one by hand, it is marked as checked. *Gaps* are things I could not verify.

Notation used throughout, unless a source line says otherwise:
- $r,\varphi,\lambda$: geocentric (spherical) radius, latitude and longitude. $\phi$: geodetic (ellipsoidal) latitude. $\theta = 90^\circ-\varphi$: colatitude.
- $t=\sin\varphi=\cos\theta$, $u=\cos\varphi=\sin\theta$.
- $n$ (or $\ell$, $l$): degree. $m$: order.
- $GM$ and $R$ (or $a_e$): the scaling parameters of the model.
- $\bar C_{nm},\bar S_{nm}$: fully (4π) normalised Stokes coefficients. $\bar P_{nm}$: 4π fully normalised associated Legendre functions (ALFs), without the Condon–Shortley phase.

---

## 1. The geopotential SH expansion, 4π normalisation, and geodetic vs geocentric coordinates

### Takeaway
The conventional model is $V=\frac{GM}{r}\sum_n (a_e/r)^n\sum_m(\bar C_{nm}\cos m\lambda+\bar S_{nm}\sin m\lambda)\bar P_{nm}(\sin\varphi)$. The ALFs are 4π-normalised, $\bar P_{nm}=\sqrt{(2-\delta_{m0})(2n+1)(n-m)!/(n+m)!}\,P_{nm}$, and geodesy omits the Condon–Shortley phase. The expansion must be evaluated in **geocentric spherical** coordinates. Any geodetic input $(\phi,\lambda,h)$ must therefore go through ECEF first, and the coefficients must be rescaled whenever $GM$ or $R$ differ between two models.

### Cited Findings
- **The expansion (IERS 2010, eq. 6.1).**
  $$V(r,\varphi,\lambda)=\frac{GM}{r}\sum_{n=0}^{N}\left(\frac{a_e}{r}\right)^{n}\sum_{m=0}^{n}\big[\bar C_{nm}\cos(m\lambda)+\bar S_{nm}\sin(m\lambda)\big]\bar P_{nm}(\sin\varphi),\qquad \bar S_{n0}=0$$
  Here $\bar C_{nm},\bar S_{nm}$ are the normalised geopotential coefficients and $\bar P_{nm}$ the normalised ALFs — [IERS Conventions 2010, ch. 6](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)
- **Normalisation (IERS 2010, eqs. 6.2a–b, 6.3).** $\bar P_{nm}=N_{nm}P_{nm}$ with $N_{nm}=\sqrt{\frac{(n-m)!\,(2n+1)\,(2-\delta_{0m})}{(n+m)!}}$. Unnormalised coefficients are $C_{nm}=N_{nm}\bar C_{nm}$ and $S_{nm}=N_{nm}\bar S_{nm}$ — [IERS Conventions 2010, ch. 6](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)
- **Unnormalised ALFs without the Condon–Shortley phase.** $P_{lm}(\mu)=(1-\mu^2)^{m/2}\frac{d^m}{d\mu^m}P_l(\mu)$, with $P_l(\mu)=\frac{1}{2^l l!}\frac{d^l}{d\mu^l}(\mu^2-1)^l$. In SHTOOLS the default `csphase=1` excludes the $(-1)^m$ phase, and `csphase=-1` appends it. The phase does not change power spectra — [SHTOOLS real spherical harmonics](https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html)
- **4π orthogonality.** $\int_\Omega Y_{lm}Y_{l'm'}\,d\Omega=4\pi\,\delta_{ll'}\delta_{mm'}$. Coefficients follow from $f_{lm}=\frac{1}{4\pi}\int_\Omega f\,Y_{lm}\,d\Omega$, and the power per degree is $S_{ff}(l)=\sum_m f_{lm}^2$. The same page gives conversions to Schmidt and orthonormalised forms: Schmidt drops the factor $(2l+1)$, and orthonormalised divides the 4π form by $\sqrt{4\pi}$ — [SHTOOLS real spherical harmonics](https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html)
- **Integral of the squared real ALF.** For 4π-normalised real functions, $\int_{-1}^{1}\bar P_{lm}^2\,d\mu = 2(2-\delta_{0m})$ — [SHTOOLS PlmBar doc](https://github.com/SHTOOLS/SHTOOLS/blob/v4.11.7/src/fdoc/plmbar.md)
- **What the expansion represents.** $\bar C_{00}\cdot GM$ is the model's $GM$, so $\bar C_{00}$ is usually set to 1. The degree-1 terms are zero when the coordinate origin is the geocentre. $\bar C_{21},\bar S_{21}$ relate to the mean rotation pole. A model up to $\ell_{max}$ has $(\ell_{max}+1)^2$ coefficients plus $GM$ and $R$. The expansion uses spherical coordinates, so the spherical latitude $\varphi$ must be distinguished from the ellipsoidal latitude $\phi$ — [Barthelmes STR09/02 rev. 2013, §4.1](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Rescaling between $GM$/$R$ conventions (eq. 111).** If the normal-field coefficients $\hat C^U_{\ell m}$ are given with $GM^U, R^U$, then on the model's $GM, R$:
  $$C^U_{\ell m}=\hat C^U_{\ell m}\,\frac{GM^U}{GM}\left(\frac{R^U}{R}\right)^{\ell}$$
  Apply this before differencing with eq. 109, $C^T_{\ell m}=C^W_{\ell m}-C^U_{\ell m}$ — [Barthelmes STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Centrifugal potential (eqs. 123–124).**
  $$\Phi=\tfrac12\omega^2r^2\cos^2\varphi,\quad \Phi_r=\omega^2r\cos^2\varphi,\quad \Phi_\lambda=0,\quad \Phi_\varphi=-\omega^2r^2\cos\varphi\sin\varphi$$
  — [Barthelmes STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Resolution (eqs. 112, 114).** The shortest resolvable half-wavelength is $\psi_{min}\approx\pi R/\ell_{max}$, or with the equal-area argument $\psi_{min}=4\arcsin(1/(\ell_{max}+1))$. Tabulated values for $\ell_{max}$ = 180, 360, 2000 are 111.1, 55.6, 10.0 km (eq. 112) and 140.7, 70.5, 12.7 km (eq. 114) — [Barthelmes STR09/02 rev. 2013, Table 1](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Geodetic to ECEF.**
  $$X=(N_\phi+h)\cos\phi\cos\lambda,\quad Y=(N_\phi+h)\cos\phi\sin\lambda,\quad Z=((1-e^2)N_\phi+h)\sin\phi,\quad N_\phi=\frac{a}{\sqrt{1-e^2\sin^2\phi}}$$
  Geocentric and geodetic longitudes are equal. For ECEF to geodetic, the same source lists Bowring (Newton–Raphson; one iteration is usually enough), Ferrari/Vermeille closed forms, and the Zhu/Heikkinen closed form, with references (Bowring 1976; Vermeille 2002, J. Geod. 76(8); Fukushima 1999, J. Geod. 73(11)) — [Wikipedia: Geographic coordinate conversion](https://en.wikipedia.org/wiki/Geographic_coordinate_conversion)

### Inferences
- **Implementation recipe for a point given as $(\phi,\lambda,h)$ on GRS80/WGS84.**
  1. Compute $X,Y,Z$.
  2. Set $r=\sqrt{X^2+Y^2+Z^2}$, $\varphi=\operatorname{atan2}(Z,\sqrt{X^2+Y^2})$, and keep $\lambda$.
  3. Evaluate the SH series with $t=\sin\varphi$, $u=\cos\varphi$.
  4. To express vectors in the local ellipsoidal (geodetic) frame, rotate the spherical-frame gradient by the angle $\phi-\varphi$ about the east axis.
  
  This is direct geometry built on the cited ECEF formulas.
- **Use $\operatorname{atan2}$, never $\arcsin(Z/r)$, for the geocentric latitude.** It keeps precision near the poles. Store $t$ and $u$ directly and avoid recomputing trigonometric functions.
- **Index layout.** SHTOOLS packs $(l,m)$ at `l*(l+1)/2+m` (0-based, one triangle per coefficient type) — [PlmBar doc](https://github.com/SHTOOLS/SHTOOLS/blob/v4.11.7/src/fdoc/plmbar.md). A Rust implementation can use the same triangular packing with separate $C$ and $S$ vectors, which keeps the inner loops cache-friendly.

### Gaps
- I did not retrieve Heiskanen & Moritz or Hofmann-Wellenhof & Moritz text. The coordinate formulas above come from Wikipedia, not from those primary textbooks.

---

## 2. Fully normalised ALFs: stable high-degree recursions, derivatives, Clenshaw/Horner summation, underflow

### Takeaway
Up to degree about 2700–2800, use the **modified forward column recursion** of Holmes & Featherstone (2002). It recurses on $\bar P_{nm}/u^m$, seeds the sectorals with a scale factor around $10^{-280}$, and applies $u^m$ by Horner's scheme over $m$. That gives values and first derivatives at all latitudes; at the poles, the first derivatives need a separate formula. Above that degree, or at the poles for derivatives, use extended-exponent arithmetic: Fukushima's (2012) X-numbers, or Wittwer et al. (2008).

### Cited Findings
- **Holmes & Featherstone (2002), J. Geod. 76:279–299.** Near the poles, ultra-high-degree ALFs span thousands of orders of magnitude, so standard recursions fail. Their modified recursions give scaled ALFs and first derivatives, combined by Horner's scheme. Partial sums are complete to degree and order 2700 at all latitudes, except at the poles for first derivatives. The method is at least as precise and efficient as Clenshaw's method and simpler to implement — [Holmes & Featherstone 2002 (DOI)](https://www.doi.org/10.1007/S00190-002-0216-2); [Curtin espace record](https://espace.curtin.edu.au/handle/20.500.11937/22940)
- **SHTOOLS.** It uses a standard three-term recursion with Holmes & Featherstone scaling to prevent overflow, and is accurate to about degree 2800 — [SHTOOLS PlmBar doc](https://github.com/SHTOOLS/SHTOOLS/blob/v4.11.7/src/fdoc/plmbar.md)
- **Orekit.** It implements the Holmes–Featherstone *modified forward row* method in its `HolmesFeatherstoneAttractionModel`, a Java reference implementation for accelerations — [Orekit source (via search)](https://www.orekit.org/static/jacoco/org.orekit.forces.gravity/HolmesFeatherstoneAttractionModel.java.html)
- **Recursion as implemented in gravity-toolkit `plm_holmes`/`plm_colombo`.** With $x=t$ and $u=\sqrt{1-x^2}$ (where $u=0$ is replaced by machine epsilon):
  - Seeds: $\bar P_{00}=1$, $\bar P_{10}=\sqrt3\,t$, $\bar P_{11}=\sqrt3\,u$.
  - Sectoral: $\bar P_{mm}=u\sqrt{\frac{2m+1}{2m}}\,\bar P_{m-1,m-1}$.
  - First off-diagonal: $\bar P_{m+1,m}=\sqrt{2m+3}\;t\,\bar P_{mm}$.
  - General: $\bar P_{nm}=a_{nm}\,t\,\bar P_{n-1,m}-b_{nm}\,\bar P_{n-2,m}$, with
  $$a_{nm}=\sqrt{\frac{(2n-1)(2n+1)}{(n-m)(n+m)}},\qquad b_{nm}=\sqrt{\frac{(2n+1)(n+m-1)(n-m-1)}{(n-m)(n+m)(2n-3)}}$$
  - Holmes variant: the scale factor is `scalef = 1.0e-280`. The sectoral seed `pmm = sqrt(2)*scalef` is multiplied by $\sqrt{2m+1}/\sqrt{2m}$ at each order *without* the $u$ factor. Columns are later rescaled by `rescalem = (1/scalef)·u^m`.
  
  — [gravity-toolkit associated_legendre.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/associated_legendre.py)
- **Derivative as implemented in gravity-toolkit (with respect to colatitude $\theta$).** $\frac{d\bar P_{nm}}{d\theta}=\frac{1}{u}\big(n\,t\,\bar P_{nm}-f_{nm}\bar P_{n-1,m}\big)$, with $f_{nm}=\sqrt{\frac{(n^2-m^2)(2n+1)}{2n-1}}$. For sectorals, $\frac{d\bar P_{mm}}{d\theta}=m\frac{t}{u}\bar P_{mm}$. The only guard against the polar singularity is replacing $u=0$ by machine epsilon — [gravity-toolkit associated_legendre.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/associated_legendre.py)
- **Fukushima (2012), J. Geod. 86:271–285, doi:10.1007/s00190-011-0519-2.** "Numerical computation of spherical harmonics of arbitrary degree and order by extending exponent of floating point numbers". A number is stored as a pair: an IEEE-754 double plus a signed integer exponent of a large power-of-2 radix (an "X-number"), building on Smith et al. (1981). A follow-up in J. Geod. 86:1019–1028 (2012) covers first to third derivatives — [science.ecosyste.ms summary](https://science.ecosyste.ms/projects/28565); [DPRI Kyoto page](https://www.dpac.dpri.kyoto-u.ac.jp/en/?p=343)
- **A degree limit reported in another field.** One source (atmospheric modelling, not geodesy) reports that plain double precision with a three-term recurrence can fail to give accurate ALFs above wavenumber about 1700, and that X-numbers fix this — [WCRP Blue Book 2016, Juang](https://wcrp-climate.org/WGNE/BlueBook/2016/individual-articles/06_Juang_Hann_Ming_Henry_Numerics.pdf). This conflicts in magnitude with Holmes & Featherstone's 2700 for the scaled recursion. The difference is probably scaled versus unscaled recursion and which latitudes were tested.
- **Wittwer, Klees, Seitz & Heck (2008).** "Ultra-high degree spherical harmonic analysis and synthesis using extended-range arithmetic", J. Geod. 82(4–5):223–229 — [KIT library record](https://publikationen.bibliothek.kit.edu/1000011597). Rexer & Hirt describe extended-range (X-number) arithmetic as the way to avoid under- and overflow at ultra-high degree — [Rexer & Hirt (TUM mediaTUM)](https://mediatum.ub.tum.de/doc/1369255/document.pdf)
- **Other recursions.** gravity-toolkit also ships a Mohlenkamp/Jacobi-polynomial (Szegő) recursion, `plm_mohlenkamp`, as an alternative method — [gravity-toolkit associated_legendre.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/associated_legendre.py)

### Inferences
- **Clenshaw over degree for fixed order $m$.** This derives from the cited three-term recurrence. To evaluate $S_m(t)=\sum_{n=m}^{N}c_n\bar P_{nm}(t)$, with $c_n=(R/r)^n\bar C_{nm}$ (and likewise for $\bar S$):
  $$y_{N+1}=y_{N+2}=0,\qquad y_n=c_n+a_{n+1,m}\,t\,y_{n+1}-b_{n+2,m}\,y_{n+2}\quad(n=N,\dots,m),\qquad S_m=y_m\,\bar P_{mm}$$
  The result reduces to $y_m\bar P_{mm}$ because $b_{m+1,m}=0$ and $a_{m+1,m}=\sqrt{2m+3}$, which matches the off-diagonal seed. Checked by hand.
- **Horner over order (the Holmes–Featherstone idea).** Write $\bar P_{mm}=u^m\,\sigma_m$ with $\sigma_m=\sqrt3\prod_{k=2}^{m}\sqrt{(2k+1)/(2k)}$ for $m\ge1$, computed in scaled form as $10^{-280}\sigma_m$. Then
  $$V=\frac{GM}{r}\sum_{m=0}^{M}u^m\big[A_m\cos m\lambda+B_m\sin m\lambda\big]$$
  where $A_m,B_m$ are the Clenshaw/forward-column sums divided by $u^m$, in scaled units. Evaluate it as $(\dots((Z_M u+Z_{M-1})u+Z_{M-2})\dots)u+Z_0$, multiplying by $10^{280}$ as late as possible. The tiny $u^m$ is never formed explicitly, which is the source of the stability.
- **Pole-safe colatitude derivative.** This is a standard identity, not from a fetched source; I checked it by hand for $n\le2$.
  $$\frac{d\bar P_{nm}}{d\theta}=\tfrac12\Big[k_m\sqrt{(n+m)(n-m+1)}\,\bar P_{n,m-1}-\sqrt{(n-m)(n+m+1)}\,\bar P_{n,m+1}\Big]$$
  Here $k_1=\sqrt2$ and $k_m=1$ for $m\ge2$. For $m=0$ use $\frac{d\bar P_{n0}}{d\theta}=-\sqrt{n(n+1)/2}\,\bar P_{n1}$. The latitude derivative is $d/d\varphi=-d/d\theta$. This form avoids $1/u$ and is preferable for points near the poles, such as polar-orbiting satellites.
- **Self-checks to unit-test.**
  - The addition-theorem identity $\sum_{m=0}^{n}\bar P_{nm}(t)^2=2n+1$ holds for every $t$ in 4π normalisation.
  - The quadrature check $\int_{-1}^1\bar P_{nm}^2=2(2-\delta_{m0})$, which is cited.
  
  The first identity is very sensitive to underflow failures near the poles.
- **Degree thresholds for Rust.**
  - Up to about 2190 (EGM2008): plain `f64` with Holmes–Featherstone scaling is enough.
  - Above about 2700, or for derivatives exactly at the poles: implement an X-number type `(f64, i32)` with renormalisation after each multiply/add, used only in the sectoral seeds and the first few recursion steps of each column.
  - The design keeps the hot loop in plain `f64` once values are back in range. This follows the published approach at the concept level.

### Gaps
- I could not read the Holmes & Featherstone (2002) full text. Their exact equation numbers, and their recommended treatment of the first derivative *at* the pole, are unverified.
- Fukushima's X-number radix (commonly quoted as $B=2^{960}$) and its normalisation thresholds are from memory and **unverified**.
- I found no source text for Wittwer et al. (2008) beyond the citation.

---

## 3. Gravity-field functionals per Barthelmes (ICGEM STR09/02), normal field, gradients, deflections, tide systems

### Takeaway
ICGEM computes every functional from $T=W-U$. The coefficients are $C^T=C^W-C^U$, after rescaling $C^U$ to the model's $GM$ and $R$, and the normal field has only even zonal terms. Main results:
- Height anomaly on the ellipsoid: Bruns' formula, $\zeta=T/\gamma$.
- Geoid: the same, minus a topographic term $2\pi G\rho H^2/\gamma$.
- Gravity disturbance and anomaly: either exact, from $|\nabla W|-|\nabla U|$, or with the spherical-approximation factors $(\ell+1)$ and $(\ell-1)$.

The tide-system conversion of $\bar C_{20}$ is a constant offset, $-4.1736\times10^{-9}$ (zero-tide minus tide-free) for EGM2008.

### Cited Findings
**Definitions** — all from [Barthelmes STR09/02 (2009)](https://icgem.gfz.de/docs/str-0902.pdf) and the [rev. 2013 edition](https://icgem.gfz.de/docs/str-0902-revised.pdf):
- **Potentials.** $W=W_a+\Phi$; $U=U_a+\Phi$; $T=W_a-U_a$ (harmonic outside the masses). $\vec g=\nabla W$, $g=|\nabla W|$; $\vec\gamma=\nabla U$, $\gamma=|\nabla U|$.
- **Geoid.** $W(h=N,\lambda,\phi)=U_0$ (eq. 11). Topographic height: $h_t=N+H$ (eq. 12).
- **Height anomaly.** $W(h_t,\lambda,\phi)=U(h_t-\zeta,\phi)$ (eq. 20). Generalised form: $W(h,\lambda,\phi)=U(h-\zeta_g,\phi)$ (eq. 25).
- **Gravity disturbance.** $\delta g=|\nabla W(h,\lambda,\phi)|-|\nabla U(h,\phi)|$ (eq. 27).
- **Classical gravity anomaly.** $\Delta g_{cl}=|\nabla W^c(N,\lambda,\phi)|-|\nabla U(0,\phi)|$ (eq. 30).
- **Modern (Molodensky) anomaly.** $\Delta g(h,\lambda,\phi)=|\nabla W(h,\lambda,\phi)|-|\nabla U(h-\zeta_g,\phi)|$ (eq. 31).

**Iterative geoid and height anomaly:**
- Geoid: $N_{i+1}=N_i+\frac{1}{\gamma(0,\phi)}[W(N_i,\lambda,\phi)-U_0]$ (eq. 41). With $N_0=0$ this gives Bruns, $N_1=T(0)/\gamma(0)$ (eq. 43). The second step is $N_2\approx T(N_1)/\gamma(0)$ (eq. 50).
- Height anomaly: $\zeta_{i+1}=\zeta_i+\frac{1}{\gamma(0)}[W(h_t)-U(h_t-\zeta_i)]$ (eq. 82).
- In 64-bit arithmetic both converge in 2–3 steps (geoid) or 3–4 steps (height anomaly).

**Topography and the geoid (spherical-shell approximation, eqs. 68–70, 84):**
- $V_t(N)-V_t^c(N)\approx-2\pi G\rho H^2$.
- $N\approx N_1^c-\frac{2\pi G\rho H^2}{\gamma(0,\phi)}=\frac{T^c(0,\lambda,\phi)-2\pi G\rho H^2}{\gamma(0,\phi)}$.
- $N-\zeta_e\approx\frac{2\pi G\rho H^2}{\gamma(0,\phi)}$.

— [STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf)

**SH formulas as used by ICGEM** — [STR09/02 rev. 2013 §4](https://icgem.gfz.de/docs/str-0902-revised.pdf). In all of these, $r_e=r_e(\varphi)$ is the geocentric radius of the ellipsoid at the computation point.
- **Height anomaly on the ellipsoid, i.e. the pseudo-height anomaly (eq. 118):**
  $$\zeta_{e1}(\lambda,\varphi)=\frac{GM}{r_e\,\gamma(r_e,\varphi)}\sum_{\ell=0}^{\ell_{max}}\left(\frac{R}{r_e}\right)^{\ell}\sum_{m=0}^{\ell}\bar P_{\ell m}(\sin\varphi)\big(C^T_{\ell m}\cos m\lambda+S^T_{\ell m}\sin m\lambda\big)$$
- **Geoid with the topographic correction (eq. 116):**
  $$N^s_1=\zeta_{e1}-\frac{2\pi G\rho}{\gamma(r_e,\varphi)}\Big[R\sum_{\ell,m}\bar P_{\ell m}(\sin\varphi)\big(C^{topo}_{\ell m}\cos m\lambda+S^{topo}_{\ell m}\sin m\lambda\big)\Big]^2$$
  The bracket is $H$ from a topography SH model (eq. 115). Eq. 117 is the same with $r=r_\zeta=r_e+\zeta_e$ in the $T$ term.
- **Height anomaly at the surface (eq. 119):**
  $$\tilde\zeta_1=\zeta_{e1}+\frac{H+N}{\gamma(0,\varphi)}\left.\frac{\partial T^c}{\partial r}\right|_{r=r_e}$$
- **Exact gravity (eqs. 120–122):**
  $$|\nabla W|=\sqrt{[W_{a,r}+\Phi_r]^2+\Big[\tfrac{1}{r\cos\varphi}(W_{a,\lambda}+\Phi_\lambda)\Big]^2+\Big[\tfrac1r(W_{a,\varphi}+\Phi_\varphi)\Big]^2}$$
  with
  $$W_{a,r}=-\frac{GM}{r^2}\sum_\ell\left(\frac Rr\right)^\ell(\ell+1)\sum_m\bar P_{\ell m}\big(C^W_{\ell m}\cos m\lambda+S^W_{\ell m}\sin m\lambda\big)$$
  $$W_{a,\lambda}=\frac{GM}{r}\sum_\ell\left(\frac Rr\right)^\ell\sum_m m\,\bar P_{\ell m}\big(S^W_{\ell m}\cos m\lambda-C^W_{\ell m}\sin m\lambda\big)$$
  $$W_{a,\varphi}=\frac{GM}{r}\sum_\ell\left(\frac Rr\right)^\ell\sum_m\frac{\partial\bar P_{\ell m}(\sin\varphi)}{\partial\varphi}\big(C^W_{\ell m}\cos m\lambda+S^W_{\ell m}\sin m\lambda\big)$$
  $|\nabla U|$ uses the same formulas with $C^U,S^U$. So $\delta g$ and the anomalies are computed **exactly** at arbitrary points.
- **Spherical approximations (rev. 2013 eqs. 125–126; note the signs):**
  $$\delta g_{sa}(r,\lambda,\varphi)=\frac{GM}{r^2}\sum_{\ell}\left(\frac Rr\right)^{\ell}(\ell+1)\sum_m\bar P_{\ell m}\big(C^T_{\ell m}\cos m\lambda+S^T_{\ell m}\sin m\lambda\big)$$
  $$\Delta g_{sa}(r,\lambda,\varphi)=\frac{GM}{r^2}\sum_{\ell}\left(\frac Rr\right)^{\ell}(\ell-1)\sum_m\bar P_{\ell m}\big(C^T_{\ell m}\cos m\lambda+S^T_{\ell m}\sin m\lambda\big)$$
  The underlying definitions are $\delta g_{sa}=-\partial T^c/\partial r$ (eq. 92) and $\Delta g_{sa}=-\partial T^c/\partial r-\frac{2}{r}T^c$ (eqs. 100, 104). The **2009 edition prints these with a leading minus sign**, contradicting its own eq. 92 definition. The rev. 2013 edition prints them positive, consistent with $\delta g=-\partial T/\partial r$ — compare the [2009 edition](https://icgem.gfz.de/docs/str-0902.pdf) with [rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Accuracy of the approximations.** Barthelmes states that "Nowadays this formula ... is not accurate enough for most practical applications". On EIGEN-GL04C, the wrms/min/max differences are:

  | Difference | wrms | min | max |
  |---|---|---|---|
  | $\Delta g_{cl}-\Delta g$ | 0.73 mGal | −12 mGal | 26 mGal |
  | $\Delta g_{sa}-\Delta g$ | 0.69 mGal | −13 mGal | 23 mGal |
  | $\Delta g_{sa}-\Delta g_{cl}$ | 0.11 mGal | −1.48 mGal | 0.83 mGal |
  | $N_1^c-N^c$ | 1.3 mm | −15 mm | 22 mm |
  | $N_2^c-N^c$ | 1.4×10⁻⁷ m | | |
  | $\zeta_{e1}-\zeta$ | 3.8 cm | | max 1.1 m |
  | $\zeta_{e1}-\tilde N_2^s$ | | | 3.67 m (in continental areas) |

  — [STR09/02 (2009), Tables 2–3](https://icgem.gfz.de/docs/str-0902.pdf)
- **Normal-field coefficients.** Only even zonal terms are present. $C^U_{00}, C^U_{20}, C^U_{40}, C^U_{60}$ (sometimes $C^U_{80}$) suffice in practice, and $S^T_{\ell m}=S^W_{\ell m}$ — [STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf). In fully normalised form, $\bar C_{n0}=-J_n/\sqrt{2n+1}$, with $J_n=0$ for odd $n$ — [GeographicLib NormalGravity](https://geographiclib.sourceforge.io/C++/doc/classGeographicLib_1_1NormalGravity.html)

**GRS80 closed formulas** — [Moritz, GRS80 (corrected Geodesist's Handbook text)](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf). This version corrects a $\theta$/$\phi$ interchange that appeared in the 1992 printing.
- **Second-eccentricity function and $e^2$ (solved iteratively):**
  $$2q_0=\Big(1+\frac{3}{e'^2}\Big)\arctan e'-\frac{3}{e'},\qquad e^2=3J_2+\frac{4}{15}\frac{\omega^2a^3}{GM}\frac{e^3}{2q_0}$$
  Also $e'=e/\sqrt{1-e^2}$ and $m=\frac{\omega^2a^2b}{GM}$ — [Moritz](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf); [Wikipedia GRS80](https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980)
- **Normal potential on the ellipsoid:**
  $$U_0=\frac{GM}{E}\arctan e'+\frac13\omega^2a^2,\qquad E=\sqrt{a^2-b^2}$$
- **Zonal series and coefficients:**
  $$V=\frac{GM}{r}\Big[1-\sum_{n\ge1}\Big(\frac ar\Big)^{2n}J_{2n}P_{2n}(\cos\theta)\Big],\qquad J_{2n}=(-1)^{n+1}\frac{3e^{2n}}{(2n+1)(2n+3)}\Big(1-n+5n\frac{J_2}{e^2}\Big)$$
- **Equatorial and polar normal gravity:**
  $$\gamma_e=\frac{GM}{ab}\Big(1-m-\frac m6\frac{e'q_0'}{q_0}\Big),\qquad \gamma_p=\frac{GM}{a^2}\Big(1+\frac m3\frac{e'q_0'}{q_0}\Big),\qquad q_0'=3\Big(1+\frac{1}{e'^2}\Big)\Big(1-\frac{1}{e'}\arctan e'\Big)-1$$
- **Somigliana's closed form:**
  $$\gamma=\frac{a\gamma_e\cos^2\phi+b\gamma_p\sin^2\phi}{\sqrt{a^2\cos^2\phi+b^2\sin^2\phi}}=\gamma_e\frac{1+k\sin^2\phi}{\sqrt{1-e^2\sin^2\phi}},\qquad k=\frac{b\gamma_p}{a\gamma_e}-1$$
- **Series form.** $\gamma=\gamma_e(1+0.0052790414\sin^2\phi+0.0000232718\sin^4\phi+0.0000001262\sin^6\phi+0.0000000007\sin^8\phi)$, with relative error about $10^{-10}$. The older 1930-style series is accurate only to 0.1 mGal.
- **Clairaut check.** The gravity flattening is $f^*=(\gamma_p-\gamma_e)/\gamma_e$, and $f+f^*=\frac{\omega^2b}{\gamma_e}\big(1+\frac{e'q_0'}{2q_0}\big)$.

**Reference constants:**

| Constant | GRS80 (Moritz) | WGS84 (NGA) | EGM2008 (IERS 2010) |
|---|---|---|---|
| $a$ (m) | 6378137 | 6378137.0 | $a_e$ = 6378136.3 |
| $GM$ (m³ s⁻²) | 3986005×10⁸ | 3.986004418×10¹⁴ | 3986004.415×10⁸ (TT-compatible) |
| $J_2$ | 108263×10⁻⁸ | — | — |
| $1/f$ | 298.257222101 | 298.257223563 | — |
| $\omega$ (rad s⁻¹) | 7292115×10⁻¹¹ | 7292115×10⁻¹¹ | — |

- **GRS80 defining constants** are $a$, $GM$, $J_2$ and $\omega$ (row values above).
- **GRS80 derived geometric constants:**
  - $b$ = 6356752.3141 m
  - $E$ = 521854.0097 m
  - $c=a^2/b$ = 6399593.6259 m
  - $e^2$ = 0.00669438002290
  - $e'^2$ = 0.00673949677548
  - $f$ = 0.00335281068118
  - $1/f$ = 298.257222101
  - $R_1$ = 6371008.7714 m, $R_2$ = 6371007.1810 m, $R_3$ = 6371000.7900 m
- **GRS80 derived physical constants:**
  - $U_0$ = 62636860.850 m² s⁻²
  - $J_4$ = −0.00000237091222, $J_6$ = 0.00000000608347, $J_8$ = −0.00000000001427
  - $m$ = 0.00344978600308
  - $\gamma_e$ = 9.7803267715 m s⁻², $\gamma_p$ = 9.8321863685 m s⁻²
  - $f^*$ = 0.005302440112
  - $k$ = 0.001931851353
  
  Moritz advises computing these from the closed formulas when higher accuracy is needed — [Moritz GRS80](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf)
- **WGS84.** The defining parameters are $a$, $1/f$, $\omega$ and $GM$ as in the table. The GPS-user $GM$ is 3.9860050×10¹⁴ m³ s⁻². EGM2008 is approved for DoD use under NGA.STND.0036_1.0 (2014-07-08) — [NGA WGS84 page](https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84). Secondary sources give WGS84 $\gamma_e$ = 9.7803253359 m s⁻² — [hpmuseum forum](https://hpmuseum.org/forum/post-152584.html); [ahrs IGF docs](https://ahrs.readthedocs.io/en/latest/geodesy/igf.html)
- **EGM2008 scaling.** The parameters to use with EGM2008 coefficients are $GM$ = 398600.4415 km³ s⁻² and $a_e$ = 6378136.3 m. The TCG-compatible value 398600.4418 is used with the two-body term — [IERS 2010 ch. 6](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf). GeographicLib's EGM2008 file uses ModelRadius 6378136.3, ModelMass 3986004.415e8, ReferenceMass 3986004.418e8, and **HeightOffset −0.41 m** — [GeographicLib gravity docs](https://geographiclib.sourceforge.io/C++/doc/gravity.html)

**Zero-degree term:**
- EGM2008 geoid heights relative to WGS84 include a zero-degree term of about **−41 cm** — [Argentina Newton report](https://www.isgeoid.polimi.it/Newton/Newton_4/Report_A4_Argentina.pdf); [ACRS 2012 paper using hsynth_WGS84](https://acrs-aars.org/proceeding/ACRS2012/Proceeding%20ACRS%202012/Technical%20Sessions/G1%20Photogrammetry%20(4)/G1-2.pdf); [GeographicLib](https://geographiclib.sourceforge.io/C++/doc/gravity.html)
- GeographicLib notes that the different model and reference masses leave a $1/r$ term in $T$ at large distances — [GeographicLib gravity docs](https://geographiclib.sourceforge.io/C++/doc/gravity.html)
- ICGEM's calculation service offers a "Zero Degree Term" option and a choice of tide system: model's own, tide-free, zero-tide or mean-tide. The form page does not document the formulas — [ICGEM calcgrid](https://icgem.gfz.de/calcgrid)
- **$W_0$ values.** $W_0=(62636855.611\pm0.5)$ m² s⁻² (Bursa et al. 1998), and a local Fennoscandian value of $62636855.8\pm0.36$ m² s⁻². $W_0$ and $R_0=GM/W_0$ are independent of the tide system — [IAG SC3 report, in the same Geodesist's Handbook PDF](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf)

**Deflections of the vertical and anomaly in GeographicLib's spherical approximation.** Here $\vec\delta=\nabla T$ with $x$ east, $y$ north and $z$ up:
- $\Delta g=\delta_z-2T/R$
- $\xi=-\delta_y/\gamma$ (north–south component)
- $\eta=-\delta_x/\gamma$ (east–west component)

The Bruns-based geoid has RMS error about 1.5 mm (max about 10 mm). The zeta-to-N correction lies in [−5.05, 0.05] m — [GeographicLib gravity docs](https://geographiclib.sourceforge.io/C++/doc/gravity.html). The sign of $\delta_z$ in the anomaly expression is as the fetched summary gave it. Check it against GeographicLib's convention, $\delta=\nabla T$ with $z$ up, so $\delta_z=\partial T/\partial r$. Barthelmes' definition is $\Delta g_{sa}=-\partial T/\partial r-2T/r$.

**Gravity gradients.** Conventional local-north-oriented-frame expressions have factors at the ALFs and their derivatives that "tend to infinity when approaching the poles". Petrovskaya & Vershkov built alternative non-singular expressions for the local north-oriented frame and the orbital frame, valid for any orbit, including polar ones (J. Geod. 81:171, 2007, online 2006, doi:10.1007/s00190-006-0087-z) — [Springer record](https://link.springer.com/article/10.1007/s00190-006-0087-z); [Artificial Satellites summary](https://reference-global.com/article/10.2478/v10018-007-0015-3)

**Tide systems:**
- Zero-tide includes the permanent indirect (deformation) contribution in $\bar C_{20}$; "conventional tide-free" excludes it. The conversion is
  $$\Delta\bar C^{perm}_{20}=A_0H_0k_{20}=(4.4228\times10^{-8})(-0.31460)\,k_{20},\qquad A_0=\frac{1}{R_e\sqrt{4\pi}}$$
  For EGM2008, $\bar C_{20}^{zt}-\bar C_{20}^{tf}=-4.1736\times10^{-9}$. The conventional zero-tide value $\bar C_{20}=-0.48416948\times10^{-3}$ (rate $11.6\times10^{-12}$/yr, epoch J2000) becomes $-0.48416531\times10^{-3}$ tide-free. With a zero-tide model, only $\Delta\bar C_{20}-\Delta\bar C_{20}^{perm}$ must be added for solid tides (eq. 6.13) — [IERS 2010 §6.2.2, Table 6.2](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)
- **Definitions.** Tide-free removes all tidal effects. Zero-frequency (zero-tide) includes the indirect distortion but not the direct one. Mean-tide includes both the direct and the indirect permanent distortion. For tide-free conversions, $k_2=0.3$ is adopted. $J_2$ (tide-free, unnormalised) is $1082626.7\times10^{-9}$ versus $1082635.9\times10^{-9}$ zero-tide. Also $dJ_2=-d\bar C_{20}=1.37532\times10^{-8}\,k_s$ (fully normalised), with $k_s=0.9383$ the secular fluid Love number — [IAG SC3 report in Geodesist's Handbook PDF](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf)
- **Solid-tide coefficient changes (IERS eq. 6.6), for completeness:**
  $$\Delta\bar C_{nm}-i\Delta\bar S_{nm}=\frac{k_{nm}}{2n+1}\sum_{j=2}^{3}\frac{GM_j}{GM_\oplus}\Big(\frac{R_e}{r_j}\Big)^{n+1}\bar P_{nm}(\sin\Phi_j)e^{-im\lambda_j}$$
  — [IERS 2010 ch. 6](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)

### Inferences
- **Zero-degree term formula** (standard in Heiskanen & Moritz / Hofmann-Wellenhof & Moritz; not verified in fetched text):
  $$N_0=\frac{GM-GM_{ref}}{r\,\gamma}-\frac{W_0-U_0}{\gamma}$$
  Add it to $\zeta$ or $N$ when the model's $GM$ (with $\bar C_{00}$) differs from the reference ellipsoid's $GM$, and/or when $W_0\ne U_0$. The −0.41 m used for EGM2008/WGS84 is this term evaluated once as a constant.
- **Mean-tide $\bar C_{20}$.** It follows from the SC3 definitions and the IERS formula structure by adding the *direct* permanent term ($k\to1$): $\bar C_{20}^{mt}=\bar C_{20}^{zt}+A_0H_0\approx\bar C_{20}^{zt}-1.3914\times10^{-8}$. The SC3 fully normalised factor ($1.37532\times10^{-8}$) differs from $A_0H_0=1.3914\times10^{-8}$ by about 1%. That reflects different $H_0$ conventions, so pick one and document it. Strictly, the direct permanent tide potential is not harmonic in the exterior like a $\bar C_{20}$ term ($r^2$ versus $r^{-3}$ growth), so "mean-tide $\bar C_{20}$" is a surface convention only.
- **Conventional gravity-gradient tensor in the local north-east-up frame.** I derived and checked this against the Laplacian; it is not from a fetched source. With $V_r,V_\varphi,V_\lambda,V_{rr},\dots$ the spherical partials:
  $$V_{zz}=V_{rr},\qquad V_{xx}=\frac{V_r}{r}+\frac{V_{\varphi\varphi}}{r^2},\qquad V_{yy}=\frac{V_r}{r}-\frac{\tan\varphi}{r^2}V_\varphi+\frac{V_{\lambda\lambda}}{r^2\cos^2\varphi}$$
  $$V_{xy}=\frac{V_{\varphi\lambda}}{r^2\cos\varphi}+\frac{\sin\varphi}{r^2\cos^2\varphi}V_\lambda,\qquad V_{xz}=\frac{V_{r\varphi}}{r}-\frac{V_\varphi}{r^2},\qquad V_{yz}=\frac{V_{r\lambda}}{r\cos\varphi}-\frac{V_\lambda}{r^2\cos\varphi}$$
  Here $x$ is north, $y$ east, $z$ radial up. If $y$ is west (the GOCE LNOF convention), flip the signs of $V_{xy}$ and $V_{yz}$. The trace equals the spherical Laplacian, so $V_{xx}+V_{yy}+V_{zz}=0$ outside the masses, which is a free unit test. The radial terms in SH form are $V_{rr}=\frac{GM}{r^3}\sum(n+1)(n+2)(R/r)^n\sum_m\bar P_{nm}(\cdot)$ and $V_{r}=-\frac{GM}{r^2}\sum(n+1)(R/r)^n\sum_m\bar P_{nm}(\cdot)$, from differentiating eq. 108. The $1/\cos\varphi$ factors are what Petrovskaya & Vershkov remove. For a pole-robust alternative, evaluate the Cartesian Hessian with a Cunningham/Pines-type formulation and rotate.
- **Deflections in SH form** (from the GeographicLib definitions and the $W_{a,\varphi},W_{a,\lambda}$ series above):
  $$\xi=-\frac{1}{\gamma r}\frac{\partial T}{\partial\varphi},\qquad \eta=-\frac{1}{\gamma r\cos\varphi}\frac{\partial T}{\partial\lambda}$$
  These are spherical approximations. Exact (Helmert) deflections need $\nabla W$ and $\nabla U$ expressed in the ellipsoidal-normal frame.
- **Recommended implementation order:**
  1. $T$ coefficients (rescale with eq. 111, subtract the even zonals).
  2. Geocentric coordinates.
  3. One pass of a Holmes–Featherstone summation returning $T, T_r, T_\varphi, T_\lambda$ (and second derivatives if needed).
  4. Exact $\gamma$ from the closed Somigliana formula at $h=0$, or from the normal coefficients via eq. 121 at $h\ne0$.
  5. Functionals.
  
  This reuses one summation for every quantity.

### Gaps
- I did not obtain NGA.STND.0036 or TR8350.2 text. WGS84 derived constants ($\gamma_e$, $\gamma_p$, $U_0$, $\bar C_{20}$ of the normal field, $k$, $e^2$) are therefore not primary-sourced here. Compute them from the four defining constants with the Moritz closed formulas, and cross-check against GRS80's published values.
- The normal-gravity height series $\gamma_h=\gamma[1-\frac2a(1+f+m-2f\sin^2\phi)h+\frac{3}{a^2}h^2]$ is from memory (Heiskanen & Moritz 2-215 / TR8350.2) and not verified. Prefer the exact eq. 121 using the normal coefficients.
- The crustal density ICGEM uses in $2\pi G\rho H^2$ (commonly 2670 kg/m³) was not confirmed. Neither was the exact ICGEM zero-degree formula.
- I did not fetch the full Petrovskaya–Vershkov non-singular formulas.

---

## 4. Time-variable gravity for hydrology: ΔC/ΔS to surface mass and EWH, Love numbers, filtering, C20/C30, degree-1, GIA, leakage, mascons, basin statistics, uncertainty

### Takeaway
The core is Wahr et al. (1998), eq. 14: $\Delta\sigma=\frac{a\rho_{ave}}{3}\sum\frac{2l+1}{1+k_l}\bar P_{lm}(\Delta C\cos+\Delta S\sin)$, with $\rho_{ave}$ = 5517 kg/m³, $\rho_w$ = 1000 kg/m³, and PREM load Love numbers ($k_2$ = −0.303). Then:
- Smoothing: Jekeli Gaussian weights $W_l$, by recursion, inserted as $2\pi W_l$.
- Destriping or decorrelation: the Swenson–Wahr moving quadratic fit per order and parity, or the DDK $(E^{-1}+aS^{-1})^{-1}E^{-1}$ filters.
- Replacements and corrections: $C_{20}$ always, and $C_{30}$ after Aug 2016, from TN-14; degree-1 from TN-13; subtract a GIA rate model.
- Leakage: restore signal with least-squares gain factors, and propagate uncertainty with the covariance-weighted sums below.

### Cited Findings
**Wahr, Molenaar & Bryan (1998), JGR 103(B12):30205–30230** — [WHOI PDF](https://geodynamicsprogram.whoi.edu/wp-content/uploads/sites/14/2018/10/wahr_et_al.1998_117046.pdf):
- **Geoid change (eq. 3).** $\Delta N(\theta,\phi)=a\sum_{l}\sum_m\tilde P_{lm}(\cos\theta)(\Delta C_{lm}\cos m\phi+\Delta S_{lm}\sin m\phi)$. Here $a$ is the Earth radius, $\theta$ colatitude and $\phi$ east longitude. Normalisation (eq. 2) is $\tilde P_{lm}=\sqrt{(2-\delta_{m0})(2l+1)\frac{(l-m)!}{(l+m)!}}P_{lm}$, i.e. the 4π normalisation.
- **Volumetric density to coefficients (eq. 4):**
  $$\begin{Bmatrix}\Delta C_{lm}\\\Delta S_{lm}\end{Bmatrix}=\frac{3}{4\pi a\rho_{ave}(2l+1)}\int\Delta\rho(r,\theta,\phi)\tilde P_{lm}(\cos\theta)\Big(\frac ra\Big)^{l+2}\begin{Bmatrix}\cos m\phi\\\sin m\phi\end{Bmatrix}\sin\theta\,d\theta\,d\phi\,dr$$
  with $\rho_{ave}$ = **5517 kg/m³**.
- **Thin layer and load deformation.** Surface density is $\Delta\sigma=\int\Delta\rho\,dr$ over a thin layer (eq. 5), valid when $(l_{max}+2)H/a\ll1$. Load deformation adds a factor $k_l$ (eq. 7), so the total response scales by $(1+k_l)$ (eq. 8).
- **Surface mass to coefficients (eq. 12):**
  $$\begin{Bmatrix}\Delta C_{lm}\\\Delta S_{lm}\end{Bmatrix}=\frac{3\rho_w}{\rho_{ave}}\frac{1+k_l}{2l+1}\begin{Bmatrix}\Delta\hat C_{lm}\\\Delta\hat S_{lm}\end{Bmatrix}$$
  with
  $$\begin{Bmatrix}\Delta\hat C_{lm}\\\Delta\hat S_{lm}\end{Bmatrix}=\frac{1}{4\pi a\rho_w}\int_0^{2\pi}\!\!d\phi\int_0^\pi\!\!\sin\theta\,d\theta\,\Delta\sigma(\theta,\phi)\tilde P_{lm}(\cos\theta)\begin{Bmatrix}\cos m\phi\\\sin m\phi\end{Bmatrix}\quad\text{(eq. 11)}$$
- **Coefficients to surface mass (eqs. 13–14):**
  $$\Delta\sigma(\theta,\phi)=\frac{a\rho_{ave}}{3}\sum_{l=0}^{\infty}\sum_{m=0}^{l}\tilde P_{lm}(\cos\theta)\frac{2l+1}{1+k_l}\big(\Delta C_{lm}\cos m\phi+\Delta S_{lm}\sin m\phi\big)$$
  $\Delta\sigma/\rho_w$ is the equivalent water thickness, with $\rho_w$ = **1000 kg/m³**. The inverse (eq. 15) is $\Delta N=\frac{3a\rho_w}{\rho_{ave}}\sum\tilde P_{lm}\frac{1+k_l}{2l+1}(\Delta\hat C\cos+\Delta\hat S\sin)$.
- **Elastic load Love numbers $k_l$** (Table 1; PREM, computed by D. Han as in Han & Wahr 1995):

  | $l$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | $k_l$ | +0.000 | +0.027 | −0.303 | −0.194 | −0.132 | −0.104 | −0.089 | −0.081 | −0.076 | −0.072 | −0.069 |

  | $l$ | 12 | 15 | 20 | 30 | 40 | 50 | 70 | 100 | 150 | 200 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | $k_l$ | −0.064 | −0.058 | −0.051 | −0.040 | −0.033 | −0.027 | −0.020 | −0.014 | −0.010 | −0.007 |

  Linear interpolation for other $l<200$ gives errors below 0.05%. The $l=1$ value assumes a centre-of-figure (CF) origin, $k_1=-(h_1+2l_1)/3$. In a centre-of-mass (CM) frame $k_1=-1$, and the GRACE degree-1 terms vanish. Anelastic effects perturb $(1+k_l)$ by at most a few percent. $k_0=0$.
- **Gaussian averaging (Jekeli 1981, normalised; eqs. 32–34):**
  $$W(\alpha)=\frac{b}{2\pi}\frac{\exp[-b(1-\cos\alpha)]}{1-e^{-2b}},\qquad b=\frac{\ln2}{1-\cos(r/a)}$$
  $$W_0=\frac{1}{2\pi},\quad W_1=\frac{1}{2\pi}\Big[\frac{1+e^{-2b}}{1-e^{-2b}}-\frac1b\Big],\quad W_{l+1}=-\frac{2l+1}{b}W_l+W_{l-1}$$
  Here $r$ is the averaging radius, the distance at which $W$ falls to half its central value.
- **Smoothed surface mass (form of eq. 35, without the land mask):**
  $$\bar{\Delta\sigma}(\theta,\phi)=\frac{2a\rho_{ave}\pi}{3}\sum_{l,m}\frac{2l+1}{1+k_l}W_l\,\tilde P_{lm}(\cos\theta)\big[\Delta C_{lm}\cos m\phi+\Delta S_{lm}\sin m\phi\big]$$
- **Error degree amplitudes (eqs. 23–26):**
  $$\delta N_l^{sp/t}=a\Big[\frac1T\sum_m\int_0^T(\delta C_{lm}^2+\delta S_{lm}^2)\,dt\Big]^{1/2}$$
  The space-and-time error variance is $\sum_l(\delta N_l)^2$. A per-coefficient error is $\delta N_l/\sqrt{2l+1}$. For an annual cosine/sine fit to $N$ samples per year, each coefficient's amplitude error is the per-sample error times $\sqrt{2/N}$. Averaging errors scale as $\sqrt{12/30}$ when moving from 12-day to 30-day solutions.

**gravity-toolkit implementation (Sutterley):**
- **Constants.** $G$ = 6.67430e-8 cgs; $GM_e$ = 3986004.418e14 cm³ s⁻² (WGS84); $GM$ without atmosphere = $GM_e-3.5\times10^{14}$; $a$ = 6.378137e8 cm; $f$ = 1/298.257223563.
- **Derived radius and density.** $R_e=a(1-f)^{1/3}$ and $\rho_e=\frac{3}{4}\frac{GM}{G\pi R_e^3}$ (about 5.51 g/cm³).
- **Degree factors from SH to spatial units:**
  - cm EWH: $\rho_eR_e(2l+1)/(3(1+k_l))$
  - mm geoid: $10R_e$
  - µGal: $10^6GM(l+1)/R_e^2$
  - mm elastic uplift: $10R_eh_l/(1+k_l)$
  - mbar: $g_{wmo}\rho_eR_e(2l+1)/(3\cdot10^3(1+k_l))$
  
  An optional ellipsoidal correction divides by $(1-f)$. All are cited to Wahr et al. (1998).
  
  — [gravity-toolkit units.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/units.py)
- **Gaussian weights in gravity-toolkit.** Uses $R_e$ = 6371.0 km, the same $W_0,W_1$ and recursion written as $W_l=\frac{1-2l}{b}W_{l-1}+W_{l-2}$, and is used as $2\pi W_l$. Once $W_l$ drops below `CUTOFF=1e-10`, the remaining tail is set to `CUTOFF`; this guards against the recursion's instability and underflow at high $l$ — [gravity-toolkit gauss_weights.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/gauss_weights.py)
- **Love numbers and frames in gravity-toolkit.** PREM load Love numbers (Dziewonski & Anderson 1981; Farrell 1972) are the default, with alternative sets selectable for uncertainty estimation. GRACE fields are in a CM frame. A CF option adjusts the Love numbers. SLR can replace $C_{20}, C_{21}, S_{21}, C_{22}, S_{22}, C_{30}, C_{40}, C_{50}$. Pole-tide and atmospheric corrections are needed before RL06 but not for RL06 — [gravity-toolkit Spatial Maps background](https://gravity-toolkit.readthedocs.io/en/latest/background/Spatial-Maps.html)

**Destriping (Swenson & Wahr 2006, GRL 33, L08402):**
- Stripes come from correlated coefficients. The method fits and removes smooth polynomials across degree for fixed order, separately for even and odd degrees — [UCAR record](https://impacts.ucar.edu/en/publications/post-processing-removal-of-correlated-errors-in-grace-data/)
- Limitations: signal distortion is poorly understood, and the method becomes inapplicable as $m\to l_{max}$ because too few coefficients remain to fit — [Springer ESIN 2018](https://link.springer.com/article/10.1007/s12145-018-0368-0)
- gravity-toolkit implementation: a moving-window **quadratic** least-squares fit over $2n_s+1$ points centred on each degree, with $n_s=\mathrm{round}(15\,e^{-m/10})$ (minimum 2). Even and odd degrees are smoothed separately, and the smoothed series is subtracted only for $m\ge5$. Defaults are $l_{max}$=60 — [gravity-toolkit destripe_harmonics.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/destripe_harmonics.py)

**DDK decorrelation (Kusche 2007; Kusche et al. 2009, J. Geod. 83:903–913):**
- Filter: $W_a=(E^{-1}+aS^{-1})^{-1}E^{-1}$, with $E$ the error covariance (GRACE normal matrix $E^{-1}$) and $S$ the signal covariance from a power law $l^{-p}$. It is equivalent to quadratic regularisation and gives a dense, anisotropic filter. An order-convolution (block-diagonal) simplification makes it practical.
- Parameters, all with $p=4$ (DDK1–3 from Kusche et al. 2009 Table 1; DDK4–8 are ICGEM's values). Equivalent Gaussian radii for DDK1–3 are given under two definitions.

  | Filter | $a$ | Gaussian radius (Kusche 2007 def.) | Gaussian radius ($\omega_l=\tfrac12$ def.) |
  |---|---|---|---|
  | DDK1 | 1e14 | 1350 km | 530 km |
  | DDK2 | 1e13 | 900 km | 340 km |
  | DDK3 | 1e12 | 660 km | 240 km |
  | DDK4 | 5e11 | | |
  | DDK5 | 1e11 | | |
  | DDK6 | 5e10 | | |
  | DDK7 | 1e10 | | |
  | DDK8 | 5e9 | | |

- Kusche et al. (2009) write the surface mass coefficients as $\sigma_{lmq}=\frac{2\rho_e}{3\rho_w(1+k_l)}(l+\tfrac12)\delta x_{lmq}$.

— [ICGEM Note on the DDK filters](https://icgem.gfz.de/docs/Note_on_the_DDK_filters.pdf). Filter files and code are at [github.com/strawpants/GRACE-filter](https://github.com/strawpants/GRACE-filter) (linked from the note).

**C20 and C30 (TN-14):**
- The GSFC file `TN-14_C30_C20_GSFC_SLR.txt` (Loomis et al. 2020) recommends replacing $C_{20}$ for the whole GRACE and GRACE-FO span, and $C_{30}$ **after August 2016**. That is when single-accelerometer and transplant data begin for GRACE; they cover all of GRACE-FO.
- GRACE $C_{20}$ shows a spurious 161-day signal. TN-14 replaces TN-11. A long-term SLR-only $C_{20}$ series, June 1976 to July 2026, is also published (Loomis et al. 2025).

— [NASA GSFC SLR page](https://earth.gsfc.nasa.gov/geo/data/slr)

**Degree 1 (TN-13):**
- The method is Swenson, Chambers & Wahr (2008), JGR 113, B08410. It combines GRACE degree ≥2 coefficients with ocean-model degree-1 terms in a matrix equation. Sun et al. (2016) extended it with a sea-level-equation-consistent ocean. The result is distributed as TN-13, one per centre (JPL, CSR, GFZ) from RL06 onward — [JPL geocenter page](https://grace.jpl.nasa.gov/data/get_data/geocenter)

**GIA:**
- GFZ processing subtracts the ICE-6G_D (VM5a) model (Peltier et al. 2018) from Level-2 Stokes coefficients — [EGU22 presentation](https://presentations.copernicus.org/EGU22/EGU22-6961_presentation-h141992.pdf)
- GIA rates of Stokes coefficients can be derived from GIA uplift models — [Remote Sensing 2020](https://www.mdpi.com/2072-4292/12/7/1209)
- The choice of GIA model is the primary driver of differences in mass trends — [Peltier et al. comment on Chambers et al. 2010](https://www.atmosp.physics.utoronto.ca/%7epeltier/pubs_recent/Peltier,%20W.R.%20%20et%20al.,%20Comment%20on%20Chambers%20et%20al%20(2010),%20JGR-Solid%20Earth%20117,%20B11403,%202012.pdf)
- Published ICE-6G_C coefficients have excess power at degree ≥90 — [ANU assessment](https://researchportalplus.anu.edu.au/en/publications/an-assessment-of-the-ice6gcvm5a-glacial-isostatic-adjustment-mode/)

**Leakage and gain factors (Landerer & Swenson 2012, WRR 48, W04531):**
- Filtering is destriping plus a 300 km Gaussian, at $l_{max}\approx60$ (≈330 km wavelength).
- Gain factor: $k=\arg\min_k M$ with $M=\sum(S_T-kS_F)^2$ (eq. 1). Here $S_T$ is the unfiltered model time series and $S_F$ the same model filtered identically to GRACE. For example, $k=1.44$ cut leakage variance about 85% for one basin.
- Gridded leakage error (eq. 2): $E^l_g=\mathrm{RMS}(S_T-kS_F)\cdot\frac{\mathrm{RMS}_{GRACE}}{\mathrm{RMS}_{model}}$. The total error is leakage and measurement errors in quadrature.
- Spatial error covariance (eq. 3): $\mathrm{Cov}(x_i,x_j)=\sigma_i\sigma_j\exp\!\big(-\frac{d_{ij}^2}{2d_0^2}\big)$, with $d_0$ = 300 km for measurement errors and 100 km for leakage.
- Regional variance (eq. 4): $\mathrm{var}=\sum_i\sum_jw_iw_j\mathrm{Cov}(x_i,x_j)$, with $w_i$ area weights.
- Gridded gain factors give accuracy comparable to basin gain factors. Frequency-dependent factors may be needed where seasonal and trend signals have different spatial scales.

— [Landerer & Swenson 2012 (JPL PDF)](https://grace.jpl.nasa.gov/system/downloadable_items/3_Landerer.2011WR011453.pdf)

**Mascons:**
- JPL RL05M/RL06 (Watkins et al. 2015) parameterises monthly gravity as **3° equal-area spherical-cap mascons**, with constraints from near-global geophysical models. A Coastline Resolution Improvement filter halves coastal leakage in synthetic tests — [Caltech record](https://resolver.caltech.edu/CaltechAUTHORS:20230307-21998000.14); [Wiese et al. 2016 WRR](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/grace/open/L1B/GFZ/AOD1B/RL04/docs/Wiese_2016_WRR.pdf)
- CSR (Save et al. 2016, JGR Solid Earth) uses Tikhonov regularisation on a 1° geodesic grid with GRACE-only, time-variable constraints, choosing the regularisation parameter by the L-ribbon method — [CSR RL05 mascons page](https://www2.csr.utexas.edu/grace/RL05_mascons.html); [GFZ GRACE refs](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2016JGRB..121.7547S.html)
- All these constraints are variants of Tikhonov regularisation — [TU Delft paper](https://repository.tudelft.nl/file/File_835c10e5-c109-4cd9-8093-1d58b9981819)

### Inferences
- **Monthly processing pipeline** (synthesis of the cited sources):
  1. Read GSM coefficients. Rescale to a common $GM$ and $R$ if centres differ (Barthelmes eq. 111).
  2. Replace $\bar C_{20}$ (and $\bar C_{30}$ after Aug 2016) with TN-14. Insert degree-1 from TN-13.
  3. Subtract a static mean field over a reference period, e.g. 2004.0–2009.999.
  4. Subtract GIA as $\dot C^{GIA}_{lm}(t-t_0)$.
  5. Optionally destripe or apply DDK.
  6. Convert with $\frac{a\rho_{ave}}{3\rho_w}\frac{2l+1}{1+k_l}\cdot2\pi W_l$.
  7. Synthesise on a grid.
  8. Apply gain factors.
  
  Keep $l=0$ out (set it to 0) unless you model total mass change explicitly; Wahr et al. drop $l=0$.
- **Basin average directly from coefficients.** I derived this from Wahr eq. 14 and the 4π orthogonality. Let $\vartheta(\theta,\phi)$ be the basin indicator function, with 4π coefficients $\vartheta^{c,s}_{lm}=\frac{1}{4\pi}\int\vartheta\,\bar P_{lm}\{\cos,\sin\}m\phi\,d\Omega$; then $\vartheta_{00}$ = basin area / $4\pi$. The basin-mean EWH is
  $$\overline{\mathrm{EWH}}=\frac{a\rho_{ave}}{3\rho_w\,\vartheta_{00}}\sum_{l,m}\frac{2l+1}{1+k_l}\tilde W_l\big(\vartheta^c_{lm}\Delta C_{lm}+\vartheta^s_{lm}\Delta S_{lm}\big)$$
  Here $\tilde W_l=2\pi W_l$, or 1 if unsmoothed. Its variance is $\mathbf g^\top\Sigma_{CS}\mathbf g$, with $\mathbf g$ the vector of the multipliers above and $\Sigma_{CS}$ the coefficient covariance. With only calibrated diagonal errors this reduces to $\sum g^2\sigma^2$. Validate against brute-force grid integration with $\cos\varphi$ area weights.
- **Trends and seasonal amplitudes.** Least-squares fit $x(t)=c_0+c_1(t-t_0)+a_1\cos\omega t+b_1\sin\omega t+a_2\cos2\omega t+b_2\sin2\omega t$, plus optional 161-day S2 or other alias terms. Amplitude is $A=\sqrt{a_1^2+b_1^2}$ and phase $\arctan2(b_1,a_1)$. The parameter covariance is $\hat\sigma^2(A^\top A)^{-1}$, or $(A^\top C^{-1}A)^{-1}$ with monthly errors. The amplitude error is $\sigma_A=\sqrt{(a_1^2\sigma_{a}^2+b_1^2\sigma_b^2+2a_1b_1\sigma_{ab})}/A$. This matches the annual-fit treatment in Wahr eqs. 25–26.
- **Spherical-cap mascon coefficients.** I derived this from Wahr eq. 12 plus the addition theorem; it is not from a fetched source, so validate it numerically. A cap of angular radius $\alpha$, centred at $(\theta_0,\lambda_0)$, with uniform EWH $h_w$ has
  $$\Delta C_{lm}=\frac{3\rho_w}{a\rho_{ave}}\frac{1+k_l}{2l+1}\cdot h_w\cdot\frac{P_{l-1}(\cos\alpha)-P_{l+1}(\cos\alpha)}{2(2l+1)}\,\tilde P_{lm}(\cos\theta_0)\cos m\lambda_0$$
  $\Delta S_{lm}$ is the same with $\sin m\lambda_0$. Here $P_l$ are unnormalised Legendre polynomials and $P_{-1}\equiv1$, so for $l=0$ the cap factor is $(1-\cos\alpha)/2$. The design matrix maps the mascon vector $\mathbf x$ to coefficients or to KBR range-rates. The regularised estimate is $\hat{\mathbf x}=(A^\top WA+\mu R)^{-1}A^\top W\mathbf y$, i.e. Tikhonov.
- **Gaussian weights in Rust.** The recursion is unstable as $l$ grows for small radii, because $W_l$ decays and cancellation amplifies error. Clamp as gravity-toolkit does, or compute $W_l$ by direct quadrature of $W(\alpha)$ with Gauss–Legendre for robustness.
- **gravity-toolkit's $\rho_e$ differs from Wahr's.** It is derived from $GM/(G\cdot\frac43\pi R_e^3)$, about 5.51 g/cm³, rather than Wahr's 5517 kg/m³. That is a ≈0.1% systematic difference in EWH; document the choice.

### Gaps
- Wang et al. (2012) load Love numbers (PREM/iasp91/ak135 with Crust 2.0) were not fetched; no table is given here. Han & Wahr (1995) values come via Wahr (1998) Table 1 only, to degree 200. Higher degrees need a downloadable table, e.g. gravity-toolkit's Love number files.
- I did not fetch TN-14/TN-13 file headers. The column layout, units, and whether TN-14 values are anomalies or full values, with which mean removed, are unverified. Check the file headers before coding.
- I did not retrieve the A et al. (2013) GIA paper or ICE-6G_D coefficient files. Units and reference epoch are unverified.
- Watkins et al. (2015) and Save et al. (2016) full texts were not read. Their mascon partials and regularisation matrices are not reproduced.
- The original Swenson & Wahr (2006) parameters were not read from the paper. The window $\approx30e^{-m/10}+1$ (gravity-toolkit's $2\cdot\mathrm{round}(15e^{-m/10})+1$) is consistent with commonly quoted $A=30, K=10$ values, but I did not confirm this from the original.

---

## 5. Mission error analysis: degree variances, degree RMS, error amplitudes, cumulative geoid error, Kaula's rule, commission vs omission, transfer factors

### Takeaway
Everything is built from per-degree sums of squared coefficients, multiplied by degree-dependent transfer factors. Using the cited Barthelmes, Wahr and gravity-toolkit formulas:

| Quantity | Transfer factor |
|---|---|
| Geoid | $R$ |
| Gravity disturbance | $\frac{GM}{R^2}(n+1)$ |
| Gravity anomaly | $\frac{GM}{R^2}(n-1)$ |
| Radial gradient $T_{rr}$ | $\frac{GM}{R^3}(n+1)(n+2)$ |
| EWH | $\frac{R\rho_{ave}}{3\rho_w}\frac{2n+1}{1+k_n}$ |

Kaula's rule, $\approx10^{-5}/n^2$ per coefficient, gives the signal model used for omission error.

### Cited Findings
- **Error degree amplitude (Wahr 1998, eq. 24).** $\delta N_l=a\big[\sum_m(\delta C^2_{lm}+\delta S^2_{lm})\big]^{1/2}$, time-averaged. The global error variance is $\sum_l(\delta N_l)^2$ (eq. 23), and the per-coefficient error is $\delta N_l/\sqrt{2l+1}$ — [Wahr et al. 1998](https://geodynamicsprogram.whoi.edu/wp-content/uploads/sites/14/2018/10/wahr_et_al.1998_117046.pdf)
- **Degree power and power per coefficient.** Degree power is $S(l)=\sum_{m}f_{lm}^2$ in 4π normalisation; power per coefficient is $S(l)/(2l+1)$; total energy is $4\pi S(l)$ — [SHTOOLS real SH](https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html)
- **Kaula's rule.** It bounds the RMS of the coefficients of each degree as $K/n^2$, with $K=10^{-5}$ for Earth — [AGU FM15 abstract](https://agu.confex.com/agu/fm15/webprogram/Paper80878.html). The corresponding degree variance is $(2n+1)(10^{-5}/n^2)^2$. It is too high below about degree 10, and altimetry-derived spectra deviate from it, falling below to about $n\sim60$ and above from 60 to 400 — [NASA NTRS Kaula 1964 review](https://ntrs.nasa.gov/api/citations/19640005753/downloads/19640005753.pdf); [NASA NTRS 1978](https://ntrs.nasa.gov/api/citations/19780023690/downloads/19780023690.pdf). The search summaries could not pin which NTRS document carries each statement.
- **Degree-dependent factors** are given by Barthelmes (rev. 2013, eqs. 118, 125, 126) and gravity-toolkit `units.py` (mmGH $=10R_e$; µGal $=10^6GM(l+1)/R_e^2$; cmwe $=\rho_eR_e(2l+1)/(3(1+k_l))$) — [STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf); [gravity-toolkit units.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/units.py)
- **Commission vs omission in practice.** Barthelmes notes that truncation at $\ell_{max}$ is a low-pass filter. The difference between the true potential and the series depends on both coefficient accuracy (commission) and $\ell_{max}$ (omission) — [STR09/02 rev. 2013 §4.1](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Measurement-error magnitudes.** RMS measurement error for 300 km-filtered monthly GRACE TWS reaches up to 36 mm EWH at low latitudes and falls below 15 mm poleward — [Landerer & Swenson 2012](https://grace.jpl.nasa.gov/system/downloadable_items/3_Landerer.2011WR011453.pdf)

### Inferences
- **Definitions to implement**, for degree $n$ with $\Delta$ standing for either signal or error coefficients:
  - Degree variance: $c_n=\sum_{m=0}^n(\bar C_{nm}^2+\bar S_{nm}^2)$
  - Degree RMS (per coefficient): $\sqrt{c_n/(2n+1)}$
  - Degree amplitude of functional $F$: $\sigma_n^F=\kappa_n^F\sqrt{c_n}$
  - Cumulative: $\Sigma^F(N)=\big[\sum_{n=2}^{N}(\kappa^F_n)^2c_n\big]^{1/2}$
  
  The transfer factors $\kappa_n^F$ at the reference sphere $r=R$ are:

  | $F$ | $\kappa_n^F$ |
  |---|---|
  | Geoid | $R$ |
  | $\delta g$ | $\frac{GM}{R^2}(n+1)$ |
  | $\Delta g$ | $\frac{GM}{R^2}(n-1)$ |
  | $T_{rr}$ | $\frac{GM}{R^3}(n+1)(n+2)$ |
  | EWH | $\frac{R\rho_{ave}}{3\rho_w}\frac{2n+1}{1+k_n}$ |
  | Surface density | $\frac{R\rho_{ave}}{3}\frac{2n+1}{1+k_n}$ |
  | Elastic uplift | $R\frac{h_n}{1+k_n}$ |

  At satellite altitude, multiply by $(R/r)^{n+1}$ for potential-type quantities, and include the corresponding powers of $R/r$ in the prefactors. Smoothing multiplies by $2\pi W_n$.
- **Cross-functional conversion.** Converting an error spectrum between functionals is the ratio of $\kappa$'s, e.g. $\sigma^{\Delta g}_n=\frac{GM}{R^3}(n-1)\sigma^{N}_n\approx\frac{\gamma}{R}(n-1)\sigma^N_n$.
- **Omission error from Kaula.** $\Sigma^N_{om}(N)\approx R\big[\sum_{n>N}(2n+1)10^{-10}n^{-4}\big]^{1/2}\approx R\cdot10^{-5}/N$, i.e. about 63.7 m$/N$. That gives about 1 m at $N=60$ and about 3 cm at $N=2190$. This is an order-of-magnitude estimate; Kaula over-predicts high-degree power for some bands per the cited NTRS note. A better approach is to fit a degree-variance model, e.g. Tscherning–Rapp, to a high-degree model such as EGM2008.
- **Commission error.** Use $\delta\bar C,\delta\bar S$ (calibrated formal errors, or month-to-month scatter after removing signal) in the same cumulative formula. Total error is $\sqrt{\text{commission}^2+\text{omission}^2}$ under independence.

### Gaps
- I did not verify the exact form of the Tscherning–Rapp (1974) degree-variance model, nor any normalisation-specific constants for Kaula's rule. Treat the factor $10^{-5}$ as approximate.

---

## 6. Fast spherical harmonic transforms (gridded synthesis and analysis)

### Takeaway
Separate variables. For each latitude ring, compute $A_m(\varphi)=\sum_n c_{nm}\bar P_{nm}(\sin\varphi)$ with the Section 2 recursions, which costs $O(L^3)$ for the whole grid. Then do one real FFT of length $\ge2L+1$ along longitude per ring, which costs $O(L^2\log L)$. Gauss–Legendre grids, $(L+1)\times(2L+1)$, and Driscoll–Healy grids, $N\times N$ or $N\times2N$ with $L=N/2-1$, both give **exact** quadrature for analysis. SHTns shows that well-vectorised $O(L^3)$ Gauss–Legendre transforms with on-the-fly ALFs beat asymptotically fast algorithms up to at least $L=1023$.

### Cited Findings
- **Grid formats** — [SHTOOLS grid formats](https://shtools.github.io/SHTOOLS/grid-formats.html):
  - **Gauss–Legendre (GLQ).** Latitudes are the $L+1$ zeros of the Legendre polynomial; quadrature is exact; longitudes are regular. Grid size is $(L+1)\times(2L+1)$, for the FFT in longitude.
  - **Driscoll–Healy.** The DH1 grid is $N\times N$ and DH2 is $N\times2N$, with $N$ equally spaced latitudes, an exact quadrature, and $L_{max}=N/2-1$. The 90°N band is down-weighted to zero, and the grid includes 90°N but not 90°S.
  - Reference: Driscoll & Healy 1994, Adv. Appl. Math. 15:202–250.
- **SHTns (Schaeffer 2013).** Explicitly vectorised Gauss–Legendre transforms that compute ALFs on the fly. This reduces memory and time, and the method is $O(N^3)$ yet outperforms fast algorithms up to $N=1023$ — [arXiv:1202.6522](https://arxiv.org/abs/1202.6522); [HAL](https://hal-insu.archives-ouvertes.fr/insu-00675145)
- **Analysis integral.** $f_{lm}=\frac{1}{4\pi}\int_\Omega f\,Y_{lm}\,d\Omega$ in 4π normalisation — [SHTOOLS real SH](https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html)
- **ICGEM approach.** ICGEM/Barthelmes transform gridded topography to SH using Sneeuw (1994), Geophys. J. Int. 118:707–716, which covers least-squares and quadrature analysis — [STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf)
- **Timing reference.** GeographicLib evaluates EGM2008 (degree 2190) at about 78 ms per point. With `GravityCircle`, which reuses work along a latitude circle, the cost drops to about 92 µs per point — [GeographicLib gravity docs](https://geographiclib.sourceforge.io/C++/doc/gravity.html)

### Inferences
- **Synthesis on a ring at latitude $\varphi_i$:**
  1. Compute $a_m=\sum_n(R/r)^n\bar C_{nm}\bar P_{nm}(t_i)$ and $b_m=\sum_n(R/r)^n\bar S_{nm}\bar P_{nm}(t_i)$.
  2. Form the complex spectrum $F_m=a_m-ib_m$ for $m\le L$, zero-padded to $N_\lambda$.
  3. Take a real inverse FFT, giving $f(\lambda_j)=\mathrm{Re}\sum_mF_me^{im\lambda_j}$.
  
  Use equatorial symmetry, $\bar P_{nm}(-t)=(-1)^{n+m}\bar P_{nm}(t)$, to compute both hemispheres from one recursion, which halves the cost.
- **Analysis on a GLQ grid.** With Gauss–Legendre nodes $t_i$ and weights $w_i$ on $[-1,1]$ ($\sum w_i=2$):
  $$\bar C_{nm}=\frac{1}{4\pi}\sum_iw_i\bar P_{nm}(t_i)\Big[\frac{2\pi}{N_\lambda}\sum_jf_{ij}\cos m\lambda_j\Big]$$
  $\bar S_{nm}$ is the same with $\sin$. The inner bracket is a forward real FFT. This is exact for band-limited $f$ up to degree $L$ when there are $L+1$ nodes and $N_\lambda\ge2L+1$.
- **Complexity.** Point-wise evaluation costs $O(L^2)$ per point, so a full grid point by point is $O(L^4)$. Ring-based evaluation is $O(L^3)$, which is the reason to reuse per-ring work, as GeographicLib's `GravityCircle` does. Fast Legendre transforms with $O(L^2\log^2L)$ complexity exist, but per SHTns are not worth it below $L\approx1000$.
- **For a Rust crate.** Compute GL nodes with Newton iteration on $P_{L+1}$, or the Bogaert (2014) $O(1)$ method. Use `rustfft`/`realfft` for the longitude step. Parallelise over latitude rings.

### Gaps
- I did not research the ducc0 internals or algorithm and found no citable source on it. I also did not verify the exact Driscoll–Healy quadrature weight formula, and I found no published benchmark comparing Rust SHT implementations.

---

## 7. Validation: test values and test suites

### Takeaway
There is no single official "EGM2008 test-point" table among the sources I retrieved. The practical validation set is a combination of:
- **Closed-form constants**: the GRS80 derived values, the IERS $\bar C_{20}$ zero-tide to tide-free value, and the Wahr Love numbers.
- **Analytic identities**: the ALF sum rule, orthogonality, the Laplace trace, and Gaussian $W_0=1/2\pi$.
- **Reference implementations**: GeographicLib, which agrees with NGA's harmonic synthesis code to 78 pm for EGM2008; SHTOOLS; gravity-toolkit tests; and ICGEM grids.
- **Barthelmes' published approximation-error statistics**, as order-of-magnitude checks.

### Cited Findings
- **GeographicLib vs NGA `hsynth_WGS84.f`.** `GravityModel::GeoidHeight` matches to 1.1 mm (EGM84), 23 nm (EGM96), and 78 pm for EGM2008 (23 pm after fixing issues in the NGA code). The Bruns approximation error is about 1.5 mm RMS (10 mm max) — [GeographicLib gravity docs](https://geographiclib.sourceforge.io/C++/doc/gravity.html)
- **Grid-interpolated EGM96 example.** `echo 16:46:33N 3:00:34W | GeoidEval` returns 28.7068 m, the EGM96 geoid at Timbuktu from the default `egm96-5` grid. Interpolation error is 0.140 m max (4.6 mm RMS) for bilinear and 3.2 mm max (0.7 mm RMS) for cubic — [GeoidEval man page](https://geographiclib.sourceforge.io/C++/doc/GeoidEval.1.html)
- **IERS $\bar C_{20}$ check.** Zero-tide $-0.48416948\times10^{-3}$ corresponds to tide-free $-0.48416531\times10^{-3}$. For EGM2008, the difference is $-4.1736\times10^{-9}$ — [IERS 2010 ch. 6](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)
- **GRS80 derived constants** to 10–12 significant digits ($U_0$, $J_4$, $J_6$, $J_8$, $\gamma_e$, $\gamma_p$, $e^2$, $f$, $m$, $k$; see Section 3) — [Moritz GRS80](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf); [Wikipedia GRS80](https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980)
- **gravity-toolkit tests.** The `test/` directory contains `test_legendre.py`, `test_harmonics.py`, `test_love_numbers.py`, `test_gaussian.py`, `test_units.py`, `test_geocenter.py`, `test_gia.py`, `test_point_masses.py`, `test_sea_level.py` and others. Reference outputs include `out.recurs.csv.gz`, `out.geoid.green_ice.0.5.2008.60.gz` and `land.fcn.1_deg.gz` — [gravity-toolkit test dir](https://github.com/tsutterley/gravity-toolkit/tree/main/test)
- **ICGEM calculation service.** It produces grids of functionals for any model, with a choice of reference system (radius, $GM$, flattening, $\omega$), tide system, zero-degree term, Gaussian filter and height. The form page currently says "Due to high load, no new calculations are accepted" — [ICGEM calcgrid](https://icgem.gfz.de/calcgrid); [ICGEM theory](https://icgem.gfz.de/theory)
- **Barthelmes' magnitude checks.** On EIGEN-GL04C and EIGEN-6C2:
  - $\Delta g_{cl}-\Delta g$: wrms 0.73 mGal (−12/+26)
  - $\Delta g_{sa}-\Delta g_{cl}$: wrms 0.11 mGal
  - $N_1^c-N^c$: wrms 1.3 mm
  - $\zeta_e$ vs geoid: up to 3.67 m in mountains
  
  — [STR09/02 (2009)](https://icgem.gfz.de/docs/str-0902.pdf)
- **Wahr 1998 Table 1 Love numbers** for unit-testing the interpolation of $k_l$ (Section 4) — [Wahr et al. 1998](https://geodynamicsprogram.whoi.edu/wp-content/uploads/sites/14/2018/10/wahr_et_al.1998_117046.pdf)

### Inferences
- **Recommended Rust test suite:**
  1. **ALF tests.** $\bar P_{nm}$ at a few $t$ values against closed forms: $\bar P_{10}=\sqrt3t$, $\bar P_{11}=\sqrt3u$, $\bar P_{20}=\frac{\sqrt5}{2}(3t^2-1)$, $\bar P_{21}=\sqrt{15}tu$, $\bar P_{22}=\frac{\sqrt{15}}{2}u^2$. Also the sum rule $\sum_m\bar P_{nm}^2=2n+1$ up to $n=2700$ at $\varphi$ = 0°, 45°, 89.9°, 89.9999°.
  2. **Derivative tests.** Finite differences, and the two derivative formulas against each other away from the poles.
  3. **GRS80.** Recompute every derived constant from the 4 defining ones and match Moritz's values.
  4. **Normal field consistency.** Evaluate $U$ via the SH series with $-J_{2n}/\sqrt{4n+1}$ at the equator and pole, and match $U_0$ and $\gamma_e$/$\gamma_p$.
  5. **Laplace trace.** The gradient tensor trace is ≈0 at random exterior points.
  6. **Round trip.** GLQ synthesis then analysis recovers random coefficients to ~1e-13.
  7. **EGM2008 geoid.** Heights at random points against GeographicLib (`GeoidHeight` or the `Gravity` tool with `-n egm2008`). This is the closest available stand-in for NGA test values.
  8. **Gravity-toolkit outputs.** EWH maps against `out.green_ice...` reference files, or against gravity-toolkit run on the same GSM file.
  9. **IERS $\bar C_{20}$** tide conversion.
- Of these, items 1–6 need no external data.

### Gaps
- I could not locate NGA's official EGM2008 test-point input and output files (`hsynth_WGS84` test data) or their numeric values; the NGA page only lists the zip archive.
- The contents of the pyshtools/SHTOOLS test suites were not inspected.
- The contents of gravity-toolkit's test assertions were not opened; the names only suggest what each tests.
