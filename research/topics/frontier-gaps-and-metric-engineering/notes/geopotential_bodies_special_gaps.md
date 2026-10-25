# Three gaps in `geopotential`: long thin prisms, the polyhedron's far field and mesh checks, and Carlson's principal values and incomplete integrals

Scope: `crates/geopotential/src/bodies.rs` (`prism`, `cubature`, `polyhedron`) and `crates/geopotential/src/special.rs` (`carlson_rf/rd/rj/rc`), as of commit 551424b. Background is the Phase 1 issue (FusionVault/physica#2, Items 1, 2 and 6).

**Where things are.** Everything is under `geopotential_gaps/` next to this file. This folder is a symlink to `research/topics/frontier-gaps-and-metric-engineering/notes/`.
- `scripts/` holds my own Python and the `.out` files they printed. They need numpy and scipy only.
- `downloads/<source>/` holds untrusted downloads, one folder per source. They were read, never run:
  - Carlson's arXiv TeX;
  - DLMF chapter 19 and §5.14 pages;
  - A&S page scans 613, 616, 618 and 625;
  - Choclo, Harmonica, polyhedral-gravity-model, Boost.Math and GSL sources;
  - the Eros mesh.

**How to read the values.** "Published" means quoted from a source I opened, with the equation or table named. "Recomputed" means my scripts produced it and it is not published; each such value is cross-checked by at least two independent methods.

## Gap 1 — Long thin prisms: what is the most accurate implementable method, and what tests it?

### Takeaway

**What is wrong.** The corner sums lose about 10–25 ε·r³/(abc) far away. Next to a very long prism they also lose about ε·L²/V, where L is the long side, even close in. Splitting a prism into near-cubes does not cure the far-field loss: each piece loses about ε·r² in absolute terms whatever its size.

**The fix.** The cure is a Gauss–Legendre cubature, with the order on each axis set from the Bernstein ellipse of 1/r:

  ρ_a = (d + √(d² + h_a²))/h_a,  n_a = ⌈ln(1/ε)/(2 ln ρ_a)⌉

Here d is the distance from p to the box and h_a is the half-width along axis a. The rule is:
- use the cubature when it is affordable;
- otherwise bisect an elongated prism along its longest side, only when p lies outside the box;
- otherwise use the corner sums.

**How well it works.** In a Python port this holds every case tested to ≤3×10⁻¹⁴ in V, ∇V and the Hessian, from 0.01 off the surface out to 100 Brillouin radii:
- a cube;
- a 1×1×100 needle;
- a 1×1×10⁴ needle;
- a 100×100×1 plate.

Today's code loses up to 5×10⁻⁸. This is the same idea as Fukushima's 2019 Taylor series, without its coefficient tables.

### Cited Findings

**The diagnosis is published.** Fukushima (2019) says:
- the exact prism formulae "suffer from round-off errors when the evaluation point is outside the Brillouin sphere";
- the error "grows cubically with respect to the distance";
- rewriting with addition theorems does not settle it.

He builds a 3-D Taylor series up to 16th order for a uniform prism. Its 8th-order truncation gives 8 digits beyond 2.7 Brillouin radii, at 11% of the cost of the exact formula. He also gives an adaptive procedure that combines the series orders with the exact formula. — [Fukushima 2019, GJI 220, 610, abstract via Crossref](https://doi.org/10.1093/gji/ggz449)

**The series and its stopping rule** (full text read through WebFetch; OUP blocks curl):
- The expansion is about the prism centre: V = GM Σₙ P₂ₙ/R^(4n+1), with P₂ₙ homogeneous of degree 2n and coefficients C_ijk polynomial in the half-sides a, b, c.
- The sum is truncated when |P₂ₙ|/R^(4n) < δ.
- If order 16 is still too large, the exact formula is used instead.
- For half-sides (2.5, 3.1, 5.0), the distances for 8 digits in V are R/R_B > 6400 at order 0, 17 at order 2, 10 at order 4, 4.2 at order 6, 2.7 at order 8, 2.2 at order 10 and 1.7 at orders 12–16.
- The double-precision exact formula gives 8 digits for R/R_B < 170.

— [Fukushima 2019, article page](https://academic.oup.com/gji/article/220/1/610/5583624)

**Published test values: Fukushima 2019, Table 4.** The prism has half-sides (2.5, 3.1, 5.0), the point is (30, 40, 50) from the centre, and G = ρ = 1. — [Fukushima 2019, Table 4](https://academic.oup.com/gji/article/220/1/610/5583624)
- V rigorous: +4.385413035403429E+00
- V at order 0: +4.384062043356595E+00
- V at order 8: +4.385413035403303E+00
- g_x: −2.638376163531892E−02
- g_y: −3.515460560153808E−02
- g_z: −4.380838400289370E−02
- Γ_xx: −4.019417444564296E−04
- Coefficients: C₂₀₀ = −3.685, C₀₂₀ = −2.005, C₀₀₂ = +5.69, C₄₀₀ = +9.5305741666666666667

**Recomputed: Table 4 reproduces.** — [prism_checks.out](geopotential_gaps/scripts/prism_checks.out)
- V: 4.3854130354034266 by 10-point split Gauss–Legendre, and 4.3854130354034275 by extended-precision corner sums.
- g and Γ_xx: all agree with the table to ≤1×10⁻¹⁵ relative.
- Order 0 is M/R = 310/√5000 = 4.3840620433565949.
- The quadrupole coefficients are exactly (2a²−b²−c²)/6 = −3.685, (2b²−a²−c²)/6 = −2.005 and (2c²−a²−b²)/6 = +5.69.
- Double-precision corner sums give 4.3854130354011431. That is 5×10⁻¹³ off at only R/R_B ≈ 11.

**Recomputed: the issue's numbers reproduce.** Errors are the maximum over 12 directions, at 10 body diagonals. — [prism_checks.out](geopotential_gaps/scripts/prism_checks.out)

| Prism | V | ∇V | Hessian | ε·r³/abc |
|---|---|---|---|---|
| 1×1×10 | 2.2×10⁻¹⁰ | 6.6×10⁻¹⁰ | 1.7×10⁻¹⁰ | 2.3×10⁻¹¹ |
| 1×1×100 needle | 2.0×10⁻⁸ | 4.3×10⁻⁸ | 1.0×10⁻⁸ | 2.2×10⁻⁹ |

So the loss is about 9–20 ε·r³/(abc).

**Recomputed: the loss also bites near the body.** At a point 2 units off the axis at mid-length, the corner sums lose more as the prism gets longer. — [prism_checks.out](geopotential_gaps/scripts/prism_checks.out)

| Length L | V | ∇V |
|---|---|---|
| 10 | 1×10⁻¹⁵ | 5×10⁻¹⁵ |
| 100 | 1.3×10⁻¹⁴ | 7×10⁻¹⁴ |
| 10⁴ | 8.5×10⁻¹³ | 3.1×10⁻¹¹ |

**Recomputed: comparing strategies on the 1×1×100 needle.** Each cell is the maximum relative error over V, ∇V and the Hessian and 12 directions, with the kernel evaluations in brackets. — [prism_checks.out](geopotential_gaps/scripts/prism_checks.out)

| Distance (diagonals) | Crate today | Near-cubes × crate | Split GL, ratio 8, 4³ | Per-axis Bernstein GL |
|---|---|---|---|---|
| 1 | 4×10⁻¹¹ (8) | ≤7×10⁻¹⁶ (3200) | 7×10⁻¹² | 7×10⁻¹⁵ (336) |
| 2 | 5×10⁻¹⁰ (8) | ≤7×10⁻¹⁶ (3200) | 8×10⁻¹² | 4×10⁻¹⁵ (90) |
| 5 | 7×10⁻⁹ (8) | ≤7×10⁻¹⁶ (3200) | 9×10⁻¹² | 4×10⁻¹⁶ (63) |
| 9.99 | 5×10⁻⁸ (8) | ≤7×10⁻¹⁶ (3200) | 8×10⁻¹² | 6×10⁻¹⁵ (45) |
| 10.01 | 8×10⁻¹² (64) | ≤7×10⁻¹⁶ (3200) | 8×10⁻¹² | 6×10⁻¹⁵ (45) |
| 50 | 4×10⁻¹⁶ (64) | ≤7×10⁻¹⁶ (3200) | 4×10⁻¹⁶ | 3×10⁻¹⁶ (16) |

**Recomputed: the proposed recursive algorithm** (tolerance 10⁻¹⁴). Accuracy ranges from 0.01 off the surface to 100 R_B. — [prism_recursive.out](geopotential_gaps/scripts/prism_recursive.out)

| Body | Worst error | Evaluations |
|---|---|---|
| Cube | ≤3×10⁻¹⁴ | ≤216 |
| 1×1×100 needle | ≤8×10⁻¹⁵ | ≤3.7×10³ near the surface, 352 at 1.5 R_B, 16 at 100 R_B |
| 100×100×1 plate | ≤3×10⁻¹⁴ | ≤4.6×10⁴ near the face |
| 1×1×10⁴ needle | ≤2×10⁻¹⁴ | ≤6.3×10³ near, 16–88 beyond 1.5 R_B |

For comparison, today's whole-box corner sums for the 1×1×10⁴ needle lose 1×10⁻¹¹ at a point 2 units off its surface.

**Recomputed: cost against tolerance.** — [prism_tol_cost.out](geopotential_gaps/scripts/prism_tol_cost.out)

| Case | Tolerance 10⁻¹⁴ | Tolerance 10⁻¹³ |
|---|---|---|
| Plate, 0.3 above the face | 1.0×10⁻¹⁴, 34031 evaluations | 2.4×10⁻¹⁴, 8 evaluations (no split) |
| 1×1×10⁴ needle, 2 off its side | 4.8×10⁻¹⁵, 6386 evaluations | 1.9×10⁻¹⁴, 2753 evaluations |
| 1×1×100 needle, 0.3 off its side | 1.4×10⁻¹⁵, 3419 evaluations | 3.6×10⁻¹⁵, 1123 evaluations |

**Recomputed regression values** (not published). Each is from 8- and 10-point split GL, which agree to ≤7×10⁻¹⁶; the crate's current error is in brackets. — [prism_test_values.out](geopotential_gaps/scripts/prism_test_values.out)
- **1×1×100 needle centred at the origin, p = (300, 400, 1200):**
  - V = 7.6952599930725910e-02 (1×10⁻¹⁴)
  - ∇V = (−1.3687958249365563e-05, −1.8250610999154124e-05, −5.4670865023521322e-05) (1×10⁻¹³)
  - T_zz = 7.0929393171032028e-08 and T_xy = 9.7553576758747929e-09 (6×10⁻¹³)
- **Same needle, p = (30, 40, 120):**
  - V = 7.9989392258983916e-01 (8×10⁻¹²)
  - ∇V = (−1.7475322889634304e-03, −2.3300430521503370e-03, −5.9813860121975910e-03) (4×10⁻¹¹)
  - T_zz = 7.9410558597455347e-05 and T_xy = 1.7803965137863311e-05
- **100×100×1 plate, p = (120, −90, 300):**
  - V = 2.9660346942256535e+01 (4×10⁻¹²)
  - ∇V = (−3.0763825064585165e-02, 2.3071490847281481e-02, −7.8619525049547487e-02)
  - T_zz = 3.6670482722705658e-04 and T_xy = −7.0133912185296821e-05
- **Identity check:** monopole plus quadrupole, V ≈ M/R + M Σ C_i X_i²/R⁵, leaves residuals of O((h/R)⁴), as it should: 1.6×10⁻⁷ for the needle at R = 1300.

**Fukushima 2020** (GJI 222, 1898) rewrites Bessel's exact prism formula with the addition theorems of atan and ln. It is 2.3–3.7× faster and gives "a slight increase in the computing precision". — [Fukushima 2020 abstract, Crossref](https://doi.org/10.1093/gji/ggaa240). The full text returned HTTP 403 and I could not open it.

**Harmonica and Choclo** (Fatiando a Terra) evaluate the Nagy et al. kernels at every distance.
- A "safe" ln is credited to Fukushima 2020:
  - ln((y²+z²)/(r−x)) for x < 0;
  - −ln(−2x) on the axis;
  - 0 at r = 0.
- A "safe" atan returns ±π/2 or 0 when the denominator is 0. The crate already does the same.
- There is no far-field switch.

— [choclo prism/_kernels.py](https://github.com/fatiando/choclo/blob/main/src/choclo/prism/_kernels.py); [choclo references.rst](https://github.com/fatiando/choclo/blob/main/doc/references.rst)

Harmonica's `prism_gravity` only checks singular points and warns. I found no issue on far-field precision in either repo (`gh search issues`). — [harmonica prisms/gravity.py](https://github.com/fatiando/harmonica/blob/main/src/harmonica/_forward/prisms/gravity.py)

**Adaptive discretisation, as Harmonica does it for tesseroids:**
- A 2-point Gauss–Legendre rule on each axis (`GLQ_DEGREES = (2, 2, 2)`).
- A stack-based split into halves along every axis where distance/L_axis < D.
- `DISTANCE_SIZE_RATII = {"potential": 1, "g_z": 2.5}`.

— [harmonica tesseroid_gravity.py](https://github.com/fatiando/harmonica/blob/main/src/harmonica/_forward/tesseroid_gravity.py); [_tesseroid_utils.py](https://github.com/fatiando/harmonica/blob/main/src/harmonica/_forward/_tesseroid_utils.py)

Uieda et al. (2016) found that a 0.1% maximum relative error needs D = 1 for the potential, 1.5 for the acceleration and 8 for the gradients. — [Uieda, Barbosa & Braitenberg 2016, Geophysics 81(5) F41, abstract](https://livrepository.liverpool.ac.uk/id/eprint/3057019). Harmonica's 2.5 for g_z differs from the paper's 1.5.

### Inferences

**Why splitting into near-cubes alone is wrong for the far field.** The absolute corner-sum error is ε times the largest corner term (≈ r²·ln r), whatever the piece's size. So N pieces lose N·ε·r² against one piece's ε·r². Splitting helps only for points near a long side, where it shrinks the corner distances.

**Why the cubature is the right far field.** An n-point Gauss–Legendre rule on each axis integrates every moment up to degree 2n−1 exactly. A tensor cubature therefore *is* a multipole expansion to that degree plus a convergent remainder, i.e. Fukushima's series without the C_ijk tables. Its error goes like ρ^(−2n), with ρ the Bernstein-ellipse parameter of the nearest complex singularity of 1/r along each axis.

**Recommended algorithm** for `bodies::prism`, about 60–80 new lines that replace `FAR` and the fixed 4×4×4 `cubature`:
1. Compute d, the distance from p to the closed box. If d > 0, set n_a = max(1, ⌈ln(1/ε)/(2 ln ρ_a)⌉) with ρ_a = (d + √(d²+h_a²))/h_a.
2. If every n_a ≤ 16 and n_x·n_y·n_z ≤ 512, return a tensor Gauss–Legendre cubature with those orders. Get the nodes from the existing `crate::gauss_legendre(n)`, which is O(n²) per call; cache them only if profiling asks.
3. Otherwise compute the corner sums. If d > 0, the box is elongated (L_max > 2 L_min), and ε·R_max²(π/2 + |ln R_max|) > tol·|V| (R_max being the farthest corner), bisect the longest side and recurse.
4. Otherwise return the corner sums.

Use tol = 10⁻¹³ by default: plates then cost 8 evaluations and needles about 10³, at ≤2×10⁻¹⁴.

**Risks.**
- **Split planes through p.** Bisecting when p lies inside or on the box would put p on an internal face. Each half would then return its *outside* Hessian limit, and the sum would be off by 4π in T_nn. The d > 0 guard avoids this. Inside points keep today's whole-box corner sums.
- **The Bernstein estimate is heuristic.** It is conservative because it uses the distance to the box, not to the axis. It was validated only on the shapes above.
- **The Hessian's kernel 1/r³ has a larger constant** than V's, so check T separately. The tests here found T ≤ 8×10⁻¹⁵.
- **The switch is not perfectly continuous:** the field can jump at the ≤10⁻¹⁴ level where the method changes.

**Tests to change or add.**
- Change `far_field`: it pins the 10-diagonal hand-over at 10⁻⁹.
- Change `the_cubature_is_gauss_legendre`: it pins the 4-node rule.
- Add Fukushima's Table 4.
- Add the needle and plate regression values.
- Add additivity: prism(box) = prism(lower half) + prism(upper half), at points off the split plane.
- Add a check that the corner-sum and cubature branches agree at the switch distance.

### Gaps

- I could not open the full text of Fukushima 2020 (403). Its abstract is all I have.
- I did not extract the Fukushima 2019 coefficient tables beyond C₂₀₀, C₀₂₀, C₀₀₂ and C₄₀₀. The values above come from WebFetch's reading of the page. I re-derived and checked all of Table 4 except C₄₀₀ and the order-8 value.
- Fukushima 2018 (GJI 215, 864; a recursive formula for vertically polynomial density) and Fukushima 2017 (AJ 154, 145; a general finite body and Eros) were not opened; Crossref gives no abstracts.
- The Bernstein-ellipse order rule rests on standard Gauss quadrature theory. I did not open a primary source for its constant. It is validated numerically here only.

## Gap 2 — The polyhedron: far-field method, mesh checks and singular points

### Takeaway

**Far field.** Hand over to a cubature over the signed tetrahedra. Each tetrahedron joins the vertex centroid to one face and is collapsed onto a cube (the Duffy map). The rule is Gauss–Jacobi in u and v and Gauss–Legendre in w. The hand-over is at 5 circumradii with n = 5 (125 points per face), or at 10 with n = 4 (64 per face):
- On a cube mesh and on the 14744-face Eros mesh this gives ≤6×10⁻¹⁴ in V.
- The analytic sums lose up to 3×10⁻⁸ at 200 radii and 3×10⁻⁶ at 1000.
- Gauss–Jacobi with n points integrates every polynomial of degree ≤2n−1 over a tetrahedron exactly. The same rule therefore yields exact multipole or spherical-harmonic coefficients, which is the content of Werner 1997, without his recursion.

**Mesh checks.** Use an O(F log F) directed-edge test plus the sign of the volume, not polyhedral-gravity-model's O(F²) ray casting:
- every directed edge appears exactly once, and so does its reverse;
- the signed volume is positive.

This catches a flipped face, an open mesh and an inverted mesh.

### Cited Findings

**polyhedral-gravity-model's orientation check.** — [mesh_integrity.rst](https://github.com/esa/polyhedral-gravity-model/blob/main/docs/background/mesh_integrity.rst); [Polyhedron.cpp](https://github.com/esa/polyhedral-gravity-model/blob/main/src/polyhedralGravity/model/Polyhedron.cpp)
- For each face it casts a ray from the centroid along the normal, starting 10⁻¹⁴ off the face.
- It counts the faces the ray crosses (Möller–Trumbore), de-duplicating crossing points through a `std::set`.
- An odd count means the normal points inwards. The docs say the naive version is O(n²).
- The majority orientation wins.
- Its options are `AUTOMATIC` (the default; it warns about the O(n²) cost), `VERIFY` (throws), `HEAL` (swaps vertices 0 and 1 of the violating faces) and `DISABLE`.
- It also rejects any triangle with zero area and any out-of-range index.

**Its singularities and tolerances** follow Tsoulis 2012. — [GravityModelDetail.cpp](https://github.com/esa/polyhedral-gravity-model/blob/main/src/polyhedralGravity/model/GravityModelDetail.cpp); [UtilityFloatArithmetic.h](https://github.com/esa/polyhedral-gravity-model/blob/main/src/polyhedralGravity/util/UtilityFloatArithmetic.h)
- **Singularity terms**, by where P′ (p projected on the face's plane) lies:
  - inside the face: −2π·h_p;
  - on a segment: −π·h_p;
  - at a vertex: −θ·h_p, with θ the interior angle;
  - outside the face: 0.
- **LN** is set to 0 when P′ sits on a vertex of the segment's line, or when |s₁+s₂| and |l₁+l₂| are below the tolerance.
- **AN** is set to 0 when h_p or h_pq is below the tolerance.
- The tolerance is an absolute `EPSILON_ZERO_OFFSET = 1e-14`, not scaled to the mesh.
- There is no far-field treatment.

Tsoulis (2012) obtains the full tensor by turning volume integrals into surface and then line integrals. Where p projects onto a face's plane or a segment's line, extra "correction terms" are needed, and a Fortran program is supplied. — [Tsoulis 2012, Geophysics 77(2) F1, abstract](https://doi.org/10.1190/geo2010-0334.1)

Tsoulis & Gavriilidou (2021) review the line-integral approach: it "involves no approximations" but carries singularities at certain positions of the computation point, and comes with Matlab code. — [Geophysical Prospecting 69, 1745, abstract](https://doi.org/10.1111/1365-2478.13134)

**The library's published tests:**
- **Cube [−1,1]³** at ρ = 1 and G = 6.67430e-11, against an analytic prism solution:
  - at (−5,−5,−5): V = 6.16562723600017E-11 and accelerations of 4.11088282102974e-12 on each axis;
  - at (−5,−4.5,−5): V = 6.37068727767111E-11 and accelerations of (4.53491850099316e-12, 4.08129150394744e-12, 4.53491850099315e-12).
  — [analytic_cube_solution_density1.txt](https://github.com/esa/polyhedral-gravity-model/blob/main/test/resources/analytic_cube_solution_density1.txt); generator [analytical_cube_gravity.py](https://github.com/esa/polyhedral-gravity-model/blob/main/script/analytical_cube_gravity.py)
- **Tsoulis's test prism** (x ∈ [−20,0], y ∈ [0,10], z ∈ [15,25]) at the origin, checked against intermediate quantities from Tsoulis's Fortran to 10⁻⁶. — [GravityModelTest.cpp](https://github.com/esa/polyhedral-gravity-model/blob/main/test/model/GravityModelTest.cpp)
- **Eros:** the JOSS paper reports only timing: 253 µs per point for the 14744-face Eros mesh on an M1 Pro. — [paper.md](https://github.com/esa/polyhedral-gravity-model/blob/main/paper/paper.md); [JOSS 9(98) 6384](https://doi.org/10.21105/joss.06384)

**Recomputed: the library's cube rows reproduce.** — [polyhedron_checks.out](geopotential_gaps/scripts/polyhedron_checks.out)
- The crate's polyhedron gives G·V = 6.165627236000088e-11 and the prism 6.165627236000387e-11.
- A 10-point split GL reference gives 6.165627236000158e-11, so the file is right to 1.8×10⁻¹⁵ (and to −2.0×10⁻¹⁵ at the second point).
- The accelerations agree to all printed digits.

**Recomputed: the mesh checks work.** — [polyhedron_checks.out](geopotential_gaps/scripts/polyhedron_checks.out)
- The cube passes, with volume +1.
- The inverted cube is reported as volume −1, "flip all faces". It silently gives −V, and its trace inside is +4π.
- One flipped face is reported as "3 directed edges used twice". Silently, it gives a V that is 13.7% wrong at an inside point.
- One missing face is reported as "3 edges without a partner".
- The Eros mesh (polyhedral-gravity-model's copy) passes: 7374 of the 24235 listed vertices are used, there are 14744 faces, and the volume is +0.291330570919549. Its trace at the vertex centroid is −4π to 1.3×10⁻¹⁴.

**Recomputed: today's far-field loss.** Each figure is the maximum relative error over V, ∇V and T, in body sizes. — [polyhedron_checks.out](geopotential_gaps/scripts/polyhedron_checks.out)

| Distance | Cube mesh, crate | Cube mesh, log1p | Eros, crate | Eros, log1p |
|---|---|---|---|---|
| 10 | 3.0×10⁻¹³ | 1.1×10⁻¹³ | 8.0×10⁻¹¹ | 5.5×10⁻¹¹ |
| 100 | 2.2×10⁻¹⁰ | 1.4×10⁻¹¹ | 8.3×10⁻⁹ | 1.5×10⁻⁸ |
| 1000 | 2.7×10⁻⁷ | 1.6×10⁻⁸ | 3.2×10⁻⁶ | 1.0×10⁻⁶ |

The log1p column replaces L_e = ln((s+e)/gap) with log1p(2e/gap). It helps the cube about tenfold but does not help Eros consistently.

**Recomputed: the tetrahedral cubature.** Relative error in V; sums use `math.fsum` so that the test's own summation adds no error. — [polyhedron_handover.out](geopotential_gaps/scripts/polyhedron_handover.out)

| Mesh, distance in circumradii | Analytic | GJ n=4 | GJ n=5 | GJ n=6 |
|---|---|---|---|---|
| Cube, 5 | 7×10⁻¹⁵ | 2×10⁻¹¹ | 6×10⁻¹⁴ | 7×10⁻¹⁶ |
| Cube, 10 | 6×10⁻¹⁴ | 5×10⁻¹⁴ | 4×10⁻¹⁶ | 4×10⁻¹⁶ |
| Eros, 5 | 2×10⁻¹³ | 3×10⁻¹³ | 4×10⁻¹⁶ | 7×10⁻¹⁶ |
| Eros, 10 | 8×10⁻¹³ | 7×10⁻¹⁶ | ≤7×10⁻¹⁶ | ≤7×10⁻¹⁶ |
| Eros, 200 | 3×10⁻⁸ | ≤4×10⁻¹⁶ | ≤4×10⁻¹⁶ | ≤4×10⁻¹⁶ |

With naive (einsum) summation over the 10⁶ point terms, V, g and T plateau at about 2×10⁻¹³ instead ([polyhedron_checks.out](geopotential_gaps/scripts/polyhedron_checks.out)). Compensated summation across faces matters.

**Recomputed: exactness of the collapsed rule** — [duffy_exactness.out](geopotential_gaps/scripts/duffy_exactness.out)
- The map is x = A + u[(B−A) + v((C−B) + w(D−C))], with Jacobian 6·Vol·u²·v.
- With Gauss–Jacobi weights u² and v, n = 2, 3 and 4 are exact through total degree 3, 5 and 7, i.e. 2n−1.
- With plain Gauss–Legendre in u and v they are exact only through 2n−3.
- The exact moments come from the Dirichlet integral over the unit simplex, ∫x^a y^b z^c = a!b!c!/(a+b+c+3)!. — [DLMF 5.14.2](https://dlmf.nist.gov/5.14)

**Recomputed Eros regression values** (mesh units, unit Gρ; not published). Analytic sums against GJ n = 10. — [eros_values.out](geopotential_gaps/scripts/eros_values.out)

| p | V, analytic | V, GJ n=10 | Difference |
|---|---|---|---|
| (2, 0, 0) | 1.511930505573376e-01 | 1.511930505573354e-01 | 1×10⁻¹⁴ |
| (0, 0, 5) | 5.809882197502869e-02 | 5.809882197501222e-02 | 3×10⁻¹³ |
| (20, 30, 40) | 5.409795408354512e-03 | 5.409795405316109e-03 | 6×10⁻¹⁰ (the analytic loss) |

Werner (1997, Computers & Geosciences 23(10), 1071–1077) exists, per Crossref, but has no abstract there. — [Crossref record, DOI 10.1016/S0098-3004(97)00110-6](https://doi.org/10.1016/S0098-3004(97)00110-6)

### Inferences

**Far-field method for `bodies::polyhedron`.**
- **Switch rule:** if |p − c| > 5·R_c, with c the vertex centroid and R_c the largest vertex distance from c, use the tetrahedral Gauss–Jacobi rule with n = 5. Otherwise use today's edge/face sums.
  - c and R_c cost one O(V) pass, which is cheap next to the O(F) sums.
  - The same loop gives ∇V and the Hessian from the point kernels, as the prism's `cubature` already does.
- **Size:** about 40 lines, plus Gauss–Jacobi nodes for weights u² and v. Either:
  - hard-code the n = 5 tables, generating them once with a Newton iteration on Jacobi polynomials (about 25 lines; Golub–Welsch would need an eigen-solver), or
  - use plain Gauss–Legendre with n + 1 points (6³ = 216 per face) and skip Jacobi altogether.
- **Summation:** sum each face's points, then add across faces with Neumaier compensation, which the crate's tests already use.
- **Cost:** about 125 kernel evaluations per face against one log/atan face-and-edge evaluation, so roughly 5–10× today's cost per far point, still O(F).
- **Many far points:** for many points around one body (an orbit around Eros), a precomputed solid-harmonic expansion is far cheaper. Its exact coefficients are the same cubature applied to r^l Y_lm, a polynomial of degree l, with n ≥ (l+1)/2. That would need a stateful type (`Polyhedron::new(verts, faces)` caching c, R_c, the mesh check and the coefficients). Werner 1997 is the classical reference for this, but I could not open it.

**Mesh checks.** Add `bodies::check_mesh(verts, faces) -> Result<f64, MeshError>` returning the volume. It takes about 30 lines:
1. Reject repeated indices and zero-area faces.
2. Make 3F directed-edge keys i·V + j, then sort, or use a HashMap.
3. Any key that occurs twice means inconsistent orientation or a non-manifold edge.
4. Any key whose reverse is absent means an open mesh.
5. The volume is Σ v₀·(v₁×v₂)/6; negative means inward, so flip every face or negate the result.

**Wiring the mesh check in.**
- `polyhedron` is stateless, so calling `check_mesh` on every evaluation doubles the cost. Expose it and document it, with a `debug_assert!` in `polyhedron`. Or make it the constructor of a validated mesh type.
- A free run-time hint already exists: the trace is −Σω. It equals 0 or −4π, and on the surface −2π or a vertex's solid angle. A value far from these flags an open mesh at almost any p. An inverted closed mesh is visible only from inside (+4π).
- polyhedral-gravity-model's per-face ray parity costs O(F²): for Eros, about 2×10⁸ ray–triangle tests. It can also miscount when a ray grazes an edge or vertex; it de-duplicates crossing points only, and uses an absolute 10⁻¹⁴ tolerance.

**The Hessian on edges.** On a true edge, entries along the edge dyad E_e = n̂_A n̂_eᴬᵀ + n̂_B n̂_eᴮᵀ diverge like L_e, which is physical.
- **Inner diagonals** of a flat, triangulated quad are a defect, not physics. The two coplanar triangles carry equal and opposite dyads, so E_e ≡ 0 and the true Hessian is finite (the outside limit, as for the prism). The crate sums face by face, so on the diagonal itself it forms ∞ − ∞ = NaN.
- **Fix:** pair the two faces of each edge, which the directed-edge map supplies, and form E_e before multiplying by L_e. Skip the edge when E_e = 0, and return ±∞ only where E_e ≠ 0. This needs the edge map, i.e. the mesh-check pass. Otherwise keep the docs as they are.
- polyhedral-gravity-model instead sets LN = 0 at such points: a finite convention, not the limit.

**Small win:** computing L_e as log1p(2e/gap) is optional. It helps the cube far field tenfold but is no substitute for the cubature.

### Gaps

- **Werner 1997** (exact spherical-harmonic coefficients of a polyhedron): not opened. ScienceDirect was not tried, and the HAL line-integral paper on the same topic was blocked by Anubis.
- **No published final Eros numbers found:**
  - polyhedral-gravity-model tests only intermediate quantities on its "BigTest" mesh;
  - the JOSS paper gives only timings;
  - Fukushima 2017 (AJ 154, 145, a general body and Eros) was not opened.
- I did not test the Tsoulis/Fortran values for the test prism, since its final V, g and T are not in the test file.
- Cost ratios are estimates from operation counts. Nothing was timed in Rust.

## Gap 3 — Carlson's principal values and the incomplete Legendre integrals

### Takeaway

**The principal values, as one-line pre-transformations onto the existing positive-argument codes.**
- **R_C, y < 0:** use Carlson's (2.13a), R_C(x,−y) = √(x/(x+y))·R_C(x+y, y).
- **R_J, p < 0, either of two equivalent forms:**
  - Carlson's (2.23a), with y the middle value;
  - DLMF 19.20.14, with z the largest value, which Boost uses.

Both reproduce Carlson's check values to their 14 printed digits, and agree with scipy and with QUADPACK's Cauchy-weight quadrature to ≤7×10⁻¹⁶.

**The incomplete integrals F(φ|m), E(φ|m) and Π(n;φ|m) through DLMF 19.25:**
- Reduce φ by periodicity (DLMF 19.2.10).
- Compute Δ² = cos²φ + (1−m)sin²φ, not 1 − m sin²φ. Otherwise E by DLMF 19.25.10 loses 4×10⁻⁹ near m → 1, φ → π/2.
- Pick the E formula by the range of m.
- Take Π beyond its pole as a principal value, as DLMF does; Boost returns a domain error there instead.

All reproduce the A&S tables to their last digit.

**One correction to the literature.** Carlson's text puts the zero of R_J(2,3,4,p) at "approximately −1.2552". Three independent methods put it at p = −1.2521959347.

### Cited Findings

**Definitions and principal values in Carlson 1995** — [Carlson, arXiv:math/9409227 TeX](https://arxiv.org/abs/math/9409227)
- **R_J**, (1.2): (3/2)∫₀^∞ [(t+x)(t+y)(t+z)]^(−1/2)(t+p)^(−1) dt, with "the Cauchy principal value … taken if p is real and negative".
- **R_C**, (1.3): the same convention for y < 0.
- **R_C principal value**, (2.13a): R_C(x,−y) = (x/(x+y))^(1/2)·R_C(x+y, y), y > 0. It "vanishes if x=0".
- **R_J principal value**, (2.23a):
  - (y+q)R_J(x,y,z,−q) = (p−y)R_J(x,y,z,p) − 3R_F(x,y,z) + 3(xyz/(xz+pq))^(1/2)·R_C(xz+pq, pq), q > 0;
  - here p−y = (z−y)(y−x)/(y+q);
  - "if we permute the values of x,y,z so that (z−y)(y−x) ≥ 0, then p ≥ y > 0".
- **Check values** (§3):
  - R_C(1/4,−2) = ln2/3 = 0.23104 90601 8665;
  - R_J(2,3,4,−0.5) = 0.24723 81970 3052;
  - R_J(2,3,4,−5) = −0.12711 23004 2964;
  - R_C(9/4,2) = ln 2 = 0.69314 71805 5995;
  - R_C(0,1/4) = π.
- **The zero crossing:** "R_J vanishes when the fourth variable is approximately −1.2552, and cancellation between terms on the right side of (2.23a) will lead to loss of significant figures".
- **Legendre's integrals in Carlson's own sign convention for Π**, (4.5)–(4.7), with c = csc²φ:
  - F(φ,k) = sin φ·R_F(cos²φ, 1−k²sin²φ, 1);
  - E = R_F(c−1, c−k², c) − (k²/3)R_D(c−1, c−k², c);
  - Π(φ,k,n) = ∫(1 + n sin²θ)^(−1)(…)^(−1/2) dθ = R_F − (n/3)R_J(c−1, c−k², c, c+n).

**DLMF 19.20.14–17** — [DLMF §19.20](https://dlmf.nist.gov/19.20)
- With 0 ≤ x < y < z: (q+z)R_J(x,y,z,−q) = (p−z)R_J(x,y,z,p) − 3R_F + 3(xyz/(xy+pq))^(1/2)·R_C(xy+pq, pq), where p = (z(x+y+q) − xy)/(z+q), so that x < y < p < z.
- Complete case (x = 0): (q+z)R_J(0,y,z,−q) = (p−z)R_J(0,y,z,p) − 3R_F(0,y,z), with p = z(y+q)/(z+q).
- **Closed forms (19.20.8):**
  - R_J(0,y,y,−q) = −3π/(2√y(y+q));
  - R_J(x,y,y,p) = 3(R_C(x,y) − R_C(x,p))/(p−y).
- **Sign rule:** for real, distinct variables, R_J is a "circular" case if (p−x)(p−y)(p−z) > 0 and a "hyperbolic" case if it is < 0. The hyperbolic cases include the principal values.

**DLMF 19.6.15** gives R_C(0,y) = 0 for y < 0, and R_C(0,y) = ½πy^(−1/2) otherwise. — [DLMF §19.6](https://dlmf.nist.gov/19.6)

**DLMF 19.25** (corrected in v1.1.7 to 0 ≤ φ ≤ π/2, after Karney) — [DLMF §19.25](https://dlmf.nist.gov/19.25)
- **First kind**, 19.25.5: F(φ,k) = R_F(c−1, c−k², c).
- **Second kind, three forms:**
  - 19.25.9: E = R_F − (k²/3)R_D(c−1, c−k², c);
  - 19.25.10: E = k′²R_F(c−1, c−k², c) + (k²k′²/3)R_D(c−1, c, c−k²) + k²√(c−1)/(√c·√(c−k²)), for c > k²;
  - 19.25.11: E = −(k′²/3)R_D(c−k², c, c−1) + √(c−k²)/(√c·√(c−1)), for φ ≠ π/2.
  - "All terms … are nonnegative when k² ≤ 0, 0 ≤ k² ≤ 1, or 1 ≤ k² ≤ c, respectively."
- **Third kind**, 19.25.14: Π(φ,α²,k) − F = (α²/3)R_J(c−1, c−k², c, c−α²).
- **Third kind beyond the pole**, 19.25.16: "If α² > c, then the Cauchy principal value is" −(ω²/3)R_J(c−1, c−k², c, c−ω²) + √((c−1)(c−k²)/((α²−1)(1−ω²)))·R_C(c(α²−1)(1−ω²), (α²−c)(c−ω²)), with ω² = k²/α². This route needs no negative R_J argument.
- **Complete principal value**, 19.25.4: Π(α²,k) = −(k²/3α²)R_J(0, 1−k², 1, 1−k²/α²), for k² < 1 < α².

**DLMF 19.36(i)** — [DLMF §19.36](https://dlmf.nist.gov/19.36)
- "If (19.25.9) is used when 0 ≤ k² ≤ 1, cancellations may lead to loss of significant figures when k² is close to 1 and φ > π/4 … The cancellations can be eliminated, however, by using (19.25.10)."
- Principal values are transformed "by means of (19.20.14) and (19.2.20)" before duplication.

**Definitions and ranges** — [DLMF §19.2](https://dlmf.nist.gov/19.2); [DLMF §19.7](https://dlmf.nist.gov/19.7)
- **Π's convention**, 19.2.7: Π(φ,α²,k) = ∫₀^φ dθ/(√(1−k²sin²θ)(1−α²sin²θ)). The Cauchy principal value is taken where 1−α²sin²φ vanishes inside the path. This is the crate's sign, (n/3)R_J(…, 1−n).
- **Periodicity**, 19.2.10: F(mπ ± φ) = 2mK ± F(φ), E(mπ ± φ) = 2mE ± E(φ), D likewise.
- **Reciprocal modulus**, 19.7.4: F(φ,1/k) = kF(β,k), with sin β = sin φ/k ≤ 1.
- **Change of parameter**, 19.7.8: Π(φ,α²,k) + Π(φ,ω²,k) = F + √c·R_C((c−1)(c−k²), (c−α²)(c−ω²)), α²ω² = k². It "gives the Cauchy principal value … when α² > csc²φ".

**Special values** — [DLMF §19.6](https://dlmf.nist.gov/19.6)
- F(φ,1) = gd⁻¹(φ).
- E(φ,1) = sin φ.
- F(φ,0) = E(φ,0) = φ.
- Π(φ,1,0) = tan φ.
- Π(φ,α²,0) = R_C(c−1, c−α²).
- Π(α²,k) = K − Π(k²/α², k) for 1 < α² (19.6.5).

**A&S tables** (modular angle α, so m = sin²α). — [A&S p. 613](https://www.convertit.com/Go/ConvertIt/Reference/AMS55.ASP?Res=150&Page=613); [p. 616](https://www.convertit.com/Go/ConvertIt/Reference/AMS55.ASP?Res=150&Page=616); [p. 618](https://www.convertit.com/Go/ConvertIt/Reference/AMS55.ASP?Res=150&Page=618); [p. 625](https://www.convertit.com/Go/ConvertIt/Reference/AMS55.ASP?Res=150&Page=625)
- **Table 17.5, F(φ\α), p. 613:**
  - F(30°\45°) = 0.53562 273
  - F(30°\60°) = 0.54222 911
  - F(30°\90°) = 0.54930 614
  - F(20°\40°) = 0.35198 869
- **Table 17.6, E(φ\α), pp. 616 and 618:**
  - E(30°\45°) = 0.51204 932
  - E(30°\60°) = 0.50609 207
  - E(90°\45°) = 1.35064 388
  - E(70°\60°) = 1.02663 689
  - E(75°\85°) = 0.96992 212
  - E(80°\88°) = 0.98568 915
  - E(85°\88°) = 0.99748 392
- **Table 17.9, Π(n;φ\α), p. 625**, defined with (1 − n sin²θ)^(−1), the same sign as DLMF and the crate:
  - Π(0.5;60°\45°) = 1.38218
  - Π(0.3;75°\60°) = 1.93879
  - Π(0.6;45°\30°) = 0.91527
  - Π(0.5;90°\45°) = 2.70129
  - Π(0.1;90°\0°) = 1.65576
  - Π(0;75°\90°) = 2.02759

**Recomputed: Carlson's check values.** — [carlson_checks.out](geopotential_gaps/scripts/carlson_checks.out)

| Value | Ours | Error |
|---|---|---|
| R_C(1/4,−2) | 0.231049060186648 | 1.1×10⁻¹⁶ against ln2/3 |
| R_J(2,3,4,−0.5) | 0.247238197030516 | within the table's 14-digit rounding |
| R_J(2,3,4,−5) | −0.127112300429639 | within the table's 14-digit rounding |

**Recomputed: the two R_J forms agree with each other and with independent codes.** Carlson's (2.23a), DLMF 19.20.14, `scipy.special.elliprj` and QUADPACK `weight='cauchy'` agree to ≤6.7×10⁻¹⁶ at:
- (2,3,4,−0.5);
- (2,3,4,−5);
- (0.5,1,7,−3);
- (0.1,1,2,−0.7).

The DLMF 19.20.8 closed forms reproduce to ≈4×10⁻¹⁶.

**Recomputed: the zero of R_J(2,3,4,p)** is at p = −1.2521959347, by our root, scipy and QUADPACK. Carlson's −1.2552 does not match: at p = −1.2552, R_J = −5.011×10⁻⁴ by both our code and quadrature. Near the root our code and scipy differ by 5×10⁻¹¹ relative, which is the cancellation Carlson warns of.

**Recomputed: the incomplete integrals against A&S.**
- F and E match the 8-decimal tables to within table rounding (≤4.3×10⁻⁹), and scipy's Cephes `ellipkinc`/`ellipeinc` to round-off.
- Π matches the 5-decimal Table 17.9 within rounding (≤4.7×10⁻⁶), and adaptive quadrature to 10⁻¹².

**Recomputed: the special values.**
- F(φ,1) = atanh(sin φ), E(φ,1) = sin φ, F(φ,0) = E(φ,0) = φ and Π(φ,1,0) = tan φ.
- Π(φ,3,0) = R_C(c−1, c−3), a principal value.
- Periodicity at φ = π+0.4, 2π−0.3 and −π+1.0 agrees with 2jK + F and with scipy to round-off, and F and E are odd.
- m > 1 against reciprocal-modulus 19.7.4:
  - F(0.5|2.5) = 0.570176528700536 and E(0.5|2.5) = 0.444842124053142;
  - scipy returns NaN for both, although they are real.
- m = −3 and −50 agree with quadrature to 2×10⁻¹⁶.

**Recomputed: Π beyond its pole** (n sin²φ > 1). Our result is compared with DLMF 19.25.16; the spread is given in brackets.
- Π(2;1.2|0.5) = 0.253604270701507 (3×10⁻¹⁵)
- Π(4;0.9|0.3) = 0.259293626045840 (4×10⁻¹⁶)
- Π(1.5;1.4|0.9) = −0.973823998688881 (4×10⁻¹⁶)
- Π(10;0.5|0.2) = 0.233655250908165 (3×10⁻¹⁶)

QUADPACK's Cauchy weight failed for (1.5;1.4|0.9): it returned −0.6368 and warned "Extremely bad integrand behavior". Subtracting the singularity analytically instead gives −0.97382399868806 ([pv_check.py](geopotential_gaps/scripts/pv_check.py)).

The complete principal value matches DLMF 19.6.5 (Π(2|0.5) = −0.313544683465184 and Π(5|0.9) = −0.346655296948390), and periodicity holds with n > 1.

**Recomputed: E near m = 1.** Relative error against quadrature, from 19.25.10 with three ways of computing Δ²:

| Case | Δ² = 1 − m sin²φ | Δ² = cos²φ + (1−m)sin²φ | 19.25.9 |
|---|---|---|---|
| m = 1−10⁻¹⁰, φ = 1.5707 | 4.3×10⁻⁹ | 3×10⁻¹⁶ | 1.6×10⁻¹⁵ |

— [e_near_one.out](geopotential_gaps/scripts/e_near_one.out)

### Inferences

**Where it goes in `special.rs`, about 80 lines.**

*Principal values.*
- **`carlson_rc`:** add at the top `if y < 0.0 { return (x / (x - y)).sqrt() * carlson_rc(x - y, -y); }`.
- **`carlson_rj`:** replace the `p > 0` debug assert with `if p < 0.0` and this body:
  - sort x ≤ y ≤ z;
  - set q = −p and p′ = (z(x+y+q) − xy)/(z+q);
  - return ((p′−z)R_J(x,y,z,p′) − 3R_F + 3√(xyz/(xy+p′q))·R_C(xy+p′q, p′q))/(z+q).
  This is Boost's DLMF form. Carlson's middle-value form is equally exact here.
- **`cel` comment:** its doc says p > 0. With R_J for p < 0, `cel` could gain p < 0 by DLMF 19.25.19, cel = aR_F + ((b−pa)/3)·R_J(0, k_c², 1, p), but Bulirsch's own iteration does not cover it. Leave `cel` alone unless a caller needs it.

*Incomplete integrals,* new `pub fn` `elliptic_f(phi, m)`, `elliptic_e(phi, m)` and `elliptic_pi(n, phi, m)`, in parameter m to match `complete_elliptic`:
1. **Reduce φ:** j = round(φ/π) and r = φ − jπ, so |r| ≤ π/2. Then s = sin r, c = cos r.
2. **Δ²:** use c² + (1−m)s² for 0 ≤ m ≤ 1, and 1 − m s² otherwise.
3. **Where to return NaN:**
   - Δ² < 0, where the result is complex;
   - m ≥ 1 with j ≠ 0;
   - n s² = 1 exactly, where the result is infinite.
4. **F** = s·R_F(c², Δ², 1) + 2jK.
5. **E**, chosen by m:
   - m ≤ 0: s·R_F − (m/3)s³R_D(c², Δ², 1);
   - 0 < m < 1: (1−m)s·R_F(c², Δ², 1) + (m(1−m)/3)s³R_D(c², 1, Δ²) + m s c/Δ;
   - m = 1: 2j + sin r;
   - 1 < m: −((1−m)/3)s³R_D(Δ², 1, c²) + sΔ/c;
   - then add 2jE.
6. **Π** = s·R_F + (n/3)s³R_J(c², Δ², 1, 1 − n s²) + 2jΠ(n|m), with R_J's principal value when 1 − n s² < 0.
7. **For |φ| > 1/ε**, return 2φK/π, as Boost does.

The existing inline F in `gravitics::spacetime`'s Darwin test (`s2.sqrt() * carlson_rf(1.0 - s2, 1.0 - k2 * s2, 1.0)`) is a caller to switch to `elliptic_f`.

**Pitfalls to document.**
- **The sign of n.** Carlson's (4.7) and GSL use (1 + n sin²θ)^(−1). DLMF, A&S, Boost and the crate use (1 − n sin²θ)^(−1).
- **R_J changes sign** for p < 0 between circular and hyperbolic regions, and loses relative accuracy near its zero (5×10⁻¹¹ at 10⁻⁶ from the root).
- **m > 1** is real only for |sin φ| ≤ 1/√m, and periodicity then fails.
- **Π diverges logarithmically** at φ = asin(1/√n), and the complete Π at n = 1.

**Test values to pin.**
- Carlson's three principal values and R_C(0,−2) = 0.
- DLMF 19.20.8, both closed forms.
- The A&S entries above, at 10⁻⁸ for F and E and 10⁻⁵ for Π.
- The DLMF 19.6 special values.
- Periodicity, and oddness in φ.
- m > 1 against 19.7.4.
- E at m = 1−10⁻¹⁰, φ = 1.5707 against quadrature (≤3×10⁻¹⁶ with the cancellation-free Δ²), and E(φ,1) = sin φ at exactly m = 1.
- Π with n > 1 against DLMF 19.25.16 and 19.6.5.

### Gaps

- **Carlson's −1.2552** disagrees with three independent computations (−1.25220). I treat it as a typo in the arXiv text, but did not check the published version (Numer. Algorithms 10, 13, 1995).
- **The A&S values were read by eye** from 200-dpi page scans (NIST, public domain). The digits were cross-checked against computation, but a misread digit would show up only as a mismatch beyond table rounding. None did.
- **The scipy docs for elliprj/elliprc were not opened.** I used scipy 1.18.0's outputs as an independent code.

## Mature implementations compared (Harmonica/Choclo, polyhedral-gravity-model, Boost.Math, GSL, scipy)

### Takeaway

None of the mature codes treats the far-field cancellation of the analytic prism or polyhedron formulae:
- Harmonica/Choclo add only safe log/atan;
- polyhedral-gravity-model adds singularity terms and O(F²) orientation checks.

For elliptic integrals:
- Boost returns principal values for R_C and R_J, uses DLMF 19.25.10 with cancellation-free Δ², but refuses Legendre Π beyond its pole;
- GSL has no principal values and uses the cancelling 19.25.9 for E.

### Cited Findings

**Boost.Math `ellint_rc`** (y < 0) sets prefix = √(x/(x−y)), x ← x−y and y ← −y. It then uses closed forms (atan, or log1p differences), not duplication. — [Boost ellint_rc.hpp](https://github.com/boostorg/math/blob/develop/include/boost/math/special_functions/ellint_rc.hpp)

**Boost.Math `ellint_rj`** (p < 0) sorts x ≤ y ≤ z, then sets q = −p and p = (z(x+y+q) − xy)/(z+q). It returns ((p−z)R_J − 3R_F + 3√(xyz/(xy+pq))·R_C(xy+pq, pq))/(z+q), which is DLMF 19.20.14. — [Boost ellint_rj.hpp](https://github.com/boostorg/math/blob/develop/include/boost/math/special_functions/ellint_rj.hpp)

**Boost `ellint_1`** — [Boost ellint_1.hpp](https://github.com/boostorg/math/blob/develop/include/boost/math/special_functions/ellint_1.hpp)
- φ < 0 is handled by oddness.
- φ > 1/ε returns 2φK/π.
- Otherwise φ is reduced modulo π/2, with a sign flip on odd multiples.
- It applies DLMF 19.25.5 with c−1 = cot²φ, which is cancellation-free.
- c−k² is computed as (c−1) + (1−k²) when c/k² ∈ (0.9, 1.1).
- k²sin²φ ≥ 1 is a domain error.

**Boost `ellint_2`** — [Boost ellint_2.hpp](https://github.com/boostorg/math/blob/develop/include/boost/math/special_functions/ellint_2.hpp)
- It applies DLMF 19.25.10, with 1−k² = (1−|k|)(1+|k|) and c−k² = (c−1) + (1−k²), "calculated without cancellation".
- |k| = 1 is a special case: 2m + sin(φ − mπ).
- Small φ uses a series shortcut.

**Boost `ellint_3`** — [Boost ellint_3.hpp](https://github.com/boostorg/math/blob/develop/include/boost/math/special_functions/ellint_3.hpp)
- It raises a domain error "result is complex for v > 1 / sin^2(phi)".
- The complete Π requires v < 1.
- An A&S 17.7.7/8 route for v > 1 is kept in `#if 0`, "disabled … appears not to produce better error rates".

**GSL** (`specfunc/ellint.c`, read through the ampl/gsl mirror) — [GSL ellint.c](https://github.com/ampl/gsl/blob/master/specfunc/ellint.c)
- It reduces φ with nc = floor(φ/π + 0.5) and adds 2·nc times the complete integral.
- E uses sin φ·R_F − (k²/3)sin³φ·R_D, the form DLMF 19.36 warns about.
- P(φ,k,n) uses R_J(x, y, 1, 1 + n sin²φ), Carlson's sign of n.
- R_J raises a domain error for p below a tiny positive bound, so there is no principal value.

**scipy 1.18.0** — [carlson_checks.out](geopotential_gaps/scripts/carlson_checks.out)
- `elliprj` and `elliprc` return the principal values, matching Carlson's table.
- Cephes `ellipkinc`/`ellipeinc` return NaN for m > 1 even where the integral is real (φ = 0.5, m = 2.5).

### Inferences

**Copy Boost's choices with one change.**
- R_C and R_J principal values: DLMF 19.20.14 for R_J, Carlson 2.13a for R_C.
- E by 19.25.10, with Boost's cancellation-free Δ².
- The change: return the principal value for Legendre Π beyond its pole, as DLMF 19.2.7 defines it, rather than an error. Π is genuinely real there and DLMF 19.25.16 computes it cleanly.

**For the bodies** the crate would be ahead of every code surveyed:
- a far-field hand-over;
- an O(F log F) mesh check;
- finite Hessians on inner diagonals.

### Gaps

- **Harmonica's documentation pages** were not read; only the source code was. No statement on far-field precision was found in either.
- **The polyhedral-gravity-model student report** (TUM mediatum 1695208) was not opened.

---

**Not checked:**
- No Rust was compiled or run. All numbers come from Python ports of the crate's algorithms.
- Fukushima 2020 and Werner 1997 full texts were not opened.
- The Bernstein-ellipse order rule is validated numerically, not proven.
