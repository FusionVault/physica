"""Gap 1b: self-demagnetisation of a non-ellipsoidal body by cells, using only prism Hessians.

With magnetics::anomaly's convention B = (mu0/4pi)(T - I tr T) M, the H field of a uniformly
magnetised cell is H = T M / 4pi inside and out (B/mu0 - M inside, B/mu0 outside). A body of cubes
with susceptibility chi in H0 then satisfies M_i = chi (H0 + sum_j T_j(c_i) M_j / 4pi), a 3n x 3n
linear system (collocation at cell centres). Checks: the prism Hessian (sign fixed by trace -4pi
inside, 0 outside, and the point-mass far field); a staircase sphere of cubes against the exact
sphere M = chi H0 / (1 + chi/3); a single cube (N = 1/3 at its centre).
"""
import numpy as np

def prism_hessian(p, lo, hi):
    """Hessian of U = Int dV/r of the box [lo, hi] at points p (n, 3): shape (n, 3, 3). Principal arctan: valid off the face planes (collocation points here never lie on them)."""
    p = np.atleast_2d(p)
    T = np.zeros((len(p), 3, 3))
    for i, xi in enumerate((lo[0], hi[0])):
        for j, yj in enumerate((lo[1], hi[1])):
            for k, zk in enumerate((lo[2], hi[2])):
                sgn = (-1) ** (i + j + k + 1)  # +1 at the upper corner
                x, y, z = xi - p[:, 0], yj - p[:, 1], zk - p[:, 2]
                r = np.sqrt(x * x + y * y + z * z)
                T[:, 0, 0] += -sgn * np.arctan(y * z / (x * r))
                T[:, 1, 1] += -sgn * np.arctan(x * z / (y * r))
                T[:, 2, 2] += -sgn * np.arctan(x * y / (z * r))
                T[:, 0, 1] += sgn * np.log(z + r)
                T[:, 0, 2] += sgn * np.log(y + r)
                T[:, 1, 2] += sgn * np.log(x + r)
    T[:, 1, 0], T[:, 2, 0], T[:, 2, 1] = T[:, 0, 1], T[:, 0, 2], T[:, 1, 2]
    return T

# sign and normalisation checks
lo, hi = np.array([-0.5, -0.5, -0.5]), np.array([0.5, 0.5, 0.5])
c = prism_hessian(np.zeros(3), lo, hi)[0]
far = np.array([[7.0, -5.0, 11.0]]); Tf = prism_hessian(far, lo, hi)[0]
r = np.linalg.norm(far); pm = (3 * np.outer(far[0], far[0]) - r * r * np.eye(3)) / r ** 5
inside_trace = np.trace(prism_hessian(np.array([[0.2, -0.1, 0.3]]), lo, hi)[0])
print(f"cube centre -T/4pi diag {np.diag(-c / 4 / np.pi)}; trace inside {inside_trace:.12f} (-4pi = {-4*np.pi:.12f}); far field rel. err {np.abs(Tf - pm).max() / np.abs(pm).max():.1e}")

def solve(cells, s, chi, H0):
    n = len(cells)
    K = np.zeros((3 * n, 3 * n))
    for j, cj in enumerate(cells):
        T = prism_hessian(cells, cj - s / 2, cj + s / 2) / (4 * np.pi)
        K[:, 3 * j:3 * j + 3] = T.reshape(3 * n, 3)
    M = np.linalg.solve(np.eye(3 * n) - chi * K, np.tile(chi * H0, n))
    return M.reshape(n, 3)

H0 = np.array([0.0, 0.0, 1.0])
for chi in (0.1, 1.0, 10.0):
    M = solve(np.zeros((1, 3)), 1.0, chi, H0)
    print(f"single cube chi={chi}: M_z/(chi H0) = {M[0,2]/chi:.6f}  vs 1/(1+chi/3) = {1/(1+chi/3):.6f} (one-cell collocation = magnetometric N = 1/3)")
for R in (3, 5, 7):
    g = np.arange(-R, R) + 0.5
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    cells = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    cells = cells[np.linalg.norm(cells, axis=1) <= R]
    a_eq = (3 * len(cells) / (4 * np.pi)) ** (1 / 3)
    for chi in (1.0, 10.0):
        M = solve(cells, 1.0, chi, H0)
        print(f"staircase sphere R={R} ({len(cells)} cubes, equal-volume radius {a_eq:.3f}), chi={chi}: mean M_z/(chi H0) = {M[:,2].mean()/chi:.5f}; exact sphere {1/(1+chi/3):.5f}; spread {M[:,2].std()/M[:,2].mean():.3f}")
