"""Exact (observer-free) energy conditions from the Hawking-Ellis type I eigen-decomposition, checked against sampled NEC."""
import numpy as np
import metric_to_energy_conditions as m

eta = np.diag([-1.0, 1, 1, 1])


def type_one(T):
    """Eigenvalues of T^a_b; for type I one eigenvector is timelike with eigenvalue -rho, the others give p_i."""
    w, V = np.linalg.eig(eta @ T)
    if np.max(np.abs(w.imag)) > 1e-9 * np.max(np.abs(w)):
        return None
    w, V = w.real, V.real
    norms = np.einsum("ai,ab,bi->i", V, eta, V)
    i = np.argmin(norms)
    if norms[i] >= 0:
        return None
    rho, p = -w[i], np.delete(w, i)
    return rho, p


for name, g, X in (("Alcubierre side wall", m.alcubierre(1.0, 1.0, 8.0), np.array([0, 0, 1.0, 0])),
                   ("Alcubierre front axis", m.alcubierre(1.0, 1.0, 8.0), np.array([0, 1.0, 1e-3, 0])),
                   ("Natario", m.shift_metric(m.natario_X(1.0, 1.0, 4.0)), np.array([0, 0.6, 0.85, 0.0])),
                   ("Morris-Thorne r=1.2", m.morris_thorne(1.0), np.array([0, 1.2 / np.sqrt(2), 1.2 / np.sqrt(2), 0.0]))):
    c = m.conditions(g, X)
    t1 = type_one(c["T"])
    if t1 is None:
        print(f"{name:24s} not type I; sampled NEC min {c['nec']:.4e}")
        continue
    rho, p = t1
    print(f"{name:24s} type I: rho_HE = {rho: .4e} (Eulerian {c['rho']: .4e}); min(rho+p_i) = {np.min(rho + p): .4e} "
          f"vs sampled NEC min {c['nec']: .4e}; DEC: rho-max|p| = {rho - np.max(np.abs(p)): .3e}; SEC: rho+sum p = {rho + p.sum(): .3e}")
