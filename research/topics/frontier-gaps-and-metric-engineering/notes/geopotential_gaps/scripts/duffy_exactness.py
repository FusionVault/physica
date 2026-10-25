"""The collapsed (Duffy) rule over a tetrahedron, x = A + u[(B-A) + v((C-B) + w(D-C))]: with Gauss-Jacobi
nodes (weights u^2, v) and n points per direction it integrates every monomial of total degree <= 2n-1
exactly; with plain Gauss-Legendre only degree <= 2n-3. Exact values from the Dirichlet integral
(DLMF 5.14.2): over the unit simplex, int x^a y^b z^c = a! b! c! / (a+b+c+3)!."""
import math, itertools
import numpy as np
from polyhedron_checks import duffy_rule

A, B, C, D = (np.array(v, float) for v in ([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]))
vol6 = np.dot(B - A, np.cross(C - A, D - A))
for jacobi in (True, False):
    for n in (2, 3, 4):
        U, V, W, Wt = duffy_rule(n, jacobi)
        X = A + U[:, None] * ((B - A) + V[:, None] * ((C - B) + W[:, None] * (D - C)))
        exact_to = -1
        for deg in range(0, 2 * n + 2):
            ok = True
            for a, b in itertools.product(range(deg + 1), repeat=2):
                c = deg - a - b
                if c < 0:
                    continue
                q = abs(vol6) * np.sum(Wt * X[:, 0] ** a * X[:, 1] ** b * X[:, 2] ** c)
                ex = math.factorial(a) * math.factorial(b) * math.factorial(c) / math.factorial(a + b + c + 3)
                ok &= abs(q / ex - 1) < 1e-13
            if ok and exact_to == deg - 1:
                exact_to = deg
        print(f"{'Gauss-Jacobi' if jacobi else 'Gauss-Legendre'} n = {n}: exact through total degree {exact_to} (2n-1 = {2 * n - 1})")
