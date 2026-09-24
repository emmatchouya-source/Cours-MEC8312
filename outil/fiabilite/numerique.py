"""Outils numériques en Python pur (aucune dépendance externe)."""
import math
from statistics import NormalDist

_N01 = NormalDist()


def phi(u):
    """F(u) de la loi normale centrée réduite (tableau N1)."""
    return _N01.cdf(u)


def phi_inv(p):
    """u tel que F(u) = p."""
    return _N01.inv_cdf(p)


def simpson(f, a, b, n=2000):
    """Intégrale de f sur [a, b] (Simpson composite)."""
    if a == b:
        return 0.0
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3


def integrale_infinie(f, a=0.0, echelle=1.0):
    """Intégrale de f sur [a, +inf) par le changement t = a + echelle*u/(1-u)."""
    def g(u):
        if u >= 1:
            return 0.0
        t = a + echelle * u / (1 - u)
        return f(t) * echelle / (1 - u) ** 2
    return simpson(g, 0.0, 1.0 - 1e-12, 20000)


def bissection(f, a, b, tol=1e-12, iter_max=500):
    """Racine de f sur [a, b] (f(a) et f(b) de signes opposés)."""
    fa, fb = f(a), f(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ValueError("Pas de changement de signe sur l'intervalle")
    for _ in range(iter_max):
        m = (a + b) / 2
        fm = f(m)
        if fm == 0 or (b - a) / 2 < tol * max(1.0, abs(m)):
            return m
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return (a + b) / 2


def minimum(f, a, b, n_grille=400, tol=1e-10):
    """Minimum global approché de f sur [a, b] : grille puis section dorée."""
    xs = [a + (b - a) * i / n_grille for i in range(n_grille + 1)]
    vals = [f(x) for x in xs]
    i = min(range(len(xs)), key=lambda k: vals[k])
    lo = xs[max(i - 1, 0)]
    hi = xs[min(i + 1, n_grille)]
    g = (math.sqrt(5) - 1) / 2
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    fc, fd = f(c), f(d)
    while hi - lo > tol * max(1.0, abs(c)):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - g * (hi - lo)
            fc = f(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + g * (hi - lo)
            fd = f(d)
    x = (lo + hi) / 2
    return x, f(x)


# ---------- algèbre linéaire minimale ----------

def resoudre(A, b):
    """Résout A x = b (élimination de Gauss avec pivot partiel)."""
    n = len(A)
    M = [list(map(float, A[i])) + [float(b[i])] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-300:
            raise ValueError("Système singulier")
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                k = M[r][c] / M[c][c]
                if k:
                    for j in range(c, n + 1):
                        M[r][j] -= k * M[c][j]
    return [M[i][n] / M[i][i] for i in range(n)]


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def expm(A):
    """Exponentielle de matrice (mise à l'échelle et élévation au carré + Taylor)."""
    n = len(A)
    norme = max(sum(abs(x) for x in ligne) for ligne in A) if n else 0
    s = max(0, int(math.ceil(math.log2(norme))) + 1) if norme > 0.5 else 0
    As = [[x / 2 ** s for x in ligne] for ligne in A]
    E = [[float(i == j) for j in range(n)] for i in range(n)]
    terme = [ligne[:] for ligne in E]
    for k in range(1, 30):
        terme = matmul(terme, As)
        terme = [[x / k for x in ligne] for ligne in terme]
        E = [[E[i][j] + terme[i][j] for j in range(n)] for i in range(n)]
    for _ in range(s):
        E = matmul(E, E)
    return E
