"""M2 — Lois de probabilité en fiabilité.

Notations du cours (Z. Klim) :
  exponentielle : lambda ; MTTF = 1/lambda
  normale       : N(m, sigma) ; u = (t - m)/sigma
  Weibull       : b (forme, aussi noté beta), theta (échelle, aussi noté eta), t0 (position)
  log-normale   : m, sigma de T  ->  sigma_y^2 = ln[(sigma/m)^2 + 1], mu_y = ln m - sigma_y^2/2
"""
import math

from .numerique import bissection, integrale_infinie, phi, phi_inv
from .rapport import Solution, fmt


# ---------------------------------------------------------------- objets lois

class Exponentielle:
    def __init__(self, lam):
        self.lam = lam
        self.nom = f"Exponentielle(λ={fmt(lam)})"

    def R(self, t):
        return math.exp(-self.lam * t)

    def f(self, t):
        return self.lam * math.exp(-self.lam * t)

    def h(self, t):
        return self.lam

    def H(self, t):
        return self.lam * t

    def mttf(self):
        return 1 / self.lam


class Weibull:
    def __init__(self, b, theta, t0=0.0):
        self.b, self.theta, self.t0 = b, theta, t0
        self.nom = f"Weibull(b={fmt(b)}, θ={fmt(theta)}" + (f", t0={fmt(t0)})" if t0 else ")")

    def _x(self, t):
        return max(t - self.t0, 0.0) / (self.theta - self.t0)

    def R(self, t):
        return math.exp(-self._x(t) ** self.b)

    def f(self, t):
        if t <= self.t0:
            return 0.0
        x = self._x(t)
        return self.b / (self.theta - self.t0) * x ** (self.b - 1) * math.exp(-x ** self.b)

    def h(self, t):
        if t <= self.t0:
            return 0.0
        return self.b / (self.theta - self.t0) * self._x(t) ** (self.b - 1)

    def H(self, t):
        return self._x(t) ** self.b

    def mttf(self):
        return self.t0 + (self.theta - self.t0) * math.gamma(1 + 1 / self.b)


class Normale:
    def __init__(self, m, sigma):
        self.m, self.sigma = m, sigma
        self.nom = f"Normale N({fmt(m)}; {fmt(sigma)})"

    def u(self, t):
        return (t - self.m) / self.sigma

    def R(self, t):
        return 1 - phi(self.u(t))

    def f(self, t):
        u = self.u(t)
        return math.exp(-u * u / 2) / (self.sigma * math.sqrt(2 * math.pi))

    def h(self, t):
        return self.f(t) / self.R(t)

    def H(self, t):
        return -math.log(self.R(t))

    def mttf(self):
        return self.m


class LogNormale:
    """Définie par la moyenne m et l'écart-type sigma de T (comme au cours),
    ou directement par mu_y et sigma_y (paramètres de Y = ln T)."""

    def __init__(self, m=None, sigma=None, mu_y=None, sigma_y=None):
        if mu_y is None:
            sigma_y = math.sqrt(math.log((sigma / m) ** 2 + 1))
            mu_y = math.log(m) - 0.5 * sigma_y ** 2
        self.m, self.sigma, self.mu_y, self.sigma_y = m, sigma, mu_y, sigma_y
        self.nom = f"Log-normale(μy={fmt(mu_y)}, σy={fmt(sigma_y)})"

    def u(self, t):
        return (math.log(t) - self.mu_y) / self.sigma_y

    def R(self, t):
        return 1 - phi(self.u(t)) if t > 0 else 1.0

    def f(self, t):
        if t <= 0:
            return 0.0
        u = self.u(t)
        return math.exp(-u * u / 2) / (self.sigma_y * t * math.sqrt(2 * math.pi))

    def h(self, t):
        return self.f(t) / self.R(t)

    def H(self, t):
        return -math.log(self.R(t))

    def mttf(self):
        return math.exp(self.mu_y + self.sigma_y ** 2 / 2)


# ---------------------------------------------------------------- utilitaires

def _liste(x):
    if x is None:
        return []
    return list(x) if isinstance(x, (list, tuple)) else [x]


def _temps_pour_R(loi, R_cible, t_max):
    return bissection(lambda t: loi.R(t) - R_cible, 0.0, t_max)


def _intervalle(sol, loi, t1, t2):
    sol.section(f"Probabilité de défaillance dans Δt = [{sol.q(t1, 'T')} ; {sol.q(t2, 'T')}]")
    R1, R2 = loi.R(t1), loi.R(t2)
    sol.etape(f"R(t1={fmt(t1)})", R1)
    sol.etape(f"R(t2={fmt(t2)})", R2)
    sol.etape("R(Δt) = R(t2)/R(t1)", R2 / R1, "R(Δt)")
    sol.etape("F(Δt) = 1 - R(t2)/R(t1)", 1 - R2 / R1, "F(Δt)")


# ---------------------------------------------------------------- solveurs

def exponentielle(lam=None, mttf=None, t=None, t1=None, t2=None, R_cible=None, unite_temps="h"):
    """Loi exponentielle. Donner lam OU mttf. t peut être une valeur ou une liste.
    unite_temps : unité des durées (h, j, an, km, cycle…) ; λ est alors en unite_temps⁻¹."""
    sol = Solution("LOI EXPONENTIELLE", unite_temps)
    if lam is None:
        sol.donnee(f"MTTF = {sol.q(mttf, 'T')}")
        lam = 1 / mttf
        sol.etape("λ = 1/MTTF", lam, "λ", "1/T")
    else:
        sol.donnee(f"λ = {sol.q(lam, '1/T')}")
        sol.resultats["λ"] = lam
        sol.unites["λ"] = sol.u("1/T")
    loi = Exponentielle(lam)
    sol.etape("MTTF = 1/λ", loi.mttf(), "MTTF", "T")
    for ti in _liste(t):
        sol.section(f"t = {sol.q(ti, 'T')}")
        sol.etape(f"λt = {fmt(lam)}·{fmt(ti)}", lam * ti)
        sol.etape("R(t) = exp(-λt)", loi.R(ti), f"R({fmt(ti)})")
        sol.etape("F(t) = 1 - exp(-λt)", 1 - loi.R(ti), f"F({fmt(ti)})")
        sol.etape("f(t) = λ·exp(-λt)", loi.f(ti), f"f({fmt(ti)})", "1/T")
        sol.etape("λ(t) = λ (constant)", lam, unite="1/T")
    if t1 is not None and t2 is not None:
        sol.section(f"Probabilité de défaillance dans Δt = [{sol.q(t1, 'T')} ; {sol.q(t2, 'T')}]")
        sol.etape(f"R(Δt) = exp[-λ(t2 - t1)] = exp[-{fmt(lam)}·{fmt(t2 - t1)}]",
                  math.exp(-lam * (t2 - t1)), "R(Δt)")
        sol.etape("F(Δt) = 1 - R(Δt)", 1 - math.exp(-lam * (t2 - t1)), "F(Δt)")
    for Rc in _liste(R_cible):
        sol.section(f"Temps pour lequel R(t) = {fmt(Rc)}")
        sol.etape(f"t = -ln(R)/λ = -ln({fmt(Rc)})/{fmt(lam)}", -math.log(Rc) / lam, f"t(R={fmt(Rc)})", "T")
    return sol


def normale(m, sigma, t=None, t1=None, t2=None, R_cible=None, unite_temps="h"):
    """Loi normale N(m, sigma). Utilise u = (t-m)/sigma et le tableau N1.
    unite_temps : unité de m, sigma et t (h, j, an, km, cycle…)."""
    sol = Solution(f"LOI NORMALE N({fmt(m)}; {fmt(sigma)})", unite_temps)
    sol.donnee(f"m = {sol.q(m, 'T')} ; σ = {sol.q(sigma, 'T')}")
    loi = Normale(m, sigma)
    sol.etape("MTTF = m", m, "MTTF", "T")
    for ti in _liste(t):
        u = loi.u(ti)
        sol.section(f"t = {sol.q(ti, 'T')}")
        sol.etape(f"u = (t - m)/σ = ({fmt(ti)} - {fmt(m)})/{fmt(sigma)}", u)
        if u < 0:
            sol.note(f"u < 0 : F(u) = 1 - F(|u|) = 1 - {fmt(phi(-u), 4)}")
        sol.etape(f"F(t) = F(u) [tableau N1]", phi(u), f"F({fmt(ti)})")
        sol.etape("R(t) = 1 - F(t)", loi.R(ti), f"R({fmt(ti)})")
        sol.etape("f(t) = exp(-u²/2)/(σ√(2π))", loi.f(ti), unite="1/T")
        sol.etape("λ(t) = f(t)/R(t)", loi.h(ti), f"λ({fmt(ti)})", "1/T")
    if t1 is not None and t2 is not None:
        _intervalle(sol, loi, t1, t2)
    for Rc in _liste(R_cible):
        sol.section(f"Temps pour lequel R(t) = {fmt(Rc)}")
        u = phi_inv(1 - Rc)
        sol.etape(f"F(u) = 1 - R = {fmt(1 - Rc)}  =>  u", u)
        sol.etape("t = m + u·σ", m + u * sigma, f"t(R={fmt(Rc)})", "T")
    return sol


def weibull(b=None, theta=None, t0=0.0, t=None, t1=None, t2=None, R_cible=None,
            beta=None, eta=None, unite_temps="h"):
    """Loi de Weibull (b = beta = forme, theta = eta = échelle, t0 = position).
    unite_temps : unité de θ, t0 et t (h, j, an, km, cycle…) ; b est sans unité."""
    b = b if b is not None else beta
    theta = theta if theta is not None else eta
    sol = Solution(f"LOI DE WEIBULL b={fmt(b)}, θ={fmt(theta)}" + (f", t0={fmt(t0)}" if t0 else ""), unite_temps)
    sol.donnee(f"b = {fmt(b)} (sans unité) ; θ = {sol.q(theta, 'T')}" + (f" ; t0 = {sol.q(t0, 'T')}" if t0 else ""))
    loi = Weibull(b, theta, t0)
    g = math.gamma(1 + 1 / b)
    sol.etape(f"Γ(1 + 1/b) = Γ({fmt(1 + 1 / b)})", g, "Γ(1+1/b)")
    if t0:
        sol.etape("MTTF = t0 + (θ - t0)·Γ(1 + 1/b)", loi.mttf(), "MTTF", "T")
    else:
        sol.etape("MTTF = θ·Γ(1 + 1/b)", loi.mttf(), "MTTF", "T")
    if b < 1:
        sol.note("b < 1 : taux de défaillance décroissant (jeunesse) -> maintenance préventive NUISIBLE")
    elif b == 1:
        sol.note("b = 1 : taux constant (loi exponentielle, λ = 1/θ) -> maintenance préventive INUTILE")
    else:
        sol.note("b > 1 : taux croissant (usure) -> maintenance préventive RAISONNABLE")
    for ti in _liste(t):
        sol.section(f"t = {sol.q(ti, 'T')}")
        x = loi._x(ti)
        sol.etape(f"(t - t0)/(θ - t0) = {fmt(ti - t0)}/{fmt(theta - t0)}", x)
        sol.etape("((t - t0)/(θ - t0))^b", x ** b)
        sol.etape("R(t) = exp[-((t - t0)/(θ - t0))^b]", loi.R(ti), f"R({fmt(ti)})")
        sol.etape("F(t) = 1 - R(t)", 1 - loi.R(ti), f"F({fmt(ti)})")
        sol.etape("f(t) = (b/θ)(t/θ)^(b-1)·exp[-(t/θ)^b]", loi.f(ti), unite="1/T")
        sol.etape("λ(t) = (b/θ)(t/θ)^(b-1)", loi.h(ti), f"λ({fmt(ti)})", "1/T")
    if t1 is not None and t2 is not None:
        _intervalle(sol, loi, t1, t2)
    for Rc in _liste(R_cible):
        sol.section(f"Temps pour lequel R(t) = {fmt(Rc)}")
        tR = t0 + (theta - t0) * (-math.log(Rc)) ** (1 / b)
        sol.etape(f"t = t0 + θ·[-ln R]^(1/b) = t0 + θ·[-ln {fmt(Rc)}]^(1/{fmt(b)})", tR, f"t(R={fmt(Rc)})", "T")
    return sol


def lognormale(m=None, sigma=None, t=None, t1=None, t2=None, R_cible=None,
               mu_y=None, sigma_y=None, unite_temps="h"):
    """Loi log-normale. Donner (m, sigma) de T comme au cours, ou (mu_y, sigma_y).
    unite_temps : unité de m, sigma et t ; μy = moyenne de ln t (t dans cette unité)."""
    sol = Solution("LOI LOG-NORMALE", unite_temps)
    if mu_y is None:
        sol.donnee(f"m = {sol.q(m, 'T')}, σ = {sol.q(sigma, 'T')} (moyenne et écart-type de T)")
        s2 = math.log((sigma / m) ** 2 + 1)
        sol.etape(f"σy² = ln[(σ/m)² + 1] = ln[({fmt(sigma)}/{fmt(m)})² + 1]", s2)
        sol.etape("σy", math.sqrt(s2), "σy")
        sol.etape(f"μy = ln m - σy²/2 = ln {fmt(m)} - 0,5·{fmt(s2)}", math.log(m) - s2 / 2, "μy")
    loi = LogNormale(m, sigma, mu_y, sigma_y)
    sol.etape("MTTF = exp(μy + σy²/2)", loi.mttf(), "MTTF", "T")
    for ti in _liste(t):
        u = loi.u(ti)
        sol.section(f"t = {sol.q(ti, 'T')}")
        sol.etape(f"u = (ln t - μy)/σy = (ln {fmt(ti)} - {fmt(loi.mu_y)})/{fmt(loi.sigma_y)}", u)
        sol.etape("F(t) = F(u) [tableau N1]", phi(u), f"F({fmt(ti)})")
        sol.etape("R(t) = 1 - F(t)", loi.R(ti), f"R({fmt(ti)})")
        sol.etape("λ(t) = f(t)/R(t)", loi.h(ti), f"λ({fmt(ti)})", "1/T")
    if t1 is not None and t2 is not None:
        _intervalle(sol, loi, t1, t2)
    for Rc in _liste(R_cible):
        sol.section(f"Temps pour lequel R(t) = {fmt(Rc)}")
        u = phi_inv(1 - Rc)
        sol.etape(f"u tel que F(u) = {fmt(1 - Rc)}", u)
        sol.etape("t = exp(μy + u·σy)", math.exp(loi.mu_y + u * loi.sigma_y), f"t(R={fmt(Rc)})", "T")
    return sol


def gamma(x):
    """Fonction Gamma (remplace l'Annexe G) avec la réduction Γ(x+1) = x·Γ(x)."""
    sol = Solution(f"FONCTION GAMMA Γ({fmt(x)})")
    y, facteurs = x, []
    while y > 2:
        y -= 1
        facteurs.append(y)
    if facteurs:
        sol.etape("Γ(x) = " + "·".join(fmt(f) for f in facteurs) + f"·Γ({fmt(y)})", math.gamma(x))
        sol.etape(f"Γ({fmt(y)}) [Annexe G]", math.gamma(y))
    sol.etape(f"Γ({fmt(x)})", math.gamma(x), "Γ")
    return sol


def mttf_numerique(R, echelle=1.0):
    """MTTF = ∫0^∞ R(t) dt pour une fonction R quelconque (callable)."""
    return integrale_infinie(R, 0.0, echelle)
