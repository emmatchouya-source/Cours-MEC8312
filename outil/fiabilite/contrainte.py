"""Théorie contrainte–résistance (chapitre à venir : formules classiques, à
confronter aux notes du professeur quand le chapitre sera couvert)."""
import math

from .numerique import phi, phi_inv, simpson
from .rapport import Solution, fmt


def normale_normale(mu_R, sigma_R, mu_S, sigma_S, R_cible=None):
    """Résistance R ~ N(μR, σR) et contrainte S ~ N(μS, σS) indépendantes :
    z = (μR - μS)/√(σR² + σS²),  fiabilité = Φ(z)."""
    sol = Solution("CONTRAINTE–RÉSISTANCE : deux lois normales")
    sol.donnee(f"Résistance N({fmt(mu_R)}; {fmt(sigma_R)}), contrainte N({fmt(mu_S)}; {fmt(sigma_S)})")
    s = math.sqrt(sigma_R ** 2 + sigma_S ** 2)
    sol.etape("√(σR² + σS²)", s)
    z = (mu_R - mu_S) / s
    sol.etape("z = (μR - μS)/√(σR² + σS²)", z, "z")
    sol.etape("Fiabilité = Φ(z) [tableau N1]", phi(z), "R")
    sol.etape("Probabilité de défaillance = 1 - Φ(z)", 1 - phi(z), "F")
    sol.etape("Coefficient de sécurité central μR/μS", mu_R / mu_S, "n")
    if R_cible is not None:
        zc = phi_inv(R_cible)
        sol.section(f"Pour une fiabilité de {fmt(R_cible)}")
        sol.etape("z requis", zc)
        sol.etape("μR requis = μS + z·√(σR² + σS²)", mu_S + zc * s, "μR requis")
    return sol


def general(loi_resistance, loi_contrainte, borne_inf, borne_sup, n=4000):
    """Cas de deux lois quelconques : R = ∫ f_S(s)·R_R(s) ds
    (R_R(s) = P(résistance > s)). Les lois sont des objets de fiabilite.lois."""
    sol = Solution("CONTRAINTE–RÉSISTANCE : lois quelconques (intégration numérique)")
    sol.donnee(f"Résistance : {loi_resistance.nom} ; contrainte : {loi_contrainte.nom}")
    val = simpson(lambda s: loi_contrainte.f(s) * loi_resistance.R(s), borne_inf, borne_sup, n)
    sol.etape("R = ∫ fS(s)·P(Résistance > s) ds", val, "R")
    sol.etape("F = 1 - R", 1 - val, "F")
    return sol
