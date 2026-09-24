"""M3 — Maintenance préventive, remplacement, maintenance corrective,
stock de pièces de rechange, analyse ABC."""
import math

from .numerique import minimum, phi_inv, simpson
from .rapport import Solution, fmt


# ------------------------------------------------------- taux de défaillance

class Taux:
    """Taux de défaillance λ(t) avec son intégrale Λ(t) = ∫0^t λ(x)dx."""

    def __init__(self, f, Lambda=None, desc="λ(t)"):
        self.f = f
        self._L = Lambda
        self.desc = desc

    def __call__(self, t):
        return self.f(t)

    def Lambda(self, t):
        return self._L(t) if self._L else simpson(self.f, 0.0, t, 2000)


def taux_constant(lam):
    return Taux(lambda t: lam, lambda t: lam * t, f"λ(t) = {fmt(lam)}")


def taux_lineaire(a, b):
    """λ(t) = a + b·t  ->  Λ(t) = a·t + b·t²/2."""
    return Taux(lambda t: a + b * t, lambda t: a * t + b * t * t / 2,
                f"λ(t) = {fmt(a)} + {fmt(b)}·t")


def taux_weibull(beta, eta):
    """λ(t) = (β/η)(t/η)^(β-1)  ->  Λ(t) = (t/η)^β."""
    return Taux(lambda t: beta / eta * (t / eta) ** (beta - 1) if t > 0 else (0.0 if beta > 1 else float("inf")),
                lambda t: (t / eta) ** beta, f"Weibull β={fmt(beta)}, η={fmt(eta)}")


def taux_expression(expr):
    """λ(t) donné comme texte, ex. '5e-4 + 5e-5*t' ou '2e-3*t**0.5'."""
    code = compile(expr, "<λ(t)>", "eval")
    env = {k: getattr(math, k) for k in ("exp", "log", "sqrt", "sin", "cos", "pi", "e")}
    return Taux(lambda t: eval(code, {"__builtins__": {}}, dict(env, t=t)), None, f"λ(t) = {expr}")


# ------------------------------------------------------- préventive systématique

def preventive_systematique(taux, T, n=None, R_cible=None, n_max=1000):
    """Maintenance préventive idéale : n interventions uniformes sur (0, T),
    période θ = T/(n+1), Rm(T) = exp[-(n+1)·Λ(T/(n+1))].

    taux    : objet Taux (taux_lineaire(a, b), taux_weibull(β, η), taux_constant(λ),
              taux_expression("...")) ou fonction Python λ(t)
    n       : nombre de maintenances (entier ou liste)
    R_cible : fiabilité minimale à garantir -> n minimal
    """
    if not isinstance(taux, Taux):
        taux = Taux(taux)
    sol = Solution("MAINTENANCE PRÉVENTIVE SYSTÉMATIQUE (idéale)")
    sol.donnee(f"{taux.desc} ; T = {fmt(T)}")
    L = taux.Lambda(T)
    sol.etape(f"Sans maintenance : R(T) = exp[-∫0^T λ(t)dt] = exp(-{fmt(L)})", math.exp(-L), "R(T) sans maintenance")
    ns = [] if n is None else (list(n) if isinstance(n, (list, tuple)) else [n])
    for k in ns:
        th = T / (k + 1)
        Lk = taux.Lambda(th)
        sol.section(f"n = {k} maintenances")
        sol.etape("θ = T/(n+1)", th)
        sol.etape(f"∫0^θ λ(t)dt", Lk)
        sol.etape(f"Rm(T) = exp[-(n+1)·∫0^θ λ(t)dt] = exp(-{k + 1}·{fmt(Lk)})",
                  math.exp(-(k + 1) * Lk), f"Rm(T) n={k}")
    if R_cible is not None:
        sol.section(f"Nombre minimal de maintenances pour R ≥ {fmt(R_cible)}")
        lim = math.exp(-taux(0.0) * T) if math.isfinite(taux(0.0)) else 0.0
        sol.etape("Limite n→∞ : exp[-λ(0)·T]", lim)
        if lim < R_cible:
            sol.note("Impossible : même avec n → ∞ la fiabilité n'atteint pas la cible.")
        else:
            for k in range(0, n_max + 1):
                Rk = math.exp(-(k + 1) * taux.Lambda(T / (k + 1)))
                if Rk >= R_cible:
                    if k > 0:
                        Rp = math.exp(-k * taux.Lambda(T / k))
                        sol.etape(f"n = {k - 1} : Rm(T)", Rp)
                    sol.etape(f"n = {k} : Rm(T)", Rk)
                    sol.etape("n minimal", k, "n_min")
                    sol.etape("θ = T/(n+1)", T / (k + 1), "θ")
                    break
    return sol


# ------------------------------------------------------- remplacement

def remplacement_age(loi, Cp, Cd, theta_max=None, theta=None):
    """Remplacement préventif à âge fixé θ :
    C(θ) = [Cp·R(θ) + Cd·(1 - R(θ))] / ∫0^θ R(t)dt, minimisé en θ*.

    loi : objet loi (ex. lois.Weibull(b, theta)) ayant R(t) et mttf().
    """
    sol = Solution("REMPLACEMENT PRÉVENTIF À ÂGE FIXÉ")
    sol.donnee(f"{loi.nom} ; Cp = {fmt(Cp)} ; Cd = {fmt(Cd)}")

    def C(th):
        if th <= 0:
            return float("inf")
        return (Cp * loi.R(th) + Cd * (1 - loi.R(th))) / simpson(loi.R, 0.0, th, 400)

    mttf = loi.mttf()
    sol.etape("Politique corrective seule : C(∞) = Cd/MTTF", Cd / mttf, "C_corrective")
    for th in ([] if theta is None else (theta if isinstance(theta, (list, tuple)) else [theta])):
        sol.etape(f"C(θ = {fmt(th)})", C(th), f"C({fmt(th)})")
    th_max = theta_max or 3 * mttf
    th_opt, c_opt = minimum(C, th_max * 1e-4, th_max, 300)
    sol.etape("θ* (dC/dθ = 0)", th_opt, "θ*")
    sol.etape("C(θ*)", c_opt, "C(θ*)")
    if c_opt < Cd / mttf * 0.999 and th_opt < th_max * 0.99:
        sol.etape("Gain C(θ*)/C(∞)", c_opt / (Cd / mttf), "gain")
    else:
        sol.note("Pas de minimum intéressant : la maintenance corrective est préférable.")
    return sol


def periode_optimale_remplacement(CI, ice, icm, k1=None, k2=None):
    """Période optimale t* = sqrt(2·CI/(ice + icm)) ; CT* = k1 + k2 - (ice+icm)/2 + sqrt(2·CI·(ice+icm))."""
    sol = Solution("PÉRIODE OPTIMALE DE REMPLACEMENT (coût total annuel minimal)")
    s = ice + icm
    sol.etape(f"t* = √[2·CI/(ice + icm)] = √[2·{fmt(CI)}/{fmt(s)}]", math.sqrt(2 * CI / s), "t*", "ans")
    if k1 is not None and k2 is not None:
        sol.etape("CT* = k1 + k2 - (ice + icm)/2 + √[2·CI·(ice + icm)]",
                  k1 + k2 - s / 2 + math.sqrt(2 * CI * s), "CT*")
    return sol


def cout_total_annuel(CI, ice, icm, t, k1=0.0, k2=0.0):
    """CT(t) = k1 + k2 + CI/t + (t-1)/2·(ice + icm)."""
    return k1 + k2 + CI / t + (t - 1) / 2 * (ice + icm)


def _ratio_weibull(x, beta, r):
    if x <= 0:
        return float("inf")
    integ = simpson(lambda u: math.exp(-u ** beta), 0.0, x, 400)
    return (1 + (1 - math.exp(-x ** beta)) * r) / integ * math.gamma(1 + 1 / beta) / (1 + r)


def weibull_systematique(beta, eta, r=None, p=None, P=None):
    """Optimisation de la période d'intervention systématique (modèle de Weibull)
    — remplace l'ABAQUE : calcul exact du minimum de C2(x)/C1.

    beta, eta : paramètres de Weibull ; r = P/p (criticité économique) ou p et P.
    """
    if r is None:
        r = P / p
    sol = Solution("OPTIMISATION DE LA PÉRIODE SYSTÉMATIQUE — MODÈLE DE WEIBULL")
    sol.donnee(f"β = {fmt(beta)} ; η = {fmt(eta)} ; r = P/p = {fmt(r)}")
    g = math.gamma(1 + 1 / beta)
    sol.etape("m∞ = MTBF = η·Γ(1 + 1/β)", eta * g, "MTBF")
    sol.note("C2(x)/C1 = [1 + (1 - e^(-x^β))·r] / ∫0^x e^(-t^β)dt · Γ(1+1/β)/(1 + r),  x = θ/η")
    x0, ratio = minimum(lambda x: _ratio_weibull(x, beta, r), 1e-4, 3.0, 300)
    sol.etape("x0 (minimum de C2/C1)", x0, "x0")
    sol.etape("C2(x0)/C1", ratio, "C2/C1")
    if ratio < 1:
        sol.etape(f"θ0 = η·x0 = {fmt(eta)}·{fmt(x0)}", eta * x0, "θ0")
        sol.note(f"C2/C1 < 1 : la maintenance systématique est rentable (gain de {fmt((1 - ratio) * 100, 3)} %).")
    else:
        sol.note("C2/C1 ≥ 1 : pas de solution, rester en maintenance corrective.")
    return sol


# ------------------------------------------------------- corrective

def corrective(A_actuel, Cmv_actuel, MTTF, c, Cpr=0.0, decimales_A=None):
    """Maintenance corrective : disponibilité optimale et coût de main-d'œuvre.

    A_actuel    : A'∞, disponibilité actuelle
    Cmv_actuel  : C'mv, coût moyen actuel de la main-d'œuvre par réparation
    MTTF        : temps (ou km) moyen de bon fonctionnement
    c           : facteur de proportionnalité des pertes (coût par unité de temps d'arrêt)
    Cpr         : coût moyen des pièces de rechange (reste constant)
    decimales_A : arrondir Å∞ comme dans les corrigés (ex. 3 -> 0,997)
    """
    sol = Solution("MAINTENANCE CORRECTIVE — DISPONIBILITÉ OPTIMALE")
    sol.donnee(f"A'∞ = {fmt(A_actuel)} ; C'mv = {fmt(Cmv_actuel)} ; MTTF = {fmt(MTTF)} ; c = {fmt(c)}")
    K1 = Cmv_actuel * (1 - A_actuel) / A_actuel
    sol.etape(f"K1 = C'mv·(1 - A'∞)/A'∞ = {fmt(Cmv_actuel)}·(1 - {fmt(A_actuel)})/{fmt(A_actuel)}", K1, "K1")
    K2 = c * MTTF
    sol.etape(f"K2 = c·MTTF = {fmt(c)}·{fmt(MTTF)}", K2, "K2")
    kappa = K2 / K1
    sol.etape("κ = K2/K1", kappa, "κ")
    A_opt = math.sqrt(kappa) / (1 + math.sqrt(kappa))
    sol.etape("Å∞ = √κ/(1 + √κ)", A_opt, "Å∞")
    if decimales_A is not None:
        A_opt = round(A_opt, decimales_A)
        sol.etape(f"Å∞ arrondi ({decimales_A} décimales)", A_opt)
    Cmv = K1 * A_opt / (1 - A_opt)
    sol.etape("Cmv = K1·Å∞/(1 - Å∞)", Cmv, "Cmv")
    if A_opt > A_actuel:
        sol.note("Å∞ > A'∞ : INTENSIFIER la maintenance (augmenter Cmv).")
    elif A_opt < A_actuel:
        sol.note("Å∞ < A'∞ : RÉDUIRE l'effort de maintenance (diminuer Cmv).")
    else:
        sol.note("Å∞ = A'∞ : la disponibilité actuelle est déjà optimale.")
    sol.section("Coûts totaux des charges d'exploitation (par réparation)")
    Cp_act = c * MTTF * (1 / A_actuel - 1)
    Cp_opt = c * MTTF * (1 / A_opt - 1)
    CT_act = Cmv_actuel + Cpr + Cp_act
    CT_opt = Cmv + Cpr + Cp_opt
    sol.etape("Pertes actuelles Cp = c·MTTF·(1/A'∞ - 1)", Cp_act, "Cp actuel")
    suffixe = "" if Cpr else "  (+ Cpr, non donné)"
    sol.etape("CT actuel = C'mv + Cp" + (" + Cpr" if Cpr else "") + suffixe, CT_act, "CT actuel")
    sol.etape("Pertes après modification = c·MTTF·(1/Å∞ - 1)", Cp_opt, "Cp optimal")
    sol.etape("CT après modification = Cmv + Cp" + (" + Cpr" if Cpr else "") + suffixe, CT_opt, "CT optimal")
    sol.etape("ΔCI = Cmv - C'mv (investissement)", Cmv - Cmv_actuel, "ΔCI")
    sol.etape("ΔCEX = CT actuel - CT optimal (diminution des charges)", CT_act - CT_opt, "ΔCEX")
    sol.etape("ΔCP = diminution des pertes", Cp_act - Cp_opt, "ΔCP")
    return sol


# ------------------------------------------------------- stock

def stock(lam, T, p0, nb_elements=1):
    """Nombre de pièces de rechange pour avoir un stock suffisant avec la probabilité p0.

    lam : taux de défaillance d'un élément ; T : durée de la mission ;
    nb_elements : nombre d'éléments identiques en service (système série).
    """
    sol = Solution("STOCK DE PIÈCES DE RECHANGE")
    lT = lam * T * nb_elements
    sol.etape("λT" + (f" (× {nb_elements} éléments)" if nb_elements != 1 else ""), lT, "λT")
    sol.section("Formule approximative (loi normale)")
    u = phi_inv(p0)
    sol.etape(f"u tel que F(u) = {fmt(p0)} [tableau N1]", u, "u")
    SN = lT + u * math.sqrt(lT)
    sol.etape("SN = λT + u·√(λT)", SN, "SN")
    sol.etape("Stock (arrondi supérieur)", math.ceil(SN - 1e-9), "stock approx")
    sol.section("Calcul exact (loi de Poisson)")
    cumul, k, terme = 0.0, 0, math.exp(-lT)
    while True:
        cumul += terme
        sol.etape(f"k = {k} : exp(-λT)·Σ(λT)^j/j!", cumul)
        if cumul >= p0 or k > 10000:
            break
        k += 1
        terme *= lT / k
    sol.etape("k minimal (Poisson)", k, "stock exact")
    return sol


# ------------------------------------------------------- analyse ABC

def analyse_abc(donnees, seuil_A=80.0, seuil_B=95.0):
    """Analyse ABC (Pareto).

    donnees : dict {machine: (Ci, Npi)} avec Ci = coût (ou heures d'arrêt) et Npi = nb de pannes.
    """
    sol = Solution("ANALYSE ABC")
    lignes = sorted(donnees.items(), key=lambda kv: -kv[1][0])
    CT = sum(v[0] for _, v in lignes)
    NPT = sum(v[1] for _, v in lignes)
    sol.etape("Coût total CT", CT, "CT")
    sol.etape("Nombre total de pannes NPT", NPT, "NPT")
    sol.section("N° | Ci | ΣCi | ΣCi/CT % | Npi | ΣNpi | ΣNpi/NPT % | classe")
    sC = sN = 0
    table = []
    for nom, (Ci, Ni) in lignes:
        sC += Ci
        sN += Ni
        pc, pn = 100 * sC / CT, 100 * sN / NPT
        classe = "A" if pc <= seuil_A + 1e-9 else ("B" if pc <= seuil_B + 1e-9 else "C")
        table.append((nom, Ci, sC, pc, Ni, sN, pn, classe))
        sol.etape(f"{nom} | {fmt(Ci)} | {fmt(sC)} | {pc:.1f} | {fmt(Ni)} | {fmt(sN)} | {pn:.1f} | {classe}")
    sol.resultats["table"] = table
    for cl in "ABC":
        membres = [str(r[0]) for r in table if r[7] == cl]
        sol.note(f"Classe {cl} : " + (", ".join(membres) if membres else "—"))
    return sol
