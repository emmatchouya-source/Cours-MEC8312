"""M3 — Systèmes réparables : chaînes de Markov, disponibilité, maintenabilité."""
import math

from .numerique import expm, resoudre
from .rapport import Solution, fmt


def _taux(sol, lam, mttf, mu, mttr):
    if lam is None:
        lam = 1 / mttf
        sol.etape(f"λ = 1/MTTF = 1/{fmt(mttf)}", lam, "λ")
    if mu is None:
        mu = 1 / mttr
        sol.etape(f"μ = 1/MTTR = 1/{fmt(mttr)}", mu, "μ")
    return lam, mu


def _ts(t):
    if t is None:
        return []
    return list(t) if isinstance(t, (list, tuple)) else [t]


def element_reparable(lam=None, mu=None, t=None, mttf=None, mttr=None):
    """Élément réparable (2 états, λ et μ constants)."""
    sol = Solution("ÉLÉMENT RÉPARABLE — Markov 2 états (1 marche, 0 panne)")
    lam, mu = _taux(sol, lam, mttf, mu, mttr)
    s = lam + mu
    sol.etape("A∞ = μ/(λ + μ) = MTTF/(MTTF + MTTR)", mu / s, "A∞")
    sol.etape("Indisponibilité ID∞ = λ/(λ + μ)", lam / s, "ID∞")
    for ti in _ts(t):
        sol.section(f"t = {fmt(ti)}")
        A = mu / s + lam / s * math.exp(-s * ti)
        sol.etape("A(t) = P1(t) = μ/(λ+μ) + λ/(λ+μ)·exp[-(λ+μ)t]", A, f"A({fmt(ti)})")
        sol.etape("P0(t) = λ/(λ+μ) - λ/(λ+μ)·exp[-(λ+μ)t]", 1 - A)
        sol.etape("R(t) = exp(-λt)", math.exp(-lam * ti), f"R({fmt(ti)})")
        sol.etape("M(t) = 1 - exp(-μt)  (maintenabilité)", 1 - math.exp(-mu * ti), f"M({fmt(ti)})")
    return sol


def temps_moyens(MUT=None, MDT=None, MTBF=None):
    """A∞ = MUT/(MDT + MUT) = MUT/MTBF ; ID∞ = MDT/MTBF."""
    sol = Solution("TEMPS MOYENS — disponibilité asymptotique")
    if MTBF is None:
        MTBF = MUT + MDT
        sol.etape("MTBF = MUT + MDT", MTBF, "MTBF")
    elif MUT is None:
        MUT = MTBF - MDT
        sol.etape("MUT = MTBF - MDT", MUT, "MUT")
    elif MDT is None:
        MDT = MTBF - MUT
        sol.etape("MDT = MTBF - MUT", MDT, "MDT")
    sol.etape("A∞ = MUT/MTBF", MUT / MTBF, "A∞")
    sol.etape("ID∞ = MDT/MTBF", MDT / MTBF, "ID∞")
    return sol


def serie_reparable(lams, mus, t=None):
    """Série d'éléments réparables (un réparateur par élément)."""
    sol = Solution("SYSTÈME RÉPARABLE SÉRIE")
    A_inf = 1.0
    for i, (l, m) in enumerate(zip(lams, mus), 1):
        a = m / (l + m)
        sol.etape(f"A{i}∞ = μ{i}/(λ{i}+μ{i})", a)
        A_inf *= a
    sol.etape("Ass(∞) = Π μi/(λi+μi)", A_inf, "Ass∞")
    for ti in _ts(t):
        sol.section(f"t = {fmt(ti)}")
        A = 1.0
        for l, m in zip(lams, mus):
            A *= m / (l + m) + l / (l + m) * math.exp(-(l + m) * ti)
        sol.etape("Ass(t) = Π Ai(t)", A, f"Ass({fmt(ti)})")
        sol.etape("Rss(t) = exp(-Σλi·t)", math.exp(-sum(lams) * ti), f"Rss({fmt(ti)})")
    return sol


def parallele_reparable(lams, mus, t=None):
    """Parallèle d'éléments réparables (un réparateur par élément, redondance active) :
    disponibilité uniquement."""
    sol = Solution("SYSTÈME RÉPARABLE PARALLÈLE — disponibilité")
    ind = 1.0
    for i, (l, m) in enumerate(zip(lams, mus), 1):
        sol.etape(f"λ{i}/(λ{i}+μ{i})", l / (l + m))
        ind *= l / (l + m)
    sol.etape("Aps(∞) = 1 - Π λi/(λi+μi)", 1 - ind, "Aps∞")
    for ti in _ts(t):
        prod = 1.0
        for l, m in zip(lams, mus):
            prod *= l / (l + m) - l / (l + m) * math.exp(-(l + m) * ti)
        sol.etape(f"Aps({fmt(ti)}) = 1 - Π[1 - Ai(t)]", 1 - prod, f"Aps({fmt(ti)})")
    return sol


def parallele2_reparable(lam=None, mu=None, t=None, mttf=None, mttr=None):
    """Deux éléments identiques réparables en parallèle (un réparateur par élément)
    — exemple des deux pompes (diapos 31-35 + document 'Fiabilité exemple')."""
    sol = Solution("DEUX ÉLÉMENTS IDENTIQUES RÉPARABLES EN PARALLÈLE")
    lam, mu = _taux(sol, lam, mttf, mu, mttr)
    sol.note("États : 2 (deux marchent), 1 (un en réparation), 0 (deux en panne)")
    sol.note("dP2/dt = -2λP2 + μP1 ; dP1/dt = -(λ+μ)P1 + 2λP2")
    mttf_s = (mu + 3 * lam) / (2 * lam ** 2)
    sol.etape(f"MTTFps = (μ + 3λ)/(2λ²)", mttf_s, "MTTFps", "")
    A_inf = 1 - (lam / (lam + mu)) ** 2
    sol.etape("Aps∞ = 1 - [λ/(λ+μ)]²", A_inf, "Aps∞")
    # racines de s² + (3λ+μ)s + 2λ² = 0
    B, C = 3 * lam + mu, 2 * lam ** 2
    d = math.sqrt(B * B - 4 * C)
    s1, s2 = (-B + d) / 2, (-B - d) / 2
    sol.etape("s1 (racine de s² + (3λ+μ)s + 2λ² = 0)", s1)
    sol.etape("s2", s2)
    for ti in _ts(t):
        sol.section(f"t = {fmt(ti)}")
        R_exact = ((s1 + B) / (s1 - s2)) * math.exp(s1 * ti) + ((s2 + B) / (s2 - s1)) * math.exp(s2 * ti)
        R_app = math.exp(-2 * lam ** 2 / (mu + 3 * lam) * ti)
        sol.etape("Rps(t) ≈ exp[-2λ²t/(μ + 3λ)]  (λ << μ, formule du cours)", R_app, f"Rps({fmt(ti)})")
        sol.etape("Rps(t) exact (Laplace, s1 et s2)", R_exact, f"Rps_exact({fmt(ti)})")
        Ai = lam / (lam + mu) - lam / (lam + mu) * math.exp(-(lam + mu) * ti)
        sol.etape("Aps(t) = 1 - [λ/(λ+μ) - λ/(λ+μ)·exp(-(λ+μ)t)]²", 1 - Ai ** 2, f"Aps({fmt(ti)})")
        sol.section(f"Comparaison : même système NON réparable, t = {fmt(ti)}")
        sol.etape("Rps(t) = 1 - [1 - exp(-λt)]²", 1 - (1 - math.exp(-lam * ti)) ** 2, f"Rps_non_rep({fmt(ti)})")
    sol.etape("MTTFps non réparable = MTTF·(1 + 1/2)", 1.5 / lam, "MTTFps_non_rep")
    return sol


def chaine_markov(transitions, initial, pannes, t=None, etats=None):
    """Chaîne de Markov quelconque (ex. 3 processeurs / 2 sources d'énergie).

    transitions : dict {(de, vers): taux}  ex. {("3,2","2,2"): 3*l1, ...}
    initial     : état initial (probabilité 1 à t=0)
    pannes      : liste des états de panne
    t           : instant(s) de calcul

    Calcule :
      - les équations différentielles dPi/dt ;
      - P(t) pour chaque état et la disponibilité A(t) = Σ P(états de marche) ;
      - la fiabilité R(t) (états de panne rendus absorbants) et la MTTF ;
      - la disponibilité stationnaire (si la chaîne n'a pas d'état absorbant).
    """
    if etats is None:
        etats = []
        for (a, b) in transitions:
            for e in (a, b):
                if e not in etats:
                    etats.append(e)
    n = len(etats)
    idx = {e: i for i, e in enumerate(etats)}
    marche = [e for e in etats if e not in pannes]
    sol = Solution(f"CHAÎNE DE MARKOV — {n} états, départ en {initial}")

    def generateur(absorbant):
        Q = [[0.0] * n for _ in range(n)]
        for (a, b), r in transitions.items():
            if absorbant and a in pannes:
                continue
            Q[idx[a]][idx[b]] += r
            Q[idx[a]][idx[a]] -= r
        return Q

    Q = generateur(False)
    sol.section("Équations différentielles")
    for e in etats:
        j = idx[e]
        termes = []
        if Q[j][j]:
            termes.append(f"{fmt(Q[j][j])}·P[{e}]")
        for (a, b), r in transitions.items():
            if b == e:
                termes.append(f"+{fmt(r)}·P[{a}]")
        sol.etape(f"dP[{e}]/dt = " + " ".join(termes) if termes else f"dP[{e}]/dt = 0")

    p0 = [0.0] * n
    p0[idx[initial]] = 1.0
    Qa = generateur(True)
    for ti in _ts(t):
        sol.section(f"t = {fmt(ti)}")
        E = expm([[x * ti for x in ligne] for ligne in Q])
        P = [sum(p0[i] * E[i][j] for i in range(n)) for j in range(n)]
        for e in etats:
            sol.etape(f"P[{e}]({fmt(ti)})", P[idx[e]])
        sol.etape("A(t) = Σ P(états de marche)", sum(P[idx[e]] for e in marche), f"A({fmt(ti)})")
        Ea = expm([[x * ti for x in ligne] for ligne in Qa])
        Pa = [sum(p0[i] * Ea[i][j] for i in range(n)) for j in range(n)]
        sol.etape("R(t) (pannes absorbantes)", sum(Pa[idx[e]] for e in marche), f"R({fmt(ti)})")
        for e in pannes:
            sol.etape(f"P(panne {e}) absorbante", Pa[idx[e]])
        sol.etape("P(F) = Σ P(pannes) = 1 - R(t)", sum(Pa[idx[e]] for e in pannes), f"P(F)({fmt(ti)})")

    # MTTF : temps moyen passé dans les états de marche avant absorption
    if initial in marche:
        m = len(marche)
        A = [[-Qa[idx[marche[j]]][idx[marche[i]]] for j in range(m)] for i in range(m)]
        b = [1.0 if e == initial else 0.0 for e in marche]
        tau = resoudre(A, b)
        sol.section("Temps moyen avant la première défaillance")
        sol.etape("MTTF = Σ temps moyens passés dans les états de marche", sum(tau), "MTTF")

    # stationnaire : π Q = 0, Σ π = 1
    absorbants = [e for e in etats if all(Q[idx[e]][j] == 0 for j in range(n))]
    if not absorbants:
        A = [[Q[j][i] for j in range(n)] for i in range(n)]
        A[-1] = [1.0] * n
        b = [0.0] * (n - 1) + [1.0]
        pi = resoudre(A, b)
        sol.section("Régime stationnaire")
        for e in etats:
            sol.etape(f"π[{e}]", pi[idx[e]])
        sol.etape("A∞ = Σ π(états de marche)", sum(pi[idx[e]] for e in marche), "A∞")
    return sol
