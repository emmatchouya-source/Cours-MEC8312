"""M2 — Fiabilité structurale : série, parallèle, r/m, stand-by, systèmes complexes."""
import itertools
import math

from .rapport import Solution, fmt


def _prod(xs):
    p = 1.0
    for x in xs:
        p *= x
    return p


def serie(R):
    """Système série : Rss = R1·R2·…·Rk (liste de fiabilités)."""
    sol = Solution("SYSTÈME EN SÉRIE")
    sol.etape("Rss = " + " · ".join(fmt(r) for r in R), _prod(R), "Rss")
    return sol


def parallele(R):
    """Système parallèle : Rps = 1 - (1-R1)(1-R2)…(1-Rk)."""
    sol = Solution("SYSTÈME EN PARALLÈLE")
    F = _prod(1 - r for r in R)
    sol.etape("Π(1 - Ri) = " + " · ".join(fmt(1 - r) for r in R), F)
    sol.etape("Rps = 1 - Π(1 - Ri)", 1 - F, "Rps")
    return sol


def r_sur_m(R, r, m):
    """Redondance active r/m (éléments identiques) :
    R = Σ_{y=r}^{m} C(m,y)·R^y·(1-R)^(m-y)."""
    sol = Solution(f"REDONDANCE ACTIVE {r}/{m} (R élément = {fmt(R)})")
    total = 0.0
    for y in range(r, m + 1):
        terme = math.comb(m, y) * R ** y * (1 - R) ** (m - y)
        sol.etape(f"y={y} : C({m},{y})·R^{y}·(1-R)^{m - y} = {math.comb(m, y)}·{fmt(R ** y)}·{fmt((1 - R) ** (m - y))}", terme)
        total += terme
    sol.etape(f"R_{r}/{m} = Σ", total, f"R_{r}/{m}")
    return sol


def serie_exp(lams, t=None):
    """Série, lois exponentielles : λss = Σλi, Rss(t) = exp(-λss·t), MTTFss = 1/λss."""
    sol = Solution("SYSTÈME EN SÉRIE — LOI EXPONENTIELLE")
    lss = sum(lams)
    sol.etape("λss = Σλi = " + " + ".join(fmt(l) for l in lams), lss, "λss")
    sol.etape("MTTFss = 1/λss", 1 / lss, "MTTFss")
    for ti in ([t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])):
        sol.etape(f"Rss({fmt(ti)}) = exp(-λss·t)", math.exp(-lss * ti), f"Rss({fmt(ti)})")
    return sol


def parallele_exp(lam, k, t=None):
    """Parallèle de k éléments identiques exponentiels :
    Rps(t) = 1 - [1 - exp(-λt)]^k ; MTTFps = (1/λ)·Σ 1/i."""
    sol = Solution(f"SYSTÈME EN PARALLÈLE — {k} ÉLÉMENTS IDENTIQUES EXPONENTIELS")
    H = sum(1 / i for i in range(1, k + 1))
    sol.etape(f"Σ(1/i), i=1..{k}", H)
    sol.etape(f"MTTFps = (1/λ)·Σ1/i = {fmt(1 / lam)}·{fmt(H)}", H / lam, "MTTFps")
    for ti in ([t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])):
        Ri = math.exp(-lam * ti)
        sol.etape(f"R élément({fmt(ti)}) = exp(-λt)", Ri)
        sol.etape(f"Rps({fmt(ti)}) = 1 - [1 - exp(-λt)]^{k}", 1 - (1 - Ri) ** k, f"Rps({fmt(ti)})")
        if lam * ti < 0.2 and k >= 2:
            sol.etape(f"λps(t) ≈ k·t^(k-1)·Πλi (λt < 0,2)", k * ti ** (k - 1) * lam ** k)
    return sol


def r_sur_m_exp(lam, r, m, t=None):
    """r/m avec éléments exponentiels identiques : MTTF = (1/λ)·Σ_{y=r}^{m} 1/y."""
    sol = Solution(f"REDONDANCE {r}/{m} — LOI EXPONENTIELLE")
    H = sum(1 / y for y in range(r, m + 1))
    sol.etape(f"Σ(1/y), y={r}..{m}", H)
    sol.etape("MTTF r/m = (1/λ)·Σ1/y", H / lam, "MTTF")
    for ti in ([t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])):
        R = math.exp(-lam * ti)
        sol.etape(f"R élément({fmt(ti)})", R)
        val = sum(math.comb(m, y) * R ** y * (1 - R) ** (m - y) for y in range(r, m + 1))
        sol.etape(f"R_{r}/{m}({fmt(ti)}) = ΣC(m,y)R^y(1-R)^(m-y)", val, f"R({fmt(ti)})")
    return sol


def standby(t=None, lamA=None, lamB=None, Rs=1.0, n=1, lamBsb=0.0):
    """Système avec réserve (stand-by), modèles 1 à 6 du cours.

    - Éléments identiques : donner lamA seulement (n = nombre d'éléments en réserve).
        Modèle 1 : n=1, Rs=1 ; Modèle 2 : n quelconque, Rs=1 ; Modèle 3 : n=1, Rs<1.
    - Éléments différents : lamA (actif), lamB (réserve une fois en marche).
        Modèle 4 : Rs=1 ; Modèle 5 : Rs<1 ; Modèle 6 : lamBsb>0 (défaillance en attente), Rs=1.
    """
    ts = [t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])
    identiques = lamB is None or (lamB == lamA and not lamBsb)
    if identiques:
        lam = lamA
        if Rs == 1.0:
            modele = 1 if n == 1 else 2
            sol = Solution(f"STAND-BY — MODÈLE {modele} (identiques, 1 actif + {n} en réserve, commutateur idéal)")
            sol.etape(f"MTTFsb = (n+1)/λ = {n + 1}/{fmt(lam)}", (n + 1) / lam, "MTTFsb")
            for ti in ts:
                lt = lam * ti
                sol.etape(f"λt", lt)
                val = sum(lt ** x * math.exp(-lt) / math.factorial(x) for x in range(n + 1))
                sol.etape(f"Rsb({fmt(ti)}) = Σ_(x=0..{n}) (λt)^x·e^(-λt)/x!", val, f"Rsb({fmt(ti)})")
        else:
            if n != 1:
                raise ValueError("Modèle 3 (commutateur non idéal) : n = 1 seulement")
            sol = Solution("STAND-BY — MODÈLE 3 (identiques, commutateur non idéal)")
            sol.etape(f"MTTFsb = (1 + Rs)/λ = (1 + {fmt(Rs)})/{fmt(lam)}", (1 + Rs) / lam, "MTTFsb")
            for ti in ts:
                lt = lam * ti
                sol.etape(f"Rsb({fmt(ti)}) = e^(-λt)·(1 + λt·Rs)", math.exp(-lt) * (1 + lt * Rs), f"Rsb({fmt(ti)})")
        return sol

    if lamBsb:
        sol = Solution("STAND-BY — MODÈLE 6 (différents, réserve pouvant défaillir en attente)")
        a, b, c = lamA, lamB, lamBsb
        sol.etape("MTTFsb = 1/λA + λA/[λB(λA + λBsb)]", 1 / a + a / (b * (a + c)), "MTTFsb")
        for ti in ts:
            val = math.exp(-a * ti) + a / (a + c - b) * (math.exp(-b * ti) - math.exp(-(a + c) * ti))
            sol.etape(f"Rsb({fmt(ti)}) = e^(-λA t) + λA/(λA+λBsb-λB)·[e^(-λB t) - e^(-(λA+λBsb)t)]", val, f"Rsb({fmt(ti)})")
        return sol

    modele = 4 if Rs == 1.0 else 5
    sol = Solution(f"STAND-BY — MODÈLE {modele} (éléments différents" + (", commutateur non idéal)" if modele == 5 else ")"))
    a, b = lamA, lamB
    sol.etape(f"MTTFsb = 1/λA + Rs/λB = 1/{fmt(a)} + {fmt(Rs)}/{fmt(b)}", 1 / a + Rs / b, "MTTFsb")
    for ti in ts:
        val = math.exp(-a * ti) + a * Rs / (b - a) * (math.exp(-a * ti) - math.exp(-b * ti))
        sol.etape(f"Rsb({fmt(ti)}) = e^(-λA t) + λA·Rs/(λB-λA)·[e^(-λA t) - e^(-λB t)]", val, f"Rsb({fmt(ti)})")
    return sol


def pont(RA, RB, RC, RD, RE):
    """Système en pont (A-C en haut, B-D en bas, E au milieu) par la méthode de
    la probabilité conditionnelle sur E (diapos 78-80)."""
    sol = Solution("SYSTÈME COMPLEXE EN PONT — probabilité conditionnelle sur E")
    R1 = (1 - (1 - RA) * (1 - RB)) * (1 - (1 - RC) * (1 - RD))
    R2 = 1 - (1 - RA * RC) * (1 - RB * RD)
    sol.section("E en marche : (A // B) en série avec (C // D)")
    sol.etape("R1 = [1-(1-RA)(1-RB)]·[1-(1-RC)(1-RD)]", R1, "R1")
    sol.section("E en panne : (A-C) en parallèle avec (B-D)")
    sol.etape("R2 = 1 - (1 - RA·RC)(1 - RB·RD)", R2, "R2")
    sol.etape("Rs = R1·RE + R2·(1 - RE)", R1 * RE + R2 * (1 - RE), "Rs")
    return sol


def systeme_chemins(chemins, R, afficher_etats=False):
    """Système complexe quelconque par énumération des états (2^n états).

    chemins : liste des chemins minimaux (listes de noms de composants) ; le système
              fonctionne si tous les composants d'au moins un chemin fonctionnent.
              Ex. pont : [["A","C"], ["B","D"], ["A","E","D"], ["B","E","C"]]
    R       : dict {nom: fiabilité}
    """
    noms = sorted(R)
    sol = Solution(f"SYSTÈME COMPLEXE — énumération des {2 ** len(noms)} états")
    sol.donnee("chemins minimaux : " + ", ".join("-".join(c) for c in chemins))
    total, n_ok = 0.0, 0
    for etat in itertools.product([1, 0], repeat=len(noms)):
        e = dict(zip(noms, etat))
        if any(all(e[c] for c in ch) for ch in chemins):
            p = _prod(R[n] if e[n] else 1 - R[n] for n in noms)
            total += p
            n_ok += 1
            if afficher_etats:
                sol.etape("".join(n if e[n] else n.lower() + "̄" for n in noms), p)
    sol.etape(f"Nombre d'états de marche", n_ok)
    sol.etape("Rs = Σ P(états de marche)", total, "Rs")
    return sol


def systeme_coupes(coupes, R):
    """Système défini par ses coupes minimales : il tombe en panne si tous les
    composants d'au moins une coupe sont en panne."""
    noms = sorted(R)
    sol = Solution(f"SYSTÈME — coupes minimales, énumération des {2 ** len(noms)} états")
    total = 0.0
    for etat in itertools.product([1, 0], repeat=len(noms)):
        e = dict(zip(noms, etat))
        if not any(all(not e[c] for c in cp) for cp in coupes):
            total += _prod(R[n] if e[n] else 1 - R[n] for n in noms)
    sol.etape("Rs", total, "Rs")
    return sol
