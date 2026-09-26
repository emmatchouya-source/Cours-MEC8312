"""M2 — Fiabilité structurale : série, parallèle, r/m, stand-by, systèmes complexes."""
import itertools
import math
import re

from .lois import Exponentielle
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


def serie_exp(lams, t=None, unite_temps="h"):
    """Série, lois exponentielles : λss = Σλi, Rss(t) = exp(-λss·t), MTTFss = 1/λss."""
    sol = Solution("SYSTÈME EN SÉRIE — LOI EXPONENTIELLE", unite_temps)
    lss = sum(lams)
    sol.etape("λss = Σλi = " + " + ".join(fmt(l) for l in lams), lss, "λss", "1/T")
    sol.etape("MTTFss = 1/λss", 1 / lss, "MTTFss", "T")
    for ti in ([t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])):
        sol.etape(f"Rss({fmt(ti)}) = exp(-λss·t)", math.exp(-lss * ti), f"Rss({fmt(ti)})")
    return sol


def parallele_exp(lam, k, t=None, unite_temps="h"):
    """Parallèle de k éléments identiques exponentiels :
    Rps(t) = 1 - [1 - exp(-λt)]^k ; MTTFps = (1/λ)·Σ 1/i."""
    sol = Solution(f"SYSTÈME EN PARALLÈLE — {k} ÉLÉMENTS IDENTIQUES EXPONENTIELS", unite_temps)
    H = sum(1 / i for i in range(1, k + 1))
    sol.etape(f"Σ(1/i), i=1..{k}", H)
    sol.etape(f"MTTFps = (1/λ)·Σ1/i = {fmt(1 / lam)}·{fmt(H)}", H / lam, "MTTFps", "T")
    for ti in ([t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])):
        Ri = math.exp(-lam * ti)
        sol.etape(f"R élément({fmt(ti)}) = exp(-λt)", Ri)
        sol.etape(f"Rps({fmt(ti)}) = 1 - [1 - exp(-λt)]^{k}", 1 - (1 - Ri) ** k, f"Rps({fmt(ti)})")
        if lam * ti < 0.2 and k >= 2:
            sol.etape(f"λps(t) ≈ k·t^(k-1)·Πλi (λt < 0,2)", k * ti ** (k - 1) * lam ** k, unite="1/T")
    return sol


def r_sur_m_exp(lam, r, m, t=None, unite_temps="h"):
    """r/m avec éléments exponentiels identiques : MTTF = (1/λ)·Σ_{y=r}^{m} 1/y."""
    sol = Solution(f"REDONDANCE {r}/{m} — LOI EXPONENTIELLE", unite_temps)
    H = sum(1 / y for y in range(r, m + 1))
    sol.etape(f"Σ(1/y), y={r}..{m}", H)
    sol.etape("MTTF r/m = (1/λ)·Σ1/y", H / lam, "MTTF", "T")
    for ti in ([t] if t is not None and not isinstance(t, (list, tuple)) else (t or [])):
        R = math.exp(-lam * ti)
        sol.etape(f"R élément({fmt(ti)})", R)
        val = sum(math.comb(m, y) * R ** y * (1 - R) ** (m - y) for y in range(r, m + 1))
        sol.etape(f"R_{r}/{m}({fmt(ti)}) = ΣC(m,y)R^y(1-R)^(m-y)", val, f"R({fmt(ti)})")
    return sol


def standby(t=None, lamA=None, lamB=None, Rs=1.0, n=1, lamBsb=0.0, unite_temps="h"):
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
            sol = Solution(f"STAND-BY — MODÈLE {modele} (identiques, 1 actif + {n} en réserve, commutateur idéal)", unite_temps)
            sol.etape(f"MTTFsb = (n+1)/λ = {n + 1}/{fmt(lam)}", (n + 1) / lam, "MTTFsb", "T")
            for ti in ts:
                lt = lam * ti
                sol.etape(f"λt", lt)
                val = sum(lt ** x * math.exp(-lt) / math.factorial(x) for x in range(n + 1))
                sol.etape(f"Rsb({fmt(ti)}) = Σ_(x=0..{n}) (λt)^x·e^(-λt)/x!", val, f"Rsb({fmt(ti)})")
        else:
            if n != 1:
                raise ValueError("Modèle 3 (commutateur non idéal) : n = 1 seulement")
            sol = Solution("STAND-BY — MODÈLE 3 (identiques, commutateur non idéal)", unite_temps)
            sol.etape(f"MTTFsb = (1 + Rs)/λ = (1 + {fmt(Rs)})/{fmt(lam)}", (1 + Rs) / lam, "MTTFsb", "T")
            for ti in ts:
                lt = lam * ti
                sol.etape(f"Rsb({fmt(ti)}) = e^(-λt)·(1 + λt·Rs)", math.exp(-lt) * (1 + lt * Rs), f"Rsb({fmt(ti)})")
        return sol

    if lamBsb:
        sol = Solution("STAND-BY — MODÈLE 6 (différents, réserve pouvant défaillir en attente)", unite_temps)
        a, b, c = lamA, lamB, lamBsb
        sol.etape("MTTFsb = 1/λA + λA/[λB(λA + λBsb)]", 1 / a + a / (b * (a + c)), "MTTFsb", "T")
        for ti in ts:
            val = math.exp(-a * ti) + a / (a + c - b) * (math.exp(-b * ti) - math.exp(-(a + c) * ti))
            sol.etape(f"Rsb({fmt(ti)}) = e^(-λA t) + λA/(λA+λBsb-λB)·[e^(-λB t) - e^(-(λA+λBsb)t)]", val, f"Rsb({fmt(ti)})")
        return sol

    modele = 4 if Rs == 1.0 else 5
    sol = Solution(f"STAND-BY — MODÈLE {modele} (éléments différents" + (", commutateur non idéal)" if modele == 5 else ")"), unite_temps)
    a, b = lamA, lamB
    sol.etape(f"MTTFsb = 1/λA + Rs/λB = 1/{fmt(a)} + {fmt(Rs)}/{fmt(b)}", 1 / a + Rs / b, "MTTFsb", "T")
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


# ------------------------------------------------------- diagramme de fiabilité quelconque

def _relie(s, t, aretes, marche):
    """Existe-t-il un chemin de s à t par des arêtes (a, b, composants, …) dont tous les composants marchent ?"""
    if s == t:
        return True
    adj = {}
    for e in aretes:
        if all(marche(c) for c in e[2]):
            adj.setdefault(e[0], []).append(e[1])
            adj.setdefault(e[1], []).append(e[0])
    vu, pile = {s}, [s]
    while pile:
        x = pile.pop()
        if x == t:
            return True
        for y in adj.get(x, []):
            if y not in vu:
                vu.add(y)
                pile.append(y)
    return False


def _compose(nom):
    return nom.startswith("R") and nom[1:].isdigit()


def _fiab_reseau(s, t, aretes, R, trace=None):
    """R_S exacte d'un réseau à deux bornes : réductions série/parallèle entre arêtes sans
    composant commun, puis conditionnement (composant répété d'abord, sinon l'arête la plus
    centrale). aretes : (a, b, composants, nom) ; trace reçoit les lignes du calcul."""
    compteur = [0]

    def ecrire(prof, texte):
        if trace is not None:
            trace.append("  " * prof + texte)

    def lab(nom):
        return nom if _compose(nom) else "R_" + nom

    def fus(x, y, genre, a, b, prof):
        compteur[0] += 1
        nom = f"R{compteur[0]}"
        if genre == "serie":
            v = x[4] * y[4]
            ecrire(prof, f"{nom} = {lab(x[3])}·{lab(y[3])} = {fmt(x[4])}·{fmt(y[4])} = {fmt(v)}  (série)")
        else:
            v = 1 - (1 - x[4]) * (1 - y[4])
            ecrire(prof, f"{nom} = 1 - (1 - {lab(x[3])})(1 - {lab(y[3])}) = 1 - (1 - {fmt(x[4])})(1 - {fmt(y[4])}) = {fmt(v)}  (parallèle)")
        return [a, b, x[2] + y[2], nom, v]

    def reduire(s, t, A, prof):
        """Une réduction série ou parallèle (ou un élagage) ; renvoie True si le graphe a changé."""
        deg, occ = {}, {}
        for a, b, comps, _, _ in A:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
            for c in comps:
                occ[c] = occ.get(c, 0) + 1
        for i, e in enumerate(A):
            if any(n not in (s, t) and deg[n] == 1 for n in e[:2]):
                del A[i]
                return True

        def seul(e):
            return all(occ[c] == 1 for c in e[2])
        for i in range(len(A)):
            for j in range(i + 1, len(A)):
                x, y = A[i], A[j]
                if seul(x) and seul(y) and {x[0], x[1]} == {y[0], y[1]}:
                    A[i] = fus(x, y, "par", x[0], x[1], prof)
                    del A[j]
                    return True
        for n in deg:
            if n in (s, t) or deg[n] != 2:
                continue
            i, j = [k for k, e in enumerate(A) if n in e[:2]]
            x, y = A[i], A[j]
            if seul(x) and seul(y):
                A[i] = fus(x, y, "serie", x[1] if x[0] == n else x[0], y[1] if y[0] == n else y[0], prof)
                del A[j]
                return True
        return False

    def go(s, t, A, prof):
        A = [list(e) for e in A]
        while True:
            if s == t:
                ecrire(prof, "-> entrée et sortie reliées : R = 1")
                return 1.0
            A = [e for e in A if e[0] != e[1]]
            if not _relie(s, t, A, lambda c: True):
                ecrire(prof, "-> plus aucun chemin de l'entrée à la sortie : R = 0")
                return 0.0
            if len(A) == 1:
                ecrire(prof, f"-> R = {lab(A[0][3])} = {fmt(A[0][4])}")
                return A[0][4]
            if not reduire(s, t, A, prof):
                break
        occ, deg = {}, {}
        for a, b, comps, _, _ in A:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
            for c in comps:
                occ[c] = occ.get(c, 0) + 1
        reps = sorted((c for c in occ if occ[c] > 1), key=lambda c: -occ[c])
        if reps:
            nom = reps[0]
            les, pv = [e for e in A if nom in e[2]], R[nom]
            ecrire(prof, f"Conditionnement sur {nom} (composant répété), R_{nom} = {fmt(pv)}")
        else:
            e = max(A, key=lambda e: deg[e[0]] + deg[e[1]] - (0.5 if _compose(e[3]) else 0))
            les, pv, nom = [e], e[4], e[3]
            ecrire(prof, f"Conditionnement sur {nom} (élément central), {lab(nom)} = {fmt(pv)}")
        autres = [e for e in A if e not in les]
        uf = {}

        def f(x):
            while x in uf:
                x = uf[x]
            return x
        for e in les:
            x, y = f(e[0]), f(e[1])
            if x != y:
                uf[y] = x
        ecrire(prof, f"Si {nom} marche (remplacé par un fil) :")
        vh = go(f(s), f(t), [[f(e[0]), f(e[1])] + e[2:] for e in autres], prof + 1)
        ecrire(prof, f"Si {nom} est en panne (retiré) :")
        vb = go(s, t, autres, prof + 1)
        v = pv * vh + (1 - pv) * vb
        ecrire(prof, f"R = {lab(nom)}·R|{nom} + (1 - {lab(nom)})·R|non {nom} = {fmt(pv)}·{fmt(vh)} + {fmt(1 - pv)}·{fmt(vb)} = {fmt(v)}")
        return v

    return go(s, t, [[a, b, list(c), n, R[n]] for a, b, c, n in aretes], 0)


def schema_fiabilite(blocs, R, t=None, entree="E", sortie="S", unite_temps="h", grandeur="R",
                     cible=None, inconnue=None):
    """Diagramme de fiabilité quelconque (série, parallèle, pont, composants répétés).

    blocs     : liste de (composant, nœud_a, nœud_b) ; un bloc relie deux nœuds, `entree` et
                `sortie` sont les bornes du système. Deux blocs du même nom sont le même
                composant (il marche ou tombe en panne partout à la fois).
                Ex. pont : [("A","E","1"), ("B","E","2"), ("C","1","S"), ("D","2","S"), ("X","1","2")]
    R         : dict {composant: valeur} ; une valeur peut aussi être une loi (Exponentielle,
                Weibull…) : sa fiabilité est prise au temps t. Un bloc r/m ou stand-by se
                donne par sa fiabilité, y compris un stand-by d'éléments différents (modèles
                4 à 6), à condition que ses éléments n'apparaissent pas ailleurs :
                X = standby(t=500, lamA=0.001, lamB=0.002, Rs=0.95)["Rsb(500)"]
                schema_fiabilite([("C", "E", "1"), ("X", "1", "S")], {"C": 0.99, "X": X})
    grandeur  : "R" (fiabilité) ou "A" (disponibilité : donner A_i = μ/(λ + μ), un réparateur
                par composant) ; la même structure s'applique.
    cible, inconnue : calcul inverse. `inconnue` = nom d'un composant -> sa valeur pour que
                R_S = cible (R_S est affine en chaque composant) ; `inconnue` = "t" -> temps t
                tel que R_S(t) = cible (toutes les valeurs de R doivent alors être des lois).
    Calcule R_S pas à pas (réductions série/parallèle puis conditionnement), vérifie par
    énumération des états, donne les chemins et coupes minimaux et l'importance de Birnbaum
    I_B = R_S(R_i = 1) - R_S(R_i = 0).
    """
    L = "A" if grandeur == "A" else "R"
    sol = Solution("DIAGRAMME DE " + ("DISPONIBILITÉ" if L == "A" else "FIABILITÉ"), unite_temps)
    lettre = (lambda s: re.sub(r"\bR(?=[_|0-9 ])", "A", s)) if L == "A" else (lambda s: s)
    Rv = {}
    for n, x in R.items():
        if hasattr(x, "R"):
            if t is None and inconnue != "t":
                raise ValueError(f"{n} est donné par une loi : préciser t")
            Rv[n] = x.R(t) if t is not None else None
        else:
            Rv[n] = x
    comps = []
    for c, _, _ in blocs:
        if c not in comps:
            comps.append(c)
    aretes = [(a, b, [c], c) for c, a, b in blocs]
    if t is not None:
        sol.donnee(f"t = {sol.q(t, 'T')}")
    if all(Rv[n] is not None for n in comps):
        sol.donnee(" ; ".join(f"{L}_{n} = {fmt(Rv[n])}" for n in comps))
        trace = []
        Rs = _fiab_reseau(entree, sortie, aretes, Rv, trace)
        sol.section("Calcul pas à pas (réductions série / parallèle, conditionnement)")
        for ligne in trace:
            sol.etape(lettre(ligne))
        sol.etape(f"{L}_S", Rs, f"{L}_S")
        if len(comps) <= 16:
            marche = {}
            for etat in itertools.product([1, 0], repeat=len(comps)):
                e = dict(zip(comps, etat))
                marche[etat] = _relie(entree, sortie, aretes, lambda c: e[c])
            verif = sum(_prod(Rv[n] if x else 1 - Rv[n] for n, x in zip(comps, etat)) for etat, ok in marche.items() if ok)
            sol.etape(f"Vérification par énumération des {2 ** len(comps)} états", verif, f"{L}_S (énumération)")
            chemins, coupes = [], []
            for etat, ok in marche.items():
                if ok and all(not marche[etat[:i] + (0,) + etat[i + 1:]] for i, x in enumerate(etat) if x):
                    chemins.append([n for n, x in zip(comps, etat) if x])
                if not ok and all(marche[etat[:i] + (1,) + etat[i + 1:]] for i, x in enumerate(etat) if not x):
                    coupes.append([n for n, x in zip(comps, etat) if not x])
            chemins.sort(key=lambda c: (len(c), c))
            coupes.sort(key=lambda c: (len(c), c))
            sol.resultats["chemins"], sol.resultats["coupes"] = chemins, coupes
            sep = lambda c: ("–" if any(len(n) > 1 for n in c) else "").join(c)
            sol.note("Chemins de succès minimaux : " + ", ".join(sep(c) for c in chemins))
            sol.note("Coupes minimales : " + ", ".join(sep(c) for c in coupes))
        sol.section(f"Importance de Birnbaum I_B = {L}_S({L}_i = 1) - {L}_S({L}_i = 0)")
        IB = {}
        for n in comps:
            IB[n] = (_fiab_reseau(entree, sortie, aretes, dict(Rv, **{n: 1.0}))
                     - _fiab_reseau(entree, sortie, aretes, dict(Rv, **{n: 0.0})))
            sol.etape(f"I_B({n})", IB[n])
        sol.resultats["I_B"] = IB
    if cible is None or inconnue is None:
        return sol
    sol.section(f"Calcul inverse : {L}_S visée = {fmt(cible)}")
    if inconnue == "t":
        lois = {n: x for n, x in R.items() if hasattr(x, "R")}

        def Rs_t(x):
            return _fiab_reseau(entree, sortie, aretes, {n: (R[n].R(x) if n in lois else R[n]) for n in comps})
        if Rs_t(0.0) < cible:
            sol.note(f"Impossible : même à t = 0, {L}_S = {fmt(Rs_t(0.0))}.")
            return sol
        hi = max((getattr(x, "mttf", lambda: 1.0)() for x in lois.values()), default=1.0)
        for _ in range(80):
            if Rs_t(hi) <= cible:
                break
            hi *= 2
        else:
            sol.note(f"Impossible : {L}_S ne descend pas sous {fmt(Rs_t(hi))}.")
            return sol
        lo = 0.0
        for _ in range(200):
            m = (lo + hi) / 2
            lo, hi = (m, hi) if Rs_t(m) > cible else (lo, m)
        sol.etape(f"t tel que {L}_S(t) = {fmt(cible)} (dichotomie, {L}_S décroît avec t)", (lo + hi) / 2, "t visé", "T")
        return sol
    h = _fiab_reseau(entree, sortie, aretes, dict(Rv, **{inconnue: 1.0}))
    b = _fiab_reseau(entree, sortie, aretes, dict(Rv, **{inconnue: 0.0}))
    sol.etape(f"{L}_S({L}_{inconnue} = 0)", b)
    sol.etape(f"{L}_S({L}_{inconnue} = 1)", h)
    if cible > h + 1e-15:
        sol.note(f"Impossible : même avec {L}_{inconnue} = 1, {L}_S = {fmt(h)} < {fmt(cible)}.")
        return sol
    if cible < b - 1e-15:
        sol.note(f"Déjà atteint : même avec {L}_{inconnue} = 0, {L}_S = {fmt(b)}.")
        return sol
    x = (cible - b) / (h - b)
    sol.etape(f"{L}_S est affine en {L}_{inconnue} : {L}_{inconnue} = ({fmt(cible)} - {fmt(b)})/{fmt(h - b)}", x, f"{L}_{inconnue} visée")
    loi = R.get(inconnue)
    if isinstance(loi, Exponentielle) and t and 0 < x < 1:
        sol.etape(f"λ_{inconnue} = -ln({fmt(x)})/t", -math.log(x) / t, f"λ_{inconnue} visé", "1/T")
    return sol
