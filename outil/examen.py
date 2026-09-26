#!/usr/bin/env python3
"""Menu interactif pour résoudre les exercices MEC8312 à l'examen.

    python3 outil/examen.py            (menu)
    python3 outil/examen.py 7          (ouvre directement le solveur n° 7)

Saisie des nombres : 0,994  0.994  5e-4  1/500  2*3  sqrt(2)  (virgule ou point).
Listes : séparer par des points-virgules, ex. 100; 200; 500.
Laisser vide = donnée non fournie (valeur par défaut).
"""
import ast
import math
import operator
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fiabilite as F  # noqa: E402

# ------------------------------------------------------------ saisie sûre

_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.Pow: operator.pow, ast.USub: operator.neg,
        ast.UAdd: operator.pos}
_FUNCS = {k: getattr(math, k) for k in ("sqrt", "exp", "log", "log10", "gamma")}
_CONSTS = {"pi": math.pi, "e": math.e}


def _eval(noeud):
    if isinstance(noeud, ast.Expression):
        return _eval(noeud.body)
    if isinstance(noeud, ast.Constant) and isinstance(noeud.value, (int, float)):
        return noeud.value
    if isinstance(noeud, ast.BinOp) and type(noeud.op) in _OPS:
        return _OPS[type(noeud.op)](_eval(noeud.left), _eval(noeud.right))
    if isinstance(noeud, ast.UnaryOp) and type(noeud.op) in _OPS:
        return _OPS[type(noeud.op)](_eval(noeud.operand))
    if isinstance(noeud, ast.Call) and isinstance(noeud.func, ast.Name) and noeud.func.id in _FUNCS:
        return _FUNCS[noeud.func.id](*[_eval(a) for a in noeud.args])
    if isinstance(noeud, ast.Name) and noeud.id in _CONSTS:
        return _CONSTS[noeud.id]
    raise ValueError("expression non reconnue")


def nombre(texte):
    texte = texte.strip().replace(",", ".").replace("^", "**").replace("×", "*")
    return _eval(ast.parse(texte, mode="eval"))


def demander(invite, genre="nombre", defaut=None):
    """genre : nombre | entier | liste | texte ; defaut renvoyé si vide."""
    suffixe = f" [{defaut}]" if defaut is not None else " [vide = non donné]"
    while True:
        try:
            brut = input(f"  {invite}{suffixe} : ").strip()
        except EOFError:
            print()
            sys.exit(0)
        if not brut:
            return defaut
        try:
            if genre == "texte":
                return brut
            if genre == "entier":
                return int(round(nombre(brut)))
            if genre == "liste":
                valeurs = [nombre(x) for x in brut.split(";") if x.strip()]
                return valeurs if len(valeurs) > 1 else valeurs[0]
            return nombre(brut)
        except Exception as err:  # noqa: BLE001
            print(f"    ! saisie invalide ({err}), recommencez")


# ------------------------------------------------------------ solveurs du menu

def s_exponentielle():
    lam = demander("λ (vide -> saisir MTTF)")
    mttf = None if lam is not None else demander("MTTF")
    return F.exponentielle(lam=lam, mttf=mttf, t=demander("t (un ou plusieurs ; séparés par ;)", "liste"),
                           t1=demander("t1 (intervalle Δt)"), t2=demander("t2"),
                           R_cible=demander("R cible (trouver t)"))


def s_normale():
    return F.normale(demander("m (moyenne)"), demander("σ (écart-type)"),
                     t=demander("t (liste avec ;)", "liste"), t1=demander("t1"), t2=demander("t2"),
                     R_cible=demander("R cible"))


def s_weibull():
    return F.weibull(b=demander("b ou β (forme)"), theta=demander("θ ou η (échelle)"),
                     t0=demander("t0 (position)", defaut=0.0), t=demander("t (liste avec ;)", "liste"),
                     t1=demander("t1"), t2=demander("t2"), R_cible=demander("R cible"))


def s_lognormale():
    m = demander("m (moyenne de T)")
    if m is None:
        return F.lognormale(mu_y=demander("μy"), sigma_y=demander("σy"),
                            t=demander("t (liste)", "liste"), R_cible=demander("R cible"))
    return F.lognormale(m=m, sigma=demander("σ (écart-type de T)"), t=demander("t (liste)", "liste"),
                        t1=demander("t1"), t2=demander("t2"), R_cible=demander("R cible"))


def s_gamma():
    return F.gamma(demander("x"))


def s_serie_parallele():
    R = demander("fiabilités des éléments (séparées par ;)", "liste")
    R = R if isinstance(R, list) else [R]
    print(F.serie(R))
    return F.parallele(R)


def s_rm():
    R = demander("R d'un élément (vide si loi exponentielle)")
    r, m = demander("r (nombre minimal qui doit fonctionner)", "entier"), demander("m (nombre total)", "entier")
    if R is not None:
        return F.r_sur_m(R, r, m)
    return F.r_sur_m_exp(demander("λ"), r, m, t=demander("t", "liste"))


def s_exp_systemes():
    choix = demander("1 = série (λ différents), 2 = parallèle k identiques", "entier", 1)
    if choix == 1:
        l = demander("λ des éléments (séparés par ;)", "liste")
        return F.serie_exp(l if isinstance(l, list) else [l], t=demander("t", "liste"))
    return F.parallele_exp(demander("λ"), demander("k", "entier"), t=demander("t", "liste"))


def s_standby():
    return F.standby(t=demander("t (liste)", "liste"), lamA=demander("λA (actif ; ou λ si identiques)"),
                     lamB=demander("λB (réserve en marche ; vide si identiques)"),
                     Rs=demander("Rs commutateur", defaut=1.0), n=demander("n en réserve", "entier", 1),
                     lamBsb=demander("λBsb (défaillance en attente, modèle 6)", defaut=0.0))


def s_complexe():
    print("  Composants : lettres (ex. A B C D E). Chemins minimaux séparés par ; ex. AC; BD; AED; BEC")
    noms = demander("noms des composants (séparés par des espaces)", "texte", "A B C D E").split()
    chemins = [list(c.strip().replace("-", "")) for c in demander("chemins", "texte", "AC; BD; AED; BEC").split(";")]
    R = {}
    commun = demander("R commun à tous (vide pour saisir chacun)")
    for n in noms:
        R[n] = commun if commun is not None else demander(f"R{n}")
    return F.systeme_chemins(chemins, R, afficher_etats=demander("afficher les états (1/0)", "entier", 0) == 1)


def s_element_reparable():
    lam = demander("λ (vide -> saisir MTTF)")
    mttf = None if lam is not None else demander("MTTF (ou MTBF)")
    mu = demander("μ (vide -> saisir MTTR)")
    mttr = None if mu is not None else demander("MTTR")
    return F.element_reparable(lam, mu, t=demander("t (liste)", "liste"), mttf=mttf, mttr=mttr)


def s_temps_moyens():
    return F.temps_moyens(MUT=demander("MUT"), MDT=demander("MDT"), MTBF=demander("MTBF"))


def s_serie_par_reparable():
    choix = demander("1 = série, 2 = parallèle (disponibilité)", "entier", 1)
    l = demander("λ des éléments (;)", "liste")
    m = demander("μ des éléments (;)", "liste")
    l, m = (l if isinstance(l, list) else [l]), (m if isinstance(m, list) else [m])
    if len(m) == 1 and len(l) > 1:
        m = m * len(l)
    f = F.serie_reparable if choix == 1 else F.parallele_reparable
    return f(l, m, t=demander("t (liste)", "liste"))


def s_deux_paralleles():
    lam = demander("λ (vide -> saisir MTTF)")
    mttf = None if lam is not None else demander("MTTF d'un élément")
    mu = demander("μ (vide -> saisir MTTR)")
    mttr = None if mu is not None else demander("MTTR")
    return F.parallele2_reparable(lam, mu, t=demander("t mission (liste)", "liste"), mttf=mttf, mttr=mttr)


def s_markov():
    print("  Saisir les transitions une par ligne : DE VERS TAUX   (ex. 3,2 2,2 3*0.001)")
    print("  Ligne vide pour terminer.")
    tr = {}
    while True:
        ligne = input("    > ").strip()
        if not ligne:
            break
        try:
            a, b, taux = ligne.split(maxsplit=2)
            tr[(a, b)] = nombre(taux)
        except Exception as err:  # noqa: BLE001
            print(f"    ! ligne invalide ({err})")
    initial = demander("état initial", "texte")
    pannes = demander("états de panne (séparés par des espaces)", "texte").split()
    return F.chaine_markov(tr, initial, pannes, t=demander("t (liste)", "liste"))


def s_preventive():
    print("  λ(t) : 1 = a + b·t ; 2 = Weibull (β, η) ; 3 = constant ; 4 = expression en t")
    c = demander("type", "entier", 1)
    if c == 1:
        taux = F.taux_lineaire(demander("a"), demander("b"))
    elif c == 2:
        taux = F.taux_weibull(demander("β"), demander("η"))
    elif c == 3:
        taux = F.taux_constant(demander("λ"))
    else:
        taux = F.taux_expression(demander("λ(t) =", "texte").replace("^", "**"))
    n = demander("n maintenances (liste ;)", "liste")
    if n is not None:
        n = [int(round(x)) for x in (n if isinstance(n, list) else [n])]
    return F.preventive_systematique(taux, demander("T (durée)"), n=n, R_cible=demander("R cible (n minimal)"))


def s_remplacement():
    return F.periode_optimale_remplacement(demander("CI (investissement)"), demander("ice (croissance coût exploitation)"),
                                           demander("icm (croissance coût maintenance)"),
                                           k1=demander("k1"), k2=demander("k2"))


def s_remplacement_age():
    loi = F.Weibull(demander("β"), demander("η"))
    return F.remplacement_age(loi, demander("Cp (préventif)"), demander("Cd (défaillance)"))


def s_weibull_systematique():
    beta, eta = demander("β"), demander("η")
    r = demander("r = P/p (liste ;)", "liste")
    if r is None:
        return F.weibull_systematique(beta, eta, p=demander("p (coût direct)"), P=demander("P (coût indirect)"))
    for ri in (r if isinstance(r, list) else [r])[:-1]:
        print(F.weibull_systematique(beta, eta, ri))
    return F.weibull_systematique(beta, eta, (r if isinstance(r, list) else [r])[-1])


def s_corrective():
    return F.corrective(demander("A'∞ (disponibilité actuelle)"), demander("C'mv (coût main-d'œuvre actuel)"),
                        demander("MTTF (ou MTBF)"), demander("c (perte par unité de temps)"),
                        Cpr=demander("Cpr (pièces)", defaut=0.0),
                        decimales_A=demander("arrondir Å∞ à n décimales (comme le corrigé)", "entier"))


def s_stock():
    return F.stock(demander("λ"), demander("T (mission)"), demander("p0 (probabilité)"),
                   nb_elements=demander("nombre d'éléments en service", "entier", 1))


def s_abc():
    print("  Une machine par ligne : NOM COUT NB_PANNES   (ligne vide pour terminer)")
    d = {}
    while True:
        ligne = input("    > ").strip()
        if not ligne:
            break
        try:
            nom, c, n = ligne.split()
            d[nom] = (nombre(c), nombre(n))
        except Exception as err:  # noqa: BLE001
            print(f"    ! ligne invalide ({err})")
    return F.analyse_abc(d, demander("seuil A %", defaut=80.0), demander("seuil B %", defaut=95.0))


def s_contrainte():
    return F.contrainte_resistance_normales(demander("μR résistance"), demander("σR"), demander("μS contrainte"),
                                            demander("σS"), R_cible=demander("fiabilité cible"))


def s_normale_table():
    u = demander("u (vide pour chercher u à partir de p)")
    if u is not None:
        print(f"  F({F.fmt(u)}) = {F.fmt(F.phi(u))}")
    else:
        p = demander("p = F(u)")
        print(f"  u = {F.fmt(F.phi_inv(p))}")
    return None


MENU = [
    ("M2 — LOIS", None),
    ("Loi exponentielle (R, F, λ, MTTF, Δt, t pour R cible)", s_exponentielle),
    ("Loi normale", s_normale),
    ("Loi de Weibull (2 ou 3 paramètres)", s_weibull),
    ("Loi log-normale", s_lognormale),
    ("Fonction Gamma", s_gamma),
    ("Table loi normale : F(u) ou u", s_normale_table),
    ("M2 — SYSTÈMES", None),
    ("Série / parallèle (fiabilités données)", s_serie_parallele),
    ("Redondance r/m", s_rm),
    ("Série ou parallèle exponentiels (MTTF, R(t))", s_exp_systemes),
    ("Stand-by (modèles 1 à 6)", s_standby),
    ("Système complexe / pont (chemins minimaux)", s_complexe),
    ("M3 — SYSTÈMES RÉPARABLES", None),
    ("Élément réparable : A(t), A∞, M(t)", s_element_reparable),
    ("Temps moyens : MUT, MDT, MTBF -> A∞", s_temps_moyens),
    ("Série / parallèle réparables (disponibilité)", s_serie_par_reparable),
    ("Deux éléments identiques réparables en parallèle", s_deux_paralleles),
    ("Chaîne de Markov quelconque", s_markov),
    ("M3 — MAINTENANCE", None),
    ("Préventive systématique : Rm(T), n minimal", s_preventive),
    ("Période optimale de remplacement t* (CI, ice, icm)", s_remplacement),
    ("Remplacement préventif à âge fixé (Weibull, Cp, Cd)", s_remplacement_age),
    ("Période systématique Weibull (remplace l'ABAQUE)", s_weibull_systematique),
    ("Maintenance corrective : disponibilité optimale", s_corrective),
    ("Stock de pièces de rechange", s_stock),
    ("Analyse ABC", s_abc),
    ("CONTRAINTE–RÉSISTANCE", None),
    ("Deux lois normales", s_contrainte),
]


def solveurs():
    return [(titre, f) for titre, f in MENU if f is not None]


def afficher_menu():
    print("\n" + "=" * 70 + "\n  MEC8312 — RÉSOLUTION RAPIDE\n" + "=" * 70)
    k = 0
    for titre, f in MENU:
        if f is None:
            print(f"\n  {titre}")
        else:
            k += 1
            print(f"   {k:2d}. {titre}")
    print("\n    q. quitter")


def lancer(numero):
    titre, f = solveurs()[numero - 1]
    print(f"\n>>> {titre}")
    try:
        sol = f()
        if sol is not None:
            print()
            print(sol)
    except Exception as err:  # noqa: BLE001
        print(f"\n  ERREUR : {err}")


def main():
    if len(sys.argv) > 1:
        lancer(int(sys.argv[1]))
        return
    while True:
        afficher_menu()
        try:
            choix = input("\nChoix : ").strip().lower()
        except EOFError:
            break
        if choix in ("q", "quit", "exit"):
            break
        if choix.isdigit() and 1 <= int(choix) <= len(solveurs()):
            lancer(int(choix))
            input("\n(Entrée pour revenir au menu)")


if __name__ == "__main__":
    main()
