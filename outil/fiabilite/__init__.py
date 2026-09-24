"""Boîte à outils MEC8312 — Fiabilité et sécurité des systèmes.

Utilisation rapide (Python, Jupyter ou `python -i`) :

    from fiabilite import *
    weibull(b=2, theta=1000, t=500)
    corrective(A_actuel=0.994, Cmv_actuel=200, MTTF=27000, c=5)

Chaque fonction renvoie une Solution qui s'affiche comme un corrigé détaillé ;
les résultats sont accessibles par sol["nom"].
"""
from .contrainte import general as contrainte_resistance_generale
from .contrainte import normale_normale as contrainte_resistance_normales
from .lois import (Exponentielle, LogNormale, Normale, Weibull, exponentielle, gamma,
                   lognormale, mttf_numerique, normale, weibull)
from .maintenance import (Taux, analyse_abc, corrective, cout_total_annuel,
                          periode_optimale_remplacement, preventive_systematique,
                          remplacement_age, stock, taux_constant, taux_expression,
                          taux_lineaire, taux_weibull, weibull_systematique)
from .markov import (chaine_markov, element_reparable, parallele2_reparable,
                     parallele_reparable, serie_reparable, temps_moyens)
from .numerique import phi, phi_inv
from .rapport import Solution, fmt
from .systemes import (parallele, parallele_exp, pont, r_sur_m, r_sur_m_exp, serie,
                       serie_exp, standby, systeme_chemins, systeme_coupes)
