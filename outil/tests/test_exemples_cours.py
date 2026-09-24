"""Chaque test reproduit un exemple résolu des notes de cours (M2, M3, exemples Topic 3).
Lancer :  python -m unittest discover -s outil/tests   (depuis la racine du dépôt)"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fiabilite import *  # noqa: E402,F403


class LoisM2(unittest.TestCase):
    def test_exponentielle_diapo27(self):
        s = exponentielle(mttf=100, t=[100, 200, 10, 1])
        self.assertAlmostEqual(s["F(100)"], 0.6321, 4)
        self.assertAlmostEqual(s["R(200)"], 0.135, 3)
        self.assertAlmostEqual(s["F(10)"], 0.095, 3)
        self.assertAlmostEqual(s["R(1)"], 0.99005, 5)  # la diapo écrit 0,9905 (coquille)

    def test_normale_diapo34(self):
        s = normale(100, 50, t=[100, 160, 40])
        self.assertAlmostEqual(s["F(100)"], 0.5, 6)
        self.assertAlmostEqual(s["F(160)"], 0.8849, 4)
        self.assertAlmostEqual(s["F(40)"], 0.1151, 4)

    def test_weibull_diapos38_40(self):
        s = weibull(b=0.5, theta=1, t=2)
        self.assertAlmostEqual(s["MTTF"], 2, 9)
        self.assertAlmostEqual(s["F(2)"], 0.7569, 4)
        self.assertAlmostEqual(weibull(b=1, theta=1, t=1)["F(1)"], 0.6321, 4)
        s = weibull(b=3.5, theta=1)
        self.assertAlmostEqual(s["MTTF"], 0.9, 2)
        # la diapo 40 écrit F = 0,5008 : c'est en fait R(0,9) ; F(0,9) = 0,4992
        self.assertAlmostEqual(weibull(b=3.5, theta=1, t=0.9)["R(0,9)"], 0.5008, 4)

    def test_gamma_diapo38(self):
        self.assertAlmostEqual(gamma(2.5)["Γ"], 1.329, 3)

    def test_lognormale_diapo47(self):
        s = lognormale(m=10, sigma=2, t=5)
        self.assertAlmostEqual(s["σy"], 0.198, 3)
        self.assertAlmostEqual(s["μy"], 2.283, 3)
        self.assertAlmostEqual(s["R(5)"], 0.99966, 4)
        self.assertAlmostEqual(normale(10, 2, t=5)["R(5)"], 0.9938, 4)

    def test_intervalle_diapos48_49(self):
        self.assertAlmostEqual(exponentielle(lam=1e-4, t1=500, t2=700)["F(Δt)"], 0.0198, 4)
        self.assertAlmostEqual(weibull(beta=1.7, eta=2000, t1=1000, t2=1100)["F(Δt)"], 0.0527, 3)


class SystemesM2(unittest.TestCase):
    def test_pont_diapos80_83(self):
        self.assertAlmostEqual(pont(0.99, 0.99, 0.99, 0.99, 0.99)["Rs"], 0.9997981, 6)
        R = {k: 0.99 for k in "ABCDE"}
        chemins = [["A", "C"], ["B", "D"], ["A", "E", "D"], ["B", "E", "C"]]
        self.assertAlmostEqual(systeme_chemins(chemins, R)["Rs"], 0.9997981, 6)
        coupes = [["A", "B"], ["C", "D"], ["A", "E", "D"], ["B", "E", "C"]]
        self.assertAlmostEqual(systeme_coupes(coupes, R)["Rs"], 0.9997981, 6)

    def test_r_sur_m_limites(self):
        self.assertAlmostEqual(r_sur_m(0.9, 3, 3)["R_3/3"], serie([0.9] * 3)["Rss"], 12)
        self.assertAlmostEqual(r_sur_m(0.9, 1, 3)["R_1/3"], parallele([0.9] * 3)["Rps"], 12)

    def test_standby(self):
        lam, t = 0.001, 500
        self.assertAlmostEqual(standby(t, lam)["Rsb(500)"], math.exp(-0.5) * 1.5, 12)
        self.assertAlmostEqual(standby(t, lam)["MTTFsb"], 2000, 9)
        self.assertAlmostEqual(standby(t, lam, n=2)["MTTFsb"], 3000, 9)
        self.assertAlmostEqual(standby(t, lam, Rs=0.9)["MTTFsb"], 1900, 9)
        s = standby(t, 0.001, 0.002)
        self.assertAlmostEqual(s["MTTFsb"], 1500, 9)
        # Modèle 6 : MTTF analytique = intégrale numérique de R(t)
        s6 = standby(None, 0.001, 0.002, lamBsb=0.0005)
        R6 = lambda x: standby(x, 0.001, 0.002, lamBsb=0.0005)[f"Rsb({fmt(x)})"]
        self.assertAlmostEqual(mttf_numerique(R6, 1000) / s6["MTTFsb"], 1, 4)


class ReparablesM3(unittest.TestCase):
    def test_deux_pompes_parallele(self):
        s = parallele2_reparable(mttf=500, mttr=8, t=150)
        self.assertAlmostEqual(s["Rps(150)"], 0.9909, 4)
        self.assertAlmostEqual(s["Rps_exact(150)"], s["Rps(150)"], 3)
        self.assertAlmostEqual(s["Aps(150)"], 0.9997, 3)
        self.assertAlmostEqual(s["Aps∞"], 0.9997, 3)
        self.assertAlmostEqual(s["MTTFps"], 16375, 6)
        self.assertAlmostEqual(s["Rps_non_rep(150)"], 0.9328, 4)
        self.assertAlmostEqual(s["MTTFps_non_rep"], 750, 9)

    def test_markov_generique_equivaut_deux_pompes(self):
        l, m = 0.002, 0.125
        tr = {("2", "1"): 2 * l, ("1", "2"): m, ("1", "0"): l, ("0", "1"): 2 * m}
        s = chaine_markov(tr, "2", ["0"], t=150)
        ref = parallele2_reparable(l, m, t=150)
        self.assertAlmostEqual(s["R(150)"], ref["Rps_exact(150)"], 9)
        self.assertAlmostEqual(s["MTTF"], ref["MTTFps"], 6)
        # avec 1 réparateur par pompe, A∞ = 1 - (λ/(λ+μ))²
        self.assertAlmostEqual(s["A∞"], ref["Aps∞"], 9)
        self.assertAlmostEqual(s["A(150)"], ref["Aps(150)"], 9)

    def test_markov_3_processeurs_2_sources(self):
        l1, l2 = 1e-3, 2e-3
        tr = {("3,2", "2,2"): 3 * l1, ("2,2", "1,2"): 2 * l1, ("1,2", "F1"): l1,
              ("3,1", "2,1"): 3 * l1, ("2,1", "1,1"): 2 * l1, ("1,1", "F1"): l1,
              ("3,2", "3,1"): 2 * l2, ("2,2", "2,1"): 2 * l2, ("1,2", "1,1"): 2 * l2,
              ("3,1", "F2"): l2, ("2,1", "F2"): l2, ("1,1", "F2"): l2}
        t = 100
        s = chaine_markov(tr, "3,2", ["F1", "F2"], t=t)
        Rp = 1 - (1 - math.exp(-l1 * t)) ** 3
        Rse = 1 - (1 - math.exp(-l2 * t)) ** 2
        self.assertAlmostEqual(s[f"R({t})"], Rp * Rse, 10)

    def test_element_reparable(self):
        s = element_reparable(mttf=27400, mttr=160, t=0)
        self.assertAlmostEqual(s["A∞"], 27400 / 27560, 12)
        self.assertAlmostEqual(s["A(0)"], 1, 12)


class MaintenanceM3(unittest.TestCase):
    def test_preventive_systematique_manuscrit(self):
        s = preventive_systematique(taux_lineaire(5e-4, 5e-5), 100, n=[2], R_cible=0.9)
        self.assertAlmostEqual(s["R(T) sans maintenance"], 0.7408, 4)
        self.assertAlmostEqual(s["Rm(T) n=2"], 0.8751, 3)
        self.assertEqual(s["n_min"], 4)
        s2 = preventive_systematique(taux_expression("5e-4 + 5e-5*t"), 100, n=4)
        # le manuscrit écrit 0,9512 x 0,9512 = 0,9045 ; le produit exact vaut 0,9048
        self.assertAlmostEqual(s2["Rm(T) n=4"], 0.9048, 4)

    def test_periode_remplacement(self):
        self.assertAlmostEqual(periode_optimale_remplacement(10000, 400, 200)["t*"], 5.8, 1)

    def test_corrective_freins(self):
        s = corrective(0.994, 200, 27000, 5, decimales_A=3)
        self.assertAlmostEqual(s["K1"], 1.2072, 4)
        self.assertAlmostEqual(s["K2"], 135000, 6)
        self.assertAlmostEqual(s["κ"], 111829, -1)
        self.assertAlmostEqual(s["Cmv"], 401, 0)
        self.assertAlmostEqual(s["CT actuel"], 1014.9, 1)
        self.assertAlmostEqual(s["CT optimal"], 807.2, 0)
        self.assertAlmostEqual(s["ΔCP"], 408.7, 0)

    def test_stock_pompe(self):
        s = stock(0.0007, 3000, 0.975)
        self.assertAlmostEqual(s["SN"], 4.94, 2)
        self.assertEqual(s["stock approx"], 5)

    def test_abaque_weibull(self):
        # diapo 73 : β = 3, r = 10 -> x0 ≈ 0,35, C2/C1 ≈ 0,36
        s = weibull_systematique(3, 10, 10)
        self.assertAlmostEqual(s["x0"], 0.35, 1)
        self.assertAlmostEqual(s["C2/C1"], 0.36, 1)
        # diapo 77 : β = 2, η = 15000 h ; r = 20 -> 0,38 à 0,23 ; r = 5 -> 0,67 à 0,46
        s = weibull_systematique(2, 15000, 20)
        self.assertAlmostEqual(s["x0"], 0.23, 1)
        self.assertAlmostEqual(s["C2/C1"], 0.38, 1)
        s = weibull_systematique(2, 15000, 5)
        self.assertAlmostEqual(s["x0"], 0.46, 1)
        self.assertAlmostEqual(s["C2/C1"], 0.67, 1)

    def test_analyse_abc(self):
        d = {1: (100, 4), 2: (32, 14), 3: (50, 4), 4: (19, 15), 5: (4, 3), 6: (30, 8),
             7: (40, 12), 8: (80, 2), 9: (55, 3), 10: (150, 5), 11: (160, 4), 12: (5, 3),
             13: (10, 8), 14: (20, 8)}
        table = analyse_abc(d)["table"]
        self.assertEqual([r[0] for r in table][:3], [11, 10, 1])
        self.assertAlmostEqual(table[0][3], 21.2, 1)
        self.assertAlmostEqual(table[5][3], 78.8, 1)
        self.assertAlmostEqual(table[5][6], 23.7, 1)  # le cours arrondit à 23,6
        self.assertEqual(table[-1][2], 755)

    def test_remplacement_age_weibull(self):
        s = remplacement_age(Weibull(3, 1000), Cp=100, Cd=1000)
        self.assertLess(s["C(θ*)"], s["C_corrective"])


class ContrainteResistance(unittest.TestCase):
    def test_normales(self):
        s = contrainte_resistance_normales(500, 40, 400, 30)
        self.assertAlmostEqual(s["z"], 2.0, 9)
        self.assertAlmostEqual(s["R"], 0.97725, 5)
        g = contrainte_resistance_generale(Normale(500, 40), Normale(400, 30), 200, 600)
        self.assertAlmostEqual(g["R"], s["R"], 6)


if __name__ == "__main__":
    unittest.main()
