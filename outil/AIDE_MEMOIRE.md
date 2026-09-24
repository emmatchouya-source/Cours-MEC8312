# Aide-mémoire : reconnaître l'exercice → choisir le solveur

Numéro = option du menu `python3 outil/examen.py`.

| L'énoncé contient… | Type | Menu | Fonction Python |
|---|---|---|---|
| λ constant, « MTTF », R(t) | Loi exponentielle | 1 | `exponentielle(lam=, mttf=, t=)` |
| m, σ, « loi normale », usure | Loi normale (u = (t−m)/σ) | 2 | `normale(m, sigma, t=)` |
| b/β (forme), θ/η (échelle), t0 | Weibull | 3 | `weibull(b=, theta=, t0=, t=)` |
| « log-normale », m et σ de T | Log-normale | 4 | `lognormale(m=, sigma=, t=)` |
| « entre t1 et t2 », « sachant qu'il a fonctionné jusqu'à t1 » | F(Δt) = 1 − R(t2)/R(t1) | 1-4 | param. `t1=, t2=` |
| « durée de vie garantie », « pour R = 0,9 » | t pour R cible | 1-4 | `R_cible=` |
| R1, R2… en série / parallèle | Structure | 7 | `serie([...])`, `parallele([...])` |
| « au moins r parmi m » | Redondance r/m | 8 | `r_sur_m(R, r, m)` |
| k éléments identiques exponentiels, MTTF système | Série/parallèle expo | 9 | `serie_exp`, `parallele_exp` |
| « en réserve », « stand-by », commutateur Rs | Stand-by modèles 1-6 | 10 | `standby(t, lamA, lamB, Rs, n, lamBsb)` |
| schéma en pont, « système complexe » | Chemins minimaux | 11 | `pont(...)`, `systeme_chemins(...)` |
| λ et μ (ou MTTF et MTTR), A(t), A∞ | Élément réparable | 12 | `element_reparable` |
| MUT, MDT, MTBF | Temps moyens | 13 | `temps_moyens` |
| plusieurs éléments réparables, un réparateur chacun | Série/parallèle réparable | 14 | `serie_reparable`, `parallele_reparable` |
| **deux** pompes identiques réparables en parallèle | Markov 3 états | 15 | `parallele2_reparable(mttf=, mttr=, t=)` |
| graphe d'états, « processeurs / sources d'énergie » | Markov quelconque | 16 | `chaine_markov(transitions, initial, pannes, t)` |
| λ(t) = a + b·t, « n maintenances sur (0,T) » | Préventive systématique | 17 | `preventive_systematique(taux_lineaire(a,b), T, n=, R_cible=)` |
| CI, ice, icm, « période optimale de remplacement » | t* = √(2CI/(ice+icm)) | 18 | `periode_optimale_remplacement` |
| Cp, Cd, remplacement à âge θ | C(θ) minimal | 19 | `remplacement_age(Weibull(b,θ), Cp, Cd)` |
| β, η, r = P/p, « abaque » | Période systématique Weibull | 20 | `weibull_systematique(beta, eta, r)` |
| A'∞, C'mv, MTTF, c | Corrective / dispo. optimale | 21 | `corrective(A, Cmv, MTTF, c)` |
| λ, durée T, probabilité p0, « pièces de rechange » | Stock | 22 | `stock(lam, T, p0)` |
| tableau machines / heures d'arrêt / pannes | Analyse ABC | 23 | `analyse_abc({...})` |
| résistance et contrainte normales | Contrainte–résistance | 24 | `contrainte_resistance_normales` |

## Pièges fréquents

- **Unités** : λ en pannes/h avec t en h (ou km/km). MTTR = 1/μ.
- **Γ(x)** : Γ(x+1) = x·Γ(x) ; si x entier, Γ(x) = (x−1)!.
- **Loi normale, u < 0** : F(−u) = 1 − F(u).
- **Maintenance préventive** utile seulement si λ(t) croissant (b > 1) ; inutile si constant ; nuisible si décroissant.
- **Corrective** : si Å∞ > A'∞ → intensifier (Cmv augmente) ; sinon réduire. Le corrigé arrondit Å∞ à 3 décimales
  (option `decimales_A=3`) — sans arrondi Cmv = 404 $ au lieu de 401 $.
- **Stock** : arrondir SN à l'entier supérieur ; le calcul Poisson exact confirme.
- **Parallèle réparable** : la formule du cours R ≈ exp[−2λ²t/(μ+3λ)] suppose λ ≪ μ ; l'outil donne aussi la valeur exacte.
