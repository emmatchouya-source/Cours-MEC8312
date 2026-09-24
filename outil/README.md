# Outil de résolution rapide — MEC8312

Boîte à outils Python **sans aucune dépendance** (Python 3.8+ standard, rien à installer,
fonctionne hors ligne) qui résout les exercices du cours et affiche un **corrigé détaillé
étape par étape**, avec les notations du professeur.

## Trois façons de l'utiliser

### 1. Menu interactif (le plus rapide à l'examen)

```bash
python3 outil/examen.py        # menu
python3 outil/examen.py 21     # ouvre directement le solveur n° 21 (maintenance corrective)
```

Saisie : `0,994`, `0.994`, `5e-4`, `1/500`, `sqrt(2)` sont acceptés. Listes : `100; 200; 500`.
Laisser vide = donnée absente.

### 2. Python / Jupyter (le plus souple)

```python
import sys; sys.path.insert(0, "outil")
from fiabilite import *

weibull(b=2, theta=1000, t=[200, 500], R_cible=0.9)
corrective(A_actuel=0.994, Cmv_actuel=200, MTTF=27000, c=5)
sol = stock(lam=0.0007, T=3000, p0=0.975)
sol["SN"]            # -> 4.94
```

### 3. Avec Claude

Ouvrir le dépôt dans Claude Code (ou coller l'énoncé) : le fichier `CLAUDE.md` à la racine lui
explique comment identifier le type d'exercice et appeler le bon solveur.

## Contenu

| Module | Couvre |
|---|---|
| `lois.py` | exponentielle, normale, Weibull (2-3 param.), log-normale, Γ, F(Δt), t pour R cible |
| `systemes.py` | série, parallèle, r/m, stand-by (modèles 1 à 6), pont, système quelconque par chemins ou coupes minimaux |
| `markov.py` | élément réparable A(t)/A∞/M(t), temps moyens, série/parallèle réparables, 2 éléments en parallèle (exact + approx.), **chaîne de Markov quelconque** (équations, P(t), R(t), MTTF, A∞) |
| `maintenance.py` | préventive systématique Rm(T) + n minimal, remplacement à âge fixé, t* optimal, **optimisation Weibull sans abaque**, corrective (K1, K2, κ, Å∞, Cmv, coûts), stock (normal + Poisson exact), analyse ABC |
| `contrainte.py` | contrainte–résistance (normales, lois quelconques) — à valider quand le chapitre sera vu |

Voir `AIDE_MEMOIRE.md` pour reconnaître rapidement chaque type d'exercice.

## Validation

Chaque exemple résolu des notes est rejoué comme test :

```bash
python3 -m unittest discover -s outil/tests -v
```

Au passage, les tests ont relevé trois coquilles dans les notes :
- diapo 27 (M2) : R(0,01·MTTF) = 0,99005 (et non 0,9905) ;
- diapo 40 (M2) : 0,5008 est R(MTTF), pas F(MTTF) (F = 0,4992) ;
- note manuscrite (maintenance systématique) : 0,9512 × 0,9512 = 0,9048 (et non 0,9045) —
  la conclusion n = 4 reste vraie.

## Checklist avant l'examen

1. Copier le dossier `outil/` sur l'ordinateur de l'examen et lancer une fois `python3 outil/examen.py`.
2. Lancer les tests (doivent afficher `OK`).
3. Refaire 2-3 exemples du cours avec le menu pour prendre les réflexes de saisie.
4. Ajouter les nouveaux chapitres (sécurité, risques, allocation de redondance) au fur et à mesure.
