# Dépôt de préparation — MEC8312 Fiabilité et sécurité des systèmes (Polytechnique Montréal)

Notes de cours dans `Topic 1/`, `Topic 2/`, `Topic 3/` ; boîte à outils de résolution dans `outil/`.

## Résoudre un exercice

1. Identifier le type d'exercice avec `outil/AIDE_MEMOIRE.md`.
2. Appeler le solveur correspondant du paquet `outil/fiabilite` plutôt que de calculer à la main :
   ```bash
   cd outil && python3 -c "from fiabilite import *; print(corrective(0.994, 200, 27000, 5))"
   ```
3. Présenter la réponse en français, dans l'ordre des étapes affichées (données → formule →
   application numérique → résultat → conclusion), avec les notations du cours (λ, μ, b/θ, A∞, Å∞, MTTF…).
4. Si un exercice sort du périmètre des solveurs, le résoudre en suivant les formules des notes
   (voir les PDF) et proposer d'ajouter le solveur manquant dans `outil/fiabilite` avec un test.

## Conventions

- Python standard uniquement (pas de numpy/scipy) : l'outil doit tourner hors ligne sans installation.
- Chaque solveur renvoie une `Solution` (`rapport.py`) : `sol.etape(texte, valeur, nom)`.
- Tout nouvel exemple résolu du cours doit être ajouté à `outil/tests/test_exemples_cours.py`.
- Tests : `python3 -m unittest discover -s outil/tests`.
