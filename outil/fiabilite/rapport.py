"""Présentation des solutions : chaque solveur renvoie un objet Solution
qui contient les étapes de calcul (style corrigé) et les résultats."""


def fmt(x, chiffres=6):
    """Formate un nombre à la française (virgule décimale)."""
    if isinstance(x, bool):
        return "oui" if x else "non"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.{chiffres}g}".replace(".", ",")
    return str(x)


class Solution:
    def __init__(self, titre):
        self.titre = titre
        self.etapes = []
        self.resultats = {}

    def donnee(self, texte):
        self.etapes.append(("donnee", texte))

    def etape(self, texte, valeur=None, nom=None, unite=""):
        """Ajoute une ligne de calcul ; si `nom` est donné, la valeur est
        aussi enregistrée comme résultat."""
        if valeur is not None:
            texte = f"{texte} = {fmt(valeur)}{(' ' + unite) if unite else ''}"
        self.etapes.append(("calcul", texte))
        if nom is not None:
            self.resultats[nom] = valeur

    def section(self, texte):
        self.etapes.append(("section", texte))

    def note(self, texte):
        self.etapes.append(("note", texte))

    def __getitem__(self, cle):
        return self.resultats[cle]

    def __contains__(self, cle):
        return cle in self.resultats

    def texte(self):
        lignes = ["=" * 70, self.titre, "=" * 70]
        for genre, t in self.etapes:
            if genre == "section":
                lignes += ["", f"--- {t} ---"]
            elif genre == "donnee":
                lignes.append(f"  Donnée : {t}")
            elif genre == "note":
                lignes.append(f"  >> {t}")
            else:
                lignes.append(f"  {t}")
        if self.resultats:
            lignes += ["", "RÉSULTATS :"]
            for k, v in self.resultats.items():
                if isinstance(v, (list, tuple, dict)):
                    continue
                lignes.append(f"  {k:<28} {fmt(v)}")
        return "\n".join(lignes)

    __str__ = texte

    def __repr__(self):
        return self.texte()
