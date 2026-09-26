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


# Dimensions (comme dans la calculatrice) : T durée, 1/T taux, 1/T2, $ coût, $/T, $/T2,
# P probabilité, - sans dimension, #mot nombre entier (#pièce, #élément…).
_PLURIEL = {"an": "ans", "cycle": "cycles"}
_DIMS_TEMPS = {"T", "1/T", "1/T2", "$/T", "$/T2"}


def unite(dim, unite_temps="h", valeur=None):
    """Texte de l'unité d'une grandeur de dimension `dim` quand le temps est
    exprimé en `unite_temps` (h, j, an, km, cycle…) ; `valeur` accorde le
    pluriel (1,5 an, 5,8 ans). Un texte qui n'est pas une dimension est rendu tel quel."""
    if not dim or dim in ("P", "-"):
        return ""
    pluriel = dim.startswith("#") if valeur is None else not abs(valeur) < 2
    if dim.startswith("#"):
        mot = dim[1:]
        return mot + "s" if pluriel and not mot.endswith(("s", "x")) else mot
    ut = _PLURIEL.get(unite_temps, unite_temps) if pluriel else unite_temps
    return {"T": ut, "1/T": unite_temps + "⁻¹", "1/T2": unite_temps + "⁻²", "$": "$",
            "$/T": "$/" + unite_temps, "$/T2": "$/" + unite_temps + "²"}.get(dim, dim)


class Solution:
    def __init__(self, titre, unite_temps="h"):
        self.titre = titre
        self.unite_temps = unite_temps
        self.etapes = []
        self.resultats = {}
        self.unites = {}

    def donnee(self, texte):
        self.etapes.append(("donnee", texte))

    def u(self, dim, valeur=None):
        """Unité d'une dimension dans l'unité de temps de la solution."""
        return unite(dim, self.unite_temps, valeur)

    def q(self, valeur, dim):
        """Valeur avec son unité, ex. sol.q(500, "T") -> '500 h'."""
        u = self.u(dim, valeur)
        return fmt(valeur) + (" " + u if u else "")

    def etape(self, texte, valeur=None, nom=None, unite=""):
        """Ajoute une ligne de calcul ; si `nom` est donné, la valeur est
        aussi enregistrée comme résultat. `unite` : dimension ("T", "1/T",
        "$", "$/T"…) convertie selon l'unité de temps, ou texte libre."""
        u = self.u(unite, valeur if isinstance(valeur, (int, float)) else None) if unite else ""
        if valeur is not None:
            texte = f"{texte} = {fmt(valeur)}{(' ' + u) if u else ''}"
        self.etapes.append(("calcul", texte))
        if nom is not None:
            self.resultats[nom] = valeur
            if u:
                self.unites[nom] = u

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
                u = self.unites.get(k, "")
                lignes.append(f"  {k:<28} {fmt(v)}{(' ' + u) if u else ''}")
        return "\n".join(lignes)

    __str__ = texte

    def __repr__(self):
        return self.texte()
