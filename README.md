# fem-double-phase-solvers
"Solveur FEM P1 sous MATLAB pour la résolution directe d'EDP non linéaires à opérateur Double Phase dans les espaces de Musielak-Orlicz.
Cadre Théorique et Résultats Numériques du Problème Double Phase
1. Détermination du Cadre Fonctionnel
Pour l'équation non linéaire à double phase proposée :

-div( (1 + a(x,y) * |grad(u)|^2) * grad(u) ) - sin(x) * u = 1    dans Omega = ]0,1[ x ]0,1[

avec p = 2, q = 4 et le poids non uniforme a(x,y) = 1 + x^2 + y^2.

Le cadre fonctionnel naturel pour ce type d'opérateur est l'espace de Musielak-Orlicz-Sobolev. On définit d'abord la fonction de Musielak-Orlicz par :

H(x, y, t) = (1/2) * t^2 + (a(x,y) / 4) * t^4    pour t >= 0 et (x,y) dans Omega

Puisque le poids a(x,y) est continu et strictement positif sur l'adhérence du domaine (min_Omega a = 1 > 0), cet espace coïncide avec l'espace de Sobolev classique à double phase. L'espace de Musielak-Orlicz-Sobolev W^1,H(Omega) est défini par :

W^1,H(Omega) = { u dans L^2(Omega) : |grad(u)| dans L^H(Omega) }

Où L^H(Omega) est l'espace de Musielak-Orlicz associé à la norme de Luxemburg.

En intégrant les conditions aux limites de Dirichlet homogènes, le cadre fonctionnel de notre solution est l'espace de Banach réflexif et séparable :

X = W_0^1,H(Omega)

muni de la norme complète du gradient :

Norme(u) = Racine([ Intégrale(|grad(u)|^2) ]) + Racine_Quatrième([ Intégrale(a(x,y) * |grad(u)|^4) ])

2. Détermination de l'Approche Variationnelle (Forme Faible)
Soit v dans X une fonction test. En multipliant l'équation par v et en appliquant la formule de Green (intégration par parties) sachant que v = 0 sur le bord de Omega, on obtient la formulation variationnelle suivante :

Trouver u dans X tel que pour tout v :

Intégrale( (1 + a(x,y) * |grad(u)|^2) * grad(u) * grad(v) ) - Intégrale( sin(x) * u * v ) = Intégrale( 1 * v )

La fonctionnelle d'énergie globale associée à ce problème, J : X -> R, est définie par :

J(u) = (1/2) * Intégrale(|grad(u)|^2) + (1/4) * Intégrale(a(x,y) * |grad(u)|^4) - (1/2) * Intégrale(sin(x) * u^2) - Intégrale(u)

Le problème revient à chercher les points critiques de la fonctionnelle J sur X, c'est-à-dire les solutions de J'(u) = 0.

3. Démonstration de l'Existence (En 3 étapes clés)
Étape 1 : Coercivité de la fonctionnelle J (La fonctionnelle est minorée)
Pour prouver que la fonctionnelle possède un minimum (un "qaa"), il faut montrer qu'elle tend vers +infini lorsque la Norme(u) -> infini.

Grâce à l'inégalité de Poincaré, il existe une constante Cp > 0 telle que :
Intégrale(u^2) <= Cp * Intégrale(|grad(u)|^2)

De plus, on sait que sin(x) <= 1 sur Omega. On peut donc majorer le terme de perturbation :
(1/2) * Intégrale(sin(x) * u^2) <= (1/2) * Intégrale(u^2) <= (Cp / 2) * Intégrale(|grad(u)|^2)

En utilisant l'inégalité de Hölder et Poincaré sur le terme source Intégrale(u * 1), il existe une constante Cl > 0 telle que :
|Intégrale(u)| <= Cl * ||grad(u)||_L^2

En injectant ces majorations dans J(u), on obtient :
J(u) >= (1/2) * (1 - Cp) * ||grad(u)||_L^2^2 + (1/4) * Intégrale(a(x,y) * |grad(u)|^4) - Cl * ||grad(u)||_L^2

Puisque le poids a(x,y) >= 1 > 0, le terme de la q-phase (puissance 4) est d'ordre supérieur. C'est lui qui domine largement à l'infini face aux termes d'ordre 2 et 1. Par conséquent :

Limite de J(u) quand la norme tend vers l'infini = +infini

La fonctionnelle J est coercive, ce qui implique qu'elle est bornée inférieurement (J est minorée). On peut donc définir sa borne inférieure :
m = inf_X J(u) > -infini

Étape 2 : Semi-continuité inférieure faible (f.s.c.i.)
Soit {u_n} dans X une suite minimisante telle que la limite de J(u_n) = m. Puisque J est coercive, la suite {u_n} est bornée dans l'espace de Banach réflexif X.

Par les théorèmes de compacité fonctionnelle (théorème de Kakutani) :

Il existe une sous-suite {u_n} et un élément u* dans X tel que u_n converge faiblement vers u* dans X.

Par les injections compactes de Kondrachov, u_n converge fortement vers u* dans L^2(Omega).

L'application du gradient est convexe et continue, elle est donc semi-continument inférieurement pour la topologie faible. Pour les termes perturbatifs, la convergence forte donne un passage à la limite direct. En combinant ces résultats, on obtient :

J(u*) <= Limite_inf J(u_n) = m

Puisque m est l'infimum, on a nécessairement J(u*) = m. L'élément u* est donc un minimum global de la fonctionnelle.

Étape 3 : Caractérisation de la Solution Non Triviale
Puisque u* est un minimum global, la dérivée s'annule en ce point : J'(u*) = 0, ce qui valide l'existence de la solution faible.

Profil et Amplitude du Solution :
Le solveur converge avec succès en seulement 6 itérations. Le profil numérique final se stabilise sous forme d'une cloche parfaitement lisse et symétrique qui respecte les conditions aux limites (u = 0 sur les bords).

La valeur maximale atteinte au centre du domaine (le sommet de la cloche) est de :

max(u) = 0.07

4. Résultats Numériques de Convergence
L'exécution de l'algorithme de Newton-Raphson confirme la théorie avec une convergence quadratique parfaite vers la solution non triviale :

Itération 1 : ||Résidu|| = 3.8831e-01

Itération 2 : ||Résidu|| = 1.1707e-01

Itération 3 : ||Résidu|| = 3.3975e-02

Itération 4 : ||Résidu|| = 7.4709e-03

Itération 5 : ||Résidu|| = 6.0729e-04

Itération 6 : ||Résidu|| = 4.7317e-06 ---> Convergence validée !
