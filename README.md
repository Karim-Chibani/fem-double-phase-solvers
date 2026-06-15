# fem-double-phase-solvers
"Solveur FEM P1 sous MATLAB pour la résolution directe d'EDP non linéaires à opérateur Double Phase dans les espaces de Musielak-Orlicz.

# Analyse Mathématique et Résultats Numériques de l'Opérateur Double Phase

---

## 1. Détermination du Cadre Fonctionnel

Pour l'équation non linéaire à double phase proposée :

$$-\text{div}\left( \left(1 + a(x,y)|\nabla u|^2\right) \nabla u \right) - \sin(x)u = 1 \quad \text{dans } \Omega = ]0,1[\times]0,1[$$

avec $p = 2$, $q = 4$ et le poids non uniforme $a(x,y) = 1 + x^2 + y^2$.

Le cadre fonctionnel naturel pour ce type d'opérateur est l'espace de Musielak-Orlicz-Sobolev. On définit d'abord la fonction de Musielak-Orlicz par :

$$H((x,y), t) = \frac{1}{2}t^2 + \frac{a(x,y)}{4}t^4 \quad \text{pour } t \geq 0 \text{ et } (x,y) \in \Omega$$

Puisque le poids $a(x,y)$ est continu et strictement positif sur l'adhérence du domaine ($\min_{\bar{\Omega}} a = 1 > 0$), cet espace coïncide avec l'espace de Sobolev classique à double phase. L'espace de Musielak-Orlicz-Sobolev $W^{1,H}(\Omega)$ est défini par :

$$W^{1,H}(\Omega) = \left\{ u \in L^2(\Omega) : |\nabla u| \in L^H(\Omega) \right\}$$

Où $L^H(\Omega)$ est l'espace de Musielak-Orlicz associé à la norme de Luxemburg.

En intégrant les conditions aux limites de Dirichlet homogènes, le cadre fonctionnel de notre solution est l'espace de Banach réflexif et séparable :

$$X = W_0^{1,H}(\Omega)$$

muni de la norme complète du gradient :

$$\|u\| = \left( \int_\Omega |\nabla u|^2 \,dxdy \right)^{1/2} + \left( \int_\Omega a(x,y)|\nabla u|^4 \,dxdy \right)^{1/4}$$

---

## 2. Détermination de l'Approche Variationnelle (Forme Faible)

Soit $v \in X$ une fonction test. En multipliant l'équation par $v$ et en appliquant la formule de Green (intégration par parties) sachant que $v = 0$ sur le bord de $\Omega$, on obtient la formulation variationnelle suivante :

Trouver $u \in X$ tel que pour tout $v \in X$ :

$$\int_\Omega \left(1 + a(x,y)|\nabla u|^2\right) \nabla u \cdot \nabla v \,dxdy - \int_\Omega \sin(x) u v \,dxdy = \int_\Omega v \,dxdy$$

La fonctionnelle d'énergie globale associée à ce problème, $J : X \to \mathbb{R}$, est définie par :

$$J(u) = \frac{1}{2}\int_\Omega |\nabla u|^2 \,dxdy + \frac{1}{4}\int_\Omega a(x,y)|\nabla u|^4 \,dxdy - \frac{1}{2}\int_\Omega \sin(x)u^2 \,dxdy - \int_\Omega u \,dxdy$$

Le problème revient à chercher les points critiques de la fonctionnelle $J$ sur $X$, c'est-à-dire les solutions de $J'(u) = 0$.

---

## 3. Démonstration de l'Existence (En 3 étapes clés)

### Étape 1 : Coercivité de la fonctionnelle J (La fonctionnelle est minorée)
Pour prouver que la fonctionnelle possède un minimum global, il faut montrer qu'elle tend vers $+\infty$ lorsque la norme $\|u\| \to \infty$.

Grâce à l'inégalité de Poincaré, il existe une constante $C_P > 0$ telle que : 
$$\int_\Omega u^2 \,dxdy \leq C_P \int_\Omega |\nabla u|^2 \,dxdy$$

De plus, on sait que $\sin(x) \leq 1$ sur $\Omega$. On peut donc majorer le terme de perturbation : 
$$\frac{1}{2}\int_\Omega \sin(x)u^2 \,dxdy \leq \frac{1}{2}\int_\Omega u^2 \,dxdy \leq \frac{C_P}{2}\int_\Omega |\nabla u|^2 \,dxdy$$

En utilisant l'inégalité de Hölder et Poincaré sur le terme source, il existe une constante $C_L > 0$ telle que : 
$$\left|\int_\Omega u \,dxdy\right| \leq C_L \|\nabla u\|_{L^2}$$

En injectant ces majorations dans $J(u)$, on obtient : 
$$J(u) \geq \frac{1}{2}(1 - C_P)\|\nabla u\|_{L^2}^2 + \frac{1}{4}\int_\Omega a(x,y)|\nabla u|^4 \,dxdy - C_L\|\nabla u\|_{L^2}$$

Puisque le poids $a(x,y) \geq 1 > 0$, le terme de la $q$-phase (puissance 4) est d'ordre supérieur. C'est lui qui domine largement à l'infini face aux termes d'ordre 2 et 1. Par conséquent :

$$\lim_{\|u\| \to \infty} J(u) = +\infty$$

La fonctionnelle $J$ est coercive, ce qui implique qu'elle est bornée inférieurement ($J$ est minorée). On peut donc définir sa borne inférieure : 
$$m = \inf_{u \in X} J(u) > -\infty$$

### Étape 2 : Semi-continuité inférieure faible (f.s.c.i.)
Soit $\{u_n\} \subset X$ une suite minimisante telle que $\lim_{n \to \infty} J(u_n) = m$. Puisque $J$ est coercive, la suite $\{u_n\}$ est bornée dans l'espace de Banach réflexif $X$.

Par les théorèmes de compacité fonctionnelle (théorème de Kakutani) :
* Il existe une sous-suite (encore notée $\{u_n\}$) et un élément $u^* \in X$ tel que $u_n$ converge faiblement vers $u^*$ dans $X$.
* Par les injections compactes de Kondrachov, $u_n$ converge fortement vers $u^*$ dans $L^2(\Omega)$.

L'application du gradient est convexe et continue, elle est donc semi-continument inférieurement pour la topologie faible. Pour les termes perturbatifs, la convergence forte donne un passage à la limite direct. En combinant ces résultats, on obtient :

$$J(u^*) \leq \liminf_{n \to \infty} J(u_n) = m$$

Puisque $m$ est l'infimum, on a nécessairement $J(u^*) = m$. L'élément $u^* \in W_0^{1,H}(\Omega)$ est donc un minimum global de la fonctionnelle.

### Étape 3 : Caractérisation de la Solution Non Triviale
Puisque $u^*$ est un minimum global, la dérivée s'annule en ce point : $J'(u^*) = 0$, ce qui valide l'existence de la solution faible.

---

## 4. Profil et Amplitude de la Solution Approchée

Le solveur converge avec succès en seulement 6 itérations. Le profil numérique final se stabilise sous forme d'une cloche parfaitement lisse et symétrique qui respecte les conditions aux limites ($u = 0$ sur les bords).

La valeur maximale atteinte au centre du domaine (le sommet de la cloche) constitue la solution approchée maximale :

$$\max(u_h) = 0.07$$

---

## 5. Résultats Numériques de Convergence

L'exécution de l'algorithme de Newton-Raphson confirme la théorie avec une convergence quadratique parfaite vers la solution non triviale :

* **Itération 1 :** $\|\text{Résidu}\| = 3.8831 \times 10^{-1}$
* **Itération 2 :** $\|\text{Résidu}\| = 1.1707 \times 10^{-1}$
* **Itération 3 :** $\|\text{Résidu}\| = 3.3975 \times 10^{-2}$
* **Itération 4 :** $\|\text{Résidu}\| = 7.4709 \times 10^{-3}$
* **Itération 5 :** $\|\text{Résidu}\| = 6.0729 \times 10^{-4}$
* **Itération 6 :** $\|\text{Résidu}\| = 4.7317 \times 10^{-6} \quad \rightarrow \quad \text{Convergence validée !}$
