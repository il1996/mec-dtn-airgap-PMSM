# R5 — note de livraison

> Sortie : `MEC_BLDC/R5_nsh_out.txt`
> Script : `MEC_BLDC/RUN_R5_NSH.m`
> Exécuté le 9 août 2026. PMSM 15/14, $M_s = 540$, $n_{st}=4$, $n_y=3$,
> $\phi = 0$, base P0, solveur linéaire $\mu_r = 3000$, $k_f = 0{,}3250$.

## Tranché : $n_{sh} = 2$. La Table 12 est juste, le texte ne l'est pas

**Aucun chiffre n'a été transcrit dans un sens ni dans l'autre.** Chaque
valeur de comparaison est **relue du `.tex` par expression régulière**, chaque
valeur de référence est **recalculée par la chaîne**. Le script confronte les
deux ; il ne les recopie pas.

### La garde de la v4, exécutée

Le compte d'inconnues tranche seul, et il tranche sans ambiguïté :

| $n_{sh}$ | $L_s$ | rangées | inconnues | $B_{g1}$ (T) |
|---|---|---|---|---|
| 1 | 9 | 10 | **5400** | 1,07787 |
| **2** | **11** | **12** | **6480** | **1,07908** |
| 3 | 13 | 14 | 7560 | 1,07937 |
| 4 | 15 | 16 | 8640 | 1,07947 |

$L_s = 2n_{sh} + n_{st} + n_y$. Les deux comptes déclarés par la note de la
Table 12 — 5400 et 6480 — **sont exactement ceux que `mesh_bldc` produit à
$n_{sh} = 1$ et $n_{sh} = 2$**. À $n_{sh} = 4$ le compte serait 8640.

### La garde indépendante : $B_{g1}$ désigne la même configuration

Le compte pourrait être juste et la colonne venir d'ailleurs. Second test,
sur la grandeur elle-même, relue de la Table 12 :

| $n_{sh}$ | $B_{g1}$ chaîne | écart / col. 1 | écart / col. 2 |
|---|---|---|---|
| 1 | 1,07787 | **+1,75e−06** | −1,21e−03 |
| 2 | 1,07908 | +1,21e−03 | **+4,96e−06** |
| 4 | 1,07947 | +1,60e−03 | +3,94e−04 |

**Colonne 1 ← $n_{sh}=1$, colonne 2 ← $n_{sh}=2$**, aux arrondis de la
cinquième décimale. À $n_{sh}=4$ l'écart à la colonne 2 vaut 3,94e−04, **quatre-vingts fois la tolérance**.

**GARDE PASSÉE.** Les deux tests concordent.

*À déclarer sans l'arrondir : l'écart de la colonne 2 vaut 4,96e−06 contre une
tolérance de 5e−06 — il passe de justesse.* Il est cohérent avec ce
qu'imprime `A1_table7_out.txt`, mais c'est le compte d'inconnues, exact et
entier, qui porte la conclusion ; $B_{g1}$ la confirme, il ne la fonde pas.

## Le manuscrit est juste partout sauf à une ligne

Les six occurrences de $n_{sh}$, relevées automatiquement :

| ligne | valeur | statut |
|---|---|---|
| 1656 | 4 | ✅ note de la Table 10 — « fifteen radial layers » ⇒ $L_s = 15$ ⇒ $n_{sh}=4$ |
| 1663 | 4 | ✅ identification, résidu 0,06 pt sur six cellules |
| 1690 | 1 | ✅ note de la Table 11 — neuf couches |
| 2014–2015 | 1 et 2 | ✅ **5400 et 6480, vérifiés** |
| **2358** | **1 et 4** | ❌ **le défaut** |

Les lignes « $M_s$ & inconnues » des Tables 10 et 11, relues elles aussi,
impliquent $n_{sh} \in \{1, 4\}$ : 5760/8640/14400 à seize rangées
($n_{sh}=4$), 1800/3600/5400/9000 à dix ($n_{sh}=1$). **Les trois tableaux
sont exacts et mutuellement cohérents.**

## Correction à appliquer — aucun `.tex` n'a été modifié

**Fichier** `article/ArticleI_DtN_PMSM.tex`, **ligne 2358** :

```
avant : between $n_{\mathrm{sh}}=1$ and $n_{\mathrm{sh}}=4$ and both columns are
après : between $n_{\mathrm{sh}}=1$ and $n_{\mathrm{sh}}=2$ and both columns are
```

**Il y a un choix à faire, et il vous revient.** La phrase porte deux
affirmations distinctes, et **elles n'ont pas le même statut** :

- *« the shoe layering is varied between 1 and 4 »* — **vraie**. Le balayage
  d'identification couvre bien 1 à 4 (`X1_table5b_reconcile_out.txt`), et la
  Table 10 est publiée à $n_{sh}=4$.
- *« both columns are reported throughout »* — **fausse à 4**. Les deux
  colonnes portées d'un bout à l'autre sont celles de la Table 12, donc 1 et 2.

Le remplacement mécanique ci-dessus rend la phrase vraie en **rétrécissant**
l'intervalle de sensibilité annoncé, alors que le papier en explore un plus
large. Une seconde rédaction les garde toutes deux — *« varied from
$n_{sh}=1$ to $n_{sh}=4$ (Table 10), the columns at $n_{sh}=1$ and
$n_{sh}=2$ being reported side by side throughout »* — mais **c'est une
réécriture, pas une correction de chiffre**, et je ne la fais pas à votre
place.

## Rien à régénérer

La v4 prévoyait *« régénérer le tableau si nécessaire »*. **Ce n'est pas
nécessaire** : la chaîne reproduit les deux colonnes publiées de la Table 12,
ce que la garde établit. Le défaut est éditorial, pas numérique.

## Ce que R5 ne remet pas en cause

**La ligne 3318** — *« 41,7 points de pourcentage à une couche par
sous-région de bec et 42,6 à quatre »* — **est juste** et ne doit pas être
touchée : à $M_s = 900$, $n_{sh}=1$ donne −9,7 % / +32,0 % (41,7 pt) et
$n_{sh}=4$ donne −18,8 % / +23,8 % (42,6 pt), les deux sortant du même
balayage. C'est un des rares endroits où $n_{sh}=4$ est cité à bon droit.
