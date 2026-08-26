# Bloc G — clôture

**Sorties : `outputs/MEC_BLDC/G1_gamma_out.txt` … `G5_count_out.txt`, 26 août 2026.**

Cinq blocs écrits pour clore les points que la relecture avait laissés ouverts et
qui ne pouvaient pas être tranchés depuis le manuscrit seul. Chacun porte sa
garde ; **aucune n'a échoué**.

| bloc | ce qu'il tranche | garde | issue |
|---|---|---|---|
| `RUN_G1_GAMMA` | la convention sur γ et la croissance du crochet | γ est additif, donc la **pente par décade doit être identique** dans les deux conventions | \|Δ\| = 4,4×10⁻¹⁶ — **passée** |
| `RUN_G1_GAMMA` | idem | à N_h = M_s/2 la chaîne déverrouillée redonne la verrouillée | écart max sur B_r = 0 T — **passée** |
| `RUN_G2_CLOTURE` | les deux références du facteur de réluctance | k_r relu du `.mat` contre le projet ANSYS | identiques à 15 décimales — **passée** |
| `RUN_G3_DUMP` | la matrice complète de k_r, les quatre largeurs Newton | idem | **passée** |
| `RUN_G4_POINTS` | les rendements en points et tous les recomptages | l'écart relatif recomposé doit redonner le vecteur `E` sauvegardé par R7 | écart max 7,1×10⁻¹⁵ — **passée** |
| `RUN_G5_COUNT` | le sens du mouvement sous raffinement du bec | écartées + rapprochées + inchangées = 17 | 13 + 4 + 0 = 17 — **passée** |

---

## G1 — la convention sur γ

**La question.** Deux quantités circulaient sous le même intitulé « bracket » :
`ln N + ln|2 sin(d/2)|` au panneau bas de la Table 6, et
`ln N + γ + ln|2 sin(d/2)|` au panneau (a) de la Table 7. La différence est
exactement γ.

**Ce qui tranche.** L'équation (28) est la forme exacte. Le développement

> Σ<sub>n≤N</sub> (1/n) sin²(nd/2) = ½ [ ln N + γ + ln\|2 sin(d/2)\| ] + o(1)

contient γ, qui vient de la somme harmonique et ne peut pas en être retiré.
**γ est retenu partout.**

**Configuration déclarée.** PMSM 15/14 750 W (`machine_bldc`), chaîne
`cogging_mec` + `airgap_magnet`, base P0, solveur linéaire µ_i = 3000,
k_fringe = 0,325, **M_s = 1080 fixe**, N_h = 540, 1080, 2160, 4320,
γ = 0,57721566490153287.

**La valeur.**

```
    SANS gamma : 1.144728475583 -> 3.224170017263
                 croissance = 181.6536922103 %   -> arrondi : +182 %
    AVEC gamma : 1.721944140485 -> 3.801385682165
                 croissance = 120.7612658733 %   -> arrondi : +120.8 %
```

**Ce qui ne bouge pas, et c'est l'essentiel.** γ est une constante additive :

```
    pente SANS gamma : 2.302585092994 par decade
    pente AVEC gamma : 2.302585092994 par decade
    ln 10            : 2.302585092994
    ecart AVEC gamma a ln10 : +0.0000 %
```

**Substitutions au manuscrit.** Croissance **+182 % → +120,8 %**. Équation (30)
réécrite avec γ des deux côtés. En-tête et quatre valeurs du panneau bas de la
Table 6 portés à **1,7219 / 2,4151 / 3,1082 / 3,8014** ; les colonnes B_g1 et
l'écart sont inchangés. Légende de la Table 7 : la phrase disant que les deux
tables n'étaient pas comparables est retirée, puisqu'elles portent désormais la
même forme.

---

## G3 — les quatre largeurs d'encadrement, à pleine précision

**Ce que le manuscrit disait.** Que le crochet de 41,7 points à n_sh = 1 venait
d'une exécution dont la colonne Newton n'était pas reproduite, donc qu'il n'était
pas auditable. **C'est faux.**

Le balayage complet existe et il est dans cette archive :
`outputs/MEC_BLDC/X1_table5b_reconcile_out.txt`, **5 août 2026**, produit par
`RUN_X1_TABLE5B_RECONCILE.m`, données dans `code/MEC_BLDC/X1_table5b.mat`.

**Configuration.** ν = 8, n_sh = 1 à 4, M_s = 900, deux solveurs (linéaire
µ_r = 3000, et Newton exact sur `bh_curve` M350), base P0, numax = ⌊M_s/2⌋,
k_fringe = 0,325 (inopérant, contrôlé au §0 du bloc), φ = 0, référence EF
ν = 8 à 0,0196594150 T.

| n_sh | ν = 8 linéaire | ν = 8 Newton | largeur |
|---|---|---|---|
| 1 | −9,717218 % | +31,982123 % | **41,699342 points** |
| 2 | −16,018999 % | +26,183413 % | **42,202412 points** |
| 3 | −17,987571 % | +24,526287 % | **42,513858 points** |
| 4 | −18,809286 % | +23,829427 % | **42,638713 points** |

Les quatre valeurs publiées — 41,7 / 42,2 / 42,5 / 42,6 — sont exactes, et
l'argument de non-fermeture repose sur le balayage entier, non sur une paire.

---

## G2 et G3 — la matrice du facteur de réluctance

**Le défaut.** Le §5.2 du manuscrit interdit de rapporter une grandeur dont la
référence offre plusieurs déterminations sans qu'aucune soit nommée. Or l'écart
**défavorable** était donné pour une colonne et l'écart **favorable** pour
l'autre, contre deux références différentes.

**Les six valeurs**, toutes de la chaîne du Tableau 13 (`A1_table7.mat`,
M_s = 540, N_p = 61, base P0, solveur linéaire) et des `.tab` relus à
l'exécution :

| | **k_r** | / 1,086928844128470 | / 1,066446589284564 |
|---|---|---|---|
| maillage n_sh = 1 | 1,070690086039355 | **−1,4940 %** | **+0,3979 %** |
| maillage n_sh = 2 | 1,066652526241467 | **−1,8655 %** | **+0,0193 %** |
| réseau à une dent | 1,052198648049849 | **−3,1953 %** | **−1,3360 %** |

Le point écarté est le **premier** des quatorze : 531,1580 A contre 704,5 à
791,7 A pour les treize autres.

**Deux conséquences.** Le « +0,03 % » publié ne se régénère pas : à pleine
précision c'est **+0,02 %**. Et la référence n'offre pas deux déterminations :
elle en offre **une**, 1,086929, lue dans `Output Variables Table 2.tab` ; la
seconde est **reconstruite** par `RUN_KRKL_MESH.m:19` en retirant un point. La
note du manuscrit le dit désormais — c'est un grief plus fort, pas plus faible.

---

## G4 — les rendements en points, et les recomptages

Le §5.3 impose que les grandeurs qui **sont** des pourcentages soient comparées
par une différence en points. Deux lignes du Tableau 20 y contrevenaient.

| ligne | modèle | référence | relatif | **en points** |
|---|---|---|---|---|
| 25, rendement en charge | 86,3340330680 | 86,5147786624 | −0,208919 % | **−0,180746 p.p.** |
| 29, rendement de crête | 85,4455763512 | 89,4319268475 | −4,457413 % | **−3,986350 p.p.** |

| | convention actuelle | lignes 25 et 29 en points |
|---|---|---|
| dans les 5 | **21 / 29** | **21 / 29** |
| hors les trois non validées | **20 / 26** | **20 / 26** |
| à moins d'un point du seuil | **3** | **2** |

Les deux comptes annoncés ne bougent pas. **Le troisième bouge** : la ligne 29
passe de 0,54 point du seuil à 1,01 point et en sort. La légende de la Fig. 10
passe de « Three quantities … all three » à « Two quantities … both ».

**Et un défaut trouvé au passage.** Le **−3,98 p.p.** du Tableau 16 avait été
formé sur une exécution **périmée** : sa note cite 85,4476, valeur du run
préfixe, alors que le run final imprime 85,4456. À pleine précision,
85,4455763512 − 89,4319268475 = **−3,986350** → **−3,99 p.p.**

---

## G5 — le sens du mouvement sous raffinement du bec

La conclusion affirmait que raffiner le bec écarte **douze** des dix-sept
grandeurs de la Table 13. Recompté à pleine précision depuis `A1_table7.mat` :

```
  ECARTE    : 13
  rapproche : 4
  inchange  : 0
  GARDE : 13 + 4 + 0 = 17 (doit valoir 17)   GARDE PASSEE
```

Les quatre qui se rapprochent sont B_t rms, la FMM moyenne par entrefer, la
mutuelle et k_l. **Douze → treize.**

---

## Ce que le bloc G ne tranche pas

- **La bande d'incertitude sur les grandeurs de champ.** Elle exigerait de
  re-résoudre le projet ANSYS à deux densités de maillage. Non fait ; voir
  `docs/OPEN_POINTS.md`.
- **Le DOI.** Aucun n'est frappé, et aucun n'est inventé. Le manuscrit ne porte
  plus de gabarit non résolvable.
