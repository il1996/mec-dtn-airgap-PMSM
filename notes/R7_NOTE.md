# R7 — note de livraison

> Sortie : `MEC_BLDC/R7_scorecard_out.txt`
> Script : `MEC_BLDC/RUN_R7_SCORECARD.m` (exécute `BLDC_MEC_COMPLET`, 155 s)
> Exécuté le 9 août 2026. PMSM 15/14, cartes `mec_map.mat` / `mec_map_field.mat`,
> vérifiées identiques à un rebuild en R6. Seuil |écart| < 5 %, strict, comme le
> code. Aucun chiffre transcrit : la liste est celle que le programme construit.

## Le 23/29 est exact

**GARDE PASSÉE** : le recomptage indépendant donne **23 sur 29**, identique au
compteur `nok` que le programme calcule lui-même (`BLDC_MEC_COMPLET.m:748-753`).

Le manuscrit annonce *« twenty-three fall within 5 % »* — phrase relue du `.tex`,
non transcrite. **L'annonce est vérifiée.** Contrairement au précédent du
« five of the seventeen », qui s'était révélé faux au recomptage, celui-ci tient.

## Pourquoi le rapporteur en compte au moins sept

La figure trace des barres ; trois grandeurs sont à moins d'un point du seuil,
et **les trois sont à l'intérieur** :

| grandeur | écart |
|---|---|
| Couple à l'arrêt (balayage) | **+4,9520 %** |
| Rendement max (balayage) | −4,4552 % |
| Mutuelle $M$ | +4,1538 % |

Le couple à l'arrêt à **+4,95 %** est à cinq centièmes de point du seuil : sur
un graphique à barres il est indiscernable de la ligne. Le comptage visuel du
rapporteur est donc explicable sans que le manuscrit soit en faute — mais il
signale que **la figure ne permet pas de vérifier son propre titre**.

## Hors grandeurs déclarées non validées : 21 / 26

Trois exclusions, chacune justifiée par une phrase **du manuscrit** :

| # | grandeur | écart | motif |
|---|---|---|---|
| 20 | Pertes aimant à vide | +2,4450 % | §6.6 — *« The reference cannot resolve the quantity on which it is being compared. »* |
| 24 | Pertes aimant en charge | −49,9506 % | §6.6 + Table 19 — écart −1,6 W sous le plancher de résolution de 52,7 W |
| 27 | Couple à l'arrêt (balayage) | +4,9520 % | §5 — *« cannot be claimed as a validation … excluded from the comparison »* |

**Compte hors grandeurs non validées : 21 sur 26**, soit cinq dépassements au
lieu de six. La v4 demande que ce soit **ce compte qui figure au résumé**.

Les cinq dépassements qui subsistent :

| grandeur | écart |
|---|---|
| Puissance max (balayage) | **−25,4270 %** |
| Courant de bus moyen | −13,8796 % |
| Pertes fer à vide | +8,1401 % |
| Bt tangentiel rms | +7,4505 % |
| Vitesse établie en charge | −6,1797 % |

L'exclusion retire le plus gros dépassement du lot (−49,95 % sur les pertes
aimant en charge). **Le compte s'améliore à peine — 79,3 % contre 80,8 % — et
c'est en soi rassurant** : la synthèse ne repose pas sur ce que le manuscrit
écarte.

## Un défaut découvert en recomptant, et il est sérieux

L'entrée #20 vaut **+2,4450 %** : `Ppm` = 0,341619 W contre 0,333466 W de
référence. C'est **exactement la valeur publiée puis retirée** — 0,342 W,
+2,4 %.

Or le §6.6 de l'Article I écrit :

> *« The no-load value on which the 2.4 % rested is not produced by any chain
> that can be reconstructed. »*

**Cette phrase est fausse.** La chaîne existe, elle est dans le programme
maître, et elle a tourné il y a 155 secondes :

```
BLDC_MEC_COMPLET.m:338-340
  Rnl = cogging_mec(M,1680,0,181,M.muI,kfr,2*pi/Ns);
  Ppm = pm_loss_load(M,Rnl,zeros(1680,3),tnl,zeros(3,721),om_nl*tnl,41);
```

C'est `pm_loss_load` — la routine de la charge, appelée à courants nuls sur un
pas d'encoche — et non `pm_loss`, la chaîne que le §6.5 formalise et que A3 a
gelée. Sur le même point de fonctionnement et contre la même référence, les
deux donnent :

| chaîne | $P_{pm}$ à vide | écart |
|---|---|---|
| `pm_loss_load`, réseau à une dent (programme maître) | **0,341619 W** | **+2,45 %** |
| `pm_loss`, réseau à une dent (A3, gelée) | 0,2036 W | −39,1 % |
| `pm_loss`, maillage polaire (A3, gelée) | 0,2516 W | −24,7 % |

**Ce n'est donc pas un problème de reproductibilité, c'est un problème de
choix de chaîne** — et les deux chaînes diffèrent d'un facteur 1,7 sur la même
grandeur. *Je n'isole pas ici la cause de l'écart entre les deux routines : les
deux exécutions diffèrent aussi par l'échantillonnage (181 contre 241
positions). Ce qui est établi est que la valeur publiée se régénère.*

### Ce qu'il faut corriger, et ce qui vous revient

**À corriger** : la phrase *« is not produced by any chain that can be
reconstructed »*. Elle doit dire que la valeur n'est pas produite par la
**chaîne gelée du §6.5**, ce qui est vrai et vérifié, plutôt que par aucune
chaîne, ce qui est faux.

**Ce que je ne tranche pas** : laquelle des deux routines est *le modèle*.
C'est une décision d'auteur, et elle a des conséquences — si `pm_loss_load`
est le modèle à vide, l'argument du §6.6 sur la disparition de la dichotomie
doit être réexaminé ; si c'est `pm_loss`, alors le programme maître et la
Fig. 10 doivent être changés, pas seulement le texte.

**Incohérence interne à signaler** : la Fig. 10 telle que publiée porte une
barre à +2,4 % pour une grandeur que le §6.6 déclare retirée. Le tableau et le
texte ne disent pas la même chose.

## Ce que R7 ne remet pas en cause

Les vingt-neuf grandeurs sont bien vingt-neuf : `BLDC_MEC_COMPLET` contient
exactement vingt-neuf appels `SC=add(...)`, vérifié par recherche. Le
dénominateur est juste, le seuil appliqué est bien celui du code (`abs(e)<5`,
strict), et le recomptage reproduit le compteur du programme.
