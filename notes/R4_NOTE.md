# R4 — note de livraison

> Sortie : `MEC_BLDC/R4_unlock_out.txt`
> Script : `MEC_BLDC/RUN_R4_UNLOCK.m`
> Exécuté le 9 août 2026. PMSM 15/14, base P0, solveur linéaire μ_i = 3000,
> k_fringe = 0,3250.

## La prédiction du §3.5 est confirmée

**Le verrou, d'abord.** Sur un facteur **16** en $M_s$, le crochet
$\ln N + \ln|2\sin(d/2)|$ vaut :

| $M_s$ | 540 | 1080 | 2160 | 4320 | 8640 |
|---|---|---|---|---|---|
| crochet | 1,144724 | 1,144728 | 1,144730 | 1,144730 | 1,144730 |

soit **+0,0005 %**. Et $\ln\pi = 1{,}144730$ : **le crochet vaut $\ln\pi$ à
six décimales**. La stationnarité n'est pas approchée, elle est exacte, et
elle est celle que le mécanisme prédit.

**Le déverrouillage ensuite.** À pavage **fixe** $M_s = 1080$, en libérant
$N_h$ :

| $N_h$ | $N/$verrou | crochet | $B_{g1}$ (T) | écart |
|---|---|---|---|---|
| 540 | 1 | 1,144728 | 1,078648 | — |
| 1080 | 2 | 1,837876 | 1,079783 | +0,1052 % |
| 2160 | 4 | 2,531023 | 1,080847 | +0,2038 % |
| 4320 | 8 | 3,224170 | 1,081745 | +0,2871 % |

**Le crochet reprend sa croissance : +182 % sur ce balayage.** Et sa pente
vaut **2,302585 par décade**, contre $\ln 10 = 2{,}302585$ prédit — **écart
−0,0000 %**. La prédiction est vérifiée à la précision machine.

**$B_{g1}$ dérive de façon monotone** et ne se stabilise pas : +0,29 % sur un
facteur 8 en $N$. **La dérive revient**, exactement comme le §3.5 l'annonce.

## La GARDE passe

À $N_h = M_s/2$ exactement, la chaîne déverrouillée redonne la chaîne
verrouillée :

```
Bg1 verrouille   : 1.078648294812032 T
Bg1 deverrouille : 1.078648294812032 T
ecart            : 0.000e+00 T
ecart max sur Br : 0.000e+00 T
```

**Identité au bit près.** Le portage n'a donc changé **que** la troncature.

## Ce qui n'est pas publiable

**Le point $N_h = 8640$ est à écarter.** La matrice devient singulière
(`RCOND = NaN`) et $B_{g1}$ sort NaN. Le mécanisme est connu et doit être
énoncé plutôt que masqué : la projection $W(n) = (2/n\pi)\sin(n d/2)\cos(n\theta)$
**s'annule aux multiples de $M_s$**. À $N_h = 8\,M_s$, huit familles
d'harmoniques ne contribuent rien, et l'opérateur perd progressivement son
rang. Le balayage exploitable s'arrête à $N_h = 4\,M_s$.

C'est une limite de la voie retenue, pas du résultat : la croissance du
crochet et la dérive de $B_{g1}$ sont établies sur un facteur 8, ce qui suffit
à falsifier l'immunité.

## Un point sur lequel je ne conclus pas

**La dérive revient, mais elle est plus faible que sur la MAS** : +0,29 % sur
$B_{g1}$ ici, contre +5,99 % sur $X_m$ pour un balayage comparable côté MAS.

**Ce rapport n'explique pas cette différence d'amplitude**, et n'en propose
aucune. En particulier, il ne l'attribue **pas** à l'épaisseur d'anneau : cette
explication est réfutée par l'Article I §3.5 — la queue s'installe dès
$n \approx 7$ sur le PMSM contre $n \approx 337$ sur la MAS, donc **plus tôt**,
et si l'épaisseur était le mécanisme cette machine dériverait **davantage**.

Ce qui est établi est **l'existence** de la dérive, non sa magnitude.

## Choix de mise en œuvre, à déclarer

La v4 demandait de **porter la condensation de Schur** sur le bore du PMSM.
Je ne l'ai pas fait, et voici pourquoi : la condensation exigerait de
reconstruire le réseau de réluctances sur $N_s = 15$ nœuds au lieu de
$N_{\text{surf}}$, donc de **réécrire la chaîne** — ce que la règle du dossier
proscrit et ce qui a déjà coûté cinq itérations sur l'ondulation.

La contrainte de `cogging_mec` porte sur $n_{\max}$ **inférieur** à
$N_{\text{surf}}/2$ (rang déficient). **Rien n'interdit de monter au-dessus.**
Le troisième argument, jusqu'ici ignoré, est désormais honoré s'il atteint ce
seuil, avec avertissement sinon. Le déverrouillage est obtenu sans toucher au
reste, et **la garde le prouve** : identité au bit près au point de recouvrement.

Si la condensation reste souhaitée pour d'autres raisons — un panneau PMSM à
la Table 5, par exemple — elle reste à faire ; mais elle n'est pas nécessaire
à R4.

## Critère d'acceptation v4

*« Montrer que le crochet reprend sa croissance, avec la pente prédite »* —
**satisfait** : +182 %, pente $\ln 10$ à −0,0000 %.

*Garde exécutée* — **oui, et elle passe** : 0,000e+00 T d'écart.

**Bénéfice attendu, obtenu** : le meilleur raisonnement du manuscrit devient
une démonstration. Le panneau LaTeX pour la Table 8 est dans la sortie.
