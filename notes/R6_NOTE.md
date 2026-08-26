# R6 — note de livraison

> Sorties : `MEC_BLDC/R6_satmap_out.txt` (blocs R6b + R6c), `MEC_BLDC/R6a_grid_out.txt`
> Scripts : `RUN_R6A_GRID.m`, `RUN_R6B_FEA.m`, `RUN_R6C_CORR.m`
> Exécuté le 9 août 2026. PMSM 15/14, $N_{surf} = 720$, $N_\theta = 61$, base P0,
> Newton sur B(H) M350-50A, $k_f = 0{,}75$ (rôle induit) **et** $k_f = 0$ (diagnostic).
> Aucun cache : `mec_map.mat` datait du 28/07 alors que `machine_bldc.m` et
> `airgap_magnet.m` ont changé le 04/08 — le cache de `BLDC_MEC_COMPLET` n'est
> indexé que sur le nom de fichier. Tout est reconstruit.
>
> **Vérifié après coup, et rassurant** : le cache relu donne
> 104,318239 / 82,244042 / 62,309367 / 12,057481 mH à 0 / 5 / 10 / 25 A,
> **identiques au rebuild de R6a à la même grille**. Les modifications du 04/08
> n'ont donc pas touché la carte. La précaution restait justifiée — rien dans
> le code ne l'aurait signalé — mais elle n'a rien changé au résultat.

## Les gardes, d'abord

| garde | résultat |
|---|---|
| deux routes EF indépendantes (flux totalisé / tensions de nœud) | ✅ écart moyen **3,78 %**, max 9,37 % |
| la correction cale-t-elle $L(0)$ sur sa cible ? | ✅ écart **9,7 × 10⁻¹⁷** |
| la soustraction est-elle constante en courant ? | ✅ exacte, par construction |
| la chaîne reconstruite reproduit-elle la Fig. 12 publiée ? | ✅ **+0,02 %** à 0 A, **−0,35 %** à 25 A |

L'extraction côté référence n'est donc pas en cause, et la chaîne côté réseau
est bien celle du manuscrit.

## Ce que R6 demandait, et ce qui est mesurable

**2, 5 et 10 A : mesurés. 15 et 20 A : impossibles, et rien n'est extrapolé.**

Le projet EF ne contient **aucun balayage en courant** — `magnetostique(Armature-Field)`
n'a qu'un point, $i_a = 1$ A. La seule source de $\lambda(i)$ est la montée de
courant du transitoire en charge, à rotor quasi immobile : elle plafonne à
**12,269 A**. Au-delà, le rotor accélère de façon monotone et le terme
$\omega\,\partial\lambda/\partial\theta$ atteint **250 %** de la variation de
flux observée — ce n'est plus une correction, c'est le signal.

## Le résultat

À la position rotor **de la fenêtre EF elle-même** ($\theta_e \approx -0{,}14°$) :

| $i$ (A) | EF route 1 | EF route 2 | EF moy. | réseau à $\theta_W$ | écart |
|---|---|---|---|---|---|
| 2 | 64,10 | 65,33 | **64,71** | 89,70 | **+39 %** |
| 5 | 39,56 | 40,81 | **40,19** | 48,49 | **+21 %** |
| 10 | 13,61 | 14,81 | **14,21** | 33,98 | **+139 %** |

**Le réseau est au-dessus de la référence à tous les courants mesurés.**

### La spécification se trompe de signe

La v4 pose que « la carte **sur-sature** à courant intermédiaire ». Les données
disent l'inverse : $L_{MEC} > L_{EF}$ partout. Le réseau **sous-sature**. Le
sens de l'effet est cohérent avec le déficit observé — `drive_mec.m:148` donne
$\mathrm{d}i/\mathrm{d}t = (V - 2Ri - e)/L_{eq}$, donc un $L$ trop grand
ralentit le transitoire et abaisse le courant, ce qui est bien ce qu'on
constate (5,59 A contre 9,14 A à 872 tr/min).

### Le piège nommé par la spécification, chiffré

Comparer la **sécante** de l'EF à l'**incrémentale** du réseau donnerait
+43 / +40 / +52 % au lieu de +55 / +105 / +348 % (colonnes moyennées). À
12,27 A la sécante EF vaut 36,19 mH et l'incrémentale 11,69 mH — **facteur
3,1**. Les deux sont affichées côte à côte dans la sortie.

## Deux hypothèses testées, deux hypothèses réfutées

### La grille de courant : réfutée

**Le signe est inversé.** À 2, 5 et 10 A la grille de production *sous*-estime
$L$ de 1,75 / 0,01 / 1,16 mH. La raison est que `gradient` sur grille non
uniforme n'est pas une différence centrée mais **la corde brute entre nœuds
encadrants** — premier ordre — dont le milieu est décalé de $28/81 = 0{,}3457$ A
à tout nœud intérieur de la loi quadratique. Le terme $L'\cdot 0{,}3457$,
négatif, domine celui de Jensen. Et $L(i)$ n'est **pas convexe** : elle est
concave sur $[0\,;\,3{,}5]$ A, où tombe justement le point 2 A.

**L'amplitude est dérisoire.** Passer de la grille de production à la grille
fine ajoute en moyenne +0,98 mH sur 0–10 A ; reporté sur la constante de temps,
cela fait passer le réseau de 5,59 A à **5,52 A**, quand il faut **+3,55 A**.
Soit **2 % de l'écart, et dans le mauvais sens**.

**Le défaut est anti-corrélé au mécanisme invoqué.** À 297 tr/min le réseau
travaille à 17,76 A, là où la corde est la plus large, et n'est qu'à −8,7 % ;
à 1446 tr/min il travaille à 0,82 A, sous le deuxième nœud de grille, là où
$L$ est indépendante de la grille à sept chiffres — et le déficit y vaut
−25,5 %, trois fois pire.

**Mais la grille n'est pas convergée en haut de plage** : dispersion **26,6 %
à 15 A** et **41,6 % à 25 A** sur $n_a$ = 19 → 43. C'est un défaut réel, et la
section suivante en donne la cause.

### La calibration effaçant la pré-saturation par les aimants : réfutée

`mec_map.m:215` cale $L(0)$ sur `inductance_mec(...)`, réseau **linéaire et
aimants éteints**. J'ai supposé que cette calibration détruisait la
pré-saturation de la denture par les aimants. **Test direct**, carte rejouée à
$k_f = 0$ où il n'y a plus de double comptage à corriger :

| | mH |
|---|---|
| $L$ brute à $i=0$, aimants ON, non linéaire | 103,7563 |
| cible du code : $2L_d$ linéaire, aimants OFF | 104,3182 |
| EF magnétostatique, aimants OFF | 104,6831 |
| **EF transitoire, aimants ON, $i \to 0$** | **72,2471** |

**Les aimants coûtent au réseau −0,56 mH (−0,5 %) et à l'EF −32,44 mH
(−31,0 %).** La calibration ne cache rien : le réseau ne produit pas l'effet.
Mon hypothèse tombe.

## Ce qui survit — et c'est le résultat

### L'effondrement publié n'est produit par aucun réseau

| $i$ (A) | publiée ($k_f{=}0{,}75$ corrigée) | **brute** du même réseau | brute sans pont tangentiel |
|---|---|---|---|
| 0 | 104,3 | 133,7 | 103,8 |
| 10 | 63,6 | 93,0 | 90,9 |
| 25 | **9,5** | **38,9** | 58,2 |
| **facteur** | **11,0** | **3,4** | **1,8** |

La courbe publiée de la Fig. 12 tombe d'un facteur 11 ; la physique brute du
même réseau, d'un facteur 3,4 ; le réseau sans pont tangentiel, d'un facteur
1,8. **Le facteur publié vient de ce qu'on retranche une constante de
29,36 mH à une grandeur qui, elle, ne descend qu'à 38,9 mH.**

### Et la prémisse de cette soustraction est fausse

`mec_map.m:203-213` justifie la soustraction ainsi : *« L'excès est un chemin
dans l'AIR : il ne sature pas, la correction est donc SOUSTRACTIVE. »* Si
c'était vrai, l'écart entre le réseau avec pont tangentiel et le réseau sans
serait constant. Mesure : **+29,921 mH à 0 A, −19,303 mH à 25 A**. Il **change
de signe**, sur une amplitude de 49,2 mH.

Retrancher une constante calibrée à courant nul est donc exact à courant nul,
et de plus en plus faux ensuite.

### Cela explique la non-convergence en haut de plage

À 25 A la valeur publiée est la **différence de 38,9 et 29,4 mH**. Tout écart
relatif sur la brute s'y trouve multiplié par **4,1**. La dispersion de 41,6 %
que R6a mesure à 25 A, contre 0,3 % à 5 A, est cohérente avec ce facteur.
*Statut : explication cohérente avec les nombres, non mesurée indépendamment —
R6a n'a pas archivé $\alpha$ pour chaque grille.*

### La Fig. 12 publie une moyenne d'une grandeur qui varie d'un facteur 3

| $i$ (A) | moyenne sur $\theta$ | à $\theta_W$ | min sur le tour | max sur le tour |
|---|---|---|---|---|
| 0 | 104,32 | 103,39 | 103,28 | 106,00 |
| 5 | 82,34 | 48,49 | **35,97** | **106,02** |
| 10 | 63,64 | 33,98 | **31,56** | **105,93** |

**Et le modèle d'entraînement n'utilise pas cette moyenne** : `drive_mec.m:117`
interpole `La`/`Lb` à la position instantanée. La courbe publiée n'est donc pas
la grandeur que le modèle intègre.

## Ce sur quoi je ne conclus pas

**La référence se contredit elle-même à courant nul.** Son étude
magnétostatique donne 104,68 mH (aimants éteints), son transitoire 72,25 mH
(aimants présents), soit **−31,0 %**, à 116 At où aucune saturation d'induit
n'est possible. Le réseau, lui, s'accorde avec la magnétostatique à −0,35 %.

**Le transitoire EF est écranté par des courants induits.** `Plot 1_loss.tab` :
les pertes solides valent **52,3 % des pertes Joule au premier échantillon**,
16,7 % au deuxième, 5,1 % au quatrième — le poids relatif est **maximal
précisément à bas courant**, là où l'on voudrait mesurer. Une part du $L$ EF
bas est donc un effet dynamique, absent d'une carte quasi statique.

**Conséquence sur la portée du résultat : les écarts +39 / +21 / +139 % sont
une borne SUPÉRIEURE de l'erreur quasi statique du modèle**, pas une mesure de
celle-ci. Les séparer exigerait un balayage magnétostatique en courant à rotor
figé — la même dépendance ANSYS que R3.

## Écart relevé en passant

**$R_{ph}$.** Le circuit EF utilise **10,220000 Ω exactement** — vérifié à
10⁻¹³ près sur les échantillons $|i| > 1$ A via $(I_{v} - I_{e})/i$ — contre
`M.Rph = 10.222304` dans `machine_bldc.m:76` et 10,222 Ω au manuscrit.
Sans effet ici (0,023 %), mais c'est la valeur à utiliser pour dépouiller l'EF.

## Affirmation publiée à corriger

Le manuscrit écrit que placer la perméance d'air en série avec la réluctance du
fer *« restores the physical collapse »*. La saturation de la fuite d'encoche
mesurée dans la carte ne descend qu'à $\lambda_{ef}/\lambda_{tot} = 0{,}999$ à
0 A, **0,991 à 5 A et 0,810 à 23 A** : elle contribue peu. L'effondrement
publié est porté par la soustraction constante, pas par ce mécanisme.

**Signalé pour correction. La chaîne n'a été ajustée dans aucun sens.**

## Un diagnostic réfuté est déclaré

La v4 demande : *« Si les deux cartes coïncident, l'hypothèse tombe et l'erreur
est ailleurs ; il faudra alors le dire dans le manuscrit. »* Elles ne
coïncident pas — mais **ni la sur-saturation annoncée, ni la grille, ni la
calibration des aimants n'en sont la cause**. Ce qui reste, mesuré, est que la
courbe $L(i)$ publiée doit une part majeure de sa forme à une soustraction
constante dont la prémisse est réfutée.
