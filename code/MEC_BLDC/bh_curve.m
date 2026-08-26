function BH = bh_curve(which)
%BH_CURVE  Courbe d'aimantation B(H) et derivees pour le solveur non lineaire.
%
%   BH = bh_curve()         -> M350-50A REELLE (defaut)
%   BH = bh_curve('M800')   -> ancien proxy M800-50A (comparaison seulement)
%
%   *** SOURCE : le PROJET ANSYS LUI-MEME ***
%   La table M350-50A est extraite mot pour mot du materiau
%     'Cogent Power - M350-50A, B-H at 50Hz'  (bloc 'BHCoordinates', 81 points)
%   du fichier C:\Users\hp\Desktop\ANSYS-\BLDC\BLDC.aedt, affecte aux objets
%   'stator' ET 'rotor'. Ce n'est donc plus un proxy : c'est EXACTEMENT la
%   courbe que resout la reference EF.
%
%   Difference majeure avec l'ancien proxy M800-50A : la M350 sature bien plus
%   TOT (1,827 T a 32 kA/m, contre 2,364 T a 270 kA/m) et sa permeabilite de
%   coude est bien PLUS ELEVEE (mu_r max ~ 6700 a 0,9 T contre ~4400) - deux
%   effets qui comptent pour la denture (a8 croit avec mu_r) et pour la
%   saturation en charge.
%
%   Sorties : BH.Bof(H), BH.Hof(B), BH.dHdB(B), BH.mur(B)   (vectorisees)
%   Au-dela du dernier point : extrapolation a la pente mu0 (monotone,
%   physique en saturation profonde) - meme convention qu'ANSYS.

if nargin<1 || isempty(which), which='M350'; end
mu0=4*pi*1e-7;

switch upper(which)
case 'M350'
% --- Cogent Power M350-50A, B-H at 50Hz  (ANSYS BLDC.aedt, 81 points) -----
T=[      0        0
      9.03568  0.025
     17.6854   0.05
     25.5633   0.075
     32.2833   0.1
     37.5786   0.125
     41.6583   0.15
     44.8505   0.175
     47.4833   0.2
     49.8375   0.225
     52.0042   0.25
     54.0271   0.275
     55.95     0.3
     57.8115   0.325
     59.6292   0.35
     61.4156   0.375
     63.1833   0.4
     64.9438   0.425
     66.7042   0.45
     68.4708   0.475
     70.25     0.5
     72.0482   0.525
     73.8729   0.55
     75.732    0.575
     77.6333   0.6
     79.5841   0.625
     81.5896   0.65
     83.6544   0.675
     85.7833   0.7
     87.9826   0.725
     90.2646   0.75
     92.6435   0.775
     95.1333   0.8
     97.7531   0.825
    100.542    0.85
    103.543    0.875
    106.8      0.9
    110.358    0.925
    114.267    0.95
    118.575    0.975
    123.333    1.0
    128.602    1.025
    134.479    1.05
    141.076    1.075
    148.5      1.1
    156.909    1.125
    166.646    1.15
    178.102    1.175
    191.667    1.2
    207.922    1.225
    228.208    1.25
    254.057    1.275
    287.0      1.3
    329.26     1.325
    385.833    1.35
    462.406    1.375
    564.667    1.4
    698.927    1.425
    874.0      1.45
   1099.32     1.475
   1384.33     1.5
   1737.49     1.525
   2163.38     1.55
   2665.57     1.575
   3247.67     1.6
   3912.76     1.625
   4661.94     1.65
   5495.77     1.675
   6414.83     1.7
   7407.02     1.72474
   8409.48     1.74794
   9346.65     1.76806
  10143.0      1.78356
  10794.8      1.79346
  11585.9      1.79912
  12871.7      1.80246
  15008.0      1.80539
  18240.8      1.80945
  22378.0      1.81465
  27118.3      1.82061
  32160.0      1.82694];
case 'M800'
% --- ancien proxy (conserve pour la comparaison documentee) ---------------
T=[      0      0
     43.02 0.0398
     69.73 0.0835
     88.59 0.1323
    103.35 0.1888
    116.25 0.2583
    130.00 0.3592
    159.56 0.6343
    170.84 0.7347
    186.72 0.8601
    204.00 0.9724
    217.92 1.0460
    232.85 1.1108
    249.52 1.1693
    268.78 1.2232
    291.75 1.2733
    320.11 1.3204
    356.45 1.3650
    405.04 1.4073
    473.42 1.4478
    575.46 1.4866
    737.25 1.5240
   1003.32 1.5600
   1424.27 1.5949
   2011.97 1.6287
   2732.41 1.6615
   3553.37 1.6934
   4465.91 1.7245
   5480.28 1.7548
   6621.08 1.7844
   7926.87 1.8133
   9453.39 1.8417
  11280.56 1.8694
  13523.54 1.8966
  16348.96 1.9232
  19994.67 1.9494
  24783.33 1.9751
  31101.82 2.0004
  39304.40 2.0252
  49546.65 2.0497
  61679.26 2.0737
  75322.75 2.0975
  90041.78 2.1208
 105466.11 2.1439
 121323.48 2.1666
 137426.62 2.1890
 157713.96 2.2166
 182093.34 2.2491
 206364.10 2.2811
 230437.74 2.3125
 254265.30 2.3433
 270000.00 2.3636];
otherwise
    error('bh_curve:mat','Nuance inconnue : %s',which);
end

H=T(:,1); B=T(:,2);
Hend=H(end); Bend=B(end);
slope=diff(H)./diff(B);              % dH/dB par segment

    function Bo=Bof(Hin)
        Hin=Hin(:); s=sign(Hin); s(s==0)=1; a=abs(Hin);
        Bo=zeros(size(a)); in=a<=Hend;
        Bo(in)=interp1(H,B,a(in),'linear');
        Bo(~in)=Bend+(a(~in)-Hend)*mu0;
        Bo=s.*Bo;
    end
    function Ho=Hof(Bin)
        Bin=Bin(:); s=sign(Bin); s(s==0)=1; a=abs(Bin);
        Ho=zeros(size(a)); in=a<=Bend;
        Ho(in)=interp1(B,H,a(in),'linear');
        Ho(~in)=Hend+(a(~in)-Bend)/mu0;
        Ho=s.*Ho;
    end
    function d=dHdB(Bin)
        a=abs(Bin(:)); d=zeros(size(a)); in=a<=Bend;
        idx=discretize(a(in),B); idx(isnan(idx))=1;
        d(in)=slope(min(idx,numel(slope)));
        d(~in)=1/mu0;
        d=max(d,1e-6);
    end
    function m=mur(Bin)
        a=abs(Bin(:)); h=Hof(a); m=ones(size(a));
        nz=a>1e-9; m(nz)=a(nz)./(mu0*max(h(nz),1e-12));
        m(~nz)=B(2)/(mu0*H(2));
    end
BH.Bof=@Bof; BH.Hof=@Hof; BH.dHdB=@dHdB; BH.mur=@mur;
BH.Bsat=Bend; BH.Hsat=Hend; BH.mu0=mu0; BH.table=T; BH.name=upper(which);
end
