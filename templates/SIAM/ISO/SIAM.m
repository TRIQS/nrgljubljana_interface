def1ch[nrimp=1];
Heps = (-U1 / 2) number[ d[] ];
Hint = U1 hubbard[ d[] ];
Himp = Heps + Hint;
Hhyb = gammaPolCh[1] hop[f[0], d[]];
H = H0 + Himp + Hhyb;
isoselfopd = ( Chop @ Expand @ komutator[ (2 U1) pow[ isospinz[ d[] ], 2 ], d[#1, #2] ] )&;
