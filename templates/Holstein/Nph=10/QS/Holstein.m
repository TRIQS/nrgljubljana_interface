def1ch[nrimp=1];
Heps = eps1 number[d[]];
HU = U1 hubbard[d[]];
nph = ToExpression @ optionvalue["Nph"];
Hph = omega phononnumber[nph] + g1 nc[number[d[]]-n1, phononplus[nph] + phononminus[nph]];
Hint = HU + Hph;
MAKEPHONON = 1; (* One phonon mode *)
Himp = Heps + Hint;
Hhyb = gammaPolCh[1] hop[f[0], d[]];
H = H0 + Himp + Hhyb;
selfopd = ( Chop @ Expand @ komutator[Hint, d[#1, #2]] )&;

SigmaHartree = Expand @ antikomutator[ selfopd[CR, #1], d[AN, #1] ] &;
SigmaHartreeAvg := Expand @ (SigmaHartree[UP] + SigmaHartree[DO]) / 2;

Print["SigmaHartree[UP]=", SigmaHartree[UP] ];
Print["SigmaHartree[DO]=", SigmaHartree[DO] ];
Print["SigmaHartree=", SigmaHartreeAvg ];