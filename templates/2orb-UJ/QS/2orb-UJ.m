def2ch[nrimp=2];
impuritybasis = {d[1],d[2]};

Heps1 = eps1 number[d[1]];
Hint1 = U1 hubbard[d[1]];
Hhyb1 = gammaPolCh[1] hop[f[0], d[1]]; (* f[0] is the 0-th site of the Wilson chain for channel 1 *)

Heps2 = eps2 number[d[2]];
Hint2 = U2 hubbard[d[2]];
Hhyb2 = gammaPolCh[2] hop[f[1], d[2]]; (* f[1] is the 0-th site of the Wilson chain for channel 2 *)

H12 = U12 nc[number[d[1]], number[d[2]]]+ J12 spinspin[d[1], d[2]];

Heps = Heps1 + Heps2;
Hint = Hint1 + Hint2 + H12;
Hhyb = Hhyb1 + Hhyb2;
H = H0 + Heps + Hint + Hhyb;

selfopd1 = ( Chop @ Expand @ komutator[Hint, d[#1, 1, #2]] )&;
selfopd2 = ( Chop @ Expand @ komutator[Hint, d[#1, 2, #2]] )&;

SigmaHartree1 = Expand @ antikomutator[ selfopd1[CR, #1], d[AN, 1, #1] ] &;
SigmaHartree2 = Expand @ antikomutator[ selfopd2[CR, #1], d[AN, 2, #1] ] &;

SigmaHartreeAvg1 := Expand @ (SigmaHartree1[UP] + SigmaHartree1[DO]) / 2;
SigmaHartreeAvg2 := Expand @ (SigmaHartree2[UP] + SigmaHartree2[DO]) / 2;

Print["SigmaHartree1[UP]=", SigmaHartree1[UP] ];
Print["SigmaHartree1[DO]=", SigmaHartree1[DO] ];
Print["SigmaHartree1=", SigmaHartreeAvg1 ];

Print["SigmaHartree2[UP]=", SigmaHartree2[UP] ];
Print["SigmaHartree2[DO]=", SigmaHartree2[DO] ];
Print["SigmaHartree2=", SigmaHartreeAvg2 ];