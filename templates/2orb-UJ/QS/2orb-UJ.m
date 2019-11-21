def2ch[nrimp=2];

Heps1 = eps1 number[d[]];
Hint1 = U1 hubbard[d[]];
Hhyb1 = gammaPol1 hop[f[0], d[]];
H1 = Heps1 + Hint1;

Heps2 = eps2 number[a[]];
Hint2 = U2 hubbard[a[]];
Hhyb2 = gammaPol2 hop[f[1], a[]];
H2 = Heps2 + Hint2;

H12 = U12 nc[number[d[]], number[a[]]]+ J12 spinspin[d[], a[]];

Heps = Heps1 + Heps2;
Hint = Hint1 + Hint2 + H12;
Hhyb = Hhyb1 + Hhyb2;
H = H0 + Heps + Hint + Hhyb;

params = Join[params, {
  gammaPol1 -> Sqrt[(1/Pi) thetaCh[1] (extraGamma1) gammaA],
  gammaPol2 -> Sqrt[(1/Pi) thetaCh[2] (extraGamma2) gammaA]
}];

selfopd = ( Chop @ Expand @ komutator[Hint, d[#1, #2]] )&;
selfopa = ( Chop @ Expand @ komutator[Hint, a[#1, #2]] )&;
