def1ch[nrimp=1];
Heps = eps1 number[d[]] + B1 spinz[d[]];
Hint = U1 hubbard[d[]];
Himp = Heps + Hint;
Hhyb = gammaPolCh[1] hop[f[0], d[]];
H = H0 + Himp + Hhyb;
selfopd = ( Chop @ Expand @ komutator[Hint /. params, d[#1, #2]] )&;
