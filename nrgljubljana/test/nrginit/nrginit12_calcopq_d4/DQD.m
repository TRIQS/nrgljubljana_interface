def1ch[nrimp=2];
Heps = eps1 number[d] + eps2 number[a];
Hpot = Hint = U1 hubbard[d[]] + U2 hubbard[a[]];
Himp = Heps + Hint + t12 hop[d[],a[]];
Hhyb = gammaPolCh[1] hop[f[0], d[]];
H = H0 + Himp + Hhyb;
selfopd = ( Chop @ Expand @ komutator[Hint /. params, d[#1, #2]] )&;
selfopcd = ( Chop @ Expand @ komutator[Hint /. params, ((-1)^#2 d[1-#1, 1-#2]) ] )&;
selfopa = ( Chop @ Expand @ komutator[Hint /. params, a[#1, #2]] )&;
selfopca = ( Chop @ Expand @ komutator[Hint /. params, ((-1)^#2 a[1-#1, 1-#2]) ] )&;
Print[selfopd[CR,UP]];
Print[selfopcd[CR,UP]];
Print[selfopa[CR,UP]];
Print[selfopca[CR,UP]];
