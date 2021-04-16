tt={};

tt = Join[tt, mtSingletOp["n_d1",     number[d[1]] ]];
tt = Join[tt, mtSingletOp["n_d2",     number[d[2]] ]];
tt = Join[tt, mtSingletOp["hop12",    hop[d[1], d[2]] ]];
tt = Join[tt, mtSingletOp["n_d1^2",   pow[number[d[1]],2] ]];
tt = Join[tt, mtSingletOp["n_d2^2",   pow[number[d[2]],2] ]];
tt = Join[tt, mtSingletOp["S_d1S_d2", spinspin[d[1], d[2]] ]];
tt = Join[tt, mtSingletOp["n_d1n_d2", nc[number[d[1]], number[d[2]]] ]];

tt = Join[tt, mtDoubletOp["A_d1", d[1] ]];
tt = Join[tt, mtDoubletOp["A_d2", d[2] ]];
tt = Join[tt, mtDoubletOp["self_d1", selfopd1 ]];
tt = Join[tt, mtDoubletOp["self_d2", selfopd2 ]];

tt = Join[tt, mtTripletOp["sigma_d1", d[1] ]];
tt = Join[tt, mtTripletOp["sigma_d2", d[2] ]];

tt
