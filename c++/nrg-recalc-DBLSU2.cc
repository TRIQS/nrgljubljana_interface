// *** WARNING!!! Modify nrg-recalc-DBLSU2.cc.m4, not nrg-recalc-DBLSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Dec 2009
// This file pertains to (I1,I2) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers

namespace DBLSU2 {
#include "dblsu2/dblsu2-2ch-def.dat"
}

// Recalculate matrix elements of a doublet tenzor operator
void SymmetryDBLSU2::recalc_doublet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) {
  LOOP(diag, is1) {
    Invar I1   = INVAR(is1);
    Ispin ii11 = I1.get("II1");
    Ispin ii21 = I1.get("II2");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21);
    {
      nrglog('f',
             "RECALC(fn="
                << "dblsu2/dblsu2-2ch-doubletm0.dat"
                << ", len=" << DBLSU2::LENGTH_D1_2CH << ", Iop=" << Invar(2, 0) << ")");
      if (diag.count(Ip)) {
        struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doubletm0.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D1_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, recalc_table, DBLSU2::LENGTH_D1_2CH, Invar(2, 0));
      }
    };

    Ip = Invar(ii11 + 1, ii21);
    {
      nrglog('f',
             "RECALC(fn="
                << "dblsu2/dblsu2-2ch-doubletp0.dat"
                << ", len=" << DBLSU2::LENGTH_D1_2CH << ", Iop=" << Invar(2, 0) << ")");
      if (diag.count(Ip)) {
        struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doubletp0.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D1_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, recalc_table, DBLSU2::LENGTH_D1_2CH, Invar(2, 0));
      }
    };

    Ip = Invar(ii11, ii21 - 1);
    {
      nrglog('f',
             "RECALC(fn="
                << "dblsu2/dblsu2-2ch-doublet0m.dat"
                << ", len=" << DBLSU2::LENGTH_D2_2CH << ", Iop=" << Invar(0, 2) << ")");
      if (diag.count(Ip)) {
        struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doublet0m.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D2_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, recalc_table, DBLSU2::LENGTH_D2_2CH, Invar(0, 2));
      }
    };

    Ip = Invar(ii11, ii21 + 1);
    {
      nrglog('f',
             "RECALC(fn="
                << "dblsu2/dblsu2-2ch-doublet0p.dat"
                << ", len=" << DBLSU2::LENGTH_D2_2CH << ", Iop=" << Invar(0, 2) << ")");
      if (diag.count(Ip)) {
        struct Recalc recalc_table[] = {
#include "dblsu2/dblsu2-2ch-doublet0p.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_D2_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, recalc_table, DBLSU2::LENGTH_D2_2CH, Invar(0, 2));
      }
    };
  }
}

// Override the recalc_f definition: we need to track the type (1 or 2) of
// the f-matrices.

// Driver routine for recalc_f()
void SymmetryDBLSU2::recalc_irreduc(const DiagInfo &diag) {
  // Convention: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1;

    Ispin ii1p = Ip.get("II1");
    Ispin ii2p = Ip.get("II2");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = getnn();

    // RECALC_F_TAB_... (filename, channel_number, matrix_number, array_length)

    // type 1: [f^dag_UP, f_DO]
    // type 2: [f^dag_DO, f_UP]

    I1 = Invar(ii1p + 1, ii2p);
    nrglog('f', "isoup-a I1=" << I1);
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isoup-a.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isoup-a.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][1], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };

    I1 = Invar(ii1p, ii2p + 1);
    nrglog('f', "isoup-b I1=" << I1);
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isoup-b.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isoup-b.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][1], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };

    I1 = Invar(ii1p - 1, ii2p);
    nrglog('f', "isodown-a I1=" << I1);
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isodown-a.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isodown-a.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][1], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };

    I1 = Invar(ii1p, ii2p - 1);
    nrglog('f', "isodown-b I1=" << I1);
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type1-isodown-b.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };
    {
      if (diag.count(I1)) {
        struct Recalc_f recalc_table[] = {
#include "dblsu2/dblsu2-2ch-type2-isodown-b.dat"
        };
        BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLSU2::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][1], Ip, I1, recalc_table, DBLSU2::LENGTH_I_2CH);
      }
    };
  }
}

#undef SPINX
#define SPINX(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#ifdef NRG_COMPLEX
#undef SPINY
#define SPINY(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)
#undef Complex
#define Complex(x, y) cmpl(x, y)
#endif // NRG_COMPLEX

void SymmetryDBLSU2::recalc_global(DiagInfo &diag, string name, MatrixElements &cnew) {
  if (name == "SZtot") {
    LOOP(diag, is1) {
      Invar I1          = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    } // LOOP
  }

#ifdef NRG_COMPLEX
  if (name == "SYtot") {
    LOOP(diag, is1) {
      Invar I1          = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spiny.dat"
          break;
        default: my_assert_not_reached();
      }
    } // LOOP
  }
#endif // NRG_COMPLEX

  if (name == "SXtot") {
    LOOP(diag, is1) {
      Invar I1          = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "dblsu2/dblsu2-2ch-spinx.dat"
          break;
        default: my_assert_not_reached();
      }
    } // LOOP
  }
}
