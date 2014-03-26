/* dweb.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal aa, ee, gg, bb, dprey, dpred;
} ppar1_;

#define ppar1_1 ppar1_

struct {
    integer np, ns;
    doublereal ax, ay, acoef[4]	/* was [2][2] */, bcoef[2], dx, dy, alph, 
	    beta, fpi, diff[2], cox[2], coy[2];
    integer mx, my, mxns;
} ppar2_;

#define ppar2_1 ppar2_

/* Table of constant values */

static integer c__1 = 1;
static integer c__40 = 40;
static integer c_b58 = 104883;
static integer c__1640 = 1640;

/* ***BEGIN PROLOGUE  DWEB */
/* ***REFER TO  DDASKR */
/* ***DATE WRITTEN   020813   (YYMMDD) */
/* ***REVISION DATE  021217   Added JROOT output value. */

/* ***AUTHORS  A. C. Hindmarsh, P. N. Brown */
/*            Lawrence Livermore National Laboratory */
/*            Livermore, CA 94551, USA */

/* ***DESCRIPTION */

/* ----------------------------------------------------------------------- */
/* Example program for DDASKR. */
/* DAE system derived from ns-species interaction PDE in 2 dimensions. */

/* This is the double precision version. */
/* ----------------------------------------------------------------------- */

/* This program solves a DAE system that arises from a system */
/* of partial differential equations.  The PDE system is a food web */
/* population model, with predator-prey interaction and diffusion on */
/* the unit square in two dimensions.  The dependent variable vector is */

/*         1   2        ns */
/*   c = (c , c , ..., c  ) */

/* and the PDEs are as follows.. */

/*     i               i      i */
/*   dc /dt  =  d(i)*(c    + c   )  +  R (x,y,c)  (i=1,...,ns/2) */
/*                     xx     yy        i */

/*                     i      i */
/*   0       =  d(i)*(c    + c   )  +  R (x,y,c)  (i=(ns/2)+1,...,ns) */
/*                     xx     yy        i */

/* where */
/*                  i          ns         j */
/*   R (x,y,c)  =  c *(b(i) + sum a(i,j)*c ) */
/*    i                       j=1 */

/* The number of species is ns = 2*np, with the first np being prey and */
/* the last np being predators.  The coefficients a(i,j), b(i), d(i) are */

/*   a(i,i) = -a  (all i) */
/*   a(i,j) = -g  (i .le. np, j .gt. np) */
/*   a(i,j) =  e  (i .gt. np, j .le. np) */
/*   b(i) =  b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .le. np) */
/*   b(i) = -b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .gt. np) */
/*   d(i) = dprey  (i .le. np) */
/*   d(i) = dpred  (i .gt. np) */

/* The various scalar parameters are set in subroutine SETPAR. */

/* The boundary conditions are.. normal derivative = 0. */
/* A polynomial in x and y is used to set the initial conditions. */

/* The PDEs are discretized by central differencing on a MX by MY mesh. */

/* The root function is R(T,Y,Y') = average(c1) - 20. */

/* The DAE system is solved by DDASKR with three different method options: */
/* (1) direct band method for the linear systems (internal Jacobian), */
/* (2) preconditioned Krylov method for the linear systems, without */
/*     block-grouping in the reaction-based factor, and */
/* (3) preconditioned Krylov method for the linear systems, with */
/*     block-grouping in the reaction-based factor. */

/* In the Krylov cases, the preconditioner is the product of two factors: */
/* (a) The spatial factor uses a fixed number of Gauss-Seidel iterations */
/* based on the diffusion terms only. */
/* (b) The reaction-based factor is a block-diagonal matrix based on */
/* the partial derivatives of the interaction terms R only. */
/* With block-grouping, only a subset of the ns by ns blocks are computed. */
/* An integer flag, JPRE, is set in the main program to specify whether */
/* the preconditioner is to use only one of the two factors or both, */
/* and in which order. */

/* The reaction-based preconditioner factor is set up and solved in */
/* seven subroutines -- */
/*   DMSET2, DRBDJA, DRBDPS  in the case of no block-grouping, and */
/*   DGSET2, GSET1, DRBGJA, DRBGPS  in the case of block-grouping. */
/* These routines are provided separately for general use on problems */
/* arising from reaction-transport systems. */

/* Two output files are written.. one with the problem description and */
/* performance statistics on unit LOUT = 9, and one with solution */
/* profiles at selected output times on unit LCOUT = 10. */
/* The solution file is written only in the case of the direct method. */
/* ----------------------------------------------------------------------- */
/* Note.. in addition to the main program and subroutines given below, */
/* this program requires the BLAS routine DAXPY. */
/* ----------------------------------------------------------------------- */
/* References */
/* [1] Peter N. Brown and Alan C. Hindmarsh, */
/*     Reduced Storage Matrix Methods in Stiff ODE Systems, */
/*     J. Appl. Math. & Comp., 31 (1989), pp. 40-91. */
/* [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, */
/*     Using Krylov Methods in the Solution of Large-Scale Differential- */
/*     Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488. */
/* [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, */
/*     Consistent Initial Condition Calculation for Differential- */
/*     Algebraic Systems, LLNL Report UCRL-JC-122175, August 1995; */
/*     submitted to SIAM J. Sci. Comp. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   SETPAR, DGSET2, CINIT, DDASKR, OUTWEB */

/* ***END PROLOGUE  DWEB */

/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static integer lout = 9;
    static integer lcout = 10;

    /* Format strings */
    static char fmt_20[] = "(\002 DWEB: Example program for DDASKR packag"
	    "e\002//\002 Food web problem with NS species, NS =\002,i4/\002 P"
	    "redator-prey interaction and diffusion on a 2-D square\002/)";
    static char fmt_25[] = "(\002 Matrix parameters..  a =\002,e12.4,\002   "
	    "e =\002,e12.4,\002   g =\002,e12.4/21x,\002 b parameter =\002,e1"
	    "2.4//\002 Diffusion coefficients.. dprey =\002,e12.4,\002   dpre"
	    "d =\002,e12.4//\002 Rate parameters alpha =\002,e12.4,\002 and b"
	    "eta =\002,e12.4/)";
    static char fmt_30[] = "(\002 Mesh dimensions (MX,MY) =\002,2i4,5x,\002 "
	    "Total system size is NEQ =\002,i7/)";
    static char fmt_35[] = "(\002 Root function is R(Y) = average(c1) - 2"
	    "0\002/)";
    static char fmt_70[] = "(\002 Tolerance parameters.. RTOL =\002,e10.2"
	    ",\002   ATOL =\002,e10.2//\002 Internal I.C. calculation flag IN"
	    "FO(11) =\002,i2,\002   (0 = off, 1 = on)\002/\002 Predator I.C. "
	    "guess =\002,e10.2//\002 Alternate error test flag INFO(16) =\002"
	    ",i2,\002  (0 = off, 1 = on)\002)";
    static char fmt_80[] = "(//80(\002.\002)//\002 Linear solver method flag"
	    " INFO(12) =\002,i2,\002   (0 = direct, 1 = Krylov)\002/)";
    static char fmt_90[] = "(\002 Difference-quotient banded Jacobian,\002"
	    ",\002 half-bandwidths =\002,i4)";
    static char fmt_100[] = "(\002 Preconditioner flag is JPRE =\002,i3/\002"
	    "  (1 = reaction factor A_R, 2 = spatial factor A_S,\002,\002 3 ="
	    " A_S*A_R, 4 = A_R*A_S )\002/)";
    static char fmt_110[] = "(\002 No block-grouping in reaction factor\002)";
    static char fmt_120[] = "(\002 Block-grouping in reaction factor\002/"
	    "\002 Number of groups =\002,i5,\002   (NGX by NGY, NGX =\002,i3"
	    ",\002,  NGY =\002,i3,\002)\002)";
    static char fmt_140[] = "(//\002   t\002,12x,\002Ave.c1  NSTEP  NRE  NNI"
	    "  NLI  NPE  NQ\002,4x,\002H\002,10x,\002AVLIN\002)";
    static char fmt_160[] = "(e13.5,f10.5,i5,i6,3i5,i4,e11.2,f9.4)";
    static char fmt_165[] = "(15x,\002*****   Root found, JROOT =\002,i3)";
    static char fmt_170[] = "(//\002 Final time reached =\002,e12.4//)";
    static char fmt_220[] = "(//\002 Final statistics for this run..\002/"
	    "\002 RWORK size =\002,i8,\002   IWORK size =\002,i6/\002 Number "
	    "of time steps            =\002,i5/\002 Number of residual evalua"
	    "tions  =\002,i5/\002 Number of root fn. evaluations  =\002,i5"
	    "/\002 Number of Jac. or prec. evals.  =\002,i5/\002 Number of pr"
	    "econditioner solves =\002,i5/\002 Number of nonlinear iterations"
	    "  =\002,i5/\002 Number of linear iterations     =\002,i5/\002 Av"
	    "erage Krylov subspace dimension =\002,f8.4/i3,\002 nonlinear con"
	    "v. failures,\002,i5,\002 linear conv. failures\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal t, cc[800];
    static integer ng;
    static doublereal hu;
    static integer jbg, nli, neq, nni, nre, npe, nxg, nyg, nps, nrt, nst, nqu;
    extern /* Subroutine */ int avc1_(doublereal *, doublereal *);
    static integer idid, ncfl, ncfn, info[20], ipar[2], meth;
    static doublereal atol;
    static integer jpre;
    static doublereal rpar[800];
    static integer nrte;
    static doublereal rtol;
    static integer iout, nout;
    static doublereal tout, c1ave;
    static integer imod3;
    extern /* Subroutine */ int jacrs_();
    extern /* Subroutine */ int cinit_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), setid_(integer *, integer *, integer *, integer 
	    *, integer *, integer *);
    static doublereal avlin;
    static integer leniw;
    extern /* Subroutine */ int rtweb_();
    static integer lenrw, iwork[1640], jroot;
    static doublereal rwork[104883];
    extern /* Subroutine */ int dgset2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), dmset2_(
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static integer nlidif;
    static doublereal predic;
    static integer nnidif;
    extern /* Subroutine */ int ddaskr_(U_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, U_fp, U_fp, U_fp, integer *, 
	    integer *);
    extern /* Subroutine */ int resweb_();
    extern /* Subroutine */ int setpar_(void), outweb_(doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */ int psolrs_();
    static doublereal ccprime[800];

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_25, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_35, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_70, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_165, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_170, 0 };
    static cilist io___60 = { 0, 0, 0, fmt_220, 0 };



/* Dimension solution arrays and work arrays. */

/* When INFO(12) = 0, with INFO(5) = 0, INFO(6) = 1: */
/*    The length required for RWORK is */
/*        60 + 3*NRT + (2*ML+MU+11)*NEQ + 2*(NEQ/(ML+MU+1) + 1) . */
/*    For MX = MY = (even number) and ML = MU = NS*MX, this length is */
/*        60 + 3*NRT + (3*NS*MX + 11)*NEQ + MY . */
/*    The length required for IWORK is  40 + 2*NEQ . */

/* When INFO(12) = 1: */
/*    The length required for RWORK is */
/*        101 + 3*NRT + 19*NEQ + LENWP = 104 + 19*NEQ + NS*NS*NGRP . */
/*    The length required for IWORK is */
/*         40 + NEQ + LENIWP = 40 + NEQ + NS*NGRP . */

/* The dimensions for the various arrays are set below using parameters */
/*   MAXN    which must be .ge. NEQ = NS*MX*MY, */
/*   MAXS    which must be .ge. NS, */
/*   MAXM    which must be .ge. MAX(MX,MY). */



/* The COMMON blocks /PPAR1/ and /PPAR2/ contain problem parameters. */


/* Set output unit numbers for main output and tabulated solution. */

/* Open output files. */
    o__1.oerr = 0;
    o__1.ounit = lout;
    o__1.ofnmlen = 5;
    o__1.ofnm = "wdout";
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = lcout;
    o__1.ofnmlen = 6;
    o__1.ofnm = "wccout";
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/* Call SETPAR to set basic problem parameters. */
    setpar_();

/* Set remaining problem parameters. */
    neq = ppar2_1.ns * ppar2_1.mx * ppar2_1.my;
    ppar2_1.mxns = ppar2_1.mx * ppar2_1.ns;
    ppar2_1.dx = ppar2_1.ax / (real) (ppar2_1.mx - 1);
    ppar2_1.dy = ppar2_1.ay / (real) (ppar2_1.my - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = ppar2_1.dx;
	ppar2_1.cox[i__ - 1] = ppar2_1.diff[i__ - 1] / (d__1 * d__1);
/* L10: */
/* Computing 2nd power */
	d__1 = ppar2_1.dy;
	ppar2_1.coy[i__ - 1] = ppar2_1.diff[i__ - 1] / (d__1 * d__1);
    }

/* Set NRT = number of root functions. */
    nrt = 1;

    io___6.ciunit = lout;
    s_wsfe(&io___6);
    do_fio(&c__1, (char *)&ppar2_1.ns, (ftnlen)sizeof(integer));
    e_wsfe();
    io___7.ciunit = lout;
    s_wsfe(&io___7);
    do_fio(&c__1, (char *)&ppar1_1.aa, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar1_1.ee, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar1_1.gg, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar1_1.bb, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar1_1.dprey, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar1_1.dpred, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar2_1.alph, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ppar2_1.beta, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___8.ciunit = lout;
    s_wsfe(&io___8);
    do_fio(&c__1, (char *)&ppar2_1.mx, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ppar2_1.my, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&neq, (ftnlen)sizeof(integer));
    e_wsfe();
    io___9.ciunit = lout;
    s_wsfe(&io___9);
    e_wsfe();

/* Here set the flat initial guess for the predators. */
    predic = 1e5;

/* Set remaining method parameters for DDASKR. */
/* These include the INFO array and tolerances. */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L50: */
	info[i__ - 1] = 0;
    }

/* Here set INFO(11) = 1, indicating I.C. calculation requested. */
    info[10] = 1;

/* Here set INFO(14) = 1 to get the computed initial values. */
    info[13] = 1;

/* Here set INFO(15) = 1 to signal that a preconditioner setup routine */
/* is to be called in the Krylov case. */
    info[14] = 1;

/* Here set INFO(16) = 1 to get alternative error test (on the */
/* differential variables only). */
    info[15] = 1;

/* Here set the tolerances. */
    rtol = 1e-5;
    atol = rtol;

    io___14.ciunit = lout;
    s_wsfe(&io___14);
    do_fio(&c__1, (char *)&rtol, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&atol, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&info[10], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&predic, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&info[15], (ftnlen)sizeof(integer));
    e_wsfe();

/* Set NOUT = number of output times. */
    nout = 18;

/* Loop over method options: */
/* METH = 0 means use INFO(12) = 0 (direct) */
/* METH = 1 means use INFO(12) = 1 (Krylov) without block-grouping in */
/*          the reaction-based factor in the preconditioner. */
/* METH = 2 means use INFO(12) = 1 (Krylov) with block-grouping in */
/*          the reaction-based factor in the preconditioner. */
/* A block-grouping flag JBG, communicated through IPAR, is set to */
/* 0 (no block-grouping) or 1 (use block-grouping) with METH = 1 or 2. */
/* Reset INFO(1) = 0 and INFO(11) = 1. */

    for (meth = 0; meth <= 2; ++meth) {
	info[11] = min(meth,1);
	info[0] = 0;
	info[10] = 1;
	jbg = meth - 1;
	ipar[1] = jbg;

	io___19.ciunit = lout;
	s_wsfe(&io___19);
	do_fio(&c__1, (char *)&info[11], (ftnlen)sizeof(integer));
	e_wsfe();

/* In the case of the direct method, set INFO(6) = 1 to signal a banded */
/* Jacobian, set IWORK(1) = IWORK(2) = MX*NS, the half-bandwidth, and */
/* call SETID to set the IWORK segment ID indicating the differential */
/* and algebraic components. */
	if (info[11] == 0) {
	    info[5] = 1;
	    iwork[0] = ppar2_1.mxns;
	    iwork[1] = ppar2_1.mxns;
	    setid_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &c__40,
		     iwork);
	    io___21.ciunit = lout;
	    s_wsfe(&io___21);
	    do_fio(&c__1, (char *)&ppar2_1.mxns, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/* In the case of the Krylov method, set and print various */
/* preconditioner parameters. */
	if (info[11] == 1) {
/* First set the preconditioner choice JPRE. */
/*  JPRE = 1 means reaction-only (block-diagonal) factor A_R */
/*  JPRE = 2 means spatial factor (Gauss-Seidel) A_S */
/*  JPRE = 3 means A_S * A_R */
/*  JPRE = 4 means A_R * A_S */
/* Use IPAR to communicate JPRE to the preconditioner solve routine. */
	    jpre = 3;
	    ipar[0] = jpre;
	    io___23.ciunit = lout;
	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&jpre, (ftnlen)sizeof(integer));
	    e_wsfe();
/* Here call DMSET2 if JBG = 0, or DGSET2 if JBG = 1, to set the 2-D */
/* mesh parameters and block-grouping data, and the IWORK segment ID */
/* indicating the differential and algebraic components. */
	    if (jbg == 0) {
		dmset2_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &
			c__40, iwork);
		io___24.ciunit = lout;
		s_wsfe(&io___24);
		e_wsfe();
	    }
	    if (jbg == 1) {
		nxg = 5;
		nyg = 5;
		ng = nxg * nyg;
		dgset2_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &
			nxg, &nyg, &c__40, iwork);
		io___28.ciunit = lout;
		s_wsfe(&io___28);
		do_fio(&c__1, (char *)&ng, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nxg, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nyg, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	}

/* Set the initial T and TOUT, and call CINIT to set initial values. */
	t = 0.;
	tout = 1e-8;
	cinit_(cc, ccprime, &predic, rpar);

	nli = 0;
	nni = 0;

	io___36.ciunit = lout;
	s_wsfe(&io___36);
	e_wsfe();

/* Loop over output times and call DDASKR.  At each output time, */
/* print average c1 value and performance data. */
/* The first call, with IOUT = 0, is to calculate initial values only. */
/* After the first call, reset INFO(11) = 0 and the initial TOUT. */
/* If a root was found, we flag this, and return to the DDASKR call. */

	i__1 = nout;
	for (iout = 0; iout <= i__1; ++iout) {

L150:
	    ddaskr_((U_fp)resweb_, &neq, &t, cc, ccprime, &tout, info, &rtol, 
		    &atol, &idid, rwork, &c_b58, iwork, &c__1640, rpar, ipar, 
		    (U_fp)jacrs_, (U_fp)psolrs_, (U_fp)rtweb_, &nrt, &jroot);

	    nst = iwork[10];
	    nre = iwork[11];
	    npe = iwork[12];
	    nnidif = iwork[18] - nni;
	    nni = iwork[18];
	    nlidif = iwork[19] - nli;
	    nli = iwork[19];
	    nqu = iwork[7];
	    hu = rwork[6];
	    avlin = 0.;
	    if (nnidif > 0) {
		avlin = (real) nlidif / (real) nnidif;
	    }

	    if (meth == 0) {
		imod3 = iout - iout / 3 * 3;
		if (imod3 == 0) {
		    outweb_(&t, cc, &ppar2_1.ns, &ppar2_1.mx, &ppar2_1.my, &
			    lcout);
		}
	    }

	    avc1_(cc, &c1ave);
	    io___51.ciunit = lout;
	    s_wsfe(&io___51);
	    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&c1ave, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&nst, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nre, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nni, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nli, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&npe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nqu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&hu, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&avlin, (ftnlen)sizeof(doublereal));
	    e_wsfe();

	    if (idid == 5) {
		io___52.ciunit = lout;
		s_wsfe(&io___52);
		do_fio(&c__1, (char *)&jroot, (ftnlen)sizeof(integer));
		e_wsfe();
		goto L150;
	    }
	    if (idid < 0) {
		io___53.ciunit = lout;
		s_wsfe(&io___53);
		do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
		e_wsfe();
		goto L210;
	    }

	    if (tout > .9) {
		tout += 1.;
	    }
	    if (tout < .9) {
		tout *= 10.;
	    }
	    if (iout == 0) {
		info[10] = 0;
		tout = 1e-8;
		nli = 0;
		nni = 0;
	    }
/* L200: */
	}

L210:
	lenrw = iwork[17];
	leniw = iwork[16];
	nst = iwork[10];
	nre = iwork[11];
	npe = iwork[12];
	nni = iwork[18];
	nli = iwork[19];
	nps = iwork[20];
	if (nni > 0) {
	    avlin = (real) nli / (real) nni;
	}
	ncfn = iwork[14];
	ncfl = iwork[15];
	nrte = iwork[35];
	io___60.ciunit = lout;
	s_wsfe(&io___60);
	do_fio(&c__1, (char *)&lenrw, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&leniw, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nst, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nre, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrte, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&npe, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nps, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nni, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nli, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&avlin, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ncfn, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ncfl, (ftnlen)sizeof(integer));
	e_wsfe();

/* L300: */
    }
    s_stop("", (ftnlen)0);
/* ------  End of main program for DWEB example program ------------------ */
    return 0;
} /* MAIN__ */

/* Subroutine */ int setpar_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal pi;

/* ----------------------------------------------------------------------- */
/* This routine sets the basic problem parameters, namely */
/* AX, AY, NS, MX, MY,  problem coefficients ACOEF, BCOEF, DIFF, */
/* ALPH, BETA, using parameters NP, AA, EE, GG, BB, DPREY, DPRED. */
/* ----------------------------------------------------------------------- */

    ppar2_1.ax = 1.;
    ppar2_1.ay = 1.;
    ppar2_1.np = 1;
    ppar2_1.mx = 20;
    ppar2_1.my = 20;
    ppar1_1.aa = 1.;
    ppar1_1.ee = 1e4;
    ppar1_1.gg = 5e-7;
    ppar1_1.bb = 1.;
    ppar1_1.dprey = 1.;
    ppar1_1.dpred = .05;
    ppar2_1.alph = 50.;
    ppar2_1.beta = 100.;
    ppar2_1.ns = ppar2_1.np << 1;
    i__1 = ppar2_1.np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ppar2_1.np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ppar2_1.acoef[ppar2_1.np + i__ + (j << 1) - 3] = ppar1_1.ee;
	    ppar2_1.acoef[i__ + (ppar2_1.np + j << 1) - 3] = -ppar1_1.gg;
/* L10: */
	}
	ppar2_1.acoef[j + (j << 1) - 3] = -ppar1_1.aa;
	ppar2_1.acoef[ppar2_1.np + j + (ppar2_1.np + j << 1) - 3] = 
		-ppar1_1.aa;
	ppar2_1.bcoef[j - 1] = ppar1_1.bb;
	ppar2_1.bcoef[ppar2_1.np + j - 1] = -ppar1_1.bb;
	ppar2_1.diff[j - 1] = ppar1_1.dprey;
	ppar2_1.diff[ppar2_1.np + j - 1] = ppar1_1.dpred;
/* L20: */
    }
    pi = 3.141592653589793;
    ppar2_1.fpi = pi * 4.;

    return 0;
/* ------------  End of Subroutine SETPAR  ------------------------------- */
} /* setpar_ */

/* Subroutine */ int setid_(integer *mx, integer *my, integer *ns, integer *
	nsd, integer *lid, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i0, i00, jx, jy, nsdp1;

/* ----------------------------------------------------------------------- */
/* This routine sets the ID array in IWORK, indicating which components */
/* are differential and which are algebraic. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --iwork;

    /* Function Body */
    nsdp1 = *nsd + 1;
    i__1 = *my;
    for (jy = 1; jy <= i__1; ++jy) {
	i00 = *mx * *ns * (jy - 1) + *lid;
	i__2 = *mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    i0 = i00 + *ns * (jx - 1);
	    i__3 = *nsd;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L10: */
		iwork[i0 + i__] = 1;
	    }
	    i__3 = *ns;
	    for (i__ = nsdp1; i__ <= i__3; ++i__) {
/* L20: */
		iwork[i0 + i__] = -1;
	    }
/* L30: */
	}
/* L40: */
    }

    return 0;
/* ------------  End of Subroutine SETID  -------------------------------- */
} /* setid_ */

/* Subroutine */ int cinit_(doublereal *cc, doublereal *ccprime, doublereal *
	predic, doublereal *rpar)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t, x, y;
    static integer jx, jy;
    static doublereal fac;
    static integer npp1, ioff;
    extern /* Subroutine */ int fweb_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static doublereal argx, argy;
    static integer iyoff;

/* ----------------------------------------------------------------------- */
/* This routine computes and loads the vectors of initial values. */
/* ----------------------------------------------------------------------- */

/* Load CC. */
    /* Parameter adjustments */
    --rpar;
    --ccprime;
    --cc;

    /* Function Body */
    npp1 = ppar2_1.np + 1;
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	y = (real) (jy - 1) * ppar2_1.dy;
	argy = y * 16. * y * (ppar2_1.ay - y) * (ppar2_1.ay - y);
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    x = (real) (jx - 1) * ppar2_1.dx;
	    argx = x * 16. * x * (ppar2_1.ax - x) * (ppar2_1.ax - x);
	    ioff = iyoff + ppar2_1.ns * (jx - 1);
	    fac = ppar2_1.alph * x * y + 1. + ppar2_1.beta * sin(ppar2_1.fpi *
		     x) * sin(ppar2_1.fpi * y);
	    i__3 = ppar2_1.np;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L10: */
		cc[ioff + i__] = (real) i__ * argx * argy + 10.;
	    }
	    i__3 = ppar2_1.ns;
	    for (i__ = npp1; i__ <= i__3; ++i__) {
/* L15: */
		cc[ioff + i__] = *predic;
	    }
/* L20: */
	}
/* L30: */
    }

/* Load CCPRIME. */
    t = 0.;
    fweb_(&t, &cc[1], &ccprime[1], &rpar[1]);
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ioff = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = npp1; i__ <= i__3; ++i__) {
/* L40: */
		ccprime[ioff + i__] = 0.;
	    }
/* L50: */
	}
/* L60: */
    }

    return 0;
/* ------------  End of Subroutine CINIT  -------------------------------- */
} /* cinit_ */

/* Subroutine */ int outweb_(doublereal *t, doublereal *c__, integer *ns, 
	integer *mx, integer *my, integer *lun)
{
    /* Format strings */
    static char fmt_10[] = "(/1x,79(\002-\002)/30x,\002At time t = \002,e16."
	    "8/1x,79(\002-\002))";
    static char fmt_20[] = "(\002 the species c(\002,i2,\002) values are "
	    "\002)";
    static char fmt_25[] = "(6(1x,g12.6))";
    static char fmt_35[] = "(1x,79(\002-\002),/)";

    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, jx, jy;

    /* Fortran I/O blocks */
    static cilist io___82 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___86 = { 0, 0, 0, fmt_25, 0 };
    static cilist io___88 = { 0, 0, 0, fmt_35, 0 };


/* ----------------------------------------------------------------------- */
/* This routine prints the values of the individual species densities */
/* at the current time T, to logical unit LUN. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    c_dim1 = *ns;
    c_dim2 = *mx;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    io___82.ciunit = *lun;
    s_wsfe(&io___82);
    do_fio(&c__1, (char *)&(*t), (ftnlen)sizeof(doublereal));
    e_wsfe();

    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___84.ciunit = *lun;
	s_wsfe(&io___84);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
	for (jy = *my; jy >= 1; --jy) {
	    io___86.ciunit = *lun;
	    s_wsfe(&io___86);
	    i__2 = *mx;
	    for (jx = 1; jx <= i__2; ++jx) {
		do_fio(&c__1, (char *)&c__[i__ + (jx + jy * c_dim2) * c_dim1],
			 (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
/* L30: */
	}
	io___88.ciunit = *lun;
	s_wsfe(&io___88);
	e_wsfe();
/* L40: */
    }

    return 0;
/* ------------  End of Subroutine OUTWEB  ------------------------------- */
} /* outweb_ */

/* Subroutine */ int resweb_(doublereal *t, doublereal *u, doublereal *uprime,
	 doublereal *cj, doublereal *delta, integer *ires, doublereal *rpar, 
	integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, jx, jy, ic0, ici;
    extern /* Subroutine */ int fweb_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static integer iyoff;

/* ----------------------------------------------------------------------- */
/* This routine computes the residual vector, using Subroutine FWEB */
/* for the right-hand sides. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --delta;
    --uprime;
    --u;

    /* Function Body */
    fweb_(t, &u[1], &delta[1], &rpar[1]);

    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ic0 = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic0 + i__;
		if (i__ > ppar2_1.np) {
		    delta[ici] = -delta[ici];
		} else {
		    delta[ici] = uprime[ici] - delta[ici];
		}
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }

    return 0;
/* ------------  End of Subroutine RESWEB  ------------------------------- */
} /* resweb_ */

/* Subroutine */ int fweb_(doublereal *t, doublereal *cc, doublereal *crate, 
	doublereal *rpar)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ic, jx, jy, ici;
    extern /* Subroutine */ int webr_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);
    static integer idxl, idyl, idxu, idyu;
    static doublereal dcxli, dcyli;
    static integer iyoff;
    static doublereal dcyui, dcxui;

/* ----------------------------------------------------------------------- */
/* This routine computes the right-hand sides of all the equations */
/* and returns them in the array CRATE. */
/* The interaction rates are computed by calls to WEBR, and these are */
/* saved in RPAR(1),...,RPAR(NEQ) for use later in preconditioning. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --rpar;
    --crate;
    --cc;

    /* Function Body */
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	idyu = ppar2_1.mxns;
	if (jy == ppar2_1.my) {
	    idyu = -ppar2_1.mxns;
	}
	idyl = ppar2_1.mxns;
	if (jy == 1) {
	    idyl = -ppar2_1.mxns;
	}
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1) + 1;
/* Get interaction rates at one point (X,Y). */
	    webr_(t, &jx, &jy, &cc[ic], &rpar[ic]);
	    idxu = ppar2_1.ns;
	    if (jx == ppar2_1.mx) {
		idxu = -ppar2_1.ns;
	    }
	    idxl = ppar2_1.ns;
	    if (jx == 1) {
		idxl = -ppar2_1.ns;
	    }
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__ - 1;
/* Do differencing in Y. */
		dcyli = cc[ici] - cc[ici - idyl];
		dcyui = cc[ici + idyu] - cc[ici];
/* Do differencing in X. */
		dcxli = cc[ici] - cc[ici - idxl];
		dcxui = cc[ici + idxu] - cc[ici];
/* Collect terms and load CRATE elements. */
		crate[ici] = ppar2_1.coy[i__ - 1] * (dcyui - dcyli) + 
			ppar2_1.cox[i__ - 1] * (dcxui - dcxli) + rpar[ici];
/* L20: */
	    }
/* L40: */
	}
/* L60: */
    }
    return 0;
/* ------------  End of Subroutine FWEB  --------------------------------- */
} /* fweb_ */

/* Subroutine */ int webr_(doublereal *t, integer *jx, integer *jy, 
	doublereal *c__, doublereal *crate)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal x, y, fac;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ----------------------------------------------------------------------- */
/* This routine computes one block of the interaction term R of the */
/* system, namely block (JX,JY), for use in preconditioning. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --crate;
    --c__;

    /* Function Body */
    y = (real) (*jy - 1) * ppar2_1.dy;
    x = (real) (*jx - 1) * ppar2_1.dx;
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	crate[i__] = 0.;
    }
    i__1 = ppar2_1.ns;
    for (j = 1; j <= i__1; ++j) {
	daxpy_(&ppar2_1.ns, &c__[j], &ppar2_1.acoef[(j << 1) - 2], &c__1, &
		crate[1], &c__1);
/* L15: */
    }
    fac = ppar2_1.alph * x * y + 1. + ppar2_1.beta * sin(ppar2_1.fpi * x) * 
	    sin(ppar2_1.fpi * y);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	crate[i__] = c__[i__] * (ppar2_1.bcoef[i__ - 1] * fac + crate[i__]);
    }
    return 0;
/* ------------  End of Subroutine WEBR  --------------------------------- */
} /* webr_ */

/* Subroutine */ int jacrs_(doublereal *res, integer *ires, integer *neq, 
	doublereal *t, doublereal *cc, doublereal *ccprime, doublereal *rewt, 
	doublereal *savr, doublereal *wk, doublereal *h__, doublereal *cj, 
	doublereal *wp, integer *iwp, integer *ier, doublereal *rpar, integer 
	*ipar)
{
    static integer jbg;
    extern /* Subroutine */ int webr_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *), drbdja_(doublereal *, doublereal *, 
	    doublereal *, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), drbgja_(doublereal *, 
	    doublereal *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* This routine interfaces to Subroutine DRBDJA or Subroutine DRBGJA, */
/* depending on the flag JBG = IPAR(2), to generate and preprocess the */
/* block-diagonal Jacobian corresponding to the reaction term R. */
/* If JBG = 0, we call DRBDJA, with no block-grouping. */
/* If JBG = 1, we call DRBGJA, and use block-grouping. */
/* Array RPAR, containing the current R vector, is passed to DRBDJA and */
/* DRBGJA as argument R0, consistent with the loading of RPAR in FWEB. */
/* The external name WEBR is passed, as the routine which computes the */
/* individual blocks of the R vector. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwp;
    --wp;
    --wk;
    --savr;
    --rewt;
    --ccprime;
    --cc;

    /* Function Body */
    jbg = ipar[2];
    if (jbg == 0) {
	drbdja_(t, &cc[1], &rpar[1], (S_fp)webr_, &wk[1], &rewt[1], cj, &wp[1]
		, &iwp[1], ier);
    } else {
	drbgja_(t, &cc[1], &rpar[1], (S_fp)webr_, &wk[1], &rewt[1], cj, &wp[1]
		, &iwp[1], ier);
    }

    return 0;
/* ------------  End of Subroutine JACRS  -------------------------------- */
} /* jacrs_ */

/* Subroutine */ int psolrs_(integer *neq, doublereal *t, doublereal *cc, 
	doublereal *ccprime, doublereal *savr, doublereal *wk, doublereal *cj,
	 doublereal *wt, doublereal *wp, integer *iwp, doublereal *b, 
	doublereal *eplin, integer *ier, doublereal *rpar, integer *ipar)
{
    extern /* Subroutine */ int gs_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal hl0;
    static integer jbg, jpre;
    extern /* Subroutine */ int drbdps_(doublereal *, doublereal *, integer *)
	    , drbgps_(doublereal *, doublereal *, integer *);

/* ----------------------------------------------------------------------- */
/* This routine applies the inverse of a product preconditioner matrix */
/* to the vector in the array B.  Depending on the flag JPRE, this */
/* involves a call to GS, for the inverse of the spatial factor, and/or */
/* a call to DRBDPS or DRBGPS for the inverse of the reaction-based */
/* factor (CJ*I_d - dR/dy).  The latter factor uses block-grouping */
/* (with a call to DRBGPS) if JBG = 1, and does not (with a call to */
/* DRBDPS) if JBG = 0.  JBG is communicated as IPAR(2). */
/* The array B is overwritten with the solution. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --b;
    --iwp;
    --wp;
    --wk;
    --savr;
    --ccprime;
    --cc;

    /* Function Body */
    jpre = ipar[1];
    *ier = 0;
    hl0 = 1. / *cj;

    jbg = ipar[2];

    if (jpre == 2 || jpre == 3) {
	gs_(neq, &hl0, &b[1], &wk[1]);
    }

    if (jpre != 2) {
	if (jbg == 0) {
	    drbdps_(&b[1], &wp[1], &iwp[1]);
	}
	if (jbg == 1) {
	    drbgps_(&b[1], &wp[1], &iwp[1]);
	}
    }

    if (jpre == 4) {
	gs_(neq, &hl0, &b[1], &wk[1]);
    }

    return 0;
/* ------------  End of Subroutine PSOLRS  ------------------------------- */
} /* psolrs_ */

/* Subroutine */ int gs_(integer *n, doublereal *hl0, doublereal *z__, 
	doublereal *x)
{
    /* Initialized data */

    static integer itmax = 5;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ic, ii, jx, jy, ici;
    static doublereal dinv[2];
    static integer iter;
    static doublereal beta1[2], beta2[2];
    static integer iyoff;
    static doublereal gamma1[2], gamma2[2], elamda;

/* ----------------------------------------------------------------------- */
/* This routine provides the inverse of the spatial factor for a */
/* product preconditoner in an NS-species reaction-diffusion problem. */
/* It performs ITMAX = 5 Gauss-Seidel iterations to compute an */
/* approximation to (A_S)-inverse * Z, where A_S = I - hl0*Jd, and Jd */
/* represents the diffusion contributions to the Jacobian. */
/* The solution is vector returned in Z. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --x;
    --z__;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Write matrix as A = D - L - U. */
/* Load local arrays BETA, BETA2, GAMMA1, GAMMA2, and DINV. */
/* ----------------------------------------------------------------------- */
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	elamda = 1. / (*hl0 * 2. * (ppar2_1.cox[i__ - 1] + ppar2_1.coy[i__ - 
		1]) + 1.);
	beta1[i__ - 1] = *hl0 * ppar2_1.cox[i__ - 1] * elamda;
	beta2[i__ - 1] = beta1[i__ - 1] * 2.;
	gamma1[i__ - 1] = *hl0 * ppar2_1.coy[i__ - 1] * elamda;
	gamma2[i__ - 1] = gamma1[i__ - 1] * 2.;
	dinv[i__ - 1] = elamda;
/* L10: */
    }
/* ----------------------------------------------------------------------- */
/* Begin iteration loop. */
/* Load array X with (D-inverse)*Z for first iteration. */
/* ----------------------------------------------------------------------- */
    iter = 1;
/* Zero X in all its components, since X is added to Z at the end. */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
/* L15: */
	x[ii] = 0.;
    }

    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__;
		x[ici] = dinv[i__ - 1] * z__[ici];
		z__[ici] = 0.;
/* L30: */
	    }
/* L40: */
	}
/* L50: */
    }
    goto L160;
/* ----------------------------------------------------------------------- */
/* Calculate (D-inverse)*U*X. */
/* ----------------------------------------------------------------------- */
L70:
    ++iter;
    jy = 1;
    jx = 1;
    ic = ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L75: */
	x[ici] = beta2[i__ - 1] * x[ici + ppar2_1.ns] + gamma2[i__ - 1] * x[
		ici + ppar2_1.mxns];
    }
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L80: */
	    x[ici] = beta1[i__ - 1] * x[ici + ppar2_1.ns] + gamma2[i__ - 1] * 
		    x[ici + ppar2_1.mxns];
	}
/* L85: */
    }
    jx = ppar2_1.mx;
    ic = ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L90: */
	x[ici] = gamma2[i__ - 1] * x[ici + ppar2_1.mxns];
    }
    i__1 = ppar2_1.my - 1;
    for (jy = 2; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	jx = 1;
	ic = iyoff;
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L95: */
	    x[ici] = beta2[i__ - 1] * x[ici + ppar2_1.ns] + gamma1[i__ - 1] * 
		    x[ici + ppar2_1.mxns];
	}
	i__2 = ppar2_1.mx - 1;
	for (jx = 2; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__;
/* L100: */
		x[ici] = beta1[i__ - 1] * x[ici + ppar2_1.ns] + gamma1[i__ - 
			1] * x[ici + ppar2_1.mxns];
	    }
/* L105: */
	}
	jx = ppar2_1.mx;
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L110: */
	    x[ici] = gamma1[i__ - 1] * x[ici + ppar2_1.mxns];
	}
/* L115: */
    }
    jy = ppar2_1.my;
    iyoff = ppar2_1.mxns * (jy - 1);
    jx = 1;
    ic = iyoff;
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L120: */
	x[ici] = beta2[i__ - 1] * x[ici + ppar2_1.ns];
    }
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L125: */
	    x[ici] = beta1[i__ - 1] * x[ici + ppar2_1.ns];
	}
/* L130: */
    }
    jx = ppar2_1.mx;
    ic = iyoff + ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L135: */
	x[ici] = 0.;
    }
/* ----------------------------------------------------------------------- */
/* Calculate [(I - (D-inverse)*L)]-inverse * X. */
/* ----------------------------------------------------------------------- */
L160:
    jy = 1;
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L170: */
	    x[ici] += beta1[i__ - 1] * x[ici - ppar2_1.ns];
	}
/* L175: */
    }
    jx = ppar2_1.mx;
    ic = ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L180: */
	x[ici] += beta2[i__ - 1] * x[ici - ppar2_1.ns];
    }
    i__1 = ppar2_1.my - 1;
    for (jy = 2; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	jx = 1;
	ic = iyoff;
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L185: */
	    x[ici] += gamma1[i__ - 1] * x[ici - ppar2_1.mxns];
	}
	i__2 = ppar2_1.mx - 1;
	for (jx = 2; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__;
		x[ici] = x[ici] + beta1[i__ - 1] * x[ici - ppar2_1.ns] + 
			gamma1[i__ - 1] * x[ici - ppar2_1.mxns];
/* L195: */
	    }
/* L200: */
	}
	jx = ppar2_1.mx;
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
	    x[ici] = x[ici] + beta2[i__ - 1] * x[ici - ppar2_1.ns] + gamma1[
		    i__ - 1] * x[ici - ppar2_1.mxns];
/* L205: */
	}
/* L210: */
    }
    jy = ppar2_1.my;
    iyoff = ppar2_1.mxns * (jy - 1);
    jx = 1;
    ic = iyoff;
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L215: */
	x[ici] += gamma2[i__ - 1] * x[ici - ppar2_1.mxns];
    }
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
	    x[ici] = x[ici] + beta1[i__ - 1] * x[ici - ppar2_1.ns] + gamma2[
		    i__ - 1] * x[ici - ppar2_1.mxns];
/* L220: */
	}
/* L225: */
    }
    jx = ppar2_1.mx;
    ic = iyoff + ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
	x[ici] = x[ici] + beta2[i__ - 1] * x[ici - ppar2_1.ns] + gamma2[i__ - 
		1] * x[ici - ppar2_1.mxns];
/* L230: */
    }
/* ----------------------------------------------------------------------- */
/* Add increment X to Z. */
/* ----------------------------------------------------------------------- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L300: */
	z__[i__] += x[i__];
    }

    if (iter < itmax) {
	goto L70;
    }
    return 0;
/* ------------  End of Subroutine GS  ----------------------------------- */
} /* gs_ */

/* Subroutine */ int avc1_(doublereal *cc, doublereal *c1ave)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer jx, jy;
    static doublereal sum;
    static integer npp1, ioff, iyoff;

/* ----------------------------------------------------------------------- */
/* This routine computes the average value of c1. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --cc;

    /* Function Body */
    sum = 0.;
    npp1 = ppar2_1.np + 1;
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ioff = iyoff + ppar2_1.ns * (jx - 1);
	    sum += cc[ioff + 1];
/* L20: */
	}
/* L30: */
    }

    *c1ave = sum / (ppar2_1.mx * ppar2_1.my);

    return 0;
/* ------------  End of Subroutine AVC1  --------------------------------- */
} /* avc1_ */

/* Subroutine */ int rtweb_(integer *neq, doublereal *t, doublereal *cc, 
	doublereal *cp, integer *nrt, doublereal *rval, doublereal *rpar, 
	integer *ipar)
{
    extern /* Subroutine */ int avc1_(doublereal *, doublereal *);
    static doublereal c1ave;


/* This routine sets RVAL = average(c1) - 20.0. */


    /* Parameter adjustments */
    --cp;
    --cc;

    /* Function Body */
    avc1_(&cc[1], &c1ave);
    *rval = c1ave - 20.;

    return 0;
/* ------------  End of Subroutine RTWEB  -------------------------------- */
} /* rtweb_ */

