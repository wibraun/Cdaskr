/* dheatilu.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__2594 = 2594;
static integer c__4468 = 4468;
static integer c__1 = 1;
static integer c__9 = 9;

/* ***BEGIN PROLOGUE  DHEATILU */
/* ***REFER TO  DDASKR */
/* ***DATE WRITTEN   020813   (YYMMDD) */
/* ***REVISION DATE */

/* ***DESCRIPTON */

/* ----------------------------------------------------------------------- */
/* Example program for DDASKR. */
/* DAE system derived from the discretized heat equation on a square. */

/* This is the double precision version. */
/* ----------------------------------------------------------------------- */

/* This program solves a DAE system that arises from the heat equation, */
/*   du/dt = u   + u */
/*            xx    yy */
/* posed on the 2-D unit square with zero Dirichlet boundary conditions. */
/* An M+2 by M+2 mesh is set on the square, with uniform spacing 1/(M+1). */
/* The spatial deriviatives are represented by standard central finite */
/* difference approximations.  At each interior point of the mesh, */
/* the discretized PDE becomes an ODE for the discrete value of u. */
/* At each point on the boundary, we pose the equation u = 0.  The */
/* discrete values of u form a vector U, ordered first by x, then by y. */
/* The result is a DAE system G(t,U,U') = 0 of size NEQ = (M+2)*(M+2). */

/* Initial conditions are posed as u = 16x(1-x)y(1-y) at t = 0. */
/* The problem is solved by DDASKR on the time interval t .le. 10.24. */

/* The root functions are R1(U) = max(u) - 0.1, R2(U) = max(u) - 0.01. */

/* The Krylov linear system solution method, with preconditioning, is */
/* selected.  The preconditioner is a sparse matrix with half-bandwidths */
/* equal to 1, i.e. a tridiagonal matrix.  (The true half-bandwidths */
/* are equal to M+2.)  This corresponds to ignoring the y-direction */
/* coupling in the ODEs, for purposes of preconditioning.  The extra */
/* iterations resulting from this approximation are offset by the lower */
/* storage and linear system solution costs for a tridiagonal matrix. */

/* The routines DJACILU and DPSOLILU that generate and solve the sparse */
/* preconditioner are provided in a separate file for general use. */

/* The output times are t = .01 * 2**n (n = 0,...,10).  The maximum of */
/* abs(u) over the mesh, and various performance statistics, are printed. */

/* For details and test results on this problem, see the reference: */
/*   Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, */
/*   Using Krylov Methods in the Solution of Large-Scale Differential- */
/*   Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   UINIT, DDASKR */

/* ***END PROLOGUE  DHEATILU */

/* Here are necessary declarations.  The dimension statements use a */
/* maximum value for the mesh parameter M. */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_15[] = "(\002 Error return from DSPSETUP: IERR = \002,i5)"
	    ;
    static char fmt_30[] = "(\002 DHEATILU: Heat Equation Example Program fo"
	    "r DDASKR\002//\002    M+2 by M+2 mesh, M =\002,i3,\002,  System "
	    "size NEQ =\002,i4/\002    Root functions are: R1 = max(u) - 0."
	    "1\002,\002 and R2 = max(u) - 0.01\002//\002    Linear solver met"
	    "hod flag INFO(12) =\002,i3,\002    (0 = direct, 1 = Krylov)\002"
	    "/\002    Preconditioner is a sparse approximation with ML =\002,"
	    "i3,\002  MU =\002,i3/\002    Incomplete factorization option "
	    "=\002,i2,\002    (1 = ILUT, 2 = ILUTP)\002/,\002    Tolerances a"
	    "re RTOL =\002,e10.1,\002   ATOL =\002,e10.1//)";
    static char fmt_40[] = "(5x,\002t\002,12x,\002UMAX\002,8x,\002NQ\002,8x"
	    ",\002H\002,8x,\002STEPS\002,5x,\002NNI\002,5x,\002NLI\002/)";
    static char fmt_60[] = "(e15.5,e12.4,i5,e14.3,i7,i9,i8)";
    static char fmt_61[] = "(20x,\002*****   Root found, JROOT =\002,2i3)";
    static char fmt_65[] = "(//\002 Final time reached =\002,e12.4//)";
    static char fmt_90[] = "(//\002 Final statistics for this run..\002/\002"
	    "   RWORK size =\002,i5,\002   IWORK size =\002,i4/\002   Number "
	    "of time steps ................ =\002,i5/\002   Number of residua"
	    "l evaluations ...... =\002,i5/\002   Number of res. evals. for p"
	    "recond.    =\002,i5/\002   Number of root function evaluations ."
	    " =\002,i5/\002   Number of preconditioner evaluations  =\002,i5"
	    "/\002   Number of preconditioner solves ..... =\002,i5/\002   Nu"
	    "mber of nonlinear iterations ...... =\002,i5/\002   Number of li"
	    "near iterations ......... =\002,i5/\002   Average Krylov subspac"
	    "e dimension =\002,f8.4/i5,\002 nonlinear conv. failures,\002,i5"
	    ",\002 linear conv. failures\002)";
    static char fmt_100[] = "(\002 Minimum lengths for work arrays WP and IW"
	    "P: \002,i7,1x,i7)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int dpsolilu_();
    extern /* Subroutine */ int dspsetup_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static integer i__, m;
    static doublereal t, u[144];
    static integer ml;
    static doublereal dx, hu;
    static integer mu, nli, neq, nni, npe, nre, liw, nps, nrt, lrw, nqu, nst, 
	    idid, ncfl, ncfn, info[20], ipar[34];
    static doublereal atol;
    extern /* Subroutine */ int resh_();
    static integer ierr;
    static doublereal rpar[4];
    static integer nrte;
    static doublereal umax, rtol;
    static integer iout, lout, nout;
    static doublereal tout;
    static integer mband;
    static doublereal coeff, avdim;
    extern /* Subroutine */ int uinit_(doublereal *, doublereal *, doublereal 
	    *, integer *);
    static integer iwork[4508], jroot[2];
    static doublereal rwork[5293];
    extern /* Subroutine */ int ddaskr_(U_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, U_fp, U_fp, U_fp, integer *, 
	    integer *);
    extern /* Subroutine */ int rtheat_();
    static doublereal uprime[144];
    static integer lwpmin;
    extern /* Subroutine */ int djacilu_();
    static integer liwpmin;

    /* Fortran I/O blocks */
    static cilist io___17 = { 0, 0, 0, fmt_15, 0 };
    static cilist io___18 = { 0, 0, 0, 0, 0 };
    static cilist io___19 = { 0, 0, 0, 0, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_61, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_65, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_100, 0 };


/* Load parameters for sparse preconditioner. */
/* =2 means ILUTP preconditioner used */
/* =1 means ILUT preconditioner used */
/* Load workspace lengths. */
/* Set LOUT, the unit number of the output device. */
    lout = 6;
/* Open matrix output file if JACOUT .EQ. 1. */
    if (FALSE_) {
	ipar[28] = 1;
	o__1.oerr = 0;
	o__1.ounit = 1;
	o__1.ofnmlen = 20;
	o__1.ofnm = "Heat_Test_Matrix.dat";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

/* Here set parameters for the problem being solved.  Use RPAR and IPAR */
/* to communicate these to the other routines. */

    m = 10;
    dx = 1. / (m + 1);
    neq = (m + 2) * (m + 2);
    coeff = 1. / (dx * dx);
    ipar[32] = neq;
    ipar[33] = m;
    rpar[2] = dx;
    rpar[3] = coeff;

/* Set NRT = number of root functions. */
    nrt = 2;

/* Here set the lengths of the preconditioner work arrays WP and IWP, */
/* load them into IWORK, and set the total lengths of WORK and IWORK. */
    iwork[26] = 2594;
    iwork[27] = 4468;
    lrw = 5293;
    liw = 4508;
/* Load values into IPAR and RPAR for sparse preconditioner. */
    ml = 1;
    mu = 1;
    ipar[0] = ml;
    ipar[1] = mu;
    ipar[2] = 5;
    ipar[3] = 5;
    ipar[4] = 1;
    ipar[5] = 5;
    ipar[6] = 1;
    ipar[7] = 1;
    ipar[8] = 2;
    ipar[9] = 0;
    ipar[10] = 1;
    ipar[29] = 0;
    rpar[0] = .001;
    rpar[1] = .01;
/* Check IPAR, RPAR, LENWP and LENIWP for illegal entries and long */
/* enough work array lengths. */
    dspsetup_(&neq, &c__2594, &c__4468, rpar, ipar, &ierr, &lwpmin, &liwpmin);
    if (ierr != 0) {
	io___17.ciunit = lout;
	s_wsfe(&io___17);
	do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	e_wsfe();
	if (lwpmin > 2594) {
	    io___18.ciunit = lout;
	    s_wsle(&io___18);
	    do_lio(&c__9, &c__1, " More WP work array length needed", (ftnlen)
		    33);
	    e_wsle();
	}
	if (liwpmin > 4468) {
	    io___19.ciunit = lout;
	    s_wsle(&io___19);
	    do_lio(&c__9, &c__1, " More IWP work array length needed", (
		    ftnlen)34);
	    e_wsle();
	}
	s_stop("", (ftnlen)0);
    }

/* Call subroutine UINIT to initialize U and UPRIME. */

    uinit_(u, uprime, rpar, ipar);

/* ----------------------------------------------------------------------- */
/* Here we set up the INFO array, which describes the various options */
/* in the way we want DDASKR to solve the problem. */
/* In this case, we select the iterative preconditioned Krylov method, */
/* and we supply the sparse preconditioner routines DJACILU/DPSOLILU. */

/* We first initialize the entire INFO array to zero, then set select */
/* entries to nonzero values for desired solution options. */

/* To select the Krylov iterative method for the linear systems, */
/* we set INFO(12) = 1. */

/* Since we are using a preconditioner that involves approximate */
/* Jacobian elements requiring preprocessing, we have a JAC routine, */
/* namely subroutine DJACILU, and we must set INFO(15) = 1 to indicate */
/* this to DDASKR. */

/* No other entries of INFO need to be changed for this example. */
/* ----------------------------------------------------------------------- */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L20: */
	info[i__ - 1] = 0;
    }

    info[11] = 1;
    info[14] = 1;

/* Here we set tolerances for DDASKR to indicate how much accuracy */
/* we want in the solution, in the sense of local error control. */
/* For this example, we ask for pure absolute error control with a */
/* tolerance of 1.0D-5. */
    rtol = 0.;
    atol = 1e-5;

/* Here we generate a heading with important parameter values. */

    io___26.ciunit = lout;
    s_wsfe(&io___26);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&neq, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&info[11], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ml, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mu, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&rtol, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&atol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___27.ciunit = lout;
    s_wsfe(&io___27);
    e_wsfe();

/* ----------------------------------------------------------------------- */
/* Now we solve the problem. */

/* DDASKR will be called to compute 11 intermediate solutions from */
/* tout = 0.01 to tout = 10.24 by powers of 2. */

/* We pass to DDASKR the names DJACILU and DPSOLILU for the JAC and PSOL */
/* routines to do the preconditioning. */

/* At each output time, we compute and print the max-norm of the */
/* solution (which should decay exponentially in t).  We also print */
/* some relevant statistics -- the current method order and step size, */
/* the number of time steps so far, and the numbers of nonlinear and */
/* linear iterations so far. */

/* If a root was found, we flag this, and return to the DDASKR call. */

/* If DDASKR failed in any way (IDID .lt. 0) we print a message and */
/* stop the integration. */
/* ----------------------------------------------------------------------- */

    nout = 11;
    t = 0.;
    tout = .01;
    i__1 = nout;
    for (iout = 1; iout <= i__1; ++iout) {

L45:
	ddaskr_((U_fp)resh_, &neq, &t, u, uprime, &tout, info, &rtol, &atol, &
		idid, rwork, &lrw, iwork, &liw, rpar, ipar, (U_fp)djacilu_, (
		U_fp)dpsolilu_, (U_fp)rtheat_, &nrt, jroot);

	umax = 0.;
	i__2 = neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L50: */
/* Computing MAX */
	    d__2 = umax, d__3 = (d__1 = u[i__ - 1], abs(d__1));
	    umax = max(d__2,d__3);
	}

	hu = rwork[6];
	nqu = iwork[7];
	nst = iwork[10];
	nni = iwork[18];
	nli = iwork[19];
	io___41.ciunit = lout;
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&umax, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&nqu, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&hu, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&nst, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nni, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nli, (ftnlen)sizeof(integer));
	e_wsfe();

	if (idid == 5) {
	    s_wsfe(&io___42);
	    do_fio(&c__1, (char *)&jroot[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jroot[1], (ftnlen)sizeof(integer));
	    e_wsfe();
	    goto L45;
	}

	if (idid < 0) {
	    io___43.ciunit = lout;
	    s_wsfe(&io___43);
	    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    goto L80;
	}

	tout *= 2.;
/* L70: */
    }

/* Here we display some final statistics for the problem. */
/* The ratio of NLI to NNI is the average dimension of the Krylov */
/* subspace involved in the Krylov linear iterative method. */
L80:
    nst = iwork[10];
    npe = iwork[12];
    nre = iwork[11] + npe * mband;
    liw = iwork[16];
    lrw = iwork[17];
    nni = iwork[18];
    nli = iwork[19];
    nps = iwork[20];
    if (nni != 0) {
	avdim = (real) nli / (real) nni;
    }
    ncfn = iwork[14];
    ncfl = iwork[15];
    nrte = iwork[35];
    io___52.ciunit = lout;
    s_wsfe(&io___52);
    do_fio(&c__1, (char *)&lrw, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&liw, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nst, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nre, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ipar[29], (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nrte, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npe, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nps, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nni, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nli, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&avdim, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ncfn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ncfl, (ftnlen)sizeof(integer));
    e_wsfe();
    io___53.ciunit = lout;
    s_wsfe(&io___53);
    do_fio(&c__1, (char *)&lwpmin, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&liwpmin, (ftnlen)sizeof(integer));
    e_wsfe();

/* ------  End of main program for DHEATILU example program -------------- */
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Subroutine */ int uinit_(doublereal *u, doublereal *uprime, doublereal *
	rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal dx, xj, yk;
    static integer neq, ioff;


/* This routine computes and loads the vector of initial values. */
/* The initial U values are given by the polynomial u = 16x(1-x)y(1-y). */
/* The initial UPRIME values are set to zero.  (DDASKR corrects these */
/* during the first time step.) */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --uprime;
    --u;

    /* Function Body */
    neq = ipar[33];
    m = ipar[34];
    dx = rpar[3];

    i__1 = m + 1;
    for (k = 0; k <= i__1; ++k) {
	yk = k * dx;
	ioff = (m + 2) * k;
	i__2 = m + 1;
	for (j = 0; j <= i__2; ++j) {
	    xj = j * dx;
	    i__ = ioff + j + 1;
	    u[i__] = xj * 16. * (1. - xj) * yk * (1. - yk);
/* L10: */
	}
/* L20: */
    }
    i__1 = neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	uprime[i__] = 0.;
    }
    return 0;
/* ------------  End of Subroutine UINIT  -------------------------------- */
} /* uinit_ */

/* Subroutine */ int resh_(doublereal *t, doublereal *u, doublereal *uprime, 
	doublereal *cj, doublereal *delta, integer *ires, doublereal *rpar, 
	integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m, m2, neq, ioff;
    static doublereal temx, temy, coeff;


/* This is the user-supplied RES subroutine for this example. */
/* It computes the residuals for the 2-D discretized heat equation, */
/* with zero boundary values. */


/* Set problem constants using IPAR and RPAR. */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --delta;
    --uprime;
    --u;

    /* Function Body */
    neq = ipar[33];
    m = ipar[34];
    coeff = rpar[4];
    m2 = m + 2;

/* Load U into DELTA, in order to set boundary values. */
    i__1 = neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	delta[i__] = u[i__];
    }

/* Loop over interior points, and load residual values. */
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	ioff = m2 * k;
	i__2 = m;
	for (j = 1; j <= i__2; ++j) {
	    i__ = ioff + j + 1;
	    temx = u[i__ - 1] + u[i__ + 1];
	    temy = u[i__ - m2] + u[i__ + m2];
	    delta[i__] = uprime[i__] - (temx + temy - u[i__] * 4.) * coeff;
/* L20: */
	}
/* L30: */
    }

    return 0;
/* ------------  End of Subroutine RESH  --------------------------------- */
} /* resh_ */

/* Subroutine */ int rtheat_(integer *neq, doublereal *t, doublereal *u, 
	doublereal *up, integer *nrt, doublereal *rval, doublereal *rpar, 
	integer *ipar)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal umax;


/* This routine finds the max of U, and sets RVAL(1) = max(u) - 0.1, */
/* RVAL(2) = max(u) - 0.01. */


    /* Parameter adjustments */
    --up;
    --u;
    --rval;

    /* Function Body */
    umax = 0.;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing MAX */
	d__1 = umax, d__2 = u[i__];
	umax = max(d__1,d__2);
    }
    rval[1] = umax - .1;
    rval[2] = umax - .01;

    return 0;
/* ------------  End of Subroutine RTHEAT  ------------------------------- */
} /* rtheat_ */

