/* dkrdem.f -- translated by f2c (version 20100827).
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
    integer neq;
} local_;

#define local_1 local_

/* Table of constant values */

static integer c__1 = 1;

/* DKRDEM program: DDASKR demonstration program */
/* ----------------------------------------------------------------------- */

/* ***BEGIN PROLOGUE  DKRDEM */
/* ***DATE WRITTEN   020813     (YYMMDD) */
/* ***REVISION DATE  021217   Added JROOT output value in Problem 2. */
/* ***AUTHORS  Linda R. Petzold and Alan C. Hindmarsh */
/*             LAWRENCE LIVERMORE NATIONAL LABORATORY */
/*             LIVERMORE, CA    94550 */
/* ***PURPOSE  Quick check for routine DDASKR. */
/* ***DESCRIPTION */
/*       Demonstration program for DDASKR. */
/*       This version is in double precision. */

/*       DDASKR is used to solve two simple problems, */
/*       one nonstiff and one intermittently stiff. */
/*       If the errors are too large, or other difficulty occurs, */
/*       a warning message is printed.  All output is on unit LUN. */

/*       To run the demonstration problems with full printing, */
/*       set KPRINT = 3. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DDASKR,RES1,RT1,RES2,JAC2,RT2 */
/* ***END PROLOGUE DKRDEM */

/* ----------------------------------------------------------------------- */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_110[] = "(\002 DKRDEM: Demonstration Program for DDASK"
	    "R\002///\002 Problem 1..\002//\002 Problem is  dY/dT = ((2*LOG(Y"
	    ")+8)/T - 5)*Y,  Y(1) = 1\002/\002 Solution is  Y(T) = EXP(-T**2 "
	    "+ 5*T - 4)\002/\002 Root functions are..\002/\002 R1 = dY/dT  (r"
	    "oot at T = 2.5)\002/\002 R2 = LOG(Y) - 2.2491  (roots at T = 2.4"
	    "7 and T = 2.53)\002/\002 RTOL =\002,e10.1,\002   ATOL =\002,e10."
	    "1,\002   JTYPE =\002,i3/)";
    static char fmt_130[] = "(\002 At t =\002,e15.7,5x,\002y =\002,e15.7,5x"
	    ",\002error =\002,e12.4)";
    static char fmt_135[] = "(//\002 WARNING.. Error exceeds 1000 * tolera"
	    "nce\002//)";
    static char fmt_150[] = "(/\002      Root found at t =\002,e15.7,5x,\002"
	    "JROOT =\002,2i5)";
    static char fmt_160[] = "(\002      Error in t location of root is\002,e"
	    "12.4/)";
    static char fmt_165[] = "(//\002 WARNING.. Root error exceeds 1.0D-3\002"
	    "//)";
    static char fmt_190[] = "(/\002 Final statistics for this run..\002/\002"
	    " number of steps =\002,i5/\002 number of Gs    =\002,i5/\002 (ex"
	    "cluding Js)  =\002,i5/\002 number of Js    =\002,i5/\002 number "
	    "of Rs    =\002,i5/\002 error overrun   =\002,e10.2)";
    static char fmt_200[] = "(/80(\002-\002)//\002 Problem 2.. Van Der Pol o"
	    "scillator\002//\002 Problem is dY1/dT = Y2,  dY2/dT = 100*(1-Y1*"
	    "*2)*Y2 - Y1\002/\002            Y1(0) = 2,  Y2(0) = 0\002/\002 R"
	    "oot function is  R(T,Y,YP) = Y1\002/\002 RTOL =\002,e10.1,\002  "
	    " ATOL =\002,2e10.1)";
    static char fmt_210[] = "(/80(\002.\002)//\002 Solution with JTYPE =\002"
	    ",i2/)";
    static char fmt_230[] = "(\002 At t =\002,e15.7,5x,\002y1 =\002,e15.7,"
	    "5x,\002y2 =\002,e15.7)";
    static char fmt_240[] = "(/\002      Root found at t =\002,e15.7,\002  J"
	    "ROOT =\002,i3)";
    static char fmt_250[] = "(\002      Error in t location of root is\002,e"
	    "12.4/)";
    static char fmt_255[] = "(//\002 WARNING.. Root error exceeds 1.0\002//)";
    static char fmt_280[] = "(/\002 Final statistics for this run..\002/\002"
	    " number of steps =\002,i5/\002 number of Gs    =\002,i5/\002 (ex"
	    "cluding Js)  =\002,i5/\002 number of Js    =\002,i5/\002 number "
	    "of Rs    =\002,i5)";
    static char fmt_300[] = "(/80(\002-\002)//\002 Number of errors encounte"
	    "red =\002,i3)";
    static char fmt_700[] = "(/,\002 **********DDASKR passed all tests******"
	    "****\002)";
    static char fmt_800[] = "(/,\002 **********DDASKR failed some tests*****"
	    "****\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double exp(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal t, y[2], er, yt;
    extern /* Subroutine */ int rt1_(), rt2_();
    static integer nje, nre;
    static doublereal ero;
    static integer liw, lun, nrt, lrw, nst;
    extern /* Subroutine */ int jac2_(), res1_(), res2_();
    static integer idid, nrea, info[20], ipar;
    static doublereal atol[2];
    static integer jdum;
    static doublereal rpar;
    static integer nerr, nrte;
    static doublereal errt;
    static integer iout;
    static doublereal rtol[2], tout;
    static integer leniw, ipass, lenrw;
    static doublereal psdum;
    static integer iwork[100], jtype, jroot[2], kroot;
    static doublereal tzero, rwork[100];
    extern /* Subroutine */ int ddaskr_();
    static doublereal yprime[2];
    static integer kprint;

    /* Fortran I/O blocks */
    static cilist io___18 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_135, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_150, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_160, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_165, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_190, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_200, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_210, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___50 = { 0, 0, 0, fmt_250, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_255, 0 };
    static cilist io___52 = { 0, 0, 0, fmt_280, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_300, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___55 = { 0, 0, 0, fmt_800, 0 };



    lun = 6;
    kprint = 3;
    ipass = 1;
    nerr = 0;
/* ----------------------------------------------------------------------- */
/* First problem. */
/* The initial value problem is.. */
/*   dY/dT = ((2*LOG(Y) + 8)/T - 5)*Y,  Y(1) = 1,  1 .LE. T .LE. 6 */
/* The solution is  Y(T) = EXP(-T**2 + 5*T - 4), YPRIME(1) = 3 */
/* The two root functions are.. */
/*   R1(T,Y,Y') = dY/dT  (with root at T = 2.5), */
/*   R2(T,Y,Y') = LOG(Y) - 2.2491  (with roots at T = 2.47 and 2.53) */
/* ----------------------------------------------------------------------- */
/* Set all input parameters and print heading. */
    for (i__ = 1; i__ <= 20; ++i__) {
/* L10: */
	info[i__ - 1] = 0;
    }
    local_1.neq = 1;
    t = 1.;
    y[0] = 1.;
    tout = 2.;
    rtol[0] = 0.;
    atol[0] = 1e-6;
    lrw = 100;
    liw = 100;
    idid = 0;

/* Set INFO(11) = 1 if DDASKR is to compute the initial YPRIME, and */
/* generate an initial guess for YPRIME.  Otherwise, set INFO(11) = 0 */
/* and supply the correct initial value for YPRIME. */

    info[10] = 0;
    yprime[0] = 3.;

/* Note: JTYPE indicates the Jacobian type: */
/*       JTYPE = 1 ==> Jacobian is dense and user-supplied */
/*       JTYPE = 2 ==> Jacobian is dense and computed internally */

    jtype = 2;
    info[4] = 2 - jtype;
    nrt = 2;
    if (kprint >= 2) {
	io___18.ciunit = lun;
	s_wsfe(&io___18);
	do_fio(&c__1, (char *)&rtol[0], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&atol[0], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&jtype, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/* Call DDASKR in loop over TOUT values = 2, 3, 4, 5, 6. */
    ero = 0.;
    for (iout = 1; iout <= 5; ++iout) {

L120:
	ddaskr_((U_fp)res1_, &local_1.neq, &t, y, yprime, &tout, info, rtol, 
		atol, &idid, rwork, &lrw, iwork, &liw, &rpar, &ipar, &jdum, &
		psdum, (U_fp)rt1_, &nrt, jroot);

/* Print Y and error in Y, and print warning if error too large. */
	yt = exp(-t * t + t * 5. - 4.);
	er = y[0] - yt;
	if (kprint > 2) {
	    io___30.ciunit = lun;
	    s_wsfe(&io___30);
	    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&y[0], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&er, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (idid < 0) {
	    goto L185;
	}
	er = abs(er) / atol[0];
	ero = max(ero,er);
	if (er < 1e3) {
	    goto L140;
	}
	if (kprint >= 2) {
	    ipass = 0;
	    io___31.ciunit = lun;
	    s_wsfe(&io___31);
	    e_wsfe();
	}
	++nerr;
L140:
	if (idid != 5) {
	    goto L175;
	}

/* If a root was found, write results and check root location. */
/* Then return to DDASKR to continue the integration. */
	if (kprint > 2) {
	    io___32.ciunit = lun;
	    s_wsfe(&io___32);
	    do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&jroot[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jroot[1], (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	if (jroot[0] != 0) {
	    errt = t - 2.5;
	}
	if (jroot[1] != 0 && t <= 2.5) {
	    errt = t - 2.47;
	}
	if (jroot[1] != 0 && t > 2.5) {
	    errt = t - 2.53;
	}
	if (kprint > 2) {
	    io___34.ciunit = lun;
	    s_wsfe(&io___34);
	    do_fio(&c__1, (char *)&errt, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (abs(errt) < .001) {
	    goto L170;
	}
	if (kprint >= 2) {
	    ipass = 0;
	    io___35.ciunit = lun;
	    s_wsfe(&io___35);
	    e_wsfe();
	}
	++nerr;
L170:
	goto L120;

/* If no root found, increment TOUT and loop back. */
L175:
	tout += 1.;
/* L180: */
    }

/* Problem complete.  Print final statistics. */
L185:
    if (idid < 0) {
	++nerr;
    }
    nst = iwork[10];
    nre = iwork[11];
    nje = iwork[12];
    nrte = iwork[35];
    lenrw = 0;
    leniw = 0;
    nrea = nre;
    if (jtype == 2) {
	nre += local_1.neq * nje;
    }

    if (kprint > 2) {
	io___43.ciunit = lun;
	s_wsfe(&io___43);
	do_fio(&c__1, (char *)&nst, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nre, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrea, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nje, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nrte, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ero, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* ----------------------------------------------------------------------- */
/* Second problem (Van Der Pol oscillator). */
/* The initial value problem is.. */
/*   dY1/dT = Y2,  dY2/dT = 100*(1 - Y1**2)*Y2 - Y1, */
/*   Y1(0) = 2,  Y2(0) = 0,  0 .LE. T .LE. 200 */
/*   Y1PRIME(0) = 0, Y2PRIME(0) = -2 */
/* The root function is  R(t,Y,Y') = Y1. */
/* An analytic solution is not known, but the zeros of Y1 are known */
/* to 15 figures for purposes of checking the accuracy. */
/* ----------------------------------------------------------------------- */

/* Reset INFO array */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L195: */
	info[i__ - 1] = 0;
    }

/* Set tolerance parameters and print heading. */
/* Note that INFO(2) is set to 1, indicating that RTOL and ATOL */
/* are arrays.  Each entry of RTOL and ATOL must then be defined. */

    info[1] = 1;
    rtol[0] = 1e-6;
    rtol[1] = 1e-6;
    atol[0] = 1e-6;
    atol[1] = 1e-4;
    if (kprint >= 2) {
	io___44.ciunit = lun;
	s_wsfe(&io___44);
	do_fio(&c__1, (char *)&rtol[0], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&atol[0], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&atol[1], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* Note: JTYPE indicates the Jacobian type: */
/*       JTYPE = 1 ==> Jacobian is dense and user-supplied */
/*       JTYPE = 2 ==> Jacobian is dense and computed internally */

/* Loop over JTYPE = 1, 2.  Set remaining parameters and print JTYPE. */
    for (jtype = 1; jtype <= 2; ++jtype) {

/*     Set INFO(1) = 0 to indicate start of a new problem */
/*     Set INFO(5) = 2-JTYPE to tell DDASKR the Jacobian type. */

	info[0] = 0;
	info[4] = 2 - jtype;
	local_1.neq = 2;
	t = 0.;
	y[0] = 2.;
	y[1] = 0.;
	yprime[0] = 0.;
	yprime[1] = -2.;
	tout = 20.;
	nrt = 1;
	if (kprint > 2) {
	    io___45.ciunit = lun;
	    s_wsfe(&io___45);
	    do_fio(&c__1, (char *)&jtype, (ftnlen)sizeof(integer));
	    e_wsfe();
	}

/* Call DDASKR in loop over TOUT values = 20, 40, ..., 200. */
	for (iout = 1; iout <= 10; ++iout) {

L220:
	    ddaskr_((U_fp)res2_, &local_1.neq, &t, y, yprime, &tout, info, 
		    rtol, atol, &idid, rwork, &lrw, iwork, &liw, &rpar, &ipar,
		     (U_fp)jac2_, &psdum, (U_fp)rt2_, &nrt, jroot);

/* Print Y1 and Y2. */
	    if (kprint > 2) {
		io___46.ciunit = lun;
		s_wsfe(&io___46);
		do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&y[0], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&y[1], (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (idid < 0) {
		goto L275;
	    }
	    if (idid != 5) {
		goto L265;
	    }

/* If a root was found, write results and check root location. */
/* Then return to DDASKR to continue the integration. */
	    if (kprint > 2) {
		io___47.ciunit = lun;
		s_wsfe(&io___47);
		do_fio(&c__1, (char *)&t, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&jroot[0], (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    kroot = (integer) (t / 81.2 + .5);
	    tzero = (doublereal) (kroot - 1) * 81.41853556212 + 
		    81.17237787055;
	    errt = t - tzero;
	    if (kprint > 2) {
		io___50.ciunit = lun;
		s_wsfe(&io___50);
		do_fio(&c__1, (char *)&errt, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (errt < 1.) {
		goto L260;
	    }
	    if (kprint >= 2) {
		ipass = 0;
		io___51.ciunit = lun;
		s_wsfe(&io___51);
		e_wsfe();
	    }
	    ++nerr;
L260:
	    goto L220;

/* If no root found, increment TOUT and loop back. */
L265:
	    tout += 20.;
/* L270: */
	}

/* Problem complete.  Print final statistics. */
L275:
	if (idid < 0) {
	    ++nerr;
	}
	nst = iwork[10];
	nre = iwork[11];
	nje = iwork[12];
	nrte = iwork[35];
	lenrw = 0;
	leniw = 0;
	nrea = nre;
	if (jtype == 2) {
	    nre += local_1.neq * nje;
	}
	if (kprint >= 2) {
	    io___52.ciunit = lun;
	    s_wsfe(&io___52);
	    do_fio(&c__1, (char *)&nst, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nre, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nrea, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nje, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nrte, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
/* L290: */
    }


    if (kprint >= 2) {
	io___53.ciunit = lun;
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&nerr, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    if (nerr > 0) {
	ipass = 0;
    }
    if (ipass == 1 && kprint > 1) {
	io___54.ciunit = lun;
	s_wsfe(&io___54);
	e_wsfe();
    }
    if (ipass == 0 && kprint != 0) {
	io___55.ciunit = lun;
	s_wsfe(&io___55);
	e_wsfe();
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */


/* Subroutine */ int res1_(doublereal *t, doublereal *y, doublereal *yprime, 
	doublereal *cj, doublereal *delta, integer *ires, doublereal *rpar, 
	integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int f1_(integer *, doublereal *, doublereal *, 
	    doublereal *);


/*     Check Y to make sure that it is valid input. */
/*     If Y is less than or equal to zero, this is invalid input. */

    /* Parameter adjustments */
    --delta;
    --yprime;
    --y;

    /* Function Body */
    if (y[1] <= 0.) {
	*ires = -1;
	return 0;
    } else {

/*        Call F1 to obtain F(T,Y) */

	f1_(&local_1.neq, t, &y[1], &delta[1]);

/*        Form G = Y' - F(T,Y) */

	i__1 = local_1.neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    delta[i__] = yprime[i__] - delta[i__];
/* L10: */
	}
    }

    return 0;
} /* res1_ */


/* Subroutine */ int f1_(integer *neq, doublereal *t, doublereal *y, 
	doublereal *ydot)
{
    /* Builtin functions */
    double log(doublereal);

    /* Parameter adjustments */
    --ydot;
    --y;

    /* Function Body */
    ydot[1] = ((log(y[1]) * 2. + 8.) / *t - 5.) * y[1];
    return 0;
} /* f1_ */


/* Subroutine */ int rt1_(integer *neq, doublereal *t, doublereal *y, 
	doublereal *yp, integer *nrt, doublereal *rval, doublereal *rpar, 
	integer *ipar)
{
    /* Builtin functions */
    double log(doublereal);

    /* Parameter adjustments */
    --rval;

    /* Function Body */
    rval[1] = *yp;
    rval[2] = log(*y) - 2.2491;
    return 0;
} /* rt1_ */


/* Subroutine */ int res2_(doublereal *t, doublereal *y, doublereal *yprime, 
	doublereal *cj, doublereal *delta, integer *ires, doublereal *rpar, 
	integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int f2_(integer *, doublereal *, doublereal *, 
	    doublereal *);


/*     CALL F2 TO OBTAIN F(T,Y) */

    /* Parameter adjustments */
    --delta;
    --yprime;
    --y;

    /* Function Body */
    f2_(&local_1.neq, t, &y[1], &delta[1]);

/*     FORM G = Y' - F(T,Y) */

    i__1 = local_1.neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	delta[i__] = yprime[i__] - delta[i__];
    }

    return 0;
} /* res2_ */


/* Subroutine */ int f2_(integer *neq, doublereal *t, doublereal *y, 
	doublereal *ydot)
{
    /* Parameter adjustments */
    --ydot;
    --y;

    /* Function Body */
    ydot[1] = y[2];
    ydot[2] = (1. - y[1] * y[1]) * 100. * y[2] - y[1];
    return 0;
} /* f2_ */


/* Subroutine */ int jac2_(doublereal *t, doublereal *y, doublereal *yprime, 
	doublereal *pd, doublereal *cj, doublereal *rpar, integer *ipar)
{

/* First define the Jacobian matrix for the right-hand side */
/* of the ODE: Y' = F(T,Y) , i.e. dF/dY. */

    /* Parameter adjustments */
    pd -= 3;
    --y;

    /* Function Body */
    pd[3] = 0.;
    pd[5] = 1.;
    pd[4] = y[1] * -200. * y[2] - 1.;
    pd[6] = (1. - y[1] * y[1]) * 100.;

/* Next update the Jacobian with the right-hand side to form the */
/* DAE Jacobian: CJ*dR/dY' + dR/dY = CJ*I - dF/dY. */

    pd[3] = *cj - pd[3];
    pd[5] = -pd[5];
    pd[4] = -pd[4];
    pd[6] = *cj - pd[6];

    return 0;
} /* jac2_ */


/* Subroutine */ int rt2_(integer *neq, doublereal *t, doublereal *y, 
	doublereal *yp, integer *nrt, doublereal *rval, doublereal *rpar, 
	integer *ipar)
{
    /* Parameter adjustments */
    --rval;
    --yp;
    --y;

    /* Function Body */
    rval[1] = y[1];
    return 0;
} /* rt2_ */

