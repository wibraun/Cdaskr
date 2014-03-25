/* dilupre.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__3 = 3;
static integer c__80 = 80;
static real c_b21 = 0.f;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* ----------------------------------------------------------------------- */

/*             Preconditioner Routines for Sparse Problems */
/*                          13 December 2000 */

/* The following triple of subroutines -- DSPSETUP, DJACILU and */
/* DPSOLILU -- provides a general-purpose sparse incomplete LU (ILU) */
/* preconditioner matrix for use with the DDASPK solver, with the Krylov */
/* linear system method.  When using DDASPK to solve a problem */
/* G(t,y,y') = 0, whose iteration matrix (Jacobian) */
/*    J = dG/dy + c * dG/dy'  (c = scalar) */
/* is a general sparse matrix, these routines can be used to generate */
/* an approximation to J as the preconditioner and to solve the resulting */
/* sparse linear system, in conjunction with the Krylov method option */
/* (INFO(12) = 1) in DDASPK. */

/* The incomplete LU factorization is achieved via one of two routines */
/* from the SPARSKIT library available from Yousef Saad at the */
/* University of Minnesota.  The two routines are ILUT and ILUTP. */
/* See below for detailed descriptions of these routines. */

/* Descriptions of the above routines are as follows: */

/* DSPSETUP   - Setup routine for the preconditioner.  This routine */
/*              checks the user input for illegal values, and calculates */
/*              the minimum length needed for preconditioner workspace */
/*              arrays in DDASPK. */

/* DJACILU   -  This routine is a version of JAC for DDASPK. */
/*              It uses finite-differences to calculate the Jacobian */
/*              matrix J in sparse format, and then performs an */
/*              incomplete LU factorization using either ILUT or ILUTP. */
/*              DJACILU must be declared EXTERNAL in the user's main */
/*              program and passed through to DDASPK. */

/* DPSOLILU   - This routine is a version of PSOL for DDASPK. */
/*              It uses the incomplete factorization calculated in */
/*              DJACILU.  DPSOLILU must be declared EXTERNAL in the */
/*              user's main program and passed through to DDASPK. */


/* The routines ILUT and ILUTP are part of the SPARSKIT library, */
/* and are contained in the file 'dsparsk.f'.  ILUT performs an */
/* incomplete LU factorization of a sparse matrix using a dual */
/* thresholding technique based on a drop tolerance (TOLILUT below) */
/* and a level of fill-in parameter (LFILILUT).  LFILILUT controls */
/* the amount of fill-in allowed in the factorization (limited to a */
/* maximum of 2*LFILILUT*NEQ, but normally much less).  Increasing */
/* LFILILUT will generally make the ILU factorization more accurate. */
/* TOLILUT also controls the accuracy of the ILU factorization via */
/* a drop tolerance based on element size.  Descreasing TOLILUT */
/* will increase the amount of fill-in and make for a more accurate */
/* factorization.  ILUTP is a variant of ILUT that in addition performs */
/* pivoting based on a tolerance ratio PERMTOL. */

/* An important aspect of using incomplete factorization techniques */
/* is that of reordering the rows and columns in the Jacobian matrix */
/* J before performing the ILU.  In this package, this is accomplished */
/* via the parameter IREORDER, which when equal to 1 performs */
/* a reverse Cuthill-McKee (RCM) reordering before performing the */
/* ILU factorization.  Based on the limited amount of testing done so */
/* far, RCM seems the best overall choice.  It is possible to include */
/* a different reordering technique if desired. */

/* To use these routines in conjunction with DDASPK, the user's calling */
/* program should include the following, in addition to setting the other */
/* DDASPK input parameters. */

/* (a) Dimension the array IPAR to have length at least 30, and load the */
/*     following parameters into IPAR as */

/*      IPAR(1)  = ML        - The lower bandwidth used in calculating */
/*                             the approximate Jacobian matrix. */
/*      IPAR(2)  = MU        - The upper bandwidth used in calculating */
/*                             the approximate Jacobian matrix. */
/*      IPAR(3)  = LENPFAC   - The average number of nonzeros in a */
/*                             row of the Jacobian matrix.  The */
/*                             maximum of nonzeros allowed in the */
/*                             Jacobian matrix is NNZMX = LENPFAC*NEQ. */
/*                             LENPFAC must be .GE. 2. */
/*      IPAR(4)  = LENPLUFAC - The average amount of fill-in per row */
/*                             in the factored Jacobian matrix.  The */
/*                             maximum number of nonzeros allowed */
/*                             in the factored Jacobian matrix is */
/*                             LENPLUMX = NNZMX + LENPLUFAC*NEQ. */
/*                             LENPLUFAC must be .GE. 2. */
/*      IPAR(5)  = IPREMETH  - Preconditioner type flag. */
/*                             =1 means ILUT preconditioner used */
/*                             =2 means ILUTP preconditioner used */
/*      IPAR(6)  = LFILILUT  - Fill-in parameter for ILUT and ILUTP. */
/*                             The largest LFILILUT elements per row */
/*                             of the L and U factors are kept.  Each */
/*                             row of L and U will have a maximum of */
/*                             LFILILUT elements in addition to */
/*                             their original number of nonzero */
/*                             elements. */
/*      IPAR(7)  = IREORDER  - Reordering flag. */
/*                             =0 means that no reordering of the */
/*                                rows and columns of the Jacobian */
/*                                matrix is performed before the */
/*                                incomplete factorization is performed. */
/*                             =1 means that a reverse Cuthill-McKee */
/*                                (RCM) reordering is performed. */
/*      IPAR(8)  = ISRNORM   - Row norm flag. */
/*                             =1 means that row norms of the Jacobian */
/*                                matrix are computed and used as */
/*                                row scalings when solving the */
/*                                preconditioner linear system P*x=b. */
/*                             =0 means no row norm scaling is used. */
/*      IPAR(9)  = NORMTYPE  - Type of row norm scaling for ISRNORM */
/*                             =0 means max-norm is used. */
/*                             =1 means 1-norm is used. */
/*                             =2 means 2-norm is used. */
/*      IPAR(10) = JACOUT    - Output Jacobian flag. */
/*                             =1 means write the calculated Jacobian */
/*                                matrix along with the initial value of */
/*                                the residual G to a file pointed to by */
/*                                the logical unit number in IPAR(29). */
/*                                This is done after any reordering and */
/*                                scaling is performed.  The user must */
/*                                have attached the unit number to a */
/*                                file before calling DDASPK.  The */
/*                                integration is then halted by setting */
/*                                IRES = -2 (and a false message about */
/*                                failure to initialize is printed). */
/*                             =0 means no Jacobian matrix output. */
/*                             The matrix and initial residual G are */
/*                             written in Boeing-Harwell format. */
/*      IPAR(11) = JSCALCOL  - Flag for scaling columns of the */
/*                             Jacobian matrix by the inverses of the */
/*                             elements in the EWT array. */
/*                             =0 means no scaling. */
/*                             =1 means perform scaling. */

/*      IPAR(21:28)          - Used to hold pointer information. */
/*      IPAR(29)             - Logical unit number to write matrix output */
/*                             file on.  Only needed when JACOUT = 1. */
/*      IPAR(30)             - On return from DDASPK, this value */
/*                             holds the number of calls to the */
/*                             RES routine used in the preconditioner */
/*                             evaluations. */

/* (b) Dimension the array RPAR to have length at least 2, and load the */
/*     following parameters into RPAR as */

/*      RPAR(1)  = TOLILUT   - Drop tolerance for use by ILUT and ILUTP. */
/*                             TOLILUT must be .ge. 0.  Larger values */
/*                             of TOLILUT cause less fill-in.  Good */
/*                             values range from 0.001 to 0.01. */
/*      RPAR(2)  = PERMTOL   - Tolerance ratio used in determining column */
/*                             pivoting by ILUTP.  PERMTOL must be */
/*                             .ge. 0.  Good values are from 0.1 to */
/*                             0.01.  Two columns are permuted only if */
/*                             A(i,j)*PERMTOL .GT. A(i,i). */

/*     The two paramaters TOLILUT and LFILILUT gives the user a great */
/*     deal of flexibility: one can use TOLILUT=0 to get a strategy */
/*     based on keeping the largest elements in each row of L and U. */
/*     Taking TOLILUT .NE. 0 but LFILILUT=NEQ will give the usual */
/*     threshold strategy (however, fill-in is then unpredictable). */

/* (c) Include the names DJACILU and DPSOLILU in an EXTERNAL statement. */
/*     Set INFO(15) = 1 to indicate that a JAC routine exists. */
/*     Then in the call to DDASPK, pass the names DJACILU and DPSOLILU */
/*     as the arguments JAC and PSOL, respectively. */

/* (d) The DDASPK work arrays RWORK and IWORK must include segments WP */
/*     and IWP for use by DJACILU/DPSOLILU.  The lengths of these depend */
/*     on the problem size, half-bandwidths, and other parameters */
/*     as follows: */
/*       LWP =  length of RWORK segment WP = */
/*              2*LENPFAC*NEQ + LENPLUFAC*NEQ + ISRNORM*NEQ + NEQ */
/*       LIWP = length of IWORK segment IWP = */
/*              3*NEQ+1 + 3*LENPFAC*NEQ + 2*LENPLUFAC*NEQ */
/*                 + 2*IREORDER*NEQ + 2*(IPREMETH-1)*NEQ */
/*     Load these lengths in IWORK as */
/*       IWORK(27) = LWP */
/*       IWORK(28) = LIWP */
/*     and include these values in the declared size of RWORK and IWORK. */


/* The DJACILU and DPSOLILU routines generate and solve the sparse */
/* preconditioner matrix P within the preconditioned Krylov algorithm */
/* used by DDASPK when INFO(12) = 1.  P is generated and ILU-factored */
/* periodically during the integration, and the factors are used to */
/* solve systems Px = b as needed. */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int dspsetup_(integer *neq, integer *lwp, integer *liwp, 
	doublereal *rpar, integer *ipar, integer *ierr, integer *lwp_min__, 
	integer *liwp_min__)
{
    static integer jscalcol, ireorder, ipremeth, lfililut, lenplumx, normtype,
	     lrownrms, lenplufac, lbw, lju, ubw, ljac, ljlu, lplu, neqp1, 
	    liwk1, lrwk1, ljaci, ljacj, lmask, lperm, nnzmx, jacout, lqperm, 
	    lenpfac, llevels;
    static doublereal permtol;
    static integer isrnorm;
    static doublereal tolilut;

/* ... Version of 12-12-00 */
/* ... Calculate storage needed for ILU decomposition of the Jacobian */
/*     matrix for use as a preconditioner by the DDASPK solver. */
/*     Also, check for illegal input. */
/* ... Input arguments: */
/* total number of equations */
/* user real workspace */
/* user integer workspace */
/* current length of WP for DDASPK */
/* ... Output arguments: */
/* current length of IWP for DDASPK */
/* IERR between 1 and 11, means there's an */
/* illegal value for IPAR(IERR). */
/* IERR = 12 means IPAR(29) is illegal. */
/* IERR = 21 means RPAR(1) is illegal. */
/* IERR = 22 means RPAR(2) is illegal. */
/* IERR = 30 means more WP length is needed. */
/* IERR = 31 means more IWP length is needed. */
/* error flag (0 means success, else failure) */
/* minimum WP length needed. */
/* ... Local variables: */
/* minimum IWP length needed. */
/* =1 causes row normalization of JAC. */
/* 2-norm row scaling */
/* =0,1,2 for max-norm, 1-norm, or */
/* =1 causes row and column reordering of JAC. */
/*   be written to a file and then exit with */
/*   ierr = 1 to signal a stop to DDASPK. */
/* =1 causes the Jacobian matrix and SAVR to */
/*   to be scaled by EWT-inverse */
/* =1 causes the columns of the Jacobian matrix */
/* ... Load values from IPAR and RPAR.  Check for illegal values. */
    /* Parameter adjustments */
    --ipar;
    --rpar;

    /* Function Body */
    lbw = ipar[1];
/* LBW must be .gt. 0 */
    if (lbw <= 0) {
	*ierr = 1;
	return 0;
    }
    ubw = ipar[2];
/* UBW must be .gt. 0 */
    if (ubw <= 0) {
	*ierr = 2;
	return 0;
    }
    lenpfac = ipar[3];
/* LENPFAC must be .ge. 2 */
    if (lenpfac <= 1) {
	*ierr = 3;
	return 0;
    }
    lenplufac = ipar[4];
/* LENPLUFAC must be .ge. 2 */
    if (lenplufac <= 1) {
	*ierr = 4;
	return 0;
    }
    ipremeth = ipar[5];
/* IPREMETH must be .eq. 1 or 2 currently */
    if (ipremeth != 1 && ipremeth != 2) {
	*ierr = 5;
	return 0;
    }
    lfililut = ipar[6];
/* LFILILUT must be .ge. 0 */
    if (lfililut < 0) {
	*ierr = 6;
	return 0;
    }
    ireorder = ipar[7];
/* IREORDER must be 0 or 1 */
    if (ireorder < 0 || ireorder > 1) {
	*ierr = 7;
	return 0;
    }
    isrnorm = ipar[8];
/* ISRNORM must be 0 or 1 */
    if (isrnorm < 0 || isrnorm > 1) {
	*ierr = 8;
	return 0;
    }
    normtype = ipar[9];
/* NORMTYPE must be 0, 1, or 2 */
    if (normtype < 0 || normtype > 2) {
	*ierr = 9;
	return 0;
    }
    jacout = ipar[10];
/* JACOUT must be 0 or 1 */
    if (jacout < 0 || jacout > 1) {
	*ierr = 10;
	return 0;
    }
    jscalcol = ipar[11];
/* JSCALCOL must be 0 or 1 */
    if (jscalcol < 0 || jscalcol > 1) {
	*ierr = 11;
	return 0;
    }
    if (jacout == 1) {
/* IPAR(29) must be .gt. 0 */
	if (ipar[29] <= 0) {
	    *ierr = 12;
	    return 0;
	}
    }
    tolilut = rpar[1];
/* TOLILUT must be .ge. 0.0 */
    if (tolilut < 0.f) {
	*ierr = 21;
	return 0;
    }
    if (ipremeth == 2) {
	permtol = rpar[2];
/* PERMTOL must be .ge. 0.0 */
	if (permtol < 0.f) {
	    *ierr = 22;
	    return 0;
	}
    }
/* ... Calculate minimum work lengths for WP and IWP arrays. */
    neqp1 = *neq + 1;
    nnzmx = lenpfac * *neq;
    lenplumx = nnzmx + lenplufac * *neq;
/* ... Set up pointers into WP */
    ljac = 1;
    lrownrms = nnzmx + ljac;
    if (isrnorm == 1) {
	lplu = lrownrms + *neq;
    } else {
	lplu = lrownrms;
    }
    lrwk1 = lplu + lenplumx;
    *lwp_min__ = lrwk1 + *neq - 1;
    if (*lwp < *lwp_min__) {
	*ierr = 30;
/* more WP length needed. */
	return 0;
    }
/* ... Set up pointers into IWP */
    ljaci = 1;
    ljacj = ljaci + neqp1;
    lju = ljacj + nnzmx;
    ljlu = lju + max(lenplumx,neqp1);
    if (ireorder != 0) {
	lperm = ljlu + lenplumx;
	lqperm = lperm + *neq;
	liwk1 = lqperm + *neq;
	llevels = ljlu + nnzmx;
/* assumes that LENPLUFAC >= 2. */
	lmask = llevels + *neq;
    } else {
	lperm = 0;
	lqperm = 0;
	llevels = 0;
	lmask = 0;
	liwk1 = ljlu + lenplumx;
    }
    *liwp_min__ = liwk1 + (*neq << 1) - 1;
    if (ipremeth == 2) {
	*liwp_min__ += *neq << 1;
    }
    if (*liwp < *liwp_min__) {
	*ierr = 31;
/* more IWP length needed. */
	return 0;
    }
    *ierr = 0;
    return 0;
/* ------------  End of Subroutine DSPSETUP  ----------------------------- */
} /* dspsetup_ */

/* Subroutine */ int djacilu_(U_fp res, integer *ires, integer *neq, 
	doublereal *t, doublereal *y, doublereal *yprime, doublereal *rewt, 
	doublereal *savr, doublereal *wk, doublereal *h__, doublereal *cj, 
	doublereal *wp, integer *iwp, integer *ierr, doublereal *rpar, 
	integer *ipar)
{
    /* Initialized data */

    static char pmeth[8*4] = "ILUT    " "ILUTP   " "ILU0    " "MILU0   ";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer jscalcol, ireorder, ipremeth, lfililut, lenplumx, normtype,
	     lrownrms, i__, lenplufac, nre, lbw;
    static char msg[80];
    static integer lju, ubw, ljac, ifmt, ljlu, lplu, neqp1, liwk1, lrwk1, 
	    ljaci, ljacj, lmask;
    extern /* Subroutine */ int djilu_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, char *, integer *, integer *, ftnlen);
    static integer lperm;
    static char title[72];
    static integer iunit;
    extern /* Subroutine */ int prtmt_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, char *, char *, char *, char *
	    , integer *, integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen)
	    ;
    static doublereal sqrtn;
    static integer nnzmx;
    extern /* Subroutine */ int djcalc_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, U_fp, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *), amudia_(integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    roscal_(integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *);
    static integer jacout, liperm;
    extern /* Subroutine */ int dvperm_(integer *, doublereal *, integer *);
    static integer lqperm;
    extern /* Subroutine */ int xerrwd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, real *, real *, 
	    ftnlen);
    static integer lenpfac;
    extern /* Subroutine */ int djreord_(integer *, integer *, integer *, 
	    char *, doublereal *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, ftnlen);
    static char premeth[8];
    static integer llevels;
    static doublereal permtol;
    static integer isrnorm;
    static doublereal tolilut;

/* ... Version of 12-12-00 */
/* ... Calculate ILU decomposition of the Jacobian matrix */
/*     for use as a preconditioner by the DDASPK solver. */
/* ... Input arguments: */
/* total number of equations */
/* independent variable t */
/* most recent iterate of solution vector y */
/* most recent iterate of solution vector y' */
/* current residual evaluated at (T,Y,YPRIME) */
/* scale factors for Y and YPRIME */
/* function that evaluates residuals */
/* error flag for RES routine */
/* work space available to this subroutine */
/* current step size */
/* scalar proportional to 1/H */
/* user real workspace */
/* ... Output arguments: */
/* user integer workspace */
/* matrix elements of ILU */
/* array indices for elements of ILU */
/* Jacobian NRE is returned in IPAR(30) */
/* number of RES calls needed to evaluate */
/* ... Local variables: */
/* error flag (0 means success, else failure) */
/* =1 causes row normalization of Jac. */
/* 2-norm row scaling */
/* =0,1,2 for max-norm, 1-norm, or */
/* =1 causes row and column reordering of Jac. */
/*   be written to a file and then exit with */
/*   IRES = -2 to signal a stop to DDASPK. */
/* =1 causes the Jacobian matrix and SAVR to */
/* logical unit number to use when JACOUT .EQ. 1 */
/*   to be scaled by REWT-inverse */
/* =1 causes the columns of the Jacobian matrix */
    /* Parameter adjustments */
    --wk;
    --savr;
    --rewt;
    --yprime;
    --y;
    --wp;
    --iwp;
    --rpar;
    --ipar;

    /* Function Body */
/* ... Zero out NRE counter */
    nre = 0;
/* ... Load values from IPAR and RPAR. */
    lbw = ipar[1];
    ubw = ipar[2];
    lenpfac = ipar[3];
    lenplufac = ipar[4];
    ipremeth = ipar[5];
    lfililut = ipar[6];
    ireorder = ipar[7];
    isrnorm = ipar[8];
    normtype = ipar[9];
    jacout = ipar[10];
    jscalcol = ipar[11];
    tolilut = rpar[1];
    permtol = rpar[2];
    s_copy(premeth, pmeth + (ipremeth - 1 << 3), (ftnlen)8, (ftnlen)8);
/* ... Set pointers into the WP and IWP arrays. */
    neqp1 = *neq + 1;
    nnzmx = lenpfac * *neq;
    lenplumx = nnzmx + lenplufac * *neq;
/* ... Set up pointers into WP */
    ljac = 1;
    lrownrms = nnzmx + ljac;
    if (isrnorm == 1) {
	lplu = lrownrms + *neq;
    } else {
	lplu = lrownrms;
    }
    lrwk1 = lplu + lenplumx;
/* ... Set up pointers into IWP */
    ljaci = 1;
    ljacj = ljaci + neqp1;
    lju = ljacj + nnzmx;
    ljlu = lju + lenplumx;
/* ... Calculate Jacobian matrix. */
    *ierr = 0;
    djcalc_(neq, t, &y[1], &yprime[1], &savr[1], &lbw, &ubw, &wk[1], &rewt[1],
	     (U_fp)res, h__, cj, &nnzmx, &wp[ljac], &iwp[ljacj], &iwp[ljaci], 
	    &wp[lplu], &iwp[ljlu], &iwp[lju], &ipar[1], &rpar[1], ires, &nre, 
	    ierr);
    if (*ires < 0) {
	return 0;
    }
    if (*ierr != 0) {
	return 0;
    }
/* ... Save NRE value for user output. */
    ipar[30] += nre;
/* ... Modify pointers into IWP */
    ljlu = lju + neqp1;
    if (ireorder != 0) {
	lperm = ljlu + lenplumx;
	lqperm = lperm + *neq;
	liwk1 = lqperm + *neq;
	llevels = ljlu + nnzmx;
/* assumes that LENPLUFAC >= 2. */
	lmask = llevels + *neq;
    } else {
	lperm = 0;
	lqperm = 0;
	llevels = 0;
	lmask = 0;
	liwk1 = ljlu + lenplumx;
    }
    if (s_cmp(premeth, "ILUTP", (ftnlen)8, (ftnlen)5) == 0) {
	liperm = liwk1 + (*neq << 1);
    } else {
	liperm = liwk1;
    }
/* ... Multiply Jacobian columns by inverse of scaling vector REWT. */
/*     In PSOLILU, the WGHT array equals REWT/SQRT(NEQ), so we must */
/*     be consistent here. */
    if (jscalcol == 1) {
	sqrtn = sqrt((real) (*neq));
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wk[i__] = sqrtn / rewt[i__];
/* L10: */
	}
	amudia_(neq, &c__0, &wp[ljac], &iwp[ljacj], &iwp[ljaci], &wk[1], &wp[
		ljac], &iwp[ljacj], &iwp[ljaci]);
    }
/* ... Normalize Jacobian rows, if desired. */
    if (isrnorm == 1) {
	roscal_(neq, &c__0, &normtype, &wp[ljac], &iwp[ljacj], &iwp[ljaci], &
		wp[lrownrms], &wp[ljac], &iwp[ljacj], &iwp[ljaci], ierr);
	if (*ierr != 0) {
	    return 0;
	}
    }
/* ... Reorder Jacobian rows and columns, if desired. */
    if (ireorder != 0) {
	djreord_(neq, &neqp1, &nnzmx, premeth, &wp[ljac], &iwp[ljacj], &iwp[
		ljaci], &wp[lplu], &iwp[ljlu], &iwp[lju], &iwp[lperm], &iwp[
		lqperm], &iwp[llevels], &iwp[lmask], &ireorder, (ftnlen)8);
    }
/* ... Write matrix JAC and scaled RES to file if JACOUT .eq. 1. */
    if (jacout == 1) {
	iunit = ipar[29];
	if (isrnorm == 1) {
	    i__1 = *neq;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		savr[i__] *= wp[lrownrms + i__ - 1];
/* L20: */
	    }
	}
	if (ireorder != 0) {
	    dvperm_(neq, &savr[1], &iwp[lperm]);
	}
	s_copy(title, " DDASPK Test Matrix ", (ftnlen)72, (ftnlen)20);
	ifmt = 15;
	prtmt_(neq, neq, &wp[ljac], &iwp[ljacj], &iwp[ljaci], &savr[1], "NN", 
		title, "SPARSKIT", "RUA", &ifmt, &c__3, &iunit, (ftnlen)2, (
		ftnlen)72, (ftnlen)8, (ftnlen)3);
	s_copy(msg, "DJACILU -- Jacobian Matrix written to file.", (ftnlen)80,
		 (ftnlen)43);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	*ierr = 1;
	*ires = -2;
	return 0;
    }
/* ... Compute ILU decomposition. */
    i__1 = *neq + 1;
    djilu_(neq, &i__1, &nnzmx, &wp[ljac], &iwp[ljacj], &iwp[ljaci], &iwp[lju],
	     &wp[lplu], &iwp[ljlu], &wp[lrwk1], &iwp[liwk1], &lenplumx, &
	    tolilut, &lfililut, &permtol, premeth, &iwp[liperm], ierr, (
	    ftnlen)8);
    if (*ierr == -2 || *ierr == -3) {
	*ires = -2;
/* Stop since more storage needed. */
    }
/* ... Save pointers for use in DPSOLILU into IPAR array. */
    ipar[21] = lplu;
    ipar[22] = lju;
    ipar[23] = ljlu;
    ipar[24] = lrownrms;
    ipar[25] = lperm;
    ipar[26] = lqperm;
    return 0;
/* ------------  End of Subroutine DJACILU  ------------------------------ */
} /* djacilu_ */

/* Subroutine */ int djcalc_(integer *neq, doublereal *t, doublereal *y, 
	doublereal *yprime, doublereal *r0, integer *ml, integer *mu, 
	doublereal *r1, doublereal *rewt, S_fp res, doublereal *h__, 
	doublereal *cj, integer *nnzmx, doublereal *jac, integer *ja, integer 
	*ia, doublereal *rcoo, integer *jcoo, integer *icoo, integer *ipar, 
	doublereal *rpar, integer *ires, integer *nre, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j, i1, i2, jj, mba;
    static doublereal del;
    static char msg[80];
    static integer nnz, meb1;
    static doublereal squr;
    static integer mband;
    extern doublereal d1mach_(integer *);
    static integer meband;
    static doublereal delinv;
    extern /* Subroutine */ int coocsr_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *), 
	    xerrwd_(char *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, real *, real *, ftnlen);
    static doublereal uround, jacelem;

/* ... Version of 10-6-95 */
/* ... Calculate Jacobian matrix (derivatives with respect to each */
/*     dependent variable of the right-hand side of each rate equation). */
/*     Lower and upper bandwidths are used to select for computation */
/*     only those Jacobian elements that may be nonzero. */
/*     Estimates of Jacobian elements are computed by finite differences. */
/*     The Jacobian is stored in compressed sparse row format. */
/* ... Input arguments: */
/* total number of equations */
/* independent variable t */
/* most recent iterate of solution vector y */
/* most recent iterate of solution vector y' */
/* current residual evaluated at (T,Y,YPRIME) */
/* array of scaling factors for Y and YPRIME */
/* function that evaluates residuals */
/* error flag for RES routine */
/* lower and upper bandwidths */
/* maximum number of nonzeros in Jacobian */
/* current step size */
/* scalar proportional to 1/H */
/* user real workspace */
/* ... Work-array argument: */
/* user integer workspace */
/* ... Output arguments: */
/* work space available to this subroutine */
/* nonzero Jacobian elements */
/* col indices of nonzero Jacobian elements */
/* pointers to beginning of each row in JAC,JA */
/* ... Workspace for temporary storage of Jacobian elements: */
/* nonzero Jacobian elements */
/* col indices of nonzero Jacobian elements */
/* ... Local variables: */
/* row indices of nonzero Jacobian elements */
/* ... Set band parameters. */
    /* Parameter adjustments */
    --ia;
    --rewt;
    --r1;
    --r0;
    --yprime;
    --y;
    --icoo;
    --jcoo;
    --rcoo;
    --ja;
    --jac;
    --ipar;
    --rpar;

    /* Function Body */
    nnz = 1;
    mband = *ml + *mu + 1;
    mba = min(mband,*neq);
    meband = mband + *ml;
    meb1 = meband - 1;
/* ... Set the machine unit roundoff UROUND and SQRT(UROUND), used to */
/* ... set increments in the difference quotient procedure. */
    uround = d1mach_(&c__4);
    squr = sqrt(uround);
/* ... Initial error flags. */
    *ierr = 0;
    *ires = 0;
/* ... Make MBA calls to RES to approximate the Jacobian. */
/* ... Here, R0(1),...,R0(neq) contains the base RES value, and */
/*     R1(1),...,R1(NEQ) contains the perturbed values of RES. */
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *neq;
	i__3 = mband;
	for (jj = j; i__3 < 0 ? jj >= i__2 : jj <= i__2; jj += i__3) {
	    jac[jj] = y[jj];
	    jac[jj + *neq] = yprime[jj];
/* Computing MAX */
	    d__4 = (d__1 = y[jj], abs(d__1)), d__5 = (d__2 = *h__ * yprime[jj]
		    , abs(d__2)), d__4 = max(d__4,d__5), d__5 = (d__3 = 1.f / 
		    rewt[jj], abs(d__3));
	    del = squr * max(d__4,d__5);
	    d__1 = *h__ * yprime[jj];
	    del = d_sign(&del, &d__1);
	    del = y[jj] + del - y[jj];
	    y[jj] += del;
	    yprime[jj] += *cj * del;
/* L10: */
	}
	(*res)(t, &y[1], &yprime[1], cj, &r1[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	++(*nre);
	i__3 = *neq;
	i__2 = mband;
	for (jj = j; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    y[jj] = jac[jj];
	    yprime[jj] = jac[jj + *neq];
/* Computing MAX */
	    d__4 = (d__1 = y[jj], abs(d__1)), d__5 = (d__2 = *h__ * yprime[jj]
		    , abs(d__2)), d__4 = max(d__4,d__5), d__5 = (d__3 = 1.f / 
		    rewt[jj], abs(d__3));
	    del = squr * max(d__4,d__5);
	    d__1 = *h__ * yprime[jj];
	    del = d_sign(&del, &d__1);
	    del = y[jj] + del - y[jj];
	    delinv = 1.f / del;
/* Computing MAX */
	    i__4 = 1, i__5 = jj - *mu;
	    i1 = max(i__4,i__5);
/* Computing MIN */
	    i__4 = *neq, i__5 = jj + *ml;
	    i2 = min(i__4,i__5);
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* ... Calculate possibly nonzero Jacobian elements for this variable, */
/*     and store nonzero elements in coordinate format. */
		jacelem = (r1[i__] - r0[i__]) * delinv;
		if (jacelem != 0.f) {
		    if (nnz > *nnzmx) {
			s_copy(msg, "DJCALC -- More storage needed for Jacob"
				"ian.", (ftnlen)80, (ftnlen)43);
			xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &
				c__0, &c__0, &c_b21, &c_b21, (ftnlen)80);
			s_copy(msg, "DJCALC -- Increase LENPFAC.", (ftnlen)80,
				 (ftnlen)27);
			xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &
				c__0, &c__0, &c_b21, &c_b21, (ftnlen)80);
			s_copy(msg, "DJCALC -- Storage exceeded at (I,J) = ("
				"I1,I2)", (ftnlen)80, (ftnlen)45);
			xerrwd_(msg, &c__80, &c__0, &c__0, &c__2, &i__, &jj, &
				c__0, &c_b21, &c_b21, (ftnlen)80);
			*ierr = 1;
			*ires = -2;
			return 0;
		    }
		    rcoo[nnz] = jacelem;
		    jcoo[nnz] = jj;
		    icoo[nnz] = i__;
		    ++nnz;
		}
/* L20: */
	    }
/* L30: */
	}
/* L40: */
    }
    --nnz;

/* ... Convert Jacobian from coordinate to compressed sparse row format. */
    coocsr_(neq, &nnz, &rcoo[1], &icoo[1], &jcoo[1], &jac[1], &ja[1], &ia[1]);
    return 0;
/* ------------  End of Subroutine DJCALC  ------------------------------- */
} /* djcalc_ */

/* Subroutine */ int dpsolilu_(integer *neq, doublereal *t, doublereal *y, 
	doublereal *yprime, doublereal *r0, doublereal *wk, doublereal *cj, 
	doublereal *wght, doublereal *wp, integer *iwp, doublereal *bl, 
	doublereal *eplin, integer *ierr, doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer jscalcol, ireorder, ipremeth, lrownrms, i__, lju, ljlu, 
	    lplu, lperm;
    extern /* Subroutine */ int lusol_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), dvperm_(integer *, 
	    doublereal *, integer *);
    static integer lqperm, isrnorm;

/* ... Version of 12-5-00 */
/* ... Solve the linear system P*x=c, using elements of P loaded into */
/*     arrays WP and IWP. */

/* ... Input arguments: */
/* total number of equations */
/* independent variable t */
/* most recent iterate of solution vector y */
/* most recent iterate of solution vector y' */
/* function value G(T,Y,YPRIME) */
/* scaling array for Y and YPRIME */
/* scalar proportional to 1/H */
/* not used */
/* matrix elements of ILU */
/* array indices for elements of ILU */
/* user workspace */
/* ... Work-array argument: */
/* user workspace */
/* ... In-out argument: */
/* work space available to this subroutine */
/* ... Output arguments: */
/* on input, c of P*x=c; on output, x */
/* ... Local variables: */
/* error flag (0 only possible value here) */
/* ... Load IPREMETH, IREORDER and ISRNORM values from IPAR. */
    /* Parameter adjustments */
    --bl;
    --wght;
    --wk;
    --r0;
    --yprime;
    --y;
    --wp;
    --iwp;
    --rpar;
    --ipar;

    /* Function Body */
    ipremeth = ipar[5];
    ireorder = ipar[7];
    isrnorm = ipar[8];
    jscalcol = ipar[11];
/* ... Load pointers into WP and iWP arrays. */
    lplu = ipar[21];
    lju = ipar[22];
    ljlu = ipar[23];
    lrownrms = ipar[24];
    lperm = ipar[25];
    lqperm = ipar[26];
/* ... Scale c by multiplying by row-normalization factors, if used. */
    if (isrnorm == 1) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    bl[i__] *= wp[lrownrms + i__ - 1];
/* L10: */
	}
    }
/* ... Solve P*x=c for a preconditioner stored as a sparse matrix in */
/*     compressed sparse row format. */
/*     If rows and columns of P were reordered (permuted), permute c, */
/*     then use inverse permutation on x. */
    if (ipremeth == 1 || ipremeth == 2) {
	if (ireorder == 1) {
	    dvperm_(neq, &bl[1], &iwp[lperm]);
	}
	lusol_(neq, &bl[1], &wk[1], &wp[lplu], &iwp[ljlu], &iwp[lju]);
	if (ireorder == 1) {
	    dvperm_(neq, &wk[1], &iwp[lqperm]);
	}
    }
/* ... Unscale x by dividing by column scaling vector WGHT. */
    if (jscalcol == 1) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	    bl[i__] = wk[i__] / wght[i__];
	}
    } else {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	    bl[i__] = wk[i__];
	}
    }
    *ierr = 0;
    return 0;
/* ------------  End of Subroutine DPSOLILU  ----------------------------- */
} /* dpsolilu_ */

/* Subroutine */ int djilu_(integer *neq, integer *neqp1, integer *nnzmx, 
	doublereal *jac, integer *ja, integer *ia, integer *ju, doublereal *
	plu, integer *jlu, doublereal *rwk1, integer *iwk1, integer *lenplumx,
	 doublereal *tolilut, integer *lfililut, doublereal *permtol, char *
	premeth, integer *iperm, integer *ierr, ftnlen premeth_len)
{
    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static char msg[80];
    extern /* Subroutine */ int ilut_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    static logical error;
    extern /* Subroutine */ int ilutp_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerrwd_(char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    real *, real *, ftnlen);

/* ... Version of 12-12-00 */
/* ... Compute ILU decomposition of Jacobian and return it in any one */
/*     of the storage formats. */
/* ... Input arguments: */
/* total number of equations */
/* NEQ + 1 */
/* maximum number of nonzeros in Jacobian */
/* nonzero Jacobian elements */
/* col indices of nonzero Jacobian elements */
/* pointers to beginning of each row in jac,ja */
/* ... Output arguments: */
/* matrix elements of ILU */
/* sizes and array indices for elements of ILU */
/* matrix PLU,JLU */
/* pointer to beginning of each row of U in */
/* ... Local variables: */
/* error flag */
    /* Parameter adjustments */
    --iperm;
    --iwk1;
    --rwk1;
    --ju;
    --ia;
    --ja;
    --jac;
    --jlu;
    --plu;

    /* Function Body */
    error = FALSE_;
    if (s_cmp(premeth, "ILUT", (ftnlen)8, (ftnlen)4) == 0) {
/* ... Use incomplete factorization routine ILUT from SparsKit. */
	ilut_(neq, &jac[1], &ja[1], &ia[1], lfililut, tolilut, &plu[1], &jlu[
		1], &ju[1], lenplumx, &rwk1[1], &iwk1[1], ierr);
	if (*ierr != 0) {
	    s_copy(msg, "DJILU -- Error return from ILUT: IERR = (I1)", (
		    ftnlen)80, (ftnlen)44);
	    xerrwd_(msg, &c__80, &c__0, &c__0, &c__1, ierr, &c__0, &c__0, &
		    c_b21, &c_b21, (ftnlen)80);
	    error = TRUE_;
	}
    } else if (s_cmp(premeth, "ILUTP", (ftnlen)8, (ftnlen)5) == 0) {
/* ... Use incomplete factorization routine ILUTP from SparsKit. */
	ilutp_(neq, &jac[1], &ja[1], &ia[1], lfililut, tolilut, permtol, neq, 
		&plu[1], &jlu[1], &ju[1], lenplumx, &rwk1[1], &iwk1[1], &
		iperm[1], ierr);
	if (*ierr != 0) {
	    s_copy(msg, "DJILU -- Error return from ILUTP: IERR = (I1)", (
		    ftnlen)80, (ftnlen)45);
	    xerrwd_(msg, &c__80, &c__0, &c__0, &c__1, ierr, &c__0, &c__0, &
		    c_b21, &c_b21, (ftnlen)80);
	    error = TRUE_;
	}
/* ... Put in other options here for incomplete factorizations. */
    }
    if (error) {
	s_copy(msg, "DJILU -- IERR .NE. 0 means one of the following has occ"
		"urred:", (ftnlen)80, (ftnlen)61);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    IERR >  0   --> Zero pivot encountered at step numb"
		"er IERR.", (ftnlen)80, (ftnlen)63);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    IERR = -1   --> Error. input matrix may be wrong.", (
		ftnlen)80, (ftnlen)53);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "                     (The elimination process has gener"
		"ated a", (ftnlen)80, (ftnlen)61);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "                     row in L or U with length > NEQ.)", 
		(ftnlen)80, (ftnlen)54);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    IERR = -2   --> Matrix L overflows.", (ftnlen)80, (
		ftnlen)39);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    IERR = -3   --> Matrix U overflows.", (ftnlen)80, (
		ftnlen)39);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    IERR = -4   --> Illegal value for LFILILUT.", (
		ftnlen)80, (ftnlen)47);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    IERR = -5   --> Zero row encountered.", (ftnlen)80, (
		ftnlen)41);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    ", (ftnlen)80, (ftnlen)4);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    For IERR = -2 or -3, increase the value of LENPLUFA"
		"C or", (ftnlen)80, (ftnlen)59);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    decrease the value of LFILILUT if LENPLUFAC cannot "
		"be", (ftnlen)80, (ftnlen)57);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
	s_copy(msg, "    increased.", (ftnlen)80, (ftnlen)14);
	xerrwd_(msg, &c__80, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b21,
		 &c_b21, (ftnlen)80);
    }
    return 0;
/* ------------  End of Subroutine DJILU  -------------------------------- */
} /* djilu_ */

/* Subroutine */ int djreord_(integer *neq, integer *neqp1, integer *nnzmx, 
	char *premeth, doublereal *jac, integer *ja, integer *ia, doublereal *
	awk, integer *jwk, integer *iwk, integer *perm, integer *qperm, 
	integer *levels, integer *mask, integer *ireorder, ftnlen premeth_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int bfs_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *), atob_(integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    static integer nlev;
    extern /* Subroutine */ int dperm_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *);
    static integer nfirst;
    extern /* Subroutine */ int rversp_(integer *, integer *);
    static integer maskval;

/* ... Version of 10-6-95 */
/* ... If desired, reorder the Jacobian matrix. */
/* ... Input arguments: */
/* total number of equations */
/* NEQ + 1 */
/* maximum number of nonzeroes in Jacobian */
/* nonzero Jacobian elements */
/* column indices of nonzero Jacobian elements */
/* indices of 1st nonzero element in each row */
/* ... Work-array arguments: */
/* used in reordering the rows and columns of */
/* the Jacobian matrix. */
/* Integer array containing the permutation */
/* permutation in array perm. */
/* Integer array holding the inverse of the */
/* subroutine.   See subroutine BFS for */
/* more details. */
/* Work array used by the bfs reordering */
/* subroutine.  See BFS subroutine. */
/* Work array used by the BFS reordering */
/* of the Jacobian matrix is desired. */
/* = 1 means a reverse Cuthill-McKee */
/*     reordering of the rows and columns */
/*     of the Jacobian is done. */
/* = 0 means no reordering. */
/* ... Local variables: */
/* Flag used to determine if a reordering */
/* See subroutine BFS for more details. */
/* Number of levels in levels array. */
/* Scalar used with MASK. */
    /* Parameter adjustments */
    --mask;
    --levels;
    --qperm;
    --perm;
    --iwk;
    --ia;
    --jwk;
    --awk;
    --ja;
    --jac;

    /* Function Body */
    if (*ireorder == 1) {
/* ... Copy JAC, JA, and IA to AWK, JWK, and IWK. */
	atob_(neq, &jac[1], &ja[1], &ia[1], &awk[1], &jwk[1], &iwk[1]);
/* ... Perform a Cuthill-McKee reordering of the Jacobian. */
	nfirst = 1;
	perm[1] = 0;
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mask[i__] = 1;
	}
	maskval = 1;
	qperm[1] = 1;
	bfs_(neq, &jwk[1], &iwk[1], &nfirst, &perm[1], &mask[1], &maskval, &
		qperm[1], &levels[1], &nlev);
/* ... Reverse the permutation to obtain the reverse Cuthill-McKee */
/*     reordering. */
	rversp_(neq, &qperm[1]);
/* ... Calculate the inverse of QPERM and put it in PERM. */
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    perm[qperm[i__]] = i__;
	}
/* ... Permute rows and columns of Jacobian using PERM. */
	dperm_(neq, &awk[1], &jwk[1], &iwk[1], &jac[1], &ja[1], &ia[1], &perm[
		1], &perm[1], &c__1);
/* ... End of If block */
    }
    return 0;
/* ------------  End of Subroutine DJREORD  ------------------------------ */
} /* djreord_ */

