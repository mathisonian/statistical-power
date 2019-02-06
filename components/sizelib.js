/* size_lib.js
 * 2005 April
 * REMorse, for the MGH Biostatistics Center
 * This file is based off of fortran code used by David Schoenfeld
 * to implement the original sample size program.
 * It contains a number of functions shared between the different
 * test types.
 */


/* The various procedures have a number of named constants used in them.
 * Using named constants makes the code cleaner than having the contstants
 * provided inline.
 * Rather than have the contstants have to be re-evaluated each time through
 * the function, we declare them here.  I'm trying to mimic C's "static"
 * declaration.
 */

var ZERO = 0.0, HALF = 0.5, ONE = 1.0, ONE_P_HALF = 1.5, TWO = 2.0, FOUR = 4.0, TWELVE = 12.0;

// constants for function alnorm
var ALNORM_A1 = 5.75885480458, ALNORM_A2 = 2.62433121679, ALNORM_A3 = 5.92885724438;
var ALNORM_B1 = -29.8213557807, ALNORM_B2 = 48.6959930692;
var ALNORM_C1 = -3.8052e-8, ALNORM_C2 = 3.98064794e-4, ALNORM_C3 = -.151679116635, ALNORM_C4 = 4.8385912808, ALNORM_C5 = .742380924027, ALNORM_C6 = 3.99019417011;
var ALNORM_D1 = 1.00000615302, ALNORM_D2 = 1.98615381364, ALNORM_D3 = 5.29330324926, ALNORM_D4 = -15.1508972451, ALNORM_D5 = 30.789933034;

var ALNORM_LTONE = 7.0, ALNORM_UTZERO = 18.66, ALNORM_CON = 1.28;
var ALNORM_P = 0.398942280444, ALNORM_Q = .39990348504, ALNORM_R = .398942280385;

// constants for function alngam
var ALNGAM_ALR2PI = .918938533204673;
var ALNGAM_R1 = [ -2.66685511495, -24.4387534237, -21.9698958928, 11.1667541262, 3.13060547623, .607771387771, 11.9400905721, 31.4690115749, 15.234687407 ];
var ALNGAM_R2 = [ -78.3359299449, -142.046296688, 137.519416416, 78.6994924154, 4.16438922228, 47.066876606, 313.399215894, 263.505074721, 43.3400022514 ];
var ALNGAM_R3 = [ -212159.572323, 230661.510616, 27464.7644705, -40262.1119975, -2296.6072978, -116328.495004 ,-146025.937511, -24235.7409629, -570.691009324 ];
var ALNGAM_R4 = [ .279195317918525, .4917317610505968, .0692910599291889, 3.350343815022304, 6.012459259764103 ];
var ALNGAM_XLGE = 5.1e6, ALNGAM_XLGST = 1e305; // These are architecture specific constants.  They probably need to be fixed...

// constants for the betain function
var BETAIN_ACU = 1e-15;

// constants for the tnc function
var TNC_ITRMAX = 100.1, TNC_ERRMAX = 1e-6;
var TNC_R2PI = .79788456080286535588, TNC_ALNRPI = .57236494292470008707;

// constants for the ppnd7 function
var PPND7_SPLIT1 = 0.425, PPND7_SPLIT2 = 5.0;
var PPND7_CONST1 = 0.180625, PPND7_CONST2 = 1.6;
var PPND7_A0 = 3.3871327179E+00, PPND7_A1 = 5.0434271938E+01, PPND7_A2 = 1.5929113202E+02, PPND7_A3 = 5.9109374720E+01;
var PPND7_B1 = 1.7895169469E+01, PPND7_B2 = 7.8757757664E+01, PPND7_B3 = 6.7187563600E+01;
var PPND7_C0 = 1.4234372777E+00, PPND7_C1 = 2.7568153900E+00, PPND7_C2 = 1.3067284816E+00, PPND7_C3 = 1.7023821103E-01;
var PPND7_D1 = 7.3700164250E-01, PPND7_D2 = 1.2021132975E-01;
var PPND7_E0 = 6.6579051150E+00, PPND7_E1 = 3.0812263860E+00, PPND7_E2 = 4.2868294337E-01, PPND7_E3 = 1.7337203997E-02;
var PPND7_F1 = 2.4197894225E-01, PPND7_F2 = 1.2258202635E-02;

// constants for the ppft function
var PPFT_B21 = 4.0;
var PPFT_B31 = 96.0, PPFT_B32 = 5.0, PPFT_B33 = 16.0, PPFT_B34 = 3.0;
var PPFT_B41 = 384.0, PPFT_B42 = 3.0, PPFT_B43 = 19.0, PPFT_B44 = 17.0, PPFT_B45 = -15.0;
var PPFT_B51 = 9216.0, PPFT_B52 = 79.0, PPFT_B53 = 776.0, PPFT_B54 = 1482.0, PPFT_B55 = -1920.0, PPFT_B56 = -945.0;

//constants for the PPFNML function
var PPFNML_P0 = -.322232431088, PPFNML_P1 = -1.0, PPFNML_P2 = -.342242088547, PPFNML_P3 = -.0204231210245, PPFNML_P4 = -4.53642210148e-5;
var PPFNML_Q0 = .099348462606, PPFNML_Q1 = .588581570495, PPFNML_Q2 = .531103462366, PPFNML_Q3 = .10353775285, PPFNML_Q4 = .0038560700634;

function alnorm(x, upper) {
    /* Algorithm AS66 Applied Statistics (1973) Vol. 22,  No. 3
     * Evaluates the tail area of the standardised normal curve
     * from x to infinity if upper is true or
     * from minus infinity to x if upper is false
     */

    var ret = 0.0;

    if (x < ZERO) {
        upper = !upper;
        x = -x;
    }

    if ((x <= ALNORM_LTONE) || (upper && (x <= ALNORM_UTZERO))) {
        var y = x * x * HALF;
        if (x > ALNORM_CON) {
            ret = ALNORM_R * Math.exp(-y) / (x + ALNORM_C1 + ALNORM_D1 / (x + ALNORM_C2 + ALNORM_D2 / (x + ALNORM_C3 + ALNORM_D3 / (x + ALNORM_C4 + ALNORM_D4 / (x + ALNORM_C5 + ALNORM_D5 / (x + ALNORM_C6))))));
        } else {
            ret = HALF - x * (ALNORM_P - ALNORM_Q * y / (y + ALNORM_A1 + ALNORM_B1 / (y + ALNORM_A2 + ALNORM_B2 / (y + ALNORM_A3))));
        }
    }

    if (!upper) {
        ret = ONE - ret;
    }

    return ret;
}

function alngam(x) {
    /* Algorithm AS245  Applied Statistics (1989) Vol. 38, No. 2
     * Calculation of the logarithm of the gamma function
     */

    var ret = 0.0;

    // make sure x is a valid option
    if (x >= ALNGAM_XLGST) {
        throw "alngam: x >= xlgst";
    } else if (x <= ZERO) {
        throw "alngam: x < zero";
    }

    if (x < ONE_P_HALF) {
        var y;
        if (x < HALF) {
            ret = - Math.log(x);
            y = x + ONE;
            if (y == ONE) {
                // x < machine epsilon; break out early
                return ret;
            }
        } else {
            ret = ZERO;
            y = x;
            x = (x - HALF) - HALF;
        }
        ret += x * ((((ALNGAM_R1[4] * y + ALNGAM_R1[3]) * y + ALNGAM_R1[2]) * y + ALNGAM_R1[1]) * y + ALNGAM_R1[0]) / ((((y + ALNGAM_R1[8]) * y + ALNGAM_R1[7]) * y + ALNGAM_R1[6]) * y + ALNGAM_R1[5]);
    } else if (x < FOUR) {
        var y = x - ONE - ONE;
        ret = y * ((((ALNGAM_R2[4] * x + ALNGAM_R2[3]) * x + ALNGAM_R2[2]) * x + ALNGAM_R2[1]) * x + ALNGAM_R2[0]) / ((((x + ALNGAM_R2[8]) * x + ALNGAM_R2[7]) * x + ALNGAM_R2[6]) * x + ALNGAM_R2[5]);
    } else if (x < TWELVE) {
        ret = ((((ALNGAM_R3[4] * x + ALNGAM_R3[3]) * x + ALNGAM_R3[2]) * x + ALNGAM_R3[1]) * x + ALNGAM_R3[0]) / ((((x + ALNGAM_R3[8]) * x + ALNGAM_R3[7]) * x + ALNGAM_R3[6]) * x + ALNGAM_R3[5]);
    } else {
        var y = Math.log(x);
        ret = x * (y - ONE) - HALF * y + ALNGAM_ALR2PI;
        if (x <= ALNGAM_XLGE) {
            var x1 = ONE / x;
            var x2 = x1 * x1;
            ret += x1 * ((ALNGAM_R4[2] * x2 + ALNGAM_R4[1]) * x2 + ALNGAM_R4[0]) / ((x2 + ALNGAM_R4[4]) * x2 + ALNGAM_R4[3]);
        }
    }

    return ret;
}

function betain(x, p, q, beta) {
    /* Algorithm AS63  Applied Statistics (1973), Vol. 22, No. 3
     * computes incomplete beta function ratio for arguments
     * x between zero and one, p and q positive.
     * log of complete beta function, beta, is assumed to be known
     */

    var ret = x;

    var ai, cx, pp, ns, qq, rx, xx, psq, indx, temp, term;

    if ((p <= ZERO) || (q <= ZERO)) {
        throw "betain: p or q <= zero";
    } else if ((x < ZERO) || (x > ONE)) {
        throw "betain: x < zero or x > one";
    }

    if ((x == ZERO) || (x == ONE)) {
        return ret;
    }

    psq = p + q;
    cx = ONE - x;

    if (p >= (psq * x)) {
        xx = x;
        pp = p;
        qq = q;
        indx = false;
    } else {
        xx = cx;
        cx = x;
        pp = q;
        qq = p;
        indx = true;
    }

    term = ONE;
    ai = ONE;
    ret = ONE;
    ns = Math.round(qq + cx * psq);
    rx = xx / cx;

    temp = qq - ai;
    if (ns == 0) {
        rx = xx;
    }
    term = term * temp * rx / (pp + ai);
    ret += term;
    temp = Math.abs(term);

    while(!((temp <= BETAIN_ACU) && (temp <= (BETAIN_ACU * ret)))) {
        ai += ONE;
        ns = ns - 1;
        if (ns < 0) {
            temp = psq;
            psq += ONE;
        }
        else {
            temp = qq - ai;
            if (ns == 0) {
                rx = xx;
            }
        }
        term = term * temp * rx / (pp + ai);
        ret += term;
        temp = Math.abs(term);
    }

    ret = ret * Math.exp(pp * Math.log(xx) + (qq - ONE) * Math.log(cx) - beta) / pp;

    if (indx) {
        ret = ONE - ret;
    }

    return ret;
}


function tnc(t, df, delta) {
    /* ALGORITHM AS 243  APPL. STATIST. (1989), VOL.38, NO. 1
     *
     * Cumulative probability at T of the non-central t-distribution
     * with DF degrees of freedom (may be fractional) and non-centrality
     * parameter DELTA.
     *
     * Note - requires the following auxiliary routines
     * ALOGAM (X)                         - ACM 291 or AS 245
     * BETAIN (X, A, B, ALBETA, IFAULT)   - AS 63 (updated in ASR 19)
     * ALNORM (X, UPPER)                  - AS 66
     */

    var ret = ZERO;

    var a, b, p, q, s, x, en, tt, del, rxb, godd, xodd, geven, xeven, errbd, lambda, albeta, negdel;

    if (df <= ZERO) {
        throw "tnc: df <= zero";
    }

    tt = t;
    del = delta;
    negdel = false;
    if (t < ZERO) {
        negdel = true;
        tt = -tt;
        del = -del;
    }

    en = ONE;
    x = t * t / (t * t + df);
    if (x > ZERO) {
        lambda = del * del;
        p = HALF * Math.exp(-HALF * lambda);
        q = TNC_R2PI * p * del;
        s = HALF - p;
        a = HALF;
        b = HALF * df;
        rxb = Math.pow((ONE - x), b);
        albeta = TNC_ALNRPI + alngam(b) - alngam(a + b);
        xodd = betain(x, a, b, albeta);
        godd = TWO * rxb * Math.exp(a * Math.log(x) - albeta);
        xeven = ONE - rxb;
        geven = b * x * rxb;

        ret = (p * xodd) + (q * xeven);

        do {
            a += ONE;
            xodd -= godd;
            xeven -= geven;
            godd = godd * x * (a + b - ONE) / a;
            geven = geven * x * (a + b - HALF) / (a + HALF);
            p = p * lambda / (TWO * en);
            q = q * lambda / (TWO * en + ONE);
            s -= p;
            en += ONE;
            ret += (p * xodd) + (q * xeven);
            errbd = TWO * s * (xodd - godd);
        } while ( ((errbd > TNC_ERRMAX) && (en <= TNC_ITRMAX)) );
    }

    if (en > TNC_ITRMAX) {
        return ret;
        throw "tnc: too many iterations";
    }

    var ddel = del;
    ret += alnorm(ddel, true);
    if (negdel) {
        ret = ONE - ret;
    }

    return ret;
}

function ppnd7(p) {
    /* ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-
     * 484.
     *
     * Produces the normal deviate Z corresponding to a given lower
     * tail area of P; Z is accurate to about 1 part in 10**7.
     */

    var ret = ZERO;
    var q, r;

    q = p - HALF;
    if (Math.abs(q) <= PPND7_SPLIT1) {
        r = PPND7_CONST1 - (q * q);
        ret = q * (((PPND7_A3 * r + PPND7_A2) * r + PPND7_A1) * r + PPND7_A0) / (((PPND7_B3 * r + PPND7_B2) * r + PPND7_B1) * r + ONE);
    } else {
        if (q < ZERO) {
            r = p;
        } else {
            r = ONE - p;
        }
        if (r <= ZERO) {
            //throw "ppnd7: r < zero (because p > one)";
            return ZERO;
        }
        r = Math.sqrt(-Math.log(r));
        if (r <= PPND7_SPLIT2) {
            r = r - PPND7_CONST2;
            ret = (((PPND7_C3 * r + PPND7_C2) * r + PPND7_C1) * r + PPND7_C0) / ((PPND7_D2 * r + PPND7_D1) * r + ONE);
        } else {
            r = r - PPND7_SPLIT2;
            ret = (((PPND7_E3 * r + PPND7_E2) * r + PPND7_E1) * r + PPND7_E0) / ((PPND7_F2 * r + PPND7_F1) * r + ONE);
        }
        if (q < ZERO) {
            ret = -ret;
        }
    }
    return ret;
}

function ppft(p, idf) {
    /* LATEST REVISION  -  2005 Apr 22 (REM)
     *
     * THIS FUNCTION IS A VERSION OF DATAPAC SUBROUTINE
     * TPPF, WITH MODIFICATIONS TO FACILITATE CONVERSION TO
     * DOUBLE PRECISION AUTOMATICALLY USING THE NAG, INC. CODE APT,
     * AND TO CORRESPOND TO STARPAC CONVENTIONS.
     *
     * PURPOSE--THIS SUBROUTINE COMPUTES THE PERCENT POINT
     *          FUNCTION VALUE FOR THE STUDENT"S T DISTRIBUTION
     *          WITH INTEGER DEGREES OF FREEDOM PARAMETER = IDF.
     *          THE STUDENT"S T DISTRIBUTION USED
     *          HEREIN IS DEFINED FOR ALL X,
     *          AND ITS PROBABILITY DENSITY FUNCTION IS GIVEN
     *          IN THE REFERENCES BELOW.
     *          NOTE THAT THE PERCENT POINT FUNCTION OF A DISTRIBUTION
     *          IS IDENTICALLY THE SAME AS THE INVERSE CUMULATIVE
     *          DISTRIBUTION FUNCTION OF THE DISTRIBUTION.
     * ERROR CHECKING--NONE
     * RESTRICTIONS--IDF SHOULD BE A POSITIVE INTEGER VARIABLE.
     *             --P SHOULD BE BETWEEN 0.0D0 (EXCLUSIVELY)
     *               AND 1.0D0 (EXCLUSIVELY).
     * COMMENT--FOR IDF = 1 AND IDF = 2, THE PERCENT POINT FUNCTION
     *          FOR THE T DISTRIBUTION EXISTS IN SIMPLE CLOSED FORM
     *          AND SO THE COMPUTED PERCENT POINTS ARE EXACT.
     *        --FOR OTHER SMALL VALUES OF IDF (IDF BETWEEN 3 AND 6,
     *          INCLUSIVELY), THE APPROXIMATION
     *          OF THE T PERCENT POINT BY THE FORMULA
     *          GIVEN IN THE REFERENCE BELOW IS AUGMENTED
     *          BY 3 ITERATIONS OF NEWTON"S METHOD FOR
     *          ROOT DETERMINATION.
     *          THIS IMPROVES THE ACCURACY--ESPECIALLY FOR
     *          VALUES OF P NEAR 0 OR 1.
     * REFERENCES--NATIONAL BUREAU OF STANDARDS APPLIED MATHMATICS
     *             SERIES 55, 1964, PAGE 949, FORMULA 26.7.5.
     *           --JOHNSON AND KOTZ, CONTINUOUS UNIVARIATE
     *             DISTRIBUTIONS--2, 1970, PAGE 102,
     *             FORMULA 11.
     *           --FEDERIGHI, "EXTENDED TABLES OF THE
     *             PERCENTAGE POINTS OF STUDENT"S T
     *             DISTRIBUTION, JOURNAL OF THE
     *             AMERICAN STATISTICAL ASSOCIATION,
     *             1969, PAGES 683-688.
     *           --HASTINGS AND PEACOCK, STATISTICAL
     *             DISTRIBUTIONS--A HANDBOOK FOR
     *             STUDENTS AND PRACTITIONERS, 1975,
     *             PAGES 120-123.
     * WRITTEN BY--JAMES J. FILLIBEN
     *             STATISTICAL ENGINEERING DIVISION
     *             NATIONAL BUREAU OF STANDARDS
     *             WASHINGTON, D. C. 20234
     * ORIGINAL VERSION--OCTOBER   1975.
     * UPDATED         --NOVEMBER  1975.
     *
     * MODIFIED BY     --JANET R. DONALDSON, DECEMBER 7, 1981
     *                   STATISTICAL ENGINEERING DIVISION
     *                   NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
     *
     * Translated to Javascript by Richard Morse, MGH Biostatistics Center
     *                             April 2005
     */

    var maxit = 5;
    var ret = 0.0;

    if (idf == 1) {
        var arg = Math.PI * p;
        ret = -Math.cos(arg)/Math.sin(arg);
    } else if (idf == 2) {
        ret = (Math.SQRT2/2.0) * ((2.0 * p) - 1.0) / (Math.sqrt(p * (1.0 - p)))
    } else if (idf >= 3) {
        ret = ppfnml(p);
        var d1 = ret;
        var d3 = Math.pow(d1, 3);
        var d5 = Math.pow(d1, 5);
        var d7 = Math.pow(d1, 7);
        var d9 = Math.pow(d1, 9);
        var term1 = d1;
        var term2 = (1.0 / PPFT_B21) * (d3 + d1) / idf;
        var term3 = (1.0 / PPFT_B31) * (PPFT_B32*d5 + PPFT_B33*d3 + PPFT_B34*d1) / (Math.pow(idf, 2));
        var term4 = (1.0 / PPFT_B41) * (PPFT_B42*d7 + PPFT_B43*d5 + PPFT_B44*d3 + PPFT_B45*d1) / (Math.pow(idf, 3));
        var term5 = (1.0 / PPFT_B51) * (PPFT_B52*d9 + PPFT_B53*d7 + PPFT_B54*d5 + PPFT_B55*d3 + PPFT_B56*d1) / (Math.pow(idf, 4));
        ret = term1 + term2 + term3 + term4 + term5;
        if (idf == 3) {
            var con = Math.PI * (p - 0.5);
            var arg = ret / Math.sqrt(idf);
            var z = Math.atan(arg);
            var s, c;
            for(var i = 0; i < maxit; i++) {
                s = Math.sin(z); c = Math.cos(z);
                z = z - (z + s * c - con) / (2.0 * c * c);
            }
            ret = Math.sqrt(idf) * s/c;
        } else if (idf == 4) {
            var con = 2.0 * (p - 0.5);
            var arg = ret / Math.sqrt(idf);
            var z = Math.atan(arg);
            var s, c;
            for(var i = 0; i < maxit; i++) {
                s = Math.sin(z); c = Math.cos(z);
                z = z - ((1.0 + 0.5 * c * c) * s - con) / (1.5 * c * c * c);
            }
            ret = Math.sqrt(idf) * s/c;
        } else if (idf == 5) {
            var con = Math.PI * (p - 0.5);
            var arg = ret / Math.sqrt(idf);
            var z = Math.atan(arg);
            var s, c;
            for(var i = 0; i < maxit; i++) {
                s = Math.sin(z); c = Math.cos(z);
                z = z - (z + (c + (2.0/3.0)*c*c*c)*s - con)/((8.0/3.0)*Math.pow(c, 4));
            }
            ret = Math.sqrt(idf) * s/c;
        } else if (idf == 6) {
            var con = 2.0 * (p - 0.5);
            var arg = ret / Math.sqrt(idf);
            var z = Math.atan(arg);
            var s, c;
            for (var i = 0; i < maxit; i++) {
                s = Math.sin(z); c = Math.cos(z);
                z = z - ((1.0 + 0.5*c*c + 0.375*Math.pow(c, 4))*s - con)/((15.0/8.0)*Math.pow(c, 5));
            }
            ret = Math.sqrt(idf) * s/c;
        }
    } else {
        ret = 0.0;
    }

    return ret;
}

function ppfnml(p) {
    /* LATEST REVISION  -  2005 April (REM)
     *
     * THIS FUNCTION IS A VERSION OF DATAPAC SUBROUTINE
     * NORPPF, WITH MODIFICATIONS TO FACILITATE CONVERSION TO
     * DOUBLE PRECISION AUTOMATICALLY USING THE NAG, INC. CODE APT, AND
     * TO CORRESPOND TO STARPAC CONVENTIONS.
     *
     * PURPOSE--THIS SUBROUTINE COMPUTES THE PERCENT POINT
     *          FUNCTION VALUE FOR THE NORMAL (GAUSSIAN)
     *          DISTRIBUTION WITH MEAN = 0 AND STANDARD DEVIATION = 1.
     *          THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
     *          THE PROBABILITY DENSITY FUNCTION
     *          F(X) = (1/SQRT(2*PI))*EXP(-X*X/2).
     *          NOTE THAT THE PERCENT POINT FUNCTION OF A DISTRIBUTION
     *          IS IDENTICALLY THE SAME AS THE INVERSE CUMULATIVE
     *          DISTRIBUTION FUNCTION OF THE DISTRIBUTION.
     * ERROR CHECKING--NONE
     * RESTRICTIONS--P SHOULD BE BETWEEN 0.0D0 AND 1.0D0, EXCLUSIVELY.
     * REFERENCES--ODEH AND EVANS, THE PERCENTAGE POINTS
     *             OF THE NORMAL DISTRIBUTION, ALGORTIHM 70,
     *             APPLIED STATISTICS, 1974, PAGES 96-97.
     *           --EVANS, ALGORITHMS FOR MINIMAL DEGREE
     *             POLYNOMIAL AND RATIONAL APPROXIMATION,
     *             M. SC. THESIS, 1972, UNIVERSITY
     *             OF VICTORIA, B. C., CANADA.
     *           --HASTINGS, APPROXIMATIONS FOR DIGITAL
     *             COMPUTERS, 1955, PAGES 113, 191, 192.
     *           --NATIONAL BUREAU OF STANDARDS APPLIED MATHEMATICS
     *             SERIES 55, 1964, PAGE 933, FORMULA 26.2.23.
     *           --FILLIBEN, SIMPLE AND ROBUST LINEAR ESTIMATION
     *             OF THE LOCATION PARAMETER OF A SYMMETRIC
     *             DISTRIBUTION (UNPUBLISHED PH.D. DISSERTATION,
     *             PRINCETON UNIVERSITY), 1969, PAGES 21-44, 229-231.
     *           --FILLIBEN, "THE PERCENT POINT FUNCTION",
     *             (UNPUBLISHED MANUSCRIPT), 1970, PAGES 28-31.
     *           --JOHNSON AND KOTZ, CONTINUOUS UNIVARIATE
     *             DISTRIBUTIONS--1, 1970, PAGES 40-111.
     *           --THE KELLEY STATISTICAL TABLES, 1948.
     *           --OWEN, HANDBOOK OF STATISTICAL TABLES,
     *             1962, PAGES 3-16.
     *           --PEARSON AND HARTLEY, BIOMETRIKA TABLES
     *             FOR STATISTICIANS, VOLUME 1, 1954,
     *             PAGES 104-113.
     * COMMENTS--THE CODING AS PRESENTED BELOW
     *           IS ESSENTIALLY IDENTICAL TO THAT
     *           PRESENTED BY ODEH AND EVANS
     *           AS ALGORTIHM 70 OF APPLIED STATISTICS.
     *           THE PRESENT AUTHOR HAS MODIFIED THE
     *           ORIGINAL ODEH AND EVANS CODE WITH ONLY
     *           MINOR STYLISTIC CHANGES.
     *         --AS POINTED OUT BY ODEH AND EVANS
     *           IN APPLIED STATISTICS,
     *           THEIR ALGORITHM REPRESENTES A
     *           SUBSTANTIAL IMPROVEMENT OVER THE
     *           PREVIOUSLY EMPLOYED
     *           HASTINGS APPROXIMATION FOR THE
     *           NORMAL PERCENT POINT FUNCTION--
     *           THE ACCURACY OF APPROXIMATION
     *           BEING IMPROVED FROM 4.5*(10**-4)
     *           TO 1.5*(10**-8).
     *
     * WRITTEN BY--JAMES J. FILLIBEN
     *             STATISTICAL ENGINEERING LABORATORY
     *             NATIONAL BUREAU OF STANDARDS
     *             WASHINGTON, D. C. 20234
     * ORIGINAL VERSION--JUNE      1972.
     * UPDATED         --SEPTEMBER 1975.
     * UPDATED         --NOVEMBER  1975.
     * UPDATED         --OCTOBER   1976.
     *
     * MODIFIED BY     --JANET R. DONALDSON, DECEMBER 7, 1981
     *                   STATISTICAL ENGINEERING DIVISION
     *                   NATIONAL BUREAU OF STANDARDS, BOULDER, COLORDAO
     *
     * Translated to Javascript by Richard Morse, MGH Biostatistics Center
     *                             April 2005
     */

    var ret = 0.0;

    if (p == 0.5) {
        ret = 0.0;
    } else {
        var r = p;
        if (p > 0.5) { r = 1.0 - r; }
        var t = Math.sqrt(-2.0 * Math.log(r));
        var anum = ((((((((t * PPFNML_P4) + PPFNML_P3) * t) + PPFNML_P2) * t) + PPFNML_P1) * t) + PPFNML_P0);
        var aden = ((((((((t * PPFNML_Q4) + PPFNML_Q3) * t) + PPFNML_Q2) * t )+ PPFNML_Q1) * t) + PPFNML_Q0);
        ret = t + (anum / aden);
        if (p < 0.5) { ret = -ret; }
    }

    return ret;
}


// The following functions were written by David Schoenfeld, Ph. D.
// for his original Fortran version of the program.

function cri(nn, of, al) {
    return -tin(al, nn-of);
}

function pr(x, cr, n, of) {
    return 1.0 - tnc(cr, n - of, Math.sqrt(n) * x);
}

function gg(bet, d, al) {
    return Math.pow((-ppnd7(al) + ppnd7(bet)), 2) / Math.pow(d, 2);
}

function nn(al, bet, x, of) {
    var n1 = Math.max(Math.round(gg(bet, x, al) + 0.5), of + 1);
    if (n1 > 30000) {
        throw "nn: n1 > 30000";
    }
    var cr1 = cri(n1, of, al);
    var cr2 = cr1;
    var n2 = n1;
    var ret = -1;
    var t1 = pr(x, cr1, n1, of);
    var t2 = t1;
    while(ret < 0) {
        if ((t1 < bet) && (bet <= t2)) {
            if (n2 <= (n1 + 1)) {
                ret = n2;
            } else {
                var n3 = Math.floor((n1 + n2)/2.0 + 0.5);
                var cr3 = cri(n3, of, al);
                var t3 = pr(x, cr3, n3, of);
                if (bet <= t3) {
                    n2 = n3;
                    t2 = t3;
                    cr2 = cr3;
                } else {
                    n1 = n3;
                    t1 = t3;
                    cr1 = cr3;
                }
            }
        } else {
            if (t2 < bet) {
                n1 = n2;
                t1 = t2;
                cr1 = cr2;
                n2 = n2 + 6;
                cr2 = cri(n2, of, al);
                t2 = pr(x, cr2, n2, of);
            } else {
                n2 = n1;
                t2 = t1;
                cr2 = cr1;
                if (n1 > (of + 1)) {
                    n1 = Math.max(n1 - 6, of + 1);
                    cr1 = cri(n1, of, al);
                    t1 = pr(x, cr1, n1, of);
                } else {
                    ret = of + 1;
                }
            }
        }
    }
    return ret;
}

function fdelt(cr, bet, n, of) {
    var d1 = (cr + ppnd7(bet)) / Math.sqrt(n);
    var d2 = d1;
    var t1 = pr(d1, cr, n, of);
    var t2 = t1;
    var d = 0;
    while(d == 0) {
        if ((t1 < bet) && (bet <= t2)) {
            if ((d2 <= d1) || (Math.abs(t1-t2) <= .003)) {
                d = d2;
            } else {
                var d3 = (d1 + d2) / 2.0;
                var t3 = pr(d3, cr, n, of);
                if (bet <= t3) {
                    d2 = d3;
                    t2 = t3;
                } else {
                    d1 = d3;
                    t1 = t3;
                }
            }
        } else {
            if (t2 < bet) {
                d1 = d2;
                t1 = t2;
                d2 += 0.2;
                t2 = pr(d2, cr, n, of);
            } else {
                d2 = d1;
                t2 = t1;
                d1 = Math.max((d1 - 0.2), 0.0);
                t1 = pr(d1, cr, n, of);
            }
        }
    }
    return d;
}

function tin(p, didf) {
    return ppft(p, Math.round(didf));
}

function anorin(p) {
    return ppnd7(p);
}

// Following are some helper functions written by Richard Morse


function val_or_defined(v, name) {
    // values which come from form elements are never undefined, so if the value is empty, return undefined
    // we also coerce all values going through this function to be numbers...
    var result = ((v != "") ? Number(v) : undefined);
    if (defined(result) && isNaN(result)) { throw "The value '" + v + "' supplied for " + name + " is not a number"; }
    return result;
}

function invnorm(p) {
    var
	p1 = 2.515517,
	p2 = 0.802853,
	p3 = 0.010328,
	q1 = 1.432788,
	q2 = 0.189269,
	q3 = 0.001308;

    var pp;
    if (p > 0.5) {
	pp = 1 - p;
    } else {
	pp = p;
    }

    var y = Math.sqrt(-2.0 * Math.log(pp));
    var m = p1 + (p2 * y) + (p3 * y * y);
    var n = 1 + (q1 * y) + (q2 * y * y) + (q3 * y * y * y);
    var o = y - (m / n);

    var z;
    if (p < 0.5) {
	z = (- o);
    } else {
	z = o;
    }

    return z;
}

function defined(v) {
  // basically, the equivalent to Perl's "defined"
  return ((typeof(v) != "undefined") ? true : false);
}


function doSizeAssoc(alpha, dep, ind, sampleSize, power, diff_in) {
    var n, prob, diff, t_case, looking_for;

    if (!defined(sampleSize)) {
        looking_for = "total number of patients";
        if (!defined(power)) { throw "\"Power\" must be defined"; }
        if (power <= 0) { throw "\"Power\" should not be less than zero"; }
        if (power >= 0.999999) { throw "\"Power\" is too close to 1"; }

        if (!defined(diff_in)) { throw "\"Minimal detectable difference\" must be defined"; }
        diff = diff_in;

        if (!defined(ind) && !defined(dep)) {
            n = nn(alpha, power, diff, 2);
            prob = pr(diff, cri(n, 2, alpha), n, 2);
            t_case = "!ind_!dep";
        } else if (!defined(ind) && defined(dep)) {
            var dd = diff / dep;
            n = nn(alpha, power, dd, 2);
            prob = pr(dd, cri(n, 2, alpha), n, 2);
            t_case = "!ind_dep";
        } else if (defined(ind) && !defined(dep)) {
            var dd = ind * diff;
            n = nn(alpha, power, dd, 2);
            prob = pr(dd, cri(n, 2, alpha), n, 2);
            t_case = "ind_!dep"
        } else if (defined(ind) && defined(dep)) {
            var dd = ind * diff / dep;
            n = nn(alpha, power, dd, 2);
            prob = pr(dd, cri(n, 2, alpha), n, 2);
            t_case = "ind_dep";
        } else {
            throw "this can't happen!";
        }
    } else if (!defined(power)) {
        looking_for = "detection probability (power)";
        n = sampleSize;

        if (!defined(diff_in)) { throw "\"Minimal detectable difference\" must be defined"; }
        diff = diff_in;

        var cr = cri(n, 2, alpha);
        if (!defined(ind) && !defined(dep)) {
            prob = pr(diff, cr, n, 2);
            t_case = "!ind_!dep";
        } else if (!defined(ind) && defined(dep)) {
            prob = pr((diff / dep), cr, n, 2);
            t_case = "!ind_dep";
        } else if (defined(ind) && !defined(dep)) {
            prob = pr((diff * ind), cr, n, 2);
            t_case = "ind_!dep"
        } else if (defined(ind) && defined(dep)) {
            prob = pr((ind * diff / dep), cr, n, 2);
            t_case = "ind_dep";
        } else {
            throw "this can't happen!";
        }
    } else if (!defined(diff_in)) {
        looking_for = "minimal detectable difference";
        n = sampleSize;
        prob = power;

        var cr = cri(n, 2, alpha);
        var del = fdelt(cr, power, n, 2);

        if (!defined(ind) && !defined(dep)) {
            diff = del;
            t_case = "!ind_!dep";
        } else if (!defined(ind) && defined(dep)) {
            diff = del * dep;
            t_case = "!ind_dep";
        } else if (defined(ind) && !defined(dep)) {
            diff = del / ind;
            t_case = "ind_!dep"
        } else if (defined(ind) && defined(dep)) {
            diff = del * dep / ind;
            t_case = "ind_dep";
        } else {
            throw "this can't happen!";
        }
    } else {
        throw "You need to leave one of the three inputs (number of patients, power of the study, minimal detectable difference) empty, so that it can be calculated.";
    }

    return { n: n, prob: Math.floor(prob * 100), diff: diff, t_case: t_case, looking_for: looking_for, params: { alpha: alpha, dep: dep, ind: ind, power: power, sampleSize: sampleSize, diff_in: diff_in } };
  }

module.exports = {
  doSizeAssoc
}