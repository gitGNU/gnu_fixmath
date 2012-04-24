/*  Copyright (C) 2005-2011, Axis Communications AB, LUND, SWEDEN
 *
 *  This file is part of Fixmath.
 *
 *  Fixmath is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  You can use the comments under either the terms of the GNU Lesser General
 *  Public License version 3 as published by the Free Software Foundation,
 *  either version 3 of the License or (at your option) any later version, or
 *  the GNU Free Documentation License version 1.3 or any later version
 *  published by the Free Software Foundation; with no Invariant Sections, no
 *  Front-Cover Texts, and no Back-Cover Texts.
 *  A copy of the license is included in the documentation section entitled
 *  "GNU Free Documentation License".
 *
 *  Fixmath is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License and a copy of the GNU Free Documentation License along
 *  with Fixmath. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 *  @file   fixtest.c
 *  @brief  Correctness tests for fixed-point math library.
 */

#include <stdbool.h> /* Boolean types      */
#include <stdint.h>  /* uint8_t            */
#include <stdlib.h>  /* rand()             */
#include <stdio.h>   /* printf()           */
#include <string.h>  /* strcmp(), strlen() */
#include <errno.h>   /* errno facility     */
#include <assert.h>  /* assert()           */
#include <math.h>    /* roundf(), ...      */
#include <float.h>   /* DBL_EPSILON        */
#include <check.h>   /* Unit test fwk      */
#include "fixmath.h" /* Fixed-point API    */


/*
 * -------------------------------------------------------------
 *  Macros
 * -------------------------------------------------------------
 */

/**
 *  Standard MIN().
 */
#undef  MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/**
 *  Standard MAX().
 */
#undef  MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/*
 * -------------------------------------------------------------
 *  Constants
 * -------------------------------------------------------------
 */

#define NTESTS (1 << 14)

/*
 * -------------------------------------------------------------
 *  Local functions fwd declare
 * -------------------------------------------------------------
 */

static Suite*
fx_test_fixmath(void);

static bool
fx_test_driver(uint32_t (*test)(), uint32_t (*ref)(), bool zero);

static bool
fx_test_ftox_xtof(void);

static bool
fx_test_dtox_xtod(void);

static bool
fx_test_itox_xtoi(void);

static bool
fx_test_xtox(void);

static bool
fx_test_rounding(void);

static bool
fx_test_mulx(void);

static bool
fx_test_divx(void);

static bool
fx_test_rdivx(void);

static bool
fx_test_sqrtx(void);

static bool
fx_test_isqrtx(void);

static bool
fx_test_invx(void);

static bool
fx_test_exp_driver(fixed_t (*test)(fixed_t, unsigned),
                   double (*ref)(double),
                   double (*inv)(double));
static bool
fx_test_log_driver(fixed_t (*test)(fixed_t, unsigned),
                   double (*ref)(double),
                   double (*inv)(double));
static bool
fx_test_powx(void);

static bool
fx_test_sincos_driver(fixed_t (*test)(fixed_t, unsigned),
                      double (*ref)(double));
static bool
fx_test_itoa(void);

static bool
fx_test_xtoa(void);

static uint32_t
fx_test_clz_macro(uint32_t word);

static uint32_t
fx_test_clz_func(uint32_t word);

static uint32_t
fx_test_clz_impl(uint32_t word);

static uint32_t
fx_test_ctz_macro(uint32_t word);

static uint32_t
fx_test_ctz_func(uint32_t word);

static uint32_t
fx_test_ctz_impl(uint32_t word);

static uint32_t
fx_test_bitcount_macro(uint32_t word);

static uint32_t
fx_test_bitcount_func(uint32_t word);

static uint32_t
fx_test_bitcount_impl(uint32_t word);

static uint32_t
fx_test_clz_ref(uint32_t word);

static uint32_t
fx_test_ctz_ref(uint32_t word);

static uint32_t
fx_test_bitcount_ref(uint32_t word);

static uint32_t
fx_test_urand(void);

static double
fx_test_drand(double min, double max);

static double
fx_test_exp2_ref(double x);

static double
fx_test_exp10_ref(double x);

double
log2(double x);


/*
 * -------------------------------------------------------------
 *  Test cases
 * -------------------------------------------------------------
 */

START_TEST(tcase_clz_macro)
{
    fail_unless(fx_test_driver(&fx_test_clz_macro, &fx_test_clz_ref, false));
}
END_TEST

START_TEST(tcase_clz_func)
{
    fail_unless(fx_test_driver(&fx_test_clz_func, &fx_test_clz_ref, false));
}
END_TEST

START_TEST(tcase_clz_impl)
{
    fail_unless(fx_test_driver(&fx_test_clz_impl, &fx_test_clz_ref, false));
}
END_TEST

START_TEST(tcase_ctz_macro)
{
    fail_unless(fx_test_driver(&fx_test_ctz_macro, &fx_test_ctz_ref, false));
}
END_TEST

START_TEST(tcase_ctz_func)
{
    fail_unless(fx_test_driver(&fx_test_ctz_func, &fx_test_ctz_ref, false));
}
END_TEST

START_TEST(tcase_ctz_impl)
{
    fail_unless(fx_test_driver(&fx_test_ctz_impl, &fx_test_ctz_ref, false));
}
END_TEST

START_TEST(tcase_bitcount_macro)
{
    fail_unless(fx_test_driver(&fx_test_bitcount_macro,
                               &fx_test_bitcount_ref, true));
}
END_TEST

START_TEST(tcase_bitcount_func)
{
    fail_unless(fx_test_driver(&fx_test_bitcount_func,
                               &fx_test_bitcount_ref, true));
}
END_TEST

START_TEST(tcase_bitcount_impl)
{
    fail_unless(fx_test_driver(&fx_test_bitcount_impl,
                               &fx_test_bitcount_ref, true));
}
END_TEST

START_TEST(tcase_ftox_xtof)
{
    fail_unless(fx_test_ftox_xtof());
}
END_TEST

START_TEST(tcase_dtox_xtod)
{
    fail_unless(fx_test_dtox_xtod());
}
END_TEST

START_TEST(tcase_itox_xtoi)
{
    fail_unless(fx_test_itox_xtoi());
}
END_TEST

START_TEST(tcase_xtox)
{
    fail_unless(fx_test_xtox());
}
END_TEST

START_TEST(tcase_rounding)
{
    fail_unless(fx_test_rounding());
}
END_TEST

START_TEST(tcase_mulx)
{
    fail_unless(fx_test_mulx());
}
END_TEST

START_TEST(tcase_divx)
{
    fail_unless(fx_test_divx());
}
END_TEST

START_TEST(tcase_rdivx)
{
    fail_unless(fx_test_rdivx());
}
END_TEST

START_TEST(tcase_sqrtx)
{
    fail_unless(fx_test_sqrtx());
}
END_TEST

START_TEST(tcase_isqrtx)
{
    fail_unless(fx_test_isqrtx());
}
END_TEST

START_TEST(tcase_invx)
{
    fail_unless(fx_test_invx());
}
END_TEST

START_TEST(tcase_expx)
{
    fail_unless(fx_test_exp_driver(&fx_expx, &exp, &log));
}
END_TEST

START_TEST(tcase_exp2x)
{
    fail_unless(fx_test_exp_driver(&fx_exp2x, &fx_test_exp2_ref, &log2));
}
END_TEST

START_TEST(tcase_exp10x)
{
    fail_unless(fx_test_exp_driver(&fx_exp10x, &fx_test_exp10_ref, &log10));
}
END_TEST

START_TEST(tcase_logx)
{
    fail_unless(fx_test_log_driver(&fx_logx, &log, &exp));
}
END_TEST

START_TEST(tcase_log2x)
{
    fail_unless(fx_test_log_driver(&fx_log2x, &log2, &fx_test_exp2_ref));
}
END_TEST

START_TEST(tcase_log10x)
{
    fail_unless(fx_test_log_driver(&fx_log10x, &log10, &fx_test_exp10_ref));
}
END_TEST

START_TEST(tcase_powx)
{
    fail_unless(fx_test_powx());
}
END_TEST

START_TEST(tcase_sinx)
{
    fail_unless(fx_test_sincos_driver(&fx_sinx, &sin));
}
END_TEST

START_TEST(tcase_cosx)
{
    fail_unless(fx_test_sincos_driver(&fx_cosx, &cos));
}
END_TEST

START_TEST(tcase_itoa)
{
    fail_unless(fx_test_itoa());
}
END_TEST

START_TEST(tcase_xtoa)
{
    fail_unless(fx_test_xtoa());
}
END_TEST


/*
 * -------------------------------------------------------------
 *  Exported functions
 * -------------------------------------------------------------
 */

int main(void)
{
    int num;

    /* Create the Suite Runner object */
    SRunner *sr = srunner_create(fx_test_fixmath());

    /* Run the tests */
    srunner_run_all(sr, CK_ENV);
    num = srunner_ntests_failed(sr);
    srunner_free(sr);

    return num ? EXIT_FAILURE : EXIT_SUCCESS;
}

/*
 * -------------------------------------------------------------
 *  Local functions
 * -------------------------------------------------------------
 */

static Suite*
fx_test_fixmath(void)
{
    Suite *s  = suite_create("Fixed-point math");
    TCase *tc;

    tc = tcase_create("fx_clz (inline)");
    tcase_add_test(tc, tcase_clz_macro);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_clz (call)");
    tcase_add_test(tc, tcase_clz_func);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_clz (impl)");
    tcase_add_test(tc, tcase_clz_impl);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_ctz (inline)");
    tcase_add_test(tc, tcase_ctz_macro);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_ctz (call)");
    tcase_add_test(tc, tcase_ctz_func);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_ctz (impl)");
    tcase_add_test(tc, tcase_ctz_impl);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_bitcount (inline)");
    tcase_add_test(tc, tcase_bitcount_macro);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_bitcount (call)");
    tcase_add_test(tc, tcase_bitcount_func);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_bitcount (impl)");
    tcase_add_test(tc, tcase_bitcount_impl);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_ftox and fx_xtof");
    tcase_add_test(tc, tcase_ftox_xtof);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_dtox and dx_xtof");
    tcase_add_test(tc, tcase_dtox_xtod);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_itox_xtoi");
    tcase_add_test(tc, tcase_itox_xtoi);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_xtox");
    tcase_add_test(tc, tcase_xtox);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_rounding");
    tcase_add_test(tc, tcase_rounding);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_mulx");
    tcase_add_test(tc, tcase_mulx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_divx");
    tcase_add_test(tc, tcase_divx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_rdivx");
    tcase_add_test(tc, tcase_rdivx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_sqrtx");
    tcase_add_test(tc, tcase_sqrtx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_isqrtx");
    tcase_add_test(tc, tcase_isqrtx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_invx");
    tcase_add_test(tc, tcase_invx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_expx");
    tcase_add_test(tc, tcase_expx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_exp2x");
    tcase_add_test(tc, tcase_exp2x);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_exp10x");
    tcase_add_test(tc, tcase_exp10x);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_logx");
    tcase_add_test(tc, tcase_logx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_log2x");
    tcase_add_test(tc, tcase_log2x);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_log10x");
    tcase_add_test(tc, tcase_log10x);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_powx");
    tcase_add_test(tc, tcase_powx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_sinx");
    tcase_add_test(tc, tcase_sinx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_cosx");
    tcase_add_test(tc, tcase_cosx);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_itoa");
    tcase_add_test(tc, tcase_itoa);
    suite_add_tcase(s, tc);

    tc = tcase_create("fx_xtoa");
    tcase_add_test(tc, tcase_xtoa);
    suite_add_tcase(s, tc);

    return s;
}

static bool
fx_test_driver(uint32_t (*test)(), uint32_t (*ref)(), bool zero)
{
    static const uint32_t mtab[] =
        {0xffffffffUL, 0xffff0000UL, 0x0000ffffUL, 0x00ffff00UL,
         0xff000000UL, 0x00ff0000UL, 0x0000ff00UL, 0x0000ff00UL,
         0xff00ff00UL, 0x00ff00ffUL, 0xff0000ffUL, 0xf0f0f0f0UL};

    int k;

    /* Perform explicit zero test */
    if (zero && (*test)(0) != (*ref)(0)) {
        return false;
    }

    /* Test random words */
    for (k = 0; k < NTESTS; k++) {
        uint32_t word = fx_test_urand();
        int      i;

        for (i = 0; i < (int)(sizeof mtab / sizeof mtab[0]); i++) {
            uint32_t w = word & mtab[i];

            if ((zero || w) && (*test)(w) != (*ref)(w)) {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_ftox_xtof(void)
{
    int k;

    for (k = 0; k < NTESTS;  k++) {
        float fval = fx_test_drand(-0.5, 0.5);
        int   frac;

        for (frac = 0; frac < 32; frac++) {
            float   tol = 0.5f*(1 + (frac >= 24));
	    float   err;
	    fixed_t xval;


            xval = fx_ftox(fval, frac);
            err  = fval - fx_xtof(xval, frac);

            /* Check that the numeric result is correct */
            if (fabs(err)*(1UL << frac) > tol) {
                return false;
            }

            /* Check that the inline and call versions are equal */
            if (fx_ftox(fval, frac) != (fx_ftox)(fval, frac) ||
                fx_xtof(xval, frac) != (fx_xtof)(xval, frac))
            {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_dtox_xtod(void)
{
    int k;

    for (k = 0; k < NTESTS;  k++) {
        double dval = fx_test_drand(-0.5, 0.5);
        int    frac;

        for (frac = 0; frac < 32; frac++) {
            fixed_t xval;
            double  err;

            /* Check that the numeric result is correct */
            xval = fx_dtox(dval, frac);
            err  = dval - fx_xtod(xval, frac);

            if (fabs(err)*(1UL << frac) > .5) {
                return false;
            }

            /* Check that the inline and call versions are equal */
            if (fx_dtox(dval, frac) != (fx_dtox)(dval, frac) ||
                fx_xtod(xval, frac) != (fx_xtod)(xval, frac))
            {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_itox_xtoi(void)
{
    int k;
    for (k = 0; k < NTESTS;  k++) {

        int frac;

        for (frac = 0; frac < 32; frac++) {
            int     rng  = (1UL << (31 - frac)) - 1;
            double  rval = fx_test_drand(-rng, rng);
            fixed_t xval = fx_dtox(rval, frac);
            double  dval = fx_xtod(xval, frac);
            int     ival = fx_xtoi(xval, frac);

            /* Check fx_xtoi() numeric result */
            if (ival != floor(dval)) {
                return false;
            }

            /* Check fx_itox() numeric result */
            if (fx_xtoi(fx_itox(ival, frac), frac) != ival) {
                return false;
            }

            /* Check inline vs call equivalence */
            if ((fx_xtoi)(xval, frac) != fx_xtoi(xval, frac) ||
                (fx_itox)(ival, frac) != fx_itox(ival, frac))
            {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_xtox(void)
{
    int ntests = (int)sqrt(NTESTS);
    int k;

    for (k = 0; k < ntests;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double  rval = fx_test_drand(-0.5, 0.5);
            fixed_t xval = fx_dtox(rval, frac);
            int     f2;

            for (f2 = 0; f2 < 32; f2++) {
                fixed_t x2  = fx_xtox(xval, frac, f2);
                fixed_t err = xval - fx_xtox(x2, f2, frac);

                /* Check the numeric result */
                if (fabs(err) > fx_dtox(0.5, frac)) {
                    return false;
                }

                /* Check inline vs call equivalence */
                if (fx_xtox(xval, frac, f2) != (fx_xtox)(xval, frac, f2)) {
                    return false;
                }
            }
        }
    }

    return true;
}

static bool
fx_test_rounding(void)
{
    int k;
    for (k = 0; k < NTESTS;  k++) {

        int frac;

        for (frac = 0; frac < 31; frac++) {
            int     rng  = (1UL << (31 - frac)) - 1;
            double  rval = fx_test_drand(-rng, rng);
            fixed_t xval = fx_dtox(rval, frac);
            fixed_t err;

            /* Check fx_roundx() numeric value */
            err = xval - fx_itox(fx_roundx(xval, frac), frac);
            if (abs(err) > fx_dtox(0.5, frac)) {
                return false;
            }

            /* Check fx_floorx() numeric value */
            err = xval - fx_itox(fx_floorx(xval, frac), frac);
            if (err < 0 || err > fx_itox(1, frac)) {
                return false;
            }

            /* Check fx_ceilx() numeric value */
            err = xval - fx_itox(fx_ceilx(xval, frac), frac);
            if (err > 0 || err < -fx_itox(1, frac)) {
                return false;
            }

            /* Check inline vs call equivalence */
            if (fx_roundx(xval, frac) != (fx_roundx)(xval, frac) ||
                fx_floorx(xval, frac) != (fx_floorx)(xval, frac) ||
                fx_ceilx (xval, frac) != (fx_ceilx) (xval, frac))
            {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_mulx(void)
{
    int k;
    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double  r1   = fx_test_drand(-0.9, 0.9);
            double  r2   = fx_test_drand(-0.9, 0.9);
            fixed_t x1   = fx_dtox(r1, frac);
            fixed_t x2   = fx_dtox(r2, frac);
            double  d1   = fx_xtod(x1, frac);
            double  d2   = fx_xtod(x2, frac);

            /* Check numeric result */
            if (fabs(fx_xtod(fx_mulx(x1, x2, frac), frac) - d1*d2) >
                0.5*fx_xtod(1, frac))
            {
                return false;
            }

            /* Check inline vs call equivalence */
            if (fx_mulx(x1, x2, frac) != (fx_mulx)(x1, x2, frac)) {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_divx(void)
{
    int k;
    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double  r1   = fx_test_drand(-0.5, 0.5);
            double  r2   = fx_test_drand(0.0, 1.0) > 0.5 ?
                           fx_test_drand(0.5, 0.9) :
                          -fx_test_drand(0.5, 0.9);
            fixed_t x1   = fx_dtox(r1, frac);
            fixed_t x2   = fx_dtox(r2, frac);
            double  d1   = fx_xtod(x1, frac);
            double  d2   = fx_xtod(x2, frac);

            /* Check numeric result */
            if (fabs(fx_xtod(fx_divx(x1, x2, frac), frac) - d1/d2) >
                0.5*fx_xtod(1, frac))
            {
                return false;
            }

            /* Check inline vs call equivalence */
            if (fx_divx(x1, x2, frac) != (fx_divx)(x1, x2, frac)) {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_rdivx(void)
{
    fixed_t   val  =  0x08765432;
    fx_rdiv_t rdiv = {0x12345678, 31};

    /* Check inline vs call equivalence.
     * The numeric value is tested in fx_isqrtx and fx_invx test cases.
     */
    return fx_rdivx(val, &rdiv) == (fx_rdivx)(val, &rdiv);
}

static bool
fx_test_sqrtx(void)
{
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    fx_sqrtx(-1, 16);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_sqrtx(1, 33);
    if (errno != EDOM) {
        return false;
    }

    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double  rval  = fx_test_drand(0.0, fx_xtod(INT32_MAX, frac));
            fixed_t xval  = fx_dtox(rval, frac);
            double  dval  = fx_xtod(xval, frac);

            /* Check numeric result */
            if (abs(fx_sqrtx(xval, frac) - fx_dtox(sqrt(dval), frac)) >
                2*(1 + (frac == 31)))
            {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_isqrtx(void)
{
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    fx_isqrtx(-1, 16, NULL);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_isqrtx(0, 16, NULL);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_isqrtx(1, 33, NULL);
    if (errno != EDOM) {
        return false;
    }

    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double    r1   = fx_test_drand(-0.1, 0.1);
            double    r2   = fx_test_drand(0.1, 0.9);
            fixed_t   x1   = fx_dtox(r1, frac);
            fixed_t   x2   = fx_dtox(r2, frac);

            if (x2 != 0) {
                double    d2 = fx_xtod(x2, frac);
                fixed_t   isq;
                fx_rdiv_t rdiv;

                /* Compute the inverse square root */
                isq = fx_isqrtx(x2, frac, &rdiv);

                /* Check the inverted numeric value */
                if (frac < 31 &&
                    fabs(1.0/sqrt(d2))*(double)(1UL << frac) < INT32_MAX &&
                    abs(isq - fx_rdivx(fx_itox(1, frac), &rdiv)) > 1)
                {
                    return false;
                }

                /* Check the reciprocal multiplier */
                if (abs(fx_rdivx(x1, &rdiv) -
                        fx_divx(x1, fx_sqrtx(x2, frac), frac)) >
                    2 + (frac == 31))
                {
                    return false;
                }
            }
        }
    }

    return true;
}

static bool
fx_test_invx(void)
{
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    fx_invx(0, 16, NULL);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_invx(1, 33, NULL);
    if (errno != EDOM) {
        return false;
    }

    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double    r1   = fx_test_drand(-0.1, 0.1);
            double    r2   = fx_test_drand(0.0, 1.0) > 0.5 ?
                             fx_test_drand(0.1, 0.9) :
                            -fx_test_drand(0.1, 0.9);
            fixed_t   x1   = fx_dtox(r1, frac);
            fixed_t   x2   = fx_dtox(r2, frac);

            if (x2 != 0) {
                double    d2 = fx_xtod(x2, frac);
                fixed_t   inv;
                fx_rdiv_t rdiv;

                /* Compute the inverse */
                inv = fx_invx(x2, frac, &rdiv);

                /* Check the inverted numeric value */
                if (frac < 31 &&
                    fabs(1.0/d2)*(double)(1UL << frac) < INT32_MAX &&
                    abs(inv - fx_rdivx(fx_itox(1, frac), &rdiv)) > 1)
                {
                    return false;
                }

                /* Check the reciprocal multiplier */
                if (abs(fx_rdivx(x1, &rdiv) - fx_divx(x1, x2, frac)) >
                    4 + 3*(frac == 31))
                {
                    return false;
                }
            }
        }
    }

    return true;
}

static bool
fx_test_exp_driver(fixed_t (*test)(fixed_t, unsigned),
                   double (*ref)(double),
                   double (*inv)(double))
{
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    (*test)(1, 32);
    if (errno != EDOM) {
        return false;
    }

    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double   rng  = 0.5*(1UL << (31 - frac));
            double   rval = fx_test_drand(DBL_EPSILON, (*inv)(rng));
            fixed_t  xval = fx_dtox(rval, frac);
            double   dval = fx_xtod(xval, frac);

            /* Test the exponential value */
            if (abs((*test)(xval, frac) - fx_dtox((*ref)(dval), frac)) > 2) {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_log_driver(fixed_t (*test)(fixed_t, unsigned),
                   double (*ref)(double),
                   double (*inv)(double))
{
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    (*test)(1, 32);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    (*test)(0, 0);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    (*test)(-1, 0);
    if (errno != EDOM) {
        return false;
    }

    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double   rng  = 0.5*(1UL << (31 - frac));
            double   rval = fx_test_drand((*inv)(-rng), rng);
            fixed_t  xval = fx_dtox(rval, frac);
            double   dval = fx_xtod(xval, frac);

            /* Test the logarithm value */
            if (abs((*test)(xval, frac) - fx_dtox((*ref)(dval), frac)) > 2) {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_powx(void)
{
    int ntests = sqrt(NTESTS);
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    fx_powx(1, 32, 1, 0);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_powx(1, 0, 1, 32);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_powx(-1, 0, 1, 1);
    if (errno != EDOM) {
        return false;
    }
    errno = 0;
    fx_powx(0, 0, -1, 1);
    if (errno != EDOM) {
        return false;
    }

    /* Check edge cases */
    fail_unless(fx_powx(0, 1, 5, 1) == 0, "0 ^ x != 0");

    for (k = 0; k < ntests;  k++) {
        int frac1, frac2;

        /* Check worst case error bounds: (1 + epsilon)^y for large y */
        for (frac1 = 0; frac1 < 32; frac1++) {
            int     offset = rand() > RAND_MAX / 2 ? 1 : -1;
            fixed_t xval1  = frac1 == 31 ? INT32_MAX :
                                           fx_itox(1, frac1) + offset;
            double  dval1  = fx_xtod(xval1, frac1);

            for (frac2 = 0; frac2 < 32; frac2++) {
                double  max1  = fx_xtod(INT32_MAX, frac1);
                double  max2  = fx_xtod(INT32_MAX, frac2);
                double  dval2 = MIN(log(max1) / log(xval1), max2);
                fixed_t xval2 = fx_dtox(dval2, frac2);
                double  value = pow(dval1, dval2);
                fixed_t error = abs(fx_powx(xval1, frac1, xval2, frac2) -
                                    fx_dtox(value, frac1));

                /* Check error bounds */
                if (error > fabs(dval2)*value + 5) {
                    printf("error=%d ulp\n", error);
                    return false;
                }
            }
        }

        /* Run random tests */
        for (frac1 = 0; frac1 < 32; frac1++) {
            double  rng1  = 0.5*(1UL << (31 - frac1));
            double  rval1 = fx_test_drand(0.0, rng1);
            fixed_t xval1 = MAX(fx_dtox(rval1, frac1), 1);
            double  dval1 = fx_xtod(xval1, frac1);
            double  lim2  = log(rng1) / log(dval1);

            for (frac2 = 0; frac2 < 32; frac2++) {
                double  rng2   = 0.5*(1UL << (31 - frac2));
                double  rval2  = dval1 > 1.0 ?
                                 fx_test_drand(-rng2, MIN(rng2, lim2)) :
                                 fx_test_drand(MAX(-rng2, lim2), rng2);
                fixed_t xval2  = fx_dtox(rval2, frac2);
                double  dval2  = fx_xtod(xval2, frac2);
                double  value  = pow(dval1, dval2);
                fixed_t error  = abs(fx_powx(xval1, frac1, xval2, frac2) -
                                     fx_dtox(value, frac1));

                /* Check error bounds */
                if (error > fabs(dval2)*value + 5) {
                    printf("error=%d ulp\n", error);
                    return false;
                }
            }
        }
    }

    return true;
}

static bool
fx_test_sincos_driver(fixed_t (*test)(fixed_t, unsigned),
                      double (*ref)(double))
{
    int k;

    /* Check that errno is set to EDOM on invalid input */
    errno = 0;
    (*test)(1, 32);
    if (errno != EDOM) {
        return false;
    }

    for (k = 0; k < NTESTS;  k++) {
        int frac;

        for (frac = 0; frac < 32; frac++) {
            double   rng  = 0.5*(1UL << (31 - frac));
            double   rval = fx_test_drand(0.0, 1.0) > 0.5 ?
                            fx_test_drand(0,  rng) :
                            fx_test_drand(0, -rng);
            fixed_t  xval = fx_dtox(rval, frac);
            double   dval = fx_xtod(xval, frac);

            /* Test the sin/cos value */
            if (abs((*test)(xval, frac) - fx_dtox((*ref)(dval), frac)) >
                3 + (frac == 31))
            {
                return false;
            }
        }
    }

    return true;
}

static bool
fx_test_itoa(void)
{
    int k, i;

    /* Check the error handling */
    errno = 0;
    if (fx_itoa(NULL, 123) >= 0 ||
        errno != EINVAL)
    {
        return false;
    }

    /* Check operation */
    for (k = 0; k < 32; k++) {
        uint32_t mask = ~0UL >> k;

        for (i = 0; i < NTESTS/32; i++) {
            static char buf[32];
            static char ref[32];
            int32_t     val;
            int         len;

            /* Compute the integer value */
            val = fx_test_urand() & mask;
            if (rand() < RAND_MAX/2) {
                val = -val;
            }

            /* Perform string conversion */
            len = fx_itoa(buf, val);
            sprintf(ref, "%d", val);

            /* Check the result */
            if (len != (int)strlen(buf) ||
                strcmp(buf, ref) != 0)
            {
                return false;
            }

        }
    }

    return true;
}

static bool
fx_test_xtoa(void)
{
    static char buf[32];
    int i, k;

    /* Check the error handling */
    errno = 0;
    if (fx_xtoa(NULL, 123, 0, 0) >= 0 ||
        errno != EINVAL)
    {
        return false;
    }
    errno = 0;
    if (fx_xtoa(buf, 123, 32, 0) >= 0 ||
        errno != EINVAL)
    {
        return false;
    }


    for (i = 0; i < NTESTS/32/32; i++) {
        for (k = 0; k < 32; k++) {
            fixed_t xval;
            int     frac;

            /* Compute the fixed-point value */
            xval = fx_test_urand() & (~0UL >> k);
            if (rand() < RAND_MAX/2) {
                xval = -xval;
            }

            for (frac = 0; frac < 32; frac++) {
                double dval = fx_xtod(xval, frac);
                int    digs;

                for (digs = 0; digs <= 10; digs++) {
                    static char ref[32];
                    static char fmt[8];
                    double      tval;
                    int         expn;
                    int         len;

                    /* Construct the format string */
                    sprintf(fmt, "%%.%de", digs);

                     /* Perform string conversion */
                    len  = fx_xtoa(buf, xval, frac, digs);
                    tval = atof(buf);

                    /* Perform reference string conversion */
                    sprintf(ref, fmt, dval);

                    /* Fetch the exponent value for the last digit */
                    expn = atoi(&ref[strlen(ref) - 3]) - digs;

                    /* Check the return value */
                    if (len != (int)strlen(buf)) {
                        printf("Invalid string length: %d(%d)\n",
                               len, (int)strlen(buf));
                        return false;
                    }

                    /* Check the numeric value */
                    if (fabs(tval - dval)/fx_test_exp10_ref(expn) > 1.0 &&
                        fabs(tval - dval)/fx_xtod(1, frac)        > 1.0)
                    {
                        printf("Wrote     %s\n", buf);
                        printf("Should be %s\n", ref);
                        return false;
                    }
                }
            }

        }
    }

    return true;
}

static uint32_t
fx_test_clz_macro(uint32_t word)
{
    return fx_clz(word);
}

static uint32_t
fx_test_clz_func(uint32_t word)
{
    return (fx_clz)(word);
}

static uint32_t
fx_test_clz_impl(uint32_t word)
{
    uint32_t ret;
    FX_IMPL_CLZ(ret, word);
    return ret;
}

static uint32_t
fx_test_ctz_macro(uint32_t word)
{
    return fx_ctz(word);
}

static uint32_t
fx_test_ctz_func(uint32_t word)
{
    return (fx_ctz)(word);
}

static uint32_t
fx_test_ctz_impl(uint32_t word)
{
    uint32_t ret;
    FX_IMPL_CTZ(ret, word);
    return ret;
}

static uint32_t
fx_test_bitcount_macro(uint32_t word)
{
    return fx_bitcount(word);
}

static uint32_t
fx_test_bitcount_func(uint32_t word)
{
    return (fx_bitcount)(word);
}

static uint32_t
fx_test_bitcount_impl(uint32_t word)
{
    uint32_t ret;
    FX_IMPL_BITCOUNT(ret, word);
    return ret;
}

static uint32_t
fx_test_clz_ref(uint32_t word)
{
    uint32_t bit;
    int      cnt = 0;

    for (bit = 1UL << 31; bit; bit >>= 1) {
        if (word & bit) {
            break;
        }
        cnt++;
    }

    return cnt;
}

static uint32_t
fx_test_ctz_ref(uint32_t word)
{
    uint32_t bit;
    int      cnt = 0;

    for (bit = 1; bit; bit <<= 1) {
        if (word & bit) {
            break;
        }
        cnt++;
    }

    return cnt;
}

static uint32_t
fx_test_bitcount_ref(uint32_t word)
{
    uint32_t bit;
    int      cnt = 0;

    for (bit = 1; bit; bit <<= 1) {
        if (word & bit) {
            cnt++;
        }
    }

    return cnt;
}

static uint32_t
fx_test_urand(void)
{
    return rand() | (rand() << 16);
}

static double
fx_test_drand(double min, double max)
{
    return rand()*(max - min) / (double)RAND_MAX + min;
}

static double
fx_test_exp2_ref(double x)
{
    return exp(x*log(2));
}

static double
fx_test_exp10_ref(double x)
{
    return exp(x*log(10));
}
