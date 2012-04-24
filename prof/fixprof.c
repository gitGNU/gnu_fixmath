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
 *  @file   fixprof.c
 *  @brief  Fixed-point math library profiling.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include "fixmath.h"

/*
 * -------------------------------------------------------------
 *  Constants
 * -------------------------------------------------------------
 */

/**
 *  The number of milliseconds per test.
 */
#define FX_PROF_MSECS 100

/**
 *  The length of the data buffers.
 */
#define FX_PROF_LEN 16

/*
 * -------------------------------------------------------------
 *  Macros
 * -------------------------------------------------------------
 */

#define FX_PROF_EXEC(stmt)                                \
do {                                                      \
    for (fx_prof_iter = 0; fx_prof_run; fx_prof_iter++) { \
        int i;                                            \
        for (i = 0; i < FX_PROF_LEN; i++) {               \
            stmt;                                         \
        }                                                 \
    }                                                     \
} while (0)

#define FX_PROF_ENTRY(func) \
    {"fx_" #func, fx_prof_ ## func}

#define FX_PROF_FLOAT(func, desc) \
    {desc, fx_prof_ ## func}

/*
 * -------------------------------------------------------------
 *  Local functions fwd declare
 * -------------------------------------------------------------
 */

static void
fx_prof_clz(void);

static void
fx_prof_ctz(void);

static void
fx_prof_bitcount(void);

static void
fx_prof_itox(void);

static void
fx_prof_ftox(void);

static void
fx_prof_dtox(void);

static void
fx_prof_xtox(void);

static void
fx_prof_xtoi(void);

static void
fx_prof_xtof(void);

static void
fx_prof_xtod(void);

static void
fx_prof_itoa(void);

static void
fx_prof_xtoa(void);

static void
fx_prof_roundx(void);

static void
fx_prof_ceilx(void);

static void
fx_prof_floorx(void);

static void
fx_prof_addx(void);

static void
fx_prof_subx(void);

static void
fx_prof_mulx(void);

static void
fx_prof_divx(void);

static void
fx_prof_rdivx(void);

static void
fx_prof_invx(void);

static void
fx_prof_sqrtx(void);

static void
fx_prof_isqrtx(void);

static void
fx_prof_expx(void);

static void
fx_prof_exp2x(void);

static void
fx_prof_exp10x(void);

static void
fx_prof_logx(void);

static void
fx_prof_log2x(void);

static void
fx_prof_log10x(void);

static void
fx_prof_powx(void);

static void
fx_prof_sinx(void);

static void
fx_prof_cosx(void);

static void
fx_prof_fadd(void);

static void
fx_prof_fsub(void);

static void
fx_prof_fmul(void);

static void
fx_prof_fdiv(void);

static void
fx_prof_sqrt(void);

static void
fx_prof_exp(void);

static void
fx_prof_log(void);

static void
fx_prof_pow(void);

static void
fx_prof_sin(void);

static void
fx_prof_cos(void);

static void
fx_prof_sprintfd(void);

static void
fx_prof_sprintfe(void);

static void
fx_prof_setf(float *buf, int len, float value);

static void
fx_prof_setd(double *buf, int len, double value);

static void
fx_prof_sighandler(int signum);

static void
fx_prof_print(float value);


/*
 * -------------------------------------------------------------
 *  Global variables
 * -------------------------------------------------------------
 */

struct fx_prof_table_st {
    const char  *name;
    void       (*func)(void);
} fx_prof_table[] = {
    FX_PROF_ENTRY(clz),
    FX_PROF_ENTRY(ctz),
    FX_PROF_ENTRY(bitcount),
    FX_PROF_ENTRY(itox),
    FX_PROF_ENTRY(ftox),
    FX_PROF_ENTRY(dtox),
    FX_PROF_ENTRY(xtox),
    FX_PROF_ENTRY(xtoi),
    FX_PROF_ENTRY(xtof),
    FX_PROF_ENTRY(xtod),
    FX_PROF_ENTRY(itoa),
    FX_PROF_ENTRY(xtoa),
    FX_PROF_ENTRY(roundx),
    FX_PROF_ENTRY(ceilx),
    FX_PROF_ENTRY(floorx),
    FX_PROF_ENTRY(addx),
    FX_PROF_ENTRY(subx),
    FX_PROF_ENTRY(mulx),
    FX_PROF_ENTRY(divx),
    FX_PROF_ENTRY(rdivx),
    FX_PROF_ENTRY(invx),
    FX_PROF_ENTRY(sqrtx),
    FX_PROF_ENTRY(isqrtx),
    FX_PROF_ENTRY(expx),
    FX_PROF_ENTRY(exp2x),
    FX_PROF_ENTRY(exp10x),
    FX_PROF_ENTRY(logx),
    FX_PROF_ENTRY(log2x),
    FX_PROF_ENTRY(log10x),
    FX_PROF_ENTRY(powx),
    FX_PROF_ENTRY(sinx),
    FX_PROF_ENTRY(cosx),
    FX_PROF_FLOAT(fadd,     "add"),
    FX_PROF_FLOAT(fsub,     "subtract"),
    FX_PROF_FLOAT(fmul,     "multiply"),
    FX_PROF_FLOAT(fdiv,     "divide"),
    FX_PROF_FLOAT(sqrt,     "sqrt"),
    FX_PROF_FLOAT(exp,      "exp"),
    FX_PROF_FLOAT(log,      "log"),
    FX_PROF_FLOAT(pow,      "pow"),
    FX_PROF_FLOAT(sin,      "sin"),
    FX_PROF_FLOAT(cos,      "cos"),
    FX_PROF_FLOAT(sprintfe, "sprintf %e"),
    FX_PROF_FLOAT(sprintfd, "sprintf %d")
};

static volatile long fx_prof_iter;
static volatile int  fx_prof_run;
static int           fx_prof_int[FX_PROF_LEN];
static uint32_t      fx_prof_u32[FX_PROF_LEN];
static fixed_t       fx_prof_x1[FX_PROF_LEN];
static fixed_t       fx_prof_x2[FX_PROF_LEN];
static fixed_t       fx_prof_x3[FX_PROF_LEN];
static float         fx_prof_f1[FX_PROF_LEN];
static float         fx_prof_f2[FX_PROF_LEN];
static float         fx_prof_f3[FX_PROF_LEN];
static double        fx_prof_d1[FX_PROF_LEN];
static double        fx_prof_d2[FX_PROF_LEN];
static double        fx_prof_d3[FX_PROF_LEN];

/*
 * -------------------------------------------------------------
 *  Exported functions
 * -------------------------------------------------------------
 */

int
main(void)
{
    struct itimerval itm = {{0, 0}, { FX_PROF_MSECS / 1000,
                                     (FX_PROF_MSECS % 1000)*1000}};
    struct sigaction act;
    unsigned         k;

    memset(&act, 0, sizeof act);
    act.sa_handler = &fx_prof_sighandler;
    sigaction(SIGPROF, &act, NULL);

    printf("function        ops/sec\n"
           "-----------------------\n");

    for (k = 0; k < sizeof fx_prof_table /
                    sizeof fx_prof_table[0]; k++)
    {
        fx_prof_iter = 0;
        fx_prof_run  = 1;

        setitimer(ITIMER_PROF, &itm, NULL);

        (*fx_prof_table[k].func)();

        printf("%-16s", fx_prof_table[k].name);
        fx_prof_print((float)fx_prof_iter*
                      (float)FX_PROF_LEN*
                      1000.0/FX_PROF_MSECS);
        printf("\n");
    }

    return 0;
}


/*
 * -------------------------------------------------------------
 *  Local functions
 * -------------------------------------------------------------
 */

static void
fx_prof_clz(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_clz(fx_prof_u32[i])));
}

static void
fx_prof_ctz(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_ctz(fx_prof_u32[i])));
}

static void
fx_prof_bitcount(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_bitcount(fx_prof_u32[i])));
}

static void
fx_prof_itox(void)
{
    FX_PROF_EXEC((fx_prof_x1[i] = fx_itox(fx_prof_int[i], 31)));
}

static void
fx_prof_ftox(void)
{
    FX_PROF_EXEC((fx_prof_x1[i] = fx_ftox(fx_prof_f1[i], 31)));
}

static void
fx_prof_dtox(void)
{
    FX_PROF_EXEC((fx_prof_x1[i] = fx_dtox(fx_prof_d1[i], 31)));
}

static void
fx_prof_xtox(void)
{
    FX_PROF_EXEC((fx_prof_x2[i] = fx_xtox(fx_prof_x1[i], 31, 24)));
}

static void
fx_prof_xtoi(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_xtoi(fx_prof_x1[i], 31)));
}

static void
fx_prof_xtof(void)
{
    FX_PROF_EXEC((fx_prof_f1[i] = fx_xtof(fx_prof_x1[i], 31)));
}

static void
fx_prof_xtod(void)
{
    FX_PROF_EXEC((fx_prof_d1[i] = fx_xtod(fx_prof_x1[i], 31)));
}

static void
fx_prof_itoa(void)
{
    char buf[32];
    memset(fx_prof_u32, 0xcc, sizeof fx_prof_u32);
    FX_PROF_EXEC(fx_itoa(buf, fx_prof_u32[i]));
}

static void
fx_prof_xtoa(void)
{
    char buf[32];
    memset(fx_prof_x1, 0xcc, sizeof fx_prof_x1);
    FX_PROF_EXEC(fx_xtoa(buf, fx_prof_x1[i], 31, 9));
}

static void
fx_prof_roundx(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_roundx(fx_prof_x1[i], 31)));
}

static void
fx_prof_ceilx(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_ceilx(fx_prof_x1[i], 31)));
}

static void
fx_prof_floorx(void)
{
    FX_PROF_EXEC((fx_prof_int[i] = fx_floorx(fx_prof_x1[i], 31)));
}

static void
fx_prof_addx(void)
{
    FX_PROF_EXEC((fx_prof_x1[i] = fx_addx(fx_prof_x1[i], fx_prof_x2[i])));
}

static void
fx_prof_subx(void)
{
    FX_PROF_EXEC((fx_prof_x1[i] = fx_subx(fx_prof_x1[i], fx_prof_x2[i])));
}

static void
fx_prof_mulx(void)
{
    FX_PROF_EXEC((fx_prof_x1[i] = fx_mulx(fx_prof_x1[i], fx_prof_x2[i], 31)));
}

static void
fx_prof_divx(void)
{
    memset(fx_prof_x2, 1, sizeof fx_prof_x2);
    FX_PROF_EXEC((fx_prof_x3[i] = fx_divx(fx_prof_x1[i], fx_prof_x2[i], 31)));
}

static void
fx_prof_rdivx(void)
{
    volatile fx_rdiv_t rdiv = {0xccccccccUL, 31};
    FX_PROF_EXEC((fx_prof_x2[i] = fx_rdivx(fx_prof_x1[i], &rdiv)));
}

static void
fx_prof_invx(void)
{
    FX_PROF_EXEC((fx_invx(0x7fffffff, 31, NULL)));
}

static void
fx_prof_sqrtx(void)
{
    FX_PROF_EXEC((fx_sqrtx(0x7fffffff, 31)));
}

static void
fx_prof_isqrtx(void)
{
    FX_PROF_EXEC((fx_isqrtx(0x7fffffff, 31, NULL)));
}

static void
fx_prof_expx(void)
{
    FX_PROF_EXEC((fx_expx(0x7fffffff, 31)));
}

static void
fx_prof_exp2x(void)
{
    FX_PROF_EXEC((fx_exp2x(0x7fffffff, 31)));
}

static void
fx_prof_exp10x(void)
{
    FX_PROF_EXEC((fx_exp10x(0x7fffffff, 31)));
}

static void
fx_prof_logx(void)
{
    FX_PROF_EXEC((fx_logx(0x7fffffff, 31)));
}

static void
fx_prof_log2x(void)
{
    FX_PROF_EXEC((fx_log2x(0x7fffffff, 31)));
}

static void
fx_prof_log10x(void)
{
    FX_PROF_EXEC((fx_log10x(0x7fffffff, 31)));
}

static void
fx_prof_powx(void)
{
    FX_PROF_EXEC((fx_powx(0x7fffffff, 31, 0x7fffffff, 31)));
}

static void
fx_prof_sinx(void)
{
    FX_PROF_EXEC((fx_sinx(0x7fffffff, 31)));
}

static void
fx_prof_cosx(void)
{
    FX_PROF_EXEC((fx_cosx(0x7fffffff, 31)));
}

static void
fx_prof_fadd(void)
{
    fx_prof_setf(fx_prof_f1, FX_PROF_LEN, 1.0f);
    fx_prof_setf(fx_prof_f2, FX_PROF_LEN, 1.0f);
    FX_PROF_EXEC((fx_prof_f3[i] = fx_prof_f1[i] + fx_prof_f2[i]));
}


static void
fx_prof_fsub(void)
{
    fx_prof_setf(fx_prof_f1, FX_PROF_LEN, 1.0f);
    fx_prof_setf(fx_prof_f2, FX_PROF_LEN, 1.0f);
    FX_PROF_EXEC((fx_prof_f3[i] = fx_prof_f1[i] - fx_prof_f2[i]));
}

static void
fx_prof_fmul(void)
{
    fx_prof_setf(fx_prof_f1, FX_PROF_LEN, 1.0f);
    fx_prof_setf(fx_prof_f2, FX_PROF_LEN, 1.0f);
    FX_PROF_EXEC((fx_prof_f3[i] = fx_prof_f1[i] * fx_prof_f2[i]));
}

static void
fx_prof_fdiv(void)
{
    fx_prof_setf(fx_prof_f1, FX_PROF_LEN, 1.0f);
    fx_prof_setf(fx_prof_f2, FX_PROF_LEN, 1.0f);
    FX_PROF_EXEC((fx_prof_f3[i] = fx_prof_f1[i] / fx_prof_f2[i]));
}

static void
fx_prof_sqrt(void)
{
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 0.1);
    FX_PROF_EXEC((fx_prof_d2[i] = sqrt(fx_prof_d1[i])));
}

static void
fx_prof_exp(void)
{
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 0.1);
    FX_PROF_EXEC((fx_prof_d2[i] = exp(fx_prof_d1[i])));
}

static void
fx_prof_log(void)
{
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 0.1);
    FX_PROF_EXEC((fx_prof_d2[i] = log(fx_prof_d1[i])));
}

static void
fx_prof_pow(void)
{
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 0.1);
    fx_prof_setd(fx_prof_d2, FX_PROF_LEN, 0.1);
    FX_PROF_EXEC((fx_prof_d3[i] = pow(fx_prof_d1[i], fx_prof_d2[i])));
}

static void
fx_prof_sin(void)
{
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 0.1);
    FX_PROF_EXEC((fx_prof_d2[i] = sin(fx_prof_d1[i])));
}

static void
fx_prof_cos(void)
{
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 0.1);
    FX_PROF_EXEC((fx_prof_d2[i] = cos(fx_prof_d1[i])));
}

static void
fx_prof_sprintfd(void)
{
    char buf[32];
    memset(fx_prof_u32, 0xcc, sizeof fx_prof_u32);
    FX_PROF_EXEC(sprintf(buf, "%d", fx_prof_u32[i]));
}

static void
fx_prof_sprintfe(void)
{
    char buf[32];
    fx_prof_setd(fx_prof_d1, FX_PROF_LEN, 1.23456789);
    FX_PROF_EXEC(sprintf(buf, "%.9e", fx_prof_d1[i]));
}

static void
fx_prof_setf(float *buf, int len, float value)
{
    int k;
    for (k = 0; k < len; k++) {
        buf[k] = value;
    }
}

static void
fx_prof_setd(double *buf, int len, double value)
{
    int k;
    for (k = 0; k < len; k++) {
        buf[k] = value;
    }
}

static void
fx_prof_sighandler(int signum)
{
    (void)signum;
    fx_prof_run = 0;
}

static void
fx_prof_print(float value)
{
    float aval = fabs(value);
    char  suff = ' ';

    if (aval < 1.0e3f) {
        suff = ' ';
    }
    else if (aval < 1.0e6f) {
        value *= 1.0e-3f;
        suff   = 'k';
    }
    else if (aval < 1.0e9f) {
        value *= 1.0e-6f;
        suff   = 'M';
    }
    else if (aval < 1.0e12f) {
        value *= 1.0e-9f;
        suff   = 'G';
    }
    else {
        value *= 1.0e-12f;
        suff   = 'T';
    }

    printf("%6.2f%c", value, suff);
}
