#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "so3_correlate_fftw.h"
#include "soft_fftw.h"
#include "soft20_wrap.h"

/// This struct is to preserve memory for later calculation.
typedef struct workspace_impl
{
    fftw_complex *workspace1;
    fftw_complex *workspace2;
    double *workspace3;
    double *sigCoefR;
    double *sigCoefI;
    double *patCoefR;
    double *patCoefI;
    fftw_complex *so3Coef;
    fftw_complex *so3Sig_x;
    fftw_complex *so3Sig_y;
    fftw_complex *so3Sig_z;
    fftw_plan p1_x;
    fftw_plan p1_y;
    fftw_plan p1_z;
} workspace_impl;

/**
 * @brief This function reserves memory inside the struct `workspace.
 *
 * @param bwIn Bandwidth of input signal.
 * @param bwOut Parameter to define the resolution of SO(3) rotation.
 */
void reserve(workspace * const w, const int bwIn, const int bwOut)
{
    // workspace_impl
    w->impl = (workspace_impl*) malloc( sizeof(workspace_impl) );

    // members of workspace_impl
    const int n = 2 * bwIn;
    w->impl->workspace1 = fftw_malloc( sizeof(fftw_complex) * (8 * pow(bwOut, 3)) );
    w->impl->workspace2 = fftw_malloc( sizeof(fftw_complex) * ((14 * bwIn*bwIn) + (48 * bwIn)) );
    w->impl->workspace3 = (double *) malloc( sizeof(double) * (12*n + n*bwIn) );
    w->impl->sigCoefR = (double*) malloc( sizeof(double) * bwIn * bwIn );
    w->impl->sigCoefI = (double*) malloc( sizeof(double) * bwIn * bwIn );
    w->impl->patCoefR = (double*) malloc( sizeof(double) * bwIn * bwIn );
    w->impl->patCoefI = (double*) malloc( sizeof(double) * bwIn * bwIn );
    w->impl->so3Coef = fftw_malloc( sizeof(fftw_complex) * ((4*pow(bwOut,3) - bwOut)/3) );

    const int rank = 2;
    const int howmany = n * n;
    const int na[2] = {1, 2 * bwOut};
    const int inembed[2] = {n, n * n};
    const int onembed[2] = {n, n * n};
    const int istride = 1;
    const int ostride = 1;
    const int idist = n;
    const int odist = n;
    w->impl->so3Sig_x = fftw_malloc( sizeof(fftw_complex) * (8*bwOut*bwOut*bwOut) );
    w->impl->so3Sig_y = fftw_malloc( sizeof(fftw_complex) * (8*bwOut*bwOut*bwOut) );
    w->impl->so3Sig_z = fftw_malloc( sizeof(fftw_complex) * (8*bwOut*bwOut*bwOut) );
    w->impl->p1_x = fftw_plan_many_dft(rank, na, howmany, w->impl->workspace1, inembed, istride, idist, w->impl->so3Sig_x, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    w->impl->p1_y = fftw_plan_many_dft(rank, na, howmany, w->impl->workspace1, inembed, istride, idist, w->impl->so3Sig_y, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    w->impl->p1_z = fftw_plan_many_dft(rank, na, howmany, w->impl->workspace1, inembed, istride, idist, w->impl->so3Sig_z, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);

}

/// Release memory
void clear(workspace *w)
{
    fftw_free(w->impl->workspace1);
    fftw_free(w->impl->workspace2);
    free(w->impl->workspace3);
    free(w->impl->sigCoefR);
    free(w->impl->sigCoefI);
    free(w->impl->patCoefR);
    free(w->impl->patCoefI);
    fftw_free(w->impl->so3Coef);
    fftw_free(w->impl->so3Sig_x);
    fftw_free(w->impl->so3Sig_y);
    fftw_free(w->impl->so3Sig_z);
    fftw_destroy_plan(w->impl->p1_x);
    fftw_destroy_plan(w->impl->p1_y);
    fftw_destroy_plan(w->impl->p1_z);
    free(w->impl);
}

/**
 * @brief This function is to calculate correlation for each degree of SO(3) rotation
 *
 * @param bwIn Maximum bandwidth l of spherical hamononic function of x0, x1, etc.
 * @param bwOut Parameter to define the resolution of SO(3) rotation.
 * @param x0 Target coefficients of spherical harmonic function basis. Array Length of x0 should be bwIn^2.
 * @param x1 Source coefficients of spherical harmonic function basis. Array Length of x1 should be bwIn^2.
 * @param alpha Z-axis-wise rotation.
 * @param beta Y-axis-wise rotation.
 * @param gamma Z-axis-wise rotation.
 */
void calc_rot_dist_c(const int bwIn, const int bwOut,
        const double * const x0, const double * const x1,
        const double * const y0, const double * const y1,
        const double * const z0, const double * const z1,
        double * const alpha, double * const beta, double * const gamma,
        double * const dist,
        workspace * const w)
{
    const int degLim = bwOut;

    {
        w->impl->sigCoefR[0] = 0;
        w->impl->sigCoefI[0] = 0;
        w->impl->patCoefR[0] = 0;
        w->impl->patCoefI[0] = 0;
    }
    for (int i = 1; i < bwIn; i++)
    {
        w->impl->sigCoefR[i] = x0[i];
        w->impl->sigCoefI[i] = 0;
        w->impl->patCoefR[i] = x1[i];
        w->impl->patCoefI[i] = 0;
    }
    so3CombineCoef_fftw(bwIn, bwOut, degLim,
        w->impl->sigCoefR, w->impl->sigCoefI, w->impl->patCoefR, w->impl->patCoefI, w->impl->so3Coef);
    Inverse_SO3_Naive_fftw(bwOut, w->impl->so3Coef, w->impl->so3Sig_x,
            w->impl->workspace1, w->impl->workspace2, w->impl->workspace3, &w->impl->p1_x, 0);

    for (int i = 1; i < bwIn; i++)
    {
        w->impl->sigCoefR[i] = y0[i];
        // w->impl->sigCoefI[i] = 0;
        w->impl->patCoefR[i] = y1[i];
        // w->impl->patCoefI[i] = 0;
    }
    so3CombineCoef_fftw(bwIn, bwOut, degLim,
        w->impl->sigCoefR, w->impl->sigCoefI, w->impl->patCoefR, w->impl->patCoefI, w->impl->so3Coef);
    Inverse_SO3_Naive_fftw(bwOut, w->impl->so3Coef, w->impl->so3Sig_y,
            w->impl->workspace1, w->impl->workspace2, w->impl->workspace3, &w->impl->p1_y, 0);

    for (int i = 1; i < bwIn; i++)
    {
        w->impl->sigCoefR[i] = z0[i];
        // w->impl->sigCoefI[i] = 0;
        w->impl->patCoefR[i] = z1[i];
        // w->impl->patCoefI[i] = 0;
    }
    so3CombineCoef_fftw(bwIn, bwOut, degLim,
        w->impl->sigCoefR, w->impl->sigCoefI, w->impl->patCoefR, w->impl->patCoefI, w->impl->so3Coef);
    Inverse_SO3_Naive_fftw(bwOut, w->impl->so3Coef, w->impl->so3Sig_z,
            w->impl->workspace1, w->impl->workspace2, w->impl->workspace3, &w->impl->p1_z, 0);

    /* now find max value */
    double maxval = 0.0 ;
    double maxloc = 0 ;
    for (int i = 0 ; i < 8 * pow(bwOut, 3) ; i++)
    {
        /*
        if (so3Sig[i][0] >= maxval)
        {
            maxval = so3Sig[i][0];
            maxloc = i ;
        }
        */
#define NORM( x ) ( (x[0])*(x[0]) + (x[1])*(x[1]) )
        double tmpval = NORM( w->impl->so3Sig_x[i] ) + NORM( w->impl->so3Sig_y[i] ) + NORM( w->impl->so3Sig_z[i] );
        if ( tmpval > maxval )
        {
            maxval = tmpval;
            maxloc = i ;
        }
    }

    int tmp;
    const int ii = floor( maxloc / (4.*bwOut*bwOut) );
    tmp = maxloc - (ii*4.*bwOut*bwOut);
    const int jj = floor( tmp / (2.*bwOut) );
    tmp = maxloc - (ii *4*bwOut*bwOut) - jj*(2*bwOut);
    const int kk = tmp ;

    *alpha = M_PI * jj / ((double) bwOut);
    *beta = M_PI * (2*ii + 1) / (4.*bwOut);
    *gamma = M_PI * kk / ((double) bwOut);
    *dist = maxval;
}
