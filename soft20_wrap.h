/// This struct is to preserve memory for later calculation.
struct workspae_impl;
typedef struct workspace
{
    struct workspace_impl *impl;
} workspace;

/**
 * @brief This function reserves memory inside the struct `workspace.
 *
 * @param bwIn Bandwidth of input signal.
 * @param bwOut Parameter to define the resolution of SO(3) rotation.
 */
void reserve(workspace * const w, const int bwIn, const int bwOut);

/// Release memory
void clear(workspace *w);

/**
 * @brief This function is to calculate correlation for each degree of SO(3) rotation
 *
 * @param bwIn Maximum bandwidth l of spherical hamononic function of x0, x1, etc.
 * @param bwOut Parameter to define the resolution of SO(3) rotation.
 * @param x0 Target coefficients of spherical harmonic function basis. Array Length of x0 should be bwIn^2.
 * @param x1 Source coefficients of spherical harmonic function basis. Array Length of x1 should be bwIn^2.
 * @param alpha Z-axis-wise rotation.
 * @param beta Y-axis-wise rotation.
 * @param dist L2-norm between shape0 and shape1
 * @param gamma Z-axis-wise rotation.
 */
void calc_rot_dist_c(const int bwIn, const int bwOut,
        const double * const x0, const double * const x1,
        const double * const y0, const double * const y1,
        const double * const z0, const double * const z1,
        double * const alpha, double * const beta, double * const gamma,
        double * const dist,
        workspace * const w);
