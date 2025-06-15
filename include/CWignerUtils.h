/**
 * @defgroup WignerUtils Static Wigner Utility Functions
 * @brief Numerical tools and energy computations for Wigner-based coalescence.
 * @{
 */

#ifndef CWIGNERUTILS
#define CWIGNERUTILS

#include "TF2.h"
#include "TFile.h"
#include "TH2.h"

/**
 * @class wignerUtils
 * @brief Static utility class for Wigner function and coalescence probability calculations.
 *
 * This class provides TF2-compatible static functions for computing Wigner distributions,
 * energy components, and coalescence observables used in two-particle correlation studies.
 * It also includes tools for numerical integration and access to deuteron wavefunction data.
 */
class wignerUtils
{
public:
    /**
     * @brief Compute the standard Wigner source function.
     * @param x Coordinate array: x[0] = radius, x[1] = momentum.
     * @param pm Parameter array: [norm, radius, k*].
     * @return Value of the Wigner function.
     */
    static double wignerSource(double *x, double *pm);

    /**
     * @brief Compute the Wigner source function multiplied by jacobianW2.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Value of Wigner function × Jacobian.
     */
    static double wignerSource2(double *x, double *pm);

    /**
     * @brief Compute the Wigner function weighted by a Jacobian.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Wigner × Jacobian value.
     */
    static double jacobianFun(double *x, double *pm);

    /**
     * @brief Alternative Jacobian-weighted Wigner function.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Wigner × Jacobian (alternative form).
     */
    static double jacobianW2(double *x, double *pm);

    /**
     * @brief Compute the classical kinetic energy.
     * @param x Coordinate array.
     * @param pm Parameter array: pm[3] = reduced mass.
     * @return Kinetic energy value.
     */
    static double kineticEnergy(double *x, double *pm);

    /**
     * @brief Compute square-well potential energy.
     * @param x Coordinate array.
     * @param pm Parameter array: pm[4] = width, pm[5] = depth.
     * @return Potential energy value.
     */
    static double potentialEnergy(double *x, double *pm);

    /**
     * @brief Compute total energy (Hamiltonian).
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Hamiltonian = kinetic + potential.
     */
    static double hamiltonian(double *x, double *pm);

    /**
     * @brief Compute Wigner-weighted kinetic energy.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Weighted kinetic energy.
     */
    static double wK(double *x, double *pm);

    /**
     * @brief Compute Wigner-weighted potential energy.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Weighted potential energy.
     */
    static double wV(double *x, double *pm);

    /**
     * @brief Compute Wigner-weighted total energy.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Weighted Hamiltonian.
     */
    static double wH(double *x, double *pm);

    /**
     * @brief Evaluate deuteron Wigner function from histogram derived from numerical integration, the file is in deuteronFunction/.
     * @param x Coordinate array.
     * @param pm Unused.
     * @return Interpolated value from deuteron histogram.
     */
    static double wignerDeuteron(double *x, double *pm);

    /**
     * @brief Compute Jacobian-weighted deuteron Wigner function.
     * @param x Coordinate array.
     * @param pm Unused.
     * @return Histogram value × Jacobian.
     */
    static double wignerDeuteronIntegral(double *x, double *pm);

    /**
     * @brief Compute deuteron coalescence probability.
     * @param x Coordinate array.
     * @param pm Parameter array.
     * @return Product of source and deuteron Wigner × Jacobian.
     */
    static double coalescenceProbability(double *x, double *pm);

    /**
     * @brief Compute effective source radius based on k* and R₀.
     * @param k Relative momentum (k*).
     * @param r0 Reference radius.
     * @return Effective radius.
     */
    static double radius(double k, double r0);

    /**
     * @brief Compute effective k* based on total k* and radius.
     * @param k Input relative momentum.
     * @param radius Source radius.
     * @return Effective k* value.
     */
    static double kStarEff(double k, double radius);

    /**
     * @brief Numerically integrate a TF2 over specified (x, p) range.
     * @param function Pointer to TF2 object.
     * @param minX Lower x (radius) limit.
     * @param maxX Upper x (radius) limit.
     * @param minP Lower p (momentum) limit.
     * @param maxP Upper p (momentum) limit.
     * @return Integral result.
     */
    static double integral(TF2 *function, double minX = mMinX, double maxX = mMaxX, double minP = mMinP, double maxP = mMaxP);

    /// @brief Get minimum radius used for integration.
    static double getMinX();

    /// @brief Get maximum radius used for integration.
    static double getMaxX();

    /// @brief Get minimum momentum used for integration.
    static double getMinP();

    /// @brief Get maximum momentum used for integration.
    static double getMaxP();

    /// @brief Get the h-bar * c conversion constant in GeV·fm.
    static double getHCut();

    /// @brief Set minimum radius for integration.
    static void setMinX(double minX);

    /// @brief Set maximum radius for integration.
    static void setMaxX(double maxX);

    /// @brief Set minimum momentum for integration.
    static void setMinP(double minP);

    /// @brief Set maximum momentum for integration.
    static void setMaxP(double maxP);

    /**
     * @brief Set all integration bounds in a single call.
     * @param minX Minimum radius.
     * @param maxX Maximum radius.
     * @param minP Minimum momentum.
     * @param maxP Maximum momentum.
     */
    static void setIntegrationRanges(double minX, double maxX, double minP, double maxP);

    /// @brief If true, uses TF2::Integral instead of manual integration, good for testing.
    static bool testMode;

private:
    // Constants for integration and physical conversion
    static double mHCut;   ///< ℏ·c conversion factor [GeV·fm]
    static double mMinX;   ///< Minimum radius for integration.
    static double mMaxX;   ///< Maximum radius for integration.
    static double mMinP;   ///< Minimum momentum for integration.
    static double mMaxP;   ///< Maximum momentum for integration.
    static double mDx;     ///< dx step for manual integration.
    static double mDp;     ///< dp step for manual integration.
    static double mFactor; ///< Conversion factor used in radius/k* calculations.

    // ROOT objects for deuteron Wigner function
    static TH2D *mH;             ///< 2D histogram with deuteron Wigner data.
    static TFile *mFileDeuteron; ///< ROOT file holding the histogram.
};

#endif
/// @}