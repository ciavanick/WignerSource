#ifndef CWIGNERUTILS
#define CWIGNERUTILS

#include "TF2.h"
#include "TFile.h"
#include "TH2.h"

/**
 * Static utility class for Wigner function calculations.
 * Provides TF2-compatible function implementations.
 */

class wignerUtils {
public:
    /**
     * Compute the standard Wigner source function.
     * Uses position x and parameter array pm.
     */
    static double wignerSource(double *x, double *pm);

    /**
     * Compute the Wigner source function multiplied by jacobianW2.
     */
    static double wignerSource2(double *x, double *pm);

    /**
     * Compute the Wigner source function weighted by a Jacobian.
     */
    static double jacobianFun(double *x, double *pm);

    /**
     * Compute an alternative Jacobian-weighted Wigner function.
     */
    static double jacobianW2(double *x, double *pm);

    /**
     * Compute classical kinetic energy.
     * Depends on momentum p and reduced mass pm[3].
     */
    static double kineticEnergy(double *x, double *pm);

    /**
     * Return a square-well potential value based on radius and threshold.
     */
    static double potentialEnergy(double *x, double *pm);

    /**
     * Return the Hamiltonian as the sum of kinetic and potential energy.
     */
    static double hamiltonian(double *x, double *pm);

    /**
     * Return the Wigner function for kinetic energy.
     */
    static double wK(double *x, double *pm);

    /**
     * Return the Wigner function for potential energy.
     */
    static double wV(double *x, double *pm);

    /**
     * Return the Wigner function for the Hamiltonian.
     */
    static double wH(double *x, double *pm);

    /**
     * Return the value of the deuteron Wigner function from histogram data.
     */
    static double wignerDeuteron(double *x, double *pm);

    /**
     * Return deuteron Wigner function multiplied by Jacobian.
     */
    static double wignerDeuteronIntegral(double *x, double *pm);

    /**
     * Calculate the probability for coalescence using deuteron and source.
     */
    static double coalescenceProbability(double *x, double *pm);

    /**
     * Compute the effective radius based on input momentum and initial radius.
     */
    static double radius(double k, double r0);

    /**
     * Compute the effective k* value based on total and wave momenta.
     */
    static double kStarEff(double k, double radius);

    /**
     * Compute a numeric integral over a TF2 within specified (x, p) bounds.
     * If no bounds are provided, use the internal default limits.
     */
    static double integral(TF2 *function, double minX = mMinX, double maxX = mMaxX, double minP = mMinP, double maxP = mMaxP);

    /**
     * Get the minimum value of X (radius).
     */
    static double getMinX();

    /**
     * Get the maximum value of X (radius).
     */
    static double getMaxX();

    /**
     * Get the minimum value of P (momentum).
     */
    static double getMinP();

    /**
     * Get the maximum value of P (momentum).
     */
    static double getMaxP();

    /**
     * Get the h-bar * c conversion factor in GeV*fm.
     */
    static double getHCut();

    /**
     * Set the minimum radius for numerical integration.
     */
    static void setMinX(double minX);

    /**
     * Set the maximum radius for numerical integration.
     */
    static void setMaxX(double maxX);

    /**
     * Set the minimum momentum for numerical integration.
     */
    static void setMinP(double minP);

    /**
     * Set the maximum momentum for numerical integration.
     */
    static void setMaxP(double maxP);

    /**
     * Set all (x, p) integration limits in a single call.
     */
    static void setIntegrationRanges(double minX, double maxX, double minP, double maxP);

private:
    // Internal constants for numerical integration and physics constants
    static double mHCut;
    static double mMinX;
    static double mMaxX;
    static double mMinP;
    static double mMaxP;
    static double mDx;
    static double mDp;
    static double mFactor;

    // Deuteron data histogram and ROOT file
    static TH2D* mH;
    static TFile *mFileDeuteron;
};

#endif
