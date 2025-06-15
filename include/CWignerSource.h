/**
 * @defgroup WignerSource Wigner Source Class
 * @brief Class that manages source properties and computes coalescence observables.
 * @{
 */

#ifndef CWIGNERSOURCE
#define CWIGNERSOURCE

#include "TF2.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @class wignerSource
 * @brief Class to compute deuteron coalescence probability and source properties in momentum and coordinate space.
 *
 * This class manages and evaluates TF2-based Wigner functions used to estimate the probability
 * of deuteron formation (coalescence), as well as the associated energy components of the source
 * such as kinetic energy, potential energy, and the total Hamiltonian.
 */
class wignerSource
{
public:
    /**
     * @brief Constructor with optional TF2 name suffix.
     * @param name Identifier suffix for internal TF2 functions.
     */
    wignerSource(TString name = "") : mName(name) {}

    /**
     * @brief Initialize all internal TF2 functions.
     *
     * @param testMode Enable test mode (optional). If true, the integral is computed
     * using ROOT's built-in TF2::Integral method instead of the custom implementation,
     * which increases speed and efficiency but reduce numerical precision.
     */
    void initFunctions(bool testMode = false);

    /**
     * @brief Apply current parameters to all TF2 functions.
     */
    void setFunctionsParameters();

    /**
     * @brief Set the source radius independently from the parametrization.
     * @param radius Source radius.
     */
    void setRadius(double radius);

    /**
     * @brief Set the reference radius R0.
     * @param r0 Reference radius.
     */
    void setR0(double r0);

    /**
     * @brief Set relative momentum and update radius/k* accordingly.
     * @param k Relative momentum (k*).
     */
    void setRadiusK(double k);

    /**
     * @brief Set the relative momentum k* directly.
     * @param k Relative momentum (k*).
     */
    void setKstar(double k);

    /**
     * @brief Set the internal kinetic value (input k*).
     * @param k Relative momentum (k*).
     */
    void setKIn(double k);

    /**
     * @brief Set TF2 integration ranges in radius and momentum.
     * @param xmin Minimum radius.
     * @param ymin Minimum momentum.
     * @param xmax Maximum radius.
     * @param ymax Maximum momentum.
     */
    void setRanges(double xmin, double ymin, double xmax, double ymax);

    /**
     * @brief Set the reduced mass of the system.
     * @param mu Reduced mass.
     */
    void setMu(double mu);

    /**
     * @brief Set the potential spatial width.
     * @param rWidth Width in spatial dimension.
     */
    void setRWidth(double rWidth);

    /**
     * @brief Set the depth of the potential well.
     * @param v0 Potential depth.
     */
    void setV0(double v0);

    /// @brief Get the main Wigner TF2 function.
    TF2 *getWignerFunction();

    /// @brief Get the Wigner function squared.
    TF2 *getWignerFunctionForItself();

    /// @brief Get the Wigner function times the Jacobian.
    TF2 *getWignerFunctionForJacobian();

    /// @brief Get the alternative Wigner-Jacobian TF2.
    TF2 *getWignerFunction2ForJacobian();

    /// @brief Get the TF2 function for kinetic energy.
    TF2 *getKineticEnergyFunction();

    /// @brief Get the TF2 function for potential energy.
    TF2 *getPotentialEnergyFunction();

    /// @brief Get the TF2 function for the Hamiltonian.
    TF2 *getHamiltonianFunction();

    /// @brief Get Wigner-weighted kinetic energy function.
    TF2 *getWignerKinetic();

    /// @brief Get Wigner-weighted potential energy function.
    TF2 *getWignerPotential();

    /// @brief Get Wigner-weighted Hamiltonian function.
    TF2 *getWignerHamiltonan();

    /// @brief Get the Wigner function for the deuteron.
    TF2 *getWignerDeuteron();

    /// @brief Get the integral version of the deuteron Wigner function.
    TF2 *getWignerDeuteronIntegral();

    /// @brief Get the TF2 function for coalescence probability.
    TF2 *getCoalescenceProbability();

    /// @brief Get the current normalization constant.
    double getNorm();

    /// @brief Get the integral of Wigner-weighted kinetic energy.
    double getwK();

    /// @brief Get the integral of Wigner-weighted potential energy.
    double getwV();

    /// @brief Get the integral of Wigner-weighted Hamiltonian.
    double getwH();

    /// @brief Get the current value of the source radius.
    double getRadius();

    /// @brief Get the current value of k*.
    double getKStar();

    /// @brief Get the minimum radius value used.
    double getRMin();

    /// @brief Get the maximum radius value used.
    double getRMax();

    /// @brief Get the minimum momentum value used.
    double getPMin();

    /// @brief Get the maximum momentum value used.
    double getPMax();

    /// @brief Get the current reduced mass.
    double getMu();

    /// @brief Get the current radial width of the potential.
    double getRWidth();

    /// @brief Get the current depth of the potential well.
    double getV0();

    /// @brief Get the deuteron coalescence probability.
    double getcoal();

    /// @brief Get the integral over the deuteron Wigner function.
    double getDeuteronInt();

    /**
     * @brief Check normalization of the WxW function.
     * @return Integral result scaled by h^3.
     */
    double checkWxW();

    /**
     * @brief Set parameters from an external text file.
     * @param txtfile Input file name (default: "default.txt").
     */
    void SetFromTxt(const std::string &txtfile = "default.txt");

private:
    double mR0 = 1.;        ///< Reference radius R0.
    double mRadius = 1.;    ///< Source radius.
    double mKStar = 0.050;  ///< Effective relative momentum.
    double mKin = 0.050;    ///< Stores the input k* value passed to setRadiusK(), used to compute radius and effective k*.
    double mNorm = 1.;      ///< Normalization constant.
    double mMu = 0.938 / 2; ///< Reduced mass.
    double mRWidth = 3.2;   ///< Width of the potential well.
    double mV0 = -17.4E-3;  ///< Depth of the potential well.
    TString mName = "";     ///< Suffix for TF2 naming.

    double mRMin = 0;   ///< Minimum radius.
    double mRMax = 50;  ///< Maximum radius.
    double mPMin = 0;   ///< Minimum momentum.
    double mPMax = 1.5; ///< Maximum momentum.

    TF2 *mW = nullptr;            ///< Wigner function.
    TF2 *mWxW = nullptr;          ///< Wigner function squared.
    TF2 *mWxJ = nullptr;          ///< Wigner × Jacobian.
    TF2 *mWxJforItself = nullptr; ///< Alternative Wigner × Jacobian.
    TF2 *mK = nullptr;            ///< Kinetic energy.
    TF2 *mV = nullptr;            ///< Potential energy.
    TF2 *mH = nullptr;            ///< Hamiltonian.
    TF2 *mWK = nullptr;           ///< Wigner-weighted kinetic energy.
    TF2 *mWV = nullptr;           ///< Wigner-weighted potential energy.
    TF2 *mWH = nullptr;           ///< Wigner-weighted Hamiltonian.
    TF2 *mD = nullptr;            ///< Deuteron Wigner function.
    TF2 *mDInt = nullptr;         ///< Integral over deuteron Wigner.
    TF2 *mC = nullptr;            ///< Coalescence probability.

    /**
     * @brief Set 3 parameters for a TF2.
     * @param function TF2 function to configure.
     */
    void setThreeParam(TF2 *function);

    /**
     * @brief Set 4 parameters for a TF2.
     * @param function TF2 function to configure.
     */
    void setFourParam(TF2 *function);

    /**
     * @brief Set 6 parameters for a TF2.
     * @param function TF2 function to configure.
     */
    void setSixParam(TF2 *function);

    /**
     * @brief Calculate normalization for the Wigner × Jacobian function.
     */
    void normalization();

    /// @brief Update normalization in all TF2s.
    void reSetNorm();

    /// @brief Update radius in all TF2s.
    void reSetRadius();

    /// @brief Update k* in all TF2s.
    void reSetKStar();

    /// @brief Update reduced mass in all TF2s.
    void reSetMu();

    /// @brief Update potential width in all TF2s.
    void reSetRWidth();

    /// @brief Update potential depth in all TF2s.
    void reSetV0();

    /**
     * @brief Read parameters from a file.
     * @param filename File name to read from.
     * @return Vector of parameter values.
     */
    std::vector<double> readParamsFromFile(const std::string &filename);
};

#endif
/// @}