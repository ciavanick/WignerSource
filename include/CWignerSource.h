#ifndef CWIGNERSOURCE
#define CWIGNERSOURCE

#include "TF2.h"

class wignerSource {
public:
    /**
     * Constructor with optional TF2 name suffix.
     * Used to uniquely identify all generated TF2 objects.
     */
    wignerSource(TString name = "") : mName(name) {}

    /**
     * Initialize all TF2-based Wigner-related functions.
     * Must be called before using get*Function methods.
     */
    void initFunctions();

    /**
     * Apply current object parameters to all initialized functions.
     * Used after manually setting parameters like radius or mass.
     */
    void setFunctionsParameters();

    /**
     * Set the source radius and update dependent parameters.
     */
    void setRadius(double radius);

    /**
     * Set the reference radius R0 used in certain calculations.
     */
    void setR0(double r0);

    /**
     * Set kinetic energy and compute derived quantities like radius and k*.
     * r0 is the kinetic energy input.
     */
    void setRadiusK(double r0);

    /**
     * Set the relative momentum k* directly.
     */
    void setKstar(double k);

    /**
     * Set kinetic input value used internally.
     */
    void setKIn(double k);

    /**
     * Set the spatial and momentum boundaries for TF2 functions.
     * xmin/xmax define radius range, ymin/ymax define momentum range.
     */
    void setRanges(double xmin, double ymin, double xmax, double ymax);

    /**
     * Set the reduced mass mu for the system.
     */
    void setMu(double mu);

    /**
     * Set the width of the potential function in spatial dimensions.
     */
    void setRWidth(double rWidth);

    /**
     * Set the depth of the potential well.
     */
    void setV0(double v0);

    /// Get the main Wigner TF2 function.
    TF2* getWignerFunction();

    /// Get the Wigner function calculated with itself.
    TF2* getWignerFunctionForItself();

    /// Get the Wigner function multiplied by the Jacobian.
    TF2* getWignerFunctionForJacobian();

    /// Get an alternative Wigner-Jacobian function.
    TF2* getWignerFunction2ForJacobian();

    /// Get the TF2 function representing kinetic energy.
    TF2* getKineticEnergyFunction();

    /// Get the TF2 function representing potential energy.
    TF2* getPotentialEnergyFunction();

    /// Get the TF2 function representing the full Hamiltonian.
    TF2* getHamiltonianFunction();

    /// Get the Wigner-weighted kinetic energy function.
    TF2* getWignerKinetic();

    /// Get the Wigner-weighted potential energy function.
    TF2* getWignerPotential();

    /// Get the Wigner-weighted Hamiltonian function.
    TF2* getWignerHamiltonan();

    /// Get the Wigner function used for the deuteron.
    TF2* getWignerDeuteron();

    /// Get the integral version of the deuteron Wigner function.
    TF2* getWignerDeuteronIntegral();

    /// Get the function used to compute coalescence probability.
    TF2* getCoalescenceProbability();

    /// Get the current normalization constant used in TF2s.
    double getNorm();

    /// Get the integral of the Wigner-weighted kinetic energy.
    double getwK();

    /// Get the integral of the Wigner-weighted potential energy.
    double getwV();

    /// Get the integral of the Wigner-weighted Hamiltonian.
    double getwH();

    /// Get the current value of the source radius.
    double getRadius();

    /// Get the current value of k*.
    double getKStar();

    /// Get the minimum radius value used for TF2 domains.
    double getRMin();

    /// Get the maximum radius value used for TF2 domains.
    double getRMax();

    /// Get the minimum momentum value used for TF2 domains.
    double getPMin();

    /// Get the maximum momentum value used for TF2 domains.
    double getPMax();

    /// Get the current value of reduced mass.
    double getMu();

    /// Get the current radial width of the potential.
    double getRWidth();

    /// Get the current potential well depth.
    double getV0();

    /// Get the integral of the coalescence probability function.
    double getcoal();

    /// Get the integral over the deuteron Wigner function.
    double getDeuteronInt();

    /**
     * Verify the normalization of the WxW function.
     * Returns the scaled result of its integral.
     */
    double checkWxW();

private:
    double mR0 = 1.;
    double mRadius = 1.;
    double mKStar = 0.050;
    double mKin = 0.050;
    double mNorm = 1.;
    double mMu = 0.938 / 2;
    double mRWidth = 3.2; // 2.1
    double mV0 = -17.4E-3; // -0.0337
    TString mName = "";

    double mRMin = 0;
    double mRMax = 50;
    double mPMin = 0;
    double mPMax = 1.5;

    TF2 *mW = nullptr;
    TF2 *mWxW = nullptr;
    TF2 *mWxJ = nullptr;
    TF2 *mWxJforItself = nullptr;
    TF2 *mK = nullptr;
    TF2 *mV = nullptr;
    TF2 *mH = nullptr;
    TF2 *mWK = nullptr;
    TF2 *mWV = nullptr;
    TF2 *mWH = nullptr;
    TF2 *mD = nullptr;
    TF2 *mDInt = nullptr;
    TF2 *mC = nullptr;

    /**
     * Assign 3 parameters to a TF2 function.
     */
    void setThreeParam(TF2 *function);

    /**
     * Assign 4 parameters to a TF2 function.
     */
    void setFourParam(TF2 *function);

    /**
     * Assign 6 parameters to a TF2 function.
     */
    void setSixParam(TF2 *function);

    /**
     * Compute normalization constant for the Jacobian-weighted Wigner function.
     */
    void normalization();

    void reSetNorm();
    void reSetRadius();
    void reSetKStar();
    void reSetMu();
    void reSetRWidth();
    void reSetV0();
};

#endif
