#ifndef CWIGNERUTILS
#define CWIGNERUTILS

#include "TF2.h"
#include "TFile.h"
#include "TH2.h"

class wignerUtils{
    public:
        static double wignerSource(double *x, double *pm);
        static double wignerSource2(double *x, double *pm);
        static double jacobianFun(double *x, double *pm);
        static double jacobianW2(double *x, double *pm);
        static double kineticEnergy(double *x, double *pm);
        static double potentialEnergy(double *x, double *pm);
        static double hamiltonian(double *x, double *pm);
        static double wK(double *x, double *pm);
        static double wH(double *x, double *pm);
        static double wV(double *x, double *pm);
        static double wignerDeuteron(double *x, double *pm);
        static double wignerDeuteronIntegral(double *x, double *pm);
        static double coalescenceProbability(double *x, double *pm);
        static double radius(double k, double r0);
        static double kStarEff(double k, double radius);

        static double integral(TF2 *function, double minX = mMinX, double maxX = mMaxX, double minP = mMinP, double maxP = mMaxP);

        static double getMinX();
        static double getMaxX();
        static double getMinP();
        static double getMaxP();
        static double getHCut();

        static void setMinX(double minX);
        static void setMaxX(double maxX);
        static void setMinP(double minP);
        static void setMaxP(double maxP);
        static void setIntegrationRanges(double minX, double maxX, double minP, double maxP);
    private:
        static double mHCut;
        static double mMinX;
        static double mMaxX;
        static double mMinP;
        static double mMaxP;
        static double mDx;
        static double mDp;
        static double mFactor;

        static TH2D* mH;
        static TFile *mFileDeuteron;
};

#endif