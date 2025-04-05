#ifndef CWIGNERSOURCE
#define CWIGNERSOURCE

#include "TF2.h"

class wignerSource{
    public:
        wignerSource(TString name = "") : mName(name) {}
        void initFunctions();
        void setFunctionsParameters();
        void setRadius(double r0);
        void setRadiusK(double r0);
        void setKstar(double k);
        void setRanges(double xmin, double ymin, double xmax, double ymax);
        void setMu(double mu);
        void setRWidth(double rWidth);
        void setV0(double v0);

        TF2* getWignerFunction();
        TF2* getWignerFunctionForItself();
        TF2* getWignerFunctionForJacobian();
        TF2* getWignerFunction2ForJacobian();
        TF2* getKineticEnergyFunction();
        TF2* getPotentialEnergyFunction();
        TF2* getHamiltonianFunction();
        TF2* getWignerKinetic();
        TF2* getWignerPotential();
        TF2* getWignerHamiltonan();

        double getNorm();
        double getwK();
        double getwV();
        double getwH();
        double getRadius();
        double getKStar();
        double getRMin();
        double getRMax();
        double getPMin();
        double getPMax();
        double getMu();
        double getRWidth();
        double getV0();
        
        double checkWxW();

    private:
        double mR0 = 5.;
        double mKStar = 0.050;
        double mNorm = 1.;
        double mMu = 0.938/2;
        double mRWidth = 2.1;
        double mV0 = -0.0337;
        TString mName = "";

        double mRMin = 0;
        double mRMax = 50;
        double mPMin = 0;
        double mPMax = 0.300;

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

        void setThreeParam(TF2 *function);
        void setFourParam(TF2 *function);
        void setSixParam(TF2 *function);

        void normalization();

        void reSetNorm();
        void reSetRadius();
        void reSetKStar();
        void reSetMu();
        void reSetRWidth();
        void reSetV0();

};


#endif