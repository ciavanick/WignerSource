#include "CWignerUtils.h"
#include "TMath.h"
#include "TF2.h"

double wignerUtils::mHCut = 0.1973; //GeV fm
double wignerUtils::mMinX = 0.;
double wignerUtils::mMaxX = 100;
double wignerUtils::mMinP = 0.;
double wignerUtils::mMaxP = 3.000;
double wignerUtils::mDx = 0.01;
double wignerUtils::mDp = 0.001;
//_________________________________________________________________________
double wignerUtils::wignerSource(double *x, double *pm){
    double r = x[0];
    double p = x[1];

    return pm[0] * TMath::Exp(-r * r * 0.5 / (pm[1] * pm[1]) - ((p - pm[2]) * (p - pm[2])) * 0.5 * (pm[1] * pm[1]) / (mHCut * mHCut)) / ((2 * TMath::Pi() * mHCut)*(2 * TMath::Pi() * mHCut)*(2 * TMath::Pi() * mHCut));
}
//_________________________________________________________________________
double wignerUtils::wignerSource2(double *x, double *pm){
    return wignerSource(x, pm) * jacobianW2(x, pm);
}
//_________________________________________________________________________
double wignerUtils::jacobianFun(double *x, double *pm){
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p);
    if(pm[2] > 1e-3) jacobian *= 4 * TMath::Pi() * 2 * TMath::Pi() * (1 - TMath::Exp(-2 * pm[1] * pm[1] / (mHCut * mHCut) * p * pm[2])) / p / pm[2] / (pm[1] * pm[1]) * mHCut * mHCut;
    else             jacobian *= 16 * TMath::Pi() * TMath::Pi();
    return jacobian * wignerSource(x, pm);
}
//_________________________________________________________________________
double wignerUtils::jacobianW2(double *x, double *pm){
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p);
    if(pm[2] > 1e-3) jacobian *= 4 * TMath::Pi() * 2 * TMath::Pi() * (1 - TMath::Exp(-4 * pm[1] * pm[1] / (mHCut * mHCut) * p * pm[2])) / p / pm[2] / (pm[1] * pm[1]) * mHCut * mHCut * 0.5;
    else             jacobian *= 16 * TMath::Pi() * TMath::Pi();
    return jacobian * wignerSource(x, pm);
}
//_________________________________________________________________________
double wignerUtils::kineticEnergy(double *x, double *pm){
    double r = x[0];
    double p = x[1];

    return (p * p) / (2 * pm[3]);
}
//_________________________________________________________________________
double wignerUtils::potentialEnergy(double *x, double *pm){
    double r = x[0];
    double p = x[1];

    if (r < pm[4])
        return pm[5];
    else
        return 0;
}
//_________________________________________________________________________
double wignerUtils::hamiltonian(double *x, double *pm){
    return kineticEnergy(x, pm) + potentialEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wK(double *x, double *pm){
    return jacobianFun(x, pm) * kineticEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wV(double *x, double *pm){
    return jacobianFun(x, pm) * potentialEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wH(double *x, double *pm){
    return jacobianFun(x, pm) * hamiltonian(x, pm);
}
//_________________________________________________________________________
double wignerUtils::rSource(double k){
    return TMath::Sqrt(3*0.5)*mHCut/k;
}
//_________________________________________________________________________
double wignerUtils::getMinX(){
    return mMinX;
}
//_________________________________________________________________________
double wignerUtils::getMaxX(){
    return mMaxX;
}
//_________________________________________________________________________
double wignerUtils::getMinP(){
    return mMinP;
}
//_________________________________________________________________________
double wignerUtils::getHCut(){
    return mHCut;
}
//_________________________________________________________________________
double wignerUtils::getMaxP(){
    return mMaxP;
}
//_________________________________________________________________________
void wignerUtils::setMinX(double minX){
    mMinX = minX;
}
//_________________________________________________________________________
void wignerUtils::setMaxX(double maxX){
    mMaxX = maxX;
}
//_________________________________________________________________________
void wignerUtils::setMinP(double minP){
    mMinP = minP;
}
//_________________________________________________________________________
void wignerUtils::setMaxP(double maxP){
    mMaxP = maxP;
}
//_________________________________________________________________________
void wignerUtils::setIntegrationRanges(double minX, double maxX, double minP, double maxP){
    setMinX(minX);
    setMaxX(maxX);
    setMinP(minP);
    setMaxP(maxP);
}
//_________________________________________________________________________
double wignerUtils::integral(TF2 *function){
    double res = 0;
    for (float x = mDx / 2 + mMinX; x < mMaxX; x += mDx)
    {
        for (float p = mDp / 2 + mMinP; p < mMaxP; p += mDp)
        {
            res += function->Eval(x, p);
        }
    }
    res *= mDx * mDp;
    
    return res;
}