#include "CWignerUtils.h"
#include "TMath.h"
#include "TF2.h"

double wignerUtils::mHCut = 0.1973; // GeV fm
double wignerUtils::mMinX = 0.;
double wignerUtils::mMaxX = 20.;
double wignerUtils::mMinP = 0.;
double wignerUtils::mMaxP = 0.6;
double wignerUtils::mDx = 0.01;
double wignerUtils::mDp = 0.001;
double wignerUtils::mFactor = sqrt(3. / 8) * mHCut;

TFile *wignerUtils::mFileDeuteron = new TFile("deuteronFunction/wigner2.root", "READ");
TH2D *wignerUtils::mH = (TH2D *)mFileDeuteron->Get("h");

bool wignerUtils::testMode = false;
//_________________________________________________________________________
double wignerUtils::wignerSource(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double norm = 1. / (TMath::Pi() * mHCut);
    norm *= norm * norm;

    return pm[0] * norm * TMath::Exp(-r * r * 0.25 / (pm[1] * pm[1]) - 4 * (p * p + pm[2] * pm[2] - 2 * p * pm[2]) * (pm[1] * pm[1]) / (mHCut * mHCut));
}
//_________________________________________________________________________
double wignerUtils::wignerSource2(double *x, double *pm)
{
    return wignerSource(x, pm) * jacobianW2(x, pm);
}
//_________________________________________________________________________
double wignerUtils::jacobianFun(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p) * 16 * TMath::Pi() * TMath::Pi();

    float kstarP = pm[2] * p;
    if (kstarP < 1E-16)
    {
        kstarP = 1E-16;
    }
    double alpha = 8 * kstarP * pm[1] * pm[1] / (mHCut * mHCut);

    jacobian *= 0.5 * (1 - TMath::Exp(-2 * alpha)) / alpha;

    return wignerSource(x, pm) * jacobian;
}
//_________________________________________________________________________
double wignerUtils::jacobianW2(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p) * 16 * TMath::Pi() * TMath::Pi();

    float kstarP = pm[2] * p;
    if (kstarP < 1E-16)
    {
        kstarP = 1E-16;
    }
    double alpha = 16 * kstarP * pm[1] * pm[1] / (mHCut * mHCut);
    jacobian *= 0.5 * (1 - TMath::Exp(-2 * alpha)) / alpha;
    return jacobian * wignerSource(x, pm);
}
//_________________________________________________________________________
double wignerUtils::kineticEnergy(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return (p * p) / (2 * pm[3]);
}
//_________________________________________________________________________
double wignerUtils::potentialEnergy(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double E = 0;
    if (r < pm[4])
    {
        E = pm[5];
    }
    return E;
}
//_________________________________________________________________________
double wignerUtils::hamiltonian(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return kineticEnergy(x, pm) + potentialEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wK(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return jacobianFun(x, pm) * kineticEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wV(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return jacobianFun(x, pm) * potentialEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wH(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return jacobianFun(x, pm) * hamiltonian(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wignerDeuteron(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return mH->Interpolate(r, p);
}
//_________________________________________________________________________
double wignerUtils::wignerDeuteronIntegral(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = 16 * TMath::Pi() * TMath::Pi() * r * r * p * p;
    return mH->Interpolate(r, p) * jacobian;
}
//_________________________________________________________________________
double wignerUtils::coalescenceProbability(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p) * 16 * TMath::Pi() * TMath::Pi();

    float kstarP = pm[2] * p;
    if (kstarP < 1E-16)
    {
        kstarP = 1E-16;
    }
    double alpha = 8 * kstarP * pm[1] * pm[1] / (mHCut * mHCut);

    jacobian *= 0.5 * (1 - TMath::Exp(-2 * alpha)) / alpha;

    return wignerDeuteron(x, pm) * wignerSource(x, pm) * jacobian;
}
//_________________________________________________________________________
double wignerUtils::radius(double k, double r0)
{
    double radiusWave = mFactor / k;
    return sqrt(radiusWave * radiusWave + r0 * r0);
}
//_________________________________________________________________________
double wignerUtils::kStarEff(double k, double radius)
{
    double kstarWave = mFactor / radius;
    return radius = sqrt(k * k - kstarWave * kstarWave);
}
//_________________________________________________________________________
double wignerUtils::getMinX()
{
    return mMinX;
}
//_________________________________________________________________________
double wignerUtils::getMaxX()
{
    return mMaxX;
}
//_________________________________________________________________________
double wignerUtils::getMinP()
{
    return mMinP;
}
//_________________________________________________________________________
double wignerUtils::getHCut()
{
    return mHCut;
}
//_________________________________________________________________________
double wignerUtils::getMaxP()
{
    return mMaxP;
}
//_________________________________________________________________________
void wignerUtils::setMinX(double minX)
{
    mMinX = minX;
}
//_________________________________________________________________________
void wignerUtils::setMaxX(double maxX)
{
    mMaxX = maxX;
}
//_________________________________________________________________________
void wignerUtils::setMinP(double minP)
{
    mMinP = minP;
}
//_________________________________________________________________________
void wignerUtils::setMaxP(double maxP)
{
    mMaxP = maxP;
}
//_________________________________________________________________________
void wignerUtils::setIntegrationRanges(double minX, double maxX, double minP, double maxP)
{
    setMinX(minX);
    setMaxX(maxX);
    setMinP(minP);
    setMaxP(maxP);
}
//_________________________________________________________________________
double wignerUtils::integral(TF2 *function, double minX, double maxX, double minP, double maxP)
{
    double res = 0;
    if (testMode == false)
    {
        for (float x = mDx / 2 + minX; x < maxX; x += mDx)
        {
            for (float p = mDp / 2 + minP; p < maxP; p += mDp)
            {
                res += function->Eval(x, p);
            }
        }
        res *= mDx * mDp;
    }
    else
    {
        res = function->Integral(minX, maxX, minP, maxP);
    }
    return res;
}