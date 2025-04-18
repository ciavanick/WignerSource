#include "CWignerSource.h"
#include "CWignerUtils.h"

void wignerSource::initFunctions()
{
    mW = new TF2("w" + mName, wignerUtils::wignerSource, mRMin, mRMax, mPMin, mPMax, 3);
    mWxJ = new TF2("wxj" + mName, wignerUtils::jacobianFun, mRMin, mRMax, mPMin, mPMax, 3);
    mWxJforItself = new TF2("mWxJforItself" + mName, wignerUtils::jacobianW2, mRMin, mRMax, mPMin, mPMax, 3);
    mWxW = new TF2("wxw" + mName, wignerUtils::wignerSource2, mRMin, mRMax, mPMin, mPMax, 3);
    mK = new TF2("K" + mName, wignerUtils::kineticEnergy, mRMin, mRMax, mPMin, mPMax, 4);
    mV = new TF2("V" + mName, wignerUtils::potentialEnergy, mRMin, mRMax, mPMin, mPMax, 6);
    mH = new TF2("H" + mName, wignerUtils::hamiltonian, mRMin, mRMax, mPMin, mPMax, 6);
    mWK = new TF2("WxK" + mName, wignerUtils::wK, mRMin, mRMax, mPMin, mPMax, 4);
    mWV = new TF2("WxV" + mName, wignerUtils::wV, mRMin, mRMax, mPMin, mPMax, 6);
    mWH = new TF2("WxH" + mName, wignerUtils::wH, mRMin, mRMax, mPMin, mPMax, 6);
    mD = new TF2("WD" + mName, wignerUtils::wignerDeuteron, mRMin, mRMax, mPMin, mPMax, 0);
    mDInt = new TF2("WDInt" + mName, wignerUtils::wignerDeuteronIntegral, mRMin, mRMax, mPMin, mPMax, 0);
    mC = new TF2("CoalescenceProb" + mName, wignerUtils::coalescenceProbability, mRMin, mRMax, mPMin, mPMax, 3);
    setFunctionsParameters();
}
//_________________________________________________________________________
void wignerSource::setThreeParam(TF2 *function)
{
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mRadius);
    function->SetParameter(2, mKStar);
}
//_________________________________________________________________________
void wignerSource::setFourParam(TF2 *function)
{
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mRadius);
    function->SetParameter(2, mKStar);
    function->SetParameter(3, mMu);
}
//_________________________________________________________________________
void wignerSource::setSixParam(TF2 *function)
{
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mRadius);
    function->SetParameter(2, mKStar);
    function->SetParameter(3, mMu);
    function->SetParameter(4, mRWidth);
    function->SetParameter(5, mV0);
}
//_________________________________________________________________________
void wignerSource::setFunctionsParameters()
{
    setThreeParam(mWxJ);
    normalization();
    setThreeParam(mW);
    setThreeParam(mWxJforItself);
    setThreeParam(mWxW);
    setThreeParam(mC);
    setFourParam(mK);
    setSixParam(mV);
    setSixParam(mH);
    setFourParam(mWK);
    setSixParam(mWV);
    setSixParam(mWH);
}
//_________________________________________________________________________
void wignerSource::setRadius(double radius)
{
    mRadius = radius;
    mWxJ->SetParameter(0, 1.);
    reSetRadius();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
void wignerSource::setR0(double r0)
{
    mR0 = r0;
}
//_________________________________________________________________________
void wignerSource::setRadiusK(double k)
{
    mKin = k;
    mRadius = wignerUtils::radius(k, mR0);
    mKStar = wignerUtils::kStarEff(k, mRadius);
    mWxJ->SetParameter(0, 1.);
    reSetRadius();
    reSetKStar();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
// to check if it is useful
void wignerSource::setKstar(double k)
{
    mKStar = wignerUtils::kStarEff(k, mRadius);
    mWxJ->SetParameter(0, 1.);
    reSetKStar();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
void wignerSource::setKIn(double k)
{
    mKin = k;
}
//_________________________________________________________________________
void wignerSource::setRanges(double xmin, double ymin, double xmax, double ymax)
{
    if (xmin < wignerUtils::getMinX() || xmax > wignerUtils::getMaxX() || ymin < wignerUtils::getMinP() || ymax > wignerUtils::getMaxP())
    {
        std::cout << "Invalid ranges specified; previous ranges kept.\n";
    }
    else
    {
        mW->SetRange(xmin, ymin, xmax, ymax);
        mRMin = xmin;
        mRMax = xmax;
        mPMin = ymin;
        mPMax = ymax;
    }
}
//_________________________________________________________________________
void wignerSource::setMu(double mu)
{
    mMu = mu;
    reSetMu();
}
//_________________________________________________________________________
void wignerSource::setRWidth(double rWidth)
{
    mRWidth = rWidth;
    reSetRWidth();
}
void wignerSource::setV0(double v0)
{
    mV0 = v0;
    reSetV0();
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunction()
{
    return mW;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunctionForItself()
{
    return mWxW;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunctionForJacobian()
{
    return mWxJ;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerFunction2ForJacobian()
{
    return mWxJforItself;
}
//_________________________________________________________________________
TF2 *wignerSource::getKineticEnergyFunction()
{
    return mK;
}
//_________________________________________________________________________
TF2 *wignerSource::getPotentialEnergyFunction()
{
    return mV;
}
//_________________________________________________________________________
TF2 *wignerSource::getHamiltonianFunction()
{
    return mH;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerKinetic()
{
    return mWK;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerPotential()
{
    return mWV;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerHamiltonan()
{
    return mWH;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerDeuteron()
{
    return mD;
}
//_________________________________________________________________________
TF2 *wignerSource::getWignerDeuteronIntegral()
{
    return mDInt;
}
//_________________________________________________________________________
TF2 *wignerSource::getCoalescenceProbability()
{
    return mC;
}
//_________________________________________________________________________
void wignerSource::normalization()
{
    double norm = 1. / wignerUtils::integral(mWxJ, 0., TMath::Max(5. * mRadius, 20.), 0., 0.6);
    mNorm = norm;
}
//_________________________________________________________________________
double wignerSource::getNorm()
{
    return mNorm;
}
//_________________________________________________________________________
double wignerSource::getwK()
{
    return wignerUtils::integral(mWK);
}
//_________________________________________________________________________
double wignerSource::getwV()
{
    return wignerUtils::integral(mWV);
}
//_________________________________________________________________________
double wignerSource::getwH()
{
    return wignerUtils::integral(mWH);
}
//_________________________________________________________________________
double wignerSource::checkWxW()
{
    return wignerUtils::integral(mWxW) * (wignerUtils::getHCut() * 2 * TMath::Pi()) * (wignerUtils::getHCut() * 2 * TMath::Pi()) * (wignerUtils::getHCut() * 2 * TMath::Pi());
}
//_________________________________________________________________________
double wignerSource::getRadius()
{
    return mRadius;
}
//_________________________________________________________________________
double wignerSource::getKStar()
{
    return mKStar;
}
//_________________________________________________________________________
double wignerSource::getRMin()
{
    return mRMin;
}
//_________________________________________________________________________
double wignerSource::getRMax()
{
    return mRMax;
}
//_________________________________________________________________________
double wignerSource::getPMin()
{
    return mPMin;
}
//_________________________________________________________________________
double wignerSource::getPMax()
{
    return mPMax;
}
//_________________________________________________________________________
double wignerSource::getMu()
{
    return mMu;
}
//_________________________________________________________________________
double wignerSource::getRWidth()
{
    return mRWidth;
}
//_________________________________________________________________________
double wignerSource::getV0()
{
    return mV0;
}
//_________________________________________________________________________
double wignerSource::getcoal()
{
    return wignerUtils::integral(mC)*(wignerUtils::getHCut()*2*TMath::Pi())*(wignerUtils::getHCut()*2*TMath::Pi())*(wignerUtils::getHCut()*2*TMath::Pi());
}
//_________________________________________________________________________
double wignerSource::getDeuteronInt()
{
    return wignerUtils::integral(mDInt);
}
//_________________________________________________________________________
void wignerSource::reSetNorm()
{
    mW->SetParameter(0, mNorm);
    mWxJ->SetParameter(0, mNorm);
    mWxJforItself->SetParameter(0, mNorm);
    mWxW->SetParameter(0, mNorm);
    mK->SetParameter(0, mNorm);
    mV->SetParameter(0, mNorm);
    mH->SetParameter(0, mNorm);
    mWK->SetParameter(0, mNorm);
    mWV->SetParameter(0, mNorm);
    mWH->SetParameter(0, mNorm);
    mC->SetParameter(0, mNorm);
}
//_________________________________________________________________________
void wignerSource::reSetRadius()
{
    mW->SetParameter(1, mRadius);
    mWxJ->SetParameter(1, mRadius);
    mWxJforItself->SetParameter(1, mRadius);
    mWxW->SetParameter(1, mRadius);
    mK->SetParameter(1, mRadius);
    mV->SetParameter(1, mRadius);
    mH->SetParameter(1, mRadius);
    mWK->SetParameter(1, mRadius);
    mWV->SetParameter(1, mRadius);
    mWH->SetParameter(1, mRadius);
    mC->SetParameter(1, mRadius);
}
//_________________________________________________________________________
void wignerSource::reSetKStar()
{
    mW->SetParameter(2, mKStar);
    mWxJ->SetParameter(2, mKStar);
    mWxJforItself->SetParameter(2, mKStar);
    mWxW->SetParameter(2, mKStar);
    mK->SetParameter(2, mKStar);
    mV->SetParameter(2, mKStar);
    mH->SetParameter(2, mKStar);
    mWK->SetParameter(2, mKStar);
    mWV->SetParameter(2, mKStar);
    mWH->SetParameter(2, mKStar);
    mC->SetParameter(2, mKStar);
}
//_________________________________________________________________________
void wignerSource::reSetMu()
{
    mK->SetParameter(3, mMu);
    mV->SetParameter(3, mMu);
    mH->SetParameter(3, mMu);
    mWK->SetParameter(3, mMu);
    mWV->SetParameter(3, mMu);
    mWH->SetParameter(3, mMu);
}
//_________________________________________________________________________
void wignerSource::reSetRWidth()
{
    mK->SetParameter(4, mRWidth);
    mV->SetParameter(4, mRWidth);
    mH->SetParameter(4, mRWidth);
    mWK->SetParameter(4, mRWidth);
    mWV->SetParameter(4, mRWidth);
    mWH->SetParameter(4, mRWidth);
}
//_________________________________________________________________________
void wignerSource::reSetV0()
{
    mK->SetParameter(5, mV0);
    mV->SetParameter(5, mV0);
    mH->SetParameter(5, mV0);
    mWK->SetParameter(5, mV0);
    mWV->SetParameter(5, mV0);
    mWH->SetParameter(5, mV0);
}
