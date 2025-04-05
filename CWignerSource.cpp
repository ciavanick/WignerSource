#include "CWignerSource.h"
#include "CWignerUtils.h"

void wignerSource::initFunctions(){
    mW = new TF2("w" + mName, wignerUtils::wignerSource, mRMin, mRMax, mPMin, mPMax, 3);
    mWxJ = new TF2("wxj" + mName, wignerUtils::jacobianFun, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 3);
    mWxJforItself = new TF2("mWxJforItself" + mName, wignerUtils::jacobianW2, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 3);
    mWxW = new TF2("wxw" + mName, wignerUtils::wignerSource2, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 3);
    mK = new TF2("K" + mName, wignerUtils::kineticEnergy, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 4);
    mV = new TF2("V" + mName, wignerUtils::potentialEnergy, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 6);
    mH = new TF2("H" + mName, wignerUtils::hamiltonian, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 6);
    mWK = new TF2("WxK" + mName, wignerUtils::wK, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 4);
    mWV = new TF2("WxV" + mName, wignerUtils::wV, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 6);
    mWH = new TF2("WxH" + mName, wignerUtils::wH, wignerUtils::getMinX(), wignerUtils::getMaxX(), wignerUtils::getMinP(), wignerUtils::getMaxP(), 6);
    setFunctionsParameters();
}
//_________________________________________________________________________
void wignerSource::setThreeParam(TF2 *function){
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mR0);
    function->SetParameter(2, mKStar);
}
//_________________________________________________________________________
void wignerSource::setFourParam(TF2 *function){
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mR0);
    function->SetParameter(2, mKStar);
    function->SetParameter(3, mMu);
}
//_________________________________________________________________________
void wignerSource::setSixParam(TF2 *function){
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mR0);
    function->SetParameter(2, mKStar);
    function->SetParameter(3, mMu);
    function->SetParameter(4, mRWidth);
    function->SetParameter(5, mV0);
}
//_________________________________________________________________________
void wignerSource::setFunctionsParameters(){
    setThreeParam(mWxJ);
    normalization();
    setThreeParam(mW);
    setThreeParam(mWxJforItself);
    setThreeParam(mWxW);
    setFourParam(mK);
    setSixParam(mV);
    setSixParam(mH);
    setFourParam(mWK);
    setSixParam(mWV);
    setSixParam(mWH);
}
//_________________________________________________________________________
void wignerSource::setRadius(double r0){
    mR0 = r0;
    mWxJ->SetParameter(0, 1.);
    reSetRadius();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
void wignerSource::setRadiusK(double k){
    mR0 = wignerUtils::rSource(k);
    mKStar = 0.;
    mWxJ->SetParameter(0, 1.);
    reSetRadius();
    reSetKStar();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
void wignerSource::setKstar(double k){
    mKStar = k;
    mWxJ->SetParameter(0, 1.);
    reSetKStar();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
void wignerSource::setRanges(double xmin, double ymin, double xmax, double ymax){
    if(xmin < wignerUtils::getMinX() || xmax > wignerUtils::getMaxX() || ymin < wignerUtils::getMinP() || ymax > wignerUtils::getMaxP()){
        std::cout << "Invalid ranges specified; previous ranges kept.\n";
    }else{
        mW->SetRange(xmin, ymin, xmax, ymax);
        mRMin = xmin;
        mRMax = xmax;
        mPMin = ymin;
        mPMax = ymax;
    }

}
//_________________________________________________________________________
void wignerSource::setMu(double mu){
    mMu = mu;
    reSetMu();
}
//_________________________________________________________________________
void wignerSource::setRWidth(double rWidth){
    mRWidth = rWidth;
    reSetRWidth();
}
void wignerSource::setV0(double v0){
    mV0 = v0;
    reSetV0();
}
//_________________________________________________________________________
TF2* wignerSource::getWignerFunction(){
    return mW;
}
//_________________________________________________________________________
TF2* wignerSource::getWignerFunctionForItself(){
    return mWxW;
}
//_________________________________________________________________________
TF2* wignerSource::getWignerFunctionForJacobian(){
    return mWxJ;
}
//_________________________________________________________________________
TF2* wignerSource::getWignerFunction2ForJacobian(){
    return mWxJforItself;
}
//_________________________________________________________________________
TF2* wignerSource::getKineticEnergyFunction(){
    return mK;
}
//_________________________________________________________________________
TF2* wignerSource::getPotentialEnergyFunction(){
    return mV;
}
//_________________________________________________________________________
TF2* wignerSource::getHamiltonianFunction(){
    return mH;
}
//_________________________________________________________________________
TF2* wignerSource::getWignerKinetic(){
    return mWK;
}
//_________________________________________________________________________
TF2* wignerSource::getWignerPotential(){
    return mWV;
}
//_________________________________________________________________________
TF2* wignerSource::getWignerHamiltonan(){
    return mWH;
}
//_________________________________________________________________________
void wignerSource::normalization(){
    double norm = 1./wignerUtils::integral(mWxJ);
    mNorm = norm;
}
//_________________________________________________________________________
double wignerSource::getNorm(){
    return mNorm;
}
//_________________________________________________________________________
double wignerSource::getwK(){
    return wignerUtils::integral(mWK);
}
//_________________________________________________________________________
double wignerSource::getwV(){
    return wignerUtils::integral(mWV);
}
//_________________________________________________________________________
double wignerSource::getwH(){
    return wignerUtils::integral(mWH);
}
//_________________________________________________________________________
double wignerSource::checkWxW(){
    return wignerUtils::integral(mWxW) * (4 * wignerUtils::getHCut() * TMath::Pi()) * (4 * wignerUtils::getHCut() * TMath::Pi()) * (4 * wignerUtils::getHCut() * TMath::Pi());
}
//_________________________________________________________________________
double wignerSource::getRadius(){
    return mR0;
}
//_________________________________________________________________________
double wignerSource::getKStar(){
    return mKStar;
}
//_________________________________________________________________________
double wignerSource::getRMin(){
    return mRMin;
}
//_________________________________________________________________________
double wignerSource::getRMax(){
    return mRMax;
}
//_________________________________________________________________________
double wignerSource::getPMin(){
    return mPMin;
}
//_________________________________________________________________________
double wignerSource::getPMax(){
    return mPMax;
}
//_________________________________________________________________________
double wignerSource::getMu(){
    return mMu;
}
//_________________________________________________________________________
double wignerSource::getRWidth(){
    return mRWidth;
}
//_________________________________________________________________________
double wignerSource::getV0(){
    return mV0;
}
//_________________________________________________________________________
void wignerSource::reSetNorm(){
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
}
//_________________________________________________________________________
void wignerSource::reSetRadius(){
    mW->SetParameter(1, mR0);
    mWxJ->SetParameter(1, mR0);
    mWxJforItself->SetParameter(1, mR0);
    mWxW->SetParameter(1, mR0);
    mK->SetParameter(1, mR0);
    mV->SetParameter(1, mR0);
    mH->SetParameter(1, mR0);
    mWK->SetParameter(1, mR0);
    mWV->SetParameter(1, mR0);
    mWH->SetParameter(1, mR0);
}
//_________________________________________________________________________
void wignerSource::reSetKStar(){
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
}
//_________________________________________________________________________
void wignerSource::reSetMu(){
    mK->SetParameter(3, mMu);
    mV->SetParameter(3, mMu);
    mH->SetParameter(3, mMu);
    mWK->SetParameter(3, mMu);
    mWV->SetParameter(3, mMu);
    mWH->SetParameter(3, mMu);
}
//_________________________________________________________________________
void wignerSource::reSetRWidth(){
    mK->SetParameter(4, mRWidth);
    mV->SetParameter(4, mRWidth);
    mH->SetParameter(4, mRWidth);
    mWK->SetParameter(4, mRWidth);
    mWV->SetParameter(4, mRWidth);
    mWH->SetParameter(4, mRWidth);
}
//_________________________________________________________________________
void wignerSource::reSetV0(){
    mK->SetParameter(5, mV0);
    mV->SetParameter(5, mV0);
    mH->SetParameter(5, mV0);
    mWK->SetParameter(5, mV0);
    mWV->SetParameter(5, mV0);
    mWH->SetParameter(5, mV0);
}
