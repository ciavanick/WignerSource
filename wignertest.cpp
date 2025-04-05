#include "CWignerSource.h"

void wignertest(){
    wignerSource *fw = new wignerSource;
    fw->initFunctions();
    fw->setRanges(0.,50, 0., 0.5);

    TString s = "data4/data.root";

    TFile *file = new TFile(s, "RECREATE");
    std::cout<<"Creating "<<s<<"\n";

    TTree *tree = new TTree("tree", "W x W");

    double r0, WW, norm, wH, wK, wV, k;
    tree->Branch("r0", &r0, "r0/D");
    tree->Branch("WxW", &WW, "WW/D");
    tree->Branch("norm", &norm, "norm/D");
    tree->Branch("wH", &wH, "wH/D");
    tree->Branch("wK", &wK, "wK/D");
    tree->Branch("wV", &wV, "wV/D");
    tree->Branch("k", &k, "k/D");

    auto w = fw->getWignerFunction();
    for(double i = 0.01; i < 0.3; i = i + 0.01){
        fw->setRadiusK(i);
        r0 = fw->getRadius();
        k = fw->getKStar();
        norm = fw->getNorm();
        WW = fw->checkWxW();
        wK = fw->getwK();
        wV = fw->getwV();
        wH = fw->getwH();

        std::cout << "i : "<< i << " r0:  " << r0 << " k*: " << k << " Norm: " << norm <<" Check: "<< WW << " K: " << wK << " V: "<< wV << " H: " << wH << "\n";
        tree->Fill();
    }
    // Write and close the file
    tree->Write();
    file->Close();

}