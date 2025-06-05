void wignersim(double range_start, double range_end, double increment, TString outfile)
{
    wignerSource *fw = new wignerSource;
    fw->initFunctions();
    fw->setRanges(0., 20, 0., 0.6);
    fw->setR0(1.);

    TFile *file = new TFile(outfile, "RECREATE");
    std::cout << "Creating " << outfile << "\n";

    TTree *tree = new TTree("tree", "W x W");

    double r0, WW, norm, wH, wK, wV, k, coal;
    tree->Branch("r0", &r0, "r0/D");
    tree->Branch("WxW", &WW, "WW/D");
    tree->Branch("coal", &coal, "coal/D");
    tree->Branch("norm", &norm, "norm/D");
    tree->Branch("wH", &wH, "wH/D");
    tree->Branch("wK", &wK, "wK/D");
    tree->Branch("wV", &wV, "wV/D");
    tree->Branch("k", &k, "k/D");

    auto w = fw->getWignerFunction();
    std::cout << "deuteron int : " << fw->getDeuteronInt() << "\n";

    for (double i = range_start; i < range_end; i += increment)
    {
        fw->setRadiusK(i);
        r0 = fw->getRadius();
        k = i;
        norm = fw->getNorm();
        WW = fw->checkWxW();
        wK = fw->getwK();
        wV = fw->getwV();
        wH = fw->getwH();
        coal = fw->getcoal();

        std::cout << "i : " << i 
                  << " coal: " << fw->getcoal() 
                  << " r0:  " << r0 
                  << " k*: " << k 
                  << " Norm: " << norm 
                  << " Check: " << WW 
                  << " K: " << wK 
                  << " V: " << wV 
                  << " H: " << wH << "\n";

        tree->Fill();
    }

    tree->Write();
    file->Close();
}
