
void wignertest2()
{
    wignerSource *fw = new wignerSource;
    fw->initFunctions();
    fw->setRanges(0., 20, 0., 0.6);
    fw->setR0(1.);

    TString s = "data4/data.root";

    TFile *file = new TFile(s, "RECREATE");
    std::cout << "Creating " << s << "\n";

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
    std::cout << "deuteron int : " << fw->getDeuteronInt() << "\n";
    float xk[1000], xP[1000], xK[1000], xE[1000], xV[1000], xkf[1000];
    int np = 0;
    double mRed = 0.938 / 2;
    double mRdeu = 3.2;
    for (double i = 1E-3 / 2; i < 0.1; i += 2E-3)
    {
        xk[np] = i;
        fw->setRadiusK(i);
        r0 = fw->getRadius();
        k = i;
        norm = fw->getNorm();
        WW = fw->checkWxW();
        double intW2 = fw->checkWxW();
        wK = fw->getwK();
        double intWkin = fw->getwK();
        wV = fw->getwV();
        double intWpot = fw->getwV();
        wH = fw->getwH();

        xK[np] = i * i * 0.5 / mRed;
        xV[np] = intWpot;

        if (xK[np] > xV[np])
        {
            xE[np] = xK[np] + xV[np];
            xkf[np] = sqrt(xE[np] * 2 * mRed);
        }
        else
        {
            xE[np] = 0;
            xkf[np] = 0;
        }

        if (intW2 < 0.00001)
        {
            intW2 = 1;
        }

        xP[np] = fw->getcoal();

        std::cout << "i : " << i << " coal: " << fw->getcoal() << " r0:  " << r0 << " k*: " << k << " Norm: " << norm << " Check: " << WW << " K: " << wK << " V: " << wV << " H: " << wH << "\n";
        tree->Fill();
        ++np;
    }
    // Write and close the file
    tree->Write();
    file->Close();

    int clean_np = 0;
    double xk_clean[1000], xP_clean[1000], xV_clean[1000]; // and others

    for (int i = 0; i < np; ++i)
    {
        if (xkf[i] == 0 || xP[i] == 0)
            continue;
        xk_clean[clean_np] = xk[i];
        xP_clean[clean_np] = xP[i];
        xV_clean[clean_np] = xV[i];
        clean_np++;
    }

    TCanvas *c = new TCanvas("wignertest","wignertest");
    c->Divide(2, 2);
    c->cd(1);
    TGraph *g = new TGraph(np, xk_clean, xP_clean);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    g->SetTitle(Form("Coal prob R_{deuteron} = %.1f fm;k*_{0} (GeV/c);P(k*_{0})", mRdeu));
    g->SetMinimum(0);
    //g->SetMaximum(0.8);
    c->cd(2);
    g = new TGraph(np, xk, xkf);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    g->SetTitle("Effect of interaction;k*_{0} (GeV/c);k*_{f} (GeV/c)");
    TLine *l = new TLine(0, 0, 0.5, 0.5);
    l->SetLineColor(2);
    l->Draw("SAME");
    c->cd(3);
    g = new TGraph(np, xK, xE);
    g->SetMarkerStyle(20);
    g->SetTitle("Effect of interaction;K.E._{0} (GeV);E_{f} (GeV)");
    g->Draw("AP");
    l = new TLine(0, 0, 0.25, 0.25);
    l->SetLineColor(2);
    l->Draw("SAME");
    c->cd(4);
    g = new TGraph(np, xk_clean, xV_clean);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    g->SetTitle("V interaction;k*_{0} (GeV/c);V (GeV)");
}