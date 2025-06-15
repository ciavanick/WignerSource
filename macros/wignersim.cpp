/**
 * @defgroup WignerMacros ROOT Macros
 * @brief Functions for Wigner simulation and plotting.
 *
 * This group contains standalone utilities for computing Wigner observables
 * and visualizing simulation results using ROOT.
 * @{
 */
/**
 * @brief Run Wigner simulation over a k* range and store results in a ROOT file.
 *
 * This function performs a scan over a range of relative momenta (k*) and computes:
 * - Effective source radius
 * - Normalization factor
 * - Wigner function normalization check (WxW)
 * - Wigner-weighted kinetic, potential, and total energy
 * - Deuteron coalescence probability
 *
 * All computed values are written to a ROOT TTree stored in the specified output file.
 *
 * @param range_start Starting value of k* (must be >= 0).
 * @param range_end   Ending value of k* (must be >= range_start).
 * @param increment   Step size in k*.
 * @param outfile     Output ROOT file name to store the TTree.
 * @param txtinput    Path to the input text file with Wigner source parameters.
 */

void wignersim(double range_start, double range_end, double increment, TString outfile, const std::string &txtinput)
{
    wignerSource *fw = new wignerSource;
    fw->initFunctions();
    fw->SetFromTxt(txtinput);

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

    if (range_start < 0 || range_end < 0 || range_start > range_end)
    {
        std::cerr << "invalid range of k\n";
        return;
    }
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
/** @} */