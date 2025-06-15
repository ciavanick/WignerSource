/// @addtogroup WignerMacros
/// @{
/**
 * @brief Generate and save plots from Wigner simulation output stored in a ROOT file.
 *
 * This macro opens a ROOT file containing a TTree named `"tree"` and generates
 * TGraph plots of coalescence probability, kinetic/potential/Hamiltonian energy vs. k* and r₀.
 * Each plot is saved as both a `.pdf` and `.root` file in the given output folder.
 *
 * Required tree branches:
 * - "coal": deuteron coalescence probability
 * - "k":    relative momentum (k*)
 * - "r0":   effective source radius
 * - "wH":   Wigner-weighted Hamiltonian
 * - "wK":   Wigner-weighted kinetic energy
 * - "wV":   Wigner-weighted potential energy
 *
 * @param folder    Path to the output folder where plots will be saved (e.g., "simtest4").
 * @param filename  Name of the input ROOT file inside the folder (e.g., "data_merged.root").
 *
 * @note If any required branch or the TTree is missing, an error is printed and the function returns.
 * @note A cut on r₀ is applied (0 ≤ r₀ ≤ 10 fm) when plotting energies vs. r₀.
 */
void makeplots(const char* folder, const char* filename)
{
    // Construct the full file path
    TString filepath = TString::Format("%s/%s", folder, filename);

    TFile *file = TFile::Open(filepath);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: could not open file " << filepath << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: TTree 'tree' not found in file " << filepath << std::endl;
        return;
    }

    // Check required branches
    const char* required_branches[] = {"coal", "k", "r0", "wH", "wK", "wV"};
    for (auto br : required_branches) {
        if (!tree->GetBranch(br)) {
            std::cerr << "Error: '" << br << "' branch not found in the tree." << std::endl;
            return;
        }
    }

    auto draw_and_save = [&](const char* yvar, const char* xvar, const char* title,
                             const char* xaxis, const char* yaxis, const char* base_filename,
                             const char* cut = "")
    {
        TCanvas *c = new TCanvas(Form("c_%s_vs_%s", yvar, xvar), title, 800, 600);
        tree->Draw(Form("%s:%s", yvar, xvar), cut, "AP");
        gPad->SetGrid();
        TGraph *gr = (TGraph*)gPad->GetPrimitive("Graph");
        if (gr) {
            gr->SetTitle(Form("%s;%s;%s", title, xaxis, yaxis));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(1.2);

            TString pdf_path = TString::Format("%s/%s.pdf", folder, base_filename);
            c->SaveAs(pdf_path);

            TString root_path = TString::Format("%s/%s.root", folder, base_filename);
            TFile *outfile = new TFile(root_path, "RECREATE");
            gr->Write("graph");
            outfile->Close();
        } else {
            std::cerr << "Warning: Could not retrieve graph for " << yvar << " vs " << xvar << std::endl;
        }

        delete c;
    };

    draw_and_save("coal", "k",  "Coalescence Probability vs. k*", "k* (GeV/c)", "P_{coal}", "coal_vs_k");
    draw_and_save("wH",   "k",  "Hamiltonian vs. k*", "k* (GeV/c)", "H (GeV)", "hamiltonian_vs_k");
    draw_and_save("wK",   "k",  "Kinetic Energy vs. k*", "k* (GeV/c)", "Kinetic Energy (GeV)", "kinetic_vs_k");
    draw_and_save("wV",   "k",  "Potential Energy vs. k*", "k* (GeV/c)", "Potential Energy (GeV)", "potential_vs_k");

    const char* r0_cut = "r0 >= 0 && r0 <= 10";
    draw_and_save("wH",   "r0", "Hamiltonian vs. r_{0}", "r_{0} (fm)", "H (GeV)", "hamiltonian_vs_r0", r0_cut);
    draw_and_save("wK",   "r0", "Kinetic Energy vs. r_{0}", "r_{0} (fm)", "Kinetic Energy (GeV)", "kinetic_vs_r0", r0_cut);
    draw_and_save("wV",   "r0", "Potential Energy vs. r_{0}", "r_{0} (fm)", "Potential Energy (GeV)", "potential_vs_r0", r0_cut);

    file->Close();
}

/** @} */  // End of WignerMacros group