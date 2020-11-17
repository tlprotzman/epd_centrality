/**
 * \brief plotting the multiplicty found in the tpc against that of the tof detector
 * 
 * We hope to be able to use this as a selector of events to avoid pile up in the
 * tpc detector
 * 
 * \author Tristan Protzman
 * \date September 24, 2020
 * \email tlprotzman@gmail.com
 * 
 */

#include <stdlib.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVectorD.h>

void tpcVsTofSelection(const char *inFileName = "data/detector_data.root") {
    TFile inFile(inFileName);
    TVectorD *tpc;
    TVectorD *tof;

    inFile.GetObject("tpc_multiplicity", tpc);
    inFile.GetObject("tof_multiplicity", tof);
    inFile.Close();


    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetStatX(0.38);
    gStyle->SetStatY(0.85);

    uint32_t tpcBins = 150;
    int32_t tpcMin = 0;
    int32_t tpcMax = 300;

    uint32_t tofBins = 150;
    int32_t tofMin = 0;
    int32_t tofMax = 450; 

    TH2D *tofVsTpc = new TH2D("TPC vs TOF", "2D Histo;b TOF Tray Multiplicity;TPC RefMult",
                                  tpcBins, tpcMin, tpcMax,
                                  tofBins, tofMin, tofMax);

    int selectionWindow = 50;
    TH2D *windowTofVsTpc = new TH2D("Windowed", "2D Histo;TPC RefMult;b TOF Tray Multiplicity",
                                        tofBins, tofMin, tofMax,
                                        tpcBins, tpcMin, tpcMax);

    float tolerance1 = 0.8;
    TH2D *toleranceTofVsTpc = new TH2D("tolerance", "2D Histo;b TOF Tray Multiplicity;TPC RefMult",
                                        tofBins, tofMin, tofMax,
                                        tpcBins, tpcMin, tpcMax);

    float percentDifference = 0.8;
    TH2D *tolerance2TofVsTpc = new TH2D("percent difference", "2D Histo;b TOF Tray Multiplicity;TPC RefMult",
                                    tofBins, tofMin, tofMax,
                                    tpcBins, tpcMin, tpcMax);



    for (uint32_t i = 0; i < tpc->GetNrows(); i++) {
        double tpcVal, tofVal;
        tpcVal = (*tpc)[i];
        tofVal = (*tof)[i] / 2;
        tofVsTpc->Fill(tpcVal, tofVal);
        // Windowed
        if (abs((2 * tpcVal) - tofVal) < selectionWindow) {
            windowTofVsTpc->Fill(tofVal, tpcVal);
        }
        // Tolerance
        if (tofVal * (1 - tolerance1) < 2 * tpcVal && tofVal * (1 + tolerance1) > 2 * tpcVal) {
            toleranceTofVsTpc->Fill(tofVal, tpcVal);
        }
        // Percent Difference
        if (abs(((2 * tpcVal) - tofVal) / (((2 * tpcVal) + tofVal) / 2)) < percentDifference) {
            tolerance2TofVsTpc->Fill(tofVal, tpcVal);
        }
    }

    TCanvas *canvas = new TCanvas("canvas", "canvas");
    // canvas->Divide(2, 2);
    
    canvas->cd(1);
    gPad->SetLogz();
    tofVsTpc->SetTitle("TofMult vs RefMult1 , 7.7 GeV");
    tofVsTpc->Draw("Colz");

    // canvas->cd(2);
    // gPad->SetLogz();
    // windowTofVsTpc->SetTitle(Form("TPC Multiplicity vs TOF Multiplicity, 7.7 GeV, Selection window +-%d", selectionWindow));
    // windowTofVsTpc->Draw("Colz");

    // canvas->cd(3);
    // gPad->SetLogz();
    // toleranceTofVsTpc->SetTitle(Form("TPC Multiplicity vs TOF Multiplicity, 7.7 GeV, Tolerance +-%.f%%", tolerance1 * 100));
    // toleranceTofVsTpc->Draw("Colz");

    // canvas->cd(4);
    // gPad->SetLogz();
    // tolerance2TofVsTpc->SetTitle(Form("TPC Multiplicity vs TOF Multiplicity, 7.7 GeV, Max Percent Difference +-%.f%%", percentDifference * 100));
    // tolerance2TofVsTpc->Draw("Colz");
    
    canvas->Draw();
}
