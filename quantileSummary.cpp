/**
 * \brief Plots the percetiles of the global as well as the
 *        predicted percentiles as created by quantiles.cpp
 * 
 * \author Tristan Protzman
 * \date 15 October 2020
 * \email tlprotzman@gmail.com
 * \affiliation Lehigh University
 * 
 */

#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TGraph.h>

#include <iostream> 


// Recreates figure 11, comparing tpc quantiles to epd quantiles
void quantileComparison(TDirectory *dir, const char *name) {
    // dir->pwd();
    // dir->ls();
    // We want to compare the 0-5% (95-100%), 20-30% (70-80%), and 90-100% (0-10%) ranges
    TH1D *tpc_95_100, *epd_95_100;
    
    TH1D *tpc_70_80, *epd_70_80;
    TH1D *tpc_70_75, *epd_70_75;
    TH1D *tpc_75_80, *epd_75_80;
    
    TH1D *tpc_0_10, *epd_0_10;
    TH1D *tpc_0_5, *epd_0_5;
    TH1D *tpc_5_10, *epd_5_10;

    dir->GetObject("tpc_95%-100%;1", tpc_95_100);
    dir->GetObject("epd_95%-100%;1", epd_95_100);

    dir->GetObject("tpc_70%-75%;1", tpc_70_75);
    dir->GetObject("epd_70%-75%;1", epd_70_75);
    dir->GetObject("tpc_75%-80%;1", tpc_75_80);
    dir->GetObject("epd_75%-80%;1", epd_75_80);

    dir->GetObject("tpc_0%-5%;1", tpc_0_5);
    dir->GetObject("epd_0%-5%;1", epd_0_5);
    dir->GetObject("tpc_5%-10%;1", tpc_5_10);
    dir->GetObject("epd_5%-10%;1", epd_5_10);
    
    // Merge 70-75 and 75-80 ranges into 70-80 range
    tpc_70_80 = (TH1D*)tpc_70_75->Clone("tpc_70_80");
    tpc_70_80->Add(tpc_75_80);

    epd_70_80 = (TH1D*)epd_70_75->Clone("epd_70_80");
    epd_70_80->Add(epd_75_80);

    // Merge 0-5 and 5-10 ranges into 70-80 range
    tpc_0_10 = (TH1D*)tpc_0_5->Clone("tpc_0_10");
    tpc_0_10->Add(tpc_5_10);

    epd_0_10 = (TH1D*)epd_0_5->Clone("epd_0_10");
    epd_0_10->Add(epd_5_10);

    TH1D *quantileTPCProjects[3];
    TH1D *quantileEPDProjects[3];

    quantileTPCProjects[0] = tpc_0_10;
    quantileTPCProjects[1] = tpc_70_80;
    quantileTPCProjects[2] = tpc_95_100;

    quantileEPDProjects[0] = epd_0_10;
    quantileEPDProjects[1] = epd_70_80;
    quantileEPDProjects[2] = epd_95_100;

    // Plot results
    TCanvas canvas(name, name, 1000, 1000);
    char *formatA = "hist l p";
    char *formatB = "same hist l p";
    char *format = formatA;
    quantileTPCProjects[0]->SetTitle(name);
    quantileTPCProjects[0]->SetStats(0);

    for (uint32_t i = 0; i < 3; i++) {
        gPad->SetLogy();
        quantileTPCProjects[i]->SetLineColor(kRed);
        quantileTPCProjects[i]->SetMarkerColor(kRed);
        quantileTPCProjects[i]->SetMarkerStyle(kOpenTriangleUp);
        quantileTPCProjects[i]->SetMarkerSize(0.5);
        quantileTPCProjects[i]->Draw(format);
        format = formatB;

        quantileEPDProjects[i]->SetLineColor(kBlue);
        quantileEPDProjects[i]->SetMarkerColor(kBlue);
        quantileEPDProjects[i]->SetMarkerStyle(kStar);
        quantileEPDProjects[i]->SetMarkerSize(0.5);
        quantileEPDProjects[i]->Draw(format);

        TLegend *legend = new TLegend(0.65, 0.78, 0.75, 0.85);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.03);
        legend->AddEntry(quantileTPCProjects[i]->GetName(), "RefMult", "l");
        legend->AddEntry(quantileEPDProjects[i]->GetName(), name, "l");
        legend->Draw();
    }
    canvas.Draw();
    canvas.SaveAs(Form("histograms/%s.png", name));

}

// Compare the variance of 10% quantiles, as in figure 12 in the paper I am referencing
TGraph **quantileVarianceComparison(TDirectory *dir, const char *name) {
    TH1D *epdQuantiles[20];
    TH1D *tpcQuantiles[20];

    TList *keys = dir->GetListOfKeys();
    uint32_t i = 0;
    for (TIter itr = keys->begin(); itr != keys->end(); ++itr) {
        const char *key = (*itr)->GetName();
        dir->GetObject(key, tpcQuantiles[i]);   // Quick and dirty and bad but maybe it'll work for now?
        ++itr;
        key = (*itr)->GetName();
        dir->GetObject(key, epdQuantiles[i]);
        if (epdQuantiles[i] == nullptr || tpcQuantiles[i] == nullptr) {
            std::cerr << "it didn't work" << std::endl;
            return nullptr;
        }
        i++;
    }

    double *tpcVariance = (double*)malloc(20 * sizeof(double));
    double *epdVariance = (double*)malloc(20 * sizeof(double));
    double *count = (double*)malloc(20 * sizeof(double)); //memory leaks here
    for (uint32_t i = 0; i < 20; i++) {
        tpcVariance[i] = tpcQuantiles[i]->GetRMS();
        if (tpcQuantiles[i]->GetRMS() < 0.0001) {
            epdVariance[i] = 0;
        }
        else {
            epdVariance[i] = epdQuantiles[i]->GetRMS() / tpcQuantiles[i]->GetRMS();
        }
        count[i] = i + 1;
    }
    TGraph **graphs = (TGraph**)malloc(2 * sizeof(TGraph*));
    graphs[0] = new TGraph(20, count, tpcVariance);
    graphs[1] = new TGraph(20, count, epdVariance);
    graphs[1]->SetTitle(name);

    return graphs;
}

void quantileSummary(char *infile="data/epd_tpc_relations.root") {
    TFile rootFile(infile);
    TDirectory *quantile_directory = rootFile.GetDirectory("quantiles");

    TList *methods = quantile_directory->GetListOfKeys();

    TGraph **varianceGraphs[5];    // WILL NOT SCALE AS IS

    int i = 0;
    for (TIter method = methods->begin(); method != methods->end(); ++method) {
        const char *methodName = (*method)->GetName();
        std::cout << methodName << std::endl;
        // quantileComparison(quantile_directory->GetDirectory(methodName), methodName);
        varianceGraphs[i] = quantileVarianceComparison(quantile_directory->GetDirectory(methodName), methodName);
        i++;
    }

    TCanvas c2("variances", "variances", 1000, 1000);
    for (uint32_t i = 0; i < 5; i++) {
        if (i == 0) {
            varianceGraphs[i][1]->Draw("AC*");
        }
        else {
            varianceGraphs[i][1]->Draw("C*");
        }
    }
    // varianceGraphs[0][0]->Draw();
    c2.Draw();
    c2.SaveAs("histograms/variances.png");
    rootFile.Close();
}
