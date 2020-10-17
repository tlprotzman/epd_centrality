/**
 * \brief Plots the percentiles of the global as well as the
 *        predicted percentiles as created by quantiles.cpp
 * 
 *      As a note there is not a chance this doesn't have memory leakes in its current form
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
#include <TColor.h>
#include <TStyle.h>

#include <iostream> 

// Since data is stored as 5 percent ranges, it is convenient to be able to ask
// for a larger range in dynamic way, combinding them as needed.
TH1D *getQuantileRange(int min, int max, TDirectory *dir, const char *mode) {
    if (min % 5 != 0 || max % 5 != 0) {
        std::cerr << "Quantiles are stored in 5%% increments, cannot get desired range: " << min << "-" << max << std::endl;
        return nullptr;
    }

    TH1D *range = nullptr;

    int numQuantiles = (max - min) / 5;
    std::cout << "Loading " << numQuantiles << "bins" << std::endl;
    dir->GetObject(Form("%s_%d%%-%d%%;1", mode, min, min + 5), range);
    // std::cout << "form: " << Form("%s_%d%%-%d%%;1", mode, min, min + 5) << std::endl;
    if (range == nullptr) {
        std::cerr << "WTF" << std::endl;
    }
    
    for (uint32_t i = 1; i < numQuantiles; i++) {
        // std::cout << "This better not be running" << std::endl;
        TH1D *temp = nullptr;
        dir->GetObject(Form("%s_%d%%-%d%%;1", mode, min + 5 * i, min + 5 * (i + 1)), temp);
        // std::cout << "Form 2: " << Form("%s_%d%%-%d%%;1", mode, min + 5 * i, min + 5 * (i + 1)) << std::endl;
        if (temp == nullptr) {
            std::cerr << "WTF" << std::endl;
        }
        range->Add(temp);
        delete temp;
    }
    TH1D *retVal = (TH1D*)range->Clone(Form("%s_%d_%d", mode, min, max));
    return retVal;
}

// Recreates figure 11, comparing tpc quantiles to epd quantiles
void quantileComparison(TDirectory *dir, const char *name) {
    // dir->pwd();
    // dir->ls();
    // We want to compare the 0-5% (95-100%), 20-30% (70-80%), and 90-100% (0-10%) ranges
    TH1D *tpc_95_100, *epd_95_100;
    TH1D *tpc_65_70,  *epd_65_70;
    TH1D *tpc_15_20,  *epd_15_20;

    tpc_15_20 = (TH1D*)getQuantileRange(50, 60, dir, "tpc")->Clone();
    epd_15_20 = (TH1D*)getQuantileRange(50, 60, dir, "epd")->Clone();

    tpc_65_70 = (TH1D*)getQuantileRange(75, 85, dir, "tpc")->Clone();
    epd_65_70 = (TH1D*)getQuantileRange(75, 85, dir, "epd")->Clone();

    tpc_95_100 = (TH1D*)getQuantileRange(95, 100, dir, "tpc")->Clone();
    epd_95_100 = (TH1D*)getQuantileRange(95, 100, dir, "epd")->Clone();

    TH1D **quantileTPCProjects = (TH1D**)malloc(3 * sizeof(TH1D*));
    TH1D **quantileEPDProjects = (TH1D**)malloc(3 * sizeof(TH1D*));


    quantileTPCProjects[0] = tpc_15_20;
    quantileTPCProjects[1] = tpc_65_70;
    quantileTPCProjects[2] = tpc_95_100;

    quantileEPDProjects[0] = epd_15_20;
    quantileEPDProjects[1] = epd_65_70;
    quantileEPDProjects[2] = epd_95_100;

    // Plot results
    TCanvas *canvas = new TCanvas(name, name, 1000, 1000);
    canvas->SetLogy();
    tpc_15_20->SetMinimum(1);

    TLegend *legend = new TLegend(0.65, 0.78, 0.75, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(0.03);

    char *formatA = "hist l p";
    char *formatB = "same hist l p";
    char *format = formatB;
    quantileTPCProjects[0]->SetTitle(name);
    quantileTPCProjects[0]->SetStats(0);

    for (int32_t i = 0; i < 3; i++) {

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

        
    }
    legend->AddEntry(quantileTPCProjects[0], "RefMult", "l");
    legend->AddEntry(quantileEPDProjects[0], name, "l");
    legend->Draw();

    canvas->Draw();
    canvas->SaveAs(Form("histograms/%s.png", name));

}

// Compare the variance of 10% quantiles, as in figure 12 in the paper I am referencing
TGraph *quantileVarianceComparison(TDirectory *dir, const char *name) {
    const int numQuantiles = 10;
    int quantilesRange = 100 / numQuantiles;
    TH1D *epdQuantiles[numQuantiles];
    TH1D *tpcQuantiles[numQuantiles];

    TList *keys = dir->GetListOfKeys();
    for (uint32_t i = 0; i < numQuantiles; i++) {
        tpcQuantiles[i] = getQuantileRange(i * quantilesRange, (i + 1) * quantilesRange, dir, "tpc");
        epdQuantiles[i] = getQuantileRange(i * quantilesRange, (i + 1) * quantilesRange, dir, "epd");
    }

    double *tpcVariance = (double*)malloc(numQuantiles * sizeof(double));
    double *epdVariance = (double*)malloc(numQuantiles * sizeof(double));
    double *count = (double*)malloc(numQuantiles * sizeof(double)); //memory leaks here
    for (uint32_t i = 0; i < numQuantiles; i++) {
        tpcVariance[i] = tpcQuantiles[i]->GetRMS();
        if (tpcQuantiles[i]->GetRMS() < 0.0001) {
            epdVariance[i] = 1;
        }
        else {
            epdVariance[i] = epdQuantiles[i]->GetRMS() / tpcQuantiles[i]->GetRMS();
        }
        count[i] = i + 1;
    }
    TGraph *graphs;
    graphs = new TGraph(numQuantiles, count, epdVariance);
    graphs->SetTitle(name);

    return graphs;
}

void quantileSummary(char *infile="data/epd_tpc_relations.root") {
    TFile rootFile(infile);
    TDirectory *quantile_directory = rootFile.GetDirectory("quantiles");

    TList *methods = quantile_directory->GetListOfKeys();

    uint32_t numMethods = 3;
    TGraph *varianceGraphs[numMethods];    // WILL NOT SCALE AS IS

    int i = 0;
    for (TIter method = methods->begin(); method != methods->end(); ++method) {
        const char *methodName = (*method)->GetName();
        std::cout << methodName << std::endl;
        quantileComparison(quantile_directory->GetDirectory(methodName), methodName);
        varianceGraphs[i] = quantileVarianceComparison(quantile_directory->GetDirectory(methodName), methodName);
        i++;
    }

    // Plotting variance graph

    TCanvas c2("variances", "variances", 1000, 1000);
    TLegend *legend = new TLegend(0.65, 0.25, 0.75, 0.45);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(0.03);
    gPad->SetLogy();

    for (uint32_t i = 0; i < numMethods; i++) {
        varianceGraphs[i]->SetLineColor(i + 1);
        if (i == 0) {
            varianceGraphs[i]->Draw("AC*");
        }
        else {
            varianceGraphs[i]->Draw("C*");
        }
        legend->AddEntry(varianceGraphs[i], varianceGraphs[i]->GetTitle(), "l");
    }
    c2.Draw();
    legend->Draw();
    c2.SaveAs("histograms/variances.png");
    rootFile.Close();
}
