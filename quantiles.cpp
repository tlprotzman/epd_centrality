/**
 * \brief Plots the percetiles of the global as well as the
 *        predicted percentiles
 * 
 * \author Tristan Protzman
 * \date 6 October 2020
 * \email tlprotzman@gmail.com
 * \affiliation Lehigh University
 * 
 */


#include <TROOT.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphQQ.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <stdint.h>

#include "quantiles.h"

const int32_t numberQuantiles = 100;

// getQuantileRange creates a histogram of the data in the specified quantile range
TH2D *getQuantileRange(int lowerQuantile, int upperQuantile, TH2D *histogram, double *quantiles) {
    TH2D *quantileHistogram = new TH2D(Form("quantiles%d_%d", lowerQuantile, upperQuantile), Form("%d%% - %d%% Multiplicity", lowerQuantile, upperQuantile),
            histogram->GetNbinsX(), histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax(),
            histogram->GetNbinsY(), histogram->GetYaxis()->GetXmin(), histogram->GetYaxis()->GetXmax());


    quantileHistogram->SetXTitle("Multiplicity");
    quantileHistogram->SetYTitle("Counts");

    double min, max;
    lowerQuantile -= 1;      // Since the quantile array gives the upper bound of the quantile, the lower bound should be the one below
    if (lowerQuantile < 0) { // the requested quantile
        min = histogram->GetXaxis()->GetXmin(); // If we want from 0, set lower bound to be the minimum
    }
    else {
        min = quantiles[lowerQuantile];
    }
    if (upperQuantile > numberQuantiles - 1) {
        max = histogram->GetXaxis()->GetXmax();
    }
    else{
        max = quantiles[upperQuantile - 1];
    }
    double *centersX = (double*)malloc(sizeof(double) * histogram->GetNbinsX());
    double *centersY = (double*)malloc(sizeof(double) * histogram->GetNbinsY());
    histogram->GetXaxis()->GetCenter(centersX);
    histogram->GetYaxis()->GetCenter(centersY);
    std::cout << "Selecting events between " << min << "  and " << max << std::endl;
    for (uint32_t i = 0; i < histogram->GetNbinsX(); i++) {
        if (min < centersX[i] && centersX[i] < max) {
            for (uint32_t j = 0; j < histogram->GetNbinsY(); j++) {
                for (uint32_t k = 0; k < histogram->GetBinContent(i, j); k++) {
                    quantileHistogram->Fill(centersX[i], centersY[j]);
                }
            }
        }
    }
    free(centersX);
    free(centersY);
    return quantileHistogram;
}

void quantileAnalysis(TH2D *histogram) {
    // Storing the quantiles
    double xxQuantiles[numberQuantiles];
    double yxQuantiles[numberQuantiles];
    double xyQuantiles[numberQuantiles];
    double yyQuantiles[numberQuantiles];

    TH1D *xProjection = histogram->ProjectionX();
    for (uint32_t i = 0; i < numberQuantiles; i++) {
        xxQuantiles[i] = double(i + 1) / numberQuantiles;
        yxQuantiles[i] = double(i + 1) / numberQuantiles;
        xyQuantiles[i] = 0;
        xyQuantiles[i] = 0;
    }
    TH1D *yProjection = histogram->ProjectionY();
    xProjection->GetQuantiles(numberQuantiles, xyQuantiles, xxQuantiles);
    yProjection->GetQuantiles(numberQuantiles, yyQuantiles, yxQuantiles);

    for (uint32_t i = 0; i < numberQuantiles; i++) {
        std::cout << "Quantile " << i << ": " << xyQuantiles[i] << std::endl;
    }

    TH2D *quantileHistogram = getQuantileRange(50, 75, histogram, xyQuantiles);
    TH1D *quantileXProject = quantileHistogram->ProjectionX();
    TH1D *quantileYProject = quantileHistogram->ProjectionY();

    TGraph *xgraph = new TGraph(numberQuantiles, xxQuantiles, xyQuantiles);
    TGraph *ygraph = new TGraph(numberQuantiles, yxQuantiles, yyQuantiles);
    TGraphQQ *qq = new TGraphQQ(numberQuantiles, yyQuantiles, numberQuantiles, xyQuantiles);
    TCanvas *canvas = new TCanvas("canvas", "testing");
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(11);
    canvas->Divide(2, 4);

    {
    canvas->cd(1);
    gPad->SetLogz();
    histogram->Draw("Colz");

    // canvas->cd(2);
    // gPad->SetGrid();
    // qq->SetTitle("Quantile-Quantile Graph");
    // qq->SetMarkerStyle(21);
    // qq->SetMarkerColor(kRed);
    // qq->SetMarkerSize(0.5);
    // qq->Draw("ap");
    
    canvas->cd(3);
    xProjection->SetTitle("TPC Projection");
    xProjection->SetXTitle("TPC RefMult");
    xProjection->SetYTitle("Counts");
    xProjection->Draw();
    
    canvas->cd(4);
    yProjection->SetTitle("Linear Weights Projection");
    yProjection->SetXTitle("Linear Weight Multiplicity");
    yProjection->SetYTitle("Counts");
    yProjection->Draw();
    
    canvas->cd(5);
    xgraph->SetMarkerStyle(21);
    xgraph->SetMarkerColor(kRed);
    xgraph->SetMarkerSize(0.5);
    xgraph->SetTitle("TPC Quantiles");
    xgraph->GetXaxis()->SetTitle("Quantile");
    xgraph->GetYaxis()->SetTitle("TPC RefMult");
    xgraph->Draw("ap");
    
    canvas->cd(6);
    ygraph->SetMarkerStyle(21);
    ygraph->SetMarkerColor(kRed);
    ygraph->SetMarkerSize(0.5);
    ygraph->SetTitle("Linear Weights Quantiles");
    ygraph->GetXaxis()->SetTitle("Quantile");
    ygraph->GetYaxis()->SetTitle("Linear Weight Multiplicity");
    ygraph->Draw("ap");

    canvas->cd(2); // too many
    gPad->SetLogz();
    quantileHistogram->Draw("Colz");

    canvas->cd(7);
    quantileXProject->SetLineColor(kRed);
    quantileXProject->SetMarkerColor(kRed);
    quantileXProject->SetMarkerStyle(kOpenTriangleUp);
    quantileXProject->SetMarkerSize(0.5);
    quantileXProject->Draw("hist l p");
    
    canvas->cd(8);
    quantileYProject->SetLineColor(kBlue);
    quantileYProject->SetMarkerColor(kBlue);
    quantileYProject->SetMarkerStyle(kStar);
    quantileYProject->SetMarkerSize(0.5);
    quantileYProject->Draw("hist l p");
    }


    TH2D **quantileHistograms = (TH2D**)malloc(sizeof(TH2D*) * 4);
    TH1D **quantileXProjects = (TH1D**)malloc(sizeof(TH1D*) * 4);
    TH1D **quantileYProjects = (TH1D**)malloc(sizeof(TH1D*) * 4);

    quantileHistograms[0] = getQuantileRange(20, 25, histogram, xyQuantiles);
    quantileHistograms[1] = getQuantileRange(45, 50, histogram, xyQuantiles);
    quantileHistograms[2] = getQuantileRange(70, 75, histogram, xyQuantiles);
    quantileHistograms[3] = getQuantileRange(95, 100, histogram, xyQuantiles);

    for (uint32_t i = 0; i < 4; i++) {
        quantileXProjects[i] = quantileHistograms[i]->ProjectionX();
        std::cout << quantileXProjects[i]->ComputeIntegral() << std::endl;
        quantileYProjects[i] = quantileHistograms[i]->ProjectionY();
        std::cout << quantileYProjects[i]->ComputeIntegral() << std::endl << std::endl;
    }


    TCanvas *canvas2 = new TCanvas("canvas2", "Quantile Ranges", 1000, 1000);
    canvas2->Divide(2, 2);
    for (uint32_t i = 0; i < 4; i++) {
        canvas2->cd(i + 1);
        quantileXProjects[i]->SetLineColor(kRed);
        quantileXProjects[i]->SetMarkerColor(kRed);
        quantileXProjects[i]->SetMarkerStyle(kOpenTriangleUp);
        quantileXProjects[i]->SetMarkerSize(0.5);
        quantileXProjects[i]->Draw("hist l p");
        
        quantileYProjects[i]->SetLineColor(kBlue);
        quantileYProjects[i]->SetMarkerColor(kBlue);
        quantileYProjects[i]->SetMarkerStyle(kStar);
        quantileYProjects[i]->SetMarkerSize(0.5);
        quantileYProjects[i]->Draw("same hist l p");

        TLegend *legend = new TLegend(0.65, 0.78, 0.75, 0.85);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.03);
        legend->AddEntry(quantileXProjects[i]->GetName(), "RefMult", "l");
        legend->AddEntry(quantileYProjects[i]->GetName(), "Linear Weights", "l");
        legend->Draw();

    }
    canvas2->Draw();

}

void quantiles(const char *inHistName="data/hist_data.root") {
    TFile inHist(inHistName);
    TH2D *histogram;
    inHist.GetObject("histogram", histogram);
    histogram->SetDirectory(0); // what is this black magic
    inHist.Close();
    quantileAnalysis(histogram);
    // TH2D *testingHistogram = new TH2D("test", "test",
    //                                     100, 0, 100*100,
    //                                     100, 0, 100);
    // for (double i = 0; i < 100; i+=0.01) {
    //     testingHistogram->Fill(i*i, i);
    // }
    // quantileAnalysis(testingHistogram);
}