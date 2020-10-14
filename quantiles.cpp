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
TH2D **getQuantileRange(int lowerQuantile, int upperQuantile, TH2D *histogram, double *quantileX, double* quantileY, const char *method) {
    TH2D **quantileHistogram = (TH2D**)malloc(2 * sizeof(TH2D*));
    quantileHistogram[0] = new TH2D(Form("quantilesX%d_%d", lowerQuantile, upperQuantile), Form("%s", method),
            histogram->GetNbinsX(), histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax(),
            histogram->GetNbinsY(), histogram->GetYaxis()->GetXmin(), histogram->GetYaxis()->GetXmax());

    quantileHistogram[1] = new TH2D(Form("quantilesY%d_%d", lowerQuantile, upperQuantile), Form("%s", method),
            histogram->GetNbinsX(), histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax(),
            histogram->GetNbinsY(), histogram->GetYaxis()->GetXmin(), histogram->GetYaxis()->GetXmax());


    quantileHistogram[0]->SetXTitle("Multiplicity");
    quantileHistogram[0]->SetYTitle("Counts");
 
    quantileHistogram[1]->SetXTitle("Multiplicity");
    quantileHistogram[1]->SetYTitle("Counts");

    double minX, maxX;
    double minY, maxY;
    lowerQuantile -= 1;      // Since the quantile array gives the upper bound of the quantile, the lower bound should be the one below
    if (lowerQuantile < 0) { // the requested quantile
        minX = histogram->GetXaxis()->GetXmin(); // If we want from 0, set lower bound to be the minimum
        minY = histogram->GetYaxis()->GetXmin(); // If we want from 0, set lower bound to be the minimum
    }
    else {
        minX = quantileX[lowerQuantile];
        minY = quantileY[lowerQuantile];
    }
    if (upperQuantile > numberQuantiles - 1) {
        maxX = histogram->GetXaxis()->GetXmax();
        maxY = histogram->GetYaxis()->GetXmax();
    }
    else{
        maxX = quantileX[upperQuantile - 1];
        maxY = quantileY[upperQuantile - 1];
    }
    double *centersX = (double*)malloc(sizeof(double) * histogram->GetNbinsX());
    double *centersY = (double*)malloc(sizeof(double) * histogram->GetNbinsY());
    histogram->GetXaxis()->GetCenter(centersX);
    histogram->GetYaxis()->GetCenter(centersY);
    std::cout << "Selecting X events between " << minX << "  and " << maxX << std::endl;
    for (uint32_t i = 0; i < histogram->GetNbinsX(); i++) {
        if (minX < centersX[i] && centersX[i] < maxX) {
            for (uint32_t j = 0; j < histogram->GetNbinsY(); j++) {
                for (uint32_t k = 0; k < histogram->GetBinContent(i, j); k++) {
                    quantileHistogram[0]->Fill(centersX[i], centersY[j]);
                }
            }
        }
    }
    
    std::cout << "Selecting Y events between " << minY << "  and " << maxY << std::endl;
    for (uint32_t i = 0; i < histogram->GetNbinsX(); i++) {
        for (uint32_t j = 0; j < histogram->GetNbinsY(); j++) {
            if (minY < centersY[j] && centersY[j] < maxY) {
                for (uint32_t k = 0; k < histogram->GetBinContent(i, j); k++) {
                    quantileHistogram[1]->Fill(centersX[i], centersY[j]);
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
        yyQuantiles[i] = 0;
    }
    TH1D *yProjection = histogram->ProjectionY();
    xProjection->GetQuantiles(numberQuantiles, xyQuantiles, xxQuantiles);
    yProjection->GetQuantiles(numberQuantiles, yyQuantiles, yxQuantiles);

    // for (uint32_t i = 0; i < numberQuantiles; i++) {
    //     std::cout << "Quantile " << i << ": " << xyQuantiles[i] << std::endl;
    // }

    TH2D **quantileHistogram = getQuantileRange(95, 100, histogram, xyQuantiles, yyQuantiles, "");
    TH1D *quantileXProject = quantileHistogram[0]->ProjectionX();
    TH1D *quantileYProject = quantileHistogram[1]->ProjectionX();

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

    canvas->cd(7); // too many
    gPad->SetLogz();
    quantileHistogram[0]->Draw("Colz");

    canvas->cd(8); // too many
    gPad->SetLogz();
    quantileHistogram[1]->Draw("Colz");

    canvas->cd(2);
    quantileXProject->SetLineColor(kRed);
    quantileXProject->SetMarkerColor(kRed);
    quantileXProject->SetMarkerStyle(kOpenTriangleUp);
    quantileXProject->SetMarkerSize(0.5);
    quantileXProject->Draw("hist l p");
    
    // canvas->cd(8);
    quantileYProject->SetLineColor(kBlue);
    quantileYProject->SetMarkerColor(kBlue);
    quantileYProject->SetMarkerStyle(kStar);
    quantileYProject->SetMarkerSize(0.5);
    quantileYProject->Draw("same hist l p");
    }


    TH2D ***quantileHistograms = (TH2D***)malloc(sizeof(TH2D**) * 4);
    TH1D **quantileXProjects = (TH1D**)malloc(sizeof(TH1D*) * 4);
    TH1D **quantileYProjects = (TH1D**)malloc(sizeof(TH1D*) * 4);

    const char *method = "Ridge, -1e5";
    quantileHistograms[0] = getQuantileRange(0, 10, histogram, xyQuantiles, yyQuantiles, method);
    quantileHistograms[1] = getQuantileRange(70, 80, histogram, xyQuantiles, yyQuantiles, method);
    quantileHistograms[2] = getQuantileRange(95, 100, histogram, xyQuantiles, yyQuantiles, method);
    quantileHistograms[3] = getQuantileRange(94, 100, histogram, xyQuantiles, yyQuantiles, method);

    for (uint32_t i = 0; i < 4; i++) {
        quantileXProjects[i] = quantileHistograms[i][0]->ProjectionX();
        // std::cout << quantileXProjects[i]->ComputeIntegral() << std::endl;
        quantileYProjects[i] = quantileHistograms[i][1]->ProjectionX();
        // std::cout << quantileYProjects[i]->ComputeIntegral() << std::endl << std::endl;
    }


    TCanvas *canvas2 = new TCanvas("canvas2", "Quantile Ranges", 1000, 1000);
    // canvas2->Divide(2, 2);
    char *formatA = "hist l p";
    char *formatB = "same hist l p";
    char *format = formatA;
    for (uint32_t i = 0; i < 3; i++) {
        // canvas2->cd(i + 1);
        gPad->SetLogy();
        quantileXProjects[i]->SetLineColor(kRed);
        quantileXProjects[i]->SetMarkerColor(kRed);
        quantileXProjects[i]->SetMarkerStyle(kOpenTriangleUp);
        quantileXProjects[i]->SetMarkerSize(0.5);
        quantileXProjects[i]->Draw(format);
        format = formatB;

        quantileYProjects[i]->SetLineColor(kBlue);
        quantileYProjects[i]->SetMarkerColor(kBlue);
        quantileYProjects[i]->SetMarkerStyle(kStar);
        quantileYProjects[i]->SetMarkerSize(0.5);
        quantileYProjects[i]->Draw(format);

        TLegend *legend = new TLegend(0.65, 0.78, 0.75, 0.85);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextSize(0.03);
        legend->AddEntry(quantileXProjects[i]->GetName(), "RefMult", "l");
        legend->AddEntry(quantileYProjects[i]->GetName(), "Linear Weights", "l");
        legend->Draw();

    }
    canvas2->Draw();
    TFile output("data/quantiles.root", "UPDATE");
    canvas2->Write("ridge_-1e5");
    output.Close();

}

void quantiles(const char *inHistName="data/linear_hist_data.root") {
    TFile inHist(inHistName);
    TH2D *histogram;
    inHist.GetObject("histogram", histogram);
    histogram->SetDirectory(0); // what is this black magic
    inHist.Close();
    quantileAnalysis(histogram);
}