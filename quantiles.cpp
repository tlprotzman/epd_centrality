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

#include <stdint.h>

#include "quantiles.h"

const int32_t numberQuantiles = 100;

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
    // for (uint32_t i = 0; i < numberQuantiles; i++){
    //     std::cout << xQuantiles[i] << "\t" << yQuantiles[i] << std::endl;
    // }
    TGraph *xgraph = new TGraph(numberQuantiles, xxQuantiles, xyQuantiles);
    TGraph *ygraph = new TGraph(numberQuantiles, yxQuantiles, yyQuantiles);
    TCanvas *canvas2 = new TCanvas("canvas2", "testing");
    TGraphQQ *qq = new TGraphQQ(numberQuantiles, yyQuantiles, numberQuantiles, xyQuantiles);
    canvas2->Divide(2, 3);

    canvas2->cd(1);
    gPad->SetLogz();
    histogram->Draw("Colz");
    
    canvas2->cd(2);
    gPad->SetGrid();
    qq->SetTitle("Quantile-Quantile Graph");
    qq->SetMarkerStyle(21);
    qq->SetMarkerColor(kRed);
    qq->SetMarkerSize(0.5);
    qq->Draw("ap");
    
    canvas2->cd(3);
    xProjection->SetTitle("TPC Projection");
    xProjection->SetXTitle("TPC RefMult");
    xProjection->SetYTitle("Counts");
    xProjection->Draw();
    
    canvas2->cd(4);
    yProjection->SetTitle("Linear Weights Projection");
    yProjection->SetXTitle("Linear Weight Multiplicity");
    yProjection->SetYTitle("Counts");
    yProjection->Draw();
    
    canvas2->cd(5);
    xgraph->SetMarkerStyle(21);
    xgraph->SetMarkerColor(kRed);
    xgraph->SetMarkerSize(0.5);
    xgraph->SetTitle("TPC Quantiles");
    xgraph->Draw("ap");
    
    canvas2->cd(6);
    ygraph->SetMarkerStyle(21);
    ygraph->SetMarkerColor(kRed);
    ygraph->SetMarkerSize(0.5);
    ygraph->SetTitle("Linear Weights Quantiles");
    ygraph->Draw("ap");

}