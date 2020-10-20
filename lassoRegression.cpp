/**
 * \brief Generates the vector W using lasso 
 * regression to correlate the EPD nMIP data to the TPC multiplicity
 * 
 * \author Tristan Protzman
 * \date September 30, 2020
 * \email tlprotzman@gmail.com
 * \affiliation Lehigh University
 * 
 */

// #define DEBUG

#include <iostream>

// Root headers
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVectorD.h"

const uint32_t dim = 17;

// Takes the data matrix c and global vector g and generates 
// weights relating the two using ridge regression
TMatrixD* generateWeights(const TMatrixD *c, const TVectorD * g, float alpha) {
    // c is the data matrix
    // g is the response vector
    TMatrixD *data = new TMatrixD(TMatrixD::kTransposed, *c);
    TMatrixD *ones = new TMatrixD(data->GetNcols(), 1);
    TMatrixD *mean = new TMatrixD(data->GetNrows(), 1);
    
    // Calculate the average of the data matrix;
    for (uint32_t i = 0; i < ones->GetNrows(); i++) {
        (*ones)[i][0] = 1;
    }
    mean->Mult(*data, *ones);
    TMatrixD *centeringData = new TMatrixD(*mean, TMatrixD::kMult, TMatrixD(TMatrixD::kTransposed, *ones));
    (*data) -= (*centeringData);    // Centering the data matrix

    // Calculating the mean of the response vector
    double avgResponse = 0;
    for (uint32_t i = 0; i < g->GetNoElements(); i++) {
        avgResponse += (*g)[i];
    }
    avgResponse /= g->GetNoElements();
    std::cout << "Average Response: " << avgResponse << std::endl;
    TMatrixD *response = new TMatrixD(g->GetNoElements(), 1);
    for (uint32_t i = 0; i < g->GetNoElements(); i++) {
        (*response)[i][0] = (*g)[i] - avgResponse;
    }

    // Set up initial random weights vector
    TMatrixD *weights = new TMatrixD(data->GetNcols(), 1);
    for (uint32_t i = 0; i < weights->GetNrows(), i++) {
        (*weights)[i][0] = 0.1 * i; // MAKE RANDOM
    }

    // main weights loop
    double weightDelta = 100000;
    double epsilon = 0.5;
    while 

    
    return nullptr;
}


// Uses the generated weights and the truncated nMIPs data to predict TPC multiplicity
// X_t = sum_r W_r * C_{r, t}
TVectorD* predictTPCMultiplicity(TMatrixD *weights, TMatrixD *epdData) {
    uint32_t numEvents = epdData->GetNcols();
    TVectorD *predictedTCPMultiplicity = new TVectorD(numEvents);   // Store our guesses
    for (uint32_t i = 0; i < numEvents; i++) {
        if (i%100 == 0) {
            // std::cout << "Working on event " << i << "/" << numEvents << std::endl;
        }
        (*predictedTCPMultiplicity)[i] = 0;//(*weights)[0][0]; // Add the bias
        for (uint32_t j = 0; j < dim - 1; j++) {
            (*predictedTCPMultiplicity)[i] += (*epdData)[j][i] * (*weights)[j + 1][0]; // Add the weighted input
        }
    }
    // predictedTCPMultiplicity->Print();
    return predictedTCPMultiplicity;
}

void lassoRegression(const char *inFileName = "data/detector_data.root", float alpha=1e7) {
    std::cout << "Running..." <<std::endl;
    
    TFile inFile(inFileName);
    TMatrixD *c;
    TVectorD *g;
    inFile.GetObject("ring_sums", c);           // should probably check for existance but...
    inFile.GetObject("tpc_multiplicity", g);

    std::cout << "Generating Weights.." << std::endl;
    TMatrixD *weights = generateWeights(c, g, alpha);
    return;
    std::cout << "Weights Generated" << std::endl;
    std::cout << "Weights:\n";
    weights->Print();

    std::cout << "Applying linear weights..." << std::endl;
    TVectorD *predictions = predictTPCMultiplicity(weights, c);
    inFile.Close();


    uint32_t predictBins = 200;
    int32_t predictMin = -100;
    int32_t predictMax = 300;
    
    uint32_t realBins = 175;
    int32_t realMin = 0;
    int32_t realMax = 350; 
    
    TH2D *ridge_histogram = new TH2D(Form("alpha=%f", alpha), Form("alpha=%f;TPC RefMult;Ridge Regression Prediction", alpha),
                                  realBins, realMin, realMax,
                                  predictBins, predictMin, predictMax);
    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        ridge_histogram->Fill((*g)[i], (*predictions)[i]);
    }
    
    TFile outFile("data/epd_tpc_relations.root", "UPDATE");
    outFile.mkdir("methods", "methods", true);
    outFile.cd("methods");
    ridge_histogram->Write(Form("ridge_%f_detector", alpha));
    outFile.Close();

    bool draw = false;
    if (!draw) {
        return;
    }

    // Everything from here down is plotting
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(11);
    gStyle->SetStatX(0.38);
    gStyle->SetStatY(0.85);


    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    gPad->SetLogz();
    ridge_histogram->Draw("Colz");

    std::cout << "Plotted " << g->GetNrows() << " events\n";
}