/**
 * \brief Generates the vector W giving the ridge 
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
float alpha = 1e8;

// Takes the data matrix c and global vector g and generates 
// weights relating the two using ridge regression
TMatrixD* generateWeights(const TMatrixD *c, const TVectorD * g) {
    TMatrixD* data = new TMatrixD(c->GetNrows() + 1, c->GetNcols());
    TMatrixD* expected = new TMatrixD(dim, 1);
    TMatrixD* identity = new TMatrixD(dim, dim);
    int32_t numEvents = c->GetNcols();


    // Add row of ones to the data matrix
    for (uint32_t i = 0; i < c->GetNrows(); i++) {
        for (uint32_t j = 0; j < c->GetNcols(); j++) {
            (*data)[i+1][j] = (*c)[i][j];
        }
    }
    for (uint32_t i = 0; i < c->GetNrows(); i++) {
        (*data)[i][0] = 1;
    }

    // Make the identity matrix
    for (uint32_t q = 0; q < dim; q++) {
        for (uint32_t t = 0; t < dim; t++) {
            (*identity)[q][t] = 0;
        }
        (*identity)[q][q] = alpha;
    }

    // Make the true values vector
    for (uint32_t t = 0; t < dim - 1; t++) {
        (*expected)[t+1][0] = 0;
        for (uint32_t j = 0; j < numEvents; j++) {
            (*expected)[t+1][0] += (*g)[j] * (*c)[t][j];
        }
    }
    (*expected)[0][0] = numEvents;

    TMatrixD *data_t = new TMatrixD(TMatrixD::kTransposed, *data);
    TMatrixD *first = new TMatrixD(dim, dim);
    first->Mult(*data, *data_t);
    first->Print();
    *first += *identity;
    first->Invert();
    *first *= *data;
    TMatrixD *weights = new TMatrixD(dim, 1);
    weights->Mult(*first, *expected);
    
    delete expected;
    delete identity;
    delete data_t;
    delete first;

    return weights;
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

void ridgeRegression(const char *inFileName = "data/simRingSums.root") {
    std::cout << "Running..." <<std::endl;
    
    TFile inFile(inFileName);
    TMatrixD *c;
    TVectorD *g;
    inFile.GetObject("ring_sums", c);           // should probably check for existance but...
    inFile.GetObject("tpc_multiplicity", g);

    std::cout << "Generating Weights.." << std::endl;
    TMatrixD *weights = generateWeights(c, g);
    std::cout << "Weights Generated" << std::endl;
    std::cout << "Weights:\n";
    weights->Print();

    std::cout << "Applying linear weights..." << std::endl;
    TVectorD *predictions = predictTPCMultiplicity(weights, c);
    inFile.Close();


    uint32_t bins = 162;
    int32_t predictMin = -500;
    int32_t predictMax = -1400;
    int32_t realMin = -26;
    int32_t realMax = 300; 
    
    alpha = 0;
    weights = generateWeights(c, g);
    predictions = predictTPCMultiplicity(weights, c);
    TH2D *a_0 = new TH2D("alpha=0", "alpha=0;TPC RefMult;Ridge Regression Prediction",
                                  bins, realMin, realMax,
                                  bins, predictMin, predictMax);
    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        a_0->Fill((*g)[i], (*predictions)[i]);
    }


    alpha = -1e5;
    weights = generateWeights(c, g);
    predictions = predictTPCMultiplicity(weights, c);
    TH2D *a_1e5 = new TH2D("alpha=1e5", "alpha=1e5;TPC RefMult;Ridge Regression Prediction",
                                  bins, realMin, realMax,
                                  bins, predictMin, predictMax);
    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        a_1e5->Fill((*g)[i], (*predictions)[i]);
    }


    alpha = -1e6;
    weights = generateWeights(c, g);
    predictions = predictTPCMultiplicity(weights, c);
    TH2D *a_1e6 = new TH2D("alpha=1e6", "alpha=1e6;TPC RefMult;Ridge Regression Prediction",
                                  bins, realMin, realMax,
                                  bins, predictMin, predictMax);
    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        a_1e6->Fill((*g)[i], (*predictions)[i]);
    }


    alpha = -1e7;
    weights = generateWeights(c, g);
    predictions = predictTPCMultiplicity(weights, c);
    TH2D *a_1e7 = new TH2D("alpha=1e7", "alpha=1e7;TPC RefMult;Ridge Regression Prediction",
                                  bins, realMin, realMax,
                                  bins, predictMin, predictMax);
    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        a_1e7->Fill((*g)[i], (*predictions)[i]);
    }
    

    // Everything from here down is plotting
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(11);
    gStyle->SetStatX(0.38);
    gStyle->SetStatY(0.85);




    TCanvas *canvas = new TCanvas("canvas", "canvas", 1900, 1000);
    canvas->Divide(2, 2);

    canvas->cd(1);
    gPad->SetLogz();
    a_0->Draw("Colz");

    canvas->cd(2);
    gPad->SetLogz();
    a_1e5->Draw("Colz");

    canvas->cd(3);
    gPad->SetLogz();
    a_1e6->Draw("Colz");

    canvas->cd(4);
    gPad->SetLogz();
    a_1e7->Draw("Colz");

    std::cout << "Plotted " << g->GetNrows() << " events\n";
}