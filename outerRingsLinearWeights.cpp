/**
 * \brief Generates the vector W giving the linear 
 * weights to correlate the EPD nMIP data to the TPC multiplicity
 * using only the outer rings
 * 
 * \author Tristan Protzman
 * \date November 17, 2020
 * \email tlprotzman@gmail.com
 * 
 */

// #define DEBUG


#include <iostream>

// Root headers
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphQQ.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVectorD.h"
#include "TMatrixDUtils.h"
#include "TVectorDfwd.h"


const uint32_t real_dim = 17;
const uint32_t inner_ring = 7;
const uint32_t dim = real_dim - inner_ring;

// (Step 1) A_{q, t}   = \sum_j=1^Nevents C_{q, j} C{t, j}
// (Step 2) A_{17, t}  = \sum_j=1^Nevents C_{t, j}
// (Step 3) A_{17, 17} = Nevents

// (Step 4) B_t  = \sum_j=1^Nevents G_j C{t, j}
// (Step 5) B_17 = \sum_j=1^Nevents G_j


// Takes the matrix C and the vector G and generates the weight vector W
TMatrixD* generateWeights (const TMatrixD *c, const TVectorD *g) {
    TMatrixD *a = new TMatrixD(dim, dim);    // Creates a double precision matrix
    TMatrixD *b = new TMatrixD(dim, 1);
    int32_t numEvents = c->GetNcols();

    std::cerr << "Processings " << numEvents << " events.\n";

    // Generating A
    // Step 1
    std::cout << "Generating A..." << std::endl;
    for (uint32_t q = 0; q < dim - 1; q++) {
        for (uint32_t t = 0; t < dim - 1; t++) {
            (*a)[q][t] = 0; // likely not needed
            for (uint32_t j = 0; j < numEvents; j++) {
                (*a)[q][t] += (*c)[q + inner_ring][j] * (*c)[t + inner_ring][j];
            }
        }
    }

    // Step 2
    for (uint32_t t = 0; t < dim - 1; t++) {
        for (uint32_t j = 0; j < numEvents; j++) {
            (*a)[dim - 1][t] += (*c)[t+inner_ring][j];
            (*a)[t][dim - 1] += (*c)[t+inner_ring][j];
        }
    }

    // Step 3
    (*a)[dim -1][dim -1] = numEvents;

    std::cout << "Generated A" << std::endl;

    // Step 4
    for (uint32_t t = 0; t < dim - 1; t++) {
        (*b)[t][0] = 0;
        for (uint32_t j = 0; j < numEvents; j++) {
            (*b)[t][0] += (*g)[j] * (*c)[t+inner_ring][j];
        }
    }

    // Step 5
    for (uint32_t j = 0; j < numEvents; j++) {
        (*b)[dim - 1][0] += (*g)[j];
    }

    std::cout << "Generated B" << std::endl;
    
    // Debug printing
    #ifdef DEBUG
    std::cout << "A:\n";
    a->Print();
    std::cout << "\n\nB:\n";
    b->Print();
    #endif // DEBUG

    // std::cout << "\n\n\nA:\n";
    // for (uint32_t i = 0; i < 17; i++) {
    //     for (uint32_t j = 0; j < 17; j++) {
    //         std::cout << (*a)[i][j] << ", ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << (*b)[i][0] << ", ";
    // }
    // std::cout << "\n\n\n";

    // Inverting A to solve Ax=B for x
    std::cout << "Inverting A" << std::endl;
    a->Invert();
    TMatrixD *weights = new TMatrixD(17, 1);
    weights->Mult(*a, *b);
    

    #ifdef DEBUG
    std::cout << "A Inverted:\n";
    a->Print();
    #endif // DEBUG

    weights->Print();
    for (uint32_t i = real_dim - 1; i >= inner_ring; i--) {
        (*weights)[i][0] = (*weights)[i - inner_ring][0]; 
    }
    for(uint32_t i = 0; i < inner_ring; i++) {
        (*weights)[i][0] = 0;
    }
    return weights;
}

// Uses the generated weights and the truncated nMIPs data to predict TPC multiplicity
// X_t = sum_r W_r * C_{r, t} + W_17
TVectorD* predictTPCMultiplicity(TMatrixD *weights, TMatrixD *epdData) {
    uint32_t numEvents = epdData->GetNcols();
    TVectorD *predictedTCPMultiplicity = new TVectorD(numEvents);   // Store our guesses
    for (uint32_t i = 0; i < numEvents; i++) {
        // if (i%100 == 0) {
        //     std::cout << "Working on event " << i << "/" << numEvents << std::endl;
        // }
        (*predictedTCPMultiplicity)[i] = (*weights)[real_dim - 1][0]; // Add the bias
        for (uint32_t j = 0; j < real_dim - 1; j++) {
            (*predictedTCPMultiplicity)[i] += (*epdData)[j][i] * (*weights)[j][0]; // Add the weighted input
        }
    }
    // predictedTCPMultiplicity->Print();
    return predictedTCPMultiplicity;
}

void outerRingsLinearWeights(const char *inFileName = "data/detector_data.root") {
    std::cout << "Running..." <<std::endl;
    
    TFile inFile(inFileName);
    TMatrixD *c;
    TVectorD *g;
    inFile.GetObject("ring_sums", c);           // should probably check for existance but...
    inFile.GetObject("tpc_multiplicity", g);
    TMatrixD *detector_sums;
    TVectorD *detector_refmult;
    TFile detector("data/detector_data.root");
    detector.GetObject("ring_sums", detector_sums);


    std::cout << "Generating Weights.." << std::endl;
    TMatrixD *weights = generateWeights(c, g);
    std::cout << "Weights Generated" << std::endl;
    std::cout << "Weights:\n";
    weights->Print();

    std::cout << "Applying linear weights..." << std::endl;
    TVectorD *predictions = predictTPCMultiplicity(weights, detector_sums);
    inFile.Close();
    detector.GetObject("tpc_multiplicity", g);
    detector.Close();
    

    // Everything from here down is plotting

    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetStatX(0.38);
    gStyle->SetStatY(0.85);

    uint32_t predictBins = 200;
    int32_t predictMin = -100;
    int32_t predictMax = 300;
    
    uint32_t realBins = 175;
    int32_t realMin = 0;
    int32_t realMax = 350;

    TH2D *predictVsReal = new TH2D("linear_simulated", "2D Histo;RefMult1;X_{#zeta'}",
                                  realBins, realMin, realMax,
                                  predictBins, predictMin, predictMax);
    predictVsReal->SetTitle("X_{#zeta'} vs RefMult1, 7.7 GeV, TOF Selected, Outer 9 Rings");


    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        predictVsReal->Fill((*g)[i], (*predictions)[i]);
    }


    bool draw = true;
    if (draw) {
        TCanvas *canvas = new TCanvas("canvas", "canvas", 700, 500);
        gPad->SetLogz();
        predictVsReal->Draw("Colz");
    }

    std::cout << "Plotted " << g->GetNrows() << " events\n";

    TFile outFile("data/epd_tpc_relations.root", "UPDATE");
    outFile.mkdir("methods", "methods", true);
    outFile.cd("methods");
    weights->Write("linear_weights_outer");
    predictVsReal->Write("linear");
    outFile.Close();
}