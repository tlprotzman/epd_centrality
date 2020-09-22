/**
 * \brief Generates the vector W giving the linear 
 * weights to correlate the EPD nMIP data to the TPC multiplicity
 * 
 * \author Tristan Protzman
 * \date September 17, 2020
 * \email tlprotzman@gmail.com
 * 
 */

#define DEBUG

#include <iostream>

// Root headers
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVectorD.h"

const uint32_t dim = 17;

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
    for (uint32_t q = 0; q < dim - 1; q++) {
        for (uint32_t t = 0; t < dim - 1; t++) {
            for (uint32_t j = 0; j < numEvents; j++) {
                (*a)[q][t] += (*c)[q][j] * (*c)[t][j];
            }
        }
    }

    // Step 2
    for (uint32_t t = 0; t < dim - 1; t++) {
        for (uint32_t j = 0; j < numEvents; j++) {
            (*a)[dim - 1][t] += (*c)[t][j];
            (*a)[t][dim - 1] += (*c)[t][j];
        }
    }

    // Step 3
    (*a)[dim -1][dim -1] = numEvents;

    // Step 4
    for (uint32_t t = 0; t < dim - 1; t++) {
        for (uint32_t j = 0; j < numEvents; j++) {
            (*b)[t][0] += (*g)[j] * (*c)[t][j];
        }
    }

    // Step 5
    for (uint32_t j = 0; j < numEvents; j++) {
        (*b)[dim - 1][0] += (*g)[j];
    }

    // std::cerr << "A and B generated\n";
    
    // Debug printing
    #ifdef DEBUG
    std::cout << "A:\n";
    a->Print();
    std::cout << "\n\nB:\n";
    b->Print();
    #endif // DEBUG

    // Inverting A to solve Ax=B for x
    a->Invert();
    TMatrixD *weights = new TMatrixD(17, 1);
    weights->Mult(*a, *b);
    

    #ifdef DEBUG
    std::cout << "A Inverted:\n";
    a->Print();
    #endif // DEBUG

    
    return weights;
}

TVectorD* predictTPCMultiplicity(TMatrixD *weights, TMatrixD *epdData) {
    uint32_t numEvents = epdData->GetNcols();
    TVectorD *predictedTCPMultiplicity = new TVectorD(numEvents);
    for (uint32_t i = 0; i < numEvents; i++) {
        (*predictedTCPMultiplicity)[i] = (*weights)[16][0];
        for (uint32_t j = 0; j < 16; j++) {
            (*predictedTCPMultiplicity) += (*epdData)[j][i] * (*weights)[j][0];
        }
    }
    // predictedTCPMultiplicity->Print();
    return predictedTCPMultiplicity;
}


void linearWeights() {
    std::cout << "Running..." <<std::endl;
    
    TFile inFile("ringSums.root");
    TMatrixD *c;
    TVectorD *g;
    inFile.GetObject("ring_sums", c);
    inFile.GetObject("tpc_multiplicity", g);
    TMatrixD *weights = generateWeights(c, g);
    std::cout << "Weights:\n";
    weights->Print();
    TVectorD *predictions = predictTPCMultiplicity(weights, c);
    inFile.Close();

    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    uint32_t bins = 100;
    int32_t predictMin = 0;
    int32_t predictMax = 35000;
    int32_t realMin = 0;
    int32_t realMax = 200; 

    TH2D *predictVsReal = new TH2D("predictVsReal", "2D Histo;Real;Predicted",
                                  bins, realMin, realMax,
                                  bins, predictMin, predictMax);


    for (uint32_t i = 0; i < g->GetNrows(); i++) {
        predictVsReal->Fill((*g)[i], (*predictions)[i]);
    }

    TCanvas *canvas = new TCanvas("Canvas", "Canvas", 800, 800);
    predictVsReal->Draw("Colz");
}