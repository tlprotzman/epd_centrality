/**
 * \brief Converts the data from the simulation format to what is
 * expected in the linear weighting program.  Probably should be 
 * phased out shortly once I understand ROOT's tuples better.
 * 
 * \author Tristan Protzman
 * \email tlprotzman@gmail.com
 * \date September 23, 2020
 */

#include <iostream>



// ROOT Headers
#include "TROOT.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TLeafF.h"
#include "TTree.h"
#include "TBranch.h"
#include "TNtuple.h"
#include "TTreeReader.h"

const uint8_t RINGS = 16;

void simulationDataPreprocessor(const char *inFileName = "CentralityNtupleout06212020_7.7.root" ) {
    std::cout << "Converting data format..." << std::endl;
    TFile *inFile = TFile::Open(inFileName);

    TTreeReader eventReader("Rings", inFile);
    Long64_t numEvents = eventReader.GetEntries();
    std::cout << "Processing " << numEvents << " events" << std::endl;

    TMatrixD *ringSums = new TMatrixD(RINGS, numEvents);
    TVectorD *tpcMultiplicity = new TVectorD(numEvents);

    TTreeReaderValue<Float_t> *rings[RINGS];
    for (uint8_t i = 1; i <= RINGS; i++) {
        rings[i] = new TTreeReaderValue<Float_t>(eventReader, Form("r%02d", i));
        // std::cout << Form("r%02d", i) << std::endl;;
    }
    TTreeReaderValue<Float_t> refMul(eventReader, "RefMult1");

    Long64_t i = 0;
    while(eventReader.Next()) {
        for (uint8_t j = 0; j < RINGS; j++) {
            (*ringSums)[j][i] = 5;
        }
        (*tpcMultiplicity)[i] = *refMul;
        i++;
    }
    // ringSums->Print();
    tpcMultiplicity->Print();
}