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
#include <vector>



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

    std::vector<TTreeReaderValue<Float_t>> ringReaders;
    for (uint32_t i = 1; i <= RINGS; i++) {
        TTreeReaderValue<Float_t> reader(eventReader, Form("r%02d", i));
        ringReaders.push_back(reader);
    }
    TTreeReaderValue<Float_t> refMul(eventReader, "RefMult1");

    Long64_t i = 0;
    while(eventReader.Next()) {
        for (uint32_t j = 0; j < RINGS; j++) {
            (*ringSums)[j][i] = *(ringReaders[j]);
        }
        (*tpcMultiplicity)[i] = *refMul;
        i++;
    }
    // ringSums->Print();
    // tpcMultiplicity->Print();

    TFile outFile("simRingSums.root", "RECREATE");
    ringSums->Write("ring_sums");
    tpcMultiplicity->Write("tpc_multiplicity");
    outFile.Close();

}