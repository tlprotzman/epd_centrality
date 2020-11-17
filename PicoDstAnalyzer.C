/**
 * \brief Example of how to read a file (list of files) using StPicoEvent classes
 *
 * RunPicoDstAnalyzer.C is an example of reading STAR picoDst format.
 * One can use either picoDst file or a list of picoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * When using ROOT5, one needs to run RunAnalyzer.C macro when run processing.
 * This will handle libraries loading, etc.
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 * \email nigmatkulov@gmail.com
 */

// This is needed for calling standalone classes (not needed on RACF)
#define _VANILLA_ROOT_

// C++ headers
#include <iostream>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"

// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) 
R__LOAD_LIBRARY(StRoot/StPicoEvent/libStPicoDst)
#endif

// inFile - is a name of name.picoDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.picoDst.root files

//_________________
void PicoDstAnalyzer(const Char_t *inFile = "data/files.list") {
    
    std::cout << "Hi! Lets do some physics, Master!" << std::endl;
    
    StPicoDstReader* picoReader = new StPicoDstReader(inFile);
    picoReader->Init();
    
    
    // This is a way if you want to spead up IO
    std::cout << "Explicit read status for some branches" << std::endl;
    picoReader->SetStatus("*",0);
    picoReader->SetStatus("Event",1);
    picoReader->SetStatus("EpdHit",1);
    std::cout << "Status has been set" << std::endl;
    
    if( !picoReader->chain() ) {
        std::cout << "No chain has been found." << std::endl;
    }
    Long64_t eventsInTree = picoReader->tree()->GetEntries();
    std::cout << "eventsInTree: "  << eventsInTree << std::endl;
    Long64_t events2read = picoReader->chain()->GetEntries();
    
    std::cout << "Number of events to read: " << events2read
    << std::endl;
    
    TString OutFileName = "data/out.root";
    TFile *file1 = TFile::Open(OutFileName.Data(),"RECREATE");
    
    // Histogramming
    // Event
    TH1F *hRefMult = new TH1F("hRefMult",
                              "Reference multiplicity;refMult",
                              500, -0.5, 499.5);
    TH2F *hVtxXvsY = new TH2F("hVtxXvsY",
                              "hVtxXvsY",
                              200,-10.,10.,200,-10.,10.);
    TH1F *hVtxZ = new TH1F("hVtxZ","hVtxZ",
                           140, -70., 70.);
    
    // EPD
    
    TH1* mNmipDists[2][12][31];
    TH1* mAdcDists[2][12][31];
    TH2* hRingvsRegMult[2][16];
    TH1F* hSizeEpd = new TH1F("hSizeEpd","",1000,0,1000);
    
    // 1D histograms for nMIP distributions
    for (int ew=0; ew<2; ew++){
        for (int pp=1; pp<13; pp++){
            for (int tt=1; tt<32; tt++){
                mNmipDists[ew][pp-1][tt-1] = new TH1D(Form("NmipEW%dPP%dTT%d",ew,pp,tt),Form("NmipEW%dPP%dTT%d",ew,pp,tt),1000,0,50);
                mAdcDists[ew][pp-1][tt-1] = new TH1D(Form("ADCEW%dPP%dTT%d",ew,pp,tt),Form("ADCEW%dPP%dTT%d",ew,pp,tt),4095,0,4095);
            }
        }
        for (int r = 0;r<16;r++){
            hRingvsRegMult[ew][r] = new TH2F(Form("hRingvsRegMultEW%iRing%i",ew,r+1),Form("hRingvsRegMultEW%iRing%i",ew,r+1),500,-0.5,499.5,500,0,500);
            hRingvsRegMult[ew][r]->GetXaxis()->SetTitle("refMult");
            hRingvsRegMult[ew][r]->GetYaxis()->SetTitle(Form("sum TrNmip ring %i",r+1));
            
        }
    }

    TMatrixD *sums = new TMatrixD(16, 4);
    TVectorD *tpcMultiplicity = new TVectorD(4);
    TVectorD *tofMultiplicity = new TVectorD(4);
    uint32_t numEvents = 0;
    
    
    // Loop over events
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
        
        if (iEvent % 1000 == 0) {
            std::cout << "Working on event #[" << (iEvent+1)
            << "/" << events2read << "]" << std::endl;
        }
        
        Bool_t readEvent = picoReader->readPicoEvent(iEvent);
        if( !readEvent ) {
            std::cout << "Something went wrong" << std::endl;
            break;
        }
        
        // Retrieve picoDst
        StPicoDst *dst = picoReader->picoDst();
        
        // Retrieve event information
        StPicoEvent *event = dst->event();
        if( !event ) {
            std::cout << "Something went wrong" << std::endl;
            break;
        }
        
        //Select good events
        if ( abs(event->primaryVertex().Z()) > 70)  // what is 70? collision happened within 70 cm on z axis
            continue;
        float Vr = sqrt(event->primaryVertex().X()*event->primaryVertex().X()+event->primaryVertex().Y()*event->primaryVertex().X());   // xy distance of collision
        if (Vr > 2) // within 2 cm of beamline
            continue;

        // Selection on tof vs tpc multiplicity
        // float tolerance = 0.8;
        // UShort_t tofMult = event->btofTrayMultiplicity();
        // Int_t tpcMult = event->refMult();
        // if (!(tofMult * (1 - tolerance) < 2 * tpcMult && tofMult * (1 + tolerance) > 2 * tpcMult)) {
        //     continue;
        // }

        
        //Fill eventwise distributions
        hRefMult->Fill( event->refMult() );
        
        TVector3 pVtx = event->primaryVertex();
        hVtxXvsY->Fill( event->primaryVertex().X(), event->primaryVertex().Y() );
        hVtxZ->Fill( event->primaryVertex().Z() );
        
        UInt_t Nepd = dst->numberOfEpdHits();
        hSizeEpd->Fill(Nepd);
        
        float ringsum[2][16];
        for (int i = 0;i<2;i++){
            for (int j = 0;j<16;j++){
                ringsum[i][j] = 0;
            }
        }
        for (UInt_t iepd = 0;iepd<Nepd;iepd++){
            StPicoEpdHit* epdhit = dst->epdHit(iepd);
            // epdhit->Print();
            
            int ew;
            if (epdhit->side() < 0)
                ew = 0;
            else
                ew = 1;
            float nMip;
            if (epdhit->nMIP() < 0.2)   // some kind of clamping between 0.2 and 2?
                nMip = 0;
            else if (epdhit->nMIP() > 3)
                nMip = 3;
            else nMip = epdhit->nMIP();
            ringsum[ew][(int)epdhit->tile()/2]+=nMip;
            mNmipDists[ew][epdhit->position()-1][epdhit->tile()-1]->Fill(nMip);
            mAdcDists[ew][epdhit->position()-1][epdhit->tile()-1]->Fill(epdhit->adc());
        } 
        
        for (int i = 0;i<2;i++){
            for (int j = 0;j<16;j++){
                hRingvsRegMult[i][j]->Fill(event->refMult(),ringsum[i][j]);
            }
        }

        // Saving events to file to generate weights
        for (uint32_t i = 0; i < 16; i++) {
            (*sums)[i][numEvents] = ringsum[0][i] + ringsum[1][i];
        }
        (*tpcMultiplicity)[numEvents] = event->refMult();
        (*tofMultiplicity)[numEvents] = event->btofTrayMultiplicity();
        numEvents++;
        if (numEvents >= sums->GetNcols()) {
            sums->ResizeTo(sums->GetNrows(), sums->GetNcols() * 2);
            tpcMultiplicity->ResizeTo(tpcMultiplicity->GetNrows() * 2);
            tofMultiplicity->ResizeTo(tofMultiplicity->GetNrows() * 2);
        }
        
        
    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
    
    file1->Write();
    file1->Close();
    
    
    picoReader->Finish();

    // Trim off empty space
    sums->ResizeTo(sums->GetNrows(), numEvents);
    tpcMultiplicity->ResizeTo(numEvents);
    tofMultiplicity->ResizeTo(numEvents);

    // std::cout << "SUMS:\n";
    // sums->Print();

    TFile outFile("data/detector_data.root", "RECREATE");
    sums->Write("ring_sums");
    tpcMultiplicity->Write("tpc_multiplicity");
    tofMultiplicity->Print();
    tofMultiplicity->Write("tof_multiplicity");
    outFile.Close();

    delete sums;
    delete tpcMultiplicity;
    delete tofMultiplicity;
    
    std::cout << "Analysis complete" << std::endl;
    
    
}
