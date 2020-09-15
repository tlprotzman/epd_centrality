/*
* This seems to as the name implies just draw the histogram from the preprocessed
* data.
*/

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <string>

void SetHist(TH1* h, std::string xt ="", std::string yt ="",int color = 1, int marker = 20,int width = 3, float size = 1.0);
void SetHist(TH1* h, int color = 1);
void SetLeg(TLegend* l,float txtsize=0.03);

void DrawEpdHist(std::string infile = "out.root"){
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *_file0 = TFile::Open(infile.c_str());
    TH1* mNmipDists[2][12][31];
    TH1* mAdcDists[2][12][31];
    TH2* hRingvsRegMult[2][16]; // The distribution that's getting plotted, ring multiplicity vs TPC multiplicity?
    
    // 1D histograms for nMIP distributions  
    for (int east_west=0; east_west<2; east_west++){    // 0 = east, 1 = west
        // for (int pp=1; pp<13; pp++){
        //     for (int tt=1; tt<32; tt++){
        //         mNmipDists[east_west][pp-1][tt-1] = (TH1*)gROOT->FindObject(Form("NmipEW%dPP%dTT%d",east_west,pp,tt));
        //         mAdcDists[east_west][pp-1][tt-1] = (TH1*)gROOT->FindObject(Form("ADCEW%dPP%dTT%d",east_west,pp,tt));
        //     }
        // }
        for (int r = 0;r<16;r++){
            hRingvsRegMult[east_west][r] = (TH2*)gROOT->FindObject(Form("hRingvsRegMultEW%iRing%i",east_west,r+1));
        }
    }
    
    TCanvas* can1 = new TCanvas("can1","",4000,4000);
    can1->Divide(4,4);
    const int east_west = 1;
    for (int ring = 1;ring < 17; ring++){
        can1->cd(ring);
        gPad->SetLogz();
        hRingvsRegMult[east_west][ring-1]->Draw("COLZ");
        hRingvsRegMult[east_west][ring-1]->GetYaxis()->SetRangeUser(0,100);
        hRingvsRegMult[east_west][ring-1]->GetXaxis()->SetRangeUser(0,300);
        TLegend* leg1 = new TLegend(0.5,0.5,0.85,0.85);
        SetLeg(leg1,0.07);
        leg1->AddEntry(hRingvsRegMult[east_west][ring-1],(east_west == 0 ?"East" : "West"), "");
        leg1->AddEntry(hRingvsRegMult[east_west][ring-1],Form("Ring %i",ring),"");
        leg1->Draw();
    }
    can1->Draw();
    
}

//*******************************************************************************************

void SetHist(TH1* h, std::string xt ="", std::string yt ="",int color = 1, int marker = 20,int width = 3, float size = 1.0)
{
    h->SetLineWidth(width);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerSize(size);
    h->SetMarkerStyle(marker);
    h->GetYaxis()->SetTitle(yt.c_str());
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetXaxis()->SetTitle(xt.c_str());
}

void SetHist(TH1* h, int color = 1)
{
    h->SetLineWidth(3);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerSize(1);
    h->GetYaxis()->SetTitleOffset(1.6);
    
}


void SetLeg(TLegend* l,float txtsize=0.03){
    l->SetBorderSize(0);
    l->SetFillColor(0);
    l->SetTextSize(txtsize);
}

