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
    TH2* hRingvsRegMult[2][16];
    
    // 1D histograms for nMIP distributions
    for (int ew=0; ew<2; ew++){
        for (int pp=1; pp<13; pp++){
            for (int tt=1; tt<32; tt++){
                mNmipDists[ew][pp-1][tt-1] = (TH1*)gROOT->FindObject(Form("NmipEW%dPP%dTT%d",ew,pp,tt));
                mAdcDists[ew][pp-1][tt-1] = (TH1*)gROOT->FindObject(Form("ADCEW%dPP%dTT%d",ew,pp,tt));
            }
        }
        for (int r = 0;r<16;r++){
            hRingvsRegMult[ew][r] = (TH2*)gROOT->FindObject(Form("hRingvsRegMultEW%iRing%i",ew,r+1));
        }
    }
    
    TCanvas* can1 = new TCanvas("can1","",4000,4000);
    can1->Divide(4,4);
    for (int ring = 1;ring < 17; ring++){
        can1->cd(ring);
        gPad->SetLogz();
        hRingvsRegMult[0][ring-1]->Draw("COLZ");
        hRingvsRegMult[0][ring-1]->GetYaxis()->SetRangeUser(0,100);
        hRingvsRegMult[0][ring-1]->GetXaxis()->SetRangeUser(0,300);
        TLegend* leg1 = new TLegend(0.5,0.5,0.85,0.85);
        SetLeg(leg1,0.07);
        leg1->AddEntry(hRingvsRegMult[0][ring-1],"East","");
        leg1->AddEntry(hRingvsRegMult[0][ring-1],Form("Ring %i",ring),"");
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

