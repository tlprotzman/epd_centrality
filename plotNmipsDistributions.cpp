/**
 * \brief Plots the nmip distribution for each ring in the EPD detector
 *        for simulated and real data   
 * 
 * \author Tristan Protzman
 * \date October 20, 2020
 * \email tlprotzman@gmail.com
 * \affiliation Lehigh University
 * 
 */

#include <TROOT.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <THStack.h>
#include <TLegend.h>
#include <TText.h>

#include <iostream>

const int RINGS = 16;

void plotNmipsDistributions() {
    // Open data files
    TFile simulated_data("data/simulated_data.root");
    TFile detector_data("data/detector_data.root");

    // Load simulated data
    TMatrixD *sim_nmips;
    TVectorD *sim_impact_parameter;
    TVectorD *sim_refmult1;
    simulated_data.GetObject("ring_sums", sim_nmips);
    simulated_data.GetObject("tpc_multiplicity", sim_refmult1);
    simulated_data.GetObject("impact_parameter", sim_impact_parameter);
    simulated_data.Close();

    // Load detector data
    TMatrixD *det_nmips;
    TVectorD *det_refmult1;
    detector_data.GetObject("ring_sums", det_nmips);
    detector_data.GetObject("tpc_multiplicity", det_refmult1);
    detector_data.Close();

    uint32_t num_bins = 60;
    int32_t lower_bin = 0;
    int32_t upper_bin = 60;

    // Create and fill histograms
    TH1D **sim_histograms = (TH1D**)malloc(RINGS * sizeof(TH1D*));
    TH1D **sim_histograms_bFiltered = (TH1D**)malloc(RINGS * sizeof(TH1D*));
    TH1D **det_histograms = (TH1D**)malloc(RINGS * sizeof(TH1D*));
    for (uint32_t i = 0; i < RINGS; i++) {
        sim_histograms[i] = new TH1D(Form("sim_nmips_ring_%d", i + 1),
                                     Form("UrQMD nmips distribution, ring %d", i + 1),
                                     num_bins, lower_bin, upper_bin);
        sim_histograms_bFiltered[i] = new TH1D(Form("sim_nmips_ring_%d", i + 1),
                                     Form("UrQMD nmips distribution, ring %d, b<7.5", i + 1),
                                     num_bins, lower_bin, upper_bin);
        det_histograms[i] = new TH1D(Form("det_nmips_ring_%d", i + 1),
                                     Form("Detector nmips distribution, ring %d", i + 1),
                                     num_bins, lower_bin, upper_bin);

        sim_histograms[i]->SetXTitle("nmips");
        sim_histograms[i]->SetYTitle("Count");
        sim_histograms_bFiltered[i]->SetXTitle("nmips");
        sim_histograms_bFiltered[i]->SetYTitle("Count");
        det_histograms[i]->SetXTitle("nmips");
        det_histograms[i]->SetYTitle("Count");

        // Fill data    
        // Simulation    
        for (uint32_t j = 0; j < sim_nmips->GetNcols(); j++) {
            // if ((*sim_impact_parameter)[j] < 7.5) {
                sim_histograms[i]->Fill((*sim_nmips)[i][j]);//, 1. / sim_nmips->GetNcols());
            // }
        }
        for (uint32_t j = 0; j < sim_nmips->GetNcols(); j++) {
            if ((*sim_impact_parameter)[j] < 7.5) {
                sim_histograms_bFiltered[i]->Fill((*sim_nmips)[i][j]);//, 1. / sim_nmips->GetNcols());
            }
        }
        // Detector Data
        for (uint32_t j = 0; j < det_nmips->GetNcols(); j++) {
            // if ((*det_refmult1)[j] > 25) {
                det_histograms[i]->Fill((*det_nmips)[i][j]);//, 1. / det_nmips->GetNcols());
            // }
        }
    }

    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    TCanvas *canvas = new TCanvas("canvas", "Nmips Comparison", 1000, 1000);
    canvas->Divide(4, 4);
    for (uint32_t i = 0; i < RINGS; i++) {
        TLegend *legend = new TLegend(0.75,0.7,0.95,0.85);
        THStack *stack = new THStack(Form("sim_nmips_ring_%d", i + 1), Form("nmips distribution, ring %d", i + 1));
        canvas->cd(i + 1);

        std::cout << "Ring " << i + 1 << " Max Bin: " << det_histograms[i]->GetMaximumBin() << "\t" << sim_histograms[i]->GetMaximumBin() << std::endl;
        sim_histograms[i]->Scale(1. / sim_histograms[i]->Integral());        
        sim_histograms_bFiltered[i]->Scale(1. / sim_histograms_bFiltered[i]->Integral());        
        det_histograms[i]->Scale(1. / det_histograms[i]->Integral());        

        sim_histograms[i]->SetLineColor(1);
        sim_histograms[i]->SetMarkerColor(1);
        sim_histograms[i]->SetMarkerStyle(24);
        sim_histograms[i]->SetMarkerSize(0.5);
        stack->Add(sim_histograms[i]);
        sim_histograms_bFiltered[i]->SetLineColor(kGreen + 2);
        sim_histograms_bFiltered[i]->SetMarkerColor(kGreen + 2);
        sim_histograms_bFiltered[i]->SetMarkerStyle(32);
        sim_histograms_bFiltered[i]->SetMarkerSize(0.5);
        stack->Add(sim_histograms_bFiltered[i]);
        det_histograms[i]->SetLineColor(2);
        det_histograms[i]->SetMarkerColor(2);
        det_histograms[i]->SetMarkerStyle(30);
        det_histograms[i]->SetMarkerSize(0.7);
        stack->Add(det_histograms[i]);
        
        stack->Draw("nostack hist p");
        stack->GetXaxis()->SetTitle("nmips");
        stack->GetYaxis()->SetTitle("count");
        stack->GetXaxis()->SetTickLength(0);
        stack->GetYaxis()->SetTickLength(0);
        
        legend->AddEntry(sim_histograms[i], "UrQMD");
        legend->AddEntry(sim_histograms_bFiltered[i], "UrQMD, b<7.5");
        legend->AddEntry(det_histograms[i], "Detector");
        legend->Draw();
    }
    canvas->Draw();
    canvas->SaveAs("histograms/nmips_distributions.png");

    // Plotting nmips vs refmult1 for detector data
    int32_t refmult1_bins, refmult1_min, refmult1_max;
    int32_t nmips_bins, nmips_min, nmips_max;
    refmult1_bins = 50;
    refmult1_min = 0;
    refmult1_max = 300;
    nmips_bins = 50;
    nmips_min = 0;
    nmips_max = 60;

    TH2D **det_nmips_refmult1 = (TH2D**)malloc(RINGS * sizeof(TH2D*));
    for (uint32_t i = 0; i < RINGS; i++) {
        det_nmips_refmult1[i] = new TH2D(Form("det_nmips_refmult_%d)", i+1), Form("nMIPs vs RefMult1, Detector, ring %d", i + 1),
                                     refmult1_bins, refmult1_min, refmult1_max,
                                     nmips_bins, nmips_min, nmips_max);
        det_nmips_refmult1[i]->SetXTitle("refmult1");
        det_nmips_refmult1[i]->SetYTitle("nmips");
        for (uint32_t j = 0; j < det_nmips->GetNcols(); j++) {
            det_nmips_refmult1[i]->Fill((*det_refmult1)[j], (*det_nmips)[i][j]);
        }
    }

    TCanvas *canvas2 = new TCanvas("canvas2", "RefMult1 vs nMIPs", 1000, 1000);
    canvas2->Divide(4, 4);
    for (uint32_t i = 0; i < RINGS; i++) {
        canvas2->cd(i + 1);
        gPad->SetLogz();
        det_nmips_refmult1[i]->Draw("colz");
    }
    canvas2->Draw();
    canvas2->SaveAs("histograms/det_nmips_refmult1.png");


    // Plotting nmips vs refmult1 for simulation
    TH2D **sim_nmips_refmult1 = (TH2D**)malloc(RINGS * sizeof(TH2D*));
    for (uint32_t i = 0; i < RINGS; i++) {
        sim_nmips_refmult1[i] = new TH2D(Form("sim_nmips_refmult_%d)", i+1), Form("nMIPs vs RefMult1, UrQMD, ring %d", i + 1),
                                     refmult1_bins, refmult1_min, refmult1_max,
                                     nmips_bins, nmips_min, nmips_max);
        sim_nmips_refmult1[i]->SetXTitle("refmult1");
        sim_nmips_refmult1[i]->SetYTitle("nmips");
        for (uint32_t j = 0; j < sim_nmips->GetNcols(); j++) {
            sim_nmips_refmult1[i]->Fill((*sim_refmult1)[j], (*sim_nmips)[i][j]);
        }
    }

    TCanvas *canvas3 = new TCanvas("canvas3", "RefMult1 vs nMIPs", 1000, 1000);
    canvas3->Divide(4, 4);
    for (uint32_t i = 0; i < RINGS; i++) {
        canvas3->cd(i + 1);
        gPad->SetLogz();
        sim_nmips_refmult1[i]->Draw("colz");
    }
    canvas3->Draw();
    canvas3->SaveAs("histograms/sim_nmips_refmult1.png");

    // Plotting 3 rings
    TCanvas *canvas4 = new TCanvas("Canvas4", "RefMult1 vs nMIPs", 1000, 1000);
    canvas4->Divide(2, 3);
    TText **pearson_coefficients = (TText**)malloc(6 * sizeof(TText*));

    int ring_selection[] = {0, 2, 4};
    for (uint32_t i = 0; i < 3; i++) {
        canvas4->cd(2 * i + 1);
        gPad->SetLogz();
        det_nmips_refmult1[ring_selection[i]]->Draw("colz");
        double corr = det_nmips_refmult1[ring_selection[i]]->GetCorrelationFactor();
        pearson_coefficients[2 * i] = new TText();
        pearson_coefficients[2 * i]->DrawText(120, 50, Form("Pearsons Coefficient: %f", corr));


        canvas4->cd(2 * i + 2);
        gPad->SetLogz();
        sim_nmips_refmult1[ring_selection[i]]->Draw("colz");
        corr = sim_nmips_refmult1[ring_selection[i]]->GetCorrelationFactor();
        pearson_coefficients[2 * i + 1] = new TText();
        pearson_coefficients[2 * i + 1]->DrawText(120, 50, Form("Pearsons Coefficient: %f", corr));
    }

    canvas4->Draw();
    canvas4->SaveAs("histograms/det_sim_comparison.png");

}