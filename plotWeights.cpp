/**
 * \brief Compares the weights generated from detector data and UrQMD
 *        data
 * 
 * \author Tristan Protzman
 * \date November 17, 2020
 * \email tlprotzman@gmail.com
 * \affiliation Lehigh University
 * 
 */

#include <TROOT.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <TLegend.h>

void plotWeights() {
    TFile detector_data("data/epd_tpc_relations.root");
    TFile simulator_data("data/epd_tpc_relations_simulated.root");
    TMatrixD *detector_weights, *detector_weights_outer, *simulator_weights, *simulator_weights_outer;

    detector_data.GetDirectory("methods")->GetObject("linear_weights", detector_weights);
    if (detector_weights == nullptr) {
        printf("Could not read detector linear weights, why?\n");
        exit(1);
    }

    detector_data.GetDirectory("methods")->GetObject("linear_weights_outer", detector_weights_outer);
    if (detector_weights_outer == nullptr) {
        printf("Could not read detector linear weights, why?\n");
        exit(1);
    }

    simulator_data.GetDirectory("methods")->GetObject("linear_weights", simulator_weights);
    if (simulator_weights == nullptr) {
        printf("Could not read simulator linear weights, why?\n");
        exit(1);
    }

    simulator_data.GetDirectory("methods")->GetObject("linear_weights_outer", simulator_weights_outer);
    if (simulator_weights == nullptr) {
        printf("Could not read simulator outer ring linear weights, why?\n");
        exit(1);
    }

    TGraph *detector_weight_graph = new TGraph(16);
    TGraph *detector_weight_outer_graph = new TGraph(16);
    TGraph *simulator_weight_graph = new TGraph(16);
    TGraph *simulator_weight_outer_graph = new TGraph(16);
    TGraph *zero = new TGraph(2);


    for (uint32_t i = 0; i < 16; i++) {
        detector_weight_graph->SetPoint(i, i + 1, (*detector_weights)[i][0]);
        simulator_weight_graph->SetPoint(i, i + 1, (*simulator_weights)[i][0]);
        // printf("%d\t%f\n", i + 1, (*detector_weights)[i][0]);
    }
    uint32_t outer_sim_rings_used = 0;
    for (uint32_t i = 0; i < 16; i++) {
        if ((*simulator_weights_outer)[i][0] != 0) {
            simulator_weight_outer_graph->SetPoint(outer_sim_rings_used, i + 1, (*simulator_weights_outer)[i][0]);
            outer_sim_rings_used++;
            printf("%d\n", i);
        }
    }
    simulator_weight_outer_graph->Set(outer_sim_rings_used);

    uint32_t outer_det_rings_used = 0;
    for (uint32_t i = 0; i < 16; i++) {
        if ((*detector_weights_outer)[i][0] != 0) {
            detector_weight_outer_graph->SetPoint(outer_det_rings_used, i + 1, (*detector_weights_outer)[i][0]);
            outer_det_rings_used++;
            printf("%d\n", i);
        }
    }
    detector_weight_outer_graph->Set(outer_det_rings_used);

    zero->SetPoint(0, 0, 0);
    zero->SetPoint(1, 16, 0);

    detector_weight_graph->SetTitle(Form("Detector Weights, Bias=%f", (*detector_weights)[16][0]));
    simulator_weight_graph->SetTitle(Form("Simulator Weights, Bias=%f", (*simulator_weights)[16][0]));
    simulator_weight_outer_graph->SetTitle(Form("Simulator Weights, Outer %d Rings, Bias=%f", outer_sim_rings_used, (*simulator_weights_outer)[16][0]));
    detector_weight_outer_graph->SetTitle(Form("Detector Weights, Outer %d Rings, Bias=%f", outer_sim_rings_used, (*detector_weights_outer)[16][0]));

    TCanvas *canvas = new TCanvas("Weights");
    TMultiGraph *graph = new TMultiGraph();
    TLegend *legend;

    detector_weight_graph->SetLineColor(kBlue);
    simulator_weight_graph->SetLineColor(kRed);
    simulator_weight_outer_graph->SetLineColor(kGreen + 2);
    detector_weight_outer_graph->SetLineColor(kOrange);

    detector_weight_graph->SetMarkerStyle(24);
    simulator_weight_graph->SetMarkerStyle(25);
    simulator_weight_outer_graph->SetMarkerStyle(26);
    detector_weight_outer_graph->SetMarkerStyle(27);

    detector_weight_graph->SetMarkerSize(1.5);
    simulator_weight_graph->SetMarkerSize(1.5);
    simulator_weight_outer_graph->SetMarkerSize(1.5);
    detector_weight_outer_graph->SetMarkerSize(1.5);

    graph->Add(detector_weight_graph);
    graph->Add(detector_weight_outer_graph);
    graph->Add(simulator_weight_graph);
    graph->Add(simulator_weight_outer_graph);

    graph->SetTitle("Linear Weighting Weights");
    graph->GetXaxis()->SetTitle("Ring");
    graph->GetYaxis()->SetTitle("Weight");

    graph->Draw("alp");
    legend = canvas->BuildLegend(0.13, 0.6, 0.43, 0.8);
    legend->Draw();
    zero->Draw("l");
    canvas->Draw();

    detector_data.Close();
    simulator_data.Close();
}