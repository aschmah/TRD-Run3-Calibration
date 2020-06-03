#include <TFile.h>
// #include <iostream>
// #include <fstream>
// #include "TString.h"
// #include "TChain.h"
// #include "TObject.h"


TFile* calib_file;
TFile* outputfile;

vector<TGraph*> vec_tgraph_rms;

// without v_fit check
// Int_t Defect_TRD_detectors[167] = {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 36, 40, 43, 47, 49, 50, 51, 55, 59, 64, 88, 92, 113, 116, 119, 123, 125, 129, 130,
// 131, 132, 136, 137, 140, 141, 142, 143, 144, 146, 148, 149, 150, 156, 157, 159, 162, 163, 164, 169, 171, 175, 180, 181, 190, 191,
// 194, 197, 207, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 236, 238, 241, 245, 249, 255, 265, 274, 277, 287,
// 302, 304, 305, 308, 309, 310, 311, 317, 318, 319, 320, 323, 326, 328, 335, 348, 365, 368, 371, 377, 386, 389, 391, 395, 397, 400,
// 401, 402, 403, 404, 405, 406, 407, 413, 417, 419, 425, 427, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 442, 443, 445,
// 447, 448, 449, 452, 455, 456, 461, 462, 463, 464, 465, 466, 467, 470, 476, 482, 483, 484, 485, 490, 491, 493, 497, 500, 502, 504,
// 520, 526, 533, 536, 537, 538};

// with v_fit check
Int_t Defect_TRD_detectors[271] =  {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 47, 48,
49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 64, 88, 92, 113, 114, 115, 116, 117, 118, 119, 123, 125, 129, 130, 131, 132, 136, 137,
140, 141, 142, 143, 144, 146, 148, 149, 150, 156, 157, 159, 162, 163, 164, 169, 171, 175, 180, 181, 182, 183, 184, 185, 190, 191, 
194, 197, 207, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 241, 245, 249, 255, 
265, 274, 277, 287, 302, 304, 305, 306, 307, 308, 309, 310, 311, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 
335, 348, 365, 368, 371, 377, 384, 385, 386, 387, 388, 389, 391, 395, 397, 400, 401, 402, 403, 404, 405, 406, 407, 413, 417, 419, 
425, 427, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 442, 443, 445, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 
461, 462, 463, 464, 465, 466, 467, 470, 476, 480, 481, 482, 483, 484, 485, 490, 491, 493, 497, 498, 499, 500, 501, 502, 503, 504, 
520, 526, 533, 536, 537, 538};


void get_vec_tgraph_rms()
{
    TH2D* h2d_delta_alpha_det;
    TH1D* y_projection_bin;

    vec_tgraph_rms.resize(540);
    for (int i_det=0; i_det<540; i_det++)
    {
        vec_tgraph_rms[i_det] = new TGraph();
    }

    TDirectory *delta_impact_circle_dir = (TDirectory *)calib_file->Get("Delta_impact_circle;1");
    delta_impact_circle_dir->cd();

    Double_t x;
    Double_t y;

    int i_det = 0;
    for (auto&& keyAsObj : *delta_impact_circle_dir->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        if (key->ReadObj()->IsA()->InheritsFrom("TH2")) 
        {
            h2d_delta_alpha_det = (TH2D *)key->ReadObj();
            // cout << endl << key->GetName() << endl;

            for (int i_bin=1; i_bin<14; i_bin++)
            {
                Double_t center = h2d_delta_alpha_det->GetXaxis()->GetBinCenter(i_bin);

                y_projection_bin = (TH1D *)h2d_delta_alpha_det->ProjectionY("y_projection_bin", i_bin, i_bin); 
                Double_t rms = y_projection_bin->GetRMS();

                vec_tgraph_rms[i_det]->SetPoint(i_bin, center, rms);
            }
            i_det ++;
        }
    }
}


void draw_vec_tgraph_rms()
{
    TString HistName;
    char NoP[50];

    vector<TCanvas*> vec_tgraph_can_rms;
    vec_tgraph_can_rms.resize(36); // 12 sector blocks with 3 sectors each (18)

    TH1D* h2D_dummy_Delta_vs_impact_circle = new TH1D("h2D_dummy_Delta_vs_impact_circle","h2D_dummy_Delta_vs_impact_circle",90,50,140);
    h2D_dummy_Delta_vs_impact_circle->SetStats(0);
    h2D_dummy_Delta_vs_impact_circle->SetTitle("");
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleOffset(0.85);
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleOffset(0.78);
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelOffset(0.0);
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelOffset(0.01);
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelSize(0.08);
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelSize(0.08);
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleSize(0.08);
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleSize(0.08);
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetNdivisions(505,'N');
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetNdivisions(505,'N');
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->CenterTitle();
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->CenterTitle();
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitle("impact angle");
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitle("#Delta #alpha");
    h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetRangeUser(70,110);
    h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetRangeUser(0,15);

    for(Int_t i_sec_block = 0; i_sec_block < 36; i_sec_block++)
    {
        HistName = "vec_tgraph_can_rms";
        HistName += i_sec_block;
        vec_tgraph_can_rms[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_tgraph_can_rms[i_sec_block] ->Divide(5,3);


        for(Int_t i_det = i_sec_block*15; i_det < i_sec_block*15 + 15; i_det++)
        {
            Int_t det_colour = kBlack;

            bool is_defective = find(begin(Defect_TRD_detectors), end(Defect_TRD_detectors), i_det) != end(Defect_TRD_detectors);
            if (is_defective){
                if (i_det != 0)
                {
                    det_colour = kRed;
                }
            }
            
            Int_t iPad = i_det % 15 + 1;


            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetTicks(1,1);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetGrid(0,0);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetFillColor(10);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
            vec_tgraph_can_rms[i_sec_block] ->cd(iPad);

            h2D_dummy_Delta_vs_impact_circle->Draw("h");

            vec_tgraph_rms[i_det] ->Draw("same C*");


            HistName = "Det: ";
            HistName += i_det;

            TLatex* text = new TLatex(0.75, 0.9, HistName.Data());
            text->SetTextFont(42);
            text->SetTextSize(0.055);
            text->SetNDC();
            text->SetTextColor(det_colour);
            text->Draw();
        }

        HistName = "";
        sprintf(NoP,"%s%d","half_sector_", (Int_t)i_sec_block);
        HistName += NoP;
        HistName += ".png";
        vec_tgraph_can_rms[i_sec_block]->Print((char*)HistName.Data());
    }

    outputfile ->cd();
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
       vec_tgraph_rms[i_det]->Write();
    }
}


void h2d_error_analysis()
{
    calib_file = TFile::Open("../Data/TRD_Calib_circle_56.root");
    outputfile = new TFile("./delta_alpha_rms.root", "RECREATE");

    get_vec_tgraph_rms();
    draw_vec_tgraph_rms();
}
