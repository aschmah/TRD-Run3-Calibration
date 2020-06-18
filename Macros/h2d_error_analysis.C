#include <TFile.h>
// #include <iostream>
// #include <fstream>
// #include "TString.h"
// #include "TChain.h"
// #include "TObject.h"


TFile* TBase_calibrate_hist_output;
TFile* rms_output;

TGraph *tg_anode_hv;
TGraph *tg_drift_hv;
TGraph *vdrift_fit;

TString HistName;

vector<TGraph*> vec_tgraph_rms;

TGraph* tg_mean_rms_good_chamber;

int det_good_suspicious[540] = {};

// without v_fit check
// Int_t Defect_TRD_detectors[167] = {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 36, 40, 43, 47, 49, 50, 51, 55, 59, 64, 88, 92, 113, 116, 119, 123, 125, 129, 130,
// 131, 132, 136, 137, 140, 141, 142, 143, 144, 146, 148, 149, 150, 156, 157, 159, 162, 163, 164, 169, 171, 175, 180, 181, 190, 191,
// 194, 197, 207, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 236, 238, 241, 245, 249, 255, 265, 274, 277, 287,
// 302, 304, 305, 308, 309, 310, 311, 317, 318, 319, 320, 323, 326, 328, 335, 348, 365, 368, 371, 377, 386, 389, 391, 395, 397, 400,
// 401, 402, 403, 404, 405, 406, 407, 413, 417, 419, 425, 427, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 442, 443, 445,
// 447, 448, 449, 452, 455, 456, 461, 462, 463, 464, 465, 466, 467, 470, 476, 482, 483, 484, 485, 490, 491, 493, 497, 500, 502, 504,
// 520, 526, 533, 536, 537, 538};

// bad ADC
// Int_t bad_adc[137] = {2, 5, 15, 17, 27, 31, 32, 34, 36, 40, 43, 49, 50, 55, 59, 62, 64, 88, 92, 113, 116, 119, 125, 129, 130, 
// 132, 136, 137, 140, 141, 142, 143, 146, 148, 149, 156, 159, 162, 163, 164, 168, 171, 175, 180, 181, 190, 191, 194, 206, 207, 
// 213, 214, 219, 221, 226, 228, 230, 231, 236, 238, 241, 249, 255, 265, 273, 274, 277, 302, 308, 309, 310, 311, 317, 318, 319, 
// 320, 326, 328, 335, 348, 365, 368, 371, 377, 386, 389, 400, 401, 402, 403, 404, 405, 406, 407, 413, 417, 425, 429, 431, 432, 
// 433, 434, 435, 436, 437, 442, 443, 445, 447, 448, 449, 452, 455, 456, 461, 462, 463, 464, 465, 466, 467, 470, 482, 483, 484, 
// 485, 490, 491, 493, 497, 500, 502, 504, 526, 536, 537, 538};

Int_t bad_adc[140] = {2, 4, 5, 15, 17, 27, 31, 32, 34, 36, 40, 43, 49, 50, 55, 59, 62, 64, 88, 92, 113, 116, 119, 125, 129, 130, 
132, 136, 137, 140, 141, 142, 143, 146, 148, 149, 156, 159, 162, 163, 164, 168, 171, 175, 180, 181, 184, 190, 191, 194, 206, 207, 
213, 214, 219, 221, 226, 228, 230, 231, 236, 238, 241, 242, 249, 255, 265, 273, 274, 277, 302, 308, 309, 310, 311, 317, 318, 319, 
320, 326, 328, 335, 348, 365, 368, 371, 377, 386, 389, 400, 401, 402, 403, 404, 405, 406, 407, 413, 417, 425, 429, 431, 432, 433, 
434, 435, 436, 437, 442, 443, 445, 447, 448, 449, 452, 455, 456, 461, 462, 463, 464, 465, 466, 467, 470, 482, 483, 484, 485, 490, 
491, 493, 497, 500, 502, 504, 526, 536, 537, 538};

// with v_fit check
Int_t Defect_TRD_detectors[215] =  {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 47, 48,
49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 64, 88, 92, 113, 114, 115, 116, 117, 118, 119, 123, 125, 129, 130, 131, 132, 136, 137,
140, 141, 142, 143, 144, 146, 148, 149, 150, 156, 157, 159, 162, 163, 164, 169, 171, 175, 180, 181, 182, 183, 184, 185, 190, 191, 
194, 197, 207, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 241, 245, 249, 255, 
265, 274, 277, 287, 302, 304, 305, 306, 307, 308, 309, 310, 311, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 
335, 348, 365, 368, 371, 377, 384, 385, 386, 387, 388, 389, 391, 395, 397, 400, 401, 402, 403, 404, 405, 406, 407, 413, 417, 419, 
425, 427, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 442, 443, 445, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 
461, 462, 463, 464, 465, 466, 467, 470, 476, 480, 481, 482, 483, 484, 485, 490, 491, 493, 497, 498, 499, 500, 501, 502, 503, 504, 
520, 526, 533, 536, 537, 538};

Int_t official_qa[91] = {2, 15, 17, 27, 31, 32, 36, 40, 43, 49, 50, 55, 59, 64, 88, 92, 113, 116, 119, 132, 180, 181, 190,
191, 194, 207, 215, 219, 221, 226, 227, 228, 230, 231, 233, 236, 238, 241, 249, 255, 277, 287, 302, 308, 310, 311, 317, 318, 319, 
320, 326, 328, 335, 348, 368, 377, 386, 389, 402, 403, 404, 405, 406, 407, 432, 433, 434, 435, 436, 437, 452, 455, 456, 462, 463, 
464, 465, 466, 467, 470, 482, 483, 484, 485, 490, 491, 493, 500, 502, 504, 538};


void draw_2d_delta_alpha()
{
    TLine* line_delta_equal_zero = new TLine(70,0,110,0);
    line_delta_equal_zero -> SetLineWidth(2);
    line_delta_equal_zero -> SetLineColor(kGreen);

    vector<TH2D*> vec_TH2D_Delta_vs_impact_circle;
    vec_TH2D_Delta_vs_impact_circle.resize(540);

    TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_impact_circle;1");
    delta_impact_circle_dir->cd();

    Double_t x;
    Double_t y;

    int i_det = 0;
    for (auto&& keyAsObj : *delta_impact_circle_dir->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        if (key->ReadObj()->IsA()->InheritsFrom("TH2")) 
        {
            vec_TH2D_Delta_vs_impact_circle[i_det] = (TH2D *)key->ReadObj();
            i_det ++;
        }
    }

    vector<TCanvas*> vec_can_2D_Delta_vs_impact_circle;
    vec_can_2D_Delta_vs_impact_circle.resize(30);

    int rows = 3;
    int columns = 6;
    int n_pads = rows*columns;
    for(Int_t i_canvas = 0; i_canvas < 30; i_canvas++)
    {
        HistName = "vec_can_2D_Delta_vs_impact_circle";
        HistName += i_canvas;
        vec_can_2D_Delta_vs_impact_circle[i_canvas] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_can_2D_Delta_vs_impact_circle[i_canvas] ->Divide(columns,rows); 

        for(Int_t i_det = i_canvas*n_pads; i_det < i_canvas*n_pads + n_pads; i_det++)
        {
            Double_t top_left = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,4,110,124);
            Double_t top_right = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(8,11,110,124);
            double x_pos_text = 0.75;
            if (top_left < top_right)
            {
                x_pos_text = 0.25;
            }

            Double_t total_inner;
            Double_t total_outer;
            
            bool draw_not = 0;
            bool draw_outer = 0;
            total_inner = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,75,125);
            // cout << "det: " << i_det << " | " << "total: " << total << endl;
            if (total_inner < 20)
            {
                total_outer = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,125,200);
                total_outer += vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,1,75);
                if (total_outer > total_inner && total_outer > 10)
                {
                    draw_outer = 1;
                    x_pos_text = 0.75;
                }
                else if (total_inner < 1)
                {
                    draw_not = 1;
                }
            }

            // get HV
            double x;
            double y;
            double anode_hv;
            Int_t text_anode_hv_colour = kBlack;

            x = -1;
            for (int point = 0; point < tg_anode_hv->GetN(); point++)
            {
                tg_anode_hv->GetPoint(point, x, y);
                if (x == i_det)
                {
                    anode_hv = y;
                    if (anode_hv < 1400)
                    {
                        text_anode_hv_colour = kOrange-3;
                    }
                    if (anode_hv < 10)
                    {
                        text_anode_hv_colour = kRed;
                    }
                    goto next;
                }
            }
            text_anode_hv_colour = kRed;
            next:

            double drift_hv;
            Int_t text_drift_hv_colour = kBlack;
            x = -1;
            for (int point = 0; point < tg_drift_hv->GetN(); point++)
            {
                tg_drift_hv->GetPoint(point, x, y);
                if (x == i_det)
                {
                    drift_hv = y;
                    if (drift_hv < 1920)
                    {
                        text_drift_hv_colour = kOrange-3;
                    }
                    if (drift_hv < 10)
                    {
                        text_drift_hv_colour = kRed;
                    }
                    goto next2;
                }
            }
            text_drift_hv_colour = kRed;
            next2:

            // get ADC
            if (i_det == 184)
            {
                cout << "here" << endl;
            }
            bool has_bad_adc = find(begin(bad_adc), end(bad_adc), i_det) != end(bad_adc);

            if (i_det == 184)
            {
                cout << has_bad_adc << endl;
            }
            bool has_been_calibrated = 0;
            for (Int_t vdrift_det=0; vdrift_det<540; vdrift_det++)
            {
                vdrift_fit->GetPoint(vdrift_det, x, y);
                if ((Int_t)x == i_det)
                {
                    has_been_calibrated = 1;
                }
            }

            Int_t det_colour = kBlack;
            if (!has_been_calibrated)
            {
                det_colour = kRed;
            }

            Int_t iPad = i_det % n_pads + 1;

            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetTicks(1,1);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetGrid(0,0);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetFillColor(10);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetRightMargin(0.01);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetTopMargin(0.01);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetBottomMargin(0.2);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad)->SetLeftMargin(0.2);
            vec_can_2D_Delta_vs_impact_circle[i_canvas] ->cd(iPad);
           
            vec_TH2D_Delta_vs_impact_circle[i_det]->SetStats(0);
            vec_TH2D_Delta_vs_impact_circle[i_det]->SetTitle("");
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetTitleOffset(0.85);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetTitleOffset(0.78);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetLabelOffset(0.0);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetLabelOffset(0.01);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetLabelSize(0.08);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetLabelSize(0.08);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetTitleSize(0.08);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetTitleSize(0.08);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetNdivisions(505,'N');
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetNdivisions(505,'N');
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->CenterTitle();
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->CenterTitle();
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetTitle("impact angle");
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetTitle("#Delta #alpha");
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetRangeUser(70,110);

            if (draw_outer)
            {
                vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetRangeUser(-10,90);
            }
            else
            {
                vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetRangeUser(-24,24);
            }
            
            bool is_defective = find(begin(official_qa), end(official_qa), i_det) != end(official_qa);
            if (is_defective)
            {
                vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetAxisColor(kRed);
                vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetAxisColor(kRed);
            }

            if (!draw_not) {vec_TH2D_Delta_vs_impact_circle[i_det] ->Draw("colz");}
            else {vec_TH2D_Delta_vs_impact_circle[i_det] ->Draw("AXIS");}

            line_delta_equal_zero->Draw();

            TLatex* det_number = new TLatex(x_pos_text, 0.925, Form("Det: %d", (Int_t)i_det));
            det_number->SetTextFont(42);
            det_number->SetTextSize(0.06);
            det_number->SetNDC();
            det_number->SetTextColor(det_colour);
            det_number->Draw();

            // draw HV
            TLatex* text_anode_hv = new TLatex(x_pos_text, 0.875, Form("A: %d", (Int_t)anode_hv));
            text_anode_hv->SetTextFont(42);
            text_anode_hv->SetTextSize(0.06);
            text_anode_hv->SetNDC();
            text_anode_hv->SetTextColor(text_anode_hv_colour);
            text_anode_hv->Draw();

            TLatex* text_drift_hv = new TLatex(x_pos_text, 0.825, Form("D: %d", (Int_t)drift_hv));
            text_drift_hv->SetTextFont(42);
            text_drift_hv->SetTextSize(0.06);
            text_drift_hv->SetNDC();
            text_drift_hv->SetTextColor(text_drift_hv_colour);
            text_drift_hv->Draw();

            if (has_bad_adc)
            {
                TLatex* text_adc = new TLatex(x_pos_text, 0.775, "ADC");
                text_adc->SetTextFont(42);
                text_adc->SetTextSize(0.06);
                text_adc->SetNDC();
                text_adc->SetTextColor(kRed);
                text_adc->Draw();
            }
        }
        vec_can_2D_Delta_vs_impact_circle[i_canvas]->Print(Form("canvas%d_%d_to_%d.png", i_canvas, i_canvas*n_pads, i_canvas*n_pads + n_pads));
    }
}


void get_vec_tgraph_rms()
{
    TH2D* h2d_delta_alpha_det;
    TH1D* y_projection_bin;

    vec_tgraph_rms.resize(540);
    for (int i_det=0; i_det<540; i_det++)
    {
        vec_tgraph_rms[i_det] = new TGraph();
    }

    TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_impact_circle;1");
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

                // for (int bin=0; bin<200; bin++)
                // {
                //     cout << "bin: " << bin << " | " << "fills: " << y_projection_bin->GetBinContent(bin) << endl;;
                // }
                Double_t rms = y_projection_bin->GetRMS();

                vec_tgraph_rms[i_det]->SetPoint(i_bin, center, rms);
                // break;
            }
            i_det ++;
            // break;
        }
    }
}


void draw_vec_tgraph_rms()
{
    TString HistName;

    vector<TCanvas*> vec_tgraph_can_rms;
    vec_tgraph_can_rms.resize(30); // 12 sector blocks with 3 sectors each (18)


    // rms_output ->cd();
    // for(Int_t i_det = 0; i_det < 540; i_det++)
    // {
    //    vec_tgraph_rms[i_det]->Write();
    // }


    int rows = 3;
    int columns = 6;
    int n_pads = rows*columns;
    for(Int_t i_canvas = 0; i_canvas < 2; i_canvas++)
    {
        HistName = "vec_tgraph_can_rms";
        HistName += i_canvas;
        vec_tgraph_can_rms[i_canvas] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_tgraph_can_rms[i_canvas] ->Divide(columns,rows); 

        for(Int_t i_det = i_canvas*n_pads; i_det < i_canvas*n_pads + n_pads; i_det++)
        {
            double x_pos_text = 0.75;

            bool draw_not = 0;
            bool draw_outer = 0;

            // get HV
            double x;
            double y;
            double anode_hv;
            Int_t text_anode_hv_colour = kBlack;

            x = -1;
            for (int point = 0; point < tg_anode_hv->GetN(); point++)
            {
                tg_anode_hv->GetPoint(point, x, y);
                if (x == i_det)
                {
                    anode_hv = y;
                    if (anode_hv < 1400)
                    {
                        text_anode_hv_colour = kOrange-3;
                    }
                    if (anode_hv < 10)
                    {
                        text_anode_hv_colour = kRed;
                    }
                    goto next;
                }
            }
            text_anode_hv_colour = kRed;
            next:

            double drift_hv;
            Int_t text_drift_hv_colour = kBlack;
            x = -1;
            for (int point = 0; point < tg_drift_hv->GetN(); point++)
            {
                tg_drift_hv->GetPoint(point, x, y);
                if (x == i_det)
                {
                    drift_hv = y;
                    if (drift_hv < 1920)
                    {
                        text_drift_hv_colour = kOrange-3;
                    }
                    if (drift_hv < 10)
                    {
                        text_drift_hv_colour = kRed;
                    }
                    goto next2;
                }
            }
            text_drift_hv_colour = kRed;
            next2:

            // get ADC
            bool has_bad_adc = find(begin(bad_adc), end(bad_adc), i_det) != end(bad_adc);

            bool has_been_calibrated = 0;
            for (Int_t vdrift_det=0; vdrift_det<540; vdrift_det++)
            {
                vdrift_fit->GetPoint(vdrift_det, x, y);
                if ((Int_t)x == i_det)
                {
                    has_been_calibrated = 1;
                }
            }

            Int_t det_colour = kBlack;
            if (!has_been_calibrated)
            {
                det_colour = kRed;
            }

            Int_t iPad = i_det % n_pads + 1;

            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetTicks(1,1);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetGrid(0,0);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetFillColor(10);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetRightMargin(0.01);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetTopMargin(0.01);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetBottomMargin(0.2);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad)->SetLeftMargin(0.2);
            vec_tgraph_can_rms[i_canvas] ->cd(iPad);
           
            vec_tgraph_rms[i_det]->SetTitle("");
            vec_tgraph_rms[i_det]->GetXaxis()->SetTitleOffset(0.85);
            vec_tgraph_rms[i_det]->GetYaxis()->SetTitleOffset(1.1);
            vec_tgraph_rms[i_det]->GetXaxis()->SetLabelOffset(0.0);
            vec_tgraph_rms[i_det]->GetYaxis()->SetLabelOffset(0.01);
            vec_tgraph_rms[i_det]->GetXaxis()->SetLabelSize(0.08);
            vec_tgraph_rms[i_det]->GetYaxis()->SetLabelSize(0.08);
            vec_tgraph_rms[i_det]->GetXaxis()->SetTitleSize(0.08);
            vec_tgraph_rms[i_det]->GetYaxis()->SetTitleSize(0.08);
            vec_tgraph_rms[i_det]->GetXaxis()->SetNdivisions(505,'N');
            vec_tgraph_rms[i_det]->GetYaxis()->SetNdivisions(505,'N');
            vec_tgraph_rms[i_det]->GetXaxis()->CenterTitle();
            vec_tgraph_rms[i_det]->GetYaxis()->CenterTitle();
            vec_tgraph_rms[i_det]->GetXaxis()->SetTitle("impact angle [deg]");
            vec_tgraph_rms[i_det]->GetYaxis()->SetTitle("#sigma_{#Delta#alpha} [deg]");
            vec_tgraph_rms[i_det]->GetXaxis()->SetRangeUser(70,110);
            vec_tgraph_rms[i_det]->GetYaxis()->SetRangeUser(0,24);
            
            
            bool is_defective = find(begin(official_qa), end(official_qa), i_det) != end(official_qa);
            if (is_defective)
            {
                vec_tgraph_rms[i_det]->GetYaxis()->SetAxisColor(kRed);
                vec_tgraph_rms[i_det]->GetXaxis()->SetAxisColor(kRed);
            }

            vec_tgraph_rms[i_det]->Draw("AC*");
            tg_mean_rms_good_chamber->SetLineStyle(2);
            tg_mean_rms_good_chamber->SetLineColor(kRed);
            tg_mean_rms_good_chamber->SetLineWidth(2);
            tg_mean_rms_good_chamber->Draw("C");


            TLatex* det_number = new TLatex(x_pos_text, 0.925, Form("Det: %d", (Int_t)i_det));
            det_number->SetTextFont(42);
            det_number->SetTextSize(0.06);
            det_number->SetNDC();
            det_number->SetTextColor(det_colour);
            det_number->Draw();

            // draw HV
            TLatex* text_anode_hv = new TLatex(x_pos_text, 0.875, Form("A: %d", (Int_t)anode_hv));
            text_anode_hv->SetTextFont(42);
            text_anode_hv->SetTextSize(0.06);
            text_anode_hv->SetNDC();
            text_anode_hv->SetTextColor(text_anode_hv_colour);
            text_anode_hv->Draw();

            TLatex* text_drift_hv = new TLatex(x_pos_text, 0.825, Form("D: %d", (Int_t)drift_hv));
            text_drift_hv->SetTextFont(42);
            text_drift_hv->SetTextSize(0.06);
            text_drift_hv->SetNDC();
            text_drift_hv->SetTextColor(text_drift_hv_colour);
            text_drift_hv->Draw();

            if (has_bad_adc)
            {
                TLatex* text_adc = new TLatex(x_pos_text, 0.775, "ADC");
                text_adc->SetTextFont(42);
                text_adc->SetTextSize(0.06);
                text_adc->SetNDC();
                text_adc->SetTextColor(kRed);
                text_adc->Draw();
            }
        }
        // vec_tgraph_can_rms[i_canvas]->Print(Form("canvas%d_%d_to_%d.png", i_canvas, i_canvas*n_pads, i_canvas*n_pads + n_pads));
    }
}


void get_mean_rms_vs_layer()
{
    TGraph* tg_rms_vs_layer = new TGraph();

    TFile* flag_file = TFile::Open("./chamber_QC_flags.root");

    vector<int> *flag_pointers;
    vector<int> flags;
    flag_file->GetObject("QC_flags", flag_pointers);

    for(vector<int>::iterator it = flag_pointers->begin(); it != flag_pointers->end(); ++it) 
    {
        // cout << *it << endl;
        flags.push_back(*it);
    }

    int layer_count[6] = {0};
    double layer_sum[6] = {0.0};

    for (int i_det=0; i_det<540; i_det++)
    {
        bool det_suspicious = (flags[i_det] > 0);
        double mean = vec_tgraph_rms[i_det]->GetMean(2);

        // cout << endl << "det: " << i_det << endl;
        // cout << "mean: " << mean << endl;

        int layer_number = i_det % 6;

        // if (flags[i_det] == 8 || flags[i_det] == 24)
        if (flags[i_det] == 8)
        // if (!det_suspicious)
        {
            layer_count[layer_number] ++;
            layer_sum[layer_number] += mean;
        }
    }


    for (int i_layer=0; i_layer<6; i_layer++) 
    {
        cout << layer_count[i_layer] << endl;
        if (layer_count[i_layer] == 0) {continue;}
        double layer_mean = layer_sum[i_layer] / layer_count[i_layer];
    
        // cout << layer_count[i_layer] << endl;

        tg_rms_vs_layer->SetPoint(i_layer, i_layer, layer_mean);
    }

    // TCanvas* can_tg_rms_vs_layer = new TCanvas("can_tg_rms_vs_layer","can_tg_rms_vs_layer",10,10,400,400);
    // can_tg_rms_vs_layer->SetRightMargin(0.01);
    // can_tg_rms_vs_layer->SetTopMargin(0.1);
    // can_tg_rms_vs_layer->SetBottomMargin(0.15);
    // can_tg_rms_vs_layer->SetLeftMargin(0.15);

    // tg_rms_vs_layer->SetTitle("Drift HV OR Drift HV + ADC");
    // tg_rms_vs_layer->SetMarkerColor(kBlue);

    // tg_rms_vs_layer->GetXaxis()->SetTitle("Layer");
    // tg_rms_vs_layer->GetYaxis()->SetTitle("Mean RMS");
    // tg_rms_vs_layer->GetXaxis()->CenterTitle();
    // tg_rms_vs_layer->GetYaxis()->CenterTitle();
    // tg_rms_vs_layer->GetXaxis()->SetLabelSize(0.05);
    // tg_rms_vs_layer->GetYaxis()->SetLabelSize(0.05);
    // tg_rms_vs_layer->GetXaxis()->SetTitleSize(0.05);
    // tg_rms_vs_layer->GetYaxis()->SetTitleSize(0.05);
    // tg_rms_vs_layer->GetYaxis()->SetRangeUser(0,15);
    // tg_rms_vs_layer->SetFillColor(kBlue);
    // tg_rms_vs_layer->Draw("AB");
}


void get_mean_rms_vs_impact_angle()
{
    TFile* flag_file = TFile::Open("./chamber_QC_flags.root");

    vector<int> *flag_pointers;
    vector<int> flags;
    flag_file->GetObject("QC_flags", flag_pointers);

    for(vector<int>::iterator it = flag_pointers->begin(); it != flag_pointers->end(); ++it) 
    {
        // cout << *it << endl;
        flags.push_back(*it);
    }
    
    double mean_rms_good_chamber_vs_impact_angle[13] = {0.0};

    double impact_angle_rms_sum[13] = {0.0};
    int impact_angle_rms_count[13] = {0};

    for (int i_det=0; i_det<540; i_det++)
    {
        bool det_suspicious = (flags[i_det] > 0);

        double point_impact_angle;
        double point_rms;
        
        if (!det_suspicious)
        {
            for (int i_point=0; i_point<13; i_point++)
            {
                vec_tgraph_rms[i_det]->GetPoint(i_point+1, point_impact_angle, point_rms);

                impact_angle_rms_sum[i_point] += point_rms;
                impact_angle_rms_count[i_point] ++;
            }
        }
        // cout << endl << "det: " << i_det << endl;
        // cout << "mean: " << mean << endl;
    }

    tg_mean_rms_good_chamber = new TGraph();

    for (int i_impact_angle=0; i_impact_angle<13; i_impact_angle++)
    {
        double impact_angle = 70 + (i_impact_angle+1) * 40/14;
        double mean_rms = impact_angle_rms_sum[i_impact_angle] / impact_angle_rms_count[i_impact_angle];

        tg_mean_rms_good_chamber->SetPoint(i_impact_angle, impact_angle, mean_rms);
    }

    // TCanvas* tg_mean_rms_good_chamber_can = new TCanvas("tg_mean_rms_good_chamber_can","tg_mean_rms_good_chamber_can",10,10,400,400);
    // tg_mean_rms_good_chamber->GetYaxis()->SetRangeUser(0,15);
    // tg_mean_rms_good_chamber->Draw("AC*");
}


void draw_ole_error_comp()
{
    TFile *ole_transformed_sigma_dy_file = TFile::Open("./ole_transformed_sigma_dy.root");
    // TCanvas *anode_hv_can = (TCanvas *)anode_hv_file->Get("tg_HV_anode_vs_det_can;1");
    TGraph* tg_ole_transformed_sigma_dy = (TGraph *)ole_transformed_sigma_dy_file->Get("Graph");
    

    TCanvas* tg_ole_error_comp_can = new TCanvas("tg_ole_error_comp_can","tg_ole_error_comp_can",10,10,400,400);
    

    tg_mean_rms_good_chamber->GetYaxis()->SetRangeUser(0,10);
    tg_mean_rms_good_chamber->GetXaxis()->SetTitle("Impact Angle [deg]");
    tg_mean_rms_good_chamber->GetYaxis()->SetTitle("#sigma_{#Delta#alpha} [deg]");
    tg_mean_rms_good_chamber->GetXaxis()->CenterTitle();
    tg_mean_rms_good_chamber->GetYaxis()->CenterTitle();
    tg_mean_rms_good_chamber->SetLineColor(kAzure+1);
    tg_mean_rms_good_chamber->SetMarkerColor(kAzure+1);
    tg_mean_rms_good_chamber->SetLineWidth(2);
    tg_mean_rms_good_chamber->SetMarkerStyle(8);

    tg_ole_transformed_sigma_dy->SetMarkerStyle(8);
    tg_ole_transformed_sigma_dy->SetMarkerColor(kTeal+1);
    tg_ole_transformed_sigma_dy->SetLineWidth(2);
    tg_ole_transformed_sigma_dy->SetLineColor(kTeal+1);

    tg_mean_rms_good_chamber->Draw("ACP");
    tg_ole_transformed_sigma_dy->Draw("CP");

    auto legend = new TLegend(0.1,0.77,0.48,0.9);
    gStyle->SetLegendTextSize(0.04);
    legend->AddEntry(tg_mean_rms_good_chamber,"Our result","l");
    legend->AddEntry(tg_ole_transformed_sigma_dy,"Ole's result","l");
    legend->Draw();
}


void h2d_error_analysis()
{
    TBase_calibrate_hist_output = TFile::Open("../Data/TRD_Calib_circle_3456_binfix.root");

    TFile *anode_hv_file = TFile::Open("../Data/HV_anode_vs_det_265338.root");
    TCanvas *anode_hv_can = (TCanvas *)anode_hv_file->Get("tg_HV_anode_vs_det_can;1");
    tg_anode_hv = (TGraph *)anode_hv_can->FindObject("tg_HV_anode_vs_det");

    TFile *drift_hv_file = TFile::Open("../Data/HV_drift_vs_det_265338.root");
    TCanvas *drift_hv_can = (TCanvas *)drift_hv_file->Get("tg_HV_drift_vs_det_can;1");
    tg_drift_hv = (TGraph *)drift_hv_can->FindObject("tg_HV_drift_vs_det");

    TFile* calibration_params = TFile::Open("../Data/TRD_Calib_vDfit_and_LAfit_3456.root");
    vdrift_fit = (TGraph*)calibration_params->Get("tg_v_fit_vs_det;1");
    
    rms_output = new TFile("./delta_alpha_rms.root", "RECREATE");


    // draw_2d_delta_alpha();


    get_vec_tgraph_rms();

    get_mean_rms_vs_impact_angle();

    // draw_vec_tgraph_rms();

    get_mean_rms_vs_layer();

    draw_ole_error_comp();
}
