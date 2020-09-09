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
TGraph *tg_vdrift_fit;
TGraph *tg_la_fit;

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

    // alpha or offset
    TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_impact_circle;1");
    // TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_offset");
    delta_impact_circle_dir->cd();

    Double_t x;
    Double_t y;

    int i_det = 0;
    int object_counter = 0;
    for (auto&& keyAsObj : *delta_impact_circle_dir->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        if (key->ReadObj()->IsA()->InheritsFrom("TH2")) 
        {
            vec_TH2D_Delta_vs_impact_circle[i_det] = (TH2D *)key->ReadObj();
            // vec_TH2D_Delta_vs_impact_circle[i_det]->RebinY(4);
            // object_counter ++;
            // if (object_counter % 2 == 0)
            // {
            //     i_det ++;
            // }
            i_det ++;
        }
    }

    vector<TCanvas*> vec_can_2D_Delta_vs_impact_circle;
    vec_can_2D_Delta_vs_impact_circle.resize(30);

    int rows = 3;
    int columns = 6;
    int n_pads = rows*columns;
    for(Int_t i_canvas = 10; i_canvas < 20; i_canvas++)
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

            // rebin(10)
            // Double_t total_offset = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,5,55);
            // rebin(2)
            // Double_t total_offset = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,100,200);
            // rebin(4)
            // Double_t total_offset = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,50,100);
            // no rebin
            // Double_t total_offset = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,200,400);

            // alpha or offset
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

            // if (total_offset < 1)
            // {
            //     draw_not = 1;
            // }

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
                tg_vdrift_fit->GetPoint(vdrift_det, x, y);
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
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetNdivisions(505, 0);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->CenterTitle();
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->CenterTitle();
            // vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetNdivisions(4);
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetTitle("impact angle [deg]");

            // alpha or offset
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetTitle("  #Delta#alpha [deg]");
            // vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetTitle("#Deltay [mm]");
            vec_TH2D_Delta_vs_impact_circle[i_det]->GetXaxis()->SetRangeUser(70,110);

            // alpha or offset
            if (draw_outer)
            {
                vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetRangeUser(-10,90);
            }
            else
            {
                vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetRangeUser(-25,25);
                // vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetRangeUser(-20,21);
                // vec_TH2D_Delta_vs_impact_circle[i_det]->GetYaxis()->SetRangeUser(-10,10);
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
        // vec_can_2D_Delta_vs_impact_circle[i_canvas]->Print(Form("canvas%d_%d_to_%d.png", i_canvas, i_canvas*n_pads, i_canvas*n_pads + n_pads));
    }
}


Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
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

    // alpha or offset
    TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_impact_circle;1");
    // TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_offset");
    delta_impact_circle_dir->cd();

    Double_t x;
    Double_t y;

    int i_det = 0;
    int object_counter = 0;
    for (auto&& keyAsObj : *delta_impact_circle_dir->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        if (key->ReadObj()->IsA()->InheritsFrom("TH2")) 
        {
            
            // weird bug in root file where each TH2D is duplicated once
            // if (object_counter < 2)
            // {
            //     object_counter ++;
            //     continue;
            // }

            h2d_delta_alpha_det = (TH2D *)key->ReadObj();
            // h2d_delta_alpha_det->RebinY(4);

            cout << endl << key->GetName() << endl;

            // 168 - 173
            // if (i_det < 6)
            // {
            //     i_det ++;
            //     continue;
            // }

            for (int i_bin=1; i_bin<14; i_bin++)
            {
                Double_t center = h2d_delta_alpha_det->GetXaxis()->GetBinCenter(i_bin);

                y_projection_bin = (TH1D *)h2d_delta_alpha_det->ProjectionY("y_projection_bin", i_bin, i_bin); 
                int n_bins = y_projection_bin->GetEntries();

                // for (int bin=75; bin<125; bin++)
                // {
                //     cout << "bin: " << bin << " | " << "fills: " << y_projection_bin->GetBinContent(bin) << endl;;
                // }


                // TCanvas* gauss_fit_can = new TCanvas("gauss_fit_can","gauss_fit_can",10,10,1600,1000);
                // gauss_fit_can->cd();

                // alpha or offset
                TF1* func_Gauss_fit = new TF1("func_Gauss_fit",GaussFitFunc,-50,50,3);
                // TF1* func_Gauss_fit = new TF1("func_Gauss_fit",GaussFitFunc,-25,25,3);

                // y_projection_bin->Draw();
                // func_Gauss_fit->Draw();

                Double_t rms = y_projection_bin->GetRMS();

                double amplitude = 20;
                double mean = 0;
                double sigma = rms;

                for(Int_t i = 0; i < 3; i++)
                {
                    func_Gauss_fit ->ReleaseParameter(i);
                    func_Gauss_fit ->SetParameter(i,0.0);
                    func_Gauss_fit ->SetParError(i,0.0);
                }

                func_Gauss_fit ->SetParameter(0,amplitude);
                func_Gauss_fit ->SetParameter(1,mean);
                func_Gauss_fit ->SetParameter(2,sigma);

                // alpha or offset
                y_projection_bin->GetXaxis()->SetRangeUser(-50,50);
                // y_projection_bin->GetXaxis()->SetRangeUser(-15,15);

                //Fit
                y_projection_bin ->Fit("func_Gauss_fit","QMN","",mean-1.5*sigma,mean+1.5*sigma);

                amplitude = func_Gauss_fit ->GetParameter(0);
                mean      = func_Gauss_fit ->GetParameter(1);
                sigma     = fabs(func_Gauss_fit ->GetParameter(2));
                
                double chi2 = func_Gauss_fit->GetChisquare();

                double err = func_Gauss_fit->GetParError(2);

                cout << err << endl;

                // y_projection_bin->Draw();
                // func_Gauss_fit->Draw("same");

                for(Int_t i = 0; i < 3; i++)
                {
                    func_Gauss_fit ->ReleaseParameter(i);
                    func_Gauss_fit ->SetParameter(i,0.0);
                    func_Gauss_fit ->SetParError(i,0.0);
                }

                func_Gauss_fit ->SetParameter(0,amplitude);
                func_Gauss_fit ->SetParameter(1,mean);
                func_Gauss_fit ->SetParameter(2,sigma);

                //fitten
                y_projection_bin ->Fit("func_Gauss_fit","QMN","",mean-1.5*sigma,mean+1.5*sigma);

                // y_projection_bin ->GetXaxis()->SetRangeUser(mean-4*sigma,mean+4*sigma);

                // Double_t bin_cent_max_bin = vec_ADC_spec_cat[det][i_cat][i_run_ids] ->GetBinCenter(vec_ADC_spec_cat[det][i_cat][i_run_ids] ->GetMaximumBin());

                //Parameter auslesen
                amplitude = func_Gauss_fit ->GetParameter(0);
                mean      = func_Gauss_fit ->GetParameter(1);
                // mean_err  = func_Gauss_fit ->GetParError(1);
                sigma     = fabs(func_Gauss_fit ->GetParameter(2));


                // Double_t rms = y_projection_bin->GetRMS();
 
                if (n_bins < 100 || err > 1.5)
                {
                    sigma = 0;
                }

                cout << "BIN: " << i_bin << " | " << "rms: " << rms << " | " << "sigma: " << sigma << endl;

                // cout << center << endl;
                // rms or gaussian sigma
                vec_tgraph_rms[i_det]->SetPoint(i_bin, center, sigma);
                // if (n_bins > 20)
                // {
                //     vec_tgraph_rms[i_det]->SetPoint(i_bin, center, sigma);
                // }

                // gauss_fit_can->Print(Form("det%d_bin%d_rebin%d.png", i_det, i_bin, 1));
                // break;

            }
            // object_counter ++;
            // if (object_counter % 2 == 0)
            // {
            //     i_det ++;
            // }
            i_det ++;
            // cout << i_det << " | " << vec_tgraph_rms[i_det]->GetN() << endl;
            // if (i_det > 173)
            // {
            //     break;
            // }
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
    for(Int_t i_canvas = 9; i_canvas < 10; i_canvas++)
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
                tg_vdrift_fit->GetPoint(vdrift_det, x, y);
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

            // alpha or offset
            vec_tgraph_rms[i_det]->GetYaxis()->SetTitle("#sigma_{#Delta#alpha} [deg]");
            // vec_tgraph_rms[i_det]->GetYaxis()->SetTitle("#sigma_{#Deltay} [mm]");
            vec_tgraph_rms[i_det]->GetXaxis()->SetRangeUser(70,110);

            // alpha or offset
            vec_tgraph_rms[i_det]->GetYaxis()->SetRangeUser(0,24);
            // vec_tgraph_rms[i_det]->GetYaxis()->SetRangeUser(0,5);
            
            
            bool is_defective = find(begin(official_qa), end(official_qa), i_det) != end(official_qa);
            if (is_defective)
            {
                cout << i_det << endl;
                vec_tgraph_rms[i_det]->GetYaxis()->SetAxisColor(kRed);
                vec_tgraph_rms[i_det]->GetXaxis()->SetAxisColor(kRed);
                // vec_tgraph_rms[i_det]->SetPoint(90, 90, 15);
                // vec_tgraph_rms[i_det]->Draw("AXIS");
                // tg_mean_rms_good_chamber->Draw("AC");
            }

            // if (vec_tgraph_rms[i_det]->GetN() == 0)
            // {
            //     vec_tgraph_rms[i_det]->Draw("AXIS");
            // }
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
        // if (flags[i_det] == 8)
        if (!det_suspicious)
        {
            // if (point_rms < 20 && point_rms > 0.1)
            {
            layer_count[layer_number] ++;
            layer_sum[layer_number] += mean;
            }
        }
    }


    for (int i_layer=0; i_layer<6; i_layer++) 
    {
        // cout << layer_count[i_layer] << endl;
        if (layer_count[i_layer] == 0) {continue;}
        double layer_mean = layer_sum[i_layer] / layer_count[i_layer];
    
        cout << layer_mean << endl;

        tg_rms_vs_layer->SetPoint(i_layer, i_layer, layer_mean);
    }

    TCanvas* can_tg_rms_vs_layer = new TCanvas("can_tg_rms_vs_layer","can_tg_rms_vs_layer",10,10,400,400);
    can_tg_rms_vs_layer->SetRightMargin(0.01);
    can_tg_rms_vs_layer->SetTopMargin(0.1);
    can_tg_rms_vs_layer->SetBottomMargin(0.15);
    can_tg_rms_vs_layer->SetLeftMargin(0.15);

    tg_rms_vs_layer->SetTitle("Drift HV OR Drift HV + ADC");
    tg_rms_vs_layer->SetMarkerColor(kBlue);

    tg_rms_vs_layer->GetXaxis()->SetTitle("Layer");
    tg_rms_vs_layer->GetYaxis()->SetTitle("Mean RMS");
    tg_rms_vs_layer->GetXaxis()->CenterTitle();
    tg_rms_vs_layer->GetYaxis()->CenterTitle();
    tg_rms_vs_layer->GetXaxis()->SetLabelSize(0.05);
    tg_rms_vs_layer->GetYaxis()->SetLabelSize(0.05);
    tg_rms_vs_layer->GetXaxis()->SetTitleSize(0.05);
    tg_rms_vs_layer->GetYaxis()->SetTitleSize(0.05);
    // tg_rms_vs_layer->GetYaxis()->SetRangeUser(0,15);
    tg_rms_vs_layer->SetFillColor(kBlue);
    tg_rms_vs_layer->Draw("AB");
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
            // cout << i_det << " | " << vec_tgraph_rms[i_det]->GetN() << endl;
            for (int i_point=0; i_point<13; i_point++)
            {
                vec_tgraph_rms[i_det]->GetPoint(i_point+1, point_impact_angle, point_rms);
                if (point_rms < 1000 && point_rms > 0.01)
                {
                    impact_angle_rms_sum[i_point] += point_rms;
                    impact_angle_rms_count[i_point] ++;
                }
            }
        }
        // cout << endl << "det: " << i_det << endl;
        // cout << "mean: " << mean << endl;
        // break;
    }

    tg_mean_rms_good_chamber = new TGraph();

    for (int i_impact_angle=0; i_impact_angle<13; i_impact_angle++)
    {
        double impact_angle = 70 + (i_impact_angle+1) * 40/14;
        // cout << impact_angle_rms_count[i_impact_angle] << endl;
        double mean_rms = impact_angle_rms_sum[i_impact_angle] / impact_angle_rms_count[i_impact_angle];

        cout << impact_angle << " | " << mean_rms << endl;
        // cout << mean_rms << "+";
        
        tg_mean_rms_good_chamber->SetPoint(i_impact_angle, impact_angle, mean_rms);
    }
    // cout << endl;

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


void draw_2d_average_offset()
{
    TDirectory *delta_impact_circle_dir = (TDirectory *)TBase_calibrate_hist_output->Get("Delta_offset");
    delta_impact_circle_dir->cd();

    int rebin = 4;
    // int bin = 1;

    TH2D* TH2D_delta_offset_sum = new TH2D("TH2D_delta_offset_sum", "TH2D_delta_offset_sum",13,70,110,600/rebin,-30,30);

    Double_t x;
    Double_t y;

    int i_det = 0;
    for (auto&& keyAsObj : *delta_impact_circle_dir->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        if (key->ReadObj()->IsA()->InheritsFrom("TH2")) 
        {
            if (i_det % 6 == 0 || i_det % 6 == 5)
            {
                i_det ++;
                continue;
            }
            TH2D* h = (TH2D *)key->ReadObj();
            h->RebinY(rebin);
            TH2D_delta_offset_sum->Add(h);

            i_det ++;
        }
    }

    TCanvas* TH2D_delta_offset_sum_can = new TCanvas("TH2D_delta_offset_sum_can","TH2D_delta_offset_sum_can",10,10,1000,800);
    TH2D_delta_offset_sum_can->cd();
    TH2D_delta_offset_sum->GetXaxis()->CenterTitle();
    TH2D_delta_offset_sum->GetYaxis()->CenterTitle();
    TH2D_delta_offset_sum->GetXaxis()->SetTitle("impact angle [deg]");
    TH2D_delta_offset_sum->GetYaxis()->SetTitle("#Deltay [mm]");
    TH2D_delta_offset_sum->GetXaxis()->SetRangeUser(70,110);
    TH2D_delta_offset_sum->GetYaxis()->SetRangeUser(-11,11);
    TH2D_delta_offset_sum_can->SetRightMargin(0.2);
    gStyle->SetStatX(0.3);
    gStyle->SetStatY(0.9);                
    gStyle->SetStatW(0.2);                
    gStyle->SetStatH(0.1);  
    // TH2D_delta_offset_sum_can->SetLogZ();
    gPad->SetLogz(1);
    TH2D_delta_offset_sum->SetTitle("");
    TH2D_delta_offset_sum->Draw("colz");

    for (int i_bin=1; i_bin<14; i_bin++)
    {
        TH1D* y_projection_bin = (TH1D *)TH2D_delta_offset_sum->ProjectionY("y_projection_bin", i_bin, i_bin); 
        int n_bins = y_projection_bin->GetEntries();

        // TCanvas* gauss_fit_can = new TCanvas("gauss_fit_can","gauss_fit_can",10,10,800,600);
        // gauss_fit_can->cd();

        TF1* func_Gauss_fit = new TF1("func_Gauss_fit",GaussFitFunc,-5,5,3);


        Double_t rms = y_projection_bin->GetRMS();

        double amplitude = 3000;
        double mean = 0;
        double sigma = rms;

        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);

        y_projection_bin->GetXaxis()->SetRangeUser(-5,5);

        //Fit
        y_projection_bin ->Fit("func_Gauss_fit","QMN","",mean-2.*sigma,mean+2.*sigma);

        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));


        // y_projection_bin->Draw();
        // func_Gauss_fit->Draw("same");

        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);

        //fitten
        y_projection_bin ->Fit("func_Gauss_fit","QMN","",mean-1.2*sigma,mean+1.2*sigma);

        //Parameter auslesen
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        // mean_err  = func_Gauss_fit ->GetParError(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));


        // Double_t rms = y_projection_bin->GetRMS();

        if (n_bins < 20)
        {
            sigma = 1000;
        }

        cout << "BIN: " << i_bin << " | " << "rms: " << rms << " | " << "sigma: " << sigma << endl;

        // gauss_fit_can->Print(Form("bin%d_rebin%d.png", i_bin, rebin));
        // break;
    }
    

    // y_projection_bin->Draw();
}


void systematic_fit_uncertainty()
{
    TH1D* h1d_vd_fits = new TH1D("vd_fits", "Drift Velocity", 100, 0.6, 1.6);
    TH1D* h1d_la_fits = new TH1D("la_fits", "Lorentz Angle", 100, -0.3, 0.1);

    TCanvas *vd_fit_can = new TCanvas("vd_fits","vd_fits",10,10,800,600);
    TCanvas *la_fit_can = new TCanvas("la_fits","la_fits",10,10,800,600);

    double anode_det;
    double anode_hv;
    double drift_det;
    double drift_hv;

    double vd_det;
    double vd_fit;
    double la_det;
    double la_fit;


    for (int ipoint = 0; ipoint < tg_drift_hv->GetN(); ipoint++)
    {
        // tg_anode_hv->GetPoint(ipoint, anode_det, anode_hv);
        tg_drift_hv->GetPoint(ipoint, drift_det, drift_hv);
        
        // if (anode_det == drift_det)
        {
            // cout << "anode_hv: " << anode_hv << " | " << "drift_hv: " << drift_hv << endl;
            // if ((1500 < anode_hv && anode_hv < 1600) && (1950 < drift_hv && drift_hv < 2250))
            if (2150 < drift_hv && drift_hv < 2170)
            {
                for (int ipoint_drift = 0; ipoint_drift < tg_vdrift_fit->GetN(); ipoint_drift++)
                {
                    tg_vdrift_fit->GetPoint(ipoint_drift, vd_det, vd_fit);
                    if (vd_det == drift_det)
                    {
                        // cout << vd_det <<  endl;
                        h1d_vd_fits->Fill(vd_fit);
                    }
                }

                for (int ipoint_la = 0; ipoint_la < tg_la_fit->GetN(); ipoint_la++)
                {
                    tg_la_fit->GetPoint(ipoint_la, la_det, la_fit);
                    if (la_det == drift_det)
                    {
                        // cout << vd_det <<  endl;
                        h1d_la_fits->Fill(la_fit);
                    }
                }
            }
        }
    }

    h1d_vd_fits->GetXaxis()->CenterTitle();
    h1d_vd_fits->GetYaxis()->CenterTitle();
    h1d_vd_fits->GetXaxis()->SetTitle("vD fit");
    h1d_vd_fits->GetYaxis()->SetTitle("number of chambers");

    h1d_la_fits->GetXaxis()->CenterTitle();
    h1d_la_fits->GetYaxis()->CenterTitle();
    h1d_la_fits->GetXaxis()->SetTitle("LA fit [rad.]");
    h1d_la_fits->GetYaxis()->SetTitle("number of chambers");

    vd_fit_can->cd();
    h1d_vd_fits->Draw();

    la_fit_can->cd();
    h1d_la_fits->Draw();
}


void plot_vd_vs_drifthv()
{
    TGraph* tg_vd_vs_drifthv = new TGraph();
    TCanvas* tg_vd_vs_drifthv_can = new TCanvas("tg_vd_vs_drifthv","tg_vd_vs_drifthv",10,10,800,600);

    double drift_det;
    double drift_hv;

    double vd_det;
    double vd_fit;

    for (int ipoint = 0; ipoint < tg_drift_hv->GetN(); ipoint++)
    {
        tg_drift_hv->GetPoint(ipoint, drift_det, drift_hv);

        for (int ipoint2 = 0; ipoint2 < tg_vdrift_fit->GetN(); ipoint2++)
        {
            tg_vdrift_fit->GetPoint(ipoint2, vd_det, vd_fit);
            
            if (vd_det == drift_det)
            {
                tg_vd_vs_drifthv->SetPoint(tg_vd_vs_drifthv->GetN(), drift_hv, vd_fit);
            }
        }
    }

    tg_vd_vs_drifthv->GetXaxis()->CenterTitle();
    tg_vd_vs_drifthv->GetYaxis()->CenterTitle();
    tg_vd_vs_drifthv->GetXaxis()->SetTitle("Drift HV [V]");
    tg_vd_vs_drifthv->GetYaxis()->SetTitle("vD fit");

    tg_vd_vs_drifthv_can->cd();
    tg_vd_vs_drifthv->Draw("A*");
        
}


void h2d_error_analysis()
{
    // TBase_calibrate_hist_output = TFile::Open("../Data/TRD_Calib_circle_3456_minos_postfix_5plus_aminlt4.root");
    // TBase_calibrate_hist_output = TFile::Open("../Data/TRD_Calib_circle_3456_minos_postfix_5plus_aminlt1.root");
    TBase_calibrate_hist_output = TFile::Open("../Data/TRD_Calib_circle_3456_minos_postfix_3plus_aminlt1_50k.root");
    // TBase_calibrate_hist_output = TFile::Open("../Data/TRD_Calib_circle_3456_minos_ymath_50k.root");

    TFile *anode_hv_file = TFile::Open("../Data/HV_anode_vs_det_265338.root");
    TCanvas *anode_hv_can = (TCanvas *)anode_hv_file->Get("tg_HV_anode_vs_det_can;1");
    tg_anode_hv = (TGraph *)anode_hv_can->FindObject("tg_HV_anode_vs_det");

    TFile *drift_hv_file = TFile::Open("../Data/HV_drift_vs_det_265338.root");
    TCanvas *drift_hv_can = (TCanvas *)drift_hv_file->Get("tg_HV_drift_vs_det_can;1");
    tg_drift_hv = (TGraph *)drift_hv_can->FindObject("tg_HV_drift_vs_det");

    // TFile* calibration_params = TFile::Open("../Data/TRD_Calib_vDfit_and_LAfit_3456.root");
    TFile* calibration_params = TFile::Open("../Data/TRD_Calib_vDfit_and_LAfit_3456_minos_5plus_aminlt1_50k.root");
    tg_vdrift_fit = (TGraph*)calibration_params->Get("tg_v_fit_vs_det;1");
    tg_la_fit = (TGraph*)calibration_params->Get("tg_LA_factor_fit_vs_det;1");
    
    rms_output = new TFile("./delta_alpha_rms.root", "RECREATE");


    // draw_2d_delta_alpha();


    // get_vec_tgraph_rms();

    // get_mean_rms_vs_impact_angle();

    // draw_vec_tgraph_rms();

    // get_mean_rms_vs_layer();

    // draw_ole_error_comp();

    // draw_2d_average_offset();

    // systematic_fit_uncertainty();

    plot_vd_vs_drifthv();
}
