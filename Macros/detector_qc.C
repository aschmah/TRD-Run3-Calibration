#include <TFile.h>
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TChain.h"

#include "TObject.h"

#include "../Ali_AS_Event.h"
#include "../Ali_AS_EventLinkDef.h"

ClassImp(Ali_AS_Event)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)

Ali_AS_Event*     AS_Event;
Ali_AS_Tracklet*  AS_Tracklet;
Ali_AS_Track*     AS_Track;
Ali_AS_TRD_digit* AS_Digit;

TList *detector_vs_adc_list;
TH1D *max_hist;
TH1D *plateu_hist;
TH1D *ratio_hist;
TH1D *all_defects_hist;
vector<int> all_defects;
vector<int> adc_defects;
vector<int> hv_anode_defects;
vector<int> hv_drift_defects;
TFile* outputfile;
TFile* calib_file;
TChain* input_SE;

TString JPsi_TREE   = "Tree_AS_Event";
TString JPsi_BRANCH = "Tree_AS_Event_branch";

Int_t offical_qa[91] = {2, 15, 17, 27, 31, 32, 36, 40, 43, 49, 50, 55, 59, 64, 88, 92, 113, 116, 119, 132, 180, 181, 190,
    191, 194, 207, 215, 219, 221, 226, 227, 228, 230, 231, 233, 236, 238, 241, 249, 255, 277, 287, 302, 308, 310, 311, 317, 318, 319, 
    320, 326, 328, 335, 348, 368, 377, 386, 389, 402, 403, 404, 405, 406, 407, 432, 433, 434, 435, 436, 437, 452, 455, 456, 462, 463, 
    464, 465, 466, 467, 470, 482, 483, 484, 485, 490, 491, 493, 500, 502, 504, 538};

Int_t bad_hv[86] = {2, 5, 8, 12, 17, 26, 27, 29, 30, 31, 32, 36, 40, 43, 49, 50, 51, 59, 64, 88, 92, 113, 116, 119, 132, 
    181, 184, 190, 191, 197, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 236, 241, 249, 255, 265, 274, 277, 287, 302, 
    308, 309, 310, 311, 316, 317, 318, 319, 320, 326, 328, 335, 348, 368, 377, 386, 389, 452, 455, 456, 461, 470, 483, 484, 485, 490, 491, 
    493, 494, 498, 500, 502, 504, 538};

Int_t offical_and_hv[111] = {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 36, 40, 43, 49, 50, 51, 55, 59, 64, 88, 92, 113, 116, 
    119, 132, 180, 181, 184, 190, 191, 194, 197, 207, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 236, 
    238, 241, 249, 255, 265, 274, 277, 287, 302, 308, 309, 310, 311, 316, 317, 318, 319, 320, 326, 328, 335, 348, 368, 377, 
    386, 389, 402, 403, 404, 405, 406, 407, 432, 433, 434, 435, 436, 437, 452, 455, 456, 461, 462, 463, 464, 465, 466, 467, 
    470, 482, 483, 484, 485, 490, 491, 493, 494, 498, 500, 502, 504, 538};


void Init_tree(TString SEList)
{
    cout << "Initialize tree" << endl;

    // TString JPsi_TREE   = "Tree_AS_Event";
    // TString JPsi_BRANCH = "Tree_AS_Event_branch";

    // TString pinputdir = "/misc/alidata120/alice_u/schmah/TRD_offline_calib/Data/";
    TString pinputdir = "/home/jasonb/Documents/masters/calibration/";
    //TString pinputdir = "/home/ceres/berdnikova/TRD-Run3-Calibration/";

    AS_Event = new Ali_AS_Event();
    // AS_Track = new Ali_AS_Track();
    AS_Tracklet = new Ali_AS_Tracklet();
    // AS_Digit = new Ali_AS_TRD_digit();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( JPsi_TREE.Data(), JPsi_TREE.Data() );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = pinputdir+addfile;
                    input_SE ->AddFile(addfile.Data(),-1, JPsi_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( JPsi_BRANCH, &AS_Event );
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
        }
    }
}


tuple<Double_t, Double_t, Double_t> getADCParams(const char *profile_name)
{
    TProfile *h = (TProfile *)detector_vs_adc_list->FindObject(profile_name);

    Double_t max;
    Double_t plateu;
    Double_t ratio;

    max = max + h->GetBinContent(2) + h->GetBinContent(3);
    max = max / 2;

    for (int i = 7; i < 21; ++i)
    {
        plateu += h->GetBinContent(i);
    }
    plateu = plateu / 14;

    ratio = max / plateu;

    //look for ADC plots with ion tail
    // if (h->GetBinContent(24) < h->GetBinContent(15) * 0.75 && plateu > 1)
    // {
    //     string detector_str(profile_name, 14, 3);
    //     cout << detector_str << ", ";
    // }

    //get total sum of ADC over 24 tb for good chambers
    // if (max > 40)
    // {
    //     Double_t ADC_sum = 0;
    //     for (Int_t i=0; i <24; i++)
    //     {
    //         ADC_sum += h->GetBinContent(i);
    //     }
    //     cout << ADC_sum << endl;
    // }

    return {max, plateu, ratio};
}


void get_ADC_defects()
{
    detector_vs_adc_list = (TList *)calib_file->Get("TRD_Digits_output;1");

    // defects = new TGraph("adc_defects_hist", "All Defects", 100, 0, 540);

    max_hist = new TH1D("max_hist", "", 200, 10, 100);
    plateu_hist = new TH1D("plateu_hist", "", 200, 0, 60);
    ratio_hist = new TH1D("ratio_hist", "", 200, 0, 8);

    TIter next(detector_vs_adc_list);
    TObject *object = 0;
    while ((object = next()))
    {
        if (strcmp(object->ClassName(), "TProfile") == 0)
        {
            auto [max, plateu, ratio] = getADCParams(object->GetName());
            string detector_str(object->GetName(), 14, 3);
            int detector = stoi(detector_str);

            bool exists = std::find(std::begin(offical_and_hv), std::end(offical_and_hv), detector) != std::end(offical_and_hv);
            if (exists)
            {
                continue;
            }

            if (max < 30 || max > 72)
            {
                adc_defects.push_back(detector);
                all_defects.push_back(detector);
                // cout << detector << endl;
                // continue;
            }

            max_hist->Fill(max);
            plateu_hist->Fill(plateu);
            ratio_hist->Fill(ratio);
        }
    }
}


void draw_ADC_hists()
{
    TLine *upper = new TLine(72,0,72,25);
    upper -> SetLineWidth(2);
    upper -> SetLineColor(kRed);

    TLine *lower = new TLine(30,0,30,25);
    lower -> SetLineWidth(2);
    lower -> SetLineColor(kRed);

    max_hist->GetXaxis()->SetTitle("ADC Value");
    max_hist->GetYaxis()->SetTitle("Number of Detectors");
    max_hist->GetXaxis()->CenterTitle();
    max_hist->GetYaxis()->CenterTitle();

    TCanvas *max_hist_can = new TCanvas("max_hist_can", "Max");
    max_hist_can->cd();
    max_hist->Draw();
    upper->Draw();
    lower->Draw();

    plateu_hist->GetXaxis()->SetTitle("ADC Value");
    plateu_hist->GetYaxis()->SetTitle("Number of Detectors");
    plateu_hist->GetXaxis()->CenterTitle();
    plateu_hist->GetYaxis()->CenterTitle();

    TCanvas *plateu_hist_can = new TCanvas("plateu_hist_can", "Plateu");
    plateu_hist_can->cd();
    plateu_hist->Draw();

    ratio_hist->GetXaxis()->SetTitle("Ratio ADC[2,3] / ADC[7,20]");
    ratio_hist->GetYaxis()->SetTitle("Number of Detectors");
    ratio_hist->GetXaxis()->CenterTitle();
    ratio_hist->GetYaxis()->CenterTitle();

    TCanvas *ratio_hist_can = new TCanvas("ratio_hist_can", "Max/Plateu");
    ratio_hist_can->cd();
    ratio_hist->Draw();
}


void get_HV_defects()
{
    cout << endl;

    TFile *hv_file = TFile::Open("../Data/HV_anode_vs_det_265338.root");
    TCanvas *anode_hv_can = (TCanvas *)hv_file->Get("tg_HV_anode_vs_det_can;1");
    TGraph *anode_hv = (TGraph *)anode_hv_can->FindObject("tg_HV_anode_vs_det");

    TCanvas * new_anode_hv = new TCanvas("new_anode_hv","new_anode_hv",0,540,1000,540);
    TLine* anode_hv_line = new TLine(0,1400,590,1400);
    anode_hv_line -> SetLineWidth(2);
    anode_hv_line -> SetLineColor(kRed);

    TGraph *anode_hv_unique = new TGraph();
    anode_hv_unique->SetMarkerColor(kRed);

    anode_hv->GetXaxis()->SetTitle("Detector");
    anode_hv->GetYaxis()->SetTitle("Anode HV [V]");
    anode_hv->GetXaxis()->CenterTitle();
    anode_hv->GetYaxis()->CenterTitle();

    for (int i = 0; i < 540; i++)
    {
        Double_t hv = anode_hv->Eval(i);
        if (hv < 1400)
        {
            hv_anode_defects.push_back(i);
            all_defects.push_back(i);
        }
    }

    for (int i =0; i<hv_anode_defects.size(); i++)
    {
        Double_t x;
        Double_t y;

        bool exists = std::find(std::begin(offical_qa), std::end(offical_qa), hv_anode_defects[i]) != std::end(offical_qa);
        if (!exists)
        {
            anode_hv->GetPoint(hv_anode_defects[i], x, y);
            cout << "x: " << x << " | " << "y: " << y << endl;
            anode_hv_unique->SetPoint(anode_hv_unique->GetN(), x, y);
        }  
    }

    // new_anode_hv->cd();
    anode_hv->Draw("AP");
    anode_hv_unique->Draw("*");
    anode_hv_line->Draw();

    for (int i = 0; i < hv_anode_defects.size(); i++)
    {
        cout << hv_anode_defects[i] << ", ";
    }



    TFile *drift_hv_file = TFile::Open("../Data/HV_drift_vs_det_265338.root");
    TCanvas *drift_hv_can = (TCanvas *)drift_hv_file->Get("tg_HV_drift_vs_det_can;1");
    TGraph *drift_hv = (TGraph *)drift_hv_can->FindObject("tg_HV_drift_vs_det");

    TCanvas * new_drift_hv = new TCanvas("new_drift_hv","new_drift_hv",0,540,1000,540);
    TLine* drift_hv_line = new TLine(0,1920,590,1920);
    drift_hv_line -> SetLineWidth(2);
    drift_hv_line -> SetLineColor(kRed);

    TGraph *drift_hv_unique = new TGraph();
    drift_hv_unique->SetMarkerColor(kRed);

    drift_hv->GetXaxis()->SetTitle("Detector");
    drift_hv->GetYaxis()->SetTitle("Drift HV [V]");
    drift_hv->GetXaxis()->CenterTitle();
    drift_hv->GetYaxis()->CenterTitle();

    for (int i = 0; i < 540; i++)
    {
        Double_t hv = drift_hv->Eval(i);
        if (hv < 1920)
        {
            hv_drift_defects.push_back(i);
            all_defects.push_back(i);
        }
    }

    for (int i =0; i<hv_drift_defects.size(); i++)
    {
        Double_t x;
        Double_t y;

        bool exists = std::find(std::begin(offical_qa), std::end(offical_qa), hv_drift_defects[i]) != std::end(offical_qa);
        if (!exists)
        {
            drift_hv->GetPoint(hv_drift_defects[i], x, y);
            cout << "x: " << x << " | " << "y: " << y << endl;
            drift_hv_unique->SetPoint(drift_hv_unique->GetN(), x, y);
        }  
    }

    new_drift_hv->cd();
    drift_hv->Draw("AP");
    drift_hv_unique->Draw("*");
    drift_hv_line->Draw();


    // for (int i = 0; i < hv_drift_defects.size(); i++)
    // {
    //     cout << hv_drift_defects[i] << ", ";
    // }

    sort( all_defects.begin(), all_defects.end() );
    all_defects.erase( unique( all_defects.begin(), all_defects.end() ), all_defects.end() );

    for (int i = 0; i < all_defects.size(); i++)
    {
        cout << all_defects[i] << ", ";
    }
}


void get_noise_defects()
{
    TTree* tree = (TTree *)calib_file->Get("Tree_AS_Event;1");
    // TBranch* branch2 = (TBranch *)tree->GetBranch("fTracklets");
    // TLeaf* leaf = (TLeaf *)tree->GetLeaf("fTracklets.detector");

    TH1D* num_tracklets_hist = new TH1D("num_tracklets_hist", "", 540, 0, 540);
    // TH1D* num_tracklets_hist = new TH1D("num_tracklets_hist", "", 100, -400, 400);

    TH2D* ADC_digits = new TH2D("ADC_digits","ADC_digits",540,0,540,70,0,1000);


    // tree->Draw("fTracklets.detector>>myh");
    // TH1F* detectors_tracklets = (TH1F *)gDirectory->Get("myh");


    //TRACKLETS
    // for (int i=0; i<540; i++)
    // {
    //     cout << detectors_tracklets->GetBinContent(500) << endl;
    //     num_tracklets_hist->Fill(detectors_tracklets->GetBinContent(i));
    // }

    // for(Long64_t i_event = 0; i_event < 3000; i_event++)
    // {
    //     if(i_event % 20 == 0) printf("i_event: %lld out of %lld \n",i_event);
    //     if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

    //     UShort_t NumTracks = AS_Event ->getNumTracklets(); // number of tracks in this event
    //     for(Int_t i_tracklet = 0; i_tracklet < NumTracks; i_tracklet++)
    //     {
    //         AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
    //         Short_t  i_det_tracklet = AS_Tracklet ->get_detector();
    //         num_tracklets_hist->Fill(i_det_tracklet);
    //     }
    // }


    //DIGITS
    for(Long64_t i_event = 0; i_event < 200; i_event++)
    {
        if(i_event % 20 == 0) printf("i_event: %lld out of %lld \n",i_event);
        if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

        UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event
        for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        {
            AS_Track              = AS_Event ->getTrack( i_track ); // take the track

            UShort_t  fNumTRDdigits        = AS_Track ->getNumTRD_digits();
            for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
            {
                AS_Digit              = AS_Track ->getTRD_digit(i_digits);

                Int_t    layer        = AS_Digit ->get_layer();
                Int_t    sector       = AS_Digit ->get_sector();
                Int_t    stack        = AS_Digit ->get_stack();
                Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);

                // bool exists = std::find(std::begin(all_defects), std::end(all_defects), detector) != std::end(all_defects);
                // if (exists)
                // {
                //     continue;
                // }

                Float_t  dca_z        = AS_Digit ->getdca_z();   
                Float_t  zpos         = AS_Digit ->get_pos(0, 2);

                Double_t ADC_sum = 0;

                for (Int_t i=0; i <24; i++)
                {
                    ADC_sum += AS_Digit->getADC_time_value(i);
                }

                num_tracklets_hist->Fill(detector, ADC_sum);
                ADC_digits->Fill(detector, ADC_sum, 1);
            }
        }
    }


    // num_tracklets_hist->GetXaxis()->SetTitle("Detector");
    // num_tracklets_hist->GetYaxis()->SetTitle("Number of Tracklets");

    num_tracklets_hist->GetXaxis()->SetTitle("Detector");
    num_tracklets_hist->GetYaxis()->SetTitle("Number of Digits");

    num_tracklets_hist->GetXaxis()->CenterTitle();
    num_tracklets_hist->GetYaxis()->CenterTitle();


    ADC_digits->GetXaxis()->SetTitle("dectector");
    ADC_digits->GetYaxis()->SetTitle("ADC sum over 24 tb");

    ADC_digits->GetXaxis()->CenterTitle();
    ADC_digits->GetYaxis()->CenterTitle();
    ADC_digits->SetStats(0);

    TCanvas *num_tracklets_hist_can = new TCanvas("num_tracklets_hist_can", "Max");
    num_tracklets_hist_can->cd();
    // num_tracklets_hist_can->SetLogz();
    num_tracklets_hist->Draw();
    // ADC_digits->Draw("colz");
    
}


void get_no_calibration()
{
    TFile* calibration_params = TFile::Open("../Data/TRD_Calib_vDfit_and_LAfit_digits.root");
    TGraph* vdrift_fit = (TGraph*)calibration_params->Get(";1");

    Double_t vdrift;
    Double_t x;
    Double_t y;    
    
    for (Int_t i_det=0; i_det<540; i_det++)
    {
        bool found = 0;

        for (Int_t vdrift_det=0; vdrift_det<540; vdrift_det++)
        {

            vdrift = vdrift_fit->GetPoint(vdrift_det, x, y);
            if ((Int_t)x == i_det)
            {
                found = 1;
            }
        }

        if (!found)
        {
            cout << i_det << ", ";
        }
    }

    int count = 0;
    TFile* h1d_detector_hit_file = TFile::Open("../Data/h_detector_hit.root");
    TH1D* h1d_detector_hit = (TH1D*)h1d_detector_hit_file->Get("h_detector_hit;1");

    // draw detector hits hist with axis labels
    // TCanvas *h1d_detector_hit_can = new TCanvas("h1d_detector_hit_can", "h1d_detector_hit_can");
    // h1d_detector_hit->GetXaxis()->SetTitle("Detector");
    // h1d_detector_hit->GetYaxis()->SetTitle("Number of Tracklets per 10k Events");

    // h1d_detector_hit->GetXaxis()->CenterTitle();
    // h1d_detector_hit->GetYaxis()->CenterTitle();

    // h1d_detector_hit_can->cd();
    // h1d_detector_hit->Draw();

    cout << endl << endl;

    // number of hits per chameber with no fit
    TH1D* h1d_hits_no_fit = new TH1D("h1d_hits_no_fit", "h1d_hits_no_fit", 56, 0, 10);

    TCanvas *h1d_hits_no_fit_can = new TCanvas("h1d_hits_no_fit_can", "h1d_hits_no_fit_can");
    h1d_hits_no_fit->GetXaxis()->SetTitle("Detector");
    h1d_hits_no_fit->GetYaxis()->SetTitle("Number of Tracklets per 10k Events");

    h1d_hits_no_fit->GetXaxis()->CenterTitle();
    h1d_hits_no_fit->GetYaxis()->CenterTitle();

    // indexed from 1
    for (Int_t i_det=0; i_det<540; i_det++)
    {
        int hits = h1d_detector_hit->GetBinContent(i_det+1);
        // cout << "idet: " << i_det << " | " << "hits: " << hits << endl;
        if (hits < 2)
        {
            count ++;
            // cout << i_det << ", ";
        }
    }
    cout << endl << count << endl;

    int det_not_in_badhv_plus_official[49] = {33, 34, 35, 37, 38, 39, 41, 48, 52, 53, 54, 56, 57, 58, 114, 115, 117, 118, 
    182, 183, 185, 234, 235, 237, 239, 306, 307, 321, 322, 323, 324, 325, 327, 329, 384, 385, 387, 388, 450, 451, 453, 454, 
    480, 481, 499, 501, 503, 533, 537};

    int det_not_in_official[56] = {30, 33, 34, 35, 37, 38, 39, 41, 48, 51, 52, 53, 54, 56, 57, 58, 114, 115, 117, 
    118, 182, 183, 184, 185, 213, 234, 235, 237, 239, 306, 307, 309, 321, 322, 323, 324, 325, 327, 329, 384, 385, 387, 388, 
    450, 451, 453, 454, 461, 480, 481, 498, 499, 501, 503, 533, 537};

    count = 0;
    for (int i_det=0; i_det<56  ; i_det++)
    {
        int hits = h1d_detector_hit->GetBinContent(det_not_in_official[i_det]+1);
        h1d_hits_no_fit->Fill(hits);
        cout << "det: " << det_not_in_official[i_det] << " | " << "hits: " << hits << endl;
        if (hits < 3) {count ++;}
    }
    cout << count << endl;

    h1d_hits_no_fit_can->cd();
    h1d_hits_no_fit->Draw();
}


void draw_all_defects_hist()
{
    all_defects_hist = new TH1D("all_defects_hist", "All Defects", 540, 0, 540);
    
    cout << all_defects.size() << endl;
    for (int i = 0; i < all_defects.size(); i++)
    {
        cout << all_defects[i] << ", ";
    }

    for (int i=0; i<540; i++)
    {
        bool exists = std::find(std::begin(all_defects), std::end(all_defects), i) != std::end(all_defects);
        if (!exists)
        {
            all_defects_hist->SetBinContent(i,1);
        }
    }
    all_defects_hist->GetXaxis()->SetTitle("Chamber");
    all_defects_hist->GetYaxis()->SetTitle("Good=1, Bad=0");
    all_defects_hist->GetXaxis()->CenterTitle();
    all_defects_hist->GetYaxis()->CenterTitle();

    TCanvas *all_defects_hist_can = new TCanvas("all_defects_hist_can", "All_defects");
    all_defects_hist_can->cd();
    all_defects_hist->Draw();

    outputfile = new TFile("./chamber_QC.root","RECREATE");
    outputfile ->cd();
    all_defects_hist->Write();
}


void write_all_defects()
{
    vector<int> v;

    for (int i=0; i<540; i++)
    {
        int flag = 0;
        if (count(adc_defects.begin(), adc_defects.end(), i))
        {
            flag += 1;
        }
	    if (count(hv_anode_defects.begin(), hv_anode_defects.end(), i))
        {
            flag += 2;
        }
        if (count(hv_drift_defects.begin(), hv_drift_defects.end(), i))
        {
            flag += 4;
        }
        v.push_back(flag);
    }

    // for (int i=0; i<v.size(); i++)
    // {
    //     cout << v[i] << endl;
    // }

    // sort( v.begin(), v.end() );
    // v.erase( unique( v.begin(), v.end() ), v.end() );


    TFile* flag_file = new TFile("./chamber_QC_flags.root","RECREATE");
    flag_file->WriteObject(&v,"QC_flags");
}


void detector_qc()
{
    calib_file = TFile::Open("../../Merge_TRD_Calib_V4.root");

    // TFile* hist = TFile::Open("./chamber_QC.root");
    // TH1D *histo = (TH1D *)hist->Get("all_defects_hist;1");

    
    // Init_tree("../list_tree_off_trkl_V1.txt");

    // get_ADC_defects();

    // draw_ADC_hists();

    // get_HV_defects();

    // get_noise_defects();

    // get_no_calibration();

    // draw_all_defects_hist();

    // write_all_defects();
}
