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

    //get sum of ADC over 24 tb for good chambers
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

            if (max < 40)
            {
                adc_defects.push_back(detector);
                all_defects.push_back(detector);
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
    max_hist->GetXaxis()->SetTitle("ADC Value");
    max_hist->GetYaxis()->SetTitle("Number of Detectors");
    max_hist->GetXaxis()->CenterTitle();
    max_hist->GetYaxis()->CenterTitle();

    TCanvas *max_hist_can = new TCanvas("max_hist_can", "Max");
    max_hist_can->cd();
    max_hist->Draw();

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

    for (int i = 0; i < 540; i++)
    {
        Double_t hv = anode_hv->Eval(i);
        if (hv < 1400)
        {
            hv_anode_defects.push_back(i);
            all_defects.push_back(i);
        }
    }

    // for (int i = 0; i < hv_anode_defects.size(); i++)
    // {
    //     cout << hv_anode_defects[i] << ",";
    // }

    TFile *drift_hv_file = TFile::Open("../Data/HV_drift_vs_det_2016.root");
    TCanvas *drift_hv_can = (TCanvas *)drift_hv_file->Get("tg_HV_drift_vs_det_can;1");
    TGraph *drift_hv = (TGraph *)drift_hv_can->FindObject("tg_HV_drift_vs_det");

    for (int i = 0; i < 540; i++)
    {
        Double_t hv = drift_hv->Eval(i);
        if (hv < 1920)
        {
            hv_drift_defects.push_back(i);
            all_defects.push_back(i);
        }
    }

    // for (int i = 0; i < hv_drift_defects.size(); i++)
    // {
    //     cout << hv_drift_defects[i] << ",";
    // }


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

void draw_all_defects_hist()
{
    all_defects_hist = new TH1D("all_defects_hist", "All Defects", 540, 0, 540);
    
    cout << all_defects.size() << endl;
    for (int i = 0; i < all_defects.size(); i++)
    {
        cout << all_defects[i] << ",";
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

    // cout << histo->GetBinContent(0) << endl;
    // cout << histo->GetBinContent(1) << endl;
    // cout << histo->GetBinContent(2) << endl;
    // cout << histo->GetBinContent(3) << endl;
    // cout << histo->GetBinContent(4) << endl;
    // cout << histo->GetBinContent(542) << endl;

    // TFile* test = TFile::Open("./chamber_QC_flags.root");
    // vector<int> *t;
    // test->GetObject("QC_flags", t);

    // for(vector<int>::iterator it = t->begin(); it != t->end(); ++it) {
    //   cout << *it << endl;
    // }

    // Init_tree("../list_tree_off_trkl_V1.txt");

    get_ADC_defects();

    // draw_ADC_hists();

    get_HV_defects();

    // get_noise_defects();

    draw_all_defects_hist();

    // write_all_defects();
}
