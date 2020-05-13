#include <TFile.h>

TList *detector_vs_adc_list;
TH1D *max_hist;
TH1D *plateu_hist;
TH1D *ratio_hist;
TH1D *all_defects_hist;
vector<int> all_defects;
TFile* outputfile;


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

    return {max, plateu, ratio};
}


void draw_defects_hist()
{
    all_defects_hist = new TH1D("all_defects_hist", "All Defects", 540, 0, 539);
    
    cout << all_defects.size() << endl;
    // for (int i = 0; i < all_defects.size(); i++)
    // {
    //     cout << all_defects[i] << ",";
    // }

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
    // for(Int_t i = 0; i < 540; i++)
    // {
    //     all_defects_hist[i]          ->Write();
    // }
}


void draw_hists()
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

    ratio_hist->GetXaxis()->SetTitle("");
    ratio_hist->GetYaxis()->SetTitle("Number of Detectors");
    ratio_hist->GetXaxis()->CenterTitle();
    ratio_hist->GetYaxis()->CenterTitle();

    TCanvas *ratio_hist_can = new TCanvas("ratio_hist_can", "Max/Plateu");
    ratio_hist_can->cd();
    ratio_hist->Draw();
}


void detector_qc()
{

    TFile *calib_file = TFile::Open("../../Merge_TRD_Calib_V4.root");
    detector_vs_adc_list = (TList *)calib_file->Get("TRD_Digits_output;1");

    // defects = new TGraph("defectives_hist", "All Defects", 100, 0, 540);


    max_hist = new TH1D("max_hist", "Average ADC Over Timebins [2,3]", 200, 10, 100);
    plateu_hist = new TH1D("plateu_hist", "Average ADC Over Timebins [7,20]", 200, 0, 60);
    ratio_hist = new TH1D("ratio_hist", "Ratio of Max/Plateu", 200, 0, 8);

    vector<int> defectives;

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
                defectives.push_back(detector);
                all_defects.push_back(detector);
                continue;
            }

            max_hist->Fill(max);
            plateu_hist->Fill(plateu);
            ratio_hist->Fill(ratio);
        }
    }

    // draw_hists();


    //HV stuff

    cout << endl;

    vector<int> hv_anode_defectives;
    vector<int> hv_drift_defectives;

    TFile *hv_file = TFile::Open("../Data/HV_anode_vs_det_265338.root");
    TCanvas *anode_hv_can = (TCanvas *)hv_file->Get("tg_HV_anode_vs_det_can;1");
    TGraph *anode_hv = (TGraph *)anode_hv_can->FindObject("tg_HV_anode_vs_det");

    for (int i = 0; i < 540; i++)
    {
        Double_t hv = anode_hv->Eval(i);
        if (hv < 1400)
        {
            hv_anode_defectives.push_back(i);
            all_defects.push_back(i);
        }
    }

    // for (int i = 0; i < hv_anode_defectives.size(); i++)
    // {
    //     cout << hv_anode_defectives[i] << ",";
    // }

    TFile *drift_hv_file = TFile::Open("../Data/HV_drift_vs_det_2016.root");
    TCanvas *drift_hv_can = (TCanvas *)drift_hv_file->Get("tg_HV_drift_vs_det_can;1");
    TGraph *drift_hv = (TGraph *)drift_hv_can->FindObject("tg_HV_drift_vs_det");

    for (int i = 0; i < 540; i++)
    {
        Double_t hv = drift_hv->Eval(i);
        if (hv < 1920)
        {
            hv_drift_defectives.push_back(i);
            all_defects.push_back(i);
        }
    }

    // for (int i = 0; i < hv_drift_defectives.size(); i++)
    // {
    //     cout << hv_drift_defectives[i] << ",";
    // }

    draw_defects_hist();
}
