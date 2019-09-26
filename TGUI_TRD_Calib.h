


#include "Ana_Digits_functions.h"

#include "Ali_AS_Event.h"
#include "TBase_TRD_Calib.h"
#include "Ali_AS_EventLinkDef.h"

#include "TBase_TRD_Calib.h"
#include "TBase_TRD_CalibLinkDef.h"
ClassImp(TBase_TRD_Calib)


ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Event)


//---------------------------------------------------------------------------------
//class TGUI_TRD_Calib : public TGMainFrame, public TBase_TRD_Calib
class TGUI_TRD_Calib : public TGMainFrame
{
private:
    TRootEmbeddedCanvas *fCanvas_HV_vs_time        = NULL;
    TGMainFrame* Frame_Main;
    TGHorizontalFrame *hframe_Main[4];
    TGVerticalFrame   *vframe_Main[4];
    TGVerticalFrame   *vframe_stat_Main[4];

    TGNumberEntry*     arr_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_stat[3];

    TGTextButton *Button_exit;
    TGTextButton *Button_load;
    TGTextButton *Button_save;
    TGTextButton *Button_draw3D;
    TGTextButton *Button_draw3D_track;

    TBase_TRD_Calib *Base_TRD_Calib;
    vector<TPolyMarker3D*> vec_TPM3D_digits;
    TPolyMarker3D* TPM3D_cluster;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;

    vector< vector< vector<Float_t> > > vec_digit_track_info; // [track number][digit number][info ->] x,y,z,time,ADC,sector,stack,layer,row,column,dca
    vector< vector<Float_t> >           vec_track_info; // [track number][info ->] dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p

    TCanvas* c_3D        = NULL;
    TCanvas* c_3D_track  = NULL;

    Int_t color_layer[6] = {kBlack,kGreen,kBlue,kMagenta,kCyan,kYellow};
    vector<TPolyMarker3D*> vec_TPM3D_single_track_digit_layer;
    vector<TPolyMarker3D*> vec_TPM3D_single_track_digit;
    TPolyMarker3D* TPM3D_single;

public:
    TGUI_TRD_Calib();
    virtual ~TGUI_TRD_Calib();
    Int_t LoadData();
    Int_t Draw3D();
    Int_t Draw3D_track();
    ClassDef(TGUI_TRD_Calib, 0)
};
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
TGUI_TRD_Calib::TGUI_TRD_Calib() : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    //-------------------------------------
    cout << "TGUI_TRD_Calib started" << endl;
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    SetCleanup(kDeepCleanup);
    TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
    //-------------------------------------



    //-------------------------------------
    vec_TPM3D_single_track_digit_layer.resize(6); // layer
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_single_track_digit_layer[i_layer] = new TPolyMarker3D();
    }
    TPM3D_single = new TPolyMarker3D();
    TPM3D_cluster = new TPolyMarker3D();
    //-------------------------------------



    //-------------------------------------
    // Create horizontal splitter
    Frame_Main = new TGMainFrame(gClient->GetRoot(), 400, 100);
    Frame_Main ->SetWindowName("Buttons");

    //--------------
    // A horizontal frame
    hframe_Main[0]  = new TGHorizontalFrame(Frame_Main,200,100);

    // exit button
    Button_exit = new TGTextButton(hframe_Main[0], "&Exit ","gApplication->Terminate(0)");
    hframe_Main[0]->AddFrame(Button_exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // load button
    Button_load = new TGTextButton(hframe_Main[0], "&Load ",10);
    Button_load->Connect("Clicked()", "TGUI_TRD_Calib", this, "LoadData()");
    hframe_Main[0]->AddFrame(Button_load, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_save = new TGTextButton(hframe_Main[0], "&Save ","gApplication->Terminate(0)");
    hframe_Main[0]->AddFrame(Button_save, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw3D = new TGTextButton(hframe_Main[0], "&Draw 3D ",10);
    Button_draw3D->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw3D()");
    hframe_Main[0]->AddFrame(Button_draw3D, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw3D_track = new TGTextButton(hframe_Main[0], "&Draw 3D track ",10);
    Button_draw3D_track->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw3D_track()");
    hframe_Main[0]->AddFrame(Button_draw3D_track, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Frame_Main ->AddFrame(hframe_Main[0], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    TString arr_label_stat[3] = {"# events:          ","# tracks:       ","# digits:           "};
    hframe_Main[1]  = new TGHorizontalFrame(Frame_Main,200,100);
    for(Int_t i_param = 0; i_param < 3; i_param++)
    {
        vframe_stat_Main[i_param] = new TGVerticalFrame(hframe_Main[1], 400,200);
        TString label_entry = arr_label_stat[i_param];
        arr_Label_NEntry_stat[i_param] = new TGLabel(vframe_stat_Main[i_param], label_entry.Data(), myGC(), myfont->GetFontStruct());
        vframe_stat_Main[i_param]->AddFrame(arr_Label_NEntry_stat[i_param], new TGLayoutHints(kLHintsCenterX,20,20,2,2));
        hframe_Main[1]->AddFrame(vframe_stat_Main[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    }
    Frame_Main ->AddFrame(hframe_Main[1], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    TString arr_label_params_ana[4] = {"Event","Track","blubb B","blubb C"};
    hframe_Main[2]  = new TGHorizontalFrame(Frame_Main,200,100);
    for(Int_t i_param = 0; i_param < 4; i_param++)
    {
        vframe_Main[i_param] = new TGVerticalFrame(hframe_Main[2], 200,200);
        arr_NEntry_ana_params[i_param] = new TGNumberEntry(vframe_Main[i_param], 0.0, 12,(TGNumberFormat::EStyle) 0);
        arr_NEntry_ana_params[i_param] ->SetNumStyle( TGNumberFormat::kNESInteger); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
        vframe_Main[i_param]->AddFrame(arr_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TString label_entry = arr_label_params_ana[i_param];
        arr_Label_NEntry_ana_params[i_param] = new TGLabel(vframe_Main[i_param], label_entry.Data(), myGC(), myfont->GetFontStruct());
        vframe_Main[i_param]->AddFrame(arr_Label_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        hframe_Main[2]->AddFrame(vframe_Main[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    }
    Frame_Main ->AddFrame(hframe_Main[2], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    Frame_Main ->Resize(450,300); // size of frame
    Frame_Main ->MapSubwindows();
    Frame_Main ->MapWindow();
    Frame_Main ->Move(1350,850); // position of frame
    //-------------------------------------

    Base_TRD_Calib = new TBase_TRD_Calib();
    Base_TRD_Calib ->Init_tree("list_tree.txt");
    Base_TRD_Calib ->Loop_event(0);
    vec_TPM3D_digits = Base_TRD_Calib ->get_PM3D_digits();

    LoadData();

}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
TGUI_TRD_Calib::~TGUI_TRD_Calib()
{
    // Clean up
    Cleanup();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::LoadData()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);

    Int_t Event = arr_NEntry_ana_params[0]->GetNumberEntry()->GetNumber();
    Base_TRD_Calib ->Loop_event(Event);
    N_Events = Base_TRD_Calib ->get_N_Events();
    N_Tracks = Base_TRD_Calib ->get_N_Tracks();
    N_Digits = Base_TRD_Calib ->get_N_Digits();
    vec_TPM3D_digits = Base_TRD_Calib ->get_PM3D_digits();

    arr_Label_NEntry_stat[0] ->SetText(Form("# events: %lld",N_Events));
    arr_Label_NEntry_stat[1] ->SetText(Form("# tracks: %lld",N_Tracks));
    arr_Label_NEntry_stat[2] ->SetText(Form("# digits: %lld",N_Digits));

   vec_track_info       = Base_TRD_Calib ->get_track_info();
   vec_digit_track_info = Base_TRD_Calib ->get_digit_track_info();

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw3D()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw3D->ChangeBackground(green);

    //if(!c_3D) c_3D    = new TCanvas("c_3D","c_3D",10,10,800,800);

    Base_TRD_Calib ->Draw_TRD();

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_digits[i_layer] ->SetMarkerColor(color_layer[i_layer]);
        vec_TPM3D_digits[i_layer] ->SetMarkerSize(5.0);
        vec_TPM3D_digits[i_layer] ->SetMarkerStyle(20);
        vec_TPM3D_digits[i_layer] ->DrawClone("");
    }

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw3D_track()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw3D_track->ChangeBackground(green);


    //if(!c_3D_track) c_3D_track    = new TCanvas("c_3D_track","c_3D_track",10,10,800,800);

    Int_t i_track = arr_NEntry_ana_params[1]->GetNumberEntry()->GetNumber();

    Base_TRD_Calib ->Draw_TRD();
    Base_TRD_Calib ->Draw_track(i_track);
    Base_TRD_Calib ->Draw_neighbor_tracks(i_track);

    vector<Int_t> vec_merge_time_bins;
    vec_merge_time_bins.resize(4);
    vec_merge_time_bins[0] = 0;
    vec_merge_time_bins[1] = 5;
    vec_merge_time_bins[2] = 12;
    vec_merge_time_bins[3] = 23;

    Base_TRD_Calib ->set_merged_time_bins(vec_merge_time_bins);
    vector< vector<TVector3> > vec_TV3_digit_pos_cluster = Base_TRD_Calib ->make_clusters(i_track); // layer, merged time bin

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_TV3_digit_pos_cluster[i_layer].size(); i_time_merge++)
        {
            TPM3D_cluster ->SetNextPoint(vec_TV3_digit_pos_cluster[i_layer][i_time_merge][0],vec_TV3_digit_pos_cluster[i_layer][i_time_merge][1],vec_TV3_digit_pos_cluster[i_layer][i_time_merge][2]);
        }
    }


    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_single_track_digit_layer[i_layer] ->Clear();
    }

    vec_TPM3D_single_track_digit.clear();


    Float_t dca            = vec_track_info[i_track][0];
    Float_t TPCdEdx        = vec_track_info[i_track][1];
    Float_t momentum       = vec_track_info[i_track][2];
    Float_t eta_track      = vec_track_info[i_track][3];
    Float_t pT_track       = vec_track_info[i_track][4];
    Float_t TOFsignal      = vec_track_info[i_track][5];
    Float_t Track_length   = vec_track_info[i_track][6];
    Float_t TRDsumADC      = vec_track_info[i_track][7];
    Float_t TRD_signal     = vec_track_info[i_track][8];
    Float_t nsigma_TPC_e   = vec_track_info[i_track][9];
    Float_t nsigma_TPC_pi  = vec_track_info[i_track][10];
    Float_t nsigma_TPC_p   = vec_track_info[i_track][11];

    Int_t n_digits_track = (Int_t)vec_digit_track_info[i_track].size();

    printf("TGUI_TRD_Calib::Draw3D_track(), i_track: %d, n_digits_track: %d \n",i_track,n_digits_track);

    for(Int_t i_digit = 0; i_digit < n_digits_track; i_digit++)
    {
        Float_t x_pos   = vec_digit_track_info[i_track][i_digit][0];
        Float_t y_pos   = vec_digit_track_info[i_track][i_digit][1];
        Float_t z_pos   = vec_digit_track_info[i_track][i_digit][2];
        Int_t time_bin  = (Int_t)vec_digit_track_info[i_track][i_digit][3];
        Float_t raw_ADC = vec_digit_track_info[i_track][i_digit][4];
        Int_t sector    = (Int_t)vec_digit_track_info[i_track][i_digit][5];
        Int_t stack     = (Int_t)vec_digit_track_info[i_track][i_digit][6];
        Int_t layer     = (Int_t)vec_digit_track_info[i_track][i_digit][7];
        Int_t row       = (Int_t)vec_digit_track_info[i_track][i_digit][8];
        Int_t column    = (Int_t)vec_digit_track_info[i_track][i_digit][9];
        Float_t dca     = vec_digit_track_info[i_track][i_digit][10];

        Double_t digit_size = raw_ADC/50.0;

        //if(!(row == 14 && column == 120)) continue;
        TPM3D_single ->SetNextPoint(x_pos,y_pos,z_pos);
        vec_TPM3D_single_track_digit_layer[layer]   ->SetNextPoint(x_pos,y_pos,z_pos);
        vec_TPM3D_single_track_digit.push_back((TPolyMarker3D*)TPM3D_single ->Clone());
        Int_t N_digits = (Int_t)vec_TPM3D_single_track_digit.size();
        vec_TPM3D_single_track_digit[N_digits - 1] ->SetMarkerSize(1.0);
        vec_TPM3D_single_track_digit[N_digits - 1] ->SetMarkerColor(kGreen+2);
        vec_TPM3D_single_track_digit[N_digits - 1] ->SetMarkerStyle(20);

        printf("  -> i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, row: %d, column: %d, time_bin: %d, sector: %d, stack: %d, layer: %d, raw_ADC: %4.3f \n",i_digit,x_pos,y_pos,z_pos,row,column,time_bin,sector,stack,layer,raw_ADC);
    }


    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_single_track_digit_layer[i_layer] ->SetMarkerColor(color_layer[i_layer]);
        vec_TPM3D_single_track_digit_layer[i_layer] ->SetMarkerSize(1.0);
        vec_TPM3D_single_track_digit_layer[i_layer] ->SetMarkerStyle(20);
        vec_TPM3D_single_track_digit_layer[i_layer] ->DrawClone("");
    }

    TPM3D_cluster ->SetMarkerColor(kRed);
    TPM3D_cluster ->SetMarkerSize(1.5);
    TPM3D_cluster ->SetMarkerStyle(20);
    TPM3D_cluster ->DrawClone("");


#if 0
    for(Int_t i_digit = 0; i_digit < (Int_t)vec_TPM3D_single_track_digit.size(); i_digit++)
    {
        printf("i_digit: %d \n",i_digit);
        vec_TPM3D_single_track_digit[i_digit] ->DrawClone("");
    }
#endif


    return 1;
}
//---------------------------------------------------------------------------------


