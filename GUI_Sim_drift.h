
static Int_t i_detector_global = 0;
static vector<TProfile*> vec_tp_Delta_vs_impact;
static const Double_t l_drift = 0.0335; // m
static const Double_t B_field = 0.5; // T = kg * s^-2 * A^-1
static const Double_t delta_t = 0.1; // mus
static const Int_t    N_time  = 24;



//----------------------------------------------------------------------------------------
TGraph* calc_Delta_alpha(Double_t B_field_use, Double_t E_field, Double_t v_drift_use, Double_t vD_use)
{
    TGraph* tg_delta_vs_angle_single = new TGraph();

    Double_t v_drift = 1.0*v_drift_use*1E06/100.0;
    Double_t alpha_L = TMath::ATan(v_drift*B_field_use/E_field);
    Double_t v_drift_vD = 1.0*vD_use*1E06/100.0; // m/s

    Int_t phi_counter = 0;
    for(Double_t track_phi_deg = 135.0; track_phi_deg >= 45.0; track_phi_deg -= 1.0)
    {
        Double_t track_phi   = track_phi_deg*TMath::DegToRad();
        Double_t track_slope = TMath::Tan(track_phi);

        vector< vector<Double_t> > vec_xy_pos_vD;
        vec_xy_pos_vD.resize(2); // start, stop
        vec_xy_pos_vD[0].resize(2); // x, y
        vec_xy_pos_vD[1].resize(2); // x, y

        Int_t flag_first_time = 0;
        for(Int_t i_time = 0; i_time < N_time; i_time++)
        {
            Double_t time  = i_time*delta_t*1E-06; // s
            Double_t y_pos = time*v_drift;     // original track cluster positions
            Double_t x_pos = y_pos/track_slope;

            for(Double_t i_drift_time = 0; i_drift_time < 0.00003; i_drift_time += 0.00000003)
            {
                Double_t x_pos_drift = x_pos;  // drifted cluster position without Lorentz-angle effect
                Double_t y_pos_drift = y_pos + v_drift*i_drift_time;

                Double_t x_pos_drift_Lorentz = x_pos - TMath::Tan(alpha_L) * v_drift*i_drift_time; // drifted cluster position with Lorentz-angle effect
                Double_t y_pos_drift_Lorentz = y_pos + v_drift*i_drift_time;

                Double_t x_pos_Lorentz    = x_pos_drift_Lorentz; // Lorentz-angle shifted cluster position
                Double_t y_pos_Lorentz    = y_pos_drift_Lorentz - v_drift*i_drift_time; // same as y_pos

                Double_t x_pos_Lorentz_vD = x_pos_drift_Lorentz; // Lorentz-angle shifted cluster position with wrong drift velocity vD
                Double_t y_pos_Lorentz_vD = y_pos_drift_Lorentz - v_drift_vD*i_drift_time;

                //printf("v_drift: %4.3f, v_drift_vD: %4.3f \n",v_drift,v_drift_vD);

                if(y_pos_drift > l_drift && y_pos_drift > 0.0) // at anonde wire
                {
                    if(i_drift_time == 0) break;
                    if(!flag_first_time)
                    {
                        vec_xy_pos_vD[0][0] = x_pos_Lorentz_vD;
                        vec_xy_pos_vD[0][1] = y_pos_Lorentz_vD;
                        flag_first_time = 1;
                    }
                    else
                    {
                        //if(y_pos_Lorentz_vD >= 0.0 && y_pos_Lorentz_vD <= l_drift)
                        {
                            vec_xy_pos_vD[1][0] = x_pos_Lorentz_vD;
                            vec_xy_pos_vD[1][1] = y_pos_Lorentz_vD;
                        }
                    }
                    break;
                }
            } // end of drift time loop
        } // end of cluster time loop


        Double_t track_phi_vD  = TMath::ATan((vec_xy_pos_vD[1][1] - vec_xy_pos_vD[0][1])/(vec_xy_pos_vD[1][0] - vec_xy_pos_vD[0][0]));
        if(track_phi_vD < 0.0) track_phi_vD += TMath::Pi();
        Double_t delta_phi_deg = -(track_phi - track_phi_vD)*TMath::RadToDeg();

        tg_delta_vs_angle_single ->SetPoint(phi_counter,track_phi_deg,delta_phi_deg);
        phi_counter++;
    } // end of phi loop

    return tg_delta_vs_angle_single;
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
// function to be minimized
void Chi2_TRD_vDrift(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t )
{
    sum = 0;

    Double_t B_field_use = par[0];
    Double_t E_field_use = par[1];
    Double_t v_drift_use = par[2];
    Double_t vD_use      = par[3];

    TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(B_field_use,E_field_use,v_drift_use,vD_use);

    for(Int_t i_bin = 1; i_bin <= vec_tp_Delta_vs_impact[i_detector_global]->GetNbinsX(); i_bin++)
    {
        Double_t impact_angle     = vec_tp_Delta_vs_impact[i_detector_global] ->GetBinCenter(i_bin);
        Double_t Delta_alpha      = vec_tp_Delta_vs_impact[i_detector_global] ->GetBinContent(i_bin);
        Double_t Delta_alpha_err  = vec_tp_Delta_vs_impact[i_detector_global] ->GetBinError(i_bin);
        Delta_alpha_err = 1.0;

        if(Delta_alpha_err == 0.0) continue;
        if(impact_angle < 78.0) continue;
        if(impact_angle > 100.0) continue;

        Double_t Delta_alpha_sim = tg_Delta_vs_impact_single ->Eval(impact_angle);

        sum += TMath::Power(Delta_alpha_sim - Delta_alpha,2)/TMath::Power(Delta_alpha_err,2);
        //sum += TMath::Power(Delta_alpha_sim - Delta_alpha,2);

    }

    delete tg_Delta_vs_impact_single;
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    // Text alignment: https://root.cern.ch/doc/master/classTAttText.html#T1
    // align = 10*HorizontalAlign + VerticalAlign
    // horizontal: 1=left adjusted, 2=centered, 3=right adjusted
    // vertical: 1=bottom adjusted, 2=centered, 3=top adjusted


    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLine* PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;

    return Zero_line;
}
//----------------------------------------------------------------------------------------


//#include "Ana_Digits_functions.h"




//---------------------------------------------------------------------------------
class GUI_Sim_drift : public TGMainFrame
{
private:

    TString HistName;
    char NoP[50];

    TGraph* tg_vdrift_vs_det;
    TGraph* tg_HV_drift_vs_det;
    TGraph* tg_v_fit_vs_det;
    TGraph* tg_vD_fit_vs_det;

    TCanvas* can_vdrift;


    TPolyMarker* tpm_track;
    TPolyMarker* tpm_track_Lorentz;
    TPolyMarker* tpm_track_Lorentz_vD;
    TPolyMarker* tpm_track_Lorentz_vD_used;

    vector<TPolyMarker*> tpm_track_drift;
    vector<TPolyMarker*> tpm_track_drift_Lorentz;

    vector< vector<TGraph*> > tg_delta_vs_angle;

    Int_t phi_counter_plot   = 48;
    Int_t vD_counter_plot    = 0;
    Int_t v_counter_plot     = 0;

    vector< Double_t> vec_v_fit;
    vector< Double_t> vec_vD_fit;

    Double_t v_fit  = 1.56;
    Double_t vD_fit = 1.56;
    Double_t E_fit  = 2000.0/0.0335;


    TRootEmbeddedCanvas *fCanvas_HV_vs_time        = NULL;
    TGMainFrame* Frame_Main;
    TGHorizontalFrame *hframe_Main[4];
    TGVerticalFrame   *vframe_Main[4];
    TGVerticalFrame   *vframe_stat_Main[4];

    TGNumberEntry*     arr_NEntry_det;
    TGLabel*           TGLabel_det;


    TGGroupFrame*  fGroupFrames[7];
    TGCheckButton* fCheckBox_sel[3];


    TGNumberEntry*     arr_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_stat[3];

    TGTextButton *Button_exit;
    TGTextButton *Button_load;
    TGTextButton *Button_save;
    TGTextButton *Button_minimize;
    TGTextButton *Button_draw_data;
    TGTextButton *Button_draw3D_track;
    TGTextButton *Button_Calibrate;

    vector<TGHorizontalFrame*>  vec_Hframe;
    vector<TGHSlider*>          vec_slider;
    vector<TGLayoutHints*>      vec_LayoutHints;
    vector<TGTextEntry*>        vec_TextEntry;
    vector<TGTextBuffer*>       vec_TextBuffer;
    vector<TGLabel*>            vec_TGLabel;

    vector<TPolyMarker3D*> vec_TPM3D_digits;
    vector<TPolyLine3D*>   vec_TPL3D_tracklets;
    TPolyMarker3D* TPM3D_cluster;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;

    vector< vector< vector<Float_t> > > vec_digit_track_info; // [track number][digit number][info ->] x,y,z,time,ADC,sector,stack,layer,row,column,dca
    vector< vector<Float_t> >           vec_track_info; // [track number][info ->] dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p

    TCanvas* c_3D        = NULL;
    TCanvas* c_3D_track  = NULL;

    Int_t color_layer[6] = {kRed,kGreen,kBlue,kMagenta,kCyan,kYellow};
    vector<TPolyMarker3D*> vec_TPM3D_single_track_digit_layer;
    vector<TPolyMarker3D*> vec_TPM3D_single_track_digit;
    TPolyMarker3D* TPM3D_single;

public:
    GUI_Sim_drift();
    virtual ~GUI_Sim_drift();
    void DoText(const char *text);
    Int_t LoadData();
    Int_t Do_Minimize();
    Int_t Draw_data();
    Int_t Draw3D_track();
    Int_t Calibrate();
    ClassDef(GUI_Sim_drift, 0)
};
//---------------------------------------------------------------------------------




//---------------------------------------------------------------------------------
GUI_Sim_drift::GUI_Sim_drift() : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    //-------------------------------------
    cout << "GUI_Sim_drift started" << endl;
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    SetCleanup(kDeepCleanup);
    TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
    //-------------------------------------

    vec_tp_Delta_vs_impact.resize(540);
    tpm_track_drift.resize(N_time);
    tpm_track_drift_Lorentz.resize(N_time);

    tpm_track                 = new TPolyMarker();
    tpm_track_Lorentz         = new TPolyMarker();
    tpm_track_Lorentz_vD      = new TPolyMarker();
    tpm_track_Lorentz_vD_used = new TPolyMarker();

    vec_Hframe.resize(5);
    vec_slider.resize(5);
    vec_LayoutHints.resize(5);
    vec_TextEntry.resize(5);
    vec_TextBuffer.resize(5);
    vec_TGLabel.resize(5);


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
    Button_load->Connect("Clicked()", "GUI_Sim_drift", this, "LoadData()");
    hframe_Main[0]->AddFrame(Button_load, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_save = new TGTextButton(hframe_Main[0], "&Save ","gApplication->Terminate(0)");
    hframe_Main[0]->AddFrame(Button_save, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // minimize button
    Button_minimize = new TGTextButton(hframe_Main[0], "&Minimize ",10);
    Button_minimize->Connect("Clicked()", "GUI_Sim_drift", this, "Do_Minimize()");
    hframe_Main[0]->AddFrame(Button_minimize, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    // draw 3D button
    Button_draw_data = new TGTextButton(hframe_Main[0], "&Draw data ",10);
    Button_draw_data->Connect("Clicked()", "GUI_Sim_drift", this, "Draw_data()");
    hframe_Main[0]->AddFrame(Button_draw_data, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    arr_NEntry_det = new TGNumberEntry(hframe_Main[0], 0.0, 12,(TGNumberFormat::EStyle) 0);
    arr_NEntry_det ->SetNumStyle( TGNumberFormat::kNESInteger); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
    arr_NEntry_det ->Connect("ValueSet(Long_t)", "GUI_Sim_drift", this, "Draw_data()");
    hframe_Main[0] ->AddFrame(arr_NEntry_det, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

    TGLabel_det = new TGLabel(hframe_Main[0], "det", myGC(), myfont->GetFontStruct());
    hframe_Main[0] ->AddFrame(TGLabel_det, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

    Frame_Main ->AddFrame(hframe_Main[0], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //------------------------------------------------------------



    //------------------------------------------------------------
    fGroupFrames[0] = new TGGroupFrame(Frame_Main, new TGString("Input"),kHorizontalFrame|kRaisedFrame);
    TString label_checkbox[3] = {"Use slider","Use fit","C"};
    for(Int_t i_cb = 0; i_cb < 3; i_cb++)
    {
        fCheckBox_sel[i_cb]  = new TGCheckButton(fGroupFrames[0], new TGHotString(label_checkbox[i_cb].Data()), -1);
        //fCheckBox_sel[i_cb] ->SetState(kButtonDown);
        fCheckBox_sel[i_cb] ->SetState(kButtonUp);
        fCheckBox_sel[i_cb] ->Connect("Clicked()", "GUI_Sim_drift", this, "Draw_data()");
        fGroupFrames[0]->AddFrame(fCheckBox_sel[i_cb], new TGLayoutHints(kLHintsTop | kLHintsLeft,0, 0, 5, 0));
    }
    Frame_Main ->AddFrame(fGroupFrames[0], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //------------------------------------------------------------



    //------------------------------------------------------------
    // Add a slider
    Int_t start_pos_slider[5] = {9,15,9,0,0};
    TString vec_label[5] = {"v","vD","UD","Rx","fboost"};
    for(Int_t i_param = 0; i_param < 3; i_param++)
    {
        vec_Hframe[i_param] = new TGHorizontalFrame(Frame_Main, 200, 200);
        vec_slider[i_param] = new TGHSlider(vec_Hframe[i_param],250,kSlider1|kScaleDownRight,1);
        vec_slider[i_param]->Connect("PositionChanged(Int_t)", "GUI_Sim_drift",this, "Draw_data()");
        vec_slider[i_param]->SetRange(0,20);
        vec_slider[i_param]->SetPosition(start_pos_slider[i_param]);
        vec_LayoutHints[i_param] = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 1, 1, 1, 1); // handles size of slider and number box
        vec_Hframe[i_param]->AddFrame(vec_slider[i_param], vec_LayoutHints[i_param]);
        Frame_Main->AddFrame(vec_Hframe[i_param], vec_LayoutHints[i_param]);


        // Add number text box
        vec_TextEntry[i_param] = new TGTextEntry(vec_Hframe[i_param], vec_TextBuffer[i_param] = new TGTextBuffer(5));
        //vec_TextEntry[i_param] ->SetDefaultSize(0,0);
        vec_TextEntry[i_param] ->SetTextColor(kRed);
        vec_TextEntry[i_param]->SetToolTipText("Minimum (left) Value of Slider");
        vec_TextBuffer[i_param]->AddText(0, "0.0");
        //vec_TextEntry[i_param]->Connect("TextChanged(char*)", "GUI_Sim_drift", this,"DoText(char*)");
        //vec_Hframe[i_param]->Resize(100, 25);
        vec_Hframe[i_param]->AddFrame(vec_TextEntry[i_param], vec_LayoutHints[i_param]);
        Frame_Main->AddFrame(vec_Hframe[i_param], vec_LayoutHints[i_param]);
        vec_TGLabel[i_param] = new TGLabel(vec_Hframe[i_param], vec_label[i_param].Data());
        vec_Hframe[i_param]->AddFrame(vec_TGLabel[i_param], new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

        vec_Hframe[i_param] ->SetHeight(1);
        //vec_Hframe[i_param] ->DrawBorder();

    }
    //----------------------------------------------------



    Frame_Main ->Resize(550,400); // size of frame
    Frame_Main ->MapSubwindows();
    Frame_Main ->MapWindow();
    Frame_Main ->Move(850,750); // position of frame
    //-------------------------------------


    LoadData();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
GUI_Sim_drift::~GUI_Sim_drift()
{
    // Clean up
    Cleanup();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t GUI_Sim_drift::LoadData()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);

    //---------------------------------------------------------
    // Read the data
    //TFile* input_data = TFile::Open("./Data/TRD_Calib.root");
    TFile* input_data = TFile::Open("./Data/TRD_Calib_Zero.root");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det] = (TProfile*)input_data->Get(Form("vec_th1d_Delta_vs_impact_%d",i_det));
    }

    TFile* input_vdrift = TFile::Open("./Data/vdrift_vs_det_2016_B.root");
    //TCanvas* can_vdrift = (TCanvas*)input_vdrift->Get("tg_vdrift_vs_det_can");
    can_vdrift = (TCanvas*)input_vdrift->Get("tg_vdrift_vs_det_can");
    can_vdrift ->Draw();
    tg_vdrift_vs_det = (TGraph*)can_vdrift->FindObject("tg_vdrift_vs_det");

    TFile* input_HVdrift = TFile::Open("./Data/HV_drift_vs_det_2016_B.root");
    TCanvas* can_HVdrift = (TCanvas*)input_HVdrift->Get("tg_HV_drift_vs_det_can");
    can_HVdrift ->Draw();
    tg_HV_drift_vs_det = (TGraph*)can_HVdrift->FindObject("tg_HV_drift_vs_det");
    //---------------------------------------------------------

    return 1;
}
//---------------------------------------------------------------------------------

#if 1

//single chamber - if we need it later
//---------------------------------------------------------------------------------
Int_t GUI_Sim_drift::Do_Minimize()
{

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_minimize->ChangeBackground(green);

    Int_t i_detector = arr_NEntry_det->GetNumberEntry()->GetNumber();

    Double_t det = -1.0;
    Double_t v_drift_in = -1.0;
    tg_vdrift_vs_det ->GetPoint(i_detector,det,v_drift_in);
    Double_t HV_drift_in = -1.0;
    tg_HV_drift_vs_det ->GetPoint(i_detector,det,HV_drift_in);
    Double_t E_field = 1.0*(HV_drift_in/l_drift); // V/cm = kg * m * s^-3 * A^-1

    //-------------------------------------------------
    // Minimization
    i_detector_global = i_detector;
    TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
    min->SetFCN(Chi2_TRD_vDrift);
    Double_t pStart[4] = {B_field,HV_drift_in/l_drift,v_drift_in,1.56}; // B-field, E-field, v_drift, vD_drift (1.7 instead of 1.56)
    //Double_t pStart[4] = {B_field,HV_drift_in/l_drift,1.48,1.56}; // B-field, E-field, v_drift, vD_drift (1.7 instead of 1.56)
   
    // 73: 1.48  after fit: 1.269

    min->SetParameter(0,"B_field",pStart[0],0.01,0,0);
    min->SetParameter(1,"E_field",pStart[1]*0.6,0.01,0,0);
    min->SetParameter(2,"v_drift",pStart[2],0.01,0,0);
    min->SetParameter(3,"vD_drift",pStart[3],0.01,0,0);

    min->FixParameter(0);
    min->FixParameter(1);
    min->FixParameter(3);

    Double_t arglist[2];
    arglist[0] = 1000; // number of function calls
    arglist[1] = 0.001; // tolerance
    min->ExecuteCommand("MIGRAD",arglist,2);
    // get fit parameters
    Double_t parFit[4];
    for(int i = 0; i < 4; ++i)
    {
        parFit[i] = min->GetParameter(i);
    }

    v_fit  = parFit[2];
    vD_fit = parFit[3];
    E_fit  = parFit[1];

    printf("v_drift: %4.3f, vD_drift: %4.3f, E_field: %4.3f \n",parFit[2],parFit[3],parFit[1]);
    //-------------------------------------------------

    Draw_data();
    return 1;
}

#endif
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
// minimize all chambers
//---------------------------------------------------------------------------------
#if 0
Int_t GUI_Sim_drift::Do_Minimize()
{

    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_minimize->ChangeBackground(green);

    //Int_t i_detector = arr_NEntry_det->GetNumberEntry()->GetNumber();

    Double_t det;
    Double_t v_drift_in;
    Double_t HV_drift_in;
    Double_t E_field;

    printf("test 1 \n");

    vec_v_fit.resize(540);
    vec_vD_fit.resize(540);

    TGraph* tg_v_fit_vs_det  = new TGraph();
    TGraph* tg_vD_fit_vs_det = new TGraph();

    for (Int_t i_detector = 0; i_detector < 540; ++i_detector)
    {
        det         = -1.0;
        v_drift_in  = -1.0;
        HV_drift_in = -1.0;
        E_field     = -1.0;

        tg_vdrift_vs_det ->GetPoint(i_detector,det,v_drift_in);
        tg_HV_drift_vs_det ->GetPoint(i_detector,det,HV_drift_in);
        E_field = 1.0*(HV_drift_in/l_drift); // V/cm = kg * m * s^-3 * A^-1

        vec_v_fit[i_detector] = 1.56;
        vec_vD_fit[i_detector] = 1.56;

        printf("test 2 \n");
        printf("i_detector = %d \n",i_detector);

//#if 1
        //-------------------------------------------------
        // Minimization
        i_detector_global = i_detector;

        Int_t N = vec_tp_Delta_vs_impact[i_detector_global]->GetEntries();
        printf("N_points = %d \n",N);

        if (N >= 10 && E_field > 0)

        {

            TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
            min->SetFCN(Chi2_TRD_vDrift);
            Double_t pStart[4] = {B_field,HV_drift_in/l_drift,v_drift_in,1.56}; // B-field, E-field, v_drift, vD_drift (1.7 insread of 1.56)
            //Double_t pStart[4] = {B_field,HV_drift_in/l_drift,v_drift_in-0.1,1.56}; // B-field, E-field, v_drift, vD_drift (1.7 insread of 1.56)
           
            min->SetParameter(0,"B_field",pStart[0],0.01,0,0);
            min->SetParameter(1,"E_field",pStart[1],0.01,0,0);
            min->SetParameter(2,"v_drift",pStart[2],0.01,0,0);
            min->SetParameter(3,"vD_drift",pStart[3],0.01,0,0);

            printf("test 2.1 \n");

            min->FixParameter(0);
            min->FixParameter(1);
            min->FixParameter(3);

            Double_t arglist[2];
            arglist[0] = 1000; // number of function calls
            arglist[1] = 0.001; // tolerance

            printf("test 2.2  \n");

            //min->ExecuteCommand("MINIMIZE",arglist,2);
            min->ExecuteCommand("MIGRAD",arglist,2);

            printf("test 2.3  \n");
            // get fit parameters
            Double_t parFit[4];

            printf("test 2.4  \n");
            for(int i = 0; i < 4; ++i)
            {
                parFit[i] = min->GetParameter(i);
            }

            //min->SetParameter(0,"B_field",parFit[0],0.01,0,0);
            //min->SetParameter(1,"E_field",parFit[1],0.01,0,0);
            //min->SetParameter(2,"v_drift",parFit[2],0.01,0,0);
            //min->SetParameter(3,"vD_drift",parFit[3],0.01,0,0);

            //printf("test 2.4 new SetParam  \n");

            //min->ExecuteCommand("MIGRAD",arglist,2);

            //for(int i = 0; i < 4; ++i)
            //{
            //    parFit[i] = min->GetParameter(i);
            //}

            vec_v_fit[i_detector]  = parFit[2];
            vec_vD_fit[i_detector] = parFit[3];
            E_fit                  = parFit[1];

            printf("v_drift: %4.3f, vD_drift: %4.3f, E_field: %4.3f \n",parFit[2],parFit[3],parFit[1]);
            //-------------------------------------------------
            //#endif
            printf("test 3 \n");

            

            delete min;
            printf("test 4 \n");
            
        }

        tg_v_fit_vs_det  ->SetPoint(i_detector,i_detector,vec_v_fit[i_detector]);
        tg_vD_fit_vs_det ->SetPoint(i_detector,i_detector,vec_vD_fit[i_detector]);

        //tg_vD_fit_vs_det ->SetPoint(i_detector,i_detector,i_detector);

        printf("i_detector = %d \n",i_detector);
        printf("tg_v_fit_vs_det -> getbincont: %4.3f \n",tg_v_fit_vs_det ->Eval(i_detector));
        printf("tg_vD_fit_vs_det -> getbincont: %4.3f \n",tg_vD_fit_vs_det ->Eval(i_detector));
        //printf("v_drift: %4.3f, vD_drift: %4.3f, E_field: %4.3f \n",parFit[2],parFit[3],parFit[1]);


    }

    printf("test 6 \n");

    can_vdrift       ->cd();

    tg_v_fit_vs_det ->SetMarkerColor(kRed);
    tg_v_fit_vs_det ->SetLineColor(kRed);
    tg_v_fit_vs_det ->SetMarkerSize(1.0);
    tg_v_fit_vs_det ->Draw();
    printf("test 7 \n");

    //can_vdrift       ->Modify();
    can_vdrift       ->Update();
    printf("test 8 \n");

    //Draw_data();

    return 1;

}
#endif
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t GUI_Sim_drift::Draw_data()
{

    Double_t v_slider   = 0.1   + 0.1*vec_slider[0]->GetPosition();
    Double_t vD_slider  = 0.06  + 0.1*vec_slider[1]->GetPosition();
    Double_t UD_slider  = 500.0 + 100*vec_slider[2]->GetPosition();

    Double_t arr_param_val[5] = {v_slider,vD_slider,UD_slider,0.0,0.0};

    // Handle slider widgets.
    for(Int_t i_param = 0; i_param < 3; i_param++)
    {
        char buf[32];
        //sprintf(buf, "%d", vec_slider[i_param]->GetPosition());
        sprintf(buf,"%2.4f",arr_param_val[i_param]);
        vec_TextBuffer[i_param]->Clear();
        vec_TextBuffer[i_param]->AddText(0, buf);
        vec_TextEntry[i_param]->SetCursorPosition(vec_TextEntry[i_param]->GetCursorPosition());
        vec_TextEntry[i_param]->Deselect();
        gClient->NeedRedraw(vec_TextEntry[i_param]);
    }

    if(tpm_track)
    {
        delete tpm_track;
        tpm_track                 = new TPolyMarker();
    }
    if(tpm_track_Lorentz)
    {
        delete tpm_track_Lorentz;
        tpm_track_Lorentz         = new TPolyMarker();
    }
    if(tpm_track_Lorentz_vD)
    {
        delete tpm_track_Lorentz_vD;
        tpm_track_Lorentz_vD      = new TPolyMarker();
    }
    if(tpm_track_Lorentz_vD_used)
    {
        delete tpm_track_Lorentz_vD_used;
        tpm_track_Lorentz_vD_used = new TPolyMarker();
    }

    Int_t i_detector = arr_NEntry_det->GetNumberEntry()->GetNumber();
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw_data->ChangeBackground(green);

    Double_t det = -1.0;
    Double_t v_drift_in = -1.0;
    tg_vdrift_vs_det ->GetPoint(i_detector,det,v_drift_in);
    Double_t HV_drift_in = -1.0;
    tg_HV_drift_vs_det ->GetPoint(i_detector,det,HV_drift_in);
    Double_t E_field = 1.0*(HV_drift_in/l_drift); // V/cm = kg * m * s^-3 * A^-1

    printf("detector: %d, v_drift: %4.3f, HV_drift: %4.3f \n",i_detector,v_drift_in,HV_drift_in);

    Double_t vD_use = 1.56;
    Double_t v_drift_use = v_drift_in;
    if(fCheckBox_sel[0]->GetState() == kButtonDown) // slider
    {
        printf("Use slider data \n");
        v_drift_use = v_slider;
        vD_use      = vD_slider;
        E_field     = 1.0*(UD_slider/l_drift); // V/cm = kg * m * s^-3 * A^-1;
    }

    if(fCheckBox_sel[1]->GetState() == kButtonDown) // fit
    {
        printf("Use fit data \n");
        v_drift_use = v_fit;
        vD_use      = vD_fit;
        E_field     = E_fit;
    }


    //------------------------------------------------------------------------
    Int_t v_counter = 0;
    vector<Double_t> vec_v_drift_vals;
    vec_v_drift_vals.clear();
    tg_delta_vs_angle.clear();
    //for(Double_t fac_v_drift = 1.2; fac_v_drift >= 0.8; fac_v_drift -= 0.05) // true drift velocity
    for(Double_t fac_v_drift = 1.0; fac_v_drift >= 1.0; fac_v_drift -= 0.05) // true drift velocity
    {
        //Double_t v_drift = 1.0*1.65*1E06/100.0*fac_v_drift; // m/s
        //Double_t v_drift = 1.0*0.612*1E06/100.0*fac_v_drift; // m/s
        //Double_t v_drift = 0.7*1E06/100.0*fac_v_drift; // m/s
        Double_t v_drift = 1.0*v_drift_use*1E06/100.0;
        Double_t alpha_L = TMath::ATan(v_drift*B_field/E_field);
        //printf("alpha_L: %4.3f, v_drift: %4.6f, ref: %4.6f \n",alpha_L*TMath::RadToDeg(),v_drift,1.56*1E06/100.0);
        tg_delta_vs_angle.resize(v_counter+1);
        Int_t vD_counter = 0;
        //for(Double_t fac_v_drift_vD = 1.1; fac_v_drift_vD >= 0.74; fac_v_drift_vD -= 0.05)
        for(Double_t fac_v_drift_vD = 1.0; fac_v_drift_vD >= 1.0; fac_v_drift_vD -= 0.05)
        {
            //Double_t v_drift_vD = v_drift*fac_v_drift_vD*0.8;
            Double_t v_drift_vD = 1.0*vD_use*1E06/100.0; // m/s
            if(v_counter == 0) vec_v_drift_vals.push_back(v_drift_vD);

            printf("v_drift: %4.6f, v_drift_vD: %4.6f, E_field: %4.3f \n",v_drift,v_drift_vD,E_field);

            tg_delta_vs_angle[v_counter].push_back(new TGraph());
            Int_t phi_counter = 0;
            for(Double_t track_phi_deg = 135.0; track_phi_deg >= 45.0; track_phi_deg -= 1.0)
            {
                Double_t track_phi   = track_phi_deg*TMath::DegToRad();
                Double_t track_slope = TMath::Tan(track_phi);

                vector< vector<Double_t> > vec_xy_pos_vD;
                vec_xy_pos_vD.resize(2); // start, stop
                vec_xy_pos_vD[0].resize(2); // x, y
                vec_xy_pos_vD[1].resize(2); // x, y

                //printf("v_drift_vD: %4.3f \n",v_drift_vD);
                Int_t flag_first_time = 0;
                for(Int_t i_time = 0; i_time < N_time; i_time++)
                {
                    Double_t time  = i_time*delta_t*1E-06; // s
                    Double_t y_pos = time*v_drift;     // original track cluster positions
                    //Double_t y_pos = time*1.56*1E06/100.0;     // original track cluster positions
                    Double_t x_pos = y_pos/track_slope;
                    //printf("phi_counter: %d, phi_counter_plot: %d, vD_counter: %d, vD_counter_plot: %d  \n",phi_counter,phi_counter_plot,vD_counter,vD_counter_plot);
                    //printf("i_time: %d, pos: {%4.3f, %4.3f} \n",i_time,x_pos,y_pos);
                    if(phi_counter == phi_counter_plot && vD_counter == vD_counter_plot && v_counter == v_counter_plot)
                    {
                        //printf("track_phi_deg: %4.3f \n",track_phi_deg);
                        tpm_track ->SetNextPoint(x_pos,y_pos);
                    }


                    if(phi_counter == phi_counter_plot && vD_counter == vD_counter_plot && v_counter == v_counter_plot)
                    {
                        tpm_track_drift[i_time]         = new TPolyMarker();
                        tpm_track_drift_Lorentz[i_time] = new TPolyMarker();
                    }
                    for(Double_t i_drift_time = 0; i_drift_time < 0.00003; i_drift_time += 0.00000003)
                    {
                        Double_t x_pos_drift = x_pos;  // drifted cluster position without Lorentz-angle effect
                        Double_t y_pos_drift = y_pos + v_drift*i_drift_time;

                        Double_t x_pos_drift_Lorentz = x_pos - TMath::Tan(alpha_L) * v_drift*i_drift_time; // drifted cluster position with Lorentz-angle effect
                        Double_t y_pos_drift_Lorentz = y_pos + v_drift*i_drift_time;

                        Double_t x_pos_Lorentz    = x_pos_drift_Lorentz; // Lorentz-angle shifted cluster position
                        Double_t y_pos_Lorentz    = y_pos_drift_Lorentz - v_drift*i_drift_time; // same as y_pos

                        Double_t x_pos_Lorentz_vD = x_pos_drift_Lorentz; // Lorentz-angle shifted cluster position with wrong drift velocity vD
                        Double_t y_pos_Lorentz_vD = y_pos_drift_Lorentz - v_drift_vD*i_drift_time;

                        //printf("v_drift: %4.3f, v_drift_vD: %4.3f \n",v_drift,v_drift_vD);

                        if(y_pos_drift > l_drift && y_pos_drift > 0.0) // at anonde wire
                        {
                            if(i_drift_time == 0) break;
                            if(phi_counter == phi_counter_plot && vD_counter == vD_counter_plot && v_counter == v_counter_plot)
                            {
                                tpm_track_Lorentz    ->SetNextPoint(x_pos_Lorentz,y_pos_Lorentz);
                                tpm_track_Lorentz_vD ->SetNextPoint(x_pos_Lorentz_vD,y_pos_Lorentz_vD);

                                //printf("i_time: %d, pos_orig: {%4.3f, %4.3f}, pos: {%4.3f, %4.3f}, pos_vD: {%4.3f, %4.3f}, i_drift_time: %4.10f, Dx: %4.4f, Dx_vD: %4.3f \n",i_time,x_pos_drift_Lorentz,y_pos_drift_Lorentz,x_pos_Lorentz,y_pos_Lorentz,x_pos_Lorentz_vD,y_pos_Lorentz_vD,i_drift_time,v_drift*i_drift_time,v_drift_vD*i_drift_time);
                            }
                            //if(y_pos_Lorentz_vD >= 0.0 && y_pos_Lorentz_vD <= l_drift && !flag_first_time)
                            if(!flag_first_time)
                            {
                                vec_xy_pos_vD[0][0] = x_pos_Lorentz_vD;
                                vec_xy_pos_vD[0][1] = y_pos_Lorentz_vD;
                                flag_first_time = 1;
                            }
                            else
                            {
                                //if(y_pos_Lorentz_vD >= 0.0 && y_pos_Lorentz_vD <= l_drift)
                                {
                                    vec_xy_pos_vD[1][0] = x_pos_Lorentz_vD;
                                    vec_xy_pos_vD[1][1] = y_pos_Lorentz_vD;
                                }
                            }
                            break;
                        }
                        if(phi_counter == phi_counter_plot && vD_counter == vD_counter_plot && v_counter == v_counter_plot)
                        {
                            tpm_track_drift[i_time]         ->SetNextPoint(x_pos_drift,y_pos_drift);
                            tpm_track_drift_Lorentz[i_time] ->SetNextPoint(x_pos_drift_Lorentz,y_pos_drift_Lorentz);
                        }
                    } // end of drift time loop
                } // end of cluster time loop


                if(phi_counter == phi_counter_plot && vD_counter == vD_counter_plot && v_counter == v_counter_plot)
                {
                    tpm_track_Lorentz_vD_used ->SetNextPoint(vec_xy_pos_vD[0][0],vec_xy_pos_vD[0][1]);
                    tpm_track_Lorentz_vD_used ->SetNextPoint(vec_xy_pos_vD[1][0],vec_xy_pos_vD[1][1]);
                }

                Double_t track_phi_vD  = TMath::ATan((vec_xy_pos_vD[1][1] - vec_xy_pos_vD[0][1])/(vec_xy_pos_vD[1][0] - vec_xy_pos_vD[0][0]));
                if(track_phi_vD < 0.0) track_phi_vD += TMath::Pi();
                //Double_t delta_phi_deg = fabs(track_phi - track_phi_vD)*TMath::RadToDeg();
                Double_t delta_phi_deg = -(track_phi - track_phi_vD)*TMath::RadToDeg();

                if(phi_counter == phi_counter_plot && vD_counter == vD_counter_plot && v_counter == v_counter_plot) printf(" ---> ");
                //printf("Lorentz angle: %4.3f degree, init angle: %4.3f, reconstructed angle: %4.3f \n",alpha_L*TMath::RadToDeg(),track_phi*TMath::RadToDeg(),track_phi_vD*TMath::RadToDeg());

                tg_delta_vs_angle[v_counter][vD_counter] ->SetPoint(phi_counter,track_phi_deg,delta_phi_deg);
                phi_counter++;
            } // end of phi loop
            vD_counter++;
        } // end of vD loop
        v_counter++;
    } // end of v loop
    //------------------------------------------------------------------------


    //------------------------------------------------------------------------
    TCanvas* can_track = new TCanvas("can_track","can_track",600,50,600,600);
    can_track->cd()->SetTicks(1,1);
    can_track->cd()->SetGrid(0,0);
    can_track->cd()->SetFillColor(10);
    can_track->cd()->SetRightMargin(0.01);
    can_track->cd()->SetLeftMargin(0.18);
    can_track->cd()->SetBottomMargin(0.12);
    can_track->cd()->SetTopMargin(0.01);
    TH1F* h_frame = can_track->cd()->DrawFrame(-0.025,-0.007,0.025,0.043,"h_frame");
    h_frame->SetStats(0);
    h_frame->SetTitle("");
    h_frame->GetXaxis()->SetTitleOffset(1.0);
    h_frame->GetYaxis()->SetTitleOffset(1.7);
    h_frame->GetXaxis()->SetLabelSize(0.05);
    h_frame->GetYaxis()->SetLabelSize(0.05);
    h_frame->GetXaxis()->SetTitleSize(0.05);
    h_frame->GetYaxis()->SetTitleSize(0.05);
    h_frame->GetXaxis()->SetNdivisions(505,'N');
    h_frame->GetYaxis()->SetNdivisions(505,'N');
    h_frame->GetXaxis()->CenterTitle();
    h_frame->GetYaxis()->CenterTitle();
    h_frame->GetXaxis()->SetTitle("x (m)");
    h_frame->GetYaxis()->SetTitle("y (m)");

    PlotLine(-0.025,0.025,0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(-0.025,0.025,0.03,0.03,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)


    tpm_track ->SetMarkerStyle(20);
    tpm_track ->SetMarkerSize(0.5);
    tpm_track ->SetMarkerColor(kRed);
    tpm_track ->Draw("same");


    tpm_track_Lorentz ->SetMarkerStyle(20);
    tpm_track_Lorentz ->SetMarkerSize(0.5);
    tpm_track_Lorentz ->SetMarkerColor(kMagenta);
    tpm_track_Lorentz ->Draw("same");

    tpm_track_Lorentz_vD ->SetMarkerStyle(20);
    tpm_track_Lorentz_vD ->SetMarkerSize(0.5);
    tpm_track_Lorentz_vD ->SetMarkerColor(kCyan);
    tpm_track_Lorentz_vD ->Draw("same");

    tpm_track_Lorentz_vD_used ->SetMarkerStyle(20);
    tpm_track_Lorentz_vD_used ->SetMarkerSize(1.0);
    tpm_track_Lorentz_vD_used ->SetMarkerColor(kOrange+1);
    tpm_track_Lorentz_vD_used ->Draw("same");

    for(Int_t i_time = 0; i_time < N_time; i_time++)
    {
        if(!(i_time % 5 == 0)) continue;
        tpm_track_drift[i_time] ->SetMarkerStyle(24);
        tpm_track_drift[i_time] ->SetMarkerSize(0.1);
        tpm_track_drift[i_time] ->SetMarkerColor(kBlue);
        tpm_track_drift[i_time] ->Draw("same");

        tpm_track_drift_Lorentz[i_time] ->SetMarkerStyle(24);
        tpm_track_drift_Lorentz[i_time] ->SetMarkerSize(0.1);
        tpm_track_drift_Lorentz[i_time] ->SetMarkerColor(kGreen);
        tpm_track_drift_Lorentz[i_time] ->Draw("same");
    }

    plotTopLegend((char*)"track",0.72,0.5,0.045,kRed,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    plotTopLegend((char*)"Lorentz angle",0.25,0.265,0.045,kMagenta,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    plotTopLegend((char*)"wrong v_{D}",0.32,0.36,0.045,kCyan,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    plotTopLegend((char*)"drift",0.44,0.6,0.045,kGreen,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    can_track ->GetCanvas()->cd(1)->Modified();
    can_track ->GetCanvas()->cd(1)->Update();
    //------------------------------------------------------------------------




    //------------------------------------------------------------------------
    Int_t arr_color_A[6] = {kBlack,kGray+2,kGray+1,kGray,kBlue-9,kBlue};
    Int_t arr_color_B[6] = {kGreen+3,kGreen-3,kRed,kRed+3,kMagenta,kMagenta-9};

    TCanvas* can_delta_vs_angle = new TCanvas("can_delta_vs_angle","can_delta_vs_angle",50,50,600,600);
    can_delta_vs_angle->cd()->SetTicks(1,1);
    can_delta_vs_angle->cd()->SetGrid(0,0);
    can_delta_vs_angle->cd()->SetFillColor(10);
    can_delta_vs_angle->cd()->SetRightMargin(0.01);
    can_delta_vs_angle->cd()->SetLeftMargin(0.18);
    can_delta_vs_angle->cd()->SetBottomMargin(0.12);
    can_delta_vs_angle->cd()->SetTopMargin(0.01);
    TH1F* h_frame_delta_vs_angle = can_delta_vs_angle->cd()->DrawFrame(70.0,-16.5,110.0,16.5,"h_frame_delta_vs_angle");
    h_frame_delta_vs_angle->SetStats(0);
    h_frame_delta_vs_angle->SetTitle("");
    h_frame_delta_vs_angle->GetXaxis()->SetTitleOffset(1.0);
    h_frame_delta_vs_angle->GetYaxis()->SetTitleOffset(1.7);
    h_frame_delta_vs_angle->GetXaxis()->SetLabelSize(0.05);
    h_frame_delta_vs_angle->GetYaxis()->SetLabelSize(0.05);
    h_frame_delta_vs_angle->GetXaxis()->SetTitleSize(0.05);
    h_frame_delta_vs_angle->GetYaxis()->SetTitleSize(0.05);
    h_frame_delta_vs_angle->GetXaxis()->SetNdivisions(505,'N');
    h_frame_delta_vs_angle->GetYaxis()->SetNdivisions(505,'N');
    h_frame_delta_vs_angle->GetXaxis()->CenterTitle();
    h_frame_delta_vs_angle->GetYaxis()->CenterTitle();
    h_frame_delta_vs_angle->GetXaxis()->SetTitle("impact angle (deg.)");
    h_frame_delta_vs_angle->GetYaxis()->SetTitle("#Delta#alpha (deg.)");

    for(Int_t i_vD = 0; i_vD < (Int_t)tg_delta_vs_angle[0].size(); i_vD++)
    //for(Int_t i_vD = 0; i_vD < (Int_t)1; i_vD++)
    {
        //tg_delta_vs_angle[0][i_vD] ->SetLineColor(arr_color_A[i_vD]);
        tg_delta_vs_angle[0][i_vD] ->SetLineColor(kRed);
        tg_delta_vs_angle[0][i_vD] ->SetLineWidth(2);
        tg_delta_vs_angle[0][i_vD] ->Draw("same");

        HistName = "";
        sprintf(NoP,"%4.2f",vec_v_drift_vals[i_vD]/10000.0);
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.22,0.15+i_vD*0.086,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        //tg_delta_vs_angle[5][i_vD] ->SetLineColor(arr_color_B[i_vD]);
        //tg_delta_vs_angle[5][i_vD] ->SetLineColor(kRed);
        //tg_delta_vs_angle[5][i_vD] ->SetLineWidth(2);
        //tg_delta_vs_angle[5][i_vD] ->Draw("same");
    }


#if 0
    TGraph* tg_Delta_vs_impact_single = calc_Delta_alpha(B_field,E_field,v_drift_use,vD_use);
    tg_Delta_vs_impact_single ->SetLineColor(kGreen+1);
    tg_Delta_vs_impact_single ->SetLineWidth(2);
    tg_Delta_vs_impact_single ->Draw("same");
#endif

    vec_tp_Delta_vs_impact[i_detector] ->SetLineColor(kBlack);
    vec_tp_Delta_vs_impact[i_detector] ->SetLineWidth(2);
    vec_tp_Delta_vs_impact[i_detector] ->SetLineStyle(1);
    vec_tp_Delta_vs_impact[i_detector] ->Draw("same hl");


    HistName = "";
    sprintf(NoP,"%4.0f",(Double_t)i_detector);
    HistName += NoP;
    HistName += ", v_{D} = ";
    sprintf(NoP,"%4.3f",v_drift_in);
    HistName += NoP;
    HistName += ", HV = ";
    sprintf(NoP,"%4.1f",HV_drift_in);
    HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.24,0.91,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    if(fCheckBox_sel[1]->GetState() == kButtonDown) // fit
    {
        HistName = "";
        sprintf(NoP,"%4.0f",(Double_t)i_detector);
        HistName += NoP;
        HistName += ", vf_{D} = ";
        sprintf(NoP,"%4.3f",v_fit*1.35);
        HistName += NoP;
        HistName += ", HVf = ";
        sprintf(NoP,"%4.1f",E_fit*l_drift);
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.24,0.85,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }


    printf("detector: %d, v_drift: %4.3f, HV_drift: %4.3f \n",i_detector,v_drift_in,HV_drift_in);

    PlotLine(70.0,110.0,0.0,0.0,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(90.0,90.0,-13.5,13.5,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

    can_delta_vs_angle ->GetCanvas()->cd(1)->Modified();
    can_delta_vs_angle ->GetCanvas()->cd(1)->Update();
    //------------------------------------------------------------------------

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t GUI_Sim_drift::Calibrate()
{
    printf("GUI_Sim_drift::Calibrate() \n");
    Pixel_t green;
    gClient->GetColorByName("green", green);
    //Button_Calibrate->ChangeBackground(green);


    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t GUI_Sim_drift::Draw3D_track()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    //Button_draw3D_track->ChangeBackground(green);



    return 1;
}
//---------------------------------------------------------------------------------


