

static const Double_t TRD_drift_start = 0.0;
static const Double_t TRD_drift_stop  = 0.03;
static const Double_t TRD_anode_plane = 0.0335;
static const Double_t TRD_ampl_stop   = 0.037;
static const Int_t    N_clusters      = 100;
static const Double_t step_clusters   = 0.001; // cm

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



//----------------------------------------------------------------------------------------
void Draw_TRD_detector_2D()
{
    TCanvas* can_track = new TCanvas("can_track","can_track",600,50,800,600);
    can_track->cd()->SetTicks(1,1);
    can_track->cd()->SetGrid(0,0);
    can_track->cd()->SetFillColor(10);
    can_track->cd()->SetRightMargin(0.01);
    can_track->cd()->SetLeftMargin(0.18);
    can_track->cd()->SetBottomMargin(0.12);
    can_track->cd()->SetTopMargin(0.01);

    Double_t x_range = 0.035;

    TH1F* h_frame = can_track->cd()->DrawFrame(-x_range,-0.007,x_range,0.043,"h_frame");
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

    PlotLine(-x_range,x_range,TRD_drift_start,TRD_drift_start,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(-x_range,x_range,TRD_drift_stop,TRD_drift_stop,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(-x_range,x_range,TRD_ampl_stop,TRD_ampl_stop,kBlack,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    PlotLine(-x_range,x_range,TRD_anode_plane,TRD_anode_plane,kBlue+1,3,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    TBox* box_amp = new TBox(-x_range,TRD_drift_stop,x_range,TRD_ampl_stop);
    box_amp->SetFillColor(kRed-8);
    box_amp->SetFillStyle(3001);
    box_amp->Draw("same");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector2* Create_TRD_track(Double_t impact_angle, Double_t Lorentz_angle, Double_t drift_vel_ratio, Double_t Lorentz_angle_pre_corr)
{
    // Lorentz_angle_pre_corr: The Lorentz angle pre correction which was applied for the data used for calibration

    Double_t x_dir = TMath::Cos(impact_angle);
    Double_t y_dir = TMath::Sin(impact_angle);
    Double_t slope = 1.0;
    if(x_dir != 0.0) slope = y_dir/x_dir;

    Double_t Lorentz_tan   = TMath::Tan(Lorentz_angle);
    Double_t Lorentz_slope = 1.0;
    if(Lorentz_tan != 0.0) Lorentz_slope = 1.0/Lorentz_tan;

    TPolyMarker* TPM_trd_track_clusters      = new TPolyMarker();
    TPolyMarker* TPM_trd_anode_plane         = new TPolyMarker();
    TPolyMarker* TPM_trd_Lorentz_anode_plane = new TPolyMarker();
    TPolyMarker* TPM_trd_Lorentz_drift       = new TPolyMarker();
    TPolyMarker* TPM_trd_Lorentz_drift_pre_corr       = new TPolyMarker();

    for(Int_t i_cluster = 0; i_cluster < N_clusters; i_cluster++)
    {
        Double_t x_val = x_dir*step_clusters*i_cluster;
        Double_t y_val = TRD_drift_start + y_dir*step_clusters*i_cluster;

        if(y_val > TRD_ampl_stop) break;
        TPM_trd_track_clusters ->SetNextPoint(x_val,y_val);
    }

    TPM_trd_track_clusters ->SetMarkerSize(0.5);
    TPM_trd_track_clusters ->SetMarkerColor(kRed);
    TPM_trd_track_clusters ->SetMarkerStyle(20);
    TPM_trd_track_clusters ->Draw("");

    Double_t x_anode_hit = TRD_anode_plane/slope;
    Double_t y_anode_hit = TRD_anode_plane;
    TPM_trd_anode_plane ->SetNextPoint(x_anode_hit,y_anode_hit);
    TPM_trd_anode_plane ->SetMarkerSize(2.0);
    TPM_trd_anode_plane ->SetMarkerColor(kGreen+1);
    TPM_trd_anode_plane ->SetMarkerStyle(29);
    TPM_trd_anode_plane ->Draw("");

    // Physical hit of first cluster at anode wire after Lorentz shift
    Double_t x_Lorentz_anode_hit = TRD_anode_plane/Lorentz_slope;
    Double_t y_Lorentz_anode_hit = TRD_anode_plane;
    TPM_trd_Lorentz_anode_plane ->SetNextPoint(x_Lorentz_anode_hit,y_Lorentz_anode_hit);
    TPM_trd_Lorentz_anode_plane ->SetMarkerSize(2.0);
    TPM_trd_Lorentz_anode_plane ->SetMarkerColor(kBlack);
    TPM_trd_Lorentz_anode_plane ->SetMarkerStyle(29);
    TPM_trd_Lorentz_anode_plane ->Draw("");

    // Reconstructed hit of first cluster at chamber entrance before pre Lorentz angle correction applied
    Double_t x_Lorentz_drift_hit = x_Lorentz_anode_hit;
    Double_t y_Lorentz_drift_hit = TRD_anode_plane - TRD_anode_plane*drift_vel_ratio;
    TPM_trd_Lorentz_drift ->SetNextPoint(x_Lorentz_drift_hit,y_Lorentz_drift_hit);
    TPM_trd_Lorentz_drift ->SetMarkerSize(2.0);
    TPM_trd_Lorentz_drift ->SetMarkerColor(kCyan+1);
    TPM_trd_Lorentz_drift ->SetMarkerStyle(29);
    TPM_trd_Lorentz_drift ->Draw("");

    // Reconstructed hit of first cluster at chamber entrance after pre Lorentz angle correction
    Double_t x_Lorentz_drift_hit_pre_corr = x_Lorentz_anode_hit - (TRD_anode_plane - y_Lorentz_drift_hit)*TMath::Tan(Lorentz_angle_pre_corr);
    Double_t y_Lorentz_drift_hit_pre_corr = y_Lorentz_drift_hit;
    TPM_trd_Lorentz_drift_pre_corr ->SetNextPoint(x_Lorentz_drift_hit_pre_corr,y_Lorentz_drift_hit_pre_corr);
    TPM_trd_Lorentz_drift_pre_corr ->SetMarkerSize(2.0);
    TPM_trd_Lorentz_drift_pre_corr ->SetMarkerColor(kOrange+1);
    TPM_trd_Lorentz_drift_pre_corr ->SetMarkerStyle(29);
    TPM_trd_Lorentz_drift_pre_corr ->Draw("");

    PlotLine(x_Lorentz_drift_hit_pre_corr,x_anode_hit,y_Lorentz_drift_hit_pre_corr,y_anode_hit,kGreen+2,2,2); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

    TVector2* TV2_trd_track = new TVector2(x_dir,y_dir);
    return TV2_trd_track;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TGraph* Create_TRD_Delta_alpha_vs_impact_angle(Double_t Lorentz_angle, Double_t drift_vel_ratio, Double_t Lorentz_angle_pre_corr)
{
    TGraph* TG_Delta_alpha_vs_impact_angle = new TGraph();

    Int_t i_point = 0;
    for(Double_t impact_angle = 65.0*TMath::DegToRad(); impact_angle < 115.0*TMath::DegToRad(); impact_angle += 1.0*TMath::DegToRad())
    {

        // Direction vector of incoming track
        Double_t x_dir = TMath::Cos(impact_angle);
        Double_t y_dir = TMath::Sin(impact_angle);

        // Slope of incoming track
        Double_t slope = 10000000.0;
        if(x_dir != 0.0) slope = y_dir/x_dir;

        // Slope of Lorentz angle
        Double_t Lorentz_tan   = TMath::Tan(Lorentz_angle);
        Double_t Lorentz_slope = 10000000.0;
        if(Lorentz_tan != 0.0) Lorentz_slope = 1.0/Lorentz_tan;

        // Hit point of incoming track with anode plane
        Double_t x_anode_hit = TRD_anode_plane/slope;
        Double_t y_anode_hit = TRD_anode_plane;

        // Hit point at anode plane of Lorentz angle shifted cluster from the entrance
        Double_t x_Lorentz_anode_hit = TRD_anode_plane/Lorentz_slope;
        Double_t y_Lorentz_anode_hit = TRD_anode_plane;

        // Cluster location within drift cell of cluster from entrance after drift velocity ratio is applied
        Double_t x_Lorentz_drift_hit = x_Lorentz_anode_hit;
        Double_t y_Lorentz_drift_hit = TRD_anode_plane - TRD_anode_plane*drift_vel_ratio;

        // Reconstructed hit of first cluster at chamber entrance after pre Lorentz angle correction
        Double_t x_Lorentz_drift_hit_pre_corr = x_Lorentz_anode_hit - y_Lorentz_drift_hit*TMath::Tan(Lorentz_angle_pre_corr);
        Double_t y_Lorentz_drift_hit_pre_corr = TRD_anode_plane - TRD_anode_plane*drift_vel_ratio;

        Double_t impact_angle_track = TMath::ATan2(y_anode_hit,x_anode_hit);

        Double_t Delta_x_Lorentz_drift_hit = x_anode_hit - x_Lorentz_drift_hit_pre_corr;
        Double_t Delta_y_Lorentz_drift_hit = y_anode_hit - y_Lorentz_drift_hit_pre_corr;
        Double_t impact_angle_rec   = TMath::ATan2(Delta_y_Lorentz_drift_hit,Delta_x_Lorentz_drift_hit);

        Double_t Delta_angle = -(impact_angle_track - impact_angle_rec);
        TG_Delta_alpha_vs_impact_angle ->SetPoint(i_point,impact_angle*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
        i_point++;
        //printf("impact_angle: %4.3f, impact_angle_track: %4.3f, impact_angle_rec: %4.3f, Delta: {%4.3f, %4.3f}, x_anode_hit: %4.3f, x_Lorentz_drift_hit: %4.3f \n",impact_angle*TMath::RadToDeg(),impact_angle_track*TMath::RadToDeg(),impact_angle_rec*TMath::RadToDeg(),Delta_x_Lorentz_drift_hit,Delta_y_Lorentz_drift_hit,x_anode_hit,x_Lorentz_drift_hit);
    }

    return TG_Delta_alpha_vs_impact_angle;
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
void TRD_drift_sim()
{
    printf("TRD_drift_sim started \n");

    Double_t impact_angle             = 75.0;
    Double_t Lorentz_angle            = -10;
    Double_t Drift_vel_ratio          = 1.1; // 0.8
    Double_t Lorentz_angle_pre_corr   = -16; // -7.5

    Draw_TRD_detector_2D();
    Create_TRD_track(impact_angle*TMath::DegToRad(),Lorentz_angle*TMath::DegToRad(),Drift_vel_ratio,Lorentz_angle_pre_corr*TMath::DegToRad());

    TCanvas* can_Delta_alpha_vs_impact = new TCanvas("can_Delta_alpha_vs_impact","can_Delta_alpha_vs_impact",10,10,600,600);
    can_Delta_alpha_vs_impact ->cd();
    TGraph* TG_Delta_alpha_vs_impact = Create_TRD_Delta_alpha_vs_impact_angle(Lorentz_angle*TMath::DegToRad(),Drift_vel_ratio,Lorentz_angle_pre_corr*TMath::DegToRad());
    TG_Delta_alpha_vs_impact ->Draw("Al");

}
//----------------------------------------------------------------------------------------



