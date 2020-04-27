#include <TFile.h>

static TGraph* tg_pad_response;
static TGraph* tg_time_response;

static Int_t N_pads = 16;
static Double_t pad_width = 0.00725; // m
static const Double_t TRD_drift_start = 0.0;
static const Double_t TRD_drift_stop  = 0.03;
static const Double_t TRD_anode_plane = 0.0335;
static const Double_t TRD_ampl_stop   = 0.037;
static const Int_t    N_clusters      = 100;
static const Double_t step_clusters   = 0.00125; // cm

vector<Double_t> recon_points;
vector<Double_t> primary_points;

TLinearFitter *lf_recon = new TLinearFitter(2);
TLinearFitter *lf_primary = new TLinearFitter(2);

static vector< vector<Double_t> > vec_pad_pos;
static vector< vector<Double_t> > vec_cluster_pos;

static TGraph* tg_rec_cluster;

static TH2D* h2D_ADC_xy;


//----------------------------------------------------------------------------------------
void Load_data()
{
    tg_pad_response  = new TGraph(58);
    tg_time_response = new TGraph();



    Double_t arr_pad_response[58][2] =
    {
        {-1.975729998 ,  0.002407932},
        {-1.865301489 ,  0.001416431},
        {-1.81645811  ,   0.000566572},
        {-1.735760353 ,  0.001416431},
        {-1.64019722  ,   0.000566572},
        {-1.565870338 ,  0.000566572},
        {-1.474554456 ,  0.000424929},
        {-1.385362198 ,  0.004390935},
        {-1.304664441 ,  0.009348442},
        {-1.22609031  ,   0.019263456},
        {-1.166628805 ,  0.033144476},
        {-1.098672799 ,  0.050991501},
        {-1.039211294 ,  0.068838527},
        {-0.990367915 ,  0.090651558},
        {-0.939400911 ,  0.112464589},
        {-0.879939406 ,  0.151133144},
        {-0.818354275 ,  0.192776204},
        {-0.765263646 ,  0.233427762},
        {-0.712173016 ,  0.28796034},
        {-0.642093385 ,  0.346458924},
        {-0.595373632 ,  0.399008499},
        {-0.535912127 ,  0.465439093},
        {-0.484945122 ,  0.521954674},
        {-0.421236367 ,  0.588385269},
        {-0.342662235 ,  0.666713881},
        {-0.287447981 ,  0.717280453},
        {-0.247099102 ,  0.741076487},
        {-0.185513972 ,  0.780736544},
        {-0.138794218 ,  0.800566572},
        {-0.100568965 ,  0.815439093},
        {-0.051725586 ,  0.822379603},
        {0.001365044  ,   0.828328612},
        {0.054455673  ,   0.823371105},
        {0.113917178  ,   0.811473088},
        {0.175502308  ,   0.787677054},
        {0.215851187  ,   0.769830028},
        {0.256200065  ,   0.735127479},
        {0.313537945  ,   0.692492918},
        {0.37299945   ,   0.63796034},
        {0.411224703  ,   0.604249292},
        {0.462191708  ,   0.544759207},
        {0.528024088  ,   0.47733711},
        {0.595980094  ,   0.4},
        {0.657565224  ,   0.332577904},
        {0.717026729  ,   0.276062323},
        {0.78710636   ,   0.214589235},
        {0.889040369  ,   0.143201133},
        {0.954872749  ,   0.104532578},
        {1.022828755  ,   0.073796034},
        {1.109897387  ,   0.044050992},
        {1.20546052   ,   0.023229462},
        {1.322259905  ,   0.005382436},
        {1.419946663  ,   0.001416431},
        {1.490026294  ,   0.000566572},
        {1.596207553  ,   0.000424929},
        {1.706636063  ,   0.002407932},
        {1.836177199  ,   0.000424929},
        {1.959347459  ,   0.000424929}
    };


    Double_t arr_time_response[73][2] =
    {
        {0.009707273,	0.288065844},
        {0.050055426 ,   0.576131687},
        {0.097838453 ,   2.77E-12   },
        {0.146681107 ,   0.288065844},
        {0.167916632 ,   0.576131687},
        {0.177466392 ,   3.16872428 },
        {0.182765259 ,   7.201646091},
        {0.18805466  ,   14.97942387},
        {0.192277877 ,   24.48559671},
        {0.197558538 ,   35.72016461},
        {0.202832645 ,   49.5473251 },
        {0.205978756 ,   65.10288066},
        {0.210186681 ,   80.65843621},
        {0.215459331 ,   95.0617284 },
        {0.220718872 ,   114.6502058},
        {0.223857701 ,   133.0864198},
        {0.228061984 ,   150.0823045},
        {0.233325895 ,   167.9423868},
        {0.23754183  ,   180.3292181},
        {0.241769417 ,   188.1069959},
        {0.244929366 ,   198.1893004},
        {0.249143116 ,   211.4403292},
        {0.256545218 ,   223.5390947},
        {0.260761153 ,   235.9259259},
        {0.264988012 ,   243.9917695},
        {0.268155244 ,   251.1934156},
        {0.274512282 ,   256.6666667},
        {0.27768971  ,   259.8353909},
        {0.276622799 ,   261.8518519},
        {0.28086568  ,   263.5802469},
        {0.286169646 ,   265.5967078},
        {0.290416168 ,   265.8847737},
        {0.294665603 ,   265.0205761},
        {0.298916495 ,   263.5802469},
        {0.304229928 ,   261.8518519},
        {0.306357195 ,   260.4115226},
        {0.313798623 ,   256.9547325},
        {0.315931716 ,   253.2098765},
        {0.318065537 ,   249.1769547},
        {0.323384796 ,   245.1440329},
        {0.329770237 ,   239.382716 },
        {0.335102605 ,   230.1646091},
        {0.341500427 ,   219.5061728},
        {0.354284418 ,   202.7983539},
        {0.366002227 ,   187.81893  },
        {0.376664778 ,   170.2469136},
        {0.389443672 ,   155.5555556},
        {0.402228392 ,   138.5596708},
        {0.416078565 ,   120.1234568},
        {0.440554877 ,   98.51851852},
        {0.465020264 ,   81.2345679 },
        {0.485218009 ,   72.01646091},
        {0.511791728 ,   60.781893  },
        {0.540484703 ,   51.27572016},
        {0.581916517 ,   42.9218107 },
        {0.614849459 ,   36.2962963 },
        {0.655204895 ,   33.7037037 },
        {0.692377078 ,   30.24691358},
        {0.743354281 ,   26.21399177},
        {0.803885614 ,   23.04526749},
        {0.873973989 ,   19.58847737},
        {0.946183805 ,   16.99588477},
        {1.033258998 ,   14.40329218},
        {1.122457088 ,   12.09876543},
        {1.205281389 ,   10.94650206},
        {1.305095421 ,   9.506172839},
        {1.436762377 ,   8.641975309},
        {1.552503598 ,   7.201646091},
        {1.659749591 ,   6.049382716},
        {1.758498169 ,   6.049382716},
        {1.833888324 ,   5.473251029},
        {1.916711169 ,   4.897119342},
        {1.997408203 ,   5.185185185}
    };



    for(Int_t i_point = 0; i_point < 58; i_point++)
    {
        tg_pad_response ->SetPoint(i_point,arr_pad_response[i_point][0],arr_pad_response[i_point][1]);
    }

    for(Int_t i_point = 0; i_point < 73; i_point++)
    {
        tg_time_response ->SetPoint(i_point,arr_time_response[i_point][0],arr_time_response[i_point][1]);
    }


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



//----------------------------------------------------------------------------------------
// Calculates integral of graph between x_start and x_stop
Double_t calc_integral(TGraph *g, Double_t x_start, Double_t x_stop)
{
     Double_t area = 0.0;
     Double_t x_point1 = 0.0;
     Double_t x_point2 = 0.0;    
     Double_t y_point1 = 0.0;
     Double_t y_point2 = 0.0;
     Int_t n_points = 0;

     for (Int_t i_point = 0; i_point < g ->GetN() -1; i_point++)
     {
        g ->GetPoint(i_point,x_point1,y_point1);

        if (x_point1 < x_start) continue;
        
        g ->GetPoint(i_point+1,x_point2,y_point2);

        if (x_point2 > x_stop) break;

        area += (x_point2-x_point1)*((y_point1+y_point2)/2);
        n_points++;
        
        //printf("\n");
        //printf("x_start: %4.3f, x_stop: %4.3f \n",x_start,x_stop);
        //printf("i_point: %d, x_point1: %4.3f, x_point2: %4.3f, y_point1: %4.3f, y_point2: %4.3f, area: %4.3f \n",i_point,x_point1,x_point2,y_point1,y_point2,area);
        //printf("\n");

     }

     if (n_points == 0) 
     {
        area = g ->Eval(x_point1);
     }
     
     return area;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TRD_detector_2D()
{
    tg_rec_cluster = new TGraph();

    vector<Double_t> vec_start_stop;
    vec_start_stop.resize(2);

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

    Double_t pad_pos = -x_range;
    Int_t i_pad = 0;
    while(pad_pos < x_range)
    {
        TBox* box_pad = new TBox(pad_pos+0.0003,TRD_ampl_stop,pad_pos+pad_width-0.0003,TRD_ampl_stop+0.001);
        box_pad->SetFillColor(kBlue);
        box_pad->SetFillStyle(3001);
        box_pad->Draw("same");

        vec_start_stop[0] = pad_pos;
        vec_start_stop[1] = pad_pos + pad_width;
        vec_pad_pos.push_back(vec_start_stop);

        pad_pos += pad_width;

        printf("i_pad: %d, pos: {%4.5f, %4.5f} \n",i_pad,vec_start_stop[0],vec_start_stop[1]);
        i_pad++;
    }

    h2D_ADC_xy = new TH2D("h2D_ADC_xy","h2D_ADC_xy",(Int_t)vec_pad_pos.size(),vec_pad_pos[0][0],vec_pad_pos[(Int_t)vec_pad_pos.size()-1][1],24,0.0,0.03);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector2* Create_TRD_track(Double_t impact_angle, Double_t Lorentz_angle, Double_t drift_vel_ratio, Double_t Lorentz_angle_pre_corr)
{
    // Lorentz_angle_pre_corr: The Lorentz angle pre correction which was applied for the data used for calibration
    lf_primary->SetFormula("x ++ 1");

    vec_cluster_pos.clear();
    vector<Double_t> vec_pos;
    vec_pos.resize(2);

    Double_t x_dir = TMath::Cos(impact_angle);
    Double_t y_dir = TMath::Sin(impact_angle);
    Double_t slope = 1.0;
    if(x_dir != 0.0) slope = y_dir/x_dir;
    x_dir /= y_dir;
    y_dir /= y_dir;

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
        Double_t y_val = 0.5*step_clusters + TRD_drift_start + y_dir*step_clusters*i_cluster;

        //Lorentz angle correction

        Double_t shift = TMath::Tan(Lorentz_angle)*(0.03 - y_val);
        x_val += shift;

        lf_primary->AddPoint(&x_val,y_val);

        if(y_val > TRD_ampl_stop) break;
        TPM_trd_track_clusters ->SetNextPoint(x_val,y_val);

        //printf("i_cluster: %d, step_clusters: %4.3f, y_dir*step_clusters*i_cluster: %4.3f, x_val: %4.3f, y_val: %4.3f \n \n",i_cluster,step_clusters,y_dir*step_clusters*i_cluster,x_val,y_val);
        vec_pos[0] = x_val;
        vec_pos[1] = y_val;

        

        vec_cluster_pos.push_back(vec_pos);
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
    Double_t y_Lorentz_drift_hit_pre_corr = TRD_anode_plane - TRD_anode_plane*drift_vel_ratio;
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
Int_t Find_pad_index(Double_t x_pos, Double_t &rel_pos_center)
{
    for(Int_t i_pad = 0; i_pad < (Int_t)vec_pad_pos.size(); i_pad++)
    {
        Double_t pad_start = vec_pad_pos[i_pad][0];
        Double_t pad_stop  = vec_pad_pos[i_pad][1];

        if(x_pos > pad_start && x_pos <= pad_stop)
        {
            Double_t center_pos = (pad_start + pad_stop)/2.0;
            rel_pos_center = (x_pos - center_pos)/pad_width;
            //printf("x_pos: %4.5f, center_pos: %4.5f, rel_pos_center: %4.5f \n",x_pos,center_pos,rel_pos_center);
            return i_pad;
        }
    }

    return -1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Make_ion_tail_convolution()
{
    printf("Make_ion_tail_convolution \n");

    h2D_ADC_xy ->Reset();
    Double_t prf = 0.0;
    Double_t trf = 0.0;

    for(Int_t i_cluster = 0; i_cluster < (Int_t)vec_cluster_pos.size(); i_cluster++)
    {
        Double_t x_pos = vec_cluster_pos[i_cluster][0];
        Double_t y_pos = vec_cluster_pos[i_cluster][1];

        if(y_pos > TRD_drift_stop) break;

        Double_t rel_pos_center = 0.0;
        Double_t rel_pos_x_start = 0.0;
        Double_t rel_pos_x_stop = 0.0;

        Int_t i_pad = Find_pad_index(x_pos,rel_pos_center);
        // printf(" \n");
        // printf("i_cluster: %d, pos: {%4.5f, %4.5f}, i_pad: %d, rel_pos_center: %4.5f \n",i_cluster,x_pos,y_pos,i_pad,rel_pos_center);

        for(Int_t i_ion_tail = i_cluster; i_ion_tail >= 0.0; i_ion_tail--)
        {
            Double_t i_time       = i_ion_tail*0.1; // in mus
            Double_t y_pos_charge = y_pos - i_time*1.25*0.01; // 1.56 cm/mus, one time bin = 0.1 mus


            for(Int_t i_charge_share = -2; i_charge_share <= 2; i_charge_share++) // two pads on the left and right for charge sharing
            {
                Double_t x_pos_charge = x_pos + i_charge_share*pad_width;
                Double_t rel_pos_center_offset = (Double_t)i_charge_share; // in pad units
                Double_t rel_pos_center_pad_response = rel_pos_center;

                //printf("x_pos: %4.3f, x_pos_charge: %4.3f, rel_pos_center_offset: %4.3f, rel_pos_center_pad_response: %4.3f \n",x_pos,x_pos_charge,rel_pos_center_offset,rel_pos_center_pad_response);
                //printf("vec_pad_pos[Find_pad_index(x_pos_charge,rel_pos_x_start)][0]: %4.3f, vec_pad_pos[Find_pad_index(x_pos_charge,rel_pos_x_start)][1]: %4.3f \n",vec_pad_pos[Find_pad_index(x_pos_charge,rel_pos_x_start)][0],vec_pad_pos[Find_pad_index(x_pos_charge,rel_pos_x_start)][1]);
                trf = calc_integral(tg_time_response,0.3 + i_time - 0.05, 0.3 + i_time + 0.05);
                //prf = calc_integral(tg_pad_response, floor(rel_pos_center_pad_response),ceil(rel_pos_center_pad_response));

                Double_t low_int_prf = (-0.5 + rel_pos_center_offset) - rel_pos_center_pad_response;
                Double_t up_int_prf  = (0.5  + rel_pos_center_offset) - rel_pos_center_pad_response;
                prf = calc_integral(tg_pad_response,low_int_prf,up_int_prf);

                //printf("tg_pad_response->Eval(rel_pos_center_pad_response): %4.3f \n",tg_pad_response->Eval(rel_pos_center_pad_response));
                //Double_t ADC          = tg_time_response ->Eval(0.3 + i_time)*tg_pad_response->Eval(rel_pos_center_pad_response);
                //Double_t ADC          = tg_time_response ->Eval(0.3 + i_time)*prf;

                Double_t ADC = trf*prf;

                //compensate for aplification region effect
                //I think this is correct?
                if (i_cluster >= 21)
                {
                    ADC *= 2;
                }
                
                // printf("   i_charge_share: %d, i_ion_tail: %d, i_time: %4.3f, y_pos_charge: %4.5f, ADC: %4.3f \n",i_charge_share,i_ion_tail,i_time,y_pos_charge,ADC);

                h2D_ADC_xy ->Fill(x_pos_charge,y_pos_charge,ADC);

                // printf("   - i_charge_share: %d, x_pos_charge: %4.3f, rel_pos_center: %4.3f, rel_pos_center_pad_response: %4.3f, low_int_prf: %4.3f, up_int_prf: %4.3f, prf: %4.3f, trf: %4.3f, ADC: %4.3f \n",i_charge_share,x_pos_charge,rel_pos_center,rel_pos_center_pad_response,low_int_prf,up_int_prf,prf,trf,ADC);
            }
        }
    }

    TCanvas* can_ADC_xy = new TCanvas("can_ADC_xy","can_ADC_xy",600,50,800,600);
    can_ADC_xy->cd()->SetTicks(1,1);
    can_ADC_xy->cd()->SetGrid(0,0);
    can_ADC_xy->cd()->SetFillColor(10);
    can_ADC_xy->cd()->SetRightMargin(0.15);
    can_ADC_xy->cd()->SetLeftMargin(0.18);
    can_ADC_xy->cd()->SetBottomMargin(0.12);
    can_ADC_xy->cd()->SetTopMargin(0.01);

    h2D_ADC_xy->SetStats(0);
    h2D_ADC_xy->SetTitle("");
    h2D_ADC_xy->GetXaxis()->SetTitleOffset(1.0);
    h2D_ADC_xy->GetYaxis()->SetTitleOffset(1.7);
    h2D_ADC_xy->GetXaxis()->SetLabelSize(0.05);
    h2D_ADC_xy->GetYaxis()->SetLabelSize(0.05);
    h2D_ADC_xy->GetXaxis()->SetTitleSize(0.05);
    h2D_ADC_xy->GetYaxis()->SetTitleSize(0.05);
    h2D_ADC_xy->GetXaxis()->SetNdivisions(505,'N');
    h2D_ADC_xy->GetYaxis()->SetNdivisions(505,'N');
    h2D_ADC_xy->GetXaxis()->CenterTitle();
    h2D_ADC_xy->GetYaxis()->CenterTitle();
    h2D_ADC_xy->GetXaxis()->SetTitle("x (m)");
    h2D_ADC_xy->GetYaxis()->SetTitle("y (m)");
    h2D_ADC_xy->DrawCopy("colz");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Reconstruct_tracklet(Double_t Lorentz_angle)
{
    printf("Reconstruct_tracklet \n");

    lf_recon->SetFormula("x ++ 1");

    tg_rec_cluster ->Clear();

    for(Int_t i_biny = 1; i_biny <= h2D_ADC_xy ->GetNbinsY(); i_biny++)
    {
        Double_t y_pos = h2D_ADC_xy ->GetYaxis()->GetBinCenter(i_biny);
        Double_t x_pos_weighted = 0.0;
        Double_t sum_ADC        = 0.0;
        for(Int_t i_binx = 1; i_binx <= h2D_ADC_xy ->GetNbinsX(); i_binx++)
        {
            Double_t x_pos = h2D_ADC_xy ->GetXaxis()->GetBinCenter(i_binx);
            Double_t ADC   = h2D_ADC_xy ->GetBinContent(i_binx,i_biny);

            x_pos_weighted += x_pos*ADC;
            sum_ADC        += ADC;
        }

        if(sum_ADC > 0.0)
        {
            x_pos_weighted /= sum_ADC;
        }
        else x_pos_weighted = -999.0;

        tg_rec_cluster ->SetPoint(i_biny - 1,x_pos_weighted,y_pos);
        lf_recon->AddPoint(&x_pos_weighted, y_pos);
    }

    tg_rec_cluster ->SetMarkerSize(0.8);
    tg_rec_cluster ->SetMarkerStyle(20);
    tg_rec_cluster ->SetMarkerColor(kGreen+2);
    tg_rec_cluster ->Draw("same p");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void straight_line_fits()
{
    lf_recon->Eval();
    lf_primary->Eval();

    Double_t parFit_recon[2];
    for(int i = 0; i < 2; ++i)
    {
        parFit_recon[i] = lf_recon->GetParameter(i);
        cout <<  "recon line: " << parFit_recon[i] << endl;
    }

    Double_t parFit_primary[2];
    for(int i = 0; i < 2; ++i)
    {
        parFit_primary[i] = lf_primary->GetParameter(i);
        cout << "primary line: " << parFit_primary[i] << endl;
    }

    // range [0, 0.01]
    for(int i = 0; i < 100; ++i)
    {
        Double_t x = 0.01/100 * i;
        Double_t y_recon = parFit_recon[0]*x + parFit_recon[1];
        Double_t y_primary = parFit_primary[0]*x + parFit_primary[1];
        // cout << "x val: " << x << endl;
        // cout << "y vals: " << y_recon << "    " << y_primary << endl;

        cout << "ratio at x = " << x << ": " << y_recon/y_primary << endl;
    } 
}
//----------------------------------------------------------------------------------------
