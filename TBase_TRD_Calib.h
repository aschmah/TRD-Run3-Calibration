#ifndef __TBASE_TRD_CALIB_H__
#define __TBASE_TRD_CALIB_H__

using namespace std;
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TChain.h"

#include "TObject.h"

#include "Ali_AS_Event.h"
#include "Ali_AS_EventLinkDef.h"
#include "Ana_Digits_functions.h"

ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_Event)



//----------------------------------------------------------------------------------------
Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
									 TVector3 &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.Mag()>0))
    {
      return -1000000.;
    }
  
  TVector3 diff = base-point;

  TVector3 cross = dir.Cross(diff);
  
  return cross.Mag()/dir.Mag();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class TBase_TRD_Calib
{
private:
    Ali_AS_Event*     AS_Event;
    Ali_AS_Track*     AS_Track;
    Ali_AS_Tracklet*  AS_Tracklet;
    Ali_AS_TRD_digit* AS_Digit;


    TFile* outputfile;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;

    TChain* input_SE;
    TString JPsi_TREE   = "Tree_AS_Event";
    TString JPsi_BRANCH = "Tree_AS_Event_branch";
    Long64_t file_entries_total;

    Float_t digit_pos[3];

    TPolyMarker3D* TPM3D_single;
    vector<TPolyMarker3D*> vec_TPM3D_digits;
    vector<TPolyLine3D*> vec_TPL3D_tracklets;
    vector<TPolyLine3D*> vec_TPL3D_tracklets_match;
    vector<Float_t> vec_digit_single_info;
    vector< vector<Float_t> > vec_digit_info;
    vector< vector< vector<Float_t> > > vec_digit_track_info;
    vector<Float_t> vec_track_single_info;
    vector< vector<Float_t> > vec_track_info;

    AliHelix aliHelix;
    TPolyLine3D* TPL3D_helix;
    TPolyLine3D* fit_line;
    TPolyLine3D* vec_tracklets_line;

    TString HistName;
    char NoP[50];

    // TVector3 vec_datapoints;
    vector <TCanvas*> ADC_vs_time;

    vector<TPolyLine3D*> vec_TPL3D_helix_neighbor;
    vector< vector<TPolyLine3D*> > vec_TPL3D_TRD_center;
    vector< vector<TVector3> >     vec_TV3_TRD_center; // 540 chambers, 3 axes

    Int_t arr_layer_detector[6];

    TGLViewer *TGL_viewer;

    TH1D* h_delta_angle_perp_impact;
    TH1D* h_detector_hit;


    // TRD offline Tracklets
    vector< vector<TVector3> > vec_TV3_Tracklet_pos;
    vector< vector<TVector3> > vec_TV3_Tracklet_dir;


    // TRD 3D graphics
    AliTRDgeometry* fGeo;
    TGeoManager     *geom;
    TGeoMaterial    *vacuum;
    TGeoMedium      *Air;
    TGeoVolume      *top;
    TGeoMaterial *Fe;
    TGeoMaterial *M_outer_tube;
    TGeoMaterial *Material_TRD_box;
    TGeoMedium *Iron;
    TGeoMedium *Me_outer_tube;
    TGeoMedium *Medium_TRD_box;

    TGeoVolume *inner_field_tube;
    TGeoVolume *outer_field_tube;


    TGeoCombiTrans* combitrans[540];
    TGeoVolume *TRD_boxes[540];

    TPolyLine3D* z_BeamLine;
    TPolyLine3D* x_BeamLine;
    TPolyLine3D* y_BeamLine;

    Int_t Not_installed_TRD_detectors[19] = {402,403,404,405,406,407,432,433,434,435,436,437,462,463,464,465,466,467,538};
    Int_t Defect_TRD_detectors[84]        = {2,5,17,24,26,30,31,32,36,40,41,43,49,50,59,62,64,78,88,92,93,107,110,111,113,116,119,131,161,
    165,182,184,188,190,191,215,219,221,223,226,227,233,236,239,241,249,255,265,277,287,302,308,310,311,318,319,320,326,328,335,348,354,368,377,380,
    386,389,452,455,456,470,474,476,483,484,485,490,491,493,494,500,502,504,506};

    Double_t max_dca_z_to_track = 8.0; // in cm
    Double_t max_dca_r_to_track = 1.0; // in cm
    vector<Int_t> vec_merge_time_bins;

    vector< vector<TVector3> > vec_TV3_digit_pos_cluster;    // layer, merged time bin
    vector< TVector3> vec_TV3_digit_pos_cluster_t0; // layer, x, y, z
    vector< vector<TH1F*> > th1f_ADC_vs_time;
    Int_t color_layer[7] = {kOrange+2,kGreen,kBlue,kMagenta,kCyan,kPink+5,kYellow};
    Int_t line_width_layer[7] = {3,3,3,3,3,3,10};
    vector<Int_t> vec_layer_in_fit;
    vector< vector< vector<Double_t> > > vec_tracklet_fit_points;
    vector<TProfile*> vec_tp_Delta_vs_impact;



public:
    TBase_TRD_Calib();
    ~TBase_TRD_Calib();
    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t event);
    Long64_t get_N_Events() {return N_Events;}
    Long64_t get_N_Tracks() {return N_Tracks;}
    Long64_t get_N_Digits() {return N_Digits;}
    vector< vector<Float_t> > get_track_info() {return vec_track_info;}
    vector< vector< vector<Float_t> > > get_digit_track_info() {return vec_digit_track_info;}
    vector<TPolyMarker3D*> get_PM3D_digits() {return vec_TPM3D_digits;}
    vector<TPolyLine3D*> get_PL3D_tracklets() {return vec_TPL3D_tracklets;}
    vector<TPolyLine3D*> get_PL3D_tracklets_match() {return vec_TPL3D_tracklets_match;}
    void Draw_track(Int_t i_track);
    void Draw_line(Int_t i_track);
    void Draw_tracklets_line(Int_t i_track);
    void Draw_neighbor_tracks(Int_t i_track);
    TGLViewer* Draw_TRD();
    void set_dca_to_track(Double_t dca_r, Double_t dca_z) {max_dca_r_to_track = dca_r; max_dca_z_to_track = dca_z;}
    void set_merged_time_bins(vector<Int_t> vec_merge_time_bins_in) {vec_merge_time_bins = vec_merge_time_bins_in;}
    TPolyLine3D* get_helix_polyline(Int_t i_track);
    TPolyLine3D* get_straight_line_fit(Int_t i_track);
    void get_tracklets_fit(Int_t i_track);
    vector< vector<TVector3> >  make_clusters(Int_t i_track);
    vector< vector< vector<Double_t> > > get_tracklet_fit_points() {return vec_tracklet_fit_points;}
    void make_plots_ADC(Int_t i_track);
    void Calibrate();
    void Track_Tracklets();

    ClassDef(TBase_TRD_Calib, 1)
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TBase_TRD_Calib::TBase_TRD_Calib()
{
    outputfile = new TFile("./TRD_Calib.root","RECREATE");

    vec_TV3_Tracklet_pos.resize(540);
    vec_TV3_Tracklet_dir.resize(540);

    // Standard time bins
    vec_merge_time_bins.resize(24+1);
    for(Int_t i_time = 0; i_time < (24+1); i_time++)
    {
        vec_merge_time_bins[i_time] = i_time;
    }

    vec_TPL3D_TRD_center.resize(540);
    vec_TV3_TRD_center.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TPL3D_TRD_center[i_det].resize(3);
        vec_TV3_TRD_center[i_det].resize(3);
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_TPL3D_TRD_center[i_det][i_xyz] = new TPolyLine3D();
        }
    }

    h_delta_angle_perp_impact = new TH1D("h_delta_angle_perp_impact","h_delta_angle_perp_impact",240,-30,30);
    h_detector_hit            = new TH1D("h_detector_hit","h_detector_hit",540,0,540);

    vec_layer_in_fit.resize(7);
    vec_tracklets_line = new TPolyLine3D();


    vec_digit_single_info.resize(14); // x,y,z,time,ADC,sector,stack,layer,row,column,dca,dca_x,dca_y,dca_z
    vec_track_single_info.resize(12); // dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p

    TPM3D_single = new TPolyMarker3D();
    vec_TPM3D_digits.resize(6); // layers
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_digits[i_layer] = new TPolyMarker3D();
    }

    vec_tracklet_fit_points.resize(7); // layers 0-5, 6 = fit through first time cluster points of all layers
    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        vec_tracklet_fit_points[i_layer].resize(2); // start, stop point
        for(Int_t i_start_stop = 0; i_start_stop < 2; i_start_stop++)
        {
            vec_tracklet_fit_points[i_layer][i_start_stop].resize(3); // x,y,z
        }
    }

    TPL3D_helix = new TPolyLine3D();
    fit_line    = new TPolyLine3D();

    //vec_tracklets_line.resize(6);
    //for (Int_t i_layer = 0; i_layer < 6; i_layer++)
    //{
    //    vec_tracklets_line[i_layer] = new TPolyLine3D();
    //}

    // TRD 3D graphics
    fGeo = new AliTRDgeometry;
    geom         = new TGeoManager("geom","My 3D Project");
    vacuum       = new TGeoMaterial("vacuum",0,0,0);
    vacuum       ->SetTransparency(0);
    Air          = new TGeoMedium("Air",0,vacuum);
    top          = geom->MakeBox("top",Air,500,500,500);
    geom         ->SetTopVolume(top);
    geom         ->SetTopVisible(0);
    geom         ->SetVisLevel(4);

    Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
    Fe->SetTransparency(80); // higher value means more transparent, 100 is maximum

    M_outer_tube = new TGeoMaterial("M_outer_tube",55.84,26.7,7.87);
    M_outer_tube->SetTransparency(60); // higher value means more transparent, 100 is maximum

    Material_TRD_box = new TGeoMaterial("Material_TRD_box",55.84,26.7,7.87);
    Material_TRD_box->SetTransparency(70); // higher value means more transparent, 100 is maximum

    Iron                = new TGeoMedium("Iron",1,Fe);
    Me_outer_tube       = new TGeoMedium("Me_outer_tube",1,M_outer_tube);
    Medium_TRD_box      = new TGeoMedium("Medium_TRD_box",1,Material_TRD_box);

    inner_field_tube    = geom->MakeTube("inner_field_tube",Iron,49.5,50.0,510.0/2.0);  // r_min, r_max, dz (half of total length)
    outer_field_tube    = geom->MakeTube("outer_field_tube",Me_outer_tube,556.0/2.0,556.0/2.0+0.5,510.0/2.0);  // r_min, r_max, dz (half of total length)

    inner_field_tube       ->SetLineColor(4);
    outer_field_tube       ->SetLineColor(kRed-8);

    top->AddNodeOverlap(inner_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_tube,1,new TGeoTranslation(0,0,0));

    z_BeamLine    = new TPolyLine3D(2);
    z_BeamLine    ->SetPoint(0,0,0,-550);
    z_BeamLine    ->SetPoint(1,0,0,550);
    z_BeamLine    ->SetLineStyle(0);
    z_BeamLine    ->SetLineColor(kBlue);
    z_BeamLine    ->SetLineWidth(2);

    x_BeamLine    = new TPolyLine3D(2);
    x_BeamLine    ->SetPoint(0,0,0,0);
    x_BeamLine    ->SetPoint(1,50.0,0,0);
    x_BeamLine    ->SetLineStyle(0);
    x_BeamLine    ->SetLineColor(kGreen);
    x_BeamLine    ->SetLineWidth(2);

    y_BeamLine    = new TPolyLine3D(2);
    y_BeamLine    ->SetPoint(0,0,0,0);
    y_BeamLine    ->SetPoint(1,0,50.0,0);
    y_BeamLine    ->SetLineStyle(0);
    y_BeamLine    ->SetLineColor(kRed);
    y_BeamLine    ->SetLineWidth(2);



    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        Int_t flag_not_installed_TRD = 0;
        for(Int_t i_not_installed = 0; i_not_installed < 19; i_not_installed++)
        {
            if(TRD_detector == Not_installed_TRD_detectors[i_not_installed])
            {
                flag_not_installed_TRD = 1;
                break;
            }
        }
        if(flag_not_installed_TRD) continue;


        Int_t TRD_color = kGray+1;
        Int_t flag_defect_TRD = 0;
        for(Int_t i_defect = 0; i_defect < 84; i_defect++)
        {
            if(TRD_detector == Defect_TRD_detectors[i_defect])
            {
                flag_defect_TRD = 1;
                break;
            }
        }
        if(flag_defect_TRD)
        {
            TRD_color = kGreen+1;
            //continue;
        }


        Int_t                TRD_sector         = fGeo         ->GetSector(TRD_detector);
        Int_t                TRD_stack          = fGeo         ->GetStack(TRD_detector);
        Int_t                TRD_layer          = fGeo         ->GetLayer(TRD_detector);
        Float_t              TRD_time0          = fGeo         ->GetTime0(TRD_layer);

        Float_t              TRD_chamber_length = fGeo->GetChamberLength(TRD_layer,TRD_stack);
        Float_t              TRD_chamber_width  = fGeo->GetChamberWidth(TRD_layer);
        Float_t              TRD_chamber_height = 8.4;

        AliTRDpadPlane*      padplane           = fGeo         ->GetPadPlane(TRD_detector);
        Double_t             TRD_col_end        = padplane     ->GetColEnd();
        Double_t             TRD_row_end        = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
        Double_t             TRD_col_start      = padplane     ->GetCol0();
        Double_t             TRD_row_start      = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
        Double_t             TRD_row_end_ROC    = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
        Double_t             TRD_col_spacing    = padplane     ->GetColSpacing();
        Double_t             TRD_row_spacing    = padplane     ->GetRowSpacing();

        Double_t Rotation_angle     = ((360.0/18.0)/2.0) + ((Double_t)TRD_sector)*(360.0/18.0);

        TString HistName = "TRD_box_";
        HistName += TRD_detector;
        TRD_boxes[TRD_detector] = geom->MakeBox(HistName.Data(),Medium_TRD_box,TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0); // dx, dy, dz
        TRD_boxes[TRD_detector]->SetLineColor(TRD_color);
        //TRD_boxes[TRD_detector] = geom->MakeBox(HistName.Data(),Medium_TRD_box,TRD_chamber_width/2.0,1.0,TRD_chamber_length/2.0); // dx, dy, dzd

        Double_t             loc[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glb[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,loc,glb);

        Double_t             locZ[3]           = {TRD_time0-50.0,0.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbZ[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locZ,glbZ);

        Double_t             locX[3]           = {TRD_time0,50.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbX[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locX,glbX);

        Double_t             locY[3]           = {TRD_time0,0.0,50.0+(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbY[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locY,glbY);

        vec_TPL3D_TRD_center[TRD_detector][0] ->SetPoint(0,glb[0],glb[1],glb[2]);
        vec_TPL3D_TRD_center[TRD_detector][0] ->SetPoint(1,glbX[0],glbX[1],glbX[2]);
        vec_TPL3D_TRD_center[TRD_detector][0] ->SetLineColor(kRed);
        vec_TPL3D_TRD_center[TRD_detector][0] ->SetLineWidth(2);
        vec_TPL3D_TRD_center[TRD_detector][0] ->SetLineStyle(1);

        vec_TPL3D_TRD_center[TRD_detector][1] ->SetPoint(0,glb[0],glb[1],glb[2]);
        vec_TPL3D_TRD_center[TRD_detector][1] ->SetPoint(1,glbY[0],glbY[1],glbY[2]);
        vec_TPL3D_TRD_center[TRD_detector][1] ->SetLineColor(kGreen);
        vec_TPL3D_TRD_center[TRD_detector][1] ->SetLineWidth(2);
        vec_TPL3D_TRD_center[TRD_detector][1] ->SetLineStyle(1);

        vec_TPL3D_TRD_center[TRD_detector][2] ->SetPoint(0,glb[0],glb[1],glb[2]);
        vec_TPL3D_TRD_center[TRD_detector][2] ->SetPoint(1,glbZ[0],glbZ[1],glbZ[2]);
        vec_TPL3D_TRD_center[TRD_detector][2] ->SetLineColor(kBlue);
        vec_TPL3D_TRD_center[TRD_detector][2] ->SetLineWidth(2);
        vec_TPL3D_TRD_center[TRD_detector][2] ->SetLineStyle(1);

        vec_TV3_TRD_center[TRD_detector][0].SetXYZ(glbX[0]-glb[0],glbX[1]-glb[1],glbX[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][1].SetXYZ(glbY[0]-glb[0],glbY[1]-glb[1],glbY[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][2].SetXYZ(glbZ[0]-glb[0],glbZ[1]-glb[1],glbZ[2]-glb[2]);

        // Normalize
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_TV3_TRD_center[TRD_detector][i_xyz] *= 1.0/vec_TV3_TRD_center[TRD_detector][i_xyz].Mag();
        }


        combitrans[TRD_detector] = new TGeoCombiTrans();
        combitrans[TRD_detector] ->RotateZ(Rotation_angle + 90.0);
        combitrans[TRD_detector] ->SetTranslation(glb[0],glb[1],glb[2]);
        top->AddNodeOverlap(TRD_boxes[TRD_detector],1,combitrans[TRD_detector]);

    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
vector< vector<TVector3> >  TBase_TRD_Calib::make_clusters(Int_t i_track)
{
    //printf("TBase_TRD_Calib::make_clusters(%d) \n",i_track);
    // Merge the ADC digits according to time bins set in set_merged_time_bins
    // Digit space points are merged using the ADC values as weight

    Int_t N_merged_time_bis = (Int_t)vec_merge_time_bins.size();

    AS_Track     = AS_Event ->getTrack( i_track ); // take the track

    //----------------------------------------------
    // TRD digit information
    UShort_t  fNumTRDdigits = AS_Track ->getNumTRD_digits();

    //printf("i_track: %d, fNumTRDdigits: %d \n",i_track,fNumTRDdigits);

    TVector3 TV3_digit_pos;

    //for old TV3
    vec_TV3_digit_pos_cluster.clear();
    vec_TV3_digit_pos_cluster.resize(6); // 6 layers
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        arr_layer_detector[i_layer] = -1;
        vec_TV3_digit_pos_cluster[i_layer].resize(N_merged_time_bis);
    }

    //for new vector<Double_t>  vector<vector<vector<Double_t>>> vec_Dt_digit_pos_cluster[i_layer][i_merged_time_bin][i_xyzADC]
    vec_Dt_digit_pos_cluster.resize(7); // 6 layers + global line (first time points close to anode wire)
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_Dt_digit_pos_cluster[i_layer].resize(N_merged_time_bis);
    }
    vec_Dt_digit_pos_cluster[6].resize(6); // at maximum 6 layer points for global fit
    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        for(Int_t i_merged_time_bis = 0; i_merged_time_bis < (Int_t)vec_Dt_digit_pos_cluster[i_layer].size(); i_merged_time_bis++)
        {
            vec_Dt_digit_pos_cluster[i_layer][i_merged_time_bis].resize(4);
            for(Int_t i_val = 0; i_val < 4; i_val++)
            {
                vec_Dt_digit_pos_cluster[i_layer][i_merged_time_bis][i_val] = -999.0;
            }
        }
    }


    vector< vector<Double_t> > vec_weight_digits_merged;
    vector< vector< vector<Double_t> > > vec_pos_merge;
    vec_weight_digits_merged.resize(6); // layer
    vec_pos_merge.resize(6); // layer
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_weight_digits_merged[i_layer].resize(N_merged_time_bis);
        vec_pos_merge[i_layer].resize(N_merged_time_bis);
        for(Int_t i_time_merge = 0; i_time_merge < (N_merged_time_bis); i_time_merge++)
        {
            vec_pos_merge[i_layer][i_time_merge].resize(3); // x,y,z
        }
    }

    // Loop over all digits matched to the track
    for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
    {
        //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
        AS_Digit              = AS_Track ->getTRD_digit(i_digits);
        Int_t    layer        = AS_Digit ->get_layer();
        Int_t    sector       = AS_Digit ->get_sector();
        Int_t    column       = AS_Digit ->get_column();
        Int_t    stack        = AS_Digit ->get_stack();
        Int_t    row          = AS_Digit ->get_row();
        Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);
        Float_t  dca_to_track = AS_Digit ->getdca_to_track();
        Float_t  dca_x        = AS_Digit ->getdca_x();
        Float_t  dca_y        = AS_Digit ->getdca_y();
        Float_t  dca_z        = AS_Digit ->getdca_z();
        Float_t  ImpactAngle  = AS_Digit ->getImpactAngle();


        //printf(" \n");
        //printf("----------------------------------------------- \n");
        //printf("i_digits: %d, detector: %d, layer: %d, sector: %d, stack: %d, row: %d, column: %d \n",i_digits,detector,layer,sector,stack,row,column);

        Float_t dca_phi = TMath::Sqrt(dca_x*dca_x + dca_y*dca_y);
        if(dca_phi > 3.0)  continue;

        arr_layer_detector[layer] = detector;

        for(Int_t i_time_merge = 0; i_time_merge < (N_merged_time_bis); i_time_merge++)
        {
            Int_t i_time_start = vec_merge_time_bins[i_time_merge];
            Int_t i_time_stop  = vec_merge_time_bins[i_time_merge + 1];
            
            for(Int_t i_time = i_time_start; i_time < i_time_stop; i_time++)
            {
                Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time) - 10.0;  // baseline correction

                //printf("i_time_merge: %d, i_time: %d, ADC: %4.3f \n",i_time_merge,i_time,ADC);
                if(ADC <= 0.0) continue; // Don't use negative ADC values
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    TV3_digit_pos[i_xyz] = AS_Digit ->get_pos(i_time,i_xyz); // get original digit space point
                    vec_pos_merge[layer][i_time_merge][i_xyz] += ADC*TV3_digit_pos[i_xyz]; // use ADC value as weight
                }
                vec_weight_digits_merged[layer][i_time_merge] += ADC; // keep track of the weights used

                //printf("layer, timebin, ADC: %d, %d, %4.3f \n",layer,i_time_merge,ADC);
                //printf("layer, timebin, ADC+: %d, %d, %4.3f \n",layer,i_time_merge,vec_weight_digits_merged[layer][i_time_merge]);
                //printf("pos: {%4.3f, %4.3f, %4.3f} \n,",digit_pos[0],digit_pos[1],digit_pos[2]);
            }
        }
    }

    // Calculate the average cluster vector for each merged time bin and layer
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        //printf("i_layer: %d \n",i_layer);
        for(Int_t i_time_merge = 0; i_time_merge < (N_merged_time_bis); i_time_merge++)
        {
            //printf("   i_time_merge: %d \n",i_time_merge);
            if(vec_weight_digits_merged[i_layer][i_time_merge] > 0.0)
            {
                for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
                {
                    if(i_xyzADC < 3)  //coordinates
                    {
                        vec_TV3_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC] = vec_pos_merge[i_layer][i_time_merge][i_xyzADC]/vec_weight_digits_merged[i_layer][i_time_merge];
                        vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC]  = vec_pos_merge[i_layer][i_time_merge][i_xyzADC]/vec_weight_digits_merged[i_layer][i_time_merge];
                    }
                    if(i_xyzADC == 3)  //ADC value
                    {
                        vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC] = vec_weight_digits_merged[i_layer][i_time_merge];
                    }
                    //printf("i_xyzADC: %d, value: %4.3f \n",i_xyzADC,vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC]);
                }
                //printf("       pos+ADC: {%4.3f} \n",vec_weight_digits_merged[i_layer][i_time_merge]);
            }
        }
    }

    // Fill the points for the global fit
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
        {
            vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] = vec_Dt_digit_pos_cluster[i_layer][0][i_xyzADC];
        }
    }

    return vec_TV3_digit_pos_cluster;
    //return vec_Dt_digit_pos_cluster;  //to be activated
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//vector<TCanvas*>  TBase_TRD_Calib::make_plots_ADC(Int_t i_track)
void TBase_TRD_Calib::make_plots_ADC(Int_t i_track)
{
    printf("TBase_TRD_Calib::make_plots_ADC(%d) \n",i_track);
    AS_Track               = AS_Event ->getTrack( i_track ); // take the track
    UShort_t fNumTRDdigits = AS_Track ->getNumTRD_digits();   
    Int_t N_time_bin = 24;

    th1f_ADC_vs_time.resize(6);
    ADC_vs_time.resize(6);

    //loop over digits

    Int_t n_digits_per_layer[6] = {0,0,0,0,0,0};

    for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
    {
        //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
        AS_Digit              = AS_Track ->getTRD_digit(i_digits);
        Int_t    layer        = AS_Digit ->get_layer();
        n_digits_per_layer[layer]++;
    }
    
        //cout << "n_digits_per_layer[0]f: " << n_digits_per_layer[0] <<endl;
        //cout << "n_digits_per_layer[1]f: " << n_digits_per_layer[1] <<endl;
        //cout << "n_digits_per_layer[2]f: " << n_digits_per_layer[2] <<endl;
        //cout << "n_digits_per_layer[3]f: " << n_digits_per_layer[3] <<endl;
        //cout << "n_digits_per_layer[4]f: " << n_digits_per_layer[4] <<endl;
        //cout << "n_digits_per_layer[5]f: " << n_digits_per_layer[5] <<endl;

    for (Int_t i_layer = 0; i_layer< 6; i_layer++)
    {
        th1f_ADC_vs_time[i_layer].resize(n_digits_per_layer[i_layer]);
        for(UShort_t i_digits_per_layer = 0; i_digits_per_layer < n_digits_per_layer[i_layer]; i_digits_per_layer++)
        {
            th1f_ADC_vs_time[i_layer][i_digits_per_layer] = new TH1F(Form("th1f_ADC_vs_time_%d_%d",i_layer,i_digits_per_layer),Form("th1f_ADC_vs_time_%d_%d",i_layer,i_digits_per_layer),24,0,23);
        }
    }


    Int_t layer_check  = 0;
    Int_t i_digits_loc = 0;

    for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
    {

        AS_Digit      = AS_Track ->getTRD_digit(i_digits);
        Int_t layer   = AS_Digit ->get_layer();

        if (layer != layer_check)
        {
            i_digits_loc = 0;
        }

        if (layer != layer_check && layer != layer_check+1 && layer_check != 0)
        {
            cout << "!!! layers are not organised properly !!! " << endl;
        }

        //cout << "Test 2" << << endl;

        //cout << "i_digits: " << i_digits << endl;
        //cout << "i_digits_loc: " << i_digits_loc << endl;
        //cout << "layer: " << layer << endl;
        //cout << "layer check: " << layer_check << endl;

        for(Int_t i_time_bin = 0; i_time_bin < N_time_bin; i_time_bin++)
        {
            //cout << "Test 2.2" << endl;

            Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time_bin) - 10.0;  // baseline correction
            //cout << "ADC: " << ADC <<  endl;

            if(ADC <= 0.0) continue;//{th1f_ADC_vs_time[layer][i_digits_loc]->AddBinContent(i_time_bin, 0);} // Don't use negative ADC values
            //cout << "i_digits" << i_digits <<  endl;
            //cout << "i_time_bin" << i_time_bin <<  endl;
            th1f_ADC_vs_time[layer][i_digits_loc]->AddBinContent(i_time_bin+1, ADC);
            //cout << "th1f: " << th1f_ADC_vs_time[layer][i_digits_loc]->GetBinContent(i_time_bin) <<  endl;
            //cout << "i_digits: " << i_digits << endl;
            //cout << "i_digits_loc: " << i_digits_loc << endl;
        }

        layer_check = layer;
        i_digits_loc++; // = i_digits;
        //cout << "Test 4" << endl;
    }

    //cout << "Test 5" << endl;
    Int_t N_pads_y = 2;
    Int_t N_pads_x = 0;
    Int_t N_digits_per_layer[6];

    for (Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        N_digits_per_layer[i_layer] = th1f_ADC_vs_time[i_layer].size();
        N_pads_x = ceil(N_digits_per_layer[i_layer]/N_pads_y)+1;

        ADC_vs_time[i_layer] = new TCanvas(Form("ADC_vs_time_%d",i_layer),Form("ADC_vs_time_%d",i_layer),100,200,1500,620);
        //cout << "Test 7" << endl;
        ADC_vs_time[i_layer] ->SetTopMargin(0.02);
        ADC_vs_time[i_layer] ->SetBottomMargin(0.18);
        ADC_vs_time[i_layer] ->SetRightMargin(0.2);
        ADC_vs_time[i_layer] ->SetLeftMargin(0.2);
        ADC_vs_time[i_layer] ->Divide(N_pads_x,N_pads_y,0.01,0.01);
        //cout << "N_digits_per_layer[i_layer]: " << N_digits_per_layer[i_layer] << endl;
        //cout << "N_pads_x: " << N_pads_x << endl;
        //cout << "N_pads_y: " << N_pads_y << endl;

        for (Int_t i_pad = 0; i_pad < N_digits_per_layer[i_layer]; i_pad++)
        {
            ADC_vs_time[i_layer]->cd(i_pad+1)->SetTicks(1,1);
            //cout << "Test 11" << endl;
            ADC_vs_time[i_layer]->cd(i_pad+1)->SetGrid(0,0);
            ADC_vs_time[i_layer]->cd(i_pad+1)->SetFillColor(10);
            ADC_vs_time[i_layer]->cd(i_pad+1)->SetRightMargin(0.01);
            ADC_vs_time[i_layer]->cd(i_pad+1)->SetTopMargin(0.01);
            //HistName = Form("ADC_vs_time_%d_",i_layer);
            //HistName += i_pad;


            ADC_vs_time[i_layer]->cd(i_pad+1);
            //th1f_ADC_vs_time[i_layer][i_pad]->SetLineColor(color_layer[i_layer]-i_pad);
            th1f_ADC_vs_time[i_layer][i_pad]->SetLineColor(kBlack);
            th1f_ADC_vs_time[i_layer][i_pad]->GetXaxis()->SetTitle("Time bin");
            th1f_ADC_vs_time[i_layer][i_pad]->GetYaxis()->SetTitle("ADC counts");
            th1f_ADC_vs_time[i_layer][i_pad]->Draw();

            //cout << "Test 12" << endl;
            //cout << "th1f: " << th1f_ADC_vs_time[i_layer][i_pad]->GetBinContent(0) <<  endl;
        }
    }
}

//----------------------------------------------------------------------------------------
TPolyLine3D* TBase_TRD_Calib::get_straight_line_fit(Int_t i_track)
{
    printf("TBase_TRD_Calib::get_straight_line_fit((%d) \n",i_track);

    //fit merged digits with a straight line

    TGraph2D*    tg_cluster_points              = new TGraph2D();
    TPolyLine3D* digits_fit_line = new TPolyLine3D();

    // Fill the 2D graph
    Double_t p0[4] = {10,20,1,2};

    // generate graph with the 3d points

    Int_t i_layer_notempty = 0;

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        //printf("i_layer: %d \n",i_layer);
        if(vec_Dt_digit_pos_cluster[i_layer][0][0] != 0 && vec_Dt_digit_pos_cluster[i_layer][0][1] != 0 && vec_Dt_digit_pos_cluster[i_layer][0][2] != 0)
        {
            tg_cluster_points->SetPoint(i_layer_notempty,vec_Dt_digit_pos_cluster[i_layer][0][0],vec_Dt_digit_pos_cluster[i_layer][0][1],vec_Dt_digit_pos_cluster[i_layer][0][2]);
            //dt->SetPointError(N,0,0,err);
            //Double_t* point = tg_cluster_points->GetX();
            //cout << "layer: " <<  i_layer  << endl;
            //cout << "layer not empty : " <<  i_layer_notempty  << endl;
            //cout << "point: " <<  point[i_layer_notempty]  << endl;
            i_layer_notempty++;
        }
    }

    if(i_layer_notempty == 0)
    {
        printf("No digits found for this track \n");
        return digits_fit_line;
    }

    // fit the graph now

    TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
    min->SetObjectFit(tg_cluster_points);
    //min->SetFCN(SumDistance2);
    min->SetFCN(SumDistance2_X);


    Double_t arglist[10];
    arglist[0] = 3;
    //min->ExecuteCommand("SET PRINT",arglist,1);

    Double_t a0[3] = {0,0,0};
    Double_t a1[3] = {0,0,0};

    tg_cluster_points -> GetPoint(0,a0[0],a0[1],a0[2]);
    tg_cluster_points -> GetPoint(i_layer_notempty-1,a1[0],a1[1],a1[2]);

    printf("point start: {%4.3f, %4.3f, %4.3f} \n",a0[0],a0[1],a0[2]);
    printf("point end: {%4.3f, %4.3f, %4.3f} \n",a1[0],a1[1],a1[2]);

#if 0
    TPolyLine3D* digits_fit_line_init = new TPolyLine3D();
    digits_fit_line_init ->SetNextPoint(a0[0],a0[1],a0[2]);
    digits_fit_line_init ->SetNextPoint(a1[0],a1[1],a1[2]);
    digits_fit_line_init ->SetLineColor(kMagenta);
    digits_fit_line_init ->SetLineWidth(3);
    digits_fit_line_init ->DrawClone("ogl");
#endif


    TVector3 vec_a0;
    vec_a0.SetXYZ(a0[0],a0[1],a0[2]);
    TVector3 vec_a1;
    vec_a1.SetXYZ(a1[0],a1[1],a1[2]);
    TVector3 vec_u = vec_a1 - vec_a0;
    TVector3 vec_u_perp;
    Double_t pStart[4]; //= {1,1,1,1};
#if 0
    //---------------------------------
    // z-plane
    TVector3 vec_x0 = vec_a0 - vec_u*(vec_a0.Z()/vec_u.Z());
    vec_u_perp.SetXYZ(vec_u[0],vec_u[1],vec_u[2]);
    vec_u_perp *= 1.0/vec_u[2];
    TVector3 vec_x1 = vec_x0 + vec_u_perp;
    pStart[0] = vec_x0.X();
    pStart[1] = vec_x1.X() - pStart[0];
    pStart[2] = vec_x0.Y();
    pStart[3] = vec_x1.Y() - pStart[2];
    //---------------------------------
#endif


#if 1
    //---------------------------------
    // x-plane
    TVector3 vec_x0 = vec_a0 - vec_u*(vec_a0.X()/vec_u.X());
    vec_u_perp.SetXYZ(vec_u[0],vec_u[1],vec_u[2]);
    vec_u_perp *= 1.0/vec_u[0];
    TVector3 vec_x1 = vec_x0 + vec_u_perp;
    pStart[0] = vec_x0.Z();
    pStart[1] = vec_x1.Z() - pStart[0];
    pStart[2] = vec_x0.Y();
    pStart[3] = vec_x1.Y() - pStart[2];
    //---------------------------------
#endif


    //TVector3 a0 =

    //TVector3 u = a1-a0;
    //x0 = a0  (a0z/uz)*u;
    //x1 = x0 + u / uz;


   

    cout << "pStart[0]" << pStart[0] << endl;
    cout << "pStart[1]" << pStart[1] << endl;
    cout << "pStart[2]" << pStart[2] << endl;
    cout << "pStart[3]" << pStart[3] << endl;

    min->SetParameter(0,"x0",pStart[0],0.01,0,0);
    min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
    min->SetParameter(2,"y0",pStart[2],0.01,0,0);
    min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

    arglist[0] = 1000; // number of function calls
    arglist[1] = 0.001; // tolerance
    min->ExecuteCommand("MIGRAD",arglist,2);

    //if (minos) min->ExecuteCommand("MINOS",arglist,0);
    int nvpar,nparx;
    Double_t amin,edm, errdef;
    min->GetStats(amin,edm,errdef,nvpar,nparx);
    min->PrintResults(1,amin);
    //tg_cluster_points->Draw("p0");

    // get fit parameters
    Double_t parFit[4];
    for(int i = 0; i <4; ++i)
    {
        parFit[i] = min->GetParameter(i);
        //parFit[i] = pStart[i];
    }


    // draw the fitted line
    int n = 1000;
    Double_t t0 = -500.0;
    Double_t dt = 1;
    TVector3 TV3_line_point;
    Int_t i_point = 0;
    for(int i = 0; i <n; ++i)
    {
        Double_t t = t0 + dt*i;
        Double_t x,y,z;
        //line(t,parFit,x,y,z);
        line_X(t,parFit,x,y,z);
        TV3_line_point.SetXYZ(x,y,z);
        //printf("point: {%4.3f, %4.3f, %4.3f} \n",x,y,z);

        Double_t distance = 1000.0;
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            TVector3 vec_TV3_digit_pos_cluster;
            vec_TV3_digit_pos_cluster.SetXYZ(vec_Dt_digit_pos_cluster[i_layer][0][0],vec_Dt_digit_pos_cluster[i_layer][0][1],vec_Dt_digit_pos_cluster[i_layer][0][2]);
            TVector3 TV3_diff = vec_TV3_digit_pos_cluster - TV3_line_point;
            Double_t distance_layer = TV3_diff.Mag();
            if(distance_layer < distance) distance = distance_layer;
        }
        if(TV3_line_point.Perp() > 300.0 && TV3_line_point.Perp() < 380.0 && distance < 10.0)
        {
            digits_fit_line->SetPoint(i_point,x,y,z);
            //printf("point line: {%4.3f, %4.3f, %4.3f} \n",x,y,z);
            i_point++;
        }
    }
    digits_fit_line->SetLineColor(kRed);
    //digits_fit_line->Draw("same");

    return digits_fit_line;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_line(Int_t i_track)
{
    fit_line = get_straight_line_fit(i_track);

    fit_line    ->SetLineStyle(0);
    fit_line    ->SetLineColor(kRed);
    fit_line    ->SetLineWidth(2);
    //fit_line    ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::get_tracklets_fit(Int_t i_track)
{
    // Is fitting the tracklets and doing the global fit through all first cluster points of all available layers
    //printf("TBase_TRD_Calib::get_tracklets_fit((%d) \n",i_track);

    //fit merged digits with a straight line
    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        for(Int_t i_start_stop = 0; i_start_stop < 2; i_start_stop++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                vec_tracklet_fit_points[i_layer][i_start_stop][i_xyz] = -999.0;
            }
        }
    }

    Double_t p0[4] = {10,20,1,2};
    Int_t i_layer_notempty = 0;

    //declare things that will be needed later
    Double_t arglist[10];
    Double_t a0[3] = {0,0,0};
    Double_t a1[3] = {0,0,0};

    TVector3 vec_a0;
    TVector3 vec_a1;
    TVector3 vec_u;
    TVector3 vec_x0;
    TVector3 vec_u_perp;
    TVector3 vec_x1;

    Double_t pStart[4]; //= {1,1,1,1};

    int nvpar,nparx;
    Double_t amin,edm, errdef;
    
    Double_t parFit[4];

    Int_t    n;
    Double_t t0;
    Double_t dt;
    TVector3 TV3_line_point;
    Double_t i_point;

    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        vec_layer_in_fit[i_layer] = 0;
    }

    //loop over layers
    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        global_layer = i_layer;

#if 0
        for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_Dt_digit_pos_cluster[i_layer].size(); i_time_merge++)
        {
              printf("layer: %d, i_time_merge: %d, point: {%4.3f, %4.3f, %4.3f} \n",i_layer,i_time_merge,vec_Dt_digit_pos_cluster[i_layer][i_time_merge][0],vec_Dt_digit_pos_cluster[i_layer][i_time_merge][1],vec_Dt_digit_pos_cluster[i_layer][i_time_merge][2]);
        }
#endif

        Int_t i_time_merge_AB[2] = {-1,-1};
        for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_Dt_digit_pos_cluster[i_layer].size(); i_time_merge++)
        {
            if(vec_Dt_digit_pos_cluster[i_layer][i_time_merge][0] == -999.0 && vec_Dt_digit_pos_cluster[i_layer][i_time_merge][1] == -999.0 && vec_Dt_digit_pos_cluster[i_layer][i_time_merge][2] == -999.0) continue;
            else i_time_merge_AB[0] = i_time_merge;
        }
        if(i_time_merge_AB[0] == -1) continue; // all values are 0

        for(Int_t i_time_merge = ((Int_t)vec_Dt_digit_pos_cluster[i_layer].size() - 1); i_time_merge >= 0; i_time_merge--)
        {
            if(vec_Dt_digit_pos_cluster[i_layer][i_time_merge][0] == -999.0 && vec_Dt_digit_pos_cluster[i_layer][i_time_merge][1] == -999.0 && vec_Dt_digit_pos_cluster[i_layer][i_time_merge][2] == -999.0) continue;
            else i_time_merge_AB[1] = i_time_merge;
        }

        if(i_time_merge_AB[0] == i_time_merge_AB[1]) continue; // no fit possible with just one point

        //printf("TBase_TRD_Calib::get_tracklets_fit(%d), i_layer: %d \n",i_track,i_layer);


        TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
        //min->SetObjectFit(tracklets_gr);

        arglist[0] = 3;
        //min->ExecuteCommand("SET PRINT",arglist,1);
        Double_t arglist_B[1] = {-1};
        min->ExecuteCommand("SET PRIntout",arglist_B,1);
        min->ExecuteCommand("SET NOWarnings",arglist_B,1);

        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            a0[i_xyz] = vec_Dt_digit_pos_cluster[i_layer][i_time_merge_AB[0]][i_xyz];
            a1[i_xyz] = vec_Dt_digit_pos_cluster[i_layer][i_time_merge_AB[1]][i_xyz];
        }

        //printf("point start: {%4.3f, %4.3f, %4.3f} \n",a0[0],a0[1],a0[2]);
        //printf("point end: {%4.3f, %4.3f, %4.3f} \n",a1[0],a1[1],a1[2]);

#if 0
        TPolyLine3D* digits_fit_line_init = new TPolyLine3D();
        digits_fit_line_init ->SetNextPoint(a0[0],a0[1],a0[2]);
        digits_fit_line_init ->SetNextPoint(a1[0],a1[1],a1[2]);
        digits_fit_line_init ->SetLineColor(kMagenta);
        digits_fit_line_init ->SetLineWidth(3);
        digits_fit_line_init ->DrawClone("ogl");
#endif

        Int_t flag_XZ = 1;
        if(fabs(a0[0] - a1[0]) > fabs(a0[2] - a1[2]))
        {
            flag_XZ = 0; // x is used
        }

        vec_a0.SetXYZ(a0[0],a0[1],a0[2]);
        vec_a1.SetXYZ(a1[0],a1[1],a1[2]);
        vec_u = vec_a1 - vec_a0;
        if(flag_XZ == 0)
        {
            min->SetFCN(SumDistance2_X_tr);
            vec_x0 = vec_a0 - vec_u*(vec_a0.X()/vec_u.X());
            vec_u_perp.SetXYZ(vec_u[0],vec_u[1],vec_u[2]);
            vec_u_perp *= 1.0/vec_u[0];
        }
        else
        {
            min->SetFCN(SumDistance2_tr);
            vec_x0 = vec_a0 - vec_u*(vec_a0.Z()/vec_u.Z());
            vec_u_perp.SetXYZ(vec_u[0],vec_u[1],vec_u[2]);
            vec_u_perp *= 1.0/vec_u[2];
        }
        vec_x1 = vec_x0 + vec_u_perp;

        //TVector3 u = a1-a0;
        //x0 = a0  (a0z/uz)*u;
        //x1 = x0 + u / uz;

        if(flag_XZ == 0)
        {
            pStart[0] = vec_x0.Z();
            pStart[1] = vec_x1.Z() - pStart[0];
            pStart[2] = vec_x0.Y();
            pStart[3] = vec_x1.Y() - pStart[2];
        }
        else
        {
            pStart[0] = vec_x0.X();
            pStart[1] = vec_x1.X() - pStart[0];
            pStart[2] = vec_x0.Y();
            pStart[3] = vec_x1.Y() - pStart[2];
        }

        //cout << "pStart[0]" << pStart[0] << endl;
        //cout << "pStart[1]" << pStart[1] << endl;
        //cout << "pStart[2]" << pStart[2] << endl;
        //cout << "pStart[3]" << pStart[3] << endl;

        min->SetParameter(0,"x0",pStart[0],0.01,0,0);
        min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
        min->SetParameter(2,"y0",pStart[2],0.01,0,0);
        min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

        arglist[0] = 1000; // number of function calls
        arglist[1] = 0.001; // tolerance
        //printf("     ------------------ MIGRAD ------------------ \n");
        min->ExecuteCommand("MIGRAD",arglist,2);
        //printf("     ------------------ END ------------------ \n");

        //if (minos) min->ExecuteCommand("MINOS",arglist,0);

        min->GetStats(amin,edm,errdef,nvpar,nparx);
        //min->PrintResults(1,amin);

        // get fit parameters
        for(int i = 0; i < 4; ++i)
        {
            parFit[i] = min->GetParameter(i);
            //parFit[i] = pStart[i];
        }


        // draw the fitted line
        n = 1000;
        t0 = -500.0;
        dt = 1;

        Double_t x_A, y_A, z_A;
        t0 = ((a0[1]-parFit[2])/parFit[3])-500;
        if(flag_XZ == 0) line_X(t0,parFit,x_A,y_A,z_A);
        else line(t0,parFit,x_A,y_A,z_A);
        //printf("start point: {%4.3f, %4.3f, %4.3f} \n",x_A,y_A,z_A);

        i_point = 0;
        TVector3 vec_AB[2];
        TVector3 vec_AB_2D[2];
        vec_AB[0].SetXYZ(-999.0,-999.0,-999.0);
        vec_AB[1].SetXYZ(-999.0,-999.0,-999.0);
        vec_AB_2D[0].SetXYZ(-999.0,-999.0,-999.0);
        vec_AB_2D[1].SetXYZ(-999.0,-999.0,-999.0);
        for(int i = 0; i <n; ++i)
        {
            Double_t t = t0 + dt*i;
            Double_t x,y,z;
            if(flag_XZ == 0) line_X(t,parFit,x,y,z);
            else line(t,parFit,x,y,z);
            TV3_line_point.SetXYZ(x,y,z);
            //printf("point: {%4.3f, %4.3f, %4.3f} \n",x,y,z);

            Double_t distance = 1000.0;
            for(Int_t i_time_bin = 0; i_time_bin < (Int_t)vec_Dt_digit_pos_cluster[i_layer].size(); i_time_bin++)
            {
                TVector3 vec_TV3_digit_pos_cluster;
                vec_TV3_digit_pos_cluster.SetXYZ(vec_Dt_digit_pos_cluster[i_layer][i_time_bin][0],vec_Dt_digit_pos_cluster[i_layer][i_time_bin][1],vec_Dt_digit_pos_cluster[i_layer][i_time_bin][2]);
                TVector3 TV3_diff = vec_TV3_digit_pos_cluster - TV3_line_point;
                Double_t distance_layer = TV3_diff.Mag();
                if(distance_layer < distance) distance = distance_layer;
            }
            //printf("perp: %4.3f, distance: %4.3f \n",TV3_line_point.Perp(),distance);
            if(TV3_line_point.Perp() > 270.0 && TV3_line_point.Perp() < 380.0 && distance < 100.0)  // 100 just in case
            {
                if(i_point == 0)
                {
                    vec_AB[0].SetXYZ(x,y,z);
                    vec_AB_2D[0].SetXYZ(x,y,0.0);
                }
                else
                {
                    vec_AB[1].SetXYZ(x,y,z);
                    vec_AB_2D[1].SetXYZ(x,y,0.0);
                }
                //printf("i_layer: %d, i: %d, point line: {%4.3f, %4.3f, %4.3f} \n",i_layer,i,x,y,z);
                i_point++;
            }
        }

        Int_t order_AB[2] = {0,1};
        if(vec_AB_2D[0].Mag() > vec_AB_2D[1].Mag())
        {
            order_AB[0] = 1;
            order_AB[1] = 0;
        }

        for(Int_t i_AB = 0; i_AB < 2; i_AB++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                vec_tracklet_fit_points[i_layer][i_AB][i_xyz] = vec_AB[order_AB[i_AB]][i_xyz];
            }
        }


        vec_layer_in_fit[i_layer_notempty] = i_layer;
        i_layer_notempty++;

        delete min;

    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_tracklets_line(Int_t i_track)
{
    get_tracklets_fit(i_track);

    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        if(vec_tracklet_fit_points[i_layer][0][0] > -999.0 && vec_tracklet_fit_points[i_layer][1][0] > -999.0)
        {
            vec_tracklets_line ->SetPoint(0,vec_tracklet_fit_points[i_layer][0][0],vec_tracklet_fit_points[i_layer][0][1],vec_tracklet_fit_points[i_layer][0][2]);
            vec_tracklets_line ->SetPoint(1,vec_tracklet_fit_points[i_layer][1][0],vec_tracklet_fit_points[i_layer][1][1],vec_tracklet_fit_points[i_layer][1][2]);
            //printf("i_layer_notemtpy: %d, layer: %d \n",i_layer_notempty,vec_layer_in_fit[i_layer_notempty]);
            vec_tracklets_line ->SetLineStyle(1);
            vec_tracklets_line ->SetLineColor(color_layer[i_layer]);
            vec_tracklets_line ->SetLineWidth(line_width_layer[i_layer]);
            vec_tracklets_line ->DrawClone("ogl");
            //printf("i_layer: %d \n", i_layer);
        }
        else
        {
            printf("TBase_TRD_Calib::Draw_tracklets_line(%d), i_layer: %d has no entry \n",i_track,i_layer);
        }
    }
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
TPolyLine3D* TBase_TRD_Calib::get_helix_polyline(Int_t i_track)
{
    TPolyLine3D* TPL3D_helix_track = new TPolyLine3D();
    AS_Track      = AS_Event ->getTrack( i_track ); // take the track
    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
        aliHelix.fHelix[i_param] = AS_Track ->getHelix_param(i_param);
    }

    Double_t helix_point[3];
    Double_t pathA = 0.0;
    Double_t radius = 0.0;
    for(Int_t i_step = 0; i_step < 400; i_step++)
    {
        pathA = i_step*3.0;
        aliHelix.Evaluate(pathA,helix_point);

        radius = TMath::Sqrt(TMath::Power(helix_point[0],2.0) + TMath::Power(helix_point[1],2.0));
        TPL3D_helix_track ->SetPoint(i_step,helix_point[0],helix_point[1],helix_point[2]);
        //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_step,helix_point[0],helix_point[1],helix_point[2]);
        if(radius > 368.0) break;
    }

    return TPL3D_helix_track;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_track(Int_t i_track)
{
    TPL3D_helix = get_helix_polyline(i_track);

    TPL3D_helix    ->SetLineStyle(0);
    TPL3D_helix    ->SetLineColor(kRed);
    TPL3D_helix    ->SetLineWidth(5);
    TPL3D_helix    ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_neighbor_tracks(Int_t i_track_sel)
{
    // Draws all tracks with hits in the same TRD module(s)

    UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event
    AS_Track      = AS_Event ->getTrack( i_track_sel ); // take the track
    UShort_t  fNumTRDdigits = AS_Track ->getNumTRD_digits();
    vector<Int_t> vec_detectors_hit;

    // Get the detectors which were hit from the track of interest
    for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
    {
        //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
        AS_Digit              = AS_Track ->getTRD_digit(i_digits);
        Int_t    layer        = AS_Digit ->get_layer();
        Int_t    sector       = AS_Digit ->get_sector();
        Int_t    stack        = AS_Digit ->get_stack();
        Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);

        Int_t flag_exist = 0;
        for(Int_t i_ele = 0; i_ele < (Int_t)vec_detectors_hit.size(); i_ele++)
        {
            if(vec_detectors_hit[i_ele] == detector)
            {
                flag_exist = 1;
                break;
            }
        }
        if(!flag_exist) vec_detectors_hit.push_back(detector);
    }


    for(Int_t i_ele = 0; i_ele < (Int_t)vec_detectors_hit.size(); i_ele++)
    {
        //printf("detector hit: %d \n",vec_detectors_hit[i_ele]);
    }


    // Loop over all tracks and find those with hits in vec_detectors_hit
    for(UShort_t i_track = 0; i_track < NumTracks; ++i_track) // loop over all tracks of the actual event
    {
        if(i_track == i_track_sel) continue; // Don't use the selected track twice
        AS_Track      = AS_Event ->getTrack( i_track ); // take the track
        fNumTRDdigits = AS_Track ->getNumTRD_digits();

        for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
        {
            //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
            AS_Digit              = AS_Track ->getTRD_digit(i_digits);
            Int_t    layer        = AS_Digit ->get_layer();
            Int_t    sector       = AS_Digit ->get_sector();
            Int_t    stack        = AS_Digit ->get_stack();
            Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);

            Int_t flag_exist = 0;
            for(Int_t i_ele = 0; i_ele < (Int_t)vec_detectors_hit.size(); i_ele++)
            {
                if(vec_detectors_hit[i_ele] == detector)
                {
                    flag_exist = 1;
                    break;
                }
            }
            if(flag_exist)
            {
                printf("Added track: %d to neighbor tracks \n",i_track);
                vec_TPL3D_helix_neighbor.push_back(get_helix_polyline(i_track));
                break;
            }
        }
    }

    for(Int_t i_track_neighbor = 0; i_track_neighbor < (Int_t)vec_TPL3D_helix_neighbor.size(); i_track_neighbor++)
    {
        vec_TPL3D_helix_neighbor[i_track_neighbor]  ->SetLineStyle(0);
        vec_TPL3D_helix_neighbor[i_track_neighbor]  ->SetLineColor(kGreen+2);
        vec_TPL3D_helix_neighbor[i_track_neighbor]  ->SetLineWidth(2);
        vec_TPL3D_helix_neighbor[i_track_neighbor]  ->DrawClone("ogl");
    }

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TGLViewer* TBase_TRD_Calib::Draw_TRD()
{
    top->DrawClone("ogl");
    TGL_viewer = (TGLViewer *)gPad->GetViewer3D();
    TGL_viewer ->SetClearColor(kBlack);

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        if(!(i_det % 30 == 0)) continue;
        //if(!(i_det == 0)) continue;
        vec_TPL3D_TRD_center[i_det][0] ->DrawClone("ogl");
        vec_TPL3D_TRD_center[i_det][1] ->DrawClone("ogl");
        vec_TPL3D_TRD_center[i_det][2] ->DrawClone("ogl");
    }

#if 0
    x_BeamLine    ->DrawClone("ogl");
    y_BeamLine    ->DrawClone("ogl");
    z_BeamLine    ->DrawClone("ogl");
#endif

    return TGL_viewer;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Init_tree(TString SEList)
{
    cout << "Initialize tree" << endl;
    TString pinputdir = "/misc/alidata120/alice_u/schmah/TRD_offline_calib/Data/";
    //TString pinputdir = "/home/ceres/berdnikova/TRD-Run3-Calibration/";

    AS_Event = new Ali_AS_Event();
    AS_Track = new Ali_AS_Track();
    AS_Tracklet = new Ali_AS_Tracklet();
    AS_Digit = new Ali_AS_TRD_digit();

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

    file_entries_total = input_SE->GetEntries();
    N_Events = file_entries_total;
    cout << "Total number of events in tree: " << file_entries_total << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t TBase_TRD_Calib::Loop_event(Long64_t event)
{
    printf("Loop event number: %lld \n",event);

    if (!input_SE->GetEntry( event )) return 0; // take the event -> information is stored in event

    N_Digits = 0;

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        delete vec_TPM3D_digits[i_layer];
        vec_TPM3D_digits[i_layer] = new TPolyMarker3D();
    }

    vec_TPL3D_tracklets.clear();

    //---------------------------------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event
    Double_t EventVertexX         = AS_Event ->getx();
    Double_t EventVertexY         = AS_Event ->gety();
    Double_t EventVertexZ         = AS_Event ->getz();
    Int_t    N_tracks_event       = AS_Event ->getN_tracks();
    Int_t    N_TRD_tracklets      = AS_Event ->getN_TRD_tracklets();
    Int_t    N_TRD_tracklets_offline = AS_Event ->getNumTracklets();
    Float_t  V0MEq                = AS_Event ->getcent_class_V0MEq();

    N_Tracks = NumTracks;

    vec_track_info.clear();
    vec_digit_track_info.clear();
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_Tracklet_pos.clear();
        vec_TV3_Tracklet_dir.clear();
    }

    // Loop over all tracklets
    Double_t scale_factor_length = 2.0;
    for(UShort_t i_tracklet = 0; i_tracklet < N_TRD_tracklets_offline; ++i_tracklet) // loop over all tracklets of the actual event
    {
        AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
        TVector3 TV3_offset     = AS_Tracklet ->get_TV3_offset();
        TVector3 TV3_dir        = AS_Tracklet ->get_TV3_dir();
        Short_t  i_det_tracklet = AS_Tracklet ->get_detector();

        vec_TV3_Tracklet_pos[i_det_tracklet].push_back(TV3_offset);
        vec_TV3_Tracklet_dir[i_det_tracklet].push_back(TV3_dir);

        Double_t impact_angle = TV3_dir.Angle(vec_TV3_TRD_center[i_det_tracklet][2]);
        if(impact_angle > TMath::Pi()*0.5) impact_angle -= TMath::Pi();

        //if(!(i_det_tracklet >= 114 && i_det_tracklet <= 119)) continue;
        //if(!(i_det_tracklet >= 378 && i_det_tracklet <= 382)) continue;

        Float_t points[6] =
        {
            (Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),
            (Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2])
        };
        printf("i_tracklet: %d, out of %d, impact_angle: %4.3f, offset: {%4.3f, %4.3f, %4.3f}, end: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,N_TRD_tracklets_offline,impact_angle*TMath::RadToDeg(),TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + TV3_dir[0],TV3_offset[1] + TV3_dir[1],TV3_offset[2] + TV3_dir[2]);

        vec_TPL3D_tracklets.push_back(new TPolyLine3D(2,points));
    }

    // Loop over all tracks
    for(UShort_t i_track = 0; i_track < NumTracks; ++i_track) // loop over all tracks of the actual event
    {
        //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
        AS_Track      = AS_Event ->getTrack( i_track ); // take the track
        Double_t nsigma_TPC_e   = AS_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = AS_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = AS_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = AS_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = AS_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = AS_Track ->getTRDSignal();
        Double_t TRDsumADC      = AS_Track ->getTRDsumADC();
        Double_t dca            = AS_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = AS_Track ->get_TLV_part();
        UShort_t NTPCcls        = AS_Track ->getNTPCcls();
        UShort_t NTRDcls        = AS_Track ->getNTRDcls();
        UShort_t NITScls        = AS_Track ->getNITScls();
        Float_t TPCchi2         = AS_Track ->getTPCchi2();
        Float_t TPCdEdx         = AS_Track ->getTPCdEdx();
        Float_t TOFsignal       = AS_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = AS_Track ->getTrack_length();

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();


        vec_track_single_info[0]  = dca;
        vec_track_single_info[1]  = TPCdEdx;
        vec_track_single_info[2]  = momentum;
        vec_track_single_info[3]  = eta_track;
        vec_track_single_info[4]  = pT_track;
        vec_track_single_info[5]  = TOFsignal;
        vec_track_single_info[6]  = Track_length;
        vec_track_single_info[7]  = TRDsumADC;
        vec_track_single_info[8]  = TRD_signal;
        vec_track_single_info[9]  = nsigma_TPC_e;
        vec_track_single_info[10] = nsigma_TPC_pi;
        vec_track_single_info[11] = nsigma_TPC_p;
        vec_track_info.push_back(vec_track_single_info);

        //----------------------------------------------
        // TRD digit information
        UShort_t  fNumTRDdigits = AS_Track ->getNumTRD_digits();

        //printf("i_track: %d, fNumTRDdigits: %d \n",i_track,fNumTRDdigits);

        vec_digit_info.clear();
        for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
        {
            //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
            AS_Digit              = AS_Track ->getTRD_digit(i_digits);
            Int_t    layer        = AS_Digit ->get_layer();
            Int_t    sector       = AS_Digit ->get_sector();
            Int_t    column       = AS_Digit ->get_column();
            Int_t    stack        = AS_Digit ->get_stack();
            Int_t    row          = AS_Digit ->get_row();
            Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);
            Float_t  dca_to_track = AS_Digit ->getdca_to_track();
            Float_t  dca_x        = AS_Digit ->getdca_x();
            Float_t  dca_y        = AS_Digit ->getdca_y();
            Float_t  dca_z        = AS_Digit ->getdca_z();
            Float_t  ImpactAngle  = AS_Digit ->getImpactAngle();

            for(Int_t i_time = 0; i_time < 24; i_time++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    digit_pos[i_xyz] = AS_Digit ->get_pos(i_time,i_xyz);
                }
                //printf("track: %d/%d, digit: %d/%d \n",i_track,NumTracks,i_digits,fNumTRDdigits);
                //printf("pos: {%4.3f, %4.3f, %4.3f} \n,",digit_pos[0],digit_pos[1],digit_pos[2]);
                vec_TPM3D_digits[layer] ->SetNextPoint(digit_pos[0],digit_pos[1],digit_pos[2]);

                Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time);


                // x,y,z,time,ADC,sector,stack,layer,row,column,dca
                vec_digit_single_info[0]  = digit_pos[0];
                vec_digit_single_info[1]  = digit_pos[1];
                vec_digit_single_info[2]  = digit_pos[2];
                vec_digit_single_info[3]  = i_time;
                vec_digit_single_info[4]  = ADC;
                vec_digit_single_info[5]  = sector;
                vec_digit_single_info[6]  = stack;
                vec_digit_single_info[7]  = layer;
                vec_digit_single_info[8]  = row;
                vec_digit_single_info[9]  = column;
                vec_digit_single_info[10] = dca_to_track;
                vec_digit_single_info[11] = dca_x;
                vec_digit_single_info[12] = dca_y;
                vec_digit_single_info[13] = dca_z;

                //printf("dca_full: %4.3f, dca: {%4.3f, %4.3f, %4.3f} \n",dca_to_track,dca_x,dca_y,dca_z);
                vec_digit_info.push_back(vec_digit_single_info);

                N_Digits ++;
            }

        } // end of digits loop

        vec_digit_track_info.push_back(vec_digit_info);

    } // end of track loop

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Track_Tracklets()
{
    printf("TBase_TRD_Calib::Track_Tracklets() \n");

    Double_t scale_factor_length = 8.0;
    vec_TPL3D_tracklets_match.clear();

    for(Int_t i_sector = 0; i_sector < 18; i_sector++)
    {
        for(Int_t i_stack = 0; i_stack < 5; i_stack++)
        {
            for(Int_t i_layer = 5; i_layer >= 1; i_layer--)
            {
                Int_t i_detector = i_layer + 6*i_stack + 6*5*i_sector;

                if(!(i_detector >= 378 && i_detector <= 382)) continue;

                Int_t N_matches = 0;
                for(Int_t i_tracklet = 0; i_tracklet < (Int_t)vec_TV3_Tracklet_pos[i_detector].size(); i_tracklet++)
                {
                    TVector3 TV3_pos = vec_TV3_Tracklet_pos[i_detector][i_tracklet];
                    TVector3 TV3_dir = vec_TV3_Tracklet_dir[i_detector][i_tracklet];
                    Double_t dir_length = TV3_dir.Mag();
                    if(dir_length > 0.0) TV3_dir *= 1.0/dir_length;

                    for(Int_t i_layer_B = (i_layer - 1); i_layer_B >= 0; i_layer_B--)
                    {
                        Int_t i_detector_B = i_layer_B + 6*i_stack + 6*5*i_sector;
                        for(Int_t i_tracklet_B = 0; i_tracklet_B < (Int_t)vec_TV3_Tracklet_pos[i_detector_B].size(); i_tracklet_B++)
                        {
                            TVector3 TV3_pos_B = vec_TV3_Tracklet_pos[i_detector_B][i_tracklet_B];
                            TVector3 TV3_dir_B = vec_TV3_Tracklet_dir[i_detector_B][i_tracklet_B];
                            Double_t dir_length_B = TV3_dir_B.Mag();
                            if(dir_length_B > 0.0) TV3_dir_B *= 1.0/dir_length_B;
                            Double_t distance = calculateMinimumDistanceStraightToPoint(TV3_pos,TV3_dir,TV3_pos_B);
                            Double_t projection = TV3_dir.Dot(TV3_dir_B);

                            if(distance < 15.0 && fabs(1.0-projection) < 0.1)
                            {
                                N_matches++;
                                Float_t points[6] =
                                {
                                    (Float_t)(TV3_pos_B[0]),(Float_t)(TV3_pos_B[1]),(Float_t)(TV3_pos_B[2]),
                                    (Float_t)(TV3_pos_B[0] + scale_factor_length*TV3_dir_B[0]),(Float_t)(TV3_pos_B[1] + scale_factor_length*TV3_dir_B[1]),(Float_t)(TV3_pos_B[2] + scale_factor_length*TV3_dir_B[2])
                                };

                                vec_TPL3D_tracklets_match.push_back(new TPolyLine3D(2,points));
                                printf("Match, detector: %d, detector B: %d, distance: %4.3f, projection: %4.3f, length: %4.3f \n",i_detector,i_detector_B,distance,projection,TV3_dir_B.Mag());
                            }
                        }
                    }
                }
                printf("     ----> N_matches: %d \n",N_matches);
            }
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Calibrate()
{
    printf("TBase_TRD_Calib::Calibrate() \n");

    vec_tp_Delta_vs_impact.resize(540);

    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_%d",i_det),Form("vec_th1d_Delta_vs_impact_%d",i_det),360,-360,360);
    }

    //for(Long64_t i_event = 0; i_event < file_entries_total; i_event++)
    for(Long64_t i_event = 0; i_event < 50000; i_event++)
    {
        printf("i_event: %lld out of %lld \n",i_event,file_entries_total);
        if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

        UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event

        for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        {
            //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
            AS_Track      = AS_Event ->getTrack( i_track ); // take the track
            Double_t nsigma_TPC_e   = AS_Track ->getnsigma_e_TPC();
            Double_t nsigma_TPC_pi  = AS_Track ->getnsigma_pi_TPC();
            Double_t nsigma_TPC_p   = AS_Track ->getnsigma_p_TPC();
            Double_t nsigma_TOF_e   = AS_Track ->getnsigma_e_TOF();
            Double_t nsigma_TOF_pi  = AS_Track ->getnsigma_pi_TOF();
            Double_t TRD_signal     = AS_Track ->getTRDSignal();
            Double_t TRDsumADC      = AS_Track ->getTRDsumADC();
            Double_t dca            = AS_Track ->getdca();  // charge * distance of closest approach to the primary vertex
            TLorentzVector TLV_part = AS_Track ->get_TLV_part();
            UShort_t NTPCcls        = AS_Track ->getNTPCcls();
            UShort_t NTRDcls        = AS_Track ->getNTRDcls();
            UShort_t NITScls        = AS_Track ->getNITScls();
            Float_t TPCchi2         = AS_Track ->getTPCchi2();
            Float_t TPCdEdx         = AS_Track ->getTPCdEdx();
            Float_t TOFsignal       = AS_Track ->getTOFsignal(); // in ps (1E-12 s)
            Float_t Track_length    = AS_Track ->getTrack_length();

            Float_t momentum        = TLV_part.P();
            Float_t eta_track       = TLV_part.Eta();
            Float_t pT_track        = TLV_part.Pt();
            Float_t theta_track     = TLV_part.Theta();

            if(pT_track < 3.0) continue;

            make_clusters(i_track);
            //get_straight_line_fit(i_track);
            get_tracklets_fit(i_track);

            vector<TVector3> vec_TV3_tracklet_vectors;
            vec_TV3_tracklet_vectors.resize(7);


            if(vec_tracklet_fit_points[6][0][0] == -999.0  && vec_tracklet_fit_points[6][1][0] == -999.0) continue; // no global fit available

            Double_t impact_angle[6] = {0.0};
            Double_t delta_x_local_global[6] = {0.0}; // local chamber coordinate system, global fit
            for(Int_t i_layer = 6; i_layer >= 0; i_layer--)
            {
                if(vec_tracklet_fit_points[i_layer][0][0] > -999.0 && vec_tracklet_fit_points[i_layer][1][0] > -999.0)
                {
                    Double_t delta_x = vec_tracklet_fit_points[i_layer][1][0] - vec_tracklet_fit_points[i_layer][0][0];
                    Double_t delta_y = vec_tracklet_fit_points[i_layer][1][1] - vec_tracklet_fit_points[i_layer][0][1];
                    Double_t delta_z = vec_tracklet_fit_points[i_layer][1][2] - vec_tracklet_fit_points[i_layer][0][2];

                    //printf("i_layer: %d, delta_z: %4.3f \n",i_layer,delta_z);
                    TVector3 TV3_delta(delta_x,delta_y,delta_z);


                    vec_TV3_tracklet_vectors[i_layer].SetXYZ(delta_x,delta_y,0.0);
                    if(vec_TV3_tracklet_vectors[i_layer].Mag() <= 0.0)
                    {
                        printf("A: {%4.3f, %4.3f}, B: {%4.3f, %4.3f} \n",vec_tracklet_fit_points[i_layer][1][0],vec_tracklet_fit_points[i_layer][1][1],vec_tracklet_fit_points[i_layer][0][0],vec_tracklet_fit_points[i_layer][0][1]);
                        continue;
                    }
                    vec_TV3_tracklet_vectors[i_layer] *= 1.0/vec_TV3_tracklet_vectors[i_layer].Mag(); // Normalize


                    // Calculate impact angles of global track to every single layer
                    if(i_layer == 6) // Global track
                    {
                        for(Int_t i_layerB = 0; i_layerB < 6; i_layerB++)
                        {
                            if(arr_layer_detector[i_layerB] < 0) continue;
                            delta_x_local_global[i_layerB] = vec_TV3_tracklet_vectors[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][0]);
                            Double_t proj_radial = TV3_delta.Dot(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2]);
                            //if(arr_layer_detector[i_layerB] == 0)
                            //{
                            //    printf("detector: %d, proj_radial: %4.3f \n",arr_layer_detector[i_layerB],proj_radial);
                                //vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2].Print();
                                //TV3_delta.Print();
                            //}
                            Double_t sign_direction_impactB = TMath::Sign(1.0,delta_x_local_global[i_layerB]);
                            impact_angle[i_layerB] = vec_TV3_tracklet_vectors[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2]);
                            if(impact_angle[i_layerB] > TMath::Pi()*0.5) impact_angle[i_layerB] -= TMath::Pi();
                            //printf("i_layerB: %d, impact_angle: %4.3f \n",i_layerB,impact_angle[i_layerB]*TMath::RadToDeg());
                            impact_angle[i_layerB] = 0.5*TMath::Pi() - sign_direction_impactB*impact_angle[i_layerB];
                            //if(impact_angle[i_layerB]*TMath::RadToDeg() > 89.0 && impact_angle[i_layerB]*TMath::RadToDeg() < 91.0)
                            //{
                            //    printf("  --> i_layerB: %d, impact_angle: %4.3f, vec: {%4.3f, %4.3f} \n",i_layerB,impact_angle[i_layerB]*TMath::RadToDeg(),vec_tracklet_fit_points[i_layerB][0][0],vec_tracklet_fit_points[i_layerB][1][0]);
                            //}
                        }
                    }


                    if(i_layer < 6) // tracklets
                    {
                        if(arr_layer_detector[i_layer] < 0) continue;
                        Int_t detector = arr_layer_detector[i_layer];
                        Double_t delta_x_local_tracklet = vec_TV3_tracklet_vectors[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);
                        Double_t sign_angle = 1.0;
                        if(delta_x_local_tracklet < delta_x_local_global[i_layer]) sign_angle = -1.0;
                        Double_t Delta_angle  = sign_angle*vec_TV3_tracklet_vectors[6].Angle(vec_TV3_tracklet_vectors[i_layer]);
                        if(Delta_angle > TMath::Pi()*0.5)  Delta_angle -= TMath::Pi();
                        if(Delta_angle < -TMath::Pi()*0.5) Delta_angle += TMath::Pi();

                        h_detector_hit ->Fill(detector);
                        //if(detector == 69) printf("impact_angle: %4.3f, Delta_angle: %4.3f \n",impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());

                        //vec_tp_Delta_vs_impact[arr_layer_detector[i_layer]] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                        //vec_tp_Delta_vs_impact[detector] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),fabs(Delta_angle*TMath::RadToDeg()));
                        vec_tp_Delta_vs_impact[detector] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                        if(impact_angle[i_layer]*TMath::RadToDeg() > 89.0 && impact_angle[i_layer]*TMath::RadToDeg() < 91.0)
                        //if(impact_angle[i_layer]*TMath::RadToDeg() > 80 && impact_angle[i_layer]*TMath::RadToDeg() < 83)
                        {
                        //    printf("detector: %d, track: %d, impact angle: %4.3f, Delta angle: %4.3f, delta_x_local_tracklet: %4.3f \n",detector,i_track,impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg(),delta_x_local_tracklet);
                            h_delta_angle_perp_impact ->Fill(Delta_angle*TMath::RadToDeg());
                        }
                        //printf("detector: %d, track: %d, impact angle: %4.3f, Delta angle: %4.3f \n",detector,i_track,impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                    }

                }
            }
        }
    }


    vector<TCanvas*> vec_can_Delta_vs_impact;
    vec_can_Delta_vs_impact.resize(6); // 6 sector blocks with 3 sectors each (18)

    TH1D* h_dummy_Delta_vs_impact = new TH1D("h_dummy_Delta_vs_impact","h_dummy_Delta_vs_impact",90,50,140);
    h_dummy_Delta_vs_impact->SetStats(0);
    h_dummy_Delta_vs_impact->SetTitle("");
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_Delta_vs_impact->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_Delta_vs_impact->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_Delta_vs_impact->GetXaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact->GetYaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact->GetXaxis()->CenterTitle();
    h_dummy_Delta_vs_impact->GetYaxis()->CenterTitle();
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitle("impact angle");
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitle("#Delta #alpha");
    h_dummy_Delta_vs_impact->GetXaxis()->SetRangeUser(70,110);
    h_dummy_Delta_vs_impact->GetYaxis()->SetRangeUser(0,24);

    Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    for(Int_t i_sec_block = 0; i_sec_block < 6; i_sec_block++)
    {
        HistName = "vec_can_Delta_vs_impact_";
        HistName += i_sec_block;
        vec_can_Delta_vs_impact[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_can_Delta_vs_impact[i_sec_block] ->Divide(5,3); // x = stack, y = sector

        for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
        {
            Int_t i_sector = i_sec_block + 6*i_sec_sub;
            for(Int_t i_stack = 0; i_stack < 5; i_stack++)
            {
                Int_t iPad = i_sec_sub*5 + i_stack + 1;
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTicks(1,1);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetGrid(0,0);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetFillColor(10);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad);
                h_dummy_Delta_vs_impact->Draw("h");

                for(Int_t i_layer = 0; i_layer < 6; i_layer++)
                {
                    Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                    printf("detector: %d \n",i_detector);
                    vec_tp_Delta_vs_impact[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
                    vec_tp_Delta_vs_impact[i_detector] ->SetLineWidth(2);
                    vec_tp_Delta_vs_impact[i_detector] ->SetLineStyle(1);
                    vec_tp_Delta_vs_impact[i_detector] ->Draw("same hl");

                    HistName = "";
                    sprintf(NoP,"%4.0f",(Double_t)i_detector);
                    HistName += NoP;
                    plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                }
            }
        }
    }


    TCanvas* can_delta_angle_perp_impact = new TCanvas("can_delta_angle_perp_impact","can_delta_angle_perp_impact",500,10,500,500);
    can_delta_angle_perp_impact ->cd();
    h_delta_angle_perp_impact  ->Draw();

    TCanvas* can_detector_hit = new TCanvas("can_detector_hit","can_detector_hit",500,500,500,500);
    can_detector_hit ->cd();
    h_detector_hit  ->Draw();

    outputfile ->cd();
    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det] ->Write();
    }

}
//----------------------------------------------------------------------------------------



#endif // __TBASE_TRD_CALIB_H__