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
ClassImp(Ali_AS_offline_Tracklet)
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
    Ali_AS_offline_Tracklet*  AS_offline_Tracklet;
    Ali_AS_TRD_digit* AS_Digit;


    TFile* outputfile;
    TFile* outputfile_trkl;
    TFile* calibration_params;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;

    TGraph* tg_HV_drift_vs_det;
    TGraph* tg_HV_anode_vs_det;
    TGraph* tg_vdrift_vs_det;

    TGraph* vdrift_fit;
    TGraph* LA_fit;
    Double_t tv3_dirs_online[6][2] = {0.0};

    Long64_t Event_active = 0;
    Double_t tracklets_min[7] = {-1.0};
    Double_t tracklets_min_circle = -1.0;
    Double_t tracklets_min_online_circle = -1.0;

    TChain* input_SE;
    TString JPsi_TREE   = "Tree_AS_Event";
    TString JPsi_BRANCH = "Tree_AS_Event_branch";
    Long64_t file_entries_total;

    Float_t digit_pos[3];

    TEvePointSet* TPM3D_single;
    vector<TEvePointSet*> vec_TPM3D_digits;
    vector<TEveLine*> vec_TPL3D_tracklets; // online tracklets
    vector<TEveLine*> vec_TPL3D_offline_tracklets; // offline tracklets
    vector<TEveLine*> vec_TPL3D_tracklets_match;
    vector<Float_t> vec_digit_single_info;
    vector< vector<Float_t> > vec_digit_info;
    vector< vector< vector<Float_t> > > vec_digit_track_info;
    vector<Float_t> vec_track_single_info;
    vector< vector<Float_t> > vec_track_info;
    vector<TPolyLine*> vec_TPL_online_tracklets;
    vector<TPolyLine*> vec_TPL_offline_tracklets;

    AliHelix aliHelix;
    TEveLine* TPL3D_helix;
    TEveLine* fit_line;
    TEveLine* vec_tracklets_line = NULL;
    vector<TPolyLine*> vec_tracklets_line_2D;
    vector<TPolyLine*> vec_TPL_online_tracklets_selected;
    vector<TPolyLine*> vec_TPL_online_tracklets_selected_corrected;

    TPolyLine* TPL_helix;

    TString HistName;
    char NoP[50];

    // TVector3 vec_datapoints;
    vector <TCanvas*> ADC_vs_time;

    vector<TEveLine*> vec_TPL3D_helix_neighbor;
    vector< vector<TEveLine*> > vec_TPL3D_TRD_center;
    vector<TVector3> vec_TV3_TRD_center_offset; // 540 chambers
    vector< vector<TVector3> >     vec_TV3_TRD_center; // 540 chambers, 3 axes

    Int_t arr_layer_detector[6];

    TGLViewer *TGL_viewer;

    TH1D* h_delta_angle_perp_impact;
    TH1D* h_delta_angle_perp_impact_circle;

    TH1D* h_detector_hit;
    vector< vector<TH1D*> > vec_h_diff_helix_line_impact_angle; // [all,-,+][540]
    vector< vector<TH1D*> > vec_h_diff_ref_loc_off_trkl; // [all,-,+][540]


    // TRD offline Tracklets
    vector< vector<TVector3> > vec_TV3_Tracklet_pos;
    vector< vector<TVector3> > vec_TV3_Tracklet_dir;


    // TRD 3D graphics
    vector<TEveBox*> vec_eve_TRD_detector_box;
    vector<TVector3> vec_TV3_local_pos;


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

    TEveLine* z_BeamLine;
    TEveLine* x_BeamLine;
    TEveLine* y_BeamLine;

    Int_t Not_installed_TRD_detectors[19] = {402,403,404,405,406,407,432,433,434,435,436,437,462,463,464,465,466,467,538};

    // list from Yvonne
    // Int_t Defect_TRD_detectors[84] = {2, 5, 17, 24, 26, 30, 31, 32, 36, 40, 41, 43, 49, 50, 59, 62, 64, 78, 88, 92, 93, 107, 110, 
    // 111, 113, 116, 119, 131, 161, 165, 182, 184, 188, 190, 191, 215, 219, 221, 223, 226, 227, 233, 236, 239, 241, 249, 255, 265, 277, 
    // 287, 302, 308, 310, 311, 318, 319, 320, 326, 328, 335, 348, 354, 368, 377, 380, 386, 389, 452, 455, 456, 470, 474, 476, 483, 484, 
    // 485, 490, 491, 493, 494, 500, 502, 504, 506};

    // list from official QA http://aliqatrd.web.cern.ch/aliqatrd/data/2016/LHC16q/pass1_CENT_wSDD/000265338/
    // Int_t Defect_TRD_detectors[91] = {2, 15, 17, 27, 31, 32, 36, 40, 43, 49, 50, 55, 59, 64, 88, 92, 113, 116, 119, 132, 180, 181, 190,
    // 191, 194, 207, 215, 219, 221, 226, 227, 228, 230, 231, 233, 236, 238, 241, 249, 255, 277, 287, 302, 308, 310, 311, 317, 318, 319, 
    // 320, 326, 328, 335, 348, 368, 377, 386, 389, 402, 403, 404, 405, 406, 407, 432, 433, 434, 435, 436, 437, 452, 455, 456, 462, 463, 
    // 464, 465, 466, 467, 470, 482, 483, 484, 485, 490, 491, 493, 500, 502, 504, 538};

    // from Jason's QC (cut max adc [2,3] > 40)
    // bad adc only
    // Int_t Defect_TRD_detectors[148] = {2, 15, 17, 27, 31, 32, 36, 40, 43, 47, 49, 50, 51, 55, 59, 64, 88, 92, 116, 119, 123, 125, 129,
    // 130, 131, 132, 136, 137, 140, 141, 142, 143, 144, 146, 148, 149, 150, 156, 157, 159, 162, 163, 164, 169, 171, 175, 180, 181, 190,
    // 191, 194, 207, 213, 219, 221, 226, 228, 230, 231, 236, 238, 241, 245, 249, 255, 277, 302, 304, 305, 308, 309, 310, 311, 317, 318,
    // 319, 320, 323, 326, 328, 335, 348, 365, 368, 371, 377, 386, 389, 391, 395, 397, 400, 401, 402, 403, 404, 405, 406, 407, 413, 417,
    // 419, 425, 427, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 442, 443, 445, 447, 448, 449, 452, 455, 456, 461, 462, 463,
    // 464, 465, 466, 467, 470, 482, 483, 484, 485, 490, 491, 493, 497, 500, 502, 504, 526, 533, 536, 537, 538};

    // bad HV only
    // Int_t Defect_TRD_detectors[86] = {2, 5, 8, 12, 17, 26, 27, 29, 30, 31, 32, 36, 40, 43, 49, 50, 51, 59, 64, 88, 92, 113, 116, 119, 132, 
    // 181, 184, 190, 191, 197, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 236, 241, 249, 255, 265, 274, 277, 287, 302, 
    // 308, 309, 310, 311, 316, 317, 318, 319, 320, 326, 328, 335, 348, 368, 377, 386, 389, 452, 455, 456, 461, 470, 483, 484, 485, 490, 491, 
    // 493, 494, 498, 500, 502, 504, 538};

    // no v_fit
    // Int_t Defect_TRD_detectors[135] = {2, 15, 17, 27, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 48, 49, 50, 51, 52, 53, 54, 
    // 55, 56, 57, 58, 59, 64, 88, 92, 113, 114, 115, 116, 117, 118, 119, 132, 180, 181, 182, 183, 184, 185, 190, 191, 194, 207, 213, 215, 
    // 219, 221, 226, 227, 228, 230, 231, 233, 234, 235, 236, 237, 238, 239, 241, 249, 255, 277, 287, 302, 306, 307, 308, 309, 310, 311, 
    // 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 335, 348, 368, 377, 384, 385, 386, 387, 388, 389, 402, 403, 404, 
    // 405, 406, 407, 432, 433, 434, 435, 436, 437, 450, 451, 452, 453, 454, 455, 456, 461, 462, 463, 464, 465, 466, 467, 470, 480, 481, 
    // 482, 483, 484, 485, 490};

    // bad adc and bad hv
    Int_t Defect_TRD_detectors[167] = {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 36, 40, 43, 47, 49, 50, 51, 55, 59, 64, 88, 92, 113, 
    116, 119, 123, 125, 129, 130, 131, 132, 136, 137, 140, 141, 142, 143, 144, 146, 148, 149, 150, 156, 157, 159, 162, 163, 164, 169, 
    171, 175, 180, 181, 190, 191, 194, 197, 207, 213, 214, 215, 219, 220, 221, 226, 227, 228, 230, 231, 232, 233, 236, 238, 241, 245, 
    249, 255, 265, 274, 277, 287, 302, 304, 305, 308, 309, 310, 311, 317, 318, 319, 320, 323, 326, 328, 335, 348, 365, 368, 371, 377, 
    386, 389, 391, 395, 397, 400, 401, 402, 403, 404, 405, 406, 407, 413, 417, 419, 425, 427, 429, 430, 431, 432, 433, 434, 435, 436, 
    437, 438, 439, 442, 443, 445, 447, 448, 449, 452, 455, 456, 461, 462, 463, 464, 465, 466, 467, 470, 476, 482, 483, 484, 485, 490, 
    491, 493, 497, 500, 502, 504, 520, 526, 533, 536, 537, 538};

    // with v_fit check
    Int_t Defect_TRD_detectors_wfit[271] =  {2, 5, 8, 12, 15, 17, 26, 27, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 47, 48,
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

    Double_t max_dca_z_to_track = 8.0; // in cm
    Double_t max_dca_r_to_track = 1.0; // in cm
    vector<Int_t> vec_merge_time_bins;

    vector< vector<TVector3> > vec_TV3_digit_pos_cluster;    // layer, merged time bin
    vector< TVector3> vec_TV3_digit_pos_cluster_t0; // layer, x, y, z
    vector< vector<TH1F*> > th1f_ADC_vs_time;
    //Int_t color_layer[7] = {kOrange+2,kGreen,kBlue,kMagenta,kCyan,kPink+5,kYellow};
    Int_t color_layer[7] = {kGray,kGreen,kBlue,kMagenta,kCyan,kYellow,kOrange};
    Int_t line_width_layer[7] = {3,3,3,3,3,3,2};
    vector<Int_t> vec_layer_in_fit;
    vector< vector< vector<Double_t> > > vec_tracklet_fit_points;
    vector< vector< vector<Double_t> > > vec_online_tracklet_points;
    vector< vector< vector<Double_t> > > vec_corrected_online_tracklet_points;
    vector<Double_t> online_dy;
    Int_t draw_corrected_online_det_numbers[6] = {0};
    vector<TProfile*> vec_tp_Delta_vs_impact;
    vector<TH2D*> vec_TH2D_Delta_vs_impact;
    vector<TProfile*> vec_tp_Delta_vs_impact_circle;
    vector<TH2D*> vec_TH2D_Delta_vs_impact_circle;

    vector<TH2D*> vec_TH2D_delta_offset_vs_impact_circle;
    Double_t vec_delta_offset[6] = {0.0};

    TH2D* h2D_delta_alpha;
    TH2D* h2D_delta_offset;

    vector<TVector2> vec_TV2_points;
    vector<TVector3> vec_dir_vec_circle;
    vector<TVector3> vec_dir_vec_circle_notempty;


    vector<TPolyLine*> vec_TPL_circle_tracklets;

    //things for Draw_2D_circle

    Double_t par_circle[3] = {0.0,0.0,0.0};
    vector<Double_t> vec_phi;
    Double_t par_circle_online[3] = {0.0,0.0,0.0};


public:
    TBase_TRD_Calib();
    ~TBase_TRD_Calib();
    void Init_tree(TString SEList);
    void Init_QA();
    Int_t Loop_event(Long64_t event);
    Long64_t get_N_Events() {return N_Events;}
    Long64_t get_N_Tracks() {return N_Tracks;}
    Long64_t get_N_Digits() {return N_Digits;}
    vector< vector<Float_t> > get_track_info() {return vec_track_info;}
    vector< vector< vector<Float_t> > > get_digit_track_info() {return vec_digit_track_info;}
    vector<TEvePointSet*> get_PM3D_digits() {return vec_TPM3D_digits;}
    vector<TEveLine*> get_PL3D_tracklets() {return vec_TPL3D_tracklets;}
    vector<TEveLine*> get_PL3D_tracklets_match() {return vec_TPL3D_tracklets_match;}
    void Draw_track(Int_t i_track);
    void Draw_2D_track(Int_t i_track);
    void Draw_line(Int_t i_track);
    void Draw_tracklets_line(Int_t i_track);
    void Draw_tracklets_line_2D(Int_t i_track);
    void Draw_neighbor_tracks(Int_t i_track);
    void Draw_online_tracklets();
    void Draw_corrected_online_tracklets();
    void calc_impact_hist();
    void Draw_offline_tracklets();
    TGLViewer* Draw_TRD();
    void set_dca_to_track(Double_t dca_r, Double_t dca_z) {max_dca_r_to_track = dca_r; max_dca_z_to_track = dca_z;}
    void set_merged_time_bins(vector<Int_t> vec_merge_time_bins_in) {vec_merge_time_bins = vec_merge_time_bins_in;}
    TEveLine* get_helix_polyline(Int_t i_track);
    TPolyLine* get_helix_polyline_2D(Int_t i_track);
    TEveLine* get_straight_line_fit(Int_t i_track);
    void get_tracklets_fit(Int_t i_track);
    Int_t get_2D_global_circle_fit();
    vector<TPolyLine*> get_online_tracklets(Int_t i_track);

    Int_t select_online_tracklets();

    vector<TPolyLine*> get_offline_tracklets(Int_t i_track);
    vector< vector<TVector3> >  make_clusters(Int_t i_track);
    vector< vector< vector<Double_t> > > get_tracklet_fit_points() {return vec_tracklet_fit_points;}
    void make_plots_ADC(Int_t i_track);
    void Calibrate();
    void Calibrate_on_trkl();
    void Track_Tracklets(); // for online tracklets
    void Draw_2D_circle();
    void Draw_2D_online_circle();

    Int_t drawing_online_tracklets_flag = 0;

    ClassDef(TBase_TRD_Calib, 1)
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TBase_TRD_Calib::TBase_TRD_Calib()
{
    //outputfile = new TFile("./TRD_Calib.root","RECREATE");

    //outputfile = new TFile("./Data/TRD_Calib_circle_3456_minos.root","RECREATE");
    //outputfile_trkl = new TFile("./Data/TRD_Calib_on_trkl_3456_minos.root","RECREATE");

    outputfile = new TFile("./Data/TRD_Calib_circle_vD_1.2_LA_0.16133.root","RECREATE");
    outputfile_trkl = new TFile("./Data/TRD_Calib_on_trkl_vD_1.2_LA_0.16133.root","RECREATE");


    Init_QA();

    vec_TV3_Tracklet_pos.resize(540);
    vec_TV3_Tracklet_dir.resize(540);
    vec_h_diff_helix_line_impact_angle.resize(3); // [all,-,+]
    vec_h_diff_ref_loc_off_trkl.resize(3); // [all,-,+]
    for(Int_t i_charge = 0; i_charge < 3; i_charge++)
    {
        vec_h_diff_helix_line_impact_angle[i_charge].resize(547);
        vec_h_diff_ref_loc_off_trkl[i_charge].resize(547);
    }

    // Standard time bins
    vec_merge_time_bins.resize(24+1);
    for(Int_t i_time = 0; i_time < (24+1); i_time++)
    {
        vec_merge_time_bins[i_time] = i_time;
    }

    vec_TPL3D_TRD_center.resize(540);
    vec_TV3_TRD_center.resize(540);
    vec_TV3_TRD_center_offset.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TPL3D_TRD_center[i_det].resize(3);
        vec_TV3_TRD_center[i_det].resize(3);
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_TPL3D_TRD_center[i_det][i_xyz] = new TEveLine();
        }
    }

    h_delta_angle_perp_impact = new TH1D("h_delta_angle_perp_impact","h_delta_angle_perp_impact",240,-30,30);
    h_delta_angle_perp_impact_circle = new TH1D("h_delta_angle_perp_impact_circle","h_delta_angle_perp_impact_circle",240,-30,30);

    h_detector_hit            = new TH1D("h_detector_hit","h_detector_hit",540,0,540);

    vec_layer_in_fit.resize(7);


    vec_digit_single_info.resize(14); // x,y,z,time,ADC,sector,stack,layer,row,column,dca,dca_x,dca_y,dca_z
    vec_track_single_info.resize(12); // dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p

    TPM3D_single = new TEvePointSet();
    vec_TPM3D_digits.resize(6); // layers
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_digits[i_layer] = new TEvePointSet();
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

    TPL3D_helix = new TEveLine();
    fit_line    = new TEveLine();


    fGeo = new AliTRDgeometry;
    vec_eve_TRD_detector_box.resize(540);
    vec_TV3_local_pos.resize(8);
    TEveManager::Create();

    for(Int_t i_charge = 0; i_charge < 3; i_charge++)
    {
        for(Int_t TRD_detector = 0; TRD_detector < 547; TRD_detector++)
        {
            HistName = "vec_h_diff_helix_line_impact_angle_c_";
            HistName += i_charge;
            HistName += "_det_";
            HistName += TRD_detector;
            vec_h_diff_helix_line_impact_angle[i_charge][TRD_detector] = new TH1D(HistName.Data(),HistName.Data(),100,-10.0,10.0);

            HistName = "vec_h_diff_ref_loc_off_trkl_";
            HistName += i_charge;
            HistName += "_det_";
            HistName += TRD_detector;
            vec_h_diff_ref_loc_off_trkl[i_charge][TRD_detector] = new TH1D(HistName.Data(),HistName.Data(),400,-2.0,2.0);
        }
    }

    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        vec_eve_TRD_detector_box[TRD_detector] = new TEveBox;

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

        HistName = "TRD_box_";
        HistName += TRD_detector;
        vec_eve_TRD_detector_box[TRD_detector] ->SetName(HistName.Data());

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


        combitrans[TRD_detector] = new TGeoCombiTrans();
        combitrans[TRD_detector] ->RotateZ(Rotation_angle + 90.0);
        combitrans[TRD_detector] ->SetTranslation(glb[0],glb[1],glb[2]);

        vec_TV3_local_pos[0].SetXYZ(-TRD_chamber_width/2.0,-TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[1].SetXYZ(TRD_chamber_width/2.0,-TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[2].SetXYZ(TRD_chamber_width/2.0,TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[3].SetXYZ(-TRD_chamber_width/2.0,TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[4].SetXYZ(-TRD_chamber_width/2.0,-TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[5].SetXYZ(TRD_chamber_width/2.0,-TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[6].SetXYZ(TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[7].SetXYZ(-TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0);

        vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(kCyan);
        vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(80);
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            Double_t arr_pos_loc[3] = {vec_TV3_local_pos[i_vertex][0],vec_TV3_local_pos[i_vertex][1],vec_TV3_local_pos[i_vertex][2]};
            Double_t arr_pos_glb[3] = {0.0,0.0,0.0};
            combitrans[TRD_detector] ->LocalToMaster(arr_pos_loc,arr_pos_glb);
            /*
            arr_pos_glb[0] += glb[0];
            arr_pos_glb[1] += glb[1];
            arr_pos_glb[2] += glb[2];
            */
            if(TRD_detector == 106)
            {
                printf("i_vertex: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);

                cout << "TRD_detector: " << TRD_detector << ", Rotation_angle: " << Rotation_angle << ", sector: " << TRD_sector
                    << ", length: " << TRD_chamber_length << ", width: " << TRD_chamber_width << ", height: " << TRD_chamber_height
                    << ", layer: " << TRD_layer << ", stack: " << TRD_stack << ", time0: " << TRD_time0
                    << ", loc = {" << loc[0] << ", " << loc[1] << ", " << loc[2] << "}"
                    << ", glb = {" << glb[0] << ", " << glb[1] << ", " << glb[2] << "}" << ", TRD_col_end: " << TRD_col_end << endl;

            }

            vec_eve_TRD_detector_box[TRD_detector]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
        }


        vec_TV3_TRD_center_offset[TRD_detector].SetXYZ(glb[0],glb[1],glb[2]);

        vec_TV3_TRD_center[TRD_detector][0].SetXYZ(glbX[0]-glb[0],glbX[1]-glb[1],glbX[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][1].SetXYZ(glbY[0]-glb[0],glbY[1]-glb[1],glbY[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][2].SetXYZ(glbZ[0]-glb[0],glbZ[1]-glb[1],glbZ[2]-glb[2]);

        printf("det: %d, vec_TV3_TRD_center[%d][2]: {%4.3f, %4.3f, %4.3f} \n",TRD_detector,TRD_detector,vec_TV3_TRD_center[TRD_detector][2].X(),vec_TV3_TRD_center[TRD_detector][2].Y(),vec_TV3_TRD_center[TRD_detector][2].Z());

    }

    for(Int_t i_det = 0; i_det < 6; i_det++)
    {
        printf("layer: %i_det, radius: %4.3f cm \n",i_det,vec_TV3_TRD_center_offset[i_det].Perp());
    }

    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
    }
    gEve->Redraw3D(kTRUE);


    vec_tracklets_line_2D.resize(7);

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Init_QA()
{
    printf("TBase_TRD_Calib::Init_QA() \n");
    TFile* inputfile_QA = TFile::Open("/home/ceres/schmah/ALICE/TRD_Run3_calib/QA_out_year_2016_V2.root");
    // TFile* inputfile_QA = TFile::Open("/home/jasonb/Documents/masters/calibration/QA_out_year_2016_V2.root");
    tg_HV_drift_vs_det = (TGraph*)inputfile_QA->Get("tg_HV_drift_vs_det");
    tg_HV_anode_vs_det = (TGraph*)inputfile_QA->Get("tg_HV_anode_vs_det");;
    tg_vdrift_vs_det   = (TGraph*)inputfile_QA->Get("tg_vdrift_vs_det");;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
vector< vector<TVector3> >  TBase_TRD_Calib::make_clusters(Int_t i_track)
{
    //printf("TBase_TRD_Calib::make_clusters(%d) \n",i_track);
    // Merge the ADC digits according to time bins set in set_merged_time_bins
    // Digit space points are merged using the ADC values as weight

    Int_t N_merged_time_bins = (Int_t)vec_merge_time_bins.size();

    AS_Track     = AS_Event ->getTrack( i_track ); // take the track

    //----------------------------------------------
    // TRD digit information
    UShort_t  fNumTRDdigits = AS_Track ->getNumTRD_digits();

    //printf("i_track: %d, fNumTRDdigits: %d \n",i_track,fNumTRDdigits);

    TVector3 TV3_digit_pos;

    //for old TV3
    vec_TV3_digit_pos_cluster.clear();
    vec_TV3_digit_pos_cluster.resize(7); // 6 layers
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        arr_layer_detector[i_layer] = -1;
        vec_TV3_digit_pos_cluster[i_layer].resize(N_merged_time_bins);
        for(Int_t i_time = 0; i_time < N_merged_time_bins; i_time++)
        {
            vec_TV3_digit_pos_cluster[i_layer][i_time].SetXYZ(-999.0,-999.0,-999.0);
        }
    }

    //for new vector<Double_t>  vector<vector<vector<Double_t>>> vec_Dt_digit_pos_cluster[i_layer][i_merged_time_bin][i_xyzADC]
    vec_Dt_digit_pos_cluster.clear();
    vec_Dt_digit_pos_cluster.resize(7); // 6 layers + global line (first time points close to anode wire)
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_Dt_digit_pos_cluster[i_layer].resize(N_merged_time_bins);
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
        vec_weight_digits_merged[i_layer].resize(N_merged_time_bins);
        vec_pos_merge[i_layer].resize(N_merged_time_bins);
        for(Int_t i_time_merge = 0; i_time_merge < (N_merged_time_bins); i_time_merge++)
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

        //printf("test 1 \n");

        ///FOR CLEAN CALIB: IGNORE ALL SUSPICIOUS CHAMBERS
        Int_t flag_defect = 0;
        for (Int_t i_defect = 0; i_defect < 91; i_defect++)
        {
            if (detector == official_qa[i_defect])
            {
                flag_defect = 1;
                break;
            }
        }
        ////////////

        if (flag_defect == 1) continue;

        for(Int_t i_time_merge = 0; i_time_merge < (N_merged_time_bins); i_time_merge++)
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
        for(Int_t i_time_merge = 0; i_time_merge < (N_merged_time_bins); i_time_merge++)
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
                        //vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC] = vec_weight_digits_merged[i_layer][i_time_merge];

                    }
                    //printf("i_xyzADC: %d, value: %4.3f \n",i_xyzADC,vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC]);
                }
                Double_t radius = vec_TV3_digit_pos_cluster[i_layer][i_time_merge].Perp();
                //printf("i_layer: %d, i_time: %d, pos: {%4.3f, %4.3f, %4.3f}, radius: %4.3f, ADC: %4.3f \n",i_layer,i_time_merge,vec_Dt_digit_pos_cluster[i_layer][i_time_merge][0],vec_Dt_digit_pos_cluster[i_layer][i_time_merge][1],vec_Dt_digit_pos_cluster[i_layer][i_time_merge][2],radius,vec_Dt_digit_pos_cluster[i_layer][i_time_merge][3]);
            }
        }
    }

    // Fill the points for the global fit
     vec_TV3_digit_pos_cluster[6].resize(6);
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_Dt_digit_pos_cluster[6][i_layer].resize(4);
        vec_TV3_digit_pos_cluster[6][i_layer].SetXYZ(0.0,0.0,0.0);

        for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
        {
            vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] = 0.0;
            //vec_TV3_digit_pos_cluster[6][i_layer][i_xyzADC] = 0.0;
            Double_t n_points_used = 0.0;
            for(Int_t i_time_merge = 0; i_time_merge < 5; i_time_merge++)
            {
                // Attention, [layer][time][xyz,ADC] -> [6][layer][xyz,ADC] on purpose to imitate the tracklet structure for the global fit
                if(vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC] <= -999.0)
                {
                    //printf("WARNING in TBase_TRD_Calib::make_clusters, odd value: %4.3f, layer: %d, xyzADC: %d, time: %d \n",vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC],i_layer,i_xyzADC,i_time_merge);
                    continue;
                }
                vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] += vec_Dt_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC]*vec_Dt_digit_pos_cluster[i_layer][i_time_merge][3];
                //vec_TV3_digit_pos_cluster[6][i_layer][i_xyzADC] += vec_TV3_digit_pos_cluster[i_layer][i_time_merge][i_xyzADC]*vec_TV3_digit_pos_cluster[i_layer][i_time_merge][3];

                n_points_used += vec_Dt_digit_pos_cluster[i_layer][i_time_merge][3];
            }
            if(n_points_used > 0.0)
            {
                vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] /= (Double_t)n_points_used;
                //vec_TV3_digit_pos_cluster[6][i_layer][i_xyzADC] /= (Double_t)n_points_used;
            }
            else
            {
                vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] = -999.0;
                //vec_TV3_digit_pos_cluster[6][i_layer][i_xyzADC] = -999.0;
            }

            for (Int_t i_layer = 0; i_layer < 6; i_layer++)
            {
                vec_TV3_digit_pos_cluster[6][i_layer].SetXYZ(vec_Dt_digit_pos_cluster[6][i_layer][0],vec_Dt_digit_pos_cluster[6][i_layer][1],vec_Dt_digit_pos_cluster[6][i_layer][2]);
                //printf("vec_Dt_digit_pos_cluster[6][i_layer][3]: %4.3f \n",vec_Dt_digit_pos_cluster[6][i_layer][3]);

            }
            //printf("vec_TV3_digit_pos_cluster[6][i_layer][0]: %4.3f \n",vec_TV3_digit_pos_cluster[6][i_layer][0]);


            //vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] = vec_Dt_digit_pos_cluster[i_layer][0][i_xyzADC]; // old version
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


        //cout << "i_digits: " << i_digits << endl;
        //cout << "i_digits_loc: " << i_digits_loc << endl;
        //cout << "layer: " << layer << endl;
        //cout << "layer check: " << layer_check << endl;

        for(Int_t i_time_bin = 0; i_time_bin < N_time_bin; i_time_bin++)
        {

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
        ADC_vs_time[i_layer] ->SetTopMargin(0.02);
        ADC_vs_time[i_layer] ->SetBottomMargin(0.18);
        ADC_vs_time[i_layer] ->SetRightMargin(0.2);
        ADC_vs_time[i_layer] ->SetLeftMargin(0.2);
        ADC_vs_time[i_layer] ->Divide(N_pads_x,N_pads_y,0.01,0.01);
        //cout << "N_digits_per_layer[i_layer]: " << N_digits_per_layer[i_layer] << endl;
        //cout << "N_pads_x: " << N_pads_x << endl;
        //cout << "N_pads_y: " << N_pads_y << endl;

        for(Int_t i_pad = 0; i_pad < N_digits_per_layer[i_layer]; i_pad++)
        {
            ADC_vs_time[i_layer]->cd(i_pad+1)->SetTicks(1,1);
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

            //cout << "th1f: " << th1f_ADC_vs_time[i_layer][i_pad]->GetBinContent(0) <<  endl;
        }
    }
}

//----------------------------------------------------------------------------------------

TEveLine* TBase_TRD_Calib::get_straight_line_fit(Int_t i_track)
{
    printf("TBase_TRD_Calib::get_straight_line_fit((%d) \n",i_track);

    //fit merged digits with a straight line

    TGraph2D*    tg_cluster_points              = new TGraph2D();
    TEveLine* digits_fit_line = new TEveLine();

    // Fill the 2D graph
    Double_t p0[4] = {10,20,1,2};

    // generate graph with the 3d points

    Int_t i_layer_notempty = 0;

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        //printf("i_layer: %d \n",i_layer);
        if(vec_Dt_digit_pos_cluster[i_layer][0][0] != -999.0 && vec_Dt_digit_pos_cluster[i_layer][0][1] != -999.0 && vec_Dt_digit_pos_cluster[i_layer][0][2] != -999.0)
        {
            tg_cluster_points->SetPoint(i_layer_notempty,vec_Dt_digit_pos_cluster[i_layer][0][0],vec_Dt_digit_pos_cluster[i_layer][0][1],vec_Dt_digit_pos_cluster[i_layer][0][2]);
            //printf("i_layer: %d, point: {%4.3f, %4.3f, %4.3f} \n",i_layer,vec_Dt_digit_pos_cluster[i_layer][0][0],vec_Dt_digit_pos_cluster[i_layer][0][1],vec_Dt_digit_pos_cluster[i_layer][0][2]);
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
    TEveLine* digits_fit_line_init = new TEveLine();
    digits_fit_line_init ->SetNextPoint(a0[0],a0[1],a0[2]);
    digits_fit_line_init ->SetNextPoint(a1[0],a1[1],a1[2]);
    digits_fit_line_init ->SetLineColor(kMagenta);
    digits_fit_line_init ->SetLineWidth(3);
    //digits_fit_line_init ->DrawClone("ogl");
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
            digits_fit_line->SetNextPoint(x,y,z);
            printf("point line: {%4.3f, %4.3f, %4.3f} \n",x,y,z);
            i_point++;
        }
    }

    digits_fit_line ->SetLineStyle(1);
    digits_fit_line ->SetLineWidth(3);
    digits_fit_line ->SetLineColor(kGreen);
    gEve ->AddElement(digits_fit_line);
    //digits_fit_line->Draw("same");

    return digits_fit_line;

}

//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_line(Int_t i_track)
{
    fit_line = get_straight_line_fit(i_track);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t TBase_TRD_Calib::get_2D_global_circle_fit()
{

    //printf("start get_2D_global_circle_fit \n");

    // Is fitting  through all first cluster points of all available layers with a 2D circle
    // First the parameters are estimated by calculating them with three points
    //vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC]

    // vec_Dt_digit_pos_cluster[6][i_layer][i_xyzADC] // use that one

    //vector<TVector2> vec_TV2_points;
    vec_TV2_points.resize(6);

    Int_t i_layer_notempty = 0;

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        if(vec_Dt_digit_pos_cluster[6][i_layer][0] != -999.0 && vec_Dt_digit_pos_cluster[6][i_layer][1] != -999.0)
        {
            vec_TV2_points[i_layer_notempty].SetX(vec_Dt_digit_pos_cluster[6][i_layer][0]);
            vec_TV2_points[i_layer_notempty].SetY(vec_Dt_digit_pos_cluster[6][i_layer][1]);
            i_layer_notempty++;
        }
    }

    vec_phi.resize(i_layer_notempty);

    for (Int_t i_layer = 0; i_layer < i_layer_notempty; i_layer++)
    {
        vec_phi[i_layer] = 0.0;
    }

    //n_layer_notempty = i_layer_notempty-1;  ?? check later

    if (i_layer_notempty-1 < 2)  //<2
    {
        //printf("!!! less than 3 points to fit circle - you got bad circle !!! \n");
        return 0;
    }

    for (Int_t layer = 0; layer < 6; layer++)
    {
        if (arr_layer_detector[layer] == 229)
        {printf("SECTOR OF 229!! i_layer_notempty: %d \n",i_layer_notempty);
            for (Int_t i_lay = 0; i_lay < 6; i_lay++)
            {
                printf("arr_layer_det[layer] = %d \n",arr_layer_detector[i_lay]);
            }
        }
    }

    Double_t p0[3] = {10,20,1};
    tracklets_min_circle = -1.0;

    //printf(" -------- MINOS Fit get_2D_global_circle_fit -------- \n");
    TVirtualFitter *min = TVirtualFitter::Fitter(0,3);

    Double_t arglist_A[1] = {-1};
    Double_t arglist_B[1] = {0};
    Double_t arglist_C[1] = {-1};
    //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
    //min->ExecuteCommand("SHOw FCNvalue",arglist_A,1);
    min->ExecuteCommand("SET PRINT",arglist_C,1);
    min->ExecuteCommand("SET NOWarnings",arglist_B,1);

    min->SetFCN(sum_distance_circ_point_2D);

    //Double_t arglist_A[1] = {-1};
    //Double_t arglist_B[1] = {0};
    //Double_t arglist_C[1] = {-1};
    //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
    //min->ExecuteCommand("SHOw FCNvalue",arglist_A,1);
    //min->ExecuteCommand("SET NOWarnings",arglist_B,1);
    //min->ExecuteCommand("SET PRINT",arglist_C,1);



    Double_t arglist[10];
    arglist[0] = 3;

    Double_t pStart[3] = {0.0,0.0,0.0}; //= {1,1,1,1};

    // old circle fit possibly causing layer 1 and 4 bias and radius underestimation
    // Double_t x1 = vec_TV2_points[0].X();
    // Double_t y1 = vec_TV2_points[0].Y();
    // Double_t x2 = vec_TV2_points[1].X();
    // Double_t y2 = vec_TV2_points[1].Y();
    // Double_t x3 = vec_TV2_points[2].X();
    // Double_t y3 = vec_TV2_points[2].Y();

    Double_t x1 = vec_TV2_points[0].X();
    Double_t y1 = vec_TV2_points[0].Y();
    Double_t x2 = vec_TV2_points[vec_TV2_points.size()-2].X();
    Double_t y2 = vec_TV2_points[vec_TV2_points.size()-2].Y();
    Double_t x3 = vec_TV2_points[vec_TV2_points.size()-1].X();
    Double_t y3 = vec_TV2_points[vec_TV2_points.size()-1].Y();


    Double_t A_help = x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2;
    Double_t B_help = (x1*x1 + y1*y1)*(y3 - y2) + (x2*x2 + y2*y2)*(y1 - y3) + (x3*x3 + y3*y3)*(y2 - y1);
    Double_t C_help = (x1*x1 + y1*y1)*(x2 - x3) + (x2*x2 + y2*y2)*(x3 - x1) + (x3*x3 + y3*y3)*(x1 - x2);
    Double_t D_help = (x1*x1 + y1*y1)*(x3*y2 - x2*y3) + (x2*x2 + y2*y2)*(x1*y3 - x3*y1) + (x3*x3 + y3*y3)*(x2*y1 - x1*y2);

    //printf("A_help: {%4.3f, %4.3f}, B_help: {%4.3f, %4.3f}, C_help: {%4.3f, %4.3f}  \n",A_help,B_help,C_help);

    Double_t a_param = 0.0;
    Double_t b_param = 0.0;
    Double_t R_param = 0.0;
  

    if(A_help != 0.0)
    {
        a_param = -B_help/(2*A_help);
        b_param = -C_help/(2*A_help);
        R_param = TMath::Sqrt(TMath::Power(x1 - a_param,2.0) + TMath::Power(y1 - b_param,2.0));
    }

    //printf("pointA: {%4.3f, %4.3f}, pointB: {%4.3f, %4.3f}, pointC: {%4.3f, %4.3f}, circle(a,b,R): {%4.3f, %4.3f, %4.3f} \n",x1,y1,x2,y2,x3,y3,a_param,b_param,R_param);

    pStart[0] = a_param;
    pStart[1] = b_param;
    pStart[2] = R_param;


    //cout << "pStart[0]" << pStart[0] << endl;
    //cout << "pStart[1]" << pStart[1] << endl;
    //cout << "pStart[2]" << pStart[2] << endl;

    min->SetParameter(0,"a_param",pStart[0],0.01,0,0);
    min->SetParameter(1,"b_param",pStart[1],0.01,0,0);
    min->SetParameter(2,"R_param",pStart[2],0.01,0,0);

    arglist[0] = 1000; // number of function calls
    arglist[1] = 0.001; // tolerance

    //min->ExecuteCommand("MIGRAD",arglist,2);


    //min->ExecuteCommand("MINOS",arglist,2);
    min->ExecuteCommand("MINIMIZE",arglist,2);

    //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
    //min->ExecuteCommand("SET NOWarnings",arglist_B,1);
    //min->ExecuteCommand("SET PRINT",arglist_C,1);
    //printf(" -------- END MINOS Fit get_2D_global_circle_fit -------- \n");

    //if (minos) min->ExecuteCommand("MINOS",arglist,0);
    Int_t nvpar,nparx;
    Double_t amin,edm, errdef;
    min->GetStats(amin,edm,errdef,nvpar,nparx);
    //min->PrintResults(1,amin);

    Double_t parFit_circ[3];

    for(Int_t i_par = 0; i_par <3; ++i_par)
    {
        parFit_circ[i_par] = min->GetParameter(i_par);
        //parFit[i] = pStart[i];
    }

    for(Int_t i_par = 0; i_par <3; ++i_par)
    {
        par_circle[i_par] = parFit_circ[i_par];
    }

    //printf("amin: %4.3f \n",amin);
    tracklets_min_circle = amin;


    // Calculate direction vectors from 2D circle fit
    TVector2 dir_vector;
    TVector2 circle_center_vector;
    circle_center_vector.SetX(parFit_circ[0]);
    circle_center_vector.SetY(parFit_circ[1]);

        for (Int_t layer = 0; layer < 6; layer++)
    {
        if (arr_layer_detector[layer] == 229) printf("parFit_circ[0]: %4.3f, parFit_circ[1]: %4.3f \n",parFit_circ[0],parFit_circ[1]);
    }


    //printf("parFit_circ[0]: %4.3f, parFit_circ[1]: %4.3f \n",parFit_circ[0],parFit_circ[1]);

    //printf("i_layer_notempty: %d,  \n",i_layer_notempty);

    for(Int_t i_point = 0; i_point < i_layer_notempty; i_point++)
    {
        dir_vector = vec_TV2_points[i_point];

        //printf("dir_vector.X(): %4.3f, dir_vector.Y(): %4.3f \n",dir_vector.X(),dir_vector.Y());

        dir_vector -= circle_center_vector;
        Double_t phi_val = dir_vector.Phi();

        //printf("phi_val: %4.3f \n",phi_val);

        //vec_phi.push_back(phi_val);
        vec_phi[i_point] = phi_val;

        //printf("i_point: %d, vec_phi[i_point]: %4.3f, \n",i_point,vec_phi[i_point]);
    }

    Double_t delta_phi = -(vec_phi[1] - vec_phi[0])/1000.0;

    //printf("delta_phi: %4.3f, vec_phi[1]: %4.3f, vec_phi[0]: %4.3f  \n",delta_phi,vec_phi[1],vec_phi[0]);

    vec_dir_vec_circle.resize(6);

    for (Int_t i_layer = 0; i_layer <6; ++i_layer)
    {
        vec_dir_vec_circle[i_layer].SetXYZ(-999.0,-999.0,-999.0);   //because i maybe need it to have all 6 dimensions in Calibrate
    }

    //this one will have only n_notempty dimensions
    vec_dir_vec_circle_notempty.resize(i_layer_notempty);

    for(Int_t i_point = 0; i_point < i_layer_notempty; i_point++)
    //for(Int_t i_point = 0; i_point < 6; i_point++)

    {
        //printf("from here \n");
        Double_t phi_val = vec_phi[i_point];

        Double_t x_val   =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val);
        Double_t y_val   =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val);

        Double_t x_valB  =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val + delta_phi);
        Double_t y_valB  =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val + delta_phi);
        //printf("x_val: %4.3f, x_valB: %4.3f, y_val: %4.3f, y_valB: %4.3f \n",x_val,x_valB,y_val,y_valB);


        //TVector3 dir_vec_circle;               //maybe i can use this one

        vec_dir_vec_circle_notempty[i_point].SetX(x_valB - x_val);
        vec_dir_vec_circle_notempty[i_point].SetY(y_valB - y_val);
        vec_dir_vec_circle_notempty[i_point].SetZ(0.0);

        //printf("vec_dir_vec_circle_notempty[i_point].X: %4.3f, vec_dir_vec_circle_notempty[i_point].Y: %4.3f, vec_dir_vec_circle_notempty[i_point].Z: %4.3f \n",vec_dir_vec_circle_notempty[i_point].X(),vec_dir_vec_circle_notempty[i_point].Y(),vec_dir_vec_circle_notempty[i_point].Z());
        //printf("to here; what happened?? \n");


        //printf("x_val: %4.3f, x_valB: %4.3f, y_val: %4.3f, y_valB: %4.3f \n",x_val,x_valB,y_val,y_valB);  probably correct

        if(vec_dir_vec_circle_notempty[i_point].Mag() > 0.0)
        {
            vec_dir_vec_circle_notempty[i_point] *= 6.0/vec_dir_vec_circle_notempty[i_point].Mag();
        }
    }

    Int_t i_layer_notempty_check = 0;

    //hope it's gonna work (if it's needed at all)
    // cout << endl << "CIRCLE_FIT():" << endl;

    for (Int_t i_layer = 0; i_layer < 6; ++i_layer)
    {
        if(vec_Dt_digit_pos_cluster[6][i_layer][0] != -999.0 && vec_Dt_digit_pos_cluster[6][i_layer][1] != -999.0)
        {
            vec_dir_vec_circle[i_layer].SetXYZ(vec_dir_vec_circle_notempty[i_layer_notempty_check].X(),vec_dir_vec_circle_notempty[i_layer_notempty_check].Y(),vec_dir_vec_circle_notempty[i_layer_notempty_check].Z());
            Double_t length_vec = vec_dir_vec_circle[i_layer].Mag();
            if(length_vec > 0.0)
            {
                vec_dir_vec_circle[i_layer] *= -1.0/length_vec;
            }
            i_layer_notempty_check++;

            //printf("vec_dir_vec_circle[i_layer].X: %4.3f, vec_dir_vec_circle[i_layer].Y: %4.3f, vec_dir_vec_circle[i_layer].Z: %4.3f \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());



            // calculate distance along local y coord for residual information/plots
            TVector3 given_offset = {vec_Dt_digit_pos_cluster[6][i_layer][0], vec_Dt_digit_pos_cluster[6][i_layer][1], vec_Dt_digit_pos_cluster[6][i_layer][2]};

            int sector_number = arr_layer_detector[i_layer]/30;
            Double_t transform_angle = -(90 -(sector_number*20 + 10))*TMath::DegToRad();

            TVector3 local_y_unit_vec = (TVector3){TMath::Cos(transform_angle), TMath::Sin(transform_angle), 0}.Unit();

            // b_offset = tracklet_offset[b_det][b_trkl];
            // b_dir = tracklet_dir[b_det][b_trkl];

            TVector3 pos;

            Double_t radius_offset = sqrt(pow(given_offset[0] - parFit_circ[0], 2) + pow(given_offset[1] - parFit_circ[1], 2));
            Double_t delta_radius = radius_offset -  parFit_circ[2];

            TVector3 radial_unit_vec = (TVector3){given_offset[0] - parFit_circ[0], given_offset[1] - parFit_circ[1], 0}.Unit();

            if (local_y_unit_vec.Mag() != 0)
            {
                pos = given_offset - (delta_radius) * local_y_unit_vec * (1 / (local_y_unit_vec * radial_unit_vec));
            }

            Double_t offset_diff = (given_offset - pos).Mag()*10;

            Double_t offset_local_y = given_offset[0]*cos(transform_angle) - given_offset[1]*sin(transform_angle);
            Double_t pos_local_y = pos[0]*cos(transform_angle) - pos[1]*sin(transform_angle);

            // if (radius_offset < parFit_circ[2])
            if (pos_local_y > offset_local_y)
            {
                offset_diff = -offset_diff;
            }

            // cout << endl <<  "SECTOR: " << sector_number << " | " << "local y: " << local_y_unit_vec[0] << ", " << local_y_unit_vec[1] << endl;
            // cout << "offset point: " << given_offset[0] << ", " << given_offset[1] << ", " << given_offset[2] << endl;
            // cout << "intersect point: " << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;
            // cout << "offset radius: " << radius_offset << " | " << "circle radius: " << parFit_circ[2] << endl;
            // cout << "radial dist: " << delta_radius << endl;
            // cout << "layer: " << i_layer << " | " << "offset diff: " << offset_diff << endl;

            // cout << parFit_circ[0] << ", " << parFit_circ[1] << ", " << parFit_circ[2] << endl;

            vec_delta_offset[i_layer] = offset_diff;




        }
        //printf("vec_dir_vec_circle[i_layer].X: %4.3f, vec_dir_vec_circle[i_layer].Y: %4.3f, vec_dir_vec_circle[i_layer].Z: %4.3f \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());
        for (Int_t layer = 0; layer < 6; layer++)
        {
        if (arr_layer_detector[layer] == 229) printf("vec_dir_vec_circle[i_layer].X: %4.3f, vec_dir_vec_circle[i_layer].Y: %4.3f, vec_dir_vec_circle[i_layer].Z: %4.3f \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());
        }

    }




    delete min;
    return 1;
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

    Double_t layer_dist_min_max[7][2] =
    {
        {280.0,290.0},
        {290.0,300.0},
        {300.0,310.0},
        {310.0,320.0},
        {320.0,330.0},
        {330.0,340.0},
        {290.0,380.0},
    };


    Double_t delta_layer = 12.5;

    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        tracklets_min[i_layer] = -1.0;

        if(i_layer < 6)
        {
            layer_dist_min_max[i_layer][0] = 295.0 + i_layer*delta_layer;
            layer_dist_min_max[i_layer][1] = 307.0 + i_layer*delta_layer;
        }

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

        //printf("     ------------------ MIGRAD get_tracklets_fit ------------------ \n");
        Double_t arglist_A[1] = {-1};
        Double_t arglist_B[1] = {0};
        Double_t arglist_C[1] = {-1};
        //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
        min->ExecuteCommand("SET PRINT",arglist_C,1);
        min->ExecuteCommand("SET NOWarnings",arglist_B,1);

        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            a0[i_xyz] = vec_Dt_digit_pos_cluster[i_layer][i_time_merge_AB[0]][i_xyz];
            a1[i_xyz] = vec_Dt_digit_pos_cluster[i_layer][i_time_merge_AB[1]][i_xyz];
        }

        //printf("point start: {%4.3f, %4.3f, %4.3f} \n",a0[0],a0[1],a0[2]);
        //printf("point end: {%4.3f, %4.3f, %4.3f} \n",a1[0],a1[1],a1[2]);

#if 0
        TEveLine* digits_fit_line_init = new TEveLine();
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

       
        min->ExecuteCommand("MIGRAD",arglist,2);
        //printf("     ------------------ END get_tracklets_fit ------------------ \n");

        //if (minos) min->ExecuteCommand("MINOS",arglist,0);

        min->GetStats(amin,edm,errdef,nvpar,nparx);
        //min->PrintResults(1,amin);

        // get fit parameters
        for(int i = 0; i < 4; ++i)
        {
            parFit[i] = min->GetParameter(i);
            //parFit[i] = pStart[i];
        }

        //printf("i_layer: %d, amin: %4.3f, par: {%4.3f, %4.3f, %4.3f, %4.3f} \n",i_layer,amin,parFit[0],parFit[1],parFit[2],parFit[3]);

        tracklets_min[i_layer] = amin;

        // draw the fitted line
        n = 5000;   // 1000
        t0 = -500.0; // -500
        dt = 0.2;

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

        for(int i = 0; i < n; ++i)
        {
            Double_t t = t0 + dt*i;
            Double_t x,y,z;
            if(flag_XZ == 0) line_X(t,parFit,x,y,z);
            else line(t,parFit,x,y,z);
            TV3_line_point.SetXYZ(x,y,z);
            //if(i_layer == 5) printf("point: {%4.3f, %4.3f, %4.3f} \n",x,y,z);

            Double_t distance = 1000.0;
            for(Int_t i_time_bin = 0; i_time_bin < (Int_t)vec_Dt_digit_pos_cluster[i_layer].size(); i_time_bin++)
            {
                TVector3 vec_TV3_digit_pos_cluster;
                vec_TV3_digit_pos_cluster.SetXYZ(vec_Dt_digit_pos_cluster[i_layer][i_time_bin][0],vec_Dt_digit_pos_cluster[i_layer][i_time_bin][1],vec_Dt_digit_pos_cluster[i_layer][i_time_bin][2]);
                TVector3 TV3_diff = vec_TV3_digit_pos_cluster - TV3_line_point;
                Double_t distance_layer = TV3_diff.Mag();
                if(distance_layer < distance) distance = distance_layer;
            }
            //if(i_layer == 5) printf("point: {%4.3f, %4.3f, %4.3f}, perp: %4.3f, layer min/max: {%4.3f, %4.3f}, distance: %4.3f \n",x,y,z,TV3_line_point.Perp(),layer_dist_min_max[i_layer][0],layer_dist_min_max[i_layer][1],distance);
            //if(TV3_line_point.Perp() > 270.0 && TV3_line_point.Perp() < 380.0 && distance < 100.0)  // 100 just in case
            if(TV3_line_point.Perp() > layer_dist_min_max[i_layer][0] && TV3_line_point.Perp() < layer_dist_min_max[i_layer][1] && distance < 100.0)  // 100 just in case
            {
                //if(i_layer == 5) printf("Accepted \n");
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
vector<TPolyLine*> TBase_TRD_Calib::get_online_tracklets(Int_t i_track)
{
    // Loop over all tracklets
    Double_t scale_factor_length = 3.0;
    Int_t    N_TRD_tracklets_offline = AS_Event ->getNumTracklets();
    vec_TPL_online_tracklets.clear();
    for(UShort_t i_tracklet = 0; i_tracklet < N_TRD_tracklets_offline; ++i_tracklet) // loop over all tracklets of the actual event
    {
        AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
        TVector3 TV3_offset     = AS_Tracklet ->get_TV3_offset(); // online tracklets
        TVector3 TV3_dir        = AS_Tracklet ->get_TV3_dir();    // online tracklets
        Short_t  i_det_tracklet = AS_Tracklet ->get_detector();

        Int_t flag_match = 0;
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            if(i_det_tracklet == arr_layer_detector[i_layer])
            {
                flag_match = 1;
                break;
            }
        }

        if(!flag_match) continue;

        Double_t impact_angle = TV3_dir.Angle(vec_TV3_TRD_center[i_det_tracklet][2]);
        if(impact_angle > TMath::Pi()*0.5) impact_angle -= TMath::Pi();

        Float_t points[6] =
        {
            (Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),
            (Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2])
        };
        //printf("i_tracklet: %d, out of %d, impact_angle: %4.3f, offset: {%4.3f, %4.3f, %4.3f}, end: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,N_TRD_tracklets_offline,impact_angle*TMath::RadToDeg(),TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + TV3_dir[0],TV3_offset[1] + TV3_dir[1],TV3_offset[2] + TV3_dir[2]);

        vec_TPL_online_tracklets.push_back(new TPolyLine());
        vec_TPL_online_tracklets[(Int_t)vec_TPL_online_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]));
        vec_TPL_online_tracklets[(Int_t)vec_TPL_online_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]));
    }

    return vec_TPL_online_tracklets;
}
//----------------------------------------------------------------------------------------

Int_t TBase_TRD_Calib::select_online_tracklets()
{

    //take parameters of circle fit and chose tracklets which are not more than 3 cm away from circle fit

    //prepare vec_online_tracklet_points[i_layer][i_start_stop][i_XY] where we will keep tracklets points

    Int_t N_trkl_per_layer = 0;

    Double_t scale_factor_length = 3.0;
    Int_t    N_TRD_tracklets_online = AS_Event ->getNumTracklets();

    vec_online_tracklet_points.resize(6); // layer
    vec_corrected_online_tracklet_points.resize(6);
    online_dy.resize(6);
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_online_tracklet_points[i_layer].resize(2);
        vec_corrected_online_tracklet_points[i_layer].resize(2);
        for(Int_t i_startstop = 0; i_startstop < 2; i_startstop++)
        {
            vec_online_tracklet_points[i_layer][i_startstop].resize(2); // x,y,z
            vec_corrected_online_tracklet_points[i_layer][i_startstop].resize(2);
            for(Int_t i_xy = 0; i_xy < 2; i_xy++)
            {
                vec_online_tracklet_points[i_layer][i_startstop][i_xy] = -999.0;
            }
        }
    }

#if 0
    for (Int_t i_layer = 0; i_layer < 6; ++i_layer)
    {
        printf("vec_online_tracklet_points[%d][0][0] before everything: %4.3f \n",i_layer,vec_online_tracklet_points[i_layer][0][0]);
    }
#endif


    //loop over each online tracklet and select those who are not more that 3cm away from circle fit
    for(UShort_t i_tracklet = 0; i_tracklet < N_TRD_tracklets_online; ++i_tracklet) // loop over all tracklets of the actual event
    {
        //Int_t debug = 1;
        //if(debug==1) break;

        AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
        TVector3 TV3_offset     = AS_Tracklet ->get_TV3_offset(); // online tracklets
        TVector3 TV3_dir        = AS_Tracklet ->get_TV3_dir();    // online tracklets
        Short_t  i_det_tracklet = AS_Tracklet ->get_detector();

        Int_t    TRD_layer      = fGeo        ->GetLayer(i_det_tracklet);

        Double_t dy = AS_Tracklet -> get_online_dy();

        //printf("TRD_layer: %d \n",TRD_layer);

        Int_t flag_match = 0;
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            //printf("weird loop. i_layer: %d, i_det_tracklet: %d, arr_layer_detector[i_layer]: %d \n",i_layer,i_det_tracklet,arr_layer_detector[i_layer]);

            if(i_det_tracklet == arr_layer_detector[i_layer])
            {
                flag_match = 1;
                break;
            }
        }

        if(!flag_match) continue;

        Double_t impact_angle = TV3_dir.Angle(vec_TV3_TRD_center[i_det_tracklet][2]);
        if(impact_angle > TMath::Pi()*0.5) impact_angle -= TMath::Pi();
        
        Double_t dist = distance_circ_point_2D(TV3_offset[0],TV3_offset[1],par_circle);
        
        if (dist<3)

        // X,Y of beginning; X,Y of end

        if(dist>=1) continue;

        if (dist < 1.0)
        {
            //printf("we have a match! dist: %4.3f \n",dist);

            if (vec_online_tracklet_points[TRD_layer][0][0] != -999.0)
            {
                //printf("tracklet already taken.. \n");
            }
            if (vec_online_tracklet_points[TRD_layer][0][0] == -999.0)
            {
                //printf("tracklet is not taken \n");

                vec_online_tracklet_points[TRD_layer][0][0] = TV3_offset[0];
                vec_online_tracklet_points[TRD_layer][0][1] = TV3_offset[1];
                vec_online_tracklet_points[TRD_layer][1][0] = TV3_offset[0] + scale_factor_length*TV3_dir[0];
                vec_online_tracklet_points[TRD_layer][1][1] = TV3_offset[1] + scale_factor_length*TV3_dir[1];

                draw_corrected_online_det_numbers[TRD_layer] = i_det_tracklet;
                tv3_dirs_online[TRD_layer][0] = TV3_dir[0];
                tv3_dirs_online[TRD_layer][1] = TV3_dir[1];


                vec_corrected_online_tracklet_points[TRD_layer][0][0] = TV3_offset[0];
                vec_corrected_online_tracklet_points[TRD_layer][0][1] = TV3_offset[1];
                vec_corrected_online_tracklet_points[TRD_layer][1][0] = TV3_offset[0];
                vec_corrected_online_tracklet_points[TRD_layer][1][1] = TV3_offset[1];

                online_dy[TRD_layer] = dy;
                // vec_corrected_online_tracklet_points[TRD_layer][1][0] = TV3_offset[0] + scale_factor_length*TV3_dir[0];
                // vec_corrected_online_tracklet_points[TRD_layer][1][1] = TV3_offset[1] + scale_factor_length*TV3_dir[1];

                //Float_t points[6] =
                //{
                //    (Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),
                //    (Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2])
                //};
                //printf("i_tracklet: %d, out of %d, impact_angle: %4.3f, offset: {%4.3f, %4.3f, %4.3f}, end: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,N_TRD_tracklets_online,impact_angle*TMath::RadToDeg(),TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + TV3_dir[0],TV3_offset[1] + TV3_dir[1],TV3_offset[2] + TV3_dir[2]);
            }
            //printf("vec_online_tracklet_points[%d][0][0] after: %4.3f \n",TRD_layer,vec_online_tracklet_points[TRD_layer][0][0]);
        }
    }
    // cout << "this is the track: " << TV3_offset.X() << ' ' << TV3_offset.Y() << endl;
    // we will be fitting to vec_online_tracklet_points[layer_i][0][xy_i]

    vec_TV2_points.resize(6);

    Int_t i_layer_notempty = 0;

    for (Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        Double_t online_tracklet_X = vec_online_tracklet_points[i_layer][0][0];
        Double_t online_tracklet_Y = vec_online_tracklet_points[i_layer][0][1];
        if (online_tracklet_X != -999.0 && online_tracklet_Y != -999.0)
        {
            vec_TV2_points[i_layer].SetX(online_tracklet_X);
            vec_TV2_points[i_layer].SetY(online_tracklet_Y);
            i_layer_notempty++;
        }

        // cout << "tracklet points layer: " << i_layer << endl;
        // cout << "X: " << vec_online_tracklet_points[i_layer][0][0] << endl;
        // cout << "Y: " << vec_online_tracklet_points[i_layer][0][1] << endl;
    } 


    //n_layer_notempty = i_layer_notempty-1;  ?? check later

    if (i_layer_notempty < 3)
    {
        // printf("!!! less than 3 points to fit online circle - you got bad circle !!! \n");
        return 0;
    }
    
    Double_t p0[3] = {10,20,1};
    tracklets_min_online_circle = -1.0;

    TVirtualFitter *min = TVirtualFitter::Fitter(0,3);

    min->SetFCN(sum_distance_circ_point_2D);

    Double_t arglist_A[1] = {-1};
    Double_t arglist_B[1] = {0};
    Double_t arglist_C[1] = {-1};

    //print all information of minimization procedure
    Double_t arglist_D[1] = {3};


    printf("     ------------------ MIGRAD  select_online_tracklets ------------------ \n");
    //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
    //min->ExecuteCommand("SHOw FCNvalue",arglist_A,1);
    //min->ExecuteCommand("SET NOWarnings",arglist_B,1);
    //min->ExecuteCommand("SET PRINT",arglist_C,1);
    printf("     ------------------ END  select_online_tracklets ------------------ \n");


    Double_t arglist[10];
    arglist[0] = 3;

    Double_t pStart[3]; //= {1,1,1,1};

    Double_t x1 = vec_TV2_points[0].X();
    Double_t y1 = vec_TV2_points[0].Y();
    Double_t x2 = vec_TV2_points[1].X();
    Double_t y2 = vec_TV2_points[1].Y();
    Double_t x3 = vec_TV2_points[2].X();
    Double_t y3 = vec_TV2_points[2].Y();


    Double_t A_help = x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2;
    Double_t B_help = (x1*x1 + y1*y1)*(y3 - y2) + (x2*x2 + y2*y2)*(y1 - y3) + (x3*x3 + y3*y3)*(y2 - y1);
    Double_t C_help = (x1*x1 + y1*y1)*(x2 - x3) + (x2*x2 + y2*y2)*(x3 - x1) + (x3*x3 + y3*y3)*(x1 - x2);
    Double_t D_help = (x1*x1 + y1*y1)*(x3*y2 - x2*y3) + (x2*x2 + y2*y2)*(x1*y3 - x3*y1) + (x3*x3 + y3*y3)*(x2*y1 - x1*y2);

    //printf("A_help: {%4.3f, %4.3f}, B_help: {%4.3f, %4.3f}, C_help: {%4.3f, %4.3f}  \n",A_help,B_help,C_help);

    Double_t a_param = 0.0;
    Double_t b_param = 0.0;
    Double_t R_param = 0.0;


    if(A_help != 0.0)
    {
        a_param = -B_help/(2*A_help);
        b_param = -C_help/(2*A_help);
        R_param = TMath::Sqrt(TMath::Power(x1 - a_param,2.0) + TMath::Power(y1 - b_param,2.0));
    }

    //printf("pointA: {%4.3f, %4.3f}, pointB: {%4.3f, %4.3f}, pointC: {%4.3f, %4.3f}, circle(a,b,R): {%4.3f, %4.3f, %4.3f} \n",x1,y1,x2,y2,x3,y3,a_param,b_param,R_param);

    pStart[0] = a_param;
    pStart[1] = b_param;
    pStart[2] = R_param;


    //cout << "pStart[0]" << pStart[0] << endl;
    //cout << "pStart[1]" << pStart[1] << endl;
    //cout << "pStart[2]" << pStart[2] << endl;

    min->SetParameter(0,"a_param",pStart[0],0.01,0,0);
    min->SetParameter(1,"b_param",pStart[1],0.01,0,0);
    min->SetParameter(2,"R_param",pStart[2],0.01,0,0);

    arglist[0] = 1000; // number of function calls
    arglist[1] = 0.001; // tolerance

    printf("     ------------------ MIGRAD  select_online_tracklets B ------------------ \n");
    //min->ExecuteCommand("MIGRAD",arglist,2);
    //min->ExecuteCommand("MINOS",arglist,2);
    min->ExecuteCommand("MINIMIZE",arglist,2);
    printf("     ------------------ END  select_online_tracklets B ------------------ \n");

    //if (minos) min->ExecuteCommand("MINOS",arglist,0);
    Int_t nvpar,nparx;
    Double_t amin,edm, errdef;
    min->GetStats(amin,edm,errdef,nvpar,nparx);
    tracklets_min_online_circle = amin;
    //min->PrintResults(1,amin);

    Double_t parFit_circ[3];

    for(Int_t i_par = 0; i_par <3; ++i_par)
    {
        parFit_circ[i_par] = min->GetParameter(i_par);
        //parFit[i] = pStart[i];
    }

    for(Int_t i_par = 0; i_par <3; ++i_par)
    {
        par_circle_online[i_par] = parFit_circ[i_par];
    }


    //tangents
    TVector2 dir_vector;
    TVector2 circle_center_vector;
    circle_center_vector.SetX(parFit_circ[0]);
    circle_center_vector.SetY(parFit_circ[1]);

    //printf("parFit_circ[0]: %4.3f, parFit_circ[1]: %4.3f \n",parFit_circ[0],parFit_circ[1]);

    //printf("i_layer_notempty: %d,  \n",i_layer_notempty);

    for(Int_t i_point = 0; i_point < i_layer_notempty; i_point++)
    {
        dir_vector = vec_TV2_points[i_point];

        //printf("dir_vector.X(): %4.3f, dir_vector.Y(): %4.3f \n",dir_vector.X(),dir_vector.Y());

        dir_vector -= circle_center_vector;
        Double_t phi_val = dir_vector.Phi();

        //printf("phi_val: %4.3f \n",phi_val);

        //vec_phi.push_back(phi_val);
        vec_phi[i_point] = phi_val;

        //printf("i_point: %d, vec_phi[i_point]: %4.3f, \n",i_point,vec_phi[i_point]);
    }

    Double_t delta_phi = -(vec_phi[1] - vec_phi[0])/1000.0;

    //printf("delta_phi: %4.3f, vec_phi[1]: %4.3f, vec_phi[0]: %4.3f  \n",delta_phi,vec_phi[1],vec_phi[0]);

    vec_dir_vec_circle.resize(6);

    for (Int_t i_layer = 0; i_layer <6; ++i_layer)
    {
        vec_dir_vec_circle[i_layer].SetXYZ(-999.0,-999.0,-999.0);   //because i maybe need it to have all 6 dimensions in Calibrate
    }

    //this one will have only n_notempty dimensions
    vec_dir_vec_circle_notempty.resize(i_layer_notempty);

    for(Int_t i_point = 0; i_point < i_layer_notempty; i_point++)
    //for(Int_t i_point = 0; i_point < 6; i_point++)

    {
        //printf("from here \n");
        Double_t phi_val = vec_phi[i_point];

        Double_t x_val   =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val);
        Double_t y_val   =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val);

        Double_t x_valB  =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val + delta_phi);
        Double_t y_valB  =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val + delta_phi);
        //printf("x_val: %4.3f, x_valB: %4.3f, y_val: %4.3f, y_valB: %4.3f \n",x_val,x_valB,y_val,y_valB);


        //TVector3 dir_vec_circle;               //maybe i can use this one

        vec_dir_vec_circle_notempty[i_point].SetX(x_valB - x_val);
        vec_dir_vec_circle_notempty[i_point].SetY(y_valB - y_val);
        vec_dir_vec_circle_notempty[i_point].SetZ(0.0);

        //printf("vec_dir_vec_circle_notempty[i_point].X: %4.3f, vec_dir_vec_circle_notempty[i_point].Y: %4.3f, vec_dir_vec_circle_notempty[i_point].Z: %4.3f \n",vec_dir_vec_circle_notempty[i_point].X(),vec_dir_vec_circle_notempty[i_point].Y(),vec_dir_vec_circle_notempty[i_point].Z());
        //printf("to here; what happened?? \n");


        // printf("x_val: %4.3f, x_valB: %4.3f, y_val: %4.3f, y_valB: %4.3f \n",x_val,x_valB,y_val,y_valB);  probably correct

        if(vec_dir_vec_circle_notempty[i_point].Mag() > 0.0)
        {
            vec_dir_vec_circle_notempty[i_point] *= 6.0/vec_dir_vec_circle_notempty[i_point].Mag();
        }
    }

    Int_t i_layer_notempty_check = 0;

    //hope it's gonna work (if it's needed at all)
    for (Int_t i_layer = 0; i_layer < 6; ++i_layer)
    {
        if(vec_Dt_digit_pos_cluster[6][i_layer][0] != -999.0 && vec_Dt_digit_pos_cluster[6][i_layer][1] != -999.0)
        {
            vec_dir_vec_circle[i_layer].SetXYZ(vec_dir_vec_circle_notempty[i_layer_notempty_check].X(),vec_dir_vec_circle_notempty[i_layer_notempty_check].Y(),vec_dir_vec_circle_notempty[i_layer_notempty_check].Z());
            Double_t length_vec = vec_dir_vec_circle[i_layer].Mag();
            if(length_vec > 0.0)
            {
                vec_dir_vec_circle[i_layer] *= -1.0/length_vec;
            }
            i_layer_notempty_check++;

            // printf("vec_dir_vec_circle[i_layer].X: %4.3f, vec_dir_vec_circle[i_layer].Y: %4.3f, vec_dir_vec_circle[i_layer].Z: %4.3f \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());

        }
        // printf("vec_dir_vec_circle[i_layer].X: %4.3f, vec_dir_vec_circle[i_layer].Y: %4.3f, vec_dir_vec_circle[i_layer].Z: %4.3f \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());

    // Double_t impact_angle = TV3_dir.Angle(vec_TV3_TRD_center[i_det_tracklet][2]);
    // impact_angle_circle[i_layer] = vec_dir_vec_circle[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layer]][2]);  //i need something like that for circles

    // //printf("test 3.4 \n");

    // if(impact_angle_circle[i_layer] > TMath::Pi()*0.5) impact_angle_circle[i_layer] -= TMath::Pi();
    // //printf("i_layerB: %d, 3D_Angle_on_TRD(TPC track): %4.3f, 2D_impact_angle: %4.3f \n",i_layerB,Angle_on_TRD*TMath::RadToDeg(),impact_angle[i_layerB]*TMath::RadToDeg());
    // impact_angle_circle[i_layer] = 0.5*TMath::Pi() - sign_direction_impact_circle*impact_angle_circle[i_layer];
    // //printf("impact_angle_circle: %4.3f \n",impact_angle_circle[i_layer]);

    }

    delete min;
    return 1;
}

//----------------------------------------------------------------------------------------
vector<TPolyLine*> TBase_TRD_Calib::get_offline_tracklets(Int_t i_track)
{
    // Loop over all tracklets
    Double_t scale_factor_length = 3.0;

    vec_TPL_offline_tracklets.clear();

    UShort_t N_tracks_off = AS_Event ->getNumTracks(); // number of tracks in this event
    for(UShort_t i_track = 0; i_track < N_tracks_off; ++i_track) // loop over all tracks of the actual event
    {
        AS_Track      = AS_Event ->getTrack( i_track ); // take the track
        UShort_t  fNumOfflineTracklets = AS_Track ->getNumOfflineTracklets();
        if(fNumOfflineTracklets > 0)
        {
            for(Int_t i_tracklet = 0; i_tracklet < fNumOfflineTracklets; i_tracklet++) // layers
            {
                AS_offline_Tracklet     = AS_Track            ->getOfflineTracklet( i_tracklet ); // take the track
                TVector3 TV3_offset     = AS_offline_Tracklet ->get_TV3_offset(); // offline tracklets
                TVector3 TV3_dir        = AS_offline_Tracklet ->get_TV3_dir();    // offline tracklets
                Short_t  i_det_tracklet = AS_offline_Tracklet ->get_detector();

                //TV3_offset.Print();
                //TV3_dir.Print();

                Float_t points[6] =
                {
                    (Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),
                    (Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2])
                };

                //printf("  offline tracklets: i_tracklet: %d, posA: {%4.3f, %4.3f}, posB: {%4.3f, %4.3f} \n",i_tracklet,(Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]));
                vec_TPL_offline_tracklets.push_back(new TPolyLine());
                vec_TPL_offline_tracklets[(Int_t)vec_TPL_offline_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]));
                vec_TPL_offline_tracklets[(Int_t)vec_TPL_offline_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]));
            }
        }
    }

    return vec_TPL_offline_tracklets;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_tracklets_line_2D(Int_t i_track)
{
    if (drawing_online_tracklets_flag != 1) get_tracklets_fit(i_track);

    vector<TPolyLine*> vec_TPL_on_trkl  =  get_online_tracklets(i_track);
    vector<TPolyLine*> vec_TPL_off_trkl =  get_offline_tracklets(i_track);

    for(Int_t i_layer = 0; i_layer < 6; i_layer++) // 6 layers + global fit
    {
        if(vec_tracklet_fit_points[i_layer][0][0] > -999.0 && vec_tracklet_fit_points[i_layer][1][0] > -999.0)
        {
            vec_tracklets_line_2D[i_layer] = new TPolyLine();
            vec_tracklets_line_2D[i_layer] ->SetNextPoint(vec_tracklet_fit_points[i_layer][0][0],vec_tracklet_fit_points[i_layer][0][1]);
            vec_tracklets_line_2D[i_layer] ->SetNextPoint(vec_tracklet_fit_points[i_layer][1][0],vec_tracklet_fit_points[i_layer][1][1]);
            vec_tracklets_line_2D[i_layer] ->SetLineStyle(1);
            vec_tracklets_line_2D[i_layer] ->SetLineColor(color_layer[i_layer]);
            vec_tracklets_line_2D[i_layer] ->SetLineWidth(line_width_layer[i_layer]);
            vec_tracklets_line_2D[i_layer] ->DrawClone("l");
            printf("i_layer tracklets line: %d \n", i_layer);
        }
        else
        {
            printf("TBase_TRD_Calib::Draw_tracklets_line(%d), i_layer: %d has no entry \n",i_track,i_layer);
        }
    }

#if 0
    for(Int_t i_onl_trkl = 0; i_onl_trkl < (Int_t)vec_TPL_on_trkl.size(); i_onl_trkl++)
    {
        vec_TPL_on_trkl[i_onl_trkl] ->SetLineStyle(9);
        vec_TPL_on_trkl[i_onl_trkl] ->SetLineWidth(3);
        vec_TPL_on_trkl[i_onl_trkl] ->SetLineColor(kBlack);
        vec_TPL_on_trkl[i_onl_trkl] ->DrawClone("l");
    }

    for(Int_t i_offl_trkl = 0; i_offl_trkl < (Int_t)vec_TPL_off_trkl.size(); i_offl_trkl++)
    {
        //printf(" drawing i_offl_trkl: %d \n",i_offl_trkl);
        vec_TPL_off_trkl[i_offl_trkl] ->SetLineStyle(9);
        vec_TPL_off_trkl[i_offl_trkl] ->SetLineWidth(3);
        vec_TPL_off_trkl[i_offl_trkl] ->SetLineColor(kCyan+1);
        vec_TPL_off_trkl[i_offl_trkl] ->DrawClone("l");
    }
#endif
    ;
    select_online_tracklets();

    vec_TPL_online_tracklets_selected.clear();
    vec_TPL_online_tracklets_selected.resize(6);

    //create TPL for selected online tracklets and draw them
    for(Int_t i_layer = 0; i_layer < 6; i_layer++) // 6 layers of online tracklets
    {
        if(vec_online_tracklet_points[i_layer][0][0] > -999.0 && vec_online_tracklet_points[i_layer][1][0] > -999.0)
        {
            vec_TPL_online_tracklets_selected[i_layer] = new TPolyLine();
            vec_TPL_online_tracklets_selected[i_layer] ->SetNextPoint(vec_online_tracklet_points[i_layer][0][0],vec_online_tracklet_points[i_layer][0][1]);
            vec_TPL_online_tracklets_selected[i_layer] ->SetNextPoint(vec_online_tracklet_points[i_layer][1][0],vec_online_tracklet_points[i_layer][1][1]);
            vec_TPL_online_tracklets_selected[i_layer] ->SetLineStyle(9);
            vec_TPL_online_tracklets_selected[i_layer] ->SetLineColor(kBlack);
            vec_TPL_online_tracklets_selected[i_layer] ->SetLineWidth(3);
            vec_TPL_online_tracklets_selected[i_layer] ->DrawClone("l");
            printf("i_layer tracklets online: %d \n", i_layer);
        }
        else
        {
            printf("TBase_TRD_Calib::Draw_tracklets_line(%d), i_layer: %d has no entry \n",i_track,i_layer);
        }
    }


    HistName = "Ev. ";
    HistName += Event_active;
    HistName += ", tr. ";
    HistName += i_track;
    HistName += Form(", min: (%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f, gl. %4.3f, circle_gl. %4.3f, online_circle_gl. %4.3f)",
        tracklets_min[0],tracklets_min[1],tracklets_min[2],tracklets_min[3],tracklets_min[4],tracklets_min[5],tracklets_min[6],tracklets_min_circle, tracklets_min_online_circle);

    //sprintf(NoP,"%4.0f",(Double_t)i_detector);
    //HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.1,0.96,0.025,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    HistName = Form("det: (%d, %d, %d, %d, %d, %d)",arr_layer_detector[0],arr_layer_detector[1],arr_layer_detector[2],arr_layer_detector[3],arr_layer_detector[4],arr_layer_detector[5]);
    plotTopLegend((char*)HistName.Data(),0.24,0.92,0.03,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    Double_t arr_HV_anode[6] = {-1.0};
    Double_t arr_HV_drift[6] = {-1.0};
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        Double_t det;
        tg_HV_anode_vs_det ->GetPoint(arr_layer_detector[i_layer],det,arr_HV_anode[i_layer]);
        tg_HV_drift_vs_det ->GetPoint(arr_layer_detector[i_layer],det,arr_HV_drift[i_layer]);
    }
    HistName = Form("HVA: (%4.0f, %4.0f, %4.0f, %4.0f, %4.0f, %4.0f)",arr_HV_anode[0],arr_HV_anode[1],arr_HV_anode[2],arr_HV_anode[3],arr_HV_anode[4],arr_HV_anode[5]);
    plotTopLegend((char*)HistName.Data(),0.24,0.89,0.03,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    HistName = Form("HVD: (%4.0f, %4.0f, %4.0f, %4.0f, %4.0f, %4.0f)",arr_HV_drift[0],arr_HV_drift[1],arr_HV_drift[2],arr_HV_drift[3],arr_HV_drift[4],arr_HV_drift[5]);
    plotTopLegend((char*)HistName.Data(),0.24,0.86,0.03,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    //tg_HV_drift_vs_det = (TGraph*)inputfile_QA->Get("tg_HV_drift_vs_det");
    //tg_HV_anode_vs_det = (TGraph*)inputfile_QA->Get("tg_HV_anode_vs_det");;
    //tg_vdrift_vs_det   = (TGraph*)inputfile_QA->Get("tg_vdrift_vs_det");;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_tracklets_line(Int_t i_track)
{
    get_tracklets_fit(i_track);

    for(Int_t i_layer = 0; i_layer < 7; i_layer++) // 6 layers + global fit
    {
        if(vec_tracklet_fit_points[i_layer][0][0] > -999.0 && vec_tracklet_fit_points[i_layer][1][0] > -999.0)
        {
            if(vec_tracklets_line) delete vec_tracklets_line;
            vec_tracklets_line = new TEveLine();
            //vec_tracklets_line ->Clear();
            HistName = "tracklet_track";
            HistName += i_track;
            HistName += "_layer";
            HistName += i_layer;
            vec_tracklets_line ->SetName(HistName.Data());
            vec_tracklets_line ->SetNextPoint(vec_tracklet_fit_points[i_layer][0][0],vec_tracklet_fit_points[i_layer][0][1],vec_tracklet_fit_points[i_layer][0][2]);
            vec_tracklets_line ->SetNextPoint(vec_tracklet_fit_points[i_layer][1][0],vec_tracklet_fit_points[i_layer][1][1],vec_tracklet_fit_points[i_layer][1][2]);
            //printf("i_layer_notemtpy: %d, layer: %d \n",i_layer_notempty,vec_layer_in_fit[i_layer_notempty]);
            vec_tracklets_line ->SetLineStyle(1);
            vec_tracklets_line ->SetLineColor(color_layer[i_layer]);
            vec_tracklets_line ->SetLineWidth(line_width_layer[i_layer]);
            //vec_tracklets_line ->DrawClone("ogl");
            gEve->AddElement((TEveLine*)vec_tracklets_line->Clone());
            printf("i_layer tracklets line: %d \n", i_layer);
        }
        else
        {
            printf("TBase_TRD_Calib::Draw_tracklets_line(%d), i_layer: %d has no entry \n",i_track,i_layer);
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TEveLine* TBase_TRD_Calib::get_helix_polyline(Int_t i_track)
{
    TEveLine* TPL3D_helix_track = new TEveLine();
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
        TPL3D_helix_track ->SetNextPoint(helix_point[0],helix_point[1],helix_point[2]);
        //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f}, radius: %4.3f \n",i_step,helix_point[0],helix_point[1],helix_point[2],radius);
        if(radius > 368.0) break;
    }

    return TPL3D_helix_track;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TPolyLine* TBase_TRD_Calib::get_helix_polyline_2D(Int_t i_track)
{
    printf("TBase_TRD_Calib::get_helix_polyline_2D \n");
    TPolyLine* TPL_helix_track = new TPolyLine();
    AS_Track      = AS_Event ->getTrack( i_track ); // take the track

    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
        aliHelix.fHelix[i_param] = AS_Track ->getHelix_param(i_param);
        //printf("i_param: %d, param: %4.3f \n",i_param,AS_Track ->getHelix_param(i_param));
    }

    Double_t helix_point[3];
    Double_t pathA = 0.0;
    Double_t radius = 0.0;
    for(Int_t i_step = 0; i_step < 400; i_step++)
    {
        pathA = i_step*3.0;
        aliHelix.Evaluate(pathA,helix_point);

        radius = TMath::Sqrt(TMath::Power(helix_point[0],2.0) + TMath::Power(helix_point[1],2.0));
        if(radius < 250.0) continue;
        TPL_helix_track ->SetNextPoint(helix_point[0],helix_point[1]);
        //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f}, radius: %4.3f \n",i_step,helix_point[0],helix_point[1],helix_point[2],radius);
        if(radius > 368.0) break;
    }

    return TPL_helix_track;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_track(Int_t i_track)
{
    printf("TBase_TRD_Calib::Draw_track \n");
    TPL3D_helix = get_helix_polyline(i_track);

    TPL3D_helix    ->SetLineStyle(1);
    //TPL3D_helix    ->SetLineColor(kRed);
    TPL3D_helix    ->SetLineWidth(5);
    TPL3D_helix    ->DrawClone("ogl");
    TPL3D_helix->SetMainColor(kRed);
    gEve->AddElement(TPL3D_helix);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_2D_track(Int_t i_track){
    printf("TBase_TRD_Calib::Draw_2D_track \n");


    TPL_helix = get_helix_polyline_2D(i_track);

    Int_t n_points   = TPL_helix->GetN();
    Double_t *x_vals = TPL_helix->GetX();
    Double_t *y_vals = TPL_helix->GetY();

    Double_t x_start, y_start, x_stop, y_stop;

    Int_t n_points_use = 0;
    for(Int_t i_point = 0; i_point < n_points; i_point++)
    {
        if(x_vals[i_point] == 0.0 && y_vals[i_point] == 0.0) continue;
        if(n_points_use == 0)
        {
            x_start = x_vals[i_point];
            x_stop  = x_vals[i_point];
            y_start = y_vals[i_point];
            y_stop  = y_vals[i_point];
        }
        else
        {
            if(x_vals[i_point] < x_start) x_start = x_vals[i_point];
            if(y_vals[i_point] < y_start) y_start = y_vals[i_point];
            if(x_vals[i_point] > x_stop)  x_stop = x_vals[i_point];
            if(y_vals[i_point] > y_stop)  y_stop = y_vals[i_point];
        }
        n_points_use++;
        //printf("i_point: %d, pos: {%4.3f, %4.3f} \n",i_point,x_vals[i_point],y_vals[i_point]);
    }

    TCanvas* can_2D_track = new TCanvas("can_2D_track","can_2D_track",10,10,1200,900);
    can_2D_track ->cd();
    can_2D_track ->cd()->SetRightMargin(0.01);
    can_2D_track ->cd()->SetTopMargin(0.05);
    can_2D_track ->cd()->SetBottomMargin(0.2);
    can_2D_track ->cd()->SetLeftMargin(0.2);
    TH1F* h_frame = can_2D_track->cd()->DrawFrame(x_start-5.0,y_start-5.0,x_stop+5.0,y_stop+5.0,"h_frame");
    h_frame->SetStats(0);
    h_frame->SetTitle("");
    h_frame->GetXaxis()->SetTitleOffset(1.1);
    h_frame->GetYaxis()->SetTitleOffset(1.8);
    h_frame->GetXaxis()->SetLabelOffset(0.0);
    h_frame->GetYaxis()->SetLabelOffset(0.01);
    h_frame->GetXaxis()->SetLabelSize(0.05);
    h_frame->GetYaxis()->SetLabelSize(0.05);
    h_frame->GetXaxis()->SetTitleSize(0.05);
    h_frame->GetYaxis()->SetTitleSize(0.05);
    h_frame->GetXaxis()->SetNdivisions(505,'N');
    h_frame->GetYaxis()->SetNdivisions(505,'N');
    h_frame->GetXaxis()->CenterTitle();
    h_frame->GetYaxis()->CenterTitle();
    h_frame->GetXaxis()->SetTitle("x (cm)");
    h_frame->GetYaxis()->SetTitle("y (cm)");

    TPL_helix ->SetLineColor(kRed);
    TPL_helix ->SetLineStyle(1);
    TPL_helix ->SetLineWidth(2);
    TPL_helix ->DrawClone("l");

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
        vec_TPL3D_helix_neighbor[i_track_neighbor]->SetLineStyle(1);
        vec_TPL3D_helix_neighbor[i_track_neighbor]->SetLineWidth(2);
        vec_TPL3D_helix_neighbor[i_track_neighbor]->SetMainColor(kGreen);
        gEve->AddElement(vec_TPL3D_helix_neighbor[i_track_neighbor]);
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
    // TString pinputdir = "/home/jasonb/Documents/masters/calibration/";
    // TString pinputdir = "/misc/alidata120/alice_u/berdnikova/TRD-Run3-Calibration/3456/";

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

#if 0
    //----------------
    // Check which event has offline tracklet information from friends
    for(Int_t i_ev = 0; i_ev < 614; i_ev++)
    {
        input_SE->GetEntry( i_ev );
        if((i_ev % 100) == 0) printf("i_ev: %d \n",i_ev);
        UShort_t N_tr            = AS_Event ->getNumTracks(); // number of tracks in this event
        for(UShort_t i_track = 0; i_track < N_tr; ++i_track) // loop over all tracks of the actual event
        {
            AS_Track      = AS_Event ->getTrack( i_track ); // take the track
            UShort_t  fNumOfflineTracklets = AS_Track ->getNumOfflineTracklets();
            if(fNumOfflineTracklets > 0)
            {
                printf("    ------> i_ev: %d, i_track: %d, N_off_trkl: %d \n",i_ev,i_track,fNumOfflineTracklets);
            }
        }
    }
    //----------------
#endif

    Event_active = event;

    if (!input_SE->GetEntry( event )) return 0; // take the event -> information is stored in event

    N_Digits = 0;

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        delete vec_TPM3D_digits[i_layer];
        vec_TPM3D_digits[i_layer] = new TEvePointSet();
    }

    vec_TPL3D_tracklets.clear();
    vec_TPL3D_offline_tracklets.clear();

    //---------------------------------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event
    Double_t EventVertexX         = AS_Event ->getx();
    Double_t EventVertexY         = AS_Event ->gety();
    Double_t EventVertexZ         = AS_Event ->getz();
    Int_t    N_tracks_event       = AS_Event ->getN_tracks();
    Int_t    N_TRD_tracklets      = AS_Event ->getN_TRD_tracklets();
    Int_t    N_TRD_tracklets_online = AS_Event ->getNumTracklets(); // online tracklet
    Float_t  V0MEq                = AS_Event ->getcent_class_V0MEq();

    printf("N_TRD_tracklets_online: %d \n",N_TRD_tracklets_online);

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
    for(UShort_t i_tracklet = 0; i_tracklet < N_TRD_tracklets_online; ++i_tracklet) // loop over all tracklets of the actual event
    {
        AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
        TVector3 TV3_offset     = AS_Tracklet ->get_TV3_offset(); // online tracklets
        TVector3 TV3_dir        = AS_Tracklet ->get_TV3_dir();    // online tracklets
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
        //printf("i_tracklet: %d, out of %d, impact_angle: %4.3f, offset: {%4.3f, %4.3f, %4.3f}, end: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,N_TRD_tracklets_online,impact_angle*TMath::RadToDeg(),TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + TV3_dir[0],TV3_offset[1] + TV3_dir[1],TV3_offset[2] + TV3_dir[2]);

        vec_TPL3D_tracklets.push_back(new TEveLine());
        vec_TPL3D_tracklets[(Int_t)vec_TPL3D_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]));
        vec_TPL3D_tracklets[(Int_t)vec_TPL3D_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2]));
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
        Float_t phi_track       = TLV_part.Phi();

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
        UShort_t  fNumTRDdigits        = AS_Track ->getNumTRD_digits();
        UShort_t  fNumOfflineTracklets = AS_Track ->getNumOfflineTracklets();
        //--------------------------


        //--------------------------
        // Offline tracklet loop
        for(Int_t i_tracklet = 0; i_tracklet < fNumOfflineTracklets; i_tracklet++) // layers
        {
            AS_offline_Tracklet     = AS_Track            ->getOfflineTracklet( i_tracklet ); // take the track
            TVector3 TV3_offset     = AS_offline_Tracklet ->get_TV3_offset(); // offline tracklets
            TVector3 TV3_dir        = AS_offline_Tracklet ->get_TV3_dir();    // offline tracklets
            Short_t  i_det_tracklet = AS_offline_Tracklet ->get_detector();

            printf("offline, i_tracklet: %d, offset: {%4.3f, %4.3f, %4.3f}, dir: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,(Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),(Float_t)TV3_dir[0],(Float_t)TV3_dir[1],(Float_t)TV3_dir[2]);

            vec_TPL3D_offline_tracklets.push_back(new TEveLine());
            vec_TPL3D_offline_tracklets[(Int_t)vec_TPL3D_offline_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]));
            vec_TPL3D_offline_tracklets[(Int_t)vec_TPL3D_offline_tracklets.size()-1] ->SetNextPoint((Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2]));
        }
        //--------------------------


        printf("i_track: %d, pT: %4.3f, phi: %4.3f, eta: %4.3f, N_off_trkl: %d \n",i_track,pT_track,phi_track,eta_track,fNumOfflineTracklets);

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
    // for online tracklets -> track reconstruction
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

                //if(!(i_detector >= 378 && i_detector <= 382)) continue;

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

                                vec_TPL3D_tracklets_match.push_back(new TEveLine());
                                vec_TPL3D_tracklets_match[(Int_t)vec_TPL3D_tracklets_match.size()-1] ->SetNextPoint((Float_t)(TV3_pos_B[0]),(Float_t)(TV3_pos_B[1]),(Float_t)(TV3_pos_B[2]));
                                vec_TPL3D_tracklets_match[(Int_t)vec_TPL3D_tracklets_match.size()-1] ->SetNextPoint((Float_t)(TV3_pos_B[0] + scale_factor_length*TV3_dir_B[0]),(Float_t)(TV3_pos_B[1] + scale_factor_length*TV3_dir_B[1]),(Float_t)(TV3_pos_B[2] + scale_factor_length*TV3_dir_B[2]));
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

    // used to check if detector is bad for 2D hists
    TFile* calibration_params = TFile::Open("./Data/TRD_Calib_vDfit_and_LAfit_digits.root");
    vdrift_fit = (TGraph*)calibration_params->Get(";1");

    vec_tp_Delta_vs_impact.resize(540);
    vec_TH2D_Delta_vs_impact.resize(540);
    vec_tp_Delta_vs_impact_circle.resize(540);
    vec_TH2D_Delta_vs_impact_circle.resize(540);

    vec_TH2D_delta_offset_vs_impact_circle.resize(540);

    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_%d",i_det),Form("vec_th1d_Delta_vs_impact_%d",i_det),360,-360,360);
        vec_TH2D_Delta_vs_impact[i_det] = new TH2D(Form("vec_th2d_Delta_vs_impact_%d",i_det),Form("vec_th2d_Delta_vs_impact_%d",i_det),10,70,110,50,-25,25);
        vec_tp_Delta_vs_impact_circle[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),360,-360,360);
        vec_TH2D_Delta_vs_impact_circle[i_det] = new TH2D(Form("vec_th2d_Delta_vs_impact_circle_%d",i_det),Form("vec_th2d_Delta_vs_impact_circle_%d",i_det),13,70,110,200,-100,100);

        vec_TH2D_delta_offset_vs_impact_circle[i_det] = new TH2D(Form("vec_TH2D_delta_offset_vs_impact_circle%d",i_det),Form("vec_TH2D_delta_offset_vs_impact_circle%d",i_det),13,70,110,600,-30,30);
    }

    for(Long64_t i_event = 0; i_event < 50000; i_event++)
    //for(Long64_t i_event = 0; i_event < file_entries_total; i_event++)
    {
        if(i_event % 100 == 0) printf("i_event: %lld out of %lld \n",i_event,file_entries_total);
        if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

        UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event
        for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        {
            // cout << "i_track: " << i_track << ", of " << NumTracks << endl;
            AS_Track      = AS_Event ->getTrack( i_track ); // take the track
            TLorentzVector TLV_part = AS_Track ->get_TLV_part();
            Float_t pT_track        = TLV_part.Pt();
            // if(pT_track < 3.5) continue; // 3.5 
            if(pT_track < 1.5) continue; // 3.5

            Double_t nsigma_TPC_e   = AS_Track ->getnsigma_e_TPC();
            Double_t nsigma_TPC_pi  = AS_Track ->getnsigma_pi_TPC();
            Double_t nsigma_TPC_p   = AS_Track ->getnsigma_p_TPC();
            Double_t nsigma_TOF_e   = AS_Track ->getnsigma_e_TOF();
            Double_t nsigma_TOF_pi  = AS_Track ->getnsigma_pi_TOF();
            Double_t TRD_signal     = AS_Track ->getTRDSignal();
            Double_t TRDsumADC      = AS_Track ->getTRDsumADC();
            Double_t dca            = AS_Track ->getdca();  // charge * distance of closest approach to the primary vertex
            UShort_t NTPCcls        = AS_Track ->getNTPCcls();
            UShort_t NTRDcls        = AS_Track ->getNTRDcls();
            UShort_t NITScls        = AS_Track ->getNITScls();
            Float_t  TPCchi2        = AS_Track ->getTPCchi2();
            Float_t  TPCdEdx        = AS_Track ->getTPCdEdx();
            Float_t  TOFsignal      = AS_Track ->getTOFsignal(); // in ps (1E-12 s)
            Float_t  Track_length   = AS_Track ->getTrack_length();
            Float_t  Angle_on_TRD   = AS_Track ->getimpact_angle_on_TRD();

            Float_t  momentum       = TLV_part.P();
            Float_t  eta_track      = TLV_part.Eta();
            Float_t  theta_track    = TLV_part.Theta();


            //--------------------------
            // Offline tracklet loop

#if 0
            UShort_t  fNumOfflineTracklets = AS_Track ->getNumOfflineTracklets();
            if(fNumOfflineTracklets > 2)
            {
                for(Int_t i_tracklet = 0; i_tracklet < fNumOfflineTracklets; i_tracklet++) // layers
                {
                    AS_offline_Tracklet     = AS_Track            ->getOfflineTracklet( i_tracklet ); // take the track
                    TVector3 TV3_offset     = AS_offline_Tracklet ->get_TV3_offset(); // offline tracklets
                    TVector3 TV3_dir        = AS_offline_Tracklet ->get_TV3_dir();    // offline tracklets
                    Short_t  i_det_tracklet = AS_offline_Tracklet ->get_detector();
                    Float_t  chi2           = AS_offline_Tracklet ->get_chi2();
                    Float_t  refZ           = AS_offline_Tracklet ->get_refZ();
                    Float_t  refY           = AS_offline_Tracklet ->get_refY();
                    Float_t  refdZdx        = AS_offline_Tracklet ->get_refdZdx();
                    Float_t  refdYdx        = AS_offline_Tracklet ->get_refdYdx();
                    Float_t  locZ           = AS_offline_Tracklet ->get_locZ();
                    Float_t  locY           = AS_offline_Tracklet ->get_locY();
                    Float_t  locdZdx        = AS_offline_Tracklet ->get_locdZdx();
                    Float_t  locdYdx        = AS_offline_Tracklet ->get_locdYdx();

                    Int_t i_det_layer = i_det_tracklet%6;

                    if(i_det_tracklet < 0 || i_det_tracklet >= 540) continue;
                    Float_t diff_dYdx = locdYdx - refdYdx;

                    vec_h_diff_ref_loc_off_trkl[0][i_det_tracklet] ->Fill(diff_dYdx);
                    vec_h_diff_ref_loc_off_trkl[0][540] ->Fill(diff_dYdx);
                    vec_h_diff_ref_loc_off_trkl[0][541+i_det_layer] ->Fill(diff_dYdx);
                    if(dca < 0.0)
                    {
                        vec_h_diff_ref_loc_off_trkl[1][i_det_tracklet] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[1][540] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[1][541+i_det_layer] ->Fill(diff_dYdx);
                    }
                    if(dca > 0.0)
                    {
                        vec_h_diff_ref_loc_off_trkl[2][i_det_tracklet] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[2][540] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[2][541+i_det_layer] ->Fill(diff_dYdx);
                    }
                }
            }
#endif
            //--------------------------

            //if(dca > 0.0) continue;


            for(Int_t i_param = 0; i_param < 9; i_param++)
            {
                aliHelix.fHelix[i_param] = AS_Track ->getHelix_param(i_param);
            }
            Double_t pathA = 0.0;
            Double_t helix_point[3];
            Double_t helix_pointB[3];

            //printf("i_track: %4.3d, i_event: %4.3d \n",i_track,i_event);

            //printf("test 0 \n");
            make_clusters(i_track);
            //printf("test 1 \n");
            //get_straight_line_fit(i_track);
            get_tracklets_fit(i_track);
            //printf("test 2 \n");
            get_2D_global_circle_fit();  // i will get dir_vec_circle (not needed?) and vec_TPL_circle_tracklets[i_point] (basically i_layer)
            //printf("test 3 \n");

            //vec_tracklet_fit_points[i_layer][i_start_stop][i_xyz]   including tracklets[0-5] and glibal line fit [6]  vector< vector< vector< Double_t> > >
            //vec_TPL_circle_tracklets[i_layer]    vector< TPolyLine*> ??
            //dir_vec_circle   just a number
            //after all use vec_dir_ver_circle[i_point] :  vector<TVector3>



            vector<TVector3> vec_TV3_tracklet_vectors;
            vec_TV3_tracklet_vectors.resize(7);


            if(vec_tracklet_fit_points[6][0][0] == -999.0  && vec_tracklet_fit_points[6][1][0] == -999.0) continue; // no global fit available
            if(tracklets_min[6] > 7.5) continue;  // 7.5

            //printf("Track with pT: %4.3f used for calibration \n",pT_track);

            Int_t N_good_layers = 0;
            for(Int_t i_layer = 5; i_layer >= 0; i_layer--)
            {
                if(vec_tracklet_fit_points[i_layer][0][0] > -999.0 && vec_tracklet_fit_points[i_layer][1][0] > -999.0) N_good_layers++;
            }
            if(N_good_layers < 3) continue;    //  <5

            for (Int_t layer = 0; layer < 6; layer++)
            {
                if (arr_layer_detector[layer] == 229)
                {

                    printf("N_good_layers: %d \n",N_good_layers);
                    printf("i_event: %lld, i_track: %d \n",i_event,i_track);
                }
            }



            Double_t impact_angle[6] = {0.0};
            Double_t impact_angle_circle[6] = {0.0};

            Double_t delta_x_local_global[6] = {0.0}; // local chamber coordinate system, global fit
            Double_t delta_x_local_global_circle[6] = {0.0}; // local chamber coordinate system, global circle fit - in case it's different


            // cout << endl << "CALIBRATE():" << endl;

            for(Int_t i_layer = 6; i_layer >= 0; i_layer--)
            {
                if(vec_tracklet_fit_points[i_layer][0][0] > -999.0 && vec_tracklet_fit_points[i_layer][1][0] > -999.0)
                {
                    Double_t delta_x = vec_tracklet_fit_points[i_layer][1][0] - vec_tracklet_fit_points[i_layer][0][0];
                    Double_t delta_y = vec_tracklet_fit_points[i_layer][1][1] - vec_tracklet_fit_points[i_layer][0][1];
                    Double_t delta_z = vec_tracklet_fit_points[i_layer][1][2] - vec_tracklet_fit_points[i_layer][0][2];

                    //printf("i_layer: %d, delta_z: %4.3f \n",i_layer,delta_z);
                    TVector3 TV3_delta(delta_x,delta_y,delta_z);   //do i need it for circle..


                    vec_TV3_tracklet_vectors[i_layer].SetXYZ(delta_x,delta_y,0.0);  //here we finally get vector: consists of xStop-xStart; yStop-yStart
                    //printf("i_layer: %d, delta_x: %4.3f, delta_y: %4.3f \n",i_layer,delta_x,delta_y);

                    if(vec_TV3_tracklet_vectors[i_layer].Mag() <= 0.0)
                    {
                        //printf("A: {%4.3f, %4.3f}, B: {%4.3f, %4.3f} \n",vec_tracklet_fit_points[i_layer][1][0],vec_tracklet_fit_points[i_layer][1][1],vec_tracklet_fit_points[i_layer][0][0],vec_tracklet_fit_points[i_layer][0][1]);
                        continue;
                    }

                    vec_TV3_tracklet_vectors[i_layer] *= 1.0/vec_TV3_tracklet_vectors[i_layer].Mag(); // Normalize


                    // Calculate impact angles of global track to every single layer
                    if(i_layer == 6) // Global track
                    {
                        continue;   //try to disable
#if 0
                        Double_t path_offset = 0.0;
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
                            impact_angle[i_layerB] = vec_TV3_tracklet_vectors[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2]);  //i need something like that for circles
                            if(impact_angle[i_layerB] > TMath::Pi()*0.5) impact_angle[i_layerB] -= TMath::Pi();
                            //printf("i_layerB: %d, 3D_Angle_on_TRD(TPC track): %4.3f, 2D_impact_angle: %4.3f \n",i_layerB,Angle_on_TRD*TMath::RadToDeg(),impact_angle[i_layerB]*TMath::RadToDeg());
                            impact_angle[i_layerB] = 0.5*TMath::Pi() - sign_direction_impactB*impact_angle[i_layerB];
                            //printf("impact_angle: %4.3f \n",impact_angle[i_layerB]);


                            //---------------------------
                            // Calculate impact angle based on track helix
                            Double_t radius_layer = vec_TV3_TRD_center_offset[i_layerB].Perp();
                            Double_t impact_angle_helix = -999.0;
                            for(Int_t i_step = 0; i_step < 400; i_step++)
                            {
                                pathA = i_step*3.0 + path_offset;
                                aliHelix.Evaluate(pathA,helix_point);
                                Double_t radius_helix_point = TMath::Sqrt(TMath::Power(helix_point[0],2) + TMath::Power(helix_point[1],2));

                                if(radius_helix_point > radius_layer)
                                {
                                    aliHelix.Evaluate(pathA+0.5,helix_pointB);
                                    TVector3 dir_helix_impact_on_layer(helix_pointB[0] - helix_point[0],helix_pointB[1] - helix_point[1],0.0);
                                    Double_t delta_x_local_global_helix = dir_helix_impact_on_layer.Dot(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][0]);
                                    Double_t sign_direction_impact_helix = TMath::Sign(1.0,delta_x_local_global_helix);
                                    impact_angle_helix = dir_helix_impact_on_layer.Angle(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2]);
                                    if(impact_angle_helix > TMath::Pi()*0.5) impact_angle_helix -= TMath::Pi();
                                    impact_angle_helix = 0.5*TMath::Pi() - sign_direction_impact_helix*impact_angle_helix;
                                    //impact_angle[i_layerB] = impact_angle_helix; // test to use TPC impact angle instead of TRD one

                                    path_offset = pathA; // start from previous layer radius
                                    break;
                                }
                            }
                            vec_h_diff_helix_line_impact_angle[0][arr_layer_detector[i_layerB]] ->Fill(impact_angle[i_layerB]*TMath::RadToDeg() - impact_angle_helix*TMath::RadToDeg());
                            if(dca < 0.0) vec_h_diff_helix_line_impact_angle[1][arr_layer_detector[i_layerB]] ->Fill(impact_angle[i_layerB]*TMath::RadToDeg() - impact_angle_helix*TMath::RadToDeg());
                            if(dca > 0.0) vec_h_diff_helix_line_impact_angle[2][arr_layer_detector[i_layerB]] ->Fill(impact_angle[i_layerB]*TMath::RadToDeg() - impact_angle_helix*TMath::RadToDeg());
                            //printf("impact_angle: %4.3f, impact_angle_helix: %4.3f \n",impact_angle[i_layerB]*TMath::RadToDeg(),impact_angle_helix*TMath::RadToDeg());
                            //---------------------------


                            //if(impact_angle[i_layerB]*TMath::RadToDeg() > 89.0 && impact_angle[i_layerB]*TMath::RadToDeg() < 91.0)
                            //{
                            //    printf("  --> i_layerB: %d, impact_angle: %4.3f, vec: {%4.3f, %4.3f} \n",i_layerB,impact_angle[i_layerB]*TMath::RadToDeg(),vec_tracklet_fit_points[i_layerB][0][0],vec_tracklet_fit_points[i_layerB][1][0]);
                            //}
                        }
#endif
                    }

                    //printf("test 3 \n");

                    if(i_layer < 6) // tracklets + impact angle for circle fit
                    {
                        if(tracklets_min[i_layer] > 0.2) continue;  
                        if(arr_layer_detector[i_layer] < 0) continue;
                        Int_t detector = arr_layer_detector[i_layer];

                        // if (detector == 16)
                        // {
                        //     cout << endl << "DET 16 | " << "event: " << i_event << " | " << "track: " << i_track << endl;
                        // }

                        //printf("test 3.1 \n");

                        //-----impact angle for circle global fit HERE---------

                        //if(tracklets_min_circle == 0.0 || tracklets_min_circle > 4.0) continue;  // 3 points excluded
                        if(tracklets_min_circle > 1.0) continue; 

                        //printf("test 3.2 \n");

                        //printf("i_layer: %d \n",i_layer);

                        //printf("Global straight line dir: {%4.3f, %4.3f, %4.3f} \n",vec_TV3_tracklet_vectors[6].X(),vec_TV3_tracklet_vectors[6].Y(),vec_TV3_tracklet_vectors[6].Z());
                        //printf("Circle dir: {%4.3f, %4.3f, %4.3f} \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());

                        delta_x_local_global_circle[i_layer] = vec_dir_vec_circle[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);
                        //Double_t proj_radial = TV3_delta.Dot(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2]);    looks like this isn't used

                        //printf("test 3.3 \n");

                        Double_t sign_direction_impact_circle = TMath::Sign(1.0,delta_x_local_global_circle[i_layer]);
                        impact_angle_circle[i_layer] = vec_dir_vec_circle[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layer]][2]);  //i need something like that for circles

                        //printf("test 3.4 \n");

                        if(impact_angle_circle[i_layer] > TMath::Pi()*0.5) impact_angle_circle[i_layer] -= TMath::Pi();
                        //printf("i_layerB: %d, 3D_Angle_on_TRD(TPC track): %4.3f, 2D_impact_angle: %4.3f \n",i_layerB,Angle_on_TRD*TMath::RadToDeg(),impact_angle[i_layerB]*TMath::RadToDeg());
                        impact_angle_circle[i_layer] = 0.5*TMath::Pi() - sign_direction_impact_circle*impact_angle_circle[i_layer];
                        //printf("impact_angle_circle: %4.3f \n",impact_angle_circle[i_layer]);


                        //-----------------------------------------------------

                        //printf("test 4 \n");


                        Double_t delta_x_local_tracklet = vec_TV3_tracklet_vectors[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);

                        Double_t sign_angle        = 1.0;
                        Double_t sign_angle_circle = 1.0;
                        if(delta_x_local_tracklet < delta_x_local_global[i_layer])        sign_angle        = -1.0;
                        if(delta_x_local_tracklet < delta_x_local_global_circle[i_layer]) sign_angle_circle = -1.0;

                        //printf("delta_x_local_tracklet: %4.3f, delta_x_local_global: %4.3f, delta_x_local_global_circle: %4.3f \n",delta_x_local_tracklet,delta_x_local_global[i_layer],delta_x_local_global_circle[i_layer]);


                        //calculate Delta angle for straight line
                        Double_t Delta_angle  = sign_angle*vec_TV3_tracklet_vectors[6].Angle(vec_TV3_tracklet_vectors[i_layer]);
                        if(Delta_angle > TMath::Pi()*0.5)  Delta_angle -= TMath::Pi();
                        if(Delta_angle < -TMath::Pi()*0.5) Delta_angle += TMath::Pi();

                        //printf("Delta_angle: %4.3f \n",Delta_angle);


                        //same for delta angle tracklet vs circle
                        Double_t Delta_angle_circle  = sign_angle_circle*vec_dir_vec_circle[i_layer].Angle(vec_TV3_tracklet_vectors[i_layer]);
                        if(Delta_angle_circle > TMath::Pi()*0.5)  Delta_angle_circle -= TMath::Pi();
                        if(Delta_angle_circle < -TMath::Pi()*0.5) Delta_angle_circle += TMath::Pi();
                        //printf("Delta_angle_circle: %4.3f \n",Delta_angle_circle);

                        for (Int_t layer = 0; layer < 6; layer++)
                        {
                            if (arr_layer_detector[layer] == 229) printf("i_layer: %d, Delta_angle_circle: %4.3f \n",i_layer,Delta_angle_circle);
                        }



                        //printf("test 5 \n");


                        h_detector_hit ->Fill(detector);
                        //if(detector == 69) printf("impact_angle: %4.3f, Delta_angle: %4.3f \n",impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());

                        //vec_tp_Delta_vs_impact[arr_layer_detector[i_layer]] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                        //vec_tp_Delta_vs_impact[detector] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),fabs(Delta_angle*TMath::RadToDeg()));

#if 0
                        //fill histos for straight line global track
                        vec_tp_Delta_vs_impact[detector]   ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                        vec_TH2D_Delta_vs_impact[detector] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
#endif
                        //fill histos for circle global track
                        vec_tp_Delta_vs_impact_circle[detector]   ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());
                        vec_TH2D_Delta_vs_impact_circle[detector] ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());

                        // cout << "layer:  " << i_layer << " | " << "offset diff: " << vec_delta_offset[i_layer] << endl;
                        
                        vec_TH2D_delta_offset_vs_impact_circle[detector] ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(),vec_delta_offset[i_layer]);



                        //printf("test 6 \n");

#if 0
                        //separate hist for perpendicular impacts - straight line
                        if(impact_angle[i_layer]*TMath::RadToDeg() > 89.0 && impact_angle[i_layer]*TMath::RadToDeg() < 91.0)
                        //if(impact_angle[i_layer]*TMath::RadToDeg() > 80 && impact_angle[i_layer]*TMath::RadToDeg() < 83)
                        {
                        //    printf("detector: %d, track: %d, impact angle: %4.3f, Delta angle: %4.3f, delta_x_local_tracklet: %4.3f \n",detector,i_track,impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg(),delta_x_local_tracklet);
                            h_delta_angle_perp_impact ->Fill(Delta_angle*TMath::RadToDeg());
                        }
#endif

                        //separate hist for perpendicular impacts - cirrcle
                        if(impact_angle_circle[i_layer]*TMath::RadToDeg() > 89.0 && impact_angle_circle[i_layer]*TMath::RadToDeg() < 91.0)
                            //if(impact_angle[i_layer]*TMath::RadToDeg() > 80 && impact_angle[i_layer]*TMath::RadToDeg() < 83)
                        {
                            h_delta_angle_perp_impact_circle ->Fill(Delta_angle_circle*TMath::RadToDeg());
                        }

#if 0
                        //if(impact_angle[i_layer]*TMath::RadToDeg() > 73.0 && impact_angle[i_layer]*TMath::RadToDeg() < 78.0 && pT_track > 3.5)
                        if(impact_angle[i_layer]*TMath::RadToDeg() > 87.5 && impact_angle[i_layer]*TMath::RadToDeg() < 92.5 && pT_track > 3.5)
                        {
                            //printf("       ======================> i_event: %lld, i_track: %d \n",i_event,i_track);
                        }
                        //printf("detector: %d, track: %d, impact angle: %4.3f, Delta angle: %4.3f \n",detector,i_track,impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
#endif
                    }

                }
            }
        }
    }

    //printf("test 7 \n");

    // copy all that for circle angle histos
    // vector<TCanvas*> vec_can_Delta_vs_impact;
    // vec_can_Delta_vs_impact.resize(6); // 6 sector blocks with 3 sectors each (18)

    // TH1D* h_dummy_Delta_vs_impact = new TH1D("h_dummy_Delta_vs_impact","h_dummy_Delta_vs_impact",90,50,140);
    // h_dummy_Delta_vs_impact->SetStats(0);
    // h_dummy_Delta_vs_impact->SetTitle("");
    // h_dummy_Delta_vs_impact->GetXaxis()->SetTitleOffset(0.85);
    // h_dummy_Delta_vs_impact->GetYaxis()->SetTitleOffset(0.78);
    // h_dummy_Delta_vs_impact->GetXaxis()->SetLabelOffset(0.0);
    // h_dummy_Delta_vs_impact->GetYaxis()->SetLabelOffset(0.01);
    // h_dummy_Delta_vs_impact->GetXaxis()->SetLabelSize(0.08);
    // h_dummy_Delta_vs_impact->GetYaxis()->SetLabelSize(0.08);
    // h_dummy_Delta_vs_impact->GetXaxis()->SetTitleSize(0.08);
    // h_dummy_Delta_vs_impact->GetYaxis()->SetTitleSize(0.08);
    // h_dummy_Delta_vs_impact->GetXaxis()->SetNdivisions(505,'N');
    // h_dummy_Delta_vs_impact->GetYaxis()->SetNdivisions(505,'N');
    // h_dummy_Delta_vs_impact->GetXaxis()->CenterTitle();
    // h_dummy_Delta_vs_impact->GetYaxis()->CenterTitle();
    // h_dummy_Delta_vs_impact->GetXaxis()->SetTitle("impact angle");
    // h_dummy_Delta_vs_impact->GetYaxis()->SetTitle("#Delta #alpha");
    // h_dummy_Delta_vs_impact->GetXaxis()->SetRangeUser(70,110);
    // h_dummy_Delta_vs_impact->GetYaxis()->SetRangeUser(-24,24);

    Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    // for(Int_t i_sec_block = 0; i_sec_block < 6; i_sec_block++)
    // {
    //     HistName = "vec_can_Delta_vs_impact_";
    //     HistName += i_sec_block;
    //     vec_can_Delta_vs_impact[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

    //     vec_can_Delta_vs_impact[i_sec_block] ->Divide(5,3); // x = stack, y = sector

    //     for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
    //     {
    //         Int_t i_sector = i_sec_block + 6*i_sec_sub;
    //         for(Int_t i_stack = 0; i_stack < 5; i_stack++)
    //         {
    //             Int_t iPad = i_sec_sub*5 + i_stack + 1;
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTicks(1,1);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetGrid(0,0);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetFillColor(10);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
    //             vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad);
    //             h_dummy_Delta_vs_impact->Draw("h");

    //             for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    //             {
    //                 Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
    //                 printf("detector: %d \n",i_detector);
    //                 vec_tp_Delta_vs_impact[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
    //                 vec_tp_Delta_vs_impact[i_detector] ->SetLineWidth(2);
    //                 vec_tp_Delta_vs_impact[i_detector] ->SetLineStyle(1);
    //                 vec_tp_Delta_vs_impact[i_detector] ->Draw("same hl");

    //                 HistName = "";
    //                 sprintf(NoP,"%4.0f",(Double_t)i_detector);
    //                 HistName += NoP;
    //                 plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    //             }
    //         }
    //     }
    // }

    // printf("test 8 \n");

    // canvas for circle fit
    vector<TCanvas*> vec_can_Delta_vs_impact_circle;
    vec_can_Delta_vs_impact_circle.resize(6); // 6 sector blocks with 3 sectors each (18)

    TH1D* h_dummy_Delta_vs_impact_circle = new TH1D("h_dummy_Delta_vs_impact_circle","h_dummy_Delta_vs_impact_circle",90,50,140);
    h_dummy_Delta_vs_impact_circle->SetStats(0);
    h_dummy_Delta_vs_impact_circle->SetTitle("");
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact_circle->GetXaxis()->CenterTitle();
    h_dummy_Delta_vs_impact_circle->GetYaxis()->CenterTitle();
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitle("impact angle");
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitle("#Delta #alpha");
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetRangeUser(70,110);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetRangeUser(-24,24);

    // Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    for(Int_t i_sec_block = 0; i_sec_block < 6; i_sec_block++)
    {
        HistName = "vec_can_Delta_vs_impact_circle_";
        HistName += i_sec_block;
        vec_can_Delta_vs_impact_circle[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_can_Delta_vs_impact_circle[i_sec_block] ->Divide(5,3); // x = stack, y = sector

        for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
        {
            Int_t i_sector = i_sec_block + 6*i_sec_sub;
            for(Int_t i_stack = 0; i_stack < 5; i_stack++)
            {
                Int_t iPad = i_sec_sub*5 + i_stack + 1;
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTicks(1,1);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetGrid(0,0);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetFillColor(10);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad);
                h_dummy_Delta_vs_impact_circle->Draw("h");

                for(Int_t i_layer = 0; i_layer < 6; i_layer++)
                {
                    Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                    // printf("detector: %d \n",i_detector);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineWidth(2);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineStyle(1);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->Draw("same hl");

                    HistName = "";
                    sprintf(NoP,"%4.0f",(Double_t)i_detector);
                    HistName += NoP;
                    plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                }
            }
        }
    }


    // 2D plot
    // vector<TCanvas*> vec_can_2D_Delta_vs_impact_circle;
    // vec_can_2D_Delta_vs_impact_circle.resize(36); // 12 sector blocks with 3 sectors each (18)

    // TH1D* h2D_dummy_Delta_vs_impact_circle = new TH1D("h2D_dummy_Delta_vs_impact_circle","h2D_dummy_Delta_vs_impact_circle",90,50,140);
    // h2D_dummy_Delta_vs_impact_circle->SetStats(0);
    // h2D_dummy_Delta_vs_impact_circle->SetTitle("");
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleOffset(0.85);
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleOffset(0.78);
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelOffset(0.0);
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelOffset(0.01);
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelSize(0.08);
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelSize(0.08);
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleSize(0.08);
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleSize(0.08);
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetNdivisions(505,'N');
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetNdivisions(505,'N');
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->CenterTitle();
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->CenterTitle();
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitle("impact angle");
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitle("#Delta #alpha");
    // h2D_dummy_Delta_vs_impact_circle->GetXaxis()->SetRangeUser(70,110);
    // h2D_dummy_Delta_vs_impact_circle->GetYaxis()->SetRangeUser(-99,99);


    // for(Int_t i_sec_block = 0; i_sec_block < 36; i_sec_block++)
    // {
    //     HistName = "vec_can_2D_Delta_vs_impact_circle";
    //     HistName += i_sec_block;
    //     vec_can_2D_Delta_vs_impact_circle[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

    //     vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->Divide(5,3); // x = stack, y = sector


    //     for(Int_t i_det = i_sec_block*15; i_det < i_sec_block*15 + 15; i_det++)
    //     {
    //         Int_t det_colour = kBlack;
    //         bool is_defective = find(begin(Defect_TRD_detectors), end(Defect_TRD_detectors), i_det) != end(Defect_TRD_detectors);
    //         if (is_defective){
    //             det_colour = kRed;
    //         }
            
    //         Double_t vdrift;
    //         Double_t x;
    //         Double_t y;

    //         bool found = 0;
    //         for (Int_t vdrift_det=0; vdrift_det<540; vdrift_det++)
    //         {
    //             vdrift = vdrift_fit->GetPoint(vdrift_det, x, y);
    //             if ((Int_t)x == i_det)
    //             {
    //                 found = 1;
    //             }
    //         }

    //         if (!found)
    //         {
    //             det_colour = kRed;
    //         }
    //         // see all defects
    //         // if (det_colour == kRed)
    //         // {
    //         //     cout << i_det << ", ";
    //         // }
    //         Double_t total = vec_TH2D_Delta_vs_impact_circle[i_det]->Integral(1,11,1,98);
    //         // cout << "det: " << i_det << " | " << "total: " << total << endl;
    //         if (total == 0)
    //         {
    //             continue;
    //         }

    //         Int_t iPad = i_det % 15 + 1;

    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTicks(1,1);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetGrid(0,0);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetFillColor(10);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
    //         vec_can_2D_Delta_vs_impact_circle[i_sec_block] ->cd(iPad);

    //         h2D_dummy_Delta_vs_impact_circle->Draw("h");

    //         vec_TH2D_Delta_vs_impact_circle[i_det] ->Draw("colz same h1");

    //         HistName = "";
    //         sprintf(NoP,"Det: %d",(Int_t)i_det);
    //         HistName += NoP;
    //         plotTopLegend((char*)HistName.Data(),0.75,0.9,0.055,det_colour,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
        
    //     }

    //     HistName = "";
    //     sprintf(NoP,"%s%d","half_sector_", (Int_t)i_sec_block);
    //     HistName += NoP;
    //     HistName += ".png";
    //     vec_can_2D_Delta_vs_impact_circle[i_sec_block]->Print((char*)HistName.Data());
    // }


    TCanvas* can_delta_angle_perp_impact = new TCanvas("can_delta_angle_perp_impact","can_delta_angle_perp_impact",500,10,500,500);
    can_delta_angle_perp_impact ->cd();
    h_delta_angle_perp_impact  ->Draw();

    //same for circle perp impact
    TCanvas* can_delta_angle_perp_impact_circle = new TCanvas("can_delta_angle_perp_impact_circle","can_delta_angle_perp_impact_circle",500,10,500,500);
    can_delta_angle_perp_impact_circle ->cd();
    h_delta_angle_perp_impact_circle  ->Draw();

    TCanvas* can_detector_hit = new TCanvas("can_detector_hit","can_detector_hit",500,500,500,500);
    can_detector_hit ->cd();

    h_detector_hit->SetFillColor(kBlue);
    h_detector_hit  ->Draw();

    TFile* h_detector_hit_outputfile = new TFile("./h_detector_hit.root","RECREATE");
    h_detector_hit_outputfile ->cd();
    h_detector_hit->Write();

    printf("Write data to output file \n");
    outputfile ->cd();
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det]          ->Write();
        vec_TH2D_Delta_vs_impact[i_det]        ->Write();
        //vec_tp_Delta_vs_impact_circle[i_det]   ->Write();
        //vec_TH2D_Delta_vs_impact_circle[i_det] ->Write();
    }

    outputfile ->mkdir("Delta_impact_circle");
    outputfile ->cd("Delta_impact_circle");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact_circle[i_det]   ->Write();
        vec_TH2D_Delta_vs_impact_circle[i_det] ->Write();
    }

    outputfile ->cd();
    outputfile ->mkdir("Delta_offset");
    outputfile ->cd("Delta_offset");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TH2D_delta_offset_vs_impact_circle[i_det] ->Write();
    }


    outputfile ->cd();
    outputfile ->mkdir("Delta_impact");
    outputfile ->cd("Delta_impact");
    for(Int_t i_charge = 0; i_charge < 3; i_charge++)
    {
        for(Int_t i_det = 0; i_det < 547; i_det++)
        {
            vec_h_diff_helix_line_impact_angle[i_charge][i_det] ->Write();
            vec_h_diff_ref_loc_off_trkl[i_charge][i_det]        ->Write();
        }
    }
    printf("All data written \n");

}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Calibrate_on_trkl()
{
    printf("TBase_TRD_Calib::Calibrate_on_trkl() \n");

    vec_tp_Delta_vs_impact.resize(540);
    vec_TH2D_Delta_vs_impact.resize(540);
    vec_tp_Delta_vs_impact_circle.resize(540);
    vec_TH2D_Delta_vs_impact_circle.resize(540);

    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_%d",i_det),Form("vec_th1d_Delta_vs_impact_%d",i_det),360,-360,360);
        vec_TH2D_Delta_vs_impact[i_det] = new TH2D(Form("vec_th2d_Delta_vs_impact_%d",i_det),Form("vec_th2d_Delta_vs_impact_%d",i_det),10,70,110,50,-25,25);
        vec_tp_Delta_vs_impact_circle[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),360,-360,360);
        vec_TH2D_Delta_vs_impact_circle[i_det] = new TH2D(Form("vec_th2d_Delta_vs_impact_circle_%d",i_det),Form("vec_th2d_Delta_vs_impact_circle_%d",i_det),10,70,110,50,-25,25);
    }

    for(Long64_t i_event = 0; i_event < 1000; i_event++)
    //for(Long64_t i_event = 0; i_event < file_entries_total; i_event++)
    {
        if(i_event % 20 == 0) printf("i_event: %lld out of %lld \n",i_event,file_entries_total);
        if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

        UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event

        for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        {
            //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
            AS_Track      = AS_Event ->getTrack( i_track ); // take the track
            TLorentzVector TLV_part = AS_Track ->get_TLV_part();
            Float_t pT_track        = TLV_part.Pt();
            //if(pT_track < 3.5) continue; // 3.5
            if(pT_track < 1.5) continue; // 3.5

            Double_t nsigma_TPC_e   = AS_Track ->getnsigma_e_TPC();
            Double_t nsigma_TPC_pi  = AS_Track ->getnsigma_pi_TPC();
            Double_t nsigma_TPC_p   = AS_Track ->getnsigma_p_TPC();
            Double_t nsigma_TOF_e   = AS_Track ->getnsigma_e_TOF();
            Double_t nsigma_TOF_pi  = AS_Track ->getnsigma_pi_TOF();
            Double_t TRD_signal     = AS_Track ->getTRDSignal();
            Double_t TRDsumADC      = AS_Track ->getTRDsumADC();
            Double_t dca            = AS_Track ->getdca();  // charge * distance of closest approach to the primary vertex
            UShort_t NTPCcls        = AS_Track ->getNTPCcls();
            UShort_t NTRDcls        = AS_Track ->getNTRDcls();
            UShort_t NITScls        = AS_Track ->getNITScls();
            Float_t  TPCchi2        = AS_Track ->getTPCchi2();
            Float_t  TPCdEdx        = AS_Track ->getTPCdEdx();
            Float_t  TOFsignal      = AS_Track ->getTOFsignal(); // in ps (1E-12 s)
            Float_t  Track_length   = AS_Track ->getTrack_length();
            Float_t  Angle_on_TRD   = AS_Track ->getimpact_angle_on_TRD();

            Float_t  momentum       = TLV_part.P();
            Float_t  eta_track      = TLV_part.Eta();
            Float_t  theta_track    = TLV_part.Theta();


#if 0
            //--------------------------
            // Offline tracklet loop
            UShort_t  fNumOnlineTracklets = AS_Track ->getNumOfflineTracklets();
            if(fNumOfflineTracklets > 2)
            {
                for(Int_t i_tracklet = 0; i_tracklet < fNumOfflineTracklets; i_tracklet++) // layers
                {
                    AS_offline_Tracklet     = AS_Track            ->getOfflineTracklet( i_tracklet ); // take the track
                    TVector3 TV3_offset     = AS_offline_Tracklet ->get_TV3_offset(); // offline tracklets
                    TVector3 TV3_dir        = AS_offline_Tracklet ->get_TV3_dir();    // offline tracklets
                    Short_t  i_det_tracklet = AS_offline_Tracklet ->get_detector();
                    Float_t  chi2           = AS_offline_Tracklet ->get_chi2();
                    Float_t  refZ           = AS_offline_Tracklet ->get_refZ();
                    Float_t  refY           = AS_offline_Tracklet ->get_refY();
                    Float_t  refdZdx        = AS_offline_Tracklet ->get_refdZdx();
                    Float_t  refdYdx        = AS_offline_Tracklet ->get_refdYdx();
                    Float_t  locZ           = AS_offline_Tracklet ->get_locZ();
                    Float_t  locY           = AS_offline_Tracklet ->get_locY();
                    Float_t  locdZdx        = AS_offline_Tracklet ->get_locdZdx();
                    Float_t  locdYdx        = AS_offline_Tracklet ->get_locdYdx();

                    Int_t i_det_layer = i_det_tracklet%6;

                    if(i_det_tracklet < 0 || i_det_tracklet >= 540) continue;
                    Float_t diff_dYdx = locdYdx - refdYdx;

                    vec_h_diff_ref_loc_off_trkl[0][i_det_tracklet] ->Fill(diff_dYdx);
                    vec_h_diff_ref_loc_off_trkl[0][540] ->Fill(diff_dYdx);
                    vec_h_diff_ref_loc_off_trkl[0][541+i_det_layer] ->Fill(diff_dYdx);
                    if(dca < 0.0)
                    {
                        vec_h_diff_ref_loc_off_trkl[1][i_det_tracklet] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[1][540] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[1][541+i_det_layer] ->Fill(diff_dYdx);
                    }
                    if(dca > 0.0)
                    {
                        vec_h_diff_ref_loc_off_trkl[2][i_det_tracklet] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[2][540] ->Fill(diff_dYdx);
                        vec_h_diff_ref_loc_off_trkl[2][541+i_det_layer] ->Fill(diff_dYdx);
                    }
                }
            }
            //--------------------------
#endif


            //if(dca > 0.0) continue;


            for(Int_t i_param = 0; i_param < 9; i_param++)
            {
                aliHelix.fHelix[i_param] = AS_Track ->getHelix_param(i_param);
            }
            Double_t pathA = 0.0;
            Double_t helix_point[3];
            Double_t helix_pointB[3];

            //printf("i_track: %4.3d, i_event: %4.3d \n",i_track,i_event);


            make_clusters(i_track);
            //get_straight_line_fit(i_track);
            //get_tracklets_fit(i_track);

            //select_online_tracklets() uses offline circle fit to determine nearby tracklets
            get_2D_global_circle_fit();  // i will get dir_vec_circle (not needed?) and vec_TPL_circle_tracklets[i_point] (basically i_layer)
            select_online_tracklets();  // function to select onl tracklets close to something
            
            // for histos with corrected tracklets
            // Draw_corrected_online_tracklets();

            //vec_online_tracklet_points[i_layer][i_start_stop][i_XY] : starts and ends of online tracklets

            vector<TVector3> vec_TV3_tracklet_vectors;
            vec_TV3_tracklet_vectors.resize(6);


            //if(vec_tracklet_fit_points[6][0][0] == -999.0  && vec_tracklet_fit_points[6][1][0] == -999.0) continue; // no global fit available
            //printf("Track with pT: %4.3f used for calibration \n",pT_track);

            Int_t N_good_layers = 0;
            for(Int_t i_layer = 5; i_layer >= 0; i_layer--)
            {
                if(vec_online_tracklet_points[i_layer][0][0] > -999.0 && vec_online_tracklet_points[i_layer][1][0] > -999.0) N_good_layers++;
            }
            if(N_good_layers < 3) continue;

            if(tracklets_min_circle == 0.0 || tracklets_min_circle > 1.0) continue;


            //printf("N_good_layers: %d \n",N_good_layers);
            //printf("i_event: %lld, i_track: %d \n",i_event,i_track);


            Double_t impact_angle[6] = {0.0};
            Double_t impact_angle_circle[6] = {0.0};

            Double_t delta_x_local_global[6] = {0.0}; // local chamber coordinate system, global fit
            Double_t delta_x_local_global_circle[6] = {0.0}; // local chamber coordinate system, global circle fit - in case it's different



            for(Int_t i_layer = 5; i_layer >= 0; i_layer--)
            {
                if(vec_online_tracklet_points[i_layer][0][0] > -999.0 && vec_online_tracklet_points[i_layer][1][0] > -999.0)
                {

                    //---------------------------
                    // Calculate impact angle based on track helix    //maybe use this

                    Double_t path_offset = 0.0;
                    Double_t radius_layer = vec_TV3_TRD_center_offset[i_layer].Perp();
                    Double_t impact_angle_helix = -999.0;
                    for(Int_t i_step = 0; i_step < 400; i_step++)
                    {
                        pathA = i_step*3.0 + path_offset;
                        aliHelix.Evaluate(pathA,helix_point);
                        Double_t radius_helix_point = TMath::Sqrt(TMath::Power(helix_point[0],2) + TMath::Power(helix_point[1],2));

                        if(radius_helix_point > radius_layer)
                        {
                            aliHelix.Evaluate(pathA+0.5,helix_pointB);
                            TVector3 dir_helix_impact_on_layer(helix_pointB[0] - helix_point[0],helix_pointB[1] - helix_point[1],0.0);
                            Double_t delta_x_local_global_helix = dir_helix_impact_on_layer.Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);
                            Double_t sign_direction_impact_helix = TMath::Sign(1.0,delta_x_local_global_helix);
                            impact_angle_helix = dir_helix_impact_on_layer.Angle(vec_TV3_TRD_center[arr_layer_detector[i_layer]][2]);
                            if(impact_angle_helix > TMath::Pi()*0.5) impact_angle_helix -= TMath::Pi();
                            impact_angle_helix = 0.5*TMath::Pi() - sign_direction_impact_helix*impact_angle_helix;
                            //impact_angle[i_layer] = impact_angle_helix; // test to use TPC impact angle instead of TRD one

                            path_offset = pathA; // start from previous layer radius
                            break;
                        }
                    }
                    vec_h_diff_helix_line_impact_angle[0][arr_layer_detector[i_layer]] ->Fill(impact_angle[i_layer]*TMath::RadToDeg() - impact_angle_helix*TMath::RadToDeg());
                    if(dca < 0.0) vec_h_diff_helix_line_impact_angle[1][arr_layer_detector[i_layer]] ->Fill(impact_angle[i_layer]*TMath::RadToDeg() - impact_angle_helix*TMath::RadToDeg());
                    if(dca > 0.0) vec_h_diff_helix_line_impact_angle[2][arr_layer_detector[i_layer]] ->Fill(impact_angle[i_layer]*TMath::RadToDeg() - impact_angle_helix*TMath::RadToDeg());
                    //printf("impact_angle: %4.3f, impact_angle_helix: %4.3f \n",impact_angle[i_layer]*TMath::RadToDeg(),impact_angle_helix*TMath::RadToDeg());
                    //---------------------------

                    //back to angles calculation

                    // for histos with corrected tracklets
                    // corrected tracklets
                    // Double_t delta_x = vec_corrected_online_tracklet_points[i_layer][1][0] - vec_corrected_online_tracklet_points[i_layer][0][0];
                    // Double_t delta_y = vec_corrected_online_tracklet_points[i_layer][1][1] - vec_corrected_online_tracklet_points[i_layer][0][1];
                    
                    //online tracklets
                    Double_t delta_x = vec_online_tracklet_points[i_layer][1][0] - vec_online_tracklet_points[i_layer][0][0];
                    Double_t delta_y = vec_online_tracklet_points[i_layer][1][1] - vec_online_tracklet_points[i_layer][0][1];


                    //Double_t delta_z = vec_online_tracklet_points[i_layer][1][2] - vec_online_tracklet_points[i_layer][0][2];

                    //printf("i_layer: %d, delta_z: %4.3f \n",i_layer,delta_z);
                    TVector3 TV3_delta(delta_x,delta_y,0.0);   //do i need it for circle..


                    vec_TV3_tracklet_vectors[i_layer].SetXYZ(delta_x,delta_y,0.0);  //here we finally get vector: consists of xStop-xStart; yStop-yStart
                    //printf("i_layer: %d, delta_x: %4.3f, delta_y: %4.3f \n",i_layer,delta_x,delta_y);

                    if(vec_TV3_tracklet_vectors[i_layer].Mag() <= 0.0)
                    {
                        // printf("A: {%4.3f, %4.3f}, B: {%4.3f, %4.3f} \n",vec_tracklet_fit_points[i_layer][1][0],vec_tracklet_fit_points[i_layer][1][1],vec_tracklet_fit_points[i_layer][0][0],vec_tracklet_fit_points[i_layer][0][1]);
                        continue;
                    }

                    //vec_TV3_tracklet_vectors[i_layer] *= 1.0/vec_TV3_tracklet_vectors[i_layer].Mag(); // Normalize

                    Double_t length_vec = vec_TV3_tracklet_vectors[i_layer].Mag();
                    if(length_vec < 0.0)
                    {
                        vec_TV3_tracklet_vectors[i_layer] *= 1.0/length_vec;
                    }
                    if(length_vec > 0.0)
                    {
                        vec_TV3_tracklet_vectors[i_layer] *= -1.0/length_vec;
                    }

                    //printf("test 3 \n");

                    //if(tracklets_min[i_layer] > 0.2) continue;
                    if(arr_layer_detector[i_layer] < 0) continue;
                    Int_t detector = arr_layer_detector[i_layer];


                    //printf("test 3.1 \n");

                    //-----impact angle for circle global fit HERE---------

                    //printf("test 3.2 \n");

                    //printf("i_layer: %d \n",i_layer);

                    //printf("Global straight line dir: {%4.3f, %4.3f, %4.3f} \n",vec_TV3_tracklet_vectors[6].X(),vec_TV3_tracklet_vectors[6].Y(),vec_TV3_tracklet_vectors[6].Z());
                    //printf("Circle dir: {%4.3f, %4.3f, %4.3f} \n",vec_dir_vec_circle[i_layer].X(),vec_dir_vec_circle[i_layer].Y(),vec_dir_vec_circle[i_layer].Z());

                    delta_x_local_global_circle[i_layer] = vec_dir_vec_circle[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);
                    //Double_t proj_radial = TV3_delta.Dot(vec_TV3_TRD_center[arr_layer_detector[i_layerB]][2]);    looks like this isn't used

                    //printf("test 3.3 \n");

                    Double_t sign_direction_impact_circle = TMath::Sign(1.0,delta_x_local_global_circle[i_layer]);
                    impact_angle_circle[i_layer] = vec_dir_vec_circle[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layer]][2]);  //i need something like that for circles

                    //printf("test 3.4 \n");

                    if(impact_angle_circle[i_layer] > TMath::Pi()*0.5) impact_angle_circle[i_layer] -= TMath::Pi();
                    //printf("i_layerB: %d, 3D_Angle_on_TRD(TPC track): %4.3f, 2D_impact_angle: %4.3f \n",i_layerB,Angle_on_TRD*TMath::RadToDeg(),impact_angle[i_layerB]*TMath::RadToDeg());
                    impact_angle_circle[i_layer] = 0.5*TMath::Pi() - sign_direction_impact_circle*impact_angle_circle[i_layer];
                    //printf("impact_angle_circle: %4.3f \n",impact_angle_circle[i_layer]);


                    //-----------------------------------------------------

                    //printf("test 4 \n");


                    Double_t delta_x_local_tracklet = vec_TV3_tracklet_vectors[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);

                    //Double_t sign_angle        = 1.0;
                    Double_t sign_angle_circle = 1.0;
                    //if(delta_x_local_tracklet < delta_x_local_global[i_layer])        sign_angle        = -1.0;
                    if(delta_x_local_tracklet < delta_x_local_global_circle[i_layer]) sign_angle_circle = -1.0;

                    //printf("i_event: %d, i_track: %d \n",i_event,i_track);
                    //printf("delta_x_local_tracklet: %4.3f, delta_x_local_global_circle: %4.3f, sign_angle_circle %4.3f \n",delta_x_local_tracklet,delta_x_local_global_circle[i_layer],sign_angle_circle);


                    //calculate Delta angle for straight line
                    //Double_t Delta_angle  = sign_angle*vec_TV3_tracklet_vectors[6].Angle(vec_TV3_tracklet_vectors[i_layer]);
                    //if(Delta_angle > TMath::Pi()*0.5)  Delta_angle -= TMath::Pi();
                    //if(Delta_angle < -TMath::Pi()*0.5) Delta_angle += TMath::Pi();

                    //printf("Delta_angle: %4.3f \n",Delta_angle);


                    //same for delta angle tracklet vs circle
                    Double_t Delta_angle_circle  = sign_angle_circle*vec_dir_vec_circle[i_layer].Angle(vec_TV3_tracklet_vectors[i_layer]);
                    if(Delta_angle_circle > TMath::Pi()*0.5)  Delta_angle_circle -= TMath::Pi();
                    if(Delta_angle_circle < -TMath::Pi()*0.5) Delta_angle_circle += TMath::Pi();
                    //printf("Delta_angle_circle: %4.3f, impact_angle_circle[%d]: %4.3f \n",Delta_angle_circle*TMath::RadToDeg(),i_layer,impact_angle_circle[i_layer]*TMath::RadToDeg());



                    //printf("test 5 \n");
                    
                    // fill impact hists
                    h2D_delta_alpha ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(), Delta_angle_circle*TMath::RadToDeg(), 1);


                    h_detector_hit ->Fill(detector);
                    //if(detector == 69) printf("impact_angle: %4.3f, Delta_angle: %4.3f \n",impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());

                    //vec_tp_Delta_vs_impact[arr_layer_detector[i_layer]] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                    //vec_tp_Delta_vs_impact[detector] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),fabs(Delta_angle*TMath::RadToDeg()));

                    //fill histos for straight line global track
                    //vec_tp_Delta_vs_impact[detector]   ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                    //vec_TH2D_Delta_vs_impact[detector] ->Fill(impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());

                    //fill histos for circle global track
                    vec_tp_Delta_vs_impact_circle[detector]   ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());
                    vec_TH2D_Delta_vs_impact_circle[detector] ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());

                    //printf("test 6 \n");

#if 0
                    //separate hist for perpendicular impacts - straight line
                    if(impact_angle[i_layer]*TMath::RadToDeg() > 89.0 && impact_angle[i_layer]*TMath::RadToDeg() < 91.0)
                        //if(impact_angle[i_layer]*TMath::RadToDeg() > 80 && impact_angle[i_layer]*TMath::RadToDeg() < 83)
                    {
                        //    printf("detector: %d, track: %d, impact angle: %4.3f, Delta angle: %4.3f, delta_x_local_tracklet: %4.3f \n",detector,i_track,impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg(),delta_x_local_tracklet);
                        h_delta_angle_perp_impact ->Fill(Delta_angle*TMath::RadToDeg());
                    }
#endif

                    //separate hist for perpendicular impacts - cirrcle
                    if(impact_angle_circle[i_layer]*TMath::RadToDeg() > 89.0 && impact_angle_circle[i_layer]*TMath::RadToDeg() < 91.0)
                        //if(impact_angle[i_layer]*TMath::RadToDeg() > 80 && impact_angle[i_layer]*TMath::RadToDeg() < 83)
                    {
                        h_delta_angle_perp_impact_circle ->Fill(Delta_angle_circle*TMath::RadToDeg());
                    }


                    //if(impact_angle[i_layer]*TMath::RadToDeg() > 73.0 && impact_angle[i_layer]*TMath::RadToDeg() < 78.0 && pT_track > 3.5)
                    if(impact_angle[i_layer]*TMath::RadToDeg() > 87.5 && impact_angle[i_layer]*TMath::RadToDeg() < 92.5 && pT_track > 3.5)
                    {
                        //printf("       ======================> i_event: %lld, i_track: %d \n",i_event,i_track);
                    }
                    //printf("detector: %d, track: %d, impact angle: %4.3f, Delta angle: %4.3f \n",detector,i_track,impact_angle[i_layer]*TMath::RadToDeg(),Delta_angle*TMath::RadToDeg());
                
                }
            }
        }
    }

    // printf("test 7 \n");

#if 0
    // copy all that for circle angle histos
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
    h_dummy_Delta_vs_impact->GetYaxis()->SetRangeUser(-24,24);

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
#endif

    //printf("test 8 \n");

    //cancas for circle fit
    vector<TCanvas*> vec_can_Delta_vs_impact_circle;
    vec_can_Delta_vs_impact_circle.resize(6); // 6 sector blocks with 3 sectors each (18)

    TH1D* h_dummy_Delta_vs_impact_circle = new TH1D("h_dummy_Delta_vs_impact_circle","h_dummy_Delta_vs_impact_circle",90,50,140);
    h_dummy_Delta_vs_impact_circle->SetStats(0);
    h_dummy_Delta_vs_impact_circle->SetTitle("");
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact_circle->GetXaxis()->CenterTitle();
    h_dummy_Delta_vs_impact_circle->GetYaxis()->CenterTitle();
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitle("impact angle");
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitle("#Delta #alpha");
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetRangeUser(70,110);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetRangeUser(-24,24);

    Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    for(Int_t i_sec_block = 0; i_sec_block < 6; i_sec_block++)
    {
        HistName = "vec_can_Delta_vs_impact_circle_";
        HistName += i_sec_block;
        vec_can_Delta_vs_impact_circle[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_can_Delta_vs_impact_circle[i_sec_block] ->Divide(5,3); // x = stack, y = sector

        for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
        {
            Int_t i_sector = i_sec_block + 6*i_sec_sub;
            for(Int_t i_stack = 0; i_stack < 5; i_stack++)
            {
                Int_t iPad = i_sec_sub*5 + i_stack + 1;
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTicks(1,1);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetGrid(0,0);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetFillColor(10);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad);
                h_dummy_Delta_vs_impact_circle->Draw("h");

                for(Int_t i_layer = 0; i_layer < 6; i_layer++)
                {
                    Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                    // printf("detector: %d \n",i_detector);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineWidth(2);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineStyle(1);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->Draw("same hl");

                    HistName = "";
                    sprintf(NoP,"%4.0f",(Double_t)i_detector);
                    HistName += NoP;
                    plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                }
            }
        }
    }

#if 0
    TCanvas* can_delta_angle_perp_impact = new TCanvas("can_delta_angle_perp_impact","can_delta_angle_perp_impact",500,10,500,500);
    can_delta_angle_perp_impact ->cd();
    h_delta_angle_perp_impact  ->Draw();
#endif

    //same for circle perp impact
    TCanvas* can_delta_angle_perp_impact_circle = new TCanvas("can_delta_angle_perp_impact_circle","can_delta_angle_perp_impact_circle",500,10,500,500);
    can_delta_angle_perp_impact_circle ->cd();
    h_delta_angle_perp_impact_circle  ->Draw();

    TCanvas* can_detector_hit = new TCanvas("can_detector_hit","can_detector_hit",500,500,500,500);
    can_detector_hit ->cd();
    h_detector_hit  ->Draw();

    printf("Write data to output file \n");
    outputfile_trkl ->cd();
#if 0
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det]          ->Write();
        vec_TH2D_Delta_vs_impact[i_det]        ->Write();
        vec_tp_Delta_vs_impact_circle[i_det]   ->Write();
        vec_TH2D_Delta_vs_impact_circle[i_det] ->Write();
    }
#endif

    outputfile_trkl ->mkdir("Delta_impact_circle");
    outputfile_trkl ->cd("Delta_impact_circle");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact_circle[i_det]   ->Write();
        vec_TH2D_Delta_vs_impact_circle[i_det] ->Write();
    }

    outputfile_trkl ->cd();
    outputfile_trkl ->mkdir("Delta_impact");
    outputfile_trkl ->cd("Delta_impact");
    for(Int_t i_charge = 0; i_charge < 3; i_charge++)
    {
        for(Int_t i_det = 0; i_det < 547; i_det++)
        {
            vec_h_diff_helix_line_impact_angle[i_charge][i_det] ->Write();
            vec_h_diff_ref_loc_off_trkl[i_charge][i_det]        ->Write();
        }
    }
    printf("All data written \n");

}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_online_tracklets()
{
    printf("TBase_TRD_Calib::Draw_online_tracklets() \n");
    for(Int_t i_tracklet = 0; i_tracklet < (Int_t)vec_TPL3D_tracklets.size(); i_tracklet++)
    {
        printf("Add tracklet: %d \n",i_tracklet);
        vec_TPL3D_tracklets[i_tracklet]->SetLineStyle(1);
        vec_TPL3D_tracklets[i_tracklet]->SetLineWidth(2);
        vec_TPL3D_tracklets[i_tracklet]->SetMainColor(kOrange);
        gEve->AddElement(vec_TPL3D_tracklets[i_tracklet]);
    }
    gEve->Redraw3D(kTRUE);
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_corrected_online_tracklets()
{
    // printf("TBase_TRD_Calib::Draw_online_tracklets() \n");

    /* To plot delta alpha histo using the corrected tracklets, comment out two sets of line in calibrate_on_trkl()
    look for comment "for histos with corrected tracklets" to find these lines. Also disable drawing of
    corrected tracklets in this function because it will be called at every iteration */
    

    // TFile* calibration_params = TFile::Open("./TRD_Calib_vDfit_and_LAfit_online_trkl_new_withpre.root");
    TFile* calibration_params = TFile::Open("./Data/TRD_Calib_vDfit_and_LAfit_online_trkl_50k.root");
    // TRD_Calib_vDfit_and_LAfit_online_trkl_50k.root
    vdrift_fit = (TGraph*)calibration_params->Get(";1");

    vec_TPL_online_tracklets_selected_corrected.clear();
    vec_TPL_online_tracklets_selected_corrected.resize(6);

    Double_t impact_angle_circle[6] = {0.0};
    Double_t delta_x_local_global_circle[6] = {0.0}; 

    // create TPL for selected online tracklets and draw them
    for(Int_t i_layer = 0; i_layer < 6; i_layer++) // 6 layers of online tracklets
    {
        if(vec_online_tracklet_points[i_layer][0][0] > -999.0 && vec_online_tracklet_points[i_layer][1][0] > -999.0)
        {
            //extract LA and drift v parameters for detector from root file
            Int_t detector = draw_corrected_online_det_numbers[i_layer];
            Double_t det_LA = LA_fit->Eval(detector);           //Eval needs to be replaced
            Double_t det_vdrift = vdrift_fit->Eval(detector);

            bool is_defective = find(begin(Defect_TRD_detectors), end(Defect_TRD_detectors), detector) != end(Defect_TRD_detectors);
            if (is_defective){
                // cout << endl << "detector " << detector << " is broken - skipping" << endl;
                continue;
            }

            //calculate impact angle (used for debugging)
            delta_x_local_global_circle[i_layer] = vec_dir_vec_circle[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);
            Double_t sign_direction_impact_circle = TMath::Sign(1.0,delta_x_local_global_circle[i_layer]);
            impact_angle_circle[i_layer] = vec_dir_vec_circle[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layer]][2]);
            if(impact_angle_circle[i_layer] > TMath::Pi()*0.5) impact_angle_circle[i_layer] -= TMath::Pi();
            impact_angle_circle[i_layer] = 0.5*TMath::Pi() - sign_direction_impact_circle*impact_angle_circle[i_layer];
            
            
            //transform local impact angle into a global angle
            Int_t sector = detector/30;
            Double_t global_rotation = (90 -(sector*20 + 10))*TMath::DegToRad();
            // impact_angle_circle[i_layer] -= global_rotation;


            //transform impact to tracklet//
            // Double_t x_dir = TMath::Cos(impact_angle_circle[i_layer]);
            // Double_t y_dir = TMath::Sin(impact_angle_circle[i_layer]);
            // Double_t drift_vel_ratio = 1.546/det_vdrift;
            // det_LA -= global_rotation;


            //transform tracklet to impact//
            Double_t scale_length = 8.0;
            Double_t TRD_anode_plane = 3.35;
            if (det_vdrift == 0) {cout << "its zero" << endl; continue;};

            Double_t drift_vel_ratio = 1.546/det_vdrift;

            Double_t x_dir = -tv3_dirs_online[i_layer][0];
            Double_t y_dir = -tv3_dirs_online[i_layer][1];

            // cout << endl;

            // TVector3 center = vec_TV3_TRD_center[detector][0];
            // center[0] = -center[0];
            // center[1] = -center[1];

            // for (int i=0; i<3; i++){
            //     cout << center[i] << endl;
            // }

            // cout << "detector: " << detector << " | " << "sector number: " << sector << \
            //     " | " << "global rotation: " << global_rotation*TMath::RadToDeg() << endl;

            // cout << "LA = " << det_LA*TMath::RadToDeg() << " | " << "driftv = " << det_vdrift << endl;


            Double_t x_dir_local = x_dir*cos(global_rotation) - y_dir*sin(global_rotation);
            Double_t y_dir_local = x_dir*sin(global_rotation) + y_dir*cos(global_rotation);

            // cout << "tracklet angle: " << TMath::ATan(y_dir_local/x_dir_local)*TMath::RadToDeg() << endl;
            
            Double_t pre_corr = -0.16133;
            Double_t slope = 10000000.0;
            if(x_dir_local != 0.0) slope = y_dir_local/x_dir_local;

            if (y_dir_local < 0) {cout << "less than 0" << endl; y_dir_local = -1*y_dir_local;};

            Double_t Lorentz_tan   = TMath::Tan(det_LA);
            // Double_t Lorentz_slope = 10000000.0;
            // if(Lorentz_tan != 0.0) Lorentz_slope = 1.0/Lorentz_tan;

            Double_t x_anode_hit = TRD_anode_plane*drift_vel_ratio/slope;
            Double_t y_anode_hit = TRD_anode_plane*drift_vel_ratio;

            // cout << "(x,y) anode hit: " << x_anode_hit << "  " << y_anode_hit << endl;

            Double_t x_Lorentz_drift_hit = TMath::Tan(pre_corr)*TRD_anode_plane*drift_vel_ratio - \
                TMath::Tan(det_LA)*TRD_anode_plane;
            Double_t y_Lorentz_drift_hit = TRD_anode_plane*drift_vel_ratio - TRD_anode_plane;

            // cout << "(x,y) lorentz drift: " << x_Lorentz_drift_hit << "  " << y_Lorentz_drift_hit << endl;

            Double_t impact_angle_track = TMath::ATan2(y_anode_hit,x_anode_hit);

            Double_t Delta_x_Lorentz_drift_hit = x_anode_hit - x_Lorentz_drift_hit;
            Double_t Delta_y_Lorentz_drift_hit = y_anode_hit - y_Lorentz_drift_hit;
            Double_t impact_angle_rec = TMath::ATan2(Delta_y_Lorentz_drift_hit,Delta_x_Lorentz_drift_hit);

            // cout << "recon angle: " << impact_angle_rec*TMath::RadToDeg() << endl;


            // // back to global coords
            impact_angle_rec -= global_rotation;

            Double_t x_deflection = scale_length*TMath::Cos(impact_angle_rec);
            Double_t y_deflection = scale_length*TMath::Sin(impact_angle_rec);

            //apply correction
            vec_corrected_online_tracklet_points[i_layer][1][0] -= x_deflection;
            vec_corrected_online_tracklet_points[i_layer][1][1] -= y_deflection;

            //draw
            vec_TPL_online_tracklets_selected_corrected[i_layer] = new TPolyLine();
            vec_TPL_online_tracklets_selected_corrected[i_layer] ->SetNextPoint(vec_corrected_online_tracklet_points[i_layer][0][0],vec_corrected_online_tracklet_points[i_layer][0][1]);
            vec_TPL_online_tracklets_selected_corrected[i_layer] ->SetNextPoint(vec_corrected_online_tracklet_points[i_layer][1][0],vec_corrected_online_tracklet_points[i_layer][1][1]);
            vec_TPL_online_tracklets_selected_corrected[i_layer] ->SetLineStyle(9);
            vec_TPL_online_tracklets_selected_corrected[i_layer] ->SetLineColor(kMagenta);
            vec_TPL_online_tracklets_selected_corrected[i_layer] ->SetLineWidth(3);
            vec_TPL_online_tracklets_selected_corrected[i_layer] ->DrawClone("l");

            // printf("i_layer tracklets online: %d \n", i_layer);
            
        }
        else
        {
            // printf("TBase_TRD_Calib::Draw_tracklets_line(%d), i_layer: %d has no entry \n",i_track,i_layer);
        }
    }

    calibration_params->TFile::Close();
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::calc_impact_hist()
{

    TH2D* h2D_impact_hist = new TH2D("h2D_impact_hist","h2D_impact_hist",70,50,130,70,-1,1);
    TCanvas* can_impact_hist = new TCanvas("can_impact_hist","can_impact_hist"); 

    TProfile* profile_impact_hist = new TProfile("profile_impact_hist","profile_impact_hist",70,50,130);
    TCanvas* profile_impact_hist_can = new TCanvas("profile_impact_hist_can","profile_impact_hist_can");   
    

    for(Long64_t i_event = 0; i_event < 500; i_event++)
    {
        if(i_event % 20 == 0) printf("i_event: %lld out of %lld \n",i_event,file_entries_total);
        if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

        UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event

        for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        {
            AS_Track      = AS_Event ->getTrack( i_track ); // take the track
            TLorentzVector TLV_part = AS_Track ->get_TLV_part();
            Float_t pT_track        = TLV_part.Pt();
            // if(pT_track < 1.5) continue; // 3.5200

            make_clusters(i_track);


            get_2D_global_circle_fit();  // i will get dir_vec_circle (not needed?) and vec_TPL_circle_tracklets[i_point] (basically i_layer)
            select_online_tracklets();  // function to select onl tracklets close to something

            for(Int_t i_layer = 5; i_layer >= 0; i_layer--)
            {
                if(vec_online_tracklet_points[i_layer][0][0] > -999.0 && vec_online_tracklet_points[i_layer][1][0] > -999.0)
                {
                    Double_t impact_angle_circle[6] = {0.0};
                    Double_t delta_x_local_global_circle[6] = {0.0}; // local chamber coordinate system, global circle fit - in case it's different


                    delta_x_local_global_circle[i_layer] = vec_dir_vec_circle[i_layer].Dot(vec_TV3_TRD_center[arr_layer_detector[i_layer]][0]);
                    Double_t sign_direction_impact_circle = TMath::Sign(1.0,delta_x_local_global_circle[i_layer]);

                    impact_angle_circle[i_layer] = vec_dir_vec_circle[i_layer].Angle(vec_TV3_TRD_center[arr_layer_detector[i_layer]][2]);  //i need something like that for circles

                    if(impact_angle_circle[i_layer] > TMath::Pi()*0.5) impact_angle_circle[i_layer] -= TMath::Pi();

                    impact_angle_circle[i_layer] = 0.5*TMath::Pi() - sign_direction_impact_circle*impact_angle_circle[i_layer];
                    
                    // printf("impact_angle_circle: %4.3f \n",impact_angle_circle[i_layer]);
                    // cout << online_dy[i_layer] << endl;

                    h2D_impact_hist ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(), online_dy[i_layer], 1);
                    profile_impact_hist ->Fill(impact_angle_circle[i_layer]*TMath::RadToDeg(), online_dy[i_layer], 1);
                }
            }
            online_dy.clear();
        }
    }
    for (int i=0; i<40; i++)
    {
        cout << "Bin: " << i << endl;
        cout << profile_impact_hist->GetBinContent(i) << endl;
        cout << profile_impact_hist->GetBinEntries(i) << endl;
    }

    h2D_impact_hist->GetXaxis()->SetTitle("Impact Angle (deg)");
    h2D_impact_hist->GetYaxis()->SetTitle("dy (cm)");
    h2D_impact_hist->GetXaxis()->CenterTitle();
    h2D_impact_hist->GetYaxis()->CenterTitle();

    can_impact_hist->cd();
    h2D_impact_hist->Draw("colz");


    profile_impact_hist->GetXaxis()->SetTitle("Impact Angle (deg)");
    profile_impact_hist->GetYaxis()->SetTitle("Mean dy (cm)");
    profile_impact_hist->GetXaxis()->CenterTitle();
    profile_impact_hist->GetYaxis()->CenterTitle();

    profile_impact_hist_can->cd();
    profile_impact_hist->Draw();
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_offline_tracklets()
{
    printf("TBase_TRD_Calib::Draw_offline_tracklets() \n");
    for(Int_t i_tracklet = 0; i_tracklet < (Int_t)vec_TPL3D_offline_tracklets.size(); i_tracklet++)
    {
        printf("Add tracklet: %d \n",i_tracklet);
        vec_TPL3D_offline_tracklets[i_tracklet]->SetLineStyle(1);
        vec_TPL3D_offline_tracklets[i_tracklet]->SetLineWidth(2);
        vec_TPL3D_offline_tracklets[i_tracklet]->SetMainColor(kMagenta-7);
        gEve->AddElement(vec_TPL3D_offline_tracklets[i_tracklet]);
    }
    gEve->Redraw3D(kTRUE);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_2D_circle()
{
    // Calculate the circle parameters based on the input points
    // http://www.ambrsoft.com/trigocalc/circle3d.htm
    // (x - a)^2 + (y - b)^2  = R^2

    printf("TBase_TRD_Calib::Draw_2D_circle() \n");

    //Draw the fitted circle

    TPolyLine* tpl_circle = new TPolyLine();

    Double_t delta_i_y = par_circle[2]/20000.0;

    for(Double_t i_sign = -1.0; i_sign <= +1.0; i_sign += 2.0)
    {
        for(Double_t i_y = ((par_circle[1] - par_circle[2]) + delta_i_y); i_y < ((par_circle[1] + par_circle[2]) - delta_i_y); (i_y += delta_i_y))
        {
            Double_t i_x = i_sign*TMath::Sqrt(TMath::Power(par_circle[2],2.0) - TMath::Power(i_y - par_circle[1],2.0)) + par_circle[0];
            tpl_circle ->SetNextPoint(i_x,i_y);
            //printf("point: {%4.3f, %4.3f} \n",i_x,i_y);
        }
    }

    tpl_circle ->SetLineStyle(1);
    // tpl_circle ->SetLineWidth(3);
    tpl_circle ->SetLineWidth(2);
    //tpl_circle ->SetLineColor(kTeal+2);
    tpl_circle ->SetLineColor(kBlue);
    tpl_circle ->DrawClone("l");

    Int_t i_layer_notempty = 0;
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        if(vec_Dt_digit_pos_cluster[6][i_layer][0] != -999.0 && vec_Dt_digit_pos_cluster[6][i_layer][1] != -999.0)
        {
            i_layer_notempty++;
        }
    }

    Double_t delta_phi = -(vec_phi[1] - vec_phi[0])/1000.0;

    for(Int_t i_line = 0; i_line < (Int_t)vec_TPL_circle_tracklets.size(); i_line++)
    {
        delete vec_TPL_circle_tracklets[i_line];
    }
    vec_TPL_circle_tracklets.clear();

    vec_dir_vec_circle_notempty.resize(i_layer_notempty);

    for (Int_t i_point = 0; i_point < i_layer_notempty; ++i_point)
    {
        vec_TPL_circle_tracklets.push_back(new TPolyLine());

        Double_t phi_val = vec_phi[i_point];

        Double_t x_val   =  par_circle[0] + par_circle[2] * TMath::Cos(phi_val);
        Double_t y_val   =  par_circle[1] + par_circle[2] * TMath::Sin(phi_val);

        Double_t x_valB  =  par_circle[0] + par_circle[2] * TMath::Cos(phi_val + delta_phi);
        Double_t y_valB  =  par_circle[1] + par_circle[2] * TMath::Sin(phi_val + delta_phi);

        vec_dir_vec_circle_notempty[i_point].SetX(x_valB - x_val);
        vec_dir_vec_circle_notempty[i_point].SetY(y_valB - y_val);
        vec_dir_vec_circle_notempty[i_point].SetZ(0.0);

        if(vec_dir_vec_circle_notempty[i_point].Mag() > 0.0)
        {
            vec_dir_vec_circle_notempty[i_point] *= 6.0/vec_dir_vec_circle_notempty[i_point].Mag();

            vec_TPL_circle_tracklets[i_point] ->SetNextPoint(x_val,y_val);
            vec_TPL_circle_tracklets[i_point] ->SetNextPoint(x_val + vec_dir_vec_circle_notempty[i_point].X(),y_val + vec_dir_vec_circle_notempty[i_point].Y());
        }

        vec_TPL_circle_tracklets[i_point] ->SetLineColor(kCyan+1);
        // vec_TPL_circle_tracklets[i_point] ->SetLineWidth(5);
        vec_TPL_circle_tracklets[i_point] ->SetLineWidth(2);
        vec_TPL_circle_tracklets[i_point] ->Draw("");

    }
}


void TBase_TRD_Calib::Draw_2D_online_circle()
{
    // Calculate the circle parameters based on the input points
    // http://www.ambrsoft.com/trigocalc/circle3d.htm
    // (x - a)^2 + (y - b)^2  = R^2

    // Draw the fitted circle

    TPolyLine* tpl_circle = new TPolyLine();

    Double_t delta_i_y = par_circle_online[2]/20000.0;

    for(Double_t i_sign = -1.0; i_sign <= +1.0; i_sign += 2.0)
    {
        for(Double_t i_y = ((par_circle_online[1] - par_circle_online[2]) + delta_i_y); i_y < ((par_circle_online[1] + par_circle_online[2]) - delta_i_y); (i_y += delta_i_y))
        {
            Double_t i_x = i_sign*TMath::Sqrt(TMath::Power(par_circle_online[2],2.0) - TMath::Power(i_y - par_circle_online[1],2.0)) + par_circle_online[0];
            tpl_circle ->SetNextPoint(i_x,i_y);
            //printf("point: {%4.3f, %4.3f} \n",i_x,i_y);
        }
    }

    tpl_circle ->SetLineStyle(1);
    tpl_circle ->SetLineWidth(2);
    //tpl_circle ->SetLineColor(kTeal+2);
    tpl_circle ->SetLineColor(kGreen);
    tpl_circle ->DrawClone("l");

    //reset online circle fit points to remove possibility of being drawn falsely
    par_circle_online[0], par_circle_online[1], par_circle_online[2] = 0.0, 0.0, 0.0;

    // Int_t i_layer_notempty = 0;
    // for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    // {
    //     if(vec_Dt_digit_pos_cluster[6][i_layer][0] != -999.0 && vec_Dt_digit_pos_cluster[6][i_layer][1] != -999.0)
    //     {
    //         i_layer_notempty++;
    //     }
    // }

    //vec_phi.resize(i_layer_notempty);

    // Double_t delta_phi = -(vec_phi[1] - vec_phi[0])/1000.0;

    // for(Int_t i_line = 0; i_line < (Int_t)vec_TPL_circle_tracklets.size(); i_line++)
    // {
    //     delete vec_TPL_circle_tracklets[i_line];
    // }
    // vec_TPL_circle_tracklets.clear();

    // vec_dir_vec_circle_notempty.resize(i_layer_notempty);

    // for (Int_t i_point = 0; i_point < i_layer_notempty; ++i_point)
    // {
    //     vec_TPL_circle_tracklets.push_back(new TPolyLine());

    //     Double_t phi_val = vec_phi[i_point];

    //     Double_t x_val   =  par_circle_online[0] + par_circle_online[2] * TMath::Cos(phi_val);
    //     Double_t y_val   =  par_circle_online[1] + par_circle_online[2] * TMath::Sin(phi_val);

    //     Double_t x_valB  =  par_circle_online[0] + par_circle_online[2] * TMath::Cos(phi_val + delta_phi);
    //     Double_t y_valB  =  par_circle_online[1] + par_circle_online[2] * TMath::Sin(phi_val + delta_phi);

    //     vec_dir_vec_circle_notempty[i_point].SetX(x_valB - x_val);
    //     vec_dir_vec_circle_notempty[i_point].SetY(y_valB - y_val);
    //     vec_dir_vec_circle_notempty[i_point].SetZ(0.0);

    //     if(vec_dir_vec_circle_notempty[i_point].Mag() > 0.0)
    //     {
    //         vec_dir_vec_circle_notempty[i_point] *= 6.0/vec_dir_vec_circle_notempty[i_point].Mag();

    //         vec_TPL_circle_tracklets[i_point] ->SetNextPoint(x_val,y_val);
    //         vec_TPL_circle_tracklets[i_point] ->SetNextPoint(x_val + vec_dir_vec_circle_notempty[i_point].X(),y_val + vec_dir_vec_circle_notempty[i_point].Y());
    //     }


    //     vec_TPL_circle_tracklets[i_point] ->SetLineColor(kCyan+1);
    //     vec_TPL_circle_tracklets[i_point] ->SetLineWidth(5);
    //     vec_TPL_circle_tracklets[i_point] ->Draw("");

    // }
}

//----------------------------------------------------------------------------------------




#endif // __TBASE_TRD_CALIB_H__