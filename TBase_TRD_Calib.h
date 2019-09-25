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

ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Event)


//----------------------------------------------------------------------------------------
class TBase_TRD_Calib
{
private:
    Ali_AS_Event*     AS_Event;
    Ali_AS_Track*     AS_Track;
    Ali_AS_TRD_digit* AS_Digit;

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
    vector<Float_t> vec_digit_single_info;
    vector< vector<Float_t> > vec_digit_info;
    vector< vector< vector<Float_t> > > vec_digit_track_info;
    vector<Float_t> vec_track_single_info;
    vector< vector<Float_t> > vec_track_info;

    AliHelix aliHelix;
    TPolyLine3D* TPL3D_helix;


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
    void Draw_track(Int_t i_track);
    void Draw_TRD();

    ClassDef(TBase_TRD_Calib, 1)
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TBase_TRD_Calib::TBase_TRD_Calib()
{
    vec_digit_single_info.resize(11); // x,y,z,time,ADC,sector,stack,layer,row,column,dca
    vec_track_single_info.resize(12); // dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p

    TPM3D_single = new TPolyMarker3D();
    vec_TPM3D_digits.resize(6); // layers
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_digits[i_layer] = new TPolyMarker3D();
    }

    TPL3D_helix = new TPolyLine3D();

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
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_track(Int_t i_track)
{
    AS_Track      = AS_Event ->getTrack( i_track ); // take the track
    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
        aliHelix.fHelix[i_param] = AS_Track ->getHelix_param(i_param);
    }

    Double_t helix_point[3];
    Double_t pathA = 0.0;
    Double_t radius = 0.0;
    for(Int_t i_step = 0; i_step < 200; i_step++)
    {
        pathA = i_step*3.0;
        aliHelix.Evaluate(pathA,helix_point);

        radius = TMath::Sqrt(TMath::Power(helix_point[0],2.0) + TMath::Power(helix_point[1],2.0) + TMath::Power(helix_point[2],2.0));
        TPL3D_helix ->SetPoint(i_step,helix_point[0],helix_point[1],helix_point[2]);
        //printf("i_step: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_step,helix_point[0],helix_point[1],helix_point[2]);
        if(radius > 420.0) break;
    }

    TPL3D_helix    ->SetLineStyle(0);
    TPL3D_helix    ->SetLineColor(kBlue);
    TPL3D_helix    ->SetLineWidth(2);
    TPL3D_helix    ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Draw_TRD()
{
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
        //Double_t             loc[3]           = {TRD_time0 - 8.4,0.0,0.0};
        Double_t             glb[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,loc,glb);

        combitrans[TRD_detector] = new TGeoCombiTrans();
        combitrans[TRD_detector] ->RotateZ(Rotation_angle + 90.0);
        combitrans[TRD_detector] ->SetTranslation(glb[0],glb[1],glb[2]);
        top->AddNodeOverlap(TRD_boxes[TRD_detector],1,combitrans[TRD_detector]);

    }

    top->DrawClone("ogl");

#if 0
    x_BeamLine    ->DrawClone("ogl");
    y_BeamLine    ->DrawClone("ogl");
    z_BeamLine    ->DrawClone("ogl");
#endif

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TBase_TRD_Calib::Init_tree(TString SEList)
{
    cout << "Initialize tree" << endl;
    TString pinputdir = "/misc/alidata120/alice_u/schmah/TRD_offline_calib/Data/";

    AS_Event = new Ali_AS_Event();
    AS_Track = new Ali_AS_Track();
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

    //---------------------------------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event
    Double_t EventVertexX         = AS_Event ->getx();
    Double_t EventVertexY         = AS_Event ->gety();
    Double_t EventVertexZ         = AS_Event ->getz();
    Int_t    N_tracks_event       = AS_Event ->getN_tracks();
    Int_t    N_TRD_tracklets      = AS_Event ->getN_TRD_tracklets();
    Float_t  V0MEq                = AS_Event ->getcent_class_V0MEq();

    N_Tracks = NumTracks;

    vec_track_info.clear();

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

        printf("i_track: %d, fNumTRDdigits: %d \n",i_track,fNumTRDdigits);

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

                vec_digit_info.push_back(vec_digit_single_info);

                N_Digits ++;
            }

        } // end of digits loop

        vec_digit_track_info.push_back(vec_digit_info);

    } // end of track loop

    return 1;
}
//----------------------------------------------------------------------------------------



#endif // __TBASE_TRD_CALIB_H__