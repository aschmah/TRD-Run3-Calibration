#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

//------------------------
#include "AliHelix.h"
#include "TLorentzVector.h"
#include "TSystem.h"
//------------------------

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliKalmanTrack.h"

#include "AliTRDpadPlane.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"

#include "AliPIDResponse.h"
#include "AliPID.h"

#include "AliESDtrackCuts.h"

#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliESDRun.h"

#include "AliMultSelection.h"

#include "AliCDBEntry.h"
#include "TClonesArray.h"
#include "TGeoMatrix.h"
#include "AliAlignObjParams.h"

#include "AliTRDdigitsParam.h"
#include "AliRunTag.h"
#include "TObjString.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalOnlineGainTable.h"
#include "AliTRDCalROC.h"
#include "TPolyMarker.h"

#include "AliTRDCommonParam.h"

//------------------------
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
//------------------------

#include "Ali_AS_analysis_TRD_digits.h"


#include <iostream>
#include <iomanip>
using namespace std;

static Int_t flag_plot_event = 0;
static TString HistName;

static const Double_t mass_array[16]    = {0.0,0.0,0.00051099892,0.00051099892,0.0,0.105658369,0.105658369,0.1349766,0.13957018,0.13957018,0.497648,0.493677,0.493677,0.93956536,0.93827203,0.93827203}; // in GeV/c^2  {empty,gamma,e+,e-,neutrino,mu+,mu-,pi0,pi+,pi-,K0long,K+,K-,n,p,anti-p}
static TFile* dfile;
static TFile* TRD_alignment_file;
static TGeoHMatrix TM_TRD_rotation_det[540];
static TGeoHMatrix TM_TRD_rotation_sector[18];
static TVector3    TV3_TRD_translation[540];
static const Double_t TRD_lorentz_angle = TMath::DegToRad()*8.8;
//static const char *pathdatabase="alien://folder=/alice/data/2016/OCDB"; // for pPb
static const char *pathdatabase="alien://folder=/alice/data/2015/OCDB"; // for PbPb
static AliTRDCalDet *ChamberVdrift;
static AliTRDCalDet *ChamberT0;
static AliTRDCalDet *ChamberExB;
static AliTRDCalDet *chambergain;
static AliTRDCalOnlineGainTable *KryptoGain;
static AliTRDCalPad *PadNoise;
static AliTRDCalPad *LocalT0_pad;
static AliTRDCalROC* CalROC;
static AliTRDCalROC  *LocalT0;            //  Pad wise T0 calibration object
static const Int_t N_pT_bins = 5;
static const Double_t pT_ranges[N_pT_bins+1] = {0.2,0.5,1.0,2.0,3.0,5.0};
static const Double_t TRD_Impact_distance_in_drift = 3.35;

static AliTRDCommonParam* fParam;
//static Class_peak_finder my_class_peak_finder;
//static TH1D* h_ADC_vs_time;

ClassImp(Ali_AS_Event)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_analysis_TRD_digits)

    //________________________________________________________________________
    Ali_AS_analysis_TRD_digits::Ali_AS_analysis_TRD_digits(const char *name)
    : AliAnalysisTaskSE(name),
    fDigitsInputFileName("TRD.FltDigits.root"), fDigitsInputFile(0),
    fDigitsOutputFileName(""), fDigitsOutputFile(0),
    fDigMan(0),fGeo(0),AS_Event(0),AS_Track(0),AS_Digit(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
    fListOfHistos(0x0),fTree(0x0),h_dca(0x0),h_dca_xyz(0x0), h2D_TPC_dEdx_vs_momentum(0x0),h_ADC_tracklet(0x0), fPIDResponse(0), EsdTrackCuts(0)
{
    // Constructor

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());

    // Output slot #0 id reserved by the base class for AOD
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());

}

//_______________________________________________________________________
TFile* Ali_AS_analysis_TRD_digits::OpenDigitsFile(TString inputfile,
				   TString digfile,
				   TString opt)
{
    // we should check if we are reading ESDs or AODs - for now, only
    // ESDs are supported

    cout << "" << endl;
    cout << "In OpenDigitsFile" << endl;
    //cout << "Digits file name: " << digfile.Data() << endl;

    if(digfile == "")
    {
	cout << "WARNING: No TRD digits file available" << endl;
	return NULL;
    }

    // TGrid::Connect("alien")
    // construct the name of the digits file from the input file
    inputfile.ReplaceAll("AliESDs.root", digfile);
    //TString inputfile_LF = "alien:///alice/data/2016/LHC16q/000265525/pass1_CENT_wSDD/16000265525037.6203/TRD.FltDigits.root";
    TString inputfile_LF = inputfile;

    // open the file
    AliInfo( "opening digits file " + inputfile_LF + " with option \"" + opt + "\"");

    cout << "inputfile: " << inputfile_LF.Data() << endl;
    //TFile* dfile = new TFile(inputfile_LF, opt);
    //if(dfile) delete dfile;
    dfile = TFile::Open(inputfile_LF);
    cout << "After TRD digits file" << endl;
    cout << "" << endl;

    if(!dfile)
    {
	AliWarning("digits file '" + inputfile + "' cannot be opened");
    }

    return dfile;
}


//_______________________________________________________________________
Bool_t Ali_AS_analysis_TRD_digits::UserNotify()
{
    cout << "" << endl;
    cout << "In UserNotify" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;

    //h_ADC_vs_time = new TH1D("h_ADC_vs_time","h_ADC_vs_time",24,0,24);
    fParam = AliTRDCommonParam::Instance();

    delete fDigitsInputFile;
    delete fDigitsOutputFile;

    cout << "Digits file pointers deleted" << endl;

    cout << "All pointers deleted" << endl;

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
	(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if ( ! esdH ) return kFALSE;
    if ( ! esdH->GetTree() ) return kFALSE;
    if ( ! esdH->GetTree()->GetCurrentFile() ) return kFALSE;


    //-----------------------------------
    TList* list = esdH->GetUserInfo();
    //list->Print();
    //list->Dump();

    //cout << "List: " << endl;
    //list->ls();
    TList* list_cdblist = (TList*)list->FindObject("cdbList");  // list contains paths of calibration data
    Int_t list_size = list_cdblist->Capacity();
    //cout << "cdblist, with size: " << list_size <<  endl;
    //list_cdblist->ls();
    TObjString* obj_string = (TObjString*)list_cdblist->At(0);
    TString string = obj_string->GetString();
    //cout << "String: " << string.Data() << endl;
    Int_t index_A = string.Index("[");
    string.Remove(0,index_A+1);
    //cout << "String: " << string.Data() << endl;
    Int_t index_B = string.Index(",");
    string.Remove(index_B,string.Sizeof());
    Int_t run_number_from_list = string.Atoi();
    cout << "String: " << string.Data() << ", run_number_from_list: " << run_number_from_list << endl;
    //-----------------------------------


    // AliCDBEntry->GetObject()->IsA()->GetName()
    //-----------------------------------
    // Pad noise
    cout << "Open pad noise calibration file from database" << endl;
    AliCDBEntry *entryB = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/PadNoise",run_number_from_list);
    PadNoise = (AliTRDCalPad*)entryB->GetObject();
    cout << "Calibration data opened" << endl;
    //-----------------------------------


    //-----------------------------------
    // ChamberVdrift
    cout << "Open ChamberVdrift calibration file from database" << endl;
    AliCDBEntry *entryC = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/ChamberVdrift",run_number_from_list);
    ChamberVdrift = (AliTRDCalDet*)entryC->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", chamber vdrift: " << ChamberVdrift->GetValue(i_det) << endl;
    //}
    //-----------------------------------


    //-----------------------------------
    // ChamberT0
    cout << "Open ChamberT0 calibration file from database" << endl;
    AliCDBEntry *entryC1 = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/ChamberT0",run_number_from_list);
    ChamberT0 = (AliTRDCalDet*)entryC1->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", chamber t0: " << ChamberT0->GetValue(i_det) << endl; // in the order of "-1.36613"
    //}
    //-----------------------------------


    //-----------------------------------
    // LocalT0
    cout << "Open LocalT0 calibration file from database" << endl;
    AliCDBEntry *entryC2 = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/LocalT0",run_number_from_list);
    LocalT0_pad = (AliTRDCalPad*)entryC2->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    Int_t col = 2;
    //    Int_t row = 2;
    //    cout << "i_det: " << i_det << ", local t0: " << LocalT0_pad->GetCalROC(i_det)->GetValue(col,row) << endl; // returns 0
    //}
    //-----------------------------------


#if 0
    //-----------------------------------
    // LocalVdrift
    // Values are all 0
    cout << "Open LocalVdrift calibration file from database" << endl;
    AliCDBEntry *entryD = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/LocalVdrift",run_number_from_list);
    AliTRDCalDet *LocalVdrift = (AliTRDCalDet*)entryD->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", local vdrift: " << LocalVdrift->GetValue(i_det) << endl;
    //}
    //-----------------------------------
#endif


    //-----------------------------------
    // ChamberExB
    // Values are all 0
    cout << "Open ChamberExB calibration file from database" << endl;
    AliCDBEntry *entryE = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/ChamberExB",run_number_from_list);
    ChamberExB = (AliTRDCalDet*)entryE->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", ExB: " << ChamberExB->GetValue(i_det) << endl;
    //}
    //-----------------------------------



    //----------------------------------------------------------------------
    // Krypton calibration -> pad gain factors
    cout << "Open Krypto calibration file from database" << endl;
    AliCDBEntry *entryF = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/Krypton_2015-02",run_number_from_list);
    KryptoGain = (AliTRDCalOnlineGainTable*)entryF->GetObject();
    cout << "Calibration data opened" << endl;
    //Float_t GainFactor = KryptoGain ->GetGainCorrectionFactor(i_det,i_row,i_column);
    //----------------------------------------------------------------------



    //----------------------------------------------------------------------
    // Chamber gain factors
    cout << "Open chamber gain file from database" << endl;
    AliCDBEntry *entryG = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/ChamberGainFactor",run_number_from_list);
    chambergain = (AliTRDCalDet*)entryG->GetObject();
    //Float_t valuegainiteration = chambergain->GetValue(det);
    //----------------------------------------------------------------------



    //TObjString* ts = (TObjString*)list->FindObject("TRD/Calib/ChamberGainFactor");
    //TString string = ts->GetString();
    //cout << "string: " << string.Data() << endl;

    TString fname = esdH->GetTree()->GetCurrentFile()->GetName();
    TString Tree_name = esdH->GetTree()->GetName();
    FileStat_t file_stat;
    Int_t PathInfo = gSystem->GetPathInfo(fname.Data(),file_stat);
    cout << "PathInfo: " << PathInfo << ", fname: " << fname << endl;
    //TFile* file = TFile::Open(fname.Data());
    //cout << "Zombie: " << file->IsZombie() << ", header size: " << file->Sizeof() << ", FileBytesRead: " << file->GetFileBytesRead() << endl;

    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!inputHandler)
    {
	printf("WARNING: Inputhandler not available \n");
    }
    else
    {
	printf("Inputhandler available \n");

	fPIDResponse = inputHandler->GetPIDResponse();

        cout << "Got PID response" << endl;
    }

    //AliAnalysisManager *manB = AliAnalysisManager::GetAnalysisManager();
    //Int_t run_id = manB->GetRunFromPath();
    //cout << "fname: " << fname << ", tree name: " << Tree_name << ", run_id: " << run_id << endl;


    fEventNoInFile = -1;
    N_good_events  = 0;


    fDigitsInputFile = OpenDigitsFile(fname,fDigitsInputFileName,""); // <-
    EsdTrackCuts = new AliESDtrackCuts();

    EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
    EsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.52);
    EsdTrackCuts->AliESDtrackCuts::SetMinNClustersTPC(50); // 60, Automatically requires TPC refitted tracks?
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXY(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetPtRange(0.15,200.0); // 0.15, 200.0
    EsdTrackCuts->AliESDtrackCuts::SetEtaRange(-1.0,1.0); // 0.85

    // create the digits manager
    cout << "" << endl;
    cout << "________________________________________________________________________" << endl;
    cout << "Created AliTRDdigitsManager" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;
    cout << "" << endl;

    // create a TRD geometry, needed for matching digits to tracks
    fGeo = new AliTRDgeometry;
    if(!fGeo)
    {
	AliFatal("cannot create geometry ");
    }

    //if(fDigMan) delete fDigMan;
    fDigMan = new AliTRDdigitsManager;
    fDigMan->CreateArrays();

    for(Int_t i_det = 0; i_det < 5; i_det++)
    {
	Int_t N_columns   = fDigMan->GetDigits(i_det)->GetNcol();
	cout << "i_det: " << i_det << ", N_columns: " << N_columns << endl;
    }

    if(fname.Contains("/home/"))
    {
        TRD_alignment_file = TFile::Open("/home/ceres/schmah/ALICE/Database/TRD_Align_2016.root");
	cout << "Local alignment file loaded" << endl;
    }
    else
    {
	TRD_alignment_file = TFile::Open("alien:///alice/data/2016/OCDB/TRD/Align/Data/Run0_999999999_v1_s0.root");
        cout << "Alignment file from database loaded" << endl;
    }

    AliCDBEntry* align_cdb = (AliCDBEntry*)TRD_alignment_file->Get("AliCDBEntry");
    TClonesArray* TRD_align_array = (TClonesArray*)align_cdb->GetObject();
    // Alignment is stored for every INSTALLED chamber. Array does NOT got to 540!
    // First do the alignment per sector, then per chamber
    for(Int_t i_entry = 0; i_entry < TRD_align_array->GetEntries(); i_entry++)
    {
	AliAlignObjParams* Align_test = (AliAlignObjParams*)TRD_align_array->At(i_entry);
	char* name = (char*)Align_test->GetSymName();
	TString sname(name);
        Int_t length = sname.Length();
	TString sector_name = sname(6,2);
	TString stack_name  = sname(11,1);
	TString layer_name  = sname(15,1);
	Int_t sector = sector_name.Atoi();
	Int_t stack  = stack_name.Atoi();
	Int_t layer  = layer_name.Atoi();
        Int_t detector = layer + stack * 6 + sector * 6 * 5;

	TGeoHMatrix TM_TRD_align;
        Align_test->GetMatrix(TM_TRD_align);

	//printf("i_entry: %d, name: %s, length: %d, sector: %s, stack: %s, layer: %s, detector: %d \n",i_entry,name,sname.Length(),sector_name.Data(),stack_name.Data(),layer_name.Data(),detector);
        //TM_TRD_align.Print();
	if(length > 8)
	{
	    Align_test->GetMatrix(TM_TRD_rotation_det[detector]);
	    Double_t* rot_matrix = TM_TRD_rotation_det[detector].GetRotationMatrix();
	    //for(Int_t i = 0; i < 9; i++)
	    //{
	    //    cout << "i: " << i << ", rot_matrix: " << rot_matrix[i] << endl;
	    //}
	    Double_t TRD_translation[3];
	    Align_test->GetTranslation(TRD_translation);
	    TV3_TRD_translation[detector].SetXYZ(TRD_translation[0],TRD_translation[1],TRD_translation[2]);
	    //cout << "i_entry: " << i_entry << ", trans = {" << TV3_TRD_translation[i_entry].X() << ", " << TV3_TRD_translation[i_entry].Y() << ", " << TV3_TRD_translation[i_entry].Z() << "}" << endl;
	    //TM_TRD_rotation_det[i_entry].Print();
	}
	else
	{
	    Align_test->GetMatrix(TM_TRD_rotation_sector[sector]);
	}
    }
    cout << "TRD_align_array entries: " << TRD_align_array->GetEntries() << endl;



    cout << "End of UserNotify" << endl;
    return kTRUE;
}


//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::UserCreateOutputObjects()
{
    cout << "" << endl;
    cout << "In UserCreateOutputObjects" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;


    OpenFile(1);
    cout << "File opened" << endl;

    fListOfHistos = new TList();
    fListOfHistos ->SetOwner();

    h_dca.resize(N_pT_bins);
    for(Int_t i_pT = 0; i_pT < N_pT_bins; i_pT++)
    {
	HistName = "h_dca_pT_";
        HistName += i_pT;
	h_dca[i_pT] = new TH1D(HistName.Data(),HistName.Data(),1000,0.0,70.0);
	h_dca[i_pT] ->GetXaxis()->SetTitle("dca");
	fListOfHistos->Add(h_dca[i_pT]);
    }
    h_ADC_tracklet.resize(2);
    for(Int_t i_ADC = 0; i_ADC < 2; i_ADC++)
    {
        HistName = "h_ADC_tracklet_";
        HistName += i_ADC;
        h_ADC_tracklet[i_ADC] = new TH1D(HistName.Data(),HistName.Data(),350,-50.0,300.0);
        fListOfHistos->Add(h_ADC_tracklet[i_ADC]);
    }

    h_dca_xyz.resize(6); // x, y, z
    for(Int_t i_xyz = 0; i_xyz < 6; i_xyz++)
    {
	h_dca_xyz[i_xyz].resize(N_pT_bins);
	for(Int_t i_pT = 0; i_pT < N_pT_bins; i_pT++)
	{
	    HistName = "h_dca_xyz_";
	    HistName += i_xyz;
            HistName += "_pT_";
	    HistName += i_pT;
	    h_dca_xyz[i_xyz][i_pT] = new TH1D(HistName.Data(),HistName.Data(),1000,-70.0,70.0);
	    h_dca_xyz[i_xyz][i_pT] ->GetXaxis()->SetTitle("dca");
	    fListOfHistos->Add(h_dca_xyz[i_xyz][i_pT]);
	}
    }

    HistName = "h2D_TPC_dEdx_vs_momentum";
    h2D_TPC_dEdx_vs_momentum = new TH2D(HistName.Data(),HistName.Data(),2000,-30.0,30.0,2000,0,1000);
    fListOfHistos->Add(h2D_TPC_dEdx_vs_momentum);

    OpenFile(2);
    cout << "File opened" << endl;

    AS_Event       = new Ali_AS_Event();
    AS_Track       = new Ali_AS_Track();
    AS_Digit       = new Ali_AS_TRD_digit();
    Tree_AS_Event  = NULL;
    Tree_AS_Event  = new TTree("Tree_AS_Event" , "AS_Events" );
    Tree_AS_Event  ->Branch("Tree_AS_Event_branch"  , "AS_Event", AS_Event );

    PostData(1,fListOfHistos);
    PostData(2,Tree_AS_Event);

    cout << "PostData called" << endl;

}



//________________________________________________________________________
Bool_t Ali_AS_analysis_TRD_digits::NextEvent(Bool_t preload)
{
    fEventNoInFile++;
    //cout << "fEventNoInFile: " << fEventNoInFile << endl;
    fDigitsLoadedFlag = kFALSE;


    if(preload)
    {
	//cout << "Preload"  << endl;
	return ReadDigits();
    }
    else
    {
	//cout << "No preload"  << endl;
	return kTRUE;
    }
}

//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::UserExec(Option_t *)
{
    //cout << "" << endl;
    //cout << "Analysis started" << endl;
    //cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;


    //-----------------------------------------------------------------
    // IMPORTANT: call NextEvent() for book-keeping
    NextEvent();
    if(fEventNoInFile > 5) return;
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    // prepare event data structures
    AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD)
    {
	printf("ERROR: fESD not available\n");
	return;
    }

    //-----------------------------------------------------------------
    // Check if TRD digits (raw data) are available for this ESD event
    if(!ReadDigits()) return;
    //-----------------------------------------------------------------


    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man)
    {
        //Int_t run_id = man->GetRunFromPath(); // doesn't work
	//cout << "Got AliAnalysisManager, run_id: " << run_id << endl;
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
	if(inputHandler)
	{
	    //cout << "Got AliInputEventHandler" << endl;
	    fPIDResponse = inputHandler->GetPIDResponse();
	}
    }
    //cout << "cent: " << fPIDResponse->GetCurrentCentrality() << endl;


    Int_t          N_tracks         = fESD ->GetNumberOfTracks();
    Int_t          N_TRD_tracks     = fESD ->GetNumberOfTrdTracks();
    Int_t          N_TRD_tracklets  = fESD ->GetNumberOfTrdTracklets();
    Float_t        magF             = fESD ->GetMagneticField();
    const AliESDVertex* PrimVertex  = fESD ->GetPrimaryVertex();
    Int_t          RunNum           = fESD ->GetRunNumber();
    Double_t       T0zVertex        = fESD ->GetT0zVertex();
    AliCentrality* Centrality       = fESD ->GetCentrality();
    Double_t       MeanBeamIntAA    = fESD ->GetESDRun()->GetMeanIntensity(0,0);

    //printf("RunNum: %d \n",RunNum);

    Double_t Sign_magnetic_field = (magF/fabs(magF));
    //cout << "Trigger: " <<  fESD->GetFiredTriggerClasses() << endl;

    // Fill event information
    AS_Event ->clearTrackList();
    AS_Event ->setTriggerWord(fESD->GetFiredTriggerClasses());
    AS_Event ->setx(PrimVertex->GetX());
    AS_Event ->sety(PrimVertex->GetY());
    AS_Event ->setz(PrimVertex->GetZ());
    AS_Event ->setid(RunNum);
    AS_Event ->setN_tracks(N_tracks);
    AS_Event ->setN_TRD_tracklets(N_TRD_tracklets);
    AS_Event ->setBeamIntAA(MeanBeamIntAA);
    AS_Event ->setT0zVertex(T0zVertex);

    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    if(MultSelection)
    {
	// V0MEq, V0AEq, V0CEq, SPDTracklets

	AS_Event ->setcent_class_ZNA(MultSelection->GetMultiplicityPercentile("ZNA"));
	AS_Event ->setcent_class_ZNC(MultSelection->GetMultiplicityPercentile("ZNC"));
	AS_Event ->setcent_class_V0A(MultSelection->GetMultiplicityPercentile("V0A"));
	AS_Event ->setcent_class_V0C(MultSelection->GetMultiplicityPercentile("V0C"));
	AS_Event ->setcent_class_V0M(MultSelection->GetMultiplicityPercentile("V0M"));
	AS_Event ->setcent_class_CL0(MultSelection->GetMultiplicityPercentile("CL0"));
	AS_Event ->setcent_class_CL1(MultSelection->GetMultiplicityPercentile("CL1"));
	AS_Event ->setcent_class_SPD(MultSelection->GetMultiplicityPercentile("SPDTracklets"));
	AS_Event ->setcent_class_V0MEq(MultSelection->GetMultiplicityPercentile("V0MEq"));
	AS_Event ->setcent_class_V0AEq(MultSelection->GetMultiplicityPercentile("V0AEq"));
	AS_Event ->setcent_class_V0CEq(MultSelection->GetMultiplicityPercentile("V0CEq"));
    }


    //-----------------------------------------------------------------
    //cout << "" << endl;
    //cout << "" << endl;
    //cout << "----------------------------------------------------------------------------------------" << endl;
    //cout << "Event number: " << fEventNoInFile << ", event number with TRD digits: " << N_good_events << endl;
    //printf("Event number: %d, N_tracks: %d, cent(V0M): %f , cent(CL0): %f \n",fEventNoInFile,N_tracks,MultSelection->GetMultiplicityPercentile("SPDTracklets"));
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------
    Double_t TRD_time_per_bin        = 0.1;   // 100 ns = 0.1 mus
    Double_t TRD_drift_velocity      = 1.56;  // 1.56 cm/mus, this value is too high for p+Pb, database values are now used
    Double_t TPC_TRD_matching_window = 10.0;  // Matching window between TRD digits (pad position) and TPC track in cm
    Double_t TPC_min_radius_plot     = 114.0/2.0; // Minimum radius for TPC track to be plotted
    Double_t TPC_max_radius_plot     = 368.0; // Maximum radius for TPC track to be plotted
    Double_t TPC_radius_scan         = 368.0 - 30.0; // Radius for TPC track to be used for first match with TRD digits

    if(N_tracks == 0)
    {
	// Skip empty event
	return;
    }

    //printf("There are %d tracks in this event\n", N_tracks);
    //printf("There are %d TRD tracks in this event\n", N_TRD_tracks);
    //printf("There are %d TRD tracklets in this event\n", N_TRD_tracklets);


    Int_t N_good_tracks = 0;


    //-----------------------------------------------------------------
    // Loop over all TRD channels
    std::vector<TVector3> TV3_TRD_hits;
    std::vector<TVector3> TV3_TRD_hits_det_angle;
    std::vector<Double_t> vec_TRD_hits_ADC_value;
    std::vector<Double_t> vec_TRD_hits_time;
    std::vector< std::vector<Int_t> >  vec_TRD_det_lay_row_col;

    std::vector<TVector3> TV3_TRD_hits_middle;
    std::vector<TVector3> TV3_TRD_hits_det_angle_middle;
    std::vector< std::vector<Double_t> > vec_TRD_hits_ADC_values_time;
    std::vector< std::vector<TVector3> > vec_TRD_hits_points_time;
   // std::vector< std::vector<TVector3> > vec_TRD_hits_points_time_uncalib;
    std::vector< std::vector<Int_t> >    vec_TRD_hits_det_lay_row_col;

    vec_TRD_det_lay_row_col.resize(4);

    Int_t max_N_rows    = 0;
    Int_t max_N_columns = 0;

    Double_t sum_full_ADC_digit_det[540];
    memset(sum_full_ADC_digit_det, 0, sizeof(sum_full_ADC_digit_det)); // for automatically-allocated arrays

    cout << "Loop over all TRD detectors" << endl;
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
	Int_t N_columns   = fDigMan->GetDigits(i_det)->GetNcol();
	Int_t N_rows      = fDigMan->GetDigits(i_det)->GetNrow();
	Int_t N_times     = fDigMan->GetDigits(i_det)->GetNtime();
	Int_t N_detectors = fDigMan->GetDigits(i_det)->GetNdet();

	//AliTRDarrayADC* trd_array_adc = fDigMan->GetDigits(i_det);
	//trd_array_adc->CreateLut();

	//cout << "TRD detector: " << i_det << endl;

	Int_t i_sector = fGeo->GetSector(i_det);
	Int_t i_stack  = fGeo->GetStack(i_det);
	Int_t i_layer  = fGeo->GetLayer(i_det);

	//printf("i_det: %d, N_columns: %d, N_rows: %d, N_times: %d, i_sector: %d, i_stack: %d, i_layer: %d \n",i_det,N_columns,N_rows,N_times,i_sector,i_stack,i_layer);

	if(N_columns > max_N_columns) max_N_columns = N_columns;
	if(N_rows    > max_N_rows)    max_N_rows    = N_rows;

	for(Int_t i_column = 0; i_column < N_columns; i_column++)
	{
	    for(Int_t i_row = 0; i_row < N_rows; i_row++)
	    {
		Double_t ADC_amplitude_sum_times = 0.0;
		std::vector<Double_t> vec_ADC_time_bins;
		vec_ADC_time_bins.resize(N_times);
		std::vector<TVector3> vec_points_time_bins;
		vec_points_time_bins.resize(N_times);
		//cout << "i_column: " << i_column << ", i_row: " << i_row << ", N_times: " << N_times << endl;
		Double_t arr[24] = {0.0};
		for(Int_t i_time = 0; i_time < N_times; i_time++)
		{
		    Double_t ADC_amplitude = fDigMan->GetDigitAmp(i_row,i_column,i_time,i_det);
		    Int_t fBaseline = fDigMan->GetDigitsParam()->GetADCbaseline(i_det); // ADC baseline to be subtracted from digit ADC values, constant value of 10
		    ADC_amplitude_sum_times += ADC_amplitude;

                    //printf("ADC_amplitude: %4.3f \n",ADC_amplitude);

		    if(ADC_amplitude > 0.0)
		    {
                        sum_full_ADC_digit_det[i_det] += ADC_amplitude - fBaseline;

                        //printf("i_time: %d, ADC: %f \n",i_time,ADC_amplitude);
                        //h_ADC_vs_time ->SetBinContent(i_time+1,ADC_amplitude);
			arr[i_time] = ADC_amplitude;
			//cout << "i_time: " << i_time << ", fill arr: " << arr[i_time] << endl;
                        //cout << "fBaseline: " << fBaseline << ", ADC: " << ADC_amplitude << endl;
			Int_t x_TRD = i_column + i_sector*144;
			Int_t y_TRD = i_row    + i_stack*16 + i_layer*16*5;

			//----------------------
			// Calculate global position for fired TRD pad
			Float_t              TRD_time0        = fGeo         ->GetTime0(i_layer); // in cm
			AliTRDpadPlane*      padplane         = fGeo         ->GetPadPlane(i_det);
			Double_t             TRD_col_end      = padplane     ->GetColEnd();
			Double_t             TRD_row_end      = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
			Double_t             TRD_row_start    = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
			Double_t             TRD_row_end_ROC  = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
			Double_t             TRD_col_spacing  = padplane     ->GetColSpacing();
			Double_t             TRD_row_spacing  = padplane     ->GetRowSpacing();
			Double_t             TRD_col_pos      = padplane     ->GetColPos(i_column);
			Double_t             TRD_row_pos      = padplane     ->GetRowPos(i_row);       // fPadRow[row] + fPadRowSMOffset;
			Double_t             TRD_row_pos_ROC  = padplane     ->GetRowPosROC(i_row);    // fPadRow[row]; = pad border position
		        Double_t             TRD_col_size     = padplane     ->GetColSize(i_column);
			Double_t             TRD_row_size     = padplane     ->GetRowSize(i_row);
			Float_t              TRD_loc_Y        = TRD_col_pos + 0.5*TRD_col_size;
			Double_t             TRD_drift_time   = ((Double_t)i_time)*TRD_time_per_bin*ChamberVdrift->GetValue(i_det); // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum

                        //printf("TRD_tim0: %4.3f, vdrift: %4.3f \n",TRD_time0,ChamberVdrift->GetValue(i_det));

			Float_t lorentz_angle_corr_y = ChamberExB->GetValue(i_det)*TRD_drift_time;
			TRD_loc_Y -= Sign_magnetic_field*lorentz_angle_corr_y;
			//cout << "i_time: " << i_time << ", magF: " << magF << ", x: " << TRD_drift_time << ", lorentz_angle_corr_y: " << lorentz_angle_corr_y << endl;

			//Float_t              TRD_loc_Z        = (Float_t)i_row*8.5+8.5/2.0;
			//Float_t              TRD_loc_Y        = (Float_t)i_column*0.725 + 0.725/2.0 - (144.0*0.725)/2.0;
			//Double_t             loc[3]           = {TRD_time0,TRD_loc_Y,TRD_loc_Z+TRD_row_end};

			// Determine global position of TRD hit
			Double_t             loc[3]           = {TRD_time0 - TRD_drift_time,TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};
			Double_t             glb[3]           = {0.0,0.0,0.0};
			fGeo ->RotateBack(i_det,loc,glb);

			// Apply alignment
			Double_t glb_align_sec[3];
			Double_t glb_align[3];

			TM_TRD_rotation_sector[i_sector].LocalToMaster(glb,glb_align_sec);
			TM_TRD_rotation_det[i_det].LocalToMaster(glb_align_sec,glb_align);

			//glb_align[0] = glb[0];
                        //glb_align[1] = glb[1];
                        //glb_align[2] = glb[2];

			// Determine global pointing vector of TRD plane (vector perpendicular to plane in global system)
			Double_t             loc_vec[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
			Double_t             glb_vec[3]           = {0.0,0.0,0.0};
			fGeo ->RotateBack(i_det,loc_vec,glb_vec);


			// Apply alignment
			Double_t glb_vec_align_sec[3];
			Double_t glb_vec_align[3];
			TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_vec,glb_vec_align_sec);
			TM_TRD_rotation_det[i_det].LocalToMaster(glb_vec_align_sec,glb_vec_align);

			//for(Int_t i = 0; i < 3; i++)
			//{
			//    glb_align[i] = glb[i];
                        //    glb_vec_align[i] = glb_vec[i];
			//}


			Double_t             TRD_radius   = TMath::Sqrt(glb_align[0]*glb_align[0] + glb_align[1]*glb_align[1]);
#if 0
			if(i_det >= 0 && i_time == 0 && i_column == 0 && i_row >= 0)
			{
			    cout << "i_det: " << i_det << ", i_row: " << i_row << ", N_rows: " << N_rows << ", row_end: " << TRD_row_end
				<< ", row_end_ROC: " << TRD_row_end_ROC << ", row_pos: " << TRD_row_pos
				<< ", row_pos_ROC: " << TRD_row_pos_ROC
				<< ", row_size: " << TRD_row_size << endl;
			}
#endif

#if 0
			cout << "ADC_amplitude: " << ADC_amplitude << ", i_row: " << i_row
			    << ", i_column: " << i_column << ", i_time: " << i_time << ", detector: " << i_det
			    << ", col_end: " << TRD_col_end << ", row_end: " << TRD_row_end
			    << ", col_spacing: " << TRD_col_spacing << ", row_spacing: " << TRD_row_spacing
			    << ", col_pos: " << TRD_col_pos << ", row_pos: " << TRD_row_pos
			    << ", col_size: " << TRD_col_size << ", row_size: " << TRD_row_size << endl;
#endif

			if(ADC_amplitude > 0.0)
			{
			    TVector3 TV3_TRD_hit, TV3_TRD_hit_det_angle;
			    TV3_TRD_hit.SetXYZ(glb_align[0],glb_align[1],glb_align[2]);
			    TV3_TRD_hit_det_angle.SetXYZ(glb_vec_align[0],glb_vec_align[1],0.0); // glb_vec is from {0,0,0} to center of detector, set z component to 0 to get vector perpendicular to detector plane
			    TV3_TRD_hits.push_back(TV3_TRD_hit);
			    TV3_TRD_hits_det_angle.push_back(TV3_TRD_hit_det_angle);
			    vec_TRD_hits_ADC_value.push_back(ADC_amplitude);
			    vec_TRD_hits_time.push_back(((Double_t)i_time)*0.1); // time in mus
			    vec_TRD_det_lay_row_col[0].push_back(i_det);
			    vec_TRD_det_lay_row_col[1].push_back(i_layer);
			    vec_TRD_det_lay_row_col[2].push_back(i_row);
			    vec_TRD_det_lay_row_col[3].push_back(i_column);

			    vec_ADC_time_bins[i_time]    = ADC_amplitude;
			    vec_points_time_bins[i_time] = TV3_TRD_hit;
			}
			//----------------------
		    }
		} // end of time loop
		//--------------------------------------------------



                //--------------------------------------------------
		Int_t    max_time_bin = -1;
                Double_t max_ADC_val  = 0.0;
		if(ADC_amplitude_sum_times > 0.0)
		{
		    for(Int_t i_time = 0; i_time < N_times; i_time++)
		    {
			if(arr[i_time] > max_ADC_val)
			{
			    max_ADC_val  = arr[i_time];
                            max_time_bin = i_time;
			}
		    }
		    //printf("   ->peak at time bin: %d, ADC: %f \n",max_time_bin,max_ADC_val);
		    //my_class_peak_finder.add_TH1D(h_ADC_vs_time);
		    //my_class_peak_finder.set_debug(0);
		    //my_class_peak_finder.set_finding_parameters(0,2.0,1.2,4,1,0.0,6.0); //(0 = peak, 1 = rim), width, min S/B, exlcusion range, N_maxima, min_radius, max_width_scan
		    //my_class_peak_finder.find_peaks();
		    //vector< vector<Double_t> > vec_peak_positions = my_class_peak_finder.get_peak_positions();  // time bin, height, density
		    //printf("   ->peak at time bin: %f, height: %f, density: %f \n",vec_peak_positions[0][0],vec_peak_positions[1][0],vec_peak_positions[2][0]);

		    //my_class_peak_finder.clear();
		    //h_ADC_vs_time ->Reset();
		}
                //--------------------------------------------------


#if 0
		//--------------------------------------------------
		if(ADC_amplitude_sum_times > 0.0)
		{
		    cout << " " << endl;
		    cout << "--------------------------------------------------" << endl;
                    Double_t arr_raw[24];
		    for(int i = 0; i < 24; i++)
		    {
                        arr_raw[i] = arr[i];
		    }
		    // Tail cancelation
		    Float_t rates[2];
		    Float_t coefficients[2];

		    // Initialization (coefficient = alpha, rates = lambda)
		    Float_t r1 = 1.0;
		    Float_t r2 = 1.0;
		    Float_t c1 = 0.5;
		    Float_t c2 = 0.5;

		    r1 = 1.156;
		    r2 = 0.130;
		    c1 = 0.066;
		    c2 = 0.000;

		    coefficients[0] = c1;
		    coefficients[1] = c2;

		    Double_t dt = 0.1;

		    rates[0] = TMath::Exp(-dt/(r1));
		    rates[1] = 0.0;

		    Float_t reminder[2] = { .0, .0 };
		    Float_t correction = 0.0;
		    Float_t result     = 0.0;

		    Int_t fBaseline = 10;

		    for(int i = 0; i < 24; i++)
		    {
			result = arr[i] - correction - fBaseline;    // No rescaling
			arr[i] = (Short_t)(result + fBaseline + 0.5f);

			correction = 0.0;
			for(int k = 0; k < 2; k++)
			{
			    correction += reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
			}
		    }
		    for(int i = 0; i < 24; i++)
		    {
			cout << "time bin: " << i << ", after tail cancelation ADC: " << arr[i] << ", before: " << arr_raw[i] << endl;
		    }
		}
		//--------------------------------------------------
#endif



		//------------------------------------------------------------
		// Determine reference position (time bin with highest ADC value) for one single pad
		// Fill vector with ADC values for all 24 time bins for this pad
		if(ADC_amplitude_sum_times > 0.0 && max_time_bin > 0)
		{
		    //----------------------
		    // Calculate global position for fired TRD pad
		    Float_t              TRD_time0        = fGeo         ->GetTime0(i_layer); // in cm
		    AliTRDpadPlane*      padplane         = fGeo         ->GetPadPlane(i_det);
		    Double_t             TRD_col_end      = padplane     ->GetColEnd();
		    Double_t             TRD_row_end      = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
		    Double_t             TRD_row_start    = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
		    Double_t             TRD_row_end_ROC  = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
		    Double_t             TRD_col_spacing  = padplane     ->GetColSpacing();
		    Double_t             TRD_row_spacing  = padplane     ->GetRowSpacing();
		    Double_t             TRD_col_pos      = padplane     ->GetColPos(i_column);
		    Double_t             TRD_row_pos      = padplane     ->GetRowPos(i_row);       // fPadRow[row] + fPadRowSMOffset;
		    Double_t             TRD_row_pos_ROC  = padplane     ->GetRowPosROC(i_row);    // fPadRow[row]; = pad border position
		    Double_t             TRD_col_size     = padplane     ->GetColSize(i_column);
		    Double_t             TRD_row_size     = padplane     ->GetRowSize(i_row);
		    Float_t              TRD_loc_Y        = TRD_col_pos + 0.5*TRD_col_size;
		    Double_t             TRD_drift_time   = ((Double_t)(max_time_bin + 0.5))*TRD_time_per_bin*ChamberVdrift->GetValue(i_det); // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum
                    //Double_t             TRD_drift_time   = ((Double_t)12)*TRD_time_per_bin*ChamberVdrift->GetValue(i_det); // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum

		    Double_t kX0shift     = AliTRDgeometry::AnodePos(); //[cm]
                    Double_t fZShiftIdeal = 0.5 * (TRD_row_start + TRD_row_end);
		    // Get the detector wise defined calibration values
		    Double_t fCalVdriftDetValue = ChamberVdrift ->GetValue(i_det);
		    Double_t fCalT0DetValue     = ChamberT0     ->GetValue(i_det);
		    Double_t fCalExBDetValue    = ChamberExB    ->GetValue(i_det);

		    // Retrieve calibration values
		    // drift velocity
		    Double_t vd  = fCalVdriftDetValue * 1.0;
		    // t0
		    Double_t t0  = fCalT0DetValue     + 0.0;
                    Double_t fSamplingFrequency = fParam->GetSamplingFrequency();
		    t0 /= fSamplingFrequency;
		    // ExB correction
		    Double_t exb = fCalExBDetValue;//AliTRDCommonParam::Instance()->GetOmegaTau(vd);

		    //Float_t x = c->GetXloc(t0, vd);
		    //printf("AnodePos: %f, sf: %f, TRD_time0: %f \n",kX0shift,fSamplingFrequency,TRD_time0);

                    // http://personalpages.to.infn.it/~puccio/htmldoc/src/AliTRDtransform.cxx.html#LoiVTE
		    // Parameters to adjust the X position of clusters in the alignable volume
                    /*
		    const Double_t kX0shift = AliTRDgeometry::AnodePos(); //[cm]

		    // Retrieve calibration values
		    Int_t col = c->GetPadCol(), row = c->GetPadRow();
		    // drift velocity
		    Double_t vd  = fCalVdriftDetValue * fCalVdriftROC->GetValue(col,row);
		    // t0
		    Double_t t0  = fCalT0DetValue     + fCalT0ROC->GetValue(col,row);
		    t0 /= fSamplingFrequency;
		    // ExB correction
		    Double_t exb = fCalExBDetValue;//AliTRDCommonParam::Instance()->GetOmegaTau(vd);

		    Float_t x = c->GetXloc(t0, vd); // x = td*vd, td = (fPadTime + .5)/fFreq - t0 - 0.189; // [us]

		    // Pad dimensions
		    Double_t rs = fPadPlane->GetRowSize(row);
		    Double_t cs = fPadPlane->GetColSize(col);

		    // cluster error with diffusion corrections
		    Double_t s2  = cs*fCalPRFROC->GetValue(col, row);
		    s2 *= s2;
		    Float_t dl, dt;
		    AliTRDCommonParam::Instance()->GetDiffCoeff(dl, dt, vd);

		    Double_t y0 = fPadPlane->GetColPos(col) + .5*cs;
		    Double_t loc[] = {
			kX0shift-x,                    // Invert the X-position,
			c->GetYloc(y0, s2, cs) - x*exb,// apply ExB correction
			fPadPlane->GetRowPos(row) - .5*rs - fZShiftIdeal // move the Z-position relative to the middle of the chamber
		    };

		    // Go to tracking coordinates
		    Double_t trk[3];
		    fMatrix->LocalToMaster(loc, trk);
                    */

		    Float_t lorentz_angle_corr_y = ChamberExB->GetValue(i_det)*TRD_drift_time;
		    TRD_loc_Y -= Sign_magnetic_field*lorentz_angle_corr_y;

		    //printf("ExB correction: %f \n",Sign_magnetic_field*lorentz_angle_corr_y);

		    std::vector<Int_t> vec_det_info;
		    vec_det_info.resize(4);
		    vec_det_info[0] = i_det;
		    vec_det_info[1] = i_layer;
		    vec_det_info[2] = i_row;
		    vec_det_info[3] = i_column;

		    // Determine global position of TRD hit
		    Double_t             loc_ref[3]           = {TRD_time0 - TRD_drift_time,TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};
		    //Double_t             loc_ref[3]           = {kX0shift - (TRD_drift_time - TRD_time0),TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};

		    //printf("fZShiftIdeal: %f, old x: %f, new x: %f \n",fZShiftIdeal,TRD_time0 - TRD_drift_time,kX0shift - (TRD_drift_time - TRD_time0));

		    Double_t             glb_ref[3]           = {0.0,0.0,0.0};
		    fGeo ->RotateBack(i_det,loc_ref,glb_ref);

                    // Apply alignment
		    Double_t glb_ref_align_sec[3];
                    Double_t glb_ref_align[3];

		    TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_ref,glb_ref_align_sec);
                    TM_TRD_rotation_det[i_det].LocalToMaster(glb_ref_align_sec,glb_ref_align);

		    //cout << "old z: " << glb_ref[2] << ", new z: " <<  glb_ref_align[2] << ", delta z: " << TV3_TRD_translation[i_det].Z() << endl;

		    //glb_ref_align[0] += 1.0;
		    //glb_ref_align[1] += 1.0;
                    //glb_ref_align[2] += 1.0;

		    // Determine global pointing vector of TRD plane (vector perpendicular to plane in global system)
		    Double_t             loc_vec_ref[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
		    Double_t             glb_vec_ref[3]           = {0.0,0.0,0.0};
		    fGeo ->RotateBack(i_det,loc_vec_ref,glb_vec_ref);

		    // Apply alignment
                    Double_t glb_vec_ref_align_sec[3];
		    Double_t glb_vec_ref_align[3];
                    TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_vec_ref,glb_vec_ref_align_sec);
		    TM_TRD_rotation_det[i_det].LocalToMaster(glb_vec_ref_align_sec,glb_vec_ref_align);

		    //for(Int_t i = 0; i < 3; i++)
		    //{
		    //    glb_ref_align[i] = glb_ref[i];
		    //    glb_vec_ref_align[i] = glb_vec_ref[i];
		    //}

		    TVector3 TV3_TRD_hit_middle, TV3_TRD_hit_det_angle_middle;
		    TV3_TRD_hit_middle.SetXYZ(glb_ref_align[0],glb_ref_align[1],glb_ref_align[2]);
		    TV3_TRD_hit_det_angle_middle.SetXYZ(glb_vec_ref_align[0],glb_vec_ref_align[1],0.0); // glb_vec is from {0,0,0} to center of detector, set z component to 0 to get vector perpendicular to detector plane
		    TV3_TRD_hits_middle.push_back(TV3_TRD_hit_middle);
		    TV3_TRD_hits_det_angle_middle.push_back(TV3_TRD_hit_det_angle_middle);
		    vec_TRD_hits_det_lay_row_col.push_back(vec_det_info);

		    vec_TRD_hits_ADC_values_time.push_back(vec_ADC_time_bins);
                    vec_TRD_hits_points_time.push_back(vec_points_time_bins);
                    //vec_TRD_hits_points_time_uncalib.push_back(vec_points_time_bins);
		}
		//------------------------------------------------------------

	    } // end of row loop
	} // end of column loop
    } // end of detector loop

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        AS_Event ->setADC_sum_det(i_det,sum_full_ADC_digit_det[i_det]);
    }
    //cout << "max_N_columns: " << max_N_columns << ", max_N_rows: " << max_N_rows << endl;
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------
    // Track loop
    //cout << "" << endl;
    //cout << "-----------------------------------------------------------------" << endl;
    cout << "Start matching " << N_tracks << " TPC tracks with " << TV3_TRD_hits_middle.size() << " TRD pads" << endl;
    N_good_tracks = 0;
    Int_t N_matched_TRD_hits_total = 0;
    for(Int_t iTracks = 0; iTracks < N_tracks; iTracks++)
    {
	//---------------------------------------------------------------
	// Gather track information

	// We always want the ESD track
	AliESDtrack* track = fESD->GetTrack(iTracks);
	if(!track)
	{
	    printf("ERROR: Could not receive track %d\n", iTracks);
	    continue;
	}


	if(!EsdTrackCuts->AcceptTrack(track)) continue;

	Double_t TRD_signal   = track ->GetTRDsignal(); // truncated mean signal?
	Double_t Track_pT     = track ->Pt();
	Double_t Track_p      = track ->P();
        Double_t p_vec[3];
	track->GetPxPyPz(p_vec);
	Int_t    charge       = track ->Charge();
	Double_t Track_phi    = track ->Phi();
	Double_t Track_theta  = track ->Theta();
	Double_t Track_eta    = track ->Eta();
	Double_t TPC_chi2     = track ->GetTPCchi2();
	Double_t TPC_signal   = track ->GetTPCsignal(); // dE/dx?
	Double_t TOF_signal   = track ->GetTOFsignal(); // time-of-flight?
        Double_t Track_length = track ->GetIntegratedLength();
	UShort_t N_TPC_cls    = track ->GetTPCNcls();

	//AliKalmanTrack* trd_track = track ->GetTRDtrack();

	h2D_TPC_dEdx_vs_momentum ->Fill(Track_p,TPC_signal);

	//if(Track_pT > 0.5) printf("TOF_signal: %f, eta: %f \n",TOF_signal,Track_eta);


        Int_t pT_bin;
	for(Int_t i_pT = 0; i_pT < N_pT_bins; i_pT++)
	{
	    if(Track_pT >= pT_ranges[i_pT] && Track_pT < pT_ranges[i_pT+1])
	    {
                pT_bin = i_pT;
                break;
	    }
	}

	TLorentzVector TLV_p_vec;
	Double_t p_vec_energy = TMath::Sqrt(p_vec[0]*p_vec[0] + p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2] + 0.938*0.938);
	TLV_p_vec.SetPxPyPzE(p_vec[0],p_vec[1],p_vec[2],p_vec_energy);
	//cout << "TLV_p_vec.P: " << TLV_p_vec.P() << ", P: " << Track_p << ", TLV_p_vec.Theta: " << TLV_p_vec.Theta() << ", Theta: " << Track_theta
	//<< ", TLV_p_vec.Phi: " << TLV_p_vec.Phi() << ", phi: " << Track_phi  << endl;


	ULong_t status = track->GetStatus();
	Int_t ITS_refit = 0;
        Int_t TPC_refit = 0;
        Int_t track_status = 0;
	if(((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit))
	{
	    ITS_refit = 1;
	    track_status |= 1 << 0; // setting bit 0 to 1
	}
	if(((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit))
	{
	    TPC_refit = 1;
	    track_status |= 1 << 1; // setting bit 1 to 1
	}


          // e = 0, muon = 1, pion = 2, kaon = 3, proton = 4
	Double_t Track_PID[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        // nSigma TPC
	Track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
	Track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kMuon);
	Track_PID[2] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
	Track_PID[3] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
	Track_PID[4] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);

        // nSigma TOF, -999 in case there is no TOF hit
	Track_PID[5] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
	Track_PID[6] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kMuon);
	Track_PID[7] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
	Track_PID[8] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
	Track_PID[9] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);


	Float_t track_xy_impact,track_z_impact;
	track->GetImpactParameters(track_xy_impact,track_z_impact);
	Double_t track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);


	Double_t TRD_ADC_bin_width = 100.0;

	//-------------------
	Int_t N_ITS_cls = 0;
	for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
	{
	    if(track ->HasPointOnITSLayer(i_ITS_layer))
	    {
		N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
	    }
	}
	//-------------------

	TLorentzVector TL_vec;
	TL_vec.SetPtEtaPhiM(Track_pT,Track_eta,Track_phi,mass_array[2]);
	AS_Track  = AS_Event ->createTrack();
	AS_Track  ->clearTRD_digit_list();
	AS_Track  ->set_TLV_part(TL_vec);
	AS_Track  ->setdca(((Double_t)charge)*track_total_impact);
	AS_Track  ->setnsigma_e_TPC(Track_PID[0]);
	AS_Track  ->setnsigma_e_TOF(Track_PID[5]);
	AS_Track  ->setnsigma_pi_TPC(Track_PID[2]);
	AS_Track  ->setnsigma_pi_TOF(Track_PID[7]);
	AS_Track  ->setnsigma_K_TPC(Track_PID[3]);
	AS_Track  ->setnsigma_K_TOF(Track_PID[8]);
	AS_Track  ->setnsigma_p_TPC(Track_PID[4]);
	AS_Track  ->setnsigma_p_TOF(Track_PID[9]);
	AS_Track  ->setTRDSignal(TRD_signal);
	AS_Track  ->setNTPCcls(N_TPC_cls);
	AS_Track  ->setNITScls(N_ITS_cls);
	AS_Track  ->setStatus(track_status);
	AS_Track  ->setTPCchi2(TPC_chi2);
	AS_Track  ->setTPCdEdx(TPC_signal);
	AS_Track  ->setTOFsignal(TOF_signal);
        AS_Track  ->setTrack_length(Track_length);

	//-------------------
	// Get TRD information
	Int_t N_TRD_cls = 0;
	Double_t TRD_sum_ADC = 0.0;
	Long64_t TRD_layer_info[6]; // each value stores all 8 time slices, Long64_t has 8 byte = 8*8 bit = 64 bit in total -> one byte = 8 bit = 256 per time bin
	memset(TRD_layer_info, 0, sizeof(TRD_layer_info)); // for automatically-allocated arrays
	for(Int_t iPl = 0; iPl < 6; iPl++) // layer
	{
	    AS_Track  ->setTRD_layer(iPl,TRD_layer_info[iPl]); // set all values to 0
	}

	//printf("track: %d \n",iTracks);
	if(TRD_signal > 0.0)
	{
	    //------------------------------
	    // Get TRD PID information from ESD track
	    for(Int_t iPl = 0; iPl < 6; iPl++) // layer
	    {
		TRD_layer_info[iPl] = 0;
		Double_t TRD_momentum = track->GetTRDmomentum(iPl);
		//cout << "" << endl;
		//cout << "--------------------------------" << endl;
                Double_t sum_ADC = 0.0;
		for(int isl = 0; isl <= 7; isl++) // time slice
		{
                    Double_t TRD_ADC_time_slice = track->GetTRDslice(iPl,isl);

                    if(isl == 0 && TPC_signal > 50.0 && TPC_signal < 70.0)
                    {
                        h_ADC_tracklet[0] ->Fill(TRD_ADC_time_slice/TPC_signal);
                        //printf("iTracks: %d, dE/dx: %4.2f, layer: %d, ADC/dEdx: %4.2f \n",iTracks,TPC_signal,iPl,TRD_ADC_time_slice/TPC_signal);
                    }

		    //if(TRD_ADC_time_slice > 20000.0) printf("layer: %d, time: %d, ADC: %f \n",iPl,isl,TRD_ADC_time_slice);
                    sum_ADC += TRD_ADC_time_slice;
		    Int_t TRD_ADC_time_slice_byte = (Int_t)round((Double_t)TRD_ADC_time_slice/TRD_ADC_bin_width);
		    if(TRD_ADC_time_slice_byte > 255) TRD_ADC_time_slice_byte = 255;

		    //number |= 1 << x; // setting bit x to 1
		    //bit = (number >> x) & 1; // check bit x

		    for(Int_t i_bit = 0; i_bit < 8; i_bit++) // One single time slice 8 bit = 256
		    {
			Int_t bit_status = (TRD_ADC_time_slice_byte >> i_bit) & 1; // check bit i_bit
			Int_t bitset = i_bit + 8*isl; // range: 0..63 = 64 bit = 8 byte = Long64_t
			if(bit_status) TRD_layer_info[iPl] |= (ULong64_t)1 << bitset; // setting bit bitset to 1
		    }

		    if(TRD_ADC_time_slice > 0.0)
		    {
			TRD_sum_ADC += TRD_ADC_time_slice;
			N_TRD_cls++;
		    }
		    //cout << "isl: " << isl << ", TRD_ADC_time_slice: " << TRD_ADC_time_slice << ", TRD_ADC_time_slice_byte: " << TRD_ADC_time_slice_byte << endl;
		}
		Double_t average_sum_ADC = sum_ADC/7.0;
		//printf("sum_ADC: %f, average_sum_ADC: %f \n",sum_ADC,average_sum_ADC);

                AS_Track  ->setTRD_layer(iPl,TRD_layer_info[iPl]);
                //for(int isl = 0; isl <= 7; isl++) // time slice
                //{
                //    Float_t rec_TRD_AdC = AS_Track->getTRD_ADC(iPl,isl);
                //    if(rec_TRD_AdC > 19500) printf("    ->layer: %d, time: %d, rec ADC: %f \n",iPl,isl,rec_TRD_AdC);
                //}

		// Check the decoding
		//cout << "" << endl;
		for(int isl = 0; isl <= 7; isl++) // time slice
		{
		    // Decode TRD_layer_info back to TRD ADC values
		    ULong64_t TRD_value = 0;
		    for(Int_t i_bit = 0; i_bit < 8; i_bit++) // One single time slice 8 bit = 256
		    {
			Int_t bitcheck = i_bit + 8*isl; // range: 0..63 = 64 bit = 8 byte = Long64_t
			Int_t bit_status = (TRD_layer_info[iPl] >> bitcheck) & 1; // check bit bitcheck
			if(bit_status) TRD_value |= (ULong64_t)1 << i_bit; // setting bit i_bit to 1
		    }
		    Double_t TRD_value_decode = (Double_t)TRD_value * TRD_ADC_bin_width;
		    //cout << "isl: " << isl << ", TRD_value_decode: " << TRD_value_decode << endl;
		}
		//cout << "--------------------------------" << endl;

	    }
	}

	AS_Track  ->setNTRDcls(N_TRD_cls);
	AS_Track  ->setTRDsumADC(TRD_sum_ADC);
	//-------------------




	//Helix
        FillHelix(track,magF);

        AS_Track ->setHelix(aliHelix.fHelix[0],aliHelix.fHelix[1],aliHelix.fHelix[2],aliHelix.fHelix[3],aliHelix.fHelix[4],aliHelix.fHelix[5],aliHelix.fHelix[6],aliHelix.fHelix[7],aliHelix.fHelix[8]);

	Int_t N_track_TRD_tracklets     = track     ->GetTRDntracklets();
	//printf("N_track_TRD_tracklets: %d \n",N_track_TRD_tracklets);

	if(N_good_tracks >= 0 && N_good_tracks < 60000)
	{
            Double_t max_radius_reached = -1.0;
            Int_t flag_reached_TRD  = 1;
	    Double_t path_initA = 0.0;
	    Double_t helix_point_search[3];
	    for(Int_t i_step = 0; i_step < 2000; i_step++)
	    {
		Double_t pathA_step  = (TPC_radius_scan-20.0) + (Double_t)i_step*5.0; // use this for GetExternalParameters
                //Double_t pathA_step  = (0.0-20.0) + (Double_t)i_step*5.0; // use this for GetExternalOuterParameters
		path_initA = pathA_step;
		aliHelix.Evaluate(pathA_step,helix_point_search);
		Double_t Track_radius = TMath::Sqrt(helix_point_search[0]*helix_point_search[0] + helix_point_search[1]*helix_point_search[1]);
		//cout << "i_step: " << i_step << ", pT: " << Track_pT << ", path: " << pathA_step << ", radius: " << Track_radius << ", charge: " << charge
		//    << ", xyz = {" << helix_point_search[0] << ", " << helix_point_search[1] << ", " << helix_point_search[2] << "}" << endl;
                if(Track_radius > max_radius_reached) max_radius_reached = Track_radius;
		if(Track_radius > TPC_radius_scan) break;
		if(Track_radius < max_radius_reached)
		{
                    flag_reached_TRD = 0;
		    break;
		}
	    }

	    //cout << "flag_reached_TRD: " << flag_reached_TRD << endl;
	    if(!flag_reached_TRD) continue; // skip tracks which didn't make it to the TRD due to low momentum



#if 0
            //--------------------------------------------------------------------------------
            // Only for test purposes so far
	    // Loop over all TRD hits in the event and match them with the TPC track
	    // WARNING: All time bins are matched!
	    for(Int_t i_TRD_hit = 0; i_TRD_hit < TV3_TRD_hits.size(); i_TRD_hit++) // All ADC hits, including time bins
	    {
		Double_t dx = helix_point_search[0] - TV3_TRD_hits[i_TRD_hit].X();
		Double_t dy = helix_point_search[1] - TV3_TRD_hits[i_TRD_hit].Y();
		Double_t dz = helix_point_search[2] - TV3_TRD_hits[i_TRD_hit].Z();
		Double_t dist_to_TRD_pad = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

		if(dist_to_TRD_pad > 65.0) continue;

		Float_t pathA, dcaAB;
		FindDCAHelixPoint(TV3_TRD_hits[i_TRD_hit],aliHelix,path_initA-25.0,path_initA+25.0,pathA,dcaAB);
		Double_t helix_point[3];
		aliHelix.Evaluate(pathA,helix_point);
		h_dca[pT_bin] ->Fill(dcaAB);
		for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
		{
		    h_dca_xyz[i_xyz][pT_bin] ->Fill(helix_point[i_xyz]-TV3_TRD_hits[i_TRD_hit][i_xyz]);
		    if(i_xyz == 0
		       && fabs(helix_point[1]-TV3_TRD_hits[i_TRD_hit][1]) < 3.0
		       && fabs(helix_point[2]-TV3_TRD_hits[i_TRD_hit][2]) < 7.0
		      )
		    {
			h_dca_xyz[3][pT_bin] ->Fill(helix_point[i_xyz]-TV3_TRD_hits[i_TRD_hit][i_xyz]);
		    }
		    if(i_xyz == 1
		       && fabs(helix_point[0]-TV3_TRD_hits[i_TRD_hit][0]) < 3.0
		       && fabs(helix_point[2]-TV3_TRD_hits[i_TRD_hit][2]) < 7.0
		      )
		    {
			h_dca_xyz[4][pT_bin] ->Fill(helix_point[i_xyz]-TV3_TRD_hits[i_TRD_hit][i_xyz]);
		    }
		    if(i_xyz == 2
		       && fabs(helix_point[0]-TV3_TRD_hits[i_TRD_hit][0]) < 3.0
		       && fabs(helix_point[1]-TV3_TRD_hits[i_TRD_hit][1]) < 3.0
		      )
		    {
			h_dca_xyz[5][pT_bin] ->Fill(helix_point[i_xyz]-TV3_TRD_hits[i_TRD_hit][i_xyz]);
		    }
		}
	    }
            //--------------------------------------------------------------------------------
#endif



	    //--------------------------------------------------------------------------------
	    // Loop over all TRD middle hits in the event and match them with the TPC track
	    Int_t TRD_layer_match[6] = {0,0,0,0,0,0};
            Double_t Impact_angle_first = -100.0;
            Double_t min_dca = 100000.0;
	    for(Int_t i_TRD_hit = 0; i_TRD_hit < TV3_TRD_hits_middle.size(); i_TRD_hit++) // ADC (average for all time bins) hit for a single pad
	    {
		Double_t dx = helix_point_search[0] - TV3_TRD_hits_middle[i_TRD_hit].X();
		Double_t dy = helix_point_search[1] - TV3_TRD_hits_middle[i_TRD_hit].Y();
		Double_t dz = helix_point_search[2] - TV3_TRD_hits_middle[i_TRD_hit].Z();
		Double_t dist_to_TRD_pad = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

		if(dist_to_TRD_pad > 65.0) continue;

		Float_t pathA, dcaAB;
		FindDCAHelixPoint(TV3_TRD_hits_middle[i_TRD_hit],aliHelix,path_initA-25.0,path_initA+25.0,pathA,dcaAB);
                if(dcaAB < min_dca) min_dca = dcaAB;

                Double_t dca_xyz[3];
		Double_t helix_point_dca[3];
		aliHelix.Evaluate(pathA,helix_point_dca);
		h_dca[pT_bin] ->Fill(dcaAB);
		for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
		{
                    dca_xyz[i_xyz] = helix_point_dca[i_xyz]-TV3_TRD_hits_middle[i_TRD_hit][i_xyz];
		    h_dca_xyz[i_xyz][pT_bin] ->Fill(helix_point_dca[i_xyz]-TV3_TRD_hits_middle[i_TRD_hit][i_xyz]);
		    if(i_xyz == 0
		       && fabs(helix_point_dca[1]-TV3_TRD_hits_middle[i_TRD_hit][1]) < 3.0
		       && fabs(helix_point_dca[2]-TV3_TRD_hits_middle[i_TRD_hit][2]) < 7.0
		      )
		    {
			h_dca_xyz[3][pT_bin] ->Fill(helix_point_dca[i_xyz]-TV3_TRD_hits_middle[i_TRD_hit][i_xyz]);
		    }
		    if(i_xyz == 1
		       && fabs(helix_point_dca[0]-TV3_TRD_hits_middle[i_TRD_hit][0]) < 3.0
		       && fabs(helix_point_dca[2]-TV3_TRD_hits_middle[i_TRD_hit][2]) < 7.0
		      )
		    {
			h_dca_xyz[4][pT_bin] ->Fill(helix_point_dca[i_xyz]-TV3_TRD_hits_middle[i_TRD_hit][i_xyz]);
		    }
		    if(i_xyz == 2
		       && fabs(helix_point_dca[0]-TV3_TRD_hits_middle[i_TRD_hit][0]) < 3.0
		       && fabs(helix_point_dca[1]-TV3_TRD_hits_middle[i_TRD_hit][1]) < 3.0
		      )
		    {
			h_dca_xyz[5][pT_bin] ->Fill(helix_point_dca[i_xyz]-TV3_TRD_hits_middle[i_TRD_hit][i_xyz]);
		    }
		}


		//Double_t dca_yz[2];
		//Bool_t b_prop_to_dca =  track->PropagateToDCA(PrimVertex,magF,3.0,dca_yz);
		//Double_t dca_yz_total = TMath::Sqrt(dca_yz[0]*dca_yz[0] + dca_yz[1]*dca_yz[1]);

		//cout << "pathA: " << pathA << ", dcaAB, " << dcaAB << endl;
		Double_t helix_point[3];
		aliHelix.Evaluate(pathA,helix_point);

		Double_t helix_point_B[3];
		aliHelix.Evaluate(pathA-1.0,helix_point_B);

		TVector3 TV3_track_impact_vec_on_TRD;
		TV3_track_impact_vec_on_TRD.SetXYZ(helix_point[0]-helix_point_B[0],helix_point[1]-helix_point_B[1],helix_point[2]-helix_point_B[2]);

                Int_t N_good_match = 0;
		if(dcaAB < TPC_TRD_matching_window) // TRD hit matched with TPC track
		{
		    Int_t    TRD_det    = vec_TRD_hits_det_lay_row_col[i_TRD_hit][0];
		    Int_t    TRD_lay    = vec_TRD_hits_det_lay_row_col[i_TRD_hit][1];
		    Int_t    TRD_row    = vec_TRD_hits_det_lay_row_col[i_TRD_hit][2];
		    Int_t    TRD_col    = vec_TRD_hits_det_lay_row_col[i_TRD_hit][3];
		    Int_t    TRD_sec    = fGeo->GetSector(TRD_det);
		    Int_t    TRD_stack  = fGeo->GetStack(TRD_det);

                    Int_t fBaseline = fDigMan->GetDigitsParam()->GetADCbaseline(TRD_det); // ADC baseline to be subtracted from digit ADC values, constant value of 10

		    Int_t x_TRD = TRD_col + TRD_sec*144;
		    Int_t y_TRD = TRD_row + TRD_stack*16 + TRD_lay*16*5;

		    //Int_t row_calc = (y_TRD%16);
		    //printf("iTracks: %d, det: %d, lay: %d, sec: %d, stack: %d, row: {%d,%d}, col: %d \n",iTracks,TRD_det,TRD_lay,TRD_sec,TRD_stack,TRD_row,row_calc,TRD_col);
                   
		    TRD_layer_match[TRD_lay] = 1;

		    Double_t Impact_angle = TV3_TRD_hits_det_angle_middle[i_TRD_hit].Angle(TV3_track_impact_vec_on_TRD); // in rad
		    if(Impact_angle > TMath::Pi()/2.0) Impact_angle -= TMath::Pi()/2.0; // Impact angle of TPC track on TRD hit
                    if(N_good_match == 0) Impact_angle_first = Impact_angle;
		    Double_t Impact_distance_in_drift = 3.0;
		    if(TMath::Cos(Impact_angle) != 0.0)
		    {
			Impact_distance_in_drift = 3.0/TMath::Cos(Impact_angle); // 3.0 cm drift length
		    }
		    Double_t Impact_angle_degree = Impact_angle*TMath::RadToDeg();


		    //-------------------------------------------------------------------
		    AS_Digit = AS_Track ->createTRD_digit();
		    AS_Digit ->sethit_ids(x_TRD,y_TRD);
                    AS_Digit ->setdca_to_track(dcaAB,dca_xyz[0],dca_xyz[1],dca_xyz[2]);
		    AS_Digit ->setImpactAngle(Impact_angle);

                    Double_t Impact_distance_in_drift_angle = TRD_Impact_distance_in_drift/TMath::Cos(Impact_angle); // 3.0 cm drift length + 0.5*0.7 cm amplification region

		    //cout << "iTrack: " << iTracks << ", i_TRD_hit: " << i_TRD_hit << endl;

		    Int_t    time_max_bin  = 0;
		    Double_t ADC_value_max = 0.0;

                    Short_t ADC_Digit_values_array_tc[24];
		    for(Int_t i_time = 0; i_time < vec_TRD_hits_ADC_values_time[i_TRD_hit].size(); i_time++)
		    {
			if(i_time < 24)
			{
			    ADC_Digit_values_array_tc[i_time] = vec_TRD_hits_ADC_values_time[i_TRD_hit][i_time];
			}
		    }
                    func_tail_cancellation(ADC_Digit_values_array_tc,1);
		    for(Int_t i_time = 0; i_time < vec_TRD_hits_ADC_values_time[i_TRD_hit].size(); i_time++)
		    {
			Double_t ADC_value  = vec_TRD_hits_ADC_values_time[i_TRD_hit][i_time];

                        CalROC = PadNoise ->GetCalROC(TRD_det);
			Double_t pad_noise = CalROC ->GetValue(TRD_col,TRD_row);
			Double_t pad_gain_factor = KryptoGain ->GetGainCorrectionFactor(TRD_det,TRD_row,TRD_col);
			Double_t detector_gain = chambergain->GetValue(TRD_det);
			Double_t ADC_value_corrected    = ((ADC_value - fBaseline + pad_noise)/(detector_gain*pad_gain_factor))/Impact_distance_in_drift_angle;
			Double_t ADC_value_corrected_tc = ((ADC_Digit_values_array_tc[i_time] - fBaseline + pad_noise)/(detector_gain*pad_gain_factor))/Impact_distance_in_drift_angle;
		        //cout << "i_time: " << i_time << ", ADC_value: " << ADC_value << ", fBaseline: " << fBaseline << ", pad_noise: " << pad_noise
			//    << ", detector_gain: " << detector_gain << ", pad_gain_factor: " << pad_gain_factor <<
			//    ", Impact_distance_in_drift_angle: " << Impact_distance_in_drift_angle << ", ADC_value_corrected: " << ADC_value_corrected << endl;
			if(ADC_value_corrected > ADC_value_max)
			{
			    ADC_value_max = ADC_value_corrected;
                            time_max_bin  = i_time;
                        }

                        AS_Digit ->setADC_time_value(i_time,(Short_t)ADC_value);
                        //printf("ADC_value: %hi \n",(Short_t)ADC_value);
                        //AS_Digit ->setADC_time_value_corrected(i_time,(Short_t)ADC_value_corrected);

                        AS_Digit ->set_pos(i_time,vec_TRD_hits_points_time[i_TRD_hit][i_time].X(), vec_TRD_hits_points_time[i_TRD_hit][i_time].Y(), vec_TRD_hits_points_time[i_TRD_hit][i_time].Z());

                        //cout << "vec_TRD_hits_points_time[i_TRD_hit][i_time].X()" << vec_TRD_hits_points_time[i_TRD_hit][i_time].X() << endl;
                        //cout << << << endl;

                        // AS_Digit ->set_pos_uncalib(i_time,vec_TRD_hits_points_time_uncalib[i_TRD_hit][i_time].X(), vec_TRD_hits_points_time_uncalib[i_TRD_hit][i_time].Y(), vec_TRD_hits_points_time_uncalib[i_TRD_hit][i_time].Z());

                        //if(fEventNoInFile < 2 && iTracks == 3) printf("iTracks: %d, layer: %d, i_time: %d, pos: {%4.3f, %4.3f, %4.3f} \n",iTracks,TRD_lay,i_time,vec_TRD_hits_points_time[i_TRD_hit][i_time].X(), vec_TRD_hits_points_time[i_TRD_hit][i_time].Y(), vec_TRD_hits_points_time[i_TRD_hit][i_time].Z());

                        //here we need to save also uncalibrated 3D points vec_TRD_hits_points_time_uncalib

                        //AS_Digit ->setADC_time_value_corrected_tc(i_time,(Short_t)ADC_value_corrected_tc);
		    }

		    Float_t pathA_max, dcaAB_max;
		    FindDCAHelixPoint(vec_TRD_hits_points_time[i_TRD_hit][time_max_bin],aliHelix,path_initA-25.0,path_initA+25.0,pathA_max,dcaAB_max);

		    //cout << "ADC_value_max: " << ADC_value_max << ", time_max_bin: " << time_max_bin << ", dcaAB_max: " << dcaAB_max << endl;

#if 0
		    if(ADC_value_max > 30.0)
		    {
			Double_t helix_point[3];
			aliHelix.Evaluate(pathA_max,helix_point);
			h_dca[pT_bin] ->Fill(dcaAB_max);
			for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
			{
			    h_dca_xyz[i_xyz][pT_bin] ->Fill(helix_point[i_xyz]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][i_xyz]);
			    if(i_xyz == 0
			       && fabs(helix_point[1]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][1]) < 3.0
			       && fabs(helix_point[2]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][2]) < 7.0
			      )
			    {
				h_dca_xyz[3][pT_bin] ->Fill(helix_point[i_xyz]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][i_xyz]);
			    }
			    if(i_xyz == 1
			       && fabs(helix_point[0]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][0]) < 3.0
			       && fabs(helix_point[2]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][2]) < 7.0
			      )
			    {
				h_dca_xyz[4][pT_bin] ->Fill(helix_point[i_xyz]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][i_xyz]);
			    }
			    if(i_xyz == 2
			       && fabs(helix_point[0]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][0]) < 3.0
			       && fabs(helix_point[1]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][1]) < 3.0
			      )
			    {
				h_dca_xyz[5][pT_bin] ->Fill(helix_point[i_xyz]-vec_TRD_hits_points_time[i_TRD_hit][time_max_bin][i_xyz]);
			    }
			}
		    }
#endif
		    //-------------------------------------------------------------------

		    N_good_match++;

		} // end of matched hit

	    } // End of TRD hit loop
	    AS_Track  ->setimpact_angle_on_TRD(Impact_angle_first);

	    //cout << "min_dca: " << min_dca << endl;
	    //--------------------------------------------------------------------------------

	}  // N Good tracks >= ...

#if 0
	cout << "Track accepted, N_good_tracks: " << N_good_tracks << ", zpos: " << zpos << ", pT: " << Track_pT
	    << ", radius: " << radius << ", radius_inner: " << radius_inner << endl;
#endif
	N_good_tracks++;

    } // End of TPC track loop



    Tree_AS_Event ->Fill();

    N_good_events++;

    //-----------------------------------------------------------------
    //cout << "" << endl;
    //cout << "" << endl;
    //cout << "----------------------------------------------------------------------------------------" << endl;
    //cout << "Event number: " << fEventNoInFile << endl;
    //-----------------------------------------------------------------

    //PostData(1,fListOfHistos);
    //PostData(0,Tree_AS_Event);

}



//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::Terminate(Option_t *)
{
    cout << "In terminate" << endl;
    // Draw result to the screen
    // Called once at the end of the query

    //  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    //  if (!fOutputList) {
    //    printf("ERROR: Output list not available\n");
    //    return;
    //  }
    //
    //  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
    //  if (!fHistPt) {
    //    printf("ERROR: fHistPt not available\n");
    //    return;
    //  }
    //
    //TCanvas *c1 = new TCanvas("Ali_AS_analysis_TRD_digits","Pt",10,10,510,510);
    //c1->cd(1)->SetLogy();
    //fHistPt->DrawCopy("E");
}


//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::FillHelix(AliESDtrack* track_in, Double_t magF_in)
{
    //-------------------
    // Get helix
    // Track parametrization:
    // https://www.physi.uni-heidelberg.de/~sma/alice/LukasLayer_bachelor.pdf
    Double_t alpha, alpha_alt, x_param, p_param[5]; // p are the track paramters, 0 = Y, 1 = Z, 2 = Snp, 3 = Tgl (p[2]/pt), 4 = Signed1Pt
    track_in->GetExternalParameters(x_param,p_param); // Those are the parameters used for vertexing
    alpha_alt=track_in->GetAlpha();

    //cout << "x_param: " << x_param << endl;
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    cout << "A i_param: " << i_param << ", param: " << p_param[i_param] << endl;
    //}

    //track_in->GetOuterExternalParameters(alpha,x_param,p_param); //
    //alpha_alt = alpha;
    Double_t fX     = track_in->GetX();

    //cout << "x_param: " << x_param << endl;
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    cout << "B i_param: " << i_param << ", param: " << p_param[i_param] << endl;
    //}

    //cout << "x_param: " << x_param << endl;

    //const AliExternalTrackParam* TPC_track_param = track_in->GetOuterParam();
    //const Double_t* p_param_b = TPC_track_param->GetParameter();
    //const Double_t* p_param_b = track_in->GetOuterParam()->GetParameter();
    //const Double_t* p_param_b = track_in->GetInnerParam()->GetParameter();
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    p_param[i_param] = p_param_b[i_param];
    //}
    //x_param = fX;

    //-------------------
    // Correct way of filling aliHelix from
    // http://personalpages.to.infn.it/~puccio/htmldoc/src/AliHelix.cxx.html#PalT1E
    // line 52

    // CONCLUSION: AliTracker::GetBz() is not identical to GetC(magF_in), magnetic field is slightly different but it doesn't matter...

    Double_t fHelix_alt[9];
    Double_t x_alt,cs_alt,sn_alt;
    //track_in->GetExternalParameters(x_alt,fHelix_alt); // Those are the parameters used for vertexing
    x_alt = x_param;

    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
	fHelix_alt[i_param] = p_param[i_param];
    }

    //cout << "alpha: " << alpha << ", alpha_alt: " << alpha_alt << endl;

    //
    //circle parameters
    //PH Sometimes fP4 and fHelix[4] are very big and the calculation
    //PH of the Sqrt cannot be done. To be investigated...

    // kB2C=-0.299792458e-3; // from /AliRoot/STEER/STEERBase/AliVParticle.h
    //Double_t kB2C_test =-0.299792458e-3;
    //Double_t GetC(Double_t b) const
    //{return fP[4]*b*kB2C;}

    //Double_t par4test = fHelix_alt[4]*magF_in*kB2C_test;
    fHelix_alt[4] = track_in->GetC(magF_in);
    //fHelix_alt[4] = par4test; // take the one with the magnetic field directly from the ESD file

    cs_alt = TMath::Cos(alpha_alt);
    sn_alt = TMath::Sin(alpha_alt);

    Double_t xc_alt, yc_alt, rc_alt;
    rc_alt  =  1/fHelix_alt[4];
    xc_alt  =  x_alt-fHelix_alt[2]*rc_alt;
    Double_t dummy = 1-(x_alt-xc_alt)*(x_alt-xc_alt)*fHelix_alt[4]*fHelix_alt[4];
    yc_alt  =  fHelix_alt[0]+TMath::Sqrt(dummy)/fHelix_alt[4];

    fHelix_alt[6] = xc_alt*cs_alt - yc_alt*sn_alt;
    fHelix_alt[7] = xc_alt*sn_alt + yc_alt*cs_alt;
    fHelix_alt[8] =  TMath::Abs(rc_alt);
    //
    //
    fHelix_alt[5]=x_alt*cs_alt - fHelix_alt[0]*sn_alt;            // x0
    fHelix_alt[0]=x_alt*sn_alt + fHelix_alt[0]*cs_alt;            // y0
    fHelix_alt[2]=TMath::ATan2(-(fHelix_alt[5]-fHelix_alt[6]),fHelix_alt[0]-fHelix_alt[7]); // phi0
    if (fHelix_alt[4]>0) fHelix_alt[2]-=TMath::Pi();
    fHelix_alt[5]   = fHelix_alt[6];
    fHelix_alt[0]   = fHelix_alt[7];

    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
	aliHelix.fHelix[i_param] = fHelix_alt[i_param];
    }
    //-------------------
}



//________________________________________________________________________
void Ali_AS_analysis_TRD_digits::FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB)
{
    // V1.0
    Float_t pA[2] = {path_initA,path_initB}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    for(Int_t r = 0; r < 2; r++)
    {
	Double_t helix_point[3];
	helixA.Evaluate(pA[r],helix_point);
	testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[r]
	distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 30.0;
    while(fabs(scale_length) > 0.1 && loopcounter < 100) // stops when the length is too small
    {
	//cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
	//    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
	//    << ", d[1] = " << distarray[1] << ", flip = " << flip
	//    << ", scale_length = " << scale_length << endl;
	if(distarray[0] > distarray[1])
	{
	    if(loopcounter != 0)
	    {
		if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
		else scale = 0.7; // go on in this direction but only by the way * 0.7
	    }
	    scale_length = (pA[1]-pA[0])*scale; // the next length interval
	    pA[0]     = pA[1] + scale_length; // the new path

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[0],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[0]
	    distarray[0] = (testA-space_vec).Mag(); // new dca
	    flip = 1.0;
	}
	else
	{
	    if(loopcounter != 0)
	    {
		if(flip == -1.0) scale = 0.4;
		else scale = 0.7;
	    }
	    scale_length = (pA[0]-pA[1])*scale;
	    pA[1]     = pA[0] + scale_length;

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[1],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[1]
	    distarray[1] = (testA-space_vec).Mag();
	    flip = -1.0;
	}
	loopcounter++;
    }

    if(loopcounter >= 100) cout << "WARNING: FindDCAHelixPoint exceeded maximum of 100 loops" << endl;

    if(distarray[0] < distarray[1])
    {
	pathA = pA[0];
	dcaAB = distarray[0];
    }
    else
    {
	pathA = pA[1];
	dcaAB = distarray[1];
    }
}


//________________________________________________________________________
Bool_t Ali_AS_analysis_TRD_digits::ReadDigits()
{
    //cout << "In ReadDigits" << endl;
    // don't do anything if the digits have already been loaded
    if (fDigitsLoadedFlag) return kTRUE;

    if(!fDigMan)
    {
	AliError("no digits manager");
	return kFALSE;
    }

    // reset digit arrays
    for(Int_t det=0; det<540; det++)
    {
	fDigMan->ClearArrays(det);
	fDigMan->ClearIndexes(det);
    }


    if(!fDigitsInputFile)
    {
	AliError("digits file not available");
	return kFALSE;
    }


    // read digits from file
    TTree* tr = (TTree*)fDigitsInputFile->Get(Form("Event%d/TreeD",
						   fEventNoInFile));

    if(!tr)
    {
	//AliWarning(Form("digits tree for event %d not found", fEventNoInFile));
	return kFALSE;
    }

    fDigMan->ReadDigits(tr);
    delete tr;

    // expand digits for use in this task
    for(Int_t det=0; det<540; det++)
    {
	if(fDigMan->GetDigits(det))
	{
	    fDigMan->GetDigits(det)->Expand();
	}
    }

    fDigitsLoadedFlag = kTRUE;
    return kTRUE;
}

//________________________________________________________________________
Bool_t Ali_AS_analysis_TRD_digits::WriteDigits()
{
    cout << "In WriteDigits" << endl;
    // check for output file
    if(!fDigitsOutputFile)
    {
	AliError("digits output file not available");
	return kFALSE;
    }

    // compress digits for storage
    for(Int_t det=0; det<540; det++)
    {
	fDigMan->GetDigits(det)->Expand();
    }

    // create directory to store digits tree
    TDirectory* evdir =
	fDigitsOutputFile->mkdir(Form("Event%d", fEventNoInFile),
				 Form("Event%d", fEventNoInFile));

    evdir->Write();
    evdir->cd();

    // save digits tree
    TTree* tr = new TTree("TreeD", "TreeD");
    fDigMan->MakeBranch(tr);
    fDigMan->WriteDigits();
    delete tr;

    return kTRUE;
}


//----------------------------------------------------------------------------------------
void Ali_AS_analysis_TRD_digits::func_tail_cancellation(Short_t *arr, Int_t nexp)
{
    // Tail cancellation by deconvolution for PASA v4 TRF
    //

    Int_t fBaseline = 10;

    Float_t rates[2];
    Float_t coefficients[2];

    // Initialization (coefficient = alpha, rates = lambda)
    Float_t r1 = 1.0;
    Float_t r2 = 1.0;
    Float_t c1 = 0.5;
    Float_t c2 = 0.5;

    r1 = 1.156;
    r2 = 0.130;
    c1 = 0.066;
    c2 = 0.000;

    if (nexp == 1) {   // 1 Exponentials
        r1 = 1.156;
        r2 = 0.130;
        c1 = 0.066;
        c2 = 0.000;
    }
    if (nexp == 2) {   // 2 Exponentials
        //Double_t par[4];
        //fReconstructor->GetRecoParam()->GetTCParams(par);
        //r1 = par[0];//1.156;
        //r2 = par[1];//0.130;
        //c1 = par[2];//0.114;
        //c2 = par[3];//0.624;

        r1 = 1.156;
        r2 = 0.130;
        c1 = 0.114;
        c2 = 0.624;
    }

    coefficients[0] = c1;
    coefficients[1] = c2;

    Double_t dt = 0.1;

    rates[0] = TMath::Exp(-dt/(r1));
    rates[1] = (nexp == 1) ? .0 : TMath::Exp(-dt/(r2));

    Float_t reminder[2] = { .0, .0 };
    Float_t correction = 0.0;
    Float_t result     = 0.0;

    for (int i = 0; i < 24; i++) {

        result = arr[i] - correction - fBaseline;    // No rescaling
        arr[i] = (Short_t)(result + fBaseline + 0.5f);

        correction = 0.0;
        for (int k = 0; k < 2; k++) {
            correction += reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
        }
    }
}
//----------------------------------------------------------------------------------------
