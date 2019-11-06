

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

//------------------------
#include "TLorentzVector.h"
#include "TSystem.h"
//------------------------

#include "TClonesArray.h"
#include "TGeoMatrix.h"

#include "TObjString.h"

#include "TPolyMarker.h"


//------------------------
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
//------------------------


#include <iomanip>

#include "TMath.h"
#include "TObject.h"




using namespace std;

#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGTripleSlider.h"
#include <TGSlider.h>

#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
//#include <iostream.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <Math/Vector3D.h>
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"
#include "TList.h"
#include "TChain.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"
//#include <GSLMultiMinimizer.h>
//#include "Math/GSLMultiMinimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "TRotation.h"
#include "TSVDUnfold.h"
#include "TSystemDirectory.h"
#include "TMatrixT.h"
#include "TVector2.h"
#include "TEllipse.h"
#include "TPaveLabel.h"
#include "TBox.h"
#include "TSystemDirectory.h"
#include "TAttImage.h"
#include "TImage.h"
#include "TPaletteAxis.h"
#include "THelix.h"
#include "TView.h"
#include "TSystem.h"
#include "TASImage.h"
#include "TImage.h"
#include "TArrayD.h"
//#include "TPython.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TNtuple.h"
#include "TDatime.h"
#include "TArrow.h"

#include "AliTRDgeometry.h"

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>


//#include "../Ali_AS_Event.h"
//#include "../Ali_AS_EventLinkDef.h"

//ClassImp(Ali_AS_TRD_digit)
//ClassImp(Ali_AS_Track)
//ClassImp(Ali_AS_Event)


//static Ali_AS_Event* AS_Event;
//static Ali_AS_Track* AS_Track;
//static Ali_AS_TRD_digit* AS_Digit;

//----------------------------------------------------------------------------------------
//static TChain* input_SE;
//static TString JPsi_TREE   = "Tree_AS_Event";
//static TString JPsi_BRANCH = "Tree_AS_Event_branch";
//static Long64_t file_entries_total;

/*
// Track information
static Float_t        nsigma_TPC_e;
static Float_t        nsigma_TPC_pi;
static Float_t        nsigma_TPC_K;
static Float_t        nsigma_TPC_p;
static Float_t        nsigma_TOF_e;
static Float_t        nsigma_TOF_pi;
static Float_t        nsigma_TOF_K;
static Float_t        nsigma_TOF_p;
static Float_t        TRD_signal;
static Float_t        TRDsumADC;
static Float_t        dca;
static TLorentzVector TLV_part;
static UShort_t       NTPCcls;
static UShort_t       NTRDcls;
static UShort_t       NITScls;
static Float_t        TPCchi2;
static Float_t        TPCdEdx;
static Float_t        TOFsignal;
static Float_t        Track_length;
static Float_t        momentum;
static Float_t        theta_track;
static Float_t        eta_track;
static Float_t        pT_track;
static Int_t          charge_bin;
static Double_t       mass_sq;
*/
//----------------------------------------------------------------------------------------

/*
//----------------------------------------------------------------------------------------
void Init_tree(TString SEList)
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
    cout << "Total number of events in tree: " << file_entries_total << endl;
}
//----------------------------------------------------------------------------------------
*/

//----------------------------THINGS FOR 3D LINE FIT--------------------------------------

// define the parameteric line equation 
void line(Double_t t, Double_t *p, Double_t &x, Double_t &y, Double_t &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
   x = p[0] + p[1]*t; 
   y = p[2] + p[3]*t;
   z = t; 
}

// calculate distance line-point 
double distance2(Double_t x,double_t y,Double_t z, Double_t *p) {
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
   //XYZVector xp(x,y,z);
   //XYZVector x0(p[0], p[2], 0. );
   //XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
   TVector3 xp(x,y,z);
   TVector3 x0(p[0], p[2], 0. ); 
   TVector3 x1(p[0] + p[1], p[2] + p[3], 1. );
   TVector3 u = (x1-x0).Unit(); 
   Double_t d2 = ((xp-x0).Cross(u)) .Mag2();
   return d2; 
}


// function to be minimized 
void SumDistance2(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t ) {
   // the TGraph must be a global variable

   bool first = true; 
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   Double_t * x = gr->GetX();
   Double_t * y = gr->GetY();
   Double_t * z = gr->GetZ();
   Int_t npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) { 
      Double_t d = distance2(x[i],y[i],z[i],par);
      sum += d;
#ifdef DEBUG
      if (first) std::cout << "point " << i << "\t" 
                           << x[i] << "\t" 
                           << y[i] << "\t" 
                           << z[i] << "\t" 
                           << std::sqrt(d) << std::endl; 
#endif
   }
   if (first) 
      //std::cout << "Total sum2 = " << sum << std::endl;
   first = false;
}


//----------------------------------------------------------------------------------------

