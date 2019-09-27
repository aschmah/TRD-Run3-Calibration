
#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
//#include <iostream.h>
#include "StMessMgr.h"
#include <fstream>
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
#include "TMarker3DBox.h"
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
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"
#include "TChain.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"
//#include <GSLMultiMinimizer.h>
//#include "Math/GSLMultiMinimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "TRotation.h"
#include <algorithm>    // std::sort

// STAR includes
#include "TVector3.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "TLorentzVector.h"
#include "StPhysicalHelixD.hh"
#include <TGeoMatrix.h>


//------------------------------------------------------------------------------------------------------------
static const Float_t Pi = TMath::Pi();
static TRandom ran;
static TString HistName, HistNameB;
static char NoP[50];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
static const Int_t N_v2_vs_pt_BW = 14;
static Double_t mt_m0_cut;
static Double_t pt_cut_BW_global;
static Int_t flag_v2_BW_use[N_v2_vs_pt_BW];
static TGraphAsymmErrors* tgae_v2_stat_BW[N_v2_vs_pt_BW];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
TVector3 Calc_Circle_Center(TVector3 v_pointA, TVector3 v_pointB, TVector3 v_pointC)
{
    TVector3 v_center_of_circle;

    // Calculate center of circle
    Double_t xy_circle[2] = {0.0,0.0};
    Double_t x12  = v_pointA.X() - v_pointB.X();
    Double_t x13  = v_pointA.X() - v_pointC.X();
    Double_t y12  = v_pointA.Y() - v_pointB.Y();
    Double_t y13  = v_pointA.Y() - v_pointC.Y();
    Double_t x12p = v_pointA.X() + v_pointB.X();
    Double_t x13p = v_pointA.X() + v_pointC.X();
    Double_t y12p = v_pointA.Y() + v_pointB.Y();
    Double_t y13p = v_pointA.Y() + v_pointC.Y();

    xy_circle[0] = (y12*(x13*x13p + y13*y13p) - y13*(x12*x12p + y12*y12p))/(2.0*(x13*y12 - x12*y13));
    xy_circle[1] = (x13*(x12*x12p + y12*y12p) - x12*(x13*x13p + y13*y13p))/(2.0*(x13*y12 - x12*y13));

    v_center_of_circle.SetXYZ(xy_circle[0],xy_circle[1],0.0);
    return v_center_of_circle;
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void Convert_NDC_to_User(TPad* iPad, Double_t x_NDC, Double_t y_NDC,
                         Double_t &x_User, Double_t &y_User,
                         Int_t flag_log_x, Int_t flag_log_y)
{
    // Convert NDC coordinates to User coordinates, if x/y axes are in log scale set flag_log to 1
    Double_t xr = iPad->GetX2()-iPad->GetX1();
    Double_t yr = iPad->GetY2()-iPad->GetY1();
    x_User = x_NDC*xr + iPad->GetX1();
    y_User = y_NDC*yr + iPad->GetY1();

    if(flag_log_x == 1) x_User = TMath::Power(10,x_User);
    if(flag_log_y == 1) y_User = TMath::Power(10,y_User);
    //cout << "X1: " << iPad->GetX1() << ", X2: " << iPad->GetX2() << ", Y1: " <<
    //    iPad->GetY1() << ", Y2: " << iPad->GetY2() << ", xr: " << xr << ", yr: " << yr << endl;
}

void Convert_User_to_NDC(TPad* iPad, Double_t x_User, Double_t y_User, Double_t &x_NDC, Double_t &y_NDC)
{
    // Convert User coordinates to NDC coordinates
    Double_t xr = iPad->GetX2()-iPad->GetX1();
    Double_t yr = iPad->GetY2()-iPad->GetY1();
    x_NDC = (x_User-iPad->GetX1())/ xr;
    y_NDC = (y_User-iPad->GetY1())/ xr;
}
//------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

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
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
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
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistLine2(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Float_t x_start, Float_t x_stop, Float_t y_min, Float_t y_max)
{
    TPolyLine* Hist_line = new TPolyLine();
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y >= y_min && y < y_max && y != 0 && x >= x_start && x <= x_stop)
        {
            //cout << "x = " << x << ", y = " << y << endl;
            Hist_line->SetNextPoint(x,y);
        }
    }
    Hist_line -> Draw();
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistErrorBand(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    // Modified March 14th
    Int_t N_points = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y != 0 && x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }
    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line = new TGraph(N_total_points);
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    Hist_line -> SetFillStyle(FillStyle);
    Hist_line -> SetFillColor(FillColor);
    Int_t N_point = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y     = Histo->GetBinContent(i);
        Double_t y_err = Histo->GetBinError(i);
        if(y != 0)
        {
            Double_t x = Histo->GetBinCenter(i);
            if(x >= x_start && x <= x_stop)
            {
                //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
                Hist_line->SetPoint(N_point,x,y-y_err);
                Hist_line->SetPoint(N_total_points-2-N_point,x,y+y_err);
                if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-y_err);
                N_point++;
            }
        }
    }
    Hist_line -> Draw("f");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t One_over_x_FitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    y = 0.0;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    //if(x != 0.0 && !((par0/x) < 0.0 && par2 < 1.0)) y = pow(par0/x,par2) + par1;
    if(x != 0.0) y = par0*TMath::Power(1.0/x,par2) + par1;
    //if(x != 0.0) y =par0/x + par1 + 0.0001*par2;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t TwoGaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + par3*TMath::Gaus(x,par4,par5,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussPolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, pol0, pol1, pol2, pol3, pol4, pol5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    pol0  = par[3];
    pol1  = par[4];
    pol2  = par[5];
    pol3  = par[6];
    pol4  = par[7];
    pol5  = par[8];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t FlowFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t phi, y, v0, v1, v2, v3, v4;
    v0  = par[0];
    v1  = par[1];
    v2  = par[2];
    v3  = par[3];
    v4  = par[4];
    phi = x_val[0];
    y = v0 * (1.0 + 2.0*v1*TMath::Cos(phi) + 2.0*v2*TMath::Cos(2.0*phi)
              + 2.0*v3*TMath::Cos(3.0*phi) + 2.0*v4*TMath::Cos(4.0*phi));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod_x(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Double_t LevyFitFunc_pT(Double_t* x_val, Double_t* par)
{
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2 vs. pT
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT, a, b, c, d, n;
    pT = x_val[0];
    n  = par[0]; // number-of-constituent quarks
    a  = par[1];
    b  = par[2];
    c  = par[3];
    d  = par[4];

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. pT/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT_ncq, a, b, c, d, n;
    pT_ncq = x_val[0];
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT_ncq - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc_Poly(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0, par6, par7;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass
    par6   = par[6];
    par7   = par[7];

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = (a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d) + par6*mT_ncq + par7*mT_ncq*mT_ncq;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_m0, a, b, c, d, n, m0;
    mT_m0  = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_m0 + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncE(Double_t* x_val, Double_t* par)
{
    // Original
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        //for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            //r = r_start + j*delta_r;
            r = 0.01;

            //r = R; //

            rho = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_radial(Double_t* x_val, Double_t* par)
{
    // Blast wave function with radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            r = r_start + j*delta_r;

            //r = R; //

            rho  = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncC(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 100;
    Int_t nbins_r      = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}
    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncD(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t NCQ[14]  = {2.0,3.0,3.0,3.0,3.0,3.0,3.0,2.0,3.0,3.0,2.0,2.0,2.0,2.0};
    pt /= NCQ[PID];
    Mass[PID] /= NCQ[PID];
    //pt *= NCQ[PID];
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 25;
    Int_t nbins_r      = 25;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}

    v2 *= NCQ[PID];
    //v2 /= NCQ[PID];
    //v2 *= (1.0/(1.0 + TMath::Exp((pt-4.0/NCQ[PID])/0.6)));

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void BlastWaveSimultaneous(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
    Int_t npfits      = 0;
    Double_t chi2     = 0.0;
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};

    for(Int_t i = 0; i < 14; i++) // loop over PIDs
    {
        p[5] = i; // PID
        //Double_t pt_cut = TMath::Sqrt((Mass[i]+mt_m0_cut)*(Mass[i]+mt_m0_cut) - Mass[i]*Mass[i]); // mt-m0 cut
        Double_t pt_cut = pt_cut_BW_global;
        if(flag_v2_BW_use[i] == 1)
        {
            for(Int_t ix = 0; ix < tgae_v2_stat_BW[i]->GetN(); ix++)
            {
                Double_t x[] = {tgae_v2_stat_BW[i]->GetX()[ix]};
                if(x[0] < pt_cut)
                {
                    Double_t y      = tgae_v2_stat_BW[i]->GetY()[ix];
                    Double_t ye     = tgae_v2_stat_BW[i]->GetErrorYhigh(ix);
                    //ye = 0.01;
                    //cout << "ix = " << ix << ", x = " << x[0] << ", y = " << y << ", ye = " << ye << endl;
                    Double_t bw_val = BlastWaveFitFunc(x,p);
                    //Double_t bw_val = 0.1;
                    //ye += 0.003;
                    Double_t diff   = (y - bw_val)/ye;
                    chi2 += diff*diff;
                    npfits++;
                }
                else break;
            }
        }
    }

    fval = chi2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Hist_interpolate_and_error(TH1F* hist, Double_t x, Double_t &Int_val, Double_t &Int_err)
{
    // Linear interpolation of a histogram
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[2]   = {0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetEntries() == 0) // no data poins to interpolate
    {
        return_val = -1;
        Int_val = 0;
        Int_err = 0;
        //cout << "No entries in histogram" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t binx = 1; binx < hist->GetNbinsX(); binx++)
        {
            Double_t bin_error_val   = hist->GetBinError(binx);
            Double_t bin_x_pos       = hist->GetBinCenter(binx);
            if(bin_error_val != 0)
            {
                err_counter++;
                bin_entries[1] = hist->GetBinContent(binx);
                bin_error[1]   = hist->GetBinError(binx);
                bin_x_val[1]   = hist->GetBinCenter(binx);
                if(bin_x_pos >= x)
                {
                    binx_high = binx;
                    flag_max  = 1;
                    break;
                }
                else flag_max = 0;
            }
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val = bin_entries[1];
            Int_err = bin_error[1];
            return return_val;
        }
        for(Int_t binx_low = binx_high; binx_low > 0; binx_low--)
        {
            bin_entries[0] = hist->GetBinContent(binx_low);
            bin_error[0]   = hist->GetBinError(binx_low);
            bin_x_val[0]   = hist->GetBinCenter(binx_low);
            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        if(bin_error[0] != 0 && bin_error[1] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;
                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[1];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;
                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val = bin_entries[0];
                Int_err = bin_error[0];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t TGraphAsymmErrors_interpolate_and_error(TGraphAsymmErrors* hist, Double_t x, Double_t &Int_val, Double_t &Int_err_low,Double_t &Int_err_high)
{
    // V2: 03.12.2012 -> bug fixed to calculate error bars
    // Linear interpolation of a TGraphAsymmErrors
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[4]   = {0,0,0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetN() == 0) // no data poins to interpolate
    {
        return_val   = -1;
        Int_val      = 0;
        Int_err_low  = 0;
        Int_err_high = 0;
        //cout << "No entries in TGraphAsymmErrors" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t epoint = 0; epoint < hist->GetN(); epoint++)
        {
            hist->GetPoint(epoint,bin_x_val[1],bin_entries[1]);
            bin_error[2] = hist->GetErrorYlow(epoint);
            bin_error[3] = hist->GetErrorYhigh(epoint);

            err_counter++;
            if(bin_x_val[1] >= x)
            {
                binx_high = epoint;
                flag_max  = 1;
                break;
            }
            else flag_max = 0;
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val      = bin_entries[1];
            Int_err_low  = bin_error[2];
            Int_err_high = bin_error[3];
            return return_val;
        }
        for(Int_t epoint = binx_high; epoint >= 0; epoint--)
        {
            hist->GetPoint(epoint,bin_x_val[0],bin_entries[0]);
            bin_error[0] = hist->GetErrorYlow(epoint);
            bin_error[1] = hist->GetErrorYhigh(epoint);

            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        //cout << "bin0 = " << bin_error[0] << ", bin2 = " << bin_error[2] << endl;

        if(bin_error[0] != 0 && bin_error[2] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;

                //cout << "x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 << ", y2 = " << y2 << endl;

                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[2];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_low = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);

                y1_err = bin_error[1];
                y2_err = bin_error[3];

                termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                termC  = -(x-x2)/(x2-x1);
                termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_high = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val      = bin_entries[0];
                Int_err_low  = bin_error[0];
                Int_err_high = bin_error[1];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Get_muB(Double_t sqrt_sNN, Int_t flag_input, Double_t &muB_err)
{
    // Returns baryon chemical potential as a function of the center-of-mass energy sqrt(sNN)
    // from http://arxiv.org/pdf/1111.2406v1.pdf
    // flag_input: 0 = parametrization, 1 = STAR fits, 2 = corrected parametrization for 0-80% (see below)

    // corrected parametrization: Centrality data taken from: http://arxiv.org/pdf/0808.2041.pdf
    Double_t centrality_vals[9] = {75,65,55,45,35,25,15,7.5,2.5};
    Double_t muB_Energy[2][2][9] = // [62.4,200][values,error][70-80%,...,0-5%]
    {
        {   // Au+Au @ 62.4 GeV
            {37.7,42.5,47.0,51.3,54.2,54.5,59.4,61.0,62.7}, // values
            {6.5,5.8,5.1,5.2,5.2,5.2,5.4,5.7,6.0}  // errors
        },
        {   // Au+Au @ 200 GeV
            {14.1,15.3,17.7,18.9,18.6,21.3,21.0,22.8,21.9}, // values
            {4.2,4.2,4.2,4.2,4.2,4.2,4.2,4.5,4.5}  // errors
        }
    };

    Double_t muB     = 0.0;
    muB_err          = 0.0;
    Double_t a       = 1.482;  // +/- 0.0037 GeV
    Double_t b       = 0.3517; // +/- 0.009 GeV^-1

    if(!(flag_input == 0 || flag_input == 1 || flag_input == 2))
    {
        cout << "WARNING: Get_muB function, flag_input is wrong. Either 0 or 1. Parametrization is used." << endl;
        flag_input = 0;
    }

    if(sqrt_sNN >= 0.0 && (flag_input == 0 || flag_input == 2))
    {
        muB = a/(1.0+b*sqrt_sNN);
    }
    if(flag_input == 2) // correct the muB values for 0-80%
    {
        Double_t scale_fac[2][2];  // [62.4,200][val,error]
        Double_t muB_mean[2][2]; // [62.4,200][val,error]
        for(Int_t ebeam = 0; ebeam < 2; ebeam++)
        {
            for(Int_t val_err = 0; val_err < 2; val_err++) // calculate mean muB + error
            {
                muB_mean[ebeam][val_err] = 0.0; // approximation for 0-80%
                Double_t total_weight = 0.0;
                for(Int_t cent = 0; cent < 9; cent++)
                {
                    Double_t weight = 1.0;
                    if(cent >= 7) weight = 0.5;
                    if(val_err == 0) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent];
                    if(val_err == 1) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent]*weight*muB_Energy[ebeam][val_err][cent];
                    total_weight += weight;
                }
                if(total_weight > 0.0)
                {
                    if(val_err == 0) muB_mean[ebeam][val_err] /= total_weight;
                    if(val_err == 1)
                    {
                        muB_mean[ebeam][val_err] = TMath::Sqrt(muB_mean[ebeam][val_err]);
                        muB_mean[ebeam][val_err] /= total_weight;
                    }
                }
            }
            scale_fac[ebeam][0] = muB_mean[ebeam][0]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0); // <muB>/0-10% = 0-80%/0-10%
            Double_t term1 = muB_mean[ebeam][1]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0);
            Double_t term2 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][8]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            Double_t term3 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][7]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            scale_fac[ebeam][1] = TMath::Sqrt(term1*term1+term2*term2+term3*term3); // error on scaling factor
        }
        Double_t mean_scale_fac[2] = {(scale_fac[0][0]+scale_fac[1][0])/2.0,TMath::Sqrt(scale_fac[0][1]*scale_fac[0][1]+scale_fac[1][1]*scale_fac[1][1])/2.0}; // [value,error]
        Double_t original_muB = muB;
        muB     = muB * mean_scale_fac[0]; // 0-80% estimation
        muB_err = muB * mean_scale_fac[1]; // 0-80% error estimation
        //cout << "muB = " << original_muB << ", muB(0-80%) = " << muB << ", muB_err(0-80%) = " << muB_err << ", <scale> = " << mean_scale_fac[0] << ", <scale_err> = " << mean_scale_fac[1] << endl;
        //cout << "scale_fac(62 GeV) = " << scale_fac[0][0] << ", scale_fac(200 GeV) = " << scale_fac[1][0] << endl;
        //cout << "scale_fac_err(62 GeV) = " << scale_fac[0][1] << ", scale_fac_err(200 GeV) = " << scale_fac[1][1] << endl;
    }

    if(flag_input == 1)
    {
        // Preliminary muB parameters from Lokesh with strange particles for 0-80% Au+Au collisions
        // 7.7 :  3.61377e-01   1.51360e-02
        // 11.5:  2.59978e-01   1.06312e-02
        // 39  :  8.81846e-02   4.73062e-03
        // 200 :  22 -> Thats a guess

        // Fittet in TmuB_fit.cc macro with One_over_x_FitFunc function
        Double_t par0 = 2.13495e+00;
        Double_t par1 = 4.35916e-04;
        Double_t par2 = 8.67894e-01;
        muB = par0*TMath::Power(1.0/sqrt_sNN,par2) + par1;
    }

    return muB;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_new_Symbol(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_Point_new_Symbol(Double_t x_val, Double_t y_val, Double_t x_min_err, Double_t x_max_err,
                                Double_t y_min_err, Double_t y_max_err,
                                Int_t style, Int_t color, Float_t size)
{
    TGraphAsymmErrors* tgae = new TGraphAsymmErrors();
    tgae->SetPoint(0,x_val,y_val);
    tgae->SetPointError(0,x_min_err,x_max_err,y_min_err,y_max_err);

    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_header(FILE* HTML_file)
{
    fprintf(HTML_file,"%s \n","<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">");
    fprintf(HTML_file,"%s \n","<html>");
    fprintf(HTML_file,"%s \n","<head>");


    fprintf(HTML_file,"%s \n","  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">");
    fprintf(HTML_file,"%s \n","  <title>STAR Physics Database</title>");
    fprintf(HTML_file,"%s \n","</head>");
    fprintf(HTML_file,"%s \n"," <body bgcolor=\"white\">");

    fprintf(HTML_file,"%s \n","<h2><span class=\"normaltext\"><font color=\"#ff8c00\" face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\">");
    fprintf(HTML_file,"%s \n","Observation of an energy-dependent difference in elliptic flow between particles and anti-particles in relativistic heavy ion collisions");
    fprintf(HTML_file,"%s \n","</font>");

    fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\"><br> <br> </font></span></h2>");

    //fprintf(HTML_file,"%s \n","<hr><h3><font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">Figure 1");
    //fprintf(HTML_file,"%s \n","</font></h3>");
    //fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">");
    //fprintf(HTML_file,"%s \n","<A HREF=\"fig1.png\"> png </A> | <A HREF=\"fig1.eps\"> eps </A> </br>");
    fprintf(HTML_file,"%s \n","</font>");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_title(TString Label, FILE* HTML_file)
{
    TString label_html = "<hr><h3><font face=\\\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\\\" color=\\\"#0000cd\\\">";
    label_html += Label.Data();
    label_html += " </font> </h3> </hr>";
    //label_html += "\"";
    fprintf(HTML_file,"%s \n",label_html.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_tables(TString labelX, TString labelY, TGraphAsymmErrors* data_stat, TGraphAsymmErrors* syst_errors, TGraphAsymmErrors* glob_syst_errors, TString Label, FILE* HTML_file)
{
    fprintf(HTML_file,"%s %s %s \n","<tr><td>",Label.Data(),"</td> </tr>");
    TString label_html = "<table border=1> <tr><th> ";
    label_html += labelX.Data();
    label_html += " </th> <th> ";
    label_html += labelY.Data();
    label_html += " </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>";

    //fprintf(HTML_file,"%s \n","<table border=1> <tr><th> p<sub>T</sub> (GeV/c) </th> <th> v<sub>2</sub> </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>");

    fprintf(HTML_file,"%s \n",label_html.Data());
    Int_t flag_glob_syst = 1;
    if(glob_syst_errors == syst_errors) flag_glob_syst = -1;

    for(Int_t epoint = 0; epoint < data_stat->GetN(); epoint++)
    {
        Double_t x_val, y_val, x_syst, stat_error, low_syst, high_syst, glob_syst;
        data_stat                     ->GetPoint(epoint,x_val,y_val);
        stat_error = data_stat        ->GetErrorYhigh(epoint);
        high_syst  = syst_errors      ->GetErrorYhigh(epoint);
        low_syst   = syst_errors      ->GetErrorYlow(epoint);
        glob_syst  = glob_syst_errors ->GetErrorYhigh(0);
        if(flag_glob_syst == -1) glob_syst = 0.0;

        fprintf(HTML_file,"%s %f %s %f %s %f %s %f %s %f %s %f %s \n","<tr> <td> ",x_val," </td> <td> ",y_val," </td> <td> ",stat_error," </td> <td> ",
                low_syst," </td> <td> ",high_syst," </td> <td> ",glob_syst,"</td></tr>");

    }
    fprintf(HTML_file,"%s \n","</table><hr><table border=1>");
}
//----------------------------------------------------------------------------------------



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
TVector3 calculateDCA_vec_StraightToPoint(TVector3 &base, TVector3 &dir, TVector3 &point)
{
  // calculates the minimum distance vector of a point to a straight given as parametric straight x = base + n * dir

    TVector3 diff = base-point;
    TVector3 dir_norm = dir;
    dir_norm *= (1.0/dir.Mag());
    Double_t proj_val = diff.Dot(dir_norm);
    TVector3 proj_dir = dir_norm;
    proj_dir *= proj_val;

    TVector3 dist_vec = proj_dir - diff;

    return dist_vec;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t get_HFT_det_index(Int_t sensor_id)
{
    // Determines detector index (see definition below) based on the sonsor id
    // 0 = inner pixel left
    // 1 = outer pixel left
    // 2 = IST left
    // 3 = inner pixel right
    // 4 = outer pixel right
    // 5 = IST right

    Int_t HFT_det_index = -1;
    if(
       (sensor_id >= 1   && sensor_id <= 10)  ||
       (sensor_id >= 41  && sensor_id <= 50)  ||
       (sensor_id >= 81  && sensor_id <= 90)  ||
       (sensor_id >= 121 && sensor_id <= 130) ||
       (sensor_id >= 161 && sensor_id <= 170)
      )
    {
        HFT_det_index = 0;
    }

    if(
       (sensor_id >= 11   && sensor_id <= 40)  ||
       (sensor_id >= 51   && sensor_id <= 80)  ||
       (sensor_id >= 91   && sensor_id <= 120) ||
       (sensor_id >= 131  && sensor_id <= 160) ||
       (sensor_id >= 171  && sensor_id <= 200)
      )
    {
        HFT_det_index = 1;
    }
    if(
       sensor_id >= 1001 && sensor_id <= 1072
      )
    {
        HFT_det_index = 2;
    }
    if(
       (sensor_id >= 1+200   && sensor_id <= 10+200)  ||
       (sensor_id >= 41+200  && sensor_id <= 50+200)  ||
       (sensor_id >= 81+200  && sensor_id <= 90+200)  ||
       (sensor_id >= 121+200 && sensor_id <= 130+200) ||
       (sensor_id >= 161+200 && sensor_id <= 170+200)
      )
    {
        HFT_det_index = 3;
    }
    if(
       (sensor_id >= 11+200   && sensor_id <= 40+200)  ||
       (sensor_id >= 51+200   && sensor_id <= 80+200)  ||
       (sensor_id >= 91+200   && sensor_id <= 120+200) ||
       (sensor_id >= 131+200  && sensor_id <= 160+200) ||
       (sensor_id >= 171+200  && sensor_id <= 200+200)
      )
    {
        HFT_det_index = 4;
    }
    if(
       sensor_id >= 1073 && sensor_id <= 1145
      )
    {
        HFT_det_index = 5;
    }
    return HFT_det_index;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_histogram(TH1D* hist, TString plot_option, TString label_x, TString label_y, Double_t x_start, Double_t x_stop,
                    Double_t y_start, Double_t y_stop, Double_t title_offset_x, Double_t title_offset_y, Double_t label_offset_x,
                    Double_t label_offset_y, Double_t label_size, Double_t title_size, Int_t N_div)
{
    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(title_offset_x);
    hist->GetYaxis()->SetTitleOffset(title_offset_y);
    hist->GetXaxis()->SetLabelOffset(label_offset_x);
    hist->GetYaxis()->SetLabelOffset(label_offset_y);
    hist->GetXaxis()->SetLabelSize(label_size);
    hist->GetYaxis()->SetLabelSize(label_size);
    hist->GetXaxis()->SetTitleSize(title_size);
    hist->GetYaxis()->SetTitleSize(title_size);
    hist->GetXaxis()->SetNdivisions(N_div,'N');
    hist->GetYaxis()->SetNdivisions(N_div,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitle(label_x.Data());
    hist->GetYaxis()->SetTitle(label_y.Data());
    if(!(x_start == 0 && x_stop == 0)) hist->GetXaxis()->SetRangeUser(x_start,x_stop);
    if(!(y_start == 0 && y_stop == 0)) hist->GetYaxis()->SetRangeUser(y_start,y_stop);
    hist->DrawCopy(plot_option.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t get_HFT_pixel_sector(Int_t sensor_id, Int_t &IST_ladder)
{
    // Determines pixel sector based on the sonsor id

    // PXL
    Int_t HFT_pxl_sector = -1;
    if(sensor_id >= 1 && sensor_id <= 400)
    {
        HFT_pxl_sector = (Int_t)((sensor_id-1)/40);
    }

    // IST -> 24 ladders, 6 sensors per ladder -> 144 sensors in total, starting from sensor_id = 1001
    if(sensor_id >= 1001 && sensor_id <= 1144)
    {
        IST_ladder = (Int_t)((sensor_id-1001)/6);
    }
    else IST_ladder = -1;

    return HFT_pxl_sector;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_STAR_3D()
{

    TGeoManager *geom = new TGeoManager("geom","My 3D Project");
    //------------------Creat materials------------------------------
    TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
    TGeoMaterial *Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
    Fe->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_outer_tube = new TGeoMaterial("M_outer_tube",55.84,26.7,7.87);
    M_outer_tube->SetTransparency(93); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_IDS = new TGeoMaterial("M_IDS",55.84,26.7,7.87);
    M_IDS       ->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_beampipe = new TGeoMaterial("M_beampipe",55.84,26.7,7.87);
    M_beampipe       ->SetTransparency(70); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_Pixel_support = new TGeoMaterial("M_Pixel_support",55.84,26.7,7.87);
    M_Pixel_support    ->SetTransparency(70); // higher value means more transparent, 100 is maximum


    //------------------Create media---------------------------------
    TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);
    TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);
    TGeoMedium *Me_outer_tube = new TGeoMedium("Me_outer_tube",1,M_outer_tube);
    TGeoMedium *Me_IDS        = new TGeoMedium("Me_IDS",1,M_IDS);
    TGeoMedium *Me_beampipe   = new TGeoMedium("Me_beampipe",1,M_beampipe);
    TGeoMedium *Me_Pixel_support   = new TGeoMedium("Me_Pixel_support",1,M_Pixel_support);

    //------------------Create TOP volume----------------------------
    TGeoVolume *top = geom->MakeBox("top",Air,500,500,500);
    geom->SetTopVolume(top);
    geom->SetTopVisible(0);
    // If you want to see the boundary, please input the number, 1 instead of 0.
    // Like this, geom->SetTopVisible(1);


    TGeoVolume *inner_field_tube       = geom->MakeTube("inner_field_tube",Iron,49.5,50.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *outer_field_tube       = geom->MakeTube("outer_field_tube",Me_outer_tube,199.5,200.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_central_part       = geom->MakeTube("IDS_central_part",Me_IDS,42.8/2.0,43.0/2.0,56.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_side_parts         = geom->MakeTube("IDS_side_parts",Me_IDS,79.3/2.0,79.5/2.0,(222.7-64.0)/2.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_connection_parts_R = geom->MakeCone("IDS_connection_parts_R",Me_IDS,(64.0-56.0)/2.0,42.8/2.0,43.0/2.0,79.3/2.0,79.5/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *IDS_connection_parts_L = geom->MakeCone("IDS_connection_parts_L",Me_IDS,(64.0-56.0)/2.0,79.3/2.0,79.5/2.0,42.8/2.0,43.0/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *beampipe_central_part       = geom->MakeTube("beampipe_central_part",Me_beampipe,4.05/2.0,4.15/2.0,141.5);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_side_parts         = geom->MakeTube("beampipe_side_parts",Me_beampipe,9.52/2.0,9.62/2.0,100.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_connection_parts_R = geom->MakeCone("beampipe_connection_parts_R",Me_beampipe,(191.5-141.5)/2.0,4.05/2.0,4.15/2.0,9.52/2.0,9.62/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *beampipe_connection_parts_L = geom->MakeCone("beampipe_connection_parts_L",Me_beampipe,(191.5-141.5)/2.0,9.52/2.0,9.62/2.0,4.05/2.0,4.15/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *Pixel_support       = geom->MakeTube("Pixel_support",Me_Pixel_support,21.8/2.0,22.0/2.0,56.0);  // r_min, r_max, dz (half of total length)

    inner_field_tube       ->SetLineColor(4);
    outer_field_tube       ->SetLineColor(kRed-8);
    IDS_central_part       ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_side_parts         ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_R ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_L ->SetLineColor(2);  // Inner Detector Support (IDS)

    beampipe_central_part       ->SetLineColor(3);  // (beampipe)
    beampipe_side_parts         ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_R ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_L ->SetLineColor(3);  // (beampipe)

    Pixel_support ->SetLineColor(kYellow-3);  // (pixel support)

    top->AddNodeOverlap(inner_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,64.0+(222.7-64.0)/2.0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,-(64.0+(222.7-64.0)/2.0)));
    top->AddNodeOverlap(IDS_connection_parts_R,1,new TGeoTranslation(0,0,56.0+(64.0-56.0)/2.0));
    top->AddNodeOverlap(IDS_connection_parts_L,1,new TGeoTranslation(0,0,-(56.0+(64.0-56.0)/2.0)));

    top->AddNodeOverlap(beampipe_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,191.5+100.0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,-(191.5+100.0)));
    top->AddNodeOverlap(beampipe_connection_parts_R,1,new TGeoTranslation(0,0,141.4+(191.5-141.5)/2.0));
    top->AddNodeOverlap(beampipe_connection_parts_L,1,new TGeoTranslation(0,0,-(141.4+(191.5-141.5)/2.0)));

    top->AddNodeOverlap(Pixel_support,1,new TGeoTranslation(0,0,0));

    top->DrawClone("ogl");



    const Int_t n_TPC_points = 50;
    TPolyLine3D   *TPC_endcaps[4];
    TPolyLine3D   *TPC_tube[4];
    TPolyLine3D   *TPC_tube_lines[n_TPC_points+1];

    Float_t radius_table[4] = {200,200,3.81,3.81};
    Float_t z_val_table[4]  = {200,-200,200,-200};

    Float_t radius_table_tube[4] = {50,50,50,50};
    Float_t z_val_table_tube[4]  = {200,-200,100,-100};

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_endcaps[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_endcaps[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        TPC_endcaps[r]->SetLineStyle(0);
        TPC_endcaps[r]->SetLineColor(28); // 28
        TPC_endcaps[r]->SetLineWidth(2);
        TPC_endcaps[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_tube[r] = new TPolyLine3D();
        Float_t radius   = radius_table_tube[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table_tube[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_tube[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
            if(r == 0 && (t%4 == 0))
            {
                TPC_tube_lines[t] = new TPolyLine3D();
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_val_table_tube[r+1]);
                TPC_tube_lines[t]->SetLineStyle(0);
                TPC_tube_lines[t]->SetLineColor(28); // 28
                TPC_tube_lines[t]->SetLineWidth(1);
                //TPC_tube_lines[t]->DrawClone("ogl");
            }
        }
        TPC_tube[r]->SetLineStyle(0);
        TPC_tube[r]->SetLineColor(28); // 28
        TPC_tube[r]->SetLineWidth(2);
        TPC_tube[r]->DrawClone("ogl");
    }

    TPolyLine3D   *BeamLine;
    BeamLine       = new TPolyLine3D(2);
    BeamLine   ->SetPoint(0,0,0,-550);
    BeamLine   ->SetPoint(1,0,0,550);
    BeamLine   ->SetLineStyle(0);
    BeamLine   ->SetLineColor(4);
    BeamLine   ->SetLineWidth(2);
    BeamLine   ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Helix_3D(StPhysicalHelixD helix, Int_t line_color, Int_t line_style, Float_t line_width,
                   Float_t helix_min_radius, Float_t helix_max_radius,const Int_t n_helix_points, Float_t step_size)
{
    TPolyLine3D *helix_line;
    helix_line = new TPolyLine3D();

    StThreeVectorF helix_point;
    for(Int_t t = 0; t < n_helix_points; t++)
    {
        Double_t path_length = t*step_size;
        helix_point = helix.at(path_length);
        Double_t x_val = helix_point.x();
        Double_t y_val = helix_point.y();
        Double_t z_val = helix_point.z();
        Double_t radius = TMath::Sqrt(x_val*x_val + y_val*y_val);

        if(radius > helix_min_radius && radius < helix_max_radius)
        {
            helix_line->SetNextPoint(x_val,y_val,z_val);
        }
    }
    helix_line->SetLineStyle(line_style);
    helix_line->SetLineColor(line_color);
    helix_line->SetLineWidth(line_width);
    helix_line->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
StThreeVectorF Get_Helix_point_at_radius(StPhysicalHelixD helix, Float_t helix_max_radius, Float_t &path_length)
{
    StThreeVectorF helix_point;
    helix_point.set(0.0,0.0,0.0);

    const Int_t n_helix_points = 4000;
    Float_t step_size = 0.1;

    for(Int_t t = 0; t < n_helix_points; t++)
    {
        path_length = t*step_size;
        helix_point = helix.at(path_length);
        Double_t x_val = helix_point.x();
        Double_t y_val = helix_point.y();
        Double_t z_val = helix_point.z();
        Double_t radius = TMath::Sqrt(x_val*x_val + y_val*y_val);

        if(radius > helix_max_radius) break;
    }

    return helix_point;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_pt_tower(StPhysicalHelixD helix, Float_t pt, Float_t tower_radius, Color_t fill_color)
{
    Float_t path_length;
    StThreeVectorF dir_vec, pos_vec;
    pos_vec = Get_Helix_point_at_radius(helix,tower_radius,path_length);
    dir_vec = helix.cat(path_length);
    dir_vec /= dir_vec.mag();
    dir_vec *= 0.7*pt;

    TMarker3DBox pTowerHit;

    pTowerHit.SetDirection(dir_vec.theta()*TMath::RadToDeg(),dir_vec.phi()*TMath::RadToDeg());
    pTowerHit.SetPosition(pos_vec.x()+dir_vec.x(),pos_vec.y()+dir_vec.y(),pos_vec.z()+dir_vec.z());
    pTowerHit.SetSize(2.0,2.0,pt);
    pTowerHit.SetFillColor(fill_color);
    pTowerHit.SetLineColor(fill_color);
    pTowerHit.DrawClone("ogl");

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fHelixAtoPointdca(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB)
{
    // V1.1
    Float_t pA[2] = {0.0,-100.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    StThreeVectorF testA;
    for(Int_t r = 0; r < 2; r++)
    {
        testA     = helixA.at(pA[r]); // 3D-vector of helixA at path pA[r]
        distarray[r] = (testA-space_vec).mag(); // dca between helixA and helixB
        //cout << "r = " << r << ", dist = " << distarray[r] << endl;
    }

    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.005 && loopcounter < 100) // stops when the length is too small
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
            testA     = helixA.at(pA[0]); // 3D-vector of helixA at path pA[0]
            distarray[0] = (testA-space_vec).mag(); // new dca
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
            testA     =  helixA.at(pA[1]); // 3D-vector of helixA at path pA[0]
            distarray[1] = (testA-space_vec).mag();
            flip = -1.0;
        }
        //cout << "pA[0] = " << pA[0] << endl;
        loopcounter++;
    }
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
    //cout << "pathA = " << pathA << ", dcaAB = " << dcaAB << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t calcDeterminant(StThreeVectorF& v1,StThreeVectorF& v2,StThreeVectorF& v3)
{
  // calculating the Determinant of a 3 x 3 Matrix 
  // with the column vectors [v1, v2, v3]
  // using the RULE of SARRUS
  //
  // | v1(0)   v2(0)   v3(0) |      | v1(0) v2(0) v3(0)|v1(0) v2(0)  .
  // |                       |      |  \\     \\     X   |  /     /    . 
  // |                       |      |   \\     \\   / \\  | /     /     . 
  // |                       |      |    \\     \\ /   \\ |/     /      . 
  // |                       |      |     \\     X     \\/     /       . 
  // |                       |      |      \\   / \\    /\\    /        .  
  // |                       |      |       \\ /   \\  / |\\  /         . 
  // | v1(1)   v2(1)   v3(1) |   =  | v1(1) v2(1) v3(1)|v1(1) v2(1)  .
  // |                       |      |       / \\    /\\  | /\\          . 
  // |                       |      |      /   \\  /  \\ |/  \\         . 
  // |                       |      |     /     \\/    \\/    \\        . 
  // |                       |      |    /      /\\    /\\     \\       . 
  // |                       |      |   /      /  \\  / |\\     \\      .  
  // |                       |      |  /      /    \\/  | \\     \\     . 
  // | v1(2)   v2(2)   v3(2) |      | v1(2) v2(2) v3(2)| v1(2) v2(2) .  
  //                                 /      /     /  \\     \\     \\   .
  //                                                                
  //                                -      -     -    +     +     +  .

  return ( v1(0) * v2(1) * v3(2) 
	   + v2(0) * v3(1) * v1(2) 
	   + v3(0) * v1(1) * v2(2) 
	   - v3(0) * v2(1) * v1(2) 
	   - v1(0) * v3(1) * v2(2) 
	   - v2(0) * v1(1) * v3(2)); 
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
StThreeVectorF calculatePointOfClosestApproach(StThreeVectorF &base1, StThreeVectorF &dir1,
								    StThreeVectorF &base2, StThreeVectorF &dir2)
{
  //  calculating point of closest approach
  //        
  //        from the equations of the straight lines of g and h 
  //        g: x1 = base1 + l * dir1 
  //        h: x2 = base2 + m * dir2 
  //        
  //        you can construct the following planes:
  //        
  //        E1: e1 = base1  +  a * dir1  +  b * (dir1 x dir2)
  //        E2: e2 = base2  +  s * dir2  +  t * (dir1 x dir2)
  //        
  //        now the intersection point of E1 with g2 = {P1} 
  //        and the intersection point of E2 with g1 = {P2}
  //        
  //        form the base points of the perpendicular to both straight lines.
  //        
  //        The point of closest approach is the middle point between P1 and P2: 
  //        
  //        vertex = (p2 - p1)/2
  // 
  //        E1 ^ g2:
  //
  //           e1 = x2
  //    -->    base1  +  a * dir1  +  b * (dir1 x dir2) = base2 + m * dir2 
  //    -->    base1 - base2 = m * dir2  -  a * dir1  -  b * (dir1 x dir2)       
  //                                          (m)
  //    -->    [ dir2, -dir1, -(dir1 x dir2)] (a) = base1 - base2        
  //                                          (b)
  //           
  //           using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D12 = det [dir2, -dir1, -(dir1 x dir2)] 
  //               = det [dir2,  dir1,  (dir1 x dir2)]
  //           
  //           Dm  = det [base1 - base2, -dir1, -(dir1 x dir2)]
  //               = det [base1 - base2,  dir1,  (dir1 x dir2)]
  //  
  //            m  = Dm/D12
  //           
  //           P1: p1 = x2(m)
  //                  = base2 + Dm/D12 * dir2
  //
  //        E2 ^ g1:
  //
  //           e2 = x1
  //    -->    base2  +  s * dir2  +  t * (dir1 x dir2) = base1 + l * dir1 
  //    -->    base2 - base1 = l * dir1  -  s * dir2  -  t * (dir1 x dir2)       
  //                                          (l)
  //    -->    [ dir1, -dir2, -(dir1 x dir2)] (s) = base2 - base1        
  //                                          (t)
  //           
  //           again using CRAMER's RULE you can find the solution for m (a,b, not used)
  //           
  //           using the rules for converting determinants:
  //           
  //           D21 =  det [dir1, -dir2, -(dir1 x dir2)] 
  //               =  det [dir1,  dir2,  (dir1 x dir2)]
  //               = -det [dir2,  dir1,  (dir1 x dir2)]
  //               = -D12
  //           
  //           Dl  =  det [base2 - base1, -dir2, -(dir1 x dir2)]
  //               =  det [base2 - base1,  dir1,  (dir1 x dir2)]
  //               = -det [base1 - base2,  dir1,  (dir1 x dir2)]
  //
  //            l  =   Dl/D21
  //               = - Dl/D12
  //           
  //           P2: p2 = x1(m)
  //                  = base1 - Dl/D12 * dir1
  //           
  //           
  //           vertex = p1 + 1/2 * (p2 - p1)
  //                  = 1/2 * (p2 + p1)
  //                  = 1/2 *( (base1 + base2) +  1/D12 * ( Dm * dir2 - Dl * dir1) )
  //                      

  StThreeVectorF cross = dir1.cross(dir2); // cross product: dir1 x dir2

  // straight lines are either skew or have a cross point
	      
  StThreeVectorF diff = base1;
  diff-=base2; // Difference of two base vectors base1 - base2
		
  Double_t D;
  D =  calcDeterminant(dir2, dir1 ,cross);

  if (!(fabs(D) > 0.))
    {
      ::Warning(":calculatePointOfClosestApproach","Dirs and cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");
      return StThreeVectorF(-20000.,-20000.,-20000.);
    }

  Double_t Dm =  calcDeterminant(diff , dir1, cross);
  Double_t Dl = -calcDeterminant(diff , dir2, cross);

  StThreeVectorF vertex;
  StThreeVectorF dm;
  StThreeVectorF dl;

  dm = dir2;
  dm *= Dm;

  dl = dir1;
  dl *= Dl;

  vertex = dm - dl;

  vertex *= ((1.)/D);

  vertex+=base1;
  vertex+=base2;
  vertex*=0.5;

  return StThreeVectorF(vertex);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
StThreeVectorF calculateCrossPoint(StThreeVectorF &base1, StThreeVectorF &dir1,
							StThreeVectorF &base2, StThreeVectorF &dir2)
{ 
  // calculating cross point 
  // taking all three equations into account solving the overdetermined set of lin. equations
  // of 
  // base1 + l * dir2 =  base1 + m * dir2 
  //
  // set of lin. equations:
  //  
  //   base1(0) + l * dir1(0) = base2(0) + m * dir2(0) 
  //   base1(1) + l * dir1(1) = base2(1) + m * dir2(1)
  //   base1(2) + l * dir1(2) = base2(2) + m * dir2(2) this line is ignored
  //
  //   written in matrix form
  //
  //        l
  //   M * |   | = base2 - base1
  //       \\ m /
  //
  //   M is a 3x2 matrix
  //     
  // to solve multiply the equation by the transposed Matrix of M from the left: M 
  //     
  //  T      /  l \\                                                               .
  // M * M * |    | = M  * (base2 - base1)
  //         \\ -m /
  // MIND THE '-' of m
  //
  //     / dir1(0) dir2(0) \\                                                      .
  //     |                 |    T   / dir1(0) dir1(1) dir1(2) \\                   .
  // M = | dir1(1) dir2(1) |,  M  = |                         |
  //     |                 |        \\ dir2(0) dir2(1) dir2(2) /                   .
  //     \\ dir1(2) dir2(2) /                                    
  //
  //  T      / (dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2))   (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))  \\ .
  // M * M = |                                                                                                                |
  //         \\ (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))   (dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2))  /                        
  //
  //  T       / d1d1 d1d2 \\                           .
  // M  * M = |           |
  //          \\ d1d2 d2d2 /
  //
  // diff = base2 - base1
  //
  //  T           /  (dir1(0)*diff(0) + dir1(1)*diff(1) + dir1(2)*diff(2)) \\         .
  // M  * diff =  |                                                        |
  //              \\  (dir2(0)*diff(0) + dir2(1)*diff(1) + dir2(2)*diff(2)) /
  //
  //  T           /  d1diff  \\                                          .
  // M  * diff =  |          |
  //              \\  d2diff  /
  // 
  // now the new Matrix set is to be solved by CRAMER'S Rule:
  // 
  // / d1d1 d1d2 \\   /  l \\   /  d1diff \\                   .
  // |           | * |    | = |          |
  // \\ d1d2 d2d2 /   \\ -m /   \\  d2diff /
  //
  //     | d1d1 d1d2 |
  // D = |           | = d1d1*d2d2 - d1d2*d1d2;
  //     | d1d2 d2d2 |
  // 
  //     | d1diff d1d2 |
  // Dl= |              | = d1diff*d2d2 - d1d2*d2diff;
  //     | d2diff d2d2 |              
  //
  // l = Dl/D = l_cross
  // 
  // vertex = base1 + l_cross * dir1
  //

  Double_t d1d1 = dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2);
  Double_t d2d2 = dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2);
  Double_t d1d2 = dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2);
  
  Double_t D = d1d1*d2d2 - (d1d2*d1d2);
  
  if (!(fabs(D) > 0.))
    {
      ::Warning("calculateCrossPoint","Error while calculating cross point ... eqns are lin. dependent:returning default Vertex (-20000,-20000,-20000)");
      return StThreeVectorF(-20000.,-20000.,-20000.);
    }

  Double_t d1diff = dir1(0)*(base2(0)-base1(0))+dir1(1)*(base2(1)-base1(1))+dir1(2)*(base2(2)-base1(2));
  Double_t d2diff = dir2(0)*(base2(0)-base1(0))+dir2(1)*(base2(1)-base1(1))+dir2(2)*(base2(2)-base1(2));

  Double_t Dlambda = d1diff*d2d2-d1d2*d2diff;
  
  Double_t lambda = Dlambda/D;
  
  StThreeVectorF vertex;
  vertex += dir1;
  vertex *= lambda;
  vertex += base1;

  //cout << "Cross point calculated" << endl;
  return StThreeVectorF(vertex);

 // return StThreeVectorF(-20000.,-20000.,-20000.);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t calculateMinimumDistanceStraightToPoint(StThreeVectorF &base, StThreeVectorF &dir,
									 StThreeVectorF &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.mag()>0))
    {
      return -1000000.;
    }
  
  StThreeVectorF diff = base-point;

  StThreeVectorF cross = dir.cross(diff);
  
  return cross.mag()/dir.mag();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
StThreeVectorF calculateDCA_vec_StraightToPoint(StThreeVectorF &base, StThreeVectorF &dir,
									 StThreeVectorF &point)
{
  // calculates the minimum distance vector of a point to a straight given as parametric straight x = base + n * dir

    StThreeVectorF diff = base-point;
    Double_t proj_val = diff.dot(dir/dir.mag());
    StThreeVectorF proj_dir = proj_val*dir/dir.mag();

    StThreeVectorF dist_vec = proj_dir - diff;

    return dist_vec;
}
//----------------------------------------------------------------------------------------



