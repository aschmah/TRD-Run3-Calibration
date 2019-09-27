
#include "/global/homes/a/aschmah/STAR/Utils/functions.h"

static TChain* NT_HFT_left_right;
static TChain* NT_PXL_IST;
static Float_t x1_HFT,y1_HFT,z1_HFT,x2_HFT,y2_HFT,z2_HFT,x3_HFT,y3_HFT,z3_HFT,x4_HFT,y4_HFT,z4_HFT;
static Float_t PXL_IA_x,PXL_IA_y,PXL_IA_z,PXL_IB_x,PXL_IB_y,PXL_IB_z,PXL_OA_x,PXL_OA_y,PXL_OA_z,PXL_OB_x,PXL_OB_y,PXL_OB_z,PXL_id_IA,PXL_id_IB,PXL_id_OA,PXL_id_OB,IST_x,IST_y,IST_z,IST_id,IST_counter;
static Float_t IST_hit;
static Float_t magfac,track_qA,track_pxA,track_pyA,track_pzA,track_oxA,track_oyA,track_ozA,track_qB,track_pxB,track_pyB,track_pzB,track_oxB,track_oyB,track_ozB;
static Float_t secA,secB,NsecA,NsecB,NnohitsecA,NnohitsecB;
static Float_t idx_OA,idx_OB,idx_IA,idx_IB,idx_IST;

static Float_t id_IL_HFT,id_OL_HFT,id_IR_HFT,id_OR_HFT;
static std::vector< std::vector<TVector3> > HFT_align_data_vectors;
static std::vector< std::vector<TVector3> > PXL_sec_align_data_vectors;
static std::vector< std::vector<TVector3> > IST_align_data_vectors;
static std::vector< std::vector<StPhysicalHelixD> > TPC_align_data_vectors;
static std::vector< std::vector<Double_t> > vec_magfac;
static std::vector< std::vector<Double_t> > vec_charge;

static Int_t flag_left_right_align = 0; // 0 means we align right to left (left is fixed), 1 means the opposite

static Int_t TPC_sector_align_active = -1;

static TRotation rot_HFT_to_TPC;
static TVector3  shift_vec_HFT_to_TPC;

static TRotation rot_HFT_to_TPC_load;
static TVector3  shift_vec_HFT_to_TPC_load;

static TRotation rot_HFT_to_TPC_loadB;
static TVector3  shift_vec_HFT_to_TPC_loadB;

static const Int_t N_TPC_sectors = 24;

static TRotation rot_TPC_sector[N_TPC_sectors];
static TVector3  shift_vec_TPC_sector[N_TPC_sectors];

static TRotation rot_TPC_sector_load[N_TPC_sectors];
static TVector3  shift_vec_TPC_sector_load[N_TPC_sectors];

static TFile* Inputfile_TPC_sector_alignment[N_TPC_sectors];

static TRotation rot_TPC_sector_load_Yuri_old_new[N_TPC_sectors];
static TVector3  shift_vec_TPC_sector_load_Yuri_old_new[N_TPC_sectors];
static TFile* Inputfile_TPC_sector_alignment_Yuri_old_new;

//------------------------------------------------------------
Double_t Align_HFT(const Double_t *align_params)
{
    // Function for left/right half alignment
    Double_t delta_x   = align_params[0];
    Double_t delta_y   = align_params[1];
    Double_t delta_z   = align_params[2];
    Double_t rot_alpha = align_params[3];
    Double_t rot_beta  = align_params[4];
    Double_t rot_gamma = align_params[5];

    Double_t total_dca = 0.0;

    TVector3 HFT_hit_IL_align, HFT_hit_OL_align, HFT_hit_IR_align, HFT_hit_OR_align, HFT_dir_L_align, HFT_dir_R_align,
        dca_vec_IR_align, dca_vec_OR_align, dca_vec_IL_align, dca_vec_OL_align;

    // define and set translation vector
    TVector3 shift_vec;
    shift_vec.SetXYZ(delta_x,delta_y,delta_z);

    for(Long64_t counter = 0; counter < HFT_align_data_vectors[0].size(); counter++)
    {
        // set HFT hits for inner/outer left/right
        HFT_hit_IL_align = HFT_align_data_vectors[0][counter];
        HFT_hit_OL_align = HFT_align_data_vectors[1][counter];
        HFT_hit_IR_align = HFT_align_data_vectors[2][counter];
        HFT_hit_OR_align = HFT_align_data_vectors[3][counter];

        // apply rotation and translation
        if(flag_left_right_align == 0)
        {
            HFT_hit_IR_align.RotateZ(rot_alpha);
            HFT_hit_IR_align.RotateX(rot_beta);
            HFT_hit_IR_align.RotateZ(rot_gamma);
            HFT_hit_IR_align += shift_vec;

            HFT_hit_OR_align.RotateZ(rot_alpha);
            HFT_hit_OR_align.RotateX(rot_beta);
            HFT_hit_OR_align.RotateZ(rot_gamma);
            HFT_hit_OR_align += shift_vec;

            // set direction vectors for left half
            HFT_dir_L_align = HFT_align_data_vectors[4][counter];

            dca_vec_IR_align = calculateDCA_vec_StraightToPoint(HFT_hit_IL_align,HFT_dir_L_align,HFT_hit_IR_align); // base,dir,point
            dca_vec_OR_align = calculateDCA_vec_StraightToPoint(HFT_hit_IL_align,HFT_dir_L_align,HFT_hit_OR_align); // base,dir,point

            Double_t dca_vec_IR_align_mag = dca_vec_IR_align.Mag();
            Double_t dca_vec_OR_align_mag = dca_vec_OR_align.Mag();

            total_dca += dca_vec_IR_align_mag;
            total_dca += dca_vec_OR_align_mag;
        }
        if(flag_left_right_align == 1)
        {
            HFT_hit_IL_align.RotateZ(rot_alpha);
            HFT_hit_IL_align.RotateX(rot_beta);
            HFT_hit_IL_align.RotateZ(rot_gamma);
            HFT_hit_IL_align += shift_vec;

            HFT_hit_OL_align.RotateZ(rot_alpha);
            HFT_hit_OL_align.RotateX(rot_beta);
            HFT_hit_OL_align.RotateZ(rot_gamma);
            HFT_hit_OL_align += shift_vec;

            // set direction vectors for left half
            HFT_dir_R_align = HFT_align_data_vectors[5][counter];

            dca_vec_IL_align = calculateDCA_vec_StraightToPoint(HFT_hit_IR_align,HFT_dir_R_align,HFT_hit_IL_align); // base,dir,point
            dca_vec_OL_align = calculateDCA_vec_StraightToPoint(HFT_hit_IR_align,HFT_dir_R_align,HFT_hit_OL_align); // base,dir,point

            Double_t dca_vec_IL_align_mag = dca_vec_IL_align.Mag();
            Double_t dca_vec_OL_align_mag = dca_vec_OL_align.Mag();

            total_dca += dca_vec_IL_align_mag;
            total_dca += dca_vec_OL_align_mag;
        }
    }

    cout << "total_dca = " << total_dca << ", delta_x = " << delta_x << ", delta_y = " << delta_y << ", delta_z ="
        <<  delta_z << ", alpha = " << rot_alpha << ", beta = " << rot_beta << ", gamma = " << rot_gamma << endl;

    return total_dca;
}
//------------------------------------------------------------



//------------------------------------------------------------
Double_t Align_HFT_STS(const Double_t *align_params)
{
    // Function for sector-to-sector alignment
    Double_t delta_x   = align_params[0];
    Double_t delta_y   = align_params[1];
    Double_t delta_z   = align_params[2];
    Double_t rot_alpha = align_params[3];
    Double_t rot_beta  = align_params[4];
    Double_t rot_gamma = align_params[5];

    Double_t total_dca = 0.0;

    TVector3 HFT_hit_IR_align, HFT_hit_OR_align, HFT_hit_IA_align, HFT_hit_OA_align, HFT_dir_R_align, HFT_dir_A_align,
        dca_vec_IA_align, dca_vec_OA_align;

    // define and set translation vector
    TVector3 shift_vec;
    shift_vec.SetXYZ(delta_x,delta_y,delta_z);

    for(Long64_t counter = 0; counter < PXL_sec_align_data_vectors[0].size(); counter++)
    {
        // set HFT hits for inner/outer left/right
        HFT_hit_IR_align = PXL_sec_align_data_vectors[0][counter]; // reference inner
        HFT_hit_OR_align = PXL_sec_align_data_vectors[1][counter]; // reference outer
        HFT_hit_IA_align = PXL_sec_align_data_vectors[2][counter]; // align inner
        HFT_hit_OA_align = PXL_sec_align_data_vectors[3][counter]; // align outer

        // apply rotation and translation
        if(flag_left_right_align == 0)
        {
            HFT_hit_IA_align.RotateZ(rot_alpha);
            HFT_hit_IA_align.RotateX(rot_beta);
            HFT_hit_IA_align.RotateZ(rot_gamma);
            HFT_hit_IA_align += shift_vec;

            HFT_hit_OA_align.RotateZ(rot_alpha);
            HFT_hit_OA_align.RotateX(rot_beta);
            HFT_hit_OA_align.RotateZ(rot_gamma);
            HFT_hit_OA_align += shift_vec;

            // set direction vectors for left half
            HFT_dir_R_align = PXL_sec_align_data_vectors[4][counter]; // dir reference

            dca_vec_IA_align = calculateDCA_vec_StraightToPoint(HFT_hit_IR_align,HFT_dir_R_align,HFT_hit_IA_align); // base,dir,point
            dca_vec_OA_align = calculateDCA_vec_StraightToPoint(HFT_hit_IR_align,HFT_dir_R_align,HFT_hit_OA_align); // base,dir,point

            Double_t dca_vec_IA_align_mag = dca_vec_IA_align.Mag();
            Double_t dca_vec_OA_align_mag = dca_vec_OA_align.Mag();

            total_dca += dca_vec_IA_align_mag;
            total_dca += dca_vec_OA_align_mag;
        }
    }

    //cout << "total_dca = " << total_dca << ", delta_x = " << delta_x << ", delta_y = " << delta_y << ", delta_z ="
    //    <<  delta_z << ", alpha = " << rot_alpha << ", beta = " << rot_beta << ", gamma = " << rot_gamma << endl;

    return total_dca;
}
//------------------------------------------------------------



//------------------------------------------------------------
Double_t Align_HFT_IST(const Double_t *align_params)
{
    // Function for sector-to-sector alignment
    Double_t delta_x   = align_params[0];
    Double_t delta_y   = align_params[1];
    Double_t delta_z   = align_params[2];
    Double_t rot_alpha = align_params[3];
    Double_t rot_beta  = align_params[4];
    Double_t rot_gamma = align_params[5];

    Double_t total_dca = 0.0;

    TVector3 HFT_hit_IA_align, HFT_hit_OA_align, HFT_hit_IB_align, HFT_hit_OB_align, HFT_dir_A_align, HFT_dir_B_align,
        dca_vec_A_align, dca_vec_B_align, IST_hit;

    // define and set translation vector
    TVector3 shift_vec_IST;
    shift_vec_IST.SetXYZ(delta_x,delta_y,delta_z);

    for(Long64_t counter = 0; counter < IST_align_data_vectors[0].size(); counter++)
    {
        // set HFT hits for inner/outer left/right
        HFT_hit_IA_align = IST_align_data_vectors[0][counter]; // PXL IA
        HFT_hit_OA_align = IST_align_data_vectors[1][counter]; // PXL OA
        HFT_hit_IB_align = IST_align_data_vectors[2][counter]; // PXL IB
        HFT_hit_OB_align = IST_align_data_vectors[3][counter]; // PXL OB
        IST_hit          = IST_align_data_vectors[6][counter]; // IST hit

        // set direction vectors for PXL
        HFT_dir_A_align = IST_align_data_vectors[4][counter]; // dir A
        HFT_dir_B_align = IST_align_data_vectors[5][counter]; // dir B

        // apply rotation and translation
        IST_hit.RotateZ(rot_alpha);
        IST_hit.RotateX(rot_beta);
        IST_hit.RotateZ(rot_gamma);
        IST_hit += shift_vec_IST;

        dca_vec_A_align = calculateDCA_vec_StraightToPoint(HFT_hit_IA_align,HFT_dir_A_align,IST_hit); // base,dir,point
        dca_vec_B_align = calculateDCA_vec_StraightToPoint(HFT_hit_IB_align,HFT_dir_B_align,IST_hit); // base,dir,point

        // 170 / 1800 // resolution in mum along x and z for IST taken from
        // http://phys.kent.edu/~margetis/STAR/HFT/Work/Talks/conferences/Margetis-final.pdf (slide 6)

        Double_t sigmaX = 0.017;
        Double_t sigmaY = 0.017;
        Double_t sigmaZ = 0.18;

        Double_t dca_vec_A_align_mag = ((dca_vec_A_align.X()*dca_vec_A_align.X())/(sigmaX*sigmaX)) +
            ((dca_vec_A_align.Y()*dca_vec_A_align.Y())/(sigmaY*sigmaY)) +
            ((dca_vec_A_align.Z()*dca_vec_A_align.Z())/(sigmaZ*sigmaZ));

        Double_t dca_vec_B_align_mag = ((dca_vec_B_align.X()*dca_vec_B_align.X())/(sigmaX*sigmaX)) +
            ((dca_vec_B_align.Y()*dca_vec_B_align.Y())/(sigmaY*sigmaY)) +
            ((dca_vec_B_align.Z()*dca_vec_B_align.Z())/(sigmaZ*sigmaZ));

        //Double_t dca_vec_A_align_mag = dca_vec_A_align.Mag();
        //Double_t dca_vec_B_align_mag = dca_vec_B_align.Mag();

        total_dca += dca_vec_A_align_mag;
        total_dca += dca_vec_B_align_mag;
    }

    //cout << "total_dca = " << total_dca << ", delta_x = " << delta_x << ", delta_y = " << delta_y << ", delta_z ="
    //    <<  delta_z << ", alpha = " << rot_alpha << ", beta = " << rot_beta << ", gamma = " << rot_gamma << endl;

    return total_dca;
}
//------------------------------------------------------------



//------------------------------------------------------------
Double_t Align_HFT_to_TPC(const Double_t *align_params)
{
    // Function for sector-to-sector alignment
    Double_t delta_x   = align_params[0];
    Double_t delta_y   = align_params[1];
    Double_t delta_z   = align_params[2];
    Double_t rot_alpha = align_params[3];
    Double_t rot_beta  = align_params[4];
    Double_t rot_gamma = align_params[5];

    Double_t total_dca = 0.0;

    TVector3 HFT_hit_IA_align, HFT_hit_OA_align, HFT_hit_IB_align, HFT_hit_OB_align, HFT_dir_A_align, HFT_dir_B_align,
        dca_vec_A_align, dca_vec_B_align, IST_hit;

    // define and set translation vector
    TVector3 shift_vec_IST;
    shift_vec_IST.SetXYZ(delta_x,delta_y,delta_z);

    // 170 / 1800 // resolution in mum along x and z for IST taken from
    // http://phys.kent.edu/~margetis/STAR/HFT/Work/Talks/conferences/Margetis-final.pdf (slide 6)

    Double_t sigmaX    = 0.017;
    Double_t sigmaY    = 0.017;
    Double_t sigmaZ    = 0.18;
    Double_t sigma_PXL = 0.002;
    Double_t sigma_TPC = 0.3;

    Double_t sigmaX_TPC   = TMath::Sqrt(sigmaX*sigmaX + sigma_TPC*sigma_TPC);
    Double_t sigmaY_TPC   = TMath::Sqrt(sigmaY*sigmaY + sigma_TPC*sigma_TPC);
    Double_t sigmaZ_TPC   = TMath::Sqrt(sigmaZ*sigmaZ + sigma_TPC*sigma_TPC);
    Double_t sigmaIST_TPC = TMath::Sqrt(sigmaX_TPC*sigmaX_TPC + sigmaY_TPC*sigmaY_TPC + sigmaZ_TPC*sigmaZ_TPC);
    Double_t sigmaPXL_TPC = TMath::Sqrt(sigma_PXL*sigma_PXL + sigma_TPC*sigma_TPC);

    for(Long64_t counter = 0; counter < IST_align_data_vectors[0].size(); counter++)
    {
        // set HFT hits for inner/outer left/right
        HFT_hit_IA_align = IST_align_data_vectors[0][counter]; // PXL IA
        HFT_hit_OA_align = IST_align_data_vectors[1][counter]; // PXL OA
        HFT_hit_IB_align = IST_align_data_vectors[2][counter]; // PXL IB
        HFT_hit_OB_align = IST_align_data_vectors[3][counter]; // PXL OB
        IST_hit          = IST_align_data_vectors[6][counter]; // IST hit

        // apply rotation and translation
        HFT_hit_IA_align.RotateZ(rot_alpha);
        HFT_hit_IA_align.RotateX(rot_beta);
        HFT_hit_IA_align.RotateZ(rot_gamma);
        HFT_hit_IA_align += shift_vec_IST;

        HFT_hit_OA_align.RotateZ(rot_alpha);
        HFT_hit_OA_align.RotateX(rot_beta);
        HFT_hit_OA_align.RotateZ(rot_gamma);
        HFT_hit_OA_align += shift_vec_IST;

        HFT_hit_IB_align.RotateZ(rot_alpha);
        HFT_hit_IB_align.RotateX(rot_beta);
        HFT_hit_IB_align.RotateZ(rot_gamma);
        HFT_hit_IB_align += shift_vec_IST;

        HFT_hit_OB_align.RotateZ(rot_alpha);
        HFT_hit_OB_align.RotateX(rot_beta);
        HFT_hit_OB_align.RotateZ(rot_gamma);
        HFT_hit_OB_align += shift_vec_IST;

        IST_hit.RotateZ(rot_alpha);
        IST_hit.RotateX(rot_beta);
        IST_hit.RotateZ(rot_gamma);
        IST_hit += shift_vec_IST;

        StThreeVectorF HFT_hit_IA_align_STV, HFT_hit_OA_align_STV, HFT_hit_IB_align_STV, HFT_hit_OB_align_STV, IST_hit_STV;
        HFT_hit_IA_align_STV.set(HFT_hit_IA_align.x(),HFT_hit_IA_align.y(),HFT_hit_IA_align.z());
        HFT_hit_OA_align_STV.set(HFT_hit_OA_align.x(),HFT_hit_OA_align.y(),HFT_hit_OA_align.z());
        HFT_hit_IB_align_STV.set(HFT_hit_IB_align.x(),HFT_hit_IB_align.y(),HFT_hit_IB_align.z());
        HFT_hit_OB_align_STV.set(HFT_hit_OB_align.x(),HFT_hit_OB_align.y(),HFT_hit_OB_align.z());
        IST_hit_STV.set(IST_hit.x(),IST_hit.y(),IST_hit.z());

        for(Int_t i_TPC_track = 0; i_TPC_track < 2; i_TPC_track++)
        {
            Float_t pathIA,dcaIA;
            fHelixAtoPointdca(HFT_hit_IA_align_STV,TPC_align_data_vectors[i_TPC_track][counter],pathIA,dcaIA);

            Float_t pathOA,dcaOA;
            fHelixAtoPointdca(HFT_hit_OA_align_STV,TPC_align_data_vectors[i_TPC_track][counter],pathOA,dcaOA);

            Float_t pathIB,dcaIB;
            fHelixAtoPointdca(HFT_hit_IB_align_STV,TPC_align_data_vectors[i_TPC_track][counter],pathIB,dcaIB);

            Float_t pathOB,dcaOB;
            fHelixAtoPointdca(HFT_hit_OB_align_STV,TPC_align_data_vectors[i_TPC_track][counter],pathOB,dcaOB);

            Float_t pathIST,dcaIST;
            fHelixAtoPointdca(IST_hit_STV,TPC_align_data_vectors[i_TPC_track][counter],pathIST,dcaIST);

            Double_t dca_track =
                (dcaIA*dcaIA)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaIB*dcaIB)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaOA*dcaOA)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaOB*dcaOB)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaIST*dcaIST)/(sigmaIST_TPC*sigmaIST_TPC);

            total_dca += dca_track;

            //cout << "counter = " << counter << ", i_TPC_track = " << i_TPC_track << ", dca_track = " << dca_track << ", total_dca = " << total_dca
            //    << ", dcaIA = " << dcaIA << ", dcaIB = " << dcaIB << ", dcaOA = " << dcaOA << ", dcaOB = " << dcaOB << ", dcaIST = " << dcaIST << endl;
        }
    }

    cout << "total_dca = " << total_dca << ", delta_x = " << delta_x << ", delta_y = " << delta_y << ", delta_z ="
        <<  delta_z << ", alpha = " << rot_alpha << ", beta = " << rot_beta << ", gamma = " << rot_gamma << endl;

    return total_dca;
}
//------------------------------------------------------------



//------------------------------------------------------------
Double_t Align_TPC_Sector(const Double_t *align_params)
{
    // Function for sector-to-sector alignment
    Double_t delta_x   = align_params[0];
    Double_t delta_y   = align_params[1];
    Double_t delta_z   = align_params[2];
    Double_t rot_alpha = align_params[3];
    Double_t rot_beta  = align_params[4];
    Double_t rot_gamma = align_params[5];

    Double_t total_dca = 0.0;

    TVector3 HFT_hit_IA_align, HFT_hit_OA_align, HFT_hit_IB_align, HFT_hit_OB_align, HFT_dir_A_align, HFT_dir_B_align,
        dca_vec_A_align, dca_vec_B_align, IST_hit;

    // define and set translation vector
    TVector3 shift_vec_IST;
    shift_vec_IST.SetXYZ(delta_x,delta_y,delta_z);
    StThreeVectorF shift_vec_TPC_track;
    shift_vec_TPC_track.set(delta_x,delta_y,delta_z);

    // 170 / 1800 // resolution in mum along x and z for IST taken from
    // http://phys.kent.edu/~margetis/STAR/HFT/Work/Talks/conferences/Margetis-final.pdf (slide 6)

    Double_t sigmaX    = 0.017;
    Double_t sigmaY    = 0.017;
    Double_t sigmaZ    = 0.18;
    Double_t sigma_PXL = 0.002;
    Double_t sigma_TPC = 0.3;

    Double_t sigmaX_TPC   = TMath::Sqrt(sigmaX*sigmaX + sigma_TPC*sigma_TPC);
    Double_t sigmaY_TPC   = TMath::Sqrt(sigmaY*sigmaY + sigma_TPC*sigma_TPC);
    Double_t sigmaZ_TPC   = TMath::Sqrt(sigmaZ*sigmaZ + sigma_TPC*sigma_TPC);
    Double_t sigmaIST_TPC = TMath::Sqrt(sigmaX_TPC*sigmaX_TPC + sigmaY_TPC*sigmaY_TPC + sigmaZ_TPC*sigmaZ_TPC);
    Double_t sigmaPXL_TPC = TMath::Sqrt(sigma_PXL*sigma_PXL + sigma_TPC*sigma_TPC);

#if 0
    cout << "TPC sector aling, tracks: " << IST_align_data_vectors[0].size() << endl;
#endif

    for(Long64_t counter = 0; counter < IST_align_data_vectors[0].size(); counter++)
    {
        // set HFT hits for inner/outer left/right
        HFT_hit_IA_align = IST_align_data_vectors[0][counter]; // PXL IA
        HFT_hit_OA_align = IST_align_data_vectors[1][counter]; // PXL OA
        HFT_hit_IB_align = IST_align_data_vectors[2][counter]; // PXL IB
        HFT_hit_OB_align = IST_align_data_vectors[3][counter]; // PXL OB
        IST_hit          = IST_align_data_vectors[6][counter]; // IST hit

        /*
        // apply rotation and translation
        HFT_hit_IA_align.RotateZ(rot_alpha);
        HFT_hit_IA_align.RotateX(rot_beta);
        HFT_hit_IA_align.RotateZ(rot_gamma);
        HFT_hit_IA_align += shift_vec_IST;

        HFT_hit_OA_align.RotateZ(rot_alpha);
        HFT_hit_OA_align.RotateX(rot_beta);
        HFT_hit_OA_align.RotateZ(rot_gamma);
        HFT_hit_OA_align += shift_vec_IST;

        HFT_hit_IB_align.RotateZ(rot_alpha);
        HFT_hit_IB_align.RotateX(rot_beta);
        HFT_hit_IB_align.RotateZ(rot_gamma);
        HFT_hit_IB_align += shift_vec_IST;

        HFT_hit_OB_align.RotateZ(rot_alpha);
        HFT_hit_OB_align.RotateX(rot_beta);
        HFT_hit_OB_align.RotateZ(rot_gamma);
        HFT_hit_OB_align += shift_vec_IST;

        IST_hit.RotateZ(rot_alpha);
        IST_hit.RotateX(rot_beta);
        IST_hit.RotateZ(rot_gamma);
        IST_hit += shift_vec_IST;
        */

        StThreeVectorF HFT_hit_IA_align_STV, HFT_hit_OA_align_STV, HFT_hit_IB_align_STV, HFT_hit_OB_align_STV, IST_hit_STV;
        HFT_hit_IA_align_STV.set(HFT_hit_IA_align.x(),HFT_hit_IA_align.y(),HFT_hit_IA_align.z());
        HFT_hit_OA_align_STV.set(HFT_hit_OA_align.x(),HFT_hit_OA_align.y(),HFT_hit_OA_align.z());
        HFT_hit_IB_align_STV.set(HFT_hit_IB_align.x(),HFT_hit_IB_align.y(),HFT_hit_IB_align.z());
        HFT_hit_OB_align_STV.set(HFT_hit_OB_align.x(),HFT_hit_OB_align.y(),HFT_hit_OB_align.z());
        IST_hit_STV.set(IST_hit.x(),IST_hit.y(),IST_hit.z());

        for(Int_t i_TPC_track = 0; i_TPC_track < 1; i_TPC_track++)
        {
            StThreeVectorF vec_TPC_Origin = TPC_align_data_vectors[i_TPC_track][counter].origin();
            StThreeVectorF vec_TPC_Mom    = TPC_align_data_vectors[i_TPC_track][counter].momentum(vec_magfac[i_TPC_track][counter]);
            Int_t          TPC_charge     = TPC_align_data_vectors[i_TPC_track][counter].charge(vec_charge[i_TPC_track][counter]);

            vec_TPC_Mom.rotateX(rot_alpha);
            vec_TPC_Mom.rotateY(rot_beta);
            vec_TPC_Mom.rotateZ(rot_gamma);

            vec_TPC_Origin += shift_vec_TPC_track;

            StPhysicalHelixD helix_align = StPhysicalHelixD(vec_TPC_Mom,vec_TPC_Origin,vec_magfac[i_TPC_track][counter],TPC_charge);

            Float_t pathIA,dcaIA;
            if(HFT_hit_IA_align_STV.x() > -998.0)
            {
                fHelixAtoPointdca(HFT_hit_IA_align_STV,helix_align,pathIA,dcaIA);
            }
            else
            {
                dcaIA = 0.0;
            }

            Float_t pathOA,dcaOA;
            if(HFT_hit_OA_align_STV.x() > -998.0)
            {
                fHelixAtoPointdca(HFT_hit_OA_align_STV,helix_align,pathOA,dcaOA);
            }
            else
            {
                dcaOA = 0.0;
            }

            Float_t pathIB,dcaIB;
            if(HFT_hit_IB_align_STV.x() > -998.0)
            {
                fHelixAtoPointdca(HFT_hit_IB_align_STV,helix_align,pathIB,dcaIB);
            }
            else
            {
                dcaIB = 0.0;
            }

            Float_t pathOB,dcaOB;
            if(HFT_hit_OB_align_STV.x() > -998.0)
            {
                fHelixAtoPointdca(HFT_hit_OB_align_STV,helix_align,pathOB,dcaOB);
            }
            else
            {
                dcaOB = 0.0;
            }

            Float_t pathIST,dcaIST;
            fHelixAtoPointdca(IST_hit_STV,helix_align,pathIST,dcaIST);
            dcaIST = 0.0;

            Double_t dca_track =
                (dcaIA*dcaIA)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaIB*dcaIB)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaOA*dcaOA)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaOB*dcaOB)/(sigmaPXL_TPC*sigmaPXL_TPC) +
                (dcaIST*dcaIST)/(sigmaIST_TPC*sigmaIST_TPC);

            total_dca += dca_track;

            //cout << "counter = " << counter << ", i_TPC_track = " << i_TPC_track << ", dca_track = " << dca_track << ", total_dca = " << total_dca
            //    << ", dcaIA = " << dcaIA << ", dcaIB = " << dcaIB << ", dcaOA = " << dcaOA << ", dcaOB = " << dcaOB << ", dcaIST = " << dcaIST << endl;
        }
    }

    cout << "TPC sector: " << TPC_sector_align_active << ", N tracks: " << IST_align_data_vectors[0].size()<< ", total_dca = " << total_dca << ", delta_x = " << delta_x << ", delta_y = " << delta_y << ", delta_z ="
        <<  delta_z << ", alpha = " << rot_alpha << ", beta = " << rot_beta << ", gamma = " << rot_gamma << endl;

    return total_dca;
}
//------------------------------------------------------------



void HFT_Alignment(Long64_t start_event_use = 0, Long64_t stop_event_use = 100, Int_t flag_min_option = 0,
                   Int_t pxl_sector_L = -1, Int_t pxl_sector_R = 5)
{

    // .x HFT_Alignment.cc++(0,1000000,26,-1,-1)
    // .x HFT_Alignment.cc++(0,1000000,30,5,-1)

    // flag_min_option

    // half to half PXL alignment
    // 0 == do minimization
    // 1 == use some fixed parameters
    // 2 == do minimization but start with new parameters
    // 3 == do full sector alignment with start parameters

    // PXL sector-to-sector
    // 10 == plot sector-to-sector PXL residuals
    // 11 == do sector-to-sector PXL alignment

    // IST
    // 20 -> no alignment
    // 21 -> global alignment (all ladders together)
    // 22 -> ladder-by-ladder alignment

    // PXL+IST to TPC
    // 23 -> use global alignment but don't do minimization
    // 25 -> plot TPC residuals
    // 26 -> do HFT to TPC alignment
    //

    // TPC alignment
    // 30 -> do TPC sector alignment using HFT PXL information
    // 27 -> do first HFT to TPC alignment (sector-by-sector)
    // 31 -> do TPC sector alignment with pre HFT to TPC alignment applied from file
    // 28 -> do TPC alignment with initial pre alignment for TPC sector-to-sector alignment
    // 32 -> do TPC sector alignment with pre HFT to TPC alignment and pre TPC sector alignment applied from file

    const Double_t min_TPC_track_mom = 3.0;
    const Int_t N_sectors = 10; // number of PXL sectors

    Int_t pxl_sector_L_orig = pxl_sector_L;
    Int_t pxl_sector_R_orig = pxl_sector_R;

    //------------------------------------------------------------------------------
    cout << "HFT Alignment started" << endl;
    SetRootGraphicStyle();
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    Double_t mean, sigma, height,p0,p1,p2,p3,p4,p5;
    Double_t mean_err, sigma_err, height_err;
    Double_t mean1, sigma1, height1;
    char NoP[50];
    TString NoP2, HistName;
    char NoPb1[50];
    char NoP1[50];
    TString NoP21, HistName1;
    Double_t pt_center;
    Double_t start_fit, stop_fit;

    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(2000);
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    cout << "Open input file" << endl;
    NT_HFT_left_right  = new TChain( "NT_HFT_zero_field_cosmic", "NT_HFT_zero_field_cosmic");
    NT_PXL_IST         = new TChain( "NT_PXL_IST_zero_field_cosmic", "NT_PXL_IST_zero_field_cosmic");
    //TString addfile = "./Data/HFT_local_hit_histograms_loop2_nouse.root";
    //TString addfile = "./Data/HFT_local_hit_histograms_loop_Survey.root";
    //TString addfile = "./Data/HFT_local_hit_histograms_loop_Survey_PXL_IST_V3.root";
    //TString addfile = "./Data/HFT_local_hit_histograms_loop_Survey_PXL_IST_V5b.root"; // with Hao's survery data
    //TString addfile = "./Data/HFT_local_hit_histograms_loop_Survey_PXL_IST_V5c.root"; // with PXL+IST alignment
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/HFT_TPC_Survey_aligned_magnet_on.root"; // with PXL+IST alignment magent on
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/HFT_TPC_Survey_aligned_magnet_on_V2.root"; // with PXL+IST alignment magent on, after TPC t0 calibration FF
    //TString addfile = "./Data/HFT_local_hit_histograms_loop_Survey_PXL_IST_V5c_sim_magnet_off.root"; // with PXL+IST alignment magent on
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/HFT_TPC_Survey_aligned_magnet_on_RFF_V4.root"; // with PXL+IST alignment magent on, after TPC t0 calibration RFF
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/HFT_TPC_Survey_aligned_magnet_on_RFF_corr_V4.root"; // with PXL+IST alignment magent on, after TPC t0 calibration RFF

    // new tree structure with TPC sector information
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/Merge_RFF_PXL_align_TPC_out.root"; //
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/Merge_RFF2_PXL_align_TPC_out.root"; //
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/Merge_RFF3_PXL_align_TPC_out_full.root"; //
    //TString addfile = "/project/projectdirs/star/aschmah/HFT/Alignment/Data/Merge_RFF4_PXL_align_TPC_out_full.root";

    TString addfile      = "/project/projectdirs/star/pwg/starbulc/aschmah/HFT/Alignment/Data/Merge_RFF4_PXL_align_TPC_out_full.root";
    TString dca_matrices = "/project/projectdirs/star/pwg/starbulc/aschmah/HFT/Alignment/Data/Align_Matrix_Xin_dca_optimization.root";

    cout << "Open file: " << addfile.Data() << endl;
    NT_HFT_left_right ->AddFile(addfile.Data(),-1, "NT_HFT_zero_field_cosmic" );
    NT_PXL_IST        ->AddFile(addfile.Data(),-1, "NT_PXL_IST_zero_field_cosmic" );
    Long64_t file_entries_HFT = NT_HFT_left_right->GetEntries();
    cout << "Entries in file HFT: " << file_entries_HFT << endl;

    Long64_t file_entries_PXL_IST = NT_PXL_IST->GetEntries();
    cout << "Entries in file PXL IST: " << file_entries_PXL_IST << endl;

    const TString fig_all_output_dir = "./Pictures/"; // for "flag_save_all_figures"
    const TString out_all_format     = ".gif";        // for "flag_save_all_figures"

    cout << "Open alignment from dca optimization (Xin, Mustafa): " << dca_matrices.Data() << endl;
    TFile* file_dca_matrices = TFile::Open(dca_matrices.Data());
    TGeoMatrix* M_dca_matrices[N_sectors];
    for(Int_t i_sec = 0; i_sec < N_sectors; i_sec++)
    {
        HistName = "matrix_";
        HistName += i_sec+1;
        M_dca_matrices[i_sec] = (TGeoMatrix*)file_dca_matrices->Get(HistName.Data());
    }
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    TF1 *GaussFit              = new TF1("GaussFit",GaussFitFunc,0.0,10,3);
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    cout << "Define NTuple HFT (only left right half combinations, no IST)" << endl;
    NT_HFT_left_right->SetBranchAddress("x1",&x1_HFT);
    NT_HFT_left_right->SetBranchAddress("y1",&y1_HFT);
    NT_HFT_left_right->SetBranchAddress("z1",&z1_HFT);
    NT_HFT_left_right->SetBranchAddress("x2",&x2_HFT);
    NT_HFT_left_right->SetBranchAddress("y2",&y2_HFT);
    NT_HFT_left_right->SetBranchAddress("z2",&z2_HFT);
    NT_HFT_left_right->SetBranchAddress("x3",&x3_HFT);
    NT_HFT_left_right->SetBranchAddress("y3",&y3_HFT);
    NT_HFT_left_right->SetBranchAddress("z3",&z3_HFT);
    NT_HFT_left_right->SetBranchAddress("x4",&x4_HFT);
    NT_HFT_left_right->SetBranchAddress("y4",&y4_HFT);
    NT_HFT_left_right->SetBranchAddress("z4",&z4_HFT);
    NT_HFT_left_right->SetBranchAddress("id_IL",&id_IL_HFT);
    NT_HFT_left_right->SetBranchAddress("id_OL",&id_OL_HFT);
    NT_HFT_left_right->SetBranchAddress("id_IR",&id_IR_HFT);
    NT_HFT_left_right->SetBranchAddress("id_OR",&id_OR_HFT);

    cout << "Define NTuple PXL IST (all sector combinations)" << endl;
    NT_PXL_IST->SetBranchAddress("PXL_IA_x",&PXL_IA_x);
    NT_PXL_IST->SetBranchAddress("PXL_IA_y",&PXL_IA_y);
    NT_PXL_IST->SetBranchAddress("PXL_IA_z",&PXL_IA_z);
    NT_PXL_IST->SetBranchAddress("PXL_IB_x",&PXL_IB_x);
    NT_PXL_IST->SetBranchAddress("PXL_IB_y",&PXL_IB_y);
    NT_PXL_IST->SetBranchAddress("PXL_IB_z",&PXL_IB_z);
    NT_PXL_IST->SetBranchAddress("PXL_OA_x",&PXL_OA_x);
    NT_PXL_IST->SetBranchAddress("PXL_OA_y",&PXL_OA_y);
    NT_PXL_IST->SetBranchAddress("PXL_OA_z",&PXL_OA_z);
    NT_PXL_IST->SetBranchAddress("PXL_OB_x",&PXL_OB_x);
    NT_PXL_IST->SetBranchAddress("PXL_OB_y",&PXL_OB_y);
    NT_PXL_IST->SetBranchAddress("PXL_OB_z",&PXL_OB_z);
    NT_PXL_IST->SetBranchAddress("PXL_id_IA",&PXL_id_IA);
    NT_PXL_IST->SetBranchAddress("PXL_id_IB",&PXL_id_IB);
    NT_PXL_IST->SetBranchAddress("PXL_id_OA",&PXL_id_OA);
    NT_PXL_IST->SetBranchAddress("PXL_id_OB",&PXL_id_OB);
    NT_PXL_IST->SetBranchAddress("IST_x",&IST_x);
    NT_PXL_IST->SetBranchAddress("IST_y",&IST_y);
    NT_PXL_IST->SetBranchAddress("IST_z",&IST_z);
    NT_PXL_IST->SetBranchAddress("IST_id",&IST_id);
    NT_PXL_IST->SetBranchAddress("IST_counter",&IST_counter);
    NT_PXL_IST->SetBranchAddress("IST_hit",&IST_hit);

    if(flag_min_option >= 25 && flag_min_option < 34)
    {
        NT_PXL_IST->SetBranchAddress("magfac",&magfac);
        NT_PXL_IST->SetBranchAddress("qA",&track_qA);
        NT_PXL_IST->SetBranchAddress("pxA",&track_pxA);
        NT_PXL_IST->SetBranchAddress("pyA",&track_pyA);
        NT_PXL_IST->SetBranchAddress("pzA",&track_pzA);
        NT_PXL_IST->SetBranchAddress("oxA",&track_oxA);
        NT_PXL_IST->SetBranchAddress("oyA",&track_oyA);
        NT_PXL_IST->SetBranchAddress("ozA",&track_ozA);
        NT_PXL_IST->SetBranchAddress("qB",&track_qB);
        NT_PXL_IST->SetBranchAddress("pxB",&track_pxB);
        NT_PXL_IST->SetBranchAddress("pyB",&track_pyB);
        NT_PXL_IST->SetBranchAddress("pzB",&track_pzB);
        NT_PXL_IST->SetBranchAddress("oxB",&track_oxB);
        NT_PXL_IST->SetBranchAddress("oyB",&track_oyB);
        NT_PXL_IST->SetBranchAddress("ozB",&track_ozB);

        NT_PXL_IST->SetBranchAddress("secA",&secA);
        NT_PXL_IST->SetBranchAddress("secB",&secB);
        NT_PXL_IST->SetBranchAddress("NsecA",&NsecA);
        NT_PXL_IST->SetBranchAddress("NsecB",&NsecB);
        NT_PXL_IST->SetBranchAddress("NnohitsecA",&NnohitsecA);
        NT_PXL_IST->SetBranchAddress("NnohitsecB",&NnohitsecB);

        NT_PXL_IST->SetBranchAddress("idx_OA",&idx_OA);
        NT_PXL_IST->SetBranchAddress("idx_OB",&idx_OB);
        NT_PXL_IST->SetBranchAddress("idx_IA",&idx_IA);
        NT_PXL_IST->SetBranchAddress("idx_IB",&idx_IB);
        NT_PXL_IST->SetBranchAddress("idx_IST",&idx_IST);
    }
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    const Int_t N_pixel_dca_cuts = 10;
    TString xyz_label[3] = {"#Delta x_{global} (cm)","#Delta y_{global} (cm)","#Delta z_{global} (cm)"};

    TString label_pxl_sector_L = "sec_{L} = ";
    label_pxl_sector_L += pxl_sector_L+1;
    if(pxl_sector_L == -1) label_pxl_sector_L = "sec_{L} = all";

    TString label_pxl_sector_R = "sec_{R} = ";
    label_pxl_sector_R += pxl_sector_R+1;
    if(pxl_sector_R == -1) label_pxl_sector_R = "sec_{R} = all";

    TVector3 HFT_hit_IL,HFT_hit_OL,HFT_hit_IR,HFT_hit_OR,HFT_dir_L, HFT_dir_R, dca_vec_IR, dca_vec_OR;
    TH1F* h_pixel_line_dca[N_pixel_dca_cuts][3];
    TH2F* h2D_pixel_line_dca_vs_z[N_pixel_dca_cuts][3];

    for(Int_t i_cut = 0; i_cut < N_pixel_dca_cuts; i_cut++)
    {
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            HistName = "h_pixel_line_dca_";
            HistName += i_cut;
            HistName += "_xyz_";
            HistName += i_xyz;
            h_pixel_line_dca[i_cut][i_xyz] = new TH1F(HistName.Data(),HistName.Data(),10000,-10,10);

            HistName = "h2D_pixel_line_dca_vs_z";
            HistName += i_cut;
            HistName += "_xyz_";
            HistName += i_xyz;
            h2D_pixel_line_dca_vs_z[i_cut][i_xyz] = new TH2F(HistName.Data(),HistName.Data(),100,-25,25,800,-0.4,0.4);
        }
    }


    const Int_t N_HFT_TPC_match = 5;
    TString label_IST[N_HFT_TPC_match] = {"PXL inner A","PXL outer A","PXL inner B","PXL outer B","IST"};
    TH1F* h_pixel_line_dca_TPC[2][N_HFT_TPC_match][3][2];
    for(Int_t i_align = 0; i_align < 2; i_align++)
    {
        for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                for(Int_t i_track = 0; i_track < 2; i_track++)
                {
                    HistName = "h_pixel_line_dca_TPC";
                    HistName += "_align_";
                    HistName += i_align;
                    HistName += "_hist_";
                    HistName += i_histo;
                    HistName += "_xyz_";
                    HistName += i_xyz;
                    HistName += "_track_";
                    HistName += i_track;
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track] = new TH1F(HistName.Data(),HistName.Data(),3000,-10,10);
                }
            }
        }
    }


    TCanvas* c_pixel_line_dca_TPC[2]; // [inner,outer][x,y,z][trackA, trackB]
    for(Int_t i_track = 0; i_track < 2; i_track++)
    {
        HistName = "c_pixel_line_dca_TPC_";
        HistName += i_track;
        c_pixel_line_dca_TPC[i_track] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
        c_pixel_line_dca_TPC[i_track]->SetFillColor(10);
        c_pixel_line_dca_TPC[i_track]->SetTopMargin(0.05);
        c_pixel_line_dca_TPC[i_track]->SetBottomMargin(0.15);
        c_pixel_line_dca_TPC[i_track]->SetRightMargin(0.05);
        c_pixel_line_dca_TPC[i_track]->SetLeftMargin(0.15);
        c_pixel_line_dca_TPC[i_track]->Divide(3,N_HFT_TPC_match);


        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
            {
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetTopMargin(0.05);
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetBottomMargin(0.22);
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetRightMargin(0.05);
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetLeftMargin(0.2);
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetTicks(1,1);
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetGrid(0,0);
                c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo)->SetFillColor(10);

                for(Int_t i_align = 0; i_align < 2; i_align++)
                {
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetStats(0);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetTitle("");
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->SetTitleOffset(1.1);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->SetTitleOffset(1.2);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->SetLabelOffset(0.01);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->SetLabelSize(0.065);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->SetLabelSize(0.065);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->SetTitleSize(0.065);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->SetTitleSize(0.065);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->SetNdivisions(505,'N');
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->SetNdivisions(505,'N');
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->CenterTitle();
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->CenterTitle();
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->SetTitle(xyz_label[i_xyz]);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetYaxis()->SetTitle("counts");
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetXaxis()->SetRangeUser(-1.9,1.9);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetFillColor(kGray+1);
                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetFillStyle(3001);
                }
            }
        }
    }
    //------------------------------------------------------------------------------



    //------------------------------------------------------------------------------
    cout << "Loop data" << endl;
    HFT_align_data_vectors.resize(6); // hit_IL, hit_OL, hit_IR, hit_OR, dir_L, dir_R  (hit_OL == PXL hit Outer Left half)

    Double_t delta_x;
    Double_t delta_y;
    Double_t delta_z;
    Double_t rot_alpha;
    Double_t rot_beta;
    Double_t rot_gamma;

    Double_t delta_x_err;
    Double_t delta_y_err;
    Double_t delta_z_err;
    Double_t rot_alpha_err;
    Double_t rot_beta_err;
    Double_t rot_gamma_err;

    Double_t delta_x_init;
    Double_t delta_y_init;
    Double_t delta_z_init;
    Double_t rot_alpha_init;
    Double_t rot_beta_init;
    Double_t rot_gamma_init;

    delta_x_init   = -0.203405;
    delta_y_init   = -0.0572426;
    delta_z_init   = 0.0129719;
    rot_alpha_init = 0.680852;
    rot_beta_init  = 0.00506325;
    rot_gamma_init = -0.679878;


    TVector3 shift_vec, shift_vec_init, shift_vec_init1, shift_vec_init2, shift_vec_total;
    shift_vec_init.SetXYZ(delta_x_init,delta_y_init,delta_z_init);

    TRotation rot1, rot2, rot3, rot4;
    rot1.RotateZ(rot_alpha_init);
    rot1.RotateX(rot_beta_init);
    rot1.RotateZ(rot_gamma_init);



    //-------------------------------------
    shift_vec_init2.SetXYZ(0.00346756,-0.0015464,-0.0183748);

    rot2.RotateZ(-0.522878);
    rot2.RotateX(-0.00328327);
    rot2.RotateZ(0.522384);

    shift_vec_init.Transform(rot2);
    shift_vec_init += shift_vec_init2;

    rot3 = rot2 * rot1;
    rot1 = rot3;

    rot1.SetToIdentity();
    shift_vec_init.SetXYZ(0.0,0.0,0.0);

    rot_HFT_to_TPC.SetToIdentity();
    shift_vec_HFT_to_TPC.SetXYZ(0.0,0.0,0.0);

    rot_HFT_to_TPC_load.SetToIdentity();
    shift_vec_HFT_to_TPC_load.SetXYZ(0.0,0.0,0.0);

    rot_HFT_to_TPC_loadB.SetToIdentity();
    shift_vec_HFT_to_TPC_loadB.SetXYZ(0.0,0.0,0.0);

    for(Int_t i_sec = 0; i_sec < N_TPC_sectors; i_sec++)
    {
        rot_TPC_sector[i_sec].SetToIdentity();
        shift_vec_TPC_sector[i_sec].SetXYZ(0.0,0.0,0.0);
    }

    // rotation matrix scheme
    //  | x' |   | xx xy xz | | x |
    //  | y' | = | yx yy yz | | y |
    //  | z' |   | zx zy zz | | z |

    cout << "" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "rotation matrix scheme" << endl;
    cout << "| x' |   | xx xy xz | | x |" << endl;
    cout << "| y' | = | yx yy yz | | y |" << endl;
    cout << "| z' |   | zx zy zz | | z |" << endl;
    cout << "|" << rot1.XX() << " " << rot1.XY() << " " << rot1.XZ() << "|" << endl;
    cout << "|" << rot1.YX() << " " << rot1.YY() << " " << rot1.YZ() << "|" << endl;
    cout << "|" << rot1.ZX() << " " << rot1.ZY() << " " << rot1.ZZ() << "|" << endl;
    cout << "shift vector = {" << shift_vec_init.X() << ", "  << shift_vec_init.Y() << ", " << shift_vec_init.Z() << "}" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "" << endl;


    //shift_vec_init.SetXYZ(delta_x_init,delta_y_init,delta_z_init);

    //delta_x_init   = shift_vec_init.X();
    //delta_y_init   = shift_vec_init.Y();
    //delta_z_init   = shift_vec_init.Z();

    //rot_alpha_init = rot3.PhiZ();
    //rot_beta_init  = rot3.PhiX();
    //rot_gamma_init = rot3.PhiZ();
    //-------------------------------------

    Int_t IST_ladder_number;

    std::vector<TRotation> vector_rotation_sectors;
    std::vector<TVector3>  vector_shift_sectors;
    vector_rotation_sectors.resize(N_sectors);
    vector_shift_sectors.resize(N_sectors);

    for(Int_t i_full_sec_align = 0; i_full_sec_align < N_sectors; i_full_sec_align++) // loop for full sector alignment, only activated if flag_min_option == 3
    {
        vector_rotation_sectors[i_full_sec_align].SetToIdentity();
        vector_shift_sectors[i_full_sec_align].SetXYZ(0.0,0.0,0.0);
    }

    Int_t sec_comb_full_sec_align[2][9] = // [reference,align][sector combination]
    {
        {0,0,6,6,2,8,3,3,9},
        {5,6,1,2,8,3,9,7,4}
    };
    cout << "All containers defined" << endl;

    TString label_IO_pixel[2] = {"inner pxl.","outer pxl."};


    //---------------------------------------------------------------------------------------------------------------
    if(flag_min_option < 10) // use NT_HFT_left_right Ntuple --> left right half alignment
    {
        Int_t N_full_sec_align = 1;
        if(flag_min_option == 3) N_full_sec_align = 9;

        cout << "Doing left right half PXL alignment" << endl;
        for(Int_t i_full_sec_align = 0; i_full_sec_align < N_full_sec_align; i_full_sec_align++) // loop for full sector alignment, only activated if flag_min_option == 3
        {
            if(flag_min_option == 3) // for full sector alignment
            {
                // clear data vectors
                for(Int_t i_data = 0; i_data < 6; i_data++)
                {
                    HFT_align_data_vectors[i_data].clear();
                }

                if(sec_comb_full_sec_align[0][i_full_sec_align] <= 4)
                {
                    flag_left_right_align = 0;
                }
                else
                {
                    flag_left_right_align = 1;
                }

                // define sector combination
                if(sec_comb_full_sec_align[0][i_full_sec_align] <= 4)
                {
                    pxl_sector_L = sec_comb_full_sec_align[0][i_full_sec_align];
                    pxl_sector_R = sec_comb_full_sec_align[1][i_full_sec_align];
                }
                else
                {
                    pxl_sector_L = sec_comb_full_sec_align[1][i_full_sec_align];
                    pxl_sector_R = sec_comb_full_sec_align[0][i_full_sec_align];
                }

                cout << "pxl_sector_L = " << pxl_sector_L << ", pxl_sector_R = " << pxl_sector_R << endl;
            }

            cout << "Fill data vectors" << endl;
            if(stop_event_use > file_entries_HFT) stop_event_use = file_entries_HFT;
            for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
            {
                if (counter != 0  &&  counter % 1000 == 0)
                    cout << "." << flush;
                if (counter != 0  &&  counter % 10000 == 0)
                {
                    if((stop_event_use-start_event_use) > 0)
                    {
                        Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                        cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (HFT Alignment) " << flush;
                    }
                }

                if (!NT_HFT_left_right->GetEntry( counter )) // take the event -> information is stored in event
                    break;  // end of data chunk

                // Detector index
                // 0 = inner pixel left
                // 1 = outer pixel left
                // 2 = IST left
                // 3 = inner pixel right
                // 4 = outer pixel right
                // 5 = IST right
                Int_t index_pxl_IR = get_HFT_det_index(id_IR_HFT);
                Int_t index_pxl_OR = get_HFT_det_index(id_OR_HFT);

                Int_t sector_pxl_IL = get_HFT_pixel_sector(id_IL_HFT,IST_ladder_number);
                Int_t sector_pxl_OL = get_HFT_pixel_sector(id_OL_HFT,IST_ladder_number);

                Int_t sector_pxl_IR = get_HFT_pixel_sector(id_IR_HFT,IST_ladder_number);
                Int_t sector_pxl_OR = get_HFT_pixel_sector(id_OR_HFT,IST_ladder_number);

                // set HFT hits for inner/outer left/right
                HFT_hit_IL.SetXYZ(x1_HFT,y1_HFT,z1_HFT);
                HFT_hit_OL.SetXYZ(x2_HFT,y2_HFT,z2_HFT);
                HFT_hit_IR.SetXYZ(x3_HFT,y3_HFT,z3_HFT);
                HFT_hit_OR.SetXYZ(x4_HFT,y4_HFT,z4_HFT);

                if(flag_min_option == 2 || flag_min_option == 3)
                {
                    // apply rotation and translation
                    HFT_hit_IR.Transform(rot1);
                    HFT_hit_IR += shift_vec_init;

                    HFT_hit_OR.Transform(rot1);
                    HFT_hit_OR += shift_vec_init;

                    if(flag_min_option == 3)
                    {
                        HFT_hit_IL.Transform(vector_rotation_sectors[sector_pxl_IL]);
                        HFT_hit_IL += vector_shift_sectors[sector_pxl_IL];

                        HFT_hit_OL.Transform(vector_rotation_sectors[sector_pxl_OL]);
                        HFT_hit_OL += vector_shift_sectors[sector_pxl_OL];

                        HFT_hit_IR.Transform(vector_rotation_sectors[sector_pxl_IR]);
                        HFT_hit_IR += vector_shift_sectors[sector_pxl_IR];

                        HFT_hit_OR.Transform(vector_rotation_sectors[sector_pxl_OR]);
                        HFT_hit_OR += vector_shift_sectors[sector_pxl_OR];
                    }
                }

                // calculate direction vectors for left and right halfs
                HFT_dir_L = HFT_hit_IL - HFT_hit_OL;
                HFT_dir_R = HFT_hit_IR - HFT_hit_OR;

                // calculate dca values of left tracklets to right HFT hits for inner/outer
                dca_vec_IR = calculateDCA_vec_StraightToPoint(HFT_hit_IL,HFT_dir_L,HFT_hit_IR); // base,dir,point
                dca_vec_OR = calculateDCA_vec_StraightToPoint(HFT_hit_IL,HFT_dir_L,HFT_hit_OR); // base,dir,point

                //if(dca_vec_IR.Mag() < 0.08 && dca_vec_OR.Mag() < 0.08)
                if(dca_vec_IR.Mag() < 0.5 && dca_vec_OR.Mag() < 0.5)
                {
                    if(
                       (pxl_sector_R >= 0 && pxl_sector_L == -1 && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R) // take all sectors from left and one specific from right
                       || (pxl_sector_L >= 0 && pxl_sector_R == -1 && sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == pxl_sector_L) // take all sectors from right and one specific from left
                       || (pxl_sector_L >= 0 && pxl_sector_R >= 0 && sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == pxl_sector_L && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R) // take one specific from left and one (other) specific sector from right
                       || (pxl_sector_R == -1 && pxl_sector_L == -1) // take all sectors from left and right
                      )
                    {
                        HFT_align_data_vectors[0].push_back(HFT_hit_IL);
                        HFT_align_data_vectors[1].push_back(HFT_hit_OL);
                        HFT_align_data_vectors[2].push_back(HFT_hit_IR);
                        HFT_align_data_vectors[3].push_back(HFT_hit_OR);
                        HFT_align_data_vectors[4].push_back(HFT_dir_L);
                        HFT_align_data_vectors[5].push_back(HFT_dir_R);
                    }

                    if(
                       (pxl_sector_R_orig >= 0 && pxl_sector_L_orig == -1 && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R_orig) // take all sectors from left and one specific from right
                       || (pxl_sector_L_orig >= 0 && pxl_sector_R_orig == -1 && sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == pxl_sector_L_orig) // take all sectors from right and one specific from left
                       || (pxl_sector_L_orig >= 0 && pxl_sector_R_orig >= 0 && sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == pxl_sector_L_orig && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R_orig) // take one specific from left and one (other) specific sector from right
                       || (pxl_sector_R_orig == -1 && pxl_sector_L_orig == -1) // take all sectors from left and right
                      )
                    {
                        if(i_full_sec_align == 0)
                        {
                            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                            {
                                h_pixel_line_dca[0][i_xyz]->Fill(dca_vec_IR[i_xyz]);
                                h_pixel_line_dca[3][i_xyz]->Fill(dca_vec_OR[i_xyz]);
                                h2D_pixel_line_dca_vs_z[0][i_xyz]->Fill(z3_HFT,dca_vec_IR[i_xyz]);
                                h2D_pixel_line_dca_vs_z[3][i_xyz]->Fill(z4_HFT,dca_vec_OR[i_xyz]);
                            }
                        }
                    }
                }
            }
            cout << "Data vectors filled" << endl;
            //------------------------------------------------------------------------------



            //------------------------------------------------------------------------------
            if(flag_min_option == 0 || flag_min_option == 2 || flag_min_option == 3)
            {
                cout << "Set minimizer method" << endl;
                // Choose method upon creation between:
                // kConjugateFR, kConjugatePR, kVectorBFGS,
                // kVectorBFGS2, kSteepestDescent


                ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS ); // <-- good results
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugateFR ); // <-- used
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugatePR );
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kSteepestDescent ); <-
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS2 );

                //ROOT::Math::GSLSimAnMinimizer min;

                //ROOT::Math::GSLMinimizer min( ROOT::Math::kMigradImproved);
                //ROOT::Math::GSLMultiMinimizer min( ROOT::Math::kVectorBFGS ); does not work, header not found
                //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);
                //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
                //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigradImproved );
                //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);

                min.SetMaxFunctionCalls(1000000);
                min.SetMaxIterations(100000);
                min.SetTolerance(0.00001);

                ROOT::Math::Functor f(&Align_HFT,6);

                // sec 5, total_dca = 117.764, delta_x = -0.203844, delta_y = -0.0569418, delta_z =0.0149812, alpha = 0.617132, beta = 0.0050859, gamma = -0.616138
                // sec 5, total_dca = 116.171, delta_x = -0.203654, delta_y = -0.0570076, delta_z =0.0131591, alpha = 0.669738, beta = 0.00503699, gamma = -0.668778
                // sec 5, total_dca = 115.847, delta_x = -0.203405, delta_y = -0.0572426, delta_z =0.0129719, alpha = 0.680852, beta = 0.00506325, gamma = -0.679878

                // all sectors, total_dca = 122.958, delta_x = 0.00347126, delta_y = -0.00156148, delta_z =-0.0185277, alpha = -0.53865, beta = -0.00329105, gamma = 0.538159

                Double_t step[6]     = {0.2,0.2,0.2,0.2,0.2,0.2};
                //Double_t variable[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
                //Double_t variable[6] = {-0.0582198,0.0179723,0.447557,0.00498511,-0.446708};
                //Double_t variable[6] = {-0.203844,-0.0569418,0.0149812,0.617132,0.0050859,-0.616138};
                //Double_t variable[6] = {-0.203654,-0.0570076,0.0131591,0.669738,0.00503699,-0.668778};
                Double_t variable[6] = { -0.203405,-0.0572426,0.0129719,0.680852,0.00506325,-0.679878};
                if(flag_min_option == 2)
                {
                    for(Int_t i_var = 0; i_var < 6; i_var++)
                    {
                        variable[i_var] = 0.0;
                        step[i_var] = 0.1;
                    }
                }

                min.SetFunction(f);

                // Set the free variables to be minimized!
                min.SetVariable(0,"delta_x",variable[0], step[0]);
                min.SetVariable(1,"delta_y",variable[1], step[1]);
                min.SetVariable(2,"delta_z",variable[2], step[2]);
                min.SetVariable(3,"alpha",variable[3], step[3]);
                min.SetVariable(4,"beta",variable[4], step[4]);
                min.SetVariable(5,"gamma",variable[5], step[5]);

                // start minimization
                min.Minimize();

                const Double_t *min_params = min.X();

                for(Int_t i = 0; i < 6; i++)
                {
                    cout << "par = " << i << ", value = " << min_params[i] << endl;
                }

                delta_x   = min_params[0];
                delta_y   = min_params[1];
                delta_z   = min_params[2];
                rot_alpha = min_params[3];
                rot_beta  = min_params[4];
                rot_gamma = min_params[5];
            }

            rot4.SetToIdentity();
            rot4.RotateZ(rot_alpha);
            rot4.RotateX(rot_beta);
            rot4.RotateZ(rot_gamma);

            shift_vec.SetXYZ(delta_x,delta_y,delta_z);

            if(flag_min_option == 1)
            {
                //delta_x   = -0.203405;
                //delta_y   = -0.0572426;
                //delta_z   = 0.0129719;
                //rot_alpha = 0.680852;
                //rot_beta  = 0.00506325;
                //rot_gamma = -0.679878;

                //delta_x   = shift_vec_init.X();
                //delta_y   = shift_vec_init.Y();
                //delta_z   = shift_vec_init.Z();
                //rot_alpha = rot3.PhiZ();
                //rot_beta  = rot3.PhiX();
                //rot_gamma = rot3.PhiZ();

                rot4 = rot1;
                shift_vec = shift_vec_init;
            }
            //------------------------------------------------------------------------------


            vector_rotation_sectors[sec_comb_full_sec_align[1][i_full_sec_align]] = rot4;
            vector_shift_sectors[sec_comb_full_sec_align[1][i_full_sec_align]]    = shift_vec;
        } // end of full alignment loop
        //---------------------------------------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        cout << "Second loop after alignment" << endl;

        const Int_t N_plot_tracks = 50;
        const Int_t N_HFT_hits_per_track = 4;
        TPolyMarker3D PM_Pixel_hit[N_plot_tracks][N_HFT_hits_per_track];
        TPolyLine3D PL_Pixel_track[N_plot_tracks];

        Int_t counter_pxl_tracks = 0;

        if(stop_event_use > file_entries_HFT) stop_event_use = file_entries_HFT;
        for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
        {
            if (counter != 0  &&  counter % 1000 == 0)
                cout << "." << flush;
            if (counter != 0  &&  counter % 10000 == 0)
            {
                if((stop_event_use-start_event_use) > 0)
                {
                    Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                    cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (HFT Alignment) " << flush;
                }
            }

            if (!NT_HFT_left_right->GetEntry( counter )) // take the event -> information is stored in event
                break;  // end of data chunk

            // Detector index
            // 0 = inner pixel left
            // 1 = outer pixel left
            // 2 = IST left
            // 3 = inner pixel right
            // 4 = outer pixel right
            // 5 = IST right
            Int_t index_pxl_IR = get_HFT_det_index(id_IR_HFT);
            Int_t index_pxl_OR = get_HFT_det_index(id_OR_HFT);

            Int_t sector_pxl_IL = get_HFT_pixel_sector(id_IL_HFT,IST_ladder_number);
            Int_t sector_pxl_OL = get_HFT_pixel_sector(id_OL_HFT,IST_ladder_number);

            Int_t sector_pxl_IR = get_HFT_pixel_sector(id_IR_HFT,IST_ladder_number);
            Int_t sector_pxl_OR = get_HFT_pixel_sector(id_OR_HFT,IST_ladder_number);

            // set HFT hits for inner/outer left/right
            HFT_hit_IL.SetXYZ(x1_HFT,y1_HFT,z1_HFT);
            HFT_hit_OL.SetXYZ(x2_HFT,y2_HFT,z2_HFT);
            HFT_hit_IR.SetXYZ(x3_HFT,y3_HFT,z3_HFT);
            HFT_hit_OR.SetXYZ(x4_HFT,y4_HFT,z4_HFT);

            if(flag_min_option == 2 || flag_min_option == 3) // first transform with init parameters
            {
                // apply rotation and translation
                //HFT_hit_IR.RotateZ(rot_alpha_init);
                //HFT_hit_IR.RotateX(rot_beta_init);
                //HFT_hit_IR.RotateZ(rot_gamma_init);
                HFT_hit_IR.Transform(rot1);
                HFT_hit_IR += shift_vec_init;

                //HFT_hit_OR.RotateZ(rot_alpha_init);
                //HFT_hit_OR.RotateX(rot_beta_init);
                //HFT_hit_OR.RotateZ(rot_gamma_init);
                HFT_hit_OR.Transform(rot1);
                HFT_hit_OR += shift_vec_init;
            }

            // apply rotation and translation
            if(flag_min_option != 3) // rot4 -> global rotation matrix
            {
                HFT_hit_IR.Transform(rot4);
                HFT_hit_IR += shift_vec;
                //cout << "shift_vec = {" << shift_vec.X() << ", " << shift_vec.Y() << ", " << shift_vec.Z() << "}" << endl;
                //cout << "|" << rot4.XX() << " " << rot4.XY() << " " << rot4.XZ() << "|" << endl;
                //cout << "|" << rot4.YX() << " " << rot4.YY() << " " << rot4.YZ() << "|" << endl;
                //cout << "|" << rot4.ZX() << " " << rot4.ZY() << " " << rot4.ZZ() << "|" << endl;

                HFT_hit_OR.Transform(rot4);
                HFT_hit_OR += shift_vec;
            }
            if(flag_min_option == 3) // full sector alignment with initial parameters
            {
                HFT_hit_IL.Transform(vector_rotation_sectors[sector_pxl_IL]);
                HFT_hit_IL += vector_shift_sectors[sector_pxl_IL];

                HFT_hit_OL.Transform(vector_rotation_sectors[sector_pxl_OL]);
                HFT_hit_OL += vector_shift_sectors[sector_pxl_OL];

                HFT_hit_IR.Transform(vector_rotation_sectors[sector_pxl_IR]);
                HFT_hit_IR += vector_shift_sectors[sector_pxl_IR];

                HFT_hit_OR.Transform(vector_rotation_sectors[sector_pxl_OR]);
                HFT_hit_OR += vector_shift_sectors[sector_pxl_OR];

                pxl_sector_L = pxl_sector_L_orig;
                pxl_sector_R = pxl_sector_R_orig;
            }

            // calculate direction vectors for left and right halfs
            HFT_dir_L = HFT_hit_IL - HFT_hit_OL;
            HFT_dir_R = HFT_hit_IR - HFT_hit_OR;

            // calculate dca values of left tracklets to right HFT hits for inner/outer
            dca_vec_IR = calculateDCA_vec_StraightToPoint(HFT_hit_IL,HFT_dir_L,HFT_hit_IR); // base,dir,point
            dca_vec_OR = calculateDCA_vec_StraightToPoint(HFT_hit_IL,HFT_dir_L,HFT_hit_OR); // base,dir,point

            if(dca_vec_IR.Mag() < 0.5 && dca_vec_OR.Mag() < 0.5)
            {
                //cout << "index_pxl_IR = " << index_pxl_IR << ", index_pxl_OR = " << index_pxl_OR <<
                //    ", sector_pxl_IR = " << sector_pxl_IR << ", sector_pxl_OR = " << sector_pxl_OR << endl;

                if(
                   //pxl_sector_R >= 0 && index_pxl_IR == 3 && index_pxl_OR == 4 && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R
                   //&& sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == 3
                   (pxl_sector_R >= 0 && pxl_sector_L == -1 && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R) // take all sectors from left and one specific from right
                   || (pxl_sector_L >= 0 && pxl_sector_R == -1 && sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == pxl_sector_L) // take all sectors from right and one specific from left
                   || (pxl_sector_L >= 0 && pxl_sector_R >= 0 && sector_pxl_IL == sector_pxl_OL && sector_pxl_IL == pxl_sector_L && sector_pxl_IR == sector_pxl_OR && sector_pxl_IR == pxl_sector_R) // take one specific from left and one (other) specific sector from right
                   || (pxl_sector_R == -1 && pxl_sector_L == -1) // take all sectors from left and right
                  )
                {
                    if(counter_pxl_tracks < N_plot_tracks)
                    {
                        //cout << "counter = " << counter << ", counter_pxl_tracks = " << counter_pxl_tracks << ", x1_HFT = " << x1_HFT << endl;
                        PM_Pixel_hit[counter_pxl_tracks][0].SetNextPoint(x1_HFT,y1_HFT,z1_HFT);
                        PM_Pixel_hit[counter_pxl_tracks][1].SetNextPoint(x2_HFT,y2_HFT,z2_HFT);
                        PM_Pixel_hit[counter_pxl_tracks][2].SetNextPoint(x3_HFT,y3_HFT,z3_HFT);
                        PM_Pixel_hit[counter_pxl_tracks][3].SetNextPoint(x4_HFT,y4_HFT,z4_HFT);

                        PL_Pixel_track[counter_pxl_tracks].SetNextPoint(x2_HFT,y2_HFT,z2_HFT);
                        PL_Pixel_track[counter_pxl_tracks].SetNextPoint(x4_HFT,y4_HFT,z4_HFT);
                        counter_pxl_tracks++;
                    }

                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        h_pixel_line_dca[1][i_xyz]->Fill(dca_vec_IR[i_xyz]);
                        h_pixel_line_dca[4][i_xyz]->Fill(dca_vec_OR[i_xyz]);

                        h2D_pixel_line_dca_vs_z[1][i_xyz]->Fill(z3_HFT,dca_vec_IR[i_xyz]);
                        h2D_pixel_line_dca_vs_z[4][i_xyz]->Fill(z4_HFT,dca_vec_OR[i_xyz]);
                    }
                }
            }
        }
        //------------------------------------------------------------------------------


#if 0
        //------------------------------------------------------------------------------
        cout << "Draw 3D image" << endl;
        TH3D* h_3D_dummy = new TH3D("h_3D_dummy","h_3D_dummy",200,-300,300,200,-300,300,200,-1000,1000);
        TCanvas* c_3D    = new TCanvas("c_3D","c_3D",10,10,800,800);
        Draw_STAR_3D();

        cout << "counter_pxl_tracks = " << counter_pxl_tracks << endl;
        Int_t pxl_marker_color[4] = {2,3,4,5};
        for(Int_t pxl_track = 0; pxl_track < counter_pxl_tracks; pxl_track++)
        {
            PL_Pixel_track[pxl_track].SetLineStyle(1);
            PL_Pixel_track[pxl_track].SetLineColor(2);
            PL_Pixel_track[pxl_track].SetLineWidth(1);
            PL_Pixel_track[pxl_track].DrawClone("ogl");
            for(Int_t pxl_hit = 0; pxl_hit < N_HFT_hits_per_track; pxl_hit++)
            {
                PM_Pixel_hit[pxl_track][pxl_hit].SetMarkerStyle(20);
                PM_Pixel_hit[pxl_track][pxl_hit].SetMarkerSize(0.8);
                PM_Pixel_hit[pxl_track][pxl_hit].SetMarkerColor(pxl_marker_color[pxl_hit]);
                PM_Pixel_hit[pxl_track][pxl_hit].DrawClone();
            }
        }
        //------------------------------------------------------------------------------
#endif


        //------------------------------------------------------------------------------
        HistName = "c_pixel_line_dca";
        TCanvas* c_pixel_line_dca = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
        c_pixel_line_dca->SetFillColor(10);
        c_pixel_line_dca->SetTopMargin(0.05);
        c_pixel_line_dca->SetBottomMargin(0.15);
        c_pixel_line_dca->SetRightMargin(0.05);
        c_pixel_line_dca->SetLeftMargin(0.15);
        c_pixel_line_dca->Divide(3,2);

        Int_t tc_pixel_line_dca[6] = {1,kAzure-1,kOrange+2,1,kAzure-1,kOrange+2};
        for(Int_t i_cut = 0; i_cut < 6; i_cut++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                h_pixel_line_dca[i_cut][i_xyz]->SetLineColor(tc_pixel_line_dca[i_cut]);
            }
        }

        for(Int_t i_IO = 0; i_IO < 2; i_IO++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetTopMargin(0.05);
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetBottomMargin(0.22);
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetRightMargin(0.05);
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetLeftMargin(0.2);
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetTicks(1,1);
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetGrid(0,0);
                c_pixel_line_dca->cd(i_xyz+1+i_IO*3)->SetFillColor(10);

                h_pixel_line_dca[1+i_IO*3][i_xyz]->SetFillColor(tc_pixel_line_dca[1]);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->SetFillStyle(3001);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->SetStats(0);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->SetTitle("");
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->SetTitleOffset(1.5);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->SetLabelSize(0.055);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->SetLabelSize(0.055);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->SetTitleSize(0.055);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->SetTitleSize(0.055);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->CenterTitle();
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->CenterTitle();
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->SetTitle(xyz_label[i_xyz]);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetYaxis()->SetTitle("counts");
                h_pixel_line_dca[1+i_IO*3][i_xyz]->GetXaxis()->SetRangeUser(-0.03,0.03);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->DrawCopy("");
                h_pixel_line_dca[0+i_IO*3][i_xyz]->SetFillColor(kGray+1);
                h_pixel_line_dca[0+i_IO*3][i_xyz]->SetFillStyle(3001);
                h_pixel_line_dca[0+i_IO*3][i_xyz]->DrawCopy("same");
                h_pixel_line_dca[1+i_IO*3][i_xyz]->DrawCopy("same");

                //h_pixel_line_dca[2][i_xyz+i_IO*3]->DrawCopy("same");

                plotTopLegend((char*)label_IO_pixel[i_IO].Data(),0.72,0.885,0.055,1,0.0,42,1,1);
                plotTopLegend((char*)label_pxl_sector_L.Data(),0.72,0.885-0.06,0.055,1,0.0,42,1,1);
                plotTopLegend((char*)label_pxl_sector_R.Data(),0.72,0.885-0.12,0.055,1,0.0,42,1,1);

                //----------------------------------
                // first fit
                start_fit = -0.02;
                stop_fit  = 0.02;
                for(Int_t x = 0; x < 3; x++)
                {
                    GaussFit->ReleaseParameter(x);
                    GaussFit->SetParError(x,0.0);
                    GaussFit->SetParameter(x,0.0);
                }
                GaussFit->SetParameter(0,h_pixel_line_dca[1+i_IO*3][i_xyz]->GetBinContent(h_pixel_line_dca[1+i_IO*3][i_xyz]->FindBin(0)));
                GaussFit->SetParameter(1,0.0);
                GaussFit->SetParameter(2,0.005);
                GaussFit->SetRange(start_fit,stop_fit);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->Fit("GaussFit","QN","",start_fit,stop_fit);
                height = GaussFit->GetParameter(0);
                mean   = GaussFit->GetParameter(1);
                sigma  = GaussFit->GetParameter(2);
                //----------------------------------



                //----------------------------------
                // second fit
                start_fit = mean-2.0*sigma;
                stop_fit  = mean+2.0*sigma;
                for(Int_t x = 0; x < 3; x++)
                {
                    GaussFit->ReleaseParameter(x);
                    GaussFit->SetParError(x,0.0);
                    GaussFit->SetParameter(x,0.0);
                }
                GaussFit->SetParameter(0,height);
                GaussFit->SetParameter(1,mean);
                GaussFit->SetParameter(2,sigma*0.8);
                GaussFit->SetRange(start_fit,stop_fit);
                h_pixel_line_dca[1+i_IO*3][i_xyz]->Fit("GaussFit","QN","",start_fit,stop_fit);
                GaussFit->SetLineWidth(1);
                GaussFit->SetLineStyle(1);
                GaussFit->SetLineColor(2);
                GaussFit->SetRange(start_fit,stop_fit);
                GaussFit->DrawCopy("same");

                height = GaussFit->GetParameter(0);
                mean   = GaussFit->GetParameter(1);
                sigma  = GaussFit->GetParameter(2);

                HistName = "#mu = ";
                sprintf(NoP,"%2.1f",mean*10000.0);
                HistName += NoP;
                HistName += " #mum";
                plotTopLegend((char*)HistName.Data(),0.25,0.88-0.13,0.055,1,0.0,42,1,1);

                HistName = "#sigma = ";
                sprintf(NoP,"%2.1f",sigma*10000.0);
                HistName += NoP;
                HistName += " #mum";
                plotTopLegend((char*)HistName.Data(),0.245,0.82-0.13,0.055,1,0.0,42,1,1);

                //label_pxl_sector_L
                //----------------------------------



                //----------------------------------
                TLegend* leg_pixel_line_dca = new TLegend(0.23,0.8,0.46,0.93); // x1,y1,x2,y2
                leg_pixel_line_dca->SetBorderSize(0);
                leg_pixel_line_dca->SetFillColor(0);
                leg_pixel_line_dca->SetTextSize(0.055);
                leg_pixel_line_dca->SetTextFont(42);

                leg_pixel_line_dca->AddEntry(h_pixel_line_dca[0+i_IO*3][i_xyz],"before align.","fl");
                leg_pixel_line_dca->AddEntry(h_pixel_line_dca[1+i_IO*3][i_xyz],"after align.","fl");
                leg_pixel_line_dca->Draw();
                //----------------------------------


            }
        }
        //------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        HistName = "c2D_pixel_line_dca_vs_z";
        TCanvas* c2D_pixel_line_dca_vs_z = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
        c2D_pixel_line_dca_vs_z->SetFillColor(10);
        c2D_pixel_line_dca_vs_z->SetTopMargin(0.05);
        c2D_pixel_line_dca_vs_z->SetBottomMargin(0.15);
        c2D_pixel_line_dca_vs_z->SetRightMargin(0.05);
        c2D_pixel_line_dca_vs_z->SetLeftMargin(0.15);
        c2D_pixel_line_dca_vs_z->Divide(3,2);

        for(Int_t i_IO = 0; i_IO < 2; i_IO++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetTopMargin(0.05);
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetBottomMargin(0.22);
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetRightMargin(0.05);
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetLeftMargin(0.2);
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetTicks(1,1);
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetGrid(0,0);
                c2D_pixel_line_dca_vs_z->cd(i_xyz+1+i_IO*3)->SetFillColor(10);

                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->SetStats(0);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->SetTitle("");
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetTitleOffset(1.5);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->SetLabelSize(0.055);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetLabelSize(0.055);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->SetTitleSize(0.055);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetTitleSize(0.055);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->CenterTitle();
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->CenterTitle();
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->SetTitle("z_{pxl} (cm)");
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetTitle(xyz_label[i_xyz]);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetXaxis()->SetRangeUser(-11,11);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->GetYaxis()->SetRangeUser(-0.12,0.12);
                h2D_pixel_line_dca_vs_z[1+i_IO*3][i_xyz]->DrawCopy("colz");

                plotTopLegend((char*)label_IO_pixel[i_IO].Data(),0.25,0.885,0.055,1,0.0,42,1,1);
                plotTopLegend((char*)label_pxl_sector_L.Data(),0.25,0.885-0.06,0.055,1,0.0,42,1,1);
                plotTopLegend((char*)label_pxl_sector_R.Data(),0.25,0.885-0.12,0.055,1,0.0,42,1,1);
            }
        }
        //------------------------------------------------------------------------------
    } // end of left right half alignment
    //---------------------------------------------------------------------------------------------------------------




    //---------------------------------------------------------------------------------------------------------------
    // Sector-to-sector PXL alignment
    if(flag_min_option >= 10 && flag_min_option < 20)
    {
        cout << "Doing sector-to-sector PXL alignment" << endl;

        TVector3 PXL_hit_IO_AB[4]; // A/B
        TVector3 PXL_hit_IO_AR[4]; // align inner, align outer, reference inner, reference outer
        TVector3 PXL_hit_IR, PXL_hit_OR, PXL_hit_IA, PXL_hit_OA, dca_vec_IA, dca_vec_OA, PXL_dir_A, PXL_dir_R;
        PXL_sec_align_data_vectors.resize(6); // hit_IR, hit_OR, hit_IA, hit_OA, dir_R, dir_A  (hit_IA == PXL hit inner to be aligned, hit_OR == PXL hit outer reference)

        if(stop_event_use > file_entries_PXL_IST) stop_event_use = file_entries_PXL_IST;

        const Int_t N_full_sec_align = 10;
        Int_t sector_ref_order_array[N_full_sec_align] = {0,4,5,6,1,9,2,7,3,8};

        std::vector<TRotation> vector_rotation_sectors_STS;
        std::vector<TVector3>  vector_shift_sectors_STS;
        vector_rotation_sectors_STS.resize(N_full_sec_align);
        vector_shift_sectors_STS.resize(N_full_sec_align);

        for(Int_t i_full_sec_align = 0; i_full_sec_align < N_full_sec_align; i_full_sec_align++) // loop for full sector alignment, only activated if flag_min_option == 3
        {
            vector_rotation_sectors_STS[i_full_sec_align].SetToIdentity();
            vector_shift_sectors_STS[i_full_sec_align].SetXYZ(0.0,0.0,0.0);
        }


        //-----------------------------------------------------
        cout << "Define histograms" << endl;
        const Int_t N_histos_STS = 4;
        TH1F* h_pixel_line_dca_STS[N_full_sec_align-1][N_histos_STS][3];
        TH2F* h2D_pixel_line_dca_vs_z_STS[N_full_sec_align-1][N_histos_STS][3];

        for(Int_t i_full_sec_align = 0; i_full_sec_align < (N_full_sec_align-1); i_full_sec_align++) // loop for full sector alignment, only activated if flag_min_option == 3
        {
            for(Int_t i_histo = 0; i_histo < N_histos_STS; i_histo++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    HistName = "h_pixel_line_dca_STS_";
                    HistName += i_full_sec_align;
                    HistName += "_hist_";
                    HistName += i_histo;
                    HistName += "_xyz_";
                    HistName += i_xyz;
                    h_pixel_line_dca_STS[i_full_sec_align][i_histo][i_xyz] = new TH1F(HistName.Data(),HistName.Data(),15000,-10,10);

                    HistName = "h2D_pixel_line_dca_vs_z_STS_";
                    HistName += i_full_sec_align;
                    HistName += "_hist_";
                    HistName += i_histo;
                    HistName += "_xyz_";
                    HistName += i_xyz;
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_histo][i_xyz] = new TH2F(HistName.Data(),HistName.Data(),100,-25,25,800,-0.4,0.4);
                }
            }
        }
        //-----------------------------------------------------



        //------------------------------------------------------------------------------
        Int_t N_loop_min_apply = 1;
        if(flag_min_option == 11) N_loop_min_apply = 2; // do alignment -> two loop, first loop is alignment, second loop applies alignment
        for(Int_t i_loop_min_apply = 0; i_loop_min_apply < N_loop_min_apply; i_loop_min_apply++)
        {
            std::vector<Double_t> ref_sector_vector;

            // Loop over all sector combinations for sector-to-sector alignment
            for(Int_t i_full_sec_align = 0; i_full_sec_align < (N_full_sec_align-1); i_full_sec_align++) // loop for full sector alignment, only activated if flag_min_option == 3
            {

                Int_t pxl_sector_align = sector_ref_order_array[i_full_sec_align+1]; // sector next to align
                cout << "Fill data vectors for sector-to-sector alignment: Sector to align = " << pxl_sector_align << endl;
                ref_sector_vector.push_back(sector_ref_order_array[i_full_sec_align]); // add sector to reference sector list

                // clear data vectors
                for(Int_t i_data = 0; i_data < 6; i_data++)
                {
                    PXL_sec_align_data_vectors[i_data].clear();
                }

                for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
                {
                    if (counter != 0  &&  counter % 1000 == 0)
                        cout << "." << flush;
                    if (counter != 0  &&  counter % 10000 == 0)
                    {
                        if((stop_event_use-start_event_use) > 0)
                        {
                            Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                            cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (PXL sector-to-sector alignment) " << flush;
                        }
                    }

                    if (!NT_PXL_IST->GetEntry( counter )) // take the event -> information is stored in event
                        break;  // end of data chunk


                    // PXL sector information
                    Int_t sector_pxl_array[2][4]; // [(A, B hits), (align, ref. hits)][(inner PXL hit A, outer PXL hit A, inner PXL hit B, outer PXL hit B), (inner align, outer align, inner ref, outer ref)]
                    sector_pxl_array[0][0] = get_HFT_pixel_sector(PXL_id_IA,IST_ladder_number); // inner PXL hit A
                    sector_pxl_array[0][1] = get_HFT_pixel_sector(PXL_id_OA,IST_ladder_number); // outer PXL hit A
                    sector_pxl_array[0][2] = get_HFT_pixel_sector(PXL_id_IB,IST_ladder_number); // inner PXL hit B
                    sector_pxl_array[0][3] = get_HFT_pixel_sector(PXL_id_OB,IST_ladder_number); // outer PXL hit B

                    // PXL index information
                    Int_t index_pxl_array[2][4]; // [(A, B hits), (align, ref. hits)][(inner PXL hit A, outer PXL hit A, inner PXL hit B, outer PXL hit B), (inner align, outer align, inner ref, outer ref)]
                    index_pxl_array[0][0] = get_HFT_det_index(PXL_id_IA);
                    index_pxl_array[0][1] = get_HFT_det_index(PXL_id_OA);
                    index_pxl_array[0][2] = get_HFT_det_index(PXL_id_IB);
                    index_pxl_array[0][3] = get_HFT_det_index(PXL_id_OB);

                    // Set PXL vectors for A/B and inner/outer
                    PXL_hit_IO_AB[0].SetXYZ(PXL_IA_x,PXL_IA_y,PXL_IA_z);
                    PXL_hit_IO_AB[1].SetXYZ(PXL_OA_x,PXL_OA_y,PXL_OA_z);
                    PXL_hit_IO_AB[2].SetXYZ(PXL_IB_x,PXL_IB_y,PXL_IB_z);
                    PXL_hit_IO_AB[3].SetXYZ(PXL_OB_x,PXL_OB_y,PXL_OB_z);

#if 0
                    for(Int_t i_hit_dca_matrix = 0; i_hit_dca_matrix < 4; i_hit_dca_matrix++)
                    {
                        cout << "Hit before dca alignment: {" << PXL_hit_IO_AB[i_hit_dca_matrix][0] << ", "
                            << PXL_hit_IO_AB[i_hit_dca_matrix][1] << ", " << PXL_hit_IO_AB[i_hit_dca_matrix][2] << "}" << endl;
                        Double_t local_hit_IO_AB[4] = {PXL_hit_IO_AB[i_hit_dca_matrix][0],PXL_hit_IO_AB[i_hit_dca_matrix][1],PXL_hit_IO_AB[i_hit_dca_matrix][2],1};
                        Double_t master_hit[4];
                        M_dca_matrices[sector_pxl_array[0][i_hit_dca_matrix]]->LocalToMaster(local_hit_IO_AB,master_hit);
                        PXL_hit_IO_AB[i_hit_dca_matrix].SetXYZ(master_hit[0],master_hit[1],master_hit[2]);
                        cout << "Hit after dca alignment: {" << PXL_hit_IO_AB[i_hit_dca_matrix][0] << ", "
                            << PXL_hit_IO_AB[i_hit_dca_matrix][1] << ", " << PXL_hit_IO_AB[i_hit_dca_matrix][2] << "}" << endl;
                    }
#endif

                    Int_t hit_combination = -1;

                    Int_t AB_inner_outer_combinations[4][4] = // [sec combination, inner/outer]
                    {
                        {0,1,2,3}, // IA, OA, IB, OB
                        {2,3,0,1}, // IB, OB, IA, OA
                        {0,3,2,1}, // IA, OB, IB, OA
                        {2,1,0,3}  // IB, OA, IA, OB
                    };

                    // check all possible combinations of hits for align and reference sectors -> check which hits belong together (align sector, reference sector(s))
                    for(Int_t sec_IO_AB_comb = 0; sec_IO_AB_comb < 4; sec_IO_AB_comb++)
                    {
                        if(
                           sector_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][0]] == sector_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][1]] // inner sector == outer sector of this combination
                           && sector_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][0]] == pxl_sector_align // check if hits belong to sector which has to be aligned
                           && std::find(ref_sector_vector.begin(), ref_sector_vector.end(), sector_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][2]]) != ref_sector_vector.end() // check if other hits belong to sectors which belong to reference sectors
                           && std::find(ref_sector_vector.begin(), ref_sector_vector.end(), sector_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][3]]) != ref_sector_vector.end() // check if other hits belong to sectors which belong to reference sectors
                          )
                        {

                            for(Int_t i_IO_AB = 0; i_IO_AB < 4; i_IO_AB++)
                            {
                                PXL_hit_IO_AR[i_IO_AB]       = PXL_hit_IO_AB[AB_inner_outer_combinations[sec_IO_AB_comb][i_IO_AB]];
                                sector_pxl_array[1][i_IO_AB] = sector_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][i_IO_AB]];
                                index_pxl_array[1][i_IO_AB]  = index_pxl_array[0][AB_inner_outer_combinations[sec_IO_AB_comb][i_IO_AB]];

                                if(sector_pxl_array[1][i_IO_AB] >= 5) // right half
                                {
                                    // apply initial rotation and translation (from left right half alignment)
                                    PXL_hit_IO_AR[i_IO_AB].Transform(rot1);
                                    PXL_hit_IO_AR[i_IO_AB] += shift_vec_init;
                                }

                                PXL_hit_IO_AR[i_IO_AB].Transform(vector_rotation_sectors_STS[sector_pxl_array[1][i_IO_AB]]);
                                PXL_hit_IO_AR[i_IO_AB] += vector_shift_sectors_STS[sector_pxl_array[1][i_IO_AB]];
                            }
                            hit_combination = sec_IO_AB_comb;
                        }
                    }

                    // calculate direction vectors for align and refrence sector(s)
                    PXL_dir_A = PXL_hit_IO_AR[0] - PXL_hit_IO_AR[1]; // inner - outer for sector to be aligned
                    PXL_dir_R = PXL_hit_IO_AR[2] - PXL_hit_IO_AR[3]; // inner - outer for reference sector(s) hits

                    // calculate dca values of left tracklets to right HFT hits for inner/outer
                    dca_vec_IA = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AR[2],PXL_dir_R,PXL_hit_IO_AR[0]); // base,dir,point
                    dca_vec_OA = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AR[2],PXL_dir_R,PXL_hit_IO_AR[1]); // base,dir,point

                    if(dca_vec_IA.Mag() < 0.05 && dca_vec_OA.Mag() < 0.05 && IST_hit == 0 && hit_combination != -1)
                    {

                        //cout << "i_full_sec_align = " << i_full_sec_align << ", hit_combination = " << hit_combination << ", counter = " << counter
                        //    << ", PXL_IA_x = " << PXL_IA_x << ", PXL_OA_x = " << PXL_OA_x << ", PXL_IB_x = " << PXL_IB_x << ", PXL_OB_x = "
                        //    << PXL_OB_x << ", IST_x = " << IST_x << ", IST_counter = " << IST_counter << endl;

#if 0
                        cout << "PXL_IA = {" << PXL_IA_x << ", " << PXL_IA_y << ", " << PXL_IA_z << "}"
                            << ", PXL_OA = {" << PXL_OA_x << ", " << PXL_OA_y << ", " << PXL_OA_z << "}"
                            << ", PXL_IB = {" << PXL_IB_x << ", " << PXL_IB_y << ", " << PXL_IB_z << "}"
                            << ", PXL_OB = {" << PXL_OB_x << ", " << PXL_OB_y << ", " << PXL_OB_z << "}"
                            << ", IST = {" << IST_x << ", " << IST_y << ", " << IST_z << "}"
                            << ", dca_IA = " << dca_vec_IA.Mag() << ", dca_OA = " << dca_vec_OA.Mag()
                            << endl;
#endif
                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            h_pixel_line_dca_STS[i_full_sec_align][0+i_loop_min_apply*2][i_xyz]->Fill(dca_vec_IA[i_xyz]);
                            h_pixel_line_dca_STS[i_full_sec_align][1+i_loop_min_apply*2][i_xyz]->Fill(dca_vec_OA[i_xyz]);

                            h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][0+i_loop_min_apply*2][i_xyz]->Fill(PXL_hit_IO_AR[0].Z(),dca_vec_IA[i_xyz]);
                            h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][1+i_loop_min_apply*2][i_xyz]->Fill(PXL_hit_IO_AR[1].Z(),dca_vec_OA[i_xyz]);
                        }

                        // Add hit points to data vector
                        PXL_sec_align_data_vectors[0].push_back(PXL_hit_IO_AR[2]); // reference inner
                        PXL_sec_align_data_vectors[1].push_back(PXL_hit_IO_AR[3]); // reference outer
                        PXL_sec_align_data_vectors[2].push_back(PXL_hit_IO_AR[0]); // align inner
                        PXL_sec_align_data_vectors[3].push_back(PXL_hit_IO_AR[1]); // align outer
                        PXL_sec_align_data_vectors[4].push_back(PXL_dir_R); // dir reference
                        PXL_sec_align_data_vectors[5].push_back(PXL_dir_A); // dir align
                    }
                }


                //------------------------------------------------------------------------------
                if(flag_min_option == 11 && i_loop_min_apply == 0)
                {
                    cout << "Set minimizer method" << endl;
                    // Choose method upon creation between:
                    // kConjugateFR, kConjugatePR, kVectorBFGS,
                    // kVectorBFGS2, kSteepestDescent


                    ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS ); // <-- good results
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugateFR ); // <-- used
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugatePR );
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kSteepestDescent ); <-
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS2 );

                    //ROOT::Math::GSLSimAnMinimizer min;

                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kMigradImproved);
                    //ROOT::Math::GSLMultiMinimizer min( ROOT::Math::kVectorBFGS ); does not work, header not found
                    //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);
                    //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
                    //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigradImproved );
                    //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);

                    min.SetMaxFunctionCalls(1000000);
                    min.SetMaxIterations(100000);
                    min.SetTolerance(0.00001);

                    ROOT::Math::Functor f(&Align_HFT_STS,6);

                    Double_t step[6]     = {0.2,0.2,0.2,0.2,0.2,0.2};
                    Double_t variable[6] = { -0.203405,-0.0572426,0.0129719,0.680852,0.00506325,-0.679878};
                    if(flag_min_option == 11)
                    {
                        for(Int_t i_var = 0; i_var < 6; i_var++)
                        {
                            variable[i_var] = 0.0;
                            step[i_var] = 0.01;
                        }
                    }

                    min.SetFunction(f);

                    // Set the free variables to be minimized!
                    min.SetVariable(0,"delta_x",variable[0], step[0]);
                    min.SetVariable(1,"delta_y",variable[1], step[1]);
                    min.SetVariable(2,"delta_z",variable[2], step[2]);
                    min.SetVariable(3,"alpha",variable[3], step[3]);
                    min.SetVariable(4,"beta",variable[4], step[4]);
                    min.SetVariable(5,"gamma",variable[5], step[5]);

                    // start minimization
                    min.Minimize();

                    const Double_t *min_params = min.X();

                    //for(Int_t i = 0; i < 6; i++)
                    //{
                    //    cout << "par = " << i << ", value = " << min_params[i] << endl;
                    //}

                    delta_x   = min_params[0];
                    delta_y   = min_params[1];
                    delta_z   = min_params[2];
                    rot_alpha = min_params[3];
                    rot_beta  = min_params[4];
                    rot_gamma = min_params[5];


                    rot4.SetToIdentity();
                    rot4.RotateZ(rot_alpha);
                    rot4.RotateX(rot_beta);
                    rot4.RotateZ(rot_gamma);

                    shift_vec.SetXYZ(delta_x,delta_y,delta_z);

                    vector_rotation_sectors_STS[pxl_sector_align] = rot4;
                    vector_shift_sectors_STS[pxl_sector_align]    = shift_vec;

                    cout << "" << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "Alignment parameters for sector = " << pxl_sector_align <<  endl;
                    cout << "rotation matrix scheme" << endl;
                    cout << "| x' |   | xx xy xz | | x |" << endl;
                    cout << "| y' | = | yx yy yz | | y |" << endl;
                    cout << "| z' |   | zx zy zz | | z |" << endl;
                    cout << "|" << rot4.XX() << " " << rot4.XY() << " " << rot4.XZ() << "|" << endl;
                    cout << "|" << rot4.YX() << " " << rot4.YY() << " " << rot4.YZ() << "|" << endl;
                    cout << "|" << rot4.ZX() << " " << rot4.ZY() << " " << rot4.ZZ() << "|" << endl;
                    cout << "shift vector = {" << shift_vec.X() << ", "  << shift_vec.Y() << ", " << shift_vec.Z() << "}" << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "" << endl;
                }
                //------------------------------------------------------------------------------


            }
            // end of sector-to-sector alignment
        }
        //------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        cout << "Plot 1D residual distributions" << endl;

        TCanvas* c_pixel_line_dca_STS[2][3]; // [inner,outer][x,y,z]
        for(Int_t i_IO = 0; i_IO < 2; i_IO++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                HistName = "c_pixel_line_dca_STS_";
                HistName += i_IO;
                HistName += "_";
                HistName += i_xyz;
                c_pixel_line_dca_STS[i_IO][i_xyz] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
                c_pixel_line_dca_STS[i_IO][i_xyz]->SetFillColor(10);
                c_pixel_line_dca_STS[i_IO][i_xyz]->SetTopMargin(0.05);
                c_pixel_line_dca_STS[i_IO][i_xyz]->SetBottomMargin(0.15);
                c_pixel_line_dca_STS[i_IO][i_xyz]->SetRightMargin(0.05);
                c_pixel_line_dca_STS[i_IO][i_xyz]->SetLeftMargin(0.15);
                c_pixel_line_dca_STS[i_IO][i_xyz]->Divide(3,3);


                for(Int_t i_full_sec_align = 0; i_full_sec_align < (N_full_sec_align-1); i_full_sec_align++)
                {
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetTopMargin(0.05);
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetBottomMargin(0.22);
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetRightMargin(0.05);
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetLeftMargin(0.2);
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetTicks(1,1);
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetGrid(0,0);
                    c_pixel_line_dca_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetFillColor(10);

                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->SetStats(0);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->SetTitle("");
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->SetTitleOffset(1.2);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->SetLabelSize(0.065);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->SetLabelSize(0.065);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->SetTitleSize(0.065);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->SetTitleSize(0.065);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->CenterTitle();
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->CenterTitle();
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->SetTitle(xyz_label[i_xyz]);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetYaxis()->SetTitle("counts");
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->GetXaxis()->SetRangeUser(-0.048,0.048);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->SetFillColor(kGray+1);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->SetFillStyle(3001);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz]->DrawCopy("");

                    // after alignment
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->SetFillColor(kAzure-1);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->SetFillStyle(3001);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->SetLineColor(kAzure-1);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->DrawCopy("same");

                    Int_t align_sector = sector_ref_order_array[i_full_sec_align+1];
                    HistName = "sec. ";
                    HistName += (align_sector+1); // hardware counting starts from 1
                    plotTopLegend((char*)HistName.Data(),0.78,0.885,0.055,1,0.0,42,1,1);
                    plotTopLegend((char*)label_IO_pixel[i_IO].Data(),0.78,0.885-0.06,0.055,1,0.0,42,1,1);

                    //----------------------------------
                    // first fit
                    start_fit = -0.02;
                    stop_fit  = 0.02;
                    for(Int_t x = 0; x < 3; x++)
                    {
                        GaussFit->ReleaseParameter(x);
                        GaussFit->SetParError(x,0.0);
                        GaussFit->SetParameter(x,0.0);
                    }
                    GaussFit->SetParameter(0,h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->GetBinContent(h_pixel_line_dca[1+i_IO*3][i_xyz]->FindBin(0)));
                    GaussFit->SetParameter(1,0.0);
                    GaussFit->SetParameter(2,0.005);
                    GaussFit->SetRange(start_fit,stop_fit);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->Fit("GaussFit","QN","",start_fit,stop_fit);
                    height = GaussFit->GetParameter(0);
                    mean   = GaussFit->GetParameter(1);
                    sigma  = GaussFit->GetParameter(2);
                    //----------------------------------



                    //----------------------------------
                    // second fit
                    start_fit = mean-1.8*sigma;
                    stop_fit  = mean+1.8*sigma;
                    for(Int_t x = 0; x < 3; x++)
                    {
                        GaussFit->ReleaseParameter(x);
                        GaussFit->SetParError(x,0.0);
                        GaussFit->SetParameter(x,0.0);
                    }
                    GaussFit->SetParameter(0,height);
                    GaussFit->SetParameter(1,mean);
                    GaussFit->SetParameter(2,sigma*0.8);
                    GaussFit->SetRange(start_fit,stop_fit);
                    h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz]->Fit("GaussFit","QN","",start_fit,stop_fit);
                    GaussFit->SetLineWidth(1);
                    GaussFit->SetLineStyle(1);
                    GaussFit->SetLineColor(2);
                    GaussFit->SetRange(start_fit,stop_fit);
                    GaussFit->DrawCopy("same");

                    height = GaussFit->GetParameter(0);
                    mean   = GaussFit->GetParameter(1);
                    sigma  = GaussFit->GetParameter(2);

                    height_err = GaussFit->GetParError(0);
                    mean_err   = GaussFit->GetParError(1);
                    sigma_err  = GaussFit->GetParError(2);

                    HistName = "#mu = ";
                    sprintf(NoP,"%2.1f",mean*10000.0);
                    HistName += NoP;
                    HistName += "#pm";
                    sprintf(NoP,"%2.1f",mean_err*10000.0);
                    HistName += NoP;
                    HistName += " #mum";
                    plotTopLegend((char*)HistName.Data(),0.25,0.88-0.13,0.055,1,0.0,42,1,1);

                    HistName = "#sigma = ";
                    sprintf(NoP,"%2.1f",sigma*10000.0);
                    HistName += NoP;
                    HistName += "#pm";
                    sprintf(NoP,"%2.1f",sigma_err*10000.0);
                    HistName += NoP;
                    HistName += " #mum";
                    plotTopLegend((char*)HistName.Data(),0.245,0.82-0.13,0.055,1,0.0,42,1,1);
                    //----------------------------------



                    //----------------------------------
                    TLegend* leg_pixel_line_dca = new TLegend(0.23,0.8,0.46,0.93); // x1,y1,x2,y2
                    leg_pixel_line_dca->SetBorderSize(0);
                    leg_pixel_line_dca->SetFillColor(0);
                    leg_pixel_line_dca->SetTextSize(0.055);
                    leg_pixel_line_dca->SetTextFont(42);

                    leg_pixel_line_dca->AddEntry(h_pixel_line_dca_STS[i_full_sec_align][i_IO][i_xyz],"before align.","fl");
                    leg_pixel_line_dca->AddEntry(h_pixel_line_dca_STS[i_full_sec_align][i_IO+2][i_xyz],"after align.","fl");
                    leg_pixel_line_dca->Draw();
                    //----------------------------------
                }



                //----------------------------------
                cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
                HistName = fig_all_output_dir.Data();
                HistName += c_pixel_line_dca_STS[i_IO][i_xyz]->GetName();
                HistName += out_all_format.Data();
                c_pixel_line_dca_STS[i_IO][i_xyz]->SaveAs(HistName.Data(),"");
                cout << "Saved canvas: " << HistName.Data() << endl;
                //----------------------------------

            }
        }
        //------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        cout << "Plot 2D residual distributions" << endl;

        TCanvas* c_pixel_line_dca_vs_z_STS[2][3]; // [inner,outer][x,y,z]
        for(Int_t i_IO = 0; i_IO < 2; i_IO++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                HistName = "c_pixel_line_dca_vs_z_STS_";
                HistName += i_IO;
                HistName += "_";
                HistName += i_xyz;
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->SetFillColor(10);
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->SetTopMargin(0.05);
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->SetBottomMargin(0.15);
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->SetRightMargin(0.05);
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->SetLeftMargin(0.15);
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->Divide(3,3);


                for(Int_t i_full_sec_align = 0; i_full_sec_align < (N_full_sec_align-1); i_full_sec_align++)
                {
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetTopMargin(0.05);
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetBottomMargin(0.22);
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetRightMargin(0.05);
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetLeftMargin(0.2);
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetTicks(1,1);
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetGrid(0,0);
                    c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->cd(i_full_sec_align+1)->SetFillColor(10);

                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->SetStats(0);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->SetTitle("");
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetTitleOffset(1.2);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->SetLabelSize(0.065);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetLabelSize(0.065);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->SetTitleSize(0.065);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetTitleSize(0.065);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->CenterTitle();
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->CenterTitle();
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetTitle(xyz_label[i_xyz]);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->SetTitle("z_{global} (cm)");
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetXaxis()->SetRangeUser(-11,11);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->GetYaxis()->SetRangeUser(-0.048,0.048);
                    h2D_pixel_line_dca_vs_z_STS[i_full_sec_align][i_IO+2][i_xyz]->DrawCopy("colz");


                    Int_t align_sector = sector_ref_order_array[i_full_sec_align+1];
                    HistName = "sec. ";
                    HistName += (align_sector+1); // hardware counting starts from 1
                    plotTopLegend((char*)HistName.Data(),0.78,0.885,0.055,1,0.0,42,1,1);
                    plotTopLegend((char*)label_IO_pixel[i_IO].Data(),0.78,0.885-0.06,0.055,1,0.0,42,1,1);
                }



                //----------------------------------
                cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
                HistName = fig_all_output_dir.Data();
                HistName += c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->GetName();
                HistName += out_all_format.Data();
                c_pixel_line_dca_vs_z_STS[i_IO][i_xyz]->SaveAs(HistName.Data(),"");
                cout << "Saved canvas: " << HistName.Data() << endl;
                //----------------------------------

            }
        }
        //------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        if(flag_min_option == 11)
        {
            cout << "Save alignment parameters to output file" << endl;
            TFile* Outputfile        = new TFile("./Data/Align_params/HFT_STS_align_V1b.root","RECREATE");

            for(Int_t i_full_sec_align = 0; i_full_sec_align < N_full_sec_align; i_full_sec_align++) // loop for full sector alignment, only activated if flag_min_option == 3
            {
                HistName = "Rot_S";
                HistName += i_full_sec_align;
                vector_rotation_sectors_STS[i_full_sec_align].Write(HistName.Data());

                HistName = "Shift_S";
                HistName += i_full_sec_align;
                vector_shift_sectors_STS[i_full_sec_align].Write(HistName.Data());
            }
        }
        //------------------------------------------------------------------------------



    }
    //---------------------------------------------------------------------------------------------------------------



    //---------------------------------------------------------------------------------------------------------------
    // IST or TPC alignment
    if(flag_min_option >= 20 && flag_min_option < 30)
    {
        // flag_min_option
        // 20 -> no alignment
        // 21 -> global alignment (all ladders together)
        // 22 -> ladder-by-ladder alignment
        // 23 -> use global alignment but don't do minimization
        // 25 -> plot TPC residuals
        // 26 -> do TPC alignment
        // 27 -> do TPC alignment (HFT to TPC) and then TPC sector to HFT alignment
        // 28 -> do TPC alignment with initial pre alignment for TPC sector-to-sector alignment

        if(stop_event_use > file_entries_PXL_IST) stop_event_use = file_entries_PXL_IST;
        cout << "Doing IST alignment with " << stop_event_use << " entries" << endl;



        //-----------------------------------------------------
        if(flag_min_option == 28) // open pre alignment
        {
            // This one acts on the HFT hits
            cout << "Open HFT to TPC global alignment parameters" << endl;
            TFile* Inputfile_HFT_to_TPC_glob  = TFile::Open("./Data/Align_params/HFT_to_TPC_global_alignment.root");  // open the file

            HistName = "rot_HFT_to_TPC";
            rot_HFT_to_TPC_load = *(TRotation*)(Inputfile_HFT_to_TPC_glob->Get(HistName.Data()));
            HistName = "shift_vec_HFT_to_TPC";
            shift_vec_HFT_to_TPC_load = *(TVector3*)(Inputfile_HFT_to_TPC_glob->Get(HistName.Data()));

            // This one acts on the TPC tracks
            for(Int_t i_sec = 0; i_sec < N_TPC_sectors; i_sec++)
            {
                HistName = "./Data/Align_params/";
                HistName += "TPC_sec";
                HistName += i_sec+1;
                HistName += "_to_HFT_align.root";
                cout << "Open TPC sector alignment file: " << HistName.Data() << endl;
                Inputfile_TPC_sector_alignment[i_sec] = TFile::Open(HistName.Data());  // open the file

                HistName = "Rot_matrix";
                rot_TPC_sector_load[i_sec] = *(TRotation*)(Inputfile_TPC_sector_alignment[i_sec]->Get(HistName.Data()));

                HistName = "Shift_vec";
                shift_vec_TPC_sector_load[i_sec] = *(TVector3*)(Inputfile_TPC_sector_alignment[i_sec]->Get(HistName.Data()));
            }
        }
        //-----------------------------------------------------



        //-----------------------------------------------------
        cout << "Define histograms" << endl;
        const Int_t N_IST_ladder = 24;
        const Int_t N_histos_IST = 4;
        TH1F* h_pixel_line_dca_IST[2][N_histos_IST][3];
        TH2F* h2D_pixel_line_dca_vs_z_IST[2][N_histos_IST][3];
        TH2F* h2D_pixel_line_dca_vs_z_IST_ladder[2][N_IST_ladder][3];
        TH1F* h_TPC_track_eta[2];
        TH1F* h_TPC_track_phi[2];

        for(Int_t i_track = 0; i_track < 2; i_track++)
        {
            HistName = "h_TPC_track_eta";
            HistName += "_track_";
            HistName += i_track;
            h_TPC_track_eta[i_track] = new TH1F(HistName.Data(),HistName.Data(),100,-1.5,1.5);

            HistName = "h_TPC_track_phi";
            HistName += "_track_";
            HistName += i_track;
            h_TPC_track_phi[i_track] = new TH1F(HistName.Data(),HistName.Data(),100,-180.0,180.0);
        }


        for(Int_t i_align = 0; i_align < 2; i_align++)
        {
            for(Int_t i_histo = 0; i_histo < N_histos_IST; i_histo++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    HistName = "h_pixel_line_dca_IST";
                    HistName += "_align_";
                    HistName += i_align;
                    HistName += "_hist_";
                    HistName += i_histo;
                    HistName += "_xyz_";
                    HistName += i_xyz;
                    h_pixel_line_dca_IST[i_align][i_histo][i_xyz] = new TH1F(HistName.Data(),HistName.Data(),15000,-10,10);

                    HistName = "h2D_pixel_line_dca_vs_z_IST";
                    HistName += "_align_";
                    HistName += i_align;
                    HistName += "_hist_";
                    HistName += i_histo;
                    HistName += "_xyz_";
                    HistName += i_xyz;
                    h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz] = new TH2F(HistName.Data(),HistName.Data(),100,-25,25,4000,-0.4,0.4);
                }
            }
        }
        for(Int_t i_align = 0; i_align < 2; i_align++)
        {
            for(Int_t i_histo = 0; i_histo < N_IST_ladder; i_histo++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    HistName = "h2D_pixel_line_dca_vs_z_IST_ladder";
                    HistName += "_align_";
                    HistName += i_align;
                    HistName += "_hist_";
                    HistName += i_histo;
                    HistName += "_xyz_";
                    HistName += i_xyz;
                    h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz] = new TH2F(HistName.Data(),HistName.Data(),50,-25,25,400,-0.4,0.4);
                }
            }
        }
        //-----------------------------------------------------



        TVector3 PXL_hit_IO_AB[5]; // A/B
        TVector3 PXL_hit_IA, PXL_hit_OA, dca_vec_IA, dca_vec_OA, PXL_dir_A, PXL_dir_B, dca_vec_IST_A, dca_vec_IST_B;
        IST_align_data_vectors.resize(7); // hit_IA, hit_OA, hit_IB, hit_OB, dir_A, dir_B, hit_IST
        TPC_align_data_vectors.resize(2); // trackA, trackB


        // PXL alignment parameters
        const Int_t N_full_sec_align = 10;
        std::vector<TRotation> vector_rotation_sectors_STS;
        std::vector<TVector3>  vector_shift_sectors_STS;
        vector_rotation_sectors_STS.resize(N_full_sec_align);
        vector_shift_sectors_STS.resize(N_full_sec_align);

        TRotation IST_glob_rot;
        TVector3  IST_glob_shift;
        IST_glob_rot.SetToIdentity();
        IST_glob_shift.SetXYZ(0.0,0.0,0.0);

        std::vector<TRotation> vector_rotation_sectors_IST;
        std::vector<TVector3>  vector_shift_sectors_IST;
        vector_rotation_sectors_IST.resize(N_IST_ladder);
        vector_shift_sectors_IST.resize(N_IST_ladder);

        for(Int_t i_IST_ladder = 0; i_IST_ladder < N_IST_ladder; i_IST_ladder++)
        {
            vector_rotation_sectors_IST[i_IST_ladder].SetToIdentity();
            vector_shift_sectors_IST[i_IST_ladder].SetXYZ(0.0,0.0,0.0);
        }

        cout << "Open PXL alignment parameters from file" << endl;
        TFile* Inputfile  = TFile::Open("./Data/Align_params/HFT_STS_align_V1.root");  // open the file

        for(Int_t i_full_sec_align = 0; i_full_sec_align < N_full_sec_align; i_full_sec_align++)
        {
            HistName = "Rot_S";
            HistName += i_full_sec_align;
            vector_rotation_sectors_STS[i_full_sec_align] = *(TRotation*)(Inputfile->Get(HistName.Data()));

            HistName = "Shift_S";
            HistName += i_full_sec_align;
            vector_shift_sectors_STS[i_full_sec_align] = *(TVector3*)(Inputfile->Get(HistName.Data()));
        }

        if(flag_min_option == 22 || flag_min_option == 23)
        {
            cout << "Open IST global alignment parameters" << endl;
            TFile* Inputfile_IST_glob  = TFile::Open("./Data/Align_params/IST_global_align_V1.root");  // open the file

            HistName = "Rot_IST_glob";
            IST_glob_rot = *(TRotation*)(Inputfile_IST_glob->Get(HistName.Data()));
            HistName = "Shift_IST_glob";
            IST_glob_shift = *(TVector3*)(Inputfile_IST_glob->Get(HistName.Data()));

            cout << "" << endl;
            cout << "------------------------------------------------" << endl;
            cout << "Initial global IST alignment" << endl;
            cout << "rotation matrix scheme" << endl;
            cout << "| x' |   | xx xy xz | | x |" << endl;
            cout << "| y' | = | yx yy yz | | y |" << endl;
            cout << "| z' |   | zx zy zz | | z |" << endl;
            cout << "|" << IST_glob_rot.XX() << " " << IST_glob_rot.XY() << " " << IST_glob_rot.XZ() << "|" << endl;
            cout << "|" << IST_glob_rot.YX() << " " << IST_glob_rot.YY() << " " << IST_glob_rot.YZ() << "|" << endl;
            cout << "|" << IST_glob_rot.ZX() << " " << IST_glob_rot.ZY() << " " << IST_glob_rot.ZZ() << "|" << endl;
            cout << "shift vector = {" << IST_glob_shift.X() << ", "  << IST_glob_shift.Y() << ", " << IST_glob_shift.Z() << "}" << endl;
            cout << "------------------------------------------------" << endl;
            cout << "" << endl;
        }

        Int_t N_align = 1;
        if(flag_min_option == 21 || flag_min_option == 22 || flag_min_option == 26 || flag_min_option == 27 || flag_min_option == 28) N_align = 2;



        cout << "" << endl;
        cout << "------------------------------------------------" << endl;
        cout << "used rotation matrix scheme for PXL (left to right)" << endl;
        cout << "| x' |   | xx xy xz | | x |" << endl;
        cout << "| y' | = | yx yy yz | | y |" << endl;
        cout << "| z' |   | zx zy zz | | z |" << endl;
        cout << "|" << rot1.XX() << " " << rot1.XY() << " " << rot1.XZ() << "|" << endl;
        cout << "|" << rot1.YX() << " " << rot1.YY() << " " << rot1.YZ() << "|" << endl;
        cout << "|" << rot1.ZX() << " " << rot1.ZY() << " " << rot1.ZZ() << "|" << endl;
        cout << "shift vector = {" << shift_vec_init.X() << ", "  << shift_vec_init.Y() << ", " << shift_vec_init.Z() << "}" << endl;
        cout << "------------------------------------------------" << endl;
        cout << "" << endl;


        for(Int_t i_PXL_sec = 0; i_PXL_sec < 10; i_PXL_sec++)
        {
            cout << "" << endl;
            cout << "------------------------------------------------" << endl;
            cout << "used rotation matrix scheme for PXL (sector-to-sector) for sector " << i_PXL_sec << endl;
            cout << "| x' |   | xx xy xz | | x |" << endl;
            cout << "| y' | = | yx yy yz | | y |" << endl;
            cout << "| z' |   | zx zy zz | | z |" << endl;
            cout << "|" << vector_rotation_sectors_STS[i_PXL_sec].XX() << " " << vector_rotation_sectors_STS[i_PXL_sec].XY() << " " << vector_rotation_sectors_STS[i_PXL_sec].XZ() << "|" << endl;
            cout << "|" << vector_rotation_sectors_STS[i_PXL_sec].YX() << " " << vector_rotation_sectors_STS[i_PXL_sec].YY() << " " << vector_rotation_sectors_STS[i_PXL_sec].YZ() << "|" << endl;
            cout << "|" << vector_rotation_sectors_STS[i_PXL_sec].ZX() << " " << vector_rotation_sectors_STS[i_PXL_sec].ZY() << " " << vector_rotation_sectors_STS[i_PXL_sec].ZZ() << "|" << endl;
            cout << "shift vector = {" << vector_shift_sectors_STS[i_PXL_sec].X() << ", "  << vector_shift_sectors_STS[i_PXL_sec].Y() << ", " << vector_shift_sectors_STS[i_PXL_sec].Z() << "}" << endl;
            cout << "------------------------------------------------" << endl;
            cout << "" << endl;
        }



        for(Int_t i_align = 0; i_align < N_align; i_align++) // first loop -> alignment, second loop -> apply alignment
        {
            Int_t track_counter = 0;
            Int_t N_IST_ladder_use = 1;
            if(flag_min_option == 22) N_IST_ladder_use = N_IST_ladder;
            for(Int_t i_IST_ladder = 0; i_IST_ladder < N_IST_ladder_use; i_IST_ladder++) // activates only for flag_min_option == 22
            {
                cout << "i_IST_ladder = " << i_IST_ladder << endl;
                // clear data vectors
                for(Int_t i_data = 0; i_data < 7; i_data++)
                {
                    IST_align_data_vectors[i_data].clear();
                }
                for(Int_t i_tpc_track = 0; i_tpc_track < 2; i_tpc_track++)
                {
                    TPC_align_data_vectors[i_tpc_track].clear();
                }

                for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
                {
                    if (counter != 0  &&  counter % 1000 == 0)
                        cout << "." << flush;
                    if (counter != 0  &&  counter % 10000 == 0)
                    {
                        if((stop_event_use-start_event_use) > 0)
                        {
                            Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                            cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (IST or TPC alignment) " << flush;
                        }
                    }

                    //cout << "counter = " << counter << ", out of " << stop_event_use << endl;
                    if (!NT_PXL_IST->GetEntry( counter )) // take the event -> information is stored in event
                        break;  // end of data chunk

                    // PXL sector information
                    Int_t sector_pxl_array[5]; // [IA,OA,IB,OB,IST]
                    sector_pxl_array[0] = get_HFT_pixel_sector(PXL_id_IA,IST_ladder_number); // inner PXL hit A
                    sector_pxl_array[1] = get_HFT_pixel_sector(PXL_id_OA,IST_ladder_number); // outer PXL hit A
                    sector_pxl_array[2] = get_HFT_pixel_sector(PXL_id_IB,IST_ladder_number); // inner PXL hit B
                    sector_pxl_array[3] = get_HFT_pixel_sector(PXL_id_OB,IST_ladder_number); // outer PXL hit B
                    sector_pxl_array[4] = -1; // IST

                    // PXL index information
                    Int_t index_pxl_array[5]; // [IA,OA,IB,OB,IST]
                    index_pxl_array[0] = get_HFT_det_index(PXL_id_IA);
                    index_pxl_array[1] = get_HFT_det_index(PXL_id_OA);
                    index_pxl_array[2] = get_HFT_det_index(PXL_id_IB);
                    index_pxl_array[3] = get_HFT_det_index(PXL_id_OB);
                    index_pxl_array[4] = -1;

                    // Set PXL vectors for A/B and inner/outer
                    PXL_hit_IO_AB[0].SetXYZ(PXL_IA_x,PXL_IA_y,PXL_IA_z);
                    PXL_hit_IO_AB[1].SetXYZ(PXL_OA_x,PXL_OA_y,PXL_OA_z);
                    PXL_hit_IO_AB[2].SetXYZ(PXL_IB_x,PXL_IB_y,PXL_IB_z);
                    PXL_hit_IO_AB[3].SetXYZ(PXL_OB_x,PXL_OB_y,PXL_OB_z);
                    PXL_hit_IO_AB[4].SetXYZ(IST_x,IST_y,IST_z);

                    if(flag_min_option == 28) // apply pre alignment for HFT hits
                    {
                        for(Int_t HFT_hit = 0; HFT_hit < 5; HFT_hit++)
                        {
                            PXL_hit_IO_AB[HFT_hit].Transform(rot_HFT_to_TPC_load);
                            PXL_hit_IO_AB[HFT_hit] += shift_vec_HFT_to_TPC_load;
                        }
                    }


                    if((flag_min_option == 26 || flag_min_option == 27 || flag_min_option == 28) && i_align == 1) // apply alignment for second loop for HFT to TPC alignment
                    {
                        for(Int_t i_hit = 0; i_hit < 5; i_hit++)
                        {
                            PXL_hit_IO_AB[i_hit].Transform(rot4);
                            PXL_hit_IO_AB[i_hit] += shift_vec;
                        }
                    }

                    // TPC track information, one cosmic creates two tracks (A,B)
                    StThreeVectorF vec_gMomA, vec_OriginA, vec_gMomB, vec_OriginB;
                    StPhysicalHelixD helixA, helixB;
                    Double_t phiA,phiB,pA,pB,etaA,etaB;

                    if(flag_min_option >= 25 && flag_min_option < 30)
                    {
                        vec_gMomA.set(track_pxA,track_pyA,track_pzA);
                        vec_OriginA.set(track_oxA,track_oyA,track_ozA);
                        helixA = StPhysicalHelixD(vec_gMomA,vec_OriginA,magfac,track_qA);

                        vec_gMomB.set(track_pxB,track_pyB,track_pzB);
                        vec_OriginB.set(track_oxB,track_oyB,track_ozB);
                        helixB = StPhysicalHelixD(vec_gMomB,vec_OriginB,magfac,track_qB);

                        phiA = helixA.cat(0).phi()*TMath::RadToDeg();
                        phiB = helixB.cat(0).phi()*TMath::RadToDeg();

                        StThreeVectorF vecA = helixA.at(0);
                        StThreeVectorF vecB = helixB.at(0);

                        pA   = helixA.momentum((Double_t)magfac).mag();
                        pB   = helixB.momentum((Double_t)magfac).mag();

                        StThreeVectorF vectoratsA  = helixA.cat(0);
                        Float_t etaA = vectoratsA.pseudoRapidity();

                        StThreeVectorF vectoratsB  = helixB.cat(0);
                        Float_t etaB = vectoratsB.pseudoRapidity();

                        h_TPC_track_eta[0] ->Fill(etaA);
                        h_TPC_track_eta[1] ->Fill(etaB);

                        h_TPC_track_phi[0] ->Fill(phiA);
                        h_TPC_track_phi[1] ->Fill(phiB);

                        //cout << "counter = " << counter << ", phiA = " << phiA << ", phiB = " << phiB << ", etaA = " << etaA << ", etaB = " << etaB << endl;
                        //cout << "vecA = {" << vecA.x() << ", " << vecA.y() << ", " << vecA.z() << "}" << endl;
                        //cout << "vecB = {" << vecB.x() << ", " << vecB.y() << ", " << vecB.z() << "}" << endl;
                    }


                    if(i_align == 1 && flag_min_option == 21) // after alignment -> apply alignment from first loop
                    {
                        PXL_hit_IO_AB[4].Transform(rot4);
                        PXL_hit_IO_AB[4] += shift_vec;
                    }

                    if(flag_min_option == 22 || flag_min_option == 23) // apply global IST alignment
                    {
                        PXL_hit_IO_AB[4].Transform(IST_glob_rot);
                        PXL_hit_IO_AB[4] += IST_glob_shift;

                        if(i_align == 1)
                        {
                            get_HFT_pixel_sector(IST_id,IST_ladder_number); // IST hit
                            if(IST_ladder_number != -1)
                            {
                                PXL_hit_IO_AB[4].Transform(vector_rotation_sectors_IST[IST_ladder_number]);
                                PXL_hit_IO_AB[4] += vector_shift_sectors_IST[IST_ladder_number];
                            }
                        }
                    }


                    if(flag_min_option < 25) // apply PXL alignment
                    {
                        // Apply alignment for PXL
                        for(Int_t i_PXL_hit = 0; i_PXL_hit < 4; i_PXL_hit++)
                        {
                            if(sector_pxl_array[i_PXL_hit] >= 5) // right half
                            {
                                // apply initial rotation and translation (from left right half alignment)
                                PXL_hit_IO_AB[i_PXL_hit].Transform(rot1);
                                PXL_hit_IO_AB[i_PXL_hit] += shift_vec_init;
                            }

                            PXL_hit_IO_AB[i_PXL_hit].Transform(vector_rotation_sectors_STS[sector_pxl_array[i_PXL_hit]]);
                            PXL_hit_IO_AB[i_PXL_hit] += vector_shift_sectors_STS[sector_pxl_array[i_PXL_hit]];
                        }
                    }


                    // calculate direction vectors for align and refrence sector(s)
                    PXL_dir_A = PXL_hit_IO_AB[0] - PXL_hit_IO_AB[1]; // inner - outer for A
                    PXL_dir_B = PXL_hit_IO_AB[2] - PXL_hit_IO_AB[3]; // inner - outer for B

                    // calculate dca values of B tracklets to A HFT hits for inner/outer
                    dca_vec_IA = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[2],PXL_dir_B,PXL_hit_IO_AB[0]); // base,dir,point
                    dca_vec_OA = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[2],PXL_dir_B,PXL_hit_IO_AB[1]); // base,dir,point

                    // calculate dca values of A/B tracklets to IST hit
                    dca_vec_IST_A = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[0],PXL_dir_A,PXL_hit_IO_AB[4]); // base,dir,point
                    dca_vec_IST_B = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[2],PXL_dir_B,PXL_hit_IO_AB[4]); // base,dir,point

                    //if(i_align == 1) cout << "dca_IST_A = " << dca_vec_IST_A.Mag() << ", dca_IST_B = " << dca_vec_IST_B.Mag() << endl;

                    //if(dca_vec_IA.Mag() < 0.05 && dca_vec_OA.Mag() < 0.05 && IST_counter > -1 && dca_vec_IST_A.Mag() < 0.4 && dca_vec_IST_B.Mag() < 0.4)
                    if(
                       dca_vec_IA.Mag() < 0.05
                       && dca_vec_OA.Mag() < 0.05
                       && IST_counter > -1
                       && dca_vec_IST_A.Mag() < 0.4
                       && dca_vec_IST_B.Mag() < 0.4
                      )
                    {

#if 0
                        cout << "PXL_IA = {" << PXL_IA_x << ", " << PXL_IA_y << ", " << PXL_IA_z << "}"
                            << ", PXL_OA = {" << PXL_OA_x << ", " << PXL_OA_y << ", " << PXL_OA_z << "}"
                            << ", PXL_IB = {" << PXL_IB_x << ", " << PXL_IB_y << ", " << PXL_IB_z << "}"
                            << ", PXL_OB = {" << PXL_OB_x << ", " << PXL_OB_y << ", " << PXL_OB_z << "}"
                            << ", IST = {" << IST_x << ", " << IST_y << ", " << IST_z << "}"
                            << ", dca_IA = " << dca_vec_IA.Mag() << ", dca_OA = " << dca_vec_OA.Mag() << ", IST_counter = " << IST_counter
                            << endl;
#endif
                        get_HFT_pixel_sector(IST_id,IST_ladder_number); // IST hit
                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            h_pixel_line_dca_IST[i_align][0][i_xyz]->Fill(dca_vec_IA[i_xyz]);
                            h_pixel_line_dca_IST[i_align][1][i_xyz]->Fill(dca_vec_OA[i_xyz]);
                            h_pixel_line_dca_IST[i_align][2][i_xyz]->Fill(dca_vec_IST_A[i_xyz]);
                            h_pixel_line_dca_IST[i_align][3][i_xyz]->Fill(dca_vec_IST_B[i_xyz]);
                            h2D_pixel_line_dca_vs_z_IST[i_align][0][i_xyz]->Fill(PXL_hit_IO_AB[0].Z(),dca_vec_IA[i_xyz]);
                            h2D_pixel_line_dca_vs_z_IST[i_align][1][i_xyz]->Fill(PXL_hit_IO_AB[1].Z(),dca_vec_OA[i_xyz]);
                            h2D_pixel_line_dca_vs_z_IST[i_align][2][i_xyz]->Fill(PXL_hit_IO_AB[4].Z(),dca_vec_IST_A[i_xyz]);
                            h2D_pixel_line_dca_vs_z_IST[i_align][3][i_xyz]->Fill(PXL_hit_IO_AB[4].Z(),dca_vec_IST_B[i_xyz]);

                            if(IST_ladder_number != -1)
                            {
                                h2D_pixel_line_dca_vs_z_IST_ladder[i_align][IST_ladder_number][i_xyz]->Fill(PXL_hit_IO_AB[4].Z(),dca_vec_IST_A[i_xyz]);
                            }
                        }

                        if(flag_min_option != 22 || (flag_min_option == 22 && IST_ladder_number == i_IST_ladder))
                        {
                            if( !(flag_min_option >= 25 && flag_min_option < 30) || ((flag_min_option >= 25 && flag_min_option < 30) && pA > min_TPC_track_mom && pB > min_TPC_track_mom))
                            {

                                Int_t flag_good_tpc_track = 1;
                                Double_t dca_array[N_HFT_TPC_match][3][2]; // [HFT hit][x,y,z][first track, second track]
                                for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                                {
                                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                                    {
                                        for(Int_t i_track = 0; i_track < 2; i_track++)
                                        {
                                            dca_array[i_histo][i_xyz][i_track] = 0.0;
                                        }
                                    }
                                }

                                if(flag_min_option >= 25 && flag_min_option < 30)
                                {
                                    for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                                    {
                                        StThreeVectorF pixel_hit, helix_pointA, helix_pointB, dca_vecA, dca_vecB;
                                        pixel_hit.set(PXL_hit_IO_AB[i_histo].X(),PXL_hit_IO_AB[i_histo].Y(),PXL_hit_IO_AB[i_histo].Z());

                                        Float_t pathA,dcaA;
                                        fHelixAtoPointdca(pixel_hit,helixA,pathA,dcaA);
                                        helix_pointA = helixA.at(pathA);

                                        Float_t pathB,dcaB;
                                        fHelixAtoPointdca(pixel_hit,helixB,pathB,dcaB);
                                        helix_pointB = helixB.at(pathB);

                                        dca_vecA = helix_pointA;
                                        dca_vecA -= pixel_hit;

                                        dca_vecB = helix_pointB;
                                        dca_vecB -= pixel_hit;
                                        //cout << "tree counter = " << counter << ", HFT point = " << i_histo << ", dcaA = " << dcaA
                                        //    << ", hit = {" << PXL_hit_IO_AB[i_histo].X() << ", " << PXL_hit_IO_AB[i_histo].Y() << ", " << PXL_hit_IO_AB[i_histo].z() << "}" << endl;
                                        if(dcaA > 0.5 || dcaB > 0.5) flag_good_tpc_track = 0;

                                        if(flag_min_option == 28) // with pre alignment -> all hit points from a track have to be from one sector
                                        {
                                            if(
                                               !(
                                                 secA >= 1 && secA <= 24
                                                 && secB >= 1 && secB <= 24
                                                 && (NsecA + NnohitsecA) == 44
                                                 && (NsecB + NnohitsecB) == 44 // total number of hits + no hit from the same sector is 44 (means 45 hits in total)
                                                 && NsecA > 30
                                                 && NsecB > 30
                                                )
                                              )
                                            {
                                                flag_good_tpc_track = 0;
                                            }
                                        }

                                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                                        {
                                            dca_array[i_histo][i_xyz][0] = dca_vecA[i_xyz];
                                            dca_array[i_histo][i_xyz][1] = dca_vecB[i_xyz];
                                        }
                                    }
                                }

                                if(flag_good_tpc_track)
                                {
                                    IST_align_data_vectors[0].push_back(PXL_hit_IO_AB[0]);
                                    IST_align_data_vectors[1].push_back(PXL_hit_IO_AB[1]);
                                    IST_align_data_vectors[2].push_back(PXL_hit_IO_AB[2]);
                                    IST_align_data_vectors[3].push_back(PXL_hit_IO_AB[3]);
                                    IST_align_data_vectors[4].push_back(PXL_dir_A);
                                    IST_align_data_vectors[5].push_back(PXL_dir_B);
                                    IST_align_data_vectors[6].push_back(PXL_hit_IO_AB[4]);

                                    for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                                    {
                                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                                        {
                                            for(Int_t i_track = 0; i_track < 2; i_track++)
                                            {
                                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Fill(dca_array[i_histo][i_xyz][i_track]);
                                            }
                                        }
                                    }
                                    if(flag_min_option >= 25 && flag_min_option < 30)
                                    {

                                        StPhysicalHelixD helixAB[2] = {helixA,helixB};
                                        Double_t mag_helix[2]       = {magfac,magfac};
                                        Double_t q_helix[2]         = {track_qA,track_qB};
                                        Int_t TPC_sector_helix[2]   = {secA-1,secB-1};
                                        if(flag_min_option == 28) // apply pre alignment for TPC
                                        {
                                            // Apply TPC sector alignment
                                            for(Int_t i_helix = 0; i_helix < 2; i_helix++)
                                            {
                                                StThreeVectorD vec_TPC_Origin = helixAB[i_helix].origin();
                                                StThreeVectorD vec_TPC_Mom    = helixAB[i_helix].momentum(mag_helix[i_helix]);
                                                Int_t          TPC_charge     = helixAB[i_helix].charge(q_helix[i_helix]);

                                                TVector3 TV3_vec_TPC_Mom, TV3_vec_TPC_Origin;
                                                TV3_vec_TPC_Mom.SetXYZ(vec_TPC_Mom.x(),vec_TPC_Mom.y(),vec_TPC_Mom.z());
                                                TV3_vec_TPC_Mom.Transform(rot_TPC_sector_load[TPC_sector_helix[i_helix]]);
                                                vec_TPC_Mom.set(TV3_vec_TPC_Mom.x(),TV3_vec_TPC_Mom.y(),TV3_vec_TPC_Mom.z());
                                                TV3_vec_TPC_Origin.SetXYZ(vec_TPC_Origin.x(),vec_TPC_Origin.y(),vec_TPC_Origin.z());
                                                TV3_vec_TPC_Origin += shift_vec_TPC_sector_load[TPC_sector_helix[i_helix]];
                                                vec_TPC_Origin.set(TV3_vec_TPC_Origin.x(),TV3_vec_TPC_Origin.y(),TV3_vec_TPC_Origin.z());

                                                helixAB[i_helix] = StPhysicalHelixD(vec_TPC_Mom,vec_TPC_Origin,mag_helix[i_helix],TPC_charge);
                                            }
                                        }

                                        TPC_align_data_vectors[0].push_back(helixAB[0]);
                                        TPC_align_data_vectors[1].push_back(helixAB[1]);

                                        //cout << "accepted track_counter = " << track_counter << ", tree counter = " << counter << ", phiA = " << phiA << ", phiB = " << phiB << ", pA = " << pA << ", pB = " << pB << endl;

                                        for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                                        {
                                            Double_t dca_TPC[2];
                                            for(Int_t i_track = 0; i_track < 2; i_track++)
                                            {
                                                dca_TPC[i_track] = TMath::Sqrt(dca_array[i_histo][0][i_track]*dca_array[i_histo][0][i_track] + dca_array[i_histo][1][i_track]*dca_array[i_histo][1][i_track] + dca_array[i_histo][2][i_track]*dca_array[i_histo][2][i_track]);
                                                //cout << "i_track  = " << i_track << ", PXL = " << i_histo << ", dca_TPC = " << dca_TPC[i_track] << endl;
                                            }
                                        }

                                        track_counter++;
                                    }
                                }
                            }
                        }
                    }
                }


                //------------------------------------------------------------------------------
                if((flag_min_option == 21 || flag_min_option == 22) && i_align == 0)
                {
                    cout << "Set minimizer method" << endl;
                    // Choose method upon creation between:
                    // kConjugateFR, kConjugatePR, kVectorBFGS,
                    // kVectorBFGS2, kSteepestDescent


                    ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS ); // <-- good results
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugateFR ); // <-- used
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugatePR );
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kSteepestDescent ); <-
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS2 );

                    //ROOT::Math::GSLSimAnMinimizer min;

                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kMigradImproved);
                    //ROOT::Math::GSLMultiMinimizer min( ROOT::Math::kVectorBFGS ); does not work, header not found
                    //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);
                    //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
                    //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigradImproved );
                    //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);

                    min.SetMaxFunctionCalls(1000000);
                    min.SetMaxIterations(100000);
                    min.SetTolerance(0.0001);

                    ROOT::Math::Functor f(&Align_HFT_IST,6);

                    Double_t step[6]     = {0.2,0.2,0.2,0.2,0.2,0.2};
                    Double_t variable[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
                    for(Int_t i_var = 0; i_var < 6; i_var++)
                    {
                        variable[i_var] = 0.0;
                        step[i_var] = 0.01;
                    }

                    min.SetFunction(f);

                    // Set the free variables to be minimized!
                    min.SetVariable(0,"delta_x",variable[0], step[0]);
                    min.SetVariable(1,"delta_y",variable[1], step[1]);
                    min.SetVariable(2,"delta_z",variable[2], step[2]);
                    min.SetVariable(3,"alpha",variable[3], step[3]);
                    min.SetVariable(4,"beta",variable[4], step[4]);
                    min.SetVariable(5,"gamma",variable[5], step[5]);

                    // start minimization
                    min.Minimize();

                    const Double_t *min_params = min.X();

                    //for(Int_t i = 0; i < 6; i++)
                    //{
                    //    cout << "par = " << i << ", value = " << min_params[i] << endl;
                    //}

                    delta_x   = min_params[0];
                    delta_y   = min_params[1];
                    delta_z   = min_params[2];
                    rot_alpha = min_params[3];
                    rot_beta  = min_params[4];
                    rot_gamma = min_params[5];


                    rot4.SetToIdentity();
                    rot4.RotateZ(rot_alpha);
                    rot4.RotateX(rot_beta);
                    rot4.RotateZ(rot_gamma);

                    shift_vec.SetXYZ(delta_x,delta_y,delta_z);

                    cout << "" << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "rotation matrix scheme for function Align_HFT_IST" << endl;
                    cout << "| x' |   | xx xy xz | | x |" << endl;
                    cout << "| y' | = | yx yy yz | | y |" << endl;
                    cout << "| z' |   | zx zy zz | | z |" << endl;
                    cout << "|" << rot4.XX() << " " << rot4.XY() << " " << rot4.XZ() << "|" << endl;
                    cout << "|" << rot4.YX() << " " << rot4.YY() << " " << rot4.YZ() << "|" << endl;
                    cout << "|" << rot4.ZX() << " " << rot4.ZY() << " " << rot4.ZZ() << "|" << endl;
                    cout << "shift vector = {" << shift_vec.X() << ", "  << shift_vec.Y() << ", " << shift_vec.Z() << "}" << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "" << endl;


                    if(flag_min_option == 21 && i_align == 0)
                    {
                        cout << "Save IST global alignment parameters to file" << endl;
                        TFile* Outputfile_IST        = new TFile("./Data/Align_params/IST_global_align_V1b.root","RECREATE");

                        HistName = "Rot_IST_glob";
                        rot4.Write(HistName.Data());

                        HistName = "Shift_IST_glob";
                        shift_vec.Write(HistName.Data());
                    }

                    if(flag_min_option == 22 && i_align == 0)
                    {
                        vector_rotation_sectors_IST[i_IST_ladder] = rot4;
                        vector_shift_sectors_IST[i_IST_ladder]    = shift_vec;
                    }

                }
                //------------------------------------------------------------------------------



                //------------------------------------------------------------------------------
                if((flag_min_option == 26 || flag_min_option == 27 || flag_min_option == 28) && i_align == 0)
                {
                    cout << "Set minimizer method" << endl;
                    cout << "Number of tracks: " << IST_align_data_vectors[0].size() <<  endl;
                    // Choose method upon creation between:
                    // kConjugateFR, kConjugatePR, kVectorBFGS,
                    // kVectorBFGS2, kSteepestDescent


                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS ); // <-- good results
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugateFR ); // <-- used
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugatePR );
                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kSteepestDescent ); <-
                    ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS2 );

                    //ROOT::Math::GSLSimAnMinimizer min;

                    //ROOT::Math::GSLMinimizer min( ROOT::Math::kMigradImproved);
                    //ROOT::Math::GSLMultiMinimizer min( ROOT::Math::kVectorBFGS ); does not work, header not found
                    //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);
                    //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
                    //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigradImproved );
                    //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);

                    min.SetMaxFunctionCalls(1000000);
                    min.SetMaxIterations(100000);
                    min.SetTolerance(0.0001);

                    ROOT::Math::Functor f(&Align_HFT_to_TPC,6);

                    Double_t step[6]     = {0.2,0.2,0.2,0.2,0.2,0.2};
                    Double_t variable[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
                    for(Int_t i_var = 0; i_var < 6; i_var++)
                    {
                        variable[i_var] = 0.01; // 0.01
                        step[i_var] = 0.1; // 0.1
                    }

                    min.SetFunction(f);

                    // Set the free variables to be minimized!
                    min.SetVariable(0,"delta_x",variable[0], step[0]);
                    min.SetVariable(1,"delta_y",variable[1], step[1]);
                    min.SetVariable(2,"delta_z",variable[2], step[2]);
                    min.SetVariable(3,"alpha",variable[3], step[3]);
                    min.SetVariable(4,"beta",variable[4], step[4]);
                    min.SetVariable(5,"gamma",variable[5], step[5]);

                    // start minimization
                    min.Minimize();

                    const Double_t *min_params = min.X();
                    const Double_t *err_params = min.Errors();

                    //for(Int_t i = 0; i < 6; i++)
                    //{
                    //    cout << "par = " << i << ", value = " << min_params[i] << endl;
                    //}

                    delta_x   = min_params[0];
                    delta_y   = min_params[1];
                    delta_z   = min_params[2];
                    rot_alpha = min_params[3];
                    rot_beta  = min_params[4];
                    rot_gamma = min_params[5];

                    rot4.SetToIdentity();
                    rot4.RotateZ(rot_alpha);
                    rot4.RotateX(rot_beta);
                    rot4.RotateZ(rot_gamma);

                    shift_vec.SetXYZ(delta_x,delta_y,delta_z);

                    rot_HFT_to_TPC       = rot4;
                    shift_vec_HFT_to_TPC = shift_vec;

                    cout << "" << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "rotation matrix scheme for function Align_HFT_to_TPC" << endl;
                    cout << "| x' |   | xx xy xz | | x |" << endl;
                    cout << "| y' | = | yx yy yz | | y |" << endl;
                    cout << "| z' |   | zx zy zz | | z |" << endl;
                    cout << "|" << rot4.XX() << " " << rot4.XY() << " " << rot4.XZ() << "|" << endl;
                    cout << "|" << rot4.YX() << " " << rot4.YY() << " " << rot4.YZ() << "|" << endl;
                    cout << "|" << rot4.ZX() << " " << rot4.ZY() << " " << rot4.ZZ() << "|" << endl;
                    cout << "shift vector = {" << shift_vec.X() << ", "  << shift_vec.Y() << ", " << shift_vec.Z() << "}" << endl;
                    cout << "------------------------------------------------" << endl;
                    cout << "" << endl;



                    //------------------------------------------------------------------------------
                    if(flag_min_option == 26 || flag_min_option == 27 || flag_min_option == 28)
                    {
                        cout << "Save alignment parameters to output file" << endl;
                        TFile* Outputfile        = new TFile("./Data/Align_params/HFT_to_TPC_global_alignment_Itt2_RFF4.root","RECREATE");

                        HistName = "rot_HFT_to_TPC";
                        rot_HFT_to_TPC.Write(HistName.Data());

                        HistName = "shift_vec_HFT_to_TPC";
                        shift_vec_HFT_to_TPC.Write(HistName.Data());
                    }
                    //------------------------------------------------------------------------------


                }
                //------------------------------------------------------------------------------

            }
        }



        //------------------------------------------------------------------------------
        if(flag_min_option == 22)
        {
            cout << "Save IST ladder-by-ladder alignment parameters to file" << endl;

            TFile* Outputfile_IST_ladder = new TFile("./Data/Align_params/IST_ladder_align_V1.root","RECREATE");

            for(Int_t i_IST_ladder = 0; i_IST_ladder < N_IST_ladder; i_IST_ladder++)
            {
                HistName = "Rot_IST_ladder_";
                HistName += i_IST_ladder;
                vector_rotation_sectors_IST[i_IST_ladder].Write(HistName.Data());

                HistName = "Shift_IST_ladder_";
                HistName += i_IST_ladder;
                vector_shift_sectors_IST[i_IST_ladder].Write(HistName.Data());
            }
        }
        //------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        if(flag_min_option < 30 && flag_min_option >= 25)
        {
            cout << "Plot eta and phi of A and B tracks (two for each cosmic)" << endl;
            TCanvas* c_TPC_track_eta_phi = new TCanvas(HistName.Data(),HistName.Data(),10,10,900,600);
            c_TPC_track_eta_phi->SetFillColor(10);
            c_TPC_track_eta_phi->SetTopMargin(0.05);
            c_TPC_track_eta_phi->SetBottomMargin(0.15);
            c_TPC_track_eta_phi->SetRightMargin(0.05);
            c_TPC_track_eta_phi->SetLeftMargin(0.15);
            c_TPC_track_eta_phi->Divide(2,1);
            for(Int_t i_eta_phi = 0; i_eta_phi < 2; i_eta_phi++)
            {
                c_TPC_track_eta_phi->cd(i_eta_phi+1);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetTopMargin(0.05);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetBottomMargin(0.22);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetRightMargin(0.05);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetLeftMargin(0.2);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetTicks(1,1);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetGrid(0,0);
                c_TPC_track_eta_phi->cd(i_eta_phi+1)->SetFillColor(10);

                for(Int_t i_track = 0; i_track < 2; i_track++)
                {
                    if(i_eta_phi == 0)
                    {
                        h_TPC_track_eta[i_track]->SetStats(0);
                        h_TPC_track_eta[i_track]->SetTitle("");
                        h_TPC_track_eta[i_track]->GetXaxis()->SetTitleOffset(1.1);
                        h_TPC_track_eta[i_track]->GetYaxis()->SetTitleOffset(1.2);
                        h_TPC_track_eta[i_track]->GetYaxis()->SetLabelOffset(0.01);
                        h_TPC_track_eta[i_track]->GetXaxis()->SetLabelSize(0.065);
                        h_TPC_track_eta[i_track]->GetYaxis()->SetLabelSize(0.065);
                        h_TPC_track_eta[i_track]->GetXaxis()->SetTitleSize(0.065);
                        h_TPC_track_eta[i_track]->GetYaxis()->SetTitleSize(0.065);
                        h_TPC_track_eta[i_track]->GetXaxis()->SetNdivisions(505,'N');
                        h_TPC_track_eta[i_track]->GetYaxis()->SetNdivisions(505,'N');
                        h_TPC_track_eta[i_track]->GetXaxis()->CenterTitle();
                        h_TPC_track_eta[i_track]->GetYaxis()->CenterTitle();
                        h_TPC_track_eta[i_track]->GetXaxis()->SetTitle("#eta");
                        h_TPC_track_eta[i_track]->GetYaxis()->SetTitle("counts");
                        h_TPC_track_eta[i_track]->GetXaxis()->SetRangeUser(-1.5,1.5);
                        h_TPC_track_eta[i_track]->SetLineColor(i_track+1);
                        if(i_track == 0) h_TPC_track_eta[i_track]->DrawCopy("h");
                        else h_TPC_track_eta[i_track]->DrawCopy("same h");
                    }
                    if(i_eta_phi == 1)
                    {
                        h_TPC_track_phi[i_track]->SetStats(0);
                        h_TPC_track_phi[i_track]->SetTitle("");
                        h_TPC_track_phi[i_track]->GetXaxis()->SetTitleOffset(1.1);
                        h_TPC_track_phi[i_track]->GetYaxis()->SetTitleOffset(1.2);
                        h_TPC_track_phi[i_track]->GetYaxis()->SetLabelOffset(0.01);
                        h_TPC_track_phi[i_track]->GetXaxis()->SetLabelSize(0.065);
                        h_TPC_track_phi[i_track]->GetYaxis()->SetLabelSize(0.065);
                        h_TPC_track_phi[i_track]->GetXaxis()->SetTitleSize(0.065);
                        h_TPC_track_phi[i_track]->GetYaxis()->SetTitleSize(0.065);
                        h_TPC_track_phi[i_track]->GetXaxis()->SetNdivisions(505,'N');
                        h_TPC_track_phi[i_track]->GetYaxis()->SetNdivisions(505,'N');
                        h_TPC_track_phi[i_track]->GetXaxis()->CenterTitle();
                        h_TPC_track_phi[i_track]->GetYaxis()->CenterTitle();
                        h_TPC_track_phi[i_track]->GetXaxis()->SetTitle("#phi (deg)");
                        h_TPC_track_phi[i_track]->GetYaxis()->SetTitle("counts");
                        h_TPC_track_phi[i_track]->GetXaxis()->SetRangeUser(-180,180);
                        h_TPC_track_phi[i_track]->SetLineColor(i_track+1);
                        if(i_track == 0) h_TPC_track_phi[i_track]->DrawCopy("h");
                        else h_TPC_track_phi[i_track]->DrawCopy("same h");
                    }
                }
            }





            cout << "Plot 1D residual distributions for IST alignment" << endl;
            for(Int_t i_track = 0; i_track < 2; i_track++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                    {
                        c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo);

                        for(Int_t i_align = 0; i_align < 2; i_align++)
                        {
                            if(i_xyz < 2)
                            {
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Rebin(4);
                            }
                            if(i_xyz == 2)
                            {
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Rebin(8);
                            }
                        }


                        for(Int_t i_align = 0; i_align < 2; i_align++)
                        {
                            if(i_align == 0)
                            {
                                if(flag_min_option == 26 || flag_min_option == 27 || flag_min_option == 28)
                                {
                                    h_pixel_line_dca_TPC[1][i_histo][i_xyz][i_track]->DrawCopy("");
                                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->DrawCopy("same");
                                }
                                else
                                {
                                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->DrawCopy("");
                                }
                            }
                            if(i_align == 1 && i_histo >= 0)
                            {
                                //cout << "Plotting aligned distribution" << endl;
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetFillColor(kAzure-1);
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetLineColor(kAzure-1);
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetFillStyle(3001);
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->DrawCopy("same");

                                if(i_xyz < 3)
                                {
                                    //----------------------------------
                                    // first fit
                                    start_fit = -0.4;
                                    stop_fit  = 0.4;
                                    for(Int_t x = 0; x < 3; x++)
                                    {
                                        GaussFit->ReleaseParameter(x);
                                        GaussFit->SetParError(x,0.0);
                                        GaussFit->SetParameter(x,0.0);
                                    }
                                    GaussFit->SetParameter(0,h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetBinContent(h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->FindBin(0)));
                                    GaussFit->SetParameter(1,0.0);
                                    GaussFit->SetParameter(2,0.1);
                                    GaussFit->SetRange(start_fit,stop_fit);
                                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Fit("GaussFit","QN","",start_fit,stop_fit);
                                    height = GaussFit->GetParameter(0);
                                    mean   = GaussFit->GetParameter(1);
                                    sigma  = GaussFit->GetParameter(2);
                                    //----------------------------------



                                    //----------------------------------
                                    // second fit
                                    start_fit = mean-2.0*sigma;
                                    stop_fit  = mean+2.0*sigma;
                                    for(Int_t x = 0; x < 3; x++)
                                    {
                                        GaussFit->ReleaseParameter(x);
                                        GaussFit->SetParError(x,0.0);
                                        GaussFit->SetParameter(x,0.0);
                                    }
                                    GaussFit->SetParameter(0,height);
                                    GaussFit->SetParameter(1,mean);
                                    GaussFit->SetParameter(2,sigma*0.8);
                                    GaussFit->SetRange(start_fit,stop_fit);
                                    h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Fit("GaussFit","QN","",start_fit,stop_fit);
                                    GaussFit->SetLineWidth(1);
                                    GaussFit->SetLineStyle(1);
                                    GaussFit->SetLineColor(2);
                                    GaussFit->SetRange(start_fit,stop_fit);
                                    GaussFit->DrawCopy("same");

                                    height = GaussFit->GetParameter(0);
                                    mean   = GaussFit->GetParameter(1);
                                    sigma  = GaussFit->GetParameter(2);

                                    height_err = GaussFit->GetParError(0);
                                    mean_err   = GaussFit->GetParError(1);
                                    sigma_err  = GaussFit->GetParError(2);

                                    HistName = "#mu = ";
                                    sprintf(NoP,"%2.2f",mean*10.0);
                                    HistName += NoP;
                                    HistName += "#pm";
                                    sprintf(NoP,"%2.2f",mean_err*10.0);
                                    HistName += NoP;
                                    HistName += " mm";
                                    plotTopLegend((char*)HistName.Data(),0.25,0.88-0.13,0.055,1,0.0,42,1,1);

                                    HistName = "#sigma = ";
                                    sprintf(NoP,"%2.2f",sigma*10.0);
                                    HistName += NoP;
                                    HistName += "#pm";
                                    sprintf(NoP,"%2.2f",sigma_err*10.0);
                                    HistName += NoP;
                                    HistName += " mm";
                                    plotTopLegend((char*)HistName.Data(),0.245,0.82-0.13,0.055,1,0.0,42,1,1);
                                    //----------------------------------
                                }

                            }
                        }

                        if(flag_min_option == 26 || flag_min_option == 27 || flag_min_option == 28)
                        {
                            //----------------------------------
                            TLegend* leg_pixel_line_dca = new TLegend(0.23,0.8,0.46,0.93); // x1,y1,x2,y2
                            leg_pixel_line_dca->SetBorderSize(0);
                            leg_pixel_line_dca->SetFillColor(0);
                            leg_pixel_line_dca->SetTextSize(0.055);
                            leg_pixel_line_dca->SetTextFont(42);

                            leg_pixel_line_dca->AddEntry(h_pixel_line_dca_TPC[0][i_histo][i_xyz][i_track],"before align.","fl");
                            leg_pixel_line_dca->AddEntry(h_pixel_line_dca_TPC[1][i_histo][i_xyz][i_track],"after align.","fl");
                            leg_pixel_line_dca->Draw();
                            //----------------------------------
                        }

                        plotTopLegend((char*)label_IST[i_histo].Data(),0.78,0.885,0.055,1,0.0,42,1,1);

                    }
                }


                //----------------------------------
                cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
                HistName = fig_all_output_dir.Data();
                HistName += c_pixel_line_dca_TPC[i_track]->GetName();
                HistName += out_all_format.Data();
                c_pixel_line_dca_TPC[i_track]->SaveAs(HistName.Data(),"");
                cout << "Saved canvas: " << HistName.Data() << endl;
                //----------------------------------
                //------------------------------------------------------------------------------
            }
        }
        //------------------------------------------------------------------------------



        //------------------------------------------------------------------------------
        if(flag_min_option < 29)
        {
            cout << "Plot 1D residual distributions for IST alignment" << endl;

            TCanvas* c_pixel_line_dca_IST; // [inner,outer][x,y,z]
            HistName = "c_pixel_line_dca_IST";
            c_pixel_line_dca_IST = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
            c_pixel_line_dca_IST->SetFillColor(10);
            c_pixel_line_dca_IST->SetTopMargin(0.05);
            c_pixel_line_dca_IST->SetBottomMargin(0.15);
            c_pixel_line_dca_IST->SetRightMargin(0.05);
            c_pixel_line_dca_IST->SetLeftMargin(0.15);
            c_pixel_line_dca_IST->Divide(3,4);

            TString label_IST[4] = {"PXL inner","PXL outer","IST A","IST B"};

            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                for(Int_t i_histo = 0; i_histo < N_histos_IST; i_histo++)
                {
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetTopMargin(0.05);
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetBottomMargin(0.22);
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetRightMargin(0.05);
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetLeftMargin(0.2);
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetTicks(1,1);
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetGrid(0,0);
                    c_pixel_line_dca_IST->cd(i_xyz+1+3*i_histo)->SetFillColor(10);

                    for(Int_t i_align = 0; i_align < 2; i_align++)
                    {
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetStats(0);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetTitle("");
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetTitleOffset(1.2);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetLabelSize(0.065);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetLabelSize(0.065);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetTitleSize(0.065);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetTitleSize(0.065);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->CenterTitle();
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->CenterTitle();
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetTitle(xyz_label[i_xyz]);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetTitle("counts");
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetRangeUser(-0.35,0.35);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetFillColor(kGray+1);
                        h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetFillStyle(3001);

                        if(i_histo < 2) h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetRangeUser(-0.065,0.065);
                        if(i_histo > 1)
                        {
                            h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetRangeUser(-0.35,0.35);
                            if(i_xyz == 2)
                            {
                                h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetRangeUser(-0.52,0.52);
                            }
                        }
                    }


                    for(Int_t i_align = 0; i_align < 2; i_align++)
                    {
                        if(i_align == 0)
                        {
                            if((flag_min_option == 21 || flag_min_option == 22) && i_histo >= 2)
                            {
                                h_pixel_line_dca_IST[1][i_histo][i_xyz]->DrawCopy("");
                                h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->DrawCopy("same");
                            }
                            else
                            {
                                h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->DrawCopy("");
                            }
                        }
                        if(i_align == 1 && i_histo >= 2)
                        {
                            //cout << "Plotting aligned distribution" << endl;
                            h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetFillColor(kAzure-1);
                            h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetLineColor(kAzure-1);
                            h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->SetFillStyle(3001);
                            h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->DrawCopy("same");

                            if(i_xyz < 2)
                            {
                                //----------------------------------
                                // first fit
                                start_fit = -0.06;
                                stop_fit  = 0.06;
                                for(Int_t x = 0; x < 3; x++)
                                {
                                    GaussFit->ReleaseParameter(x);
                                    GaussFit->SetParError(x,0.0);
                                    GaussFit->SetParameter(x,0.0);
                                }
                                GaussFit->SetParameter(0,h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->GetBinContent(h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->FindBin(0)));
                                GaussFit->SetParameter(1,0.0);
                                GaussFit->SetParameter(2,0.005);
                                GaussFit->SetRange(start_fit,stop_fit);
                                h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->Fit("GaussFit","QN","",start_fit,stop_fit);
                                height = GaussFit->GetParameter(0);
                                mean   = GaussFit->GetParameter(1);
                                sigma  = GaussFit->GetParameter(2);
                                //----------------------------------



                                //----------------------------------
                                // second fit
                                start_fit = mean-1.8*sigma;
                                stop_fit  = mean+1.8*sigma;
                                for(Int_t x = 0; x < 3; x++)
                                {
                                    GaussFit->ReleaseParameter(x);
                                    GaussFit->SetParError(x,0.0);
                                    GaussFit->SetParameter(x,0.0);
                                }
                                GaussFit->SetParameter(0,height);
                                GaussFit->SetParameter(1,mean);
                                GaussFit->SetParameter(2,sigma*0.8);
                                GaussFit->SetRange(start_fit,stop_fit);
                                h_pixel_line_dca_IST[i_align][i_histo][i_xyz]->Fit("GaussFit","QN","",start_fit,stop_fit);
                                GaussFit->SetLineWidth(1);
                                GaussFit->SetLineStyle(1);
                                GaussFit->SetLineColor(2);
                                GaussFit->SetRange(start_fit,stop_fit);
                                GaussFit->DrawCopy("same");

                                height = GaussFit->GetParameter(0);
                                mean   = GaussFit->GetParameter(1);
                                sigma  = GaussFit->GetParameter(2);

                                height_err = GaussFit->GetParError(0);
                                mean_err   = GaussFit->GetParError(1);
                                sigma_err  = GaussFit->GetParError(2);

                                HistName = "#mu = ";
                                sprintf(NoP,"%2.1f",mean*10000.0);
                                HistName += NoP;
                                HistName += "#pm";
                                sprintf(NoP,"%2.1f",mean_err*10000.0);
                                HistName += NoP;
                                HistName += " #mum";
                                plotTopLegend((char*)HistName.Data(),0.25,0.88-0.13,0.055,1,0.0,42,1,1);

                                HistName = "#sigma = ";
                                sprintf(NoP,"%2.1f",sigma*10000.0);
                                HistName += NoP;
                                HistName += "#pm";
                                sprintf(NoP,"%2.1f",sigma_err*10000.0);
                                HistName += NoP;
                                HistName += " #mum";
                                plotTopLegend((char*)HistName.Data(),0.245,0.82-0.13,0.055,1,0.0,42,1,1);
                                //----------------------------------
                            }

                        }
                    }

                    if((flag_min_option == 21 || flag_min_option == 22) && i_histo >= 2)
                    {
                        //----------------------------------
                        TLegend* leg_pixel_line_dca = new TLegend(0.23,0.8,0.46,0.93); // x1,y1,x2,y2
                        leg_pixel_line_dca->SetBorderSize(0);
                        leg_pixel_line_dca->SetFillColor(0);
                        leg_pixel_line_dca->SetTextSize(0.055);
                        leg_pixel_line_dca->SetTextFont(42);

                        leg_pixel_line_dca->AddEntry(h_pixel_line_dca_IST[0][i_histo][i_xyz],"before align.","fl");
                        leg_pixel_line_dca->AddEntry(h_pixel_line_dca_IST[1][i_histo][i_xyz],"after align.","fl");
                        leg_pixel_line_dca->Draw();
                        //----------------------------------
                    }

                    plotTopLegend((char*)label_IST[i_histo].Data(),0.78,0.885,0.055,1,0.0,42,1,1);

                }
            }


            //----------------------------------
            cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
            HistName = fig_all_output_dir.Data();
            HistName += c_pixel_line_dca_IST->GetName();
            HistName += out_all_format.Data();
            c_pixel_line_dca_IST->SaveAs(HistName.Data(),"");
            cout << "Saved canvas: " << HistName.Data() << endl;
            //----------------------------------
            //------------------------------------------------------------------------------



            //------------------------------------------------------------------------------
            cout << "Plot 2D residual distributions for IST alignment" << endl;

            TCanvas* c2D_pixel_line_dca_IST[2]; // [ialign] [inner,outer][x,y,z]
            for(Int_t i_align = 0; i_align < 2; i_align++)
            {
                HistName = "c2D_pixel_line_dca_IST_";
                HistName += i_align;
                c2D_pixel_line_dca_IST[i_align] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
                c2D_pixel_line_dca_IST[i_align]->SetFillColor(10);
                c2D_pixel_line_dca_IST[i_align]->SetTopMargin(0.05);
                c2D_pixel_line_dca_IST[i_align]->SetBottomMargin(0.15);
                c2D_pixel_line_dca_IST[i_align]->SetRightMargin(0.05);
                c2D_pixel_line_dca_IST[i_align]->SetLeftMargin(0.15);
                c2D_pixel_line_dca_IST[i_align]->Divide(3,4);

                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    for(Int_t i_histo = 0; i_histo < N_histos_IST; i_histo++)
                    {
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetTopMargin(0.05);
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetBottomMargin(0.22);
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetRightMargin(0.05);
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetLeftMargin(0.2);
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetTicks(1,1);
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetGrid(0,0);
                        c2D_pixel_line_dca_IST[i_align]->cd(i_xyz+1+3*i_histo)->SetFillColor(10);

                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->SetStats(0);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->SetTitle("");
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetTitleOffset(1.2);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetLabelSize(0.065);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetLabelSize(0.065);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetTitleSize(0.065);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetTitleSize(0.065);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->CenterTitle();
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->CenterTitle();
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetTitle(xyz_label[i_xyz]);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetTitle("z_{global} (cm)");
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetRangeUser(-0.35,0.35);
                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetXaxis()->SetRangeUser(-22,22);

                        if(i_histo < 2) h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetRangeUser(-0.065,0.065);
                        if(i_histo > 1)
                        {
                            h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetRangeUser(-0.35,0.35);
                            if(i_xyz == 2)
                            {
                                h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->GetYaxis()->SetRangeUser(-0.52,0.52);
                            }
                        }

                        h2D_pixel_line_dca_vs_z_IST[i_align][i_histo][i_xyz]->DrawCopy("colz");

                        PlotLine(-22,22,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle


                        plotTopLegend((char*)label_IST[i_histo].Data(),0.78,0.885,0.055,1,0.0,42,1,1);

                    }
                }


                //----------------------------------
                cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
                HistName = fig_all_output_dir.Data();
                HistName += c2D_pixel_line_dca_IST[i_align]->GetName();
                HistName += out_all_format.Data();
                c2D_pixel_line_dca_IST[i_align]->SaveAs(HistName.Data(),"");
                cout << "Saved canvas: " << HistName.Data() << endl;
                //----------------------------------
            }
            //------------------------------------------------------------------------------



            //------------------------------------------------------------------------------
            cout << "Plot 2D residual distributions for IST alignment for every ladder" << endl;

            TCanvas* c2D_pixel_line_dca_IST_ladder[2][3]; // [before,after][x,y,z]
            for(Int_t i_align = 0; i_align < 2; i_align++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    HistName = "c2D_pixel_line_dca_IST_ladder_";
                    HistName += i_align;
                    HistName += "_";
                    HistName += i_xyz;
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1400,900);
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->SetFillColor(10);
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->SetTopMargin(0.02);
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->SetBottomMargin(0.19);
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->SetRightMargin(0.05);
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->SetLeftMargin(0.15);
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->Divide(4,6); // 24 ladders

                    for(Int_t i_histo = 0; i_histo < N_IST_ladder; i_histo++)
                    {
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetTopMargin(0.03);
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetBottomMargin(0.2);
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetRightMargin(0.05);
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetLeftMargin(0.15);
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetTicks(1,1);
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetGrid(0,0);
                        c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->cd(i_histo+1)->SetFillColor(10);

                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->SetStats(0);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->SetTitle("");
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->SetTitleOffset(1.1);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetTitleOffset(0.9);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetLabelOffset(0.01);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->SetLabelSize(0.085);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetLabelSize(0.085);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->SetTitleSize(0.085);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetTitleSize(0.085);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->SetNdivisions(505,'N');
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetNdivisions(505,'N');
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->CenterTitle();
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->CenterTitle();
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetTitle(xyz_label[i_xyz]);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->SetTitle("z_{global} (cm)");
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetYaxis()->SetRangeUser(-0.35,0.35);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->GetXaxis()->SetRangeUser(-22,22);
                        h2D_pixel_line_dca_vs_z_IST_ladder[i_align][i_histo][i_xyz]->DrawCopy("colz");

                        PlotLine(-22,22,0,0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle

                        HistName = "IST, ladder: ";
                        HistName += i_histo+1;
                        plotTopLegend((char*)HistName.Data(),0.73,0.85,0.085,1,0.0,42,1,1);
                    }

                    //----------------------------------
                    cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
                    HistName = fig_all_output_dir.Data();
                    HistName += c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->GetName();
                    HistName += out_all_format.Data();
                    c2D_pixel_line_dca_IST_ladder[i_align][i_xyz]->SaveAs(HistName.Data(),"");
                    cout << "Saved canvas: " << HistName.Data() << endl;
                    //----------------------------------
                }
            }
        }
        //------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------------------------



    //---------------------------------------------------------------------------------------------------------------
    // TPC sector alignment using PXL hits
    if(flag_min_option == 30 || flag_min_option == 27 || flag_min_option == 31 || flag_min_option == 32 || flag_min_option == 33)
    {
        // flag_min_option
        // 30 -> do TPC sector alignment
        // 27 -> do first HFT to TPC alignment
        // 31 -> do TPC sector alignment with pre HFT to TPC alignment applied from file
        // 32 -> do TPC sector alignment with pre HFT to TPC alignment and pre TPC sector alignment applied from file
        // 33 -> do TPC sector alignment with pre HFT to TPC alignment and pre TPC sector alignment applied from file and use pre alignment from BNL

        //-----------------------------------------------------
        if(flag_min_option == 31)
        {
            cout << "Open HFT to TPC global alignment parameters" << endl;
            TFile* Inputfile_HFT_to_TPC_glob  = TFile::Open("./Data/Align_params/HFT_to_TPC_global_alignment.root");  // open the file

            HistName = "rot_HFT_to_TPC";
            rot_HFT_to_TPC = *(TRotation*)(Inputfile_HFT_to_TPC_glob->Get(HistName.Data()));
            HistName = "shift_vec_HFT_to_TPC";
            shift_vec_HFT_to_TPC = *(TVector3*)(Inputfile_HFT_to_TPC_glob->Get(HistName.Data()));
        }
        //-----------------------------------------------------



        //-----------------------------------------------------
        if(flag_min_option == 32 || flag_min_option == 33) // open pre alignment
        {
            // This one acts on the HFT hits
            cout << "Open HFT to TPC global alignment parameters" << endl;
            TFile* Inputfile_HFT_to_TPC_glob  = TFile::Open("./Data/Align_params/HFT_to_TPC_global_alignment.root");  // open the file

            HistName = "rot_HFT_to_TPC";
            rot_HFT_to_TPC_load = *(TRotation*)(Inputfile_HFT_to_TPC_glob->Get(HistName.Data()));
            HistName = "shift_vec_HFT_to_TPC";
            shift_vec_HFT_to_TPC_load = *(TVector3*)(Inputfile_HFT_to_TPC_glob->Get(HistName.Data()));


            // This one acts on the HFT hits
            cout << "Open second itteration HFT to TPC global alignment parameters" << endl;
            TFile* Inputfile_HFT_to_TPC_globB  = TFile::Open("./Data/Align_params/HFT_to_TPC_global_alignment_Itt2.root");  // open the file

            HistName = "rot_HFT_to_TPC";
            rot_HFT_to_TPC_loadB = *(TRotation*)(Inputfile_HFT_to_TPC_globB->Get(HistName.Data()));
            HistName = "shift_vec_HFT_to_TPC";
            shift_vec_HFT_to_TPC_loadB = *(TVector3*)(Inputfile_HFT_to_TPC_globB->Get(HistName.Data()));

            // This one acts on the TPC tracks
            for(Int_t i_sec = 0; i_sec < N_TPC_sectors; i_sec++)
            {
                HistName = "./Data/Align_params/";
                HistName += "TPC_sec";
                HistName += i_sec+1;
                HistName += "_to_HFT_align.root";
                cout << "Open TPC sector alignment file: " << HistName.Data() << endl;
                Inputfile_TPC_sector_alignment[i_sec] = TFile::Open(HistName.Data());  // open the file

                HistName = "Rot_matrix";
                rot_TPC_sector_load[i_sec] = *(TRotation*)(Inputfile_TPC_sector_alignment[i_sec]->Get(HistName.Data()));

                HistName = "Shift_vec";
                shift_vec_TPC_sector_load[i_sec] = *(TVector3*)(Inputfile_TPC_sector_alignment[i_sec]->Get(HistName.Data()));
            }

            if(flag_min_option == 33)
            {
                // This one acts on the TPC tracks --> Yuri's old to new alignment
                HistName = "./Data/Align_params/Hybrid_Yuri_old_to_new.root";
                cout << "Open TPC sector alignment file (Yuri old to new): " << HistName.Data() << endl;
                Inputfile_TPC_sector_alignment_Yuri_old_new = TFile::Open(HistName.Data());  // open the file

                for(Int_t i_sec = 0; i_sec < N_TPC_sectors; i_sec++)
                {
                    HistName = "Hybrid_Rot_S";
                    HistName += i_sec+1;
                    cout << "Open rotation matrix: " << HistName.Data() << endl;
                    rot_TPC_sector_load_Yuri_old_new[i_sec] = *(TRotation*)(Inputfile_TPC_sector_alignment_Yuri_old_new->Get(HistName.Data()));

                    HistName = "Hybrid_Shift_S";
                    HistName += i_sec+1;
                    cout << "Open shift vector: " << HistName.Data() << endl;
                    shift_vec_TPC_sector_load_Yuri_old_new[i_sec] = *(TVector3*)(Inputfile_TPC_sector_alignment_Yuri_old_new->Get(HistName.Data()));
                }
            }
        }
        //-----------------------------------------------------



        Int_t TPC_sec_align = pxl_sector_L;
        TPC_sector_align_active = TPC_sec_align;

        if(stop_event_use > file_entries_PXL_IST) stop_event_use = file_entries_PXL_IST;
        cout << "Doing TPC sector alignment using PXL hits with " << stop_event_use << " entries" << endl;


        TVector3 PXL_hit_IO_AB[5]; // A/B
        TVector3 PXL_hit_IA, PXL_hit_OA, dca_vec_IA, dca_vec_IB, dca_vec_OA, dca_vec_OB, PXL_dir_A, PXL_dir_B, dca_vec_IST_A, dca_vec_IST_B;
        IST_align_data_vectors.resize(7); // hit_IA, hit_OA, hit_IB, hit_OB, dir_A, dir_B, hit_IST
        TPC_align_data_vectors.resize(2); // trackA, trackB
        vec_magfac.resize(2); // trackA, trackB
        vec_charge.resize(2); // trackA, trackB

        //-----------------------------------------------------
        cout << "Define histograms" << endl;
        const Int_t N_IST_ladder = 24;
        const Int_t N_histos_IST = 4;
        TH1F* h_pixel_line_dca_IST[2][N_histos_IST][3];
        TH2F* h2D_pixel_line_dca_vs_z_IST[2][N_histos_IST][3];
        TH2F* h2D_pixel_line_dca_vs_z_IST_ladder[2][N_IST_ladder][3];
        TH1F* h_TPC_track_eta[2];
        TH1F* h_TPC_track_phi[2];
        TH1D* h_Good_TPC_tracks_sector = new TH1D("h_Good_TPC_tracks_sector","h_Good_TPC_tracks_sector",24,1,24);

        const Int_t N_HFT_TPC_match = 5;

        const Int_t N_TPC_sectors = 24;
        Int_t TPC_good_tracks_sectors[N_TPC_sectors];
        for(Int_t iTPC_sec = 0; iTPC_sec < N_TPC_sectors; iTPC_sec++)
        {
            TPC_good_tracks_sectors[iTPC_sec] = 0;
        }

        // PXL alignment parameters
        const Int_t N_full_sec_align = 10;
        std::vector<TRotation> vector_rotation_sectors_STS;
        std::vector<TVector3>  vector_shift_sectors_STS;
        vector_rotation_sectors_STS.resize(N_full_sec_align);
        vector_shift_sectors_STS.resize(N_full_sec_align);

        TRotation IST_glob_rot;
        TVector3  IST_glob_shift;
        IST_glob_rot.SetToIdentity();
        IST_glob_shift.SetXYZ(0.0,0.0,0.0);

        std::vector<TRotation> vector_rotation_sectors_IST;
        std::vector<TVector3>  vector_shift_sectors_IST;
        vector_rotation_sectors_IST.resize(N_IST_ladder);
        vector_shift_sectors_IST.resize(N_IST_ladder);

        for(Int_t i_IST_ladder = 0; i_IST_ladder < N_IST_ladder; i_IST_ladder++)
        {
            vector_rotation_sectors_IST[i_IST_ladder].SetToIdentity();
            vector_shift_sectors_IST[i_IST_ladder].SetXYZ(0.0,0.0,0.0);
        }


        Int_t N_align = 2;

        delta_x = 0.0;
        delta_y = 0.0;
        delta_z = 0.0;

        for(Int_t i_align = 0; i_align < N_align; i_align++) // first loop -> alignment, second loop -> apply alignment
        {
            cout << "i_align: " << i_align << endl;

            StThreeVectorF shift_vec_TPC_track;
            shift_vec_TPC_track.set(delta_x,delta_y,delta_z);

            Int_t track_counter = 0;

            for(Int_t i_data = 0; i_data < 7; i_data++)
            {
                IST_align_data_vectors[i_data].clear();
            }
            for(Int_t i_tpc_track = 0; i_tpc_track < 2; i_tpc_track++)
            {
                TPC_align_data_vectors[i_tpc_track].clear();
                vec_magfac[i_tpc_track].clear();
                vec_charge[i_tpc_track].clear();
            }

            Double_t track_pxA_old = 0.0;
            Double_t track_pyA_old = 0.0;
            Double_t track_pzA_old = 0.0;
            Double_t track_pxB_old = 0.0;
            Double_t track_pyB_old = 0.0;
            Double_t track_pzB_old = 0.0;

            for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
            {
                if (counter != 0  &&  counter % 1000 == 0)
                    cout << "." << flush;
                if (counter != 0  &&  counter % 10000 == 0)
                {
                    if((stop_event_use-start_event_use) > 0)
                    {
                        Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
                        cout << " " << counter << " (" << event_percent << "%) " << "\n" << "==> Processing data (IST or TPC alignment) " << flush;
                    }
                }

#if 0
                cout << "" << endl;
                cout << "------------------------------------------------------" << endl;
                cout << "counter = " << counter << ", out of " << stop_event_use << endl;
#endif

                if (!NT_PXL_IST->GetEntry( counter )) // take the event -> information is stored in event
                    break;  // end of data chunk

                // PXL sector information
                Int_t sector_pxl_array[5]; // [IA,OA,IB,OB,IST]
                sector_pxl_array[0] = get_HFT_pixel_sector(PXL_id_IA,IST_ladder_number); // inner PXL hit A
                sector_pxl_array[1] = get_HFT_pixel_sector(PXL_id_OA,IST_ladder_number); // outer PXL hit A
                sector_pxl_array[2] = get_HFT_pixel_sector(PXL_id_IB,IST_ladder_number); // inner PXL hit B
                sector_pxl_array[3] = get_HFT_pixel_sector(PXL_id_OB,IST_ladder_number); // outer PXL hit B
                sector_pxl_array[4] = -1; // IST

                // PXL index information
                Int_t index_pxl_array[5]; // [IA,OA,IB,OB,IST]
                index_pxl_array[0] = get_HFT_det_index(PXL_id_IA);
                index_pxl_array[1] = get_HFT_det_index(PXL_id_OA);
                index_pxl_array[2] = get_HFT_det_index(PXL_id_IB);
                index_pxl_array[3] = get_HFT_det_index(PXL_id_OB);
                index_pxl_array[4] = -1;

                // Set PXL vectors for A/B and inner/outer
                PXL_hit_IO_AB[0].SetXYZ(PXL_IA_x,PXL_IA_y,PXL_IA_z);
                PXL_hit_IO_AB[1].SetXYZ(PXL_OA_x,PXL_OA_y,PXL_OA_z);
                PXL_hit_IO_AB[2].SetXYZ(PXL_IB_x,PXL_IB_y,PXL_IB_z);
                PXL_hit_IO_AB[3].SetXYZ(PXL_OB_x,PXL_OB_y,PXL_OB_z);
                PXL_hit_IO_AB[4].SetXYZ(IST_x,IST_y,IST_z);

                if(flag_min_option == 27 || flag_min_option == 31) // first HFT to TPC alignment
                {
                    for(Int_t HFT_hit = 0; HFT_hit < 5; HFT_hit++)
                    {
                        PXL_hit_IO_AB[HFT_hit].Transform(rot_HFT_to_TPC);
                        PXL_hit_IO_AB[HFT_hit] += shift_vec_HFT_to_TPC;
                    }
                }

                if(flag_min_option == 32 || flag_min_option == 33) // first two itterations of HFT to TPC alignment
                {
                    for(Int_t HFT_hit = 0; HFT_hit < 5; HFT_hit++)
                    {
                        PXL_hit_IO_AB[HFT_hit].Transform(rot_HFT_to_TPC_load);
                        PXL_hit_IO_AB[HFT_hit] += shift_vec_HFT_to_TPC_load;

                        PXL_hit_IO_AB[HFT_hit].Transform(rot_HFT_to_TPC_loadB);
                        PXL_hit_IO_AB[HFT_hit] += shift_vec_HFT_to_TPC_loadB;
                    }
                }

                if(idx_IA == idx_IB)
                {
                    PXL_hit_IO_AB[2].SetXYZ(-999.0,-999.0,-999.0);
                }
                if(idx_OA == idx_OB)
                {
                    PXL_hit_IO_AB[3].SetXYZ(-999.0,-999.0,-999.0);
                }


                // TPC track information, one cosmic creates two tracks (A,B)
                StThreeVectorF vec_gMomA, vec_OriginA, vec_gMomB, vec_OriginB;
                StPhysicalHelixD helixA, helixB;
                Double_t phiA,phiB,pA,pB,etaA,etaB;

#if 0
                cout << "track A p: {" <<  track_pxA << ", " << track_pyA << ", " << track_pzA << "}" << endl;
#endif

                vec_gMomA.set(track_pxA,track_pyA,track_pzA);
                vec_OriginA.set(track_oxA,track_oyA,track_ozA);
                helixA = StPhysicalHelixD(vec_gMomA,vec_OriginA,magfac,track_qA);

                vec_gMomB.set(track_pxB,track_pyB,track_pzB);
                vec_OriginB.set(track_oxB,track_oyB,track_ozB);
                helixB = StPhysicalHelixD(vec_gMomB,vec_OriginB,magfac,track_qB);

                StThreeVectorF vecA = helixA.at(0);
                StThreeVectorF vecB = helixB.at(0);

                pA   = helixA.momentum((Double_t)magfac).mag();
                pB   = helixB.momentum((Double_t)magfac).mag();

                StThreeVectorF vectoratsA  = helixA.cat(0);
                etaA = vectoratsA.pseudoRapidity();

                StThreeVectorF vectoratsB  = helixB.cat(0);
                etaB = vectoratsB.pseudoRapidity();

                //cout << "counter = " << counter << ", phiA = " << phiA << ", phiB = " << phiB << ", etaA = " << etaA << ", etaB = " << etaB << endl;
                //cout << "vecA = {" << vecA.x() << ", " << vecA.y() << ", " << vecA.z() << "}" << endl;
                //cout << "vecB = {" << vecB.x() << ", " << vecB.y() << ", " << vecB.z() << "}" << endl;

                if(idx_IA == idx_IB && idx_OA == idx_OB) continue;

                // calculate direction vectors for align and refrence sector(s)
                PXL_dir_A = PXL_hit_IO_AB[0] - PXL_hit_IO_AB[1]; // inner - outer for A
                PXL_dir_B = PXL_hit_IO_AB[2] - PXL_hit_IO_AB[3]; // inner - outer for B

                // Don't use the same pixel hit points
                dca_vec_IA = 0.0;
                if(idx_IA != idx_IB) dca_vec_IA = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[0],PXL_dir_A,PXL_hit_IO_AB[2]); // base,dir,point

                dca_vec_IB = 0.0;
                if(idx_OA != idx_OB) dca_vec_IA = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[0],PXL_dir_A,PXL_hit_IO_AB[3]); // base,dir,point

                // calculate dca values of A/B tracklets to IST hit
                dca_vec_IST_A = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[0],PXL_dir_A,PXL_hit_IO_AB[4]); // base,dir,point
                dca_vec_IST_B = calculateDCA_vec_StraightToPoint(PXL_hit_IO_AB[2],PXL_dir_B,PXL_hit_IO_AB[4]); // base,dir,point

                if(
                   dca_vec_IA.Mag()       < 0.4
                   && dca_vec_OA.Mag()    < 0.4
                   //&& dca_vec_IST_A.Mag() < 0.4
                   //&& dca_vec_IST_B.Mag() < 0.4
                   && IST_counter > -1
                   && ((secA == TPC_sec_align && pA > min_TPC_track_mom) || (secB == TPC_sec_align && pB > min_TPC_track_mom))
                   && ((secA == TPC_sec_align && (NsecA + NnohitsecA) == 44) || (secB == TPC_sec_align &&  (NsecB + NnohitsecB) == 44)) // total number of hits + no hit from the same sector is 44 (means 45 hits in total)
                   && ((secA == TPC_sec_align && NsecA > 30) || (secB == TPC_sec_align && NsecB > 30))
                   && secA >= 1 && secA <= 24
                   && secB >= 1 && secB <= 24
                  )
                {
                    Int_t flag_good_tpc_track = 1;
                    Double_t dca_array[N_HFT_TPC_match][3][2]; // [HFT hit][x,y,z][first track, second track]
                    for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                    {
                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            for(Int_t i_track = 0; i_track < 2; i_track++)
                            {
                                dca_array[i_histo][i_xyz][i_track] = 0.0;
                            }
                        }
                    }


                    StPhysicalHelixD helixAB[2];
                    helixAB[0] = helixA;
                    helixAB[1] = helixB;
                    Double_t vec_magfac_helix[2] = {magfac,magfac};
                    Double_t charge_helix[2]     = {track_qA,track_qB};



                    //---------------------------------------------------
                    if(flag_min_option == 32 || flag_min_option == 33)
                    {
                        Double_t mag_helix[2]       = {magfac,magfac};
                        Double_t q_helix[2]         = {track_qA,track_qB};
                        Int_t TPC_sector_helix[2]   = {secA-1,secB-1};
                        // Apply TPC sector pre alignment from file
                        for(Int_t i_helix = 0; i_helix < 2; i_helix++)
                        {
                            StThreeVectorD vec_TPC_Origin = helixAB[i_helix].origin();
                            StThreeVectorD vec_TPC_Mom    = helixAB[i_helix].momentum(mag_helix[i_helix]);
                            Int_t          TPC_charge     = helixAB[i_helix].charge(q_helix[i_helix]);

                            TVector3 TV3_vec_TPC_Mom, TV3_vec_TPC_Origin;
                            TV3_vec_TPC_Mom.SetXYZ(vec_TPC_Mom.x(),vec_TPC_Mom.y(),vec_TPC_Mom.z()); // set TVector3 from original momentum vector
                            TV3_vec_TPC_Mom.Transform(rot_TPC_sector_load[TPC_sector_helix[i_helix]]); // rotate direction momentum vector

                            //rot_TPC_sector_load_Yuri_old_new[i_sec]
                            //shift_vec_TPC_sector_load_Yuri_old_new[i_sec]

                            if(flag_min_option == 33)
                            {
#if 0
                                cout << "i_helix: " << i_helix << ", sector: " << TPC_sector_helix[i_helix] << endl;
                                cout << "XX: " << rot_TPC_sector_load_Yuri_old_new[TPC_sector_helix[i_helix]].XX() << endl;
                                cout << "XY: " << rot_TPC_sector_load_Yuri_old_new[TPC_sector_helix[i_helix]].XY() << endl;
                                cout << "XZ: " << rot_TPC_sector_load_Yuri_old_new[TPC_sector_helix[i_helix]].XZ() << endl;
#endif
                                TV3_vec_TPC_Mom.Transform(rot_TPC_sector_load_Yuri_old_new[TPC_sector_helix[i_helix]]); // rotate direction momentum vector
                            }
                            vec_TPC_Mom.set(TV3_vec_TPC_Mom.x(),TV3_vec_TPC_Mom.y(),TV3_vec_TPC_Mom.z()); // set rotated momentum vector in StThreeVectorD format
                            TV3_vec_TPC_Origin.SetXYZ(vec_TPC_Origin.x(),vec_TPC_Origin.y(),vec_TPC_Origin.z()); // set TVector3 from oginal Origina vector
                            TV3_vec_TPC_Origin += shift_vec_TPC_sector_load[TPC_sector_helix[i_helix]]; // shift Origin vector
                            if(flag_min_option == 33) TV3_vec_TPC_Origin += shift_vec_TPC_sector_load_Yuri_old_new[TPC_sector_helix[i_helix]]; // shift Origin vector
                            vec_TPC_Origin.set(TV3_vec_TPC_Origin.x(),TV3_vec_TPC_Origin.y(),TV3_vec_TPC_Origin.z()); // set shifted

                            helixAB[i_helix] = StPhysicalHelixD(vec_TPC_Mom,vec_TPC_Origin,mag_helix[i_helix],TPC_charge);
                        }
                    }
                    //---------------------------------------------------



                    //---------------------------------------------------
                    if(i_align == 1)
                    {
                        // Apply TPC sector alignment
                        for(Int_t i_helix = 0; i_helix < 2; i_helix++)
                        {
                            StThreeVectorF vec_TPC_Origin = helixAB[i_helix].origin();
                            StThreeVectorF vec_TPC_Mom    = helixAB[i_helix].momentum(vec_magfac_helix[i_helix]);
                            Int_t          TPC_charge     = helixAB[i_helix].charge(charge_helix[i_helix]);

                            vec_TPC_Mom.rotateX(rot_alpha);
                            vec_TPC_Mom.rotateY(rot_beta);
                            vec_TPC_Mom.rotateZ(rot_gamma);

                            vec_TPC_Origin += shift_vec_TPC_track;

                            helixAB[i_helix] = StPhysicalHelixD(vec_TPC_Mom,vec_TPC_Origin,vec_magfac_helix[i_helix],TPC_charge);
                        }
                    }
                    //---------------------------------------------------



                    for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                    {
                        StThreeVectorF pixel_hit, helix_pointA, helix_pointB, dca_vecA, dca_vecB;

                        if(PXL_hit_IO_AB[i_histo].X() < -998.0) continue;

                        pixel_hit.set(PXL_hit_IO_AB[i_histo].X(),PXL_hit_IO_AB[i_histo].Y(),PXL_hit_IO_AB[i_histo].Z());


                        Float_t pathA,dcaA;
                        fHelixAtoPointdca(pixel_hit,helixAB[0],pathA,dcaA);
                        helix_pointA = helixAB[0].at(pathA);

                        Float_t pathB,dcaB;
                        fHelixAtoPointdca(pixel_hit,helixAB[1],pathB,dcaB);
                        helix_pointB = helixAB[1].at(pathB);

                        dca_vecA = helix_pointA;
                        dca_vecA -= pixel_hit;

                        dca_vecB = helix_pointB;
                        dca_vecB -= pixel_hit;

                        //cout << "tree counter = " << counter << ", HFT point = " << i_histo << ", dcaA = " << dcaA
                        //    << ", hit = {" << PXL_hit_IO_AB[i_histo].X() << ", " << PXL_hit_IO_AB[i_histo].Y() << ", " << PXL_hit_IO_AB[i_histo].z() << "}" << endl;
                        if(dcaA > 0.5 || dcaB > 0.5) flag_good_tpc_track = 0;
                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            dca_array[i_histo][i_xyz][0] = dca_vecA[i_xyz];
                            dca_array[i_histo][i_xyz][1] = dca_vecB[i_xyz];

                            if(flag_good_tpc_track)
                            {
                                for(Int_t i_track = 0; i_track < 2; i_track++)
                                {
                                    if((i_track == 0 && secA == TPC_sec_align) || (i_track == 1 && secB == TPC_sec_align))
                                    {
                                        h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Fill(dca_array[i_histo][i_xyz][i_track]);
                                    }
                                }
                            }
                        }
                    }

#if 0
                    //if(!(track_pxA == track_pxA_old &&  track_pyA == track_pyA_old && track_pzA == track_pzA_old &&
                    //     track_pxB == track_pxB_old &&  track_pyB == track_pyB_old && track_pzB == track_pzB_old))
                    if(secA == 6 || secB == 6)
                    {
                        cout << "PXL IA: {" << PXL_IA_x << ", " << PXL_IA_y  << ", " << PXL_IA_z << "}, dca A: {" << dca_array[0][0][0] << ", " << dca_array[0][1][0] << ", " << dca_array[0][2][0] << "}, dca B: " << dca_array[0][0][1] << ", " << dca_array[0][1][1] << ", " << dca_array[0][2][1] << endl;
                        cout << "PXL OA: {" << PXL_OA_x << ", " << PXL_OA_y  << ", " << PXL_OA_z << "}, dca A: {" << dca_array[1][0][0] << ", " << dca_array[1][1][0] << ", " << dca_array[1][2][0] << "}, dca B: " << dca_array[1][0][1] << ", " << dca_array[1][1][1] << ", " << dca_array[1][2][1] << endl;
                        cout << "PXL IB: {" << PXL_IB_x << ", " << PXL_IB_y  << ", " << PXL_IB_z << "}, dca A: {" << dca_array[2][0][0] << ", " << dca_array[2][1][0] << ", " << dca_array[2][2][0] << "}, dca B: " << dca_array[2][0][1] << ", " << dca_array[2][1][1] << ", " << dca_array[2][2][1] << endl;
                        cout << "PXL OB: {" << PXL_OB_x << ", " << PXL_OB_y  << ", " << PXL_OB_z << "}, dca A: {" << dca_array[3][0][0] << ", " << dca_array[3][1][0] << ", " << dca_array[3][2][0] << "}, dca B: " << dca_array[3][0][1] << ", " << dca_array[3][1][1] << ", " << dca_array[3][2][1] << endl;
                        cout << "IST : {" << IST_x << ", " << IST_y  << ", " << IST_z << "}, IST_counter: " << IST_counter << ", dca A: {" << dca_array[4][0][0] << ", " << dca_array[4][1][0] << ", " << dca_array[4][2][0] << "}, dca B: " << dca_array[4][0][1] << ", " << dca_array[4][1][1] << ", " << dca_array[4][2][1] << endl;
                        cout << "flag_good_tpc_track: " << flag_good_tpc_track << ", secA: " << secA << endl;
                    }
#endif

                    if(flag_good_tpc_track
                       //&& (secA == 6 || secB == 6)
                       //&& (secB == 6)
                      )
                    {
                        if(!(track_pxA == track_pxA_old &&  track_pyA == track_pyA_old && track_pzA == track_pzA_old))
                        {
                            TPC_good_tracks_sectors[((Int_t)secA)-1]++;
                        }
                        if(!(track_pxB == track_pxB_old &&  track_pyB == track_pyB_old && track_pzB == track_pzB_old))
                        {
                            TPC_good_tracks_sectors[((Int_t)secB)-1]++;
                        }

                        if(secA == TPC_sec_align || secB == TPC_sec_align)
                        {
                            IST_align_data_vectors[0].push_back(PXL_hit_IO_AB[0]);
                            IST_align_data_vectors[1].push_back(PXL_hit_IO_AB[1]);
                            IST_align_data_vectors[2].push_back(PXL_hit_IO_AB[2]);
                            IST_align_data_vectors[3].push_back(PXL_hit_IO_AB[3]);
                            IST_align_data_vectors[4].push_back(PXL_dir_A);
                            IST_align_data_vectors[5].push_back(PXL_dir_B);
                            IST_align_data_vectors[6].push_back(PXL_hit_IO_AB[4]);

                            if(secA == TPC_sec_align)
                            {
                                TPC_align_data_vectors[0].push_back(helixAB[0]);
                            }
                            else
                            {
                                TPC_align_data_vectors[0].push_back(helixAB[1]);
                            }
                            TPC_align_data_vectors[1].push_back(helixAB[1]);

                            vec_magfac[0].push_back(magfac);
                            vec_magfac[1].push_back(magfac);
                            vec_charge[0].push_back(track_qA);
                            vec_charge[1].push_back(track_qB);

                            if(i_align == 0)
                            {
                                cout << "accepted track_counter = " << track_counter << ", tree counter = " << counter << ", pA = " << pA << ", pB = " << pB <<
                                    ", sec: " << TPC_good_tracks_sectors[((Int_t)TPC_sec_align)-1] << endl;
                            }
                        }

                        for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                        {
                            Double_t dca_TPC[2];
                            for(Int_t i_track = 0; i_track < 2; i_track++)
                            {
                                dca_TPC[i_track] = TMath::Sqrt(dca_array[i_histo][0][i_track]*dca_array[i_histo][0][i_track] + dca_array[i_histo][1][i_track]*dca_array[i_histo][1][i_track] + dca_array[i_histo][2][i_track]*dca_array[i_histo][2][i_track]);
                                //cout << "i_track  = " << i_track << ", PXL = " << i_histo << ", dca_TPC = " << dca_TPC[i_track] << endl;
                            }
                        }

                        track_counter++;
                    }

                    track_pxA_old = track_pxA;
                    track_pyA_old = track_pyA;
                    track_pzA_old = track_pzA;
                    track_pxB_old = track_pxB;
                    track_pyB_old = track_pyB;
                    track_pzB_old = track_pzB;

                }

            }

            if(i_align == 0)
            {
                cout << "Good tracks (all HFT combinations): " << track_counter << endl;
                Int_t Good_TPC_tracks_real = 0;
                for(Int_t iTPC_sec = 0; iTPC_sec < N_TPC_sectors; iTPC_sec++)
                {
                    cout << "Sector: " << iTPC_sec << ", good (real) tracks: " << TPC_good_tracks_sectors[iTPC_sec] << endl;
                    h_Good_TPC_tracks_sector ->SetBinContent(iTPC_sec+1,TPC_good_tracks_sectors[iTPC_sec]);
                    Good_TPC_tracks_real += TPC_good_tracks_sectors[iTPC_sec];
                }
                cout << "Good tracks (real): " << Good_TPC_tracks_real << endl;
                cout << "" << endl;
            }


#if 1

            //------------------------------------------------------------------------------
            if(i_align == 0)
            {
                cout << "Start alignment, set minimizer method" << endl;
                // Choose method upon creation between:
                // kConjugateFR, kConjugatePR, kVectorBFGS,
                // kVectorBFGS2, kSteepestDescent


                //ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS ); // <-- good results
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugateFR ); // <-- used
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugatePR );
                //ROOT::Math::GSLMinimizer min( ROOT::Math::kSteepestDescent );
                ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS2 ); // <-- best results

                //ROOT::Math::GSLSimAnMinimizer min;

                //ROOT::Math::GSLMinimizer min( ROOT::Math::kMigradImproved);
                //ROOT::Math::GSLMultiMinimizer min( ROOT::Math::kVectorBFGS ); does not work, header not found
                //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);
                //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
                //ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigradImproved );
                //ROOT::Math::Minimizer min( ROOT::Minuit::kMigradImproved);

                min.SetMaxFunctionCalls(1000000); // 1000000
                min.SetMaxIterations(100000); // 100000
                min.SetTolerance(0.001); // 0.0001

                ROOT::Math::Functor f(&Align_TPC_Sector,6);

                Double_t step[6]     = {0.2,0.2,0.2,0.01,0.01,0.01};
                Double_t variable[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
                for(Int_t i_var = 0; i_var < 6; i_var++)
                {
                    variable[i_var] = 0.005; // 0.005
                    step[i_var] = 0.2; // 0.2
                }

                min.SetFunction(f);

                // Set the free variables to be minimized!
                min.SetVariable(0,"delta_x",variable[0], step[0]);
                min.SetVariable(1,"delta_y",variable[1], step[1]);
                min.SetVariable(2,"delta_z",variable[2], step[2]);
                min.SetVariable(3,"alpha",variable[3], step[3]);
                min.SetVariable(4,"beta",variable[4], step[4]);
                min.SetVariable(5,"gamma",variable[5], step[5]);

                //min.SetVariableLimits(3,0.0,0.0);
                //min.SetVariableLimits(4,0.0,0.0);

                min.SetFixedVariable(3,"",0.0);
                min.SetFixedVariable(4,"",0.0);

                // start minimization
                min.Minimize();

                const Double_t *min_params = min.X();
                if(min.ProvidesError()) cout << "Minimizer errors are provided" << endl;
                else cout << "No minimizer errors are provided" << endl;
                //const Double_t *err_params = min.Errors();

                Int_t status = min.Status();
                Int_t status_1 = status % 10;
                stringstream s;
                cout << "MINUIT returned status: " << status << ", status_1: " << status_1 << endl;
                switch(status_1)
                {
                case 1: cout << " (Covariance was made pos defined)" << endl; break;
                case 2: cout << " (Hesse is invalid)" << endl; break;
                case 3: cout << " (Edm is above max)" << endl; break;
                case 4: cout << " (Reached call limit)" << endl; break;
                case 5: cout << " (Some other failure)" << endl; break;
                default:
                    cout << " [unexpected status code]" << endl;
                }


                Double_t err_params[6][2];

                for(Int_t i_var = 0; i_var < 6; i_var++)
                {
                    min.GetMinosError(i_var,err_params[i_var][0],err_params[i_var][1],0);
                }


                //for(Int_t i = 0; i < 6; i++)
                //{
                //    cout << "par = " << i << ", value = " << min_params[i] << endl;
                //}

                delta_x       = min_params[0];
                delta_y       = min_params[1];
                delta_z       = min_params[2];
                rot_alpha     = min_params[3];
                rot_beta      = min_params[4];
                rot_gamma     = min_params[5];

                delta_x_err   = err_params[0][0];
                delta_y_err   = err_params[1][0];
                delta_z_err   = err_params[2][0];
                rot_alpha_err = err_params[3][0];
                rot_beta_err  = err_params[4][0];
                rot_gamma_err = err_params[5][0];


                rot4.SetToIdentity();
                rot4.RotateZ(rot_alpha);
                rot4.RotateX(rot_beta);
                rot4.RotateZ(rot_gamma);

                shift_vec.SetXYZ(delta_x,delta_y,delta_z);

                cout << "" << endl;
                cout << "------------------------------------------------" << endl;
                cout << "TPC sector: " << TPC_sec_align << endl;
                cout << "rotation matrix scheme for function Align_HFT_to_TPC" << endl;
                cout << "| x' |   | xx xy xz | | x |" << endl;
                cout << "| y' | = | yx yy yz | | y |" << endl;
                cout << "| z' |   | zx zy zz | | z |" << endl;
                cout << "|" << rot4.XX() << " " << rot4.XY() << " " << rot4.XZ() << "|" << endl;
                cout << "|" << rot4.YX() << " " << rot4.YY() << " " << rot4.YZ() << "|" << endl;
                cout << "|" << rot4.ZX() << " " << rot4.ZY() << " " << rot4.ZZ() << "|" << endl;
                cout << "shift vector = {" << shift_vec.X() << ", "  << shift_vec.Y() << ", " << shift_vec.Z() << "}" << endl;
                cout << "shift vector error = {" << delta_x_err << ", "  << delta_y_err << ", " << delta_z_err << "}" << endl;
                cout << "------------------------------------------------" << endl;
                cout << "" << endl;


                cout << "" << endl;
                cout << "------------------------------------------------" << endl;
                cout << "Save alignment parameters to output file" << endl;
                HistName = "./Data/Align_params/";
                HistName += "TPC_sec";
                HistName += TPC_sec_align;
                HistName += "_to_HFT_align_Itt2.root";
                TFile* Outputfile        = new TFile(HistName.Data(),"RECREATE");

                HistName = "Rot_matrix";
                rot4.Write(HistName.Data());

                HistName = "Shift_vec";
                shift_vec.Write(HistName.Data());
                cout << "------------------------------------------------" << endl;
                cout << "" << endl;



                //------------------------------------------------------------------------------------------------------------------------------------
                // Write alignment parameters to output file
                TString Table_format_file_name = "Params_TPC_align_sector_Itt2_";
                Table_format_file_name += TPC_sec_align;
                Table_format_file_name += ".txt";
                FILE* Table_data_file = fopen(Table_format_file_name.Data(),"w");
                fprintf(Table_data_file,"%s  \n"," ");
                fprintf(Table_data_file,"%s  \n","------------------------------------------------");
                fprintf(Table_data_file,"%s \t %i \n","TPC sector: ",TPC_sec_align);
                fprintf(Table_data_file,"%s  \n","rotation matrix scheme for function Align_HFT_to_TPC");
                fprintf(Table_data_file,"%s  \n","| x' |   | xx xy xz | | x |");
                fprintf(Table_data_file,"%s  \n","| y' | = | yx yy yz | | y |");
                fprintf(Table_data_file,"%s  \n","| z' |   | zx zy zz | | z |");
                fprintf(Table_data_file,"%s \t %E %s \t %E %s \t %E %s \n","|",rot4.XX()," ",rot4.XY()," ",rot4.XZ(),"|");
                fprintf(Table_data_file,"%s \t %E %s \t %E %s \t %E %s \n","|",rot4.YX()," ",rot4.YY()," ",rot4.YZ(),"|");
                fprintf(Table_data_file,"%s \t %E %s \t %E %s \t %E %s \n","|",rot4.ZX()," ",rot4.ZY()," ",rot4.ZZ(),"|");
                fprintf(Table_data_file,"%s \t %E %s \t %E %s \t %E %s \n","angles = {",rot_alpha,", ",rot_beta,", ",rot_gamma,"}");
                fprintf(Table_data_file,"%s \t %E %s \t %E %s \t %E %s \n","shift vector (cm) = {",shift_vec.X(),", ",shift_vec.Y(),", ",shift_vec.Z(),"}");
                fprintf(Table_data_file,"%s  \n","------------------------------------------------");
                fprintf(Table_data_file,"%s  \n"," ");
                //------------------------------------------------------------------------------------------------------------------------------------


            }
            //------------------------------------------------------------------------------
#endif


        }


        HistName = "c_Good_TPC_tracks_sector";
        TCanvas* c_Good_TPC_tracks_sector = new TCanvas(HistName.Data(),HistName.Data(),10,10,700,700);
        c_Good_TPC_tracks_sector->SetFillColor(10);
        c_Good_TPC_tracks_sector->SetTopMargin(0.03);
        c_Good_TPC_tracks_sector->SetBottomMargin(0.2);
        c_Good_TPC_tracks_sector->SetRightMargin(0.05);
        c_Good_TPC_tracks_sector->SetLeftMargin(0.15);

        c_Good_TPC_tracks_sector->cd(1);
        c_Good_TPC_tracks_sector->cd(1);
        c_Good_TPC_tracks_sector->cd(1)->SetTicks(1,1);
        c_Good_TPC_tracks_sector->cd(1)->SetGrid(0,0);
        c_Good_TPC_tracks_sector->cd(1)->SetFillColor(10);

        Draw_histogram(h_Good_TPC_tracks_sector,"h","TPC sector","counts",0,0,0,0,1.0,1.4,0.01,0.01,0.06,0.06,505);
        HistName = "p > ";
        sprintf(NoP,"%2.1f",min_TPC_track_mom);
        HistName += NoP;
        HistName += " GeV/c";
        plotTopLegend((char*)HistName.Data(),0.5,0.88,0.05,1,0.0,42,1,1);

        cout << "Plot 1D residual distributions for IST alignment" << endl;
        for(Int_t i_track = 0; i_track < 2; i_track++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                for(Int_t i_histo = 0; i_histo < N_HFT_TPC_match; i_histo++)
                {
                    c_pixel_line_dca_TPC[i_track]->cd(i_xyz+1+3*i_histo);

                    for(Int_t i_align = 0; i_align < 2; i_align++)
                    {
                        if(i_xyz < 2)
                        {
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Rebin(4);
                        }
                        if(i_xyz == 2)
                        {
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Rebin(8);
                        }
                    }


                    for(Int_t i_align = 0; i_align < 2; i_align++)
                    {
                        if(i_align == 0)
                        {
                            h_pixel_line_dca_TPC[1][i_histo][i_xyz][i_track]->GetXaxis()->SetRangeUser(-0.9,0.9);
                            h_pixel_line_dca_TPC[1][i_histo][i_xyz][i_track]->DrawCopy("");
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->DrawCopy("same");
                        }
                        if(i_align == 1 && i_histo >= 0)
                        {
                            //cout << "Plotting aligned distribution" << endl;
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetFillColor(kAzure-1);
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetLineColor(kAzure-1);
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->SetFillStyle(3001);
                            h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->DrawCopy("same");

                            if(i_xyz < 3)
                            {
                                //----------------------------------
                                // first fit
                                start_fit = -0.4;
                                stop_fit  = 0.4;
                                for(Int_t x = 0; x < 3; x++)
                                {
                                    GaussFit->ReleaseParameter(x);
                                    GaussFit->SetParError(x,0.0);
                                    GaussFit->SetParameter(x,0.0);
                                }
                                GaussFit->SetParameter(0,h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->GetBinContent(h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->FindBin(0)));
                                GaussFit->SetParameter(1,0.0);
                                GaussFit->SetParameter(2,0.1);
                                GaussFit->SetRange(start_fit,stop_fit);
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Fit("GaussFit","QN","",start_fit,stop_fit);
                                height = GaussFit->GetParameter(0);
                                mean   = GaussFit->GetParameter(1);
                                sigma  = GaussFit->GetParameter(2);
                                //----------------------------------



                                //----------------------------------
                                // second fit
                                start_fit = mean-2.0*sigma;
                                stop_fit  = mean+2.0*sigma;
                                for(Int_t x = 0; x < 3; x++)
                                {
                                    GaussFit->ReleaseParameter(x);
                                    GaussFit->SetParError(x,0.0);
                                    GaussFit->SetParameter(x,0.0);
                                }
                                GaussFit->SetParameter(0,height);
                                GaussFit->SetParameter(1,mean);
                                GaussFit->SetParameter(2,sigma*0.8);
                                GaussFit->SetRange(start_fit,stop_fit);
                                h_pixel_line_dca_TPC[i_align][i_histo][i_xyz][i_track]->Fit("GaussFit","QN","",start_fit,stop_fit);
                                GaussFit->SetLineWidth(2);
                                GaussFit->SetLineStyle(1);
                                GaussFit->SetLineColor(2);
                                GaussFit->SetRange(start_fit,stop_fit);
                                GaussFit->DrawCopy("same");

                                height = GaussFit->GetParameter(0);
                                mean   = GaussFit->GetParameter(1);
                                sigma  = fabs(GaussFit->GetParameter(2));

                                height_err = GaussFit->GetParError(0);
                                mean_err   = GaussFit->GetParError(1);
                                sigma_err  = GaussFit->GetParError(2);

                                HistName = "#mu = ";
                                sprintf(NoP,"%2.2f",mean*10.0);
                                HistName += NoP;
                                HistName += "#pm";
                                sprintf(NoP,"%2.2f",mean_err*10.0);
                                HistName += NoP;
                                HistName += " mm";
                                plotTopLegend((char*)HistName.Data(),0.25,0.88-0.13,0.055,1,0.0,42,1,1);

                                HistName = "#sigma = ";
                                sprintf(NoP,"%2.2f",sigma*10.0);
                                HistName += NoP;
                                HistName += "#pm";
                                sprintf(NoP,"%2.2f",sigma_err*10.0);
                                HistName += NoP;
                                HistName += " mm";
                                plotTopLegend((char*)HistName.Data(),0.245,0.82-0.13,0.055,1,0.0,42,1,1);

                                HistName = "TPC sector = ";
                                sprintf(NoP,"%i",TPC_sec_align);
                                HistName += NoP;
                                plotTopLegend((char*)HistName.Data(),0.245,0.76-0.13,0.055,1,0.0,42,1,1);
                                //----------------------------------
                            }

                        }
                    }

                    //----------------------------------
                    TLegend* leg_pixel_line_dca = new TLegend(0.23,0.8,0.46,0.93); // x1,y1,x2,y2
                    leg_pixel_line_dca->SetBorderSize(0);
                    leg_pixel_line_dca->SetFillColor(0);
                    leg_pixel_line_dca->SetTextSize(0.055);
                    leg_pixel_line_dca->SetTextFont(42);

                    leg_pixel_line_dca->AddEntry(h_pixel_line_dca_TPC[0][i_histo][i_xyz][i_track],"before align.","fl");
                    leg_pixel_line_dca->AddEntry(h_pixel_line_dca_TPC[1][i_histo][i_xyz][i_track],"after align.","fl");
                    leg_pixel_line_dca->Draw();
                    //----------------------------------

                    plotTopLegend((char*)label_IST[i_histo].Data(),0.78,0.885,0.055,1,0.0,42,1,1);

                }
            }


            //----------------------------------
            cout << "Save figures to output directory: " << fig_all_output_dir.Data() << endl;
            HistName = fig_all_output_dir.Data();
            HistName += "Align_sec_";
            HistName += TPC_sec_align;
            HistName += "_";
            HistName += c_pixel_line_dca_TPC[i_track]->GetName();
            HistName += out_all_format.Data();
            c_pixel_line_dca_TPC[i_track]->SaveAs(HistName.Data(),"");
            cout << "Saved canvas: " << HistName.Data() << endl;
            //----------------------------------
            //------------------------------------------------------------------------------
        }


    }
    //---------------------------------------------------------------------------------------------------------------

}