
#include "TRD_ion_tail.h"

void TRD_ion_tail()
{
    printf("TRD_ion_tail started, all systems nominal. \n");

    Load_data();

    tg_pad_response ->Draw();

    Double_t impact_angle             = 107.0;//108.0;
    Double_t Lorentz_angle            = 0.1;//-7.5;
    Double_t Drift_vel_ratio          = 1; // 0.8
    Double_t Lorentz_angle_pre_corr   = -7.5; // -7.5

    Draw_TRD_detector_2D();

    for(Double_t impact_angle = 65.0; impact_angle < 115.0; impact_angle += 1.0)
    {

    Create_TRD_track(impact_angle*TMath::DegToRad(),Lorentz_angle*TMath::DegToRad(),Drift_vel_ratio,Lorentz_angle_pre_corr*TMath::DegToRad());
    Make_ion_tail_convolution();
    Create_TRD_track(impact_angle*TMath::DegToRad(),Lorentz_angle*TMath::DegToRad(),Drift_vel_ratio,Lorentz_angle_pre_corr*TMath::DegToRad());
    Reconstruct_tracklet(Lorentz_angle*TMath::DegToRad());

    straight_line_fits(impact_angle*TMath::DegToRad());

    }

    fit_TG();

    Draw_delta_angle();
}
