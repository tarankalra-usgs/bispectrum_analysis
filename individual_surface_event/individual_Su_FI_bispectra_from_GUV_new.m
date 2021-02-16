clear all ; close all ; clc; 


load('C:\Users\tkalra\Desktop\Observationaldata\paper_figs\matfiles\filter_step1_154.mat','u','v','u_band',....
                            'v_band','u_send','v_send',.......
                            'ur_maj_rot','PUV_raw', 'PUV_detrend', 'PUV_band_noIG', 'PUV_band_withIG', 'PUV_med');
%load('C:\Users\tkalra\Desktop\Observationaldata\paper_figs\matfiles\filter_step1_154.mat','u','v','u_band',....
%                            'v_band','u_send','v_send','PUV_1',.......
%                            'ur_maj_rot', 'PUV_detrend', 'PUV_band_noIG', 'PUV_band_withIG', 'PUV_med');                       

%load('C:\Users\tkalra\Desktop\Observationaldata\paper_figs\paper_figs_old\make_waveforms\puv_proc_FI_iwaves_taran_IG.mat');

%unet=sqrt(u.*u+v.*v);
%unet_band=sqrt(u_band.*u_band+v_band.*v_band) ; 
f_1=1.0/20.0; 
f_2=1.0/4.0 ;
nt=154; 


f=(PUV_raw(nt).f);
vspec_2nd_int=PUV_raw(nt).Guv ; % mean(mean_dspec_post_MVCO(:,:),1);

save('bispectral_FI_intervals_GUV_new.mat', 'f','vspec_2nd_int')