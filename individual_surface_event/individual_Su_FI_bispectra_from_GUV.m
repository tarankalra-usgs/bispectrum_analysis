clear all ; close all ; clc ; 
wh=fullfile('C:\Users\tkalra\Desktop\Observationaldata\data_FI\9921whp-cal.nc')  ; 

load('C:\Users\tkalra\Desktop\Observationaldata\paper_figs\matfiles\filter_step1_154.mat','u','v','u_band',....
                            'v_band','u_send','v_send',.......
                            'ur_maj_rot','PUV_raw', 'PUV_detrend', 'PUV_band_noIG', 'PUV_band_withIG', 'PUV_med');
                        
load('C:\Users\tkalra\Desktop\Observationaldata\matfiles\matfiles\skewness_steve.mat',....
    'Su_skewness','dn');
Su_skewness_FI=Su_skewness; 
% Tbr..

 %[M, I]=max(Su_skewness_FI)
 
 netcdf_load(wh)
%  
%  
%  vspec_FI=(vspec*0.001).^2;
% % % p_spec=(p_spec*0.001).^2 ;
% % 
%  count_post_FI=0; count_neg_FI=0 ;
% % 
% 
%  vspec_FI(vspec_FI<0)=0.0;
%  vspec_FI(vspec_FI>1e6)=0.0;
% 
 f_FI=double(frequency(:,1));

 % Multiply with frequency to get units in m^2/s from m^/(s-Hz) ;
 
  % second highest Su, 153 -- valentines day\
 nt=153 ;
 %[M, I]=max(Su_skewness_FI);
 vspec_2nd_int=squeeze(double(PUV_raw(1,1,:,153))); 
 vspec_2nd_int=vspec_2nd_int.*f_FI; 
 
 nt=1424; 
 %[M, I]=min(Su_skewness_FI);
 
 vspec_1st_int=squeeze(double(PUV_raw(1,1,:,nt))) ; 
 vspec_1st_int=vspec_1st_int.*f_FI; 
 
% 
 band_width=0.015625  ;
 df=band_width ;
% 
% % 
  f_neg_Su_FI_1st=vspec_1st_int(:);
  f_post_Su_FI_2nd=vspec_2nd_int(:); % mean(mean_dspec_post_MVCO(:,:),1);
% %  
% % %  
%   fv_neg_FI=mean(vspec_neg_FI,2);
%   fv_post_FI=mean(vspec_post_FI,2);
% % %  
% % save('bispectral.mat','grom_all_1', 'grom_all_2','spd_all_1','spd_all_2',........
% %                       'Su_post_all','Su_neg_all','f_MVCO','f_neg_MVCO',....);
% %                       'f_post_MVCO','vspec_post);

save('bispectral_FI_intervals_GUV.mat', 'f_neg_Su_FI_1st', 'f_post_Su_FI_2nd');
%       ');

plot(f_neg_