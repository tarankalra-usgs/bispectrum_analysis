clear all ; close all ; clc ; 
wh=fullfile('C:\Users\tkalra\Desktop\Observationaldata\data_FI\9921whp-cal.nc')  ; 

load('C:\Users\tkalra\Desktop\Observationaldata\matfiles\matfiles\skewness_steve.mat',....
    'Su_skewness','dn');
Su_skewness_FI=Su_skewness; 
% Tbr..

 %[M, I]=max(Su_skewness_FI)
 
 netcdf_load(wh)
 
 
 vspec_FI=(vspec*0.001).^2;
% % p_spec=(p_spec*0.001).^2 ;
% 
 count_post_FI=0; count_neg_FI=0 ;
% 

 %vspec_FI(vspec_FI<0)=0.0;
 %vspec_FI(vspec_FI>1e6)=0.0;
% 

  % second highest Su, 153 -- valentines day\
 nt=153 ;
 %[M, I]=max(Su_skewness_FI);
 vspec_2nd_int=squeeze(double(vspec_FI(1,1,:,153))); 
 
 nt=1300; 
 %[M, I]=min(Su_skewness_FI);
 
 vspec_1st_int=squeeze(double(vspec_FI(1,1,:,nt))) ; 
% nt1=1; nt2=2044;
% 
% d_spec_FI(:,:,:)=double(dspec(1,1,:,:,:));
% v_spec_FI(:,:,:)=double(vspec(1,1,:,:));
% 
 band_width=0.015625  ;
 f_FI=double(frequency(:,1));
 df=band_width ;
% 
% isave=0 ;
% d_spec_FI=d_spec_FI*1e-6; %convert from mm^2/Hz to m^2/Hz
% ndirs=length(direction);
% d_spec_FI=d_spec_FI/ndirs; %normalize to direction to get in m^2/Hz/degree
% 
% %
% for t=nt1:nt2
% %  % dspec_trans(t,:,:)=d_spec(:,:,t);
%   if(Su_skewness_FI(t)>0.0)
%       count_post_FI=count_post_FI+1;
%       dspec_post_FI(:,:,count_post_FI)=d_spec_FI(:,:,t); % double(squeeze(d_spec(:,:,t)));
%       vspec_post_FI(:,count_post_FI)=squeeze(vspec_FI(1,1,:,t)); 
%       Su_post_FI(count_post_FI)=Su_skewness_FI(t); 
%   else
%       count_neg_FI=count_neg_FI+1;
%       dspec_neg_FI(:,:,count_neg_FI)=d_spec_FI(:,:,t);
%       vspec_neg_FI(:,count_neg_FI)=squeeze(vspec_FI(1,1,:,t)); 
%       Su_neg_FI(count_neg_FI)=Su_skewness_FI(t); 
%   end
% end
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

save('bispectral_FI_intervals.mat', 'f_neg_Su_FI_1st', 'f_post_Su_FI_2nd');

plot(f_FI, f_post_Su_FI_2nd)
hold on 
plot(f_FI, f_neg_Su_FI_1st,'ro')
%       
% save('bispectra_vspec_FI.mat','f_FI','fv_post_FI','fv_neg_FI')