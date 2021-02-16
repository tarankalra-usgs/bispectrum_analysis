clear all ; close all ; clc ; 
nt1_MI=1; nt2_MI=728;
load('C:\Users\tkalra\Desktop\Observationaldata\Matanzes_Inlet\Matanzes_Inlet\matfiles\skewness_steve_MI_depth_corrected.mat',.....
    'Su_skewness','dn');

% FETCH VSPEC DATA
wh_MI=fullfile('C:\Users\tkalra\Desktop\Observationaldata\Matanzes_Inlet\Matanzes_Inlet\MI_data\11101whVp-cal.nc')  ;
Su_skewness_MI=Su_skewness; 
% Tbr..

netcdf_load(wh_MI)

d_spec_MI(:,:,:)=double(dspec(1,1,:,:,:));
v_spec_MI(:,:,:)=double(vspec(1,1,:,:));

band_width=0.00390625  ; 

f_MI=double(frequency(:,1));
df=band_width ;

isave=0 ;
d_spec_MI=d_spec_MI*1e-6; %convert from mm^2/Hz to m^2/Hz
ndirs=length(direction);
d_spec_MI=d_spec_MI/ndirs; %normalize to direction to get in m^2/Hz/degree

vspec_MI=(vspec*0.001).^2;
% p_spec=(p_spec*0.001).^2 ;

count_post_MI=0; count_neg_MI=0 ;

d_spec_MI(d_spec_MI<0)=0.0;
d_spec_MI(d_spec_MI>1e12)=0.0;


vspec_MI(vspec_MI<0)=0.0;
vspec_MI(vspec_MI>1e6)=0.0;

%
for t=nt1_MI:nt2_MI
%  % dspec_trans(t,:,:)=d_spec(:,:,t);
  if(Su_skewness_MI(t)>0.0)
      count_post_MI=count_post_MI+1;
      dspec_post_MI(:,:,count_post_MI)=d_spec_MI(:,:,t); % double(squeeze(d_spec(:,:,t)));
      vspec_post_MI(:,count_post_MI)=squeeze(vspec_MI(1,1,:,t)); 
      Su_post_MI(count_post_MI)=Su_skewness_MI(t); 
  else
      count_neg_MI=count_neg_MI+1;
      dspec_neg_MI(:,:,count_neg_MI)=d_spec_MI(:,:,t);
      vspec_neg_MI(:,count_neg_MI)=squeeze(vspec_MI(1,1,:,t)); 
      Su_neg_MI(count_neg_MI)=Su_skewness_MI(t); 
  end
end

% 
%  f_neg_MVCO=mean(mean_dspec_neg_MVCO(:,:),1);
%  f_post_MVCO=mean(mean_dspec_post_MVCO(:,:),1);
%  
% %  
  fv_neg_MI=mean(vspec_neg_MI,2);
  fv_post_MI=mean(vspec_post_MI,2);
% %  
% save('bispectral.mat','grom_all_1', 'grom_all_2','spd_all_1','spd_all_2',........
%                       'Su_post_all','Su_neg_all','f_MVCO','f_neg_MVCO',....);
%                       'f_post_MVCO','vspec_post);
      
save('bispectra_vspec_MI.mat','f_MI','fv_post_MI','fv_neg_MI')