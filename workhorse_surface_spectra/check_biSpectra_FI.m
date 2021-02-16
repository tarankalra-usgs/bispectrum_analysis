clear all ; close all ; clc ; 
wh=fullfile('C:\Users\tkalra\Desktop\Observationaldata\data_FI\9921whp-cal.nc')  ; 

load('C:\Users\tkalra\Desktop\Observationaldata\matfiles\matfiles\skewness_steve.mat',....
    'Su_skewness','dn');
Su_skewness_FI=Su_skewness; 
% Tbr..

netcdf_load(wh)
nt1=1; nt2=2044;

d_spec_FI(:,:,:)=double(dspec(1,1,:,:,:));
v_spec_FI(:,:,:)=double(vspec(1,1,:,:));

band_width=0.015625  ;
f_FI=double(frequency(:,1));
df=band_width ;

isave=0 ;
d_spec_FI=d_spec_FI*1e-6; %convert from mm^2/Hz to m^2/Hz
ndirs=length(direction);
d_spec_FI=d_spec_FI/ndirs; %normalize to direction to get in m^2/Hz/degree

vspec_FI=(vspec*0.001).^2;
% p_spec=(p_spec*0.001).^2 ;

count_post_FI=0; count_neg_FI=0 ;

d_spec_FI(d_spec_FI<0)=0.0;
d_spec_FI(d_spec_FI>1e12)=0.0;


vspec_FI(vspec_FI<0)=0.0;
vspec_FI(vspec_FI>1e6)=0.0;

%
for t=nt1:nt2
%  % dspec_trans(t,:,:)=d_spec(:,:,t);
  if(Su_skewness_FI(t)>0.0)
      count_post_FI=count_post_FI+1;
      dspec_post_FI(:,:,count_post_FI)=d_spec_FI(:,:,t); % double(squeeze(d_spec(:,:,t)));
      vspec_post_FI(:,count_post_FI)=squeeze(vspec_FI(1,1,:,t)); 
      Su_post_FI(count_post_FI)=Su_skewness_FI(t); 
  else
      count_neg_FI=count_neg_FI+1;
      dspec_neg_FI(:,:,count_neg_FI)=d_spec_FI(:,:,t);
      vspec_neg_FI(:,count_neg_FI)=squeeze(vspec_FI(1,1,:,t)); 
      Su_neg_FI(count_neg_FI)=Su_skewness_FI(t); 
  end
end

% 
%  f_neg_MVCO=mean(mean_dspec_neg_MVCO(:,:),1);
%  f_post_MVCO=mean(mean_dspec_post_MVCO(:,:),1);
%  
% %  
  fv_neg_FI=mean(vspec_neg_FI,2);
  fv_post_FI=mean(vspec_post_FI,2);
% %  
% save('bispectral.mat','grom_all_1', 'grom_all_2','spd_all_1','spd_all_2',........
%                       'Su_post_all','Su_neg_all','f_MVCO','f_neg_MVCO',....);
%                       'f_post_MVCO','vspec_post);
      
save('bispectra_vspec_FI.mat','f_FI','fv_post_FI','fv_neg_FI')