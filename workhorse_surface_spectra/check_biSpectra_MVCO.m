clear all ; close all ; clc ; 

nt1=1; nt2_adv=492;
load('C:\Users\tkalra\Desktop\Observationaldata\MVCO\MVCO\MVCO_matfiles\skewness_steve_MVCO.mat','Su_skewness','dn');
Su_skewness_MVCO=Su_skewness; 

% Tbr..
% ADCP wave data on quadpod
fn = 'C:\Users\tkalra\Desktop\Observationaldata\MVCO\MVCO\data_netcdf\10571whVp-cal.nc'
var_list = {'time','time2','wh_4061','hght_18','wvdir','wp_peak'};
for i=1:length(var_list)
   s = sprintf('%s = squeeze(ncread(fn, ''%s''));',var_list{i},var_list{i});
   eval(s)
end
%dn_wh = j2dn(time,time2);
%initial_instrument_height = double(ncreadatt(fn,'/','initial_instrument_height'))

netcdf_load(fn);

nt1=1; nt2_wh=795;
% 
band_width=0.015625  ;
f_MVCO=double(frequency(:,1));
df=band_width ;

d_spec_MVCO(:,:,:)=double(dspec(1,1,:,:,:));
%v_spec(:,:)=double(vspec(1,1,:,:));
%p_spec(:,:)=double(pspec(1,1,:,:));
%s_spec(:,:)=double(sspec(1,1,:,:));
% 
% isave=0 ;
d_spec_MVCO=d_spec_MVCO*1e-6; %convert from mm^2/Hz to m^2/Hz
ndirs=length(direction);
d_spec_MVCO=d_spec_MVCO/ndirs; %normalize to direction to get in m^2/Hz/degree

d_spec_MVCO(d_spec_MVCO<0)=0.0; 
d_spec_MVCO(d_spec_MVCO>1e12)=0.0;
%
nt2_MVCO=492; % nt2 corresponds to shorter array 
 
Su_skewness_new_MVCO(1:nt2_wh)=NaN ;
Su_skewness_new_MVCO(121:nt2_adv+120)=Su_skewness_MVCO; 

% Tbr..

f_MVCO=double(frequency(:,1));
df=band_width ;

isave=0 ;

vspec_MVCO=(vspec*0.001).^2;
% p_spec=(p_spec*0.001).^2 ;

count_post_MVCO=0; count_neg_MVCO=0 ;

vspec_MVCO(vspec_MVCO<0)=0.0;
vspec_MVCO(vspec_MVCO>1e6)=0.0;

%
for t=nt1:nt2_wh
%  % dspec_trans(t,:,:)=d_spec(:,:,t);
  if(Su_skewness_new_MVCO(t)>0.0)
      count_post_MVCO=count_post_MVCO+1;
      dspec_post_MVCO(:,:,count_post_MVCO)=d_spec_MVCO(:,:,t); % double(squeeze(d_spec(:,:,t)));
      vspec_post_MVCO(:,count_post_MVCO)=squeeze(vspec_MVCO(1,1,:,t)); 
      Su_post_MVCO(count_post_MVCO)=Su_skewness_new_MVCO(t); 
  else
      count_neg_MVCO=count_neg_MVCO+1;
      dspec_neg_MVCO(:,:,count_neg_MVCO)=d_spec_MVCO(:,:,t);
      vspec_neg_MVCO(:,count_neg_MVCO)=squeeze(vspec_MVCO(1,1,:,t)); 
      Su_neg_MVCO(count_neg_MVCO)=Su_skewness_new_MVCO(t); 
  end
end

% 
%  f_neg_MVCO=mean(mean_dspec_neg_MVCO(:,:),1);
%  f_post_MVCO=mean(mean_dspec_post_MVCO(:,:),1);
%  
% %  
  fv_neg_MVCO=mean(vspec_neg_MVCO,2);
  fv_post_MVCO=mean(vspec_post_MVCO,2);
% %  
% save('bispectral.mat','grom_all_1', 'grom_all_2','spd_all_1','spd_all_2',........
%                       'Su_post_all','Su_neg_all','f_MVCO','f_neg_MVCO',....);
%                       'f_post_MVCO','vspec_post);
      
save('bispectra_vspec_MVCO.mat','f_MVCO','fv_post_MVCO','fv_neg_MVCO')