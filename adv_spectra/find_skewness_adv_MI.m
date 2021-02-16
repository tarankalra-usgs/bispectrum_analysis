%% puv_proc_FI - Process ADV records from Fire Island ADV
% based on puv_proc_MVCO2011.m version of Feb 13, 2009
% last revised 5/1/2019

setup_nctoolbox
 clear all ; close all ; clc ;  

advsfn = fullfile('/media/taran/DATADRIVE2/Obs_data/Matanzes_Inlet/MI_data/11109vecs-a.nc'); % burstfile name
advbfn = fullfile('/media/taran/DATADRIVE2/Obs_data/Matanzes_Inlet/MI_data/11109vecb-cal.nc'); % burstfile name
 
% ncload(advbfn); % load the statistics file

zr = 0.4; % placeholder...need to check measurement elevation here

% you can see all of the variables in the burst file with:
% ncdisp(advbfn);
% % 
 time=double(ncread(advbfn,'time'));
 time2=ncread(advbfn,'time2') ; 
% % P_4023=ncread(advsfn,'P_4023'); 
% % 
% jt = time+time2/(3600*24*1000);
% dn = j2dn(time,time2);
%  
 

brange=ncreadatt(advbfn,'/','initial_instrument_height'); %THERE IS NO BRANGE IN DATA 
zoff=ncreadatt(advbfn,'/','ADVsampleVolumeOffset');
zr=brange-zoff; 

z_init = ncreadatt(advbfn,'u_1205','initial_instrument_height');

ap = 10.13 ;% std atmos. pressure (or a time series from nearby) [dBar]
p_z = ncreadatt(advbfn,'P_1','initial_instrument_height');
pdelz = p_z-z_init; % elevation diff. between velocity and pressure
zp = zr+pdelz; % elevation of pressure measurements [m] accounting for variable brange
%  
%depth = zr+pdelz+(P_1(gb)-ap); % time series of depth [decibars ~= meters]
 

%initial_instrument_height = double(ncreadatt(fn,'/','initial_instrument_height'))
%fs = ncreadatt(advbfn,'/','User_setup_Sampling_rate'); 
% Hard wired
fs=16 ;  
 
% Initliazign to max value of zero
%u_send_array(32768,1895)=0.0; v_send_array(32768,1895)=0.0; 
isave=10 ; 
if(isave==10)
%% process bursts with no QA/QC
nt1=1; nt2=1895; 
 
% set the time period in seconds for excluding infragravity wave band
t_up=20; t_low=4 ;
ifilt=1; 
isave=2;

for n=nt1:nt2 %nt2 
   %if(~isnan(depth(n)))
     bn = ncread(advbfn,'burst',n,1);      % this burst number from beginning...might just want to go from 1 to nb
     % jtb = double(ncread(advbfn,'time',n,1)+......
     %      double(ncread(advbfn,'time2',n,1))/(3600*24*1000))  ;
     jtb = double(ncread(advbfn,'time',[1 n],[1 1]))+......
           double(ncread(advbfn,'time2',[1 n],[1 1])/(3600*24*1000));
  
     dnsb = datestr(datenum(gregorian(jtb)))      ;
     jtb_rec(n)=jtb ;
      fprintf(1,'Burst %d at %s\n',bn,dnsb);
      u = ncread(advbfn,'u_1205',[1 n],[Inf 1])/100;
      v = ncread(advbfn,'v_1206',[1 n],[Inf 1])/100;
      w = ncread(advbfn,'w_1204',[1 n],[Inf 1])/100;
      p = ncread(advbfn,'P_1',[1 n],[Inf 1]); %Dont divide by 100 because it is already in dbar
%     
      depth= mean(p)+p_z; % depth modification 
      
      % remove nans
      u=u(~isnan(u))'; 
      v=v(~isnan(v))'; 
      
      if(ifilt==1)
      % DETREND U, V
      u_detrend=detrend(u); 
      v_detrend=detrend(v); 
      
      % Save Mean current velocity
      if(isave==4)
      u_mean(:,n)=u-u_detrend ;
      v_mean(:,n)=v-v_detrend ; 
       save('mean_current_adv.mat','u_mean','v_mean'); 
      end 
      % BAND PASS FILTER ; 
      u_band=iwavesbp(u_detrend, fs, t_up, t_low); %0.04 Hz in seconds 
      v_band=iwavesbp(v_detrend, fs, t_up, t_low); 
      
      % MEDIAN FILTER
      u_med=medfilt(u_band, 9) ; %smooth
      v_med=medfilt(v_band, 9) ;
      
      u_send=u_med;
      v_send=v_med; 
      else
       u_send=u; v_send=v; % that is no filtering
      end 
%       % PCA STATS ;
      [sd1 az1 sd2 az2]=pcastats(u_send*100,v_send*100,50,1);

      UBS(n) = ubstatsr( u_send, v_send, fs );
      ur_maj_rot=UBS(n).ur; % major rotated
      vr_min_rot=UBS(n).vr; % minor 
      CRS_skew(n)=UBS(n).ur_sk; 
      ang_rot(n)=UBS(n).maj_az; 
      % 
%      FIND SKEWNESS 
      Su_skewness(n)=mean(ur_maj_rot.^3)/(std(ur_maj_rot)).^3;  % Eqn. 5 in the skewness calculation 
       
      hilbert_asym=imag(hilbert(ur_maj_rot)) ; 
      Au_skewness(n)=mean(hilbert_asym.^3)/(std(ur_maj_rot)).^3 ; 
      
      ur_bar(n)=mean(ur_maj_rot); 
      ur_cube(n)=mean(ur_maj_rot.^3); 
     % get ubr, Hrms
      
      PUV(n) = puvq(p, (u), (v), depth, zp, zr, fs, 1050, 1030., 0.04, 1/6);
      
      kh = qkhfs( 2*pi/PUV(n).Tr, depth );
      Tr(n)=PUV(n).Tr;
      omega_br(n)=PUV(n).omegar ;  
      Ubr(n)=PUV(n).ubr ; 
      Hrmsu(n)=PUV(n).Hrmsu; 
     
      k(n) = kh./depth;
      Ursell(n) = 0.75*0.5*PUV(n).Hrmsu*k(n)./(kh.^3); % RRvR Eqn. 6.
     
      depth_nominal(n)=depth;
     % jtb_rec(n)=jtb ;
       if(isave==1)
         save('MI_matfiles/skewness_steve_MI_depth_nominal.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Tr','k','omega_br','Ubr','depth_nominal',....
            'Ursell','ur_bar','ur_cube','ang_rot','Au_skewness')
       end
       
       if(isave==2)
         if(length(u_send)==32768)
          u_send_array(:,n)=u_send;
          v_send_array(:,n)=v_send;
         else
          u_send_array(32768,n)=0.0; 
          v_send_array(32768,n)=0.0; 
         end 
       % ur_maj_rot_array(:,n)=ur_maj_rot; 
       % vr_min_rot_array(:,n)=vr_min_rot; 
       % save('MI_matfiles/skewness_orbital_array_MVCO.mat','Su_skewness','ur_maj_rot_array','vr_min_rot_array',........
%        %    'Hrmsu','Ubr','Ursell','dn','jtb_rec')
% 	    save('/media/taran/DATADRIVE2/Obs_data/matfiles/matfiles_post_OSM/skewness_MI_netvel.mat',.....
% 	      'Su_skewness','u_send_array','v_send_array','Hrmsu','Ubr','Ursell','jtb_rec')
       end 
       if(isave==3)
         save('MI_matfiles/mean_current.mat','u_mean','v_mean')
       end 
       if(ifilt==0)
         save('MI_matfiles/skewness_steve_nofilters.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
            'Hrmsu','Tr','k','omega_br','Ubr','depth','Ursell','dn','jtb_rec','ur_bar','ur_cube','ang_rot','Au_skewness')
       end 
     
  % end
  end
end
%  