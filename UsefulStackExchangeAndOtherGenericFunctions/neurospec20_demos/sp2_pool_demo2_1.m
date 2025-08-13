% Script sp2_pool_demo2_1.m
% Script to demonstrate use of NeuroSpec2.0
%  Pooled spectra and pooled coherence analysis
%
% Includes second level pooling on two populations.
%
%
% Copyright 2008 2018, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%

%--------------------------------------------------------------------------
% User parameters

% Number of sets of data to pool in each population
pool_tot1=10;
pool_tot2=10;

% Parameters of sine wave to include in first population: Frequency & Amplitude
% Set amplitude to zero to have both populations the same
s_freq=10;
s_amp=0.1;
%s_amp=0.0;

%--------------------------------------------------------------------------

% Number of samples in each set.
samp_set_tot=11000;
% For purposes of analysis: assume sampling rate is 1000/sec
samp_rate=1000;

% Define filter parameters
% Filtering done using 11 point moving average window.
a=1;
b=(ones(1,11)/11)';

% Set power for segment length as 10.
% T=2^10 = 1024.
seg_pwr=10;

% No options used here
opt_str='';

%------------------------------------------------------------------------------
% Set 1
% Loop round all data sets
for ind=1:pool_tot1
  % Generate data set
  dat=randn(samp_set_tot,2);
  % Filter column 1.
   dat_filt=filter(b,a,dat(:,1));
  % Column 3 is average of filtered version of column 1 and column 2.
  dat(:,3)=0.5*(dat_filt+dat(:,2));

  % Add in additional sine component: Set parameters above.
    t=1/samp_rate*(1:samp_set_tot)';
    sin_dat=s_amp*sin(2*pi*t*s_freq);
     dat(:,1)=dat(:,1)+sin_dat;
    dat(:,3)=dat(:,3)+sin_dat;

  % Normalise variance
  % Could also do this using 'n' option with sp2a2_m1
  dat(:,1)=dat(:,1)/std(dat(:,1));
  dat(:,3)=dat(:,3)/std(dat(:,3));

  % Process columns 1 and 3, these are Correlated & Filtered
  [f1(:,:,ind),t1(:,:,ind),cl1(ind),sc1(:,:,ind)] = sp2a2_m1(0,dat(:,1),dat(:,3),samp_rate,seg_pwr,opt_str);
  cl1(ind).what=['Set: ',num2str(ind)];

  % Pooled analysis
  if (ind==1)
    % Separate call for first set, creates new pooled analysis.
    [plf1,plv1]=pool_scf(sc1(:,:,ind),cl1(ind)); 
  else
    % Pass pooled variables as arguments for further sets.
    [plf1,plv1]=pool_scf(sc1(:,:,ind),cl1(ind),plf1,plv1);
  end  
end


%------------------------------------------------------------------------------
% Set 2
% Loop round all data sets
for ind=1:pool_tot2
  % Generate data set
  dat=randn(samp_set_tot,2);
  % Filter column 1.
   dat_filt=filter(b,a,dat(:,1));
  % Column 3 is average of filtered version of column 1 and column 2.
  dat(:,3)=0.5*(dat_filt+dat(:,2));

  % Normalise variance
  % Could also do this using 'n' option with sp2a2_m1
  dat(:,1)=dat(:,1)/std(dat(:,1));
  dat(:,3)=dat(:,3)/std(dat(:,3));

  % Process columns 1 and 3, these are Correlated & Filtered
  [f2(:,:,ind),t2(:,:,ind),cl2(ind),sc2(:,:,ind)] = sp2a2_m1(0,dat(:,1),dat(:,3),samp_rate,seg_pwr,opt_str);
  cl2(ind).what=['Set: ',num2str(ind)];

  % Pooled analysis
  if (ind==1)
    % Separate call for first set, creates new pooled analysis.
    [plf2,plv2]=pool_scf(sc2(:,:,ind),cl2(ind)); 
  else
    % Pass pooled variables as arguments for further sets.
    [plf2,plv2]=pool_scf(sc2(:,:,ind),cl2(ind),plf2,plv2);
  end  
end

%------------------------------------------------------------------------------
% Plotting first level pooling
freq=75;
ch_max=1;
lag_tot=100;
lag_neg=50;
chi_max=0;% Will auto scale

% Process pooled spectral coefficients & plot NB includes 4 outputs, sc  for 2nd level pooling.
[f1a,t1a,cl1a,sc1a]=pool_scf_out(plf1,plv1);
figure
cl1a.what=['Pooled analysis, sets: ',num2str(pool_tot1)];
psp2_pool6(f1a,t1a,cl1a,freq,lag_tot,lag_neg,ch_max,chi_max)


% Process pooled spectral coefficients & plot
[f2a,t2a,cl2a,sc2a]=pool_scf_out(plf2,plv2);
figure
cl2a.what=['Pooled analysis, sets: ',num2str(pool_tot2)];
psp2_pool6(f2a,t2a,cl2a,freq,lag_tot,lag_neg,ch_max,chi_max)

%------------------------------------------------------------------------------
% Second level pooling
% Now we are using spectral coefficients in sc1a and sc2a	(which are pooled)

% Pooling of the two sets
[plf3,plv3]=pool_scf(sc1a,cl1a); 
[plf3,plv3]=pool_scf(sc2a,cl2a,plf3,plv3); 

% Process pooled spectral coefficients & plot
[f3a,t3a,cl3a]=pool_scf_out(plf3,plv3);
cl3a.what=['Comparison of two populations'];

%KAT ADDITION
figure
[f5,cl5]=sp2_compf(sc1a,cl1a,1,sc2a,cl2a,1)
psp_compf1(f5,cl5,freq)


% Plotting
figure
psp2_pool6(f3a,t3a,cl3a,freq,lag_tot,lag_neg,ch_max,chi_max)

% Include chi-square tests on comparison of pooled spectra
figure
psp2_pool8_chi3(f3a,cl3a,freq,ch_max,chi_max)

