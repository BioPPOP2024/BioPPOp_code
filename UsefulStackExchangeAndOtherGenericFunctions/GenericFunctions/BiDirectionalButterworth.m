function [ filteredData ] = BiDirectionalButterworth( data, fs, fc, order )
% Lowpass filter column data using a two pass bi-directional butterworth filter. 
%
% NOTE: The specified cut-off frequency corresponds to the 3dB point for singe pass
% IIR filter. The bidirectional 2 pass nature of the filtering means this
% will correspond to 6dB for the combined filter.
% The true 3dB cut-off frequency for the bi-directional filter will be
% slightly less than this.
%
% Usage:  [ filteredData ] = BiDirectionalButterworth( data, fs, fc, order )
%   data       - array of data of size which must not contain any NaNs
%   fs         - sampling frequency of data
%   fc         - cutOff frequency for 1 pass filter
%   order      - filter order for 1 pass butterworth (default is 2)
%
% Returns:
%   filteredData - same size as data
%
% This uses the butter and filtfilt commands from the signal processing
% toolbox.
%
% David Redmill
% Email: David.Redmill@bristol.ac.uk
% University of Bristol, UK. 
% 25-11-2011
%

if(nargin ~= 4)
    error('Invalid arguments : Usage:  [ filteredData ] = BiDirectionalButterworth( data, fs, fc, order ) ');
end

if (isnan(sum(sum(sum(data)))))
    error('Invalid data invalid : contains at least 1 NaN ');
end


[b, a] = butter(order, 2*fc/fs);

filteredData = filtfilt(b, a, data);

%
% h1 = dfilt.df2(b,a);
% h2 = cascade(h1,h1);
% hfvt = fvtool(h1,h2);
%


end

