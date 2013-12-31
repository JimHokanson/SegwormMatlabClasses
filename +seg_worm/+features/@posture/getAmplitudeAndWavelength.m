function getAmplitudeAndWavelength(obj, theta_d,sx,sy,worm_lengths)
%
%
%   Inputs
%   =======================================================================
%   theta_d      : worm orientation based on fitting to an ellipse, in
%                   degrees
%   xx           : [49 x n_frames]
%   yy           : [49 x n_frames]
%   worm_lengths : [1 x n_frames], total length of each worm
%
%
%   Outputs
%   =======================================================================
%   amplitude    :
%       .max       - [1 x n_frames] max y deviation after rotating major axis to x-axis
%       .ratio     - [1 x n_frames] ratio of y-deviations (+y and -y) with worm centered
%                    on the y-axis, ratio is computed to be less than 1
%   wavelength   :
%       .primary   - [1 x n_frames]
%       .secondary - [1 x n_frames] this might not always be valid, even 
%                     when the primary wavelength is defined
%   trackLength  : [1 x n_frames]
%
%   
%   Old Name: getAmpWavelength.m
%   TODO: This function was missing from some of the original code
%   distributions. I need to make sure I upload it.
%
%
%   Nature Methods Description
%   =======================================================================
%   Amplitude. 
%   ------------------
%   Worm amplitude is expressed in two forms: a) the maximum
%   amplitude found along the worm body and, b) the ratio of the maximum
%   amplitudes found on opposing sides of the worm body (wherein the smaller of
%   these two amplitudes is used as the numerator). The formula and code originate
%   from the publication “An automated system for measuring parameters of
%   nematode sinusoidal movement”6.
%   The worm skeleton is rotated to the horizontal axis using the orientation of the
%   equivalent ellipse and the skeleton’s centroid is positioned at the origin. The
%   maximum amplitude is defined as the maximum y coordinate minus the minimum
%   y coordinate. The amplitude ratio is defined as the maximum positive y coordinate
%   divided by the absolute value of the minimum negative y coordinate. If the
%   amplitude ratio is greater than 1, we use its reciprocal.
%
%   Wavelength
%   ------------------------
%   Wavelength. The worm’s primary and secondary wavelength are computed by
%   treating the worm’s skeleton as a periodic signal. The formula and code
%   originate from the publication “An automated system for measuring
%   parameters of nematode sinusoidal movement”6. The worm’s skeleton is
%   rotated as described above for the amplitude. If there are any
%   overlapping skeleton points (the skeleton’s x coordinates are not
%   monotonically increasing or decreasing in sequence -- e.g., the worm is
%   in an S shape) then the shape is rejected, otherwise the Fourier
%   transform computed. The primary wavelength is the wavelength associated
%   with the largest peak in the transformed data. The secondary wavelength
%   is computed as the wavelength associated with the second largest
%   amplitude (as long as it exceeds half the amplitude of the primary
%   wavelength). The wavelength is capped at twice the value of the worm’s
%   length. In other words, a worm can never achieve a wavelength more than
%   double its size.
%
%   Tracklength
%   -----------------------------
%   Track Length. The worm’s track length is the range of the skeleton’s
%   horizontal projection (as opposed to the skeleton’s arc length) after
%   rotating the worm to align it with the horizontal axis. The formula and
%   code originate from the publication “An automated system for measuring
%   parameters of nematode sinusoidal movement”.


% Code based on:
% ------------------------------------------------
% BMC Genetics, 2005
% C.J. Cronin, J.E. Mendel, S. Mukhtar, Young-Mee Kim, R.C. Stirb, J. Bruck,
% P.W. Sternberg
% "An automated system for measuring parameters of nematode
% sinusoidal movement" BMC Genetics 2005, 6:5


N_POINTS_FFT   = 512;
HALF_N_FFT     = N_POINTS_FFT/2;
MIN_DIST_PEAKS = 5; %NOTE: Unfortunately the distance is in normalized
%frequency units (indices really), not in real frequency units
WAVELENGTH_PCT_MAX_CUTOFF = 0.5; %TODO: Describe
WAVELENGTH_PCT_CUTOFF = 2; %TODO: Describe

assert(size(sx,1) <= N_POINTS_FFT,'# of points used in the FFT must be more than the # of points in the skeleton')

theta_r = theta_d*pi/180;

%Unrotate worm
%-----------------------------------------------------------------
wwx = bsxfun(@times,sx,cos(theta_r)) + bsxfun(@times,sy,sin(theta_r));
wwy = bsxfun(@times,sx,-sin(theta_r)) + bsxfun(@times,sy,cos(theta_r));

%Subtract mean
%-----------------------------------------------------------------
wwx = bsxfun(@minus,wwx,mean(wwx,1));
wwy = bsxfun(@minus,wwy,mean(wwy,1));

%Amplitude Calculations
%--------------------------------------------------------------------------
% We have rotated the coordinates around the origin at the axes, not at the
% centroid of the worm. Now center at the origin by subtracting the
% centroid coordinates
%
%   ???? - It is surprising that you don't subtract the mean, then rotate,
%   this is a bit surprising ..., this will be identical here but it is
%   unclear if the theta calculation is identical ... (see eccentricity
%   function)   

% Calculate track amplitude
%--------------------------------------------------------------------------
amp1 = max(wwy,[],1);
amp2 = min(wwy,[],1);
amplitude.max   = amp1 - amp2;
amp2 = abs(amp2);
amplitude.ratio = min(amp1,amp2)./max(amp1,amp2);

% Calculate track length
%--------------------------------------------------------------------------
%NOTE: This is the x distance after rotation, and is different from the worm
%length which follows the skeleton. This will always be smaller than the
%worm length.
trackLength = max(wwx,[],1) - min(wwx,[],1);

%Wavelength Calculation
%--------------------------------------------------------------------------
dwwx = diff(wwx,1,1);

%Does the sign change? This is a check to make sure that the change in x is
%always going one way or the other. Is sign of all differences the same as
%the sign of the first, or rather, are any of the signs not the same as the
%first sign, indicating a "bad worm orientation".
bad_worm_orientation = any(bsxfun(@ne,sign(dwwx),sign(dwwx(1,:))),1);

n_frames = length(bad_worm_orientation);

p_wavelength = NaN(1,n_frames);
s_wavelength = NaN(1,n_frames);

%NOTE: Right now this varies from worm to worm which means the spectral
%resolution varies as well from worm to worm
spatial_sampling_frequency = (size(wwx,1)-1)./trackLength;

ds = 1./spatial_sampling_frequency;

frames_to_calculate = find(~bad_worm_orientation);

for iFrame = 1:length(frames_to_calculate)
    
    cur_frame = frames_to_calculate(iFrame);
    
    %Create an evenly sampled x-axis, note that ds varies
    x1 = wwx(1,cur_frame);
    x2 = wwx(end,cur_frame);
    if x1 > x2
        iwwx = x1:-ds(cur_frame):x2;
    else
        iwwx = x1:ds(cur_frame):x2;
    end
    
    iwwy = interp1(wwx(:,cur_frame), wwy(:,cur_frame), iwwx);
    
    temp = fft(iwwy, N_POINTS_FFT);
    
    iY   = abs(temp(1:HALF_N_FFT)); %NOTE: This is magnitude, not power ...
    %This is what the supplemental says, not what was done in the previous 
    %code. I'm not sure what was done for the actual paper, but I would
    %guess they used power.
    %
    %This gets used when determining the secondary wavelength, as it must
    %be greater than half the maximum to be considered a secondary
    %wavelength.
    
    %NOTE: True Amplitude = 2*abs(fft)/(length_real_data i.e. 48 or 49, not 512)
    %
    %i.e. for a sinusoid of a given amplitude, the above formula would give
    %you the amplitude of the sinusoid
    
    %Find peaks that are greater than the cutoff
    [peaks,indx] = seg_worm.util.maxPeaksDist(iY, MIN_DIST_PEAKS,true,WAVELENGTH_PCT_MAX_CUTOFF*max(iY));

    %We sort the peaks so that the largest is at the first index and will
    %be primary, this was not done in the previous version of the code
    [~,I] = sort(-1*peaks); %Sort descending by multiplying by -1
    indx  = indx(I);

    frequency_values = (indx-1)/N_POINTS_FFT*spatial_sampling_frequency(cur_frame);
    
    all_wavelengths = 1./frequency_values;
    
    p_temp = all_wavelengths(1);
    
    if length(indx) > 1
        s_temp = all_wavelengths(2);
    else
        s_temp = NaN;
    end
    
    worm_wavelength_max = WAVELENGTH_PCT_CUTOFF*worm_lengths(cur_frame);
    
    %Cap wavelengths ...
    if p_temp > worm_wavelength_max
        p_temp = worm_wavelength_max;
    end
    
    %??? Do we really want to keep this as well if p_temp == worm_2x?
    %i.e., should the secondary wavelength be valid if the primary is also
    %limited in this way ?????
    if s_temp > worm_wavelength_max
        s_temp = worm_wavelength_max;
    end
    
    p_wavelength(cur_frame) = p_temp;
    s_wavelength(cur_frame) = s_temp;
end


wavelength.primary   = p_wavelength;
wavelength.secondary = s_wavelength;

obj.wavelength = wavelength;
obj.trackLength = trackLength;
obj.amplitude = amplitude;

end