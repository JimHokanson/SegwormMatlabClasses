function [amplitude,wavelength,trackLength] = getAmplitudeAndWavelength(theta_d,xx,yy,worm_lengths)
%
%   
%   TODO: Finish documentation
%
%
%   [amplitude,wavelength,trackLength] = ...
%       seg_worm.feature_helpers.posture.getAmplitudeAndWavelength(theta_d,xx,yy,wormLengths)
%
%   Inputs
%   =======================================================================
%   theta_d     : 
%   xx          :
%   yy          :
%   wormLengths :
%
%
%   Outputs
%   =======================================================================
%   amplitude    :
%       .max       -
%       .ratio     -
%   wavelength   :
%       .primary   -
%       .secondary - this might not always be valid, even when the primary
%                    wavelength is defined
%   trackLength  :
%
%   IMPROVEMENTS:
%   ----------------------------------------------
%   1) Add error check (see code)
%   2) Finish documentation ...
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
%frequency units, not in real frequency units


%TODO: Add check that N_POINTS_FFT is not less than the # of samples ...

theta_r = theta_d*pi/180;

% xx = bsxfun(@minus,xx,mean(xx,1));
% yy = bsxfun(@minus,yy,mean(yy,1));

%Unrotate worm
%-----------------------------------------------------------------
wwx = bsxfun(@times,xx,cos(theta_r)) + bsxfun(@times,yy,sin(theta_r));
wwy = bsxfun(@times,xx,-sin(theta_r)) + bsxfun(@times,yy,cos(theta_r));

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
%   

wwx = bsxfun(@minus,wwx,mean(wwx,1));
wwy = bsxfun(@minus,wwy,mean(wwy,1));

% Calculate track amplitude
amp1 = max(wwy);
amp2 = abs(min(wwy));
amplitude.max   = amp1 - amp2;
amplitude.ratio = min(amp1,amp2)./max(amp1,amp2);

% Calculate track length
trackLength = max(wwx)-min(wwx);

%Wavelength Calculation
%--------------------------------------------------------------------------

dwwx = diff(wwx,1,1);

%Does the sign change? This is a check to make sure that the change in x is
%always going one way or the other
badWormOrient = any(bsxfun(@ne,sign(dwwx),sign(dwwx(1,:))),1);

n_frames = length(badWormOrient);

p_wavelength = NaN(1,n_frames);
s_wavelength = NaN(1,n_frames);

%Normalize worms, hold onto length ...

minx = min(wwx,[],1);
maxx = max(wwx,[],1);
lengthx = maxx - minx;

% norm_wwx = bsxfun(@rdivide,bsxfun(@minus,wwx,minx),lengthx);
% 
% n_samples = size(norm_wwx,1);
% iwwx = repmat(linspace(0,1,n_samples)',1,n_frames); %transpose to column vector
% 
% iwwy = interp1(norm_wwx,wwy,iwwx);

%NOTE: Right now this varies from worm to worm which means the spectral
%resolution varies as well from worm to worm
spatial_sampling_frequency = (size(wwx,1)-1)./lengthx;

ds = 1./spatial_sampling_frequency;

frames_to_calculate = find(~badWormOrient);

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
    
    iY  = abs(temp(1:HALF_N_FFT)); %NOTE: This is magnitude, not power ...
    %This is what the supplemental says, not what was done in code ...
    %Not sure what was done in the reference paper
    %
    %NOTE: Amplitude = 2*abs(fft)/(length_real_data i.e. 48 or 49)
    
    [peaks,indx] = seg_worm.util.maxPeaksDist(iY, MIN_DIST_PEAKS,true,0.5*max(iY));

    [~,I] = sort(-1*peaks); %Sort descending, don't like dim option ...
    indx      = indx(I);

    f = (indx-1)/512*spatial_sampling_frequency(cur_frame);
    
    all_wavelengths = 1./f;
    
    p_temp = all_wavelengths(1);
    
    if length(indx) > 1
        s_temp = all_wavelengths(2);
    else
        s_temp = NaN;
    end
    
    worm_2x = 2*worm_lengths(cur_frame);
    
    %Cap wavelengths ...
    if p_temp > worm_2x
        p_temp = worm_2x;
    end
    
    %??? Do we really want to keep this as well if p_temp == worm_2x?
    %i.e., should the secondary wavelength be valid
    if s_temp > worm_2x
        s_temp = worm_2x;
    end
    
    p_wavelength(cur_frame) = p_temp;
    s_wavelength(cur_frame) = s_temp;
end


wavelength.primary   = p_wavelength;
wavelength.secondary = s_wavelength;

% plot(p_wavelength,'bo')
% hold on
% plot(p_wavelength2,'k+')
% plot(s_wavelength,'go')
% plot(s_wavelength2,'r+')
% hold off
% keyboard


end


function [p_wavelength,s_wavelength] = h__getWavelengths(badWormOrientAll,wwxa,wwya,wormLengths)
%
%
%   This is the old code 
%

n_frames = length(badWormOrientAll);
nintervals = size(wwxa,1)-1;

p_wavelength = NaN(1,n_frames);
s_wavelength = NaN(1,n_frames);

for iFrame = 1:n_frames

   wwx = wwxa(:,iFrame);
   wwy = wwya(:,iFrame);
   badWormOrient = badWormOrientAll(iFrame);
   wormLen = wormLengths(iFrame);
    
% Wavelength
if (badWormOrient>0)
    wavelengths = [nan, nan];
else
    % for non-curled worms, can measure wavelength
    % Interpolate signal to position vertices equally
    %   distributed along X-axis...
    % (NOTE: Using signal as a factor of X-position, NOT TIME!
    %   Hence, references of "frequency" are "spatial frequency,"
    %   not temporal frequency)
    
    % Reference vector of equally distributed X-positions
    iwwx = [wwx(1) : (wwx(end)-wwx(1))/nintervals : wwx(end)];
    % Signal interpolated to reference vector
    try
        iwwy = interp1(wwx, wwy, iwwx);
    catch ME1
        msgString = getReport(ME1, 'extended','hyperlinks','off');
        msgbox(msgString);
    end
    % Calculate spatial frequency of signal (cycles/PIXEL)
    % TJ: 512 - number of sampling points
    iY = fft(iwwy, 512);
    % Power of constitutive "frequency" components (From Matlab online
    % documentation for _fft [1]_)
    iPyy = iY.* conj(iY) / 512;
    % Worm length (pixels) of the bounding box
    
    xlength = max(iwwx) - min(iwwx);
    
%     if iFrame == 322
%         keyboard
%     end
    
%     keyboard
    
%     spatial_frequency = 1./(iwwx(2) - iwwx(1));
%     
%     f4 = linspace(0,spatial_frequency,512);
%     
%     
%     f3 = spatial_frquency*(0:256/512);
%     
% 
%     f3 = 0:spatial_frequency/512:511/512spatial_frequency;
    
    % Vector of "frequencies" (cycles/pixel), from 0 (steady state factor)
    % to Nyquist frequency (i.e. 0.5*sampling frequency).
    % (In our case sampling frequency is typically 12 samples per worm).
    % (Nyquist frequency is the theoretical highest frequency that can be
    % accurately detected for a given sampling frequency.)
    f = (nintervals/xlength) * (0:256)/512;
    
    % Define search distance
    distance = 5; % round(257/48); -- we can try 10 if this gives poor
    % Peak discrimination
    % MaxPeaksDist is a function written by Ev Yemini to find peaks
    [peaks, indx] = seg_worm.util.maxPeaksDist(iPyy(1:257), distance,true,0);
    
    % We will filter the peaks that are smaller than 1/2 of the maximum
    % peak value
    indxFilt = indx(peaks>max(peaks)/2);
    
    wavelnthArray = 1./f(indxFilt);
    
    wavelnth2 = NaN;
    
    wavelnth1 = wavelnthArray(1);
    if length(indxFilt) > 1
        wavelnth2 = wavelnthArray(2);
    end
    % We will cap the value of wavelength to be at most 2x the length of
    % the worm
    if wavelnth1 > 2*wormLen
        wavelnth1 = 2*wormLen;
    end
    if wavelnth2 > wormLen
        wavelnth2 = 2*wormLen;
    end
	% Save wavelength
    wavelengths = [wavelnth1, wavelnth2];
    
end

p_wavelength(iFrame) = wavelengths(1);
s_wavelength(iFrame) = wavelengths(2);

end





end