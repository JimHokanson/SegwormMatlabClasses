%{ Finding_Stage_Movement_Algorithm

% The algorithm is as follows:
%
% 1. video2Diff differentiates a video frame by frame and outputs the
% differential variance. We load these frame differences.
%
% 2. The info file contains the tracking delay. This delay represents the
% minimum time between stage movements and, conversely, the maximum time it
% takes for a stage movement to complete. If the delay is too small, the
% stage movements become chaotic. We load the value for the delay.
%
% 3. The log file contains the initial stage location at media time 0 as
% well as the subsequent media times and locations per stage movement. Our
% algorithm attempts to match the frame differences in the video (see step
% 1) to the media times in this log file. Therefore, we load these media
% times and stage locations.
% 
% 4. If there are no stage movements, we're done.
%
% 5. The log file sometimes contains more than one entry at 0 media time.
% These represent stage movements that may have completed before the video
% begins. Therefore, we don't use them but, instead, store them in case we
% find their stage movement in the video frame differences.
%
% 6. The frame differences need to be aligned to the video frames.
% Therefore, we copy the first one and shift the rest over so that the
% frame differences are one-to-one with the video frames. Note that video
% indexing starts at 0 while Matlab indexing begins at 1. Also note, due
% to differentiation, large frame differences that occur at the beginning
% of a stage movement now represent the first frame of movement. Large
% frame differences that occur at the end of a stage movement now represent
% the first non-movement frame.
%
% 7. Compute the global Otsu threshold for the frame-differences to
% distinguish stage-movement peaks. Then compute a global non-movement
% threshold by taking all values less than the Otsu, and computing 3
% standard deviations from the median (approximately 99% of the small
% values). Please note, stage movements ramp up/down to
% accelerate/decelerate to/from the target speed. Therefore, the values
% below the Otsu threshold are contaminated with stage movement
% acceleration and decelaration. Fortunately, non-movement frames account
% for more than 50% of the frame differences. Therefore, to avoid the stage
% movement contamination, we use the median of the small values rather than
% their mean when computing the global non-movement (small) threshold. If
% this small threshold is greater than the Otsu, we've made a poor choice
% and ignore both thresholds. Otherwise, these 2 global thresholds serve as
% a backup to the local ones when distinguishing stage movements.
%
% 8. Ocasionally, computing the global Otsu threshold fails. This occurs
% when a long video has few stage movements. In this case, stage movements
% appear to be rare outliers and the Otsu method minimizes in-group
% variance by splitting the non-stage movement modality into two groups
% (most likely periods of worm activity and inactivity). Under these
% circumstances we attempt to use a global threshold at half the maximum
% frame-difference variance. As detailed above, we test this threshold to
% see whether it is sufficiently larger than 99% of the smaller movements.
%
% 9. When searching for stage movements, we use the same algorithm as the
% one above(see step 7), over a smaller well-defined, search window, to
% determine the local Otsu threshold. The local small threshold is computed
% slightly differently (we use the mean instead of the median -- see step
% 12 for an explanation). If the local Otsu threshold fails (it's smaller
% than 99% of the values below it and smaller than the global Otsu), we
% attempt to use the global one to see if it pulls out a peak.
%
% 10. A stage movement peak is defined as the largest value that exceeds
% the minimum of the global and local Otsu thresholds. To avoid a situation
% in which 2 stage movements occur within the same search window, we scan
% the window for the first value exceeding the Otsu threshold and, if any
% subsequent values drop below the minimum of global and local small
% thresholds, we cut off the remainder of the window and ignore it.
%
% 11. Once a stage-movement peak is identified, we search for a temporary
% back and front end to the movement. The stage movement must complete
% within one delay time window (see step 2). Therefore, we search for the
% minimum values, within one delay time window, before and after the peak.
% The locations of these minimum values are the temporary back and front
% ends for the stage movement. If the either value is below the small
% threshold, we may have overshot the movement and, therefore, attempt to
% find a location closer to the peak. Similarly, if either value is greater
% than the maximum of the global and local small thresholds and the
% remaining values till the end of the window are NaNs or, if either value
% is greater than the Otsu threshold, we may have undershot the movement
% and, therefore, attempt to find a location further from the peak.
%
% 12. Using one stage movement's temporary front end and the subsequent
% movement's temporary back end, we compute the small threshold. This
% interval is assumed to have no stage motion and, therefore, represents
% frame-difference variance from a non-movement interval. The local small
% threshold is defined as 3 deviations from the mean of this non-movement
% interval (99% confidence). With this small threshold, we start at the
% first peak and search forward for its real front end. Similarly, we start
% at the subsequent peak and search backward for its real back end.
%
% 13. Conservatively, the beginning and end of the video are treated as the
% end and begining of stage movements, respectively. We search for a
% temporary front end and a temporary back end, respectively, using the
% global small and Otsu thresholds.
%
% 14. After the final, logged stage motion is found in the frame
% differences, we look to see if there are any subsequent, extra peaks.
% An Otsu threshold is computed, as detailed earlier, using the interval
% spanning from the final stage movement's temporary front end till the
% final frame difference. If this threshold is unsuitable, we use the
% global Otsu threshold. If we find any extra peaks, the first peak's back
% end is located and its frame as well as the remainder of the frame
% differences (including any other extra peaks) are marked as a single
% large stage movement that extends till the end of the video. This is
% necessary since Worm Tracker users can terminate logging prior to
% terminating the video (e.g., this may occur automatically if the worm is
% lost).
%
% 15. To find a stage movement, we compute its offset media time. The first
% stage movement is offset by the delay time (see step 2). Subsequent media
% times are offset by the difference between the previous media time and
% its stage-movement peak. Therefore, each stage movement provides the
% offset for the next one. The longer the wait till the next stage
% movement, the less accurate this offset becomes. Therefore, we search for
% each stage movement using a window that begins at the last stage
% movement's temporary front end and ends at the offset media time plus
% this distance (a window with its center at the offset media time). If the
% window is too small (the offset media time is too close to the temporary
% front end of the previous stage movement), we extend its end to be the
% offset media time plus the delay time.
%
% 16. If a stage movement peak is absent, we attempt to shift the media
% times backward or forward, relative to the stage movement peaks,
% depending on whether the current peak is closer to the next or previous
% media time, respectively. When shifting the media times backward, we
% assume the first stage movement occurred prior to video recording and,
% therefore, throw away its media time and location. When shifting the
% media times forward, we look for a spare 0 media time and location (see
% step 5). If there are no spares, we swallow up all the frames prior to
% the end of the first stage movement and label them as an unusable
% movement that began prior to recording and bled into the beginning of the
% video.
% 
% 17. If we find a stage-movement peak closer to the previous offset media
% time than its own supposed offset media time, we assume we have a
% misalignment and attempt to shift the media times forward relative to the
% stage movement peaks. There are some restrictions to this check since
% small-scale, frame jitter can misrepresent the reported media time.
%
% 18. The final logged stage motion may occur after the video ends and, as
% a result, may have no representation in the frame-difference variance.
% Therefore, for the last stage movement, we check our threshold (see step
% 10) to ensure that it separates 99% of the smaller values and, thereby,
% picks up stage movement rather than splitting the non-movement modality.


%}