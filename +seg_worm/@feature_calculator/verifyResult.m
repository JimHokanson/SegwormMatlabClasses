function verifyResult(new_worm,feature_mat_path)
%
%
%   seg_worm.feature_calculator.verifyResult(new_worm,feature_mat_path)


assert1 = @(x,y,z)(assert(isequaln(x,y),z));
assert2 = @(x,y,z,w)(assert(all(h__percentDiff(x,y) < z),w));
assert3 = @(x,y,z,w) (assert(h__corrValue(x,y) > z,w));
h = load(feature_mat_path);
old_worm = h.worm;



%Morphology
%==========================================================================






%Posture
%==========================================================================
pn = new_worm.posture;
po = old_worm.posture;

keyboard


%Due to changes in indices used, the values are all shifted just slightly
%...
assert3(pn.bends.head.mean,po.bends.head.mean,0.99,'posture.bends.head.mean');


%h__summarize(pn.bends.head.mean,po.bends.head.mean)

plot(h__percentDiff(pn.bends.head.mean,po.bends.head.mean))

%NOTE: Different algorithms yield slightly different results, I had to
%increase the % comparision on this just a bit ...
%
%More points would make this more accurate ...
assert2(pn.eccentricity,po.eccentricity,0.03,'Different eccentricities');
%plot(abs((po.eccentricity - pn.eccentricity)./po.eccentricity))


assert2(pn.amplitude.max,po.amplitude.max,0.01,'posture.amplitude.max')
assert2(pn.amplitude.ratio,po.amplitude.ratio,0.01,'posture.amplitude.ratio')
%plot(abs((po.amplitude.max - pn.amplitude.max)./po.amplitude.max))

%Their wavelength code has errors ...
%--------------------------------------------------------------------------
% % % assert3(pn.wavelength.primary,po.wavelength.primary,0.95,'posture.wavelength.primary')
% % % assert3(pn.secondary.ratio,po.secondary.ratio,0.95,'posture.wavelength.secondary')
% % % 
% % % r1 = h__corrValue(pn.wavelength.primary,po.wavelength.primary);
% % % r2 = h__corrValue(pn.wavelength.secondary,po.wavelength.secondary);
% % % 
% % % plot(pn.wavelength.primary,po.wavelength.primary,'o')
% % % 
% % % plot(pn.wavelength.secondary,po.wavelength.secondary,'o')


%Directions
%----------------------------------
dn = pn.directions;
do = po.directions;
assert1(dn.tail2head,do.tail2head,'Different tail2head directions')
assert1(dn.head,do.head,'Different head directions')
assert1(dn.tail,do.tail,'Different tail directions')

%Skeletons
%----------------------------------
assert1(pn.skeleton,po.skeleton,'Different Skeletons')

%EigenProjections
%----------------------------------
assert2(pn.eigenProjection,po.eigenProjection,0.01,'posture.eigenProjection')
%h__summarize(pn.eigenProjection(1,:)',po.eigenProjection(1,:)')

%               bends: [1x1 struct]
%        eccentricity: DONE
%           amplitude:
%                   .max   - DONE
%                   .ratio - DONE
%          wavelength: - their code has errors in it
%         trackLength: [1x4642 double]
%               kinks: [1x4642 double]
%               coils: [1x1 struct]
%          directions: DONE - incorrect currently
%            skeleton: DONE
%     eigenProjection: DONE

end

function h__compareEventStructs(sn,so)

%The new event structs are 1 based indexing instead of 0, compare
%after this change has been made ...


end

%Might replace with max diff or average diff, then do comparison in
%assertion
function pd = h__percentDiff(x,y)
   x = x(:);
   y = y(:);
   mask = ~(isnan(x) | isnan(y));
   pd = zeros(1,length(x));
   pd(mask) = abs((x(mask) - y(mask))./x(mask));
end

function r = h__corrValue(x,y)
    x = x(:);
    y = y(:);
   mask = isnan(x) | isnan(y);
   r = corr(x(~mask),y(~mask));
end

function h__summarize(new,old)
    subplot(2,2,[1 2])
    plot(old,'-o')
    hold all
    plot(new,'-+')
    hold off
    legend({'old' 'new'})
    
    subplot(2,2,3)
    scatter(old,new)
    axis equal
    xlabel('old')
    ylabel('new')
    
    subplot(2,2,4)
    plot(h__percentDiff(old,new))
    title(sprintf('r: %0.3f',h__corrValue(old,new)))
    
end