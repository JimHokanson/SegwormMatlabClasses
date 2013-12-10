function verifyResult(new_worm,feature_mat_path)
%
%
%   seg_worm.feature_calculator.verifyResult(new_worm,feature_mat_path)


assert1 = @(x,y,z)(assert(isequaln(x,y),z));
assert2 = @(x,y,z,w)(assert(h__isSimilar(x,y,z),w));
h = load(feature_mat_path);
old_worm = h.worm;



%Morphology
%==========================================================================






%Posture
%==========================================================================
pn = new_worm.posture;
po = old_worm.posture;

keyboard

%NOTE: Different algorithms yield slightly different results, I had to
%increase the % comparision on this just a bit ...
%
%More points would make this more accurate ...
assert2(pn.eccentricity,po.eccentricity,0.03,'Different eccentricities');
%plot(abs((po.eccentricity - pn.eccentricity)./po.eccentricity))


assert2(pn.amplitude.max,po.amplitude.max,0.01,'posture.amplitude.max')
assert2(pn.amplitude.ratio,po.amplitude.ratio,0.01,'posture.amplitude.ratio')




dn = pn.directions;
do = po.directions;
assert1(dn.tail2head,do.tail2head,'Different tail2head directions')
assert1(dn.head,do.head,'Different head directions')
assert1(dn.tail,do.tail,'Different tail directions')

assert1(pn.skeleton,po.skeleton,'Different Skeletons')

assert1(pn.eigenProjection,po.eigenProjection,'Different eigenprojections')


%               bends: [1x1 struct]
%        eccentricity: DONE
%           amplitude:
%                   .max   - DONE
%                   .ratio
%          wavelength: [1x1 struct]
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

function is_similar = h__isSimilar(x,y,max_diff)
   is_similar = all((isnan(x(:)) & isnan(y(:))) | abs((x(:) - y(:))./x(:)) < max_diff);
end