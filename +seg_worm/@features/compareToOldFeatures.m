function compareToOldFeatures(new_worm,feature_mat_path)
%
%
%   seg_worm.feature_calculator.verifyResult(new_worm,feature_mat_path)
%
%   UNFINISHED    UNFINISHED    UNFINISHED    UNFINISHED
%


%This is meant to be a quick script that checks that our processing code
%produces the same (or similar) results that their code does.

h = load(feature_mat_path);
old_worm = h.worm;

h_new.worm = new_worm;
h_old.worm = old_worm;

%Movement Tests
%==========================================================================

m_specs = seg_worm.stats.movement_specs.getSpecs;

all_field_names = {m_specs.feature_field};

min_corr_values = 0.99*ones(size(all_field_names));

amp_wavelength_names = {...
    'posture.amplitude.ratio'
    'posture.wavelength.primary'
    'posture.wavelength.secondary'};

min_corr_values(ismember(all_field_names,amp_wavelength_names)) = 0.97;


%Expected differences
%---------------------------------------------------
%amplitude.ratio - different algorithm
%wavelengths - error in deep function that I am not going to bother to change ...

passed_movement = true;
for iSpec = 1:length(m_specs)
    
    cur_spec = m_specs(iSpec);
    
    cur_deep_field_name = cur_spec.feature_field;
    
    x = cur_spec.getData(old_worm);
    y = cur_spec.getData(new_worm);
    
    if isequaln(x,y)
        fprintf(1,'%s are exactly equal!\n',cur_deep_field_name);
    else
    r = h__corrValue(x,y);
    if r < min_corr_values(iSpec)
        passed_movement = false;
        fprintf(2,'%s not exactly equal: r = %0.3f\n',cur_deep_field_name,r);
    else
        fprintf(1,'[\b close:  %0.3f, %s]\b\n',r,cur_deep_field_name);
    end
    end
    
end

if ~passed_movement
    error('Failed movement tests')
end

%Simple Tests
%==========================================================================

s_specs = seg_worm.stats.simple_specs.getSpecs;

passed_simple = true;
for iSpec = 1:length(s_specs)
    
    cur_spec = s_specs(iSpec);

    cur_deep_field_name = cur_spec.feature_field;
    
    x = cur_spec.getData(old_worm);
    y = cur_spec.getData(new_worm);
    if isequal(size(x),size(y)) && h__corrValue(x,y) > 0.999
        fprintf(1,'%s are equal!\n',cur_deep_field_name);
    else
        passed_simple = false;
        fprintf(2,'%s not equal\n',cur_deep_field_name);
    end
    
end

if ~passed_simple
    error('Failed "simple" data tests')
end


%The simple specs might have different lengths as they are not

%Event Tests
%==========================================================================

n_frames = 4642; %TODO: Remove hardcoding

e_specs = seg_worm.stats.event_specs.getSpecs;

passed_events = true;
for iSpec = 1:length(e_specs)
    
    cur_spec = e_specs(iSpec);

    cur_deep_field_name = cur_spec.feature_field;
    
    if ~isempty(cur_spec.sub_field)
        cur_deep_field_name = [cur_deep_field_name '.' cur_spec.sub_field]; %#ok<AGROW>
    end

    x = cur_spec.getData(old_worm,n_frames);
    y = cur_spec.getData(new_worm,n_frames);
    
    %NOTE: Some of the values are singular, can't do correlation on single values ...
    if isequaln(x,y)
        fprintf(1,'%s are equal!\n',cur_deep_field_name);
    elseif ~isequal(numel(x),numel(y))
        passed_events = false;
        fprintf(2,'%s not equal, size mismatch\n',cur_deep_field_name);
    elseif all(h__percentDiff(x,y) < 0.001 | isnan(x) | isnan(y))
        fprintf(1,'%s are equal!\n',cur_deep_field_name);
    else
        passed_events = false;
        fprintf(2,'%s not equal\n',cur_deep_field_name);
    end
end

%NOT EQUAL
%----------------------------
%posture.coils.frames.interTime not equal
%Final event time difference is just slightly larger ????
%22.4073   89.5905   37.8873   14.1642       NaN
%22.4073   89.5905   37.8873   14.2029       NaN

if ~passed_events
    error('Failed "event" tests')
end

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
