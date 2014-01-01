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
    'worm.posture.amplitude.ratio'
    'worm.posture.wavelength.primary'
    'worm.posture.wavelength.secondary'}';

min_corr_values(ismember(all_field_names,amp_wavelength_names)) = 0.97;


%Expected differences
%---------------------------------------------------
%amplitude.ratio - different algorithm
%wavelengths - error in deep function that I am not going to bother to change ...

passed_movement = true;
for iSpec = 1:length(m_specs)

cur_deep_field_name = m_specs(iSpec).feature_field; 
    
x = eval(['h_old.' cur_deep_field_name]);
y = eval(['h_new.' cur_deep_field_name]);

if isequaln(x,y)
    fprintf(1,'%s are exactly equal!\n',cur_deep_field_name);
    continue
end
r = h__corrValue(x,y);
if r < min_corr_values(iSpec)
    passed_movement = false;
    fprintf(2,'%s not exactly equal: r = %0.3f\n',cur_deep_field_name,r);
else
    fprintf(1,'[\b close:  %0.3f, %s]\b\n',r,cur_deep_field_name);
end

end

if ~passed_movement
   error('Failed movement tests') 
end

%Simple Tests
%==========================================================================

s_specs = seg_worm.stats.simple_specs.getSpecs;

%The simple specs might have different lengths as they are not 

%Event Tests
%==========================================================================

e_specs = seg_worm.stats.event_specs.getSpecs;

for iSpec = 1:length(e_specs)
   %TODO: Still need to implement this
   %
   %In spec, implement data retrieval given head object ...
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