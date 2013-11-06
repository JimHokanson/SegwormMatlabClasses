function old_code

%% Define experiment annotation
%--------------------------------------------------------
%
% Formulate the feature file save data

% The WT2 information.
trackerVersion    = '2.0.4';
hardwareVersion   = '2.0';
analysisVersion   = '2.0';
systemAnnotations = [];
wt2 = struct( ...
    'tracker',      trackerVersion, ...
    'hardware',     hardwareVersion, ...
    'analysis',     analysisVersion, ...
    'annotations',  systemAnnotations);
% The video information.

globalFrameCounter = [];
load(datNameIn, 'globalFrameCounter');
if isempty(globalFrameCounter)
    load(fileInfo.expList.segDat, 'globalFrameCounter');
    if isempty(globalFrameCounter)
        error('features:missingVariable', ['Variable globalFrameCounter ',...
            'could not be found in the segmentation info file or the norm ',...
            'info file. Please re-run segmentation for this experiment.'])
    end
end

% Define frame rate
load(datNameIn, 'myAviInfo');
fps = myAviInfo.fps;

videoFrames = globalFrameCounter;
videoTime = globalFrameCounter / fps;

videoLength = struct( ...
    'frames',   videoFrames, ...
    'time',     videoTime);

height = myAviInfo.height;
width  = myAviInfo.width;

pixel2MicronScale = [];
load(datNameIn, 'pixel2MicronScale');

pixelsPerMicrons = pixel2MicronScale;

fourcc = myAviInfo.fourcc;
micronsPerPixels = struct( ...
    'x', pixelsPerMicrons(1,1), ...
    'y', pixelsPerMicrons(1,2));

resolution = struct( ...
    'fps', fps, ...
    'height', height, ...
    'width', width, ...
    'micronsPerPixels', micronsPerPixels, ...
    'fourcc', fourcc);

% The experiment information.
videoFile    = fileInfo.expList.avi;
vignetteFile = fileInfo.expList.vignette;
infoFile     = fileInfo.expList.xml;
stageLogFile = fileInfo.expList.log;
expDirectory = fileInfo.expList.dir;
expComputer  = char(java.net.InetAddress.getLocalHost.toString);
expFiles = struct( ...
    'video',        videoFile, ...
    'vignette',     vignetteFile, ...
    'info',         infoFile, ...
    'stage',        stageLogFile, ...
    'directory',    expDirectory, ...
    'computer',     expComputer);

% Initialize the info
expTimestamp = [];
expGenotype  = [];
expGene      = [];
expAllele    = [];
expStrain    = [];
expFood      = [];
expVentralSide = [];
expAgarSide  = [];
expSex       = [];
expAge       = [];
expIllumination = [];
expTemperature  = [];
expChemicals    = [];
expArena        = [];
expTracker      = [];
expExperimenter = [];
expLab        = [];
expAddress    = [];
expChromosome = []; 
expWormAnnotations = [];
expEnvironmentAnnotations = [];
expLabAnnotations = [];
expHabitutation   = [];


if handles.preferences.useDB
    experimentId = getExperimentId(conn, fileInfo.expList.fileName);
    
    % Lets connect to the database and fill in the annotation values
    sqlString = strcat('select EA.id, strain.strainName,', {' '}, ...
        'gene.geneName, allele.alleleName, C.chromosomeName,', {' '},...
        'I.datestamp, food.foodName, ventralSide.ventralSideName,', {' '},...
        'WS.sideName, TN.trackerName, sex.sexName, age.ageName,', {' '},...
        'experimenters.name, EL.address, genotype.genotypeName,', {' '},...
        'H.habitName', {' '},...
        'from exp_annotation EA',{' '},...
        'left join gene on EA.geneid = gene.geneid',{' '},...
        'left join strain on EA.strainID = strain.strainID', {' '},...
        'left join allele on EA.alleleid = allele.alleleid', {' '},...
        'left join chromosome C on EA.chromosomeID = C.chromosomeID', {' '},...
        'left join food on EA.foodID = food.foodID', {' '},...
        'left join ventralSide on EA.ventralSideID = ventralSide.ventralSideID', {' '},...
        'left join wormSide WS on EA.agarSideID = WS.sideID', {' '},...
        'left join trackerNO TN on EA.trackerID = TN.trackerID', {' '},...
        'left join sex on EA.sexID = sex.sexID', {' '},...
        'left join age on EA.ageID = age.ageID', {' '},...
        'left join experimenters on EA.experimenterID = experimenters.expID', {' '},...
        'left join experimenterLocation EL on EA.locationID = EL.locationID', {' '},...
        'left join genotype on EA.genotypeID = genotype.genotypeID', {' '},...
        'left join info I on EA.id = I.id', {' '},...
        'left join habituation H on EA.habitId = H.habitId', {' '},...
        'where EA.id =', {' '}, num2str(experimentId),';');
    curs    = exec(conn, sqlString);
    curs    = fetch(curs);
    expInfo = curs.Data;
    close(curs);
    
    if isempty(expInfo)
        warningStr = strcat('MYSQL query failed!', sqlString);
        warning('features:DBconnectForAnnotationFailed', warningStr{1});
    elseif strcmp(expInfo, 'No Data')
        % No info
    else        
        expData.id = expInfo{1,1};
        if strcmpi(expInfo{1,2},'null') || strcmpi(expInfo{1,2},'unknown')
            expData.strainName = [];
        else
            expData.strainName = expInfo{1,2};
        end
        if strcmpi(expInfo{1,3},'null') || strcmpi(expInfo{1,3},'unknown')
            expData.geneName = [];
        else
            expData.geneName = expInfo{1,3};
        end
        if strcmpi(expInfo{1,4},'null') || strcmpi(expInfo{1,4},'unknown')
            expData.alleleName = [];
        else
            expData.alleleName = expInfo{1,4};
        end
        
        if strcmpi(expInfo{1,5},'null') || strcmpi(expInfo{1,5},'unknown')
            expData.chromosomeName = [];
        else
            expData.chromosomeName = expInfo{1,5};
        end
        
        if strcmpi(expInfo{1,6},'null') || strcmpi(expInfo{1,6},'unknown')
            expData.datestamp = [];
        else
            expData.datestamp = expInfo{1,6};
        end
        
        if strcmpi(expInfo{1,7},'null') || strcmpi(expInfo{1,7},'unknown')
            expData.foodName = [];
        else
            expData.foodName = expInfo{1,7};
        end
        
        if strcmpi(expInfo{1,8},'null') || strcmpi(expInfo{1,8},'unknown')
            expData.ventralSide = [];
        else
            expData.ventralSide = expInfo{1,8};
        end
        
        if strcmpi(expInfo{1,9},'null') || strcmpi(expInfo{1,9},'unknown')
            expData.agarSide = [];
        else
            expData.agarSide = expInfo{1,9};
        end
        
        if strcmpi(expInfo{1,10},'null') || strcmpi(expInfo{1,10},'unknown')
            expData.trackerName = [];
        else
            expData.trackerName = expInfo{1,10};
        end
        
        if strcmpi(expInfo{1,11},'null') || strcmpi(expInfo{1,11},'unknown')
            expData.sex = [];
        else
            expData.sex = expInfo{1,11};
        end
        
        if strcmpi(expInfo{1,12},'null') || strcmpi(expInfo{1,12},'unknown')
            expData.age = [];
        else
            expData.age = expInfo{1,12};
        end
        
        
        if strcmpi(expInfo{1,13},'null') || strcmpi(expInfo{1,13},'unknown')
            expData.experimenter = [];
        else
            expData.experimenter = expInfo{1,13};
        end
        
        if strcmpi(expInfo{1,14},'null') || strcmpi(expInfo{1,14},'unknown')
            expData.address = [];
        else
            expData.address = expInfo{1,14};
        end
                
        if strcmpi(expInfo{1,15},'null') || strcmpi(expInfo{1,15},'unknown')
            expData.genotype = [];
        else
            expData.genotype = expInfo{1,15};
        end
        
        if strcmpi(expInfo{1,16},'null') || strcmpi(expInfo{1,16},'unknown')
            expData.habituation = [];
        else
            expData.habituation = expInfo{1,16};
        end
        
        expTimestamp    = expData.datestamp;
        expGene         = expData.geneName;
        expAllele       = expData.alleleName;
        expStrain       = expData.strainName;
        expChromosome   = expData.chromosomeName;
        expFood         = expData.foodName;
        expVentralSide  = expData.ventralSide;
        expAgarSide     =  expData.agarSide;
        expSex          = expData.sex;
        expAge          = expData.age;
        expHabitutation = expData.habituation;
        expTracker      = strrep(expData.trackerName,'tracker_','');
        expExperimenter = expData.experimenter;
        expAddress      = expData.address;
        
        expGenotype = expData.genotype;
        
        expIllumination = '627nm';
        expTemperature  = '22C';
        expChemicals    = [];
        expArena        = 'low-peptone NGM plate';
       
        
        expLab = 'William R Schafer';
        expWormAnnotations        = [];
        expLabAnnotations         = [];
        expEnvironmentAnnotations = [];
    end
else
    % get at least the ventralSide if we are operating offline
    try
        wormNameInfo = parseWormFilename(strcat(fileInfo.expList.fileName, '.avi'));
    catch ME1 %#ok<NASGU>
        wormNameInfo.side = 'unknown';
    end
    
    if strcmpi(wormNameInfo.side, 'L')
        expVentralSide = 'clockwise';
    elseif strcmpi(wormNameInfo.side, 'R')
        expVentralSide = 'anticlockwise';
    else
        expVentralSide = [];
    end
end

% Compute ventralMode
% Set ventralMode as follows -- unknown = 0, anticlockwise = 1, clockwise = 2
ventralMode = 0;
if ~isempty(expVentralSide)
    switch expVentralSide
        case 'anticlockwise'
            ventralMode = 2;
        case 'clockwise'
            ventralMode = 1;
    end
end

expWorm = struct( ...
    'genotype', expGenotype, ...
    'gene', expGene, ...
    'allele', expAllele, ...
    'strain', expStrain, ...
    'chromosome', expChromosome, ...
    'ventralSide', expVentralSide, ...
    'agarSide', expAgarSide, ...
    'sex', expSex, ...
    'age', expAge, ...
    'habituation', expHabitutation, ...
    'annotations', expWormAnnotations);   

expEnvironment = struct( ...
    'timestamp', expTimestamp, ...
    'food', expFood, ...    
    'illumination', expIllumination, ...    
    'temperature', expTemperature, ...
    'chemicals', expChemicals, ...
    'arena', expArena, ...
    'tracker', expTracker, ...
    'annotations', expEnvironmentAnnotations);   

expLab = struct( ...
    'name', expLab, ...
    'experimenter', expExperimenter, ...
    'address', expAddress, ...
    'annotations', expLabAnnotations);

experiment = struct( ...
    'worm', expWorm, ...
    'environment', expEnvironment);