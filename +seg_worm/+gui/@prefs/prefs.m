classdef prefs < handle
    %
    %   Class:
    %   seg_worm.gui.prefs
    
    properties
        prefData.redoStageMovDiff = 0;
        prefData.redoStageMovDet = 0;
        prefData.redoSeg = 0;
        prefData.deleteVideo = 0;
        prefData.videoOut = 1;
        prefData.redo_ht = 0;
        prefData.normalize = 0;
        prefData.featureSet = 1;

        prefData.standalone = 1;

        prefData.calibrated_ht = 0;
        prefData.disableWarnings_ht = 0;

        prefData.runOnCompleted = 0;

        % prefData.newVignette = 0;

        prefData.bgsubtract = 0;
        prefData.skeldisable = 0;
        prefData.excel = 0;
        prefData.caImaging = 0;

        prefData.useDB = 0;
        prefData.dbupdate = 0;
        prefData.nas = 0;
        prefData.nasOverwrite = 0;
        prefData.experimentCollectionList = 1;
        %prefData.experimentCollectionListName = 'victoriaExperimentList';
        prefData.experimentCollectionListName = 'segmentationExperimentList';

        prefData.version = 1;
    end
    
    methods
    end
    
end

