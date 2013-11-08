classdef normalized_worm < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.normalized_worm
    %
    %   This will be an interface class between the parsed worms and the
    %   feature sets that are saved to disk. The goal is to take in the
    %   code from normWorms and to have a well described set of properties
    %   for rewriting the feature code.
    %
    %   **** IMPORTANT *********
    %   Since the feature parsing is not quite yet finished, a method
    %   should be written which takes in the normalized worm files and
    %   populates the properties as they will be populated when everything
    %   is implemented.
    %
    %
    %   SegWorm/Pipeline/normWormProcess.m
    %   - ??? Is this different than normWorms?
    %
    %   See Also:
    %   seg_worm.w.util.normWorms
    %
    %
    %
    %    createObjectFromFiles
    %
    %   IMPROVEMENTS
    %   ===================================================================
    %   1) Ensure that non-segmented frames have NaN values ...
    
    
    %{
    %
    %   OLD FORMAT: This is a reference for when I recode things
    %
    %   i.e. references to blockN{5} should be replaced with obj.angles
    %   
    %       blockN{1}  = status:
    %                    s = segmented
    %                    f = segmentation failed
    %                    m = stage movement
    %                    d = dropped frame
    %       blockN{2}  = vulvaContours
    %       blockN{3}  = nonVulvaContours
    %       blockN{4}  = skeletons
    %       blockN{5}  = angles
    %       blockN{6}  = inOutTouches
    %       blockN{7}  = lengths
    %       blockN{8}  = widths
    %       blockN{9}  = headAreas
    %       blockN{10} = tailAreas
    %       blockN{11} = vulvaAreas
    %       blockN{12} = nonVulvaAreas
    %}
    
    %{
    Memory Considerations:
        500 frames - 1 MB
        20 fps - 25 seconds per file
        2.4 MB per minute
        These files need to get REALLY long
        before we need to
    %}
    
    properties (Constant,Hidden)
       %TODO: Eventually move this to a options file ...
       EIGENWORM_PATH = 'F:\worm_data\masterEigenWorms_N2.mat'
    end
    
    properties
        segmentation_status  %[1 n],char
        %    s = segmented
        %    f = segmentation failed
        %    m = stage movement
        %    d = dropped frame
        %    n??? - there is reference in some old code to this which 
        frame_codes         %see comments in seg_worm.parsing.frame_errors
        %near the bottom, I haven't yet coded in the values as constants
        %... :/
        %ex.
        %1    = ok
        %2    = stage movement
        %101  = no worm found
        %
        %   NOTE: This is needed for processing coils ...
        %   
        vulva_contours       %[49 2 n] double
        non_vulva_contours   %[49 2 n] double
        skeletons            %[49 2 n] double
        angles               %[49 n] double, degrees
        %??? Why are angles at edge undefined ??????
        %??? - NaN values ...
        %First and last 5 of 49 values ...
        in_out_touches       %[49 n] double
        lengths              %[1  n] double
        widths               %[49 n] double
        head_areas           %[1  n] double
        tail_areas           %[1  n] double
        vulva_areas          %[1  n] double
        non_vulva_areas      %[1  n] double
    end
    
    properties
       eigen_worms 
    end
    
    properties (Dependent)
       x
       y
       contour_x
       contour_y
    end
    
    methods
        function value = get.contour_x(obj)
           value = squeeze([obj.vulva_contours(:,1,:); obj.non_vulva_contours(end-1:-1:2,1,:);]); 
        end
        function value = get.contour_y(obj)
           value = squeeze([obj.vulva_contours(:,2,:); obj.non_vulva_contours(end-1:-1:2,2,:);]);
        end
        function value = get.x(obj)
           value = squeeze(obj.skeletons(:,1,:)); 
        end
        function value = get.y(obj)
           value = squeeze(obj.skeletons(:,2,:));
        end
        function value = get.eigen_worms(obj)
           if isempty(obj.eigen_worms)
              h = load(obj.EIGENWORM_PATH);
              obj.eigen_worms = h.eigenWorms;
           end
           value = obj.eigen_worms;
        end
    end
    
    methods (Static)
        function obj = getObject(norm_folder)
            %
            %   seg_worm.normalized_worm.getObject(norm_folder)
            %   
            
           file_path = fullfile(norm_folder,'norm_obj.mat');
           
           if ~exist(file_path,'file')
               %NOTE: Assumption being made here as to format ...
              norm_partial = fullfile(norm_folder,'normBlock1.mat');
              if exist(norm_partial,'file')
                  seg_worm.normalized_worm.createObjectFromFiles(norm_folder)
              else
                  error('Normalized worm object and partial objects not found')
              end
           end
            
           h = load(file_path);
           
           obj  = seg_worm.normalized_worm;
           
           sl.struct.toObject(obj,h.s);
           
        end
        function createObjectFromFiles(norm_folder,failed_frame_path)
            %
            %
            %    seg_worm.normalized_worm.createObjectFromFiles(norm_folder)
            
            
            %Get files and order
            %--------------------------------------------------------------
            x = dir(fullfile(norm_folder,'normBlock*'));
            
            names = {x.name};
            
            numbers = str2double(regexp(names,'\d+','once','match'));
            
            [num_sorted,I] = sort(numbers);
            
            n_files = length(x);
            
            if ~isequal(num_sorted,1:n_files)
                error('Coding assumption violated, files don''t go from 1 to n')
            end
            
            %Load each file and concatenate
            %--------------------------------------------------------------
            obj = seg_worm.normalized_worm;
            
            for iFile = 1:n_files
                cur_file_name = names{I(iFile)};
                file_path = fullfile(norm_folder,cur_file_name);
                
                h = load(file_path);
                % :/ The variables are named normBlock1
                
                %We'll just extract the only field, normBlock1, normBlock2, etc
                fn = fieldnames(h);
                file_values = h.(fn{1});
                
                %NOTE: These are aligned to the order
                %in the files ...
                FIELDS = {'segmentation_status'
                        'vulva_contours'        %3rd dimension concatenate
                        'non_vulva_contours'    %3
                        'skeletons'             %3
                        'angles'           
                        'in_out_touches'      
                        'lengths'
                        'widths' 
                        'head_areas'
                        'tail_areas' 
                        'vulva_areas'
                        'non_vulva_areas'};
                
                
                n_fields = length(FIELDS);
                
                for iField = 1:n_fields
                   cur_field = FIELDS{iField};
                   if iField >= 2 && iField <= 4
                       cat_dim = 3;
                   else
                       cat_dim = 2;
                   end
                       obj.(cur_field) = cat(cat_dim,obj.(cur_field),file_values{iField});
                end
            end
            
            % :/ The failed frames are saved to disk ...
            h = load(failed_frame_path);
            failed_frames = h.failedFrames; %Format:

            obj.frame_codes = seg_worm.parsing.frame_errors.segmentationStatusToCodes(...
                obj.segmentation_status,failed_frames);
            
            s = sl.obj.toStruct(obj); %#ok<NASGU>
        
            save(fullfile(norm_folder,'norm_obj.mat'),'s');
            
        end
    end
    
end

