function [intensity_mean, varargout] = intensityMeasureCyto(traj, nucLabel, properties_allframes, cytoLabel, channel_data, trackDirection, varargin)
% for backward track, channel_data is already flipped. But nucLabel and
% cytoLabel are not
	%% argument parse
	arg.max_skip_steps = 2;
    arg.background_method = 'meanlog'; % meanlog is the fastest, simple and imageJ 
    arg.numb_brightest_pixels = 10;
    arg.bkgsubstraction = true;
    arg.manualBkg = [];
    arg.cells2keep = [];
    arg.averageFac = 10;
    if nargin > 0
		parser = inputParser;
		addParameter(parser, 'max_skip_steps', arg.max_skip_steps);
        addParameter(parser, 'background_method', arg.background_method);
        addParameter(parser, 'numb_brightest_pixels', arg.numb_brightest_pixels);
        addParameter(parser, 'bkgsubstraction', arg.bkgsubstraction);
        addParameter(parser, 'manualBkg', arg.manualBkg);
        addParameter(parser, 'averageFac', arg.averageFac);
        addParameter(parser, 'cells2keep', arg.cells2keep);

		parse(parser,varargin{:});
		arg = parser.Results;
    end 

    
	%% get sizes
	im_size = size(channel_data);
	ndim_data = ndims(channel_data);
    if  ndim_data == 3
		numb_frame = im_size(3);
	elseif ndim_nucdata == 2
		numb_frame = 1;
    end

    
    switch trackDirection
    case 'forward'
        frame_seq = 1:1:numb_frame;
    case 'backward'
        frame_seq = numb_frame:-1:1;
    end  
    
	%% get numb of cells
	numb_cells = size(traj,1);

	%% find the ones that are good, jumping step < max_skip_steps
    if isempty(arg.cells2keep)
        keep_or_not = ones(1,numb_cells);
        for cn = 1:numb_cells
            tmp_array = [0 isnan(traj(cn,:)) 0 ];
            tmp_array2 = find(diff(tmp_array)==-1) - find(diff(tmp_array)==1);
            if max(tmp_array2) > arg.max_skip_steps
                keep_or_not(cn) = 0;
            end
        end
        cells2keep = find(keep_or_not);
    %     cell_numb2keep = traj(cells2keep,1); 
        numb_cells2keep = length(cells2keep);
    else
        cells2keep = arg.cells2keep;
        numb_cells2keep = length(cells2keep);
    end

	%% initialize the results matrix
	intensity_sum = zeros(numb_cells2keep, numb_frame);
	intensity_mean = zeros(numb_cells2keep, numb_frame);
    intensity_median = zeros(numb_cells2keep, numb_frame);
    area_allgoodcells = zeros(numb_cells2keep, numb_frame);

    intensity_brightest = zeros(numb_cells2keep, numb_frame);
    background_mean = zeros(1,numb_frame);
    
    %% start measure frame by frame
    mean_filt_kernel = fspecial('average', 3);
    ball_kernel = strel('ball',arg.averageFac,arg.averageFac,2);
    for fr = 1:numb_frame       
		%% substract the background
        channel_orig = double(channel_data(:,:,fr));
        if arg.bkgsubstraction
            if isempty(arg.manualBkg) % no mauanl bkg
                switch arg.background_method
                    case 'simple'
                        background = imopen(channel_orig,strel('disk',arg.averageFac)); 
                    case 'rollingBall'
                        % this is supposed to be the default preprocessing of ImageJ
                        % not finished yet, super slow
                        max_filtered = ordfilt2(channel_orig,9,true(3));           
                        smoothed = imfilter(max_filtered, mean_filt_kernel);
                        background = imopen(smoothed, ball_kernel);
                    case 'meanlog'
                        background = exp(mean2(log(channel_orig)));
                end
                background_mean(fr) = mean2(background);
            else % manual bkg
                background_mean(fr) = arg.manualBkg(fr);
                background = arg.manualBkg(fr)*ones(size(channel_orig));
            end
            channel_subbkg = channel_orig - background;
            
        else
            background_mean(fr) = 0;
            channel_subbkg = channel_orig;
        end
        
        % get cell numbers in current frame
%         cell_numb2keep_currentFrame = traj(cells2keep,fr); 
        fr2 = frame_seq(fr);
        cyto_label_currentFrame = cytoLabel(:,:,fr2);
        cyto_label_currentFrame_nucremoved = cyto_label_currentFrame;
        cyto_label_currentFrame_nucremoved(logical(nucLabel(:,:,fr2))) = 0;
        
        for cnk = 1:numb_cells2keep          
            pixelIdxList = [];
            % get pixelIdxList of nuc
%             nuc_pixelIdxList = (nucLabel(:,:,fr2) == cell_numb2keep_currentFrame(cnk));       % a logical
            
           
            % get cyto label numbers in these nuc area
%             cyto_label_currentFrame(~nuc_pixelIdxList) = 0;
%             cyto_label_numb = unique(cyto_label_currentFrame(cyto_label_currentFrame~=0));
%             cyto_label_numb(cyto_label_numb==0) = [];
            
            cell2measure = cells2keep(cnk);
            nuc_pixelIdxList = properties_allframes{fr}(cell2measure).PixelIdxList;
            cyto_label_numb = unique(cyto_label_currentFrame(nuc_pixelIdxList));
            cyto_label_numb(cyto_label_numb==0) = [];

            % nuc with no cyto
            if isempty(cyto_label_numb)
                if fr==1
                    intensity_sum(cnk,fr) = 0;
                    intensity_mean(cnk,fr) = 0;
                    intensity_median(cnk,fr) = 0;
                    intensity_brightest(cnk,fr) = 0;
                else
                    intensity_sum(cnk,fr) = intensity_sum(cnk,fr-1);
                    intensity_mean(cnk,fr) = intensity_mean(cnk,fr-1);
                    intensity_median(cnk,fr) = intensity_median(cnk,fr-1);
                    intensity_brightest(cnk,fr) = intensity_brightest(cnk,fr-1);
                end
                continue;
                
            end
            
            % get cyto pixelIdxList
            % pixelIdxList = ismember(cytoLabel(:,:,fr2), cyto_label_numb); % a logical
            % get rid of the nuc region
%             pixelIdxList = pixelIdxList & (~nuc_pixelIdxList);
            
%             pixelIdxList = ismember(cyto_label_currentFrame_nucremoved, cyto_label_numb); % a logical 

            for nl = 1:numel(cyto_label_numb)
                pixelIdxList_tmp = find(cyto_label_currentFrame_nucremoved == cyto_label_numb(nl));
                pixelIdxList = [pixelIdxList pixelIdxList_tmp'];
            end

            % if cyto were removed due to the last line
            if isempty(pixelIdxList)
                if fr==1
                    intensity_sum(cnk,fr) = 0;
                    intensity_mean(cnk,fr) = 0;
                    intensity_median(cnk,fr) = 0;
                    intensity_brightest(cnk,fr) = 0;
                else
                    intensity_sum(cnk,fr) = intensity_sum(cnk,fr-1);
                    intensity_mean(cnk,fr) = intensity_mean(cnk,fr-1);
                    intensity_median(cnk,fr) = intensity_median(cnk,fr-1);
                    intensity_brightest(cnk,fr) = intensity_brightest(cnk,fr-1);
                end    
                continue;
            end
            
            % cell area
%             cell_area = sum(sum(pixelIdxList));
            cell_area = length(pixelIdxList);
            
            % sum up all pixels and average them over the area
            pixel_intensities = channel_subbkg(pixelIdxList);
			intensity_sum(cnk,fr) = sum(pixel_intensities);
			intensity_mean(cnk,fr) = intensity_sum(cnk,fr)/cell_area;
            intensity_median(cnk,fr) = median(pixel_intensities);
            
            % sum up the bright pixels
            sorted_pixels = sort(pixel_intensities, 'ascend');
            if length(sorted_pixels)>=arg.numb_brightest_pixels
                intensity_brightest(cnk,fr) = sum(sorted_pixels(1:arg.numb_brightest_pixels))/arg.numb_brightest_pixels;
            else
                intensity_brightest(cnk,fr) = sum(sorted_pixels(1:length(sorted_pixels )))/ length(sorted_pixels);
            end
            area_allgoodcells(cnk,fr) = cell_area;
        end
    end
    
    
    %%  output

    switch nargout-1
        case 1
            varargout{1} = intensity_sum;
        case 2
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
        case 3
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
            varargout{3} = background_mean;
        case 4
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
            varargout{3} = background_mean;
            varargout{4} = intensity_median;
        case 5
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
            varargout{3} = background_mean;
            varargout{4} = intensity_median;
            varargout{5} = area_allgoodcells;
    end
