function [intensity_mean, varargout] = intensityMeasure(traj, properties_all, channel_data, varargin)
	%% argument parse
	arg.max_skip_steps = 3;
    arg.background_method = 'meanlog'; % meanlog is the fastest, simple and imageJ 
    arg.numb_brightest_pixels = 10;
    arg.bkgsubstraction = true;
    arg.meantype = 'mean'; % either mean or median
    arg.manualBkg = [];
    arg.averageFac = 20;

	if nargin > 0
		parser = inputParser;
		addParameter(parser, 'max_skip_steps', arg.max_skip_steps);
        addParameter(parser, 'background_method', arg.background_method);
        addParameter(parser, 'numb_brightest_pixels', arg.numb_brightest_pixels);
        addParameter(parser, 'bkgsubstraction', arg.bkgsubstraction);
        addParameter(parser, 'meantype', arg.meantype);
        addParameter(parser, 'manualBkg', arg.manualBkg);
        addParameter(parser, 'averageFac', arg.averageFac);

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
    

	%% get numb of cells
	numb_cells = size(traj,1);

	%% find the ones that are good, jumping step < max_skip_steps
	keep_or_not = ones(1,numb_cells);
	for cn = 1:numb_cells
		tmp_array = [0 isnan(traj(cn,:)) 0 ];
		tmp_array2 = find(diff(tmp_array)==-1) - find(diff(tmp_array)==1);
		if max(tmp_array2) > arg.max_skip_steps
			keep_or_not(cn) = 0;
		end
	end
	cells2keep = find(keep_or_not);
    cell_numb2keep = traj(cells2keep,1); 
	numb_cells2keep = length(cells2keep);

	%% initialize the results matrix
	intensity_sum = zeros(numb_cells2keep, numb_frame);
	intensity_mean = zeros(numb_cells2keep, numb_frame);
    intensity_median = zeros(numb_cells2keep, numb_frame);

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
        
        
		%% get pixel index and areas        
        for cnk = 1:numb_cells2keep
            cell2measure = cells2keep(cnk);
            if ~isnan(traj(cell2measure, fr))
                pixelIdxList = properties_all{fr}(cell2measure).PixelIdxList;
                cell_area = properties_all{fr}(cell2measure).Area;

                % sum up all pixels and average them over the area
                intensity_sum(cnk,fr) = sum(channel_subbkg(pixelIdxList));
                intensity_mean(cnk,fr) = intensity_sum(cnk,fr)/cell_area;
                intensity_median(cnk,fr) = median(channel_subbkg(pixelIdxList));

                % sum up the bright pixels
                sorted_pixels = sort(channel_subbkg(pixelIdxList), 'ascend');
                if length(sorted_pixels)>=arg.numb_brightest_pixels
                    intensity_brightest(cnk,fr) = sum(sorted_pixels(1:arg.numb_brightest_pixels))/arg.numb_brightest_pixels;
                else
                    intensity_brightest(cnk,fr) = sum(sorted_pixels(1:length(sorted_pixels )))/ length(sorted_pixels);
                end
            else
                 % sum up all pixels and average them over the area
                intensity_sum(cnk,fr) = intensity_sum(cnk,fr-1);
                intensity_mean(cnk,fr) = intensity_mean(cnk,fr-1) ;
                intensity_median(cnk,fr) = intensity_median(cnk,fr-1) ;
                intensity_brightest(cnk,fr) = intensity_brightest(cnk,fr-1);
                
            end

        end
    end
    
    center_pos = zeros(2,numb_cells2keep);
    if isfield(properties_all{1}(1), 'Centroid');
        for cnk = 1:numb_cells2keep
            cell2measure = cells2keep(cnk);
            center_pos(:,cnk) = properties_all{1}(cell2measure).Centroid; % only the position in the first frame or last frame if tracked backward
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
            varargout{4} = numb_cells2keep;
        case 5
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
            varargout{3} = background_mean;
            varargout{4} = numb_cells2keep;
            varargout{5} = center_pos;
        case 6
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
            varargout{3} = background_mean;
            varargout{4} = numb_cells2keep;
            varargout{5} = center_pos;
            varargout{6} = intensity_median;
        case 7
            varargout{1} = intensity_sum;
            varargout{2} = intensity_brightest;
            varargout{3} = background_mean;
            varargout{4} = numb_cells2keep;
            varargout{5} = center_pos;
            varargout{6} = intensity_median;
            varargout{7} = cell_numb2keep;

    end
