% do not change this file without Yanfei's permission, everyone else is using it
% if you really need to change it, make a copy of your own, and call it


%% version 2, parallel computation for different positions

segmentFuncHandle = @segNucWaterShed;
trackingFuncHand = @trackNucNearest;
measureFuncHandle = @intensityMeasure;
measureCytoFuncHandle = @intensityMeasureCyto;
segmentPhaseFuncHandle = @segmentPhaseYeast;

max_change_factor = 5; % for quality control, 5 means 5 sigmas

%% contruct channels struct to store information of channels
channels = struct('channelName', usedChName, 'channelNum', usedChNumb);

%% get number of the all channels
channels2measure_num = zeros(1,length(channels2measure));
nuc_ch_num = [];
for ch = 1:length(channels)
    
    % get nuc channel number
    if strcmp(nuc_ch_name, channels(ch).channelName)
		nuc_ch_num = channels(ch).channelNum;

    elseif strcmp(phase_ch_name, channels(ch).channelName)
        phase_ch_num = channels(ch).channelNum;
    else
    end
    
    % get all channels need to measure
    for ch2m = 1:length(channels2measure)
        if strcmp(channels2measure(ch2m), channels(ch).channelName)
            channels2measure_num(ch2m) = channels(ch).channelNum;
        end
    end
end
assert(~isempty(nuc_ch_num), 'wrong channel name for the nuclei');

% check if just want to do one file test
if oneFileTest
    numb_folders = 1;
else
    numb_folders = numel(Folders);
end


%% check folders names and file names, avoid get into long run before detecting file name erros
for fd = 1:numel(Folders)

	folderName = Folders{fd};
	%% file names untangle
	baseName = BaseFileNameExample{fd};
	nucfiles = dir([fullfile(folderName, baseName) '*c' num2str(nuc_ch_num) '.tif']);
    assert(~isempty(nucfiles), 'looks like some typos in the file name, or folder name');
    
end

% for all folders
for fd = 1:numel(Folders)

	folderName = Folders{fd};
    display(['working on folder ' folderName]);
	%% file names untangle
	baseName = BaseFileNameExample{fd};
	PosNameStartIdx = length(baseName) + 1;
	nucfiles_test = dir([fullfile(folderName, baseName) '*c' num2str(nuc_ch_num) '.tif']);
    if isstruct(nucfiles_test)
        nucfiles = dir([fullfile(folderName, baseName) '*c' num2str(nuc_ch_num) '.tif']);
        phasefiles = dir([fullfile(folderName, baseName) '*c' num2str(phase_ch_num) '.tif']);

    else
        nucfiles.name = dir([fullfile(folderName, baseName) '*c' num2str(nuc_ch_num) '.tif']);
        phasefiles.name = dir([fullfile(folderName, baseName) '*c' num2str(phase_ch_num) '.tif']);
    end
           
        
    
	fileNameEx = nucfiles(1).name;
	CIndex = find(fileNameEx == 'c',1, 'last');
    
    %% check if just want to do one file
    if oneFileTest
        numb_files = 1;
    else
        numb_files = numel(nucfiles);
    end
    
    % initilized some results var
    position_mean = cell(numb_files,1);
    position_median = cell(numb_files,1);
    position_sum = cell(numb_files,1);
    position_bkg = cell(numb_files,1);
    position_bright = cell(numb_files,1);
    
    position_mean_cyto = cell(numb_files,1);
    position_median_cyto = cell(numb_files,1);
    position_sum_cyto = cell(numb_files,1);
    position_bright_cyto = cell(numb_files,1);
    position_bkg_cyto = cell(numb_files,1);

    
    %% get nuc file name, numb of frames, image width and height
    fl = 1; %only one position is needed
    nucfileName = fullfile(folderName,nucfiles(fl).name);

    imageInformation =  imfinfo(nucfileName);
    numb_frame = numel(imageInformation);
    imwidth = imageInformation(1).Width;
    imheight = imageInformation(1).Height;
    
    %%
    start_frame = 28;
    end_frame = numb_frame;
    numb_frame2 = numb_frame-start_frame+1;
    
    %% get phase file name, numb of frames

    switch trackDirection
        case 'forward'
            frame_seq = start_frame:1:end_frame;
        case 'backward'
            frame_seq = end_frame:-1:start_frame;
    end        
    
    %% start to work
    nucLabel_allPos = cell(1,numb_files);
    reg_output_allPos = cell(1,numb_files);
    cytoLabel_allPos = cell(1,numb_files);
    traj_allPos = cell(1,numb_files);
    move_steps_allPos = cell(1,numb_files);
    properties_allframes_allPos = cell(1,numb_files);
    properties_allframes_allPos_cyto = cell(1,numb_files);
    
    %% for save how many cells at each position
    numb_cells2keep = zeros(1,numb_files);
    cells_centers = cell(1,numb_files);
    cells_numb2keep = cell(1, numb_files);
    cyto_area = cell(1,numb_files);
    
    parfor fl = 1:numb_files % for all positions, change to parfor if want parallel
        if fl == 1;
            profile on;
        end
        fprintf('working on repeat number %d in current folder\n', fl);
        %% get nuc file name, numb of frames, image width and height
		nucfileName = fullfile(folderName,nucfiles(fl).name);
        phasefileName = fullfile(folderName,phasefiles(fl).name);

		positionName = nucfiles(fl).name(PosNameStartIdx:CIndex-1);
        imageInformation =  imfinfo(nucfileName);
                
        %% read nuclei  
        nuc_data = zeros(imheight, imwidth, numb_frame2, 'uint16');        
		for fr = 1:numb_frame2
			nuc_data(:,:,fr) = imread(nucfileName,start_frame+fr-1);
		end

		%% segmentation and save, comment these two line if you have already done segmentation
		[nucLabel, reg_output] = segmentFuncHandle(nuc_data, 'registerOrNot', registerOrNot);        
        reg_output_allPos{fl} = reg_output;
        fprintf('finish nuc segmentation \n');
        
        %% cytoplasm segmentation, nuc are relabeled, the ones that does not have cyto are deleted
        phase_data = zeros(imheight, imwidth, numb_frame2, 'uint16');        
		for fr = 1:numb_frame2
            phase_data(:,:,fr) = imread(phasefileName,start_frame+fr-1);
		end
        [cytoLabel, nucLabel] = segmentPhaseFuncHandle(phase_data, nucLabel);
        nucLabel_allPos{fl} = nucLabel;
        cytoLabel_allPos{fl} = cytoLabel;
        fprintf('finish cyto segmentation \n');
        

        %% tracking 
        [traj, properties_allframes,move_steps] = trackingFuncHand(nucLabel, reg_output, 'trackDirection', trackDirection, 'max_move_distance', 7);
        traj_allPos{fl} = traj;
        move_steps_allPos{fl} = move_steps;
        properties_allframes_allPos{fl} = properties_allframes;
%         properties_allframes_allPos_cyto{fl} = properties_cyto_allframes;
        fprintf('finish tracking \n');
        
		%% intensity measure                  
        channel_file_name = [fullfile(folderName, baseName) positionName 'c'];
        channel_data = zeros(imheight, imwidth, numb_frame2, 'uint16');
        
        %% kind of weired, avoid broadcasting these varibles
        channels2measure_num;
        channels;
        frame_seq;
        
        %% preallocate 
        channel_mean = cell(length(channels2measure),1);
        channel_median = cell(length(channels2measure),1);
        channel_sum = cell(length(channels2measure),1);
        channel_brightest = cell(length(channels2measure),1);
        chbkg_mean = cell(length(channels2measure),1);
        
        channel_mean_cyto = cell(length(channels2measure),1);
        channel_median_cyto = cell(length(channels2measure),1);
        channel_sum_cyto = cell(length(channels2measure),1);
        channel_brightest_cyto = cell(length(channels2measure),1);
        chbkg_mean_cyto = cell(length(channels2measure),1);

        %% meausre for different channels
        for ch = 1:length(channels2measure) 
            ch2measure = channels2measure_num(ch);
            fileName = [channel_file_name num2str(channels(ch2measure).channelNum) '.tif'];
            
            % read the image, depend on track direction
            for fr = 1:numb_frame2
                fr2 = frame_seq(fr);
                channel_data(:,:,fr) = imread(fileName, fr2);
            end            
            
            % calculate bkground
%             if strcmp(trackDirection, 'backward')
%                 background = getBackground(channel_data, 'avg_master_mask','masterFrameNumb', 1);
%             else
%                 background = getBackground(channel_data, 'avg_master_mask','masterFrameNumb', numb_frame);
%             end
            
            % measure
            if ch2measure == nuc_ch_num % if measure nuc channel, get some more information
                % nuc area
                [channel_mean{ch}, channel_sum{ch}, channel_brightest{ch},chbkg_mean{ch}, numb_cells2keep(fl), cells_centers{fl}, channel_median{ch}, cells_numb2keep{fl}]= ...
                    measureFuncHandle(traj, properties_allframes, channel_data, 'max_skip_steps', 6, 'background_method', 'simple','averageFac', 50); 
                % cyto area
                [channel_mean_cyto{ch}, channel_sum_cyto{ch}, channel_brightest_cyto{ch}, chbkg_mean_cyto{ch}, channel_median_cyto{ch}, cyto_area{fl}] = ...
                    intensityMeasureCyto(traj, nucLabel, properties_allframes, cytoLabel, channel_data, trackDirection, 'max_skip_steps', 6, 'background_method', 'simple', 'averageFac', 50);
            else % if measure other channel
                % nuc area            
                [channel_mean{ch}, channel_sum{ch}, channel_brightest{ch},chbkg_mean{ch},~, ~, channel_median{ch}]= ...
                    measureFuncHandle(traj, properties_allframes, channel_data, 'max_skip_steps', 6, 'background_method', 'simple', 'averageFac', 50); 
                % cyto area
                [channel_mean_cyto{ch}, channel_sum_cyto{ch}, channel_brightest_cyto{ch}, chbkg_mean_cyto{ch}, channel_median_cyto{ch}]= ...
                    intensityMeasureCyto(traj, nucLabel, properties_allframes, cytoLabel, channel_data, trackDirection, 'max_skip_steps', 6, 'background_method', 'simple','averageFac', 50);
            end
        end
       fprintf('finish measurement \n');
       
       %% for saving results
        position_mean{fl} = channel_mean;
        position_median{fl} = channel_median;
        position_sum{fl} = channel_sum;
        position_bright{fl} = channel_brightest;
        position_bkg{fl} = chbkg_mean;
        
        position_mean_cyto{fl} = channel_mean_cyto;
        position_median_cyto{fl} = channel_median_cyto;
        position_sum_cyto{fl} = channel_sum_cyto;
        position_bright_cyto{fl} = channel_brightest_cyto;
        position_bkg_cyto{fl} = chbkg_mean_cyto;
       
        
        % 
        if fl == 1
            profile off;
        end
        
        
    end % end of all positions   
    fprintf('start saving results at %f \n', toc);
    %% save results
    if saveLabelAndTraj
        for fl = 1:numb_files
            positionName = nucfiles(fl).name(PosNameStartIdx:CIndex-1);
            nucLabel = nucLabel_allPos{fl};
            cytoLabel = cytoLabel_allPos{fl};
            reg_output = reg_output_allPos{fl};
            traj = traj_allPos{fl};
            move_steps = move_steps_allPos{fl};
            properties_allframes = properties_allframes_allPos{fl};
            properties_allframes_cyto = properties_allframes_allPos_cyto{fl};
            save(fullfile(folderName, [positionName '_nucLabel.mat']), 'nucLabel', 'cytoLabel','reg_output');
            save(fullfile(folderName, [positionName '_trackingresults.mat']), 'traj', 'move_steps', 'properties_allframes','properties_allframes_cyto' ,'trackDirection');        
        end
        save(fullfile(folderName,  'IntensityMeasure_para.mat'), 'position_mean', 'position_median', 'position_sum','position_bright','position_bkg',  'position_median_cyto','position_mean_cyto', ...
            'position_sum_cyto','position_bright_cyto','numb_cells2keep', 'cells_centers', 'cells_numb2keep', 'cyto_area'); 
    end
    
    %% combine results of different position
    
    % nuc
    channel_sum_all = cell(length(channels2measure),1);
    channel_mean_all = cell(length(channels2measure),1);
    channel_median_all = cell(length(channels2measure),1);
    channel_bright_all = cell(length(channels2measure),1);
    channel_bkg_all = cell(length(channels2measure),1);
    
    % cyto
    channel_sum_cyto_all = cell(length(channels2measure),1);
    channel_mean_cyto_all = cell(length(channels2measure),1);
    channel_median_cyto_all = cell(length(channels2measure),1);
    channel_bright_cyto_all = cell(length(channels2measure),1);
    channel_bkg_cyto_all = cell(length(channels2measure),1);

    
    %% save
    for ch = 1:length(channels2measure) 
        % convert to matrix
        for fl = 1:numb_files
            % nuc
            channel_sum_all{ch}{fl} = position_sum{fl}{ch}';
            channel_mean_all{ch}{fl} = position_mean{fl}{ch}';
            channel_median_all{ch}{fl} = position_median{fl}{ch}';
            channel_bright_all{ch}{fl} = position_bright{fl}{ch}';
            channel_bkg_all{ch}{fl} = position_bkg{fl}{ch}';
            % cyto
            channel_sum_cyto_all{ch}{fl} = position_sum_cyto{fl}{ch}';
            channel_mean_cyto_all{ch}{fl} = position_mean_cyto{fl}{ch}';
            channel_median_cyto_all{ch}{fl} = position_median_cyto{fl}{ch}';
            channel_bright_cyto_all{ch}{fl} = position_bright_cyto{fl}{ch}';
            channel_bkg_cyto_all{ch}{fl} = position_bkg_cyto{fl}{ch}';

        end
        % nuc
        channel_sum = cell2mat(channel_sum_all{ch});
        channel_mean = cell2mat(channel_mean_all{ch});
        channel_median = cell2mat(channel_median_all{ch});
        channel_bright = cell2mat(channel_bright_all{ch});
        channel_bkg = cell2mat(channel_bkg_all{ch});
        % cyto
        channel_sum_cyto = cell2mat(channel_sum_cyto_all{ch});
        channel_mean_cyto = cell2mat(channel_mean_cyto_all{ch});
        channel_median_cyto = cell2mat(channel_median_cyto_all{ch});
        channel_bright_cyto = cell2mat(channel_bright_cyto_all{ch});
        channel_bkg_cyto = cell2mat(channel_bkg_cyto_all{ch});

        
        
        % save csv
        ch2measure = channels2measure_num(ch);
        ChName = channels(ch2measure).channelName;
        if strcmp(trackDirection, 'backward')
            % nuc
            csvwrite(fullfile(folderName, [ChName '_mean.csv']), flipud(channel_mean));
            csvwrite(fullfile(folderName, [ChName '_median.csv']), flipud(channel_median));
            csvwrite(fullfile(folderName, [ChName '_sum.csv']), flipud(channel_sum));
            csvwrite(fullfile(folderName, [ChName '_bright.csv']), flipud(channel_bright));
            csvwrite(fullfile(folderName, [ChName '_background.csv']), flipud(channel_bkg));
            % cyto
            csvwrite(fullfile(folderName, [ChName '_cyto_mean.csv']), flipud(channel_mean_cyto));
            csvwrite(fullfile(folderName, [ChName '_cyto_median.csv']), flipud(channel_median_cyto));
            csvwrite(fullfile(folderName, [ChName '_cyto_sum.csv']), flipud(channel_sum_cyto));
            csvwrite(fullfile(folderName, [ChName '_cyto_bright.csv']), flipud(channel_bright_cyto));
            csvwrite(fullfile(folderName, [ChName '_cyto_background.csv']), flipud(channel_bkg_cyto));

        else
            % nuc
            csvwrite(fullfile(folderName, [ChName '_mean.csv']), channel_mean);
            csvwrite(fullfile(folderName, [ChName '_median.csv']), channel_median);
            csvwrite(fullfile(folderName, [ChName '_sum.csv']), channel_sum);
            csvwrite(fullfile(folderName, [ChName '_bright.csv']), channel_bright);
            csvwrite(fullfile(folderName, [ChName '_background.csv']), channel_bkg);
            % cyto
            csvwrite(fullfile(folderName, [ChName '_cyto_mean.csv']), channel_mean_cyto);
            csvwrite(fullfile(folderName, [ChName '_cyto_median.csv']), channel_median_cyto);
            csvwrite(fullfile(folderName, [ChName '_cyto_sum.csv']), channel_sum_cyto);
            csvwrite(fullfile(folderName, [ChName '_cyto_bright.csv']), channel_bright_cyto);
            csvwrite(fullfile(folderName, [ChName '_cyto_background.csv']), channel_bkg_cyto);
        end

	%% find some bad traj, Nan calls it quality control

        %% get nuc file name, numb of frames, image width and height

% 	if ch2measure == nuc_ch_num
%         bad_cells = cell(1,numb_frame-1);
%         for fr = 1:numb_frame-1
% 			change = channel_bright(fr+1,:) - channel_bright(fr,:);
% 			mean_increment = mean(change);
% 			std_increment = std(change);
% 			bad_cells{fr} = find( (change>mean_increment+max_change_factor*std_increment) | (change<mean_increment-max_change_factor*std_increment) );
%         end
%         all_bad_cells = unique( cell2mat(bad_cells) );
%         csvwrite(fullfile(folderName,'cells2exclude.csv'), all_bad_cells);
% 	end
%         
    end
    
    fprintf('finish working on folder %s at %f \n', folderName, toc);

end %end of all folders

