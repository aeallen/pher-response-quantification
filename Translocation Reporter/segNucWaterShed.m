function [nucLabel_orig, reg_output] = segNucWaterShed(nuc_data, varargin)    
%% use watershed method to segment nuclei

	%% argument parse
	arg.bkg_roll_para = 10;
	arg.smooth_para = [4 2];
	arg.dst_para = 0; % choose 0 or 1
	arg.max_size = 200;
	arg.min_size = 15;
	arg.erode_para = []; 
    arg.dilate_para = []; 
    arg.imadjust_tol = 0.001; % lower it if has two less cell
    arg.imadjust_gama = 1;
    arg.threshold_adj = 1;
    arg.peakMethod = 'max'; % two different methods to find peaks, use first round dst map or find max intensity of nuc, use 'max' or 'dst'
    arg.registerOrNot = false;

    if nargin > 0
        parser = inputParser;
        addParameter(parser, 'bkg_roll_para', arg.bkg_roll_para); % should be optimize when have different cells
        addParameter(parser, 'smooth_para', arg.smooth_para); % highly recommend to optimize, for different cells, different exposure time, different fluorphore, 
        addParameter(parser, 'dst_para', arg.dst_para); % only for peakMethod dst, choose 0 or 1
        addParameter(parser, 'max_size', arg.max_size); % should be optimize when have different cells
        addParameter(parser, 'min_size', arg.min_size); % should be optimize when have different cells
        addParameter(parser, 'erode_para', arg.erode_para); % only for peakMethod dst, % highly recommend to optimize, for different cells, different exposure time, different fluorphore, 
        addParameter(parser, 'dilate_para', arg.dilate_para); % only for peakMethod dst, % highly recommend to optimize, for different cells, different exposure time, different fluorphore, 
        addParameter(parser, 'imadjust_tol', arg.imadjust_tol);
        addParameter(parser, 'imadjust_gama', arg.imadjust_gama);
        addParameter(parser, 'threshold_adj', arg.threshold_adj);



        addParameter(parser, 'peakMethod', arg.peakMethod); % if one does not work well, try the other, 'dst' method should be more robust if the image quality is bad
        addParameter(parser, 'registerOrNot', arg.registerOrNot); 

        parse(parser,varargin{:});
        arg = parser.Results;
    end 
%% get size, could be just one frame
    im_size = size(nuc_data);
    nucLabel_orig = zeros(im_size);
    ndim_nucdata = ndims(nuc_data);
    if  ndim_nucdata == 3
        numb_frame = im_size(3);
    elseif ndim_nucdata == 2
        numb_frame = 1;
    end
    
    %% initialize registration output
    reg_output = zeros(numb_frame-1,4);
    
%%
    
    %% threshold_adj could be a array, varies on frame
    if length(arg.threshold_adj) == 1
        threshold_adj = arg.threshold_adj*ones(1,numb_frame);
    else
        threshold_adj = arg.threshold_adj;
    end
    
    for fr = 1:numb_frame
                nuc_orig = nuc_data(:,:,fr);
                background = imopen(nuc_orig,strel('disk',arg.bkg_roll_para)); % 10 looks very good
                nuc_subb = nuc_orig - background;


                %% smooth a little
                nuc_filtered = imfilter(nuc_subb,fspecial('gauss', arg.smooth_para(1),arg.smooth_para(2)));

                %% adjust the contrast
                Low_High = stretchlim(nuc_filtered, arg.imadjust_tol);
                nuc_adjust = imadjust(nuc_filtered, Low_High,[], arg.imadjust_gama);

                %% convert to binary
%                 bw_orig = im2bw(nuc_adjust);
% new version
                threshold = graythresh(nuc_adjust)*threshold_adj(fr);
                bw_orig = im2bw(nuc_adjust, threshold);

                %% fill holes
                bw = imfill(bw_orig,'holes');
                %% dst map
                dst = -bwdist(~bw, 'euclidean');  % dst is negative

                %% find peaks, or centers
                switch arg.peakMethod
                    case 'dst'
                        pks = imextendedmin(dst,arg.dst_para); % find pks based on dst map

                    case 'max'             
                        threshold = multithresh(nuc_filtered,1);
                        nuc_quantized = uint16(imquantize(nuc_filtered, threshold) - 1); %bkg will be zeros
                        nuc_peaks = nuc_filtered.*nuc_quantized;        
                        pks = imregionalmax(nuc_peaks); 
                    case 'shrink'
                        pks = bwmorph(bw, 'shrink', inf);
                end % end of switch method
                
                if ~isempty(arg.erode_para)
                    bw = imerode(bw,arg.erode_para);
                end
                
                if ~isempty(arg.dilate_para)
                    bw = imdilate(bw,arg.dilate_para);
                end
                %% recaculate dst map
                dst2 = imimposemin(dst,pks); % important step. If just use the dst, messy

                %% watershed seg
                Ld = watershed(dst2); % gives very beautiful Ld
                bw(Ld == 0) = 0; 
                
                %% size exclusion
                bw_final  = bw &~bwareaopen(bw,arg.max_size); 
                bw_final = bwareaopen(bw_final,arg.min_size);                

                %% label and output
                if ndim_nucdata == 3
                    nucLabel_orig(:,:,fr) = bwlabel(bw_final); 
                elseif ndim_nucdata == 2
                    nucLabel_orig = bwlabel(bw_final); 
                end
                
                %% registration
                if arg.registerOrNot
                    if fr == 1
                        imreg_before = nuc_orig; 
                    else
                        %% do registration
                        % method 1, need the "Efficient subpixel image
                        % registration by cross-correlation"
                        imreg_current = nuc_orig;
                        reg_output(fr-1,:) = dftregistration(fft2(imreg_before),fft2(imreg_current),100);
                        imreg_before = imreg_current;
                    end  
                end
                 
                          
    end % end for all frames            



