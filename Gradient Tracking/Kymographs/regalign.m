function [ bem1out,cdc3out] = regalign( bem1, cdc3, windowsize, doplot )
%REGALIGN Aligns bem1 and cdc3 to peak Bem1, and returns aligned and
%normalized matrices for each.

%    This function takes as input two identical sized matrices of
%intensity data, position in the Y, and time in the X.  It will then align
%both matrices based on the position of maximum intensity in the bem1
%matrix for each time point, adjusting the cdc3 matrix to maintain the
%relative spatial orientation between the two matrices.  The output will be
%centered at the peak with data extending on both sides to an extent
%determined by the windowsize.  Default windowsize is 41 pixels, 20 pixels
%to either side of the center pixel


if nargin <3;
    windowsize = 41;
    doplot = 0;
elseif nargin < 4;
    doplot = 0;
end

smoothfactor = 1;

timepoint = size(bem1, 2);

halfwin = floor(windowsize/2);
columnsum_bem1 = [];
columnsum_cdc3 = [];
bem1region_norm = [];
cdc3region_norm = [];

%% Align the matricies against the brightest point centered using the kymo_align script

timemax = size(bem1,2);
bem1diff = [];
curr_line_bem1 = [];
curr_size = [];
prev_size = [];
prezero = [];
post_zero = [];
bem1fixed = [];
cdc3fixed = [];
prev_line = [];
curr_line_bem1 = [];
prev_mat = [];
curr_mat_bem1 = [];
prev_mat_sized = [];
curr_mat_sized_bem1 = [];
score = [];
startmaxval = [];
startmaxind = [];
startsize = [];



bem1fixed = zeros(size(bem1, 1), size(bem1, 2));
cdc3fixed = zeros(size(bem1, 1), size(bem1, 2));

for i = 1:size(bem1fixed, 2);
    if sum(bem1(:,i),1) > 0;
        starttime = i;
        break
    else
    end
end

% set starting column
% pick max intensity and set that to center
start_preline_ind = find(bem1(:,starttime));
bem1_start_preline = bem1(start_preline_ind, starttime);
cdc3_start_preline = cdc3(start_preline_ind, starttime);
% add an average step to smooth the bem1
[startmaxval, startmaxind] = max(smooth(bem1_start_preline, smoothfactor));

prelinesize = size(bem1_start_preline, 1);
mid = floor(prelinesize/2);
bem1startline = zeros(prelinesize, 1);
cdc3startline = zeros(prelinesize, 1);

if startmaxind < mid;
    bem1startline(mid+1-startmaxind:end,1) = bem1_start_preline(1:end-(mid-startmaxind),1);
    bem1startline(1:mid-startmaxind, 1) = bem1_start_preline(end+1-(mid-startmaxind):end,1);
    cdc3startline(mid+1-startmaxind:end,1) = cdc3_start_preline(1:end-(mid-startmaxind),1);
    cdc3startline(1:mid-startmaxind, 1) = cdc3_start_preline(end+1-(mid-startmaxind):end,1);
    
elseif startmaxind > mid;
    bem1startline(1:end-(startmaxind-mid),1) = bem1_start_preline(startmaxind-mid+1:end,1);
    bem1startline(end-(startmaxind-mid)+1:end, 1) = bem1_start_preline(1:startmaxind-mid,1);
    
    cdc3startline(1:end-(startmaxind-mid),1) = cdc3_start_preline(startmaxind-mid+1:end,1);
    cdc3startline(end-(startmaxind-mid)+1:end, 1) = cdc3_start_preline(1:startmaxind-mid,1);
else
    bem1startline = bem1_start_preline;
    cdc3startline = cdc3_start_preline;
end

startprezero = zeros(floor((size(bem1,1)-prelinesize)/2),1);
startpostzero = zeros((size(bem1,1)-prelinesize-size(startprezero,1)),1);
bem1finalstart = [startprezero; bem1startline; startpostzero];
cdc3finalstart = [startprezero; cdc3startline; startpostzero];


bem1fixed(:,starttime) = bem1finalstart;
cdc3fixed(:,starttime) = cdc3finalstart;

for j = (starttime + 1): size(bem1fixed,2);
    if max(bem1(:,j))>0;
        nonzerosbem1 = find(bem1(:,j));
        curr_line_bem1 = bem1(nonzerosbem1,j);
        curr_line_cdc3 = cdc3(nonzerosbem1,j);
        curr_size = size(curr_line_bem1,1);
        prev_line = nonzeros(bem1fixed(:,j-1));
        prev_size = size(prev_line,1); 

        %determine larger line, and set that size to the current window,
        %calculate size difference

        if curr_size > prev_size;
            win_size = curr_size;
            sizediff = curr_size - prev_size;
            curr_big = 1;
        else
            win_size = prev_size;
            sizediff = prev_size - curr_size;
            curr_big = 0;
        end

        %Replicate the current line shifting by one each time and wrapping the end back
        %around  

        curr_mat_bem1 = zeros(curr_size,curr_size);
        curr_mat_cdc3 = zeros(curr_size,curr_size);

        for i = 1:curr_size;
            curr_mat_bem1(i:curr_size,i) = curr_line_bem1(1:curr_size+1-i,1);
            curr_mat_cdc3(i:curr_size,i) = curr_line_cdc3(1:curr_size+1-i,1);
            if i>1;
                curr_mat_bem1(1:i-1,i) = curr_line_bem1(curr_size+2-i:curr_size,1);
                curr_mat_cdc3(1:i-1,i) = curr_line_cdc3(curr_size+2-i:curr_size,1);
            else
            end
        end


        if curr_big == 1;
            prezero = zeros(floor(sizediff/2), size(curr_mat_bem1,2));
            postzero = zeros((sizediff-floor(sizediff/2)), size(curr_mat_bem1,2));
            prev_mat = repmat(prev_line, 1, size(curr_mat_bem1, 2));
            prev_mat_sized = [prezero; prev_mat; postzero];
            curr_mat_sized_bem1 = curr_mat_bem1;
            curr_mat_sized_cdc3 = curr_mat_cdc3;


        else
            prev_mat = repmat(prev_line, 1, size(curr_mat_bem1, 2));
            prezero = zeros(floor(sizediff/2), size(prev_mat, 2));
            postzero = zeros((sizediff-floor(sizediff/2)), size(prev_mat, 2));
            prev_mat_sized = prev_mat;
            curr_mat_sized_bem1 = [prezero; curr_mat_bem1; postzero];
            curr_mat_sized_cdc3 = [prezero; curr_mat_cdc3; postzero];
        end
        % score is calculated using bem1 and then the same shift is applied to
        % cdc3

    % score based on maximum
        [valscore, score_mat] = max(curr_mat_sized_bem1);
        midpoint = floor(size(curr_mat_sized_bem1, 1)/2);
        midpoint = ones(1,size(score_mat,2))*midpoint;
        score = abs(midpoint-score_mat);
        [val,al_c] = min(score);

        prezerofinal = zeros(floor((size(bem1fixed,1)-size(curr_mat_sized_bem1,1))/2),1);
        postzerofinal = zeros( (size(bem1fixed,1) - size(curr_mat_sized_bem1,1)-size(prezerofinal,1)) , 1);
        bem1fixed(:,j) = [prezerofinal; curr_mat_sized_bem1(:,al_c); postzerofinal];
        cdc3fixed(:,j) = [prezerofinal; curr_mat_sized_cdc3(:,al_c); postzerofinal];
    else
        bem1fixed(:,j) = zeros(size(bem1fixed,1),1);
        cdc3fixed(:,j) = zeros(size(cdc3fixed,1),1);
    end
    
end


%% Make the windowed kymograph
bem1region = zeros(windowsize, size(bem1fixed,2));
cdc3region = zeros(windowsize, size(bem1fixed,2));

bem1fixedsmooth = [];

% the "fixed" variables should already be centered for the size of the cell
maxind = floor(size(bem1fixed,1)/2);

if size(bem1fixed,1) < windowsize;
    startind = floor((windowsize-size(bem1fixed,1))/2);
    finalind = windowsize - startind - size(bem1fixed,1);
    bem1region(startind:finalind,:) = bem1fixed;
    cdc3region(startind:finalind,:) = cdc3fixed;
elseif size(bem1fixed,1) > windowsize;
    startind = floor((size(bem1fixed,1)-windowsize)/2);
    finalind = startind + windowsize-1;
    bem1region(:,:) = bem1fixed(startind:finalind,:);
    cdc3region(:,:) = cdc3fixed(startind:finalind,:);
elseif size(bem1fixed,1) == windowsize;
    bem1region = bem1fixed;
    cdc3region = cdc3fixed;
end
%% output is normalized to sum to 1
columnsum_bem1 = repmat(sum(bem1region), size(bem1region,1), 1);
columnsum_cdc3 = repmat(sum(cdc3region), size(cdc3region,1), 1);
bem1region_norm = bem1region ./ columnsum_bem1;
cdc3region_norm = cdc3region ./ columnsum_cdc3;
bem1cdc3_tcorr = timecorr(bem1region_norm, cdc3region_norm);

bem1out = bem1region_norm;
cdc3out = cdc3region_norm;


if doplot ==1;
figure();
subplot(1,3,1), imagesc(bem1region_norm);
title('Bem1');
subplot(1,3,2), imagesc(cdc3region_norm);
title('Cdc3');
subplot(1,3,3), plot(bem1cdc3_tcorr);
title('Bem1 Cdc3 Correlation over Time');
ylim([-1,1]);
xlabel('Time Point (5 min)');
ylabel('Correlation Coefficient');
else
end


end




