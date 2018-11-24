% wholecell_cap description, Joshua B Kelley, Elston Lab, UNC Chapel Hill
%
% Inputs for wholecell_cap are 3 tifs, mask, gfp, cherry
% •	maskdir, directory location of *mask.tif
% •	gfpdir, directory location of *gfp.tif
% •	cherrydir, directory location of *cherry.tif
% 
% 
% The threshold percentile 'threshper' and the lower threshhold percentile,
% 'threshperl' are used to calculate the polar cap angle of orientation.
% If the angle of orientation does not seem to be working, attempt
% different threshold values.
%
% Useful Outputs from wholecell_cap
% •	TLmask, a (time x 1) cell of sorted and labeled cell masks [x,y]
% •	gfpin, gfp images
% •	cherryin, cherry images
% •	gfpcell, a (cell# x 1) cell of matrices [x,y,t] which is just masked gfp data for each image
% •	celldata, a cell of (cell# x 1) containing matrices of [cellx,celly,peakx,peaky,angle,isreal] columns 3,4,5 are not currently used
% •	centroids, a cell of (cell# x 1) containing matrices of [highx,highy,lowx,lowy] which are the x and y values for the centroids of the high and low thresholded bem1 to determine cap angle
% •	cellangles, a cell (cell# x 1) of matrices (time x 1) containing angle data in degrees

%% Pick images to be input for analysis

maskdir = '/Users/AmyAllen/Documents/ThesisWork/pher-response-quantification/Gradient Tracking/Kymographs/Example data/wildtype_mask.tif';
gfpdir = '/Users/AmyAllen/Documents/ThesisWork/pher-response-quantification/Gradient Tracking/Kymographs/Example data/wildtype_GFP.tif';
% cherrydir = 'D:\MATLAB\wholecell\WT_0_150\9_18_12_s4_cdc3.tif';

preimin = cell(145,1);

tmin = 1;
tmax = [];
tmax = size(imfinfo(maskdir),1);
times = [tmin:1:tmax];

threshper = 99;
threshperl = 88;

maskin01 = [];
trackvar = [];

%% load images in
for c = 1:3;
    if c == 1;
        gadir = maskdir;
    elseif c==2;
        gadir = gfpdir;
%     elseif c == 3;
%         gadir = cherrydir;
    end
    
        
    for i = 1:tmax;
        preimin{i,1} = imread(gadir,i);
    end


    for i = 1:tmax
        if c ==1;
            maskin{i,1} = preimin{i,1};
        elseif c == 2;
            gfpin{i,1}= preimin{i,1};
%         elseif c == 3;
%             cherryin{i,1}=preimin{i,1};
        end
    end
end

    


% figure();
% subplot(1,3,1), imagesc(maskin{tmax,1});
% subplot(1,3,2), imagesc(gfpin{tmax,1});
% subplot(1,3,3), imagesc(cherryin{tmax,1});

tnum = tmax - tmin + 1;
max_field = zeros(tnum,1);
peak_thresh = cell(tmax,1);
currmax = [];
low_thresh_clean = cell(tnum,1);
low_thresh_count = zeros(tnum,1);
peak_thresh_count = zeros(tnum,1);
labeledmask_peak = cell(tnum,1);
labeledmask = cell(tnum,1);


% find the maximum value present in each time frame
% for i = 1:tnum;
%     max_field(i,1) = max(max(imin{i,1}));
% end

% convert the mask to a logical
for i = 1:tnum;
    maskin01{i,1} = ~maskin{i,1} > 0; %mask out of imagej is inverted
    cellmask{i,1} = uint16(maskin01{i,1});
end



figure();
% subplot(1,2,1), 
imagesc(gfpin{tmax,1}.*cellmask{tmax,1});
% subplot(1,2,2), imagesc(cherryin{tmax,1}.*cellmask{tmax,1});


%% label the mask and keep track of number of cells in each time point
for i = 1:tmax;
    labeledmask{i,1} = bwlabel(maskin01{i,1});
    cellcount(i,1) = max(max(labeledmask{i,1}));
end



%% use jtrack to sort the cells.  Use x,y,time,label
total = 1;
for i = 1:tnum;
    for j = 1:cellcount(i,1);
        currmask = [];
        currmask = labeledmask{i,1}==j;
        [ys,xs] = find(currmask); 
        trackvar(total,1) = mean(xs);
        trackvar(total,2) = mean(ys);
        trackvar(total,3) = i;
        trackvar(total,4) = j;
        total = total+1; %keeps track through the larger loop
    end
end

tracked = jtrackv3(trackvar,0,100); %had to up the max distance to 30 pixels for the tracking to work

%% relabel the mask based on cell identity. T(rack)L(abel)mask
% start with tracked cell 1, and assign all of its indices to 1 in the new
% mask.
TLmask = cell(tmax,1);
for i = 1:tmax;
    TLmask{i,1} = zeros(size(labeledmask{i,1},1),size(labeledmask{i,1},2));
end

for i = 1:size(tracked,1);
    for t = 1:tnum;
        if tracked{i,1}(t,4)==1;
        currentID = tracked{i,1}(t,3);
        currmask = labeledmask{t,1} == currentID;
        TLmask{t,1}(currmask) = i;
        else
        end
    end
end

%% Output a TIFF of the tracked mask for trouble shooting/validation
imwrite(uint8(TLmask{1,1}),'tracked_mask.tif', 'WriteMode', 'OverWrite');
for i = 2:tmax;
    imwrite(uint8(TLmask{i,1}),'tracked_mask.tif','WriteMode', 'append');
end

%% Threshold the GFP within each cell

% separate each cell's gfp image out
gfpcell = [];

for i = 1:size(tracked,1);
    for t = 1:tnum;
        gfpcell{i,1}(:,:,t)= gfpin{t,1}.*uint16(TLmask{t,1}==i);
    end
end
%% Calculate thresholded images for a high value and low value(per cell per timepoint)and
%  a single value which will work for all timepoints

threshper = 95;
threshperl = 88;

gfpthresh_hi = [];
gfpthresh_low = [];
per_threshs =[];
for i = 1:size(tracked, 1);
    for t = 1:tnum;
        if size(nonzeros(gfpcell{i,1}(:,:,t)),1) > 2;
            gfpthresh_hi{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > jthresh(gfpcell{i,1}(:,:,t), threshper); 
            gfpthresh_low{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > jthresh(gfpcell{i,1}(:,:,t), threshperl);
            per_threshs{i,1}(t,1) = jthresh(gfpcell{i,1}(:,:,t), threshper); %store the values of the jthresh output
            per_threshs{i,1}(t,2) = jthresh(gfpcell{i,1}(:,:,t), threshperl);
        else
            gfpthresh_hi{i,1}(:,:,t) = zeros(size(gfpcell{i,1}(:,:,t),1),size(gfpcell{i,1}(:,:,t),2)); 
            gfpthresh_low{i,1}(:,:,t) = zeros(size(gfpcell{i,1}(:,:,t),1),size(gfpcell{i,1}(:,:,t),2));
            per_threshs{i,1}(t,1) = 0; %store the values of the jthresh output
            per_threshs{i,1}(t,2) = 0;   
        end
        
    end
end


gfpconstantthresh = [];


%% Find Area of thresholded Bem1


for i = 1:size(tracked,1);
    for t=1:tnum;
        if sum(sum(gfpcell{i,1}(:,:,t)))>0;
        gfpconstantthresh{i,1}(:,:,t) = gfpcell{i,1}(:,:,t) > min(min(per_threshs{i,1}));
        else
        end
    end
end

        


%% find the centroid of each cell
celldata = cell(size(tracked,1),1); % cellx, celly, peakx, peaky, angle, isreal

for i = 1:size(tracked,1);
    for t=1:tnum;
        % find the centroid of the cell mask        
        if sum(find(TLmask{t,1}==i))>0;
            [ys,xs] = find(TLmask{t,1} == i);
            celldata{i,1}(t,1) = mean(xs);
            celldata{i,1}(t,2) = mean(ys);
            celldata{i,1}(t,6) = 1;
        else
           celldata{i,1}(t,1) = 0;
           celldata{i,1}(t,2) = 0; 
           celldata{i,1}(t,6) = 0;
        end
    end
end

figure();
imagesc(TLmask{end,1});
hold on
for i = 1:size(tracked,1);
text(celldata{i,1}(end,1),celldata{i,1}(end,2), ['Cell',num2str(i)], 'FontSize', 14);
end


%% calculate the angle of the polar cap in each cell based on the high and 
% low thresholds

centroids = cell(size(tracked,1),1); % centroids cell of [highx, highy, lowx, lowy];


for i = 1:size(tracked,1);
    for t = 1:tnum;
        if sum(sum(gfpthresh_hi{i,1}(:,:,t))) > 0;
        [ys,xs] = find(gfpthresh_hi{i,1}(:,:,t));
        centroids{i,1}(t,1)=mean(xs);
        centroids{i,1}(t,2)=mean(ys);
        else
        centroids{i,1}(t,1) = 0;
        centroids{i,1}(t,2) = 0;
        end
        if sum(sum(gfpthresh_low{i,1}(:,:,t))) > 0;
        [ys,xs] = find(gfpthresh_low{i,1}(:,:,t));
        centroids{i,1}(t,3) = mean(xs);
        centroids{i,1}(t,4) = mean(ys);
        else
        centroids{i,1}(t,3) = 0;
        centroids{i,1}(t,4) = 0;    
        end
        
    end
end

%% calculate the angles from the xy positions

cellangles = cell(size(tracked,1),1);

for i = 1:size(tracked,1);
    for t = 1:tnum;
        if sum(sum(centroids{i,1}))>0;
            cellangles{i,1}(t,1) = atan2d(centroids{i,1}(t,2)-centroids{i,1}(t,4),centroids{i,1}(t,1)-centroids{i,1}(t,3));
        else
            cellangles{i,1}(t,1) = NaN;
        end
    end
end

        







