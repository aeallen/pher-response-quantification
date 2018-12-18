% combokeeper combines the selected cells, collating angle and gfp/cherry kymographs
% creates "combined" which is a cell of {[angle], [GFP kymo], [cherry
% kymo], centroids[highx, highy, lowx, lowy]]}

% Pick "keepers", a variable containing the cell number of the cells you wish to analyze.
keepers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];

% pick linethickness for the width of the perimeter averaging
linethickness = 4;

kept=[];
keptfixed = [];
keptmat = [];
areas = [];

% Combined cell contains {cellangles,gfpkymo,cherrykymo}
combined = cell(size(keepers,2),4);
% profiles is an intermediate to hold raw kymographs
profiles = cell(size(keepers,2),2);

%load with excess zeros so data can fit, don't know maximum size of the
%kymograph until its been generated
for i = 1:size(keepers,2);
    profiles{i,1} = zeros(400,tmax);
    profiles{i,2} = zeros(400,tmax);
end

%% Generate Line profile from the masks and GFP and Cherry data
for i = 1:size(keepers,2);
    combined{i,1}= cellsfixed{keepers(i),1};
    combined{i,4} = centroids{keepers(i),1};
    for t = 1:size(TLmask,1);
        gfptemp = [];
%         cherrytemp =[];
        gfptemp = lineprof(gfpin{t,1},(TLmask{t,1}==keepers(i)),linethickness,0);
%         cherrytemp = lineprof(cherryin{t,1},(TLmask{t,1}==keepers(i)), linethickness,0); 
        profiles{i,1}(1:size(gfptemp,1),t) = gfptemp;
%         profiles{i,2}(1:size(cherrytemp,1),t) = cherrytemp;
    end
    
end

%% find maximum profile length for each cell, center and normalize kymograph
maxsize = [];
tempsize = [];
for i = 1:size(profiles,1);
    for t = 1:size(profiles{i,1},2);
        if sum(sum(profiles{i,1}(:,t))) > 1;
        tempsize(i,t) = find(profiles{i,1}(:,t),1,'last');
        else
            tempsize(i,t) = 0;
        end
        
    end
    maxsize(i,1) = max(tempsize(i,:));
    profiles{i,1} = profiles{i,1}(1:maxsize(i,1),:);
%     profiles{i,2} = profiles{i,2}(1:maxsize(i,1),:);
end

% Fill empy matrices with correct number of zeros to load centered data
centeredprofiles = cell(size(profiles,1),2);
for i = 1:size(centeredprofiles,1);
    centeredprofiles{i,1}=zeros(maxsize(i,1),size(profiles{i,1},2));
%     centeredprofiles{i,2}=zeros(maxsize(i,1),size(profiles{i,1},2));
end

% Center the data

for i = 1:size(profiles,1);
    for t = 1 : size(profiles{i,1},2);
    currentzeros = sum(profiles{i,1}(:,t)==0);
        if currentzeros > 1;
            prespace = floor(currentzeros/2);
            postspace = currentzeros - prespace;
            postindex = size(profiles{i,1},1)-postspace;
            centeredprofiles{i,1}(prespace:postindex,t)= profiles{i,1}(1:postindex-prespace+1,t);
%             centeredprofiles{i,2}(prespace:postindex,t)= profiles{i,2}(1:postindex-prespace+1,t);
        else
            centeredprofiles{i,1}(:,t) = profiles{i,1}(:,t);
%             centeredprofiles{i,2}(:,t) = profiles{i,2}(:,t);
        end
    end
    
   
end

%% Normalize the Data
normprofiles = [];
for i = 1:size(centeredprofiles,1);
    normprofiles{i,1} = (centeredprofiles{i,1} - min(nonzeros(centeredprofiles{i,1}))) ./ max(max(centeredprofiles{i,1})).*double(centeredprofiles{i,1}>0);
%     normprofiles{i,2} = (centeredprofiles{i,2} - min(nonzeros(centeredprofiles{i,2}))) ./ max(max(centeredprofiles{i,2})).*double(centeredprofiles{i,2}>0);
    combined{i,2} = normprofiles{i,1};
%     combined{i,3} = normprofiles{i,2};
end


%% Example of output

for i = 1:size(combined,1);
    figure();
    imagesc(combined{i,2});
    title('Bem1')
end



