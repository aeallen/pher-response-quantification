function [ cells ] = jtrack( polarcap, printplot, mindist )
%jtrack This function will sort data into cells based on x,y positions, and
%time point data.  It will co-sort a 4th data column into those cells.
%Input for "polarcap" is a matrix with columns 
%[ Xpos, Ypos, Slice(timepoint), Data]
% jtrack starts with the final time point, and sorts backwards.  Only cells
% present in the last time point will be sorted



if nargin < 2;
printplot = 0;    
mindist = 10;
else
end

if nargin < 3
    mindist = 10;
else
end

polarcapsize = size(polarcap);

timemax = max(polarcap(:,3));

time_points =[1:1:timemax];

numberofcells = find(polarcap(:,3) == timemax);

cellnum = size(numberofcells,1);


%Create a cell with an array for each yeast to be tracked with 4 columns,
%x,y,angle,null where the 4th column is 0 if there is no data at that time
%point, and is 1 if there is
for i = [1:cellnum];
    cells(i,1)={zeros(timemax,4)};

end

%% Populate the end time point cell data, x and y positions

[row,col]= find(polarcap(:,3) == timemax);
for i = [1:cellnum];
    cells{i,1}(timemax,:)= [polarcap(row(i),1), polarcap(row(i),2), polarcap(row(i),4),1];
end

%% assign x,y values based on distance from previous frame

for j = [ timemax-1 : -1 : 1]

    %initialize values
    row = 0;
    col = 0;

    %get values for current time point
    [row,col] = find(polarcap(:,3) == j);

    prevtime = j+ 1;
    currentcellnum = size(row,1);
    currentpos = zeros(cellnum, 2);

   %set current positions equal to total number of cells zeroedls,
   %otherwise set the current positions to current cell number, still
   %zeroed
    if currentcellnum == 0;
    % set values equal to zero if there is no spot detected
        for k = 1:cellnum;
            cells{k,1}(j,:) = [0,0,0,0];
        end
    else
        %make a list of the current positions
            currentpos = [zeros(currentcellnum,3)];
        for i= [1:currentcellnum]

            currentpos(i,:) = [polarcap(row(i), 1), polarcap(row(i), 2), polarcap(row(i), 4)];
        end




        %calculate differences from time j+1 to j for all of the current
        %positions, for cell #k and then find the minimum
        diff_mat = cell(cellnum,1);
        dist = zeros(cellnum,currentcellnum);
        distindex = zeros(cellnum,currentcellnum);
        for k= [1:cellnum]
            if cells{k,1}(prevtime,1) == 0;
               
                %shift the look back time to a non zero position
                for m = 2:(timemax-1);
                    prevtime = j + m;
                    if cells{k,1}(prevtime,1) == 0;
                    else
                        break;
                    end
                end
            else
            end
            %make a difference cell
                  for i = [1:currentcellnum]
                      %calculate dx and dy 
                 diff_mat{k,1}(i,1) = abs((cells{k,1}(prevtime,1)) - currentpos(i,1));
                 diff_mat{k,1}(i,2) = abs((cells{k,1}(prevtime,2)) - currentpos(i,2));
                    %square root of dx^2 + dy^2
                 dist(k,i) = sqrt((diff_mat{k,1}(i,1))^2 + ((diff_mat{k,1}(i,2))^2));
                 end
                      

        end
        
        %% compare distances to find minimum and assign object
        % dist: rows are cells, columns are unassigned objects
        
        [mincellval,mincellind] = min(dist,[],2);
        [minobjval,minobjind] = min(dist,[],1);
        
        % determine minimum values that are higher than the mindist
        % variable
        
        thresholdedcell = mincellval <= mindist;
        thresholdedobj = minobjval <= mindist;
        
        celldup = histc(mincellind,[1:1:currentcellnum]);
        objdup = histc(minobjind, [1:1:cellnum]);
        
        for p = 1:cellnum;
            if thresholdedcell(p,1) == 1 && celldup(mincellind(p,1),1) == 1; % the minimum distance is within the mindist threshold and the cell has only been assigned once
                % check the celldup histogram, looking at the row that has
                % the same value as the index given
                cells{p,1}(j,:)= [currentpos(mincellind(p,1),1:3),1];
            else if thresholdedcell(p,1) == 1 && celldup(mincellind(p,1),1) > 1;
                    if minobjind(1,mincellind(p,1)) == p;
                        cells{p,1}(j,:)= [currentpos(mincellind(p,1),1:3),1];
                    else
                        cells{p,1}(j,:)= [0,0,0,0];
                    end
                   
                else
                 cells{p,1}(j,:)= [0,0,0,0];
                end
            end
        end        
    end
end

%% Plot the resultant XY positions
if printplot == 1;
    
figure();
hold on;
% limited number of colors, if too many cells are present, graphing it will
% cause a crash, but the sorted values should be fine.  If you need more
% colors, add colors or write a script to pick a color no matter what
colors = ['b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r','c','m','k','y','b','g','r','c','m','k','y','b','g','r','c','m','k','y','b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r','c','m','k','y','b','g','r','c','m','k','y','b','g','r','c','m','k','y','b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g', 'r','c','m','k','y','b','g','r','c','m','k','y','b','g','r','c','m','k','y'];
axis ij;
for i = 1: cellnum;
    plot(cells{i,1}(:,1),cells{i,1}(:,2),'LineStyle','none','Marker','o', 'MarkerFaceColor', colors(1,i), 'MarkerEdgeColor', 'none')
    text(cells{i,1}(end,1),cells{i,1}(end,2), ['Cell',num2str(i)], 'FontSize', 14)
end

hold off

else
end


end

