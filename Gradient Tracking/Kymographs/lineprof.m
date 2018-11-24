function [ lineout ] = lineprof(image,mask,width,graph)
%lineprof uses a mask "mask" to do a linescan on the periphery of that mask in
%image "image".  The linescan will go in to the mask by "width" pixels


% create line


if sum(sum(mask)) > 0;
    %make the mask binary
    mask = mask > 0;

    masksmall = [];

    % structuring element for the erosion
    struct = strel('diamond', 1);

    masksmall{1,1} = imerode(mask,struct);

    % find the edge by subtracting the eroded mask from the mask, leaving a 1
    % pixel width line
    objedge= [];
    objedge{1,1} = double(mask)- double(masksmall{1,1});


    % determine the start point for the line, find the minimum x value in the
    % edge, then find the yvalues at that minumum x. 
    [y,x]=find(objedge{1,1});
    minx = min(x);
    [ymin,xmin] = find(objedge{1,1}(:,minx));

    %make boundary ROI (broi), a series of xy values defining lines for
    %improfile to use as the selection to measure, using minimum x and first y
    %associated with that x as a starting point
    broi = [];

    broi{1,1} = bwtraceboundary(objedge{1,1},[ymin(1),minx],'N');

    % define size of the line to use as number of points in the first (largest) line
    linesize = size(x,1);

    c1(:,1) = improfile(image,broi{1,1}(:,2),broi{1,1}(:,1),linesize);

    % figure();
    % plot(c1(:,1));

    %% to do a line of more than one pixel, repeat above steps
    if width > 1;
        for i = 2:width;
           masksmall{i,1}=imerode(masksmall{i-1,1},struct);
           objedge{i,1} = double(masksmall{i-1}) - double(masksmall{i,1});
           [yi,xi]=find(objedge{i,1});
           minxi = min(xi);
           [yimin,ximin]=find(objedge{i,1}(:,minxi));
           try
           broi{i,1}=bwtraceboundary(objedge{i,1},[yimin(1),minxi],'N');
           %makes the profile use the length of the longest line so that the
           %vectors can be averaged
           % use try statement, because if the perimeter collapses from
           % erosion, it breaks the boundary trace, which will still return a
           % boundary, but it won't work in improfile.  In situations where
           % this occurs, the script will just average the number of lines that
           % it can
           catch
               display(['Line Broken at ',num2str(i)]);
               break;
           end
           
           try
            c1(:,i)=improfile(image,broi{i,1}(:,2),broi{i,1}(:,1),linesize);
           catch
               display(['line broken at ',num2str(i)]); % will let the user know if this error is happening
               break;
           end
        end


    else
    end



    % output the average profile
    lineout = mean(c1,2);


    if graph >0;
        figure();
        plot(lineout);
        title('Line Out');
    else
    end
else
    lineout = 0;
end




end

