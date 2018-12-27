function [phaseLabel, nucRelabel] = segmentPhaseYeast(phase_data, nuc_label, varargin)

arg.segment_method = 'combine';
arg.outerMkrOrNot = true;

if nargin > 0
    parser = inputParser;   
    addParameter(parser, 'segment_method', arg.segment_method);  
    addParameter(parser, 'outerMkrOrNot', arg.outerMkrOrNot); 
    parse(parser,varargin{:});
    arg = parser.Results;
end 

%% size
im_size = size(phase_data);
ndim= ndims(nuc_label);
if  ndim == 3
    numb_frame = im_size(3);
elseif ndim == 2
    numb_frame = 1;
end

% initialize return data
phaseLabel = zeros(im_size);
nucRelabel = zeros(im_size);
for fr = 1:numb_frame
    %% get data for this frame
    phase_positive = phase_data(:,:,fr);
    nucLabel = nuc_label(:,:,fr);
    nuc_bw_final = logical(nucLabel); 

    %% use markor or not
    if arg.outerMkrOrNot
        % find ourter edges
        edges = cellEdgeFinderYeast(phase_positive);        
    end
    
    phase_fill = imfill(edges, 'holes');
    
    switch arg.segment_method
        case 'regular'
        % distance map, method1
        dst = -bwdist(~phase_fill, 'euclidean');  % dst is negative
        dst2 = imimposemin(dst, nuc_bw_final);
        
        L = watershed(dst2);
        negative_seg = phase_fill;
        negative_seg(L==0) = 0;
       

        case 'geodesic'
        % method 2, geodesic distance map
        dst = bwdistgeodesic(phase_fill,nuc_bw_final, 'quasi-euclidean');
        dst(isnan(dst)) = -Inf;
        dst(isinf(dst)) = Inf;  

        L = watershed(dst);
        L(~phase_fill | isinf(dst)) = 0;
        negative_seg = phase_fill;
        negative_seg(L==0) = 0;

    
        case 'combine'
        % combine the results of two methods
         % first
        dst = -bwdist(~phase_fill, 'euclidean');  % dst is negative
        dst2 = imimposemin(dst, nuc_bw_final);
        
        L = watershed(dst2);
        negative_seg = phase_fill;
        negative_seg(L==0) = 0;
        negative_seg(edges) = 0;
        negative_seg1 = negative_seg;
        
         % second
        dst = bwdistgeodesic(phase_fill,nuc_bw_final, 'quasi-euclidean');
        dst(isnan(dst)) = -Inf;
        dst(isinf(dst)) = Inf;  

        L = watershed(dst);
        L(~phase_fill | isinf(dst)) = 0;
        negative_seg = phase_fill;
        negative_seg(L==0) = 0;
        negative_seg2 = negative_seg;
        
         % combine
        negative_seg = negative_seg1 & negative_seg2;    
        otherwise
            error('wrong cyto segmentation method');
    end
    
    %% relabel cyto
    negative_label = bwlabel(negative_seg, 4);
    [nucRelabel(:,:,fr), phaseLabel(:,:,fr)] = cytoRelabel_v2(nucLabel, negative_label);

end




