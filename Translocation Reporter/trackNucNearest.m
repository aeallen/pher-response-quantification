function [traj, properties_allframes, varargout] = trackNucNearest(nucLabel,reg_output, varargin)

%% argument parse
% default
arg.max_move_distance = 10;
arg.trackDirection = 'forward'; % 
% parse
if nargin > 0
    parser = inputParser;
    addParameter(parser, 'max_move_distance', arg.max_move_distance);
    addParameter(parser, 'trackDirection', arg.trackDirection);
    parse(parser,varargin{:});
    arg = parser.Results;
end 

%% whether register or not
if all(all(reg_output==0)) 
    registerOrNot = false;
else
    registerOrNot = true;
end

%% get number of frames
numb_frame = size(nucLabel,3);

%% get index for different directions
switch arg.trackDirection
    case 'forward'
        track_frames = 1:1:numb_frame;
    case 'backward'
        track_frames = numb_frame:-1:1;
    otherwise
        error('wrong tracking direction');
end

%% initialize properties of all frames
properties_allframes = cell(1,numb_frame);

%% start tracking
for fr = 1:length(track_frames)-1    
    %% get all properties of current frames    
    if fr == 1 % do this for the first frame(depend on the tracking direction)        
        % get properties
        nucLabel_current = nucLabel(:,:,track_frames(fr));
        allProperties_current = regionprops(nucLabel_current,'Centroid', 'Area', 'PixelIdxList', 'PixelList');
        properties_allframes{fr} = allProperties_current(:); 
        
        % count the cells
        cellnum_current = length(allProperties_current);
        
        % initialize trajectories and move steps
        traj = zeros(cellnum_current,length(track_frames));
        traj(:,1) = 1:cellnum_current; 
        move_steps = zeros(cellnum_current,length(track_frames)-1);
        
        cumshift_x = zeros(1,cellnum_current);
        cumshift_y = zeros(1,cellnum_current);
    end  
    
    %% properties for the next frame
    nucLabel_next = nucLabel(:,:,track_frames(fr+1));
    allProperties_next = regionprops(nucLabel_next,'Centroid', 'Area', 'PixelIdxList', 'PixelList'); 
    cellnum_next = length(allProperties_next);
    
    % get the the shifted distance
    switch arg.trackDirection
        case 'forward'
            shift_correction_y = reg_output(track_frames(fr),3);
            shift_correction_x = reg_output(track_frames(fr),4);
        case 'backward'
            shift_correction_y = -reg_output(track_frames(fr)-1,3); % track_frames has one more element than reg_output
            shift_correction_x = -reg_output(track_frames(fr)-1,4);
    end
   
    %% calculate distances between cells in two frames and find nearest neighbor
    for cn = 1:cellnum_current
        coordinates_current = [ allProperties_current(cn).Centroid(1)  allProperties_current(cn).Centroid(2)];
        for cn2 = 1:cellnum_next
            coordinates_next = [ allProperties_next(cn2).Centroid(1)  allProperties_next(cn2).Centroid(2)];
            if registerOrNot
                % get current x, y, coordinates_next or current might have
                % other features, like area
                current_y = coordinates_current(2);
                current_x = coordinates_current(1);             
 
                % get the next x, y
                next_y_orig = coordinates_next(2);
                next_x_orig = coordinates_next(1);
                
                % make the correction
                next_y_shift = next_y_orig + shift_correction_y + cumshift_y(cn);
                next_x_shift = next_x_orig + shift_correction_x + cumshift_x(cn);
                
                % sum of squares of dx and dy
                sum_x2y2 = (current_x - next_x_shift)^2 + (current_y - next_y_shift)^2;
                
                if length(coordinates_next) >= 3 % if other feasures are considered
                    distance = sqrt(sum_x2y2 + sum( (coordinates_current(3:end) - coordinates_next(3:end)).^2 ) );
                else
                    distance = sqrt(sum_x2y2);
                end                    
            else
                distance = sqrt(sum( (coordinates_current - coordinates_next).^2 ) );
            end
            % find the nearest neighbor and the distance between them
            if cn2 == 1
                nearest_neighbor = 1;
                nearest_distance = distance;
            elseif distance < nearest_distance
                nearest_neighbor = cn2;
                nearest_distance = distance;
            end
        end % end of all cell in next

        if nearest_distance < arg.max_move_distance
            % store trajectory and moving steps
            traj(cn,fr+1) = nearest_neighbor;
            move_steps(cn,fr) = nearest_distance;
            
            % update cumulate shift
            cumshift_x(cn) = 0;
            cumshift_y(cn) = 0;
            
            % update properties
            try
                allProperties_current(cn) = allProperties_next(nearest_neighbor);
            catch
                length(allProperties_current)
                cn
                length(allProperties_next)
                nearest_neighbor
                fr
                error('check the segmentation of this frame');
            end
                
        else
            % if not found in the neiborhood, skip it, just save NaNs, no need to update
            % properties            
            traj(cn,fr+1) = NaN;
            move_steps(cn,fr) = NaN;
            % update cumulate shift
            cumshift_x(cn) = cumshift_x(cn) + shift_correction_x;
            cumshift_y(cn) = cumshift_y(cn) + shift_correction_y;
        end
    end % end of all cell in current
    
    properties_allframes{fr+1} = allProperties_current(:);    

end

switch nargout-2
    case 1
        varargout{1} = move_steps;
end

