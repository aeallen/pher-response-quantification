% This script will take the angle output from wholecell_cap and adjust the
% angles when the angle crosses the 180/-180 degree line.  If an angle is
% positive and it crosses 180, it will be converted to an angle larger than
% 180, rather than an angle between 0 and -180.  The opposite is true for
% an angle that starts out negative.

angle_diff = cell(size(tracked,1),1);
cellsfixed = cellangles;

for i= 2:size(cellangles,1);
    for t = 2:size(cellangles{i,1},1);
        if abs(cellsfixed{i,1}(t,1)-cellsfixed{i,1}(t-1,1)) > 180;
           anglenew = (180-abs(cellangles{i,1}(t,1))) + (180-abs(cellangles{i,1}(t-1,1))); % calculate the absolute difference of both angles from 180 
           
           if cellangles{i,1}(t-1,1) >= 0;
               cellsfixed{i,1}(t,1) = cellsfixed{i,1}(t-1,1) + anglenew; % if the angle is positive, keep it going positive
           else
               cellsfixed{i,1}(t,1) = cellsfixed{i,1}(t-1,1) - anglenew; % if the angle is negative, keep it going negative
           end
       else
        end
    end
end



