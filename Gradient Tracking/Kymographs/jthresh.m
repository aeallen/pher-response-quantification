function [ output_args ] = jthresh( image, percentile )
%jthresh takes an image and a precentile value, and returns the value that
%falls at that percentile.  This value can then be used to threhold the
%image to leave the pixels above the indicated percentile.
    immax = max(max(nonzeros(image)));
    immin = min(min(nonzeros(image)));
           
    bins = [immin:(immax-immin)/98:immax];
    imagehist = histc(nonzeros(image),bins);
    if size(imagehist,1) > 0;
        target = sum(imagehist)*(percentile/100);
        currsum = 0;
        for i = 1:100;
            currsum = currsum + imagehist(i,1);
            if currsum >= target;
                targetbin = i;
                break;
            else
                targetbin = 100;
            end
        end
        output_args = bins(1,targetbin);
    else
        output_args = inf;
    

end

