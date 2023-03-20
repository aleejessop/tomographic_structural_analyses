% function to calculate for a given b/w data set 
%    1. EDM of the 0 phase (as per bwDist)
%    2. index map (as per bwDist)
%    3. EDMValid : indicates if EDM voxel is independent of
%                  boundary influence.
%
% The array EDMValid has the same size as binary. After 
% completion of the routine it will indicate:
%    EDMValid(i,j,k)==1   : pixel i,j,k more than its EDM
%                           value from boundary of binary
%                           data set. 
%    EDMValid(i,j,k)==0   : pixel i,j,k so close to the 
%                           boundary that its value could 
%                           change if values in binary on or
%                           beyond boundary of 'binary' would
%                           change.
%    EDMValid(i,j,k)==2   : pixel is in the 'other' phase 
%                           (binary==1) for which EDM is simply
%                           not computed.
%
function [D,EDMValid,idx] = calculateBwDistAndLabelBoundary(binary)
    % calculate EDM and reference point
    [D,idx] = bwdist(binary);
    % Create validity map for EDM.
    EDMValid=zeros(size(binary))-1;
    %
    % Sample sizes for reference
    nx=size(D,1);
    ny=size(D,2);
    nz=size(D,3);
    % check whether EDM values are valid or whether they could
    % be affected by boundary of sample
    for x = 1:size(D,1)
        for y = 1:size(D,2)
            for z = 1:size(D,3)
                % check consistency
                if binary(x,y,z)==0 && D(x,y,z)<1
                    disp("ERROR in maxCoveringDistanceTransform");
                    disp("binary(x,y,z)==0 && D(x,y,z)<1");
                    disp(binary(x,y,z));
                    disp(D(x,y,z));
                    disp([x,y,z]);
                    pause;
                end
                if binary(x,y,z)==1 && D(x,y,z) ~=0
                    disp("ERROR in maxCoveringDistanceTransform");
                    disp("binary(x,y,z)==1 && D(x,y,z) ~=0");
                    pause;
                end
                % test if EDM sphere overlaps with sample boundary.
                % If this is the case, then the EDM value could 
                % be affected by hypothetical "1" pixels outside of the
                % sample
                if binary(x,y,z)==1
                    EDMValid(x,y,z)=2;  % simply means not assessed
                                      % as voxel in other phase
                elseif binary(x,y,z)==0
                    if ...
                            x - 1 >= D(x,y,z) && nx - x >= D(x,y,z) && ...
                            y - 1 >= D(x,y,z) && ny - y >= D(x,y,z) && ...
                            z - 1 >= D(x,y,z) && nz - z >= D(x,y,z) 
                       EDMValid(x,y,z)=1;
                    else
                        EDMValid(x,y,z)=0;
                    end
                else
                    error("binary must be 0 or 1");
                end
                
            end
        end
    end
    [D,EDMValid,idx];
end
    
    