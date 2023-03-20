% For a binary data set of a phase 1 ("solid") and a phase 0 ("void"), 
% dilate phase 1 by dilationR. Then erode the phase back. 
%
% Both happens after placing the data set in a large void image, so 
% that the dilated body always stays surrounded by a void layer. 
%
% 
function [dilated,erodeDilated] = ErodeDilateToLabelRelevantRegions(binary,dilationR,erosionR)
    % solid phase (binary = 1) is dilated by radius dilationR
    if dilationR < 1
        error("dilationR < 1");
    end
    if erosionR < 1
        error("erosionR < 1");
    end
    
    % 'Canvas' on which this occurs is the binary data set
    % with a layer of 0's around it. This layer is a little bigger
    % than the dilation radius
    enlargeR=round(dilationR+3);
    
    %   
    enlargedBinary=zeros(size(binary)+2*enlargeR);
    enlargedBinary(...
            (enlargeR+1):(enlargeR+size(binary,1)),...
            (enlargeR+1):(enlargeR+size(binary,2)),...
            (enlargeR+1):(enlargeR+size(binary,3)))=binary;
    enlargeD = bwdist(enlargedBinary);
    BinaryWithPhase1Dilated=ones(size(enlargedBinary));
    for i=1:size(enlargeD,1)
        for j=1:size(enlargeD,2)
            for k=1:size(enlargeD,3)
                % determine points to consider 
                if enlargedBinary(i,j,k)==0
                    if enlargeD(i,j,k)>dilationR
                        BinaryWithPhase1Dilated(i,j,k)=0;
                    else
                        BinaryWithPhase1Dilated(i,j,k)=1;
                    end
                end
            end
        end
    end
    clear enlargeD;
    
    % dilated data set is BinaryWithPhase1Dilated cut back
    % to the original size of binary dataset
    dilated=BinaryWithPhase1Dilated(...
        (enlargeR+1):(enlargeR+size(binary,1)),...
        (enlargeR+1):(enlargeR+size(binary,2)),...
        (enlargeR+1):(enlargeR+size(binary,3)));

    % now erode
    erodeDilated=zeros(size(binary));
    inverseBinaryWithPhase1Dilated=1-BinaryWithPhase1Dilated;
    clear BinaryWithPhase1Dilated;
    
    % calculate distance map of inverse phase
    % (that is, distance of any 1-voxel to the nearest 0-voxel)
    % on enlarged data set. 
    % Then cut back to size of binary. 
    D2 = bwdist(inverseBinaryWithPhase1Dilated);
    D2 = D2(...
        (enlargeR+1):(enlargeR+size(binary,1)),...
        (enlargeR+1):(enlargeR+size(binary,2)),...
        (enlargeR+1):(enlargeR+size(binary,3)));
    inverseBinaryWithPhase1Dilated=inverseBinaryWithPhase1Dilated(...
        (enlargeR+1):(enlargeR+size(binary,1)),...
        (enlargeR+1):(enlargeR+size(binary,2)),...
        (enlargeR+1):(enlargeR+size(binary,3)));
    for i=1:size(D2,1)
        for j=1:size(D2,2)
            for k=1:size(D2,3)
                if inverseBinaryWithPhase1Dilated(i,j,k)==0
                    if D2(i,j,k)<erosionR
                        erodeDilated(i,j,k)=0;
                    else
                        erodeDilated(i,j,k)=1;
                    end
                else
                    erodeDilated(i,j,k)=0;
                end
            end
        end
    end
    
    [dilated,erodeDilated];
    
end
