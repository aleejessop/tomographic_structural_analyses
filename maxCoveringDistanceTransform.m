
% RETURN VALUES:
% 
% CRTValues: Field that contains the CRT values for each voxel
%     > 1    : CRT values for voxels in phase 0 of input data set
%     == -2  : CRT not computed because voxel is in phase 1 of input data
%     == -1  : CRT not computed or not correctly computed for other 
%              reasons, see statusCRTComputation.
% 
% statusCRTComputation: field to indicate status of CRT computation
%     == 1   : successfully computed without boundary effects
%     == 5   : CRT computed, but covering sphere exceeded voxelised
%              data set (boundary effects likely)
%     == -10 : CRT not computed; voxel in other phase
%     == -25 : CRT not computed; outside of specified region
%     == -20 : CRT not computed or wrong. The EDM sphere at
%              this pixel not considered because its value exceeds
%              maximal distance set by parameter.
%     == -15 : CRT not correctly computed because a large EDM sphere
%              that would cover this voxel is not considered because
%              the parameter maximalCRT value was set.
%              
% PARAMETERS
%
% parameters.maximalCRTValueToCompute (optional): Consider only 
%         EDM spheres up to this radius for computing the CRT. 
%              
%
function [CRTValues,image1,D,statusCRTComputation] = maxCoveringDistanceTransform(binary,parameters)

    % Determine subset of data set for which to compute CRT 
    %      1. Default: entire data set
    %      2. If 'calculateCRTOnlyForRectangularSubset': use specified
    %         subset
    %
    VoxelsToAnalyseMinX=1;
    VoxelsToAnalyseMaxX=size(binary,1);
    VoxelsToAnalyseMinY=1;
    VoxelsToAnalyseMaxY=size(binary,2);
    VoxelsToAnalyseMinZ=1;
    VoxelsToAnalyseMaxZ=size(binary,3);
    if isfield(parameters,"calculateCRTOnlyForRectangularSubset")==1 
        
        if parameters.calculateCRTOnlyForRectangularSubset==1
            if isfield(parameters,"rectangularSubsetForCRTCalculation")==0
                error("Missing parameter: rectangularSubsetForCRTCalculation");
            else
                VoxelsToAnalyse=parameters.rectangularSubsetForCRTCalculation;
                if size(VoxelsToAnalyse,1) ~= 3 || size(VoxelsToAnalyse,2) ~= 2
                    error("size(VoxelsToAnalyse,1) ~= 3 || size(VoxelsToAnalyse,2) ~= 2");
                end
                if VoxelsToAnalyse(1,1)<1 || VoxelsToAnalyse(1,2)>size(binary,1)
                    error("VoxelsToAnalyse(1,1)<1 || VoxelsToAnalyse(1,2)>size(binary,1)");
                end
                if VoxelsToAnalyse(2,1)<1 || VoxelsToAnalyse(2,2)>size(binary,2)
                    error("VoxelsToAnalyse(2,1)<1 || VoxelsToAnalyse(2,2)>size(binary,2)");
                end
                if VoxelsToAnalyse(3,1)<1 || VoxelsToAnalyse(3,2)>size(binary,3)
                    error("VoxelsToAnalyse(3,1)<1 || VoxelsToAnalyse(3,2)>size(binary,3)");
                end
                VoxelsToAnalyseMinX=VoxelsToAnalyse(1,1);
                VoxelsToAnalyseMaxX=VoxelsToAnalyse(1,2);
                VoxelsToAnalyseMinY=VoxelsToAnalyse(2,1);
                VoxelsToAnalyseMaxY=VoxelsToAnalyse(2,2);
                VoxelsToAnalyseMinZ=VoxelsToAnalyse(3,1);
                VoxelsToAnalyseMaxZ=VoxelsToAnalyse(3,2);
            end
        end
    end
    CenterX=round((VoxelsToAnalyseMinX+VoxelsToAnalyseMaxX)/2);
    CenterY=round((VoxelsToAnalyseMinY+VoxelsToAnalyseMaxY)/2);
    CenterZ=round((VoxelsToAnalyseMinZ+VoxelsToAnalyseMaxZ)/2);
    
    % create a field in which status of CRT computation will be
    % stored
    % 0  = has not been determined yet.
    statusCRTComputation=zeros(size(binary));
    
    
    % calculate EDM
    % This is in essence just D = bwdist(binary), except that it 
    % also calculates information about the EDM spheres that 
    % are so close to sample boundary that they may be affected
    % by boundary influences.
    disp("CRT : calculating EDM");
    [D,EDMValid,idx] = calculateBwDistAndLabelBoundary(binary);
    
    %data=EDMValid; sliceview; pause;
    
    % check if CRT should be only computed up to a maximal value
    maximalCRTValueToCompute=max(size(binary)); % EDM < system size (practically infty)
    if isfield(parameters,"maximalCRTValueToCompute")
        % check parameter valid
        if parameters.maximalCRTValueToCompute <= 1
            error("parameters.maximalCRTValueToCompute can't be less than 1");
            % EDM can't be less than 1. hence this limit.
        end
        % extract parameter
        maximalCRTValueToCompute=parameters.maximalCRTValueToCompute;
        fprintf("CRT : Will compute CRT only up to given value (%f) in VOXEL units\n", ...
           maximalCRTValueToCompute );
        % basically perform an erosion on the 1 phase (the phase for which
        % EDM/CRT are NOT computed) by the desired distance.
        [dilated,dilatedThenEroded]=ErodeDilateToLabelRelevantRegions(binary,maximalCRTValueToCompute,maximalCRTValueToCompute+1);
        for x = 1:size(binary,1)
            for y = 1:size(binary,2)
                for z = 1:size(binary,3)
                    if binary(x,y,z)==0 && dilated(x,y,z) == 0
                        statusCRTComputation(x,y,z)=-20;  % CRT wrong, because 
                                                          %EDM sphere
                                                          % at this voxel
                                                          % not considered
                    end
                    if ... 
                            binary(x,y,z)==0 && ...
                            dilatedThenEroded(x,y,z) == 0 && ...
                            dilated(x,y,z) ~= 0
                        statusCRTComputation(x,y,z)=-15; % CRT wrong, because
                                                         % EDM spheres that
                                                         % overlap with
                                                         % this voxel were
                                                         % not considered.
                    end
                end
            end
        end
        clear dilated, dilatedThenEroded;
    end
    
        
    %data=statusCRTComputation;
    %sliceview
    %pause
    
    disp("CRT : Creating linear list of all EDM values");
    % entries: [DV,DistCenter,X,Y,Z,EDMStatus] for 
    % each voxel
    % this list will contain 1 entry for every voxel in void phase
    % (binary=0) which has status != -20. Adjust number of elements
    % to that number
    tmp=((binary==0)&(statusCRTComputation~=-20));
    LengthDistanceListToSort=length(find(tmp==1));
    DistListArray=zeros(LengthDistanceListToSort,6);
    clear tmp;
    
    %DV = [];
    %X = [];
    %Y = [];
    %Z = [];
    %EDMStatus = [];
    k = 1;
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
                    error("binary(x,y,z)==0 && D(x,y,z)<1");
                end
                if binary(x,y,z)==1 && D(x,y,z) ~=0
                    disp("ERROR in maxCoveringDistanceTransform");
                    disp("binary(x,y,z)==1 && D(x,y,z) ~=0");
                    error("binary(x,y,z)==1 && D(x,y,z) ~=0");
                end
                if binary(x,y,z)==1 && EDMValid(x,y,z)~=2
                    error("binary(x,y,z)==1 && EDMValid(x,y,z)~=2");
                end
                % include pore space voxels in list
                if binary(x,y,z)==0
                    if statusCRTComputation(x,y,z)~=-20   % max CRT value criterion
                        % distance to center
                        dist2ToCenter=(x-CenterX)^2+(y-CenterY)^2+(z-CenterZ)^2;
                        % old version (as lists)
                        %DV(k) = D(x,y,z);
                        %X(k) = x;
                        %Y(k) = y;
                        %Z(k) = z;
                        %EDMStatus(k) = EDMValid(x,y,z);
                        % new version, as matrix
                        if k>LengthDistanceListToSort
                            error("k>LengthDistanceListToSort");
                        end
                        DistListArray(k,1)=D(x,y,z);
                        DistListArray(k,2)=dist2ToCenter;
                        DistListArray(k,3)=x;
                        DistListArray(k,4)=y;
                        DistListArray(k,5)=z;
                        DistListArray(k,6)=EDMValid(x,y,z);
                        %DistListArray(k,:)=[D(x,y,z),dist2ToCenter,x,y,z,EDMValid(x,y,z)];
                        k = k+1;
                    end
                else
                    statusCRTComputation(x,y,z)=-10; % not compute CRT
                                                     % because voxel phase
                                                     % in phase 1. 
                    CRTValues(x,y,z)=-2; % indicating voxel is in phase 1 
                end
            end
        end
    end

    % old sort
    disp("CRT : sorting EDM values linearly");
    %[DVS, index] = sort(DV,'descend');
    %Xsort = X(index);
    %Ysort = Y(index);
    %Zsort = Z(index);
    %EDMStatusSort = EDMStatus(index);
    % Essential criterion is that EDM spheres are sorted in decreasing
    % size.
    % Sorting them also by their EDM status (with those not affected
    % by boundary conditions of sample) and by distance to sample perimeter
    % means that boundary effects are minimised in CRT data.
    DistListArray=sortrows(DistListArray,[1,6,2],{'descend','descend','ascend'});
    %size(DistListArray)
    %DistListArray(1:100,:)
    
    
    
    % CRTValues to contain the values of the CRT
    % image1 to contain the centers of EDM spheres that contribute to 
    %        CRT (related to medial axis)
    CRTValues = zeros(size(D,1),size(D,2),size(D,3))-1; % pre-set to impossible value                                                     
    image1 = zeros(size(D,1),size(D,2),size(D,3));

    % progress output only (not functionality)
    fprintf("CRT: Iteration (Ballsize): ")
    maxBallSize=DistListArray(1,1);
    printProgressMessageAt=maxBallSize-1;
    % loop
    %for k = 1:length(DVS)
    for k=1:length(DistListArray)
        
        % new
        kD=DistListArray(k,1);
        kX=DistListArray(k,3);
        kY=DistListArray(k,4);
        kZ=DistListArray(k,5);
        EDMStat=DistListArray(k,6);
        
        % old sort
        %kD = DVS(k);
        %kX = Xsort(k);
        %kY = Ysort(k);
        %kZ = Zsort(k);
        %EDMStat = EDMStatusSort(k);
        
        % progress output
        if kD<=printProgressMessageAt || k==1
            fprintf("%.1f, ", kD);
            printProgressMessageAt=printProgressMessageAt-0.05*maxBallSize;
        end
    
        
        % determine the range of voxels to iterate - restricted to subset 
        % for analysis
        kD2=kD^2;
        for xx = max(kX-(round(kD)+1),VoxelsToAnalyseMinX):min(kX+(round(kD)+1),VoxelsToAnalyseMaxX)
            for yy = max(kY-(round(kD)+1),VoxelsToAnalyseMinY):min(kY+(round(kD)+1),VoxelsToAnalyseMaxY)
                for zz = max(kZ-(round(kD)+1),VoxelsToAnalyseMinZ):min(kZ+(round(kD)+1),VoxelsToAnalyseMaxZ)
                    if binary(xx,yy,zz) == 0
                        if CRTValues(xx,yy,zz) == -1 && (xx-kX)^2 + (yy-kY)^2 + (zz-kZ)^2 <= kD2
                            % consistency check
                            if ... 
                                    statusCRTComputation(xx,yy,zz) ~= 0 && ...
                                    statusCRTComputation(xx,yy,zz) ~= -20 && ...
                                    statusCRTComputation(xx,yy,zz) ~= -15
                                error("statusCRTComputation(xx,yy,zz) ~= 0 or -20 or -15");
                            end
                            % set CRT value
                            if ... 
                                    statusCRTComputation(xx,yy,zz) ~= -20 && ...
                                    statusCRTComputation(xx,yy,zz) ~= -15
                                CRTValues(xx,yy,zz) = kD;
                            end
                            % set center of EDM sphere (kind of medial
                            % axis)
                            image1(kX,kY,kZ) = 1;
                            % set status of voxel
                            if EDMStat == 2
                                error("EDMStat == 2");
                            end
                            if EDMStat == 1
                                if ...
                                        statusCRTComputation(xx,yy,zz) ~= -20 && ...
                                        statusCRTComputation(xx,yy,zz) ~= -15
                                   statusCRTComputation(xx,yy,zz) = 1; % CRT successfully computed.
                                                                        % no
                                                                        % boundary
                                                                        % effect
                                end
                            elseif EDMStat == 0
                                if ...
                                        statusCRTComputation(xx,yy,zz) ~= -20 && ...
                                        statusCRTComputation(xx,yy,zz) ~= -15
                                    statusCRTComputation(xx,yy,zz) = 5; 
                                    % CRT value computed, but EDM sphere
                                    % could be subject to boundary effects
                                end
                            end
                        end 
                    end
                end
            end
        end
    end
    fprintf("done\n");
    
    % Set status of voxels that fall outside of windown of observation
    % -5 = CRT not to be computed as outside of specified box 
    for x = 1:size(statusCRTComputation,1)
        for y = 1:size(statusCRTComputation,2)
            for z = 1:size(statusCRTComputation,3)
                if ...
                        x < VoxelsToAnalyseMinX || x > VoxelsToAnalyseMaxX || ...
                        y < VoxelsToAnalyseMinY || y > VoxelsToAnalyseMaxY || ...
                        z < VoxelsToAnalyseMinZ || z > VoxelsToAnalyseMaxZ 
                    statusCRTComputation(x,y,z)=-25;
                end
            end
        end
    end
    
    
    % statistics
    possibleStatusValues=[5,1,0,-10,-15,-20,-25];
    nTotal=(size(statusCRTComputation,1)*size(statusCRTComputation,2)*size(statusCRTComputation,3));
    sum=0;
    for status=possibleStatusValues
        nStatus=length(find(statusCRTComputation==status));
        fprintf("CRT: voxels with status %d: %d (%.1f percent) \n", status, nStatus, nStatus/nTotal*100);
        sum=sum+nStatus;
    end
    if sum ~= nTotal
        error("sum ~= nTotal");
    end
      
    % consistency check on CRT
    for i=1:size(CRTValues,1)
        for j=1:size(CRTValues,2)
            for k=1:size(CRTValues,3)
                if CRTValues(i,j,k) ~= -2 && CRTValues(i,j,k) ~= -1 && CRTValues(i,j,k) < 1
                    error("CRT values need to be -2, -1 or positive and >= 1");
                end
            end
        end
    end
    
    [CRTValues,image1,D,statusCRTComputation];
end
