
function outputStructure=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,fractionEuclidDistanceToFit)

    incrementSteps=0.5;  % this is in pixel units

    distancesPhase0=[];
    distancesPhase1=[];
    numberVoxelsUpToDistancePhase0=[];
    numberVoxelsUpToDistancePhase1=[];

 
    for phase=0:1
        % calculate EDM of phase 1 (in run 0) or phase 0 (in run 1)
        if phase == 0
            bin = 1-binary;
        else
            bin = binary;
        end
    
        % calculate EDM
        EDM = bwdist(bin);
    
        % clip to relevant section
        EDM=EDM(VoxelsToAnalyse(1,1):VoxelsToAnalyse(1,2),VoxelsToAnalyse(2,1):VoxelsToAnalyse(2,2),VoxelsToAnalyse(3,1):VoxelsToAnalyse(3,2));
        bin=bin(VoxelsToAnalyse(1,1):VoxelsToAnalyse(1,2),VoxelsToAnalyse(2,1):VoxelsToAnalyse(2,2),VoxelsToAnalyse(3,1):VoxelsToAnalyse(3,2));

        % determine steps for analysis
        maxEDM=round(max(max(max(EDM)))*fractionEuclidDistanceToFit);
        if 1+incrementSteps >= maxEDM
            break;
        end
        EDMTargetSteps=(1+incrementSteps):incrementSteps:maxEDM;
        EDMStepsNVoxelsDiff=zeros(length(EDMTargetSteps),1);
        EDMStepsNVoxelsInt=zeros(length(EDMTargetSteps),1);
        EDMStepsAverageEDMValue=zeros(length(EDMTargetSteps),1);
        lowerBound=0;
        for i=1:length(EDMTargetSteps)
            upperBound=EDMTargetSteps(i);
            %disp([lowerBound,upperBound]);
            nVoxels=0;
            meanEDM=0;
            maxEDM=0;
            for xx = 1:size(bin,1)
                for yy = 1:size(bin,2)
                    for zz = 1:size(bin,3)
                        if bin(xx,yy,zz) == 0
                            if EDM(xx,yy,zz) < 1
                                disp("Error, EDM(xx,yy,zz) < 1 ");
                                pause;
                            end
                            if( EDM(xx,yy,zz) > lowerBound && EDM(xx,yy,zz) <= upperBound )
                                nVoxels=nVoxels+1;
                                meanEDM=meanEDM+EDM(xx,yy,zz);
                                if EDM(xx,yy,zz) > maxEDM
                                    maxEDM=EDM(xx,yy,zz);
                                end
                            end
                        end
                    end
                end
            end
            EDMStepsNVoxelsDiff(i)=nVoxels;
            %EDMStepsAverageEDMValue(i)=maxEDM;    %meanEDM/nVoxels;
            EDMStepsAverageEDMValue(i)=meanEDM/nVoxels;
            %EDMStepsAverageEDMValue(i)=(meanEDM/nVoxels+lowerBound)/2;
            for j=1:i
                EDMStepsNVoxelsInt(i) = EDMStepsNVoxelsInt(i)+EDMStepsNVoxelsDiff(j);
            end   
            lowerBound=upperBound;
        end
        % remove all that are zero
        zeroset=find(EDMStepsNVoxelsDiff~=0);
        EDMStepsNVoxelsDiff=EDMStepsNVoxelsDiff(zeroset);
        EDMStepsAverageEDMValue=EDMStepsAverageEDMValue(zeroset);
        EDMStepsNVoxelsInt=EDMStepsNVoxelsInt(zeroset);
        
           
        % append to phase1/phase0 vector
        if phase == 0
            distancesPhase0=EDMStepsAverageEDMValue;
            numberVoxelsUpToDistancePhase0=EDMStepsNVoxelsInt;
            %distances=-flip([0;EDMStepsAverageEDMValue]);
            %numberVoxelsUpToDistance=-flip([0;EDMStepsNVoxelsInt]);
            %zeroindex=length(numberVoxelsUpToDistance);
        else
            distancesPhase1=EDMStepsAverageEDMValue;
            numberVoxelsUpToDistancePhase1=EDMStepsNVoxelsInt;
            %distances=-flip([0;EDMStepsAverageEDMValue]);
            %numberVoxelsUpToDistance=-flip([0;EDMStepsNVoxelsInt]);
            %zeroindex=length(numberVoxelsUpToDistance);distances=[distances;0;EDMStepsAverageEDMValue];
            %numberVoxelsUpToDistance=[numberVoxelsUpToDistance;0;EDMStepsNVoxelsInt];
        end
    end
    
    % remove any shells with zero volume (there shouldn't be any though)
    indices=find(numberVoxelsUpToDistancePhase1~=0);
    numberVoxelsUpToDistancePhase1=numberVoxelsUpToDistancePhase1(indices);
    distancesPhase1=distancesPhase1(indices);
    indices=find(numberVoxelsUpToDistancePhase0~=0);
    numberVoxelsUpToDistancePhase0=numberVoxelsUpToDistancePhase0(indices);
    distancesPhase0=distancesPhase0(indices);
    
    % remove first measurement point for each phase (this point seems
    % to have the greatest discretisation errors)
    l=length(numberVoxelsUpToDistancePhase1);
    if l >= 3
        %numberVoxelsUpToDistancePhase1=numberVoxelsUpToDistancePhase1(2:l);
        %distancesPhase1=distancesPhase1(2:l);
    end
    l=length(numberVoxelsUpToDistancePhase0);
    if l >= 3
        %numberVoxelsUpToDistancePhase0=numberVoxelsUpToDistancePhase0(2:l);
        %distancesPhase0=distancesPhase0(2:l);
    end
    
        
    
    % Probably wrong:
    % The EDM black-to-white or white-to-black overestimates
    % the distance to the interface by up to half a direction through
    % a voxel, ie by distances between 1/2 and sqrt{3}/2 voxels.
    % This is here corrected.
    %distancesPhase1=distancesPhase1-0.5;
    %distancesPhase0=distancesPhase0-0.5;
    
    % convert to real units
    distancesPhase1=distancesPhase1*PixelSize;
    volUpToDistancePhase1=numberVoxelsUpToDistancePhase1*PixelSize^3;
    distancesPhase0=distancesPhase0*PixelSize;
    volUpToDistancePhase0=numberVoxelsUpToDistancePhase0*PixelSize^3;
       
    % create combined distance and vol vectors 
    distances=[-flip(distancesPhase0);0;distancesPhase1];
    volUpToDistance=[-flip(volUpToDistancePhase0);0;volUpToDistancePhase1];
    
    
     
    % fit polynomial  (but only if at least two data points incl 0)
    if length(numberVoxelsUpToDistancePhase1) >= 1 && length(numberVoxelsUpToDistancePhase0) >= 1
          
       % now fit (here done using linear regression)
       %
       %      alternatively (and obsolete) through 'fit' with 
       %      steinerFct = 'area*(x+meanC*x^2)+1/3*EulerEst*2*pi*x^3';
       %      % area estimate, set curvatures to 0
       %      areaEst=(volUpToDistancePhase1(1)-volUpToDistancePhase1(1))/(distancesPhase1(1)-distancesPhase0(1));
       %      startValuesAreaMeanCGaussC=[areaEst,0,0];
       %      fitvalues=fit(distances,volUpToDistance,steinerFct,'Start',startValuesAreaMeanCGaussC)
       %
       
       % use notation as in conventional linear regression (XIJ, YI, BETA)
       XXIJ=[distances.^1, distances.^2, distances.^3];
       YY=volUpToDistance;
       beta=inv(XXIJ'*XXIJ)*XXIJ'*YY;
       % convert back to area, chi, <H>
       area=beta(1);
       meanc=beta(2)/area;
       Euler=beta(3)*3/2/pi;
      
       % store values in output structure
       outputStructure.area=area;
       outputStructure.meanc=meanc;
       outputStructure.Euler=Euler;
       outputStructure.fitSuccessful=1;
    else
        outputStructure.area= -666;
        outputStructure.meanc= -666;
        outputStructure.Euler=-666;
        outputStructure.fitSuccessful=0;
    end
    %
    outputStructure.layerD=distances;
    outputStructure.layerV=volUpToDistance;
end


