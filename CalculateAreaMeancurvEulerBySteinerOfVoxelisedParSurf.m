
function outputStructure=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,parameters)

    % read parameter that prescribes what fraction of EDM range 
    % should be fitted to Steiner polynomial. 
    % This parameter should be set such that the dilated bodies
    % remain safely within the "reach" of the original body.
    if isfield(parameters,"fractionEuclidDistanceToFit")==0
        disp("fractionEuclidDistanceToFit is not specified. Using default value 0.5");
    else
        fractionEuclidDistanceToFit=parameters.fractionEuclidDistanceToFit;
    end

    distanceZero=0.05;   % must be kept below 1 
                         % purpose is to shift intervals, so that
                         % the integer EDM values fall within an 
                         % interval, rather than on the boundary
    if distanceZero > 1
        error("distanceZero must be below 1. Minimal EDM value possible");
    end
    incrementSteps=1.0;  % this is in pixel units. Should be kept at 1. 

    
    distancesPhase0=[];
    distancesPhase1=[];
    numberVoxelsUpToDistancePhase0=[];
    numberVoxelsUpToDistancePhase1=[];
    
    for phase=0:1
        
        % calculate EDM of phase 0 (in run 0) or phase 1 (in run 1)
        if phase == 1
            bin = 1-binary;
        else
            bin = binary;
        end
    
        % calculate EDM
        EDM = bwdist(bin);
    
        % clip to relevant section
        EDM=EDM(VoxelsToAnalyse(1,1):VoxelsToAnalyse(1,2),VoxelsToAnalyse(2,1):VoxelsToAnalyse(2,2),VoxelsToAnalyse(3,1):VoxelsToAnalyse(3,2));
        bin=bin(VoxelsToAnalyse(1,1):VoxelsToAnalyse(1,2),VoxelsToAnalyse(2,1):VoxelsToAnalyse(2,2),VoxelsToAnalyse(3,1):VoxelsToAnalyse(3,2));

        % maximal EDM value to be analysed
        maxEDMValueToAnalyse=round(max(max(max(EDM)))*fractionEuclidDistanceToFit);
        
        % store maximal EDM values 
        if phase == 1
            outputStructure.maxEDMPhase1=max(max(max(EDM)));
        else
            outputStructure.maxEDMPhase0=max(max(max(EDM)));
        end
        
        % Now calculate number of voxels in layers, plus additional
        % associated info such as average distance in layer
        EDMStatsNBins=round((maxEDMValueToAnalyse-distanceZero)/incrementSteps)+2;
        EDMStatsNVoxels=zeros(EDMStatsNBins,1);
        EDMStatsNVoxelsCumulative=zeros(EDMStatsNBins,1);
        EDMStatsAverageEDM=zeros(EDMStatsNBins,1);
        EDMStatsMinEDM=zeros(EDMStatsNBins,1)+maxEDMValueToAnalyse+1; 
        EDMStatsMaxEDM=zeros(EDMStatsNBins,1)-1; 
        EDMStatsNVoxels1Phase=0;
        EDMStatsNVoxels0Phase=0;
        EDMStatsNVoxels0PhaseGreatermaxEDMValueToAnalyse=0;  
        for xx = 1:size(bin,1)
            for yy = 1:size(bin,2)
                for zz = 1:size(bin,3)
                    if bin(xx,yy,zz) == 1
                        EDMStatsNVoxels1Phase=EDMStatsNVoxels1Phase+1;
                    else
                        % count basic statistic
                        EDMStatsNVoxels0Phase=EDMStatsNVoxels0Phase+1;
                        
                        % this means 0 phase
                        EDMValue=EDM(xx,yy,zz);
                        if EDMValue < 0
                            error("Negative EDM values. Not possible.");
                        elseif EDMValue < 1
                            error("Error, EDM(xx,yy,zz) < 1 in phase to be analysed. 1 is lowest value ");
                        end
                        % calculate which bin it falls into
                        if EDMValue-distanceZero < 0
                            error("EDM(xx,yy,zz)-distanceZero < 0");
                        end
                        
                        % only count those bins where EDM is 
                        % within percentage given by parameter
                        if EDMValue > maxEDMValueToAnalyse
                            EDMStatsNVoxels0PhaseGreatermaxEDMValueToAnalyse=EDMStatsNVoxels0PhaseGreatermaxEDMValueToAnalyse+1;
                        else
                            binIDX=floor((EDMValue-distanceZero)/incrementSteps)+1;
                            % the +1 reflects matlab arrays beginning with 1. 
                        
                            % error test only
                            if binIDX < 1
                                disp(binIDX);
                                error("binIDX < 1. Should never happen.")
                            end
                            if binIDX > length(EDMStatsNVoxels)
                                disp(binIDX);
                                disp(length(EDMStatsNVoxels));
                                error("binIDX > length(EDMStatsNVoxels");
                            end
                            
                            % update statistics
                            EDMStatsNVoxels(binIDX)=EDMStatsNVoxels(binIDX)+1;
                            EDMStatsAverageEDM(binIDX)=EDMStatsAverageEDM(binIDX)+EDMValue;
                            if EDMValue >= EDMStatsMaxEDM(binIDX)
                                EDMStatsMaxEDM(binIDX)=EDMValue;
                            end
                            if EDMValue <= EDMStatsMinEDM(binIDX)
                                EDMStatsMinEDM(binIDX)=EDMValue;
                            end
                        end
                    end
                end
            end
        end
      
        % calculate averages
        for i=1:length(EDMStatsNVoxels)
            if EDMStatsNVoxels(i) > 0
                EDMStatsAverageEDM(i)=EDMStatsAverageEDM(i)/EDMStatsNVoxels(i);
            end
        end
        
        % calculate cumulative number of voxels up to this layer
        for i=1:length(EDMStatsNVoxels)
            for j=1:i
                EDMStatsNVoxelsCumulative(i) = EDMStatsNVoxelsCumulative(i)+EDMStatsNVoxels(j);
            end  
        end
        
        % define value that we want to use as distance value for 
        % Steiner formula Vol(dist). There are a few possible choices, but
        % the following works best
        EDMStatsUpperDistanceValue=EDMStatsMaxEDM;
        %EDMStatsUpperDistanceValue=EDMStatsAverageEDM;
        
        % remove empty layers
        zeroset=find(EDMStatsNVoxels~=0);
        EDMStatsNVoxels=EDMStatsNVoxels(zeroset);
        EDMStatsUpperDistanceValue=EDMStatsUpperDistanceValue(zeroset);
        EDMStatsNVoxelsCumulative=EDMStatsNVoxelsCumulative(zeroset);
        EDMStatsMaxEDM=EDMStatsMaxEDM(zeroset);
        EDMStatsMinEDM=EDMStatsMinEDM(zeroset);
        EDMStatsAverageEDM=EDMStatsAverageEDM(zeroset);
            
        % append to phase1/phase0 vector
        if phase == 0
            distancesPhase0=EDMStatsUpperDistanceValue;
            numberVoxelsUpToDistancePhase0=EDMStatsNVoxelsCumulative;
            %distances=-flip([0;EDMStepsAverageEDMValue]);
            %numberVoxelsUpToDistance=-flip([0;EDMStepsNVoxelsInt]);
            %zeroindex=length(numberVoxelsUpToDistance);
        else
            distancesPhase1=EDMStatsUpperDistanceValue;
            numberVoxelsUpToDistancePhase1=EDMStatsNVoxelsCumulative;
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
    if isfield(parameters,"removeFirstLayerFromFit")==1
        if parameters.removeFirstLayerFromFit == 1
            l=length(numberVoxelsUpToDistancePhase1);
            if l >= 3 
                numberVoxelsUpToDistancePhase1=numberVoxelsUpToDistancePhase1(2:l);
                distancesPhase1=distancesPhase1(2:l);
            end
            l=length(numberVoxelsUpToDistancePhase0);
            if l >= 3
                numberVoxelsUpToDistancePhase0=numberVoxelsUpToDistancePhase0(2:l);
                distancesPhase0=distancesPhase0(2:l);
            end
        end
    end
    
        
    
    % Apply a sub-voxel-size shift to the data. Because this shift is
    % subvoxel, it vanishes in the limit where the voxel size goes to 
    % zero. 
    % Applying this shift does not affect the convergence properties.
    % Empirically, from analysis of known models it improves estimates
    % obtained for coarse discretisations.
    % The variable magnitude of the shift (depending on the coarseness
    % of the discretisation relative to the pore size = max EDM of the
    % sample) was determined based on structures where exact 
    % values are known. 
    % See file PlotBestShiftValuesForSomeKnownModels.pdf
    if isfield(parameters,"DoNotApplyShiftToDistanceValueToImproveConvergence") == 0
        shiftPhase1=0.12+0.13*exp(-outputStructure.maxEDMPhase1/30);
        distancesPhase1=distancesPhase1-shiftPhase1;
        shiftPhase0=0.12+0.13*exp(-outputStructure.maxEDMPhase0/30);
        distancesPhase0=distancesPhase0-shiftPhase0;
    end
    
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
          
       % extract parameter how to deal with Euler
       if isfield(parameters,"EulerMethod") == 0
           disp("parameter EulerMethod not specified. Using default fitAreaMeanCurvatureAndEuler");
       elseif parameters.EulerMethod == "fitAreaAndMeanCurvatureUseProvidedEuler"
                  disp("For method fitAreaAndMeanCurvatureUseProvidedEuler it is ");
               disp("necessary to provide the Euler number value as an argument ");
         % check that value to be used for Euler is provided. 
           if isfield(parameters, "ProvidedEulerValue") == 0
              disp("using parameters.ProvidedEulerValue. Didn't happen. Error.");
               error("Need to provide parameters.ProvidedEulerValue");
           end
       elseif parameters.EulerMethod == "fitAreaMeanCurvatureAndEuler"
            disp("Using fitAreaMeanCurvatureAndEuler method ");
       elseif parameters.EulerMethod == "fitHorizontalShiftToKnownExactDataForCodeDevelopment"
            disp("Using fitHorizontalShiftToKnownExactDataForCodeDevelopment method ");
            disp("      This method is for code development only");
            if isfield(parameters, "ProvidedEulerValue") == 0
              error("Need to provide parameters.ProvidedEulerValue");
            end
            if isfield(parameters, "ProvidedAreaValue") == 0
              error("Need to provide parameters.ProvidedAreaValue");
            end
            if isfield(parameters, "ProvidedMeanCurvatureValue") == 0
              error("Need to provide parameters.ProvidedMeanCurvatureValue");
            end
       else
           error("Name of Euler method not recognised. Error.");           
       end
      
       % calculate area, integral mean curvature and (depending on
       % specification) also Euler number
       if parameters.EulerMethod=="fitAreaMeanCurvatureAndEuler"
        
        
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
           
       elseif parameters.EulerMethod == "fitAreaAndMeanCurvatureUseProvidedEuler"
           
           Euler = parameters.ProvidedEulerValue; % specific for Gyroid
           
           % subtract Euler d^3 term
           volUpToDistanceMinusEulerTerm=volUpToDistance-Euler/3*2*pi*distances.^3;
           
            % use notation as in conventional linear regression (XIJ, YI, BETA)
           XXIJ=[distances.^1, distances.^2];
           YY=volUpToDistanceMinusEulerTerm;
           beta=inv(XXIJ'*XXIJ)*XXIJ'*YY;
           % convert back to area, chi, <H>
           area=beta(1);
           meanc=beta(2)/area;
      
           % store values in output structure
           outputStructure.area=area;
           outputStructure.meanc=meanc;
           outputStructure.Euler=Euler;
           outputStructure.fitSuccessful=1;
         
        elseif parameters.EulerMethod == "fitHorizontalShiftToKnownExactDataForCodeDevelopment"
         
            % This section is for code development only.
            % Determine value of a that minimises square deviations of
            % V(d-a) =
            % area*((d-a)+MeanCurvture*(d-a)^2)+1/2*Euler*2*pi*(d-a)^3
            % from measured data. 

            % exact values
            area=parameters.ProvidedAreaValue;
            meanc=parameters.ProvidedMeanCurvatureValue;
            Euler=parameters.ProvidedEulerValue;
           
            % try a=0:0.01:1 (in voxel units) and find 
            % best fit
            minSquareDev = 10e50;
            minA=[-1,-1];
            optimalShiftedDistance=0;
            for ap=-1:0.01:1
                for am=-1:0.01:1
                    reddistances=distances;
                    % apply shift: positive/negative depending on sign
                    for k=1:length(reddistances)
                        if abs(reddistances(k)/PixelSize) > 0.01  % leave point at d=0 unshifted
                            if reddistances(k)>0
                                reddistances(k)= reddistances(k)-ap*PixelSize;
                            elseif reddistances(k)<0
                                reddistances(k)=reddistances(k)+am*PixelSize;
                            end
                        end
                    end
             
                    YY=area*(reddistances+meanc*reddistances.^2)+1/3*Euler*2*pi*reddistances.^3;
                    squareDev=(sum((YY-volUpToDistance).^2))/length(YY);
                    if squareDev < minSquareDev
                        minSquareDev = squareDev;
                        minA=[am,ap];
                        optimalShiftedDistance=reddistances;
                    end
                end
            end
           outputStructure.shiftAInPixels=minA;
           outputStructure.pixelSize=PixelSize;
           outputStructure.shiftedDistances=optimalShiftedDistance;
           outputStructure.volumesUpToDistance=volUpToDistance;
           outputStructure.area=area;
           outputStructure.meanc=meanc;
           outputStructure.Euler=Euler;
           outputStructure.fitSuccessful=1;
           
       end
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


