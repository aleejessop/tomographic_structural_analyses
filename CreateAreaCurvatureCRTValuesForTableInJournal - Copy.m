clear all;
close all;




% p surface
N=90;
PixelSize=1/N;
surfaceType="p";
KnownEulerNumber=-4;
thresholds=[0,0.665,0.2802];
meanCurvaturesExact=[0,-1.13,-0.410];
EulerNumberExact=[KnownEulerNumber,KnownEulerNumber,KnownEulerNumber];
NOrientations=4; % this number must be < 4 (and check code)
surfaceAreaDiscretisation=0.02;%0.005;
surfaceArea=zeros(1,length(thresholds));
for i=1:length(surfaceArea)
    surfaceArea(i)=SurfaceAreaTranslationalUnitNodalSurface(surfaceType,thresholds(i),surfaceAreaDiscretisation);
end

% d surface
N=170 ;
PixelSize=1/N;
surfaceType="d";
KnownEulerNumber=-16;
thresholds=[0,0.196,0.462];
meanCurvaturesExact=[0,-0.57,-1.47];
EulerNumberExact=[KnownEulerNumber,KnownEulerNumber,KnownEulerNumber];
NOrientations=4; % this number must be < 4 (and check code)
surfaceAreaDiscretisation=0.02;%0.005;
surfaceArea=zeros(1,length(thresholds));
for i=1:length(surfaceArea)
    surfaceArea(i)=SurfaceAreaTranslationalUnitNodalSurface(surfaceType,thresholds(i),surfaceAreaDiscretisation);
end

% g surface
N=30 ;
PixelSize=1/N;
surfaceType="g";
KnownEulerNumber=-8;
thresholds=[0,0.249,0.587];
meanCurvaturesExact=[0,-0.44,-1.13];
EulerNumberExact=[KnownEulerNumber,KnownEulerNumber,KnownEulerNumber];
NOrientations=4; % this number must be < 4 (and check code)
surfaceAreaDiscretisation=0.02;%0.005;
surfaceArea=zeros(1,length(thresholds));
for i=1:length(surfaceArea)
    surfaceArea(i)=SurfaceAreaTranslationalUnitNodalSurface(surfaceType,thresholds(i),surfaceAreaDiscretisation);
end

% a single sphere
N=23;
PixelSize=1/N;
surfaceType="scs";
radius=[0.1,0.15,0.2,0.25,0.3];
thresholds=radius.^2;
KnownEulerNumber=2;
EulerNumberExact=ones(1,length(radius))*KnownEulerNumber;
meanCurvaturesExact=1./radius;
surfaceArea=4*pi*radius.^2;
NOrientations=1; % only 1 orientation possible for sphere


% make sure that fields entered above are all of the right format
% (they need to be row vectors - that's the assumption below)

% parameters for Steiner routine
parameters.fractionEuclidDistanceToFit=0.50;

% alternate this for the two different sections in 
% journal table:
% When Euler specified     : computer errors only for A and H
% When Euler not specified : computer errors for A, H and Chi
%
%parameters.EulerMethod="fitAreaAndMeanCurvatureUseProvidedEuler";
parameters.EulerMethod="fitAreaMeanCurvatureAndEuler";


% what to calculate
calculateAreaCurvature=1;
calculateEDMPhase1=1;
calculateEDMPhase2=1;
calculateCRTPhase1=1;
calculateCRTPhase2=1;

% warnings
warningsToAlertUserAtEnd=["The following are warnings collected in run"];

% export graphics? (only works in Matlab 2022+)
ExportGraphics=0;

% store properties as function of threshold, and separately for each
% orientation
SolidPhaseMaxEuclidDistance=zeros(length(thresholds),NOrientations);
VoidPhaseMaxEuclidDistance=zeros(length(thresholds),NOrientations);
SolidVolumeFraction=zeros(length(thresholds),NOrientations);
SolidPhaseMeanValueMaxCoverRadiusTrafo=zeros(length(thresholds),NOrientations);
SolidPhaseMeanSquareCoverRadiusTrafo=zeros(length(thresholds),NOrientations);
SolidPhaseStandardDevMaxCoverRadiusTrafo=zeros(length(thresholds),NOrientations);
VoidPhaseMeanValueMaxCoverRadiusTrafo=zeros(length(thresholds),NOrientations);
VoidPhaseMeanSquareCoverRadiusTrafo=zeros(length(thresholds),NOrientations);
VoidPhaseStandardDevMaxCoverRadiusTrafo=zeros(length(thresholds),NOrientations);
SurfaceAreaFromSteinerFit=zeros(length(thresholds),NOrientations);
AverageMeanCurvatureFromSteinerFit=zeros(length(thresholds),NOrientations);
AverageEulerIndexFromSteinerFit=zeros(length(thresholds),NOrientations);


% just to keep track of how long it takes
startTime = datetime(now,'ConvertFrom','datenum');

% filebase to store results
fileBase=join(["Properties_",surfaceType,"_NOri",string(NOrientations),"_N",string(N)],"");
disp("Will print results to the following file:");
disp(fileBase);

for i=1:length(thresholds)
    Threshold=thresholds(i);
    
    for j=1:NOrientations

        % create 3x3x3 of D surface in 100 orientation
        if j==1
            NVoxels=round([N,N,N]);
            if surfaceType=="s"  % treat sphere case separatedly
                binary=createNodalSurface(surfaceType,[3*N,3*N,3*N],PixelSize,1,[1,0,0],[0,1,0],[1.5,1.5,1.5],Threshold);
            else
                binary=createNodalSurface(surfaceType,[3*N,3*N,3*N],PixelSize,1,[1,0,0],[0,1,0],[0,0,0],Threshold);
            end
            MultiplesOfTUC=1;
        elseif j==2
            % create 3x3x3 of D surface in 110, -110 and 001 orientation
            NVoxels=round([sqrt(2)*N,sqrt(2)*N,N]);
            binary=createNodalSurface(surfaceType,3*NVoxels,PixelSize,1,[1,1,0],[-1 1 0],[0,0,0],Threshold);
            MultiplesOfTUC=sqrt(2)*sqrt(2)*1;
        elseif j==3
            % create 3x3x3 of D surface in 111, 1-10 and 11-2 orientation
            NVoxels=round([sqrt(3)*N,sqrt(2)*N,sqrt(6)*N]);
            binary=createNodalSurface(surfaceType,3*NVoxels,PixelSize,1,[1,1,1],[1 -1 0],[0,0,0],Threshold);
            MultiplesOfTUC=sqrt(3)*sqrt(2)*sqrt(6);
        elseif j==4
            % create 1x1x1 of D surface in 100,010,001 orientation but with
            % shifted origin
            NVoxels=round([N,N,N]);
            binary=createNodalSurface(surfaceType,[3*N,3*N,3*N],PixelSize,1,[1,0,0],[0,1,0],[PixelSize/2,1/pi*PixelSize,PixelSize*2/3],Threshold);
            MultiplesOfTUC=1;
        else
            disp("ERROR, should never get to here.");
            pause;
        end

        % adapt known Euler index
        if parameters.EulerMethod=="fitAreaAndMeanCurvatureUseProvidedEuler"
            parameters.ProvidedEulerValue=KnownEulerNumber*MultiplesOfTUC;
        end
        
        % calculate area, average mean curvature and Euler index by Steiner
        % polynomial
        if calculateAreaCurvature == 1
            VoxelsToAnalyse=[NVoxels(1)+1,2*NVoxels(1); NVoxels(2)+1,2*NVoxels(2);NVoxels(3)+1,2*NVoxels(3)];
            fitResults=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,parameters);
            area=fitResults.area;
            meanc=fitResults.meanc;
            Euler=fitResults.Euler;
            SurfaceAreaFromSteinerFit(i,j)=fitResults.area/MultiplesOfTUC;
            AverageMeanCurvatureFromSteinerFit(i,j)=fitResults.meanc;
            AverageEulerIndexFromSteinerFit(i,j)=fitResults.Euler/MultiplesOfTUC;
        end 
        
        % calculate Euclidean distance first phase
        if calculateEDMPhase1==1
            EDM = bwdist(binary); %calculates the euclidean distance for one slice
            EDM=EDM((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
            SolidPhaseMaxEuclidDistance(i,j)=max(max(max(EDM)));
        end
        
        % calculate Euclidean distance second phase
        if calculateEDMPhase2==1
            EDM = bwdist(1-binary); %calculates the euclidean distance for one slice
            EDM=EDM((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
            VoidPhaseMaxEuclidDistance(i,j)=max(max(max(EDM)));
        end
        
        % volume fraction
        nVoxelsWithValue1=length(find(binary==1));
        SolidVolumeFraction(i,j)=nVoxelsWithValue1/(3^3*NVoxels(1)*NVoxels(2)*NVoxels(3));
    
        % Max covering radius transform first phase
        if calculateCRTPhase1==1
            clear paramsCRT;
            paramsCRT.calculateCRTOnlyForRectangularSubset=1;
            paramsCRT.rectangularSubsetForCRTCalculation=[NVoxels'+1,2*NVoxels'];
            [MCRT,CRTCenterPoints,tmpEDM,CRTStatus]=maxCoveringDistanceTransform(binary,paramsCRT);
            clear tmpEDM;
            MCRT=MCRT(...
                paramsCRT.rectangularSubsetForCRTCalculation(1,1):paramsCRT.rectangularSubsetForCRTCalculation(1,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(2,1):paramsCRT.rectangularSubsetForCRTCalculation(2,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(3,1):paramsCRT.rectangularSubsetForCRTCalculation(3,2));
            CRTStatus=CRTStatus(...
                paramsCRT.rectangularSubsetForCRTCalculation(1,1):paramsCRT.rectangularSubsetForCRTCalculation(1,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(2,1):paramsCRT.rectangularSubsetForCRTCalculation(2,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(3,1):paramsCRT.rectangularSubsetForCRTCalculation(3,2));
            clear CRTCenterPoints;
            SolidPhaseMeanValueMaxCoverRadiusTrafo(i,j)=mean(MCRT(find(MCRT > 0)));
            SolidPhaseMeanSquareCoverRadiusTrafo(i,j)=mean((MCRT(find(MCRT > 0))).^2);
            % check that all CRT values correctly computed
            if length(find(CRTStatus(find(MCRT>0))~=1)) > 0
                disp(length(find(CRTStatus(find(MCRT>0))~=1)));
                warn="length(find(CRTStatus(find(MCRT>0))~=1)) > 0   (phase 1)";
                warningsToAlertUserAtEnd=[warningsToAlertUserAtEnd,warn];
                warning("length(find(CRTStatus(find(MCRT>0))~=1)) > 0");
            end
        end
        
        % Max covering radius transform second phase
        if calculateCRTPhase2==1
            clear paramsCRT;
            paramsCRT.calculateCRTOnlyForRectangularSubset=1;
            paramsCRT.rectangularSubsetForCRTCalculation=[NVoxels'+1,2*NVoxels'];
            [MCRT,CRTCenterPoints,tmpEDM,CRTStatus]=maxCoveringDistanceTransform(1-binary,paramsCRT);
            clear tmpEDM;
            MCRT=MCRT(...
                paramsCRT.rectangularSubsetForCRTCalculation(1,1):paramsCRT.rectangularSubsetForCRTCalculation(1,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(2,1):paramsCRT.rectangularSubsetForCRTCalculation(2,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(3,1):paramsCRT.rectangularSubsetForCRTCalculation(3,2));
            CRTStatus=CRTStatus(...
                paramsCRT.rectangularSubsetForCRTCalculation(1,1):paramsCRT.rectangularSubsetForCRTCalculation(1,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(2,1):paramsCRT.rectangularSubsetForCRTCalculation(2,2),...
                paramsCRT.rectangularSubsetForCRTCalculation(3,1):paramsCRT.rectangularSubsetForCRTCalculation(3,2));
            clear CRTCenterPoints;
            %CRTCenterPoints=CRTCenterPoints((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
            VoidPhaseMeanValueMaxCoverRadiusTrafo(i,j)=mean(MCRT(find(MCRT > 0)));
            VoidPhaseMeanSquareCoverRadiusTrafo(i,j)=mean((MCRT(find(MCRT > 0))).^2);
            % check that all CRT values correctly computed
            if length(find(CRTStatus(find(MCRT>0))~=1)) > 0
                warn="length(find(CRTStatus(find(MCRT>0))~=1)) > 0   (phase 0)";
                warningsToAlertUserAtEnd=[warningsToAlertUserAtEnd,warn];
                disp(length(find(CRTStatus(find(MCRT>0))~=1)));
                warning("length(find(CRTStatus(find(MCRT>0))~=1)) > 0");
            end
        end
        
        %F
        [Threshold,j]
    end
end

% just to keep track of how long it takes
endTime = datetime(now,'ConvertFrom','datenum');


% print out surface area values for different orientations
if calculateAreaCurvature == 1
    % surface areas
    disp("-------------------------------------------------------------------");
    disp("SurfaceAreas");
    disp("(Values are in real units for a=1 and correspond to a single TUC)");
    if NOrientations > 1
        maxDiffToCorrectValue=max(abs(SurfaceAreaFromSteinerFit-surfaceArea')')';
    else
        maxDiffToCorrectValue=abs(SurfaceAreaFromSteinerFit-surfaceArea');
    end
    disp("Thresholds / exact value / values for different orienations / max deviation / percentage");
    disp([thresholds' surfaceArea' SurfaceAreaFromSteinerFit maxDiffToCorrectValue round(maxDiffToCorrectValue./surfaceArea'*1000)/10]);
    
    % mean curvatures
    disp("-------------------------------------------------------------------");
    disp("Average Mean Curvature");
    disp("(Values are in real units for a=1)");
    if NOrientations > 1
        maxDiffToCorrectValue=max(abs(AverageMeanCurvatureFromSteinerFit-meanCurvaturesExact')')';
    else
        maxDiffToCorrectValue=abs(AverageMeanCurvatureFromSteinerFit-meanCurvaturesExact');
    end
    disp("Thresholds / exact value / values for different orienations / max deviation / percentage");
    disp([thresholds' meanCurvaturesExact' AverageMeanCurvatureFromSteinerFit maxDiffToCorrectValue round(maxDiffToCorrectValue./meanCurvaturesExact'*1000)/10]);

    if parameters.EulerMethod=="fitAreaMeanCurvatureAndEuler"
        % Euler number curvatures
        disp("-------------------------------------------------------------------");
        disp("Euler number");
        disp("(Values are those for a single translational unit cell.)");
        if NOrientations > 1
            maxDiffToCorrectValue=max(abs(AverageEulerIndexFromSteinerFit-EulerNumberExact')')';
        else
            maxDiffToCorrectValue=abs(AverageEulerIndexFromSteinerFit-EulerNumberExact');
        end
        disp("Thresholds / exact value / values for different orienations / max deviation / percentage");
        disp([thresholds' EulerNumberExact' AverageEulerIndexFromSteinerFit maxDiffToCorrectValue round(maxDiffToCorrectValue./EulerNumberExact'*1000)/10]);
     end
end



% calculate standard deviations and scale
VoidPhaseStandardDevMaxCoverRadiusTrafo=sqrt(max(0,VoidPhaseMeanSquareCoverRadiusTrafo-VoidPhaseMeanValueMaxCoverRadiusTrafo.^2));
SolidPhaseStandardDevMaxCoverRadiusTrafo=sqrt(max(0,SolidPhaseMeanSquareCoverRadiusTrafo-SolidPhaseMeanValueMaxCoverRadiusTrafo.^2));

% Create graph with data for each orientation
close all;
hold off;
figure(1)
hold on;
for j=1:NOrientations
    plot(thresholds,SolidVolumeFraction(:,j),"-+",thresholds,VoidPhaseMaxEuclidDistance(:,j),"--", thresholds,SolidPhaseMaxEuclidDistance(:,j),"-",thresholds,SurfaceAreaFromSteinerFit(:,j),"*");
    errorbar(thresholds,SolidPhaseMeanValueMaxCoverRadiusTrafo(:,j)*PixelSize,SolidPhaseStandardDevMaxCoverRadiusTrafo(:,j)*PixelSize,"*");
    errorbar(thresholds,VoidPhaseMeanValueMaxCoverRadiusTrafo(:,j)*PixelSize,VoidPhaseStandardDevMaxCoverRadiusTrafo(:,j)*PixelSize,"*");
end
if ExportGraphics == 1
    graphics = gca;
    filename=join([fileBase,"_AllOrientationsSeparate.pdf"],"");
    exportgraphics(graphics,filename,'ContentType','vector')
end
hold off;


if calculateEDMPhase1 == 1
    % mean
    disp("-------------------------------------------------------------------");
    disp("EDM max values for hollow phase");
    disp("(Values are in real units)");
    
    VoidPhaseMaxEuclidDistance=VoidPhaseMaxEuclidDistance*PixelSize;
    means=mean(VoidPhaseMaxEuclidDistance,2);
    if NOrientations > 1
        maxDiffToMean=max(abs(VoidPhaseMaxEuclidDistance-means)')';
    else
        maxDiffToMean=abs(VoidPhaseMaxEuclidDistance-means);
    end
    disp("Thresholds / mean / max EDM values for different orienations / max deviation / percentage");
    disp([thresholds' means VoidPhaseMaxEuclidDistance maxDiffToMean round(maxDiffToMean./means*1000)/10]);
end

if calculateEDMPhase2 == 1
    % mean
    disp("-------------------------------------------------------------------");
    disp("EDM max values for solid phase");
    disp("(Values are in real units)");
    
    SolidPhaseMaxEuclidDistance=SolidPhaseMaxEuclidDistance*PixelSize;
    means=mean(SolidPhaseMaxEuclidDistance,2);
    if NOrientations > 1
        maxDiffToMean=max(abs(SolidPhaseMaxEuclidDistance-means)')';
    else
        maxDiffToMean=abs(SolidPhaseMaxEuclidDistance-means);
    end
    disp("Thresholds / mean / max EDM values for different orienations / max deviation / percentage");
    disp([thresholds' means SolidPhaseMaxEuclidDistance maxDiffToMean round(maxDiffToMean./means*1000)/10]);
end



if calculateCRTPhase1 == 1
    % mean
    disp("-------------------------------------------------------------------");
    disp("CRT *mean* values for hollow phase");
    disp("(Values are in real units)");
    VoidPhaseMeanValueMaxCoverRadiusTrafo=VoidPhaseMeanValueMaxCoverRadiusTrafo*PixelSize;
    means=mean(VoidPhaseMeanValueMaxCoverRadiusTrafo,2);
    if NOrientations > 1
        maxDiffToMean=max(abs(VoidPhaseMeanValueMaxCoverRadiusTrafo-means)')';
    else
        maxDiffToMean=abs(VoidPhaseMeanValueMaxCoverRadiusTrafo-means);
    end
    disp("Thresholds / mean / CRT values for different orienations / max deviation / percentage");
    disp([thresholds' means VoidPhaseMeanValueMaxCoverRadiusTrafo maxDiffToMean round(maxDiffToMean./means*1000)/10]);
    % std
    disp("--------");
    disp("CRT *standard deviation* values for hollow phase");
    disp("(Values are in VOXELS, not in real units)");
    VoidPhaseStandardDevMaxCoverRadiusTrafo=VoidPhaseStandardDevMaxCoverRadiusTrafo*PixelSize;
    means=mean(VoidPhaseStandardDevMaxCoverRadiusTrafo,2);
    if NOrientations > 1
        maxDiffToMean=max(abs(VoidPhaseStandardDevMaxCoverRadiusTrafo-means)')';
    else
        maxDiffToMean=abs(VoidPhaseStandardDevMaxCoverRadiusTrafo-means);
    end
    disp("Thresholds / mean / CRT values for different orienations / max deviation / percentage");
    disp([thresholds' means VoidPhaseStandardDevMaxCoverRadiusTrafo maxDiffToMean round(maxDiffToMean./means*1000)/10]);
end

if calculateCRTPhase2 == 1
    disp("-------------------------------------------------------------------");
    % mean
    disp("CRT *mean* values for solid phase (phase 2)");
    disp("(Values are in real units. See pixel size above)");
    SolidPhaseMeanValueMaxCoverRadiusTrafo=SolidPhaseMeanValueMaxCoverRadiusTrafo*PixelSize;
    means=mean(SolidPhaseMeanValueMaxCoverRadiusTrafo,2);
    if NOrientations > 1
        maxDiffToMean=max(abs(SolidPhaseMeanValueMaxCoverRadiusTrafo-means)')';
    else
        maxDiffToMean=abs(SolidPhaseMeanValueMaxCoverRadiusTrafo-means);
    end
    disp("Thresholds / mean / CRT values for different orienations / max deviation / percentage");
    disp([thresholds' means SolidPhaseMeanValueMaxCoverRadiusTrafo maxDiffToMean round(maxDiffToMean./means*1000)/10]);
    % std
    disp("--------");
    disp("CRT *standard deviation* values for hollow phase");
    disp("(Values are in real units; see pixelSize)");
    SolidPhaseStandardDevMaxCoverRadiusTrafo=SolidPhaseStandardDevMaxCoverRadiusTrafo*PixelSize;
    means=mean(SolidPhaseStandardDevMaxCoverRadiusTrafo,2);
    if NOrientations > 1
        maxDiffToMean=max(abs(SolidPhaseStandardDevMaxCoverRadiusTrafo-means)')';
    else
        maxDiffToMean=abs(SolidPhaseStandardDevMaxCoverRadiusTrafo-means);
    end
    disp("Thresholds / mean / CRT values for different orienations / max deviation / percentage");
    disp([thresholds' means SolidPhaseStandardDevMaxCoverRadiusTrafo maxDiffToMean round(maxDiffToMean./means*1000)/10]);
end


% average over different orientations
VoidPhaseMeanValueMaxCoverRadiusTrafo=mean(VoidPhaseMeanValueMaxCoverRadiusTrafo,2);
SolidPhaseMeanValueMaxCoverRadiusTrafo=mean(SolidPhaseMeanValueMaxCoverRadiusTrafo,2);
VoidPhaseMeanSquareCoverRadiusTrafo=mean(VoidPhaseMeanSquareCoverRadiusTrafo,2);
SolidPhaseMeanSquareCoverRadiusTrafo=mean(SolidPhaseMeanSquareCoverRadiusTrafo,2);
SolidVolumeFraction=mean(SolidVolumeFraction,2);
SolidPhaseMaxEuclidDistance=mean(SolidPhaseMaxEuclidDistance,2);
VoidPhaseMaxEuclidDistance=mean(VoidPhaseMaxEuclidDistance,2);
SurfaceAreaFromSteinerFit=mean(SurfaceAreaFromSteinerFit,2);
AverageMeanCurvatureFromSteinerFit=mean(AverageMeanCurvatureFromSteinerFit,2)
AverageEulerIndexFromSteinerFit=mean(AverageEulerIndexFromSteinerFit,2);


% convert square means to std
VoidPhaseStandardDevMaxCoverRadiusTrafo=sqrt(max(0,VoidPhaseMeanSquareCoverRadiusTrafo-VoidPhaseMeanValueMaxCoverRadiusTrafo.^2));
SolidPhaseStandardDevMaxCoverRadiusTrafo=sqrt(max(0,SolidPhaseMeanSquareCoverRadiusTrafo-SolidPhaseMeanValueMaxCoverRadiusTrafo.^2));

% create multi-panel figure with results
figure(2)
hold off
subplot(2,3,1)
plot(thresholds,SolidVolumeFraction,"-+");
xlabel('threshold');
ylabel('solid volume fraction');
subplot(2,3,2)
plot(thresholds,surfaceArea,"-",thresholds,SurfaceAreaFromSteinerFit,"*");
xlabel('threshold');
ylabel('surface area (TUC, a=1)');
subplot(2,3,3)
hold on
plot(thresholds,VoidPhaseMaxEuclidDistance,"b*", thresholds,SolidPhaseMaxEuclidDistance,"bo");
errorbar(thresholds,SolidPhaseMeanValueMaxCoverRadiusTrafo*PixelSize,SolidPhaseStandardDevMaxCoverRadiusTrafo*PixelSize,'vertical','g-s');
errorbar(thresholds,VoidPhaseMeanValueMaxCoverRadiusTrafo*PixelSize,VoidPhaseStandardDevMaxCoverRadiusTrafo*PixelSize,'vertical','k-s');
xlabel('threshold');
ylabel('max of EDM, mean of CRT');
hold off
subplot(2,3,4)
plot(thresholds,AverageMeanCurvatureFromSteinerFit,"*-");
xlabel('threshold');
ylabel('average mean curvature');
subplot(2,3,5)
plot(thresholds,AverageEulerIndexFromSteinerFit,"*-");
xlabel('threshold');
ylabel('Euler index');
ylim([min(-20,1.1*min(AverageEulerIndexFromSteinerFit)) max(4,1.1*max(AverageEulerIndexFromSteinerFit))]);

if ExportGraphics == 1
    graphics = gca;
    filename=join([fileBase,".pdf"],"");
    exportgraphics(graphics,filename,'ContentType','vector')
end

% store results in a file
results.surfaceType=surfaceType;
results.NOrientations=NOrientations;
results.PixelSize=PixelSize;
results.N=N;
results.startTime=startTime;
results.endtime=endTime;
results.surfaceAreaDiscretisation=surfaceAreaDiscretisation;
results.thresholds=thresholds;
results.SolidPhaseStandardDevMaxCoverRadiusTrafo=SolidPhaseStandardDevMaxCoverRadiusTrafo;
results.VoidPhaseStandardDevMaxCoverRadiusTrafo=VoidPhaseStandardDevMaxCoverRadiusTrafo;
results.VoidPhaseMeanValueMaxCoverRadiusTrafo=VoidPhaseMeanValueMaxCoverRadiusTrafo;
results.SolidPhaseMeanValueMaxCoverRadiusTrafo=SolidPhaseMeanValueMaxCoverRadiusTrafo;
results.VoidPhaseMeanSquareCoverRadiusTrafo=VoidPhaseMeanSquareCoverRadiusTrafo;
results.SolidPhaseMeanSquareCoverRadiusTrafo=SolidPhaseMeanSquareCoverRadiusTrafo;
results.SolidVolumeFraction=SolidVolumeFraction;
results.SolidPhaseMaxEuclidDistance=SolidPhaseMaxEuclidDistance;
results.VoidPhaseMaxEuclidDistance=VoidPhaseMaxEuclidDistance;

% display warnings
disp("Please note these warnings above:");
for i=1:length(warningsToAlertUserAtEnd)
    disp(warningsToAlertUserAtEnd(i));
end           



