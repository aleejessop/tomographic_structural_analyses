clear all;
close all;
N=25;
PixelSize=1/N;
surfaceType="primitive";

thresholds=0;

% number of different orientations to use to evaluate nodal surface
% needs to match if statements in for j=1:NOrientations loop
NOrientations=1; 

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


% surface area of unit cell (not calculated from voxelised data; no 
% need to average over orientations)
surfaceArea=zeros(length(thresholds),1);

% just to keep track of how long it takes
startTime = datetime(now,'ConvertFrom','datenum');

% filebase to store results
fileBase=join(["Properties_",surfaceType,"_NOri",string(NOrientations),"_N",string(N)],"");
disp("Will print results to the following file:");
disp(fileBase);

for i=1:length(thresholds)
    Threshold=thresholds(i);
    
    %  calculate surface area of TUC (a=1)
    surfaceAreaDiscretisation=0.05;
    surfaceArea(i)=SurfaceAreaTranslationalUnitNodalSurface(surfaceType,Threshold,surfaceAreaDiscretisation);
    
    for j=1:NOrientations

        % create 3x3x3 of D surface in 100 orientation
        if j==1
            NVoxels=round([N,N,N]);
            binary=createNodalSurface(surfaceType,[3*N,3*N,3*N],PixelSize,1,[1,0,0],[0,1,0],[0,0,0],Threshold);
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

        % calculate area, average mean curvature and Euler index by Steiner
        % polynomial
        VoxelsToAnalyse=[NVoxels(1)+1,2*NVoxels(1); NVoxels(2)+1,2*NVoxels(2);NVoxels(3)+1,2*NVoxels(3)];
        fractionEuclidDistanceToFit=0.75;
        fitResults=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,fractionEuclidDistanceToFit);
        area=fitResults.area;
        meanc=fitResults.meanc;
        Euler=fitResults.Euler;
        SurfaceAreaFromSteinerFit(i,j)=fitResults.area/MultiplesOfTUC;
        AverageMeanCurvatureFromSteinerFit(i,j)=fitResults.meanc;
        AverageEulerIndexFromSteinerFit(i,j)=fitResults.Euler/MultiplesOfTUC;
        
        % calculate Euclidean distance first phase
        EDM = bwdist(binary); %calculates the euclidean distance for one slice
        EDM=EDM((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
        SolidPhaseMaxEuclidDistance(i,j)=max(max(max(EDM)))*PixelSize;

        % calculate Euclidean distance second phase
        EDM = bwdist(1-binary); %calculates the euclidean distance for one slice
        EDM=EDM((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
        VoidPhaseMaxEuclidDistance(i,j)=max(max(max(EDM)))*PixelSize;

        % volume fraction
        nVoxelsWithValue1=length(find(binary==1));
        SolidVolumeFraction(i,j)=nVoxelsWithValue1/(3^3*NVoxels(1)*NVoxels(2)*NVoxels(3));
    
        % Max covering radius transform first phase
        [MCRT,image1]=maxCoveringDistanceTransform(binary);
        MCRT=MCRT((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
        image1=image1((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
        SolidPhaseMeanValueMaxCoverRadiusTrafo(i,j)=mean(MCRT(find(MCRT ~= 0)));
        SolidPhaseMeanSquareCoverRadiusTrafo(i,j)=mean((MCRT(find(MCRT ~= 0))).^2);

        % Max covering radius transform second phase
        [MCRT,image1]=maxCoveringDistanceTransform(1-binary);
        MCRT=MCRT((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
        image1=image1((NVoxels(1)+1):2*NVoxels(1),(NVoxels(2)+1):2*NVoxels(2),(NVoxels(3)+1):2*NVoxels(3));
        VoidPhaseMeanValueMaxCoverRadiusTrafo(i,j)=mean(MCRT(find(MCRT ~= 0)));
        VoidPhaseMeanSquareCoverRadiusTrafo(i,j)=mean((MCRT(find(MCRT ~= 0))).^2);

        %F
        [Threshold,j]
    end
end

% just to keep track of how long it takes
endTime = datetime(now,'ConvertFrom','datenum');


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
graphics = gca;
filename=join([fileBase,"_AllOrientationsSeparate.pdf"],"");
exportgraphics(graphics,filename,'ContentType','vector')
hold off;

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

graphics = gca;
filename=join([fileBase,".pdf"],"");
exportgraphics(graphics,filename,'ContentType','vector')

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



