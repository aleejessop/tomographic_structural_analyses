maxEDMInPixels=[];
shiftInPixels=[];

%clear all; 
%clear;
% close all;

addpath('../');



for N=[9 13 27 53 79 111 177 12 24 48 96 192]

    disp(N);
    
    CrystA = 1;
    PixelSize=CrystA/N;
    NVoxels=round([N,N,N]);
    VoxelsToAnalyse=[1,N;1,N;1,N];

    % parameters for test case "s" (single sphere) or "c" (single cylinder)
    surfaceName="s";  % "s"
    radius=0.25;
    Threshold=(radius)^2;   % threshold is radius^2
    origin=[0.5,0.5,0.5];
    fractionEuclidDistanceToFit=0.4;
    KnownEulerNumber=1;
    parameters.ProvidedAreaValue=4*pi*radius^2;
    parameters.ProvidedMeanCurvatureValue=1/radius;
    parameters.ProvidedEulerValue=4*pi*radius^2*1/(radius)^2/2/pi;

    binary=createNodalSurface(surfaceName,[N,N,N],PixelSize,CrystA,[1,0,0],[0,1,0],origin,Threshold);

    % PixelSize = 0.7;
    % VoxelsToAnalyse = [42, 125; 50, 106; 35, 122];
    % incrementSteps=0.5;  % this is in pixel units %% parameter not needed.
    % Should always be 1

    parameters.fractionEuclidDistanceToFit=fractionEuclidDistanceToFit;
    %parameters.EulerMethod="fitAreaAndMeanCurvatureUseProvidedEuler"; 
    parameters.EulerMethod="fitHorizontalShiftToKnownExactDataForCodeDevelopment"; 
    %parameters.removeFirstLayerFromFit=0;
    fitResults=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,parameters);
    %fitResults=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,fractionEuclidDistanceToFit);
    area=fitResults.area;
    meanc=fitResults.meanc;
    Euler=fitResults.Euler;
    layerD=fitResults.layerD;
    layerV=fitResults.layerV;

    disp("Values of area, (average) mean curvature, Euler index:");
    disp([area,meanc,Euler]);
    finerD=1.05*min(layerD):0.01:1.05*max(layerD);
    exactCurve=parameters.ProvidedAreaValue*(finerD+parameters.ProvidedMeanCurvatureValue*finerD.^2)+1/3*parameters.ProvidedEulerValue*2*pi*finerD.^3;
    figure;
    plot(layerD,layerV,"+",finerD,exactCurve,"-",fitResults.shiftedDistances, fitResults.volumesUpToDistance,"o")
    % figure;
    % volshow(binary(42:125,50:106,35:122))

    disp('Exact values')
    [parameters.ProvidedAreaValue,parameters.ProvidedMeanCurvatureValue,parameters.ProvidedEulerValue]
    disp('Shift')
    [fitResults.shiftAInPixels]
    disp('max edm values');
    [fitResults.maxEDMPhase0,fitResults.maxEDMPhase1]
    
    maxEDMInPixels=[maxEDMInPixels,[fitResults.maxEDMPhase0,fitResults.maxEDMPhase1]]
    shiftInPixels=[shiftInPixels,fitResults.shiftAInPixels]
   

end

