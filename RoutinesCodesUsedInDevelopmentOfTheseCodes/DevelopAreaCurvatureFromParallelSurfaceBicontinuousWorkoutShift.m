%clear all; 
%clear;
% close all;

addpath('../');

maxEDMInPixels=[];
shiftInPixels=[];

%for N=[7 11 24 48 96 15 27 67 113 145 192]
for N=[7 11 24 48 96 15 ]
    
    disp("Working on ");
    disp(N);
    
    CrystA = 1;
    PixelSize=CrystA/N;
    NVoxels=round([N,N,N]);
    VoxelsToAnalyse=[N+1, 3*N-N; N+1, 3*N-N; N+1, 3*N-N];

    % parameters for bicontinuous surfaces
    Threshold=0.0; %0.43 for a volume fraction of 0.32
    surfaceName="g";
    origin=[0,0,0];
    fractionEuclidDistanceToFit=0.8;
    KnownEulerNumber=-8;

   

    % parameters for bicontinuous surfaces
    Threshold=0.; %0.43 for a volume fraction of 0.32
    surfaceName="g";
    origin=[0,0,0];
    fractionEuclidDistanceToFit=0.8;
    parameters.ProvidedAreaValue=3.0917;
    parameters.ProvidedMeanCurvatureValue=0;
    parameters.ProvidedEulerValue=-8;
    
     % parameters for bicontinuous surfaces
    Threshold=0.; %0.43 for a volume fraction of 0.32
    surfaceName="p";
    origin=[0,0,0];
    fractionEuclidDistanceToFit=0.5;
    parameters.ProvidedAreaValue=2.3526;
    parameters.ProvidedMeanCurvatureValue=0;
    parameters.ProvidedEulerValue=-4;
    
     % parameters for bicontinuous surfaces
    Threshold=0.; %0.43 for a volume fraction of 0.32
    surfaceName="d";
    origin=[0,0,0];
    fractionEuclidDistanceToFit=0.6;
    parameters.ProvidedAreaValue=3.84;
    parameters.ProvidedMeanCurvatureValue=0;
    parameters.ProvidedEulerValue=-16;
   
     % parameters for bicontinuous surfaces
    Threshold=0.665; %0.43 for a volume fraction of 0.32
    surfaceName="p";
    origin=[0,0,0];
    fractionEuclidDistanceToFit=0.5;
    parameters.ProvidedAreaValue=2.17;
    parameters.ProvidedMeanCurvatureValue=-1.13;
    parameters.ProvidedEulerValue=-4;
   
    binary=createNodalSurface(surfaceName,[3*N,3*N,3*N],PixelSize,CrystA,[1,0,0],[0,1,0],origin,Threshold);

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



