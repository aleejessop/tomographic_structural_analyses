clear all; 
clear;
% close all;
N=50;
CrystA = 1;
PixelSize=CrystA/N;
NVoxels=round([N,N,N]);
VoxelsToAnalyse=[N+1, 3*N-N; N+1, 3*N-N; N+1, 3*N-N];

% parameters for bicontinuous surfaces
Threshold=0.; %0.43 for a volume fraction of 0.32
surfaceName="p";
origin=[0,0,0];
fractionEuclidDistanceToFit=0.5;

% parameters for test case "s" (single sphere) or "c" (single cylinder)
%surfaceName="c";  % "s"
%Threshold=(0.25)^2;   % threshold is radius^2
%origin=[1.5,1.5,1.5];
%fractionEuclidDistanceToFit=0.4;

binary=createNodalSurface(surfaceName,[3*N,3*N,3*N],PixelSize,CrystA,[1,0,0],[0,1,0],origin,Threshold);

% PixelSize = 0.7;
% VoxelsToAnalyse = [42, 125; 50, 106; 35, 122];
incrementSteps=0.5;  % this is in pixel units

fitResults=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,fractionEuclidDistanceToFit);
area=fitResults.area;
meanc=fitResults.meanc;
Euler=fitResults.Euler;
layerD=fitResults.layerD;
layerV=fitResults.layerV;

disp("Values of area, (average) mean curvature, Euler index:");
disp([area,meanc,Euler]);
finerD=1.05*min(layerD):0.01:1.05*max(layerD);
fittedCurve=area*(finerD+meanc*finerD.^2)+1/3*Euler*2*pi*finerD.^3;
figure;
plot(layerD,layerV,"+",finerD,fittedCurve,"-")
% figure;
% volshow(binary(42:125,50:106,35:122))



