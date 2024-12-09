clear
filename = 'C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\Sea Urchin Diamond_111.tiff';
greyscale = tiffreadVolume(filename);
binary = imbinarize(greyscale,0.51);
PixelSize = 0.7;
VoxelsToAnalyse = [48,142;48,142;48,142];
parameters.fractionEuclidDistanceToFit=0.50;
parameters.EulerMethod="fitAreaMeanCurvatureAndEuler";


fitResults=CalculateAreaMeancurvEulerBySteinerOfVoxelisedParSurf(binary,PixelSize,VoxelsToAnalyse,parameters);
area=fitResults.area;
area_per_volume = (area*10^-9)/0.004;
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