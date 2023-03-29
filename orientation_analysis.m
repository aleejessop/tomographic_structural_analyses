
%calculate angle between vectors
clear
Znorm = [0 0 1];


for i = 1:size(Orientationanalysis,1)

    Ztest = Orientationanalysis(i,:);
    dd = dot(Znorm,Ztest);
    Znorm_mag = norm(Znorm); %norm calculates magnitude of a vector
    Ztest_mag = norm(Ztest);
    angle(i) = acosd(dd/(Znorm_mag*Ztest_mag))
end


quiver3(0,0,0,0,0,1)
hold on
for ii = 1:size(Orientationanalysis,1)

    quiver3(0,0,0,Orientationanalysis(ii,1),Orientationanalysis(ii,2),Orientationanalysis(ii,3))
    hold on
    axis equal

end



[x_s,y_s,z_s] = sphere(12);
figure(1);
surf(x_s, y_s, z_s,'FaceColor','w')
axis equal
axis off
map = [0.8 0.8 0.8];
colormap(map)
shading flat
hold on
xlabel('x')
ylabel('y')
zlabel('z')
marker_color = [0 0 0; 0.07 0.62 1.00; 0.72 0.27 1.00; 0.00 0.57 0.43; 0.98 0.98 0.51; 0.99 0.71 0.59;...
    0.99 0.74 0.99; 0.71 0.98 0.52; 0.67 0.98 0.98];
for i = 1:9
    scatter3(Orientationanalysis(i,1),Orientationanalysis(i,2),Orientationanalysis(i,3), 20,'MarkerEdgeColor',marker_color(i,:),'MarkerFaceColor',marker_color(i,:),'LineWidth',1)
    hold on
end
    grid off

