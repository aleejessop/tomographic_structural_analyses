
%calculate angle between vectors
clear
Zref = [1 0 0];

for i = 1:size(Orientationanalysis,1)

    Ztest = Orientationanalysis(i,:);
    dd = dot(Zref,Ztest);
    Znorm_mag = norm(Zref); %norm calculates magnitude of a vector
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

%% plot 111 direction relative to reference onto sphere

[x_s,y_s,z_s] = sphere(48);
R = 0.99;

figure(1);
surf(x_s*R, y_s*R, z_s*R,'FaceColor','w','EdgeColor','none')
axis equal
axis off
map = [0.8 0.8 0.8];
colormap(map)
hold on

%plot longtudinal lines

Ai = [0:15:180];

for i = 1:length(Ai)
    E=-pi:pi/180:pi;
    A=(Ai(i)*pi/180)*ones(size(E));
    rE=ones(size(E));
    [XE,YE,ZE]=sph2cart(A,E,rE);
    plot3(XE,YE,ZE,'Color',[0.8 0.8 0.8])
    hold on
end

%plot latitude lines

Ei = [-180:15:180];

for i = 1:length(Ei)
    A=-pi:pi/180:pi;
    E=(Ei(i)*pi/180)*ones(size(A));
    rE=ones(size(E));
    [XE,YE,ZE]=sph2cart(A,E,rE);
    plot3(XE,YE,ZE,'Color',[0.9 0.9 0.9])
    hold on
end


%plot Z reference point

scatter3(0,0,1,20,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)

%plot ZX plane

% scatter3(1,0,0, 20,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1)
R1 = 1.01;
E=-pi:pi/180:pi;
A=(0*pi/180)*ones(size(E));
rE=ones(size(E));
[XE,YE,ZE]=sph2cart(A,E,rE);
plot3(XE*R1,YE*R1,ZE*R1,'Color',[0.5 0.5 0.5])

%plot ZY plane

% scatter3(0,1,0, 20,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1)
E=-pi:pi/180:pi;
A=(90*pi/180)*ones(size(E));
rE=ones(size(E));
[XE,YE,ZE]=sph2cart(A,E,rE);
plot3(XE*R1,YE*R1,ZE*R1,'Color',[0.5 0.5 0.5])

%plot orientation points

scatter3(Orientationanalysis(:,1),Orientationanalysis(:,2),Orientationanalysis(:,3), 3,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1)




