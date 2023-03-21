
%calculate angle between vectors

Znorm = [0 0 1];

for i = 1:size(Orientationanalysis,1)

    Ztest = Orientationanalysis(i,:);
    dd = dot(Znorm,Ztest);
    Znorm_mag = norm(Znorm);
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

