function area = SurfaceAreaTranslationalUnitNodalSurface(surfaceName,threshold,discretisation)

if discretisation < 0.0001 || discretisation > 0.1
    disp("Error: discretisation < 10^-6 || discretisation > 0.1");
    pause;
end

% create grid
[x,y,z] = meshgrid([-0.5:discretisation:0.5]); 

if lower(surfaceName)=="primitive" || lower(surfaceName)=="p"
    V = cos(x*2*pi)+cos(y*2*pi)+cos(z*2*pi);
elseif lower(surfaceName)=="gyroid" || lower(surfaceName)=="g"
    V = cos(x*2*pi).*sin(y*2*pi)+cos(y*2*pi).*sin(z*2*pi)+cos(z*2*pi).*sin(x*2*pi);
elseif lower(surfaceName)=="diamond" || lower(surfaceName)=="d"
    V = sin(x*2*pi).*sin(y*2*pi).*sin(z*2*pi)+ sin(x*2*pi).*cos(y*2*pi).*cos(z*2*pi)+ cos(x*2*pi).*sin(y*2*pi).*cos(z*2*pi)+ cos(x*2*pi).*cos(y*2*pi).*sin(z*2*pi);
else
    disp("Should never get to here. Wrong surface type");
    pause;
end


s=isosurface(x,y,z,V,threshold);
p = patch(s);
%isonormals(x,y,z,V,p)
%view(3);
%set(p,'FaceColor',[0.5 1 0.5]);  
%set(p,'EdgeColor','none');
%camlight;
%lighting gouraud;

verts = get(p, 'Vertices');
faces = get(p, 'Faces');
if length(faces)>0
    if length(verts)==0
        disp("ERROR in SurfaceAreaNodalCalc");
        pause;
    end
    a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
    b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
    c = cross(a, b, 2);
    area = 1/2 * sum(sqrt(sum(c.^2, 2)));
else
    area=0;
end
%fprintf('\nThe surface area is %f\n\n', area);

end

