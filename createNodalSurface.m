

function binary = createNodalSurface(surfaceName,ImageSize, PixelSize, CrystA, LatticeVectorA, LatticeVectorB, Origin, Threshold)

ImageSize;
nx=ImageSize(1);
ny=ImageSize(2);
nz=ImageSize(3);

nodal=ones(nx,ny,nz);
binary=ones(nx,ny,nz);

% a=11.5;   % P surface repeat length (crystallographic a) in microns
% threshold=-0.13; % threshold for binarisation also sets the solid volume fraction

% define crystallographic vectors in terms of the 
% image coordinates
% origin of crystallographic reference frame in image coordinates
% crystallographic A direction (arbitrary length) in image coordinates
% crystallographic B direction (arbitrary length) in image coordinates
            % beware: need to manually make sure they are orthogonal, can
            % use the dot product to check.
if length(LatticeVectorA) ~= 3 || length(LatticeVectorB) ~=3 || length(Origin) ~= 3
    disp("Error, lattice vectors need to have three elements");
    pause;
end
if abs(LatticeVectorA*LatticeVectorB')>10^(-9)
    disp("Error, lattice vectors not perpendicular");
    pause;
end
A=LatticeVectorA;
A=A/sqrt(A*A');  % making them unit length
B=LatticeVectorB;
B=B/sqrt(B*B');
C=cross(A,B); % calculate C vector as perpendicular to A and B


for i=1:nx
    for j=1:ny
        for k=1:nz
            xyz=[(i-0.5)*PixelSize,(j-0.5)*PixelSize,(k-0.5)*PixelSize]; % pixel coordinate in microns
            T=xyz-Origin;
            X=T*A';  % position of xyz in rotated reference frame in microns
            Y=T*B';
            Z=T*C';
            if lower(surfaceName)=="primitive" || lower(surfaceName)=="p"
                nodal(i,j,k) = cos(X*2*pi/CrystA)+cos(Y*2*pi/CrystA)+cos(Z*2*pi/CrystA);
            elseif lower(surfaceName)=="gyroid" || lower(surfaceName)=="g"
                nodal(i,j,k) = cos(X*2*pi/CrystA).*sin(Y*2*pi/CrystA)+cos(Y*2*pi/CrystA).*sin(Z*2*pi/CrystA)+cos(Z*2*pi/CrystA).*sin(X*2*pi/CrystA);
            elseif lower(surfaceName)=="diamond" || lower(surfaceName)=="d"
                nodal(i,j,k)=sin(X*2*pi/CrystA)*sin(Y*2*pi/CrystA)*sin(Z*2*pi/CrystA)...
                + sin(X*2*pi/CrystA)*cos(Y*2*pi/CrystA)*cos(Z*2*pi/CrystA)...
                + cos(X*2*pi/CrystA)*sin(Y*2*pi/CrystA)*cos(Z*2*pi/CrystA)...
                + cos(X*2*pi/CrystA)*cos(Y*2*pi/CrystA)*sin(Z*2*pi/CrystA); % D surface
            elseif lower(surfaceName)=="sphere" || lower(surfaceName)=="s"
                nodal(i,j,k)=X^2+Y^2+Z^2;
            else
                disp("Should never get to here. Wrong surface type");
                pause;
            end
            % threshold
            if nodal(i,j,k)>Threshold
                binary(i,j,k)=1;
            else
                binary(i,j,k)=0;
            end
        end
    end
end
end