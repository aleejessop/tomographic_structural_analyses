

function binary = createNodalSurface(surfaceName, ImageSize, PixelSize, CrystA, LatticeVectorA, LatticeVectorB, Origin, Threshold)

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
if abs(A*B') > 1e-8
    error("A and B not orthogonal");
end

% make sure A, B, C are colum vectors
A=A';
B=B';
C=C';
if isequal(size(A),[1 3]) & isequal(size(B),[1 3]) & isequal(size(C),[1 3])  
    disp(A);
    disp(B);
    disp(C);
    error("A, B or C not column vector");
end

% determine rotation matrix
%
% angle around z-axis to bring A into xz plan
phi=-atan2(A(2),A(1));
Rz = [cos(phi) -sin(phi) 0 ; sin(phi) cos(phi) 0 ; 0 0 1] ;
AA=Rz*A;
% determine angle for rotation around y-axis, so that A falls onto 
% x-axis when Ry*(Rz*A) is applied
alpha=atan2(A(3),sqrt((A(1))^2+(A(2))^2));
Ry = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
AAA=Ry*(Rz*A);
% angle around x-axis to bring vector b onto the y-axis
BB=Ry*(Rz*B);
if abs(BB(1))/(B'*B) > 10e-8
    error("Vector B should be in yz plane after rotation");
end
gamma=-atan2(BB(3),BB(2));
Rx = [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma)];
BBB=Rx*BB;
[AAA,BBB]
devA=sqrt((AAA(1)-1)^2+(AAA(2)-0)^2+(AAA(3)-0)^2)/sqrt(A'*A);
devB=sqrt((BBB(1))^2+(BBB(2)-1)^2+(BBB(3)-0)^2)/sqrt(B'*B);
if devA > 1e-8 || devB > 1e-8
    error("Coordinate rotation wrong");
end
RotMatrix=Rx*Ry*Rz;
InverseRotMatrix=inv(RotMatrix);

for i=1:nx
    for j=1:ny
        for k=1:nz
            % pixel in image coordinates in microns
            xyz=[(i-0.5)*PixelSize,(j-0.5)*PixelSize,(k-0.5)*PixelSize]; % pixel coordinate in microns
            T=xyz-Origin;
             % position of xyz in rotated reference frame in microns
            XYZ=InverseRotMatrix*xyz';
            X=XYZ(1);
            Y=XYZ(2);
            Z=XYZ(3);
%             X=T*A';  % position of xyz in rotated reference frame in microns
%             Y=T*B';
%             Z=T*C';
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
            elseif lower(surfaceName)=="scspheres" || lower(surfaceName)=="scs"
                spherecenters=[[0,0,0];[1,0,0];[0,1,0];[1,1,0];[0,0,1];[1,0,1];[0,1,1];[1,1,1]];
                r=[mod(X,1),mod(Y,1),mod(Z,1)];
                nodal(i,j,k)=sqrt(min(sum((spherecenters-r).^2,2)));   
            elseif lower(surfaceName)=="bccspheres" || lower(surfaceName)=="bccs"
                spherecenters=[[0.5,0.5,0.5];[0,0,0];[1,0,0];[0,1,0];[1,1,0];[0,0,1];[1,0,1];[0,1,1];[1,1,1]];
                r=[mod(X,1),mod(Y,1),mod(Z,1)];
                nodal(i,j,k)=sqrt(min(sum((spherecenters-r).^2,2)));   
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