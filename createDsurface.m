

function binary = createDsurface(repeat_length, threshold, Origin,A,B)
d=1;  % voxel size in microns
n=165*2;  % microscopy data set size (in pixel)
L=n*d;  % resulting sample size in microns
nodal=ones(n,n,n);
binary=ones(n,n,n);

% define crystallographic vectors in terms of the 
% image coordinates
% O=[19 3.505 12];    % origin of crystallographic reference frame in image coordinates in microns
% A=[23.495,17,-7.505];  % crystallographic A direction (arbitrary length) in image coordinates
% B=[7.495,30.505,-7];  % crystallographic B direction (arbitrary length) in image coordinates
            % beware: need to manually make sure they are orthogonal, can
            % use the dot product to check.
A=A/sqrt(A*A');  % making them unit length
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


% a=11.5;   % P surface repeat length (crystallographic a) in microns
% threshold=-0.13; % threshold for binarisation also sets the solid volume fraction




for i=1:n
    for j=1:n
        for k=1:n
            % pixel in image coordinates in microns
            xyz=[i*d,j*d,k*d]; 
            T=xyz-Origin;
            % position of xyz in rotated reference frame in microns
            XYZ=InverseRotMatrix*xyz';
            X=XYZ(1);
            Y=XYZ(2);
            Z=XYZ(3);
            %X=T*A';  
            %Y=T*B';
            %Z=T*C';
            nodal(i,j,k)=sin(X*2*pi/repeat_length)*sin(Y*2*pi/repeat_length)*sin(Z*2*pi/repeat_length)...
                + sin(X*2*pi/repeat_length)*cos(Y*2*pi/repeat_length)*cos(Z*2*pi/repeat_length)...
                + cos(X*2*pi/repeat_length)*sin(Y*2*pi/repeat_length)*cos(Z*2*pi/repeat_length)...
                + cos(X*2*pi/repeat_length)*cos(Y*2*pi/repeat_length)*sin(Z*2*pi/repeat_length); % D surface
            if nodal(i,j,k)>threshold
                binary(i,j,k)=1;
            else
                binary(i,j,k)=0;
            end
        end
    end
end
end