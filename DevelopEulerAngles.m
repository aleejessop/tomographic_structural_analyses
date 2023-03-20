A0=[1 0 0]
B0=[0 1 0]
C0=[0 0 1]
for t=[-pi:0.05:-pi/2,-pi/2:0.05:0,0:0.05:pi/2,pi/2:0.05:pi]
    Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)]; 
    for s=[-pi:0.05:-pi/2,-pi/2:0.05:0,0:0.05:pi/2,pi/2:0.05:pi]
        Rx = [1 0 0; 0 cos(s) -sin(s); 0 sin(s) cos(s)];
        for u=[-pi:0.05:-pi/2,-pi/2:0.05:0,0:0.05:pi/2,pi/2:0.05:pi]
            Rz = [cos(u) -sin(u) 0 ; sin(u) cos(u) 0 ; 0 0 1] ;
            A=Rz*Ry*Rx*A0';
            B=Rz*Ry*Rx*B0';
            C=Rz*Ry*Rx*C0';
            
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
        end
    end
end


