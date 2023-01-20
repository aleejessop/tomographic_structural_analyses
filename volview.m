
%plot isosurface to look like volume

function volview(binary,s)
binary = double(binary);
B = ones(s+2);
B(2:end-1,2:end-1,2:end-1) = binary;
isosurface(B); axis equal
xlabel('x')
ylabel('y')
zlabel('z')
end