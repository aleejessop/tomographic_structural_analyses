
function [subvolume] = subvolume_extraction(binary, size_subvolume)
binary = double(binary);
nx = 1:size(binary,1);
ny = 1:size(binary,2);
nz = 1:size(binary,3);
x = reshape(nx,size_subvolume,(size(binary,1)/size_subvolume));
x = x';
y = reshape(ny,size_subvolume,(size(binary,2)/size_subvolume));
y = y';
z = reshape(nz,size_subvolume,(size(binary,3)/size_subvolume));
z = z';
for i = 1:size(z,1)
    for j = 1:size(y,1)
        for k = 1:size(x,1)
            subvolume(k,j,i,:,:,:) = binary(x(k,:),y(j,:),z(i,:));
        end
    end
end
end

