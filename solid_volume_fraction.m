

function [volume_fraction] = solid_volume_fraction(subvolume)
for i = 1:size(subvolume,3)
    for j = 1:size(subvolume,2)
        for k = 1:size(subvolume,1)
            volume_fraction(k,j,i) = nanmean(subvolume(k,j,i,:,:,:),'all');
        end
    end
end
end

