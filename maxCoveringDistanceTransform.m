
function [image0,image1,D] = maxCoveringDistanceTransform(binary)
    D = bwdist(binary);

    DV = [];
    X = [];
    Y = [];
    Z = [];
    k = 1;
    for x = 1:size(D,1)
        for y = 1:size(D,2)
            for z = 1:size(D,3)
                % check consistency
                if binary(x,y,z)==0 && D(x,y,z)<1
                    disp("ERROR in maxCoveringDistanceTransform");
                    disp("binary(x,y,z)==0 && D(x,y,z)<1");
                    disp(binary(x,y,z));
                    disp(D(x,y,z));
                    disp([x,y,z]);
                    pause;
                end
                if binary(x,y,z)==1 && D(x,y,z) ~=0
                    disp("ERROR in maxCoveringDistanceTransform");
                    disp("binary(x,y,z)==1 && D(x,y,z) ~=0");
                    pause;
                end
                % include pore space voxels in list
                if binary(x,y,z)==0
                    DV(k) = D(x,y,z);
                    X(k) = x;
                    Y(k) = y;
                    Z(k) = z;
                    k = k+1;
                end
            end
        end
    end

    [DVS, index] = sort(DV,'descend');
    Xsort = X(index);
    Ysort = Y(index);
    Zsort = Z(index);

    image0 = zeros(size(D,1),size(D,2),size(D,3));
    image1 = zeros(size(D,1),size(D,2),size(D,3));

    for k = 1:length(DVS)

        kD = DVS(k);
        kX = Xsort(k);
        kY = Ysort(k);
        kZ = Zsort(k);
% 
%         for x = 1:size(D,1)
%             for y = 1:size(D,1)
%                 for z = 1:size(D,1)
%                     if image0(x,y,z) == 0 && sqrt((x-kX)^2 + (y-kY)^2 + (z-kZ)^2)<= kD
%                         image0(x,y,z) = kD;
%                         image1(kX,kY,kZ) = 1;
%                     end
%                 end
%             end
%         end
    

        if kX-kD<2
            x = 1:kX+(round(kD)+1);
        elseif kX+kD>size(D,1)-1
            x = kX-(round(kD)+1):size(D,1);
        else
            x = kX-(round(kD)+1):kX+(round(kD)+1);
        end

        if kY-kD<2
            y = 1:kY+(round(kD)+1);
        elseif kY+kD>size(D,2)-1
            y = kY-(round(kD)+1):size(D,2);
        else
            y = kY-(round(kD)+1):kY+(round(kD)+1);
        end

        if kZ-kD<2
            z = 1:kZ+(round(kD)+1);
        elseif kZ+kD>size(D,3)-1
            z = kZ-(round(kD)+1):size(D,3);
        else
            z = kZ-(round(kD)+1):kZ+(round(kD)+1);
        end

    
        for xx = x(1):x(end)
            for yy = y(1):y(end)
                for zz = z(1):z(end)
                    if binary(xx,yy,zz) == 0
                        if image0(xx,yy,zz) == 0 && sqrt((xx-kX)^2 + (yy-kY)^2 + (zz-kZ)^2)<= kD
                            image0(xx,yy,zz) = kD;
                            image1(kX,kY,kZ) = 1;
                        end
                    end
                end
            end
        end
%     max_pore_size = round(max(image0(:)));
%     snip = 1+max_pore_size:length(image0)-max_pore_size;
%     image02 = image0(snip,snip,snip);
%     image12 = image1(snip,snip,snip);
%     volview(1-image12,size(image12))
%     pause
    end
end

% max_pore_size = round(max(image0(:)));
% snip = 1+max_pore_size:length(image0)-max_pore_size;
% image02 = image0(snip,snip,snip);
% image12 = image1(snip,snip,snip);

% figure;
% plot_surface_cube(image0, size(image0))
% 
% 
% % % 
% figure;
% volview(1-image0,size(image0))
% hold on
