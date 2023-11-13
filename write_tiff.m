

fname = 'diamond_30_030.tif';
for n = 1:size(binary, 3)
    
    % Make an RGB image:
    i_img = binary(:, :, n);
    
    % Check what you are writing
    imshow(i_img);
    drawnow;
    
    % Generate your tiff stack:
  
        if n == 1
            % First slice:
            imwrite(i_img,fname)
        else
            % Subsequent slices:
            imwrite(i_img,fname,'WriteMode','append');
        end 

    
    disp(n)
end