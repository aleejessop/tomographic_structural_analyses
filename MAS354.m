clear all

%% 1. First we need to load the 16-bit greyscale data into the workspace
filename = 'C:\Users\20220428\OneDrive - Murdoch University\Documents\MATLAB\Sea_Urchin\data\greyscale\Sea Urchin Diamond_111.tiff';
greyscale = tiffreadVolume(filename);
volshow(greyscale)
imshow(squeeze(greyscale(1,:,:)))

%% 2. In order to use these data for some analyses we need to convert it to a binary dataset.
% To do this we will first examine the distribution of greyscale values by
% plotting them in a histogram. We will also define the histogram object as
% a variable 'h' so that we can inspect later.

figure, h=histogram(greyscale);

% We now need to decide which values we will convert to ones and which
% to zeros. The distribution of greyscale values is clearly bimodal so 
% we can use the local minimum between the two modes as a threshold value.
% We can see from the histogram the local minimum is around 31800.
% If we inspect the histogram object we have a few different variables we
% can look at. We will be looking at h.Values, i.e. the frequency of the
% greyscale values and h.BinEdges, i.e. the greyscale values. To find the
% local minimum we will use the 'findpeaks' function in MATLAB but we'll 
% use the minus operator to begate the frequency values in h.Values.

[~, valleyIndexes] = findpeaks(-h.Values); 

%It's found a seven different local minima and we now have the index values
%of these. Using the index values we can look at what the frequency values
%of these local minima are.

h.Values(valleyIndexes)

%Three of these have a frequency of zero so we can ignore these. One has a
%frequency of three and we know the local minimum we are after has a much
%higher frequency than three. The next lowest (9457) is the local minimum we want.

%To find the greyscale value that corresponds this we can simply index the
%h.BinEdges to see what the greyscale values are that correspond to the
%local minimum bin edges.

h.BinEdges(valleyIndexes)

%The value is 31400 and the bin width is 200, so let's choose a greyscale
%value between 31400 and 31600, meaning 31500 as our threshold.

%% 3. Thresholding the greyscale values to generate a binary dataset
%At the moment our greyscale data is a uint16 data type which doesn't 
%work very well with most MATLAB functions. So we will first convert 
%it to numeric data type by using the function 'double'.

greyscale = double(greyscale);

%Then we will use the 'imbinarize' function in MATLAB to threshold the
%greyscale values and convert them to ones and zeros.

binary = imbinarize(greyscale,31500);

%% 4. Let's plot the binary data to see what it looks like
%We can look at the dataset in three-dimensions by using the 'volshow'
%function in MATLAB.

volshow(binary)

%I want to look at it slice by slice so I've made a loop that displays the
%the x/y face of the cube for each slice z.

s=size(binary);

figure; %2D slices of image data
for i = 1:s(3)
    clf;
    imagesc(squeeze(binary(:,:,i)),'AlphaData',0.5); axis equal tight
    colormap gray
    colorbar
    drawnow;
    hold on
    pause
end


%% 5. Components of the dataset
%We are really only interested in the largest connected component of this
%data, i.e. the voxels that make up the curved network-like structure that
%looks a bit like a single Diamond surface. But you can see that there are
%also lots of components that have a value of one that are not connected to
%this surface. We need to remove these otherwise our analyses will not work
%properly. To find all the connected components within the binary data we
%can use the MATLAB function 'bwconncomp'. 

cc = bwconncomp(binary);

%Now we have a MATLAB structure called 'cc' which contains information
%each connected component. 'cc.PixelIdxList' contains the indices of each 
%voxel of binary in each connected component. For example the first column 
%in 'cc.PixelIdxList' happens to contain the indices for each voxel in the 
%largest connected component which is the component we are interested in.
%Using the 'cc.PixelIdxList' we can specify all the components that are not
%connected to the largest component and set their voxel value to zero.

%compute number of pixels in each connected component
numPixels = cellfun(@numel,cc.PixelIdxList);
%find components that are less than the largest connected component
unconnect_component_ID = find(numPixels<max(numPixels));
%set all components that are not connected to the largest connected
%component to have a voxel value of 0
for i = 1:numel(unconnect_component_ID)
    binary(cc.PixelIdxList{unconnect_component_ID(i)}) = 0;
end

%check the result using volshow

volshow(binary)

%Now we should check if there are unconnected components that have a value
%of zero, i.e. we need to check if there are 'holes' in our solid
%structure. We can do this by inversing the structure and making all ones
%zeros and all zeros ones. In MATLAB you can use the tilde (~) symbol to do
%this.

binary_inv = ~binary; %the inverse of binary equals NOT binary..
volshow(binary_inv)

%We can see a few unconnected components so we will repeat the above on
%binary_inv.

cc_inv = bwconncomp(binary_inv);
numPixels_inv = cellfun(@numel,cc_inv.PixelIdxList);
unconnect_component_ID_inv = find(numPixels_inv<max(numPixels_inv));
for i = 1:numel(unconnect_component_ID_inv)
    binary_inv(cc_inv.PixelIdxList{unconnect_component_ID_inv(i)}) = 0;
end

volshow(binary_inv)

%Now we inverse binary_inv once more to get our clean original dataset.

binary_clean = ~binary_inv;
volshow(binary_clean)

%% 8. Calculating solid volume fraction
%Calculating the solid volume fraction is simply a matter of calculating
%the volume of voxels within the solid phase, i.e. voxels with a value of
%one and dividing this by the overall sample size.

sample_size = numel(binary_clean); 
solid_phase = nnz(binary_clean);

solid_volume_fraction = solid_phase/sample_size;

%% 7. Calculating a pore size
%We will just take one slice of our dataset to start with.
%The Euclidean distance transform will tell us for each pixel in a slice of our newly
%cleaned binary dataset the distance between it and the nearest nonzero
%pixel. In MATLAB you can calculate this for each pixel using the function
%'bwdist'.

slice = binary_clean(:,:,1);
figure; imshow(slice)
D = bwdist(slice);

%D is now a matrix of distance values which we can take a look at by
%plotting them in a histrogram or plotting them as an image.

figure; histogram(D)
figure; imagesc(D)

%Our distance values go from 0 to about 14. Distances equal to 0 are the
%values that correspond to the ones of the binary dataset. The larger
%values of D are now approximately equal to the radii of each pore in pixels.

%Now we are going to ask the question: what is the largest circle that will
%fit inside this pore? And make a mapping of this by assigning the largest
%distance value in each pore to all other pixels in that pore that are
%within this largest distance.

%To do this we are going to loop through the Euclidean distances of
%the whole slice in a descending order. So the first loop will use the
%largest distance across the whole slice, second loop will 
%be the second largest, etc. We first need to make a vector containing all
%the distances in a descending order but we also want to keep track of the
%corresponding row and column indices for each distance.

DV = []; %this will be our new vector containing all distance values
I = []; %this will be the corresponding row value
J = []; %this will be the correspodning column value
k = 1; %this defines the row value of the new vectors
for i = 1:size(D,1)
    for j = 1:size(D,2)
        DV(k) = D(i,j);
        I(k) = i;
        J(k) = j;
        k = k+1;
    end
end

%Now we can sort the 'DV' vector using the 'sort' function and we'll set
%the outputs to be a new vector called DVS and the correspodning index
%values which we can use to then sort our original row and column indices.

[DVS, index] = sort(DV,'descend');
Isort = I(index);
Jsort = J(index);

%Now we will create our algorithm that will loop through all distance
%values within D in descending order. In the first loop kD will be equal to
%the largest distance in our data, kI will be the row value of this pixel,
%and kJ will be the column value. What we now need is a way to specify that
%we are only interested in the pixels that surround kD and we
%will do this by ensuring that the distance between the location of kD
%(specified by [kI,kJ]) and the location of each pixel is less than or
%equal to kD: sqrt((i-kI)^2 + (j-kJ)^2)<= kD. We will preallocate a zeros
%dataset and fill this dataset with our new pixel values so we will also
%specify that pixel value must also be 0 and within the distance of kD.
%i.e. it has not already had a distance assigned to it. This will
%result in circles that fill the pore space with a radius of kD.

image0 = zeros(size(D,1),size(D,2));

for k = 1:length(DVS)

    kD = (DVS(k)); %distance value
    kI = Isort(k); %row value of that pixel
    kJ = Jsort(k); %column value of that pixel

    for i = 1:size(D,1)
        for j = 1:size(D,2)
            if image0(i,j) == 0 && sqrt((i-kI)^2 + (j-kJ)^2)<= kD % if pixel image0(i,j) equals 0 and the distance between pixel image0(i,j) and pixel [kI,kJ] is less than kD
               image0(i,j) = kD; %pixel image0(i,j) equals kD
            end
        end
    end
    imagesc(image0)
end



















































