%HW5
%GB comments:
1a 100
1b 100
1c 100	
1d 100
2yeast: 90. At the end, I used imerode to reduce the masks to get better separation of the objects. 
2worm: 100
2bacteria: 95 Could have imposed some morphological transformations to reduce the overlap of several masks. 
2phase: 100 
Overall: 99


% Note. You can use the code readIlastikFile.m provided in the repository to read the output from
% ilastik into MATLAB.

%% Problem 1. Starting with Ilastik

% Part 1. Use Ilastik to perform a segmentation of the image stemcells.tif
% in this folder. Be conservative about what you call background - i.e.
% don't mark something as background unless you are sure it is background.
% Output your mask into your repository. What is the main problem with your segmentation?  

%readIlastikFile
mask1 = readIlastikFile('stemcells_Simple Segmentation.h5');
imshow(mask1); % this mask is inversed

% inverted mask 
mask1Comp = imcomplement(mask1);
imshow(mask1Comp);
 
% Problem: mask does not seperate individual cells 

% Part 2. Read you segmentation mask from Part 1 into MATLAB and use
% whatever methods you can to try to improve it. 

% Apply watershed
CC = bwconncomp(mask1Comp);
state = regionprops(CC, 'Area');
area = [state.Area];
fusedCandidates = area > mean(area) + std(area); 
sublist = CC.PixelIdxList(fusedCandidates); 
sublist = cat(1,sublist{:});
fusedMask = false(size(mask1));
fusedMask(sublist) = 1; 
imshow(fusedMask, 'InitialMagnification', 'fit'); %show fused mask
s = round(0.8*sqrt(mean(area))/pi);
nucmin =imerode(fusedMask, strel('disk', s));
imshow(nucmin, 'InitialMagnification', 'fit'); 

outside = ~imdilate(fusedMask,strel('disk',1));
imshow(outside, 'InitialMagnification', 'fit');
basin =imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin|outside);
pcolor(basin); shading flat; 

L = watershed(basin);   
newmask = L >1 | (mask1Comp - fusedMask);
imshow(newmask, 'InitialMagnification', 'fit');
% Observe segmentation using label2rgb. 
cc = bwconncomp(newmask);
L = labelmatrix(cc);
rgb = label2rgb(L, 'jet',[.7 .7 .7],'shuffle');
imshow(rgb, 'InitialMagnification', 'fit');

%Ouput rgb into tif file. 
imwrite(rgb, 'HW5_Part1_Mask.tif')

% Part 3. Redo part 1 but now be more aggresive in defining the background.
% Try your best to use ilastik to separate cells that are touching. Output
% the resulting mask into the repository. What is the problem now?

%ReadIlastikFile
mask2 = readIlastikFile('stemcells_Simple Segmentation_aggresive2.h5');
figure; imshow(mask2);
mask2Comp = imcomplement(mask2);
imshow(mask2Comp);
% Problem: cells are well seperated, but the mask does not cover large area of the cells. 

% Part 4. Read your mask from Part 3 into MATLAB and try to improve
% it as best you can.

%I will dilate to make the mask bigger. 
mask2Comp_dil = imdilate(mask2Comp, strel('disk', 3));
imshow(mask2Comp_dil, []);

% use label2rgb to see segregation 
cc_mask2 = bwconncomp(mask2Comp_dil);
L_mask2 = labelmatrix(cc_mask2);
rgb_mask2 = label2rgb(L_mask2, 'jet',[.7 .7 .7],'shuffle');
imshow(rgb_mask2, 'InitialMagnification', 'fit');

%Output mask
imwrite(rgb_mask2, 'HW5_Part1_Mask2.tif');

%% Problem 2. Segmentation problems.

% The folder segmentationData has 4 very different images. Use
% whatever tools you like to try to segement the objects the best you can. Put your code and
% output masks in the repository. If you use Ilastik as an intermediate
% step put the output from ilastik in your repository as well as an .h5
% file. Put code here that will allow for viewing of each image together
% with your final segmentation. 

%%Read files 
yeast = imread('yeast.tif');
worm = imread('worms.tif');
bac = imread('bacteria.tif');

imshow(yeast,[])
figure; imshow(worm, [])
figure; imshow(bac, [])


%yeast
yeast_comp = imcomplement(yeast);
imshow(yeast_comp,[]);
mask = yeast_comp > 206; %206 biggest selection with a good border around the chunk of cell which I will need later on.  
imshow(mask); 
cell_prop = regionprops(mask, 'Area' ,'PixelIdxList', 'Image');
% Finding background area
hist([cell_prop.Area]);
ids = find([cell_prop.Area] > 10000); 
% turn background to zero
imshow(cell_prop(ids(1)).Image)
mask(cell_prop(ids(1)).PixelIdxList) = false; % the new mask has background removed. 
imshow(mask);
%clean-up mask
cleanMask = img_cleanup(mask, 3, 100);
imshow(cleanMask);
%dilate again
cleanMask2 = imdilate(cleanMask, strel('disk',3)); 
%cleanMask2 = imfill(cleanMask, 'fill'); imdilate looks better. 
imshow(cleanMask2); % mask for each  cells are now round and with no holes. 
%apply watershed
CC = bwconncomp(cleanMask2);
stats = regionprops(CC, 'Area');
area = [stats.Area]; % there are not many cells making statistics meaningless.
fusedCandidates = area> (mean(area) + std(area))/2 % divide by two to get all fused cells; 
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1, sublist{:}); 
fusedMask = false(size(mask));
fusedMask(sublist) = 1; 
imshow(fusedMask, 'InitialMagnification', 'fit')% this gives all the fused cells! 
title('fusedmask')
%eroding the centers
s = round(0.5*sqrt(mean(area))/pi); %0.5 works the best
nucmin = imerode(fusedMask,strel('disk',s));
figure; imshow(nucmin,'InitialMagnification', 'fit');
title('nucmin') 
%getting the regions outside
outside = ~imdilate(fusedMask, strel('disk',1));
imshow(outside, 'InitialMagnification', 'fit');
title('outside') 
%define basin for watershed
basin = imcomplement(bwdist(outside));
basin= imimposemin(basin, nucmin | outside); 
L = watershed(basin); 
imshow(L); colormap('jet'); caxis([0 20]);
%combining mask
newMask = L>1 | (cleanMask2 - fusedMask); 

%analyzing se
cc = bwconncomp(newMask);
L = labelmatrix(cc);
rgb = label2rgb(L, 'jet',[.7 .7 .7],'shuffle');
imshow(rgb, 'InitialMagnification', 'fit');

%Candid Analysis and Possible Improvment
%the image was difficult to segment because the cells are not bright and
%overlapping. 
%My final mask, however, did managed to segment most of the cells. 
%some cells are segmented as clumps(more than 1 cells); I can specfically
%target these cells and apply a watershed again only on them with optimzied
%parameters. 

%Overlaying two images
imshowpair(yeast,rgb, 'blend'); 


%Worm
worm = imread('worms.tif');
imshow(worm,[]);
%read ilastik file
immask = h5read('worms_Simple Segmentation_wormI.h5', '/exported_data');
immask = squeeze(immask);
data_trans = transpose(immask); 
data_trans(data_trans == 1) = 1;%1 is worm  
data_trans(data_trans == 2 )= 0;%2 is background
imshow(data_trans, []);


%Anlysis
cc1 = bwconncomp(data_trans);
L1 = labelmatrix(cc1);
rgb1 = label2rgb(L1, 'jet',[.7 .7 .7],'shuffle');
imshow(rgb1, 'InitialMagnification', 'fit');
imshowpair(worm,rgb1, 'blend'); 


%Bacteria
bac = imread('bacteria.tif');
imshow(bac, []); 
immask = h5read('bacteria_Simple Segmentation_BACI.h5', '/exported_data');
immask = squeeze(immask);

data_trans = transpose(immask); 
data_trans(data_trans == 1) = 1;%1 is worm  
data_trans(data_trans == 2 )= 0;%2 is background
imshow(data_trans, []);
bac_clean = img_cleanup(data_trans, 1, 100);

cc2 = bwconncomp(bac_clean);
L2 = labelmatrix(cc2);
rgb2 = label2rgb(L2, 'jet',[.7 .7 .7],'shuffle');
imshow(rgb2, 'InitialMagnification', 'fit');
imshowpair(bac,rgb2, 'blend'); 

%cellPhase
contrast = imread('cellPhaseContrast.png')
imshow(contrast, []);
immask = h5read('cellPhaseContrast_Simple Segmentation_contrastI.h5', '/exported_data');
immask = squeeze(immask);

data_trans = transpose(immask); 
data_trans(data_trans == 1) = 1;%1 is worm  
data_trans(data_trans == 2 )= 0;%2 is background
imshow(data_trans, [])
contrast_clean = img_cleanup(data_trans, 1, 100);

cc3 = bwconncomp(contrast_clean);
L3 = labelmatrix(cc3);
rgb3 = label2rgb(L3, 'jet',[.7 .7 .7],'shuffle');
imshow(rgb3, 'InitialMagnification', 'fit');
imshowpair(contrast,rgb3, 'blend'); 
