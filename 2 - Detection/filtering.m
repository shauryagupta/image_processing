% Get list of all TIF files in working directory
imagefiles = dir('*.tif');      
num_images = length(imagefiles);    % Number of files found
rawimages = cell(1,num_images);
grayimages = cell(1,num_images);

% Store each frame in "rawimages", store each frame in grayscale in 
% "grayimages"
for ii=1:num_images
   currentfilename = imagefiles(ii).name;
   rawimages{ii} = imread(currentfilename);
   grayimages{ii} = rgb2gray(rawimages{ii});
end

% Create structuring elements
se = strel('disk',10); %For top hat filtering
se2 = strel('disk',5); %For image opening

% Initialize Cell arrays
A=cell(1,num_images);
B=cell(1,num_images);
C=cell(1,num_images);
D=cell(1,num_images);
E=cell(1,num_images);
F=cell(1,num_images);

% Perform top hat filtering, increase contrast, binarization, and 
% image opening. 
for i=1:num_images
    %A{i}=imresize(images{i},2);
    C{i}=imtophat(grayimages{i},se);
    D{i}=imadjust(C{i});
    
    % imbinarize not available on my computer (MATLAB 2015b)
    %E{i}=imbinarize(D{i});
    
    level = graythresh(D{i});
    E{i} = im2bw(D{i},level);
    
    F{i}=imopen(E{i},se2);
    %T=adaptthresh(A{i},0.3,'ForegroundPolarity','bright');
    
%B{i}=imbinarize(A{i},T);
%C{i}=imclose(B{i},se);
end

%% Plays the series of images defined by the variable in imshow, 
% here it is the final opened image set

for j=1:num_images
    imshow(F{j})
    pause(0.1)
    
    %Stores each frame in struct M2
    M2(j)=getframe;
end

for idx = 1:num_images
    im{idx}=frame2im(M2(idx));
end
filename = 'out.gif'; % Specify the output file name
for idx = 1:num_images
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

%% See the gif created

[A,map]=imread('out.gif','frames','all'); 
mov=immovie(A,map); 
implay(mov)

%% This is what I used to create that montage for the progress report
%subplot(2,3,1), imshow(rawimages{1}); title('Raw', 'FontSize',16);
%subplot(2,3,2), imshow(grayimages{1}); title('Grayscale', 'FontSize',16);
%subplot(2,3,3), imshow(C{1}); title('Top Hat Filtering', 'FontSize',16);
%subplot(2,3,4), imshow(D{1});  title('Contrast Adjustment', ...
% 'FontSize',16);
%subplot(2,3,5), imshow(imbinarize(D{1})); title('Binarization', ...
% 'FontSize',16);
%subplot(2,3,6), imshow(imopen(imbinarize(D{1}),se2)); ...
% title('Morphological Opening', 'FontSize',16)
