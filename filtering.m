% Get list of all TIF files in working directory
imagefiles = dir('*.tif');      
num_images = length(imagefiles);    % Number of files found
rawimages=cell(1,num_images);
grayimages=cell(1,num_images);

for ii=1:num_images
   currentfilename = imagefiles(ii).name;
   rawimages{ii} = imread(currentfilename);
   grayimages{ii} = rgb2gray(rawimages{ii});
end
se=strel('disk',12);
se2=strel('disk',5);
A=cell(1,num_images);
B=cell(1,num_images);
C=cell(1,num_images);
D=cell(1,num_images);
E=cell(1,num_images);
F=cell(1,num_images);
for i=1:num_images
    %A{i}=imresize(images{i},2);
    C{i}=imtophat(grayimages{i},se);
    D{i}=imadjust(C{i});
    E{i}=imbinarize(D{i});
    F{i}=imopen(E{i},se2);
    %T=adaptthresh(A{i},0.3,'ForegroundPolarity','bright');
%B{i}=imbinarize(A{i},T);
%C{i}=imclose(B{i},se);
end
%Plays the series of images defined by the variable in imshow
for j=1:num_images
imshow(F{j})
pause(0.5)
end
%From here you can filter it with whatever you want. For example:
%%Top Hat Filtering
%Define a structuring element
%se=strel('disk',4);
%Filter images
%for k=1:num_images
%C{k}=imtophat(B{k},se);
%end
%for j=1:num_images
%imshow(B{j})
%pause(0.5)
%end
%subplot(2,3,1), imshow(rawimages{1}); title('Raw', 'FontSize',16);
%subplot(2,3,2), imshow(grayimages{1}); title('Grayscale', 'FontSize',16);
%subplot(2,3,3), imshow(C{1}); title('Top Hat Filtering', 'FontSize',16);
%subplot(2,3,4), imshow(D{1});  title('Contrast Adjustment', 'FontSize',16);
%subplot(2,3,5), imshow(imbinarize(D{1})); title('Binarization', 'FontSize',16);
%subplot(2,3,6), imshow(imopen(imbinarize(D{1}),se2)); title('Morphological Opening', 'FontSize',16)