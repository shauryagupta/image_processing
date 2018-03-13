% Get list of all TIF files in working directory
imagefiles = dir('*.tif');      
num_images = length(imagefiles);    % Number of files founds
grayimages=newimar(num_images);

for ii=1:num_images
   currentfilename = imagefiles(ii).name;
   a = readim(currentfilename);
   grayimages{ii} = a{2};
end
C=newimar(num_images);
D=newimar(num_images);
E=newimar(num_images);
F=newimar(num_images);
for i=1:num_images
    C{i}=tophat(grayimages{i},12,'elliptic');
    D{i}=stretch(C{i},1,99,0,255);
    E{i}=threshold(D{i},'otsu');
    F{i}=opening(E{i},7,'elliptic');
end
%Plays the series of images defined by the variable in imshow
for j=1:num_images
imshow(dip_array(F{j}),[]);
pause(0.5)
end
%% OLD CODE.
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