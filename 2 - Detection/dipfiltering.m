function points = dipfiltering()
% Get list of all TIF files in working directory and store in a structure
imagefiles = dir('*.tif');
num_images = length(imagefiles);    % Number of files founds
A=newimar(num_images); %Initializes array of images
%grayimages=newimar(num_images); %This creates an array of gray images if
%the original imageset is is colour

for ii=1:num_images
   currentfilename = imagefiles(ii).name; %Extracts the filename
   A{ii} = readim(currentfilename); %Loads in the image with the given filename to array "A"
   %grayimages{ii} = a{2}; %This would extract the green layer of the
   %images if the original imageset is in colour
end
%Initializes arrays
C=newimar(num_images);
D=newimar(num_images);
E=newimar(num_images);
F=newimar(num_images);
Fsmall=newimar(num_images);
points=cell(1,num_images);

for i=1:num_images %This loop conducts the filtering
    C{i}=tophat(A{i},9,'elliptic'); %Removes low intensity background signal for objects greater than given size
    D{i}=stretch(C{i},1,99,0,255); %Adjusts the histogram to increase contrast
    E{i}=threshold(D{i},'otsu'); %Automatically thresholds remaining image
    F{i}=opening(E{i},3,'elliptic'); %Removes objects less than given size
end

for j=1:num_images
    msr=measure(F{j},[],{'Size','Perimeter','P2A','Minimum','Maximum'},[]); %Extracts image information for each object in image
    lmsr=length(msr.size); %Measures the amount of objects detected
    %for i=1:lmsr ; %The next four lines provide an alternate measure of
    %"circularity", which dipimages measures through P2A on line 28. I
    %decided not to use it but I kept the code here just in case
    %circ(i)=(4*pi*msr.size(i))/(msr.perimeter(i))^2;
    %end
    %msr.circularity=circ;
    ind=find(msr.size>50); %Finds indices of objects greater than size 50 pixels
    msr_big=msr(ind); %Creates a new measurement array containing only the big objects
    numbig=size(msr_big.ID,2); %Counts the number of big objects
    Fsmall{j}=F{j}; %Initializes a new array that will only contain small objects
        for k=1:numbig %Sets the pixels containing the big objects to 0
            Fsmall{j}(msr_big.Minimum(1,k):msr_big.Maximum(1,k),msr_big.Minimum(2,k):msr_big.Maximum(2,k))=0;
        end

    msrFsmall{j}=measure(Fsmall{j},[],{'Minimum','Maximum'},[]); %Measures properties of the remaining small objects
    numsmall=size(msrFsmall{j}.ID,2); %Counts the number of small objects

        for h=1:numsmall %Stores the centroids of all small objects. NOTE that this assumes that all the objects are convex (in our case, they are)
            points{1,j}(h,1)=(msrFsmall{j}.Minimum(1,h)+msrFsmall{j}.Maximum(1,h))*0.5;
            points{1,j}(h,2)=(msrFsmall{j}.Minimum(2,h)+msrFsmall{j}.Maximum(2,h))*0.5;
        end
end
%Plays the series of images defined by the variable in imshow
% for j=1:num_images
% imshow(dip_array(Fsmall{j}),[]);
% pause(0.5)
% end
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
end
