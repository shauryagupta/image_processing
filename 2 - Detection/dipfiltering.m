function points=dipfiltering(name)

% Get list of all TIF files in working directory and store in a structure
currentdir=pwd;
wdir=strcat([currentdir,'\',name,'\']);
imagefiles = dir(strcat([wdir,'*.tif']));      
num_images = length(imagefiles);    % Number of files founds
A=newimar(num_images); %Initializes array of images
%grayimages=newimar(num_images); %This creates an array of gray images if
%the original imageset is is colour

for ii=1:num_images
   currentfilename = imagefiles(ii).name; %Extracts the filename
   A{ii} = readim(strcat([wdir,currentfilename])); %Loads in the image with the given filename to array "A"
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
%This if loop allows us to change parameters based on the movie in
%question. I have chosen values which seem to work best for each movie.
if strcmpi(name,'crop4')==1
    thsize=9; %Size of top hat filter
    opensize=3; %Size of morphological opening
    largesize=50; %Max size of objects to keep
    smallsize=15; %Min size of objects to keep
elseif strcmpi(name,'greatmovie')==1
    thsize=11;
    opensize=5;
    largesize=50;
    smallsize=20;
elseif strcmpi(name,'crop1')==1
    thsize=9;
    opensize=3;
    largesize=40;
    smallsize=15;
else
    error('Invalid Folder. This code is only meant to run on crop1, crop4, or greatmovie.')
end
for i=1:num_images %This loop conducts the filtering
    C{i}=tophat(A{i},thsize,'elliptic'); %Removes low intensity background signal for objects greater than given size
    D{i}=stretch(C{i},1,99,0,255); %Adjusts the histogram to increase contrast
    E{i}=threshold(D{i},'otsu'); %Automatically thresholds remaining image
    F{i}=opening(E{i},opensize,'elliptic'); %Removes objects less than given size
end
for j=1:num_images
    msr=measure(F{j},[],{'Size','Minimum','Maximum','CartesianBox'},[]); %Extracts image information for each object in image
    lmsr=length(msr.size); %Measures the amount of objects detected
    %for i=1:lmsr ; %The next four lines provide an alternate measure of
    %"circularity", which dipimages measures through P2A on line 28. I
    %decided not to use it but I kept the code here just in case
    %circ(i)=(4*pi*msr.size(i))/(msr.perimeter(i))^2;
    %end
    %msr.circularity=circ;
    Fsmall{j}=F{j};
    indlarge=find(msr.size>largesize); %Finds indices of objects greater than largesize
    indsmall=find(msr.size<smallsize); %Find indicies of objects smaller than smallsize
    indcat=horzcat(indlarge,indsmall); %Concatenates vectors
    if isempty(indcat)==0 %If any large or small objects were detected, do the following. This line helps prevent indexing errors in case the indcat array is empty.
    msr_bigandsmall=msr(indcat); %Creates a new measurement array containing big and small objects
    numbig=size(msr_bigandsmall.ID,2); %Counts the number of big/small objects
     %Initializes a new array that will only contain objects smaller than
     %largesize but larger than smallsize
        for k=1:numbig %Sets the pixels containing the big/small objects to 0
            Fsmall{j}(msr_bigandsmall.Minimum(1,k):msr_bigandsmall.Maximum(1,k),msr_bigandsmall.Minimum(2,k):msr_bigandsmall.Maximum(2,k))=0;
        end
    end
    msrFsmall{j}=measure(Fsmall{j},[],{'Center','Size'},[]); %Measures properties of the remaining objects
    numsmall=size(msrFsmall{j}.ID,2); %Counts the number of objects
        for h=1:numsmall %Stores the centroids of all objects.
    points{1,j}(h,1)=(msrFsmall{j}.Center(1,h));
    points{1,j}(h,2)=(msrFsmall{j}.Center(2,h));
        end
end
%Plays the series of images defined by the variable in imshow
for j=1:num_images
imshow(dip_array(Fsmall{j}),[]);
pause(0.5)
end
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