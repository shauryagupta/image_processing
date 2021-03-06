function points=dipfilteringcontrol_noisy2(name)

% Get list of all TIF files in working directory and store in a structure
currentdir=pwd;
wdir=strcat([currentdir,'/',name,'/']);
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
Cm=newimar(num_images);
Ct=newimar(num_images);
D=newimar(num_images);
E=newimar(num_images);
F=newimar(num_images);
Fsmall=newimar(num_images);
points=cell(num_images,1);
for i=1:num_images %This loop conducts the filtering
    Cm{i}=medif(A{i},8,'rectangular');
    Ct{i}=tophat(Cm{i},24,'rectangular'); %Removes low intensity background signal for objects greater than given size
    %D{i}=stretch(Ct{i},1,99,0,255); %Adjusts the histogram to increase contrast
    E{i}=threshold(Ct{i},'background'); %Automatically thresholds remaining image
    F{i}=opening(E{i},6,'rectangular'); %Removes objects less than given size
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
     Fsmall{j}=F{j}; %Initializes a new array that will only contain small objects
    ind=find(msr.size<118); %Finds indices of objects greater than size 50 pixels
    if isempty(ind)==0
    msr_big=msr(ind); %Creates a new measurement array containing only the big objects
    numbig=size(msr_big.ID,2); %Counts the number of big objects

        for k=1:numbig %Sets the pixels containing the big objects to 0
            Fsmall{j}(msr_big.Minimum(1,k):msr_big.Maximum(1,k),msr_big.Minimum(2,k):msr_big.Maximum(2,k))=0;
        end
    end
    msrFsmall{j}=measure(Fsmall{j},[],{'Minimum','Maximum','Size'},[]); %Measures properties of the remaining small objects
    numsmall=size(msrFsmall{j}.ID,2); %Counts the number of small objects
%     ind2=find(msrFsmall{j}.Minimum(2,:)==0);
%         if isempty(ind2)==0
%              msr_border{j}=msrFsmall{j}(ind2);
%              numborder=size(msr_border{j}.ID,2);
%                     for l=1:numborder %Sets the pixels containing the big objects to 0
%                         Fsmall{j}(msrFsmall{j}.Minimum(1,l):msrFsmall{j}.Maximum(1,l),msrFsmall{j}.Minimum(2,l):msrFsmall{j}.Maximum(2,l))=0;
%                     end
%             msrFsmall{j}=measure(Fsmall{j},[],{'Minimum','Maximum','Size'},[]); %Measures properties of the remaining small objects
%         end
    newnumsmall=size(msrFsmall{j}.ID,2);
        for h=1:newnumsmall %Stores the centroids of all small objects. NOTE that this assumes that all the objects are convex (in our case, they are)
    %points{1,j}(:,1)=msrFsmall{j}.ID;
    points{j,1}(h,1)=(msrFsmall{j}.Minimum(1,h)+msrFsmall{j}.Maximum(1,h))*0.5;
    points{j,1}(h,2)=(msrFsmall{j}.Minimum(2,h)+msrFsmall{j}.Maximum(2,h))*0.5;
        end
end

%Plays the series of images defined by the variable in imshow
  % for j=1:num_images
  % imshow(dip_array(Fsmall{j}),[]);
  % pause(0.01)
  % end
end
