function [detectedFeatures,clustersMMF,imageN3,errFlag] = ...
    detectSubResFeatures2D(image,cands,psfSigma,testAlpha,visual,...
    doMMF,bitDepth,saveResults,bgNoiseSigma)
%DETECTSUBRESFEATURES2D determines the positions and intensity amplitudes of sub-resolution features using mixture model fitting
%
%SYNOPSIS [detectedFeatures,clustersMMF,imageN3,errFlag] = ...
%    detectSubResFeatures2D(image,cands,psfSigma,testAlpha,visual,...
%    doMMF,bitDepth,saveResults,bgNoiseSigma)
%
%INPUT  image      : Image being analyzed.
%       cands      : Cands structure as output from fsmCenter.
%       psfSigma   : Standard deviation of point spread function (in pixels).
%       testAlpha  : Alpha-values for statistical tests. Structure with fields:
%             .alphaR: For the residuals test, comparing N+1-kernel fit to
%                      N-kernal fit. Optional. Default: 0.05.
%             .alphaA: For amplitude test. Optional. Default: 0.05.
%             .alphaD: For distance test. Optional. Default: 0.05.
%             .alphaF: Final residuals test, comparing residuals from final
%                      fit to estimated background noise.
%                      Optional. Default: 0.05.
%       visual     : 1 if user wants to view results; 0 otherwise. 
%                    Optional. Default: 0.
%       doMMF      : 1 if user wants to do mixture-model fitting, 0
%                    otherwise. Optional. Default: 1.
%       bitDepth   : Camera bit depth. Optional. Default: 14.
%       saveResults: 1 if results are to be saved (in file 'detectedFeatures.mat'),
%                    0 otherwise. Optional. Default: 0.
%       bgNoiseSigma:Standard deviation of background noise. Optional. If
%                    not input, the code will estimate it from the image.
%
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT detectedFeatures: Structure with fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       clustersMMF: Array of clusters of sub-resolution features. 
%                    Structure with fields:
%             .position  : Position of each feature in image coordinate 
%                          system (in pixels): [x y dx dy] = [y x dy dx] 
%                          in matrix coordinate system.
%             .amplitude : Intensity of each feature [A dA].
%             .bgAmp     : Background intensity [Abg dAbg].
%             .varCovMat : Variance-covariance matrix of estimated parameters.
%       imageN3    : Image with labeled features. Blue: those from cands;
%                    Red: those from mixture-model fitting; Magenta: those 
%                    from MMF which coincide with those from cands.
%                    Will be output only if visual = 1.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%This file is part of u-track.
%
%    u-track is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%    u-track is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with u-track.  If not, see <http://www.gnu.org/licenses/>.
%
%Copyright: Jaqaman, 08/05

%% Output

detectedFeatures = [];
clustersMMF = [];
imageN3 = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--detectSubResFeatures2D: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether statistical test alpha values were inputted
if nargin < 4 || isempty(testAlpha) %if not, assign defaults
    
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0.05);
    
else %if some were, check their values and assign default for the rest

    if ~isfield(testAlpha,'alphaR')
        testAlpha.alphaR = 0.05;
    else
        if testAlpha.alphaR < 0 || testAlpha.alphaR > 1
            disp('--detectSubResFeatures2D: testAlpha.alphaR should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaA')
        testAlpha.alphaA = 0.05;
    else
        if testAlpha.alphaA < 0 || testAlpha.alphaA > 1
            disp('--detectSubResFeatures2D: testAlpha.alphaA should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaD')
        testAlpha.alphaD = 0.05;
    else
        if testAlpha.alphaD < 0 || testAlpha.alphaD > 1
            disp('--detectSubResFeatures2D: testAlpha.alphaD should be between 0 and 1!');
            errFlag = 1;
        end
    end
    if ~isfield(testAlpha,'alphaF')
        testAlpha.alphaF = 0.05;
    else
        if testAlpha.alphaF < 0 || testAlpha.alphaF > 1
            disp('--detectSubResFeatures2D: testAlpha.alphaF should be between 0 and 1!');
            errFlag = 1;
        end
    end

end

%check visualization option
if nargin < 5 || isempty(visual)
    visual = 0;
else
    if visual ~= 0 && visual ~= 1
        disp('--detectSubResFeatures2D: Variable "visual" should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether to do MMF
if nargin < 6 || isempty(doMMF)
    doMMF = 1;
else
    if doMMF ~= 0 && doMMF ~= 1
        disp('--detectSubResFeatures2D: Variable "doMMF" should be 0 or 1!');
        errFlag = 1;
    end
end

%check the bit depth
if nargin < 7 || isempty(bitDepth)
    bitDepth = 14;
else
    if bitDepth <= 0 || bitDepth-floor(bitDepth) ~= 0
        disp('--detectSubResFeatures2D: Variable "bitDepth" should be a positive integer!');
    end
end

%check whether results are to be saved
if nargin < 8 || isempty(saveResults)
    saveResults = 0;
else
    if saveResults ~= 0 && saveResults ~= 1
        disp('--detectSubResFeatures2D: Variable "saveResults" should be 0 or 1!');
    end
end

%check whether back ground noise sigma is input or whether it should be
%calculated on the fly
if nargin < 9 || isempty(bgNoiseSigma)
    bgNoiseSigma = 0;
    estimateBgNoise = 1;
else
    estimateBgNoise = 0;
end

%exit if there are problems with input data
if errFlag
    disp('--detectSubResFeatures2D: Please fix input data!');
    return
end

%get number of pixels in each direction (in image coordinate system)
[numPixelsY,numPixelsX] = size(image);

%extract test alpha values from input
alphaR = testAlpha.alphaR;
alphaA = testAlpha.alphaA;
alphaD = testAlpha.alphaD;
alphaF = testAlpha.alphaF;

%Divide image by bit depth, to normalize it between 0 and 1
image = double(image)/(2^bitDepth-1);

%get background intensity information from cands
bgAmp = vertcat(cands.IBkg);
status = vertcat(cands.status);
bgAmp = bgAmp(status==1);
bgAmpMax = max(bgAmp);
bgAmpAve = mean(bgAmp);

%% Determine signal overlap

%determine which signals are overlapping, in which case they must
%be fitted together later on
[clusters,errFlag] = findOverlapPSFs2D(cands,numPixelsX,numPixelsY,psfSigma);
if errFlag
    disp('--detectSubResFeatures2D: Could not place signals in clusters!');
    return
end

%% Mixture Model Fitting

%initialize vector indicating whether clusters should be retained
keepCluster = ones(length(clusters),1);

%set optimization options
options = optimset('Jacobian','on','Display','off');

%reserve memory for clustersMMF
numClusters = length(clusters);
clustersMMF = repmat(struct('position',[],'amplitude',[],'bgAmp',[],...
    'numDegFree',[],'residuals',[]),numClusters,1);

%go over all clusters
for i = 1 : numClusters

    %get initial guess of positions and amplitudes
    numMaximaT = clusters(i).numMaxima;
    maximaPosT = clusters(i).maximaPos(:,1:2);
    maximaAmpT = clusters(i).maximaAmp;
    bgAmpT = bgAmpAve;
    clusterPixels = clusters(i).pixels(:,1:2);

    %crop part of image that is relevant to this fitting
    imageC = image(clusters(i).pixels(:,3));

    firstFit = 1; %logical variable indicating first fit attempt
    fit = 1; %logical variable indicating whether to attempt to fit
    while fit

        %collect initial guesses and lower and upper bounds of ...
        %feature positions
        x0 = maximaPosT; %initial guess
        lb = repmat(min(clusterPixels),numMaximaT,1); %lower bound
        ub = repmat(max(clusterPixels),numMaximaT,1); %upper bound
        %feature amplitudes
        x0 = [x0 maximaAmpT];
        lb(:,3) = 1e-5;
        ub(:,3) = 2*x0(:,3);
        %background intensity
        x0 = x0';
        x0 = [x0(:); bgAmpT];
        lb = lb';
        lb = [lb(:); 1e-5];
        ub = ub';
        ub = [ub(:); 2*bgAmpMax];

        %calculate number of degrees of freedom in system
        numDegFreeT = size(clusterPixels,1)-3*numMaximaT-1;

        %determine feature positions and amplitudes, and estimate background
        %intensity, using nonlinear least squares data fitting
        solutionT = lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
            clusterPixels,psfSigma);

        %get residuals and Jacobian matrix to calculate variance-covariance
        %of estimated parameters
        [residualsT,jacMatT] = fitNGaussians2D(solutionT,imageC,...
            clusterPixels,psfSigma);

        %check whether addition of 1 PSF has significantly improved the fit
        if firstFit %if this is the first fit

            firstFit = 0; %next one won't be
            
        else %if this is not the first fit

            %get test statistic, which is F-distributed
            testStat = (sum(residualsT.^2)/numDegFreeT)/...
                (sum(residuals.^2)/numDegFree);

            %get p-value of test statistic
            pValue = fcdf(testStat,numDegFreeT,numDegFree);

            %compare p-value to alpha
            %1-sided F-test: H0: F=1, H1: F<1
            if pValue < alphaR %if p-value is smaller, accept this fit
                fit = 1; %and attempt another one with an additional kernel
            else %if p-value is larger, do not accept this fit and exit
                fit = 0;
            end
        end

        if fit %if this fit is accepted (which is the default if it's the first fit)

            %update variables
            numMaxima = numMaximaT;
            numDegFree = numDegFreeT;
            solution = solutionT;
            residuals = residualsT;
            jacMat = jacMatT;

            %calculate the parameters' variance-covariance matrix and get their
            %uncertainties
            varCovMat = inv(jacMat'*jacMat)*sum(residuals.^2)/numDegFree;
            standDevVec = sqrt(diag(varCovMat));
            
            %if nothing weird took place in the fit...
            if all(isreal(standDevVec))

                %extract estimate and std of background intensity and
                %remove from vectors
                bgAmp = [solution(end) standDevVec(end)];
                solution = solution(1:end-1);
                standDevVec = standDevVec(1:end-1);

                %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
                solution = reshape(solution,3,numMaxima);
                solution = solution';
                standDevVec = reshape(standDevVec,3,numMaxima);
                standDevVec = standDevVec';

                %extract feature positions and amplitudes and their uncertainties
                maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
                maximaAmp = [solution(:,3) standDevVec(:,3)];

                %check amplitudes and remove maxima with insignificant amplitudes
                %1-sided t-test: H0: T=0, H1: T>0
                %calculate test statistic (t-distributed)
                testStat = maximaAmp(:,1)./maximaAmp(:,2);
                %             testStat = maximaAmp(:,1)./(sqrt(maximaAmp(:,2).^2+bgAmp(2)^2));
                %             testStat = maximaAmp(:,1)./bgNoiseSigma;
                %get its p-value
                pValue = 1-tcdf(testStat,numDegFree);
                %find the maxima whose p-value < alpha
                indx = find(pValue<alphaA);

                if isempty(indx) %if none of the maxima have a significant amplitude

                    fit = 0; %stop fitting this cluster
                    keepCluster(i) = 0; %mark it to be discarded

                else %if there are maxima with significant amplitude

                    %keep only maxima with significant amplitudes
                    maximaPos = maximaPos(indx,:);
                    maximaAmp = maximaAmp(indx,:);
                    numMaxima = size(maximaAmp,1);

                    %update variance-covariance matrix as well
                    indx2 = 3*(indx-1)+1;
                    indx2 = [indx2 indx2+1 indx2+2];
                    indx2 = reshape(indx2',3*size(indx2,1),1);
                    varCovMat = varCovMat(indx2,:);
                    varCovMat = varCovMat(:,indx2);

                    %check distances between maxima if there is more than 1
                    %maximum in cluster
                    %1-sided t-test: H0: T=0, H1: T>0
                    if numMaxima > 1

                        %get p-values of distances between maxima
                        pValue = zeros(numMaxima);
                        for k=1:numMaxima-1
                            for j=k+1:numMaxima

                                %calculate distance between the 2 maxima
                                x1_x2 = maximaPos(j,1) - maximaPos(k,1);
                                y1_y2 = maximaPos(j,2) - maximaPos(k,2);
                                distance = sqrt(x1_x2^2+y1_y2^2);

                                %get the standard deviation in the distance
                                j1 = 3*(j-1)+1;
                                k1 = 3*(k-1)+1;
                                stdDist = x1_x2^2*(varCovMat(j1,j1) + ...
                                    varCovMat(k1,k1) - 2*varCovMat(j1,k1)) ...
                                    + y1_y2^2*(varCovMat(j1+1,j1+1) + ...
                                    varCovMat(k1+1,k1+1) - 2*varCovMat(j1+1,k1+1)) ...
                                    + 2*x1_x2*y1_y2*(varCovMat(j1,j1+1) - ...
                                    varCovMat(j1,k1+1) - varCovMat(j1+1,k1) + ...
                                    varCovMat(k1,k1+1));
                                stdDist = sqrt(stdDist)/distance;

                                %calculate test statistic (t-distributed)
                                testStat = distance/stdDist;

                                %get p-value
                                pValue(j,k) = 1-tcdf(testStat,numDegFree);

                            end
                        end

                        repTest = 1; %logical variable indicating whether distance test is performed
                        while repTest

                            %compare maximum p-value to alpha
                            maxPValue = max(pValue(:));

                            if maxPValue > alphaD %if the maximum p-value is not significant

                                %find pair with maximum p-value
                                [indx,indx2] = find(pValue==maxPValue);

                                %out of this pair, identify maximum with smaller amplitude
                                indx2 = [indx;indx2];
                                ampBoth = maximaAmp(indx2,1);
                                indx = ones(numMaxima,1);
                                indx(indx2(ampBoth==min(ampBoth))) = 0;

                                %remove that maximum
                                indx = find(indx);
                                maximaPos = maximaPos(indx,:);
                                maximaAmp = maximaAmp(indx,:);
                                numMaxima = numMaxima - 1;
                                pValue = pValue(indx,:);
                                pValue = pValue(:,indx);

                                %repeat test if there still is more than 1 maximum
                                if numMaxima > 1
                                    repTest = 1;
                                else
                                    repTest = 0;
                                end

                            else

                                repTest = 0; %do not repeat test if all distances are significant

                            end %(if maxPValue > alpha)

                        end %(while repTest)

                    end %(if numMaxima > 1)

                    %since the acceptance of this fit implies that another fit with
                    %an additional kernel will be attempted, add a kernel
                    %positioned at the pixel with maximum residual
                    numMaximaT = numMaxima + 1; %update number of maxima
                    maximaAmpT = [maximaAmp(:,1); mean(maximaAmp(:,1))]; %signal amplitude
                    coord = clusterPixels(residuals==max(residuals),:); %position of new kernel
                    maximaPosT = [maximaPos(:,1:2); coord]; %signal positions
                    bgAmpT = bgAmp(1); %background amplitude

                end %(if isempty(indx))
                
            else %if things went wrong in the fit ...
                
                %discard this cluster and don't attempt any more kernels
                fit = 0;
                keepCluster(i) = 0;
                
            end %()

        end %(if fit)

        %if user does not want mixture-model fitting, then don't
        %attempt another fit
        if ~doMMF
            fit = 0;
        end

    end %(while fit)

    %if this is a good cluster,
    %repeat the fit with the final number of kernels to get a final
    %estimate of the parameters and the residuals from the fit
    if keepCluster(i)

        %collect initial guesses and lower and upper bounds of ...
        %feature positions
        x0 = maximaPos(:,1:2); %initial guess
        lb = repmat(min(clusterPixels),numMaxima,1); %lower bound
        ub = repmat(max(clusterPixels),numMaxima,1); %upper bound
        %feature amplitudes
        x0 = [x0 maximaAmp(:,1)];
        lb(:,3) = 1e-5;
        ub(:,3) = 2*x0(:,3);
        %background intensity
        x0 = x0';
        x0 = [x0(:); bgAmp(1)];
        lb = lb';
        lb = [lb(:); 1e-5];
        ub = ub';
        ub = [ub(:); 2*bgAmpMax];

        %calculate number of degrees of freedom in system
        numDegFree = size(clusterPixels,1)-3*numMaxima-1;

        %determine feature positions and amplitudes, and estimate background
        %intensity, using nonlinear least squares data fitting
        solution = lsqnonlin(@fitNGaussians2D,x0,lb,ub,options,imageC,...
            clusterPixels,psfSigma);

        %get residuals and Jacobian matrix to calculate variance-covariance
        %of estimated parameters
        [residuals,jacMat] = fitNGaussians2D(solution,imageC,...
            clusterPixels,psfSigma);

        %calculate the parameters' variance-covariance matrix and get their
        %uncertainties
        varCovMat = inv(jacMat'*jacMat)*sum(residuals.^2)/numDegFree;
        standDevVec = sqrt(diag(varCovMat));

        %extract estimate and std of background intensity and
        %remove from vectors
        bgAmp = [solution(end) standDevVec(end)];
        solution = solution(1:end-1);
        standDevVec = standDevVec(1:end-1);

        %reshape 3nx1 vectors "solution" and "standDevVec" into nx3 matrices
        solution = reshape(solution,3,numMaxima);
        solution = solution';
        standDevVec = reshape(standDevVec,3,numMaxima);
        standDevVec = standDevVec';

        %extract feature positions and amplitudes and their uncertainties
        maximaPos = [solution(:,1:2) standDevVec(:,1:2)];
        maximaAmp = [solution(:,3) standDevVec(:,3)];

        %store solution in clustersMMF
        clustersMMF(i).position = maximaPos;
        clustersMMF(i).amplitude = maximaAmp;
        clustersMMF(i).bgAmp = bgAmp;
        clustersMMF(i).numDegFree = numDegFree;
        clustersMMF(i).residuals = residuals;
        
    end %(if keepCluster(i))

end %(for i=length(clusters):-1:1)

%find the clusters to retain
indx = find(keepCluster);

 %if there are clusters with significant signal ...
 if ~isempty(indx)

    %retain only those clusters which have significant signals in them
    clustersMMF = clustersMMF(indx);
    
    %re-initialize keepClusters to use in next test
    keepCluster = ones(length(indx),1);
    
    %estimate background noise variance if not supplied
    if estimateBgNoise
    else
        bgNoiseVar = bgNoiseSigma^2;
    end
    
    %go over all clusters and check that their residuals are
    %comparable to the background noise 
    for iCluster = 1 : length(clustersMMF)
        
        %get test statistic, which is F-distributed
        testStat = (sum(clustersMMF(iCluster).residuals.^2)/numDegFree)/...
            bgNoiseVar;

        %get p-value of test statistic
        pValue = 1 - fcdf(testStat,numDegFree,numDegFree+3*size(clustersMMF(iCluster).bgAmp,1));

        %compare p-value to alpha
        %1-sided F-test: H0: F=1, H1: F>1
        if pValue < alphaF %if p-value is smaller than alpha, reject overall fit
            keepCluster(iCluster) = 0;
        end

    end %(for iCluster = 1 : length(clustersMMF))
    
    %determine clusters that pass the residuals test
    indx = find(keepCluster);
    
    if ~isempty(indx)

        %keep only clusters that pass the residuals test
        clustersMMF = clustersMMF(indx);
        
        %store information in structure "detectedFeatures"
        tmp = vertcat(clustersMMF.position);
        detectedFeatures.xCoord = [tmp(:,1) tmp(:,3)];
        detectedFeatures.yCoord = [tmp(:,2) tmp(:,4)];
        detectedFeatures.amp = vertcat(clustersMMF.amplitude);

    end
    
end %(if ~isempty(indx))

if isempty(indx)

    %store empty structures
    clustersMMF = struct('position',zeros(0,4),'amplitude',zeros(0,2),'bgAmp',zeros(0,2));
    detectedFeatures = struct('xCoord',zeros(0,2),'yCoord',zeros(0,2),'amp',zeros(0,2));
    
end

%save output if requested
if saveResults
    save('detectedFeatures','detectedFeatures','clustersMMF');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visual

    %make 3 layers out of original image (normalized)
    imageNorm = image/max(image(:));
    imageN3 = repmat(imageNorm,[1 1 3]);

    %place zeros in pixels of maxima from cands
    for i=1:length(clusters)
        for j=1:3
            imageN3(clusters(i).maximaPos(:,3)+(j-1)*numPixelsX*numPixelsY)=0;
        end
    end

    %place zeros in pixels of maxima from mixture-model fitting
    for i=1:length(clustersMMF)
        pos = (round(clustersMMF(i).position(:,1)-1))*numPixelsY ...
            + round(clustersMMF(i).position(:,2));
        for j=1:3
            imageN3(pos+(j-1)*numPixelsX*numPixelsY)=0;
        end
    end

    %label maxima from cands in blue
    for i=1:length(clusters)
        imageN3(clusters(i).maximaPos(:,3)+2*numPixelsX*numPixelsY)=1;
    end

    %label maxima from mixture-model fitting in red
    %a maximum from mixture-model fitting that falls in the same pixel 
    %as that from cands will appear in magenta
    for i=1:length(clustersMMF)
        pos = (round(clustersMMF(i).position(:,1)-1))*numPixelsY ...
            + round(clustersMMF(i).position(:,2));
        imageN3(pos)=1;
    end

    %plot image
    imtool(imageN3);

end

%%%%% ~~ the end ~~ %%%%%


%                     %since the acceptance of this fit implies that another fit with
%                     %an additional kernel will be attempted, add one kernel to the
%                     %cluster next to the signal with largest amplitude
%                     numMaximaT = numMaxima + 1; %update number of maxima
%                     maximaAmpT = [maximaAmp(:,1); mean(maximaAmp(:,1))]; %signal amplitude
%                     tmp = find(maximaAmpT==max(maximaAmpT),1); %signal with largest amplitude
%                     coord = maximaPos(tmp,1:2) + (2*rand(1,2)-1)*2; %position of new kernel
%                     maximaPosT = [maximaPos(:,1:2); coord]; %signal positions
%                     bgAmpT = bgAmp(1); %background amplitude

