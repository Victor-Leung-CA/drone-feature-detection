function [distance,thetaX,thetaY] = featureDetectionFunction(videoMatrix, depthMatrix)
%Drone detection using combination of SURF and BRISK detection algorithms

%---Input Parameters---%
% videoMatrix: 720x1280x3 rgb matrix
% depthMatrix: 720x1280 depth matrix, normalized to 0-150
%----------------------%

%---Output Parameters---%
% distance: distance between drone and camera [m], 1x1 double
% thetaX: angle [radians] between centroid of image frame to drone, positive value right of centroid, negative left of centroid
% thetaY: angle [radians] between centroid of image frame to drone, positive value is above centroid
%-----------------------%

%% Initial image filtering

%Convert video rgb matrix to greyscale
img_gray = rgb2gray(videoMatrix);

%% Feature detection

%SURF Features
I_points = detectSURFFeatures(img_gray);
[I_features, I_validPoints] = extractFeatures(img_gray,I_points);
%BRISK Features
I_points_BRISK = detectBRISKFeatures(img_gray);
[I_features_BRISK, I_validPoints_BRISK] = extractFeatures(img_gray,I_points_BRISK);

%% Import reference drone feature data

ref_features = load('ref_features.mat');
ref_features_BRISK = load('ref_features_BRISK.mat');
ref_validPoints = load('ref_validPoints.mat');

%% Match reference and video features

%Match SURF and BRISK descriptors
index_pairs = matchFeatures(ref_features, I_features,'MatchThreshold',20,'MaxRatio',0.8);
indexPairsBRISK = matchFeatures(ref_features_BRISK,I_features_BRISK,'MatchThreshold',20,'MaxRatio',0.85);

%Obtain matched points
%SURF
ref_matched_points = ref_validPoints(index_pairs(:,1));
I_matched_points = I_validPoints(index_pairs(:,2));
%BRISK
ref_matched_points_BRISK  = ref_validPoints_BRISK(indexPairsBRISK(:,1));
I_matched_points_BRISK = I_validPoints_BRISK(indexPairsBRISK(:,2));

%% SURF and BRISK Features matched
matchedOriginalXY  = [ref_matched_points.Location; ref_matched_points_BRISK.Location];
matchedVideoXY = [I_matched_points.Location; I_matched_points_BRISK.Location];

[tform_BOTH,~,I_inlier_points_BOTH] = estimateGeometricTransform(matchedOriginalXY,matchedVideoXY,'affine','MaxDistance',20,'Confidence',70);

%% Polygon outlines (BOTH)

dronePolygon_BOTH = [1, 1;...                       % top-left
        size(ref_img, 2), 1;...                     % top-right
        size(ref_img, 2), size(ref_img, 1);...      % bottom-right
        1, size(ref_img, 1);...                     % bottom-left
        1,1];                                       % top-left again to close the polygon

% Drone boundary transformed and rescaled
newDronePolygon_BOTH = transformPointsForward(tform_BOTH, dronePolygon_BOTH); 
%% Obtain centroid 

avgPoint_BOTH = mean(I_inlier_points_BOTH);

%Compute radius of circle
distSum = 0;
for i = 1:4
    distSum = distSum + ((newDronePolygon_BOTH(i,1)-newDronePolygon_BOTH(i+1,1))^2 + (newDronePolygon_BOTH(i,2)-newDronePolygon_BOTH(i+1,2))^2)^(1/2);
end
avgRadius = distSum/8;


%% Obtain depth of drone

%Need to open new figure for mask to be created
figure;
imshow(depthMatrix,[0,150]);

%Create mask
circle = drawcircle('Center', avgPoint_BOTH, 'Radius', avgRadius);
BW = createMask(circle);

%Calculate depth
filteredDroneDepth = depthMatrix .* BW;
droneDepthMask1 = filteredDroneDepth < 150;
droneDepthMask2 = filteredDroneDepth > 0;
filteredDroneDepth = filteredDroneDepth .* droneDepthMask1 .* droneDepthMask2;
depthtf =  filteredDroneDepth > 0; %Final mask of drone
distance = mean(reshape(filteredDroneDepth(depthtf),1,[]));

%% Obtain size of drone

%Obtain physical to pixel ratio
physicalDroneRadius = 0.25; %meters
physical2PixelRatio = physicalDroneRadius/avgRadius;

%Obtain physical translation parameters
physicalY = (avgPoint_BOTH(1) - depthInfoSize(1)/2) * physical2PixelRatio;
physicalX = (avgPoint_BOTH(1) - depthInfoSize(2)/2) * physical2PixelRatio;

%Obtain angles
thetaX = asin(physicalX/meanDepth);
thetaY = asin(physicalY/meanDepth);

end

