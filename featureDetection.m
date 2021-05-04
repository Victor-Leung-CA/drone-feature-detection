%% FEATURE DETECTION OPTIONS

segmentation = false; %Bool - 1: color segmentation + feature detection
mixedOptions = true; %Bool - 1: no color seg. for SURF and color seg. for BRISK
finalFigures = true; %Bool - 1: Only show final figures

%% Load image and compute features

ref_img = imread('refimg.png');
ref_img_gray = rgb2gray(ref_img);

if finalFigures ~= true
    imshow(ref_img_gray);
end

%% Reference Photo - Color filtering

%Convert image color to L*a*b color space
lab_ref_img = rgb2lab(ref_img);
ab_ref_img = lab_ref_img(:,:,2:3);
ab_ref_img = im2single(ab_ref_img);

%K-means color segmentation
nColors = 2;
[L,centers] = imsegkmeans(ref_img, nColors);
cls1 = (L == 2);
ref_img_filtered = uint8(cls1) .* (255-ref_img);
ref_img_filtered_gray = rgb2gray(ref_img_filtered);

if finalFigures ~= true
    montage({ref_img,ref_img_filtered})
end

%% Reference Photo - Detect image features

if segmentation
    if mixedOptions
        %SURF Features
        ref_points = detectSURFFeatures(ref_img_gray); %No color filter
        [ref_features, ref_validPoints] = extractFeatures(ref_img_gray,ref_points); %No color filter
        %BRISK Features
        ref_points_BRISK = detectBRISKFeatures(ref_img_filtered_gray);
        [ref_features_BRISK, ref_validPoints_BRISK] = extractFeatures(ref_img_filtered_gray,ref_points_BRISK);
    else
        %SURF Features
        ref_points = detectSURFFeatures(ref_img_filtered_gray); 
        [ref_features, ref_validPoints] = extractFeatures(ref_img_filtered_gray,ref_points); 
        %BRISK Features
        ref_points_BRISK = detectBRISKFeatures(ref_img_filtered_gray);
        [ref_features_BRISK, ref_validPoints_BRISK] = extractFeatures(ref_img_filtered_gray,ref_points_BRISK);
    end
else
    %SURF Features
    ref_points = detectSURFFeatures(ref_img_gray);
    [ref_features, ref_validPoints] = extractFeatures(ref_img_gray,ref_points);
    %BRISK Features
    ref_points_BRISK = detectBRISKFeatures(ref_img_gray);
    [ref_features_BRISK, ref_validPoints_BRISK] = extractFeatures(ref_img_filtered_gray,ref_points_BRISK);
end

% Display photo
if finalFigures ~= true
    figure; imshow(ref_img);
    hold on; plot(ref_points.selectStrongest(50));
end

%% Visualise strongest 25 SURF features

% figure;
% subplot(5,5,3); title('First 25 Features');
% for i=1:25
%     scale = ref_points(i).Scale;
%     image = imcrop(ref_img, [ref_points(i).Location-10*scale 20*scale 20*scale]);
%     subplot(5,5,i);
%     imshow(image);
%     hold on;
%     rectangle('Position',[5*scale 5*scale 10*scale 10*scale], 'Curvature', [1 1], 'EdgeColor', 'r')
% end

%% Video Frame - Capture image

%Frame number
frame = 60;

img = load('rgbinfo_v2.mat').ans.Data(:,:,:,frame);
img_gray = rgb2gray(img);

%% Video Frame - Color filtering 

%Convert image color to L*a*b color space
lab_img = rgb2lab(img);
ab_img = lab_img(:,:,2:3);
ab_img = im2single(ab_img);

%K means color segmentation
nColors2 = 3;
[L2,C2] = imsegkmeans(img, nColors2);
cls2 = (L2 == 1);
img_filtered = uint8(cls2) .* (255-img);
img_filtered_gray = rgb2gray(img_filtered);
img_bw = img_gray < 100;

if finalFigures ~= true
    montage({img,img_bw})
    title('Video Frame Color Filtered (K-means)');
end

%% Video Frame - Extract features

if segmentation
    if mixedOptions
        %SURF Features
        I_points = detectSURFFeatures(img_gray); %Edited - no color filter
        [I_features, I_validPoints] = extractFeatures(img_gray,I_points); %Edited - no color filter
        %BRISK Features
        I_points_BRISK = detectBRISKFeatures(img_filtered_gray);
        [I_features_BRISK, I_validPoints_BRISK] = extractFeatures(img_filtered_gray,I_points_BRISK);
    else
        %SURF Features
        I_points = detectSURFFeatures(img_filtered_gray); %Edited - no color filter
        [I_features, I_validPoints] = extractFeatures(img_filtered_gray,I_points); %Edited - no color filter
        %BRISK Features
        I_points_BRISK = detectBRISKFeatures(img_filtered_gray);
        [I_features_BRISK, I_validPoints_BRISK] = extractFeatures(img_filtered_gray,I_points_BRISK);
    end
else
    %SURF Features
    I_points = detectSURFFeatures(img_gray);
    [I_features, I_validPoints] = extractFeatures(img_gray,I_points);
    %BRISK Features
    I_points_BRISK = detectBRISKFeatures(img_gray);
    [I_features_BRISK, I_validPoints_BRISK] = extractFeatures(img_gray,I_points_BRISK);
end

%% Video Frame - Display SURF Features
if finalFigures ~= true
    figure; imshow(img);
    title('Video Frame - SURF Features');
    hold on; plot(I_points.selectStrongest(50));
end
%% Video Frame - Display BRISK Features
if finalFigures ~= true
    figure; imshow(img);
    title('Video Frame - BRISK Features');
    hold on; plot(I_points_BRISK.selectStrongest(50));
end
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

%% Visualise SURF matched points
if finalFigures ~= true
    figure, showMatchedFeatures(img, ref_img, I_matched_points, ref_matched_points, 'montage');
    title('Matched points - SURF');
    legend('Reference Points - SURF','Video Points - SURF')
end
%% Visualise BRISK matched points
if finalFigures ~= true
    figure
    showMatchedFeatures(ref_img,img,ref_matched_points_BRISK,I_matched_points_BRISK, 'montage')
    title('Putative matches using BRISK & FREAK')
    legend('Reference Points - BRISK','Video Points - BRISK')
end
%% Narrow down inliers (SURF)

[tform,ref_inlier_points, I_inlier_points] = estimateGeometricTransform(ref_matched_points,I_matched_points,'affine', 'MaxDistance',20,'Confidence',70 );

% Draw lines to matched points
if finalFigures ~= true
    figure
    showMatchedFeatures(img, ref_img, I_inlier_points, ref_inlier_points, 'montage');
    title('Matched inliers - SURF');
end

%% Narrow down inliers (BRISK)

[tform_BRISK,ref_inlier_points_BRISK, I_inlier_points_BRISK] = estimateGeometricTransform(ref_matched_points_BRISK,I_matched_points_BRISK,'affine', 'MaxDistance',50,'Confidence',40 );

% Draw lines to matched points
if finalFigures ~= true
    figure
    showMatchedFeatures(img, ref_img, I_inlier_points_BRISK, ref_inlier_points_BRISK, 'montage');
    title('Matched inliers - BRISK');
end

%% Polygon outlines (SURF)

dronePolygon = [1, 1;...                              % top-left
        size(ref_img, 2), 1;...                       % top-right
        size(ref_img, 2), size(ref_img, 1);...        % bottom-right
        1, size(ref_img, 1);...                       % bottom-left
        1,1];                                         % top-left again to close the polygon
newDronePolygon = transformPointsForward(tform, dronePolygon);

%% Polygon outlines (BRISK)

dronePolygon_BRISK = [1, 1;...                        % top-left
        size(ref_img, 2), 1;...                       % top-right
        size(ref_img, 2), size(ref_img, 1);...        % bottom-right
        1, size(ref_img, 1);...                       % bottom-left
        1,1];                                         % top-left again to close the polygon
newDronePolygon_BRISK = transformPointsForward(tform_BRISK, dronePolygon_BRISK);

%% SURF and BRISK Features matched

matchedOriginalXY  = [ref_matched_points.Location; ref_matched_points_BRISK.Location];
matchedVideoXY = [I_matched_points.Location; I_matched_points_BRISK.Location];

[tform_BOTH,ref_inlier_points_BOTH, I_inlier_points_BOTH] = estimateGeometricTransform(matchedOriginalXY,matchedVideoXY,'affine','MaxDistance',20,'Confidence',70);

if finalFigures ~= true
    figure
    showMatchedFeatures(ref_img,img,ref_inlier_points_BOTH,I_inlier_points_BOTH, 'montage')
    title('Matching points using SURF and BRISK (inliers only)')
    legend('ptsOriginal','ptsDistorted')
end

%% Polygon outlines (BOTH)

dronePolygon_BOTH = [1, 1;...                       % top-left
        size(ref_img, 2), 1;...                     % top-right
        size(ref_img, 2), size(ref_img, 1);...      % bottom-right
        1, size(ref_img, 1);...                     % bottom-left
        1,1];                                       % top-left again to close the polygon
newDronePolygon_BOTH = transformPointsForward(tform_BOTH, dronePolygon_BOTH);

%Display figure
figure; imshow(img);
title('Detected Drone - SURF, BRISK, BOTH comparisons')
hold on;
line(newDronePolygon(:, 1), newDronePolygon(:, 2), 'Color', 'g');
line(newDronePolygon_BRISK(:, 1), newDronePolygon_BRISK(:, 2), 'Color', 'r');
line(newDronePolygon_BOTH(:, 1), newDronePolygon_BOTH(:, 2), 'Color', 'b');
legend('SURF', 'BRISK', 'BOTH');

%% Obtain centroid (All)

avgPoint_SURF = mean(I_inlier_points.Location);
avgPoint_BRISK = mean(I_inlier_points_BRISK.Location);
avgPoint_BOTH = mean(I_inlier_points_BOTH);

%Compute radius of circle
%SURF
distSum = 0;
for i = 1:4
    distSum = distSum + ((newDronePolygon(i,1)-newDronePolygon(i+1,1))^2 + (newDronePolygon(i,2)-newDronePolygon(i+1,2))^2)^(1/2);
end
avgRadius_SURF = distSum/8;

%BRISK
distSum = 0;
for i = 1:4
    distSum = distSum + ((newDronePolygon_BRISK(i,1)-newDronePolygon_BRISK(i+1,1))^2 + (newDronePolygon_BRISK(i,2)-newDronePolygon_BRISK(i+1,2))^2)^(1/2);
end
avgRadius_BRISK = distSum/8;

%Both
distSum = 0;
for i = 1:4
    distSum = distSum + ((newDronePolygon_BOTH(i,1)-newDronePolygon_BOTH(i+1,1))^2 + (newDronePolygon_BOTH(i,2)-newDronePolygon_BOTH(i+1,2))^2)^(1/2);
end
avgRadius = distSum/8;

%Display figure
figure; imshow(img);
title('Centroid of matched points - SURF')
hold on;
circle1 = drawcircle('Center', avgPoint_SURF, 'Radius', avgRadius_SURF,'Label','SURF','Color','red','FaceAlpha',0.5);
circle2 = drawcircle('Center', avgPoint_BRISK, 'Radius', avgRadius_BRISK,'Label','BRISK','Color','green','FaceAlpha',0.5);
circle3 = drawcircle('Center', avgPoint_BOTH, 'Radius', avgRadius,'Label','BOTH','FaceAlpha',0.5);

%% Obtain depth of drone

%Import depthinfo matlab file
depthInfo = load('depthinfo_v2.mat').ans.Data(:,:,frame);
figure;
imshow(depthInfo,[0,150]);

%Rescaling - ignore for now
videoImageSize = size(img_gray);
depthInfoSize = size(depthInfo);
% scaleY = depthInfoSize(1)/videoImageSize(1);
% scaleX = depthInfoSize(2)/videoImageSize(2);
% circleScaled = drawcircle('Center', [avgPoint_BOTH(1)*scaleY, avgPoint_BOTH(2)*scaleX], 'Radius', avgRadius*(scaleX + scaleY)/2);

%Create mask
circle = drawcircle('Center', avgPoint_BOTH, 'Radius', avgRadius);
BW = createMask(circle);
imshow(BW);

%Calculate depth
filteredDroneDepth = depthInfo .* BW;
droneDepthMask1 = filteredDroneDepth < 150;
droneDepthMask2 = filteredDroneDepth > 0;
filteredDroneDepth = filteredDroneDepth .* droneDepthMask1 .* droneDepthMask2;
depthtf =  filteredDroneDepth > 0;
meanDepth = mean(reshape(filteredDroneDepth(depthtf),1,[]));

if ~finalFigures
    figure;
    imshow(filteredDroneDepth);
end
%% Obtain size of drone

%Obtain physical to pixel ratio
physicalDroneRadius = 0.305; %meters
physical2PixelRatio = physicalDroneRadius/avgRadius;

%Obtain physical translation parameters
physicalY = (avgPoint_BOTH(1) - depthInfoSize(1)/2) * physical2PixelRatio;
physicalX = (avgPoint_BOTH(1) - depthInfoSize(2)/2) * physical2PixelRatio;

%Obtain angles
thetaX = asin(physicalX/meanDepth)*180/pi;
thetaY = asin(physicalY/meanDepth)*180/pi;

fprintf("Frame %d: Measured Distance = %0.3d [m], X = %0.3d [m], Y = %0.3d [m] \n", frame, meanDepth, physicalX, physicalY);

%% Obtain world coordinates
pitch = -pi/4;
yaw = -pi/9;
offset = [0; 10; 50];

% Find point coordinates (x, y, z) in camera reference frame
physicalZ = sqrt(meanDepth^2 - physicalX^2 - physicalY^2);
point_cam = [physicalX; physicalY; physicalZ];

R_yaw = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
R_pitch = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
R = R_pitch * R_yaw;

point_world = R * point_cam + offset;

fprintf("X = %d [m], Y = %d [m], Z = %d [m]\n", point_world(1), point_world(2), point_world(3));