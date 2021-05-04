# Drone-Feature-Detection

An image feature detection script combining SURF and BRISK feature detection algorithms to estimate pose of drone from a 2D image.

## **How it works**

A reference image is provided and SURF and BRISK features of the reference image are analyzed. BRISK and SURF features are chosen due to its scale and rotation invariance. An additional K-means color segmentation filter is added, as testing showed that it provided better results with obstacles in the image.

![reference image of drone](https://github.com/Victor-Leung-CA/Drone-Feature-Detection/blob/main/refimg.png)

Snapshots are taken periodically from a video capturing device, and both SURF and BRISK features are analyzed. The data points of the reference and snapshot image is then combined and filtered out for the inliers using the M-Estimator Sample Consensus (MSAC) algorithm.

![snapshot image](https://github.com/Victor-Leung-CA/Drone-Feature-Detection/blob/main/rgb45v2.jpg)
![matched points](https://github.com/Victor-Leung-CA/Drone-Feature-Detection/blob/main/BOTH_matched_points.jpg)

After matching the inliers, a boundary around the drone is created using the transformation matrix, which is derived from the transformation of the inlier points.
To find the depth and estimate the relative coordinate of the drone, the boundary circle is used to create a mask, and is applied to the depth matrix provided by the depth sensor. This is coupled with a screening of depths smaller than 150m.

![drone boundary](https://github.com/Victor-Leung-CA/Drone-Feature-Detection/blob/main/Matched_points_comparison.jpg)

Using the filtered depth matrix and the approximation of the drone as a sphere with a radius of 0.305 meters, we can derive an estimation of the vector between the camera and the drone. A world point to pixel mapping is obtained for the plane intersecting the drone, and is used to approximate the distance between the centroid of the image to the drone on the same plane.
