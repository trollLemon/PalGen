#include "kmeans.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <string.h>

#include "defaults.h"

/*
  cosineDistance

  Computes the cosine distance between two pixels.

  The cosine distance of two vectors is defined as:

   1 - Cosine Similarity(Vector A, Vector B), which is:

   1 - (A⋅B) / (||A|| * ||B||)
 */
double cosineDistance(const Pixel *a, const Pixel *b) {

  const double aDotb = (a->r * b->r) + (a->g * b->g) + (a->b * b->b);
  const double magA = sqrt((a->r * a->r) + (a->g * a->g) + (a->b * a->b));
  const double magB = sqrt((b->r * b->r) + (b->g * b->g) + (b->b * b->b));

  return 1.0 - (aDotb / (magA * magB));
}



bool converged(const Pixel *oldCentroids, const Pixel *newCentroids, const int numClusters) {

  bool converged = true;

  for (size_t i = 0; i < numClusters; i++) {
      const Pixel oldCentroid = oldCentroids[i];
      const Pixel newCentroid = newCentroids[i];

    if (!(oldCentroid.r == newCentroid.r && oldCentroid.g == newCentroid.g && oldCentroid.b == newCentroid.b)) {
      converged = false;
      break;
    }
  }

  return converged;
}

bool isChosen(const Point *point, Point **chosenCentroids,
              const int numClusters) {

  for (size_t i = 0; i < numClusters; i++) {
    const Point *centroid = chosenCentroids[i];
    if (point == centroid) {
      return 1;
    }
  }

  return 0;
}

/*
   initializeClusters
   randomly selects n points from the points array as the initial centroids, and
   creates n Clusters. Each cluster is placed into the clusters array. Points
   selected as centroids are placed into the chosenCentroids array.
 */
int initializeClusters(Cluster **clusters, Point **points,
                       Point **chosenCentroids, const int numClusters,
                       const int numPoints) {

  srand(
      time(NULL)); // this is a basic seed and not optimal for complex tasks
                   // like cryptography, but for simple shuffling this is fine.

  // Fish Yates Shuffle

  for (size_t i = numPoints - 1; i > 0; i--) {

    const size_t index = rand() % (i + 1);

    Point *tempPoint = points[i];
    points[i] = points[index];
    points[index] = tempPoint;
  }

  for (size_t i = 0; i < numClusters; i++) {

    Point *point = points[i];

    Cluster *cluster = malloc(sizeof(Cluster));

    cluster->points = malloc(sizeof(Point *) * STARTING_ARRAY_SIZE);
    cluster->_arrSize = STARTING_ARRAY_SIZE;
    cluster->nPoints = 0;

    cluster->centroid.r = point->pixel.r;
    cluster->centroid.g = point->pixel.g;
    cluster->centroid.b = point->pixel.b;

    chosenCentroids[i] = point;
    clusters[i] = cluster;
  }

  return 0;
}

void firstStep(Cluster **clusters, Point **points, Point **chosenCentroids,
               int numClusters, int numPoints) {

  for (size_t i = 0; i < numPoints; i++) {

    Point *point = points[i];

    // check if this point is marked as a centroid for any cluster.
    // if it is, we should skip it.
    // TODO: use a set. This function is O(n) time, ideally this should be O(1).
    if (isChosen(point, chosenCentroids, numClusters)) {
      continue;
    }

    size_t bestClusterIndex = 0;
    double shortestDistance = INFINITY;
    for (size_t j = 0; j < numClusters; j++) {

      const Pixel *currPixel = &point->pixel;
      const Pixel *currCentroid = &clusters[j]->centroid;
      const double distance = cosineDistance(currPixel, currCentroid);

      if (distance < shortestDistance) {
        bestClusterIndex = j;
        shortestDistance = distance;
      }
    }

    Cluster *closestCluster = clusters[bestClusterIndex];
    int result = clusterAddPoint(closestCluster, point);
  }

  // update centroids
  for (size_t i = 0; i < numClusters; i++) {
    calculateNewMean(clusters[i]);
    clusterResetPoints(clusters[i]);
  }
}

void refineClusters(Cluster **clusters, Point **points, const int numClusters,
                    const int numPoints) {

  Pixel *oldCentroids = malloc(sizeof(Pixel) * numClusters);
  Pixel *newCentroids = malloc(sizeof(Pixel) * numClusters);

  if (!oldCentroids) {
    return;
  }
  if (!newCentroids) {
    return;
  }

  bool shouldIterate = 1;
  while (shouldIterate) {


    for (size_t i = 0; i < numClusters; i++) {
      oldCentroids[i] = points[i]->pixel;
    }

    for (size_t i = 0; i < numPoints; i++) {

      Point *point = points[i];

      size_t bestClusterIndex = 0;
      double shortestDistance = INFINITY;
      for (size_t j = 0; j < numClusters; j++) {

        const Pixel *currPixel = &point->pixel;
        const Pixel *currCentroid = &clusters[j]->centroid;
        const double distance = cosineDistance(currPixel, currCentroid);

        if (distance < shortestDistance) {
          bestClusterIndex = j;
          shortestDistance = distance;
        }
      }

      Cluster *closestCluster = clusters[bestClusterIndex];
      int result = clusterAddPoint(closestCluster, point);
    }

    // update centroids
    for (size_t i = 0; i < numClusters; i++) {
      calculateNewMean(clusters[i]);
      clusterResetPoints(clusters[i]);
    }

    for (size_t i = 0; i < numClusters; i++) {
      newCentroids[i] = points[i]->pixel;
    }

    shouldIterate = !converged(oldCentroids, newCentroids, numClusters);


  }

  free(oldCentroids);
  free(newCentroids);

}


void clusterResetPoints(Cluster *cluster) {
  //lazily reset point counter. We will overwrite the old pointer values in the array next iteration
  cluster->nPoints = 0;
}

int clusterAddPoint(Cluster *cluster, Point *point) {

  // double array size if full
  // TODO: grow this exponentially?
  if (cluster->nPoints == cluster->_arrSize) {
    const size_t newSize = cluster->_arrSize * 2;

    Point **newArray = malloc(sizeof(Point *) * newSize);

    if (newArray == NULL) {
      return 0;
    }

    memcpy(newArray, cluster->points, sizeof(Point *) * cluster->_arrSize);
    free(cluster->points);
    cluster->_arrSize = newSize;
    cluster->points = newArray;
  }

  cluster->points[cluster->nPoints] = point;
  cluster->nPoints++;

  return 1;
}

int calculateNewMean(Cluster *cluster) {

  if (cluster == NULL)
    return -1;

  Point **points = cluster->points;

  Pixel *centroid = &cluster->centroid;

  const int nPoints = cluster->nPoints;

  double sumR = 0;
  double sumG = 0;
  double sumB = 0;

  for (size_t i = 0; i < nPoints; i++) {
    const Point *point = points[i];
    sumR += point->pixel.r;
    sumG += point->pixel.g;
    sumB += point->pixel.b;
  }

  const double meanR = sumR / nPoints;
  const double meanG = sumG / nPoints;
  const double meanB = sumB / nPoints;

  centroid->r = (unsigned char)round(meanR);
  centroid->g = (unsigned char)round(meanG);
  centroid->b = (unsigned char)round(meanB);

  return 0;
}



Pixel* generatePalette(Point** points, const int numPoints, const int numClusters) {

  Pixel* palette = malloc(sizeof(Pixel) * numClusters);
  Cluster** clusters = malloc(sizeof(Cluster*) * numClusters);
  Point** chosenCentroids = malloc(sizeof(Point*) * numClusters);
  if (!palette || !clusters || !chosenCentroids) {
    return NULL;
  }

  initializeClusters(clusters,points,chosenCentroids,numClusters,numPoints);
  firstStep(clusters,points,chosenCentroids,numClusters,numPoints);
  free(chosenCentroids);

  refineClusters(clusters,points,numClusters,numPoints);

  for (size_t i = 0; i < numClusters; i++) {
    palette[i] = clusters[i]->centroid;
    free(clusters[i]->points);
    free(clusters[i]);
  }
  free(clusters);

  return palette;

}