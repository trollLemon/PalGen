#include "kmeans.h"

#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>

#include "defaults.h"


typedef struct {

  double rgbAccum[3];
  size_t numPoints;

} AccumulatedRGBSums;

typedef struct {

  AccumulatedRGBSums *rgbSums;


} ClusterAccumulations;

typedef struct {

  Point **points;
  Point **chosenCentroids; // for first step
  Cluster **clusters;
  ClusterAccumulations clusterAccumulations;
  size_t startIdx;
  size_t endIdx;
  size_t numClusters;

} ThreadPartition;




double squaredEuclidianDistance(const Pixel *a, const Pixel *b) {
  double dr = a->r - b->r;
  double dg = a->g - b->g;
  double db = a->b - b->b;
  return dr * dr + dg * dg + db * db;
}


bool converged(const Pixel *oldCentroids, const Pixel *newCentroids,
               const int numClusters) {

  bool converged = true;

  for (size_t i = 0; i < numClusters; i++) {
    const Pixel oldCentroid = oldCentroids[i];
    const Pixel newCentroid = newCentroids[i];

    if (!(oldCentroid.r == newCentroid.r && oldCentroid.g == newCentroid.g &&
          oldCentroid.b == newCentroid.b)) {
      converged = false;
      break;
    }
  }

  return converged;
}

bool isChosen(const Point *point, Point **chosenCentroids,
              const size_t numClusters) {

  for (size_t i = 0; i < numClusters; i++) {
    const Point *centroid = chosenCentroids[i];
    if (point == centroid) {
      return 1;
    }
  }

  return 0;
}

void assignClusters(const size_t numClusters,
                    const ClusterAccumulations *clusterAccumulations,
                    Cluster **clusters, const Point *point) {
  double bestDistance = DBL_MAX;
  size_t bestIdx = 0;

  for (size_t j = 0; j < numClusters; j++) {
    const Cluster *cluster = clusters[j];
    const double distance = squaredEuclidianDistance(&point->pixel, &cluster->centroid);

    if (distance < bestDistance) {
      bestDistance = distance;
      bestIdx = j;
    }
  }

  AccumulatedRGBSums *selectedSums = clusterAccumulations->rgbSums;
  selectedSums[bestIdx].rgbAccum[0] += (double)clusters[bestIdx]->centroid.r;
  selectedSums[bestIdx].rgbAccum[1] += (double)clusters[bestIdx]->centroid.g;
  selectedSums[bestIdx].rgbAccum[2] += (double)clusters[bestIdx]->centroid.b;
  selectedSums[bestIdx].numPoints += 1;
}

void *threadedFirstStep(void *args) {

  const ThreadPartition *partition = args;

  const ClusterAccumulations clusterAccumulations =
      partition->clusterAccumulations;

  Point **points = partition->points;
  Point **chosenCentroids = partition->chosenCentroids;
  Cluster **clusters = partition->clusters;

  const size_t startIdx = partition->startIdx;
  const size_t endIdx = partition->endIdx;
  const size_t numClusters = partition->numClusters;

  for (size_t i = startIdx; i < endIdx; i++) {
    Point *currPoint = points[i];
    if (isChosen(currPoint, chosenCentroids, numClusters)) {
      continue;
    }
    assignClusters(numClusters, &clusterAccumulations, clusters, currPoint);
  }
  return nullptr;
}

void *threadedRefine( void *args) {

  const ThreadPartition *partition = args;

  const ClusterAccumulations clusterAccumulations =
      partition->clusterAccumulations;

  Point **points = partition->points;
  Cluster **clusters = partition->clusters;

  const size_t startIdx = partition->startIdx;
  const size_t endIdx = partition->endIdx;
  const size_t numClusters = partition->numClusters;

  for (size_t i = startIdx; i < endIdx; i++) {
    Point *currPoint = points[i];
    assignClusters(numClusters, &clusterAccumulations, clusters, currPoint);
  }
  return nullptr;
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

  srand(time(nullptr) * numClusters *
        numPoints); // this is a basic seed and not optimal for complex tasks
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

    cluster->_arrSize = STARTING_ARRAY_SIZE;
    clusterResetPoints(cluster);
    cluster->centroid.r = point->pixel.r;
    cluster->centroid.g = point->pixel.g;
    cluster->centroid.b = point->pixel.b;

    chosenCentroids[i] = point;
    clusters[i] = cluster;
  }

  return 0;
}

void firstStep(Cluster **clusters, Point **points, Point **chosenCentroids,
               int numClusters, int numPoints, int nThreads) {

  pthread_t threads[nThreads];
  ThreadPartition *partitions = malloc(sizeof(ThreadPartition) * nThreads);

  const size_t partitionSize = numPoints / nThreads;
  size_t startIdx = 0;
  size_t endIdx = partitionSize;

  for (size_t i = 0; i < nThreads; i++) {
    partitions[i].startIdx = startIdx;
    partitions[i].endIdx = (i == nThreads - 1) ? numPoints : endIdx;
    partitions[i].numClusters = numClusters;
    partitions[i].chosenCentroids = chosenCentroids;
    partitions[i].points = points;
    partitions[i].clusters = clusters;
    // partitions[i].clusterAccumulations = {};
    partitions[i].clusterAccumulations.rgbSums =
        malloc(sizeof(AccumulatedRGBSums) * numClusters);
    startIdx = endIdx;
    endIdx += partitionSize;

    if (pthread_create(&threads[i], nullptr, threadedFirstStep, &partitions[i]) !=
        0) {
      // error handle
    }
  }
  for (size_t i = 0; i < nThreads; i++) {
    pthread_join(threads[i], nullptr);
  }

  // update centroids
  for (size_t j = 0; j < numClusters; j++) {
    for (size_t i = 0; i < nThreads; i++) {

      AccumulatedRGBSums *accumulations =
          partitions[i].clusterAccumulations.rgbSums;

      size_t nPoints = accumulations[j].numPoints;
      double *rgbSums = accumulations[j].rgbAccum;
      clusters[j]->nPoints += nPoints;
      for (size_t chan = 0; chan < 3; chan++) {
        clusters[j]->rgbCumm[chan] += rgbSums[chan];
      }


    }
    calculateNewMean(clusters[j]);
    clusterResetPoints(clusters[j]);
  }

  for (size_t i = 0; i < nThreads; i++) {
    free(partitions[i].clusterAccumulations.rgbSums);
  }
  free(partitions);
}


void refineClusters(Cluster **clusters, Point **points, const int numClusters,
                    const int numPoints,  int nThreads) {

  Pixel *oldCentroids = malloc(sizeof(Pixel) * numClusters);
  Pixel *newCentroids = malloc(sizeof(Pixel) * numClusters);

  if (!oldCentroids) {
    return;
  }
  if (!newCentroids) {
    return;
  }

  pthread_t threads[nThreads];
  ThreadPartition *partitions = malloc(sizeof(ThreadPartition) * nThreads);

  const size_t partitionSize = numPoints / nThreads;
  size_t startIdx = 0;
  size_t endIdx = partitionSize;
  bool shouldIterate = 1;
  while (shouldIterate) {
    for (size_t i = 0; i < numClusters; i++) {
      oldCentroids[i] = points[i]->pixel;
    }

    for (size_t i = 0; i < nThreads; i++) {
      partitions[i].startIdx = startIdx;
      partitions[i].endIdx = (i == nThreads - 1) ? numPoints : endIdx;
      partitions[i].numClusters = numClusters;
      partitions[i].points = points;
      partitions[i].clusters = clusters;
      // partitions[i].clusterAccumulations = {};
      partitions[i].clusterAccumulations.rgbSums =
          malloc(sizeof(AccumulatedRGBSums) * numClusters);
      startIdx = endIdx;
      endIdx += partitionSize;

      if (pthread_create(&threads[i], nullptr, threadedRefine, &partitions[i]) !=
          0) {
        // error handle
      }
    }
    for (size_t i = 0; i < nThreads; i++) {
      pthread_join(threads[i], nullptr);
    }

    // update centroids
    for (size_t j = 0; j < numClusters; j++) {
      for (size_t i = 0; i < nThreads; i++) {

        const AccumulatedRGBSums *accumulations =
            partitions[i].clusterAccumulations.rgbSums;

        const size_t nPoints = accumulations[j].numPoints;
        const double *rgbSums = accumulations[j].rgbAccum;
        clusters[j]->nPoints += nPoints;
        for (size_t chan = 0; chan < 3; chan++) {
          clusters[j]->rgbCumm[chan] += rgbSums[chan];
        }


      }
      calculateNewMean(clusters[j]);
      clusterResetPoints(clusters[j]);
    }


    for (size_t i = 0; i < numClusters; i++) {
      newCentroids[i] = points[i]->pixel;
    }

    shouldIterate = !converged(oldCentroids, newCentroids, numClusters);
  }

  for (size_t i = 0; i < nThreads; i++) {
    free(partitions[i].clusterAccumulations.rgbSums);
  }
  free(oldCentroids);
  free(newCentroids);
  free(partitions);
}

void clusterResetPoints(Cluster *cluster) {
  cluster->rgbCumm[0] = 0;
  cluster->rgbCumm[1] = 0;
  cluster->rgbCumm[2] = 0;
  cluster->nPoints = 0;
}

int clusterAddPoint(Cluster *cluster, const Point *point) {

  cluster->rgbCumm[0] += point->pixel.r;
  cluster->rgbCumm[1] += point->pixel.g;
  cluster->rgbCumm[2] += point->pixel.b;
  cluster->nPoints++;
  return 1;
}

int calculateNewMean(Cluster *cluster) {

  if (cluster == nullptr)
    return -1;

  const double meanR = cluster->rgbCumm[0] / (double)cluster->nPoints;
  const double meanG = cluster->rgbCumm[1] / (double)cluster->nPoints;
  const double meanB = cluster->rgbCumm[2] / (double)cluster->nPoints;

  cluster->centroid.r = (unsigned char)round(meanR);
  cluster->centroid.g = (unsigned char)round(meanG);
  cluster->centroid.b = (unsigned char)round(meanB);

  return 0;
}

Pixel *generatePalette(Point **points, const int numPoints,
                       const int numClusters, const int numThreads) {

  Pixel *palette = malloc(sizeof(Pixel) * numClusters);
  Cluster **clusters = malloc(sizeof(Cluster *) * numClusters);
  Point **chosenCentroids = malloc(sizeof(Point *) * numClusters);
  if (!palette || !clusters || !chosenCentroids) {
    return nullptr;
  }

  initializeClusters(clusters, points, chosenCentroids, numClusters, numPoints);
  firstStep(clusters, points, chosenCentroids, numClusters, numPoints, numThreads);
  free(chosenCentroids);

  refineClusters(clusters, points, numClusters, numPoints,numThreads);

  for (size_t i = 0; i < numClusters; i++) {
    palette[i] = clusters[i]->centroid;
    free(clusters[i]);
  }
  free(clusters);

  return palette;
}