

#ifndef KMEANS_H
#define KMEANS_H

#include "common.h"
#include <stddef.h>
typedef struct {

  Pixel pixel;
} Point;

typedef struct {
  Pixel centroid; // centroid value
  double rgbCumm[3];
  size_t nPoints;
} Cluster;

int initializeClusters(Cluster **clusters, Point **points, int numClusters,
                       int numPoints);
int calculateNewMean(Cluster *cluster);
int clusterAddPoint(Cluster *cluster, const Point *point);
void clusterResetPoints(Cluster *cluster);
void refineClusters(Cluster **clusters, Point **points, int numClusters,
                    int numPoints, int nThreads);

Pixel *generatePalette(Point **points, int numPoints, int numClusters,
                       int numThreads);

#endif // KMEANS_H
