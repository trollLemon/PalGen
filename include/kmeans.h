

#ifndef KMEANS_H
#define KMEANS_H

#include "common.h"
#include <stddef.h>
typedef struct {

    Pixel pixel;
    int clusterId;
} Point;

typedef struct {
    Pixel   centroid;       // centroid value
    double rgbCumm[3];
    int     clusterId;
    size_t     nPoints;
    size_t _arrSize;
} Cluster;

int initializeClusters(Cluster** clusters, Point** points, Point** chosenCentroids, int numClusters, int numPoints);
int calculateNewMean(Cluster* cluster);
int clusterAddPoint(Cluster* cluster, const Point* point);
void clusterResetPoints(Cluster* cluster);
void firstStep(Cluster** clusters, Point** points, Point** chosenCentroids, int numClusters, int numPoints, int nThreads);
void refineClusters(Cluster** clusters, Point** points, int numClusters, int numPoints, int nThreads);

Pixel* generatePalette(Point** points, int numPoints, int numClusters, int numThreads);

#endif //KMEANS_H
