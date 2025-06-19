

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
    Point** points;
    int     clusterId;
    int     nPoints;
    size_t _arrSize;
} Cluster;

int initializeClusters(Cluster** clusters, Point** points, Point** chosenCentroids, int numClusters, int numPoints);
int calculateNewMean(Cluster* cluster);
int clusterAddPoint(Cluster* cluster, Point* point);
void clusterResetPoints(Cluster* cluster);
void firstStep(Cluster** clusters, Point** points, Point** chosenCentroids, int numClusters, int numPoints);
void refineClusters(Cluster** clusters, Point** points, int numClusters, int numPoints);

Pixel* generatePalette(Point** points, int numPoints, int numClusters);

#endif //KMEANS_H
