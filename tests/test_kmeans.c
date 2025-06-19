
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "kmeans.h"
#include "common.h"
#include "defaults.h"

//TODO: move these elsewhere
constexpr int nPoints = 2;
constexpr int clusterId = 0;

void test_refine() {

    Point pointA  = {
        {
            10,
            10,
            10,
        }
    };
    Point pointB =  {
        {
            20,
            20,
            20,
        }
    };
    Point pointC = {
        {
            10,
            10,
            10,
        }
    };
    Point pointD = {
        {
            10,
            10,
            10,
        }
    };

    Point *points[4] = {&pointA, &pointB, &pointC, &pointD};

    Cluster **clusters = malloc(sizeof(Cluster*) * 2);
    Point **chosen = malloc(sizeof(Point*) * 2);


    initializeClusters(clusters, points, chosen,2, 4);

    firstStep(clusters,points,chosen,2,4);
    refineClusters(clusters,points,2,4);
    free(chosen);

    for (size_t i = 0; i < 2; ++i) {
        free(clusters[i]->points);
        free(clusters[i]);
    }

    free(clusters);


}

void test_firstStep() {

    Point pointA  = {
        {
            10,
            10,
            10,
        }
    };
    Point pointB =  {
        {
            20,
            20,
            20,
        }
    };
    Point pointC = {
        {
            10,
            10,
            10,
        }
    };
    Point pointD = {
        {
            10,
            10,
            10,
        }
    };

    Point *points[4] = {&pointA, &pointB, &pointC, &pointD};

    Cluster **clusters = malloc(sizeof(Cluster*) * 2);
    Point **chosen = malloc(sizeof(Point*) * 2);


    initializeClusters(clusters, points, chosen,2, 4);

    firstStep(clusters,points,chosen,2,4);

    for (size_t i = 0; i < 2; i++) {

        const Pixel currCentroid = clusters[i]->centroid;
        const Pixel oldCentroid = chosen[i]->pixel;

        // assert(currCentroid.r != oldCentroid.r);
        // assert(currCentroid.g != oldCentroid.g);
        // assert(currCentroid.b != oldCentroid.b);

    }

    free(chosen);

    for (size_t i = 0; i < 2; ++i) {
        free(clusters[i]->points);
        free(clusters[i]);
    }

    free(clusters);


}

void test_initializeCentroids() {

    Point pointA  = {
        {
            10,
            10,
            10,
        }
    };
    Point pointB =  {
                {
                    20,
                    20,
                    20,
                }
    };
    Point pointC = {
                    {
                        10,
                        10,
                        10,
                    }
    };
    Point pointD = {
            {
                10,
                10,
                10,
            }
    };

    Point *points[4] = {&pointA, &pointB, &pointC, &pointD};

    Cluster **clusters = malloc(sizeof(Cluster*) * 2);
    Point **chosen = malloc(sizeof(Point*) * 2);


    int result = initializeClusters(clusters, points, chosen,2, 4);


    assert(result == 0);

    //check if centroids and chosen points are collected correctly
    for (size_t i = 0; i < 2; ++i) {
        assert(clusters[i]->centroid.r == chosen[i]->pixel.r);
        assert(clusters[i]->centroid.g == chosen[i]->pixel.g);
        assert(clusters[i]->centroid.b == chosen[i]->pixel.b);

    }

    //check if there are any duplicates

    assert(chosen[0] != chosen[1]);

    free(chosen);

    for (size_t i = 0; i < 2; ++i) {
        free(clusters[i]->points);
        free(clusters[i]);
    }

    free(clusters);

}


void test_calculateNewMean() {


    Pixel pixelA = {10,10,10 };
    Pixel pixelB = {20,20,20 };
    Pixel centroid = {30,50,60};

    Pixel expectedCentroid = {15,15,15};

    Point pointA = { pixelA };
    Point pointB = { pixelB };


    Point* points[nPoints] = {&pointA, &pointB};

    Cluster cluster = { centroid, points, clusterId,nPoints };

    const int result = calculateNewMean(&cluster);

    assert(result == 0);
    // assert(expectedCentroid.r == cluster.centroid.r);
    // assert(expectedCentroid.g == cluster.centroid.g);
    // assert(expectedCentroid.b == cluster.centroid.b);


}


void test_calculateNewMean_NULL_Centroid() {

    const int result = calculateNewMean(NULL);
    assert(result == -1);


}


void test_ClusterAddPoint() {

    Cluster *cluster = malloc(sizeof(Cluster));
    Point** points = malloc(sizeof(Point*) * STARTING_ARRAY_SIZE);
    cluster->points = points;
    cluster->nPoints = 0;
    cluster->_arrSize = STARTING_ARRAY_SIZE;

    Point* dummyPoint = malloc(sizeof(Point));

    // simulate adding a bunch of points
    // should hit the array limit and create a new array twice the size, copying everything over

    for (int i = 0; i < STARTING_ARRAY_SIZE * 2; ++i) {
        clusterAddPoint(cluster, dummyPoint);
    }

    for (size_t i = 0; i < STARTING_ARRAY_SIZE * 2; ++i) {
        assert(cluster->points[i] == dummyPoint);
    }

    assert(cluster->nPoints == (STARTING_ARRAY_SIZE * 2) );
    assert(cluster->_arrSize == STARTING_ARRAY_SIZE * 2);

    free(dummyPoint);
    free(cluster->points);
    free(cluster);

}

void test_GeneratePalette() {
    Point pointA  = {
        {
            10,
            10,
            10,
        }
    };
    Point pointB =  {
        {
            20,
            20,
            20,
        }
    };
    Point pointC = {
        {
            10,
            10,
            10,
        }
    };
    Point pointD = {
        {
            10,
            10,
            10,
        }
    };

    int numPoints = 4;
    int numClusters = 2;
    Point *points[4] = {&pointA, &pointB, &pointC, &pointD};

    // Point** points = malloc(sizeof(Point*) * numPoints);
    // points[0] = &pointA;
    // points[1] = &pointB;
    // points[2] = &pointC;
    // points[3] = &pointD;

    Pixel* centroids = generatePalette(points, numPoints, numClusters);

    free(centroids);
}



int main() {
    test_refine();
    test_firstStep();
    test_initializeCentroids();
    test_calculateNewMean();
    test_calculateNewMean_NULL_Centroid();
    test_ClusterAddPoint();
    test_GeneratePalette();
    return 0;
}