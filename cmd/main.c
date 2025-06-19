#include <stdio.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "defaults.h"
#include "kmeans.h"
#include "common.h"
#include "printColor.h"








int main(int argc, char **argv) {

    int paletteSize = DEFAULT_PALETTE_SIZE;

    int width, height, channels;

    const char *filename = argv[1];

    unsigned char *image = stbi_load(filename, &width, &height, &channels, 0);

    if (!image) {
        fprintf(stderr, "Failed to load image: %s\n", filename);
        return 1;
    }

    const int numPixels = width * height;

    Point** points = malloc(numPixels * sizeof(Point*));

    int  pixelIdx = 0; //TODO: there is a better way to do this probably

    for (size_t x = 0; x < width; x++) {
        for (size_t y = 0; y < height; y++) {
            size_t index = channels * (y * width + x);
            unsigned char r = image[index];
            unsigned char g = image[index + 1];
            unsigned char b = image[index + 2];
            points[pixelIdx] = malloc(sizeof(Point));
            points[pixelIdx]->clusterId = -1;
            points[pixelIdx]->pixel.r = r;
            points[pixelIdx]->pixel.g = g;
            points[pixelIdx]->pixel.b = b;
            pixelIdx++;
        }
    }

   Pixel * palette = generatePalette(points, numPixels, paletteSize);

    for (size_t i = 0; i < paletteSize; i++) {
        const Pixel p = palette[i];
        char *color = asHex(p.r, p.g, p.b);
        printf("%s\n", color );
        free(color);
    }

    for (size_t i = 0; i < numPixels; i++) {
        free(points[i]);
    }

    free(image);
    free(points);
    free(palette);
    return 0;
}