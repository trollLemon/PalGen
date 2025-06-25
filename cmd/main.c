#include <stdio.h>
#include <unistd.h>
#define STB_IMAGE_IMPLEMENTATION
#include "common.h"
#include "defaults.h"
#include "kmeans.h"
#include "printColor.h"
#include "stb_image.h"
#include "hashes.h"
#include "khash.h"

KHASH_MAP_INIT_INT(PixelHash, int)

typedef struct {
  char *imagePath;
  int numColors;
  int numJobs;

} Args;

int parseArgs(int argc, char *argv[], Args *args) {

  int opt;
  while ((opt = getopt(argc, argv, "i:c:j:")) != -1) {
    switch (opt) {
    case 'i':
      args->imagePath = optarg;
      break;
    case 'c':
      args->numColors = atoi(optarg);
      break;
    case 'j':
      args->numJobs = atoi(optarg);
      break;
    default:
      return 1;
    }
  }

  return 0;
}

int main(int argc, char **argv) {

  Args args;

  args.imagePath = NULL;
  args.numColors = DEFAULT_PALETTE_SIZE;
  args.numJobs = 1;

  parseArgs(argc, argv, &args);

  if (args.imagePath == NULL) {
    printf("Error: please pass an image path with -i\n");
    return 1;
  }

  if (args.numColors < 1) {
    printf("Error: please pass a number of colors with -c greater than 0 -c\n");
    return 1;
  }
  if (args.numJobs < 1) {
    printf("Error: please pass a number of jobs with -j greater than 0\n");
    return 1;
  }

  int width, height, channels;

  const char *filename = argv[1];

  unsigned char *image =
      stbi_load(args.imagePath, &width, &height, &channels, 0);

  if (!image) {
    fprintf(stderr, "Failed to load image: %s\n", filename);
    return 1;
  }

  const int numPixels = width * height;
  khash_t(PixelHash) *hashes = kh_init(PixelHash);
  int ret;

  Point **points = malloc(numPixels * sizeof(Point *));

  int pixelIdx = 0; // TODO: there is a better way to do this probably

  for (size_t x = 0; x < width; x++) {
    for (size_t y = 0; y < height; y++) {

      const size_t index = channels * (y * width + x);
      const unsigned char r = image[index];
      const unsigned char g = image[index + 1];
      const unsigned char b = image[index + 2];

      const Pixel pixel = {r, g, b};
      const int hash = rgbHash(pixel);
      khiter_t k = kh_get(PixelHash, hashes, hash);

      // skip this pixel since we already have its color
      if (k != kh_end(hashes)) {
        continue;
      }

      k = kh_put(PixelHash, hashes, hash, &ret);
      kh_value(hashes, k) = 1;

      points[pixelIdx] = malloc(sizeof(Point));
      points[pixelIdx]->pixel = pixel;
      pixelIdx++;
    }
  }

  Pixel *palette =
      generatePalette(points, pixelIdx, args.numColors, args.numJobs);

  for (size_t i = 0; i < pixelIdx; i++) {
    free(points[i]);
  }

  free(image);
  free(points);

  for (size_t i = 0; i < args.numColors; i++) {
    const Pixel p = palette[i];
    char *color = asHex(p.r, p.g, p.b);
    printf("%s\n", color);
    free(color);
  }
  free(palette);
  kh_destroy(PixelHash, hashes);
  return 0;
}