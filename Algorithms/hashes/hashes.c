
#include "hashes.h"

int rgbEqual(const Pixel *p1, const Pixel *p2) {

  return p1->r == p2->r && p1->g == p2->g && p1->b == p2->b;
}

int rgbHash(const Pixel p) {
  const unsigned char r = p.r;
  const unsigned char g = p.g;
  const unsigned char b = p.b;

  return ((r & 0xff) << 16) + ((g & 0xff) << 8) + (b & 0xff);
}