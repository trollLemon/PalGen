#include <stdio.h>
#include "printColor.h"

#include <stdlib.h>


char *asHex(const unsigned char r, const unsigned char g, const unsigned char b) {

    char *hex = malloc(8);
    if (hex == NULL) {
        return NULL;
    }

    snprintf(hex, sizeof hex, "#%02x%02x%02x", r, g, b);

    return hex;
}

