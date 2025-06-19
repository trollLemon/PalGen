
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "printColor.h"



void test_asHex() {

    unsigned char r = 10;
    unsigned char g = 10;
    unsigned char b = 10;

   char expected[8] = "#0a0a0a";

    char* hex = asHex(r,g,b);


    assert(hex !=NULL);

    assert(strncasecmp(hex, expected, 8) == 0);

    free(hex);

}




int main() {

    test_asHex();
}