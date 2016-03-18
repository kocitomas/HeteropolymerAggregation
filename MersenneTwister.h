#include <iostream>
#include <cmath>

using namespace std;


#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL                 /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL                         /* Least significant 31 bits */


static unsigned long long mt[NN];                /* The array for the state vector */
static int mti=NN+1;                             /* mti==NN+1 means mt[NN] is not initialized */
void init_genrand64(unsigned long long seed);


double genrand64_real1(void);                    /* generates a double on the interval (0,1) */
unsigned long long genrand64_int64(void);
void init_by_array64(unsigned long long init_key[],unsigned long long key_length);
