#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>

#include <cstdint>
#include <iomanip>
#include <iostream>

#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;

#include <boost/multiprecision/cpp_bin_float.hpp>
using boost::multiprecision::cpp_bin_float_50;

#define INT512(u) ((int512_t)u)

typedef struct {
	int512_t vector[2];
} Vector;

typedef struct {
	int512_t matrix[2][2];
    cpp_bin_float_50 double_matrix[2][2];
} Matrix;

Vector vector_new(int512_t a, int512_t b) {
	return (Vector){.vector={a, b}};
}

int512_t vector_get_element(Vector *vector, int index) {
	return vector->vector[index];
}

void vector_set_element(Vector *vector, int512_t value, int index) {
	vector->vector[index] = value;
}

int512_t vector_dot(Vector *vector, Vector *other) {
	int512_t result = 0;
	for (int i = 0; i < 2; i++) {
		result += vector->vector[i] * other->vector[i];
	}
	return result;
}

Vector matrix_get_col(Matrix *matrix, int i);
Vector vector_multiply(Vector *vector, Matrix *m) {
	Vector ret = vector_new(0, 0);

	Vector v1 = matrix_get_col(m, 0);
	Vector v2 = matrix_get_col(m, 1);
	
	vector_set_element(&ret, vector_dot(vector, &v1), 0);
	vector_set_element(&ret, vector_dot(vector, &v2), 1);
	return ret;
}

Vector vector_scale(Vector *vector, int512_t scaler) {
	Vector ret = vector_new(0, 0);

	for (int i = 0; i < 2; i++) {
		vector_set_element(&ret, vector->vector[i] * scaler, i);
	}

	return ret;
}

Vector vector_add(Vector *vector, Vector *other) {
	Vector ret = vector_new(0, 0);

	for (int i = 0; i < 2; i++) {
		vector_set_element(&ret, vector->vector[i] + other->vector[i], i);
	}

	return ret;
}

bool vector_le(Vector *vector, Vector *other) {
	bool a = vector->vector[0] <= other->vector[0];
	bool b = vector->vector[1] <= other->vector[1];
	return a && b;
}

int512_t vector_norm_sq(Vector *vector) {
	return pow(vector->vector[0], 2) + pow(vector->vector[1], 2);
}

Matrix matrix_new() {
	return (Matrix){0};
}

int512_t matrix_get_element(Matrix *matrix, int row, int col) {
	return matrix->matrix[row][col];
}

void matrix_set_element(Matrix *matrix, int512_t value, int row, int col) {
	matrix->matrix[row][col] = value;	
}

Vector matrix_get_col(Matrix *matrix, int i) {
	return vector_new(matrix->matrix[0][i], matrix->matrix[1][i]);
}

Vector matrix_get_row(Matrix *matrix, int i) {
	return vector_new(matrix->matrix[i][0], matrix->matrix[i][1]);
}

void matrix_set_row(Matrix *matrix, int row, Vector *other) {
	for (int col = 0; col < 2; col++) {
		matrix->matrix[row][col] = other->vector[col];
	}
}

void matrix_set_col(Matrix *matrix, int col, Vector *other) {
	for (int row = 0; row < 2; row++) {
		matrix->matrix[row][col] = other->vector[row];
	}
}

Matrix matrix_scale(Matrix *matrix, cpp_bin_float_50 scaler, bool use_double) {
	Matrix ret = matrix_new();

	for (int row = 0; row < 2; row++) {
		for (int col = 0; col < 2; col++) {
            if (!use_double) {
			    matrix_set_element(&ret, (int512_t)((cpp_bin_float_50)matrix->matrix[row][col] * scaler), row, col);
            }
            else {
                // std::cout << (cpp_bin_float_50)matrix->matrix[row][col] * scaler << "\n";
                ret.double_matrix[row][col] = (cpp_bin_float_50)matrix->matrix[row][col] * scaler;
            }
		}
	}

	return ret;
}

int512_t matrix_det(Matrix *matrix) {
	return matrix->matrix[0][0] * matrix->matrix[1][1] - matrix->matrix[1][0] * matrix->matrix[0][1];
}

Matrix matrix_inverse(Matrix *matrix) {
	Matrix inv = matrix_new();

	int512_t det = matrix_det(matrix);
	assert(det != 0);

	matrix_set_element(&inv, matrix->matrix[1][1], 0, 0);
	matrix_set_element(&inv, matrix->matrix[0][0], 1, 1);
	matrix_set_element(&inv, -matrix->matrix[1][0], 1, 0);
	matrix_set_element(&inv, -matrix->matrix[0][1], 0, 1);

	inv = matrix_scale(&inv, ((cpp_bin_float_50)1 / (cpp_bin_float_50)det), true);
	return inv; 
}

Matrix matrix_swap_rows(Matrix *matrix) {
	Matrix result = matrix_new();

	Vector a = matrix_get_row(matrix, 1);
	Vector b = matrix_get_row(matrix, 0);
	
	matrix_set_row(&result, 0, &a);
	matrix_set_row(&result, 1, &b);

	return result;
}


int512_t gcdExtended(int512_t a, int512_t b, int512_t *x, int512_t *y) {
    if (a == 0) {
        *x = 0;
        *y = 1;
        return b;
    }

    int512_t x1, y1;
    int512_t gcd = gcdExtended(b % a, a, &x1, &y1);

    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}

// Function to find the modular inverse
int512_t modInverse(int512_t a, int512_t m) {
    int512_t x, y;
    int512_t gcd = gcdExtended(a, m, &x, &y);

    // Modular inverse exists only if gcd(a, m) = 1
    if (gcd != 1) {
        printf("Modular inverse does not exist\n");
        return 0; // No modular inverse
    }

    // Make sure x is positive
    int512_t modInverse = (x % (int512_t)m + m) % m;
    return modInverse;
}

void vector_print(Vector *vector) {
	// printf("%Lf %Lf\n", vector->vector[0], vector->vector[1]);
	printf("vector: {\n");
	std::cout << vector->vector[0] << std::endl;
	std::cout << vector->vector[1] << std::endl;
	printf("}\n");
}

void matrix_print(Matrix *matrix) {
	printf("matrix: {\n");
	std::cout << matrix->matrix[0][0] << std::endl;
	std::cout << matrix->matrix[0][1] << std::endl;
	std::cout << matrix->matrix[1][0] << std::endl;
	std::cout << matrix->matrix[1][1] << std::endl;
	printf("}\n");
}

void matrix_double_print(Matrix *matrix) {
	printf("matrix: {\n");
	std::cout << matrix->double_matrix[0][0] << std::endl;
	std::cout << matrix->double_matrix[0][1] << std::endl;
	std::cout << matrix->double_matrix[1][0] << std::endl;
	std::cout << matrix->double_matrix[1][1] << std::endl;
	printf("}\n");
}

bool lagrange_gauss_condition(Matrix *r) {
	Vector v_1 = matrix_get_row(r, 0);
	int512_t norm_sq_1 = vector_norm_sq(&v_1);

	Vector v_2 = matrix_get_row(r, 1);
	int512_t norm_sq_2 = vector_norm_sq(&v_2);

	return norm_sq_1 > norm_sq_2;
}

Matrix lagrange_gauss(Matrix *basis) {
	Matrix r = matrix_new();
	r.matrix[0][0] = basis->matrix[0][0];
	r.matrix[1][0] = basis->matrix[1][0];		
	r.matrix[0][1] = basis->matrix[0][1];		
	r.matrix[1][1] = basis->matrix[1][1];		

	do {
		if (lagrange_gauss_condition(&r)) {
			r = matrix_swap_rows(&r);
		}

		Vector r_0 = matrix_get_row(&r, 0);
		Vector r_1 = matrix_get_row(&r, 1);
		
		int512_t mu = vector_dot(&r_0, &r_1) / vector_norm_sq(&r_0);
		Vector temp = matrix_get_row(&r, 0);
		temp = vector_scale(&temp, -mu);

		Vector t = vector_add(&r_1, &temp);
		matrix_set_row(&r, 1, &t);
	} while(lagrange_gauss_condition(&r));
	
	return r;
}

typedef struct {
    int x, z;
} Pos;

Pos getFeatureChunkInRegion(uint64_t seed, int regX, int regZ) {
    /*
    // Vanilla like implementation.
    setSeed(&seed, regX*341873128712 + regZ*132897987541 + seed + config.salt);

    Pos pos;
    pos.x = nextInt(&seed, 24);
    pos.z = nextInt(&seed, 24);
    */
    Pos pos;
    const uint64_t K = 0x5deece66dULL;
    const uint64_t M = (1ULL << 48) - 1;
    const uint64_t b = 0xb;

    // set seed
    seed = seed + regX*341873128712ULL + regZ*132897987541ULL + 34222645;
    seed = (seed ^ K);
    seed = (seed * K + b) & M;

    uint64_t r = 25;
    if (r & (r-1))
    {
        pos.x = (int)(seed >> 17) % r;
        seed = (seed * K + b) & M;
        pos.z = (int)(seed >> 17) % r;
    }
    else
    {
        // Java RNG treats powers of 2 as a special case.
        pos.x = (int)((r * (seed >> 17)) >> 31);
        seed = (seed * K + b) & M;
        pos.z = (int)((r * (seed >> 17)) >> 31);
    }

    return pos;
}

Pos getFeaturePos(uint64_t seed, int regX, int regZ) {
    Pos pos = getFeatureChunkInRegion(seed, regX, regZ);

    pos.x = (int)(((uint64_t)regX*40 + pos.x) << 4);
    pos.z = (int)(((uint64_t)regZ*40 + pos.z) << 4);
    return pos;
}


bool pos_has_ruined_portal(uint64_t world_seed, int x, int z) {
    assert((x % 16 == 0) && (z % 16 == 0));
    int reg_x = floor((x>>4) / 40.0);
    int reg_z = floor((z>>4) / 40.0);

    Pos p = getFeaturePos(world_seed, reg_x, reg_z);
    return (p.x == x) && (p.z == z);
}

void find_solutions_in_box(int512_t a, int512_t b, int512_t target, int512_t mod, Vector min, Vector max, uint64_t world_seed) {
	target = target % mod;
	int512_t big_mod = mod;
	int512_t binv = modInverse(b, big_mod);
	// std::cout << binv << std::endl;

	int512_t new_z_center = (binv * target) % big_mod;
	// std::cout << new_z_center << std::endl;

	Vector initial_solution = vector_new(0, new_z_center);

	Vector temp = vector_new(0, -new_z_center);
	min = vector_add(&min, &temp);
	max = vector_add(&max, &temp);

	// vector_print(&min);
	// vector_print(&max);

	Matrix basis = matrix_new();
	Vector va = vector_new(0, big_mod);
	Vector vb = vector_new(1, binv * -a);

	matrix_set_row(&basis, 0, &va);
	matrix_set_row(&basis, 1, &vb);
	// matrix_print(&basis);

	Matrix reduced_basis = lagrange_gauss(&basis);
	// matrix_print(&reduced_basis);
	Matrix inv = matrix_inverse(&reduced_basis);
	// matrix_double_print(&inv);

    cpp_bin_float_50 transformed_mins[2] = {0.0, 0.0};
    cpp_bin_float_50 transformed_maxes[2] = {0.0, 0.0};

    for (int row = 0; row < 2; row++) { 
        for (int col = 0; col < 2; col++) {
            cpp_bin_float_50 element = inv.double_matrix[row][col];
            // std::cout << element << "\n";

            if (element >= 0) {
                transformed_maxes[col] +=  (element * (cpp_bin_float_50)vector_get_element(&max, row));
                transformed_mins[col] += (element * (cpp_bin_float_50)vector_get_element(&min, row));
            }
            else {
                transformed_maxes[col] += (element * (cpp_bin_float_50)vector_get_element(&min, row));
                transformed_mins[col] += (element * (cpp_bin_float_50)vector_get_element(&max, row));
            }
        }
    }

    for (int512_t x = (int512_t)transformed_mins[0] - 2; x < (int512_t)transformed_maxes[0] + 2; x++) {
        for (int512_t z = (int512_t)transformed_mins[1] - 2; z < (int512_t)transformed_maxes[1] + 2; z++) {
            Vector coords = vector_new(x, z);
            coords = vector_multiply(&coords, &reduced_basis);
            if (vector_le(&min, &coords) && vector_le(&coords, &max)) {
                coords = vector_add(&coords, &initial_solution);
                if ((int)coords.vector[0] % 16 == 0 && (int)coords.vector[1] % 16 == 0) {
                    printf("Testing: %d %d for a ruined portal\n", (int)coords.vector[0], (int)coords.vector[1]);
                    if (pos_has_ruined_portal(world_seed, (int)coords.vector[0], (int)coords.vector[1])) {
                        printf("------------------------\n");
                        printf("FOUND! %" PRIu64 "\n", world_seed);
                        printf("x=%d z=%d\n", (int)coords.vector[0], (int)coords.vector[1]);
                        printf("------------------------\n");
                    }
                }
            }
        }
    }

}

#define XRSR_MIX1          0xbf58476d1ce4e5b9
#define XRSR_MIX2          0x94d049bb133111eb
#define XRSR_MIX1_INVERSE  0x96de1b173f119089
#define XRSR_MIX2_INVERSE  0x319642b2d24d8ec3
#define XRSR_SILVER_RATIO  0x6a09e667f3bcc909
#define XRSR_GOLDEN_RATIO  0x9e3779b97f4a7c15

uint64_t mix64(uint64_t a) {
	a = (a ^ a >> 30) * XRSR_MIX1;
	a = (a ^ a >> 27) * XRSR_MIX2;
	return a ^ a >> 31;
}

uint64_t rotl64(uint64_t x, uint8_t b)
{
    return (x << b) | (x >> (64-b));
}

typedef struct {
    uint64_t lo, hi;
} Xoroshiro;

static void xSetSeed(Xoroshiro *xr, uint64_t value)
{
    const uint64_t XL = 0x9e3779b97f4a7c15ULL;
    const uint64_t XH = 0x6a09e667f3bcc909ULL;
    const uint64_t A = 0xbf58476d1ce4e5b9ULL;
    const uint64_t B = 0x94d049bb133111ebULL;
    uint64_t l = value ^ XH;
    uint64_t h = l + XL;
    l = (l ^ (l >> 30)) * A;
    h = (h ^ (h >> 30)) * A;
    l = (l ^ (l >> 27)) * B;
    h = (h ^ (h >> 27)) * B;
    l = l ^ (l >> 31);
    h = h ^ (h >> 31);
    xr->lo = l;
    xr->hi = h;
}

static void xSetFeatureSeed(Xoroshiro *xr, uint64_t p_190065_, int p_190066_, int p_190067_) {
    uint64_t i = p_190065_ + (long)p_190066_ + (long)(10000 * p_190067_);
    xSetSeed(xr, i);
}

static uint64_t xNextLong(Xoroshiro *xr)
{
    uint64_t l = xr->lo;
    uint64_t h = xr->hi;
    uint64_t n = rotl64(l + h, 17) + l;
    h ^= l;
    xr->lo = rotl64(l, 49) ^ h ^ (h << 21);
    xr->hi = rotl64(h, 28);
    return n;
}

static uint64_t xSetDecorationSeed(Xoroshiro *xr, uint64_t p_64691_, int p_64692_, int p_64693_) {
    // this.setSeed(p_64691_);
    xSetSeed(xr, p_64691_);
    uint64_t i = xNextLong(xr) | 1L;
    uint64_t j = xNextLong(xr) | 1L;
    uint64_t k = (uint64_t)p_64692_ * i + (uint64_t)p_64693_ * j ^ p_64691_;
    // this.setSeed(k);
    xSetSeed(xr, k);
    return k;
}

static int xNextInt(Xoroshiro *xr, uint32_t n)
{
    uint64_t r = (xNextLong(xr) & 0xFFFFFFFF) * n;
    if ((uint32_t)r < n)
    {
        while ((uint32_t)r < (~n + 1) % n)
        {
            r = (xNextLong(xr) & 0xFFFFFFFF) * n;
        }
    }
    return r >> 32;
}

static double xNextDouble(Xoroshiro *xr)
{
    return (xNextLong(xr) >> (64-53)) * 1.1102230246251565E-16;
}

static float xNextFloat(Xoroshiro *xr)
{
    return (xNextLong(xr) >> (64-24)) * 5.9604645E-8F;
}

static void xSkipN(Xoroshiro *xr, int count)
{
    while (count --> 0)
        xNextLong(xr);
}

static uint64_t xNextLongJ(Xoroshiro *xr)
{
    int32_t a = xNextLong(xr) >> 32;
    int32_t b = xNextLong(xr) >> 32;
    return ((uint64_t)a << 32) + b;
}

typedef struct {
    Xoroshiro internal;
} RNG; // Bruh I really didn't want to have to do this.

RNG rng_new() {
    return (RNG){.internal=(Xoroshiro){0}};
}

static void rng_set_seed(RNG *rng, uint64_t seed) {
    seed ^= XRSR_SILVER_RATIO;
    rng->internal.lo = mix64(seed);
    rng->internal.hi = mix64(seed + XRSR_GOLDEN_RATIO);
}

static void rng_set_internal(RNG *rng, uint64_t lo, uint64_t hi) {
    rng->internal.lo = lo;
    rng->internal.hi = hi;
}

static uint32_t rng_next(RNG *rng, int32_t bits) {
    return xNextLong(&rng->internal) >> (64 - bits);
}

static int32_t rng_next_int(RNG *rng, uint32_t bound) {
    uint32_t r = rng_next(rng, 31);
    uint32_t m = bound - 1;
    if ((bound & m) == 0) {
        // (int)((long)p_188504_ * (long)this.next(31) >> 31);
        r = (uint32_t)((uint64_t)bound * (uint64_t)r >> 31);
    }
    else {
        for (uint32_t u = r; (int32_t)(u - (r = u % bound) + m) < 0; u = rng_next(rng, 31));
    }
    return r;
}

static uint64_t rng_next_long(RNG *rng) {
    int32_t i = rng_next(rng, 32);
    int32_t j = rng_next(rng, 32);
    uint64_t k = (uint64_t)i << 32;
    return k + (uint64_t)j;
}

static uint64_t rng_set_feature_seed(RNG *rng, uint64_t p_190065_, int32_t p_190066_, int32_t p_190067_) {
    uint64_t i = p_190065_ + (uint64_t)p_190066_ + (uint64_t)(10000 * p_190067_);
    //printf("Salt = %" PRIu64 "\n", (uint64_t)p_190066_ + (uint64_t)(10000 * p_190067_));
    rng_set_seed(rng, i);
    return i;
}

uint64_t reverse_decoration_seed(uint64_t decorator_seed, int index, int step) {
    return decorator_seed - (uint64_t)index - 10000L * (uint64_t)step;
}

static uint64_t rng_set_decoration_seed(RNG *rng, uint64_t world_seed, int32_t x, int32_t z) {
    rng_set_seed(rng, world_seed);

    uint64_t a = rng_next_long(rng) | 1L;
    uint64_t b = rng_next_long(rng) | 1L;

    // printf("the k to recover = %" PRIu64 "\n", (a * (uint64_t)x + b * (uint64_t)z));
    uint64_t k = (a * (uint64_t)x + b * (uint64_t)z) ^ world_seed;
    // printf("real k = %" PRIu64 "\n", k);
    // printf("invert k = %" PRIu64 "\n", k ^ world_seed);
    rng_set_seed(rng, k);
    return k;
}

uint64_t gcdExtendedIter(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
    // Initialize previous coefficients (for a)
    int64_t x1 = 1, x2 = 0;
    // Initialize previous coefficients (for b)
    int64_t y1 = 0, y2 = 1;
    
    while (b > 0) {
        uint64_t q = a / b;  // quotient
        
        // Update GCD calculation
        uint64_t temp_b = b;
        b = a % b;
        a = temp_b;
        
        // Update coefficients
        int64_t temp_x = x1;
        x1 = x2;
        x2 = temp_x - q * x2;
        
        int64_t temp_y = y1;
        y1 = y2;
        y2 = temp_y - q * y2;
    }
    
    // Set final coefficients
    *x = x1;
    *y = y1;
    
    return a;  // a is now the GCD
}

// Function to find the modular inverse
uint64_t modInverse(uint64_t a, uint64_t m) {
    int64_t x, y;
    uint64_t gcd = gcdExtendedIter(a, m, &x, &y);

    if (gcd != 1) {
        printf("Modular inverse does not exist\n");
        return 0; // No modular inverse
    }

    uint64_t modInverse = (x % (int64_t)m + m);
    return modInverse;
}

static void find_coordinates(RNG *rng, uint64_t world_seed, uint64_t k) {
    rng_set_seed(rng, world_seed);

    uint64_t a = rng_next_long(rng) | 1L;
    uint64_t b = rng_next_long(rng) | 1L;
    k = k ^ world_seed;

    // printf("%" PRIu64 " %" PRIu64 " %" PRIu64 "\n", a, b, k);

    Vector min = vector_new(-30000000, -30000000);
	Vector max = vector_new(30000000, 30000000);
	find_solutions_in_box(INT512(a), INT512(b), INT512(k), INT512(1<<64), min, max, world_seed);

	// for (int64_t x = -4096; x < 4096; x++) {
	// 	uint64_t r = (k - a * x); //implicit mod 2^64

	// 	uint64_t b_inv = modInverse(b, 1ull << 63);
	// 	uint64_t z = (r * b_inv);

	// 	int64_t real_z = (int64_t)z;

	// 	if (real_z > -4096 && real_z < 4096) {
	// 		printf("Found a match! %" PRIu64 ": %" PRId64 " %" PRId64"\n", world_seed, x, real_z);	
	// 	}
	// }
}

void usage(void) {
    printf("USAGE: ./find <num_threads> <start_seed> <end_seed>\n");
}

int main(int argc, char **argv) {
    RNG rng = rng_new();

    uint64_t i = rng_set_decoration_seed(&rng, 11342447246ull, 28963264, -382544); //We have the population seed, we now have to bruteforce a location + world seed that gives us this population seed...
    rng_set_feature_seed(&rng, i, 10, 4); //decoration step. THIS IS THE DECORATOR SEED, its a function of the population seed above. Reversing it gives the population seed above.

    uint64_t loot_seed = rng_next_long(&rng);
    printf("%" PRIu64 "\n", loot_seed);
    return 0;

    if (argc <= 1) {
        usage();
        return 0;
    }
    int num_threads = atoi(argv[1]);
    uint64_t start_seed = std::stoull(argv[2]);
    uint64_t end_seed = std::stoull(argv[3]);

    uint64_t seeds_per_thread = (uint64_t)((float)(end_seed - start_seed) / (float)num_threads);


    // omp_set_thread_num();
    // printf("%s\n", pos_has_ruined_portal(1234ull, -288, -304) ? "true": "false");
    // return 0;
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        // printf("starting thread...\n");
        uint64_t k = reverse_decoration_seed(18789082084264ull, 10, 4);
        RNG rng = rng_new();

        uint64_t start = start_seed + omp_get_thread_num() * seeds_per_thread;
        uint64_t end = start_seed + (omp_get_thread_num() * seeds_per_thread) + seeds_per_thread;

        printf("Starting thread to search from: %" PRIu64 " to %" PRIu64 "\n", start, end);

        for (uint64_t world_seed = start; world_seed < end; world_seed++) {
            find_coordinates(&rng, world_seed, k);
        }
    }
    return 0;

    // // //test scenario:
    // uint64_t world_seed = 1231873789;

    // RNG rng = rng_new();

    // uint64_t i = rng_set_decoration_seed(&rng, world_seed, 22400, -19222); 
    // uint64_t pop = rng_set_feature_seed(&rng, i, 10, 4);

    // find_coordinates(&rng, world_seed, i);

    // printf("%" PRIu64 "\n", pop);
    // return 1;

    return 0;
}