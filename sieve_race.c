// Compile with: clang -O3 sieve_race.c -o sieve_race.o
// Run with: ./sieve_race.o

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<math.h>

#define CACHE_SIZE 64*1024

#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#endif

#if defined(__GNUC__) || defined(__clang__)
#define OPTNONE __attribute__((optnone))
#else
#endif

NOINLINE
OPTNONE
size_t eratos_base(size_t b)
{
    // Allocate array of data-type: unsigned characters
    char *is_prime = malloc((b + 1) * sizeof(unsigned char));

    // Initialize all entries to 1
    for (size_t i = 0; i <= b; i++) is_prime[i] = 1;

    // Prime-sieving loop
    for (size_t p = 2; p * p <= b; p++) {
        if (is_prime[p]) {
            for (size_t m = p * p; m <= b; m += p) is_prime[m] = 0;
        }
    }

    // Count primes
    size_t count = 0;
    for(size_t i = 2; i <= b; i++) {
        if(is_prime[i]) ++count; 
    } 

    // Free the array
    free(is_prime);  

    // Return the number of primes found
    return count; 
}

NOINLINE
OPTNONE
size_t euler_base(size_t b)
{
    // Allocate a single-link active list for storing prime gaps
    unsigned short* l = malloc((b + 1) * sizeof(unsigned short)); 
    
    // Initialize all gaps to 1 (all indices are active)
    for(size_t i = 0; i <= b; i++) l[i] = 1; 

    // Set up an iterator for primes
    unsigned long long prime = 2; 
    while(prime * prime <= b) // Iterate through primes
    {
        // Set up an iterator for prime candidates 
        unsigned long long cand = prime; 
        while(cand * prime <= b) // Iterate through candidates
        {   
            // Set up an iterator for powers of a particular prime
            unsigned long long power = cand * prime; 
            while(power <= b) // Iterate through powers
            {
                // Find adjacent active indices
                unsigned int next = power + l[power]; 
                unsigned int prev = power - l[power - 1]; 

                // Update the differences between adjacent active indices 
                l[next - 1] = l[prev] = next - prev; 

                // Raise prime to the next power
                power *= prime; 
            }

            // Advance to the next candidate
            cand += l[cand]; 
        }

        // Advance to the next prime
        prime += l[prime]; 
    } 
 
    // Print primes if needed
    prime = 2;
    size_t count = 0;
    while(prime <= b)
    {
        // Advance to the next prime
        prime += l[prime]; 

        // Increment count of primes found
        ++count;
    }

    free(l); 

    // Return number of primes found
    return count;
}

static const unsigned char period[8] = {1, 7, 11, 13, 17, 19, 23, 29}; 

// Convert real number value tolist index
static inline size_t val_to_ind(unsigned long long val) {return 8 * val / 30;}

// Convert list index to real number value
static inline unsigned long long ind_to_val(size_t ind) {return 30 * (ind >> 3) + period[ind & 7]; }

// Run euler30
NOINLINE
OPTNONE
size_t euler30noopt(size_t b)
{
    size_t len = val_to_ind(b);
    unsigned char* l = malloc(len); 
    
    for(size_t i = 0; i < len; i++) l[i] = 1; 

    size_t i = 1; 
    unsigned long long prime = 7; 
    while(prime * prime <= b)
    {
        size_t j = i; 
        unsigned long long multiple = prime * prime; 
        while(multiple <= b)
        {   
            unsigned long long power = multiple; 
            while(power <= b)
            {
                // Compute index of candidate and prime power product
                unsigned int ind = val_to_ind(power); 
                
                // Find adjacent active indices
                unsigned int next = ind + l[ind]; 
                unsigned int prev = ind - l[ind - 1]; 

                // Update the differences between adjacent active indices 
                l[next - 1] = l[prev] = next - prev; 

                power *= prime; 
            }

            j += l[j]; 
            multiple = ind_to_val(j) * prime; 
        }

        i += l[i]; 
        prime = ind_to_val(i); 
    } 

    size_t count = 3; 
    i = 1; 
    while(i < len)
    {
        i += l[i]; 
        ++count; 
    }

    free(l); 

    return count; 
}

// Perform wheel-factorized Sieve of Eratosthenes
// Return prime count
NOINLINE
OPTNONE
size_t eratos30noopt(size_t bound){
    // Calculate max number of indexes needed
    size_t end = (bound + 29) / 30;

    // Allocate & initialize unsigned character array 
    unsigned char* cand = malloc(end);
    memset(cand, UINT8_MAX, end);

    // Iterate through primes
    size_t i = 1;
    unsigned long long prime = 7;
    while(prime*prime <= bound) // Stop when prime square exceeds bound
    {
        // Primality check (using a bitmask to access individual bits)
        if((cand[i >> 3] >> (i & 7)) & 1)
        {
            // Iterate through the 8 residue classes modulo 30
            for(size_t j = i; j < i + 8; j++)
            {
                // Compute word location of the first composite 
                size_t k = val_to_ind(prime * ind_to_val(j));

                // Compute bit mask for "fetch-and" operation --> marking composites 
                unsigned char bitmask = ~(1u << (k & 7));

                k >>= 3;

                // Mark composites as 0
                // Account for scalar loop remainder after unrolling
                for(; k < end; k += prime)
                {
                    cand[k] &= bitmask;
                }
            }
        }

        // Advance real index, find the next virtual prime value
        prime = ind_to_val(++i);
    }

    // Begin prime counting 
    size_t count = 3;

    // Iterate through character array
    i = 1;
    prime = 7;
    while(prime <= bound)
    {
        // Primality check
        if((cand[i >> 3] >> (i & 7)) & 1)
        {
            ++count;
        }

        // Advance real index, compute the next virtual prime value
        prime = ind_to_val(++i);
    }

    // Deallocate character array
    free(cand);

    // Return prime number count
    return count;
}

static inline size_t euler30opt(size_t b)
{
    size_t len = val_to_ind(b);
    unsigned char* l = malloc(len); 
    
    for(size_t i = 0; i < len; i++) l[i] = 1; 

    size_t i = 1; 
    unsigned long long prime = 7; 
    while(prime * prime <= b)
    {
        size_t j = i; 
        unsigned long long multiple = prime * prime; 
        while(multiple <= b)
        {   
            unsigned long long power = multiple; 
            while(power <= b)
            {
                // Compute index of candidate and prime power product
                unsigned int ind = val_to_ind(power); 
                
                // Find adjacent active indices
                unsigned int next = ind + l[ind]; 
                unsigned int prev = ind - l[ind - 1]; 

                // Update the differences between adjacent active indices 
                l[next - 1] = l[prev] = next - prev; 

                power *= prime; 
            }

            j += l[j]; 
            multiple = ind_to_val(j) * prime; 
        }

        i += l[i]; 
        prime = ind_to_val(i); 
    } 

    size_t count = 3; 
    i = 1; 
    while(i < len)
    {
        i += l[i]; 
        ++count; 
    }

    free(l); 

    return count; 
}

static inline size_t eratos30opt(size_t bound){
    // Calculate max number of indexes needed
    size_t end = (bound + 29) / 30;

    // Allocate & initialize unsigned character array 
    unsigned char* cand = malloc(end);
    memset(cand, UINT8_MAX, end);

    // Iterate through primes
    size_t i = 1;
    unsigned long long prime = 7;
    while(prime*prime <= bound) // Stop when prime square exceeds bound
    {
        // Primality check (using a bitmask to access individual bits)
        if((cand[i >> 3] >> (i & 7)) & 1)
        {
            // Iterate through the 8 residue classes modulo 30
            for(size_t j = i; j < i + 8; j++)
            {
                // Compute word location of the first composite 
                size_t k = val_to_ind(prime * ind_to_val(j));

                // Compute bit mask for "fetch-and" operation --> marking composites 
                unsigned char bitmask = ~(1u << (k & 7));

                k >>= 3;

                // Mark composites as 0
                // Account for scalar loop remainder after unrolling
                for(; k < end; k += prime)
                {
                    cand[k] &= bitmask;
                }
            }
        }

        // Advance real index, find the next virtual prime value
        prime = ind_to_val(++i);
    }

    // Begin prime counting 
    size_t count = 3;

    // Iterate through character array
    i = 1;
    prime = 7;
    while(prime <= bound)
    {
        // Primality check
        if((cand[i >> 3] >> (i & 7)) & 1)
        {
            ++count;
        }

        // Advance real index, compute the next virtual prime value
        prime = ind_to_val(++i);
    }

    // Deallocate character array
    free(cand);

    // Return prime number count
    return count;
}

// Compute floor(log_3(bound)) using iterative multiplication
static inline size_t flog3(unsigned long long bound)
{
    unsigned long long power = 1; 
    unsigned int count = 0; 
    while(power <= bound)
    {
        power *= 3; 
        count++; 
    }

    return count - 1; 
}

static inline size_t eratos_seg(size_t b)
{
    // Compute square root of the bound for generating the initial segment
    size_t sqrt_b = (size_t)(sqrt(b));

    // Calculate max number of indexes needed
    size_t end = (sqrt_b + 29) / 30;

    // Allocate & initialize unsigned character array 
    unsigned char* cand = malloc(end);
    memset(cand, UINT8_MAX, end);

    // Iterate through primes
    size_t i = 1;
    unsigned long long prime = 7;
    while(prime*prime <= sqrt_b) // Stop when prime square exceeds bound
    {
        // Primality check (using a bitmask to access individual bits)
        if((cand[i >> 3] >> (i & 7)) & 1)
        {
            // Iterate through the 8 residue classes modulo 30
            for(size_t j = i; j < i + 8; j++)
            {
                // Compute word location of the first composite 
                size_t k = val_to_ind(prime * ind_to_val(j));

                // Compute bit mask for "fetch-and" operation --> marking composites 
                unsigned char bitmask = ~(1u << (k & 7));

                k >>= 3;

                // Mark composites as 0
                // Account for scalar loop remainder after unrolling
                for(; k < end; k += prime)
                {
                    cand[k] &= bitmask;
                }
            }
        }

        // Advance real index, find the next virtual prime value
        prime = ind_to_val(++i);
    }

    //printf("Initial segment completed\n"); 

    size_t cap = sqrt_b /(flog3(sqrt_b) - 1); 

    // Allocate array for storing primes
    size_t* primes = malloc(cap * sizeof(size_t));

    // Allocate a static array for prime-associated segment starting values
    unsigned long long* start = malloc(cap * sizeof(unsigned long long));

    size_t seg_size = 2 * CACHE_SIZE;

    // Allocate segment 
    unsigned char* segment = malloc(seg_size); 

    // Segment bounds
    unsigned long long low = sqrt_b + 1;
    unsigned long long high = low + seg_size;

    // Store the first three primes skipped by the wheel
    primes[0] = 2;
    start[0] = ((sqrt_b + 2) >> 1) * 2;

    primes[1] = 3;
    start[1] = ((sqrt_b + 3)/ 3) * 3;

    primes[2] = 5;
    start[2] = ((sqrt_b + 5)/ 5) * 5;

    size_t len = 3;

    // Begin prime counting 
    size_t count = 3;

    // Iterate through character array
    i = 1;
    prime = 7;
    while(prime <= sqrt_b)
    {
        // Primality check
        if((cand[i >> 3] >> (i & 7)) & 1)
        {
            // Store prime
            primes[len] = prime;
            // Store starting multiple of this prime in the next segment
            if(prime * prime < low) start[len] = ((low + prime -1) / prime) * prime;
            else start[len] = prime * prime;
            ++len; 
            ++count;
        }

        // Advance real index, compute the next virtual prime value
        prime = ind_to_val(++i);
    }

    // printf("Primes:\n");
    // for(size_t i = 0; i < len; i++)
    // {
    //     printf("%llu ", primes[i]); 
    // }
    // printf("\n");

    // Deallocate character array
    free(cand);

    while(low <= b)
    {
        // Reset segment
        memset(segment, 1, seg_size);

        for(unsigned int i = 0; i < len; i++)
        {
            unsigned long long p = primes[i]; 
            if(p * p > high) break; 

            unsigned long long j; 
            for(j = start[i]; j < high; j += p)
            {
                segment[j - low] = 0; 
            }
            start[i] = j;
        }

        unsigned int end; 
        if(seg_size < b - low) end = seg_size;
        else end = b - low; 

        for(unsigned int k = 0; k < end; k++)
        {
            if(segment[k] == 1) 
            {
                //printf("%llu\n", low + k);
                ++count;
            }
        }

        low = high;
        high = low + seg_size;
    }

    free(segment); 

    free(primes); 

    free(start);

    // Return prime number count
    return count;
}

// Timing algorithm
// Function pointer type matching the signature of the algorithm to be timed
typedef size_t (*func)(size_t n);

// Times `num_runs` calls to `algo`, returns the average elapsed time in seconds.
static inline double time_algorithm(func algo, size_t n,
                       size_t num_runs, size_t warm_up) {
    struct timespec start, end;
    double total_elapsed = 0.0;

    volatile size_t num_of_primes; 

    for(size_t w = 0; w < warm_up; w++) num_of_primes += algo(n); 
    for (size_t run = 0; run < num_runs; run++) {
        clock_gettime(CLOCK_MONOTONIC, &start);

        // Prevent compiler from reordering calls and optimizing the timing mechanism
        num_of_primes += algo(n);
        __asm__ volatile("" ::: "memory"); 

        clock_gettime(CLOCK_MONOTONIC, &end);

        total_elapsed += (end.tv_sec - start.tv_sec)  + (end.tv_nsec - start.tv_nsec) / 1e9;
    }

    return total_elapsed / num_runs;
}

void run_trials(func algo)
{
    printf("Avg runtime at input size %.lf: %.9lf sec.\n", 1e3, time_algorithm(algo, 1000, 100, 50));

    printf("Avg runtime at input size %.lf: %.9lf sec.\n", 1e6, time_algorithm(algo, 1000000, 20, 10));

    printf("Avg runtime at input size %.lf: %.9lf sec.\n", 1e9, time_algorithm(algo, 1000000000, 4, 2));
}

int main()
{
    printf("\nEratosthenes baseline (unoptimized): \n");
    run_trials(eratos_base);
    printf("\nEuler baseline (unoptimized): \n");
    run_trials(euler_base);
    printf("\nEratosthenes wheel 30 (no compiler optimization): \n");
    run_trials(eratos30noopt);
    printf("\nEuler wheel 30 (no compiler optimization): \n");
    run_trials(euler30noopt);
    printf("\nEratosthenes wheel 30 (with compiler optimization): \n");
    run_trials(eratos30opt);
    printf("\nEuler wheel 30 (with compiler optimization): \n");
    run_trials(euler30opt);
    printf("\nEratosthenes segmented & optimized: \n");
    run_trials(eratos_seg);

    return 0; 
}