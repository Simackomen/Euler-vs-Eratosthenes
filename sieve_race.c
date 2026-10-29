// Compile with: clang -O3 sieve_race.c -o sieve_race.o
// Run with: ./sieve_race.o

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

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

// Timing algorithm
// Function pointer type matching the signature of the algorithm to be timed
typedef size_t (*func)(unsigned int n);

// Times `num_runs` calls to `algo`, returns the average elapsed time in seconds.
static inline double time_algorithm(func algo, unsigned int n,
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


int main()
{
    printf("Number of primes found: %lu.\n", euler_base(100)); 
    return 0; 
}