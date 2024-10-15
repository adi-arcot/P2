#include <string.h>

#include "zmorton.hpp"
#include "binhash.hpp"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Spatial hashing implementation}
 * 
 * In the current implementation, we assume [[HASH_DIM]] is $2^b$,
 * so that computing a bitwise of an integer with [[HASH_DIM]] extracts
 * the $b$ lowest-order bits.  We could make [[HASH_DIM]] be something
 * other than a power of two, but we would then need to compute an integer
 * modulus or something of that sort.
 * 
 *@c*/

#define HASH_MASK (HASH_DIM-1)

unsigned particle_bucket(particle_t* p, float h)
{
    unsigned ix = p->x[0]/h;
    unsigned iy = p->x[1]/h;
    unsigned iz = p->x[2]/h;
    return zm_encode(ix & HASH_MASK, iy & HASH_MASK, iz & HASH_MASK);
}

unsigned particle_neighborhood(unsigned* buckets, particle_t* p, float h)
{
    /* BEGIN TASK */
    unsigned ix = p->x[0] / h;
    unsigned iy = p->x[1] / h;
    unsigned iz = p->x[2] / h;

    int count = 0;

    // Loop through all neighbors (including the current bin and up to 26 neighboring bins)
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                // Compute the neighbor's bin coordinates
                unsigned neighbor_ix = (ix + dx) & HASH_MASK;
                unsigned neighbor_iy = (iy + dy) & HASH_MASK;
                unsigned neighbor_iz = (iz + dz) & HASH_MASK;

                // Encode the neighbor's bin using the Z-Morton encoding
                buckets[count++] = zm_encode(neighbor_ix, neighbor_iy, neighbor_iz);
            }
        }
    }
    
    return count;  // Return the number of neighboring bins
    /* END TASK */
}

void hash_particles(sim_state_t* s, float h)
{
    /* BEGIN TASK */
    // Reset the hash table
    memset(s->hash, 0, STATE_HASH_SIZE * sizeof(particle_t*));

    // Hash each particle
    for (int i = 0; i < s->n; ++i) {
        particle_t* p = &s->part[i];
        
        // Calculate the bucket for this particle
        unsigned bucket = particle_bucket(p, h);

        // Add the particle to the front of the linked list for its bucket
        p->next = s->hash[bucket];
        s->hash[bucket] = p;
    }
    /* END TASK */
}
