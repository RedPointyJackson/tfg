#include <stdio.h>
#include <stdlib.h>
#include "../lib/mersennetlib/SFMT.h"
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>


#define GREEN_C     "\033[92m"
#define RED_C       "\033[91m"
#define END_C       "\033[0m" // Return to normal color
#define ERRORHEADER "\033[1m\033[91mError:\033[0m " // Bold red
#define WARNHEADER  "\033[1m\033[93mWarning:\033[0m " // Bold yellow
#define INFOHEADER  "\033[1m\033[94mInfo:\033[0m " // Bold blue
#define DIE(msg) do{                                                    \
                    fprintf(stderr,ERRORHEADER"%s",msg);                \
                    fprintf(stderr,"\tIn function %s\n",__func__);      \
                    fprintf(stderr,"\tIn file %s:%d\n",                 \
                              __FILE__,__LINE__);                       \
                    exit(EXIT_FAILURE);                                 \
                    } while(false);

#define IS_BIG_ENDIAN (!*(unsigned char *)&(uint16_t){1})

// 0 001 ⋯ 001 001 001
const uint64_t ONES_BY_TRIPLETS = 0x1249249249249249;

struct net {
  // Each Int64 variable contains the variables of 21 samples, one for
  // each 3 bytes.
  int64_t   L;
  uint64_t* spins;
  double    beta;
  double    probs[7];
  uint64_t* J_right;
  uint64_t* J_up;
  uint64_t* J_front;
  int64_t   mc_steps;
};

sfmt_t sfmt; // Random number generator object for the mersenne twister.
#define RANDOM_F (sfmt_genrand_res53(&sfmt)) // r ∈ [0,1) (53 bit resolution)
#define RANDOM_SPIN (sfmt_genrand_res53(&sfmt)>0.5 ? 1 : -1) // r ∈ {±1}
#define RANDOM_BIT (sfmt_genrand_res53(&sfmt)>0.5 ? 1 : 0) // r ∈ {0,1}

uint64_t random_bits(void){
  // For each triple of the uint64, assign a randomly 000, 001.
  uint64_t base = 0;
  for(int64_t i=0;i<64/3;i++){
    base = base << 3;
    base += RANDOM_BIT;
  }
  return base;
}

void seed_mersenne(uint32_t seed){
  sfmt_init_gen_rand(&sfmt,seed);
}

uint32_t hash_guides(uint64_t* J_up, uint64_t* J_right, uint64_t* J_front, int64_t L){
  uint32_t hash;
  int64_t i;
  for(hash = i = 0; i < L; ++i){
    hash += J_right[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
    hash += J_up[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
    hash += J_front[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
}

// above_spin = net[i + guide_up[i]], for example.
int64_t* guide_up     = NULL;
int64_t* guide_front  = NULL;
int64_t* guide_right  = NULL;
int64_t* guide_down   = NULL;
int64_t* guide_behind = NULL;
int64_t* guide_left   = NULL;

void init_guides_pbc(int64_t L){ // Periodic boundary conditions
  if (guide_up     != NULL ) free(guide_up);
  if (guide_front  != NULL ) free(guide_front);
  if (guide_right  != NULL ) free(guide_right);
  if (guide_down   != NULL ) free(guide_down);
  if (guide_behind != NULL ) free(guide_behind);
  if (guide_left   != NULL ) free(guide_left);

  guide_up     = malloc(sizeof(*guide_up    ) * L);
  guide_front  = malloc(sizeof(*guide_front ) * L);
  guide_right  = malloc(sizeof(*guide_right ) * L);
  guide_down   = malloc(sizeof(*guide_down  ) * L);
  guide_behind = malloc(sizeof(*guide_behind) * L);
  guide_left   = malloc(sizeof(*guide_left  ) * L);

  for(int64_t i=0;i<L;i++){
    guide_right[i]  = +1;
    guide_up[i]     = +L;
    guide_front[i]  = +L*L;
    guide_left[i]   = -1;
    guide_down[i]   = -L;
    guide_behind[i] = -L*L;
  }

  guide_right[L-1] = -(L-1);
  guide_up[L-1]    = -L*(L-1);
  guide_front[L-1] = -L*L*(L-1);

  guide_left[0]    = +(L-1);
  guide_down[0]    = +L*(L-1);
  guide_behind[0]  = +L*L*(L-1);
}

enum direction {RIGHT, LEFT, UP, DOWN, FRONT, BEHIND};

uint64_t get_shifted_spin(uint64_t* spins, size_t L, size_t idx, size_t N, enum direction dir){
  switch (dir)
    {
    case RIGHT:
      for(size_t i=0;i<N;i++){
        size_t x = idx % L;
        idx = idx+guide_right[x];
      } break;

    case LEFT:
      for(size_t i=0;i<N;i++){
        size_t x = idx % L;
        idx = idx+guide_left[x];
      } break;

    case UP:
      for(size_t i=0;i<N;i++){
        size_t y = (idx/L) % L;
        idx = idx+guide_up[y];
      } break;

    case DOWN:
      for(size_t i=0;i<N;i++){
        size_t y = (idx/L) % L;
        idx = idx+guide_down[y];
      } break;

    case FRONT:
      for(size_t i=0;i<N;i++){
        size_t z = (idx/L/L) % L;
        idx = idx+guide_front[z];
      } break;

    case BEHIND:
      for(size_t i=0;i<N;i++){
        size_t z = (idx/L/L) % L;
        idx = idx+guide_behind[z];
      } break;
    }
  return spins[idx];
}

void set_beta(struct net* SG, double beta){
  if (isnan(beta))
    DIE("Intended to use β = NaN.\n");
  if (beta < 0)
    DIE("Intended to use β < 0.\n");

  SG->beta = beta;
  SG->probs[0] = exp( +12*beta );
  SG->probs[1] = exp(  +8*beta );
  SG->probs[2] = exp(  +4*beta );
  SG->probs[3] = exp(   0*beta );
  SG->probs[4] = exp(  -4*beta );
  SG->probs[5] = exp(  -8*beta );
  SG->probs[6] = exp( -12*beta );
}

// Obtain a spin glass in a random configuration.
struct net get_spinglass(int64_t L, double beta){
  init_guides_pbc(L);
  struct net newnet;
  newnet.L = L;
  newnet.mc_steps = 0;
  set_beta(&newnet, beta);

  // Initialize all randomly.
  newnet.spins = malloc(sizeof(*newnet.spins) * L*L*L);
  for(int64_t i=0;i<L*L*L;i++){
    newnet.spins[i] = random_bits();
  }
  newnet.J_right = malloc( sizeof(*newnet.spins) * L*L*L);
  newnet.J_up    = malloc( sizeof(*newnet.spins) * L*L*L);
  newnet.J_front = malloc( sizeof(*newnet.spins) * L*L*L);
  for(int64_t i=0;i<L*L*L;i++){
    newnet.J_right[i] = random_bits();
    newnet.J_up[i] = random_bits();
    newnet.J_front[i] = random_bits();
  }

  // Check randomness
  int64_t identical_els = 0;
  for(int64_t i=0;i<L*L*L;i++){
    if (newnet.J_right[i] == newnet.J_right[0]) identical_els++;
    if (newnet.J_up[i]    == newnet.J_right[0]) identical_els++;
    if (newnet.J_front[i] == newnet.J_right[0]) identical_els++;
  }
  if (identical_els == 3*L*L*L){
    DIE("All bonds in J_guide are identical!");
  }
  return newnet;
}

void free_spinglass(struct net* spinglass){
  free(spinglass->spins);
  spinglass->spins = NULL;

  free(spinglass->J_up);
  spinglass->J_up = NULL;

  free(spinglass->J_right);
  spinglass->J_right = NULL;

  free(spinglass->J_front);
  spinglass->J_front = NULL;
}

 // ∑ⱼ σᵢ⊕Jᵢⱼ⊕σⱼ ∈ {0,1,⋯,6}
uint64_t local_energy(struct net* SG, int64_t idx, int64_t x, int64_t y, int64_t z){
  uint64_t* S     = SG->spins;
  uint64_t* Jr    = SG->J_right;
  uint64_t* Ju    = SG->J_up;
  uint64_t* Jf    = SG->J_front;

  uint64_t right  = S[idx] ^ Jr[ idx                   ] ^ S[ idx + guide_right[x]  ];
  uint64_t left   = S[idx] ^ Jr[ idx + guide_left[x]   ] ^ S[ idx + guide_left[x]   ];
  uint64_t up     = S[idx] ^ Ju[ idx                   ] ^ S[ idx + guide_up[y]     ];
  uint64_t down   = S[idx] ^ Ju[ idx + guide_down[y]   ] ^ S[ idx + guide_down[y]   ];
  uint64_t front  = S[idx] ^ Jf[ idx                   ] ^ S[ idx + guide_front[z]  ];
  uint64_t behind = S[idx] ^ Jf[ idx + guide_behind[z] ] ^ S[ idx + guide_behind[z] ];

  return right + left + up + down + front + behind;
}

void metropolis(struct net* SG){
  int64_t curr_idx = 0;
  int64_t L = SG->L;

  for(int64_t z=0;z<L;z++){
    for(int64_t y=0;y<L;y++){
      for(int64_t x=0;x<L;x++){

        uint64_t E = local_energy(SG, curr_idx, x, y, z);
        double randomnum = RANDOM_F;

        for(int_fast8_t spinidx=0; spinidx<21; spinidx++){

          uint64_t Pidx = E;
          Pidx = Pidx >> 3*spinidx;
          Pidx = Pidx &  0x7; // Select the first 3 bits.

          double P = SG->probs[Pidx];
          if (Pidx <= 3 || P > randomnum){
            uint64_t selectmask = 0x1;
            selectmask = selectmask << 3*spinidx;
            SG->spins[curr_idx] ^= selectmask; // Flip current σᵢ.
          }
        }
        curr_idx++;
      }
    }
  }
  SG->mc_steps++;
}

void linspace(int64_t* array, int64_t from, int64_t to, int64_t N){
  if(N < 1){DIE("linspace requires N>1.\n");}
  for(int64_t i=0; i<N;i++){
    double B = to;
    double A = from;
    array[i] = (B-A)/(N-1)*i + A;
  }
  array[0] = from;
  array[N-1] = to;
}

void logspace(int64_t* array, int64_t from, int64_t to, int64_t N){
  if(N < 1){DIE("Logspace requires N>1.\n");}
  for(int64_t i=0; i<N;i++){
    double B = log10(to);
    double A = log10(from);
    array[i] = pow(10, (B-A)/(N-1)*i + A);
  }
  array[0] = from;
  array[N-1] = to;
}

int int_compare(const void *a_ptr, const void *b_ptr){
  const int *a = a_ptr, *b = b_ptr;
  return (*a < *b) ? -1 : (*a > *b);
}

int64_t sort_unique(int64_t* array, size_t N){

  qsort(array, N, sizeof(array[0]), int_compare);
  size_t left = 0;
  size_t right = N-1;

  while (left != right){
    if (array[left] == array[left+1]){
      for(size_t i=left;i<right;i++){
        array[i] = array[i+1];
      }
      right--;
    } else {
      left++;
    }
  }

  int64_t newsize = right+1;
  return newsize;
}

double mean(double* array, int64_t N){
  double acc = 0;
  for(int64_t i=0;i<N;i++){
    acc += array[i];
  }
  return acc/N;
}

double variance(double* array, int64_t N){
  double mu = mean(array,N);
  double accdiff = 0;
  for(int64_t i=0;i<N;i++){
    accdiff += (array[i] - mu)*(array[i] - mu);
  }
  return accdiff/(N-1);
}

void bootstrap(double* array, int64_t N
               , double (*function)(double*,int64_t)
               , double* fmean, double* fstdev
               ){
  double fraction = 0.36787944117144233; // 1/e
  int64_t samples = 1000; // samples to take
  int64_t chunksize = floor(fraction*N);
  if(chunksize < 1){
    DIE("Bootstrap chunks have less than one element!\n");
  }
  double values[samples];
  for(int64_t i=0;i<samples;i++){
    double randchunk[chunksize];
    for(int64_t j=0;j<chunksize;j++){
      int64_t idx = floor(RANDOM_F*(chunksize+1));
      randchunk[j] = array[idx];
    }
    values[i] = function(randchunk,chunksize);
  }
  *fmean  = mean(values,samples);
  *fstdev = sqrt(variance(values,samples));
}

double mean_energy(struct net* SG){
  double E = 0;
  int64_t i = 0;
  int64_t L = SG->L;
  for(int64_t z=0;z<L;z++){
    for(int64_t y=0;y<L;y++){
      for(int64_t x=0;x<L;x++){
        uint64_t localEidx = local_energy(SG,i,x,y,z);
        for(int_fast8_t si=0;si<21;si++){
          // select the first 3 bits of the current spin.
          double thislocalE = (double)((localEidx >> 3*si) & 0x7);
          double Hterm = 6 - 2*thislocalE;
          E += Hterm;
        }
        i++;
      }
    }
  }
  return E/21/(L*L*L)/3/2;
}

double mean_magnetization(struct net* SG){
  double M = 0;
  int64_t L = SG->L;
  for(int64_t i=0;i<L*L*L;i++){
    for(int_fast8_t si=0;si<21;si++){
      // select the first 3 bits of the current spin.
      uint64_t thismagnet = (SG->spins[i] >> 3*si) & 0x1;
      M += thismagnet? +1 : -1;
    }
  }
  return M/21/L/L/L;
}

// Save mean and stdev of SD in the locations pointed by mn,std.
void schwingerdyson(struct net* SG, double* mn, double* std){
  int64_t L = SG->L;
  double* values = malloc(L*L*L*21 * sizeof (*values));
  int64_t ind = 0;
  for(int64_t z=0;z<L;z++){
    for(int64_t y=0;y<L;y++){
      for(int64_t x=0;x<L;x++){
        uint64_t localEidx = local_energy(SG,ind,x,y,z);
        for(int_fast8_t i=0;i<21;i++){
          // select the first 3 bits of the current spin.
          double thislocalE = (double)((localEidx >> 3*i) & 0x7);
          double Hterm = 6 - 2*thislocalE;
          values[ind*21+i] = exp(2.0*SG->beta*Hterm);
        }
        ++ind;
      }
    }
  }
  bootstrap(values,L*L*L*21,mean,mn,std);
  free(values);
}

struct anneal_contents {
  int64_t L;
  int64_t tws;
  int64_t* tw_list;
  int64_t ts;
  int64_t* t_list;
  int64_t Nmeas;
  double temp;
  uint64_t* Jup;
  uint64_t* Jright;
  uint64_t* Jfront;
  int64_t* mc_list;
  uint64_t* spins_list;
};

void free_anneal_contents(struct anneal_contents* AC){
  free(AC->tw_list);
  free(AC->t_list);
  free(AC->Jup);
  free(AC->Jright);
  free(AC->Jfront);
  free(AC->mc_list);
  free(AC->spins_list);
}

struct anneal_contents get_anneal_contents(FILE* f){
  struct anneal_contents AC;

  char header[9];
  fread(header, sizeof(char), 8, f);
  header[8] = '\0';
  if (strcmp(header,"ANNEALER")){
    DIE("File is not an annealer file.\n");
  }

  int64_t created_as_big_endian = -1;
  fread(&created_as_big_endian, sizeof(created_as_big_endian), 1, f);

  if (created_as_big_endian != IS_BIG_ENDIAN) DIE("Endianness of file differs from PC's.\n");

  fread(&(AC.L), sizeof(AC.L), 1, f);

  int64_t V = AC.L*AC.L*AC.L;

  fread(&(AC.tws), sizeof(AC.tws), 1, f);
  int64_t* tw_list = malloc(AC.tws* sizeof(*tw_list));
  fread(tw_list, sizeof(tw_list[0]), AC.tws, f);
  AC.tw_list = tw_list;


  fread(&(AC.ts), sizeof(AC.ts), 1, f);
  int64_t* t_list = malloc(AC.ts* sizeof(*t_list));
  fread(t_list, sizeof(t_list[0]), AC.ts, f);
  AC.t_list = t_list;

  fread(&(AC.Nmeas), sizeof(AC.Nmeas), 1, f);
  fread(&(AC.temp), sizeof(AC.temp), 1, f);

  uint64_t* Jup = malloc(V * sizeof(*Jup));
  fread(Jup, sizeof(Jup[0]), V, f);
  AC.Jup = Jup;

  uint64_t* Jright = malloc(V * sizeof(*Jright));
  fread(Jright, sizeof(Jright[0]), V, f);
  AC.Jright = Jright;

  uint64_t* Jfront = malloc(V * sizeof(*Jfront));
  fread(Jfront, sizeof(Jfront[0]), V, f);
  AC.Jfront = Jfront;

  int64_t* mc_list     = malloc(AC.Nmeas * sizeof(*mc_list));
  uint64_t* spins_list = malloc(AC.Nmeas*V * sizeof(*spins_list));

  for(int64_t i=0; i<AC.Nmeas; i++){
    if (feof(f)) DIE("Anneal file was shorter than expected.\n");
    fread(mc_list    + i  , sizeof(*mc_list)   , 1, f);
    fread(spins_list + i*V, sizeof(*spins_list), V, f);
    if (mc_list[i] < 1) DIE("Readed mc < 1.\n");
  }

  AC.mc_list = mc_list;
  AC.spins_list = spins_list;

  if (feof(f)) DIE("File ended before footer was readed.\n");

  char footer[4];
  fread(footer, sizeof(char), 3, f);
  footer[3] = '\0';
  if (strcmp(footer,"END")) DIE("File footer wasn't \"END\".");

  return AC;
}

// Get L, spinlist and J's from annealer's output.
// Use the lattice with mc ≥ step.
void continue_net_from_stdin(int64_t step, struct net* SG){

  struct anneal_contents AC = get_anneal_contents(stdin);
  SG->L = AC.L;
  int64_t V = AC.L * AC.L * AC.L;

  int64_t recoverpoint = 0;
  for(int64_t i=0;i<AC.Nmeas;i++){
    if (AC.mc_list[i] >= step){
      recoverpoint = AC.mc_list[i];
      for(int64_t j=0;j<V;j++){
        SG->spins[j]   = AC.spins_list[i*V + j];
        SG->J_up[j]    = AC.Jup[j];
        SG->J_right[j] = AC.Jright[j];
        SG->J_front[j] = AC.Jfront[j];
      }
      break;
    }
  }
  if (recoverpoint == 0) DIE("Step wasn't found in continuation net.\n");
}
