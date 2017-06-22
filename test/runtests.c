#define EPSILON 1e-7 // FP precision comparison

#include "../lib/glassy.c"
const int64_t theL = 8; // For the tests

#include <assert.h>
#include <time.h>

const char* BLUE       = "\033[94m";
const char* GREEN      = "\033[92m";
const char* YELLOW     = "\033[93m";
const char* RED        = "\033[91m";
const char* BOLD       = "\033[1m";
const char* UNDERLINE  = "\033[4m";
const char* ENDC       = "\033[0m";

////////////////////////////////////////////////
//    Define a μframework for unit testing    //
////////////////////////////////////////////////
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

#define mu_assert(message, test) do {           \
    if (!(test)){                               \
      return message                            \
        "\n\tin assertion "                     \
        #test                                   \
        "\n\tat "                               \
        AT;                                     \
    }                                           \
  } while (0)

#define mu_run_test(test) do {                  \
    char* message = test();                     \
    tests_run++;                                \
    if (message)                                \
      return message;                           \
  } while (0)

int tests_run = 0;

// Copy a net
struct net copy_net(struct net original){
  struct net copycat = get_spinglass(theL,original.beta);
  int64_t L = original.L;
  for(int64_t i=0;i<L*L*L;i++){
    copycat.J_up[i] = original.J_up[i];
    copycat.J_right[i] = original.J_right[i];
    copycat.J_front[i] = original.J_front[i];
    copycat.spins[i] = original.spins[i];
  }
  return copycat;
}

///////////////////////////////////////
//         TESTING FUNCTIONS         //
///////////////////////////////////////
char* utilities_test(){
  int64_t dummy[10] = {0};
  linspace(dummy,4,10,10);
  mu_assert("Linspace start fails", dummy[0] == 4);
  mu_assert("Linspace end fails", dummy[9] == 10);
  logspace(dummy,4,10,10);
  mu_assert("Logspace start fails", dummy[0] == 4);
  mu_assert("Logspace end fails", dummy[9] == 10);
  return 0;
}

char* random_numbers_test(){ // Sometimes needs user inspection
  printf("Non-deterministic tests:\n");
  // Check correct distribution and support
  double accf = 0;
  double acci = 0;
  int64_t iters = (int64_t)1e6;
  for(int64_t i=0;i<iters;i++){
    double randf = RANDOM_F;
    int8_t randi = RANDOM_SPIN;
    mu_assert("RANDOM_F ∉ (0,1)"
              , randf < 1 && randf > 0);
    mu_assert("RANDOM_SPIN ∉ {±1}"
              , abs(randi) == 1);
    accf += randf/iters;
    acci += (double)randi/iters;
  }
  printf("    %s⟨RANDOM_F⟩%s    = %.4lf       ≃ 0.5\n", BOLD, ENDC, fabs(accf));
  printf("    %s⟨RANDOM_SPIN⟩%s = %.4lf       ≃ 0\n"  , BOLD, ENDC, fabs(acci));
  return 0;
}

char* pbc_test(){
  init_guides_pbc(theL);
  // Check that they don't have any zero.
  for(int64_t i=0;i<theL;i++){
    mu_assert("+x guide has zeros!", guide_right[i]  != 0);
    mu_assert("-x guide has zeros!", guide_left[i]   != 0);
    mu_assert("+y guide has zeros!", guide_up[i]     != 0);
    mu_assert("-y guide has zeros!", guide_down[i]   != 0);
    mu_assert("+z guide has zeros!", guide_front[i]  != 0);
    mu_assert("-z guide has zeros!", guide_behind[i] != 0);
  }
  // If the guides are correct, if I iterate over they for the width
  // of the net I should get to the starting point (periodic boundary
  // conditions).
  int64_t spin = 0; // Start at (0,0,0)
  // Move in each direction and check at the end of every one of them.
  for(int64_t i=0;i<theL;i++){
    spin += guide_right[i];
  } mu_assert("+x guide is not periodic!", spin == 0);
  for(int64_t i=0;i<theL;i++){
    spin += guide_left[i];
  } mu_assert("-x guide is not periodic!", spin == 0);
  for(int64_t i=0;i<theL;i++){
    spin += guide_up[i];
  } mu_assert("+y guide is not periodic!", spin == 0);
  for(int64_t i=0;i<theL;i++){
    spin += guide_down[i];
  } mu_assert("-y guide is not periodic!", spin == 0);
  for(int64_t i=0;i<theL;i++){
    spin += guide_front[i];
  } mu_assert("+z guide is not periodic!", spin == 0);
  for(int64_t i=0;i<theL;i++){
    spin += guide_behind[i];
  } mu_assert("-z is not periodic!", spin == 0);
  // Move in diagonal over x,y,z
  for(int64_t i=0;i<theL;i++){
    // In each iteration, we are in x,y,z = (i,i,i)
    spin += guide_right[i];
    spin += guide_up[i];
    spin += guide_front[i];
  } mu_assert("Walking in diagonal is not periodic!", spin == 0);
  // Also, if I move the same amount in each direction and its
  // reverse, I should get to the starting point. That is, for two
  // opposite guides `a` and `b`, `a[i] = -b[i+1]`:
  for(int64_t i=0;i<theL;i++){
    mu_assert("x guides are not correctly paired!"
              , guide_right[i] == -guide_left[(i+1)%theL]   );
    mu_assert("y guides are not correctly paired!"
              , guide_up[i]    == -guide_down[(i+1)%theL]   );
    mu_assert("z guides are not correctly paired!"
              , guide_front[i] == -guide_behind[(i+1)%theL] );
  }

  return 0;
}


char* shift_spin_test(){
  init_guides_pbc(theL);
  struct net SG = get_spinglass(theL, 0);

  // Test that moving L returns us where we were. Is like pbc test but
  // more high-level.
  for(int64_t i=0;i<theL*theL*theL;i++){
    mu_assert("Moving function does not work (right).",
              SG.spins[i] == get_shifted_spin(SG.spins, theL, i, theL, RIGHT) );
    mu_assert("Moving function does not work (left).",
              SG.spins[i] == get_shifted_spin(SG.spins, theL, i, theL, LEFT) );
    mu_assert("Moving function does not work (up).",
              SG.spins[i] == get_shifted_spin(SG.spins, theL, i, theL, UP) );
    mu_assert("Moving function does not work (down).",
              SG.spins[i] == get_shifted_spin(SG.spins, theL, i, theL, DOWN) );
    mu_assert("Moving function does not work (front).",
              SG.spins[i] == get_shifted_spin(SG.spins, theL, i, theL, FRONT) );
    mu_assert("Moving function does not work (behind).",
              SG.spins[i] == get_shifted_spin(SG.spins, theL, i, theL, BEHIND) );
  }

  free_spinglass(&SG);
  return 0;
}

char* initialization_test(){
  double beta = 1/(RANDOM_F);
  struct net mynet = get_spinglass(theL,beta);
  mu_assert("New net does not have the intended β", mynet.beta == beta);
  mu_assert("New net probs[0] ≠ exp(+12β)", mynet.probs[0] == exp( +12*beta ));
  mu_assert("New net probs[1] ≠ exp(+8β)" , mynet.probs[1] == exp(  +8*beta ));
  mu_assert("New net probs[2] ≠ exp(+4β)" , mynet.probs[2] == exp(  +4*beta ));
  mu_assert("New net probs[3] ≠ exp(0β)"  , mynet.probs[3] == exp(   0*beta ));
  mu_assert("New net probs[4] ≠ exp(-4β)" , mynet.probs[4] == exp(  -4*beta ));
  mu_assert("New net probs[5] ≠ exp(-8β)" , mynet.probs[5] == exp(  -8*beta ));
  mu_assert("New net probs[6] ≠ exp(-12β)", mynet.probs[6] == exp( -12*beta ));
  mu_assert("New net has mc steps!", mynet.mc_steps == 0);
  mu_assert("New net does not have the intended L", mynet.L == theL);
  metropolis(&mynet);
  mu_assert("Net does not keep track of the mc steps!", mynet.mc_steps == 1);
  free_spinglass(&mynet);
  return 0;
}

char* local_energy_test(){
  // If all the guides/spins are one, all the local E terms are 6.
  struct net SG = get_spinglass(theL, 0);

  // 0 110 ⋯ 110 110 110
  const uint64_t SIX_BY_TRIPLETS = 0x6db6db6db6db6db6;
  const uint64_t ZERO = 0x0000000000000000;

  for(int64_t i=0;i<theL*theL*theL;i++){
    SG.spins[i]   = ONES_BY_TRIPLETS;
    SG.J_up[i]    = ONES_BY_TRIPLETS;
    SG.J_right[i] = ONES_BY_TRIPLETS;
    SG.J_front[i] = ONES_BY_TRIPLETS;
  }
  int64_t idx = 0;
  for(int64_t z=0;z<theL;z++){
    for(int64_t y=0;y<theL;y++){
      for(int64_t x=0;x<theL;x++){
        mu_assert("∑ⱼ σᵢ⊕Jᵢⱼ⊕σⱼ ≠ 6 with σᵢ,Jᵢ=1",
                  local_energy(&SG,idx,x,y,z) == SIX_BY_TRIPLETS);
        idx++;
      }
    }
  }
  // If now all the J's are -1, the spins are all in a energy maxima.
  // The local E terms will return the index 0.
  for(int64_t i=0;i<theL*theL*theL;i++){
    SG.J_up[i]    = ZERO;
    SG.J_right[i] = ZERO;
    SG.J_front[i] = ZERO;
  }
  idx = 0;
  for(int64_t z=0;z<theL;z++){
    for(int64_t y=0;y<theL;y++){
      for(int64_t x=0;x<theL;x++){
        mu_assert("∑ⱼ σᵢ⊕Jᵢⱼ⊕σⱼ ≠ 0 with σᵢ=1, Jᵢ=0",
                  local_energy(&SG,idx,x,y,z) == ZERO);
        idx++;
      }
    }
  }
  free_spinglass(&SG);
  return 0;
}

char* metrotest(){
  //
  //
  // Test 1:
  //
  // Metropolis at T=0 on a unitary ferromagnet gives you the same
  // net.
  //
  struct net mynet = get_spinglass(theL,1e+8);
  int64_t L = mynet.L;
  for(int64_t i=0;i<L*L*L;i++){
    // Build a unitary ferromagnet
    mynet.J_up[i]    = ONES_BY_TRIPLETS;
    mynet.J_right[i] = ONES_BY_TRIPLETS;
    mynet.J_front[i] = ONES_BY_TRIPLETS;
    mynet.spins[i]   = ONES_BY_TRIPLETS;
  }
  // Evolve and compare
  metropolis(&mynet);
  for(int64_t i=0;i<L*L*L;i++){
    mu_assert("Metropolis at unitary net (T=0) moved things"
              ,mynet.spins[i] == ONES_BY_TRIPLETS);
  }
  //
  //
  // Test 2
  //
  // Metropolis at T=∞ just changes all the spins, for ferromagnets,
  // spin glasses or whatever.
  //
  free_spinglass(&mynet);
  mynet = get_spinglass(theL,1e-8);
  struct net orig = copy_net(mynet); // a random spin glass
  // Evolve 5 times, so we get a global inversion + 4 noops
  metropolis(&mynet);
  metropolis(&mynet);
  metropolis(&mynet);
  metropolis(&mynet);
  metropolis(&mynet);
  for(int64_t i=0;i<L*L*L;i++){
    mu_assert("Metropolis at T=∞ does not change all the spins"
              ,mynet.spins[i] == (orig.spins[i] ^ ONES_BY_TRIPLETS) // a = -b
              );
  }
  // Go to the beggining again
  metropolis(&mynet);
  for(int64_t i=0;i<L*L*L;i++){
    mu_assert("Metropolis at T=∞ does not cycle"
              ,mynet.spins[i] == orig.spins[i]
              );
  }
  free_spinglass(&mynet);
  free_spinglass(&orig);
  return 0;
}

char* schwinger_dyson_test(){ // Needs user inspection
  double beta = 0.5;
  struct net mynet = get_spinglass(theL,beta);
  for(int i=0;i<100;i++){
    metropolis(&mynet);
  }
  double mn=-1;
  double std=-1;
  schwingerdyson(&mynet,&mn,&std);
  printf("    %sSchw.Dyson%s    = %.2lf ± %.2lf  ≃ 1\n",
         BOLD, ENDC, mn, std);
  free_spinglass(&mynet);
  return 0;
}

char* net_measure_test(){

  struct net mynet = get_spinglass(theL,1);
  int L = mynet.L;
  for(int i=0;i<L*L*L;i++){
    mynet.spins[i] = ONES_BY_TRIPLETS;
  }
  mu_assert("M of unitary spinglass is not 1!"
            , fabs(1 - mean_magnetization(&mynet)) < EPSILON);
  for(int i=0;i<L*L*L;i++){
    mynet.J_up[i]    = ONES_BY_TRIPLETS;
    mynet.J_right[i] = ONES_BY_TRIPLETS;
    mynet.J_front[i] = ONES_BY_TRIPLETS;
    mynet.spins[i]   = ONES_BY_TRIPLETS;
  }
  mu_assert("E of unitary ferromagnet is not -1!"
            , fabs(-1 - mean_energy(&mynet)) < EPSILON);
  mu_assert("M of unitary ferromagnet is not 1!"
            , fabs(1 - mean_magnetization(&mynet)) < EPSILON);
  for(int i=0;i<L*L*L;i++){
    mynet.J_up[i]    = 0x0;
    mynet.J_right[i] = 0x0;
    mynet.J_front[i] = 0x0;
    mynet.spins[i]   = 0x0;
  }
  mu_assert("E of unitary antiferromagnet is not 1!"
            , fabs(1 - mean_energy(&mynet)) < EPSILON);
  mu_assert("M of unitary antiferromagnet is not -1!"
            , fabs(-1 - mean_magnetization(&mynet)) < EPSILON);
  free_spinglass(&mynet);
  return 0;
}

char* all_tests(){
  mu_run_test(utilities_test);
  mu_run_test(random_numbers_test);
  mu_run_test(pbc_test);
  mu_run_test(shift_spin_test);
  mu_run_test(initialization_test);
  mu_run_test(local_energy_test);
  mu_run_test(metrotest);
  mu_run_test(schwinger_dyson_test);
  mu_run_test(net_measure_test);
  return 0;
}

///////////////////////////////////////
//                                   //
//          COLLECT RESULTS          //
//                                   //
///////////////////////////////////////
int main(void){


  // Very important! seed the generator.
  seed_mersenne(time(NULL));


  clock_t start = clock(), diff;
  char* result = all_tests();
  if (result != 0){
    printf(
           "\n────────────────────────────────────────────────────────────\n"
           "\n%s%s✗ Error:%s %s\n"
           "\n────────────────────────────────────────────────────────────\n"
           , BOLD, RED
           , ENDC
           , result
           );
  }
  else{
    printf("%s%sAll deterministic tests passed ✓%s\n"
           ,GREEN
           ,BOLD
           ,ENDC
           );
  }
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("%d tests run in %.3lf seconds.\n"
         , tests_run
         , msec/1000.0
         );

  return result != 0;
}
