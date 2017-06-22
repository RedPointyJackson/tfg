#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>
#include <argp.h>

#include "../lib/glassy.c"

#define IS_BIG_ENDIAN (!*(unsigned char *)&(uint16_t){1})




// Argument parsing



const char *argp_program_version =
  "annealer 5.0";

const char *argp_program_bug_address =
  "<a.clavero.alvarez@protonmail.com>";

char args_doc[] = "TEMP OUTFILE"; // Mandatory args doc

char doc[] =
  "annealer — Create a spin glass and make them evolve at a "
  "given temperature. Outputs a binary file to OUTFILE with "
  "the measured lattices. "
  "All the arguments are parsed as floating point numbers and "
  "floored, so one can specify 1e5 instead of 100000.\v"

  "For the provided T, the software will measure at every tw+t, where "
  "tw ∈ [1,N_ITERS] are M logarithmically spaced points and "
  "t ∈ [MIN_T,MAX_T] are NUM_T linearly spaced points."

  "\n\nHappy annealing! ♥"
  ;

struct arguments{
  double   temp;
  char* filename;

  int64_t  L;
  int64_t  continue_step;
  int64_t  iters;
  int64_t  measpoints;
  int64_t  mint;
  int64_t  maxt;
  int64_t  numt;

  bool     dry_run;
  bool     progress;
  bool     ferro;

  bool     t_log;

  uint32_t jseed;
};

struct argp_option options[] =
  {
    // {long name, short name, <ARGNAME>, flags, docs}
    {"side", 'l', "L", 0,
     "Side of the spin glass. It will have L³ spins in total."},

    {"continue", 'c', "TW", 0,
     "Use J couplings of the file provided in stdin, and spins of the "
     "first lattice found with mc ≥ TW. If TW is negative, the last value "
     "will be used. If TW is zero, the option will be ignored."},

    {"ferro", 'f', 0, 0,
     "If enabled, Jᵢ=1, ∀i."},

    {"iters", 'i', "N_ITERS", 0,
     "Metropolis sweeps to perform in each β."},

    {"j-seed", 'j', "SEED", 0, "Seed to use when creating the "
     "J couplings. If omited, /dev/urandom will be used."},

    {"meas-points", 'm', "M", 0, "Number of tw's to compute."},

    {"min-t", 't', "MIN_T₀", 0, "Minimum t₀ to compute."},

    {"max-t", 'T', "MAX_T₀", 0, "Maximum t₀ to compute."},

    {"num-t", 'n', "NUM_T", 0, "Number of t₀'s to compute. If "
     "NUM_T=1 will use only the minimum t₀, and if NUM_T=0 will "
     "not use t₀'s at all, measuring only at tw points. "
     "Must be ≥ 0."},

    {"t-log", 'L', 0, 0, "Use log-spaced t's."},

    {"dry-run", 'd', 0, 0,
     "Don't run anything, only parse command line arguments."},

    {"progress", 'p', 0, 0,
     "Output progress and other information to stdout."},

    {0} // The struct must end in an option with all 0.
  };


error_t parse_opt (int key, char *arg, struct argp_state *state){
  struct arguments *arguments = state->input;

  switch (key){
    case 'j':
      arguments->jseed = atoi(arg);
      break;
    case 'l':
      arguments->L = atof(arg);
      if (arguments->L<=1) DIE("-L must be > 1.\n");
      break;
    case 'L':
      arguments->t_log = true;
      break;
    case 'c':
      arguments->continue_step = atof(arg);
      if (arguments->continue_step<1) DIE("-c must be ≥ 1.\n");
      break;
    case 'f':
      arguments->ferro = true;
      break;
    case 'i':
      arguments->iters = atof(arg);
      if (arguments->iters<=1) DIE("-i must be > 1.\n");
      break;
    case 'm':
      arguments->measpoints = atof(arg);
      if (arguments->measpoints<=1) DIE("-m must be > 1.\n");
      break;
    case 't':
      arguments->mint = atof(arg);
      if (arguments->mint<1) DIE("-t must be ≥ 1.\n");
      break;
    case 'T':
      arguments->maxt = atof(arg);
      if (arguments->maxt<1) DIE("-T must be ≥ 1.\n");
      break;
    case 'n':
      arguments->numt = atof(arg);
      if (arguments->numt<0) DIE("-n must be ≥ 0.\n");
      break;
    case 'd':
      arguments->dry_run = true;
      break;
    case 'p':
      arguments->progress = true;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num == 0) arguments->temp = atof(arg);
      else if (state->arg_num == 1) arguments->filename = arg;
      else argp_usage(state);
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 2){
        argp_usage (state);
      }
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

struct argp argp = {options, parse_opt, args_doc, doc};





// Actual logic





int64_t measurements_done = 0;
void print_and_evolve(struct net* SG, FILE* sgfile){
  int64_t L = SG->L;
  int64_t mcs_written =  fwrite(&(SG->mc_steps), sizeof((SG->mc_steps)), 1, sgfile);
  int64_t spins_written = fwrite(SG->spins, sizeof((SG->spins[0])), L*L*L, sgfile);
  if (mcs_written != 1) DIE("mc steps write was corrupted.\n");
  if (spins_written != L*L*L) DIE("spins write was corrupted.\n");
  fflush(sgfile);
  metropolis(SG);
  measurements_done++;
}

uint32_t get_random_seed(void){
  uint32_t random_number;
  FILE* devurandom = fopen("/dev/urandom", "r");
  if (devurandom==NULL){
    DIE("Failed to open /dev/urandom.\n");
  }
  if(fread(&random_number, sizeof(uint32_t), 1, devurandom) != 1)
    DIE("Failed to read /dev/urandom.\n");
  fclose(devurandom);
  return random_number;
}

int main(int argc, char** argv){

  seed_mersenne(get_random_seed());


  struct arguments arg;
  arg.continue_step = 0;
  arg.ferro         = 0;
  arg.iters         = 1e4;
  arg.measpoints    = 10;
  arg.mint          = 2;
  arg.maxt          = 10;
  arg.numt          = 5;
  arg.dry_run       = 0;
  arg.progress      = 0;
  arg.L             = 8;
  arg.jseed         = get_random_seed();
  arg.t_log         = false;

  argp_parse(&argp, argc, argv, 0, 0, &arg);

  FILE* sgfile = fopen(arg.filename,"wb");
  if (!sgfile) DIE("Can't open output files.\n");

  if (arg.measpoints>arg.iters)
    DIE("-m must be ≤ -i.\n");
  if (arg.continue_step && arg.ferro)
    DIE("-c and -f are incompatible.\n");

  //
  // Build the net
  //

  struct net SG = get_spinglass(arg.L, 1.0/arg.temp);

  if (arg.continue_step){
    continue_net_from_stdin(arg.continue_step, &SG);
    if (arg.L != SG.L) DIE("Continuation file L differs from provided.\n");
  }
  else {
    for(int64_t i=0;i<pow(arg.L,3);i++){
      SG.spins[i] = random_bits();
    }
    seed_mersenne(arg.jseed);
    for(int64_t i=0;i<pow(arg.L,3);i++){
      SG.J_up[i]    = random_bits();
      SG.J_right[i] = random_bits();
      SG.J_front[i] = random_bits();
    }
  }

  seed_mersenne(get_random_seed());

  int64_t V = arg.L * arg.L * arg.L;

  if(arg.ferro){
    for(int64_t i=0;i<V;i++){
      SG.J_up[i]    = ONES_BY_TRIPLETS;
      SG.J_right[i] = ONES_BY_TRIPLETS;
      SG.J_front[i] = ONES_BY_TRIPLETS;
    }
  }


  //
  // Obtain tw's, t's
  //

  int64_t* tw_list = malloc(arg.measpoints * sizeof(*tw_list));
  logspace(tw_list, 1, arg.iters, arg.measpoints);
  int64_t tws = sort_unique(tw_list, arg.measpoints);

  int64_t* t_list = malloc(arg.numt * sizeof(*tw_list));
  int64_t ts;
  if (arg.numt == 0){
    t_list[0] = 0;
    ts = 1;
  }
  else if (arg.numt == 1){
    t_list[0] = arg.mint;
    ts = 1;
  }
  else if (arg.t_log){
    logspace(t_list, arg.mint, arg.maxt, arg.numt);
    ts = sort_unique(t_list, arg.numt);
  }
  else {
    linspace(t_list, arg.mint, arg.maxt, arg.numt);
    ts = sort_unique(t_list, arg.numt);
  }

  // m[i] = tw + t, ∀t,tw
  int64_t* meas_points = malloc(tws*(ts+1) * sizeof(*meas_points));
  size_t curridx = 0;
  for(int64_t tw=0;tw<tws;tw++){
    meas_points[curridx] = tw_list[tw];
    curridx++;
    for(int t=0;t<ts;t++){
      meas_points[curridx] = tw_list[tw] + t_list[t];
      curridx++;
    }
  }

  int64_t Nmeas = sort_unique(meas_points, tws*(ts+1));


  if (arg.dry_run){
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr ,"        Output    %s         \n" , arg.filename                 );
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr ,"    Big endian    %s         \n" , IS_BIG_ENDIAN?"true":"false" );
    fprintf(stderr ,"       dry_run    %s         \n" , arg.dry_run?"true":"false"   );
    fprintf(stderr ,"      progress    %s         \n" , arg.progress?"true":"false"  );
    fprintf(stderr ,"         ferro    %s         \n" , arg.ferro?"true":"false"     );
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr ,"          temp    %lf        \n" , arg.temp                     );
    fprintf(stderr ," continue_step    %"PRId64"  \n" , arg.continue_step            );
    fprintf(stderr ,"             L    %"PRId64"  \n" , arg.L                        );
    fprintf(stderr ,"         jseed    %"PRIu32"  \n" , arg.jseed                    );
    fprintf(stderr ,"         iters    %"PRId64"  \n" , arg.iters                    );
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr ,"          tw's    %"PRId64"  \n" , arg.measpoints               );
    fprintf(stderr ,"         mintw    %"PRId64"  \n" , tw_list[0]                   );
    fprintf(stderr ,"         maxtw    %"PRId64"  \n" , tw_list[arg.measpoints-1]    );
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr ,"           t's    %"PRId64"  \n" , arg.numt                     );
    fprintf(stderr ,"          mint    %"PRId64"  \n" , t_list[0]                    );
    fprintf(stderr ,"          maxt    %"PRId64"  \n" , t_list[arg.numt-1]           );
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr ,"         Nmeas    %"PRId64"  \n" , Nmeas                        );
    fprintf(stderr ,"\n"                                                             );
    fprintf(stderr , "Aborting due to --dry-run\n"                                   );
    return 0;
  }


  //
  // Binary header
  //

  int64_t items = 0;

  char beggining_of_file_marker[10] =  "ANNEALER";
  items += fwrite(beggining_of_file_marker, sizeof(char), 8, sgfile);

  int64_t endianness = IS_BIG_ENDIAN;
  items += fwrite(&endianness, sizeof(endianness), 1, sgfile);

  items += fwrite(&(SG.L), sizeof(SG.L), 1, sgfile);

  items += fwrite(&tws, sizeof(tws), 1, sgfile);
  items += fwrite(tw_list, sizeof(tw_list[0]), tws, sgfile);

  items += fwrite(&ts, sizeof(ts), 1, sgfile);
  items += fwrite(t_list, sizeof(t_list[0]), ts, sgfile);

  items += fwrite(&Nmeas, sizeof(Nmeas), 1, sgfile);

  items += fwrite(&(arg.temp), sizeof(arg.temp), 1, sgfile);

  items += fwrite(SG.J_up, sizeof(SG.J_up[0]), V, sgfile);
  items += fwrite(SG.J_right, sizeof(SG.J_right[0]), V ,sgfile);
  items += fwrite(SG.J_front, sizeof(SG.J_front[0]), V ,sgfile);

  if (items != 8 +1 + 1 + 1+tws + 1+ts + 1 + 1 + V+V+V)
    DIE("Header write was corrupted.\n");

  //
  // Actual measuring
  //

  int64_t itersdone = 0;
  int64_t start_time = time(NULL);

  for(int64_t midx=0; midx<Nmeas; midx++){

    while(itersdone < meas_points[midx]){
      metropolis(&SG); itersdone++;
    }

    print_and_evolve(&SG, sgfile); itersdone++;

    if (arg.progress){
      double SDmean = -1, SDstd = -1;
      schwingerdyson(&SG, &SDmean, &SDstd);
      double time_per_iter = (time(NULL)-start_time)*1.0/itersdone;
      double ETA = (meas_points[Nmeas-1]-itersdone)*time_per_iter / 3600;
      fprintf(stderr, "ETA %2.2lf h\t│\t%6.2lf%% (%"PRId64"/%"PRId64")\t│\t%s⟨SD⟩ = %.2lf ± %.2lf%s\n",
              ETA,
              100.0*(itersdone+1)/meas_points[Nmeas-1],
              midx+1, Nmeas,
              fabs(SDmean-1)/SDstd < 1.5 ? GREEN_C : RED_C,
              SDmean, SDstd, END_C);

    }

  }


  free_spinglass(&SG);
  free(meas_points);
  free(t_list);
  free(tw_list);


  char end_of_file_marker[10] =  "END";
  if ( 3 != fwrite(end_of_file_marker, sizeof(char), 3, sgfile))
    fprintf(stderr, WARNHEADER"Error writting binary footer.\n");

  fclose(sgfile);

  if (arg.progress){
    fprintf(stderr,"All done in %.2lf h! Have a lovely day.\n",
            (time(NULL)-start_time)/3600.0);
  }

  if (Nmeas != measurements_done){
    fprintf(stderr,
            WARNHEADER"Did %ld of the %ld scheduled measurings.\n",
            measurements_done, Nmeas);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
