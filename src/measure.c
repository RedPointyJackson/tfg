// Use --help after compilation to see how this program it works.
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <endian.h>
#include <argp.h>

#include "../lib/glassy.c"

#define IS_BIG_ENDIAN (!*(unsigned char *)&(uint16_t){1})

//
//    Things to parse command line arguments
//
const char *argp_program_version =
  "measure 5.0";

const char *argp_program_bug_address =
  "<a.clavero.alvarez@protonmail.com>";

// Doc for mandatory arguments.
char args_doc[] = "RA₁[, RB₁, RA₂, RB₂, RA₃, RB₃, …]";

// Program documentation.
char doc[] =
  "measure — Measures things in the output of annealer.c and prints in "
  "STDOUT a CSV with the measurements. The input, if relevant, is "
  "supposed to consist in ordered pairs of replicas.\v"

  "The output CSV has autoexplicative fields. If the observables include "
  "parameters, in C(tw,t) or C₄(r,tw), its respective parameters t,r "
  "will be included in the field \"parameter\". If it doesn't, any value "
  "will be used instead.\n\n"

  "If used without flags, it will go through the files anyway. Use "
  "--progress to use the software as a format checker for the files."
  ;

// This structure is used by main to communicate with parse_opt.
struct arguments{
  char**   files;
  int64_t n_files;
  bool    progress;
  bool    energy;
  bool    magnetization;
  bool    schwingerdyson;
  bool    correlation;
  bool    spatialcorrelation;
  bool    overlap;
};

// Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
struct argp_option options[] =
  {
    {"progress", 'p', 0, 0, "Show progress in stderr."},
    {"energy", 'e', 0, 0, "Measure ⟨-Jᵢⱼσᵢσⱼ⟩/6 ∈ [-1,1]."},
    {"magnetization", 'm', 0, 0, "Measure ⟨σᵢ⟩ ∈ [-1,1]."},
    {"schwingerdyson", 's', 0, 0, "Bootstrap the Schwinger-Dyson equation."},
    {"correlation", 'c', 0, 0, "Measure C(tw,tw+t) ∈ [-1,1]."},
    {"spatial-correlation", 'C', 0, 0, "Measure C₄(r,tw) ∈ [-1,1]."},
    {"overlap", 'q', 0, 0, "Measure ⟨σᵢᵃσⱼᵇ⟩ ∈ [-1,1]."},
    {0} // The struct must end in an option with all 0.
  };


// Order of parameters: KEY, ARG, STATE.
error_t parse_opt (int key, char *arg, struct argp_state *state){
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'p':
      arguments->progress = true;
      break;

    case 'e':
      arguments->energy = true;
      break;

    case 'm':
      arguments->magnetization = true;
      break;

    case 's':
      arguments->schwingerdyson = true;
      break;

    case 'c':
      arguments->correlation = true;
      break;

    case 'C':
      arguments->spatialcorrelation = true;
      break;

    case 'q':
      arguments->overlap = true;
      break;

    case ARGP_KEY_ARG:
      // if (state->arg_num < MAX_INPUT_FILES){
        arguments->n_files = state->arg_num + 1;
        arguments->files = realloc(arguments->files, arguments->n_files * sizeof(*arguments->files));
        if(!arguments->files) DIE("Can't realloc memory for arguments. Do you have a potato?\n");
        arguments->files[state->arg_num] = arg;
      // } else {
      //   fprintf(stderr
      //           ,"Maximum number of files allowed is %d (check source)\n"
      //           , MAX_INPUT_FILES);
      //   argp_usage(state);
      // }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1){
          argp_usage (state);
        }
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

// The ARGP structure itself.
struct argp argp = {options, parse_opt, args_doc, doc};

int64_t num_equalspins(int64_t spinpack1, int64_t spinpack2){
  // Correlation is proportional to the number of equal spins. As a
  // workaround, we use an exclusive OR (⊕=¬EQV) to get the number of
  // differing spins. Then, we just substract from the number of total
  // spins.
  int64_t different_count = 0;
  uint64_t AxorB = spinpack1 ^ spinpack2;
  for(int64_t spinidx = 0; spinidx<21; spinidx++){
    uint64_t sel_mask = (uint64_t)0x7 << 3*spinidx;
    uint64_t the_term = AxorB & sel_mask;
    if (the_term != 0x0){ different_count++; }
  }
  return 21 - different_count;
}

void measure_energy(struct anneal_contents* AC){
  int64_t V = AC->L * AC->L * AC->L;

  struct net thisSG = get_spinglass(AC->L, 1/AC->temp);
  for(int64_t spinidx=0; spinidx<V; spinidx++){
    thisSG.J_up[spinidx]    = AC->Jup[spinidx];
    thisSG.J_right[spinidx] = AC->Jright[spinidx];
    thisSG.J_front[spinidx] = AC->Jfront[spinidx];
  }

  char metadata[1024];
  for(int64_t mc_idx=0; mc_idx<AC->Nmeas; mc_idx++){

    for(int64_t spinidx=0; spinidx<V; spinidx++){
      thisSG.spins[spinidx] = AC->spins_list[mc_idx*V + spinidx];
    }

    sprintf(metadata,
            "%"PRId64",%"PRId64",%.8lf,%"PRIx32,
            AC->L,
            AC->mc_list[mc_idx],
            AC->temp,
            hash_guides(AC->Jup, AC->Jright, AC->Jfront, AC->L));

    printf("%s,%s,%"PRId64",%.8lf\n",
           metadata,
           "\"Energy\"",
           (int64_t)0,
           mean_energy(&thisSG));

  }
  free_spinglass(&thisSG);
}

void measure_magnetization(struct anneal_contents* AC){
  int64_t V = AC->L * AC->L * AC->L;

  struct net thisSG = get_spinglass(AC->L, 1/AC->temp);
  for(int64_t spinidx=0; spinidx<V; spinidx++){
    thisSG.J_up[spinidx]    = AC->Jup[spinidx];
    thisSG.J_right[spinidx] = AC->Jright[spinidx];
    thisSG.J_front[spinidx] = AC->Jfront[spinidx];
  }

  char metadata[1024];
  for(int64_t mc_idx=0; mc_idx<AC->Nmeas; mc_idx++){

    for(int64_t spinidx=0; spinidx<V; spinidx++){
      thisSG.spins[spinidx] = AC->spins_list[mc_idx*V + spinidx];
    }

    sprintf(metadata,
            "%"PRId64",%"PRId64",%.8lf,%"PRIx32,
            AC->L,
            AC->mc_list[mc_idx],
            AC->temp,
            hash_guides(AC->Jup, AC->Jright, AC->Jfront, AC->L));

    printf("%s,%s,%"PRId64",%.8lf\n",
           metadata,
           "\"Magnetization\"",
           (int64_t)0,
           mean_magnetization(&thisSG));

  }
  free_spinglass(&thisSG);
}

void measure_schwinger_dyson(struct anneal_contents* AC){
  int64_t V = AC->L * AC->L * AC->L;

  struct net thisSG = get_spinglass(AC->L, 1/AC->temp);
  for(int64_t spinidx=0; spinidx<V; spinidx++){
    thisSG.J_up[spinidx]    = AC->Jup[spinidx];
    thisSG.J_right[spinidx] = AC->Jright[spinidx];
    thisSG.J_front[spinidx] = AC->Jfront[spinidx];
  }

  char metadata[1024];
  for(int64_t mc_idx=0; mc_idx<AC->Nmeas; mc_idx++){

    for(int64_t spinidx=0; spinidx<V; spinidx++){
      thisSG.spins[spinidx] = AC->spins_list[mc_idx*V + spinidx];
    }

    sprintf(metadata,
            "%"PRId64",%"PRId64",%.8lf,%"PRIx32,
            AC->L,
            AC->mc_list[mc_idx],
            AC->temp,
            hash_guides(AC->Jup, AC->Jright, AC->Jfront, AC->L));

    double SDmean, SDstd;
    schwingerdyson(&thisSG, &SDmean, &SDstd);
    printf("%s,%s,%"PRId64",%.8lf\n", metadata, "\"SD mean\"" , (int64_t)0, SDmean);
    printf("%s,%s,%"PRId64",%.8lf\n", metadata, "\"SD stdev\"", (int64_t)0, SDstd);

  }
  free_spinglass(&thisSG);
}

void measure_correlation(struct anneal_contents* AC){

  if ((AC->ts == 1) && (AC->t_list[0]==0))
    DIE("Can't measure correlation: no t₀ points.");

  int64_t V = AC->L * AC->L * AC->L;

  int64_t tw_list_idx = 0;
  int64_t t_list_idx  = 0;
  char metadata[1024];

  for(int64_t mc_idx=0; mc_idx<AC->Nmeas; mc_idx++){
    int64_t this_tw = AC->mc_list[mc_idx];
    if (this_tw == AC->tw_list[tw_list_idx]){ // Found tw

      t_list_idx = 0;
      for(int64_t mmc_idx=mc_idx+1; mmc_idx<AC->Nmeas; mmc_idx++){
        int64_t this_t = AC->mc_list[mmc_idx] - AC->mc_list[mc_idx];
        if (this_t == AC->t_list[t_list_idx]){ // Found t

          sprintf(metadata,
                  "%"PRId64",%"PRId64",%.8lf,%"PRIx32,
                  AC->L,
                  AC->mc_list[mc_idx],
                  AC->temp,
                  hash_guides(AC->Jup, AC->Jright, AC->Jfront, AC->L));

          double corr = 0;
          for(int64_t spinidx=0; spinidx<V; spinidx++){
            uint64_t spinA = AC->spins_list[mc_idx*V + spinidx];
            uint64_t spinB = AC->spins_list[mmc_idx*V + spinidx];
            int64_t Neq = num_equalspins(spinA, spinB);
            corr += Neq - (21-Neq);
          }

          printf("%s,%s,%"PRId64",%.8lf\n", metadata,
                 "\"Correlation\"", this_t, corr / 21 / V);

          if (t_list_idx == AC->ts-1) break; else t_list_idx++;
        }

      }
      if (tw_list_idx == AC->tws-1) break; else tw_list_idx++;
    }
  }

}


void measure_overlap(struct anneal_contents* AC_A, struct anneal_contents* AC_B){
  int64_t V = AC_A->L * AC_A->L * AC_A->L;

  char metadata[1024];
  for(int64_t mc_idx=0; mc_idx<AC_A->Nmeas; mc_idx++){

    sprintf(metadata,
            "%"PRId64",%"PRId64",%.8lf,%"PRIx32,
            AC_A->L,
            AC_A->mc_list[mc_idx],
            AC_A->temp,
            hash_guides(AC_A->Jup, AC_A->Jright, AC_A->Jfront, AC_A->L));

    double overlap = 0; // Just the correlation between samples.
    for(int64_t spinidx=0; spinidx<V; spinidx++){
      uint64_t spinA = AC_A->spins_list[mc_idx*V + spinidx];
      uint64_t spinB = AC_B->spins_list[mc_idx*V + spinidx];
      int64_t Neq = num_equalspins(spinA, spinB);
      overlap += Neq - (21-Neq);
    }

    printf("%s,%s,%"PRId64",%.8lf\n", metadata,
           "\"Overlap\"" , (int64_t)0, overlap / 21 / V);

  }
}

void measure_spatial_correlation(struct anneal_contents* AC_A, struct anneal_contents* AC_B){
  int64_t V = AC_A->L * AC_A->L * AC_A->L;

  init_guides_pbc(AC_A->L);

  char metadata[1024];

  for(int64_t r=0;r<AC_A->L;r++){
    for(int64_t mc_idx=0; mc_idx<AC_A->Nmeas; mc_idx++){

      sprintf(metadata,
              "%"PRId64",%"PRId64",%.8lf,%"PRIx32,
              AC_A->L,
              AC_A->mc_list[mc_idx],
              AC_A->temp,
              hash_guides(AC_A->Jup, AC_A->Jright, AC_A->Jfront, AC_A->L));

      double C4 = 0;

      int64_t spinidx = 0;
      for(int64_t z=0;z<AC_A->L;z++){
        for(int64_t y=0;y<AC_A->L;y++){
          for(int64_t x=0;x<AC_A->L;x++){
            // First, q_i
            int64_t qi = 0;

            uint64_t spinA = AC_A->spins_list[mc_idx*V + spinidx];
            uint64_t spinB = AC_B->spins_list[mc_idx*V + spinidx];
            int64_t Neq = num_equalspins(spinA, spinB);
            qi += Neq - (21-Neq);

            int64_t qi_plus_r = 0;

            uint64_t spinA_r = get_shifted_spin(AC_A->spins_list + mc_idx*V, AC_A->L, spinidx, r, RIGHT);
            uint64_t spinB_r = get_shifted_spin(AC_B->spins_list + mc_idx*V, AC_B->L, spinidx, r, RIGHT);
            Neq = num_equalspins(spinA_r, spinB_r);
            qi_plus_r += (Neq - (21-Neq));

            uint64_t spinA_l = get_shifted_spin(AC_A->spins_list + mc_idx*V, AC_A->L, spinidx, r, LEFT);
            uint64_t spinB_l = get_shifted_spin(AC_B->spins_list + mc_idx*V, AC_B->L, spinidx, r, LEFT);
            Neq = num_equalspins(spinA_l, spinB_l);
            qi_plus_r += (Neq - (21-Neq));

            uint64_t spinA_u = get_shifted_spin(AC_A->spins_list + mc_idx*V, AC_A->L, spinidx, r, UP);
            uint64_t spinB_u = get_shifted_spin(AC_B->spins_list + mc_idx*V, AC_B->L, spinidx, r, UP);
            Neq = num_equalspins(spinA_u, spinB_u);
            qi_plus_r += (Neq - (21-Neq));

            uint64_t spinA_d = get_shifted_spin(AC_A->spins_list + mc_idx*V, AC_A->L, spinidx, r, DOWN);
            uint64_t spinB_d = get_shifted_spin(AC_B->spins_list + mc_idx*V, AC_B->L, spinidx, r, DOWN);
            Neq = num_equalspins(spinA_d, spinB_d);
            qi_plus_r += (Neq - (21-Neq));

            uint64_t spinA_f = get_shifted_spin(AC_A->spins_list + mc_idx*V, AC_A->L, spinidx, r, FRONT);
            uint64_t spinB_f = get_shifted_spin(AC_B->spins_list + mc_idx*V, AC_B->L, spinidx, r, FRONT);
            Neq = num_equalspins(spinA_f, spinB_f);
            qi_plus_r += (Neq - (21-Neq));

            uint64_t spinA_b = get_shifted_spin(AC_A->spins_list + mc_idx*V, AC_A->L, spinidx, r, BEHIND);
            uint64_t spinB_b = get_shifted_spin(AC_B->spins_list + mc_idx*V, AC_B->L, spinidx, r, BEHIND);
            Neq = num_equalspins(spinA_b, spinB_b);
            qi_plus_r += (Neq - (21-Neq));

            C4 += qi/21.0 * qi_plus_r/6.0/21.0;

            spinidx++;
          }
        }
      }

      printf("%s,%s,%"PRId64",%.8lf\n", metadata,
             "\"Spatial_correlation\"" , r, C4 / V);
    }

  }
}

void check_if_they_are_replicas(struct anneal_contents* AC_A, struct anneal_contents* AC_B){

  int64_t V = AC_A->L * AC_A->L * AC_A->L;

  if (AC_A->L != AC_B->L) DIE("Replicas don't share L.\n");
  if (AC_A->temp != AC_B->temp) DIE("Replicas don't share T.\n");

  for(int64_t i=0;i<V;i++){
    if (AC_A->Jup[i] != AC_B->Jup[i])       DIE("Replicas don't share J.\n");
    if (AC_A->Jright[i] != AC_B->Jright[i]) DIE("Replicas don't share J.\n");
    if (AC_A->Jfront[i] != AC_B->Jfront[i]) DIE("Replicas don't share J.\n");
  }

  if (AC_A->Nmeas != AC_B->Nmeas) DIE("Nmeas(A) ≠ Nmeas(B)\n");
  for(int64_t i=0;i<AC_A->Nmeas;i++){
    if (AC_A->mc_list[i] != AC_B->mc_list[i]) DIE("mcAᵢᵢ ≠ mcBᵢ");
  }

  if (AC_A->tws != AC_B->tws) DIE("tws(A) ≠ tws(B)\n");
  for(int64_t i=0;i<AC_A->tws;i++){
    if (AC_A->tw_list[i] != AC_B->tw_list[i]) DIE("twAᵢ ≠ twBᵢ");
  }

  if (AC_A->ts != AC_B->ts) DIE("tw(A) ≠ tw(B)\n");
  for(int64_t i=0;i<AC_A->ts;i++){
    if (AC_A->t_list[i] != AC_B->t_list[i]) DIE("tAᵢ ≠ tBᵢ");
  }
}

// Main routine. Go through all the β's.
int main(int argc, char** argv){

  seed_mersenne(42);

  struct arguments args = {0};

  argp_parse(&argp, argc, argv, 0, 0, &args);

  if (args.spatialcorrelation && args.n_files % 2 != 0)
    DIE("Number of files isn't even: can't match all the replicas.\n");

  if (args.overlap && args.n_files % 2 != 0)
    DIE("Number of files isn't even: can't match all the replicas.\n");

  printf("L,mc,T,Jhash,observable,parameter,value\n");

  struct anneal_contents AC_A; // replica A
  struct anneal_contents AC_B; // replica B

  for(int64_t i=0;i<args.n_files;i++){

    if (args.progress){
      fprintf(stderr, "%s (%"PRId64"/%"PRId64"):",
              args.files[i],i+1,args.n_files);
    }

    FILE* fA = fopen(args.files[i], "rb");
    if (fA==NULL) {DIE("Failed to read contents from file.\n");}
    AC_A = get_anneal_contents(fA);
    fclose(fA);

    if (args.progress){
      fprintf(stderr , "    ");
      fprintf(stderr , "L%"PRId64" | "     , AC_A.L);
      fprintf(stderr , "%"PRId64" tw's | " , AC_A.tws);
      fprintf(stderr , "%"PRId64" t's | "  , AC_A.ts);
      fprintf(stderr , "%"PRId64" meas | " , AC_A.Nmeas);
      fprintf(stderr , "T %lf\n"           , AC_A.temp);
    }

    // First, the measurements that doesn't require replicas
    if (args.energy)         measure_energy(&AC_A);
    if (args.magnetization)  measure_magnetization(&AC_A);
    if (args.schwingerdyson) measure_schwinger_dyson(&AC_A);
    if (args.correlation)    measure_correlation(&AC_A);

    // Then, the ones that does. They only are relevant every two
    // files, because they need two nets.
    if (i%2 == 0 && (args.spatialcorrelation || args.overlap)){
      FILE* fB = fopen(args.files[i+1], "rb");
      if (fB==NULL) {DIE("Failed to read contents from file.\n");}
      AC_B = get_anneal_contents(fB);
      fclose(fB);

      check_if_they_are_replicas(&AC_A, &AC_B);

      // measure
      if (args.spatialcorrelation) measure_spatial_correlation(&AC_A, &AC_B);
      if (args.overlap)            measure_overlap(&AC_A, &AC_B);

      free_anneal_contents(&AC_B);
    }


    free_anneal_contents(&AC_A);
  }

  if (args.progress){
    fprintf(stderr, "All files were analyzed.           \n");
  }

  return EXIT_SUCCESS;
}
