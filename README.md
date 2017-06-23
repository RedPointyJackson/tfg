TFG on spin glasses
===================

Utilities to simulate spin glasses, created for an undergrad
dissertation. Compile the main program with `make anneal` after
running `make test` to be sure everything's okay.
The main code is in the `/lib` folder, and is called `glassy.c`.

Raw data from the simulations is also available if wanted or needed,
but is terabytes (binary configurations) and gigabytes (parsed csv
measures) big, which makes it unsuitable for hosting in GitHub. If you
are interested, ask for it at `a.clavero.alvarez@protonmail.com`.

Requisites
------------------------
The sofware uses
- _C_ for the `annealer` program.
- _Julia 0.5_ for some scripts.
- _Python 3_ for the plotting facilities (`matplotlib` and `seaborn`
  libraries are extensively used). A custom style has been used for
  `matplotlib`, comment the `plt.style.use` lines or move the provided
  one in `/etc` to `.config/matplotlib/stylelib` to be able to use the
  same plotting code.

`annealer` has no dependencies; the Mersenne Twister library that it
employs is packaged inside the `lib/` folder and automatically linked.

Compilinig the LaTeX file (in Spanish) will require `minted` and
`lualatex`.

The code has been only tested on Linux systems.

Usage
------------------------
The main interface to the code is the `annealer` executable, in the
`exe` folder after doing a `make anneal`. Just run it with `--help` to
see the usage.

`measure` takes `annealer` files as input, and computes a csv with various
observables. Use `--help` to see a brief help message.

In the `study_cases` folder there are various scripts that can be used
to explore some situations or make the plots of the dissertation
(folder `tex`). Of
course, they need the original data to work and are displayed here
only for reference.
