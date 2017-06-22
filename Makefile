MTEXP=607 # Mersenne Twister exponent. Don't change, black magic.

WARNINGFLAGS = -Wall -Wextra -Wshadow -Werror -pedantic-errors -Wno-missing-field-initializers -fno-strict-aliasing -Wstrict-overflow -Wno-missing-braces
CCFLAGS= -std=c99 ${WARNINGFLAGS} -DSFMT_MEXP=${MTEXP}
OPTFLAGS= -O2 -ffast-math -march=native
DEBUGFLAGS = -g
MTFLAGS= -std=c99 -O3 -DSFMT_MEXP=${MTEXP}
MTDIR= "lib/mersennetlib"

CC="clang"

all: annealer measure

mt:
	@echo "Compiling mersenne twister"
	@${CC} ${MTFLAGS} -c ${MTDIR}/SFMT.c -o ${MTDIR}/SFMT.o

annealer: mt
	@echo "Compiling annealer.c"
	@${CC} ${CCFLAGS} ${OPTFLAGS} -c src/annealer.c -o src/annealer.o
	@mkdir -p exe
	@${CC} -lm -o exe/annealer src/annealer.o ${MTDIR}/SFMT.o

profile: annealer
	@echo "Compiling annealer with profile flags."
	@${CC} ${CCFLAGS} ${DEBUGFLAGS} -pg -c src/annealer.c -o src/annealer.o
	@mkdir -p exe
	@${CC} ${DEBUGFLAGS} -pg -lm -o exe/annealer src/annealer.o ${MTDIR}/SFMT.o
	@echo "Running..."
	@./exe/annealer -i 1e4 -m 10 42 /dev/null
	@gprof -p ./exe/annealer | less

benchmark: annealer
	@./src/benchmark.jl

measure: mt
	@echo "Compiling measure.c"
	@${CC} ${CCFLAGS} ${OPTFLAGS} -c src/measure.c -o src/measure.o
	@mkdir -p exe
	@${CC} -lm -o exe/measure src/measure.o ${MTDIR}/SFMT.o

test: annealer measure
	@${CC} ${CCFLAGS} ${DEBUGFLAGS} -c test/runtests.c -o test/runtests.o
	@mkdir -p exe
	@${CC} ${DEBUGFLAGS} -lm test/runtests.o ${MTDIR}/SFMT.o -o exe/runtests
	@echo -e "\n\nRunning glassy.c tests:"
	@echo    "---------------------------------------"
	@valgrind -q --leak-check=yes ./exe/runtests
	@echo -e "\n\nRunning annealer.c tests:"
	@echo    "---------------------------------------"
	@julia --color=yes test/testannealer.jl
	@echo -e "\n\nRunning measure.c tests:"
	@echo    "---------------------------------------"
	@julia --color=yes test/testmeasure.jl
	@echo -e "\n\nRunning continuation test:"
	@echo    "---------------------------------------"
	@./test/testcontinuation.sh

clean:
	rm -f observables.png
	rm -f fit_observables.png fit_observables.pdf
	rm -f fit.log
	rm -f *.{gcov,gcno,gcda}
	rm -f **/*.{gcov,gcno,gcda}
	rm -f **/*.o
	rm -f vgcore*
	rm -f *.net
	rm -f gmon.out


.PHONY: all mt annealer profile benchmark measure test clean
