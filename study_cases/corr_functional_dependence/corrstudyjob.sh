# Best served with L=8

SAMPLES="10"

DIR="study_cases/corr_functional_dependence"

cd ../../
make annealer > /dev/null
make measure > /dev/null

T="0.9"

trange="-t 1 -T 1000 -n 50 -L"
twrange="-i 1e3 -m 50"
L="32"

##################################
#                                #
#             ANNEAL             #
#                                #
##################################

rm -f ${DIR}/*.net

for i in $(seq ${SAMPLES}); do
    echo "Doing run ${i} of ${SAMPLES}"
    ./exe/annealer -l $L ${trange} ${twrange} ${T} ${DIR}/anneal_${i}.net
done

##################################
#                                #
#           ANALYSIS             #
#                                #
##################################
echo "Measuring ..."

./exe/measure -pc ${DIR}/anneal_*.net > ${DIR}/anneal.csv
rm -f ${DIR}/anneal_*.net

# Make a custom plot to show interesting effects.
cd ${DIR}
./customplot
