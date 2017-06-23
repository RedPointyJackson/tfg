SAMPLES="1500"

DIR="study_cases/vanishingmemory/"

cd ../../
make annealer > /dev/null
make measure > /dev/null

TA="0.9"
TB="0.5"

trange="-t 1 -T 1000 -n 100 -L"

##################################
#                                #
#             ANNEAL             #
#                                #
##################################

rm -f ${DIR}/*.net

for i in $(seq ${SAMPLES}); do
    echo "Doing run ${i} of ${SAMPLES}"
    ./exe/annealer ${trange} -i 1e3 -m 1e3 ${TA} ${DIR}/anneal_A_${i}.net
    cat ${DIR}/anneal_A_${i}.net | ./exe/annealer ${trange} -c 1e3 -i 1e3 -m 1e3 ${TB} ${DIR}/anneal_B_${i}.net
done

##################################
#                                #
#           ANALYSIS             #
#                                #
##################################
echo "Measuring A..."

./exe/measure -c ${DIR}/anneal_A_*.net > ${DIR}/annealA.csv
rm -f ${DIR}/anneal_A_*.net

echo "Measuring B..."

./exe/measure -c ${DIR}/anneal_B_*.net > ${DIR}/annealB.csv
rm -f ${DIR}/anneal_B_*.net

# Make a custom plot to show interesting effects.
cd ${DIR}
./customplot
