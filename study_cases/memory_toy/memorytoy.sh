SAMPLES="500"

DIR="study_cases/memory_toy/"

cd ../../
make annealer > /dev/null
make measure > /dev/null

TA="0.9"
TB="0.5"
t="-t 50 -n 1"
L="12"

##################################
#                                #
#             ANNEAL             #
#                                #
##################################

rm -f ${DIR}/*.net

for i in $(seq ${SAMPLES}); do
    echo -e "Doing run ${i} of ${SAMPLES}"
    # T1 quench
    ./exe/annealer -l $L ${t} -i 2e3 -m 1e3 ${TA} ${DIR}/anneal_A_${i}.net
    # Continue T1 quench
    cat ${DIR}/anneal_A_${i}.net | ./exe/annealer -l $L ${t} -c 1e3 -i 1e3 -m 1e3 ${TB} ${DIR}/anneal_B_${i}.net
    # Return to T1
    cat ${DIR}/anneal_B_${i}.net | ./exe/annealer -l $L ${t} -c 1e3 -i 1e3 -m 1e3 ${TA} ${DIR}/anneal_C_${i}.net
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

echo "Measuring C..."

./exe/measure -c ${DIR}/anneal_C_*.net > ${DIR}/annealC.csv
rm -f ${DIR}/anneal_C_*.net

# Make a custom plot to show interesting effects.
cd ${DIR}
./customplot
