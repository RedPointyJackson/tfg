#!/usr/bin/julia

using DataFrames

srand(46)

SAMPLES = 100
ITERS=200
ITERS_p=200

values = readdlm("SDvalues.dat") |> vec
# A = 5 + randn(10000)
# B = -5 + randn(10000)
# values = vcat(A,B)
info("Real mean: $(mean(values))")

function bootstrap(array, frac; N=1000)
    L = length(array)
    chunksize = floor(Int64,frac*L)
    0 < chunksize ≤ L || error("chunksize = $chunksize ∉ (0,L)")
    means = zeros(N)
    for i ∈ 1:N
        chunk = rand(array,chunksize)
        means[i] = mean(chunk)
    end
    return mean(means), std(means)
end

function jacknife(array)
    L = length(array)
    means = zeros(L)
    for i ∈ 1:L
        means[i] = mean(array[j] for j in 1:L if j≠i)
    end
    return mean(means), std(means)
end

jdf = DataFrame()
bdf = DataFrame()


fr = 1/e
for i in 1:ITERS
    smallsample = rand(values,SAMPLES)
    jackμ, jackσ = jacknife(smallsample)
    bootμ, bootσ = bootstrap(smallsample, fr)
    jdf = vcat(jdf, DataFrame(i=i,
                              μ=jackμ,
                              σ=jackσ))
    bdf = vcat(bdf, DataFrame(i=i,
                              fraction=fr,
                              μ=bootμ,
                              σ=bootσ))
end

writetable("data_jack.csv", jdf)
writetable("data_boot.csv", bdf)
