#!/usr/bin/julia --color=yes

using DataFrames
using ProgressMeter

const iters = Int64(1e5)
const meas = 300
const samples = 1

cd("../../")
run(`make annealer`)
run(`make measure`)
cd("study_cases/sergio_energy")

# Read all of sergio's data
df = DataFrame()

folders = ["L08","L16"]

for fol in folders
    L = parse(Int64,fol[2:end])
    for t in readdir(fol)
        T = parse(Float64,t[2:end])
        for subfol in readdir("$fol/$t")
            for file in readdir("$fol/$t/$subfol")
                data = readdlm("$fol/$t/$subfol/$file")
                mc = data[:,2]
                Ea = data[:,7] # replica A
                Eb = data[:,8] # replica B
                df_A = DataFrame(
                                mc = mc,
                                physicist = "Sergio (replica A)",
                                E = Ea/3,
                                T = T,
                                L = L,
                                )
                df_B = DataFrame(
                                mc = mc,
                                physicist = "Sergio (replica B)",
                                E = Eb/3,
                                T = T,
                                L = L,
                                )
                df = vcat(df,df_A)
                df = vcat(df,df_B)
            end
        end
    end
end

function replicate_experiment(L,T)
    datas = DataFrame()
    # Make the runs
    for i in 1:samples
        cmd = `../../exe/annealer -l $L -i $iters -m $meas $T run.net`
        run(cmd)
        # Analyze, save in the DataFrame
        run(pipeline(`../../exe/measure -e run.net`, stdout="run.meas"))
        data = readtable("run.meas")
        rm("run.meas")
        rm("run.net")
        μdf = DataFrame(
                        mc = data[:mc],
                        E = data[:value],
                        T = T,
                        L = L,
                        )
        datas = vcat(datas, μdf)
    end
    agg = aggregate(datas,[:mc,:T,:L],mean)
    rename!(agg,:E_mean,:E)
    agg[:physicist] = "Alejandro ($samples samples)"
    return agg
end

@showprogress for T in unique(df[:T]), L in unique(df[:L])
    df = vcat(df,replicate_experiment(L,T))
end

writetable("comparison.csv",df)

run(`./plot`)
