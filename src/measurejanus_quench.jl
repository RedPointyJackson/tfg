#!/usr/bin/env julia

Base.eval(:(have_color=true))

function die(msg)
    print_with_color(:red, STDERR, "ERROR: " * msg * "\n")
    exit(1)
end

if length(ARGS) != 1
    println("Usage: ", PROGRAM_FILE, "<folder>")
    println("    where <folder> is the output of a run in Janus.")
    println("    The program will output in stdout a csv with the")
    println("    temporal correlation of the configurations.")
    quit()
end

folder = ARGS[1]

if folder[end] == '/'
    thefolder = folder |> dirname |> basename
else
    thefolder = folder |> basename
end

if !isfile(folder * "/S00001/run.log")
    die("No run.log in folder.")
end

Jhash = @sprintf("%d",read(folder * "/S00001/conf_Js") |> hash)

cd(folder * "/S00001/confs")

confs = readdir(".")

isconf(f) = ismatch(r"conf_[0-9A-F]{10}", f)

if !all(map(isconf, confs))
    warn("Some files don't match conf_\\x{10}, filtering them.")
    filter!(isconf, confs)
end

function uncompact(load)
    N = length(load)
    data = BitArray(N*8)
    for ib in 1:N
        is = 8(ib-1)
        for j in 0:7
            data[1+is+j] = (load[ib]>>j) & 0x1
        end
    end
    return data
end

function conf2data(file)
    data = Dict()
    open(file,"r") do f
        data["itime"]           = read(f,Cint)
        data["Dt"]              = read(f,Clonglong)
        data["t"]               = read(f,Clonglong) # The MC step
        data["timewindow"]      = read(f,Clonglong)
        data["flagJ"]           = read(f,Cint)
        data["flagS"]           = read(f,Cint)
        data["flagh"]           = read(f,Cint)

        data["L"]               = read(f,Cint)
        data["β"]               = read(f,Cdouble)
        data["H"]               = read(f,Cdouble)
        data["tHon"]            = read(f,Culonglong)

        data["ir"]              = read(f,Cint)
        data["nfiles_dir"]      = read(f,Cint)

        data["energy_sigma"]    = read(f,Cdouble)
        data["energy_tau"]      = read(f,Cdouble)
        data["energymag_sigma"] = read(f,Cdouble)
        data["energymag_tau"]   = read(f,Cdouble)
        data["mag_sigma"]       = read(f,Cdouble)
        data["mag_tau"]         = read(f,Cdouble)

        num                     = data["L"]^3 ÷ 8
        data["spins_σ"]         = uncompact(read(f,Cchar,num))
        data["spins_τ"]         = uncompact(read(f,Cchar,num))
        @assert eof(f)
    end
    return data
end

function getC(bconfA, bconfB)
    spinsA_σ = bconfA["spins_σ"]
    spinsA_τ = bconfA["spins_τ"]
    spinsB_σ = bconfB["spins_σ"]
    spinsB_τ = bconfB["spins_τ"]
    N = length(spinsA_σ)
    differing_spins_σ = map(≠,spinsA_σ, spinsB_σ) |> sum
    differing_spins_τ = map(≠,spinsA_τ, spinsB_τ) |> sum
    # C = 1*N↑↑ + (-1)*N↑↓
    Cσ = (N-differing_spins_σ) - differing_spins_σ
    Cσ /= N
    Cτ = (N-differing_spins_τ) - differing_spins_τ
    Cτ /= N
    return (Cσ+Cτ)/2
end

function gettimes()

    function rectify_times!(array)
        L = length(array)
        for i in 2:L
            if array[i] == array[i-1]
                for j in i:L
                    array[j] += 1
                end
            end
        end
    end

    const tw1list = map(Int64,[1e4,1e6,1e8])
    tws = logspace(log10(1),log10(3e8+10),100)
    tws = map(x->floor(Int64,x), tws)
    tws = vcat(tws, tw1list) |> sort
    rectify_times!(tws)

    tws = unique(tws) |> sort
    ts = logspace(4,8,10)
    ts = map(x->floor(Int64,x), ts) |> sort

    return tws, ts
end

function measureandprint(confs, Jhash)
    confs = sort(confs)

    twlist, tlist = gettimes()
    ℓ   = length(confs)
    bconfs = map(conf2data, confs)

    L = bconfs[1]["L"]
    Ts = [1/c["β"] for c in bconfs]
    if Ts |> unique |> length != 1
        die("Multiple temperatures in quench!")
    end
    T = Ts[1]

    timeslist = [c["t"] for c in bconfs]
    realmag_σ_list = [c["mag_sigma"] for c in bconfs]
    realmag_τ_list = [c["mag_tau"] for c in bconfs]
    E_σ_list = [c["energy_sigma"] for c in bconfs]
    E_τ_list = [c["energy_tau"] for c in bconfs]

    E_list = 0.5*(E_σ_list + E_τ_list)
    M_list = 0.5*(realmag_σ_list + realmag_τ_list)

    # Consistency check
    measmag_σ_list = [(2sum(c["spins_σ"])-L^3)/L^3 for c in bconfs]
    measmag_τ_list = [(2sum(c["spins_τ"])-L^3)/L^3 for c in bconfs]
    if measmag_σ_list ≉ realmag_σ_list
        die("Readed magnetization is not ≈ calculated magnetization (σ spins)")
    end
    if measmag_τ_list ≉ realmag_τ_list
        die("Readed magnetization is not ≈ calculated magnetization (τ spins)")
    end

    for i in 1:ℓ
        tw = bconfs[i]["t"]
        tw ∉ twlist && continue
        for j in i+1:ℓ
            tw_plus_t = bconfs[j]["t"]
            t = tw_plus_t - tw
            t ∉ tlist && continue
            # Wohoo! measure!
            T = 1/bconfs[i]["β"]
            C = getC(bconfs[i],bconfs[j])
            println("$L,$tw,$T,\"$Jhash\",\"Energy\",0,$(E_list[i])")
            println("$L,$tw,$T,\"$Jhash\",\"Magnetization\",0,$(M_list[i])")
            println("$L,$tw,$T,\"$Jhash\",\"Correlation\",$t,$C")
        end
    end

end

println("\"L\",\"mc\",\"T\",\"Jhash\",\"observable\",\"parameter\",\"value\"")
measureandprint(confs, Jhash)
