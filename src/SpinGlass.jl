# Utils for both simulation and measuring of nets.
module SpinGlass

using DataFrames

export Net, metropolis!, magnetization, energy, χsg, measure

#
#
# Constructor methods
#
#

type Net
    Jup     :: Array{Int8,3}
    Jright  :: Array{Int8,3}
    Jfront  :: Array{Int8,3}
    σ       :: Array{Int8,3}
    L       :: Int64
end

function Net(L::Int64; kind=:spinglass)
    if kind == :spinglass
    Net(  rand([-1,1],L,L,L) # Jup
        , rand([-1,1],L,L,L) # Jright
        , rand([-1,1],L,L,L) # Jfront
        , rand([-1,1],L,L,L) # σ
        , L)                 # L
    elseif kind == :ferromagnet
    Net(  ones(Int8,L,L,L)   # Jup
        , ones(Int8,L,L,L)   # Jright
        , ones(Int8,L,L,L)   # Jfront
        , rand([-1,1],L,L,L) # σ
        , L)                 # L
    elseif kind == :antiferromagnet
    Net(  -ones(Int8,L,L,L)  # Jup
        , -ones(Int8,L,L,L)  # Jright
        , -ones(Int8,L,L,L)  # Jfront
        , rand([-1,1],L,L,L) # σ
        , L)                 # L
    else
        error("Unsuported initializer kind \"$kind\"")
    end
end


#
#
# Evolution
#
#

# Contribution to the hamiltonian of a spin.
function H(sample,x,y,z)
    L = sample.L
    left   = sample.Jright[mod1(x-1,L),y,z] * sample.σ[mod1(x-1,L), y, z]
    right  = sample.Jright[x,y,z]           * sample.σ[mod1(x+1,L), y, z]
    up     =    sample.Jup[x,y,z]           * sample.σ[x, mod1(y+1,L), z]
    down   =    sample.Jup[x,mod1(y-1,L),z] * sample.σ[x, mod1(y-1,L), z]
    front  = sample.Jfront[x,y,z]           * sample.σ[x, y, mod1(z+1,L)]
    behind = sample.Jfront[x,y,mod1(z-1,L)] * sample.σ[x, y, mod1(z-1,L)]

    neighbours_term = left + right + up + down + behind + front

    return -sample.σ[x,y,z] * neighbours_term
end

function metropolis!(sample,β)
    L = sample.L
    for z in 1:L, y in 1:L, x in 1:L
        # Get ΔE and probability of the tentative flip.
        ΔE = -2*H(sample,x,y,z)
        P = exp(-β*ΔE)
        if P > rand()
            sample.σ[x,y,z] *= -1
        end
    end
end


#
#
# Observables
#
#

magnetization(sample::Net)     = mean(sample.σ)
magnetization(spins::BitArray) = 2*(mean(spins)-0.5)

function energy(sample)
    L = sample.L
    acc = 0
    for x in 1:L, y in 1:L, z in 1:L
        acc += H(sample,x,y,z)
    end
    return acc / (L*L*L) / 3 / 2
end

function χsg(sample)
    L = sample.L
    SiSj = 0
    Si = 0
    for x in 1:L, y in 1:L, z in 1:L
        right  = sample.σ[mod1(x+1,L), y, z]
        up     = sample.σ[x, mod1(y+1,L), z]
        front  = sample.σ[x, y, mod1(z+1,L)]
        SiSj  += sample.σ[x,y,z] * (right + up + front)
        Si    += sample.σ[x,y,z]
    end
    SiSj /= 3*L*L*L
    Si   /= L*L*L
    return (SiSj-Si*Si)^2
end

#
#
# Measuring utilities.
#
#

# +-+-+--+- → BitArray
line_to_spins(line) = [x=='+'?true:false for x in line] |> BitArray
# Bitarray to spin list
function bits_to_spins(b,L)
    # From bits to +1 and -1
    ints = [x ? Int8(+1) : Int8(-1) for x in b]
    # Reshape, return
    return reshape(ints,(L,L,L))
end

# Get fields from a file created with annealer.c
function get_data(file)
    data = readlines(file)

    # Get the global
    Lline = split(data[1])
    L = parse(Int64,Lline[2])

    Jup_line = split(data[9])
    Jup = line_to_spins(Jup_line[2])

    Jright_line = split(data[10])
    Jright = line_to_spins(Jright_line[2])

    Jfront_line = split(data[11])
    Jfront = line_to_spins(Jfront_line[2])

    # Iterate over the rest of the lines, build a list of spins, β and
    # mc step.
    mcs   = Int64[]
    βs    = Float64[]
    spins = BitArray[]

    for line in data[12:end]
        #if is a comment, do nothing.
        if line[1] != '#'
            # Read mc, β, spin list
            fields = split(line,['\t','\n'])
            push!(mcs   ,parse(Int64  ,fields[1]) )
            push!(βs    ,parse(Float64,fields[2]) )
            push!(spins ,line_to_spins(fields[3]) )
        end
    end

    # Return everything
    return Jup,Jright,Jfront,mcs,βs,spins
end

# Parse the fields of a file created with annealer.c in a DataFrame,
# measuring the observables.
# Observables is a string that can contain E,C,χ,M.
function measure(file;T0=10,observables="")
    Jup,Jright,Jfront,mcs,βs,spins = get_data(file)
    Jhash = hash([Jup Jright Jfront])
    N = length(spins[1])
    L = floor(Int64,cbrt(N))
    Nβ = 1 + (βs |> diff |> countnz) # Number of β changes + 1
    # Build a tentative net with the J's. Spins will be changed on the
    # fly.
    net = Net( bits_to_spins(Jup,L)
               ,bits_to_spins(Jright,L)
               ,bits_to_spins(Jfront,L)
               ,zeros(Int8,L,L,L)
               ,L    )
    # Fill this.
    df = DataFrame()
    for i in 1:length(spins)
        σ   = spins[i]
        # Fill the net spins
        S = bits_to_spins(σ,L)
        net.σ = S
        # Compute observables if asked for
        M       = NaN
        E       = NaN
        χsglass = NaN
        if 'M' ∈ observables
            M       = magnetization(σ)
        end
        if 'E' ∈ observables
            E       = energy(net)
        end
        if 'χ' ∈ observables || 'X' ∈ observables
            χsglass = χsg(net)
        end
        # For the correlation, we need to compare spins[i] and
        # spins[i+T0]. But i is the measurement number, and we need
        # the measurement number at that temperature, nothing more.
        # That's easy, because the number of measurements per β is
        # just length(spins)/Nβ. Use a modulo operator.
        C = NaN
        if 'C' ∈ observables
            Nmeas_per_β = length(spins)/Nβ
            if mod1(i,Nmeas_per_β) > Nmeas_per_β-T0
                C = NaN
            else
                # Correlation is just the scalar product. Because they are
                # bit arrays with elements ∈(0,1), they must translated to
                # spins in ∈(-1,1).
                thisσ = S
                nextσ = bits_to_spins(spins[i+T0], L)
                # Now, they are cubes of spins. We can't use dot product
                # on that kind of tensor; use a workaround:
                C = map(*,thisσ,nextσ) |> mean
            end
        end
        μdf = DataFrame( mc  = mcs[i]
                        ,β   = βs[i]
                        ,J   = Jhash
                        ,M   = M
                        ,E   = E
                        ,χsg = χsglass
                        ,C   = C
                        ,C2  = βs[i]*(1-C)  )
        rename!(μdf,:C2,Symbol("β(1-C)"))
        df = vcat(df,μdf)
    end
    return df
end


end#module
