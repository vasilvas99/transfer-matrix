using Pkg;
Pkg.add("Arpack");
Pkg.add("CSV");
Pkg.add("DataFrames");
Pkg.add("Plots");
Pkg.add("LaTeXStrings");
Pkg.add("ParallelDataTransfer");
using CSV;
@everywhere using DataFrames;
@everywhere using LinearAlgebra;
@everywhere using Arpack;
using Plots;
using LaTeXStrings;
@everywhere using ParallelDataTransfer, Distributed


@everywhere function spin_state(num_spins, chain_configuration, spin_number)
    if spin_number > num_spins
         throw(DomainError(spin_number, "Spin number should be <= num_spins"))
    end
    spin_number -= 1 # start counting of bits from 1
    mask = 1 << spin_number
    state = chain_configuration & mask
    return state != 0 ? 1 : -1
end

@everywhere function inter_chain_contribution(num_spins, ch1, ch2)
    hij = 0
    for i in 1:num_spins
        hij += spin_state(num_spins, ch1, i)*spin_state(num_spins, ch2, i)
    end
    return hij
end

@everywhere function chain_contribution(num_spins, chain_conf)
    hij = 0
    for i in 1:num_spins
        hij += (spin_state(num_spins, chain_conf, i)
        *spin_state(num_spins, chain_conf, i%num_spins + 1))
    end
    return hij
end

@everywhere function hamiltonian(n, ch1, ch2, jnn)
    hij = (
        0.5*jnn*chain_contribution(n, ch1)
        + 0.5*jnn*chain_contribution(n, ch2)
        + jnn*inter_chain_contribution(n, ch1, ch2)
    )
    return hij
end

@everywhere function tm_element(n, ch1, ch2, jnn, temp)
    return exp(-hamiltonian(n, ch1, ch2, jnn)/temp)
end

@everywhere function tm_n(n, temp)
    n_states = 2^n
    jnn = -1
    
    # fastest way to calculate a matrix in Julia is with a list-comprehension
    tm_block = [tm_element(n, ch1, ch2, jnn, temp) for ch1 = 1:n_states, ch2 = 1:n_states]
    
    return tm_block
end

@everywhere function corr_len(n_spins, temp)
    tm_block = tm_n(n_spins, temp)
    if n_spins > 2
        # Arnoldi iteration can't find the smallest eigenvalue (only the n_states-1 largest),
        # thus if n = 1, and n_states = 2 => the second eigenvalue won't be found
        # ref: https://en.wikipedia.org/wiki/Arnoldi_iteration#Finding_eigenvalues_with_the_Arnoldi_iteration
        eig_val, eig_vec = eigs(tm_block, nev=2); # first and second largest eigenvalues
        cln = 1.0/(log(eig_val[1]/abs(eig_val[2])))
    else
        # the built-in Julia method finds ALL eigenvalues (not just the k largest), 
        # but is slower than ARPACK (ok for 2x2 matrices)
        eig_val = sort(eigvals(tm_block), rev=true)
        cln = 1.0/(log(eig_val[1]/abs(eig_val[2])))
    end
    
    return cln
end

@everywhere function calculate_ratios(n_spins, temp)
    r = [(corr_len(i, temp)*(i+1))/(corr_len(i+1, temp)*i) for i = 2:n_spins-1]
    return r
end

@everywhere function calculate_df_row(n_spins, temp0, tempinc, i)
    t = temp0+tempinc*i
    ratios = calculate_ratios(n_spins, t)
    row = [n_spins, t]
    append!(row, ratios)
    return row
end

function calculate_temp_range(; n_spins, temp0, tempinc, max_iter)
    column_names =  ["r$(i)$(i+1)" for i = 2:n_spins-1]
    prepend!(column_names, ["MaxSpins", "temp"])
    rows = DataFrame([Float64[] for i in 1:length(column_names)], column_names)
    
    r = pmap(i -> calculate_df_row(n_spins, temp0, tempinc, i), 0:max_iter+1)
    
    for row in r
        push!(rows, row)
    end
    
    return rows
end

# Parameters for final calculation
MAX_SPINS_CHAIN = 15;
TEMP0 = 2.0;
TEMP_INC = 0.025;
MAX_ITER = 40;

try
    mkdir("corr-len-output")
    println("Created directory corr-len-output")
catch e
    println("Directory corr-len-output already exists...Proceding")
end

@time temp_range = calculate_temp_range(n_spins=MAX_SPINS_CHAIN, temp0=TEMP0, tempinc=TEMP_INC, max_iter=MAX_ITER)
CSV.write("corr-len-output/corl-tmp-long_N$(MAX_SPINS_CHAIN)_T0$(TEMP0)_INC$(TEMP_INC).dat", temp_range, delim="\t")

p = plot(dpi=600)
xlabel!(p, "T")
ylabel!(p, "correlation length ratio")
for i in 2:Int(temp_range.MaxSpins[1]-1)
    plot!(p, temp_range.temp, temp_range[!,"r$(i)$(i+1)"], label=latexstring("\$r_{$(i), $(i+1)}\$"))
end

savefig(p, "corr-len-output/corl-tmp-long_N$(MAX_SPINS_CHAIN)_T0$(TEMP0)_INC$(TEMP_INC).png") 