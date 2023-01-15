# Plots the size of the largest symmetric class vs. the translation step for different chain sizes.
# Basic script usage:
# julia -p <number_of_cores> .\max_symblock_size.jl --start=<smallest_chain_size> --stop=<largest_chain_size> --step=<step_between_sizes>
# e.g. julia -p 3 .\max_symblock_size.jl --start=2 --stop=15 --step=3

if !@isdefined(PACKAGES_INSTALLED)
    using Pkg;
    Pkg.add("OffsetArrays");
    Pkg.add("Plots");
    Pkg.add("LaTeXStrings");
    Pkg.add("ArgParse");
    using ArgParse;
    using Plots;
    using LaTeXStrings;
    @everywhere using OffsetArrays;
    using Distributed
    PACKAGES_INSTALLED = true;
end

@everywhere function translate(vec, translation_steps, N)
    translation_steps %= N
    head_mask = ~(~0 << translation_steps) # generate the mask that is all zeroes except the last translation_steps bits
    head = vec & head_mask # get the last the last translation_steps of the state vector
    tail = vec >> translation_steps # get the first (n_conf-translation_steps) bits of the vector
    final = (head << (N-translation_steps)) | tail # combine accordingly
    return final
end

@everywhere function find_class(translated_vec, vec)
    res = findfirst(j -> j == translated_vec, 1:vec-1)
    return res === nothing ? -1 : res
end

@everywhere function add_new_class(vec, classes, num_vectors_class)
    classes[vec] = vec
    num_vectors_class[vec] = 1
end

@everywhere function set_class(vec, class, classes, num_vectors_class)
    classes[vec] = classes[class]
    num_vectors_class[classes[class]] += 1
end

@everywhere function classify(n_spins, tr_steps)    
    l = n_spins ÷ tr_steps
    n_conf = 2^n_spins

    classes = OffsetVector(zeros(Int64, n_conf+1), 0:n_conf)
    num_vecs_class = OffsetVector(zeros(Int64, n_conf+1), 0:n_conf)
    n_classes = 2 

    classes[0] = 0
    classes[1] = 1
    num_vecs_class[1] = 1
    num_vecs_class[0] = 1

    new_class_found = false
 
    for vec in 2:n_conf-2
        vec_temp = vec

        for _ in 1:l-1
            vec_temp = translate(vec_temp, tr_steps, n_spins)
            class = find_class(vec_temp, vec)
            if class != -1 # vector belongs to an already found class, update the class
                set_class(vec, class, classes, num_vecs_class)
                new_class_found = false
                break
            else # vector belongs to a new class, continue 
                new_class_found = true
            end
        end
        
        if new_class_found
            add_new_class(vec, classes, num_vecs_class)
            n_classes += 1 
        end

    end 
    
    classes[n_conf-1] = n_conf - 1
    num_vecs_class[n_conf-1] = 1
    n_classes += 1
    
    return (classes, num_vecs_class)
    
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--start"
            help = "Starting chain size"
            arg_type = Int
            default = 2
        "--stop"
            help = "Final chain size"
            arg_type = Int
            default = 21
        "--step"
            help = "Step between chain sizes"
            arg_type = Int
            default = 3
    end

    return parse_args(s)
end

@everywhere function count_nonzero(num_vecs_class)
    count = 0
    for i in 0:length(num_vecs_class)-1
        if num_vecs_class[i] > 0
            count += 1
        end
    end
    return count
end

function init_output_dir(dir_name)
    try
        mkdir(dir_name)
        println("Created directory $(dir_name)")
    catch e
        println("Directory $(dir_name) already exists...Proceeding")
    end
end

function main()
    cli = parse_commandline()
    START = cli["start"]
    STOP = cli["stop"]
    STEP = cli["step"]
    
    DIR_NAME = "classes-output"
    PLOT_PATH = "$(DIR_NAME)/trstep-vs-blocksize-START$(START)-STOP$(STOP)-STEP$(STEP).png"
    init_output_dir(DIR_NAME)

    p = plot(dpi=600)
    xlabel!(p, "Translation step")
    ylabel!(p, L"log_{2} (n_{classes})")
    println("Initialization done, starting calculations")
    for spins in range(START, stop=STOP, step=STEP)
        print("Running with $(spins) spins...")
        time_taken = @elapsed begin
            max_sizes = pmap((tr_step) -> log2(count_nonzero(classify(spins, tr_step)[2])), 1:spins)
            plot_legend=latexstring("\$n_{spins} = $(spins)\$")
            plot!(p, 1:spins, max_sizes, label=plot_legend)
        end
        println("Done ($(time_taken) s)")
    end
    println("Saving ...")
    savefig(p, PLOT_PATH) 
    print("Done ($(PLOT_PATH))")
end

main()


