using OffsetArrays;

function translate(vec, translation_steps, N)
    translation_steps %= N
    head_mask = ~(~0 << translation_steps) # generate the mask that is all zeroes except the last translation_steps bits
    head = vec & head_mask # get the last the last translation_steps of the state vector
    tail = vec >> translation_steps # get the first (n_conf-translation_steps) bits of the vector
    final = (head << (N-translation_steps)) | tail # combine accordingly
    return final
end

function find_class(translated_vec, vec)
    res = findfirst(j -> j == translated_vec, 1:vec-1)
    return res === nothing ? -1 : res
end

function add_new_class(vec, classes, num_vectors_class)
    classes[vec] = vec
    num_vectors_class[vec] = 1
end

function set_class(vec, class, classes, num_vectors_class)
    classes[vec] = classes[class]
    num_vectors_class[classes[class]] += 1
end

function classify(n_spins, tr_steps)    
    l = n_spins ÷ tr_steps
    n_conf = 2^n_spins

    classes = OffsetVector(zeros(UInt, n_conf+1), 0:n_conf)
    num_vecs_class = OffsetVector(zeros(UInt, n_conf+1), 0:n_conf)
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

function write_to_file(classes, num_vecs_class, filepath)
    f = open(filepath, "w")
    for i in 0:length(classes)-1
        if num_vecs_class[i] > 0
            write(f, string(num_vecs_class[i]) * "\t" * string(classes[i])*"\n")
        end
    end
    close(f)
end

function main()
    n_t = 28
    tr_t = 2
    classes_t, num_vecs_class_t = classify(n_t, tr_t)
    write_to_file(classes_t, num_vecs_class_t, "symblock_julia_n"*string(n_t)*"_tr"*string(tr_t)*".dat")
end

main()