
using CSV, DataFrames
using DelimitedFiles

function load_distance_matrix()
    df = DataFrame(CSV.File(joinpath("data", "distance_matrix.csv")))
    rename!(df, :Column1 => :species)
    filter!(x->length(split(x.species, "_")) >= 2, df)

    cols_to_keep = [name for name in names(df) if length(split(name, "_")) >= 2 || name == "species"]
    findall(x->contains(x, "leucopus"), cols_to_keep)
    select!(df, cols_to_keep)

    # Takes binomial name with space and fixes it if it has issues
    function adjust_names(string)
        replace(
            string,
            "Casiomys melanotis" => "Handleyomys melanotis",
            "Casiomys rostratus" => "Handleyomys rostratus",
            "Casiomys alfaroi" => "Handleyomys alfaroi",
            "Casiomys rhabdops" => "Handleyomys rhabdops",
            "Casiomys chapmani" => "Handleyomys chapmani",
            "Casiomys saturatior" => "Handleyomys saturatior",
            "Casiomys guerrerensis" => "Handleyomys guerrerensis",
            "Uroticellus" => "Urocitellus",
        )
    end

    # Create a dictionary to map base name to the fixed name

    rename_map = Dict([df.species[i] => adjust_names(x[1] * " " * x[2]) for (i,x) in enumerate(split.(df.species, "_"))])

    df.species = [rename_map[x] for x in df.species]
    for c in names(df)
        if c != "species" 
            rename!(df, c => rename_map[c]) 
        end 
    end 
    return df
end 

df = load_distance_matrix()
species_w_enough_occurrences = readlines(joinpath("data", "species_with_enough_occurrences_for_sdms.txt"))
distance_matrix = filter(x -> x.species ∈ species_w_enough_occurrences, df)
select!(distance_matrix, ["species", species_w_enough_occurrences...])
sort!(distance_matrix, :species)


function create_weight_df(; dist = :linear)


    score_df = DataFrame([x=>[] for x in names(distance_matrix)]...)
    for r in eachrow(distance_matrix)
        row = [x for x in r[2:end]]


        dist_func = Dict(
            :linear => x -> maximum(row) - x,
            :gaussian => x -> exp(-x^2/(2σ^2)),
            :exponential => x -> exp(-x)
        )[dist]

        score_row = dist_func.(row)
        s = sum(score_row)
        normalized_score = score_row ./ s
        push!(score_df, [r.species, normalized_score...])
    end

    score_df
end 


CSV.write(joinpath("data", "linear_weights.csv"), create_weight_df(;dist=:linear))

CSV.write(joinpath("data", "exponential_weights.csv"), create_weight_df(; dist = :exponential))

σ = 1
CSV.write(joinpath("data", "gaussian_weights.csv"), create_weight_df(; dist = :gaussian))


