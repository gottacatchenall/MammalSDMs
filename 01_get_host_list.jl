using CSV, DataFrames
using DelimitedFiles

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

writedlm(joinpath("data", "species_list.txt"), df.species, '\n')



