using CSV, DataFrames
using DelimitedFiles

df = DataFrame(CSV.File(joinpath("data", "distance_matrix.csv")))
species_list = [x[1] * " " * x[2] for x in split.(df.Column1, "_") if length(x) >= 2]

# fix typo
species_list = [replace(
    x,
    "Uroticellus" => "Urocitellus",
) for x in species_list]

species_list = replace(
    species_list,
    "Casiomys melanotis" => "Handleyomys melanotis",
    "Casiomys rostratus" => "Handleyomys rostratus",
    "Casiomys alfaroi" => "Handleyomys alfaroi",
    "Casiomys rhabdops" => "Handleyomys rhabdops",
    "Casiomys chapmani" => "Handleyomys chapmani",
    "Casiomys saturatior" => "Handleyomys saturatior",
    "Casiomys guerrerensis" => "Handleyomys guerrerensis"
)


writedlm(joinpath("data", "species_list.txt"), species_list, '\n')

