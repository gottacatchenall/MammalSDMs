using DelimitedFiles
using CSV
using DataFrames
using SpeciesDistributionToolkit


function load_occurrences(path)
    df = CSV.read(path, DataFrame)
    species = split(split(path, "/")[end], ".")[begin]
    return Occurrences([Occurrence(what=species, presence=true, where=(row.longitude, row.latitude)) for row in eachrow(df)])
end


# Check how many hosts have enough occurrences to not be stupid 
species_to_keep = []
for x in readdir(joinpath("data", "species_occurrence"))
    occs = load_occurrences(joinpath("data", "species_occurrence", x)) 
    if length(occs) > 50
        push!(species_to_keep, split(x, ".")[1])
    end
end

writedlm(joinpath("data", "species_with_enough_occurrences.txt"), species_to_keep, '\n')

species_to_keep

readlines(joinpath("data", "species_with_enough_occurrences.txt"))
