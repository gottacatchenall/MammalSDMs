using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using JSON
using CairoMakie
using GeoMakie

species = readdir("RodentSDMs")


function _multiply(c1::Makie.ColorTypes.RGBA, c2::Makie.ColorTypes.RGBA)
    return Makie.ColorTypes.RGBA(
        (c1.r * c2.r),
        (c1.g * c2.g),
        (c1.b * c2.b),
        (c1.alpha + c2.alpha) / 2,
    )
end
function _bivariate_grid(xbins, ybins, xcm, ycm)
    cx = Makie.cgrad(xcm, xbins; categorical = true)
    cy = Makie.cgrad(ycm, ybins; categorical = true)
    return [_multiply.(x, y) for x in cx, y in cy]
end




function read_sdms(dir)
    YEARS = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]
    Dict(
        :prediction => SDMLayer(joinpath(dir,"prediction.tif")),
        :uncertainty => SDMLayer(joinpath(dir, "uncertainty.tif")),
        :future => [SDMLayer(joinpath(dir, "future", x*".tif")) for x in YEARS],
        :metrics => JSON.parse(readline(joinpath(dir, "metrics.json")))
    )
end

function load_weights(focal_species)
    df = DataFrame(CSV.File(joinpath("data", "weights.csv")))

    @info focal_species
    row = filter(x->x.species == focal_species, df)
    Dict([n=>row[!, n][begin] for n in names(row) if n != "species"])

end


sdms = Dict([s=>read_sdms(joinpath("RodentSDMs", s)) for s in species])

weight_dict = load_weights("Peromyscus maniculatus")

weight_dict["Peromyscus maniculatus"]

dose = sum([weight_dict[s] .* quantize(v[:prediction]) for (s,v) in sdms])
unc = sum([weight_dict[s] .* quantize(v[:uncertainty]) for (s,v) in sdms])


dose, unc = quantize(dose), quantize(unc)


YEARS = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]


futures = [sum([weight_dict[s] .* quantize(v[:future][i]) for (s,v) in sdms]) for i in 1:4]

# Zero oceans
for f in futures
    f.indices .= 0
    f.indices .= dose.indices[1:size(dose,1),1:size(dose,2)-1]
end 
futures = [quantize(f) for f in futures]



numbins = 3
cscheme = ArcMapOrangeBlue()

f = Figure()
ax = Axis(f[1,1])
bv = bivariate!(ax, dose, unc ; cscheme..., xbins = numbins, ybins = numbins)


newpal = _bivariate_grid(
    bv.xbins[],
    bv.ybins[],
    bv.xcolormap[],
    bv.ycolormap[],
)

# Next, we turn each layer into a binarized version
xbin = discretize(bv.x[], bv.xbins[])
ybin = discretize(bv.y[], bv.ybins[])
idx = LinearIndices(newpal)

pal_position = similar(bv.x[], Int)
for k in keys(pal_position)
    pal_position[k] = idx[xbin[k], ybin[k]]
end

begin
f = Figure(size=(1200, 1000))
g = GridLayout(f[1,1])
ax = GeoAxis(g[1, 1]; aspect = DataAspect(), dest = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96")

gright = GridLayout(g[1,2])

grt = Axis(gright[1,1])
hidedecorations!(grt)
hidespines!(grt)

grb = Axis(gright[3,1])
hidedecorations!(grb)
hidespines!(grb)

le = Axis(gright[2,1]; aspect = 1, xlabelsize=22, ylabelsize=22, height=200, width=200, xlabel = "Dose", ylabel = "Uncertainty", xticklabelsvisible=false, yticklabelsvisible=false, xticksvisible=false, yticksvisible=false,)
hidedecorations!(ax)
bv = surface!(ax, pal_position; shading=false, colormap=vec(newpal))
heatmap!(le, newpal)

#=bivariatelegend!(le, dose, unc;
    cscheme...,
    xbins = numbins,
    ybins = numbins
)=#
colsize!(g, 1, Relative(0.8))
f
end
save("bivar.png", f)

#=
begin
f = Figure(size=(1200, 600))
g = GridLayout(f[1,1])
ax = Axis(g[1, 1]; aspect = DataAspect())
le = Axis(g[1, 2]; aspect = 1, xlabel = "Dose", ylabel = "Uncertainty")
hidespines!(ax)
bv = bivariate!(ax, dose, unc ; cscheme..., xbins = numbins, ybins = numbins)
bivariatelegend!(le, dose, unc;
    cscheme...,
    xbins = numbins,
    ybins = numbins
)
colsize!(g, 1, Relative(0.7))
f
end
=#


ax

function add_colorbar!(location, colorscheme; titlesize=12, titlealign=:left, halign=0.7, valign=0.16, width=0.2, height=0.03, title="")
    cbar_axis = Axis(
        location,
        width=Relative(width),
        height=Relative(height),
        halign=halign,
        valign=valign,
        xticksvisible=false,
        yaxisposition = :right,
        #yticksvisible=false,
        xticklabelsvisible=false,
        #yticklabelsvisible=false,
        yticks = ([1, 400], ["Loss", "Gain"]),
        yticklabelsize = 23,
        #title=title,
        #titlesize = titlesize,
        #titlefont=:regular,
        #titlealign=titlealign
    )
    cmap = [get(colorscheme, i) for i in 0:0.0025:1]
    X = collect(0:0.0025:1)
    heatmap!(cbar_axis, Matrix(X'), colormap=cmap)
end 


# Continuous diff
begin
f = Figure()
ax = Axis(f[1,1], title="Quantized Increase in Dose by 2100", titlealign=:left)
hidespines!(ax)
heatmap!(ax, futures[4] .- dose, colormap=:cork100, colorrange=(-0.7, 0.7))
add_colorbar!(f[1,1], Makie.ColorSchemes.cork100, title = "Dose Shift", width = 0.3, halign = 0.03)
#text!(ax, -143.5, 21, text="Loss", fontsize = 10)
#text!(ax, -118.5, 20.8, text="Gain", fontsize = 10)

f
end
save("shift.png", f)


begin
f = Figure(size=(1500,1500))
ax = GeoAxis(f[1,1], aspect=DataAspect(),  dest = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96", titlealign=:left)
#hidespines!(ax)
hidedecorations!(ax)
surface!(ax, futures[4] .- dose, shading=false, colormap=:cork100, colorrange=(-0.7, 0.7))
add_colorbar!(f[1,1], Makie.ColorSchemes.cork100, title = "Dose Shift", width = 0.025, height=0.3, valign=0.5, halign = 0.15)
text!(ax, -142.5, 30, text="Shift in Dose", rotation=0.5π, fontsize=28, )
#text!(ax, -143.5, 21, text="Loss", fontsize = 10)
#text!(ax, -118.5, 20.8, text="Gain", fontsize = 10)

f
end