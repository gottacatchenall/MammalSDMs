using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using JSON
using CairoMakie
using GeoMakie

species = readdir("RodentSDMs")

bigtheme = Theme(
    CairoMakie = (
        px_per_unit = 3,
    )
)
set_theme!(bigtheme)


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
    @info dir
    Dict(
        :prediction => SDMLayer(joinpath(dir,"prediction.tif")),
        :uncertainty => SDMLayer(joinpath(dir, "uncertainty.tif")),
        :future => [SDMLayer(joinpath(dir, "future", x*".tif")) for x in YEARS],
        :metrics => JSON.parse(readline(joinpath(dir, "metrics.json")))
    )
end

function load_weights(focal_species; type = "exponential")
    df = DataFrame(CSV.File(joinpath("data", "$(type)_weights.csv")))
    row = filter(x->x.species == focal_species, df)
    dict = Dict([n=>row[!, n][begin] for n in names(row) if n != "species"])
end


sdms = Dict([s=>read_sdms(joinpath("RodentSDMs", s)) for s in species])


# Peromyscus leucopus

weight_dict = load_weights("Peromyscus leucopus", type = "exponential")
weight_dict["Peromyscus leucopus"]
dose = sum([weight_dict[s] .* quantize(v[:prediction]) for (s,v) in sdms])
unc = sum([weight_dict[s] .* quantize(v[:uncertainty]) for (s,v) in sdms])
dose, unc = quantize(dose), quantize(unc)



YEARS = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]
futures = [sum([weight_dict[s] .* quantize(v[:future][i]) for (s,v) in sdms]) for i in 1:4]

# Zero out oceans
for f in futures
    f.indices .= 0
    f.indices .= dose.indices[1:size(dose,1),1:size(dose,2)-1]
end 
futures = [quantize(f) for f in futures]





# =================================================================
# Bivariate
# =================================================================

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

# =================================================================
# Dose Shift 
# =================================================================
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
save("shift.png", f)




# =================================================================
# Focal
# =================================================================

function get_layer_timeseries(sdms, species, thresholded = true)
    prediction = sdms[species][:prediction]
    futures = sdms[species][:future]
    # Zero out oceans
    for f in futures
        f.indices .= 0
        f.indices .= dose.indices[1:size(prediction,1),1:size(prediction,2)-1]
    end 

    layers = [sdms[species][:prediction], sdms[species][:future]...]
    if thresholded
        for l in layers
            l.grid = l.grid .> sdms[species][:metrics]["threshold"]
        end
    end

    return layers  
end

begin
titles = ["Baseline", YEARS...]

f = Figure(size=(1500,500))
g = GridLayout(f[1,1])

species = "Peromyscus leucopus"
layers = get_layer_timeseries(sdms, species)
ax = Axis(g[1,1])
limits!(ax, -0.1, 0.7, -0.02, 0.04)
text!(ax, 0, 0, text=replace(species, " " => "\n"), fontsize=25)
hidedecorations!(ax)
hidespines!(ax)
for i in eachindex(layers)
    #ax = Axis(g[1,i+1], aspect=DataAspect(), title=titles[i])
    #heatmap!(ax, layers[i])
    ax = GeoAxis(g[1,i+1], aspect=DataAspect(),  dest = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96", title=titles[i], titlesize=22,)
    surface!(ax, layers[i], shading=false, colormap=[:grey92, :seagreen4])
    hidedecorations!(ax)
end 

species = "Peromyscus maniculatus"
layers = get_layer_timeseries(sdms, species)
ax = Axis(g[2,1])
limits!(ax, -0.1, 0.7, -0.02, 0.04)
text!(ax, 0, 0, text=replace(species, " " => "\n"), fontsize=25)
hidedecorations!(ax)
hidespines!(ax)
for i in eachindex(layers)
    #ax = Axis(g[2,i+1], aspect=DataAspect())
    #heatmap!(ax, layers[i])
    ax = GeoAxis(g[2,i+1], aspect=DataAspect(),  dest = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96")
    surface!(ax, layers[i], shading=false, colormap=[:grey92, :seagreen4])
    hidedecorations!(ax)
end 

#ax = GeoAxis(g[1,1], aspect=DataAspect(),  dest = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96", titlealign=:left)
#hidespines!(ax)
#surface!(ax, quantize(sdms[species][:uncertainty]), shading=false, colormap=:Greens)
f
end
save("leucopus_and_manuculatus.png", f)



# 1) We would like 3 sets of multi-panel figures. 

# Figure 1: Weighted species distribution density with increasing dose/increasing uncertainty heatmap. 

# Figure 2: Parallel figure but instead, we will provide prediction maps under one climate change scenario without detailed predictions but with all species included. 

# Figure 3: Detailed climate predictions for just Peromyscus leucopus and Peromyscus maniculatus

# Figure 4: 3 different scaling functions that take phylogenetic distance and then output weight?


# ================================================================
# Alternative weighting functions
# ================================================================

weighting_functions = ["linear", "exponential", "gaussian"]
doses = []
uncs = []

for w in weighting_functions
    weight_dict = load_weights("Peromyscus leucopus", type = w)
    dose = sum([weight_dict[s] .* quantize(v[:prediction]) for (s,v) in sdms])
    unc = sum([weight_dict[s] .* quantize(v[:uncertainty]) for (s,v) in sdms])
    dose, unc = quantize(dose), quantize(unc)
    push!(doses, dose)
    push!(uncs, unc)
end


numbins = 3
cscheme = ArcMapOrangeBlue()

cidxs = CartesianIndices((1:2, 1:2))

bv = bivariate!(ax, doses[1], uncs[1] ; cscheme..., xbins = numbins, ybins = numbins)

function get_bivar_layer(bv)
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
    return newpal, pal_position
end 

begin 
f = Figure(size=(1500, 1500))
g = GridLayout(f[1,1])
for i in 1:3
    ci = cidxs[i]
    #ax = Axis(f[ci[1], ci[2]], title=weighting_functions[i])
    f2 = Figure()
    ax2 = Axis(f2[1,1])
    bv = bivariate!(ax2, doses[i], uncs[i] ; cscheme..., xbins = numbins, ybins = numbins)
    
    newpal, pal_position = get_bivar_layer(bv)

    ax = GeoAxis(g[ci[1], ci[2]]; title = weighting_functions[i], titlesize = 20, aspect = DataAspect(), dest = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96")
    bv = surface!(ax, pal_position; shading=false, colormap=vec(newpal))
    hidedecorations!(ax)
end
le = Axis(g[2,2]; aspect = 1, xlabelsize=22, ylabelsize=22, height=200, width=200, xlabel = "Dose", ylabel = "Uncertainty", xticklabelsvisible=false, yticklabelsvisible=false, xticksvisible=false, yticksvisible=false,)
heatmap!(le, newpal)
rowsize!(g, 1, Relative(0.5))
rowsize!(g, 2, Relative(0.5))
colsize!(g, 1, Relative(0.5))
colsize!(g, 2, Relative(0.5))
f
end 



save("distance_function.png", f)

