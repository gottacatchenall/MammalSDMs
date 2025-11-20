using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using EvoTrees
using JSON
using CairoMakie
const SDT = SpeciesDistributionToolkit

function load_baseline_climate_layers(worldclim_dir)
    baseline_directory = joinpath(worldclim_dir)

    filenames = readdir(baseline_directory)
    
    layer_paths = [
        joinpath(baseline_directory, filenames[findfirst(fn -> occursin(pattern, fn), filenames)]) 
        for pattern in ["_bio_$i.tif" for i in 1:19]
    ]
    
    # Load and convert to Float32
    layers = [SDMLayer(path; bottom=-60.) for path in layer_paths]
    return [Float32.(layer) for layer in layers]
end

function generate_pseudoabsences(
    presence_layer, 
    buffer_distance_km, 
    class_balance_ratio
)
    background = pseudoabsencemask(DistanceToEvent, presence_layer)
    buffer = pseudoabsencemask(WithinRadius, presence_layer; distance = buffer_distance_km)
    
    nodata_mask = nodata(buffer, true) 
    background.indices .= nodata_mask.indices
    
    num_pseudoabsences = Int(round(class_balance_ratio * sum(presence_layer)))
    return backgroundpoints(background, num_pseudoabsences)
end

function prepare_training_data(layers, presence_layer, absence_layer)
    X = Matrix(hcat([
        vcat(layer[findall(presence_layer)], layer[findall(absence_layer)]) 
        for layer in layers
    ]...)')
    
    y = Bool.(vcat(
        [1 for _ in findall(presence_layer)], 
        [0 for _ in findall(absence_layer)]
    ))
    
    return X, y
end

function train_model(X_train, y_train; max_depth = 6)
    return EvoTrees.fit(
        EvoTreeGaussian(max_depth = max_depth),
        x_train = X_train',
        y_train = y_train,
    )
end

function predict_distribution(model, feature_matrix)
    return EvoTrees.predict(model, feature_matrix')
end

function create_prediction_layer(model, environmental_layers)
    prediction_layer = deepcopy(environmental_layers[begin])
    uncertainty_layer = deepcopy(environmental_layers[begin])
    
    feature_matrix = Matrix(hcat([
        [layer[i] for layer in environmental_layers] 
        for i in eachindex(environmental_layers[1])
    ]...))
    
    predictions = predict_distribution(model, feature_matrix)
    
    prediction_layer.grid[findall(prediction_layer.indices)] .= predictions[:, 1]
    uncertainty_layer.grid[findall(prediction_layer.indices)] .= predictions[:, 2]
    
    return prediction_layer, uncertainty_layer
end

function calculate_evaluation_metrics(y_true, y_predicted, thresholds=0:0.001:1)
    # Calculate confusion matrices across all thresholds
    confusion_matrices = [ConfusionMatrix(y_predicted .> t, y_true) for t in thresholds]
    false_positive_rates, true_positive_rates = fpr.(confusion_matrices), tpr.(confusion_matrices)
    
    # Calculate ROC-AUC using trapezoidal rule
    roc_dx = [reverse(false_positive_rates)[i] - reverse(false_positive_rates)[i-1] for i in 2:length(false_positive_rates)]
    roc_dy = [reverse(true_positive_rates)[i] + reverse(true_positive_rates)[i-1] for i in 2:length(true_positive_rates)]
    roc_auc = sum(roc_dx .* (roc_dy ./ 2.0))
    
    # Calculate PR-AUC using trapezoidal rule
    precisions = ppv.(confusion_matrices)
    pr_dx = [reverse(true_positive_rates)[i] - reverse(true_positive_rates)[i-1] for i in 2:length(true_positive_rates)]
    pr_dy = [reverse(precisions)[i] + reverse(precisions)[i-1] for i in 2:length(precisions)]
    pr_auc = sum(pr_dx .* (pr_dy ./ 2.0))
    
    # Find optimal threshold using MCC
    _, threshold_index = findmax(mcc.(confusion_matrices))
    optimal_threshold = thresholds[threshold_index]

    return Dict(
        :prauc => pr_auc,
        :rocauc => roc_auc,
        :tss => trueskill(confusion_matrices[threshold_index]),
        :mcc => mcc(confusion_matrices[threshold_index]),
        :threshold => optimal_threshold
    )
end

function aggregate_fold_statistics(fold_stats)
    return Dict(
        metric => Dict("mean" => mean(values), "std" => std(values))
        for metric in keys(first(fold_stats)) 
        for values in [[fold[metric] for fold in fold_stats]]
    )
end


function fit_sdm(
    occurrences,
    environmental_layers;
    pseudoabsence_buffer_distance = 25.0,
    class_balance = 1.0,
    max_depth = 6,
    k = 4
)
    # Create presence layer from occurrence points
    presence_layer = mask(environmental_layers[begin], occurrences)
    
    @info "    |-> Generating pseudoabsences..."
    Random.seed!(123) # standardize across tuning tests
    absence_layer = generate_pseudoabsences(presence_layer, pseudoabsence_buffer_distance, class_balance)
    
    # Prepare training data
    features, labels = prepare_training_data(environmental_layers, presence_layer, absence_layer)
    
    # Create cross-validation folds
    Random.seed!(123) # standardize across tuning tests
    fold_indices = SDeMo.kfold(labels, features, k=k)
    
    # Storage for results
    true_labels = Bool[]
    out_of_fold_predictions = Float32[]
    
    @info "    |-> Training $(k) cross-validation folds..."
    
    # Train and evaluate each fold
    for (train_idx, validation_idx) in fold_indices
        model = train_model(features[:, train_idx], labels[train_idx]; max_depth = max_depth)
        
        # Evaluate on validation set
        validation_predictions = predict_distribution(model, features[:, validation_idx])[:, 1]

        true_labels = vcat(true_labels, labels[validation_idx])
        out_of_fold_predictions = vcat(out_of_fold_predictions, validation_predictions)
    end
    
    # Compute fit stats on out-of-fold predictions
    fit_stats = calculate_evaluation_metrics(true_labels, out_of_fold_predictions)
    optimal_threshold = fit_stats[:threshold]
    
    # Fit full model
    model = train_model(features, labels; max_depth = max_depth)
    prediction, uncertainty = create_prediction_layer(model, environmental_layers)    
    range_map = Int.(prediction .> optimal_threshold)
    
    return model, range_map, uncertainty, fit_stats, presence_layer, absence_layer
end


function write_sdm_artifacts(artifact_dir, species_name, results)
    mkpath(artifact_dir)
    output_dir = joinpath(artifact_dir, species_name)
    mkpath(output_dir)
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "presences.tif"),
        Float32.(results[:presences])
    )      
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "absences.tif"),
        Float32.(results[:absences])
    )      
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "uncertainty.tif"),
        results[:uncertainty]
    )        
    SDT.SimpleSDMLayers.save(
        joinpath(output_dir, "prediction.tif"),
        results[:prediction]
    )        

    open(joinpath(output_dir, "metrics.json"), "w") do f
        JSON.print(f, results[:metrics])
    end
end

function load_occurrences(path)
    df = CSV.read(path, DataFrame)
    species = split(split(path, "/")[end], ".")[begin]
    return Occurrences([Occurrence(what=species, presence=true, where=(row.decimalLongitude, row.decimalLatitude)) for row in eachrow(df)])
end


job_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
data_dir = "/project/def-tpoisot/mcatchen/MammalOccurrence"
artifact_dir = "/scratch/mcatchen/MammalSDMs"
worldclim_dir = "/project/def-tpoisot/mcatchen/JuliaEnvironments/ColoradoBees/SimpleSDMDatasets/WorldClim2/BioClim/"

species_paths = joinpath.(data_dir, sort(readdir(data_dir)))

@info "Loading occurrences from: $(species_paths[job_id])"
occs = load_occurrences(species_paths[job_id])

species_name = String(split(split(species_paths[job_id], "/")[end], ".")[1])

@info "Fitting model for species: $species_name"

environmental_layers = load_baseline_climate_layers(worldclim_dir)

model, range_map, uncertainty_map, statistics, presences, absences = fit_sdm(
    occs,
    environmental_layers
)

results = Dict(
    :baseline => Dict(
        :range => range_map,
        :uncertainty => uncertainty_map,
        :presences => presences,
        :absences => absences,
    ),
    :future => future_projections,
)

write_sdm_artifacts(artifact_dir, species_name, results)

