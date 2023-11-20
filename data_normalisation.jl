


# Function for normalising raw fragment counts of genomic features. The 
# function takes a count matrix (each column is a sample, and row is a 
# feature) and scales the samples so that their means are equal, while 
# simultaneously correcting for differences in their GC biases
using BlackBoxOptim, Statistics, StatsBase
function normalise_counts(counts, gc_frac)

    # Funcition for fitting a GC curve to logratios, and reducing the curve from the logratios
    function correct_gc_bias(gc_frac, logratio)

        # Polynomial to be fit to the data
        poly(x, coefs) = coefs[1] + coefs[2] * x + coefs[3] * x^2 + coefs[4] * x^3 + coefs[5] * x^4

        # Error function for the polynomial fit (total absolute error)
        valid = findall(isfinite.(logratio)) # features with count 0 get -Inf logarithmic values
        fit_to_lr = logratio[valid]
        fit_to_gc = gc_frac[valid]
        err_func(coefs) = sum([abs(fit_to_lr[k] - poly(fit_to_gc[k], coefs)) for k in 1:length(fit_to_lr)])

        # Fit the polynomial to the data
        res = bboptimize(err_func, SearchRange=(-10_000.0, 10_000.0), NumDimensions=5, 
            Method=:de_rand_1_bin, MaxSteps=50_000, TraceMode=:silent)

        # Reduce the polynomial from the logratios
        coefs = best_candidate(res)
        model = [poly(f, coefs) for f=gc_frac]
        corrected = logratio .- model

        return corrected
    end

    # Calculate the mean of all samples, against which logratios are calculated
    counts_mean = mean(counts, dims=2)[:]

    counts_corrected = zeros(size(counts)...)
    for s in 1:size(counts, 2)
        println("Sample $s/$(size(counts, 2))")
        # Calculate logratio against the mean of all samples
        logratio = log2.(counts[:,s] ./ counts_mean)

        # Remove any difference in GC bias between the sample and the mean of all samples
        logratio_corrected = correct_gc_bias(gc_frac, logratio)

        # Convert logratios back to normal count scale
        corrected = 2 .^ logratio_corrected .* counts_mean

        counts_corrected[:, s] .= corrected
    end

    return counts_corrected
end




# Function for normalising raw fragment counts of genomewide fragment counts in 
# 100bp windows. The function takes a count matrix (each column is a sample, 
# and row is a genomic region) and scales each sample so that the mean of each 
# sample is 1.0, while simultaneously correcting for GC biases.
using Statistics, BlackBoxOptim
function normalise_counts_with_control_peaks(counts, gc_frac; n_top_peaks=20_000, n_min_var_peaks=2_000)

    # Find top peaks based on the minimum of all samples. Each sample must have >0 count
    counts_min = minimum(counts, dims=2)[:];
    peaks = findall(counts_min .> 0)
    if length(peaks) < n_top_peaks
        if length(peaks) == 0; error("Did not find any positions with >0 count in >0 samples"); end
        println("WARNING: Found $(length(peaks)) positions with count >0 in >0 samples (which is less than given n_top_peaks: $n_top_peaks)")
    else
        peaks = peaks[sortperm(counts_min[peaks], rev=true)[1:n_top_peaks]]
    end
    
    # Select control peaks as having the smallest between-sample variance
    consensus_peak_counts = log2.(counts[peaks, :]);
    between_sample_var = var(consensus_peak_counts, dims=2)[:];
    control_i = sort(sortperm(between_sample_var)[1:n_min_var_peaks]);
    control_peaks = peaks[control_i];

    # Fit GC models and correct counts
    println("Fitting GC models")
    control_peak_gc = gc_frac[control_peaks]
    gc_models = map(1:size(counts, 2)) do s
        println("Sample $s/$(size(counts, 2))")
        
        # Calculate logarithmic counts for the sample
        sample_log_counts = log2.(counts[control_peaks, s])
        
        # Polynomial function to fit
        poly(x, c) = c[1] + c[2] * x + c[3] * x^2 + c[4] * x^3 + c[5] * x^4
        
        # Error function for the polynomial fit (mean squated error)
        err_func(coefs) = mean((sample_log_counts[k] - poly(control_peak_gc[k], coefs))^2 
            for k in 1:length(sample_log_counts))
            
        # Fit the polynomial to the data
        res = bboptimize(err_func, SearchRange=(-10_000.0, 10_000.0), NumDimensions=5, 
        Method=:de_rand_1_bin, MaxSteps=50_000, TraceMode=:silent)
        coefs = best_candidate(res)
            
        return coefs
    end;
    
    # Divide counts by the models, after converting the logarithmic model to regular scale
    println("\nCorrecting counts")
    gc_uniq = sort(unique(gc_frac))
    gc_idx = [only(searchsorted(gc_uniq, f)) for f=gc_frac]
    counts_corrected = copy(counts)
    for s in 1:size(counts, 2)
        println("Sample $s/$(size(counts, 2))")
        factor = [2^poly(f, gc_models[s]) for f=gc_uniq]
        counts_corrected[:,s] ./= factor[gc_idx]
    end

    return counts_corrected
end




