function setup_landscapes(forestCover, pvalue)

	# get generated landscape
	landscape = readdlm(string("Data/landscapes/forestcover_",forestCover,"_pvalue_",pvalue,"_edge.csv"), Int)

	# number of patches
	n_patch = size(landscape, 1)

	# landscape size
	n = Int(sqrt(n_patch))

	# grid of patch states (0=destroyed, 1=forest)
	x_state = reshape(landscape, (n,n))

	return n, n_patch, x_state
end


function get_interaction_network(network)

	# get network incidence matrix
	M_inc = readdlm(string("Data/",network,"/Minc.csv"), ' ', Int)

	# number of resources
	n_r = size(M_inc, 1)

	# number of consumers
	n_c = size(M_inc, 2)

	# get host forest occurrences
	pf_r = readdlm(string("Data/",network,"/host_forest_occurrence.csv"))

	return M_inc, n_r, n_c, pf_r
end


function get_theta(network)

	# theta of resources
	theta_r = readdlm(string("Data/",network,"/theta_r.csv"), Float32)

	# theta of consumers
	theta_c = readdlm(string("Data/",network,"/theta_c.csv"), Float32)

	return theta_r, theta_c
end

#initial occupancy equals to 100%:
#function setup_grids(n, n_r, n_c, theta_r, theta_c)
#
#	# grid of resources (1=present, 0=absent) for each species
#	x_r = ones(Int, n, n, n_r)
#
#	# grid of consumers (1=present, 0=absent) for each species
#	x_c = ones(Int, n, n, n_c)
#
#	# grid of resource trait values (initial trait = theta)
#	z_r = Array{Union{Missing, Float32}}(repeat(transpose(theta_r), outer=(n*n)))
#	z_r = reshape(z_r, (n,n,n_r))
# 
#	# grid of consumer trait values (initial trait = theta)
#	z_c = Array{Union{Missing, Float32}}(repeat(transpose(theta_c), outer=(n*n)))
#	z_c = reshape(z_c, (n,n,n_c))
#
#	return x_r, x_c, z_r, z_c
#end

#initial occupancy equals to 75%:
function setup_grids(n, n_r, n_c, theta_r, theta_c)

   # grid of resources (1=present, 0=absent) for each species
    x_r = rand(n, n, n_r) .< 0.75  # 0.75 means approximately 75% initial occupancy of every species
    x_r = Int.(x_r)

    # grid of consumers (1=present, 0=absent) for each species
    x_c = rand(n, n, n_c) .< 0.75  # 0.75 means approximately 75% initial occupancy of every species
    x_c = Int.(x_c)

    # grid of resource trait values (initial trait = theta)
    z_r = ifelse.(x_r .== 1, 
                  reshape(repeat(transpose(theta_r), outer=(n*n)), n, n, n_r),
                  missing)

    # grid of consumer trait values (initial trait = theta)
    z_c = ifelse.(x_c .== 1, 
                  reshape(repeat(transpose(theta_c), outer=(n*n)), n, n, n_c),
                  missing)

    return x_r, x_c, z_r, z_c
end

# 100% occupancy, but initial trait value sampled from uniform distribution between 10 and 20 instead of equal to theta
#function setup_grids(n, n_r, n_c, theta_r, theta_c)
#
#    # grid of resources (1=present, 0=absent) for each species
#    x_r = ones(Int, n, n, n_r)
#
#    # grid of consumers (1=present, 0=absent) for each species
#    x_c = ones(Int, n, n, n_c)
#
#    # grid of resource trait values (sampled from uniform distribution between 10 and 20)
#    z_r = reshape(Float32.(10 .+ 10 .* rand(n*n*n_r)), n, n, n_r)
#    z_r = convert(Array{Union{Missing, Float32}}, z_r)
#
#    # grid of consumer trait values (sampled from uniform distribution between 10 and 20)
#    z_c = reshape(Float32.(10 .+ 10 .* rand(n*n*n_c)), n, n, n_c)
#    z_c = convert(Array{Union{Missing, Float32}}, z_c)
#
#    return x_r, x_c, z_r, z_c
#end


function initialise_dataframes_store_results(tmax, n)

	# initialise dataframes for storing results
	df_dt = DataFrame(t = repeat(1:tmax, inner=n),
	                  species = repeat(1:n, outer=tmax),
					  abundance = Vector{Union{Missing, Float32}}(missing, tmax*n),
					  z_mean = Vector{Union{Missing, Float32}}(missing, tmax*n),
					  z_sd = Vector{Union{Missing, Float32}}(missing, tmax*n))

	return df_dt
end