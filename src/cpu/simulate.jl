function simulate(	model::String,
					m::Array{Float64, 1},
					c::Array{Float64, 1},
					ax::Float64,
					ay::Float64,
					az::Array{Float64, 1},
					Lx::Float64,
					Ly::Float64,
					Lz::Float64,
					number_of_frames::Array{Int64, 1},
					deltat::Float64,
					kmin::Int64,
					number_of_de_bins::Int64,
					ub_de::Float64)

	# Number of components in distribution.
	number_of_components::Int64 = length(m)
	fractions::Array{Float64, 1} = c / sum(c)
	cum_fractions::Array{Float64, 1} = cumsum(fractions)

	# Intensity of Poisson distribution of the number of particles. The factor
	# 1e12 takes into account that concentration is specified in particles/ml.
	lambda::Float64 = sum(c) * Lx * Ly * Lz / 1e12
	kmax::Int64 = maximum(number_of_frames)

	# Number of videos.
	number_of_videos::Int64 = length(number_of_frames)

	# Lower and upper bounds for the detection region.
	lbx::Float64 = 0.5 * (Lx - ax)
	ubx::Float64 = 0.5 * (Lx + ax)
	lby::Float64 = 0.5 * (Ly - ay)
	uby::Float64 = 0.5 * (Ly + ay)
	lbz::Array{Float64, 1} = 0.5 * (Lz - az)
	ubz::Array{Float64, 1} = 0.5 * (Lz + az)

	# Variables for storing positions, displacements, number of particles, etc
	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	deltax::Float64 = 0.0
	deltay::Float64 = 0.0
	deltaz::Float64 = 0.0
	number_of_particles::Int64 = 0
	k::Int64 = 0
	de::Float64 = 0.0
	std_dev_random_walk::Float64 = 0.0
	index::Int64 = 0

	# Histogram matrix for number of positions and estimated diffusion coefficients of all
	# recorded trajectories.
	H::Array{Int64, 2} = zeros(kmax, number_of_de_bins)
	d_de::Float64 = ub_de / convert(Float64, number_of_de_bins)
	idx::Int64 = 0

	for current_video = 1:number_of_videos
		number_of_particles = rand_poisson(lambda)

		for current_particle = 1:number_of_particles
			# Generate random diffusion coefficent from distribution, or more precisely,
			# a random standard deviation for the displacements.
			index = rand_weighted_index(cum_fractions)
			std_dev_random_walk = sqrt(2.0 * m[index] * deltat)

			# Random initial position.
			x = Lx * rand()
			y = Ly * rand()
			z = Lz * rand()

			# Reset empirical diffusion coefficient.
			de = 0.0

			# In detection region or not? Check z limits first because they are more restrictive.
			#if (lbz[index] <= z <= ubz[index]) & (lbx <= x <= ubx) & (lby <= y <= uby)
			if (abs(z - 0.5 * Lz) <= 0.5 * az[index]) & (abs(x - 0.5 * Lx) <= 0.5 * ax) & (abs(y - 0.5 * Ly) <= 0.5 * ay)
				k = 1
			else
				k = 0
			end

			# Let particle diffuse through remainder of video.
			for current_frame = 2:number_of_frames[current_video]
				deltax = std_dev_random_walk * randn()
				deltay = std_dev_random_walk * randn()
				deltaz = std_dev_random_walk * randn()

				x = x + deltax
				y = y + deltay
				z = z + deltaz

				x = position_periodic(x, Lx)
				y = position_periodic(y, Ly)
				z = position_periodic(z, Lz)

				#if (lbz[index] <= z <= ubz[index]) & (lbx <= x <= ubx) & (lby <= y <= uby)
				if (abs(z - 0.5 * Lz) <= 0.5 * az[index]) & (abs(x - 0.5 * Lx) <= 0.5 * ax) & (abs(y - 0.5 * Ly) <= 0.5 * ay)
					k = k + 1
					if k >= 2 # The first point inside the detection region does not contributed toward the estimated diffusion coefficient.
						de = de + deltax^2 + deltay^2 # Only in x-y plane.
					end
				elseif k > 0
					if k >= kmin
						de = de / (convert(Float64, k - 1) * 4.0 * deltat) # The '4' comes from the 2D observations.
						idx = convert(Int64, ceil(de / d_de))
						if 1 <= idx <= number_of_de_bins
							H[k, idx] += 1
						end
					end
					k = 0
					de = 0.0
				end
			end

			# If particle is inside detection region in last frame, we have now missed that trajectory and need to add it now.
			if k >= kmin # Assuming kmin >= 2 otherwise the division fails.
				de = de / (convert(Float64, k - 1) * 4.0 * deltat) # The '4' comes from the 2D observations.
				idx = convert(Int64, ceil(de / d_de))
				if 1 <= idx <= number_of_de_bins
					H[k, idx] += 1
				end
			end

		end
	end

	return H
end
