include("rand_poisson.jl")
include("periodic.jl")

function simulate_system(distribution_class::String, distribution_parameters::Array{Float64, 1}, c::Float64, ax::Float64, ay::Float64, az::Float64, L::Float64, number_of_frames::Array{Int64, 1}, deltat::Float64, kmin::Int64)
	
	# Intensity of Poisson distribution of the number of particles. The factor 
	# 1e12 takes into account that concentration is specified in particles/ml.
	lambda::Float64 = c * L^3 / 1e12
	
	# Number of particles in simulated system, which varies between videos.
	number_of_particles::Int64 = 0
	
	# Number of videos to be simulated.
	number_of_videos::Int64 = length(number_of_frames)
	
	# Lower and upper bounds for the detection region.
	lbx::Float64 = 0.5 * (L - ax)
	ubx::Float64 = 0.5 * (L + ax)
	lby::Float64 = 0.5 * (L - ay)
	uby::Float64 = 0.5 * (L + ay)
	lbz::Float64 = 0.5 * (L - az)
	ubz::Float64 = 0.5 * (L + az)
	
	# Variables for storing positions and displacements.
	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	deltax::Float64 = 0.0
	deltay::Float64 = 0.0
	deltaz::Float64 = 0.0
	
	# Number of positions and estimated diffusion coefficient of 
	# current trajectory.
	k::Int64 = 0
	de::Float64 = 0.0
	
	# Number of positions and estimated diffusion coefficient of all 
	# recorded trajectories.
	K::Array{Int64, 1} = zeros(0)
	DE::Array{Float64, 1} = zeros(0)
	
	# Standard deviation of random displacements.
	s::Float64 = 0.0 

	for current_video = 1:number_of_videos
		number_of_particles = rand_poisson(lambda)
		
		for current_particle = 1:number_of_particles
			# Generate random diffusion coefficent from distribution.
			if distribution_class == "monodisperse"
				s = sqrt(2 * distribution_parameters[1] * deltat)
			elseif distribution_class == "lognormal"
				s = sqrt(2 * exp(distribution_parameters[1] + distribution_parameters[2] * rand()) * deltat)
			end
			
			# Random initial position.
			x = L * rand()
			y = L * rand()
			z = L * rand()
			
			# Reset empirical diffusion coefficient.
			de = 0.0
			
			# In detection region or not? Check z limits first because they are more restrictive.
			if (lbz <= z <= ubz) & (lbx <= x <= ubx) & (lby <= y <= uby)
				k = 1
			else
				k = 0
			end
			
			# Let particle diffuse through remainder of video.
			for current_frame = 2:number_of_frames[current_video]
				deltax = s * randn()
				deltay = s * randn()
				deltaz = s * randn()
				
				x = x + deltax
				y = y + deltay
				z = z + deltaz
				
				x = periodic(x, L)
				y = periodic(y, L)
				z = periodic(z, L)
				
				if (lbz <= z <= ubz) & (lbx <= x <= ubx) & (lby <= y <= uby)
					k = k + 1
					if k >= 2 # The first point inside the detection region does not contributed toward the estimated diffusion coefficient.
						de = de + deltax^2 + deltay^2 # Only in x-y plane.
					end
				elseif k > 0
					if k >= kmin
						push!(K, k)
						if k >= 2
							de = de / (convert(Float64, k - 1) * 4.0 * deltat) # The '4' comes from the 2-D obs.
						else
							de = 0.0
						end
						push!(DE, de)
					end
					k = 0
					de = 0.0
				end
			end
		end
	end
	
	return (K, DE)
end

