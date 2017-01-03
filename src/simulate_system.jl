function simulate_system(mu::Float64, sigma::Float64, c::Float64, ax::Float64, ay::Float64, az::Float64, L::Float64, number_of_frames::Array{Int64, 1})
	
	#kmax::Int64 = maximum(nFrames)
	#H::Array{Int64,1} = zeros(kmax)
	
	# Convert concentration to a number of particles. The factor 1e12 takes into account
	# that concentration is specified in particles/ml. Also enlarge simulation domain
	# to accomodate an integer number of particle while maintaining correct concentration.
	number_of_particles_temp::Float64 = convert(Int64, c * L^3 / 1e12)
	number_of_particles::Int64 = convert(Int64, ceil(number_of_particles_temp))
	L = (convert(Float64, number_of_particles) / number_of_particles_temp)^(1.0/3.0) * L
	
	# Define lower and upper bounds for the detection region.
	lbx::Float64 = 0.5 * (L - ax)
	ubx::Float64 = 0.5 * (L + ax)
	lby::Float64 = 0.5 * (L - ay)
	uby::Float64 = 0.5 * (L + ay)
	lbz::Float64 = 0.5 * (L - az)
	ubz::Float64 = 0.5 * (L + az)
	
	# Simulate.
	count::Int64 = 0
	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	D::Float64 = 0.0
	s::Float64 = 0.0
	number_of_videos::Int64 = length(number_of_frames)
	
	for current_video = 1:number_of_videos
		for current_particle = 1:number_of_particles
			# Random initial position.
			x = L * rand()
			y = L * rand()
			z = L * rand()
			
			# In detection region or not?
			if (lbz <= z) & (z <= ubz) & (lbx <= x) & (x <= ubx) & (lby <= y) & (y <= uby)
				count = 1
			else
				count = 0
			end
			
			# Generate random diffusion coefficent from distribution.
			
			
			for current_frame = 2:number_of_frames[current_video]
				x = x + sigma*randn()
				y += sigma*randn()
				z += sigma*randn()
				if x > A_corr
					x -= A_corr
				elseif x < 0.0
					x += A_corr
				end
				if y > A_corr
					y -= A_corr
				elseif y < 0.0
					y += A_corr
				end
				if z > A_corr
					z -= A_corr
				elseif z < 0.0
					z += A_corr
				end
				if (lbz <= z) & (z <= ubz) & (lbx <= x) & (x <= ubx) & (lby <= y) & (y <= uby)
					count += 1
				elseif count > 0
					Hsim[count] += 1
					count = 0
				end
			end
		
	for currentVideo = 1:nVideos
		for currentParticle = 1:nParticles
			
			
			
			if count > 0
				Hsim[count] += 1
			end
		end
	end

end

