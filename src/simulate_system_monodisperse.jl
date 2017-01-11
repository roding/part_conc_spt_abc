function simulate_system_monodisperse(D::Float64, c::Float64, ax::Float64, ay::Float64, az::Float64, L::Float64, number_of_frames::Array{Int64, 1}, deltat::Float64, kmin::Int64)
	# Convert concentration to a number of particles. The factor 1e12 takes into ack
	# that concentration is specified in particles/ml. Also enlarge simulation domain
	# to accomodate an integer number of particle while maintaining correct concentration.
	number_of_particles_temp::Float64 = c * L^3 / 1e12
	number_of_particles::Int64 = rand_poisson(number_of_particles_temp)
	#L = (convert(Float64, number_of_particles) / number_of_particles_temp)^(1.0/3.0) * L
	#println((number_of_particles_temp, number_of_particles))
	# Define lower and upper bounds for the detection region.
	lbx::Float64 = 0.5 * (L - ax)
	ubx::Float64 = 0.5 * (L + ax)
	lby::Float64 = 0.5 * (L - ay)
	uby::Float64 = 0.5 * (L + ay)
	lbz::Float64 = 0.5 * (L - az)
	ubz::Float64 = 0.5 * (L + az)
	
	# Simulate.
	k::Int64 = 0
	rsq::Float64 = 0.0
	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	deltax::Float64 = 0.0
	deltay::Float64 = 0.0
	deltaz::Float64 = 0.0
	number_of_videos::Int64 = length(number_of_frames)
	
	s::Float64 = sqrt(2 * D * deltat) # Standard deviation of random displacements.
	
	K::Array{Int64, 1} = zeros(0)
	RSQ::Array{Float64, 1} = zeros(0)
	
	for current_video = 1:number_of_videos
		for current_particle = 1:number_of_particles
			# Random initial position.
			x = L * rand()
			y = L * rand()
			z = L * rand()
			
			# Reset rsq.
			rsq = 0.0
			
			# In detection region or not? Check z limits first because they are more restrictive.
			if (lbz <= z <= ubz) & (lbx <= x <= ubx) & (lby <= y <= uby)
				k = 1
			else
				k = 0
			end
			
			# Let particle diffuse through remainder of video.
			for current_frame = 2:number_of_frames[current_video]
				deltax  = s * randn()
				deltay  = s * randn()
				deltaz  = s * randn()
				x = x + deltax
				y = y + deltay
				z = z + deltaz
				
				if x > L
					x = x - L
				elseif x < 0.0
					x = x + L
				end
				
				if y > L
					y = y - L
				elseif y < 0.0
					y = y + L
				end
				
				if z > L
					z = z - L
				elseif z < 0.0
					z = z + L
				end
				
				if (lbz <= z <= ubz) & (lbx <= x <= ubx) & (lby <= y <= uby)
					k = k + 1
					rsq = rsq + deltax^2 + deltay^2 # Only in x-y plane.
				elseif k > 0
					if k >= kmin
						push!(K, k)
						push!(RSQ, rsq)
					end
					k = 0
					rsq = 0.0
				end
			end
		end
	end
	
	return (K, RSQ)
end

