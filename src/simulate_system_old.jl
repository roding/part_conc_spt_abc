include("rand_poisson.jl")

function simulate_system(mu::Float64, sigma::Float64, c::Float64, ax::Float64, ay::Float64, az::Float64, L::Float64, number_of_frames::Array{Int64, 1}, deltat::Float64, kmin::Int64)
	# Intensity of Poisson distribution of the number of particles. The factor 
	# 1e12 takes into account that concentration is specified in particles/ml.
	lambda::Float64 = c * L^3 / 1e12
	
	# Number of particles in simulated system, which varies between videos.
	number_of_particles::Int64 = 0
	
	# Lower and upper bounds for the detection region.
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
	D::Float64 = 0.0
	s::Float64 = 0.0
	number_of_videos::Int64 = length(number_of_frames)
	
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
			
			# Generate random diffusion coefficent from distribution.
			D = exp(mu + sigma * rand())
			s = sqrt(2 * D * deltat)
			
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

