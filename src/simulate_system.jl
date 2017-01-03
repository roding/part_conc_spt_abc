function simulate_system(mu::Float64, sigma::Float64, c::Float64, ax::Float64, ay::Float64, az::Float64, A::Float64, number_of_frames::Array{Int64, 1})
	
	#kmax::Int64 = maximum(nFrames)
	#H::Array{Int64,1} = zeros(kmax)
	
	# Convert concentration to a number of particles. The factor 1e12 takes into account
	# that concentration is specified in particles/ml. Also enlarge simulation domain
	# to accomodate an integer number of particle while maintaining correct concentration.
	number_of_particles_temp::Float64 = convert(Int64, c * A^3 / 1e12)
	number_of_particles::Int64 = convert(Int64, ceil(number_of_particles_temp))
	A = (convert(Float64, number_of_particles) / number_of_particles_temp)^(1.0/3.0) * A
	
	# Define lower and upper bounds for the detection region.
	xlo::Float64 = 0.5 * (A - ax)
	xhi::Float64 = 0.5 * (A + ax)
	ylo::Float64 = 0.5 * (A - ay)
	yhi::Float64 = 0.5 * (A + ay)
	zlo::Float64 = 0.5 * (A - az)
	zhi::Float64 = 0.5 * (A + az)
	
	
		count = 0
	for currentVideo = 1:nVideos
		for currentParticle = 1:nParticles
			x = A_corr * rand()
			y = A_corr * rand()
			z = A_corr * rand()
			if (zlo <= z) & (z <= zhi) & (xlo <= x) & (x <= xhi) & (ylo <= y) & (y <= yhi)
				count = 1
			else
				count = 0
			end
			for currentFrame = 2:nFrames[currentVideo]
				x += sigma*randn()
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
				if (zlo <= z) & (z <= zhi) & (xlo <= x) & (x <= xhi) & (ylo <= y) & (y <= yhi)
					count += 1
				elseif count > 0
					Hsim[count] += 1
					count = 0
				end
			end
			if count > 0
				Hsim[count] += 1
			end
		end
	end

end

