function simulate_system(mu::Float64, sigma::Float64, c::Float64, ax::Float64, ay::Float64, az::Float64, A::Float64, number_of_frames::Array{Int64, 1})
	
	#kmax::Int64 = maximum(nFrames)
	#H::Array{Int64,1} = zeros(kmax)
	
	# Convert concentration to a number of particles. The factor 1e12 takes into account
	# that concentration is specified in particles/ml. Also enlarge simulation domain
	# to accomodate an integer number of particle while maintaining correct concentration.
	number_of_particles_temp::Float64 = convert(Int64, c * A^3 / 1e12)
	number_of_particles::Int64 = convert(Int64, ceil(number_of_particles_temp))
	A = (convert(Float64, number_of_particles) / number_of_particles_temp)^(1.0/3.0) * A
	
	

end

