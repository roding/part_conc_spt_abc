workspace()

include("simulate_system.jl")
include("distance.jl")

function test_abc()
	#Inititalization.
	const t_start::Int64 = convert(Int64, time_ns())
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Experimental parameters.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	L::Float64 = 50.0#100.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(40)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 3
	
	# Simulate true system.
	mu_real::Float64 = 1.0 # µm^2/s.
	sigma_real::Float64 = 1.0 # µm^2/s.
	c_real::Float64 = 5e6 # part/ml.
	az_real::Float64 = 2.0 # µm.
	
	(K_real, RSQ_real) = simulate_system(mu_real, sigma_real, c_real, ax, ay, az_real, L, number_of_frames, deltat, kmin)
	
	# Parameter bounds for inference.
	lb_mu::Float64 = 0.5
	ub_mu::Float64 = 1.5
	lb_sigma::Float64 = 0.5
	ub_sigma::Float64 = 1.5
	lb_c::Float64 = 1e6
	ub_c::Float64 = 1e7
	lb_az::Float64 = 1.5
	ub_az::Float64 = 2.5
		
	# Inference parameters.
	number_of_abc_samples::Int64 = 1000000

	mu_sim::Float64 = 0.0
	sigma_sim::Float64 = 0.0
	c_sim::Float64 = 0.0
	az_sim::Float64 = 0.0
	
	#file_name_output = join(("abc_sample_", string(random_seed), ".dat"))
	file_name_output = "abc_sample.dat" 
	file_stream_output = open(file_name_output, "w")
	
	for current_abc_sample = 1:number_of_abc_samples
		println(current_abc_sample/number_of_abc_samples)
		
		mu_sim = lb_mu + (ub_mu - lb_mu) * rand()
		sigma_sim = lb_sigma + (ub_sigma - lb_sigma) * rand()
		c_sim = lb_c + (ub_c - lb_c) * rand()
		az_sim = lb_az + (ub_az - lb_az) * rand()
		
		(K_sim, RSQ_sim) = simulate_system(mu_sim, sigma_sim, c_sim, ax, ay, az_sim, L, number_of_frames, deltat, kmin)
		
		dist = distance(K_real, RSQ_real, K_sim, RSQ_sim)

		write(file_stream_output, mu_sim, sigma_sim, c_sim, az_sim, dist)
	end
	close(file_stream_output)

	nothing
end

test_abc()