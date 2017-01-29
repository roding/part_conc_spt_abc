workspace()

include("../src/simulate_system.jl")
include("../src/distance.jl")

function test_abc_rejection()
	#Inititalization.
	const t_start::Int64 = convert(Int64, time_ns())
	srand(1)
	
	# Acquisition parameters.
	ax::Float64 = 40.0 # �m.
	ay::Float64 = 40.0 # �m.
	L::Float64 = 100.0 # �m.
	number_of_frames::Array{Int64, 1} = 100 * ones(20000)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	
	# True system parameters.
	distribution_class::String = "monodisperse"
	D_real::Float64 = 2.5 # �m^2/s.
	distribution_parameters_real::Array{Float64, 1} = [D_real]
	
#	distribution_class::String = "lognormal"
#	m_real::Float64 = 2.5 # �m^2/s.
#	s_real::Float64 = 0.5 # �m^2/s.
#	mu_real::Float64 = log(m_real) - 0.5 * log(1 + s_real^2/m_real^2)
#	sigma_real::Float64 = sqrt(log(1 + s_real^2/m_real^2))
#	distribution_parameters_real::Array{Float64, 1} = [mu_real, sigma_real]
	
	c_real::Float64 = 1e10 # part/ml.
	az_real::Float64 = 2.0 # �m.
	
	# Simulate system.
	(K_real, DE_real) = simulate_system(distribution_class, distribution_parameters_real, c_real, ax, ay, az_real, L, number_of_frames, deltat, kmin)

	# Parameter bounds for inference.
	lb = Array(Float64, 1)
	ub = Array(Float64, 1)
	if distribution_class == "monodisperse"
		lb = [0.5 * D_real, 0.5 * c_real, 0.5 * az_real]
		ub = [2.0 * D_real, 2.0 * c_real, 2.0 * az_real]
	elseif distribution_class == "lognormal"
		lb = [0.5 * mu_real, 0.5 * sigma_real, 0.5 * c_real, 0.5 * az_real]
		ub = [2.0 * D_real, 2.0 * sigma_real, 2.0 * c_real, 2.0 * az_real]
	end
	
	# Inference parameters.
	number_of_abc_samples::Int64 = 1000

	mu_sim::Float64 = 0.0
	sigma_sim::Float64 = 0.0
	c_sim::Float64 = 0.0
	az_sim::Float64 = 0.0
	
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	file_name_output = join(("abc_sample_", string(random_seed), ".dat"))
	#file_name_output = "abc_sample.dat" 
	file_stream_output = open(file_name_output, "w")
	
	
	for current_abc_sample = 1:number_of_abc_samples
		if mod(current_abc_sample, 100) == 0
			println(current_abc_sample)
		end
		
		mu_sim = lb_mu + (ub_mu - lb_mu) * rand()
		sigma_sim = lb_sigma + (ub_sigma - lb_sigma) * rand()
		c_sim = lb_c + (ub_c - lb_c) * rand()
		az_sim = lb_az + (ub_az - lb_az) * rand()
		
		(K_sim, RSQ_sim) = simulate_system(mu_sim, sigma_sim, c_sim, ax, ay, az_sim, L, number_of_frames, deltat, kmin)
		#@time 
		#println(mean(K_sim))
		
		dist = distance2(K_real, RSQ_real, K_sim, RSQ_sim)
		#@time 
		write(file_stream_output, mu_sim, sigma_sim, c_sim, az_sim, dist)
	end
	close(file_stream_output)
	
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(t_exec/1e9)
	nothing
end

while true
    test_abc_rejection()
end
#@profile test_abc_rejection()
#Profile.print()