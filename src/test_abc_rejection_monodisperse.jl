workspace()

include("simulate_system_monodisperse.jl")
include("distance2.jl")

function test_abc_rejection_monodisperse()
	#Inititalization.
	const t_start::Int64 = convert(Int64, time_ns())
	srand(1)
	
	# Experimental parameters.
	ax::Float64 = 10.0 # µm.
	ay::Float64 = 10.0 # µm.
	L::Float64 = 15.0#100.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(40)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 3
	
	# Simulate true system.
	D_real::Float64 = 2.5 # µm^2/s.
	c_real::Float64 = 1e10 # part/ml.
	az_real::Float64 = 2.0 # µm.
	
	(K_real, DE_real) = simulate_system_monodisperse(D_real, c_real, ax, ay, az_real, L, number_of_frames, deltat, kmin)
	#println(mean(K_real))
	# Parameter bounds for inference.
	lb_D::Float64 = 0.0
	ub_D::Float64 = 5.0
	lb_c::Float64 = 0.5e10
	ub_c::Float64 = 1.5e10
	lb_az::Float64 = 1.0
	ub_az::Float64 = 3.0
        
	# Inference parameters.
	number_of_abc_samples::Int64 = 100000

	D_sim::Float64 = 0.0
	c_sim::Float64 = 0.0
	az_sim::Float64 = 0.0
	
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	file_name_output = join(("abc_sample_monodisperse_", string(random_seed), ".dat"))
	#file_name_output = "abc_sample.dat" 
	file_stream_output = open(file_name_output, "w")
	
	
	for current_abc_sample = 1:number_of_abc_samples
		if mod(current_abc_sample, 100) == 0
			println(current_abc_sample)
		end
		
		D_sim = lb_D + (ub_D - lb_D) * rand()
		c_sim = lb_c + (ub_c - lb_c) * rand()
		az_sim = lb_az + (ub_az - lb_az) * rand()
		
		(K_sim, DE_sim) = simulate_system_monodisperse(D_sim, c_sim, ax, ay, az_sim, L, number_of_frames, deltat, kmin)
		#@time 
		#println(mean(K_sim))
		
		dist = distance2(K_real, DE_real, K_sim, DE_sim)
		#@time 
		write(file_stream_output, D_sim, c_sim, az_sim, dist)
	end
	close(file_stream_output)
	
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(t_exec/1e9)
	nothing
end

#while true
    test_abc_rejection_monodisperse()
#end
#@profile test_abc_rejection()
#Profile.print()