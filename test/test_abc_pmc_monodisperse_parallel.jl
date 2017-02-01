workspace()

@everywhere include("../src/simulate_system.jl")
@everywhere include("../src/distance.jl")

function test_abc_pmc_lognormal_parallel()
	#Inititalization.
	srand(1)
	t_start::Int64 = convert(Int64, time_ns())

	output_dir::String = "output"
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	
	# Acquisition parameters.
	ax::Float64 = 5.0 # �m.
	ay::Float64 = 5.0 # �m.
	Lx::Float64 = 10.0 # �m.
	Ly::Float64 = 10.0 # �m.
	Lz::Float64 = 10.0 # �m.
	number_of_frames::Array{Int64, 1} = 250 * ones(40)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	
	# True system parameters.
	distribution_class::String = "monodisperse"
	D_real::Float64 = 2.5 # �m^2/s.
	c_real::Float64 = 1e10 # part/ml.
	az_real::Float64 = 2.0 # �m.
	
	# Simulate system.
	(K_real, DE_real) = simulate_system(distribution_class, [D_real], c_real, ax, ay, az_real, Lx, Ly, Lz, number_of_frames, deltat, kmin)

	# Inference starts.
	random_seed::Int64 = 1#convert(Int64, time_ns())
	srand(random_seed)
		
	# Parameter bounds for inference.
	lb_D::Float64 = 0.01
	ub_D::Float64 = 10.0
	lb_c::Float64 = 0.25 * c_real
	ub_c::Float64 = 4.0 * c_real
	lb_az::Float64 = 0.25 * az_real
	ub_az::Float64 = 4.0 * az_real
			
	# Inference parameters.
	number_of_abc_samples::Int64 = 128
	number_of_iterations::Int64 = 5000

	# Variables for population parameter values.
	D::Array{Float64, 1} = zeros(number_of_abc_samples)
	c::Array{Float64, 1} = zeros(number_of_abc_samples)
	az::Array{Float64, 1} = zeros(number_of_abc_samples)
	dist::Array{Float64, 1} = zeros(number_of_abc_samples)
	
	D_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	c_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	az_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	dist_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)

	# First iteration, assuming epsilon = inf so everything is accepted.
	D = lb_D + (ub_D - lb_D) * rand(number_of_abc_samples)
	c = lb_c + (ub_c - lb_c) * rand(number_of_abc_samples)
	az = lb_az + (ub_az - lb_az) * rand(number_of_abc_samples)
	
	w::Array{Float64, 1} = ones(number_of_abc_samples) / convert(Float64, number_of_abc_samples)
	w_star::Array{Float64, 1} = zeros(number_of_abc_samples)
	
	tau_D::Float64 = sqrt( 2.0 * var(D, corrected = false) )
	tau_c::Float64 = sqrt( 2.0 * var(c, corrected = false) )
	tau_az::Float64 = sqrt( 2.0 * var(az, corrected = false) )
	
	# The rest of the iterations.
	gamma = 5.0
	delta_gamma = 0.001
	epsilon::Float64 = 10^gamma
	trial_count::SharedArray{Int64, 1} = [0]
	trial_count_target::Int64 = 10 * number_of_abc_samples
	t_start_iteration::Int64 = 0 
	for current_iteration = 1:number_of_iterations
		t_start_iteration = convert(Int64, time_ns())
		trial_count[1] = 0	
		gamma = gamma - delta_gamma
		epsilon = 10^gamma
		#epsilon = 0.99 * epsilon
		@sync @parallel for current_abc_sample = 1:number_of_abc_samples
			idx = findfirst(cumsum(w) .>= rand())
			
			D_prim = D[idx]
			c_prim = c[idx]
			az_prim = az[idx]
			
			D_bis = 0.0
			c_bis = 0.0
			az_bis = 0.0

			dist_bis = Inf
			
			
			while dist_bis > epsilon 
				delta_D = tau_D * randn()
				D_bis = D_prim + delta_D
				if D_bis < lb_D
					D_bis = lb_D
					delta_D = D_bis - D_prim
				elseif D_bis > ub_D
					D_bis = ub_D
					delta_D = D_bis - D_prim
				end
							
				delta_c = tau_c * randn()
				c_bis = c_prim + delta_c
				if c_bis < lb_c
					c_bis = lb_c
					delta_c = c_bis - c_prim
				elseif c_bis > ub_c
					c_bis = ub_c
					delta_c = c_bis - c_prim
				end
				
				delta_az = tau_az * randn()
				az_bis = az_prim + delta_az
				if az_bis < lb_az
					az_bis = lb_az
					delta_az = az_bis - az_prim
				elseif az_bis > ub_az
					az_bis = ub_az
					delta_az = az_bis - az_prim
				end
				
				(K_sim, DE_sim) = simulate_system(distribution_class, [D_bis], c_bis, ax, ay, az_bis, Lx, Ly, Lz, number_of_frames, deltat, kmin)

				dist_bis = distance(K_real, DE_real, K_sim, DE_sim)
				#println(dist)
				
				trial_count[1] = trial_count[1] + 1

				#convert(Int64, time_ns()) - t_start_iteration < 
			end
			
			
			D_star[current_abc_sample] = D_bis
			c_star[current_abc_sample] = c_bis
			az_star[current_abc_sample] = az_bis
			dist_star[current_abc_sample] = dist_bis
		end
		
		println((current_iteration, trial_count[1], gamma, delta_gamma, tau_D, tau_c, tau_az))
		
		#if trial_count[1] > trial_count_target
		#	delta_gamma = 0.99 * delta_gamma
		#else
		#	delta_gamma = 1.01 * delta_gamma
		#end
		#println(trial_count[1])
		
		
		for current_abc_sample = 1:number_of_abc_samples
			w_star[current_abc_sample] = 0.0
			for i = 1:number_of_abc_samples
				w_star[current_abc_sample] = w_star[current_abc_sample] + w[i] *normpdf(D_star[current_abc_sample] - D[i], 0.0, tau_D) * normpdf(c_star[current_abc_sample] - c[i], 0.0, tau_c) * normpdf(az_star[current_abc_sample] - az[i], 0.0, tau_az)
			end
			w_star[current_abc_sample] = 1.0 / w_star[current_abc_sample]
		end
		
		w = w_star / sum(w_star)
		
		D = D_star
		c = c_star
		az = az_star
		dist = dist_star 
		
		tau_D = sqrt( 2.0 * var(D, corrected = false) )
		tau_c = sqrt( 2.0 * var(c, corrected = false) )
		tau_az = sqrt( 2.0 * var(az, corrected = false) )
		
		# Write intermediate result to file.
		file_name_output = join((output_dir, "/", "abc_pmc_md_par_it_", string(current_iteration), ".dat"))
		file_stream_output = open(file_name_output, "w")
		for current_abc_sample = 1:number_of_abc_samples
			write(file_stream_output, D[current_abc_sample], c[current_abc_sample], az[current_abc_sample], dist[current_abc_sample])
		end
		close(file_stream_output)
	end
	
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(t_exec/1e9)
	nothing
end

function normpdf(x, mu, sigma)
	return 1.0 / ( sqrt(2.0 * pi) * sigma) * exp( - 0.5 * (x - mu)^2 / sigma^2 )
end	

test_abc_pmc_lognormal_parallel()
