workspace()

@everywhere include("../src/simulate_system.jl")
@everywhere include("../src/distance.jl")

function test_abc_pmc_lognormal_parallel()
	#Inititalization.
	srand(1)
	t_start::Int64 = convert(Int64, time_ns())
	
	# Acquisition parameters.
	ax::Float64 = 40.0 # µm.
	ay::Float64 = 40.0 # µm.
	Lx::Float64 = 50.0#100.0 # µm.
	Ly::Float64 = 50.0#100.0 # µm.
	Lz::Float64 = 10.0#50.0 # µm.
	number_of_frames::Array{Int64, 1} = 250 * ones(40)
	deltat::Float64 = 0.05 # seconds
	kmin::Int64 = 2
	
	# True system parameters.
	distribution_class::String = "lognormal"
	m_real::Float64 = 2.5 # µm^2/s.
	s_real::Float64 = 0.5 # µm^2/s.
	c_real::Float64 = 1e9#1e10 # part/ml.
	az_real::Float64 = 2.0 # µm.
	
	# Simulate system.
	(K_real, DE_real) = simulate_system(distribution_class, [m_real, s_real], c_real, ax, ay, az_real, Lx, Ly, Lz, number_of_frames, deltat, kmin)

	# Inference starts.
	random_seed::Int64 = 1#convert(Int64, time_ns())
	srand(random_seed)
		
	# Parameter bounds for inference.
	lb_m::Float64 = 0.01
	ub_m::Float64 = 10.0
	lb_s::Float64 = 0.01
	ub_s::Float64 = 5.0
	lb_c::Float64 = 0.25 * c_real
	ub_c::Float64 = 4.0 * c_real
	lb_az::Float64 = 0.25 * az_real
	ub_az::Float64 = 4.0 * az_real
			
	# Inference parameters.
	number_of_abc_samples::Int64 = 64
	number_of_iterations::Int64 = 500

	# Variables for population parameter values.
	m::Array{Float64, 1} = zeros(number_of_abc_samples)
	s::Array{Float64, 1} = zeros(number_of_abc_samples)
	c::Array{Float64, 1} = zeros(number_of_abc_samples)
	az::Array{Float64, 1} = zeros(number_of_abc_samples)
	
	m_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	s_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	c_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	az_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	
	# First iteration, assuming epsilon = inf so everything is accepted.
	m = lb_m + (ub_m - lb_m) * rand(number_of_abc_samples)
	s = lb_s + (ub_s - lb_s) * rand(number_of_abc_samples)
	c = lb_c + (ub_c - lb_c) * rand(number_of_abc_samples)
	az = lb_az + (ub_az - lb_az) * rand(number_of_abc_samples)
	
	w::Array{Float64, 1} = ones(number_of_abc_samples) / convert(Float64, number_of_abc_samples)
	w_star::Array{Float64, 1} = zeros(number_of_abc_samples)
	
	tau_m::Float64 = sqrt( 2.0 * var(m, corrected = false) )
	tau_s::Float64 = sqrt( 2.0 * var(s, corrected = false) )
	tau_c::Float64 = sqrt( 2.0 * var(c, corrected = false) )
	tau_az::Float64 = sqrt( 2.0 * var(az, corrected = false) )
	
#	m_prim::Float64 = 0.0
#	s_prim::Float64 = 0.0
#	c_prim::Float64 = 0.0
#	az_prim::Float64 = 0.0
#	
#	m_bis::Float64 = 0.0
#	s_bis::Float64 = 0.0
#	c_bis::Float64 = 0.0
#	az_bis::Float64 = 0.0
#	
#	delta_m::Float64 = 0.0
#	delta_s::Float64 = 0.0
#	delta_c::Float64 = 0.0
#	delta_az::Float64 = 0.0
	
#	dist::Float64 = 0.0
#	idx::Int64 = 0
	
	# The rest of the iterations.
	epsilon = 1e6
	for current_iteration = 1:number_of_iterations
		epsilon = 0.99 * epsilon		
		@sync @parallel for current_abc_sample = 1:number_of_abc_samples
			println((current_iteration, epsilon, current_abc_sample, tau_m, tau_s, tau_c, tau_az))
			idx = findfirst(cumsum(w) .>= rand())
			
			m_prim = m[idx]
			s_prim = s[idx]
			c_prim = c[idx]
			az_prim = az[idx]
			
			m_bis = 0.0
			s_bis = 0.0
			c_bis = 0.0
			az_bis = 0.0

			dist = Inf
			while dist > epsilon
				delta_m = tau_m * randn()
				m_bis = m_prim + delta_m
				if m_bis < lb_m
					m_bis = lb_m
					delta_m = m_bis - m_prim
				elseif m_bis > ub_m
					m_bis = ub_m
					delta_m = m_bis - m_prim
				end
				
				delta_s = tau_s * randn()
				s_bis = s_prim + delta_s
				if s_bis < lb_s
					s_bis = lb_s
					delta_s = s_bis - s_prim
				elseif s_bis > ub_s
					s_bis = ub_s
					delta_s = s_bis - s_prim
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
				
				
				(K_sim, DE_sim) = simulate_system(distribution_class, [m_bis, s_bis], c_bis, ax, ay, az_bis, Lx, Ly, Lz, number_of_frames, deltat, kmin)

				dist = distance(K_real, DE_real, K_sim, DE_sim)
				#println(dist)
			end
			
			m_star[current_abc_sample] = m_bis
			s_star[current_abc_sample] = s_bis
			c_star[current_abc_sample] = c_bis
			az_star[current_abc_sample] = az_bis
		end
		
		
		for current_abc_sample = 1:number_of_abc_samples
			w_star[current_abc_sample] = 0.0
			for i = 1:number_of_abc_samples
				w_star[current_abc_sample] = w_star[current_abc_sample] + w[i] *normpdf(m_star[current_abc_sample] - m[i], 0.0, tau_m) * normpdf(s_star[current_abc_sample] - s[i], 0.0, tau_s) * normpdf(c_star[current_abc_sample] - c[i], 0.0, tau_c) * normpdf(az_star[current_abc_sample] - az[i], 0.0, tau_az)
			end
			w_star[current_abc_sample] = 1.0 / w_star[current_abc_sample]
		end
		
		w = w_star / sum(w_star)
		
		m = m_star
		s = s_star
		c = c_star
		az = az_star
		
		tau_m = sqrt( 2.0 * var(m, corrected = false) )
		tau_s = sqrt( 2.0 * var(s, corrected = false) )
		tau_c = sqrt( 2.0 * var(c, corrected = false) )
		tau_az = sqrt( 2.0 * var(az, corrected = false) )
		
	end
	
	file_name_output = join(("abc_pmc_sample_lognormal_", string(random_seed), ".dat"))
	file_stream_output = open(file_name_output, "w")
	for current_abc_sample = 1:number_of_abc_samples
		write(file_stream_output, m[current_abc_sample], s[current_abc_sample], c[current_abc_sample], az[current_abc_sample])
	end
	close(file_stream_output)
	
	t_exec::Int64 = convert(Int64, time_ns()) - t_start
	println(t_exec/1e9)
	nothing
end

function normpdf(x, mu, sigma)
	return 1.0 / ( sqrt(2.0 * pi) * sigma) * exp( - 0.5 * (x - mu)^2 / sigma^2 )
end	

test_abc_pmc_lognormal_parallel()