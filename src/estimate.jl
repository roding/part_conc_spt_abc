function estimate(distribution_class::String,
				number_of_components::Int64,
				Lx::Float64,
				Ly::Float64,
				Lz::Float64,
				kmin::Int64,
				number_of_de_bins::Int64,
				ub_de::Float64,
				lb_m::Float64,
				ub_m::Float64,
				lb_s::Float64,
				ub_s::Float64,
				lb_c::Float64,
				ub_c::Float64,
				lb_az::Float64,
				ub_az::Float64,
				number_of_abc_samples::Int64,
				gamma_initial::Float64,
				delta_gamma::Float64,
				ax::Float64,
				ay::Float64,
				number_of_frames::Array{Int64, 1},
				deltat::Float64,
				K::Array{Int64, 1},
				DE::Array{Float64, 1})
				
	#Inititalization.
	#t_start::Int64 = convert(Int64, time_ns())
	
	# Distance function histogram bin width.
	d_de::Float64 = ub_de / convert(Float64, number_of_de_bins)
	
	# Convert data to histogram form.
	kmax::Int64 = maximum(number_of_frames)
	n_K::Array{Int64, 1} = zeros(kmax)
	n_DE::Array{Int64, 1} = zeros(number_of_de_bins)
	idx::Int64 = 0
	for i = 1:length(K)
		idx = convert(Int64, ceil(DE[i] / d_de))
		if 1 <= idx <= number_of_de_bins
			n_K[K[i]] += 1
			n_DE[idx] += 1
		end
	end
	
	# Variables for population parameter values.
	m::Array{Float64, 2} = zeros(number_of_abc_samples, number_of_components)
	s::Array{Float64, 2} = zeros(number_of_abc_samples, number_of_components)
	c::Array{Float64, 2} = zeros(number_of_abc_samples, number_of_components)
	az::Array{Float64, 1} = zeros(number_of_abc_samples)
	dist::Array{Float64, 1} = zeros(number_of_abc_samples)
	
	m_star::SharedArray{Float64, 2} = zeros(number_of_abc_samples, number_of_components)
	s_star::SharedArray{Float64, 2} = zeros(number_of_abc_samples, number_of_components)
	c_star::SharedArray{Float64, 2} = zeros(number_of_abc_samples, number_of_components)
	az_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)
	dist_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)

	# Initialize assuming that epsilon = inf so that everything is accepted.
	m = lb_m + (ub_m - lb_m) * rand(number_of_abc_samples, number_of_components)
	s = lb_s + (ub_s - lb_s) * rand(number_of_abc_samples, number_of_components)
	c = lb_c + (ub_c - lb_c) * rand(number_of_abc_samples, number_of_components)
	az = lb_az + (ub_az - lb_az) * rand(number_of_abc_samples)
	
	w::Array{Float64, 1} = ones(number_of_abc_samples) / convert(Float64, number_of_abc_samples)
	cum_w::Array{Float64, 1} = cumsum(w)
	w_star::Array{Float64, 1} = zeros(number_of_abc_samples)
	
	# Displacement standard deviations.
	tau_m::Array{Float64, 1} = zeros(number_of_components)
	tau_s::Array{Float64, 1} = zeros(number_of_components)
	tau_c::Array{Float64, 1} = zeros(number_of_components)
	for current_component = 1:number_of_components
		tau_m[current_component] = sqrt( 2.0 * var(m[:, current_component], corrected = false) )
		tau_s[current_component] = sqrt( 2.0 * var(s[:, current_component], corrected = false) )
		tau_c[current_component] = sqrt( 2.0 * var(c[:, current_component], corrected = false) )
	end
	tau_az::Float64 = sqrt( 2.0 * var(az, corrected = false) )
	
	# The rest of the iterations.
	gamma::Float64 = gamma_initial
	epsilon::Float64 = 10^gamma
	trial_count::SharedArray{Int64, 1} = [0]
	t_start_iteration::Int64 = 0
	number_of_iterations::Int64 = 10
	for current_iteration = 1:number_of_iterations
		t_start_iteration = convert(Int64, time_ns())
		trial_count[1] = 0	
		gamma = gamma - delta_gamma
		epsilon = 10^gamma
		cum_w = cumsum(w)
		@sync @parallel for current_abc_sample = 1:number_of_abc_samples
			idx = rand_weighted_index(cum_w)
			
			m_prim = m[idx, :]
			s_prim = s[idx, :]
			c_prim = c[idx, :]
			az_prim = az[idx]
			
			m_bis = zeros(number_of_components)
			s_bis = zeros(number_of_components)
			c_bis = zeros(number_of_components)
			az_bis = 0.0

			dist_bis = Inf
			
			while dist_bis > epsilon 
				for current_component = 1:number_of_components
					m_bis[current_component] = displace(m_prim[current_component], tau_m[current_component], lb_m, ub_m)
					if distribution_class != "discrete"
						s_bis[current_component] = displace(s_prim[current_component], tau_s[current_component], lb_s, ub_s)
					end
					c_bis[current_component] = displace(c_prim[current_component], tau_c[current_component], lb_c, ub_c)
				end
				
				(n_K_sim, n_DE_sim) = simulate(distribution_class, 
											m_bis,
											s_bis, 
											c_bis, 
											ax, 
											ay, 
											az_bis, 
											Lx, 
											Ly, 
											Lz, 
											number_of_frames, 
											deltat, 
											kmin,
											number_of_de_bins, 
											ub_de)
				
				dist_bis = distance(n_K_real, n_DE_real, n_K_sim, n_DE_sim, d_de)
				
				trial_count[1] = trial_count[1] + 1
			end
			for current_component = 1:number_of_components
				m_star[current_abc_sample, current_component] = m_bis
				if distribution_class != "discrete"
					s_star[current_abc_sample, current_component] = s_bis
				end
				c_star[current_abc_sample, current_component] = c_bis
			end
			az_star[current_abc_sample] = az_bis
			dist_star[current_abc_sample] = dist_bis
		end
		
		println((current_iteration, trial_count[1], gamma, delta_gamma))
		
#		for current_abc_sample = 1:number_of_abc_samples
#			w_star[current_abc_sample] = 0.0
#			for i = 1:number_of_abc_samples
#				w_star[current_abc_sample] = w_star[current_abc_sample] + w[i] * normpdf(m_star[current_abc_sample] - m[i], 0.0, tau_m) * normpdf(s_star[current_abc_sample] - s[i], 0.0, tau_s) * normpdf(c_star[current_abc_sample] - c[i], 0.0, tau_c) * normpdf(az_star[current_abc_sample] - az[i], 0.0, tau_az)
#			end
#			w_star[current_abc_sample] = 1.0 / w_star[current_abc_sample]
#		end
#		w = w_star / sum(w_star)
		w = 1 ./ dist_star.^2
		w = w / sum(w)
		
		m = convert(Array{Float64, 2}, m_star)
		if distribution_class != "discrete"
			s = convert(Array{Float64, 2}, s_star)
		end
		c = convert(Array{Float64, 2}, c_star)
		az = convert(Array{Float64, 1}, az_star)
		dist = convert(Array{Float64, 1}, dist_star)
		
		for current_component = 1:number_of_components
			tau_m[current_component] = sqrt( 2.0 * var(m[:, current_component], corrected = false) )
			tau_s[current_component] = sqrt( 2.0 * var(s[:, current_component], corrected = false) )
			tau_c[current_component] = sqrt( 2.0 * var(c[:, current_component], corrected = false) )
		end
		tau_az = sqrt( 2.0 * var(az, corrected = false) )
	end
	
	#t_exec::Int64 = convert(Int64, time_ns()) - t_start
	#println(t_exec/1e9)
	return (m, s, c, az, dist, w, epsilon)
end