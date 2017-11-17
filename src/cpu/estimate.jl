function estimate(model::String,
				number_of_components::Int64,
				Lx::Float64,
				Ly::Float64,
				Lz::Float64,
				kmin::Int64,
				number_of_de_bins::Int64,
				ub_de::Float64,
				lb_m::Float64,
				ub_m::Float64,
				lb_c::Float64,
				ub_c::Float64,
				lb_az::Float64,
				ub_az::Float64,
				number_of_abc_samples::Int64,
				gamma_initial::Float64,
				gamma_adaptive::Bool,
				delta_gamma::Float64,
				weighting_scheme::String,
				ub_average_number_of_trials::Int64,
				ax::Float64,
				ay::Float64,
				number_of_frames::Array{Int64, 1},
				deltat::Float64,
				K::Array{Int64, 1},
				DE::Array{Float64, 1})

	# Distance function histogram bin width.
	d_de::Float64 = ub_de / convert(Float64, number_of_de_bins)

	# Convert data to histogram form.
	kmax::Int64 = maximum(number_of_frames)
	H::Array{Int64, 2} = zeros(kmax, number_of_de_bins)
	idx::Int64 = 0
	for i = 1:length(K)
		idx = convert(Int64, ceil(DE[i] / d_de))
		if 1 <= idx <= number_of_de_bins
			H[K[i], idx] += 1
		end
	end

	# Variables for population parameter values.
	m::Array{Float64, 2} = zeros(number_of_components, number_of_abc_samples)
	c::Array{Float64, 2} = zeros(number_of_components, number_of_abc_samples)
	az::Array{Float64, 2} = zeros(number_of_components, number_of_abc_samples)
	dist::Array{Float64, 1} = zeros(number_of_abc_samples)

	m_star::SharedArray{Float64, 2} = zeros(number_of_components, number_of_abc_samples)
	c_star::SharedArray{Float64, 2} = zeros(number_of_components, number_of_abc_samples)
	az_star::SharedArray{Float64, 2} = zeros(number_of_components, number_of_abc_samples)
	dist_star::SharedArray{Float64, 1} = zeros(number_of_abc_samples)

	# Initialize assuming that gamma = inf so that everything is accepted.
	m = lb_m + (ub_m - lb_m) * rand(number_of_components, number_of_abc_samples)
	c = lb_c + (ub_c - lb_c) * rand(number_of_components, number_of_abc_samples)
	if model == "discrete-fixed-depth"
		az = lb_az + (ub_az - lb_az) * repmat(rand(1, number_of_abc_samples), 2)
	elseif model == "discrete-variable-depth"
		az = lb_az + (ub_az - lb_az) * rand(number_of_components, number_of_abc_samples)
	end

	w::Array{Float64, 1} = ones(number_of_abc_samples) / convert(Float64, number_of_abc_samples)
	cum_w::Array{Float64, 1} = cumsum(w)
	w_star::Array{Float64, 1} = zeros(number_of_abc_samples)
	term::Float64 = 0.0

	# Displacement standard deviations.
	tau_m::Array{Float64, 1} = zeros(number_of_components)
	tau_c::Array{Float64, 1} = zeros(number_of_components)
	tau_az::Array{Float64, 1} = zeros(number_of_components)
	for current_component = 1:number_of_components
		tau_m[current_component] = sqrt( 2.0 * var(m[current_component, :], corrected = false) )
		tau_c[current_component] = sqrt( 2.0 * var(c[current_component, :], corrected = false) )
		tau_az[current_component] = sqrt( 2.0 * var(az[current_component, :], corrected = false) )
	end

	# Adaptive gamma selection
	gamma::Float64 = 0.0
	if gamma_adaptive
		@sync @parallel for current_abc_sample = 1:number_of_abc_samples
			H_sim = simulate(	model,
								m[:, current_abc_sample],
								c[:, current_abc_sample],
								ax,
								ay,
								az[:, current_abc_sample],
								Lx,
								Ly,
								Lz,
								number_of_frames,
								deltat,
								kmin,
								number_of_de_bins,
								ub_de)

			dist_star[current_abc_sample] = distance(H, H_sim)
		end
		gamma = median(dist_star)
	else
		gamma = gamma_initial
	end

	# The rest of the iterations.
	trial_count::SharedArray{Int64, 1} = zeros(number_of_abc_samples)
	is_converged::Bool = false
	while !is_converged
		trial_count = zeros(number_of_abc_samples)
		gamma = gamma - delta_gamma
		cum_w = cumsum(w)
		@sync @parallel for current_abc_sample = 1:number_of_abc_samples
			dist_bis = Inf

			m_bis = m[:, current_abc_sample]
			c_bis = c[:, current_abc_sample]
			az_bis = az[:, current_abc_sample]

			while dist_bis > gamma && mean(trial_count) < convert(Float64, ub_average_number_of_trials)
				idx = rand_weighted_index(cum_w)

				m_prim = m[:, idx]
				c_prim = c[:, idx]
				az_prim = az[:, idx]

				m_bis = zeros(number_of_components)
				c_bis = zeros(number_of_components)
				az_bis = zeros(number_of_components)

				if model == "discrete-fixed-depth"
					for current_component = 1:number_of_components
						m_bis[current_component] = displace(m_prim[current_component], tau_m[current_component], lb_m, ub_m)
						c_bis[current_component] = displace(c_prim[current_component], tau_c[current_component], lb_c, ub_c)
					end
					az_bis[1] = displace(az_prim[1], tau_az[1], lb_az, ub_az)
					az_bis[2:end] = az_bis[1]
				elseif model == "discrete-variable-depth"
					for current_component = 1:number_of_components
						m_bis[current_component] = displace(m_prim[current_component], tau_m[current_component], lb_m, ub_m)
						c_bis[current_component] = displace(c_prim[current_component], tau_c[current_component], lb_c, ub_c)
						az_bis[current_component] = displace(az_prim[current_component], tau_az[current_component], lb_az, ub_az)
					end
				end

				p = sortperm(m_bis) # To avoid label switching problems.
				m_bis = m_bis[p]
				c_bis = c_bis[p]
				az_bis = az_bis[p]

				H_sim = simulate(	model,
									m_bis,
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

				dist_bis = distance(H, H_sim)

				trial_count[current_abc_sample] = trial_count[current_abc_sample] + 1
			end

			m_star[:, current_abc_sample] = m_bis
			c_star[:, current_abc_sample] = c_bis
			az_star[:, current_abc_sample] = az_bis
			dist_star[current_abc_sample] = dist_bis
		end

		if mean(trial_count) >= ub_average_number_of_trials
			is_converged = true
			println("Converged.")
		end

		if !is_converged # Only compute new stuff if we are going to do one more iteration.
			if weighting_scheme == "pmc-standard"
				for current_abc_sample = 1:number_of_abc_samples
					w_star[current_abc_sample] = 0.0
					for i = 1:number_of_abc_samples
						term = 1.0
						for current_component = 1:number_of_components
							term = term * normpdf(m_star[current_component, current_abc_sample] - m[current_component, i], 0.0, tau_m[current_component])
							println((1, normpdf(m_star[current_component, current_abc_sample] - m[current_component, i], 0.0, tau_m[current_component])))
							term = term * normpdf(c_star[current_component, current_abc_sample] - c[current_component, i], 0.0, tau_c[current_component])
							println((2, normpdf(c_star[current_component, current_abc_sample] - c[current_component, i], 0.0, tau_c[current_component])))
							term = term * normpdf(az_star[current_component, current_abc_sample] - az[current_component, i], 0.0, tau_az[current_component])
							println((3, normpdf(az_star[current_component, current_abc_sample] - az[current_component, i], 0.0, tau_az[current_component])))
						end
						w_star[current_abc_sample] = w_star[current_abc_sample] + w[i] * term
					end
					println((1, w_star))
					w_star[current_abc_sample] = 1.0 / w_star[current_abc_sample]
				end
				w = w_star / sum(w_star)
			elseif weighting_scheme == "inverse-distance-squared"
				w = 1 ./ dist_star.^2
				w = w / sum(w)
			end

			m = convert(Array{Float64, 2}, m_star)
			c = convert(Array{Float64, 2}, c_star)
			az = convert(Array{Float64, 2}, az_star)
			dist = convert(Array{Float64, 1}, dist_star)

			for current_component = 1:number_of_components
				tau_m[current_component] = sqrt( 2.0 * var(m[current_component, :], corrected = false) )
				tau_c[current_component] = sqrt( 2.0 * var(c[current_component, :], corrected = false) )
				tau_az[current_component] = sqrt( 2.0 * var(az[current_component, :], corrected = false) )
			end

			if model == "discrete-fixed-depth" || model == "discrete-variable-depth"
				if number_of_components == 1
					println((round(gamma, 2), sum(trial_count), round(mean(trial_count), 2), round(mean(m), 2), round(mean(c), 2), round(mean(az), 2)))
				elseif number_of_components == 2
					println((round(gamma, 2), sum(trial_count), round(mean(trial_count), 2), round(mean(m[1, :]), 2), round(mean(m[2, :]), 2), round(mean(c[1, :]), 2), round(mean(c[2, :]), 2), round(mean(az[1, :]), 2), round(mean(az[2, :]), 2)))
					println((mean(w), std(w), minimum(w), maximum(w)))
				end
			end
		end
	end

	gamma = gamma + delta_gamma # Save the gamma value for the last complete iteration.

	return (m, c, az, dist, w, gamma)
end
