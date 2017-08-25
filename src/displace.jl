function displace(x_prim::Float64, tau_x::Float64, lb_x::Float64, ub_x::Float64)

	delta_x::Float64 = tau_x * randn()
	x_bis::Float64 = x_prim + delta_x
	while !( lb_x < x_bis < ub_x )
		if x_bis < lb_x
			x_bis = x_bis + 2.0 * ( lb_x - x_bis )
		end
		if x_bis > ub_x
			x_bis = x_bis - 2.0 * ( x_bis - ub_x )
		end
	end
	
	return x_bis
end