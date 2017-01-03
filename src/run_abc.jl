function generate_abc_samples()
	const T_START::Int64 = time_ns()
	const SEED::Int64 = time_ns()
	srand(SEED)
	
	# Load and parse data.
	data_file = "dna1.dat"
	
	fileStreamIn = open(data_file)
	data = readlines(fileStreamIn)
	nPixelsX::Int64 = parse(Int64,data[1])
	nPixelsY::Int64 = parse(Int64,data[2])
	deltax::Float64 = parse(Float64,data[3])
	deltat::Float64 = parse(Float64,data[4])
	nIncrements::Int64 = parse(Int64,data[5])
	r2::Float64 = parse(Float64,data[6])
	n_frames_str_array = split(data[7], ",")
	n_frames_str_array[end] = n_frames_str_array[end][1:end-1]
	nVideos = length(n_frames_str_array)
	nFrames::Array{Int64,1} = zeros(nVideos)
	for currentVideo = 1:nVideos
		nFrames[currentVideo] = parse(Int64,n_frames_str_array[currentVideo])
	end
	
	kmin::Int64 = 3
	kmax::Int64 = maximum(nFrames)
	
	hist_str_array = split(data[8], ",")
	hist_str_array[end] = hist_str_array[end][1:end-1]
	H::Array{Int64,1} = zeros(kmax)
	for i = 1:kmax
		H[i] = parse(Int64,hist_str_array[i])
	end
		
	ax::Float64 = nPixelsX * deltax
	ay::Float64 = nPixelsY * deltax
	A::Float64 = 100.0

	# Valid for all our data sets, refinedment of the bounds is recommended for faster inference.
	Dlim::Array{Float64,1} = [0.0,2.0]
	alim::Array{Float64,1} = [0.0,3.0]
	clim::Array{Float64,1} = [0.0,1.5e10] 
	
	# Generate a set of 100000 samples (draws) and save results directly to disk.
	nSamples::Int64 = 100000
	
	Hsim::Array{Int64,1} = zeros(size(H))
		
	a::Float64 = 0.0
	c::Float64 = 0.0
	D::Float64 = 0.0
	
	count::Int64 = 0
	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	sigma::Float64 = 0.0
	xlo::Float64 = 0.0
	xhi::Float64 = 0.0
	ylo::Float64 = 0.0
	yhi::Float64 = 0.0
	zlo::Float64 = 0.0
	zhi::Float64 = 0.0
	
	fileNameRes = join(["abc_sample_", data_file[1:end-4], "_", string(SEED), ".dat"])  
	fileStreamRes = open(fileNameRes,"w")
	
	for currentSample = 1:nSamples
		println(currentSample)
		
		D = 1.0/randgamma(nIncrements-1.0,(nIncrements*r2)/(4.0*deltat))	
		D = max(D,Dlim[1])
		D = min(D,Dlim[2])
		sigma = sqrt(2.0*D*deltat)
		
		a = alim[1] + (alim[2]-alim[1]) * rand()
		c = clim[1] + (clim[2]-clim[1]) * rand()

		Hsim = 0 * Hsim # reset histogram
		n = c*A*A*A / 1e12 # Normalizing to unit part/ml.
		nParticles = convert(Int64,ceil(n))
		A_corr = (nParticles/n)^(1.0/3.0) * A

		count = 0
		xlo = 0.5*(A_corr - ax)
		xhi = 0.5*(A_corr + ax)
		ylo = 0.5*(A_corr - ay)
		yhi = 0.5*(A_corr + ay)
		zlo = 0.5*(A_corr - a)
		zhi = 0.5*(A_corr + a)
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
		criterion1 = norm(H[kmin:end]-Hsim[kmin:end])
		criterion2 = norm(cumsum(H[kmin:end])-cumsum(Hsim[kmin:end]))
		write(fileStreamRes, D, a, c, criterion1, criterion2)
	end
	close(fileStreamRes)

	nothing
end

function randgamma(alpha,beta)
	if alpha >= 1.0
        d = alpha - 0.333333333333333
        c = 0.333333333333333/sqrt(d)
        while true
            x = randn()
            v = 1.0 + c*x
            while v <= 0.0
                x = randn()
                v = 1.0 + c*x
            end
            v = v*v*v
            u = rand()
            xsq = x*x
            if u < 1.0 -.0331*xsq*xsq || log(u) < 0.5*xsq + d*(1.0 - v + log(v))
                return d*v/beta
            end
        end
    else
        g = randgamma(alpha+1.0, 1.0)
        w = rand()
        return g * w^(1.0/alpha) / beta
    end
end

generate_abc_samples()