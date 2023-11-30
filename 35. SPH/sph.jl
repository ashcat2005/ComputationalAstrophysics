
using Random: randn, seed!
using SpecialFunctions: gamma
using NPZ: npzwrite

function W(x, y, z, h)
	"""
    3D Gausssian Smoothing kernel
	x : vector/matrix of x positions
	y : vector/matrix of y positions
	z : vector/matrix of z positions
	h : smoothing length
	w : evaluated smoothing function
	"""
	r = sqrt.(x.^2 + y.^2 + z.^2)
	w = (1. / (h*sqrt(pi)))^3 * exp.(-(r/h).^2)
	return w
end


function gradW(x, y, z, h)
    """
    Gradient of the 3D Gausssian Smoothing kernel
    x : vector/matrix of x positions
	y : vector/matrix of y positions
	z : vector/matrix of z positions
	h : smoothing length
    wx, wy, wz : evaluated gradient
    """
    r = sqrt.(x.^2 + y.^2 + z.^2)
    dw = -2 * exp.(-(r/h).^2) / (h^5 * (pi)^(3/2))
    dWx = dw.*x
    dWy = dw.*y
    dWz = dw.*z
    return dWx, dWy, dWz
end


function PairwiseSeparations(ri, rj)
	"""
	Get pairwise separations between 2 sets of coordinates
	ri : M x 3 matrix of positions
	rj : N x 3 matrix of positions
	dx, dy, dz : M x N matrices of separations
	"""
	
	M = size(ri,1)
	N = size(rj,1)

	# positions ri = (x,y,z)
	rix = reshape(ri[:,1], (M,1))
	riy = reshape(ri[:,2], (M,1))
	riz = reshape(ri[:,3], (M,1))

	# other set of points positions rj = (x,y,z)
	rjx = reshape(rj[:,1], (1,N))
	rjy = reshape(rj[:,2], (1,N))
	rjz = reshape(rj[:,3], (1,N))
	
	# matrices that store all pairwise particle separations: r_i - r_j
	dx = rix .- rjx
	dy = riy .- rjy
	dz = riz .- rjz
	
	return dx, dy, dz
end


function Density( r, pos, m, h )
	"""
	Get Density at sampling loctions from SPH particle distribution
	r   : M x 3 matrix of sampling locations
	pos : N x 3 matrix of SPH particle positions
	m   : particle mass
	h   : smoothing length
	rho : M x 1 vector of accelerations
	"""
	M = size(r, 1)
	dx, dy, dz = PairwiseSeparations( r, pos )
	rho = reshape(sum( m * W(dx, dy, dz, h), dims=2 ), (M,1))
	return rho
end


function Pressure(rho, k, n)
	"""
	Equation of State
	rho : vector of densities
	k   : equation of state constant
	n   : polytropic index
	P   : pressure
	"""
	P = k * rho.^(1 + 1/n)
	return P
end


function gravForce(lmbda, r)
	"""
	Simple Gravitational Force Model
	representing a constant gravitational potential
	"""
	return -lmbda.*r
end


function viscosForce(nu, v)
	"""
	Viscosity force
	"""
	return -nu.*v
end


function Acceleration( pos, vel, m, h, k, n, lmbda, nu )
	"""
	Calculate the acceleration on each SPH particle
	[ Euler Equation ]
	pos   : N x 3 matrix of positions
	vel   : N x 3 matrix of velocities
	m     : particle mass
	h     : smoothing length
	k     : equation of state constant
	n     : polytropic index
	lmbda : external force constant
	nu    : viscosity
	a     : N x 3 matrix of accelerations
	"""
	N = size(pos,1)
	
	# Calculate densities at the position of the particles
	rho = Density( pos, pos, m, h )
	
	# Get the pressures
	P = Pressure(rho, k, n)
	
	# Get pairwise distances and gradients
	dx, dy, dz = PairwiseSeparations( pos, pos )
	
	dWx, dWy, dWz = gradW( dx, dy, dz, h )
	
	# Add Pressure contribution to accelerations
	ax = - reshape( sum( m * ( P./rho.^2 .+ transpose(P)./transpose(rho).^2  ) .* dWx, dims=2), (N,1))
	ay = - reshape( sum( m * ( P./rho.^2 .+ transpose(P)./transpose(rho).^2  ) .* dWy, dims=2), (N,1))
	az = - reshape( sum( m * ( P./rho.^2 .+ transpose(P)./transpose(rho).^2  ) .* dWz, dims=2), (N,1))
	
	# pack together the acceleration components
	a = hcat(ax,ay,az)

	# Add gravtiational force
	a += gravForce(lmbda, pos)
	
	# Add viscosity force
	a += viscosForce(nu,vel)
	
	return a
end




function main()
	""" 
	N-body simulation 
	"""
	
	# Simulation parameters
	N     = 1000    # Number of particles
	t     = 0      # current time of the simulation
	tEnd  = 12     # time at which simulation ends
	dt    = 0.04   # timestep
	M     = 2      # star mass
	R     = 0.75   # star radius
	h     = 0.1    # smoothing length
	k     = 0.1    # equation of state constant
	n     = 1      # polytropic index
	nu    = 2      # damping
	
	# Generate Initial Conditions
	seed!(413)            # set the random number generator seed
	
	lmbda = 2*k*(1+n)*pi^(-3/(2*n)) * (M*gamma(5/2+n)/R^3/gamma(1+n))^(1/n) / R^2  # ~ 2.01
	m     = M/N                    # single particle mass
	pos   = randn(Float64, (N,3))  # randomly selected positions and velocities
	vel   = zeros(size(pos))
	#println('lambda = ', lmbda)
	#println('m = ', m)

	# calculate initial gravitational accelerations
	acc = Acceleration( pos, vel, m, h, k, n, lmbda, nu )
	
	# number of timesteps
	Nt = floor(Int,ceil(tEnd/dt))
	rr = zeros(100,3)
	rlin = LinRange(0,1,100)
	rr[:,1] = rlin
	rho_analytic = lmbda/(4*k) .* (R^2 .- rlin.^2)
	
	data = zeros(Nt, N, 3)
	rho_data = zeros(Nt, N)
	density_data = zeros(Nt, 100)
	println("Starting Simulation")
	println()
	# Simulation Main Loop
	for i in 1:Nt
		# (1/2) kick
		vel += acc .* dt/2
		# drift
		pos += vel .* dt
		# update accelerations
		acc = Acceleration( pos, vel, m, h, k, n, lmbda, nu )
		# (1/2) kick
		vel += acc .* dt/2
		
		# update time
		t += dt
		
		# get density for plottiny
		rho = Density( pos, pos, m, h )

		data[i,:,:] = pos
		rho_data[i,:] = reshape(rho, N)
		density_data[i,:] = reshape(Density( rr, pos, m, h ), 100)
		print("\u1b[1F")
		print("Step # ", i)
		print("\u1b[0K")
		println()
	end
	return data, rho_data, density_data, rho_analytic
end

@time data, rho_data, density_data, rho_analytic = main()

npzwrite("data/data.npz", data) 
npzwrite("data/rho_data.npz", rho_data)
npzwrite("data/density_data.npz", density_data)
npzwrite("data/rho_analytic.npz", rho_analytic)
