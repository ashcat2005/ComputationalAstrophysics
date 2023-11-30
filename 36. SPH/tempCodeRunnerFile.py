# (1/2) kick
		vel += acc * dt/2
		# drift
		pos += vel * dt
		# update accelerations
		acc = Acceleration( pos, vel, m, h, k, n, lmbda, nu )
		# (1/2) kick
		vel += acc * dt/2