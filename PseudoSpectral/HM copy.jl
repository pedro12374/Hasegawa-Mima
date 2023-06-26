using PyPlot, Random, FFTW
function main()
	nx, ny = 64, 64 # number of grid points
	Lx, Ly = 2.0*pi, 2.0*pi      # size of the domain in each direction

	# constructing the physical grid (x,y)
	dx, dy = Lx/nx, Ly/ny
	 x = 0:dx:Lx-dx
	 y = 0:dy:Ly-dy

	X  = zeros(nx,ny)
	Y  = zeros(nx,ny)
	for j in 1:ny, i in 1:nx
	     X[i, j] = x[i]
	     Y[i, j] = y[j]
	end

	# constructing the wavenumber grid (k,l)
	k  = 2.0*pi/Lx * [0:nx/2; -nx/2+1:-1];
	l  = 2.0*pi/Ly * [0:ny/2; -ny/2+1:-1];

	Kx = zeros(nx, ny)
	Ky = zeros(nx, ny)
	for j in 1:ny, i in 1:nx
	    Kx[i, j] = k[i] #cria o grid em Kx
	    Ky[i, j] = l[j] #cria o grid em Kx
	end
    
    ksq = @. Kx^2 + Ky^2 #faz k^2
	kqt = @. Kx^4 + Ky^4 + 2 * Kx^2 * Ky^2# faz k^4
    invksq = 1 ./(1 .+ksq) #faz 1/(1+k^2)
    
    
    #Definicoes de tempo
    dt = 0.025 # the time step
	tfin = 1000 # the final time for the integration
	nstep  = Int(tfin/dt) + 1;  # the total number of time steps
    t = 0:dt:tfin; 
    time = 0
	nc = Int(nx/32)

    #Cria as condicoes iniciais
    Random.seed!(1234)
	psi0 = randn(nx, ny) #cria as condicoes inicias aleatoriamentes
	psih = fft(psi0) #calcula a fft das condicoes iniciais
	psih[1, 1] = 0 #filtra os divergentes
	psih[@. ksq > (10)^2 ] .= 0 #faz uma filtragem inicial
	psih[@. ksq < (3)^2 ] .= 0
	psi0 = real(ifft(psih)) #recalculo usando as filtragens 
	phih = invksq.*psih           # calculate \hat{phi}
	psi =  real(ifft(psih))       # calculate ζ(x,y,t)
	phix = real(ifft(im*Kx.*phih))     # calculate phix velocity u(x,y,t)
	phiy = real(ifft(im*Ky.*phih))     # calculate phiy velocity v(x,y,t)
	psix =  real(ifft(im*Kx.*psih))    # calculate ∂psi/∂x
	psiy =  real(ifft(im*Ky.*psih))    # calculate ∂psi/∂y

    #Faz o plot das condicoes iniciais
    figure(1)
    phi = real(ifft(phih))
	pcolormesh(X, Y, psi)
	colorbar()
	quiver(X[1:nc:end,1:nc:end], Y[1:nc:end,1:nc:end], -1 .*phix[1:nc:end,1:nc:end], -1 .*phiy[1:nc:end,1:nc:end])
	xlabel(L"x")
	ylabel(L"y")
	title(string("Potencial & E at t=", time))
	axis("square")
	draw();
	savefig("HMI.png")
	
	#parametros
	nu = 3.0e-04      
	vu = 0.0


	fig = figure()

	for j = 2:nstep
	    
	    phih = invksq.*psih               # calculate \hat{ψ}
		phix = real(ifft(im*Kx.*phih))     # calculate phix velocity u(x,y,t)
		phiy = real(ifft(im*Ky.*phih))     # calculate phiy velocity v(x,y,t)
		psix =  real(ifft(im*Kx.*psih))    # calculate ∂psi/∂x
		psiy =  real(ifft(im*Ky.*psih))    # calculate ∂psi/∂y

	    # now we are ready to calculate the r.h.s. of our PDE
	    rhs = -fft(phix.*psiy - phiy.*psix) - nu.*(-ksq.*phih+ksq.*psih)+vu.*fft(phiy)
	    #println(rhs)

	    psih = psih + dt*rhs         # time-step forward using the simplest Euler scheme
		phih = invksq.*psih
	    psi = real(ifft(psih))
		if any(isnan, psi)
			println("Entrou")
			 println(j)
			 #println(phibk)
			 break
		 end

		if j % 200 == 1
			phi = real(ifft(phih))
			time = t[j]
			#phih = -invksq.*psih
			#CSV.write("dados_.csv",DataFrame(psi,:auto),header=false)
			pcolormesh(X, Y, phi)
			#contour(X, Y, phi)
			title(string(L"vorticity at t=", time))
			colorbar()
			quiver(X[1:nc:end,1:nc:end], Y[1:nc:end,1:nc:end], phix[1:nc:end,1:nc:end], phiy[1:nc:end,1:nc:end])
			#xlabel(L"x")
			#ylabel(L"y")
			title(string(L"vorticity e (v,u) at t=", time))
			axis("square")
			#sleep(0.001)
			display(fig)
			clf()
		end
	end
	time = tfin
	phih = invksq.*psih               # calculate \hat{ψ}
	phi = real(ifft(phih))
	phix = real(ifft(im*Kx.*phih))     # calculate phix velocity u(x,y,t)
	phiy = real(ifft(im*Ky.*phih))     # calculate phiy velocity v(x,y,t)
	
	clf()

	pcolormesh(X, Y, phi)
	colorbar()
	quiver(X[1:nc:end,1:nc:end], Y[1:nc:end,1:nc:end], -1 .*phix[1:nc:end,1:nc:end], -1 .*phiy[1:nc:end,1:nc:end])
	xlabel(L"x")
	ylabel(L"y")
	title(string("Potencial & E at t=", time))
	draw()
	axis("square");
	savefig("HMF.png")

end
main()
