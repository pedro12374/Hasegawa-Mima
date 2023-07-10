using PyPlot, Random, FFTW, DelimitedFiles, CSV, DataFrames, ProgressBars, NumericalIntegration

function mindist2(x, val)
    # using enumerate to avoid indexing
    min_i = 0
    min_x = Inf
    for (i, xi) in enumerate(x)
        dist = abs(xi - val)
        if dist < min_x
            min_x = dist
            min_i = i
        end
    end
    return min_i
end

function main()

    nx, ny = 64, 64 # number of grid points
	Lx, Ly = 2.0*pi, 2.0*pi      # size of the domain in each direction
    nc = Int(nx/32)


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
	    Ky[i, j] = l[j] #cria o grid em Ky
	end
    
    ksq = @. Kx^2 + Ky^2 #faz k^2
	kqt = @. Kx^4 + Ky^4 + 2 * Kx^2 * Ky^2# faz k^4
    invksq = 1 ./(1 .+ksq) #faz 1/(1+k^2)
    


    phih = Matrix(CSV.read("Poth.csv", DataFrame, header=false,delim='\t', types=Complex{Float64}))
    phix = real(ifft(im*Kx.*phih))     # calculate phix velocity u(x,y,t)
	phiy = real(ifft(im*Ky.*phih))


    fig, ax = plt.subplots()
	plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
	fig.set_size_inches(20*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """


	phi = real(ifft(phih))
    ax.pcolormesh(X, Y, phi)
	#fig.colorbar()
	ax.quiver(X[1:nc:end,1:nc:end], Y[1:nc:end,1:nc:end], -1 .*phix[1:nc:end,1:nc:end], -1 .*phiy[1:nc:end,1:nc:end])
	ax.set_xlabel(L"x")
	ax.set_ylabel(L"y")
	ax.set_title(string("Potencial & E"))
    ax.axis("square");
	
    dt = 0.1
    tt = 1000
    tc = Int(tt/dt)
    vecposicoes = zeros(tc,0)
    S = zeros(tc,0)
    for j in 1:32:64, u in 1:32:64
    	
        vel = zeros(2,0)
        pos = zeros(2,0)
        vpos = [x[Int(j)],y[Int(u)]]
        vvel = [0,0]
        #pos = hcat(pos,vpos)
        vel = hcat(vel,vvel)
        t = zeros(tc)
        t[1] = 0
        #Random.seed!(123)
        spos = zeros(2,0)
        spos = hcat(spos,vpos)
        yp = vpos[2]+20.0*phix[mindist2(x,vpos[1]),mindist2(y,vpos[2])]*dt
        vpos[2] = yp
        S = zeros(1)
        for i in 2:tc
            npx = vpos[1]- 20.0*phiy[mindist2(x,vpos[1]),mindist2(x,vpos[2])]*dt
            npy = vpos[2]+ 20.0*phix[mindist2(x,vpos[1]),mindist2(y,vpos[2])]*dt
            vvel = [-20.0*phiy[mindist2(x,vpos[1]),mindist2(x,vpos[2])],20.0*phix[mindist2(x,vpos[1]),mindist2(y,vpos[2])]]
            pos = [npx,npy]
            spos = hcat(spos,pos)
            vel = hcat(vel,vvel)
            mpx = pos[1]
            mpy = pos[2]
            nvx =0#  real(ifft(-im*Ky.*phih))[mindist2(x,npx)]	
            nvy = 0# real(ifft(im*Kx.*phih))[mindist2(y,npy)]
            
            vpos = [mod(mpx,2*pi),mod(mpy,2*pi)]
             
            tvel = transpose(vel)
            intv =@. sqrt(tvel[:,1]^2+tvel[:,2]^2) 
            s = integrate(1:i,intv)
            push!(S,s)
            
            #nv = [nvx,nvy]
            #println( real(ifft(im*Ky.*phih))  )
            #v = hcat(v,nv)
            t[i]=dt*i
        end
    
        pplot = transpose(mod.(spos,2*pi))
        #ax.plot(t,S)
        ax.plot(pplot[:,1],pplot[:,2],",")

        vecposicoes = hcat(vecposicoes,S)
    end
    #draw()
    plt.savefig("part.png")
    #=
    plt.close("all")
    #ax.clear()
    fig, ax = plt.subplots()
	plt.tight_layout() # isso aq nsei bem qq faz mas ajuda a deixar menos espaço em branco
	fig.set_size_inches(40*0.393, 20*0.393) # esse fatir 0.393 é p converter polegadas p cm """
    
    for i in 1:2:8
        ax.plot(vecposicoes[:,i])
    end
    plt.savefig("x.png")
    plt.cla()
    #plot(t,mpplot[:,2],"b")
    #draw()
    
    for i in 2:2:8
        ax.plot(vecposicoes[:,i])
    end
    plt.savefig("y.png")
=#
    writedlm("pos.csv",vecposicoes,"\t")
    #writedlm("v.csv",transpose(v),"\t")


end
main()
