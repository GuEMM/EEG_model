using Random, Distributions, DelimitedFiles

function Listof_KOut(C,idn,L,Ne)

    cy,cx = Tuple.(findall(x->x==idn, C))[1]

    Rin = 3

    cxi = (cx-Rin)*(cx-Rin>0) + (L+(cx-Rin))*(cx-Rin<0) + (cx-Rin==0)*L
    cxf = (cx+Rin)*(cx+Rin<L) + ((cx+Rin)-L)*(cx+Rin>L) + (cx+Rin==L)

    if (cxi>cx)||(cxf<cx)
        cc = C[:,[cxi:L; 1:cxf]]
    else
        cc = C[:,cxi:cxf]
    end

    cyi = (cy-Rin)*(cy-Rin>0) + (L+(cy-Rin))*(cy-Rin<0) + (cy-Rin==0)*L
    cyf = (cy+Rin)*(cy+Rin<L) + ((cy+Rin)-L)*(cy+Rin>L) + (cy+Rin==L)

    if (cyi>cy)||(cyf<cy)
        cc = cc[[cyi:L; 1:cyf],:]
    else
        cc = cc[cyi:cyf,:]
    end

    mask = ones(Bool,2*Rin+1,2*Rin+1)

    mask[1,1] = 0
    mask[1,2*Rin+1] = 0
    mask[2*Rin+1,1] = 0
    mask[2*Rin+1,2*Rin+1] = 0

    NeiOut = cc.*mask

    return Int.(filter(x -> (x>0)&(x <= Ne), NeiOut))
end

function Listof_KIn(C,idn,L,Ne)

    cy,cx = Tuple.(findall(x->x==idn, C))[1]

    Rin = 5

    cxi = (cx-Rin)*(cx-Rin>0) + (L+(cx-Rin))*(cx-Rin<0) + (cx-Rin==0)
    cxf = (cx+Rin)*(cx+Rin<L) + ((cx+Rin)-L)*(cx+Rin>L) + (cx+Rin==L)*L

    if (cxi>cx)||(cxf<cx)
        cc = C[:,[cxi:L; 1:cxf]]
    else
        cc = C[:,cxi:cxf]
    end

    cyi = (cy-Rin)*(cy-Rin>0) + (L+(cy-Rin))*(cy-Rin<0) + (cy-Rin==0)
    cyf = (cy+Rin)*(cy+Rin<L) + ((cy+Rin)-L)*(cy+Rin>L) + (cy+Rin==L)*L

    if (cyi>cy)||(cyf<cy)
        cc = cc[[cyi:L; 1:cyf],:]
    else
        cc = cc[cyi:cyf,:]
    end

    mask = ones(Int,2*Rin+1,2*Rin+1)

    mask[1,1] = 0
    mask[1,2*Rin+1] = 0
    mask[2*Rin+1,1] = 0
    mask[2*Rin+1,2*Rin+1] = 0

    NeiIn= cc.*mask

    return Int.(filter(x -> (x>0)&(x <= Ne), NeiIn))
end

function Vdyn(V,SVin,Vnoise,tau1h,tau2h,Vres)

    tauh = tau1h.*(V.>=Vres) .+ tau2h.*(V.<Vres)

    return (-V .+ SVin .+ Vnoise)./tauh
end

#Threshold dyn
function FVthr(Vth0,Vsat,t,tsp,ta,kaph)
    Vth = Vsat.*(t.>=tsp).*(t.<(tsp.+ta)) .+ (Vth0.+(Vsat.-Vth0).*exp.(-kaph.*(t.-(tsp.+ta)))).*((t.>=(tsp.+ta)).+(t.<tsp))
    return Vth
end

function tspK_act(tsp,tspK,KK,N)
    #Recover information of the firing time of each neighbor
    
    for i in 1:N

        idx = KK[i,:]
        
        l = sum(idx.!=0)
        
        id = idx[1:l]

        tspK[i,1:l] = tsp[id] 

    end

    return tspK
end

function BuildNet(ce)

    ci = Int(ce/2)

    ns = ce-1

    L = ce+ns

    C = zeros(Int,L,L)

    paso = 2

    nn = 1

    for i in 1:ce

        for j in 1:ce

            C[(i-1)*paso+1,(j-1)*paso+1] = nn

            nn += 1

        end
    end

    paso = 4

    for i in 1:ci

        for j in 1:ci

            C[(i-1)*paso+2,(j-1)*paso+2] = nn

            nn += 1 

        end

    end

    CC = zeros(Int,L+1,L+1)

    CC[2:L+1,2:L+1] = C

    C = CC;

    L = ce+ns+1

    Ne = ce^2
    Ni = ci^2

    N = Ne+Ni

    #List of E pre-synaptic neighbors of I neurons 
    Kin = 32
    Kien = zeros(Int,Ni,32) 

    #List of E post-synaptic neighbors of I neurons
    Kout = 12
    KKout = zeros(Int,Ni,Kout) 

    for i in 1:Ni

        idn = Ne + i

        KKout[i,:] .= Listof_KOut(C,idn,L,Ne)
        Kien[i,:] .= Listof_KIn(C,idn,L,Ne)

    end

    kei = length(first.(Tuple.(findall(e -> e==1, KKout))))

    #List of I pre-synaptic neighbors of E neurons
    Kein = zeros(Int,Ne,kei)

    for i in 1:Ne

        Kein[i,:] = Ne .+ first.(Tuple.(findall(e -> e==i, KKout)))

    end

    KK = zeros(Int,N,Kin)

    #Full list of neighbors
    KK[1:Ne,1:3] = Kein
    KK[Ne+1:N,1:Kin] = Kien;
    
    return KK,C,KKout
end

function LIF_Sim_TauMu_bymodules(N,Ne,A,KK,KIout,V0d,taur,mu,filename="Test.txt",Tt=100,dts=100,h=4/100,Transient=100,tau1=16,tau2=26,Vth0=6,Vsat=90,Vmin=-20,Vres=0,V0dext=5.48,U=0.5,ta=100,tmax=100,kap=2,Nn=100)
    
    filename1 = filename*"_States.txt"
    filename2 = filename*"_Vm.txt"

    Ni = N-Ne
    
    sampT = 10000
    
    Tt = Tt+Transient
    
    V = zeros(Float64,N)

    Vin = zeros(Float64,N)

    Vnoise = zeros(Float64,N)

    tsp = ones(Int,N)

    tnoiL = zeros(Int,Ne,Nn);
    SnoiL = zeros(Int,Ne,tmax);

    tau2h = tau2/h

    SVin = zeros(Float64,N);

    nn = 1
    
    V0h = -4. *V0d

    xstp = ones(Float64,N);
    x = ones(Float64,N);
    
    tau1h = tau1/h
    tau2h = tau2/h
    taurh = taur/h

    kaph = kap*h

    Vnor = ones(Float64,N)
    Vth = ones(Float64,N).*Vth0

    lambda = mu/(tmax*Nn)
    
    IIdn = [9,13,25,37,41]

    idn = Ne .+ IIdn

    IIV = [vcat(collect(Set(KK[KIout[idn[i]-Ne,:],1:3])),KK[idn[i],:]) for i in 1:length(idn)]
    
    IIdn = [25]
    
    idn = Ne.+IIdn
    
    IIS = [vcat(collect(Set(KK[KIout[idn[i]-ce^2,:],1:3])),KIout[idn[i]-Ne,:]) for i in 1:length(idn)]
    
    #IIS = [vcat(collect(Set(KK[KIout[idn[i]-ce^2,[4,5,8,9]],1:3])),KIout[idn[i]-Ne,:]) for i in 1:length(idn)]
    
    #IIS = collect(Iterators.flatten(IIS))
    
    #Ln = size(IIS)[1]
    
    Snei = zeros(Int,N)
    
    VneiE = zeros(Float32,size(IIV)[1])
    VneiI = zeros(Float32,size(IIV)[1])
    
    SI = 0
    
    io1 = open(filename1, "w")
    
    close(io1)
    
    io2 = open(filename2, "w")

    nclos = 0

    for t in 1:Tt-1

        Vnor[Ne+1:end] = ((Vsat.-V[Ne+1:end])./Vsat) 
        
        Vnor[1:Ne] = ((Vmin.-V[1:Ne])./Vmin)
        
        SVin = (A*Vin).*Vnor

        #############################################################
        
        SnoiL[:,nn] = sum(rand(Ne,Nn) .< lambda,dims=2)

        tnoiL[:,nn] = t.*(SnoiL[:,nn].>0)

        Vnoise[1:Ne] = ((Vsat.-V[1:Ne])./Vsat).*V0dext.*sum(SnoiL.*((t.>=tnoiL).-(t.>=(tnoiL.+tmax))),dims=2)
        
        ##############################################################

        tauh = tau1h.*(V.>=Vres) .+ tau2h.*(V.<Vres)

        V += (-V .+ SVin .+ Vnoise)./tauh

        ttp = (V .> Vth)
        
        tsp = t.*ttp + tsp.*(1 .-ttp)

        Vin[1:Ne] = (V0d*U.*xstp[1:Ne].*ttp[1:Ne] .+ Vin[1:Ne].*(1. .-ttp[1:Ne])).*(t.<tsp[1:Ne].+tmax)
        
        Vin[Ne+1:end] += V0h*U.*xstp[Ne+1:end].*ttp[Ne+1:end] .- (1 .- ttp[Ne+1:end]).*Vin[Ne+1:end]./tau2h
        
        if taur>0
            
            xstp = x.*ttp .+ xstp.*(1. .-ttp)

            x += ((1. .- x)./taurh .- U.*x.*ttp);
        
        end
        
        nn = (nn+1)*(nn<tmax) + (nn>=tmax)
	
        Vth = Vsat.*ttp .+ Vth.*(1 .-ttp) .- (Vth .- Vth0).*kaph.*(t.> tsp.+ta)
        
        if (t>=Transient)
            
            if t<=(Transient+10000000)
            
                VneiI[:] = [mean(V[IIV[i][1:9]]) for i in 1:size(IIV)[1]]
                VneiE[:] = [mean(V[IIV[i][10:end]]) for i in 1:size(IIV)[1]]
            
                writedlm(io2, [vcat(VneiI,VneiE)],',')
            
            else
                
		if nclos==0
	                close(io2)
        	        println("Vm_modules_Finished...#")
			nclos += 1
    		end
	    end
            
            if (t%dts == 0)
                
                Snei[:] =  Int8.(t .< (tsp .+ ta))
                
                io1 = open(filename1, "a")
                
                writedlm(io1, [Snei],',')
                
                close(io1)

            end
	
            if (t%100000 == 0)

		GC.gc(true)

	   end
    	end
    end
    
    return nothing

end

function LIF_Sim_TauMu_bymodules_Inoise(N,Ne,A,IIV,V0d,taur,mu,nu,filename="Test.txt",Tt=100,dts=100,h=4/100,Transient=100,tau1=16,tau2=26,Vth0=6,Vsat=90,Vmin=-20,Vres=0,V0dext=5.48,U=0.5,ta=100,tmax=100,kap=2,Nn=100)
    
    filename1 = filename*"_States.txt"
    filename2 = filename*"_Vm.txt"

    Ni = N-Ne
    
    sampT = 10000
    
    Tt = Tt+Transient
    
    V = zeros(Float64,N)

    Vin = zeros(Float64,N)

    Vnoise = zeros(Float64,N)

    tsp = ones(Int,N)
    
    tnoiL = zeros(Int,N,Nn);
    SnoiL = zeros(Int,N,tmax);

    tau2h = tau2/h

    SVin = zeros(Float64,N);

    nn = 1
    
    V0h = -4. *V0d

    xstp = ones(Float64,N);
    x = ones(Float64,N);
    
    tau1h = tau1/h
    tau2h = tau2/h
    taurh = taur/h

    kaph = kap*h

    Vnor = ones(Float64,N)
    Vth = ones(Float64,N).*Vth0

    lambdaE = Float32(mu/(tmax*Nn))
    lambdaI = Float32(nu/(tmax*Nn))
 
    Snei = zeros(Int,N)
    
    VneiE = zeros(Float32,size(IIV)[1])
    VneiI = zeros(Float32,size(IIV)[1])
    
    SI = 0
    
    io1 = open(filename1, "w")
    
    close(io1)
    
    io2 = open(filename2, "w")

    nclos=0
        
    for t in 1:Tt-1

        Vnor[Ne+1:end] = ((Vsat.-V[Ne+1:end])./Vsat) 
        
        Vnor[1:Ne] = ((Vmin.-V[1:Ne])./Vmin)
        
        SVin[1:end] = (A*Vin).*Vnor

        #############################################################
        
        SnoiL[1:Ne,nn] = reduce(+,rand(Ne,Nn) .< lambdaE ,dims=2)
        SnoiL[Ne+1:end,nn] = reduce(+,rand(Ni,Nn) .< lambdaI ,dims=2)
	
	tnoiL[:,nn] = t.*(SnoiL[:,nn].>0)
	
        Vnoise[:] = Float64.(reduce(+,SnoiL.*((t.>=tnoiL).-(t.>=(tnoiL.+tmax))),dims=2))

        ##############################################################

        tauh = tau1h.*(V.>=Vres) .+ tau2h.*(V.<Vres)

        V += (-V .+ SVin .+ ((Vsat.-V)./Vsat).*V0dext.*Vnoise)./tauh

        ttp = (V .> Vth)

        tsp = t.*ttp + tsp.*(1 .-ttp)

        Vin[1:Ne] = (V0d*U.*xstp[1:Ne].*ttp[1:Ne] .+ Vin[1:Ne].*(1. .-ttp[1:Ne])).*(t.<tsp[1:Ne].+tmax)

        Vin[Ne+1:end] += V0h*U.*xstp[Ne+1:end].*ttp[Ne+1:end] .- (1 .- ttp[Ne+1:end]).*Vin[Ne+1:end]./tau2h

        if taur>0
            
            xstp = x.*ttp .+ xstp.*(1. .-ttp)

            x += ((1. .- x)./taurh .- U.*x.*ttp);
        
        end
        
        nn = (nn+1)*(nn<tmax) + (nn>=tmax)
	
        Vth = Vsat.*ttp .+ Vth.*(1 .-ttp) .- (Vth .- Vth0).*kaph.*(t.> tsp.+ta)
        
        if (t>=Transient)
            
            if t<=(Transient+10000000)
            
                VneiI[:] = [mean(V[IIV[i][1:9]]) for i in 1:size(IIV)[1]]
                VneiE[:] = [mean(V[IIV[i][10:end]]) for i in 1:size(IIV)[1]]
            
                writedlm(io2, [vcat(VneiI,VneiE)],',')
            
            else
                
		if nclos==0
	                close(io2)
        	        println("Vm_modules_Finished...#")
			nclos += 1
    		end
	    end
            
            if (t%dts == 0)
                
                Snei[:] =  Int8.(t .< (tsp .+ ta))
                
                io1 = open(filename1, "a")
                
                writedlm(io1, [Snei],',')
                
                close(io1)

            end
	
            if (t%100000 == 0)

		GC.gc(true)

	   end
    	end
    end
    
    return nothing

end

function LIF_Sim_TauMu_States(N,Ne,A,IIV,V0d,taur,mu,nu,filename="Test.txt",Tt=100,dts=100,h=4/100,Transient=100,tau1=16,tau2=26,Vth0=6,Vsat=90,Vmin=-20,Vres=0,V0dext=5.48,U=0.5,ta=100,tmax=100,kap=2,Nn=100)
    
    filename1 = filename*"_States.txt"
    filename2 = filename*"_Vm.txt"

    Ni = N-Ne
    
    Tt = Tt+Transient
    
    V = zeros(Float64,N)

    Vin = zeros(Float64,N)

    Vnoise = zeros(Float64,N)

    tsp = ones(Int64,N)
    ttp = ones(Int64,N)

    tnoiL = zeros(Int64,N,Nn);
    SnoiL = zeros(Int64,N,tmax);

    tau2h = tau2/h

    SVin = zeros(Float64,N);

    nn = 1
    
    V0h = -4. *V0d

    xstp = ones(Float64,N);
    x = ones(Float64,N);
    
    tau1h = tau1/h
    tau2h = tau2/h
    taurh = taur/h

    kaph = kap*h

    Vnor = ones(Float64,N)
    Vth = ones(Float64,N).*Vth0

    tauh = ones(Float64,N)

    lambdaE = mu/(tmax*Nn)
    
    lambdaI = nu/(tmax*Nn)
    
    VneiE = zeros(Float32,size(IIV)[1])
    VneiI = zeros(Float32,size(IIV)[1])
 
    Snei = zeros(Int,N)
    
    io1 = open(filename1, "w")
    
    close(io1)
    
    io2 = open(filename2, "w")
    
    nclos=0

    for t in 1:Tt-1

        Vnor[Ne+1:end] = ((Vsat.-V[Ne+1:end])./Vsat) 
        
        Vnor[1:Ne] = ((Vmin.-V[1:Ne])./Vmin)
        
        SVin[1:end] = (A*Vin).*Vnor

        #############################################################
        
        SnoiL[1:Ne,nn] = reduce(+,rand(Ne,Nn) .< lambdaE ,dims=2)
	
	SnoiL[Ne+1:end,nn] = reduce(+,rand(Ni,Nn) .< lambdaI ,dims=2)

        tnoiL[:,nn] = t.*(SnoiL[:,nn].>0)

        Vnoise[:] = reduce(+,SnoiL.*((t.>=tnoiL).-(t.>=(tnoiL.+tmax))),dims=2)
        
        ##############################################################

        tauh[1:end] = tau1h.*(V.>=Vres) .+ tau2h.*(V.<Vres)
	
        V += (-V .+ SVin .+ ((Vsat.-V)./Vsat).*V0dext.*Vnoise)./tauh

        ttp[1:end] = (V .> Vth)
        
        tsp[1:end] = t.*ttp + tsp.*(1 .-ttp)

        Vin[1:Ne] = (V0d*U.*xstp[1:Ne].*ttp[1:Ne] .+ Vin[1:Ne].*(1. .-ttp[1:Ne])).*(t.<tsp[1:Ne].+tmax)
        
        Vin[Ne+1:end] += V0h*U.*xstp[Ne+1:end].*ttp[Ne+1:end] .- (1 .- ttp[Ne+1:end]).*Vin[Ne+1:end]./tau2h
        
        if taur>0
            
            xstp[1:end] = x.*ttp .+ xstp.*(1. .-ttp)

            x += ((1. .- x)./taurh .- U.*x.*ttp);
        
        end

        nn = (nn+1)*(nn<tmax) + (nn>=tmax)

        Vth[1:end] = Vsat.*ttp .+ Vth.*(1 .-ttp) .- (Vth .- Vth0).*kaph.*(t.> tsp.+ta)

        if (t>=Transient)
            
            if t<=(Transient+10000000)
            
                VneiI[:] = [mean(V[IIV[i][1:9]]) for i in 1:size(IIV)[1]]
                VneiE[:] = [mean(V[IIV[i][10:end]]) for i in 1:size(IIV)[1]]
            
                writedlm(io2, [vcat(VneiI,VneiE)],',')
            
            else
                
		if nclos==0
	                close(io2)
        	        println("Vm_modules_Finished...#")
			nclos += 1
    		end
	    end
            
            if (t%dts == 0)
                
                Snei[:] =  Int8.(t .< (tsp .+ ta))
                
                io1 = open(filename1, "a")
                
                writedlm(io1, [Snei],',')
                
                close(io1)

            end
	
            if (t%100000 == 0)

		GC.gc(true)

	   end
    	end
    	
    end
    
    return nothing
end


function LIF_Sim_TauMu_SinCur(N,Ne,A,IIV,V0d,taur,mu,nu,filename="Test.txt",Tt=100,dts=100,h=4/100,Transient=100,tau1=16,tau2=26,Vth0=6,Vsat=90,Vmin=-20,Vres=0,V0dext=5.48,U=0.5,ta=100,tmax=100,kap=2,Nn=100)
    
    filename1 = filename*"_States.txt"
    filename2 = filename*"_SynCur.txt"

    Ni = N-Ne
    
    Tt = Tt+Transient
    
    V = zeros(Float64,N)

    Vin = zeros(Float64,N)

    Vnoise = zeros(Float64,N)

    tsp = ones(Int64,N)
    ttp = ones(Int64,N)

    tnoiL = zeros(Int64,N,Nn);
    SnoiL = zeros(Int64,N,tmax);

    tau2h = tau2/h

    SVin = zeros(Float64,N);

    nn = 1
    
    V0h = -4. *V0d

    xstp = ones(Float64,N);
    x = ones(Float64,N);
    
    tau1h = tau1/h
    tau2h = tau2/h
    taurh = taur/h

    kaph = kap*h

    Vnor = ones(Float64,N)
    Vth = ones(Float64,N).*Vth0

    tauh = ones(Float64,N)

    lambdaE = mu/(tmax*Nn)
    
    lambdaI = nu/(tmax*Nn)
    
    Snei = zeros(Int,N)
    
    io1 = open(filename1, "w")
    
    close(io1)
    
    io2 = open(filename2, "w")
    
    nclos=0

    for t in 1:Tt-1

        Vnor[Ne+1:end] = ((Vsat.-V[Ne+1:end])./Vsat) 
        
        Vnor[1:Ne] = ((Vmin.-V[1:Ne])./Vmin)
        
        SVin[1:end] = (A*Vin).*Vnor

        #############################################################
        
        SnoiL[1:Ne,nn] = reduce(+,rand(Ne,Nn) .< lambdaE ,dims=2)
	
	SnoiL[Ne+1:end,nn] = reduce(+,rand(Ni,Nn) .< lambdaI ,dims=2)

        tnoiL[:,nn] = t.*(SnoiL[:,nn].>0)

        Vnoise[:] = reduce(+,SnoiL.*((t.>=tnoiL).-(t.>=(tnoiL.+tmax))),dims=2)
        
        ##############################################################

        tauh[1:end] = tau1h.*(V.>=Vres) .+ tau2h.*(V.<Vres)
	
        V += (-V .+ SVin .+ ((Vsat.-V)./Vsat).*V0dext.*Vnoise)./tauh

        ttp[1:end] = (V .> Vth)
        
        tsp[1:end] = t.*ttp + tsp.*(1 .-ttp)

        Vin[1:Ne] = (V0d*U.*xstp[1:Ne].*ttp[1:Ne] .+ Vin[1:Ne].*(1. .-ttp[1:Ne])).*(t.<tsp[1:Ne].+tmax)
        
        Vin[Ne+1:end] += V0h*U.*xstp[Ne+1:end].*ttp[Ne+1:end] .- (1 .- ttp[Ne+1:end]).*Vin[Ne+1:end]./tau2h
        
        if taur>0
            
            xstp[1:end] = x.*ttp .+ xstp.*(1. .-ttp)

            x += ((1. .- x)./taurh .- U.*x.*ttp);
        
        end

        nn = (nn+1)*(nn<tmax) + (nn>=tmax)

        Vth[1:end] = Vsat.*ttp .+ Vth.*(1 .-ttp) .- (Vth .- Vth0).*kaph.*(t.> tsp.+ta)

        if (t>=Transient)
            
            if t<=(Transient+10000000)
                
                writedlm(io2, [SVin],',')
            
            end
            
            if (t%dts == 0)
                
                Snei[:] =  Int8.(t .< (tsp .+ ta))
                
                io1 = open(filename1, "a")
                
                writedlm(io1, [Snei],',')
                
                close(io1)

            end
	
            if (t%100000 == 0)

		GC.gc(true)

	   end
    	end
    	
    end
    
    return nothing
end
