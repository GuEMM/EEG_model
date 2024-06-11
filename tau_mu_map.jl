include("EEG_LIF_model.jl");

using DelimitedFiles

#Initialize parameters

#Sample time window
dts = 100

#Time step
h = 4/100

#Time domain
T = 20000000

#Transient time domain
Transient = 20000

#Network parameters
ce = 14
Kin = 32

Ne = ce^2
Ni = Int(ce/2)^2

N = Ne+Ni

KK,C,KIout = BuildNet(ce);

#############################
#Obtain Adjacency matrix
A = zeros(Bool,N,N)

ll = zeros(Int,N)

for i in 1:N
    
    l = sum(KK[i,:].>0)
    
    ll[i] = l
    
    A[i,KK[i,1:l]] .= 1
    
end
############################
#Model parameters
Ntau = 30

V0d = 10
V0h = -4. *V0d

nu = 0.0

ddir = "tau_mu_map_$(Ntau)x$(Ntau)/"

isdir(ddir) || mkdir(ddir)

taur = parse(Float32,ARGS[1])
mu = parse(Float32,ARGS[2])

file = "Data_States_LIF_EEG_mu_$(mu)_taur_$(taur)_Vd_$(V0d)"

filename = ddir*file

LIF_Sim_TauMu_States(N,Ne,A,KK,KIout,V0d,taur,mu,nu,filename,T,dts,h,Transient)

println("taur=",taur," mu=",mu, ":: OK! ")

