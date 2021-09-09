using QuantumOptics

using PyPlot
using Plots

b = SpinBasis(1//2)
sx = 1/2*sigmax(b)
sy = 1/2*sigmay(b)
sz = 1/2*sigmaz(b)
XX = sx ⊗ sx 
YY = sy ⊗ sy 
ZZ = sz ⊗ sz

Jxy=1
Jz=1

function basis(N)    
    spinchain=b
    for ispin in 2:N
        spinchain = spinchain ⊗ b
        
    end
    return spinchain # same local
end

function CrtModelH(N)
    Htwoint=Jxy*(XX+YY)+Jz*ZZ
    HXXZ = embed(basis(N), [1,2] , Htwoint)
    for s in 2:N-1
        HXXZ += embed(basis(N), [s,s+1] , Htwoint)
    end
    return HXXZ
end

function solvemodl(N)
    E5 = eigenenergies(CrtModelH(N), 5)
    println(" ")
    #println("Lowest five eigenvalues");
    #println(E5);
    println(" ")
    println(" ")
    # full spectrum
    #EE, UU = eigenstates(dense(CrtModelH(N)))
    println(" ")
    #println("All Eigenvalues");
    #println(EE);
    println(" ")
    return E5
end

s5z = embed(basis(10) , [5] , sz)
s6z = embed(basis(10) , [6] , sz)

sup = spinup(b)
sdo = spindown(b)
psi0 = sup
for ispin in 2:10
    global psi0
    if ispin == 5
       psi0 = psi0 ⊗ sdo
    else 
       psi0 = psi0 ⊗ sup
    end 
end

psi0


times = 0:.1:6

tout, psi_t = timeevolution.schroedinger(times, psi0, CrtModelH(10))


m5z = expect(s5z, psi_t)
m6z = expect(s6z, psi_t)

using ITensors

pyplot()

let
  N = 10
  cutoff = 1E-8
  tau = 0.05
  ttotal = 1.5*4

  # Compute the number of steps to do
  Nsteps = Int(ttotal/tau)

  # Make an array of 'site' indices
  s = siteinds("S=1/2",N)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ITensor[]
  for j=1:N-1
    s1 = s[j]
    s2 = s[j+1]
    hj =       op("Sz",s1) * op("Sz",s2) +
               op("Sx",s1) * op("Sx",s2) +
               op("Sy",s1) * op("Sy",s2)
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
  end
  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates,reverse(gates))

  c = div(N,2) # center site

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> n!=c ? "Up" : "Dn")

  Szc=[]
  Szc2=[]
  # Compute and print initial <Sz> value on site c
  t = 0.0
  Sz  = ITensors.expect(psi,"Sz";site_range=c:c)
  Sz2 = ITensors.expect(psi,"Sz";site_range=c+1:c+1)
  println("$t $Sz $Sz2")
  append!(Szc,Sz)
  append!(Szc2,Sz2)

  # Do the time evolution by applying the gates
  # for Nsteps steps and printing <Sz> on site c
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    t += tau
    Sz  = ITensors.expect(psi,"Sz";site_range=c:c)
    Sz2 = ITensors.expect(psi,"Sz";site_range=c+1:c+1)
    println("$t $Sz $Sz2")
    append!(Szc,Sz)
    append!(Szc2,Sz2)
  end
  times = 0:tau:ttotal
  plot(times, Szc)
  plot(times, Szc2)
  plot(tout, m5z,linestyle = ":",label="line 1",linewidth = 5)
  plot(tout, m6z,linestyle = ":",linewidth = 5)
  return
end




