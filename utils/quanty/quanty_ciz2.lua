Verbosity(0)
-- complex version

-- input file
fin = io.open("quanty.inp", "r")
-- number of spin-orbitals
NF = fin:read("*n")
-- number of electrons
Ne = fin:read("*n")
-- one-body hopping terms
h1e = {}
for i = 1, NF do
    h1e[i] = {}
    for j = 1, NF do 
        r1 = fin:read("*n")
        r2 = fin:read("*n")
        h1e[i][j] = Complex.New(r1, r2)
    end
end
-- two-body interaction terms
v2e = {}
for i = 1, NF do
    v2e[i] = {}
    for j = 1, NF do
        v2e[i][j] = {}
        for k = 1, NF do
            v2e[i][j][k] = {}
            for l = 1, NF do
                r1 = fin:read("*n")
                r2 = fin:read("*n")
                v2e[i][j][k][l] = Complex.New(r1, r2)
            end
        end

    end
end
fin.close()

-- boson dimension
NB = 0

OppC = {}
OppA = {}
for i = 1, NF do 
    OppC[i] = NewOperator("Cr", NF, i-1)
    OppA[i] = NewOperator("An", NF, i-1)
end

-- Hamiltonian
H = 0.
for i = 1, NF do
    for j = 1, NF do
        -- one-body
        H = H + h1e[i][j]*OppC[i]*OppA[j]
        -- two-body
        for k = 1, NF do
            for l = 1, NF do
                H = H + 0.5*v2e[i][j][k][l]
                        *OppC[i]*OppC[k]*OppA[l]*OppA[j]
            end
        end
    end
end

-- asking one lowest state
Npsi = 1

-- initial restriction
StartRestrictions = {NF, NB}
sfill = ""
for i = 1, NF do
    sfill = sfill.."1"
end
StartRestrictions[3] = {sfill, Ne, Ne}

-- calculate ground state
psi = Eigensystem(H, StartRestrictions, Npsi)

-- record
fout = io.open("quanty_ci.out", "w")
fout:write( "energy: ", psi*H*psi, "\n" )
-- one-body density matrix
for i = 1, NF do
    for j = 1, i do
        res = psi*OppC[i]*OppA[j]*psi
        if math.abs(res) > 1e-6 then 
            fout:write(" 1bdm(", i, ", ",j ,"): ", res, "\n" )
        end
    end
end
fout:close()
