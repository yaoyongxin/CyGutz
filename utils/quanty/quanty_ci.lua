Verbosity(0)

-- input file
fin = io.open("quanty.inp", "r")
-- number of spin-orbitals
NF = fin:read("*n")
-- number of electrons
Ne = fin:read("*n")
-- one-body hopping terms
h1e = {}
for i = 1, NF/2 do
    h1e[i] = {}
    for j = 1, i do 
        h1e[i][j] = fin:read("*n")
    end
end
-- two-body interaction terms
v2e = {}
for i = 1, NF/2 do
    v2e[i] = {}
    for j = 1, i do
        v2e[i][j] = {}
        for k = 1, NF/2 do
            v2e[i][j][k] = {}
            for l = 1, NF/2 do
                v2e[i][j][k][l] = fin:read("*n")
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

-- S-square operator
up_list = {}
dn_list = {}
for i = 1, NF/2 do 
    up_list[i] = 2*(i-1)
    dn_list[i] = 2*i-1
end
OppSsqr = NewOperator("Ssqr", NF, up_list, dn_list)

-- Hamiltonian
H = 0.
for i = 1, NF/2 do
    for j = 1, i do
        for isp = 0, 1 do
            -- one-body
            i_ = 2*i-1+isp
            j_ = 2*j-1+isp
            H = H + h1e[i][j]*OppC[i_]*OppA[j_]
            if i > j then
                H = H + h1e[i][j]*OppC[j_]*OppA[i_]
            end
            -- two-body
            for k = 1, NF/2 do
                for l = 1, NF/2 do
                    for ksp = 0, 1 do
                        k_ = 2*k-1+ksp
                        l_ = 2*l-1+ksp
                        H = H + 0.5*v2e[i][j][k][l]
                                *OppC[i_]*OppC[k_]*OppA[l_]*OppA[j_]
                        if i > j then
                            H = H + 0.5*v2e[i][j][k][l]
                                    *OppC[j_]*OppC[k_]*OppA[l_]*OppA[i_]
                        end
                    end
                end
            end
        end
    end
end

-- constrain to S=0
H = H + OppSsqr

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
fout:write( "spin-square: ", psi*OppSsqr*psi, "\n")
-- one-body density matrix
for i = 1, NF/2 do
    i_ = 2*i-1
    for j = 1, i do
        j_ = 2*j-1
            fout:write(" 1bdm(", i, ", ",j ,"): ", 
                    psi*OppC[i_]*OppA[j_]*psi, "\n" )
    end
end
fout:close()
