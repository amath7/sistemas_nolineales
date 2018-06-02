
function Newton(f, J, x0, N=100)
    x = zeros(N,2) # x es un arreglo de tuplas x,y
    x[1,:] = x0 - (J(x0))^(-1)*f(x0)
    for i in 1:N-1
        if abs(det(J(x[i,:]))) > 1e-4
            x[i+1,:] = x[i,:] .- (J(x[i,:]))^(-1)*f(x[i,:])
        end
    end
    return x[end,:]
end

function PFs_estabilidad(RHS, J, ux, uy, μ, bandera, ϵ=1e-4)
    PFs = []
    estables = []
    inestables = []
    sillas = []
    centros = []
    focos = []
    eig_estables = []
    eig_inestables = []
    eig_sillas = []
    eig_centros = []
    eig_focos = []
    
    #la funcion RHS_df convierte el sist. de ecuaciones diferenciales que regresa RHS
    # a una funcion df que regresa una tupla de vectores que usa quiver
    function RHS_df(x::Float64, y::Float64)
        du = Float64[0.,0.]
        u = Float64[x, y]
        RHS(du, u, μ, 1.)
        P2(du)
    end
    
    # df_RHS() permite usar RHS pasándole como parámetros un vector, la usa Newton()
    function df_RHS(u)
        x = u[1]
        y = u[2]
        return RHS_df(x, y)
    end
    
    for k in ux, l in uy
        x0 = [k, l]
        pf = Newton(df_RHS, x->(J(x,μ)), x0)
        if (norm(df_RHS(pf))<1e-4) & (!any([norm(x.-pf)<ϵ for x in PFs]))
            push!(PFs, pf)
        end
    end
    for i in 1:length(PFs)
        λs = eig(J(PFs[i], μ))[1]
        eigvecs = eig(J(PFs[i], μ))[2]
        reλ = real.(λs)
        Imλ = imag.(λs)
            a = all([x<0 for x in reλ])
            b = all([x>0 for x in reλ])
            c = all([x==0 for x in Imλ])
            d = all([x==0 for x in reλ])
            f = any([x>0 for x in reλ])
        
        if a & c
            push!(estables, PFs[i])
            push!(eig_estables, eigvecs, λs)
            elseif b & c
            push!(inestables, PFs[i])
            push!(eig_inestables, eigvecs, λs)
            elseif c & f
            push!(sillas, PFs[i])
            push!(eig_sillas, eigvecs, λs)
            elseif d
            push!(centros, PFs[i])
            push!(eig_centros, eigvecs, λs)
            else 
            push!(focos, PFs[i])
            push!(eig_focos, eigvecs, λs)
        end
    end
    dict = Dict(cadena => Array[] for cadena in ["estable", "inestable",
                "silla", "centro", "foco"])
    dict_eig = Dict(cadena_eig => Array[] for cadena_eig in ["estable", "inestable", "silla"])
    dict["estable"] = estables
    #dict["eig_estable"] = eig_estables
    dict["inestable"] = inestables
    #dict["eig_inestable"] = eig_inestables
    dict["silla"] = sillas
    #dict["eig_silla"] = eig_sillas
    dict["centro"] = centros
    dict["foco"] = focos
    dict_eig["estable"] = eig_estables
    dict_eig["inestable"] = eig_inestables
    dict_eig["silla"] = eig_sillas
    dict_eig["centro"] = eig_centros
    dict_eig["foco"] = eig_focos
    if bandera == "eigenvector"
    return dict_eig
        else 
        return dict
    end
end

function campo(RHS,
        meshx=linspace(-3, 1, 6), meshy=linspace(-3, 1, 6),
        μ=1., k=10., tmax=4.0, h=0.1)

    #condiciones iniciales en un mesh
    ux = meshx
    uy = meshy
    cond_init = [[ux[i], uy[j]] for i=1:length(ux), j=1:length(uy)]
    for i in 1:length(cond_init)
        prob = ODEProblem(RHS, cond_init[i], (0, tmax), μ)
        sol = solve(prob, saveat=h, force_dtmin=true)
        a = [sol.u[i][1] for i = 1:length(sol.u)]
        sols = hcat(a, [sol.u[i][2] for i =1:length(sol.u)])
        solx = sols[:,1]
        soly = sols[:,2]
        plot!(solx, soly, lw=2, label="")
    end
    #current()
    
    #la funcion RHS_df convierte el sist. de ecuaciones diferenciales que regresa RHS
    # a una funcion df que regresa una tupla de vectores que usa quiver
    function RHS_df(x::Float64, y::Float64)
        du = Float64[0.,0.]
        u = Float64[x, y]
        RHS(du, u, μ, 1.)
        P2(du)/k
    end
    
    Y = meshy
    X = meshx
    pts = vec(P2[(X[i], Y[j]) for i=1:length(X), j=1:length(Y)])
    quiver!(pts, quiver=RHS_df)
end

function estabilidad_nolineal2D(RHS, J, meshx, meshy, μ, k, tmax, h)
    pts = PFs_estabilidad(RHS, J, meshx, meshy, μ, "eigenvalor")
    campo(RHS, meshx, meshy, μ, k, tmax, h)
    #grafica ptos fijos
    for key in keys(pts), j in 1:length(pts[key])
        if length(pts[key]) != 0
            scatter!([pts[key][j][1]], [pts[key][j][2]],
                label=key, ms=5)
        end
    end
    current()
end

function variedades(RHS, J, ux, uy, p, k, tmax, h, α, δ)
    pts = PFs_estabilidad(RHS, J, ux, uy, p, "eigenvalor")
    vector = PFs_estabilidad(RHS, J, ux, uy, p, "eigenvector")
    estabilidad_nolineal2D(RHS, J, ux, uy, p, k, tmax, h)
    for key in keys(vector)
        if length(vector[key]) != 0
            pf = pts[key][1]
            eigen_vec = vector[key][1]
            eigen_val = vector[key][2]
            for i in 1:length(eigen_val)
                signo = [-δ, δ]
                if eigen_val[i] < 0 #variedades estables
                    for δ in signo
                        pert = pf + δ.*eigen_vec[:,i]
                        prob = ODEProblem(RHS, pert, (0., -α*tmax), p)
                        sol = solve(prob, Tsit5(), force_dtmin=true)
                        a = [sol.u[i][1] for i = 1:length(sol.u)]
                        sols = hcat(a, [sol.u[i][2] for i =1:length(sol.u)])
                        solx = sols[:,1]
                        soly = sols[:,2]
                        plot!(solx, soly, lw=2, linestyle=:dash, color=:blue, label="")
                    end
                elseif eigen_val[i] >= 0
                    for δ in signo #variedades inestables
                        pert = pf + δ.*eigen_vec[:,i]
                        prob = ODEProblem(RHS, pert, (0., α*tmax), p)
                        sol = solve(prob, Tsit5(), force_dtmin=true)
                        a = [sol.u[i][1] for i = 1:length(sol.u)]
                        sols = hcat(a, [sol.u[i][2] for i =1:length(sol.u)])
                        solx = sols[:,1]
                        soly = sols[:,2]
                        plot!(solx, soly, lw=2, linestyle=:dash, color=:red, label="")
                    end
                end
            end
        end
    end
    current()
end
