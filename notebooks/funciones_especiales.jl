
# ESTE ES UN ARCHIVO.JL DE FUNCIONES ESPECIALES PARA SISTEMAS NO LINEALES


doc"""los parámetros por default son: 

campo(RHS, df, meshx_sol=linspace(-3, 1, 6), meshy_sol=linspace(-3, 1, 6), meshx_campo=linspace(-10, 10, 15),

meshy_campo=linspace(-1, 1, 15), xlims=(-10,10), ylims=(-1,1), tmax=4.0, h=0.1)

EJEMPLO:

Considérese el sig. sistema de ecuaciones diferenciales no lineales acopladas:


$$\dot{x} = x + e^{-y}$$

$$\dot{y} = -y$$

Las correspondientes formas para RHS (que es la función que usa por default
el integrador de julia DifferentialEquations.jl) y para df son:

function RHS1(du::Array{Float64,1}, u::Array{Float64,1}, p::Void, t::Float64)

$du[1] = u[1] + exp(-u[2])$

$du[2] = -u[2]$

end

ejemplo para df:

$df1(x, y) = P2(x+exp(-y), -y/10)$"""
function campo(RHS,
        meshx=linspace(-3, 1, 6), meshy=linspace(-3, 1, 6), k=10., tmax=4.0, h=0.1)
    
    #xlims=(-10,10), ylims=(-1,1)
    #plot(xlabel="x(t)", ylabel="y(t)", label="", title="espacio fase", 
     #   xlims=xlims, ylims=ylims)
    #condiciones iniciales en un mesh
    ux = meshx
    uy = meshy
    cond_init = [[ux[i], uy[j]] for i=1:length(ux), j=1:length(uy)]
    for i in 1:length(cond_init)
        prob = ODEProblem(RHS, cond_init[i], (0, tmax))
        sol = solve(prob, saveat=h)
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
        RHS(du, u, nothing, 1.)
        P2(du)/k
    end
    
    Y = meshy
    X = meshx
    pts = vec(P2[(X[i], Y[j]) for i=1:length(X), j=1:length(Y)])
    quiver!(pts, quiver=RHS_df, xlabel = "dx/dt", ylabel = "dy/dt")
end


function Newton(f, J, x0, N=100)
    x = zeros(N,2) # x es un arreglo de tuplas x,y
    x[1,:] = x0 - (J(x0))^(-1)*f(x0)
    for i in 1:N-1
        x[i+1,:] = x[i,:] .- (J(x[i,:]))^(-1)*f(x[i,:])
    end
    return x[end,:]
end


function PFs_estabilidad(RHS, J, ux, uy, bandera, ϵ=1e-4)
    PFs = []
    estables = []
    inestables = []
    sillas = []
    centros = []
    focos = []
    eig_estables = []
    eig_inestables = []
    eig_sillas = []
    #eig_centros = []
    #eig_focos = []
    
    #la funcion RHS_df convierte el sist. de ecuaciones diferenciales que regresa RHS
    # a una funcion df que regresa una tupla de vectores que usa quiver
    function RHS_df(x::Float64, y::Float64)
        du = Float64[0.,0.]
        u = Float64[x, y]
        RHS(du, u, nothing, 1.)
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
        pf = Newton(df_RHS, J, x0)
        if (!any([norm(x.-pf)<ϵ for x in PFs]))
            push!(PFs, pf)
        end
    end
    for i in 1:length(PFs)
        λs = eig(J(PFs[i]))[1]
        eigvecs = eig(J(PFs[i]))[2]
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
            else 
            push!(focos, PFs[i])
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
    if bandera == "eigenvector"
    return dict_eig
        else 
        return dict
    end
end

function estabilidad_nolineal2D(RHS, J, meshx, meshy, k, tmax, h)
    pts = PFs_estabilidad(RHS, J, meshx, meshy, "eigenvalor")
    pts_eig = PFs_estabilidad(RHS, J, meshx, meshy, "eigenvector")
    campo(RHS, meshx, meshy, k, tmax, h)
    #grafica ptos fijos
    for key in keys(pts), j in 1:length(pts[key])
        if length(pts[key]) != 0
            scatter!([pts[key][j][1]], [pts[key][j][2]],
                label=key, ms=5, size=(700,400))
            if (length(pts_eig) != 0) & !((key == "centro") || (key == "foco")) 
                #@show pts[key]
                v1 = pts_eig[key][1][:,1]
                λ1 = pts_eig[key][2][1]
                m1 = v1[2]/v1[1]
                v2 = pts_eig[key][1][:,2]
                λ2 = pts_eig[key][2][2]
                m2 = v2[2]/v2[1]
                v3 = pts[key][1]
                l1 = l2 = "a"
                if λ1 < 0
                    l1 = "espacio estable"
                    l2 = "espacio inestable"
                else
                    l1 = "espacio inestable"
                    l2 = "espacio estable"
                end
                plot!(meshx, x->(m1*(x-v3[1])+v3[2]), label=l1, linestyle=:dash, lw=2)
                plot!(meshx, x->(m2*(x-v3[1])+v3[2]), label=l2, linestyle=:dash ,lw=2)
            end
        end
    end
    current()
end
