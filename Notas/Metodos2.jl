
"Σ(f,x) hace la suma del arreglo producido de aplicar la funcion f al arreglo x, Σf(x)"
Σ(f::Function,x) = sum(f(x))

" ∫(f,a,b,δ) hace la integral numerica de f(x), desde x=a, hasta x=b. Lo hace en pasos de tamaño δ, predefinido
 como 10^-4"
function ∫(f::Function, a,b,δ=1e-4)
    if a<b
        I = Σ(f,a:δ:b)*δ
    else
        I = -Σ(f,b:δ:a)*δ
    end
end

"derivada(f,x0,h) Hace la derivada numerica df(x)/dx, evaluada en el punto x0, con pasos de tamao h predefinidos como 10^-2 "
derivada(f::Function,x,h=1e-2) =(1/280*f(x-4h)-4/105*f(x-3h)+1/5*f(x-2h)-4/5*f(x-h)+4/5*f(x+h)-1/5*f(x+2h)+4/105*f(x+3h)-1/280*f(x+4h))/h


function Adivinanza(f::Function,df::Function, x0)
    x0-f(x0)/df(x0)
end

"Newton1(f, df, x0, δ) Aplica el metodo de Newton para obtener los ceros de la ecuacion f(x)=0. df es la derivada de f, x0 
es una adivinanza inicial y δ es la precision del resultado, predefinido en 0.00001"
function Newton1(f::Function, df::Function, x0, δ=0.00001)
    contador = 0
    while abs(f(x0)) > δ
        contador +=1
        if contador > 1e5
            return ("no se encuentra ningún cero después de 10^5 pasos")
        end
        x0 = Adivinanza(f,df,x0)
    end
    return x0
end

function Adivinanza(f::Function, x0)  
    x0-f(x0)/derivada(f,x0)           
end

"Newton2(f, x0, δ) Aplica el metodo de Newton usando derivadas numericas (metodo de la secante) para obtener la
solucion de la ecuacion f(x)=0 con una presicion δ predefinida igual a 0.00001. x0 es una adivinanza inicial"
function Newton2(f::Function, x0, δ=0.00001) 
    contador = 0                             
    while abs(f(x0)) > δ
        contador +=1
        if contador > 1e5
            return ("no se encuentra ningún cero después de 10^5 pasos")
        end
        x0 = Adivinanza(f,x0)   
    end
    return x0
end

"deriva_vec(f,x,i,hh) obtiene la derivada vectorial numerica con respecto a la i-esima variable de la funcion f, 
sobre la variable x  con un paso hh, predefinido como 10^-2"
function deriva_vec(f::Function,x,i,hh=1e-2) 
    m = length(x)
    h = zeros(m)
    h[i] = hh    
    de =(1/280*f(x-4h)-4/105*f(x-3h)+1/5*f(x-2h)-4/5*f(x-h)+4/5*f(x+h)-1/5*f(x+2h)+4/105*f(x+3h)-1/280*f(x+4h))/hh
end

"Jacobiano(f,x0, h) regresa el Jacobiano numerico de la funcion f, evaluado en x0, con una presicion h predefinido
como 10^-2"
function Jacobiano(f,x0, h=1e-2)
    m = length(x0)  
    n = length(f(x0)) 
    df = zeros(n) 
    J = Float64[] 
    for i = 1:m
        df = deriva_vec(f,x0,i,h) 
        append!(J,df) 
    end
    return reshape(J,n,m) 
end 

Adivinanzav(x0, f:: Function, h) = x0 - *(inv(Jacobiano(f,x0,h)),f(x0))

"Newtonv(f, x0; δ, h) usa el metodo de Newton para regresar el vector solucion al sistema de ecuaciones f(x)=0, donde 
f es una funcion vectorial, 0 es un vector nulo de la misma dimension de f y de x. δ es la precision de la solucion
predefinida en 0.00001, h es el tamaño del paso de las derivadas, predefinido como 0.01"
function Newtonv(f::Function, x0; δ=0.00001, h=0.01)
    contador = 0
    while norm(f(x0))>δ 
        contador += 1
        if contador > 1e5
            return "no se encontro solucion en 10^5 pasos"
        end
        x0 = Adivinanzav(x0,f,h)
    end
    return x0
end

"Euler(f, x0,t0,tf,h) Usando el metodo de Euler regresa la solucion numerica a la ecuacion diferencial dx/dt = f(x,t). La solucion que regresa
esta entre los puntos t=t0 y t=tf, con pasos de tamaño h. "
function Euler(f::Function, x0,t0,tf,h)
    X = Float64[x0]
    T = Float64[t0]
    x = x0
    for t in t0:h:tf
        push!(T,t)
        x += h.*f(x,t)
        push!(X,x)
    end
    return T,X
end

"Euler2(f, x0,t0,tf,h) Usando el metodo de Euler multipaso de orden 2 regresa la solucion numerica a la ecuacion diferencial d^2x/dt^2 = f(x,t). 
La solucion que regresa esta entre los puntos t=t0 y t=tf, con pasos de tamaño h. "
function Euler2(f::Function,x0,v0,t0,tf,h)  
    x1 = x0+h*v0
    X = [x0,x1]
    T = [t0,t0+h]
    for t in t0+2h:h:tf
        push!(T,t)
        x2 = x1
        x1 = f(x1,t)*h^2 + 2*x1-x0
        x0 = x2
        push!(X,x1)
    end
    return T,X
end

"RK4(f,a,b,N,α) obtiene, mediante el metodo Runge-Kuta de orden 4 la solucion numerica de la ecuacion diferencial dx/dt = f(x,t) desde el punto t=a, hasta 
el punto t=b, con condicion inicial x(a)=α. Se obtienen N puntos de la solucion"
function RK4(f::Function,a,b,N,α)
    X = Float64[]
    T = Float64[]
    h = (b-a)/N
    t = a
    w = α
    push!(T,t)
    push!(X,w)
    for i in 1:N
        K1 = h*f(w,t)
        K2 = h*f(w+K1/2,t+h/2)
        K3 = h*f(w+K2/2,t+h/2)
        K4 = h*f(w+K3,t+h)
        w = w+(K1+2*K2+2*K3+K4)/6
        t = i*h+a
        push!(X,w)
        push!(T,t)
    end
    return T,X
end

"RK4v(a,b,N,α,f) obtiene, mediante el metodo Runge-Kuta de orden 4 la solucion numerica de la ecuacion diferencial dx/dt = f(x,t) desde el punto t=a, hasta 
el punto t=b, con condicion inicial x(a)=α. Se obtienen N puntos de la solucion"
function RK4(a,b,N,α,f::Function...)
    X = Float64[]
    T = Float64[]
    h = (b-a)/N
    t = a
    w = α
    push!(T,t)
    push!(X,w)
    for i in 1:N
        K1 = h*f(w,t)
        K2 = h*f(w+K1/2,t+h/2)
        K3 = h*f(w+K2/2,t+h/2)
        K4 = h*f(w+K3,t+h)
        w = w+(K1+2*K2+2*K3+K4)/6
        t = i*h+a
        push!(X,w)
        push!(T,t)
    end
    return T,X
end

"RK4v(a,b,N,α,f) obtiene, mediante el metodo Runge-Kuta de orden 4 la solucion numerica 
de la ecuacion diferencial dx/dt = f(x,t) desde el punto t=a, hasta 
el punto t=b, con condicion inicial x(a)=α. Se obtienen N puntos de la solucion"
function RK4v(a,b,N,α,f::Function...)
    F = collect(f)
    n = length(α)
    println(n)
    m = length(F)
    if m !=n
        return "el numero de variables α debe ser igual al numero de funciones f"
    end
    X = Float64[]
    T = Float64[]
    h = (b-a)/N
    t = a
    w = α
   
    push!(T,t)
    append!(X,w)
    K1 = zeros(n)
    K2 = zeros(n)
    K3 = zeros(n)
    K4 = zeros(n)
    for i in 1:N 
        for j in 1:n
            K1[j] = h*F[j](w,t)
        end
        for j in 1:n
            K2[j] = h*F[j](w+K1/2,t+h/2)
        end
        for j in 1:n
            K3[j] = h*F[j](w+K2/2,t+h/2)
        end
        for j in 1:n
            K4[j] = h*F[j](w+K3,t+h)
        end
        w = w+(K1+2*K2+2*K3+K4)/6
        t = i*h+a
        append!(X,w)
        push!(T,t)
    end
    return T,reshape(X,n,N+1)' 
end     
