
function segintsg(a, b, p, h, m)

    
        T = typeof(a)
        z = zero(T)
        ϵ = eps(T) * 1_000
        I=[0.0,0.0]

        h2, d = h^2, abs(h)
        q2 = p^2+h2
        q=sqrt(q2)
        ra, rb = sqrt(a^2+q2), sqrt(b^2+q2)
#n==-3,-4
        if d < ϵ 

            sgn = zero(T)
            if abs(p) < ϵ 
                I=[z,z]  
            else 
                I = [(a/(p*ra))-(b/(p*rb)),-(1/(4*p^2))*((p*b/rb^2)+atan(b/p)-((p*a/ra^2)+atan(a/p)))]
            end 
            

        else 
            sgn = sign(h)
        

            if abs(p) < ϵ 
            I=[z,z]  
            else 
            I=[sgn*(atan((p*b)/(q2+d*rb)) - atan((p*a)/(q2 + d*ra))),(p/(2*q*h2))*(atan(b/q)-atan(a/q))]
            end       
        #if d < ϵ 
         #   sgn = zero(T) 
        #else 
         #   sign(h)
        #end

        #if abs(p) < ϵ 
         #   I=z  
       # else 
        #    I=p*atan(b/sqrt(q2))/(2*h2*sqrt(q))-p*atan(a/sqrt(q2))/(2*h2*sqrt(q))
       # end

    
        #j = (q2 < ϵ^2) ? (b > 0 ? log(b/a) : log(a/b)) : log(b + rb) - log(a + ra)
        
        #if b < 0 && q2 < max(a2,b2) * (0.5e-3)^2
         #   j = log(a/b) + log((1-(q2/b2)/2) / (1-(q2/a2)/2))
        #else
         #   j = log(b + rb) - log(a + ra)
        #end
        #J = (j,z)
        #K1 = -j

       # j = z
        #J = (j,J[1])

        # n = -1
        #I2 = p*J[2] - h*I1
        #j = (b*rb - a*ra + q2*J[2])/2
        #J = (j,J[1])
        #K2 = j

        # n = 0
        #I3 = (b*p - a*p)/2
        #j = ((b*(b^2+q2)+2*q2*b) - (a*(a^2+q2)+2*q2*a))/3
        #J = (j,J[1])
        #K3 = j/2

    #for i in 4 : N3
     #   Ip = Symbol(:I,i-2)
      #  In = Symbol(:I,i)
       # Kn = Symbol(:K,i)
        #it = quote
         #   n = $i - 3
          #  $In = p/(n+2)*J[2] + T(n/(n+2))*h2*$Ip
           # j = (b*rb^(n+2) - a*ra^(n+2) + (n+2)*q2*J[2])/(n+3)
            #J = (j,J[1])
            #$Kn = j/(n+2)
        #end
       # append!(xp.args, it.args)
    #end
        I[1]=I[1]/h    
        end
    return I
end

function arcintsg(α, m, p, h) 

    
        T = typeof(h)
        P = typeof(m)

        h2 = h^2
        p2 = p^2
        q2 = h2 + p2
        q = sqrt(q2)
        d = sqrt(h2)
        I=[0.0,0.0]
        # n == -3,-4
        

        if norm(h)< eps(T)*1e3
            sgn = zero(T)
            I = [-α/q,-α/(q2*2)]
        else
            sgn = sign(h)
        
    
            I =[-α * (h/q - sgn),-α * (1/q2 - 1/h2)/2]
            I[1]=I[1]/h
        end
        return I
end


      
    
    
function circleintsg(σ, p, h) 

           
                T = typeof(h)
        
                d = norm(h)
                h2 = h^2
                p2 = p^2
                q2 = h2 + p2
                q = sqrt(q2)
                α = σ * 2π
                I=[zero(T),zero(T)]
                # n == -3,-4
                if norm(h)< eps(T)*1e3
                    sgn = zero(T)
                    
                else
                    sgn = sign(h)
                
                    I = [-α * (h/q - sgn),-α * (1/q2 - 1/h2)/2]
                    I[1]=I[1]/h
                end
        return I
            
end

function angle(p,q)
    cs = dot(p,q) / norm(p) / norm(q)
    cs = clamp(cs, -1, +1)
    return acos(cs)
end
        

function inttrianglenegativepowers(ctr, x)

    n = ctr.normal
    h = ctr.height

    ulps = 1000

    ξ = x - h*n

    I=[0.0,0.0]
    

    # segments contributions
    for i in eachindex(ctr.segments)
        a = ctr.segments[i][1]
        b = ctr.segments[i][2]
        t = b - a
        t /= norm(t)
        m = cross(t, n)
        p = dot(a-ξ,m)
        sa = dot(a-ξ,t)
        sb = dot(b-ξ,t)
        P = segintsg(sa, sb, p, h, m)
        I = I+P
        
    end

    # arc contributions
    α = 0
    for i in eachindex(ctr.arcs)
        a = ctr.arcs[i][1]
        b = ctr.arcs[i][2]
        σ = ctr.arcs[i][3]
        p = σ > 0 ? ctr.plane_outer_radius : ctr.plane_inner_radius
        u1 = (a - ξ) / p
        u2 = σ * (n × u1)
        ξb = b - ξ
        α = dot(ξb,u2) >= 0 ? σ*angle(ξb,u1) : σ*(angle(ξb,-u1) + π)
        m = (sin(α)*u1 + σ*(1-cos(α))*u2)
        #α < eps(typeof(α))*1000 && continue
        P = arcintsg(α, m, p, h)
        I = I+P
    end

    # circle contributions
    for i in eachindex(ctr.circles)
        σ = ctr.circles[i]
        p = σ > 0 ? ctr.plane_outer_radius : ctr.plane_inner_radius
        P = circleintsg(σ, p, h)
        I = I+P
    end

    return I
end
#builds the integral of R^-3 R^-4 between a point x and a triangle intersected with a spherical shell centered in x 
function inttrianglenegativepowers(p1,p2,p3,x,r,R)
    ws = workspace(typeof(p1))
    inttrianglenegativepowers(p1,p2,p3,x,r,R,ws)
end

function inttrianglenegativepowers(p1,p2,p3,x,r,R,ws)
    ctr = contour!(p1,p2,p3,x,r,R,ws)
    inttrianglenegativepowers(ctr,x)
end

function inttrianglenegativepowers(p1,p2,p3,x)
    ctr = contour(p1,p2,p3,x)
    inttrianglenegativepowers(ctr,x)
end

#= some test where we know the analytical result 
v4 = point(0.0, 0.0, 1.0)

v1 = point(0.0, 0.0, 0.0)
v2 = point(1.0, 0.0, 0.0)
v3 = point(0.0, 1.0, 0.0)


inttrianglenegativepowers(v1,v2,v3,v4,2.9,4.4)=#

