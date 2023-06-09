
#corregere val{7} nelle wilton 
#correggere i diviso 6
#valutare operazioni da fare in testa in inttriangtriang, ad esempio calcolare superfici o lunghezze degli edge
#riassumere espressioni composte e che si ripetono
#problema di t=0 per coinc triang (definito per t=0 o per t>10^-6, nel mezzo da nan)
function intpointtriangle(r1,r2,r3,rp,S,t1,t2)
    maxdist=norm(r1-rp)+norm(r1-r2)+norm(r2-r3)+norm(r1-r3)
    a=0.0
    b=0.0
    a=1/(2*S)
    c=[wiltonints(r1,r2,r3,rp,t1,t2,Val{7})[1][2],inttrianglenegativepowers(r1,r2,r3,rp,t1,t2)[1],inttrianglenegativepowers(r1,r2,r3,rp,t1,t2)[2]]
    d=inttrianglenegativepowers(r1,r2,r3,rp,t2,maxdist)
    b=a*(c[1]-3*(t1^2)*c[2]+2*(t1^3)*c[3])
    b+=a*(3*(t2^2-t1^2)*d[1]-2*(t2^3-t1^3)*d[2])
    return b
end

function intlinelineexact(a1,a2,a1′,a2′,temp1,temp2)

    #r=r12 r'=r12'
    #l12=norm(a1-a2)
    #l12′=norm(a1′-a2′)
    r=a1-a2
    r′=a1′-a2′
    γ=-atan(r[2]/r[1])
    β=-atan((r[3]/r[1])/(sqrt(1+(r[2]/r[1])^2)))
    Ry=[cos(β) 0 -sin(β); 0 1 0; sin(β) 0 cos(β)]
    Rz=[cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1]
    r′=Ry * Rz * r′
    α=-atan(r′[3]/r′[2])
    Rx=[1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]

    R=Rx * Ry * Rz 

    r22=R*(a2-a2′)
    r12′=R*(a1′-a2′)
    r12=R*(a1-a2)  
    r2′=R*(a2′)
    r2=r22+r2′
    #radius=sqrt(r22[3]^2+(r12′[2]*λ′+r2′[2]-r2[2])^2+(r12′[1]*λ′+r2′[1]-r12[1]*λ-r2[1])^2)
    #λ->[0,1],λ′->[0,1]
    #λ′=(t′-r2′[2]+r2[2])/r12′[2], t′->[r2′[2]-r2[2],(r12′[2]+r2′[2]-r2[2])],   dλ′=dt′(1/r12′[2]) 
    #radius=sqrt(r22[3]^2+(λ′)^2+(r12′[1]*(t′/r12′[2])-r12′[1]*(r2′[2]-r2[2])/r12′[2]+r2′[1]-r12[1]*λ-r2[1])^2)
    a=r12′[1]/r12′[2]
    b=-r22[1]+(r12′[1]*(r22[2]))/r12′[2]
    #radius=sqrt(r22[3]^2+(t′)^2+((a*t′+b)-r12[1]*λ)^2)
    #λ=p/r12[1], p->[0,r12[1]], dλ=dt/r12[1]
    #radius=sqrt(r22[3]^2+(t′)^2+((a*t′+b)-p)^2)
    #p=t+at'+b, t->[-at′-b,r12[1]-at′-b]
    #radius=sqrt(r22[3]^2+(t′)^2+(t)^2)
    t1′=-r22[2]
    t1=-a*(-r22[2])-b
    t2′=-r22[2]
    t2=-a*(-r22[2])-b+r12[1]
    t3′=-r22[2]+r12′[2]
    t3=-a*(r12′[2]-r22[2])-b+r12[1]
    t4′=-r22[2]+r12′[2]
    t4=-a*(r12′[2]-r22[2])-b
    p1=point(t1′,t1,0)
    p2=point(t2′,t2,0)
    p3=point(t3′,t3,0)
    p4=point(t4′,t4,0)
    p=point(0,0,r22[3])
    c=[wiltonints(p1,p2,p3,p,temp1,temp2,Val{7})[1][2]+wiltonints(p1,p3,p4,p,temp1,temp2,Val{7})[1][2],inttrianglenegativepowers(p1,p2,p3,p,temp1,temp2)[1]+inttrianglenegativepowers(p1,p3,p4,p,temp1,temp2)[1],inttrianglenegativepowers(p1,p2,p3,p,temp1,temp2)[2]+inttrianglenegativepowers(p1,p3,p4,p,temp1,temp2)[2]]
    maxdist=2*norm(a1-a2)+norm(a2-a2′)+norm(a2′-a1′)
    d=inttrianglenegativepowers(p1,p2,p3,p,temp2,maxdist)+inttrianglenegativepowers(p1,p3,p4,p,temp2,maxdist)
    I=(1/abs(r12′[2]*r12[1]))*(c[1]-3*(temp1^2)*c[2]+2*(temp1^3)*c[3])
    I+=(1/abs(r12′[2]*r12[1]))*(3*(temp2^2-temp1^2)*d[1]-2*(temp2^3-temp1^3)*d[2])
    return I

end

function intlinetriangle(a1,a2,b1,b2,b3,S,t1,t2)
    A= hcat(b1-b3, b2-b3, a2-a1)
    #find the intersection between the support of the edge a1 a2 (line) and the support of the triangle b1 b2 b3 (plane)
    normal = (b1-b3) × (b2-b3)
    normal /= norm(normal)
    if det(A)==zero(Float64)
        h = (b1 - a1) ⋅ normal
        if abs(h)== zero(Float64)
            A2=A[1:3,1:2]
                s2=A2\(a2-b3) #not a square matrix but for sure there is the solution for this system
                s=[s2[1], s2[2], 0.0]
                #j="edge 3 coplanar with triangle prime"
        else 
            A2=A[1:3,1:2]
            z12=A2\(a1-a2)
            I=(#(scalingpoint[1]-(z2′)[1])*intlineline(r2,r3,r1′,r2′,t1,t2)
            - (z12)[1] * intlinelinedeg1(b2,b3,a1,a2,t1,t2)+
               #(scalingpoint[2]-(z2′)[2])*intlineline(r1,r3,r1′,r2′,t1,t2)
               -(z12)[2] * intlinelinedeg1(b1,b3,a1,a2,t1,t2)+
               #(-(scalingpoint[1]-(z2′)[1])-(scalingpoint[2]-(z2′)[2]))*intlineline(r1,r2,r1′,r2′,t1,t2)
               +((z12)[1] + (z12)[2])*intlinelinedeg1(b1,b2,a1,a2,t1,t2)+
               intpointtriangle(b1,b2,b3,a1,S,t1,t2))
       
           return I
        end
    else
        s=A\(a2-b3)
    end

    I=(s[3]*intpointtriangle(b1,b2,b3,a2,S,t1,t2)+(1-s[3])*intpointtriangle(b1,b2,b3,a1,S,t1,t2)+
        s[1]*intlinelineexact(b2,b3,a1,a2,t1,t2)+
        s[2]*intlinelineexact(b1,b3,a1,a2,t1,t2)+
        (1-s[1]-s[2])*intlinelineexact(b1,b2,a1,a2,t1,t2))

    return I
end

function inttriangletriangle(a1,a2,a3,b1,b2,b3,t1,t2)
    l12=norm(a2-a1)
    S=sqrt((l12^2)*dot(a3-a1,a3-a2)-dot(a3-a1,a2-a1)dot(a3-a2,a2-a1))/2

    l12b=norm(b2-b1)
    Sb=sqrt((l12b^2)*dot(b3-b1,b3-b2)-dot(b3-b1,b2-b1)dot(b3-b2,b2-b1))/2

    A= hcat(a1-a3, a2-a3, b3-b1, b3-b2) 
    A4 = vcat(A,reshape([1,1,0,0],1,4))#find the intersection between the support of the edge 3 a (line) and the support of the triangle b (plane)
        if det(A4)==zero(Float64) 
            A4= vcat(A,reshape([1,0,0,0],1,4))
            if det(A4)==zero(Float64)
                normal = (b1-b3) × (b2-b3)
                normal /= norm(normal)
                h = (b1 - a1) ⋅ normal
                if abs(h) ==zero(Float64)
                    A2=A4[1:3,1:2]
                    s2=A2\(b3-a3) #not a square matrix but for sure there is the solution for this system
                    s=[s2[1], s2[2], 0.0, 0.0]
                    #j=" coplanar trinalges"
                else 
                     #use parallel formula
                    return "parallel trinagles"
                end
            else
                s=A4\vcat((b3-a3),0)
            end
        else
            s=A4\vcat((b3-a3),1)
        end

        I=4*S*Sb*(s[3]*intlinetriangle(b2,b3,a1,a2,a3,S,t1,t2)+s[4]*intlinetriangle(b3,b1,a1,a2,a3,S,t1,t2)+
            (1-s[3]-s[4])*intlinetriangle(b1,b2,a1,a2,a3,S,t1,t2)+
            s[1]*intlinetriangle(a2,a3,b1,b2,b3,Sb,t1,t2)+s[2]*intlinetriangle(a3,a1,b1,b2,b3,Sb,t1,t2)+
            (1-s[1]-s[2])*intlinetriangle(a1,a2,b1,b2,b3,Sb,t1,t2))
        
        return I
end


function inttriangletriangleadjacent(a1,a2,a3,b1,b2,b3,t1,t2)
    l12=norm(a2-a1)
    S=sqrt((l12^2)*dot(a3-a1,a3-a2)-dot(a3-a1,a2-a1)dot(a3-a2,a2-a1))/2

    l12b=norm(b2-b1)#corregere con l12b=l12 in questo caso particolare
    Sb=sqrt((l12b^2)*dot(b3-b1,b3-b2)-dot(b3-b1,b2-b1)dot(b3-b2,b2-b1))/2

    #@assert a1==b1 a2==b2
    I=(intlinelineexact(a3,a1,b2,b3,t1,t2)+intpointtriangle(a1,a2,a3,b3,S,t1,t2)+
    intlinelineexact(a3,a2,b1,b3,t1,t2)+intpointtriangle(b1,b2,b3,a3,Sb,t1,t2))*(4/6)*(S*Sb)

    return I
end


function intpointtrianglezerotime(r1,r2,r3,S,t2)
    maxdist=norm(r1-r2)+norm(r2-r3)+norm(r1-r3)
    a=0.0
    b=0.0
    a=1/(2*S)
    c=inttriangtimezero(r1,r2,r3,t2)
    d=inttrianglenegativepowers(r1,r2,r3,r3,t2,maxdist)
    b=a*c+a*(3*(t2^2)*d[1]-2*(t2^3)*d[2])
    return b
end

function intlinelineexactzerotime(a1,a2,a1′,a2′,temp2)

    #r=r12 r'=r12'
    #l12=norm(a1-a2)
    #l12′=norm(a1′-a2′)
    r=a1-a2
    r′=a1′-a2′
    γ=-atan(r[2]/r[1])
    β=-atan((r[3]/r[1])/(sqrt(1+(r[2]/r[1])^2)))
    Ry=[cos(β) 0 -sin(β); 0 1 0; sin(β) 0 cos(β)]
    Rz=[cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1]
    r′=Ry * Rz * r′
    α=-atan(r′[3]/r′[2])
    Rx=[1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]

    R=Rx * Ry * Rz 
    r12′=R*(a1′-a2′)
    r12=R*(a1-a2)  
    #radius=sqrt(r22[3]^2+(r12′[2]*λ′+r2′[2]-r2[2])^2+(r12′[1]*λ′+r2′[1]-r12[1]*λ-r2[1])^2)
    #λ->[0,1],λ′->[0,1]
    #λ′=(t′-r2′[2]+r2[2])/r12′[2], t′->[r2′[2]-r2[2],(r12′[2]+r2′[2]-r2[2])],   dλ′=dt′(1/r12′[2]) 
    #radius=sqrt(r22[3]^2+(λ′)^2+(r12′[1]*(t′/r12′[2])-r12′[1]*(r2′[2]-r2[2])/r12′[2]+r2′[1]-r12[1]*λ-r2[1])^2)
    #radius=sqrt(r22[3]^2+(t′)^2+((a*t′+b)-r12[1]*λ)^2)
    #λ=p/r12[1], p->[0,r12[1]], dλ=dt/r12[1]
    #radius=sqrt(r22[3]^2+(t′)^2+((a*t′+b)-p)^2)
    #p=t+at'+b, t->[-at′-b,r12[1]-at′-b]
    #radius=sqrt(r22[3]^2+(t′)^2+(t)^2)
    p1=point(0.0,0.0,0.0)
    p2=point(0.0,r12[1],0.0)
    p3=point(r12′[2],-r12′[1]+r12[1],0)
    p4=point(r12′[2],-r12′[1],0)

    c=inttriangtimezero(p2,p3,p1,temp2)+inttriangtimezero(p4,p3,p1,temp2)
    maxdist=2*norm(a1-a2)+norm(a2-a2′)+norm(a2′-a1′)
    d=inttrianglenegativepowers(p1,p2,p3,p1,temp2,maxdist)+inttrianglenegativepowers(p1,p4,p3,p1,temp2,maxdist)
    I=(1/abs(r12′[2]*r12[1]))*(c+3*(temp2^2)*d[1]-2*(temp2^3)*d[2])
    return I

end



function intcoinctriangles(a1,a2,a3,t1,t2)
    #t1=0
    l12=norm(a2-a1)
    S=sqrt((l12^2)*dot(a3-a1,a3-a2)-dot(a3-a1,a2-a1)dot(a3-a2,a2-a1))/2
    if t1==0.0
        I=(intlinelineexactzerotime(a1,a3,a2,a3,t2)+intpointtrianglezerotime(a1,a2,a3,S,t2))*(8/6)*(S^2)
    else
        I=(intlinelineexact(a1,a3,a2,a3,t1,t2)+intpointtriangle(a1,a2,a3,a3,S,t1,t2))*(8/6)*(S^2)
    end
    return I
end

intcoinctriangles(v1′,v2′,v3′,0,0.4)
quadrule4dtriang(v1′,v2′,v3′,v1′,v2′,v3′,0.0000000000001,0.4)



#test da cancellare
v3 = point( 0.0, 0.0, 0.0) #r3 #attenzione hai cambiato questi punti rispetto allaltro file

v1 = point(1.0, 0.0, 0.0) #r1

v2 = point( 0.0, 1.0, 0.0) #r2


v3′ = point(3.622579672754069, 4.478943189230791, -8.174104502769225) #r3'

v1′ =point(5.370603720688417, 3.683103988470699, -7.58889719823895) #r1'

v2′ = point( 3.851213402629311, 2.5141835154553185, -8.879585383227141) #r2'

vc=point()

intcoinctriangles(v1,v2,v3,0.00000001,0.7)
quadrule4dtriang(v1,v2,v3,v1,v2,v3,0.0,0.7)

wiltonints(v1,v2,v3,v4,Val{1})[1][2]

v4=point(0.25,0.25,0.0)


inttriangletriangleadjacent(v1,v2,v3,v1,v2,v3,0.1,30.0)

function quadrule2dtriang(r1,r2,r3,rp,t2)
    T = simplex(r1, r2, r3)


    #testing quadrature rules

    qpsT = quadpoints(T,13)

    function f(p,rc)

        x = cartesian(p)


        R = norm(x-rc)


        if   R ≤ t2

            return 1/R

        else

            return zero(R)

        end

    

        end

    #result 4d quadrule

    r = 0.0

    for qpT in qpsT

        pT, wT = qpT

            r +=wT * f(pT,rp)


    end

    return r

end


function inttriangtimezero(r1,r2,r3,t2) #definisce lintegral wiltonint di 1/R fra un triangolo e il suo vertice 3 (terzo argomento)
    l12=norm(r2-r1)
    S=sqrt((l12^2)*dot(r3-r1,r3-r2)-dot(r3-r1,r2-r1)dot(r3-r2,r2-r1))/2
    h=2*S/l12 

    p=roots(Polynomial([dot(r2-r3,r2-r3)-(t2)^2,2*dot(r1-r2,r2-r3),dot(r1-r2,r1-r2)]))
    
    if imag(p[1])!=0.0
        segint=0.0
        arcint=t2*((atan((dot(r2-r1,r2-r3)/l12)/h)+atan((l12-dot(r2-r1,r2-r3)/l12)/h)))
    else
        p=sort!(p)
        if p[1]<=0
            if p[2]<=0
                segint=0.0
                arcint=t2*((atan((dot(r2-r1,r2-r3)/l12)/h)+atan((l12-dot(r2-r1,r2-r3)/l12)/h)))
            elseif p[2]>=1
                segint=h*((asinh((dot(r2-r1,r2-r3)/l12)/h)+asinh((l12-dot(r2-r1,r2-r3)/l12)/h)))
                arcint=0.0
            else
                rseg1=(r1-r2)*p[2]+r2
                lseg12=norm(r2-rseg1)
                segint=h*((asinh((dot(r2-rseg1,r2-r3)/lseg12)/h)+asinh((lseg12-dot(r2-rseg1,r2-r3)/lseg12)/h)))
                arcint=t2*((atan((dot(rseg1-r1,rseg1-r3)/(l12-lseg12))/h)+atan(((l12-lseg12)-dot(rseg1-r1,rseg1-r3)/(l12-lseg12))/h)))
            end
        elseif p[1]<=1
            if p[2]>=1
                e=[p[1],1.0]
                rseg2=(r1-r2)*p[1]+r2
                lseg12=norm(r1-rseg2)
                segint=h*((asinh((dot(rseg2-r1,rseg2-r3)/lseg12)/h)+asinh((lseg12-dot(rseg2-r1,rseg2-r3)/lseg12)/h)))
                arcint=t2*((atan((dot(r2-rseg2,r2-r3)/(l12-lseg12))/h)+atan(((l12-lseg12)-dot(r2-rseg2,r2-r3)/(l12-lseg12))/h)))
                
            else
                rseg1=(r1-r2)*p[2]+r2
                rseg2=(r1-r2)*p[1]+r2
                lseg12=norm(rseg1-rseg2)
                l11=norm(rseg1-r1)
                l22=norm(rseg2-r2)
                segint=h*(asinh((dot(rseg2-rseg1,rseg2-r3)/lseg12)/h)+asinh((lseg12-dot(rseg2-rseg1,rseg2-r3)/lseg12)/h))
                arcint=t2*(((atan((dot(r2-rseg2,r2-r3)/(l22))/h)+atan(((l22)-dot(r2-rseg2,r2-r3)/(l22))/h)))+
                    ((atan((dot(rseg1-r1,rseg1-r3)/(l11))/h)+atan(((l11)-dot(rseg1-r1,rseg1-r3)/(l11))/h))))
            end
        else
            segint=0.0
            arcint=t2*((atan((dot(r2-r1,r2-r3)/l12)/h)+atan((l12-dot(r2-r1,r2-r3)/l12)/h)))
        end
    end
    I=arcint+segint
    return I
end


function quadrule4dtriang(r1,r2,r3,r1′,r2′,r3′,t1,t2)
    T = simplex(r1, r2, r3)

    S = simplex(r1′, r2′, r3′)

    #testing quadrature rules

    #qpsT = quadpoints(T,13)

     qpsS = quadpoints(S,13)

    #=function f(p)

        x = cartesian(p)

        y =

        R = norm(x-y)


            if   t1 ≤ R ≤ t2

                return 1/R

            else

                return zero(R)

            end

        end=#

    #result 4d quadrule

    r = 0.0

   # for qpT in qpsT

       # pT, wT = qpT

        for qpS in qpsS

            pS, wS = qpS
                x=cartesian(pS)
            r += wS* wiltonints(r1,r2,r3,x,t1,t2,Val{3})[1][2]

        end

    #end

    return r

end

quadrule4dtriang(v1,v2,v3,v1,v2,v3,0.0,0.4)


inttriangtimezero(v2,v3,v1,1.1)

quadrule4dtriang(v1,v2,v3,v1,v2,v3,1.1,1,3)

#testare integrale generale 4d zero time 

        

function quadrulelineline(r1,r2,r3,rp,t1,t2)
    T = simplex(r1, r2)
    S = simplex(r3, rp)

    #testing quadrature rules

    qpsT = quadpoints(T,50)

    qpsS = quadpoints(S,50)

    function f(p,q)

        x = cartesian(p)

        y = cartesian(q)

        R = norm(x-y)


            if   t1 ≤ R ≤ t2

                return 1/R

            else

                return zero(R)

            end

        end

    #result 4d quadrule

    r = 0.0

    for qpT in qpsT

        pT, wT = qpT

        for qpS in qpsS

            pS, wS = qpS

            r += wS * wT * f(pT,pS)

        end

    end

    return r

end

