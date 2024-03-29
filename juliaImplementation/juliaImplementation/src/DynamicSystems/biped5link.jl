# 5-link biped model in 2D, adapted from https://web.eecs.umich.edu/~grizzle/biped_book_web/
# dynamics equation: D(q)q'' + C(q, q') + G(q) = B*u


# confusing notation (pg 177, pg 465): 
# MZ = p^M 
# MYtorso = MZtorso
# XX = I 

p = (g=9.81, 
     Ltorso=0.63, Lfem=0.4, Ltib=0.4,
     Mtorso=12, Mfem=6.8, Mtib=3.2,
     MYtorso=0.24, MZtorso=0.24, 
     MZfem=0.11, MZtib=0.24, 
     XX_torso=0.63, XX_fem=1.33, XX_tib=0.2)

function D_matrix(q)
    D = zeros(5,5)
    D[1,1]=p(13)-4*p(6)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))-2*p(5)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))-2*p(7)*
           p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(14)+2*p(10)*p(4)*cos(2*q(1)-q(5)-q(3))+2*p(6)*p(3)^2-2*p(10)*p(3)+2*
           p(7)*p(4)^2-2*p(11)*p(4)+p(7)*p(3)^2+2*p(6)*p(4)^2+p(5)*p(3)^2+p(5)*p(4)^2;
    D[1,2]=p(7)*p(4)*p(3)*cos(q(2)-q(5)+q(1)-q(3))-p(7)*p(3)^2*cos(-q(2)+q(1))+p(11)*p(3)*cos(q(2)-
           q(4)+q(1)-q(5))-p(10)*p(3)*cos(-q(2)+q(1))+p(10)*p(4)*cos(q(2)-q(5)+q(1)-q(3))-p(11)*p(4)*cos(-q(2)+q(4)+
           q(1)-q(3));
    D[1,3]=2*p(6)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(7)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(5)*p(3)*p(4)*
           cos(2*q(1)-q(5)-q(3))-p(10)*p(4)*cos(2*q(1)-q(5)-q(3))-p(5)*p(4)^2+2*p(11)*p(4)-2*p(7)*p(4)^2-2*p(6)*
           p(4)^2-p(14);
    D[1,4]=-p(11)*(p(3)*cos(q(2)-q(4)+q(1)-q(5))-p(4)*cos(-q(2)+q(4)+q(1)-q(3)));
    D[1,5]=2*p(6)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))-p(7)*p(4)*p(3)*cos(q(2)-q(5)+q(1)-q(3))+p(7)*p(3)^2*
           cos(-q(2)+q(1))+p(8)*p(4)*sin(q(5)+q(1)-q(3))+p(9)*p(4)*cos(q(5)+q(1)-q(3))-p(9)*p(3)*cos(-2*q(5)+q(1))+
           p(5)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(10)*p(3)*cos(-q(2)+q(1))-p(10)*p(4)*cos(q(2)-q(5)+q(1)-q(3))+p(7)*
           p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(8)*p(3)*sin(-2*q(5)+q(1))-p(10)*p(4)*cos(2*q(1)-q(5)-q(3))-p(5)*p(3)^2-
           2*p(6)*p(3)^2+2*p(10)*p(3)-p(7)*p(3)^2-p(13);
    D[2,1]=p(7)*p(4)*p(3)*cos(q(2)-q(5)+q(1)-q(3))-p(7)*p(3)^2*cos(-q(2)+q(1))+p(11)*p(3)*cos(q(2)-
           q(4)+q(1)-q(5))-p(10)*p(3)*cos(-q(2)+q(1))+p(10)*p(4)*cos(q(2)-q(5)+q(1)-q(3))-p(11)*p(4)*cos(-q(2)+q(4)+
           q(1)-q(3));
    D[2,2]=p(7)*p(3)^2-2*p(11)*p(3)*cos(2*q(2)-q(4)-q(5))+p(14)+p(13);
    D[2,3]=-p(4)*(p(7)*p(3)*cos(q(2)-q(5)+q(1)-q(3))+p(10)*cos(q(2)-q(5)+q(1)-q(3))-p(11)*cos(-q(2)+
           q(4)+q(1)-q(3)));
    D[2,4]=-p(14)+p(11)*p(3)*cos(2*q(2)-q(4)-q(5));
    D[2,5]=-p(7)*p(3)^2-p(13)+p(7)*p(3)^2*cos(-q(2)+q(1))-p(11)*p(3)*cos(q(2)-q(4)+q(1)-q(5))+p(10)*
           p(3)*cos(-q(2)+q(1))+p(11)*p(3)*cos(2*q(2)-q(4)-q(5));
    D[3,1]=2*p(6)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(7)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(5)*p(3)*p(4)*
           cos(2*q(1)-q(5)-q(3))-p(10)*p(4)*cos(2*q(1)-q(5)-q(3))-p(5)*p(4)^2+2*p(11)*p(4)-2*p(7)*p(4)^2-2*p(6)*
           p(4)^2-p(14);
    D[3,2]=-p(4)*(p(7)*p(3)*cos(q(2)-q(5)+q(1)-q(3))+p(10)*cos(q(2)-q(5)+q(1)-q(3))-p(11)*cos(-q(2)+
           q(4)+q(1)-q(3)));
    D[3,3]=p(5)*p(4)^2-2*p(11)*p(4)+2*p(7)*p(4)^2+2*p(6)*p(4)^2+p(14);
    D[3,4]=-p(11)*p(4)*cos(-q(2)+q(4)+q(1)-q(3));
    D[3,5]=-p(4)*(2*p(6)*p(3)*cos(2*q(1)-q(5)-q(3))-p(7)*p(3)*cos(q(2)-q(5)+q(1)-q(3))+p(8)*sin(q(5)+
           q(1)-q(3))+p(9)*cos(q(5)+q(1)-q(3))+p(5)*p(3)*cos(2*q(1)-q(5)-q(3))-p(10)*cos(q(2)-q(5)+q(1)-q(3))+p(7)*
           p(3)*cos(2*q(1)-q(5)-q(3))-p(10)*cos(2*q(1)-q(5)-q(3)));
    D[4,1]=-p(11)*(p(3)*cos(q(2)-q(4)+q(1)-q(5))-p(4)*cos(-q(2)+q(4)+q(1)-q(3)));
    D[4,2]=-p(14)+p(11)*p(3)*cos(2*q(2)-q(4)-q(5));
    D[4,3]=-p(11)*p(4)*cos(-q(2)+q(4)+q(1)-q(3));
    D[4,4]=p(14);
    D[4,5]=-p(11)*p(3)*(-cos(q(2)-q(4)+q(1)-q(5))+cos(2*q(2)-q(4)-q(5)));
    D[5,1]=2*p(6)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))-p(7)*p(4)*p(3)*cos(q(2)-q(5)+q(1)-q(3))+p(7)*p(3)^2*
           cos(-q(2)+q(1))+p(8)*p(4)*sin(q(5)+q(1)-q(3))+p(9)*p(4)*cos(q(5)+q(1)-q(3))-p(9)*p(3)*cos(-2*q(5)+q(1))+
           p(5)*p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(10)*p(3)*cos(-q(2)+q(1))-p(10)*p(4)*cos(q(2)-q(5)+q(1)-q(3))+p(7)*
           p(3)*p(4)*cos(2*q(1)-q(5)-q(3))+p(8)*p(3)*sin(-2*q(5)+q(1))-p(10)*p(4)*cos(2*q(1)-q(5)-q(3))-p(5)*p(3)^2-
           2*p(6)*p(3)^2+2*p(10)*p(3)-p(7)*p(3)^2-p(13);
    D[5,2]=-p(7)*p(3)^2-p(13)+p(7)*p(3)^2*cos(-q(2)+q(1))-p(11)*p(3)*cos(q(2)-q(4)+q(1)-q(5))+p(10)*
           p(3)*cos(-q(2)+q(1))+p(11)*p(3)*cos(2*q(2)-q(4)-q(5));
    D[5,3]=-p(4)*(2*p(6)*p(3)*cos(2*q(1)-q(5)-q(3))-p(7)*p(3)*cos(q(2)-q(5)+q(1)-q(3))+p(8)*sin(q(5)+
           q(1)-q(3))+p(9)*cos(q(5)+q(1)-q(3))+p(5)*p(3)*cos(2*q(1)-q(5)-q(3))-p(10)*cos(q(2)-q(5)+q(1)-q(3))+p(7)*
           p(3)*cos(2*q(1)-q(5)-q(3))-p(10)*cos(2*q(1)-q(5)-q(3)));
    D[5,4]=-p(11)*p(3)*(-cos(q(2)-q(4)+q(1)-q(5))+cos(2*q(2)-q(4)-q(5)));
    D[5,5]=-2*p(7)*p(3)^2*cos(-q(2)+q(1))+2*p(9)*p(3)*cos(-2*q(5)+q(1))-2*p(10)*p(3)*cos(-q(2)+q(1))-
           2*p(8)*p(3)*sin(-2*q(5)+q(1))+p(12)+p(5)*p(3)^2+2*p(6)*p(3)^2-2*p(10)*p(3)+2*p(7)*p(3)^2+2*p(13);
    return D
end  

function C_matrix(q,dq)
    C = zeros(5,5)
    C[1,1]=sin(2*q(1)-q(5)-q(3))*p(4)*(2*dq(1)-dq(3)-dq(5))*(2*p(6)*p(3)+p(5)*p(3)+p(7)*p(3)-p(10));
    C[1,2]=-dq(2)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+q(1)-q(3))-dq(2)*p(7)*p(3)^2*sin(-q(2)+q(1))-dq(2)*
           p(11)*p(3)*sin(q(2)-q(4)+q(1)-q(5))-dq(2)*p(10)*p(3)*sin(-q(2)+q(1))-dq(2)*p(10)*p(4)*sin(q(2)-q(5)+q(1)-
           q(3))-p(11)*p(4)*sin(-q(2)+q(4)+q(1)-q(3))*dq(2)+dq(4)*p(11)*p(3)*sin(q(2)-q(4)+q(1)-q(5))+p(11)*p(4)*
           sin(-q(2)+q(4)+q(1)-q(3))*dq(4)+dq(5)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+q(1)-q(3))+dq(5)*p(7)*p(3)^2*sin(-
           q(2)+q(1))+dq(5)*p(10)*p(3)*sin(-q(2)+q(1))+dq(5)*p(10)*p(4)*sin(q(2)-q(5)+q(1)-q(3));
    C[1,3]=-sin(2*q(1)-q(5)-q(3))*p(4)*(dq(1)-dq(3))*(2*p(6)*p(3)+p(5)*p(3)+p(7)*p(3)-p(10));
    C[1,4]=p(11)*(dq(2)-dq(4))*(p(3)*sin(q(2)-q(4)+q(1)-q(5))+p(4)*sin(-q(2)+q(4)+q(1)-q(3)));
    C[1,5]=-2*dq(1)*p(6)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(1)*p(5)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-
           dq(1)*p(7)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))+dq(1)*p(10)*p(4)*sin(2*q(1)-q(5)-q(3))+dq(2)*p(7)*p(4)*p(3)*
           sin(q(2)-q(5)+q(1)-q(3))+dq(2)*p(7)*p(3)^2*sin(-q(2)+q(1))+dq(2)*p(10)*p(3)*sin(-q(2)+q(1))+dq(2)*p(10)*p(4)*
           sin(q(2)-q(5)+q(1)-q(3))+2*dq(5)*p(6)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(5)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+
           q(1)-q(3))+dq(5)*p(8)*p(4)*cos(q(5)+q(1)-q(3))-dq(5)*p(9)*p(4)*sin(q(5)+q(1)-q(3))-dq(5)*p(9)*p(3)*sin(-
           2*q(5)+q(1))+dq(5)*p(5)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(5)*p(10)*p(4)*sin(q(2)-q(5)+q(1)-q(3))+
           dq(5)*p(7)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(5)*p(8)*p(3)*cos(-2*q(5)+q(1))-dq(5)*p(10)*p(4)*sin(2*q(1)-
           q(5)-q(3))-dq(5)*p(7)*p(3)^2*sin(-q(2)+q(1))-dq(5)*p(10)*p(3)*sin(-q(2)+q(1));
    C[2,1]=-dq(1)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+q(1)-q(3))+dq(1)*p(7)*p(3)^2*sin(-q(2)+q(1))-p(11)*
           p(3)*sin(q(2)-q(4)+q(1)-q(5))*dq(1)+dq(1)*p(10)*p(3)*sin(-q(2)+q(1))-dq(1)*p(10)*p(4)*sin(q(2)-q(5)+q(1)-
           q(3))+p(11)*p(4)*sin(-q(2)+q(4)+q(1)-q(3))*dq(1)+dq(3)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+q(1)-q(3))+dq(3)*
           p(10)*p(4)*sin(q(2)-q(5)+q(1)-q(3))-p(11)*p(4)*sin(-q(2)+q(4)+q(1)-q(3))*dq(3)-dq(5)*p(7)*p(3)^2*sin(-
           q(2)+q(1))+p(11)*p(3)*sin(q(2)-q(4)+q(1)-q(5))*dq(5)-dq(5)*p(10)*p(3)*sin(-q(2)+q(1));
    C[2,2]=p(11)*p(3)*sin(2*q(2)-q(4)-q(5))*(2*dq(2)-dq(4)-dq(5));
    C[2,3]=-p(4)*(dq(1)-dq(3))*(-p(7)*p(3)*sin(q(2)-q(5)+q(1)-q(3))-p(10)*sin(q(2)-q(5)+q(1)-q(3))+
           p(11)*sin(-q(2)+q(4)+q(1)-q(3)));
    C[2,4]=-p(11)*p(3)*sin(2*q(2)-q(4)-q(5))*(dq(2)-dq(4));
    C[2,5]=-p(3)*(dq(1)*p(7)*p(3)*sin(-q(2)+q(1))-dq(1)*p(11)*sin(q(2)-q(4)+q(1)-q(5))+dq(1)*p(10)*
           sin(-q(2)+q(1))+p(11)*sin(2*q(2)-q(4)-q(5))*dq(2)+dq(5)*p(11)*sin(q(2)-q(4)+q(1)-q(5))-p(11)*sin(2*q(2)-
           q(4)-q(5))*dq(5)-dq(5)*p(7)*p(3)*sin(-q(2)+q(1))-dq(5)*p(10)*sin(-q(2)+q(1)));
    C[3,1]=-sin(2*q(1)-q(5)-q(3))*p(4)*(dq(1)-dq(5))*(2*p(6)*p(3)+p(5)*p(3)+p(7)*p(3)-p(10));
    C[3,2]=p(4)*(dq(2)*p(7)*p(3)*sin(q(2)-q(5)+q(1)-q(3))+dq(2)*p(10)*sin(q(2)-q(5)+q(1)-q(3))+dq(2)*
           p(11)*sin(-q(2)+q(4)+q(1)-q(3))-dq(4)*p(11)*sin(-q(2)+q(4)+q(1)-q(3))-dq(5)*p(7)*p(3)*sin(q(2)-q(5)+q(1)-
           q(3))-dq(5)*p(10)*sin(q(2)-q(5)+q(1)-q(3)));
    C[3,3]=0;
    C[3,4]=-p(11)*p(4)*sin(-q(2)+q(4)+q(1)-q(3))*(dq(2)-dq(4));
    C[3,5]=p(4)*(2*dq(1)*p(6)*p(3)*sin(2*q(1)-q(5)-q(3))+dq(1)*p(7)*p(3)*sin(2*q(1)-q(5)-q(3))+dq(1)*
           p(5)*p(3)*sin(2*q(1)-q(5)-q(3))-dq(1)*p(10)*sin(2*q(1)-q(5)-q(3))-dq(2)*p(7)*p(3)*sin(q(2)-q(5)+q(1)-
           q(3))-dq(2)*p(10)*sin(q(2)-q(5)+q(1)-q(3))-2*dq(5)*p(6)*p(3)*sin(2*q(1)-q(5)-q(3))+dq(5)*p(7)*p(3)*
           sin(q(2)-q(5)+q(1)-q(3))-dq(5)*p(8)*cos(q(5)+q(1)-q(3))+dq(5)*p(9)*sin(q(5)+q(1)-q(3))-dq(5)*p(5)*p(3)*sin(2*
           q(1)-q(5)-q(3))+dq(5)*p(10)*sin(q(2)-q(5)+q(1)-q(3))-dq(5)*p(7)*p(3)*sin(2*q(1)-q(5)-q(3))+dq(5)*p(10)*
           sin(2*q(1)-q(5)-q(3)));
    C[4,1]=p(11)*(dq(1)*p(3)*sin(q(2)-q(4)+q(1)-q(5))-dq(1)*p(4)*sin(-q(2)+q(4)+q(1)-q(3))+dq(3)*p(4)*
           sin(-q(2)+q(4)+q(1)-q(3))-dq(5)*p(3)*sin(q(2)-q(4)+q(1)-q(5)));
    C[4,2]=-p(11)*p(3)*sin(2*q(2)-q(4)-q(5))*(dq(2)-dq(5));
    C[4,3]=p(11)*p(4)*sin(-q(2)+q(4)+q(1)-q(3))*(dq(1)-dq(3));
    C[4,4]=0;
    C[4,5]=-p(11)*p(3)*(dq(1)*sin(q(2)-q(4)+q(1)-q(5))-sin(2*q(2)-q(4)-q(5))*dq(2)-dq(5)*sin(q(2)-
           q(4)+q(1)-q(5))+sin(2*q(2)-q(4)-q(5))*dq(5));
    C[5,1]=-2*dq(1)*p(6)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))+dq(1)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+q(1)-
           q(3))-dq(1)*p(7)*p(3)^2*sin(-q(2)+q(1))+dq(1)*p(8)*p(4)*cos(q(5)+q(1)-q(3))-dq(1)*p(9)*p(4)*sin(q(5)+q(1)-
           q(3))+dq(1)*p(9)*p(3)*sin(-2*q(5)+q(1))-dq(1)*p(5)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(1)*p(10)*p(3)*sin(-
           q(2)+q(1))+dq(1)*p(10)*p(4)*sin(q(2)-q(5)+q(1)-q(3))-dq(1)*p(7)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))+dq(1)*
           p(8)*p(3)*cos(-2*q(5)+q(1))+dq(1)*p(10)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(3)*p(7)*p(4)*p(3)*sin(q(2)-q(5)+
           q(1)-q(3))-dq(3)*p(8)*p(4)*cos(q(5)+q(1)-q(3))+dq(3)*p(9)*p(4)*sin(q(5)+q(1)-q(3))-dq(3)*p(10)*p(4)*
           sin(q(2)-q(5)+q(1)-q(3))+2*dq(3)*p(6)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))+dq(3)*p(5)*p(3)*p(4)*sin(2*q(1)-q(5)-
           q(3))+dq(3)*p(7)*p(3)*p(4)*sin(2*q(1)-q(5)-q(3))-dq(3)*p(10)*p(4)*sin(2*q(1)-q(5)-q(3))+dq(5)*p(7)*p(3)^2*
           sin(-q(2)+q(1))-dq(5)*p(9)*p(3)*sin(-2*q(5)+q(1))+dq(5)*p(10)*p(3)*sin(-q(2)+q(1))-dq(5)*p(8)*p(3)*cos(-
           2*q(5)+q(1));
    C[5,2]=p(3)*(dq(2)*p(7)*p(3)*sin(-q(2)+q(1))+dq(2)*p(11)*sin(q(2)-q(4)+q(1)-q(5))+dq(2)*p(10)*
           sin(-q(2)+q(1))-p(11)*sin(2*q(2)-q(4)-q(5))*dq(2)-dq(4)*p(11)*sin(q(2)-q(4)+q(1)-q(5))+p(11)*sin(2*q(2)-
           q(4)-q(5))*dq(4)-dq(5)*p(7)*p(3)*sin(-q(2)+q(1))-dq(5)*p(10)*sin(-q(2)+q(1)));
    C[5,3]=p(4)*(dq(1)-dq(3))*(2*p(6)*p(3)*sin(2*q(1)-q(5)-q(3))-p(7)*p(3)*sin(q(2)-q(5)+q(1)-q(3))-
           p(8)*cos(q(5)+q(1)-q(3))+p(9)*sin(q(5)+q(1)-q(3))+p(5)*p(3)*sin(2*q(1)-q(5)-q(3))-p(10)*sin(q(2)-q(5)+
           q(1)-q(3))+p(7)*p(3)*sin(2*q(1)-q(5)-q(3))-p(10)*sin(2*q(1)-q(5)-q(3)));
    C[5,4]=-p(11)*p(3)*(dq(2)-dq(4))*(sin(q(2)-q(4)+q(1)-q(5))-sin(2*q(2)-q(4)-q(5)));
    C[5,5]=p(3)*(dq(1)*p(7)*p(3)*sin(-q(2)+q(1))-dq(1)*p(9)*sin(-2*q(5)+q(1))+dq(1)*p(10)*sin(-q(2)+
           q(1))-dq(1)*p(8)*cos(-2*q(5)+q(1))-dq(2)*p(7)*p(3)*sin(-q(2)+q(1))-dq(2)*p(10)*sin(-q(2)+q(1))+2*dq(5)*
           p(9)*sin(-2*q(5)+q(1))+2*dq(5)*p(8)*cos(-2*q(5)+q(1)));
    return C
end

function B_matrix()
    B = zeros(5,4)
    B[1,1]=1;     B[1,2]=0;     B[1,3]=-2;    B[1,4]=0;
    B[2,1]=0;     B[2,2]=1;     B[2,3]=0;     B[2,4]=-2;
    B[3,1]=0;     B[3,2]=0;     B[3,3]=1;     B[3,4]=0;
    B[4,1]=0;     B[4,2]=0;     B[4,3]=0;     B[4,4]=1;
    B[5,1]=-2;    B[5,2]=-2;    B[5,3]=1;     B[5,4]=1;
    return B
end

function G_vector(q)
    G[1,1]=p(1)*(p(3)*sin(q(1)-q(5))*p(5)+p(4)*sin(q(1)-q(3))*p(5)+2*p(3)*sin(q(1)-q(5))*p(6)+2*p(4)*
           sin(q(1)-q(3))*p(6)-p(10)*sin(q(1)-q(5))+2*p(4)*sin(q(1)-q(3))*p(7)-p(11)*sin(q(1)-q(3))+p(3)*sin(q(1)-q(5))*p(7));
    G[2,1]=-p(1)*(p(10)*sin(q(2)-q(5))+p(3)*sin(q(2)-q(5))*p(7)+p(11)*sin(q(2)-q(4)));
    G[3,1]=-p(1)*sin(q(1)-q(3))*(p(4)*p(5)+2*p(4)*p(6)+2*p(4)*p(7)-p(11));
    G[4,1]=p(1)*p(11)*sin(q(2)-q(4));
    G[5,1]=p(1)*(-p(3)*sin(q(1)-q(5))*p(5)-sin(q(5))*p(9)+cos(q(5))*p(8)-2*p(3)*sin(q(1)-q(5))*p(6)+
           p(10)*sin(q(1)-q(5))+p(10)*sin(q(2)-q(5))-p(3)*sin(q(1)-q(5))*p(7)+p(3)*sin(q(2)-q(5))*p(7));
    return G
end
  
function dynamics(x, u, t)
    q, dq = x[1:5], x[6:10]
    D = D_matrix(q)
    C = C_matrix(q, dq)
    B = B_matrix()
    G = G_vector(q)
    ddq = D \ (B*u - C*dq - G)
    return [dq; ddq]
end