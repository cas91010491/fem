clear all
% close all
clc

alpha_init=1;
c_par=1e-4;
c_par2=0.9;
r_par=0.5;
tol=0;
tol2=1e-15;

tol_RM=1e-12;

Young=1;
nu=0.3;
My0=0.01;
Hardening_modulus=0.05;
Hardening_exponent=1;

nx=5;
ny=5;
nz=5;

cl=0.25;

reg_bc=[-1e-8 (nx-1)*cl+1e-8 -1e-8 (ny-1)*cl+1e-8 -1e-8 1e-8;
    -1e-8 (nx-1)*cl+1e-8 -1e-8 (ny-1)*cl+1e-8 (nz-1)*cl-1e-8 (nz-1)*cl+1e-8];
% reg_bc=[-1e-8 1e-8 -1e-8 1e-8 -1e-8 1e-8;
%     cl-1e-8 cl+1e-8 -1e-8 1e-8 -1e-8 1e-8];

%0/1: dis/force | 0/1/2/3/4/5/6: x/y/z/x&y/x&z/y&z/x&y&z | value x | value
%of y | value of z
bc_tot=[0 6 0 0 0;
    0 6 0.025 0.025 0.025];

% bc_tot=[0 6 0 0 0;
%     0 6 0.25 0 0];

% graph_points=[0 0;1 1;2 0];
% time_steps=[100;100];

graph_points=[0 0;1 1];
time_steps=[100];

%--------------------------------------------------------
%meshing

LEN_node=nx*ny*nz;
FE_node_coord=zeros(LEN_node,3);

LEN_el=(nx-1)*(ny-1)*(nz-1);
con_mat=zeros(LEN_el,8);


cnt=1;
cnt2=1;
cntlb=0;
for k=1:nz
    for j=1:ny
        for i=1:nx

            FE_node_coord(cnt,:)=[(i-1)*cl (j-1)*cl (k-1)*cl];
            cnt=cnt+1;

            cntlb=cntlb+1;

            if i==nx || j==ny || k==nz
            else
                con_mat(cnt2,:)=[cntlb cntlb+1 cntlb+1+nx cntlb+nx cntlb+nx*ny cntlb+nx*ny+1 cntlb+nx*ny+nx+1 cntlb+nx*ny+nx];
                cnt2=cnt2+1;
            end

        end
    end

end



% figure(1)
% hold on
% plot3(FE_node_coord(:,1),FE_node_coord(:,2),FE_node_coord(:,3),'bo')
%
%
% figure(1)
% hold on
% for i=1:length(con_mat(:,1))
%
%     nr1=con_mat(i,1);
%     nr2=con_mat(i,2);
%     nr3=con_mat(i,3);
%     nr4=con_mat(i,4);
%     nr5=con_mat(i,5);
%     nr6=con_mat(i,6);
%     nr7=con_mat(i,7);
%     nr8=con_mat(i,8);
%
%     x1=FE_node_coord(nr1,1);
%     y1=FE_node_coord(nr1,2);
%     z1=FE_node_coord(nr1,3);
%
%     x2=FE_node_coord(nr2,1);
%     y2=FE_node_coord(nr2,2);
%     z2=FE_node_coord(nr2,3);
%
%     x3=FE_node_coord(nr3,1);
%     y3=FE_node_coord(nr3,2);
%     z3=FE_node_coord(nr3,3);
%
%     x4=FE_node_coord(nr4,1);
%     y4=FE_node_coord(nr4,2);
%     z4=FE_node_coord(nr4,3);
%
%     x5=FE_node_coord(nr5,1);
%     y5=FE_node_coord(nr5,2);
%     z5=FE_node_coord(nr5,3);
%
%     x6=FE_node_coord(nr6,1);
%     y6=FE_node_coord(nr6,2);
%     z6=FE_node_coord(nr6,3);
%
%     x7=FE_node_coord(nr7,1);
%     y7=FE_node_coord(nr7,2);
%     z7=FE_node_coord(nr7,3);
%
%     x8=FE_node_coord(nr8,1);
%     y8=FE_node_coord(nr8,2);
%     z8=FE_node_coord(nr8,3);
%
%
%     plot3([x1;x2],[y1;y2],[z1;z2],'r')
%     plot3([x1;x4],[y1;y4],[z1;z4],'r')
%     plot3([x3;x2],[y3;y2],[z3;z2],'r')
%     plot3([x3;x4],[y3;y4],[z3;z4],'r')
%
%     plot3([x5;x6],[y5;y6],[z5;z6],'r')
%     plot3([x5;x8],[y5;y8],[z5;z8],'r')
%     plot3([x7;x6],[y7;y6],[z7;z6],'r')
%     plot3([x7;x8],[y7;y8],[z7;z8],'r')
%
%     plot3([x1;x5],[y1;y5],[z1;z5],'r')
%     plot3([x2;x6],[y2;y6],[z2;z6],'r')
%     plot3([x3;x7],[y3;y7],[z3;z7],'r')
%     plot3([x4;x8],[y4;y8],[z4;z8],'r')
%
% end
% axis off

%--------------------------------------------------------

indices_diri=[];
for i=1:length(reg_bc(:,1))

    x1=reg_bc(i,1);
    x2=reg_bc(i,2);
    y1=reg_bc(i,3);
    y2=reg_bc(i,4);
    z1=reg_bc(i,5);
    z2=reg_bc(i,6);


    [a,~]=find(FE_node_coord(:,1)>x1 & FE_node_coord(:,1)<x2 & FE_node_coord(:,2)>y1 & FE_node_coord(:,2)<y2 & FE_node_coord(:,3)>z1 & FE_node_coord(:,3)<z2);

    if bc_tot(i,1)==0

        if bc_tot(i,2)==0

            indices_diri=[indices_diri; a*3-2];

        elseif bc_tot(i,2)==1

            indices_diri=[indices_diri; a*3-1];

        elseif bc_tot(i,2)==2

            indices_diri=[indices_diri; a*3];

        elseif bc_tot(i,2)==3

            indices_diri=[indices_diri; a*3-2;a*3-1 ];

        elseif bc_tot(i,2)==4

            indices_diri=[indices_diri; a*3-2;a*3 ];

        elseif bc_tot(i,2)==5

            indices_diri=[indices_diri; a*3-1;a*3 ];

        elseif bc_tot(i,2)==6

            indices_diri=[indices_diri; a*3-2;a*3-1;a*3 ];

        end

    end

end
indices_diri=sort(indices_diri,'ascend');
indices_neumann=(1:LEN_node*3)';
indices_neumann(indices_diri)=[];
n_var=length(indices_neumann);

incr_tot=sum(time_steps);

if length(time_steps(:,1))>1
    time_steps2=time_steps;
    for i=2:length(time_steps(:,1))
        time_steps2(i)=time_steps2(i-1)+time_steps2(i);
    end
    time_steps=time_steps2;
end

figure(1)
hold on
u=zeros(LEN_node*3,1);
f_ex=zeros(LEN_node*3,1);

EPCUM=zeros(8,LEN_el);
FP_conv=zeros(9,8*LEN_el);
FP_conv(1,:)=ones(LEN_el*8,1);
FP_conv(5,:)=ones(LEN_el*8,1);
FP_conv(9,:)=ones(LEN_el*8,1);



tot_iter=0;
results=zeros(incr_tot,4);
for incr=1:incr_tot

    incr

    [a,~]=find(incr<=time_steps);
    a=a(1);

    if a>1
        ratio=(incr-time_steps(a-1))/(time_steps(a)-time_steps(a-1));
    else
        ratio=(incr)/(time_steps(a));
    end

    bc_cur=bc_tot;
    if a>1
        bc_cur(:,3:5)=bc_cur(:,3:5)*(ratio*(graph_points(a+1,2)-graph_points(a,2))/(graph_points(a+1,1)-graph_points(a,1))+graph_points(a,2));
    else
        bc_cur(:,3:5)=bc_cur(:,3:5)*(ratio*(graph_points(2,2)-graph_points(1,2))/(graph_points(2,1)-graph_points(1,1))+graph_points(1,2));
    end

    for i=1:length(reg_bc(:,1))

        x1=reg_bc(i,1);
        x2=reg_bc(i,2);
        y1=reg_bc(i,3);
        y2=reg_bc(i,4);
        z1=reg_bc(i,5);
        z2=reg_bc(i,6);


        [a,~]=find(FE_node_coord(:,1)>x1 & FE_node_coord(:,1)<x2 & FE_node_coord(:,2)>y1 & FE_node_coord(:,2)<y2 & FE_node_coord(:,3)>z1 & FE_node_coord(:,3)<z2);

        if bc_tot(i,1)==0

            if bc_tot(i,2)==0

                u(a*3-2)=bc_cur(i,3);

            elseif bc_tot(i,2)==1

                u(a*3-1)=bc_cur(i,4);

            elseif bc_tot(i,2)==2

                u(a*3)=bc_cur(i,5);

            elseif bc_tot(i,2)==3

                u(a*3-2)=bc_cur(i,3);
                u(a*3-1)=bc_cur(i,4);

            elseif bc_tot(i,2)==4

                u(a*3-2)=bc_cur(i,3);
                u(a*3)=bc_cur(i,5);

            elseif bc_tot(i,2)==5

                u(a*3-1)=bc_cur(i,4);
                u(a*3)=bc_cur(i,5);

            elseif bc_tot(i,2)==6

                u(a*3-2)=bc_cur(i,3);
                u(a*3-1)=bc_cur(i,4);
                u(a*3)=bc_cur(i,5);

            end

        else

            if bc_tot(i,2)==0

                f_ex(a*3-2)=bc_cur(i,3);

            elseif bc_tot(i,2)==1

                f_ex(a*3-1)=bc_cur(i,4);

            elseif bc_tot(i,2)==2

                f_ex(a*3)=bc_cur(i,5);

            elseif bc_tot(i,2)==3

                f_ex(a*3-2)=bc_cur(i,3);
                f_ex(a*3-1)=bc_cur(i,4);

            elseif bc_tot(i,2)==4

                f_ex(a*3-2)=bc_cur(i,3);
                f_ex(a*3)=bc_cur(i,5);

            elseif bc_tot(i,2)==5

                f_ex(a*3-1)=bc_cur(i,4);
                f_ex(a*3)=bc_cur(i,5);

            elseif bc_tot(i,2)==6

                f_ex(a*3-2)=bc_cur(i,3);
                f_ex(a*3-1)=bc_cur(i,4);
                f_ex(a*3)=bc_cur(i,5);

            end

        end

    end




    [m_new,f_2,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,u,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM);

    K_new_inv=inv(eye(n_var));
    f_new=zeros(n_var,1);
    m_old=1e100;
    cnt=0;
    while (m_old-m_new)>tol*abs(m_old)

        cnt=cnt+1;

        m_old=m_new;
        f_old=f_new;

        f_new=f_2;

        K_old_inv=K_new_inv;
        delta_f=f_new-f_old;

        if cnt==1
            h_new=-K_old_inv*f_new;
        else
            K_new_inv=K_old_inv+(delta_u'*delta_f+delta_f'*K_old_inv*delta_f)*delta_u*delta_u'/(delta_u'*delta_f)^2-(K_old_inv*delta_f*delta_u'+delta_u*delta_f'*K_old_inv)/(delta_u'*delta_f);
            h_new=-K_new_inv*f_new;
        end

        

        alpha3=alpha_init;
        ux=u;
        ux(indices_neumann)=u(indices_neumann)+alpha3*h_new;
        [m_3,f_3,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,ux,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM);


        signal1=0;
        if m_3<=m_new+c_par*alpha3*h_new'*f_new && h_new'*f_3>=c_par2*h_new'*f_new
            signal1=1;
        end

        cnt2a=0;
        while m_3<m_new+c_par*alpha3*h_new'*f_new && signal1==0

            alpha3=alpha3/r_par;

            ux=u;
            ux(indices_neumann)=u(indices_neumann)+alpha3*h_new;
            [m_3,f_3,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,ux,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM);



            cnt2a=cnt2a+1;

            if m_3<=m_new+c_par*alpha3*h_new'*f_new && h_new'*f_3>=c_par2*h_new'*f_new

                signal1=1;

            end

        end

        if signal1==0
            alpha1=0;
            m_1=m_new;
            f_1=f_new;

            alpha2=alpha3/2;
            ux=u;
            ux(indices_neumann)=u(indices_neumann)+alpha2*h_new;
            [m_2,f_2,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,ux,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM);





            signal2=0;
            cnt2=0;
            while signal2==0

                cnt2=cnt2+1;

                if alpha3-alpha1<tol2

                    signal2=1;
                    m_2=m_new;
                    f_2=f_new;

                elseif m_2>m_new+c_par*alpha2*h_new'*f_new

                    alpha3=alpha2;
                    m_3=m_2;
                    f_3=f_2;

                    alpha2=0.5*(alpha1+alpha2);

                    ux=u;
                    ux(indices_neumann)=u(indices_neumann)+alpha2*h_new;
                    [m_2,f_2,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,ux,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM);



                elseif h_new'*f_2<c_par2*h_new'*f_new

                    alpha1=alpha2;
                    m_1=m_2;
                    f_1=f_2;

                    alpha2=0.5*(alpha2+alpha3);

                    ux=u;
                    ux(indices_neumann)=u(indices_neumann)+alpha2*h_new;
                    [m_2,f_2,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,ux,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM);



                else

                    signal2=1;

                end

            end
        end

        delta_u=ux(indices_neumann)-u(indices_neumann);

        u=ux;

        if signal1==1
            m_new=m_3;
            f_2=f_3;
        else
            m_new=m_2;
        end

    end

    [a,~]=find(FE_node_coord(:,1)>-1e-8 & FE_node_coord(:,1)<(nx-1)*cl+1e-8 & FE_node_coord(:,2)>-1e-8 & FE_node_coord(:,2)<(ny-1)*cl+1e-8 & FE_node_coord(:,3)>(nz-1)*cl-1e-8 & FE_node_coord(:,3)<(nz-1)*cl+1e-8);

    FP_conv=FP_temp;
    EPCUM=EPCUM+DELTA_EPCUM;

    f_in = f_only(LEN_node,u,LEN_el,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM);
    %     plot(incr,norm([f_in(a*3-2);f_in(a*3-1);f_in(a*3)]),'bo')
    %     plot(norm([u(a*3-2);u(a*3-1);u(a*3)]),norm([f_in(a*3-2);f_in(a*3-1);f_in(a*3)]),'bo')
    results(incr,:)=[norm([u(a(1)*3-2);u(a(1)*3-1);u(a(1)*3)]) norm([sum(f_in(a*3-2));sum(f_in(a*3-1));sum(f_in(a*3))]) norm(f_2) cnt];





end

figure(1)
hold on
plot([0;results(:,1)],[0;results(:,2)],'b')

figure(2)
hold on
plot([0;results(:,1)],[0;results(:,3)],'k')

figure(3)
hold on
plot([0;results(:,1)],[0;results(:,4)],'m')


function [m,f,FP_temp,DELTA_EPCUM] = m_and_f(f_ex,u,LEN_el,indices_neumann,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM,tol_RM)

KSI=1/sqrt(3)*[-1 -1 -1;
    +1 -1 -1;
    +1 +1 -1;
    -1 +1 -1;
    -1 -1 +1;
    +1 -1 +1;
    +1 +1 +1;
    -1 +1 +1];

Ca=Young/(2+2*nu);
Cb=Young*nu/((1+nu)*(1-2*nu));

DELTA_EPCUM=zeros(8,LEN_el);
FP_temp=zeros(9,8*LEN_el);

m=-f_ex'*u;
f_in=-f_ex;
for ii=1:length(con_mat(:,1))

    nr1=con_mat(ii,1);
    nr2=con_mat(ii,2);
    nr3=con_mat(ii,3);
    nr4=con_mat(ii,4);
    nr5=con_mat(ii,5);
    nr6=con_mat(ii,6);
    nr7=con_mat(ii,7);
    nr8=con_mat(ii,8);

    x1=FE_node_coord(nr1,1);
    y1=FE_node_coord(nr1,2);
    z1=FE_node_coord(nr1,3);

    x2=FE_node_coord(nr2,1);
    y2=FE_node_coord(nr2,2);
    z2=FE_node_coord(nr2,3);

    x3=FE_node_coord(nr3,1);
    y3=FE_node_coord(nr3,2);
    z3=FE_node_coord(nr3,3);

    x4=FE_node_coord(nr4,1);
    y4=FE_node_coord(nr4,2);
    z4=FE_node_coord(nr4,3);

    x5=FE_node_coord(nr5,1);
    y5=FE_node_coord(nr5,2);
    z5=FE_node_coord(nr5,3);

    x6=FE_node_coord(nr6,1);
    y6=FE_node_coord(nr6,2);
    z6=FE_node_coord(nr6,3);

    x7=FE_node_coord(nr7,1);
    y7=FE_node_coord(nr7,2);
    z7=FE_node_coord(nr7,3);

    x8=FE_node_coord(nr8,1);
    y8=FE_node_coord(nr8,2);
    z8=FE_node_coord(nr8,3);


    a_mat=[nr1*3-2;nr1*3-1;nr1*3;nr2*3-2;nr2*3-1;nr2*3;nr3*3-2;nr3*3-1;nr3*3;nr4*3-2;nr4*3-1;nr4*3;nr5*3-2;nr5*3-1;nr5*3;nr6*3-2;nr6*3-1;nr6*3;nr7*3-2;nr7*3-1;nr7*3;nr8*3-2;nr8*3-1;nr8*3];

    xyz_all=[x1 y1 z1;x2 y2 z2;x3 y3 z3;x4 y4 z4;x5 y5 z5;x6 y6 z6;x7 y7 z7;x8 y8 z8];


    u_all=u(a_mat(1:3:24-2));
    v_all=u(a_mat(2:3:24-1));
    w_all=u(a_mat(3:3:24));

    dEdU=zeros(24,1);

    for jj=1:8

        ksi1=KSI(jj,1);
        ksi2=KSI(jj,2);
        ksi3=KSI(jj,3);

        dN1dksi1=-((ksi2 - 1)*(ksi3 - 1))/8;
        dN1dksi2=-((ksi1 - 1)*(ksi3 - 1))/8;
        dN1dksi3=-((ksi1 - 1)*(ksi2 - 1))/8;
        dN2dksi1=((ksi2 - 1)*(ksi3 - 1))/8;
        dN2dksi2=((ksi1 + 1)*(ksi3 - 1))/8;
        dN2dksi3=((ksi1 + 1)*(ksi2 - 1))/8;
        dN3dksi1=-((ksi2 + 1)*(ksi3 - 1))/8;
        dN3dksi2=-((ksi1 + 1)*(ksi3 - 1))/8;
        dN3dksi3=-((ksi1 + 1)*(ksi2 + 1))/8;
        dN4dksi1=((ksi2 + 1)*(ksi3 - 1))/8;
        dN4dksi2=((ksi1 - 1)*(ksi3 - 1))/8;
        dN4dksi3=((ksi1 - 1)*(ksi2 + 1))/8;
        dN5dksi1=((ksi2 - 1)*(ksi3 + 1))/8;
        dN5dksi2=((ksi1 - 1)*(ksi3 + 1))/8;
        dN5dksi3=((ksi1 - 1)*(ksi2 - 1))/8;
        dN6dksi1=-((ksi2 - 1)*(ksi3 + 1))/8;
        dN6dksi2=-((ksi1 + 1)*(ksi3 + 1))/8;
        dN6dksi3=-((ksi1 + 1)*(ksi2 - 1))/8;
        dN7dksi1=((ksi2 + 1)*(ksi3 + 1))/8;
        dN7dksi2=((ksi1 + 1)*(ksi3 + 1))/8;
        dN7dksi3=((ksi1 + 1)*(ksi2 + 1))/8;
        dN8dksi1=-((ksi2 + 1)*(ksi3 + 1))/8;
        dN8dksi2=-((ksi1 - 1)*(ksi3 + 1))/8;
        dN8dksi3=-((ksi1 - 1)*(ksi2 + 1))/8;

        dNdksi=[dN1dksi1 dN2dksi1 dN3dksi1 dN4dksi1 dN5dksi1 dN6dksi1 dN7dksi1 dN8dksi1;
            dN1dksi2 dN2dksi2 dN3dksi2 dN4dksi2 dN5dksi2 dN6dksi2 dN7dksi2 dN8dksi2;
            dN1dksi3 dN2dksi3 dN3dksi3 dN4dksi3 dN5dksi3 dN6dksi3 dN7dksi3 dN8dksi3];

        Jmat=dNdksi*xyz_all;

        detJ=det(Jmat);
        invJmat=inv(Jmat);

        dNdXYZ=[invJmat*dNdksi(:,1);
            invJmat*dNdksi(:,2);
            invJmat*dNdksi(:,3);
            invJmat*dNdksi(:,4);
            invJmat*dNdksi(:,5);
            invJmat*dNdksi(:,6);
            invJmat*dNdksi(:,7);
            invJmat*dNdksi(:,8)];

        F=[dNdXYZ(1:3:24-2)'*u_all dNdXYZ(2:3:24-1)'*u_all dNdXYZ(3:3:24)'*u_all;
            dNdXYZ(1:3:24-2)'*v_all dNdXYZ(2:3:24-1)'*v_all dNdXYZ(3:3:24)'*v_all;
            dNdXYZ(1:3:24-2)'*w_all dNdXYZ(2:3:24-1)'*w_all dNdXYZ(3:3:24)'*w_all]+eye(3);

        Fp=reshape(FP_conv(:,(ii-1)*8+jj),3,3)';

        epcum=EPCUM(jj,ii);

        Fe=F*inv(Fp);

        P=Ca*(Fe-inv(Fe)')+Cb*log(det(Fe))*inv(Fe)';

        Mdev=Fe'*P-trace(Fe'*P)/3*eye(3);

        yieldfunction=sqrt(1.5*(Mdev(1,1)^2+Mdev(2,2)^2+Mdev(3,3)^2+2*Mdev(1,2)^2+2*Mdev(1,3)^2+2*Mdev(2,3)^2))-My0-Hardening_modulus*epcum^Hardening_exponent;

        if yieldfunction<0

            depcum=0;
            depcum_last=0;
            Fp_temp=Fp;

        else

            if epcum==0 && Hardening_exponent~=1

                depcum_last=1e-10;

                dyieldfunctiondM=1.5*Mdev/sqrt(1.5*(Mdev(1,1)^2+Mdev(2,2)^2+Mdev(3,3)^2+2*Mdev(1,2)^2+2*Mdev(1,3)^2+2*Mdev(2,3)^2));

                [eigenvectors,eigenvalues]=eig(dyieldfunctiondM);
                Fp_temp=(exp(depcum_last*eigenvalues(1,1))*eigenvectors(:,1)*eigenvectors(:,1)'+exp(depcum_last*eigenvalues(2,2))*eigenvectors(:,2)*eigenvectors(:,2)'+exp(depcum_last*eigenvalues(3,3))*eigenvectors(:,3)*eigenvectors(:,3)')*Fp;

                Fe=F*inv(Fp_temp);

                P=Ca*(Fe-inv(Fe)')+Cb*log(det(Fe))*inv(Fe)';

                Mdev=Fe'*P-trace(Fe'*P)/3*eye(3);

                yieldfunction=sqrt(1.5*(Mdev(1,1)^2+Mdev(2,2)^2+Mdev(3,3)^2+2*Mdev(1,2)^2+2*Mdev(1,3)^2+2*Mdev(2,3)^2))-My0-Hardening_modulus*(epcum+depcum_last)^Hardening_exponent;

            else

                depcum_last=0;

                dyieldfunctiondM=1.5*Mdev/sqrt(1.5*(Mdev(1,1)^2+Mdev(2,2)^2+Mdev(3,3)^2+2*Mdev(1,2)^2+2*Mdev(1,3)^2+2*Mdev(2,3)^2));

                [eigenvectors,eigenvalues]=eig(dyieldfunctiondM);

                Fp_temp=Fp;

            end

            F11=F(1,1);
            F12=F(1,2);
            F13=F(1,3);
            F21=F(2,1);
            F22=F(2,2);
            F23=F(2,3);
            F31=F(3,1);
            F32=F(3,2);
            F33=F(3,3);

            iter_RM=0;
            while abs(yieldfunction)>tol_RM

                iter_RM=iter_RM+1;

                P11=P(1,1);
                P12=P(1,2);
                P13=P(1,3);
                P21=P(2,1);
                P22=P(2,2);
                P23=P(2,3);
                P31=P(3,1);
                P32=P(3,2);
                P33=P(3,3);

                Fe11=Fe(1,1);
                Fe12=Fe(1,2);
                Fe13=Fe(1,3);
                Fe21=Fe(2,1);
                Fe22=Fe(2,2);
                Fe23=Fe(2,3);
                Fe31=Fe(3,1);
                Fe32=Fe(3,2);
                Fe33=Fe(3,3);
                lnJe=log(det(Fe));

                invFe=inv(Fe);
                invFe11=invFe(1,1);
                invFe12=invFe(1,2);
                invFe13=invFe(1,3);
                invFe21=invFe(2,1);
                invFe22=invFe(2,2);
                invFe23=invFe(2,3);
                invFe31=invFe(3,1);
                invFe32=invFe(3,2);
                invFe33=invFe(3,3);

                invFp_temp=inv(Fp_temp);
                invFp_temp11=invFp_temp(1,1);
                invFp_temp12=invFp_temp(1,2);
                invFp_temp13=invFp_temp(1,3);
                invFp_temp21=invFp_temp(2,1);
                invFp_temp22=invFp_temp(2,2);
                invFp_temp23=invFp_temp(2,3);
                invFp_temp31=invFp_temp(3,1);
                invFp_temp32=invFp_temp(3,2);
                invFp_temp33=invFp_temp(3,3);

                dyieldfunctiondM11=dyieldfunctiondM(1,1);
                dyieldfunctiondM12=dyieldfunctiondM(1,2);
                dyieldfunctiondM13=dyieldfunctiondM(1,3);
                dyieldfunctiondM22=dyieldfunctiondM(2,2);
                dyieldfunctiondM23=dyieldfunctiondM(2,3);
                dyieldfunctiondM33=dyieldfunctiondM(3,3);

                C1111=Ca + Cb*invFe11^2 + invFe11^2*(Ca - Cb*lnJe);
                C1112=invFe11*invFe21*(Ca + Cb - Cb*lnJe);
                C1113=invFe11*invFe31*(Ca + Cb - Cb*lnJe);
                C1121=invFe11*invFe12*(Ca + Cb - Cb*lnJe);
                C1122=Cb*invFe11*invFe22 + invFe12*invFe21*(Ca - Cb*lnJe);
                C1123=Cb*invFe11*invFe32 + invFe12*invFe31*(Ca - Cb*lnJe);
                C1131=invFe11*invFe13*(Ca + Cb - Cb*lnJe);
                C1132=Cb*invFe11*invFe23 + invFe13*invFe21*(Ca - Cb*lnJe);
                C1133=Cb*invFe11*invFe33 + invFe13*invFe31*(Ca - Cb*lnJe);
                C1212=Cb*invFe12*invFe21 + invFe11*invFe22*(Ca - Cb*lnJe);
                C1213=Cb*invFe12*invFe31 + invFe11*invFe32*(Ca - Cb*lnJe);
                C1221=Ca + Cb*invFe12^2 + invFe12^2*(Ca - Cb*lnJe);
                C1222=invFe12*invFe22*(Ca + Cb - Cb*lnJe);
                C1223=invFe12*invFe32*(Ca + Cb - Cb*lnJe);
                C1231=invFe12*invFe13*(Ca + Cb - Cb*lnJe);
                C1232=Cb*invFe12*invFe23 + invFe13*invFe22*(Ca - Cb*lnJe);
                C1233=Cb*invFe12*invFe33 + invFe13*invFe32*(Ca - Cb*lnJe);
                C1312=Cb*invFe13*invFe21 + invFe11*invFe23*(Ca - Cb*lnJe);
                C1313=Cb*invFe13*invFe31 + invFe11*invFe33*(Ca - Cb*lnJe);
                C1322=Cb*invFe13*invFe22 + invFe12*invFe23*(Ca - Cb*lnJe);
                C1323=Cb*invFe13*invFe32 + invFe12*invFe33*(Ca - Cb*lnJe);
                C1331=Ca + Cb*invFe13^2 + invFe13^2*(Ca - Cb*lnJe);
                C1332=invFe13*invFe23*(Ca + Cb - Cb*lnJe);
                C1333=invFe13*invFe33*(Ca + Cb - Cb*lnJe);
                C2112=Ca + Cb*invFe21^2 + invFe21^2*(Ca - Cb*lnJe);
                C2113=invFe21*invFe31*(Ca + Cb - Cb*lnJe);
                C2122=invFe21*invFe22*(Ca + Cb - Cb*lnJe);
                C2123=Cb*invFe21*invFe32 + invFe22*invFe31*(Ca - Cb*lnJe);
                C2132=invFe21*invFe23*(Ca + Cb - Cb*lnJe);
                C2133=Cb*invFe21*invFe33 + invFe23*invFe31*(Ca - Cb*lnJe);
                C2213=Cb*invFe22*invFe31 + invFe21*invFe32*(Ca - Cb*lnJe);
                C2222=Ca + Cb*invFe22^2 + invFe22^2*(Ca - Cb*lnJe);
                C2223=invFe22*invFe32*(Ca + Cb - Cb*lnJe);
                C2232=invFe22*invFe23*(Ca + Cb - Cb*lnJe);
                C2233=Cb*invFe22*invFe33 + invFe23*invFe32*(Ca - Cb*lnJe);
                C2313=Cb*invFe23*invFe31 + invFe21*invFe33*(Ca - Cb*lnJe);
                C2323=Cb*invFe23*invFe32 + invFe22*invFe33*(Ca - Cb*lnJe);
                C2332=Ca + Cb*invFe23^2 + invFe23^2*(Ca - Cb*lnJe);
                C2333=invFe23*invFe33*(Ca + Cb - Cb*lnJe);
                C3113=Ca + Cb*invFe31^2 + invFe31^2*(Ca - Cb*lnJe);
                C3123=invFe31*invFe32*(Ca + Cb - Cb*lnJe);
                C3133=invFe31*invFe33*(Ca + Cb - Cb*lnJe);
                C3223=Ca + Cb*invFe32^2 + invFe32^2*(Ca - Cb*lnJe);
                C3233=invFe32*invFe33*(Ca + Cb - Cb*lnJe);
                C3333=Ca + Cb*invFe33^2 + invFe33^2*(Ca - Cb*lnJe);

                dMCdFe1111=P11 + C1111*Fe11 + C1121*Fe21 + C1131*Fe31;
                dMCdFe1112=C1112*Fe11 + C1212*Fe21 + C1312*Fe31;
                dMCdFe1113=C1113*Fe11 + C1213*Fe21 + C1313*Fe31;
                dMCdFe1121=P21 + C1121*Fe11 + C1221*Fe21 + C1231*Fe31;
                dMCdFe1122=C1122*Fe11 + C1222*Fe21 + C1322*Fe31;
                dMCdFe1123=C1123*Fe11 + C1223*Fe21 + C1323*Fe31;
                dMCdFe1131=P31 + C1131*Fe11 + C1231*Fe21 + C1331*Fe31;
                dMCdFe1132=C1132*Fe11 + C1232*Fe21 + C1332*Fe31;
                dMCdFe1133=C1133*Fe11 + C1233*Fe21 + C1333*Fe31;
                dMCdFe1211=C1111*Fe12 + C1121*Fe22 + C1131*Fe32;
                dMCdFe1212=P11 + C1112*Fe12 + C1212*Fe22 + C1312*Fe32;
                dMCdFe1213=C1113*Fe12 + C1213*Fe22 + C1313*Fe32;
                dMCdFe1221=C1121*Fe12 + C1221*Fe22 + C1231*Fe32;
                dMCdFe1222=P21 + C1122*Fe12 + C1222*Fe22 + C1322*Fe32;
                dMCdFe1223=C1123*Fe12 + C1223*Fe22 + C1323*Fe32;
                dMCdFe1231=C1131*Fe12 + C1231*Fe22 + C1331*Fe32;
                dMCdFe1232=P31 + C1132*Fe12 + C1232*Fe22 + C1332*Fe32;
                dMCdFe1233=C1133*Fe12 + C1233*Fe22 + C1333*Fe32;
                dMCdFe1311=C1111*Fe13 + C1121*Fe23 + C1131*Fe33;
                dMCdFe1312=C1112*Fe13 + C1212*Fe23 + C1312*Fe33;
                dMCdFe1313=P11 + C1113*Fe13 + C1213*Fe23 + C1313*Fe33;
                dMCdFe1321=C1121*Fe13 + C1221*Fe23 + C1231*Fe33;
                dMCdFe1322=C1122*Fe13 + C1222*Fe23 + C1322*Fe33;
                dMCdFe1323=P21 + C1123*Fe13 + C1223*Fe23 + C1323*Fe33;
                dMCdFe1331=C1131*Fe13 + C1231*Fe23 + C1331*Fe33;
                dMCdFe1332=C1132*Fe13 + C1232*Fe23 + C1332*Fe33;
                dMCdFe1333=P31 + C1133*Fe13 + C1233*Fe23 + C1333*Fe33;
                dMCdFe2111=P12 + C1112*Fe11 + C1122*Fe21 + C1132*Fe31;
                dMCdFe2112=C2112*Fe11 + C2122*Fe21 + C2132*Fe31;
                dMCdFe2113=C2113*Fe11 + C2213*Fe21 + C2313*Fe31;
                dMCdFe2121=P22 + C1212*Fe11 + C1222*Fe21 + C1232*Fe31;
                dMCdFe2122=C2122*Fe11 + C2222*Fe21 + C2232*Fe31;
                dMCdFe2123=C2123*Fe11 + C2223*Fe21 + C2323*Fe31;
                dMCdFe2131=P32 + C1312*Fe11 + C1322*Fe21 + C1332*Fe31;
                dMCdFe2132=C2132*Fe11 + C2232*Fe21 + C2332*Fe31;
                dMCdFe2133=C2133*Fe11 + C2233*Fe21 + C2333*Fe31;
                dMCdFe2211=C1112*Fe12 + C1122*Fe22 + C1132*Fe32;
                dMCdFe2212=P12 + C2112*Fe12 + C2122*Fe22 + C2132*Fe32;
                dMCdFe2213=C2113*Fe12 + C2213*Fe22 + C2313*Fe32;
                dMCdFe2221=C1212*Fe12 + C1222*Fe22 + C1232*Fe32;
                dMCdFe2222=P22 + C2122*Fe12 + C2222*Fe22 + C2232*Fe32;
                dMCdFe2223=C2123*Fe12 + C2223*Fe22 + C2323*Fe32;
                dMCdFe2231=C1312*Fe12 + C1322*Fe22 + C1332*Fe32;
                dMCdFe2232=P32 + C2132*Fe12 + C2232*Fe22 + C2332*Fe32;
                dMCdFe2233=C2133*Fe12 + C2233*Fe22 + C2333*Fe32;
                dMCdFe2311=C1112*Fe13 + C1122*Fe23 + C1132*Fe33;
                dMCdFe2312=C2112*Fe13 + C2122*Fe23 + C2132*Fe33;
                dMCdFe2313=P12 + C2113*Fe13 + C2213*Fe23 + C2313*Fe33;
                dMCdFe2321=C1212*Fe13 + C1222*Fe23 + C1232*Fe33;
                dMCdFe2322=C2122*Fe13 + C2222*Fe23 + C2232*Fe33;
                dMCdFe2323=P22 + C2123*Fe13 + C2223*Fe23 + C2323*Fe33;
                dMCdFe2331=C1312*Fe13 + C1322*Fe23 + C1332*Fe33;
                dMCdFe2332=C2132*Fe13 + C2232*Fe23 + C2332*Fe33;
                dMCdFe2333=P32 + C2133*Fe13 + C2233*Fe23 + C2333*Fe33;
                dMCdFe3111=P13 + C1113*Fe11 + C1123*Fe21 + C1133*Fe31;
                dMCdFe3112=C2113*Fe11 + C2123*Fe21 + C2133*Fe31;
                dMCdFe3113=C3113*Fe11 + C3123*Fe21 + C3133*Fe31;
                dMCdFe3121=P23 + C1213*Fe11 + C1223*Fe21 + C1233*Fe31;
                dMCdFe3122=C2213*Fe11 + C2223*Fe21 + C2233*Fe31;
                dMCdFe3123=C3123*Fe11 + C3223*Fe21 + C3233*Fe31;
                dMCdFe3131=P33 + C1313*Fe11 + C1323*Fe21 + C1333*Fe31;
                dMCdFe3132=C2313*Fe11 + C2323*Fe21 + C2333*Fe31;
                dMCdFe3133=C3133*Fe11 + C3233*Fe21 + C3333*Fe31;
                dMCdFe3211=C1113*Fe12 + C1123*Fe22 + C1133*Fe32;
                dMCdFe3212=P13 + C2113*Fe12 + C2123*Fe22 + C2133*Fe32;
                dMCdFe3213=C3113*Fe12 + C3123*Fe22 + C3133*Fe32;
                dMCdFe3221=C1213*Fe12 + C1223*Fe22 + C1233*Fe32;
                dMCdFe3222=P23 + C2213*Fe12 + C2223*Fe22 + C2233*Fe32;
                dMCdFe3223=C3123*Fe12 + C3223*Fe22 + C3233*Fe32;
                dMCdFe3231=C1313*Fe12 + C1323*Fe22 + C1333*Fe32;
                dMCdFe3232=P33 + C2313*Fe12 + C2323*Fe22 + C2333*Fe32;
                dMCdFe3233=C3133*Fe12 + C3233*Fe22 + C3333*Fe32;
                dMCdFe3311=C1113*Fe13 + C1123*Fe23 + C1133*Fe33;
                dMCdFe3312=C2113*Fe13 + C2123*Fe23 + C2133*Fe33;
                dMCdFe3313=P13 + C3113*Fe13 + C3123*Fe23 + C3133*Fe33;
                dMCdFe3321=C1213*Fe13 + C1223*Fe23 + C1233*Fe33;
                dMCdFe3322=C2213*Fe13 + C2223*Fe23 + C2233*Fe33;
                dMCdFe3323=P23 + C3123*Fe13 + C3223*Fe23 + C3233*Fe33;
                dMCdFe3331=C1313*Fe13 + C1323*Fe23 + C1333*Fe33;
                dMCdFe3332=C2313*Fe13 + C2323*Fe23 + C2333*Fe33;
                dMCdFe3333=P33 + C3133*Fe13 + C3233*Fe23 + C3333*Fe33;

                dFeCdFp_temp1111=-invFp_temp11*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp1112=-invFp_temp21*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp1113=-invFp_temp31*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp1121=-invFp_temp11*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp1122=-invFp_temp21*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp1123=-invFp_temp31*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp1131=-invFp_temp11*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp1132=-invFp_temp21*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp1133=-invFp_temp31*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp1211=-invFp_temp11*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp1212=-invFp_temp21*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp1213=-invFp_temp31*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp1221=-invFp_temp11*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp1222=-invFp_temp21*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp1223=-invFp_temp31*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp1231=-invFp_temp11*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp1232=-invFp_temp21*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp1233=-invFp_temp31*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp1311=-invFp_temp11*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp1312=-invFp_temp21*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp1313=-invFp_temp31*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp1321=-invFp_temp11*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp1322=-invFp_temp21*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp1323=-invFp_temp31*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp1331=-invFp_temp11*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp1332=-invFp_temp21*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp1333=-invFp_temp31*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp2111=-invFp_temp12*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp2112=-invFp_temp22*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp2113=-invFp_temp32*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp2121=-invFp_temp12*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp2122=-invFp_temp22*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp2123=-invFp_temp32*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp2131=-invFp_temp12*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp2132=-invFp_temp22*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp2133=-invFp_temp32*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp2211=-invFp_temp12*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp2212=-invFp_temp22*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp2213=-invFp_temp32*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp2221=-invFp_temp12*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp2222=-invFp_temp22*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp2223=-invFp_temp32*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp2231=-invFp_temp12*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp2232=-invFp_temp22*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp2233=-invFp_temp32*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp2311=-invFp_temp12*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp2312=-invFp_temp22*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp2313=-invFp_temp32*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp2321=-invFp_temp12*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp2322=-invFp_temp22*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp2323=-invFp_temp32*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp2331=-invFp_temp12*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp2332=-invFp_temp22*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp2333=-invFp_temp32*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp3111=-invFp_temp13*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp3112=-invFp_temp23*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp3113=-invFp_temp33*(F11*invFp_temp11 + F12*invFp_temp21 + F13*invFp_temp31);
                dFeCdFp_temp3121=-invFp_temp13*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp3122=-invFp_temp23*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp3123=-invFp_temp33*(F11*invFp_temp12 + F12*invFp_temp22 + F13*invFp_temp32);
                dFeCdFp_temp3131=-invFp_temp13*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp3132=-invFp_temp23*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp3133=-invFp_temp33*(F11*invFp_temp13 + F12*invFp_temp23 + F13*invFp_temp33);
                dFeCdFp_temp3211=-invFp_temp13*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp3212=-invFp_temp23*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp3213=-invFp_temp33*(F21*invFp_temp11 + F22*invFp_temp21 + F23*invFp_temp31);
                dFeCdFp_temp3221=-invFp_temp13*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp3222=-invFp_temp23*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp3223=-invFp_temp33*(F21*invFp_temp12 + F22*invFp_temp22 + F23*invFp_temp32);
                dFeCdFp_temp3231=-invFp_temp13*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp3232=-invFp_temp23*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp3233=-invFp_temp33*(F21*invFp_temp13 + F22*invFp_temp23 + F23*invFp_temp33);
                dFeCdFp_temp3311=-invFp_temp13*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp3312=-invFp_temp23*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp3313=-invFp_temp33*(F31*invFp_temp11 + F32*invFp_temp21 + F33*invFp_temp31);
                dFeCdFp_temp3321=-invFp_temp13*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp3322=-invFp_temp23*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp3323=-invFp_temp33*(F31*invFp_temp12 + F32*invFp_temp22 + F33*invFp_temp32);
                dFeCdFp_temp3331=-invFp_temp13*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp3332=-invFp_temp23*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);
                dFeCdFp_temp3333=-invFp_temp33*(F31*invFp_temp13 + F32*invFp_temp23 + F33*invFp_temp33);

                dFp_tempCddelta_epcum_last=Fp'*((eigenvalues(1,1)*exp(depcum_last*eigenvalues(1,1))*eigenvectors(:,1)*eigenvectors(:,1)'+eigenvalues(2,2)*exp(depcum_last*eigenvalues(2,2))*eigenvectors(:,2)*eigenvectors(:,2)'+eigenvalues(3,3)*exp(depcum_last*eigenvalues(3,3))*eigenvectors(:,3)*eigenvectors(:,3)'));

                dyieldfunctiondFp_temp11=dFeCdFp_temp1111*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1111*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1211*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1111*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1211*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1311*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1211*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1311*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1311*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1111*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2111*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2111*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1111*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1211*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2211*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2111*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1111*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2211*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1211*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1311*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2311*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2211*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1211*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2311*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1311*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2311*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1311*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1111*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2111*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3111*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3111*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1111*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1211*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2111*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2211*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3211*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3111*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2111*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3211*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1111*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1211*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1311*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2211*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2311*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3311*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3211*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2211*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3311*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1211*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1311*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2311*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3311*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2311*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1311*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2111*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3111*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2111*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2211*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3111*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3211*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3111*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2111*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2211*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2311*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3211*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3311*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3211*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2211*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2311*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3311*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3311*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2311*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3111*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3111*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3211*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3111*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3211*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3311*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3211*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3311*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3311*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp12=dFeCdFp_temp1112*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1112*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1212*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1112*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1212*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1312*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1212*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1312*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1312*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1112*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2112*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2112*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1112*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1212*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2212*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2112*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1112*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2212*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1212*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1312*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2312*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2212*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1212*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2312*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1312*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2312*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1312*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1112*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2112*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3112*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3112*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1112*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1212*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2112*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2212*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3212*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3112*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2112*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3212*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1112*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1212*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1312*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2212*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2312*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3312*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3212*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2212*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3312*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1212*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1312*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2312*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3312*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2312*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1312*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2112*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3112*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2112*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2212*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3112*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3212*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3112*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2112*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2212*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2312*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3212*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3312*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3212*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2212*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2312*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3312*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3312*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2312*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3112*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3112*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3212*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3112*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3212*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3312*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3212*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3312*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3312*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp13=dFeCdFp_temp1113*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1113*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1213*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1113*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1213*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1313*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1213*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1313*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1313*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1113*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2113*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2113*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1113*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1213*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2213*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2113*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1113*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2213*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1213*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1313*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2313*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2213*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1213*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2313*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1313*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2313*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1313*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1113*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2113*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3113*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3113*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1113*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1213*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2113*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2213*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3213*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3113*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2113*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3213*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1113*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1213*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1313*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2213*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2313*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3313*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3213*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2213*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3313*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1213*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1313*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2313*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3313*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2313*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1313*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2113*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3113*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2113*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2213*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3113*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3213*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3113*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2113*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2213*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2313*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3213*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3313*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3213*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2213*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2313*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3313*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3313*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2313*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3113*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3113*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3213*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3113*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3213*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3313*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3213*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3313*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3313*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp21=dFeCdFp_temp1121*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1121*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1221*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1121*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1221*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1321*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1221*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1321*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1321*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1121*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2121*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2121*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1121*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1221*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2221*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2121*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1121*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2221*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1221*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1321*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2321*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2221*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1221*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2321*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1321*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2321*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1321*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1121*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2121*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3121*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3121*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1121*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1221*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2121*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2221*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3221*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3121*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2121*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3221*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1121*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1221*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1321*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2221*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2321*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3321*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3221*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2221*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3321*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1221*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1321*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2321*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3321*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2321*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1321*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2121*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3121*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2121*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2221*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3121*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3221*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3121*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2121*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2221*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2321*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3221*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3321*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3221*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2221*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2321*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3321*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3321*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2321*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3121*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3121*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3221*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3121*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3221*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3321*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3221*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3321*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3321*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp22=dFeCdFp_temp1122*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1122*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1222*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1122*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1222*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1322*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1222*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1322*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1322*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1122*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2122*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2122*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1122*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1222*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2222*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2122*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1122*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2222*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1222*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1322*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2322*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2222*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1222*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2322*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1322*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2322*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1322*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1122*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2122*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3122*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3122*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1122*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1222*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2122*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2222*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3222*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3122*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2122*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3222*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1122*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1222*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1322*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2222*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2322*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3322*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3222*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2222*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3322*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1222*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1322*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2322*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3322*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2322*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1322*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2122*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3122*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2122*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2222*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3122*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3222*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3122*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2122*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2222*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2322*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3222*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3322*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3222*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2222*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2322*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3322*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3322*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2322*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3122*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3122*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3222*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3122*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3222*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3322*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3222*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3322*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3322*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp23=dFeCdFp_temp1123*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1123*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1223*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1123*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1223*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1323*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1223*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1323*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1323*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1123*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2123*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2123*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1123*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1223*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2223*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2123*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1123*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2223*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1223*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1323*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2323*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2223*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1223*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2323*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1323*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2323*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1323*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1123*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2123*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3123*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3123*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1123*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1223*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2123*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2223*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3223*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3123*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2123*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3223*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1123*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1223*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1323*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2223*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2323*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3323*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3223*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2223*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3323*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1223*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1323*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2323*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3323*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2323*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1323*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2123*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3123*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2123*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2223*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3123*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3223*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3123*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2123*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2223*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2323*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3223*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3323*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3223*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2223*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2323*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3323*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3323*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2323*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3123*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3123*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3223*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3123*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3223*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3323*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3223*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3323*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3323*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp31=dFeCdFp_temp1131*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1131*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1231*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1131*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1231*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1331*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1231*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1331*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1331*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1131*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2131*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2131*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1131*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1231*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2231*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2131*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1131*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2231*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1231*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1331*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2331*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2231*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1231*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2331*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1331*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2331*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1331*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1131*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2131*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3131*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3131*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1131*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1231*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2131*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2231*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3231*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3131*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2131*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3231*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1131*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1231*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1331*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2231*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2331*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3331*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3231*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2231*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3331*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1231*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1331*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2331*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3331*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2331*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1331*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2131*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3131*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2131*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2231*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3131*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3231*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3131*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2131*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2231*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2331*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3231*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3331*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3231*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2231*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2331*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3331*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3331*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2331*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3131*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3131*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3231*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3131*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3231*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3331*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3231*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3331*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3331*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp32=dFeCdFp_temp1132*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1132*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1232*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1132*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1232*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1332*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1232*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1332*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1332*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1132*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2132*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2132*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1132*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1232*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2232*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2132*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1132*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2232*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1232*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1332*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2332*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2232*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1232*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2332*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1332*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2332*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1332*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1132*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2132*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3132*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3132*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1132*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1232*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2132*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2232*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3232*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3132*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2132*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3232*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1132*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1232*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1332*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2232*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2332*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3332*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3232*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2232*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3332*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1232*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1332*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2332*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3332*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2332*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1332*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2132*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3132*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2132*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2232*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3132*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3232*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3132*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2132*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2232*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2332*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3232*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3332*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3232*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2232*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2332*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3332*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3332*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2332*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3132*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3132*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3232*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3132*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3232*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3332*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3232*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3332*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3332*dMCdFe3333*dyieldfunctiondM33;
                dyieldfunctiondFp_temp33=dFeCdFp_temp1133*dMCdFe1111*dyieldfunctiondM11 + dFeCdFp_temp1133*dMCdFe1211*dyieldfunctiondM12 + dFeCdFp_temp1233*dMCdFe1121*dyieldfunctiondM11 + dFeCdFp_temp1133*dMCdFe1311*dyieldfunctiondM13 + dFeCdFp_temp1233*dMCdFe1221*dyieldfunctiondM12 + dFeCdFp_temp1333*dMCdFe1131*dyieldfunctiondM11 + dFeCdFp_temp1233*dMCdFe1321*dyieldfunctiondM13 + dFeCdFp_temp1333*dMCdFe1231*dyieldfunctiondM12 + dFeCdFp_temp1333*dMCdFe1331*dyieldfunctiondM13 + dFeCdFp_temp1133*dMCdFe2111*dyieldfunctiondM12 + dFeCdFp_temp2133*dMCdFe1112*dyieldfunctiondM11 + dFeCdFp_temp2133*dMCdFe1212*dyieldfunctiondM12 + dFeCdFp_temp1133*dMCdFe2211*dyieldfunctiondM22 + dFeCdFp_temp1233*dMCdFe2121*dyieldfunctiondM12 + dFeCdFp_temp2233*dMCdFe1122*dyieldfunctiondM11 + dFeCdFp_temp2133*dMCdFe1312*dyieldfunctiondM13 + dFeCdFp_temp1133*dMCdFe2311*dyieldfunctiondM23 + dFeCdFp_temp2233*dMCdFe1222*dyieldfunctiondM12 + dFeCdFp_temp1233*dMCdFe2221*dyieldfunctiondM22 + dFeCdFp_temp1333*dMCdFe2131*dyieldfunctiondM12 + dFeCdFp_temp2333*dMCdFe1132*dyieldfunctiondM11 + dFeCdFp_temp2233*dMCdFe1322*dyieldfunctiondM13 + dFeCdFp_temp1233*dMCdFe2321*dyieldfunctiondM23 + dFeCdFp_temp2333*dMCdFe1232*dyieldfunctiondM12 + dFeCdFp_temp1333*dMCdFe2231*dyieldfunctiondM22 + dFeCdFp_temp2333*dMCdFe1332*dyieldfunctiondM13 + dFeCdFp_temp1333*dMCdFe2331*dyieldfunctiondM23 + dFeCdFp_temp1133*dMCdFe3111*dyieldfunctiondM13 + dFeCdFp_temp2133*dMCdFe2112*dyieldfunctiondM12 + dFeCdFp_temp3133*dMCdFe1113*dyieldfunctiondM11 + dFeCdFp_temp3133*dMCdFe1213*dyieldfunctiondM12 + dFeCdFp_temp1133*dMCdFe3211*dyieldfunctiondM23 + dFeCdFp_temp1233*dMCdFe3121*dyieldfunctiondM13 + dFeCdFp_temp2133*dMCdFe2212*dyieldfunctiondM22 + dFeCdFp_temp2233*dMCdFe2122*dyieldfunctiondM12 + dFeCdFp_temp3233*dMCdFe1123*dyieldfunctiondM11 + dFeCdFp_temp3133*dMCdFe1313*dyieldfunctiondM13 + dFeCdFp_temp2133*dMCdFe2312*dyieldfunctiondM23 + dFeCdFp_temp3233*dMCdFe1223*dyieldfunctiondM12 + dFeCdFp_temp1133*dMCdFe3311*dyieldfunctiondM33 + dFeCdFp_temp1233*dMCdFe3221*dyieldfunctiondM23 + dFeCdFp_temp1333*dMCdFe3131*dyieldfunctiondM13 + dFeCdFp_temp2233*dMCdFe2222*dyieldfunctiondM22 + dFeCdFp_temp2333*dMCdFe2132*dyieldfunctiondM12 + dFeCdFp_temp3333*dMCdFe1133*dyieldfunctiondM11 + dFeCdFp_temp3233*dMCdFe1323*dyieldfunctiondM13 + dFeCdFp_temp2233*dMCdFe2322*dyieldfunctiondM23 + dFeCdFp_temp3333*dMCdFe1233*dyieldfunctiondM12 + dFeCdFp_temp1233*dMCdFe3321*dyieldfunctiondM33 + dFeCdFp_temp1333*dMCdFe3231*dyieldfunctiondM23 + dFeCdFp_temp2333*dMCdFe2232*dyieldfunctiondM22 + dFeCdFp_temp3333*dMCdFe1333*dyieldfunctiondM13 + dFeCdFp_temp2333*dMCdFe2332*dyieldfunctiondM23 + dFeCdFp_temp1333*dMCdFe3331*dyieldfunctiondM33 + dFeCdFp_temp2133*dMCdFe3112*dyieldfunctiondM13 + dFeCdFp_temp3133*dMCdFe2113*dyieldfunctiondM12 + dFeCdFp_temp2133*dMCdFe3212*dyieldfunctiondM23 + dFeCdFp_temp2233*dMCdFe3122*dyieldfunctiondM13 + dFeCdFp_temp3133*dMCdFe2213*dyieldfunctiondM22 + dFeCdFp_temp3233*dMCdFe2123*dyieldfunctiondM12 + dFeCdFp_temp3133*dMCdFe2313*dyieldfunctiondM23 + dFeCdFp_temp2133*dMCdFe3312*dyieldfunctiondM33 + dFeCdFp_temp2233*dMCdFe3222*dyieldfunctiondM23 + dFeCdFp_temp2333*dMCdFe3132*dyieldfunctiondM13 + dFeCdFp_temp3233*dMCdFe2223*dyieldfunctiondM22 + dFeCdFp_temp3333*dMCdFe2133*dyieldfunctiondM12 + dFeCdFp_temp3233*dMCdFe2323*dyieldfunctiondM23 + dFeCdFp_temp2233*dMCdFe3322*dyieldfunctiondM33 + dFeCdFp_temp2333*dMCdFe3232*dyieldfunctiondM23 + dFeCdFp_temp3333*dMCdFe2233*dyieldfunctiondM22 + dFeCdFp_temp3333*dMCdFe2333*dyieldfunctiondM23 + dFeCdFp_temp2333*dMCdFe3332*dyieldfunctiondM33 + dFeCdFp_temp3133*dMCdFe3113*dyieldfunctiondM13 + dFeCdFp_temp3133*dMCdFe3213*dyieldfunctiondM23 + dFeCdFp_temp3233*dMCdFe3123*dyieldfunctiondM13 + dFeCdFp_temp3133*dMCdFe3313*dyieldfunctiondM33 + dFeCdFp_temp3233*dMCdFe3223*dyieldfunctiondM23 + dFeCdFp_temp3333*dMCdFe3133*dyieldfunctiondM13 + dFeCdFp_temp3233*dMCdFe3323*dyieldfunctiondM33 + dFeCdFp_temp3333*dMCdFe3233*dyieldfunctiondM23 + dFeCdFp_temp3333*dMCdFe3333*dyieldfunctiondM33;

                dyieldfunctionddepcum_last=-Hardening_modulus*Hardening_exponent*(epcum+depcum_last)^(Hardening_exponent-1);
                dyieldfunctionddepcum_last=dyieldfunctionddepcum_last+sum(sum([dyieldfunctiondFp_temp11 dyieldfunctiondFp_temp12 dyieldfunctiondFp_temp13; dyieldfunctiondFp_temp21 dyieldfunctiondFp_temp22 dyieldfunctiondFp_temp23; dyieldfunctiondFp_temp31 dyieldfunctiondFp_temp32 dyieldfunctiondFp_temp33].*dFp_tempCddelta_epcum_last'));





                delta_depcum_last=-yieldfunction/dyieldfunctionddepcum_last;

                depcum_last=depcum_last+delta_depcum_last;




                dyieldfunctiondM=1.5*Mdev/sqrt(1.5*(Mdev(1,1)^2+Mdev(2,2)^2+Mdev(3,3)^2+2*Mdev(1,2)^2+2*Mdev(1,3)^2+2*Mdev(2,3)^2));

                [eigenvectors,eigenvalues]=eig(dyieldfunctiondM);
                Fp_temp=(exp(depcum_last*eigenvalues(1,1))*eigenvectors(:,1)*eigenvectors(:,1)'+exp(depcum_last*eigenvalues(2,2))*eigenvectors(:,2)*eigenvectors(:,2)'+exp(depcum_last*eigenvalues(3,3))*eigenvectors(:,3)*eigenvectors(:,3)')*Fp;

                Fe=F*inv(Fp_temp);

                P=Ca*(Fe-inv(Fe)')+Cb*log(det(Fe))*inv(Fe)';

                Mdev=Fe'*P-trace(Fe'*P)/3*eye(3);

                yieldfunction=sqrt(1.5*(Mdev(1,1)^2+Mdev(2,2)^2+Mdev(3,3)^2+2*Mdev(1,2)^2+2*Mdev(1,3)^2+2*Mdev(2,3)^2))-My0-Hardening_modulus*(epcum+depcum_last)^Hardening_exponent;

            end

        end

        DELTA_EPCUM(jj,ii)=depcum_last;
        FP_temp(:,(ii-1)*8+jj)=reshape(Fp_temp',[],1);

        dEdu=zeros(24,1);
%         dEdu(1:3:24-2)=dNdXYZ(1:3:24-2)*P(1,1)+dNdXYZ(2:3:24-1)*P(1,2)+dNdXYZ(3:3:24)*P(1,3);
%         dEdu(2:3:24-1)=dNdXYZ(1:3:24-2)*P(2,1)+dNdXYZ(2:3:24-1)*P(2,2)+dNdXYZ(3:3:24)*P(2,3);
%         dEdu(3:3:24)=dNdXYZ(1:3:24-2)*P(3,1)+dNdXYZ(2:3:24-1)*P(3,2)+dNdXYZ(3:3:24)*P(3,3);
        A=inv(Fp_temp)*P';
        dEdu(1:3:24-2)=dNdXYZ(1:3:24-2)*A(1,1)+dNdXYZ(2:3:24-1)*A(2,1)+dNdXYZ(3:3:24)*A(3,1);
        dEdu(2:3:24-1)=dNdXYZ(1:3:24-2)*A(1,2)+dNdXYZ(2:3:24-1)*A(2,2)+dNdXYZ(3:3:24)*A(3,2);
        dEdu(3:3:24)=dNdXYZ(1:3:24-2)*A(1,3)+dNdXYZ(2:3:24-1)*A(2,3)+dNdXYZ(3:3:24)*A(3,3);

        dEdU=dEdU+detJ*dEdu;

        mx=2*Ca*(trace(Fe'*Fe)-3-2*log(det(Fe)))+2*Cb*log(det(Fe))^2;
        mx=mx+Hardening_modulus/(Hardening_exponent+1)*(epcum+depcum_last)^(Hardening_exponent+1);
        mx=mx+My0*depcum_last;

        m=m+detJ*mx;

    end

    f_in(a_mat,1)=f_in(a_mat,1)+dEdU;

end

f=f_in(indices_neumann);

end







function [f_in] = f_only(LEN_node,u,LEN_el,FE_node_coord,con_mat,Young,nu,My0,Hardening_modulus,Hardening_exponent,FP_conv,EPCUM)

KSI=1/sqrt(3)*[-1 -1 -1;
    +1 -1 -1;
    +1 +1 -1;
    -1 +1 -1;
    -1 -1 +1;
    +1 -1 +1;
    +1 +1 +1;
    -1 +1 +1];

Ca=Young/(2+2*nu);
Cb=Young*nu/((1+nu)*(1-2*nu));

f_in=zeros(LEN_node*3,1);
for ii=1:length(con_mat(:,1))

    nr1=con_mat(ii,1);
    nr2=con_mat(ii,2);
    nr3=con_mat(ii,3);
    nr4=con_mat(ii,4);
    nr5=con_mat(ii,5);
    nr6=con_mat(ii,6);
    nr7=con_mat(ii,7);
    nr8=con_mat(ii,8);

    x1=FE_node_coord(nr1,1);
    y1=FE_node_coord(nr1,2);
    z1=FE_node_coord(nr1,3);

    x2=FE_node_coord(nr2,1);
    y2=FE_node_coord(nr2,2);
    z2=FE_node_coord(nr2,3);

    x3=FE_node_coord(nr3,1);
    y3=FE_node_coord(nr3,2);
    z3=FE_node_coord(nr3,3);

    x4=FE_node_coord(nr4,1);
    y4=FE_node_coord(nr4,2);
    z4=FE_node_coord(nr4,3);

    x5=FE_node_coord(nr5,1);
    y5=FE_node_coord(nr5,2);
    z5=FE_node_coord(nr5,3);

    x6=FE_node_coord(nr6,1);
    y6=FE_node_coord(nr6,2);
    z6=FE_node_coord(nr6,3);

    x7=FE_node_coord(nr7,1);
    y7=FE_node_coord(nr7,2);
    z7=FE_node_coord(nr7,3);

    x8=FE_node_coord(nr8,1);
    y8=FE_node_coord(nr8,2);
    z8=FE_node_coord(nr8,3);


    a_mat=[nr1*3-2;nr1*3-1;nr1*3;nr2*3-2;nr2*3-1;nr2*3;nr3*3-2;nr3*3-1;nr3*3;nr4*3-2;nr4*3-1;nr4*3;nr5*3-2;nr5*3-1;nr5*3;nr6*3-2;nr6*3-1;nr6*3;nr7*3-2;nr7*3-1;nr7*3;nr8*3-2;nr8*3-1;nr8*3];

    xyz_all=[x1 y1 z1;x2 y2 z2;x3 y3 z3;x4 y4 z4;x5 y5 z5;x6 y6 z6;x7 y7 z7;x8 y8 z8];


    u_all=u(a_mat(1:3:24-2));
    v_all=u(a_mat(2:3:24-1));
    w_all=u(a_mat(3:3:24));

    dEdU=zeros(24,1);

    for jj=1:8

        ksi1=KSI(jj,1);
        ksi2=KSI(jj,2);
        ksi3=KSI(jj,3);

        dN1dksi1=-((ksi2 - 1)*(ksi3 - 1))/8;
        dN1dksi2=-((ksi1 - 1)*(ksi3 - 1))/8;
        dN1dksi3=-((ksi1 - 1)*(ksi2 - 1))/8;
        dN2dksi1=((ksi2 - 1)*(ksi3 - 1))/8;
        dN2dksi2=((ksi1 + 1)*(ksi3 - 1))/8;
        dN2dksi3=((ksi1 + 1)*(ksi2 - 1))/8;
        dN3dksi1=-((ksi2 + 1)*(ksi3 - 1))/8;
        dN3dksi2=-((ksi1 + 1)*(ksi3 - 1))/8;
        dN3dksi3=-((ksi1 + 1)*(ksi2 + 1))/8;
        dN4dksi1=((ksi2 + 1)*(ksi3 - 1))/8;
        dN4dksi2=((ksi1 - 1)*(ksi3 - 1))/8;
        dN4dksi3=((ksi1 - 1)*(ksi2 + 1))/8;
        dN5dksi1=((ksi2 - 1)*(ksi3 + 1))/8;
        dN5dksi2=((ksi1 - 1)*(ksi3 + 1))/8;
        dN5dksi3=((ksi1 - 1)*(ksi2 - 1))/8;
        dN6dksi1=-((ksi2 - 1)*(ksi3 + 1))/8;
        dN6dksi2=-((ksi1 + 1)*(ksi3 + 1))/8;
        dN6dksi3=-((ksi1 + 1)*(ksi2 - 1))/8;
        dN7dksi1=((ksi2 + 1)*(ksi3 + 1))/8;
        dN7dksi2=((ksi1 + 1)*(ksi3 + 1))/8;
        dN7dksi3=((ksi1 + 1)*(ksi2 + 1))/8;
        dN8dksi1=-((ksi2 + 1)*(ksi3 + 1))/8;
        dN8dksi2=-((ksi1 - 1)*(ksi3 + 1))/8;
        dN8dksi3=-((ksi1 - 1)*(ksi2 + 1))/8;

        dNdksi=[dN1dksi1 dN2dksi1 dN3dksi1 dN4dksi1 dN5dksi1 dN6dksi1 dN7dksi1 dN8dksi1;
            dN1dksi2 dN2dksi2 dN3dksi2 dN4dksi2 dN5dksi2 dN6dksi2 dN7dksi2 dN8dksi2;
            dN1dksi3 dN2dksi3 dN3dksi3 dN4dksi3 dN5dksi3 dN6dksi3 dN7dksi3 dN8dksi3];

        Jmat=dNdksi*xyz_all;

        detJ=det(Jmat);
        invJmat=inv(Jmat);

        dNdXYZ=[invJmat*dNdksi(:,1);
            invJmat*dNdksi(:,2);
            invJmat*dNdksi(:,3);
            invJmat*dNdksi(:,4);
            invJmat*dNdksi(:,5);
            invJmat*dNdksi(:,6);
            invJmat*dNdksi(:,7);
            invJmat*dNdksi(:,8)];

        F=[dNdXYZ(1:3:24-2)'*u_all dNdXYZ(2:3:24-1)'*u_all dNdXYZ(3:3:24)'*u_all;
            dNdXYZ(1:3:24-2)'*v_all dNdXYZ(2:3:24-1)'*v_all dNdXYZ(3:3:24)'*v_all;
            dNdXYZ(1:3:24-2)'*w_all dNdXYZ(2:3:24-1)'*w_all dNdXYZ(3:3:24)'*w_all]+eye(3);

        Fp=reshape(FP_conv(:,(ii-1)*8+jj),3,3)';

        epcum=EPCUM(jj,ii);

        Fe=F*inv(Fp);

        P=Ca*(Fe-inv(Fe)')+Cb*log(det(Fe))*inv(Fe)';
        

        dEdu=zeros(24,1);
%         dEdu(1:3:24-2)=dNdXYZ(1:3:24-2)*P(1,1)+dNdXYZ(2:3:24-1)*P(1,2)+dNdXYZ(3:3:24)*P(1,3);
%         dEdu(2:3:24-1)=dNdXYZ(1:3:24-2)*P(2,1)+dNdXYZ(2:3:24-1)*P(2,2)+dNdXYZ(3:3:24)*P(2,3);
%         dEdu(3:3:24)=dNdXYZ(1:3:24-2)*P(3,1)+dNdXYZ(2:3:24-1)*P(3,2)+dNdXYZ(3:3:24)*P(3,3);
        A=inv(Fp_temp)*P';
        dEdu(1:3:24-2)=dNdXYZ(1:3:24-2)*A(1,1)+dNdXYZ(2:3:24-1)*A(2,1)+dNdXYZ(3:3:24)*A(3,1);
        dEdu(2:3:24-1)=dNdXYZ(1:3:24-2)*A(1,2)+dNdXYZ(2:3:24-1)*A(2,2)+dNdXYZ(3:3:24)*A(3,2);
        dEdu(3:3:24)=dNdXYZ(1:3:24-2)*A(1,3)+dNdXYZ(2:3:24-1)*A(2,3)+dNdXYZ(3:3:24)*A(3,3);

        dEdU=dEdU+detJ*dEdu;

    end

    f_in(a_mat,1)=f_in(a_mat,1)+dEdU;

end

end