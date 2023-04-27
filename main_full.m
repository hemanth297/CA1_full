clear all; close all

%% set up the model
nu = 0.4999; xint = 4; yint = 16; % modifying Poisson's ratio and mesh
Model = model_setup(nu,xint,yint);  % Set up your model in "model_setup.m" to be read here.

%% discretize
Mesh = sub_discretization ( Model );

%% detect boundary nodes
BC = sub_get_boundary ( Model, Mesh );
%% assembly
fig = figure 
hold on
for i = 1:4
if i == 3
    [ K , f ]  =  sub_assembly_SRI ( Model , Mesh , BC);
    [ d ] = sub_solution ( K , f , BC );
elseif i ==4
    u_exact  =   Model.exact.displ ( Mesh.x_node(:,1), Mesh.x_node(:,2) );
    d= u_exact';
    d=d(:); 
else
int_pts = i; %% 1 for RI, 2 for FI
[ K , f ]  =  sub_assembly_FIRI ( Model , Mesh , BC, int_pts );
[ d ] = sub_solution ( K , f , BC );
end
%% postprocess
[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );


%%sigmaxx
x0_ind = find(abs(x_eval(:,1) - max(x_eval(:,1))/2) < 1e-2);
plot(x_eval(x0_ind,2),stress(x0_ind,1),"LineWidth",2)
xlabel("y-axis",'FontSize',14)
ylabel("\sigma_{xx} (Pa)",'FontSize',14)
title(sprintf('sigma_{xx} at y=0, %d * %d mesh, nu = %0.4f ',xint,yint,nu),'FontSize',14)
% ylim([-0.3e8 0.3e8])

%%sigmaxy
% x0_ind = find(abs(x_eval(:,1) - max(x_eval(:,1))/2) < 1e-2);
% plot(hax1,x_eval(x0_ind,2),stress(x0_ind,3),"LineWidth",2)
% xlabel("y-axis",'FontSize',14)
% ylabel("\sigma_{xy}",'FontSize',14)
% title("\sigma_{xy} at x =L/2 ",'FontSize',14)

%%uy
% y0_ind = find(abs(x_eval(:,2)) < 1e-2);
% plot(x_eval(y0_ind,1),displ(y0_ind,2),"LineWidth",2)
% xlabel("x-axis",'FontSize',14)
% ylabel("u_{y}",'FontSize',14)
% title("u_{y} at y=0 ",'FontSize',14)

end
legend(["Reduced Integration", "Full Integration", "Selective RI",'Exact'],'FontSize',12)





