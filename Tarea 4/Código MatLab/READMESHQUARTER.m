clear all; clc; close all;
global NODES 

%Chord_200_radius_4000_pitch_5_md_40_8% 1177 nodos; 2178 triangulos; md = 40, mdr = 8
%Chord_200_radius_4000_pitch_5_md_30_6% 1720 nodos; 3232 triangulos; md = 30, mdr = 6
%Chord_200_radius_4000_pitch_5_md_20_4% 4795 nodos; 9224 triangulos; md = 20, mdr = 4
Chord_200_radius_4000_pitch_5_md_10_2% 1984 nodos; 3733 triangulos; md = 30, mdr = 5

NODES = msh.POS;
TRI = msh.TRIANGLES;
TRI = TRI(:,1:3);
trimesh(TRI, NODES(:,1), NODES(:,2), NODES(:,3),'FaceColor','cyan')
title('Potential flow - Laplace equation');
box on
grid on
axis equal
xlabel('x');
ylabel('y');
Chord = 200;

NUMNODES = size(NODES,1);
NUMTRI   = size(TRI,1);

for i=1:1:NUMNODES
    text(NODES(i,1),NODES(i,2), int2str(i));
end
view(2)

KG = zeros( NUMNODES, NUMNODES );

for i = 1:1:NUMTRI

[Ke, Ke_CST, DET_JAC ] = computeKmatrix( TRI(i,1:3) );

if ( DET_JAC > 0)

n1 = TRI(i,1); n2 = TRI(i,2); n3 = TRI(i,3); 

KG( n1, n1 ) = KG( n1, n1 ) + Ke(1,1);
KG( n1, n2 ) = KG( n1, n2 ) + Ke(1,2);
KG( n1, n3 ) = KG( n1, n3 ) + Ke(1,3);

KG( n2, n1 ) = KG( n2, n1 ) + Ke(2,1);
KG( n2, n2 ) = KG( n2, n2 ) + Ke(2,2);
KG( n2, n3 ) = KG( n2, n3 ) + Ke(2,3);

KG( n3, n1 ) = KG( n3, n1 ) + Ke(3,1);
KG( n3, n2 ) = KG( n3, n2 ) + Ke(3,2);
KG( n3, n3 ) = KG( n3, n3 ) + Ke(3,3);

else

disp('SINGULAR JACOBIAN')    

NUMELEMENT = i

NODESELEMENT = TRI(i,1:3)

DET_JAC

end

end

XMIN = min(NODES(:,1));
XMAX = max(NODES(:,1));
%%Set the value of the velocity potential
NODEVAL = find(NODES(:,1)==XMAX);
for i=1:1:length(NODEVAL)
KG(NODEVAL(i),:) = zeros(1,NUMNODES);
KG(NODEVAL(i),NODEVAL(i)) = 1.0;
end

figure(2)
spy(KG)

KGS = sparse(KG);

BOUNDARYLINES = msh.LINES;
BOUNDARYLINES = BOUNDARYLINES(:,1:2);

NUMLINES = size(BOUNDARYLINES,1);
FV = zeros(NUMNODES,1); 

YMIN = min(NODES(:,2));
YMAX = max(NODES(:,2));

U = 30;

k = 1;

for i=1:1:NUMLINES
 n1 = BOUNDARYLINES(i,1); n2 = BOUNDARYLINES(i,2);
 FV(n1) = 0;
 FV(n2) = 0;
end

for i=1:1:NUMLINES

 n1 = BOUNDARYLINES(i,1); n2 = BOUNDARYLINES(i,2);
 LB = norm( [ NODES(n2,1) - NODES(n1,1), NODES(n2,2) - NODES(n1,2)   ] ); 

 if( NODES(n1,1) == XMIN & NODES(n2,1) == XMIN )
     FV(n1) = FV(n1) + LB/2*U;
     FV(n2) = FV(n2) + LB/2*U;
 end
 if( NODES(n1,1) == XMAX & NODES(n2,1) == XMAX )
     FV(n1) = FV(n1) + LB/2*U;
     FV(n2) = FV(n2) + LB/2*U;
 end    

 if( NODES(n1,1) > XMIN & NODES(n1,1) < XMAX  )
 if( NODES(n1,2) > YMIN & NODES(n1,2) < YMAX  )

     INTERIORB(k,1:2) = [n1, n2];
     k = k+1;

 end
 end

end

for i=1:1:length(NODEVAL)
FV(NODEVAL(i)) = 0.0;
end

figure(3)
plot(NODES(INTERIORB,1), NODES(INTERIORB,2), 'ob')
box on
grid on
axis equal
xlabel('x');
ylabel('y');

PHI = KG\FV;

% TF = isnan(PHI);
% PHI(TF) = 0;

ERROR = norm( KG*PHI - FV  );
FINT = KG*PHI;
VDIF  = KG*PHI - FV; 

CMAP = colormap;
NMAX = size(CMAP,1);

figure(4)
hold on
for i = 1:1:NUMTRI
n1 = TRI(i,1); n2 = TRI(i,2); n3 = TRI(i,3); 
x = [NODES(n1,1); NODES(n2,1); NODES(n3,1)];
y = [NODES(n1,2); NODES(n2,2); NODES(n3,2)];

ic1 = round(1 + ( NMAX -1 )*( PHI(n1) - min(PHI))/( max(PHI) - min(PHI) ));
ic2 = round(1 + ( NMAX -1 )*( PHI(n2) - min(PHI))/( max(PHI) - min(PHI) ));
ic3 = round(1 + ( NMAX -1 )*( PHI(n3) - min(PHI))/( max(PHI) - min(PHI) ));

%c = [ ic1; ic2; ic3];
c = [ PHI(n1); PHI(n2); PHI(n3)];

patch(x,y,c)

end
caxis([min(PHI) max(PHI)]);
colorbar


title('Potential flow - Laplace equation - velocity potential');
box on
grid on
axis equal
xlabel('x');
ylabel('y');
view(2)
hold off

figure(5)
TO = triangulation(TRI,NODES(:,1), NODES(:,2), PHI);
trisurf(TO)
box on
grid on
axis equal
xlabel('x');
ylabel('y');
title('Potential flow - Laplace equation - velocity potential');
view(2)
caxis([min(PHI) max(PHI)]);
colorbar


GRADMATRIX = zeros(NUMTRI,2);

for i=1:1:NUMTRI
    
GRAD  = phigrad( TRI(i,1:3), PHI(TRI(i,1:3)) );
GRADMATRIX(i,1:2) = GRAD; 

CENTROIDMATRIX(i,1:2) = [mean( NODES(TRI(i,1:3),1) ), mean(NODES(TRI(i,1:3),2))];

PHICENTRIOD(i) = mean(PHI(TRI(i,1:3)));

VELOCITY(i) = norm(GRAD);

PRESS(i) = 1 - dot(GRAD,GRAD)/(U^2);

VEL_MAX = max(VELOCITY)
P_MAX = max(PRESS)

end

figure(6)
hold on
trimesh(TRI, NODES(:,1), NODES(:,2), NODES(:,3),'FaceColor','cyan')
quiver( CENTROIDMATRIX(:,1), CENTROIDMATRIX(:,2), -GRADMATRIX(:,1), -GRADMATRIX(:,2) )
title('Potential flow - Laplace equation - velocity field');
box on
grid on
axis equal 
xlabel('x');
ylabel('y');
hold off

figure(7)
hold on
for i = 1:1:NUMTRI
n1 = TRI(i,1); n2 = TRI(i,2); n3 = TRI(i,3); 
x = [NODES(n1,1); NODES(n2,1); NODES(n3,1)];
y = [NODES(n1,2); NODES(n2,2); NODES(n3,2)];

c = [ VELOCITY(i); VELOCITY(i); VELOCITY(i)];

patch(x,y,c)

end
axis equal 
xlabel('x');
ylabel('y');
title('Potential flow - Laplace equation - velocity norm');
caxis([min(VELOCITY) max(VELOCITY)]);
colorbar


figure(8)
hold on
for i = 1:1:NUMTRI
n1 = TRI(i,1); n2 = TRI(i,2); n3 = TRI(i,3); 
x = [NODES(n1,1); NODES(n2,1); NODES(n3,1)];
y = [NODES(n1,2); NODES(n2,2); NODES(n3,2)];

c = [ PRESS(i); PRESS(i); PRESS(i)];

patch(x,y,c)

end
axis equal 
xlabel('x');
ylabel('y');
title('Potential flow - Laplace equation - nondimensional pressure'); % (P - Pinf) / (0.5*Rho*(Vinf^2));
caxis([min(PRESS) max(PRESS)]);
colorbar

%------------------------------

%% Cálculo del coeficiente de Lift

Lift = compute_lift(NODES, TRI, PRESS, INTERIORB, Chord)

%% Gráfico de Cp sobre el perfil alar

profile_nodes = unique(INTERIORB(:));   % nodos del borde del perfil

elem_nodes = TRI;  % alias para legibilidad

is_profile_node = ismember(elem_nodes, profile_nodes);
count_profile   = sum(is_profile_node, 2);

mask = (count_profile >= 2);
tri_filtered = elem_nodes(mask, :);

Cp_elem_all = mean( PRESS(elem_nodes), 2 );
Cp_elem = Cp_elem_all(mask);

X = NODES(:,1);
Y = NODES(:,2);

figure('Position',[200 200 1200 450])

trisurf(tri_filtered, X, Y, zeros(size(X)), Cp_elem, ...
        'EdgeColor','k','LineWidth',0.3);

view(2)
shading flat
axis equal tight
colormap(turbo)
colorbar
caxis([min(Cp_elem) max(Cp_elem)])

hold on

edges = INTERIORB;
profile_curve = unique(edges(:));

plot(X(profile_curve), Y(profile_curve), 'k-', 'LineWidth', 1.2)

title('Distribución de Cp alrededor del perfil alar')
xlabel('x')
ylabel('y')






