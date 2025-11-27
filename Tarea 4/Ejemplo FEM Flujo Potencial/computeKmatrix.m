
function [ Ke, Ke_CST, DET_JAC ]  = computeKmatrix( TRIE )
% Compute element stiffness matrix - one integration point
global NODES 

n1 = TRIE(1); n2 = TRIE(2); n3 = TRIE(3);  

x1 = NODES(n1,1); y1 = NODES(n1,2);
x2 = NODES(n2,1); y2 = NODES(n2,2);
x3 = NODES(n3,1); y3 = NODES(n3,2);

AREA2 = abs(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);

AREA = AREA2/2; %Nota: el área siempre es una cantidad positiva. 

DET_JAC = AREA2;

JAC_INV2 = (1/DET_JAC)*[ [(y3-y1), (x1-x3)];[(y1-y2), (x2-x1)] ];

JAC =  [ [x2 - x1, x3 - x1]; [y2 - y1, y3 - y1] ];

JAC_INV = ...
[ [-(y1 - y3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2),  (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)];...
[ (y1 - y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2), -(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)] ];
 

NGRAD = [ [-1 -1]; [1 0]; [0 1] ];

Ke = zeros(3,3);

BETA  = [ (y2-y3), (y3-y1), (y1-y2)];
GAMMA = [ -(x2-x3), -(x3-x1), -(x1-x2)];

for i=1:1:3
for j=1:1:3

    Ke_CST(i,j) = 1/(4*AREA)*( BETA(i)*BETA(j) + GAMMA(i)*GAMMA(j) );

    Ke(i,j) = (1/2)*(DET_JAC)*dot( JAC_INV'*NGRAD(i,1:2)', JAC_INV'*NGRAD(j,1:2)' );

end 
end



