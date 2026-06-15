% ------------------------------------------------------------
% Résolution par éléments finis d'un problème non linéaire
% Équation : -div( (1 + a(x,y)*|∇u|^2) ∇u ) - sin(x) u = f  dans Ω
% avec f = 1 (terme source constant) → solution non triviale
% Condition Dirichlet : u = 0 sur ∂Ω
% Ω = ]0,1[ x ]0,1[, a(x,y) = 1 + x^2 + y^2
% ------------------------------------------------------------
clear; clc; close all;

%% Paramètres du maillage et de la méthode
Nx = 50;                     % nombre de subdivisions par côté
tol = 1e-8;                  % tolérance pour Newton
max_iter = 25;
source = 1.0;                % amplitude du terme source constant
alpha0 = 0.5;                % amplitude de l'initialisation

%% 1. Génération du maillage (triangles P1) - intégrée ici
% Grille structurée de (Nx+1)x(Nx+1) points
x_lin = linspace(0, 1, Nx+1);
y_lin = linspace(0, 1, Nx+1);
[X, Y] = meshgrid(x_lin, y_lin);
x = X(:);
y = Y(:);
nnodes = length(x);

% Numérotation des triangles (deux par carré)
Tri = zeros(2*Nx^2, 3);
idx = @(i,j) (j-1)*(Nx+1) + i;
k = 0;
for j = 1:Nx
    for i = 1:Nx
        n1 = idx(i,j);
        n2 = idx(i+1,j);
        n3 = idx(i,j+1);
        n4 = idx(i+1,j+1);
        k = k+1; Tri(k,:) = [n1, n2, n3];
        k = k+1; Tri(k,:) = [n4, n2, n3];
    end
end
nelems = size(Tri, 1);

% Noeuds du bord
boundary_nodes = find(x==0 | x==1 | y==0 | y==1);
inner_nodes = setdiff(1:nnodes, boundary_nodes);

%% 2. Précalculs géométriques (gradients des fonctions de forme, aires, points de Gauss)
% Points de Gauss sur le triangle de référence (ξ,η) et poids
gauss_xi  = [1/6, 2/3, 1/6];
gauss_eta = [1/6, 1/6, 2/3];
gauss_w   = [1/6, 1/6, 1/6];
nqp = length(gauss_w);

dphi_dx = zeros(3, nelems);
dphi_dy = zeros(3, nelems);
detJ_vec = zeros(nelems, 1);  % déterminant du Jacobien pour chaque élément

for e = 1:nelems
    nodes = Tri(e,:);
    xe = x(nodes);
    ye = y(nodes);
    % Matrice Jacobienne (transformation élément de référence -> réel)
    J = [xe(2)-xe(1), xe(3)-xe(1);
         ye(2)-ye(1), ye(3)-ye(1)];
    detJ = abs(det(J));
    detJ_vec(e) = detJ;
    invJ = inv(J);
    % Dérivées des fonctions de forme dans l'élément réel
    dphi_dxi  = [-1, 1, 0];
    dphi_deta = [-1, 0, 1];
    dphi_dx(:,e) = invJ(1,1)*dphi_dxi + invJ(1,2)*dphi_deta;
    dphi_dy(:,e) = invJ(2,1)*dphi_dxi + invJ(2,2)*dphi_deta;
end

%% 3. Solution initiale (non triviale) : sin(πx) sin(πy)
u = alpha0 * sin(pi*x) .* sin(pi*y);
u(boundary_nodes) = 0;

%% 4. Boucle de Newton
iter = 0;
res_norm = inf;
fprintf('Itération  ||Résidu||\n');
while res_norm > tol && iter < max_iter
    iter = iter + 1;
    R = sparse(nnodes, 1);   % résidu
    J = sparse(nnodes, nnodes); % matrice jacobienne
    
    for e = 1:nelems
        nodes = Tri(e,:);
        ue = u(nodes);
        xe = x(nodes);
        ye = y(nodes);
        dpx = dphi_dx(:,e);
        dpy = dphi_dy(:,e);
        detJ = detJ_vec(e);
        
        % Gradient de u sur l'élément (constant car P1)
        gx = dpx' * ue;
        gy = dpy' * ue;
        norm_g2 = gx^2 + gy^2;
        
        Rlocal = zeros(3,1);
        Jlocal = zeros(3,3);
        
        for q = 1:nqp
            xi  = gauss_xi(q);
            eta = gauss_eta(q);
            wq  = gauss_w(q);
            % Coordonnées réelles du point d'intégration
            xq = xe(1) + (xe(2)-xe(1))*xi + (xe(3)-xe(1))*eta;
            yq = ye(1) + (ye(2)-ye(1))*xi + (ye(3)-ye(1))*eta;
            phi = [1-xi-eta, xi, eta];  % fonctions de forme
            weight = wq * detJ;
            a_q = 1 + xq^2 + yq^2;
            uh = phi * ue;
            coeff = 1 + a_q * norm_g2;
            
            % Contribution au résidu
            for i = 1:3
                Rlocal(i) = Rlocal(i) + weight * ( coeff * (gx*dpx(i)+gy*dpy(i)) ...
                                                 - sin(xq)*uh*phi(i) ...
                                                 - source * phi(i) );   % terme source constant
            end
            
            % Contribution au jacobien
            for i = 1:3
                for j = 1:3
                    term1 = coeff * (dpx(j)*dpx(i) + dpy(j)*dpy(i));
                    term2 = 2*a_q * (gx*dpx(j)+gy*dpy(j)) * (gx*dpx(i)+gy*dpy(i));
                    term3 = - sin(xq) * phi(j) * phi(i);
                    Jlocal(i,j) = Jlocal(i,j) + weight * (term1 + term2 + term3);
                end
            end
        end
        % Assemblage global
        R(nodes) = R(nodes) + Rlocal;
        J(nodes, nodes) = J(nodes, nodes) + Jlocal;
    end
    
    % Application des conditions de Dirichlet (u=0 sur le bord)
    R(boundary_nodes) = 0;
    J(boundary_nodes, :) = 0;
    J(boundary_nodes, boundary_nodes) = eye(length(boundary_nodes));
    
    % Résolution du système linéaire
    du = J \ (-R);
    u = u + du;
    u(boundary_nodes) = 0;
    
    res_norm = norm(R(inner_nodes));
    fprintf('%3d      %12.4e\n', iter, res_norm);
end

fprintf('Convergence atteinte en %d itérations.\n', iter);

%% 5. Post-traitement : visualisation de la solution
% Tracé 3D
figure;
trisurf(Tri, x, y, u, 'EdgeColor', 'none');
colormap(jet);
colorbar;
title(sprintf('Solution non triviale (terme source constant = %.2f)', source));
xlabel('x'); ylabel('y'); zlabel('u');
view(45,30);

% Courbes de niveau (interpolation sur grille fine pour contourf)
[Xg, Yg] = meshgrid(linspace(0,1,200), linspace(0,1,200));
Ug = griddata(x, y, u, Xg, Yg, 'linear');
figure;
contourf(Xg, Yg, Ug, 20);
colorbar;
title('Courbes de niveau de u');
axis equal; xlabel('x'); ylabel('y');

%% 6. Étude de convergence (optionnelle) : comparaison sur maillages raffinés
% Décommentez la section suivante pour une étude de convergence
% (nécessite d'exécuter la boucle Newton plusieurs fois)
%
% Nlist = [20, 40, 80];
% errL2 = zeros(size(Nlist));
% hlist = 1./Nlist;
% for idx = 1:length(Nlist)
%     N = Nlist(idx);
%     [x_c, y_c, Tri_c] = ... % recréer maillage
%     u_c = ... % résoudre sur ce maillage (répéter les étapes ci-dessus)
%     % Comparer avec une solution de référence sur N=160
%     errL2(idx) = ... 
% end
% figure; loglog(hlist, errL2, 'b-o'); xlabel('h'); ylabel('Erreur L^2');