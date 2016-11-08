function fem2D(scene)
% FEM explicit triangles

if nargin < 1
	scene = 0;
end

global video;
video = [];

switch(scene)
	case 0
        uniform = true; % Is triangle stiffness uniform throughout?
		nTiles = 1; % Number of tiles in x and y
        initTris(nTiles); % Create triangles from regular grid of nodes
        
        side = false; % Is a side pinned? true = top or bottom pinned.
        less = false; % Pin nodes above or below threshold.
        threshold = 0.9; % Threshold for which nodes to pin.
        pinNodes(side,less,threshold); % Pin appropriate nodes for scene.

		dt = 1e-3; % time step
		tEnd = 1.0; % end time
		drawHz = 100; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 2.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
	case 1
        uniform = true; % Is triangle stiffness uniform throughout?
%         nTiles = 1; % Number of tiles in x and y
%         nTiles = 2;
        nTiles = 4;
        %nTiles = 8;
        initTris(nTiles); % Create triangles from regular grid of nodes
        
        side = false; % Is a side pinned? true = top or bottom pinned.
        less = true;
        threshold = -0.9; % Threshold for which nodes to pin.
        pinNodes(side,less,threshold); % Pin appropriate nodes for scene.
        
		dt = 1e-3; % time step
		tEnd = 2.0; % end time
		drawHz = 100; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 5e0; % area density
		damping = 1.0; % viscous damping
		E = 1e3; % Young's modulus
		nu = 0.0; % Poisson's ratio
    case 2
        % 8x8 with top nodes fixed while reaching roughly y = -1.5
        uniform = true; % Is triangle stiffness uniform throughout?
        nTiles = 8; % Number of tiles in x and y
        initTris(nTiles); % Create triangles from regular grid of nodes
        
        side = false; % Is a side pinned? true = top or bottom pinned.
        less = false;
        threshold = 0.9; % Threshold for which nodes to pin.
        pinNodes(side,less,threshold); % Pin appropriate nodes for scene.
        
		dt = 1e-3; % time step
		tEnd = 1.0; % end time
		drawHz = 100; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 2.0; % viscous damping
		E = .3e2; % Young's modulus
		nu = 0.4; % Poisson's ratio
    case 3
        % 8x8 with top nodes fixed while reaching roughly y = -1.5 with two
        % different values for the Poisson ratio
        uniform = true; % Is triangle stiffness uniform throughout?
        nTiles = 8; % Number of tiles in x and y
        initTris(nTiles); % Create triangles from regular grid of nodes
        
        side = false; % Is a side pinned? true = top or bottom pinned.
        less = false;
        threshold = 0.9; % Threshold for which nodes to pin.
        pinNodes(side,less,threshold); % Pin appropriate nodes for scene.
        
		dt = 1e-3; % time step
		tEnd = 1.0; % end time
		drawHz = 100; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 1e0; % area density
		damping = 2.0; % viscous damping
		E = .3e2; % Young's modulus
		nu = -0.4; % Poisson's ratio. Negative ratio causes mesh to expand in width when stretched.
%         nu = 0.3; % Poisson's ratio 2. Lowering the value from .4 causes the mesh to stretch further without thinning out as much.
    case 4
        % 8x8 with left nodes fixed
        uniform = true; % Is triangle stiffness uniform throughout?
        nTiles = 8; % Number of tiles in x and y
        initTris(nTiles); % Create triangles from regular grid of nodes
        
        side = true; % Is a side pinned? true = top or bottom pinned.
        less = true;
        threshold = -0.9; % Threshold for which nodes to pin.
        pinNodes(side,less,threshold); % Pin appropriate nodes for scene.
        
		dt = 1e-3; % time step
		tEnd = 1.0; % end time
		drawHz = 100; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 11e-2; % area density
		damping = 2.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = -0.4; % Poisson's ratio
    case 5
        % 8x8 with different subset of nodes pinned and non uniform
        % stiffness
        uniform = false; % Is triangle stiffness uniform throughout?
        nTiles = 8; % Number of tiles in x and y
        initTris(nTiles); % Create triangles from regular grid of nodes
        
%         % Pin top nodes and left nodes
%         % Pin top nodes
% 		side = false; % Is a side pinned? true = top or bottom pinned.
%         less = false;
%         threshold = 0.9; % Threshold for which nodes to pin.
%         pinNodes(side,less,threshold); % Pin appropriate nodes for scene.
%         % Pin left nodes
%         for k = 1 : nNodes
%             if nodes(k).X(1) < -0.9
%                 nodes(k).fixed = true;
%             end
%         end
        
%         % Pin left and right nodes
%         % Pin left nodes
%         side = true; % Is a side pinned? true = top or bottom pinned.
%         less = true;
%         threshold = -0.9; % Threshold for which nodes to pin.
%         pinNodes(side,less,threshold); % Pin appropriate nodes for scene.
%         % Pin right nodes
%         for k = 1 : nNodes
%             if nodes(k).X(1) > 0.9
%                 nodes(k).fixed = true;
%             end
%         end
        
        %Pin nodes along line y=x
        for k = 1 : nNodes
            if nodes(k).X(1) == nodes(k).X(2)
                nodes(k).fixed = true;
            else
                nodes(k).fixed = false;
            end
        end
        
        
		dt = 1e-3; % time step
		tEnd = 1.0; % end time
		drawHz = 100; % refresh rate
		grav = [0 -9.81]'; % gravity
		rho = 15e-1; % area density
		damping = 2.0; % viscous damping
		E = 1e2; % Young's modulus
		nu = -0.4; % Poisson's ratio
    case 6 
        % Extra Credit
        nTiles = 1; % Number of tiles in x and y
        initTris(nTiles); % Create triangles from regular grid of nodes
        
        
        dt = 1e-3; % time step
        tEnd = 5.0; % end time
        drawHz = 100; % refresh rate
        grav = [0 -10]'; % gravity
        rho = 5e0; % area density
        damping = 0.0; % viscous damping
        E = 10000; % Young's modulus
        nu = -0.0; % Poisson's ratio
end

% Convert to lambda and mu
% lambda = (E*nu)/((1 + nu)*(1 -2*nu));
% mu = E/(2*(1 + nu));

function initTris(nTiles)
    % Creates triangles from a regular grid of nodes
    [nodes,tris] = createSquare(nTiles);
    nNodes = length(nodes);
    nTris = length(tris);
end
        
% % Find some nodes to fix
% for k = 1 : nNodes
% 	if nodes(k).X(2) > 0.9
% 		nodes(k).fixed = true;
% 	else
% 		nodes(k).fixed = false;
% 	end
% end

function pinNodes(side,less,threshold)
    if ~less
        for k = 1 : nNodes
            if ~side
                target = nodes(k).X(2);
            else
                target = nodes(k).X(1);
            end
            if target > threshold
                nodes(k).fixed = true;
            else
                nodes(k).fixed = false;
            end
        end
    else
        for k = 1 : nNodes
            if ~side
                target = nodes(k).X(2);
            else
                target = nodes(k).X(1);
            end
            if target < threshold
                nodes(k).fixed = true;
            else
                nodes(k).fixed = false;
            end
        end
    end
end
% Compute triangle mass and distribute to vertices
for k = 1 : nTris
	tri = tris(k).nodes;
	Xa = nodes(tri(1)).X;
	Xb = nodes(tri(2)).X;
	Xc = nodes(tri(3)).X;
    triArea = abs(((Xa(1)*(Xb(2)-Xc(2))) + (Xb(1)*(Xc(2)-Xa(2))) + (Xc(1)*(Xa(2) - Xb(2))))/2); %From www.mathopenref.com/coordtrianglearea.html
	triMass = (rho*triArea)/3;
	nodes(tri(1)).m = nodes(tri(1)).m + triMass;
	nodes(tri(2)).m = nodes(tri(2)).m + triMass;
	nodes(tri(3)).m = nodes(tri(3)).m + triMass;
    tris(k).defInv = inv([
        Xb(1)-Xa(1)  Xc(1)-Xa(1)
        Xb(2)-Xa(2)  Xc(2)-Xa(2)
        ]);
    
    % Rand to get either 1 or 2. 1 -> triangle gets alternate
    % stiffness, 2-> keep original stiffness.
    if ~uniform
        r = randi(2);
        if r == 1
            newE = 2e2; % triangle is assigned new stiffness
        else
            newE = E;
        end
    else
        newE = E;
    end
    
    tris(k).E = newE;
    
end

% Simulation loop
t0 = -inf;
for t = 0 : dt : tEnd;
	% Draw scene
	if t - t0 > 1 / drawHz
		draw(t,nodes,tris);
		t0 = t;
	end

	% Gravity force
	for k = 1 : nNodes
		nodes(k).f = nodes(k).m*grav;
	end
	
	% FEM force
    for k = 1 : nTris
        
        triNodes = tris(k).nodes;
        xa = nodes(triNodes(1)).x;
        xb = nodes(triNodes(2)).x;
        xc = nodes(triNodes(3)).x;
        
        % Compute Deformation Gradient
        F = [
            xb(1)-xa(1)  xc(1)-xa(1)
            xb(2)-xa(2)  xc(2)-xa(2)
        ] * tris(k).defInv;
           
        % Calculate Green Strain
        gs = .5*(F.'*F - eye(2));
        
        % Calculate Lambda and Mu
        lambda = (tris(k).E*nu)/((1 + nu)*(1 -2*nu));
        mu = tris(k).E/(2*(1 + nu));
        
        % Calculate Piola Kirchoff Stress
        
        P = F*(2*mu*gs + lambda*trace(gs)*eye(2));
        
        % Calculate and assign true stress
        tris(k).stress = P*F/det(F);
        
        % Calculate edge vectors
        Vab = xb - xa;
        Vbc = xc - xb;
        Vca = xa - xc;
        
        % Calculate lengths of vectors
        Lab = norm(Vab);
        Lbc = norm(Vbc);
        Lca = norm(Vca);
        
        % Calculate Normals
        Nab = [
            Vab(2)
            -Vab(1)
        ]/Lab;
        Nbc = [
            Vbc(2)
            -Vbc(1)
        ]/Lbc;
        Nca = [
            Vca(2)
            -Vca(1)
        ]/Lca;
        
        % Calculate force on each edge
        Fab = -Lab*tris(k).stress*Nab;
        Fbc = -Lbc*tris(k).stress*Nbc;
        Fca = -Lca*tris(k).stress*Nca;
        
        % Distribute forces to each node 
        nodes(triNodes(1)).f = nodes(triNodes(1)).f + Fab/2 + Fca/2;
        nodes(triNodes(2)).f = nodes(triNodes(2)).f + Fab/2 + Fbc/2;
        nodes(triNodes(3)).f = nodes(triNodes(3)).f + Fbc/2 + Fca/2;
    end
	
	% Integrate velocity and position
    for k = 1 : nNodes
        if ~nodes(k).fixed
            nodes(k).v = (nodes(k).m*nodes(k).v + dt*nodes(k).f)/(nodes(k).m + dt*damping*nodes(k).m);
            nodes(k).x = nodes(k).x + dt*nodes(k).v;
        end
    end
	
end

if ~isempty(video)
	video.close();
end

end

%%
function draw(t,nodes,tris)

global video;

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
	axis equal;
	axis(1.5*[-1 1 -1 1]); % Change axis limits here
	grid on;
	view(2);
	colormap jet;
	caxis([0 25]); % Change color limits here (comment out for auto)
	cb = colorbar;
	ylabel(cb, 'stress')
	video = VideoWriter('output','MPEG-4');
	video.open();
end
cla;
hold on;

x = [nodes.x]'; % flattened positions
f = reshape([tris.nodes],3,length(tris))'; % flattened indices
stress = reshape([tris.stress],4,length(tris)); % flattened stress
col = round(max(abs(stress)))'; % max stress entry per triangle
patch('Faces',f,'Vertices',x,'FaceVertexCData',col,'FaceColor','flat');

str = sprintf('t = %.4f', t);
title(str);
drawnow;

frame = getframe(gcf);
video.writeVideo(frame);

end

%%
function [nodes,tris] = createSquare(n)

% Regular grid with center points
x = linspace(-1,1,n+1);
y = linspace(-1,1,n+1);
[X,Y] = meshgrid(x,y);
x = reshape(X,(n+1)*(n+1),1);
y = reshape(Y,(n+1)*(n+1),1);
dx = 1/n;
dy = 1/n;
xc = linspace(-1+dx,1-dx,n);
yc = linspace(-1+dy,1-dy,n);
[Xc,Yc] = meshgrid(xc,yc);
xc = reshape(Xc,n*n,1);
yc = reshape(Yc,n*n,1);
x = [x;xc];
y = [y;yc];

nodes = [];
for i = 1 : length(x)
	nodes(i).X = [x(i),y(i)]'; %#ok<*AGROW>
	nodes(i).x = nodes(i).X;
	nodes(i).v = [0 0]';
	nodes(i).m = 0;
	nodes(i).f = [0 0]';
	nodes(i).fixed = false;
end

ng = (n+1)*(n+1);
tris = [];
for i = 1 : n
	% Index of the lower left node of the ith column
	ki = (i-1)*(n+1) + 1;
	% Index of the first center node of the ith column
	kci = ng + (i-1)*n + 1;
	for j = 1 : n
		% Index of the lower left node of the jth row of the ith column
		kij = ki + j - 1;
		% Index of the center node of the jth row of the ith column
		kcij = kci + j - 1;
		% Create the four triangles
		tris(end+1).nodes = [kij,kij+(n+1),kcij];
		tris(end+1).nodes = [kij+(n+1),kij+(n+1)+1,kcij];
		tris(end+1).nodes = [kij+(n+1)+1,kij+1,kcij];
		tris(end+1).nodes = [kij+1,kij,kcij];
	end
end
for k = 1 : length(tris)
	tris(k).stress = zeros(2);
end

end
