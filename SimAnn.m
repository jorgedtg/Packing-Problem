clear;clc;close all;
% Image processing toolbox required

N = 100;                                % Number of disks

r1 = 1;                                 % Radii of disks
r2 = 0.8;

groups = [r2, r1];                      

nr1 = 0.8;                               % Proportion of disks of each group
nr2 = 1-nr1;

circA = N*(nr1*pi*r1^2 + nr2*pi*r2^2);   % Area of all circles

pD = 0.80;                               % Objective packing density           

Lx = sqrt(pi*(464/5)*(1/pD));            % Box lengths
Ly = Lx;

disks = zeros(N, 4);

disks(:, 1) = (Lx*rand(1, N))';          % Initial position of disks
disks(:, 2) = (Ly*rand(1, N))';
disks(:, 3) = sqrt(disks(:, 1).^2 + disks(:, 2).^2); % Radial position of disks
disks(:, 4) = [ones(1, round(nr1*N))*r1, ones(1, round(nr2*N))*r2]';    % Radii of disks

% Visualize initial positons
figure;
viscircles([disks(:, 1), disks(:, 2)], disks(:, 4), 'Color', 'k', LineWidth=1);
rectangle('Position',[0 0 Lx Ly])
axis equal

%%
% Parameters of the Simulated Annealing
close all;
itLim = 100000;
T0 = 100;
U = 1e5;
T = T0;
P=0;
US = [];

found = false;

figure;
while ~found
    for it = 1:itLim

        if U == 0
            found = true;
            break
        end

        i = randi([1, N]);
        xi = (disks(i, 1));
        yi = (disks(i, 2));
        
        % Calculate change in disk position accordin to the elastic
        % extrusion potential
        [dx, dy] = change(disks, i, Lx, Ly, N);

        if isnan(dx) || isnan(dy)
            break
        end

        nx = xi + dx;
        ny = yi + dy;
        
        % When two disks are overlaping by a very small amount, the change
        % could be really small and taken to be 0, so I implemented a
        % minimum change in the position
        if yi+dy == yi && dy~=0
            ny = yi + sign(dy)*eps(yi);
        end
        if xi+dx == xi && dx~=0
            nx = xi + sign(dx)*eps(xi);
        end
        
        % New values of moved disk
        nr = sqrt(nx^2 + ny^2);
        nR = disks(i, 4);
        nrow = [nx, ny, nr, nR];
        
        ndisks = disks;
        ndisks(i, :) = nrow;
        
        % Calculate the potential energy of all disks
        U = energy(disks, Lx, Ly, N);
        Un = energy(ndisks, Lx, Ly, N);
        
        % Visualize change
        if mod(it, 10000) == 0
            cla
            viscircles([disks(:, 1), disks(:, 2)], disks(:, 4), 'Color', 'k', LineWidth=1);
            hold on
            viscircles([disks(i, 1), disks(i, 2)], disks(i, 4), 'Color', 'r', LineWidth=1);
            rectangle('Position',[0 0 Lx Ly]);
            hold off
            axis equal
            title(string(U))
            xlabel(string(it) + "   T = " + string(T) + "   P = " + string(round(P, 2)))
            xlim([-1, Lx+1])
            ylim([-1, Ly+1])
            drawnow
            % pause(1/60)
            % if it == 1
            %     pause(5)
            % end
        end

        % Change in energy
        dU = Un-U;
        
        % Change energy according to probability
        if dU <= 0
            disks = ndisks;
        else
            T = temperature(it, itLim, T0);
            P = exp(-dU/(T));
            rn = rand();

            if rn < P
                disks = ndisks;
            end
        end
    end
    ichange = 0;
    % The authors of the paper I used to base the algorithm found that
    % changing the position of the smaller disks first was more efficient.
    % I did this only when the energy is relatively big, because I noticed
    % that when the energy is small it increased a lot and prevented the
    % energy from going down.
    if U > 1e-2
        g = 1;
        while g <= length(groups)
            group = disks(disks(:, 4) == groups(g), :);
            ng = length(group);
            get = false;
            for gi = 1:ng
                % Find first disk of group with U ~= 0
                ichange = groupEnergy(group, Lx, Ly, ng);
                if ichange ~= 0
                    dchange = group(ichange, :);
                    tochange = find(ismember(disks, dchange, 'rows'));
                    xnew = Lx*rand();
                    ynew = Ly*rand();
                    disks(tochange, :) = [xnew, ynew, sqrt(xnew^2 + ynew^2), dchange(4)];
                    get = true;
                    break;
                end
            end
            if get
                break
            end
            g = g+1;
        end
    end
    if U<1e-2
        ichange=100;
    end
    % If no disks from the same group are overlaping each other, then
    % change the disk with the highest potential energy.
    if ichange == 0
        Uij = zeros(1, N);
        n = zeros(1, N);
        for i = 1:N
            for j = 1:N
                if i == j
                    continue
                end
                dij = embdepth(disks, i, j, Lx, Ly);
                if dij > 0
                    n(i) = n(i) + 1;
                end
                Uij(i) = Uij(i) + dij^2;
            end
        end
        ichange = find(Uij==max(Uij));
        dchange = disks(ichange, :);
        xnew = Lx*rand();
        ynew = Ly*rand();
        disks(ichange, :) = [xnew, ynew, sqrt(xnew^2 + ynew^2), dchange(4)];
    end
end

% Ratio obtained
ratio = (nr1*N*pi*r1^2 + nr2*N*pi*r2^2)/(Lx*Ly)



function i = groupEnergy(D, Lx, Ly, N)

U = 0;
for i = 1:N
    for j = 0:N
        if i == j
            continue
        end
        dij = embdepth(D, i, j, Lx, Ly);

        U = U + dij^2;
        
        if U > 0
            return
        end
    end
end

i = 0;
return
end

function U = energy(D, Lx, Ly, N)

U = 0;
for i = 1:N
    for j = 0:N
        if i == j
            continue
        end
        dij = embdepth(D, i, j, Lx, Ly);

        U = U + dij^2;

    end
end
end


function [dx, dy] = change(D, i, Lx, Ly, N)
dx = 0;
dy = 0;
for j = 0:N
    
    if i == j
        continue
    end
    
    if j ~= 0
        dij = embdepth(D, i, j, Lx, Ly);
    end
    Dij = dist(D, i, j);

    % if Dij < 1e-1
    %     Dij = 1e-1;
    % end

    if j == 0
        if D(i, 1) + D(i, 4) > Lx
            dxij = -D(i, 1) - D(i, 4) + Lx;
        elseif D(i, 1) - D(i, 4) < 0
            dxij = -D(i, 1) + D(i, 4);
        else
            dxij = 0;
        end
        if D(i, 2) + D(i, 4) > Ly
            dyij = -D(i, 2) - D(i, 4) + Ly;
        elseif D(i, 2) - D(i, 4) < 0
            dyij = -D(i, 2) + D(i, 4);
        else
            dyij = 0;
        end
        dij = sqrt(dxij^2 + dyij^2);
    else
        
        % if D(i, 1) - D(j, 1) == 0
        %     dxij = (0.01)*dij/Dij;
        % else
            dxij = (D(i, 1) - D(j, 1))*dij/Dij;
        % end
        dyij = (D(i, 2) - D(j, 2))*dij/Dij;
    end
    % fprintf("i = " + string(i)+", j = " +string(j)+ "  dx = " + string(dxij) + ", dy = " + string(dyij) + "  dij= "+string(dij)+"\n")
    dx = dx + dxij;
    dy = dy + dyij;
end
% fprintf("\n")
end

function Dij = dist(D, i, j)

if j == 0
    Dij = abs(D(i, 3));
else
    Dij = sqrt( (D(i, 1) - D(j, 1))^2 + (D(i, 2) - D(j, 2))^2   );
end
end



function dij = embdepth(D, i, j, Lx, Ly)

if j == 0
    if D(i, 1) + D(i, 4) > Lx
        dijx = D(i, 1) + D(i, 4) - Lx;
    elseif D(i, 1) - D(i, 4) < 0
        dijx = D(i, 1) - D(i, 4);
    else
        dijx = 0;
    end
    if D(i, 2) + D(i, 4) > Ly
        dijy = D(i, 2) + D(i, 4) - Ly;
    elseif D(i, 2) - D(i, 4) < 0
        dijy = D(i, 2) - D(i, 4);
    else
        dijy = 0;
    end
    dij = sqrt(dijx^2 + dijy^2);
    return
else
    if sqrt((D(j, 1)-D(i, 1)).^2 + (D(j, 2)-D(i, 2)).^2) < D(i, 4) + D(j, 4)
        dij = -sqrt((D(j, 1)-D(i, 1)).^2 + (D(j, 2)-D(i, 2)).^2) + D(i, 4) + D(j, 4);
        return
    else
        dij = 0;
        return
    end
end


end


function T = temperature(it, itMax, T0)

T = T0*0.9985^(it/10);%T0*(1-(it)/itMax);

end
