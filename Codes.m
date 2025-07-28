%% Part 1

clear; clc; close all;
load('fluidE.mat');

[m,n] = size(XX);

disp('Data:');
disp(['Number of spatial points: ', num2str(m)]);
disp(['Number of time points: ', num2str(n)]);
disp(['Space scale (\mum): ', num2str(space_scale)]);
disp(['Time scale (seconds): ', num2str(time_scale)]);

%All time steps.
figure; clf;
for i = 1:n
    plot(XX(:,i), YY(:,i));
    hold on;
end
hold off;
axis equal;
grid on;
xlabel('X  (\mum)');
ylabel('Y  (\mum)');
title('All flagellar shapes over time');

% Plot first few time steps
figure; clf;
for i = 1:5
    plot(XX(:,i), YY(:,i));
    hold on;
end
hold off;
axis equal;
grid on;
xlabel('X  (\mum)');
ylabel('Y  (\mum)');
title('First 5 time steps');

%% Part 2
dx = diff(XX, 1, 1);  
dy = diff(YY, 1, 1);  
phi = atan2(dy, dx);
phi_unwrapped = zeros(size(phi));

for i = 1:n
    phi_unwrapped(:,i) = unwrap(phi(:,i));
end

delta_s = space_scale; % Distance between two adjacent points
delta_t = time_scale;  % Time between frames

s = (0:size(phi_unwrapped,1)-1)' * delta_s; % Position where phi is defined
t = (0:n-1)' * delta_t; % When phi is calculated                   

figure; clf;
pcolor(s, t, phi_unwrapped');
shading interp;
colormap(jet);
colorbar;
xlabel('(\mum): Position along the flagellum where \phi(s,t) is measured', 'FontSize', 12);
ylabel('Time t (s): Time at which \phi(s,t) is calculated', 'FontSize', 12);
title('Kymograph of Flagellar Tangent Angles \phi(s,t)', 'FontSize', 14);
axis tight;
grid off;

%% Part 3:

phi_mean = mean(phi_unwrapped, 2);
demean_phi = phi_unwrapped - phi_mean;
[U,S,V] = svd(demean_phi, 'econ');

figure; clf; hold on;
for k = 1:4
    plot(U(:,k))
end
hold off;

figure; clf; hold on;
for k = 1:4
    plot(V(:,k))
end
hold off;

s = diag(S);
figure;
plot(cumsum(s.^2)./sum(s.^2), 'bo', 'MarkerFaceColor','b','MarkerSize',6)
xlabel("Mode k");  
ylabel("Cumulative energy "); 
title("Strength of Singular Values")

%% Part 4: 

V_1 = S(1,1) * V(:,1);
V_2 = S(2,2) * V(:,2);
figure;
plot(V_1, V_2, '.','MarkerFacecolor','r','Markersize', 10);
hold on;
hold off;
xlabel('V_1');
ylabel('V_2 ');
title('V_1 vs V_2');
axis equal;
grid on;
hold on;

theta = unwrap(atan2(V_2, V_1));
LstSqr_1 = [cos(theta) sin(theta)];
Soln = (LstSqr_1' * LstSqr_1)\ (LstSqr_1' * V_1);
A_1 = Soln(1,1);
B_1 = Soln(2,1);
LstSqr_2 = [cos(theta) sin(theta)];
Soln_2 = (LstSqr_2' * LstSqr_2) \ (LstSqr_2' * V_2);
A_2 = Soln_2(1,1);
B_2 = Soln_2(2,1);
V_1_fit = A_1 * cos(theta) + B_1 * sin(theta);
V_2_fit = A_2 * cos(theta) + B_2 * sin(theta);
plot(V_1_fit, V_2_fit, 'g-', 'LineWidth', 2);
hold off;


%% Part 5:

phi_recon = zeros(size(phi_unwrapped));
for i = 1:size(phi_unwrapped, 2)
    phi_recon(:,i) = phi_mean + V_1_fit(i) * U(:,1) + V_2_fit(i) * U(:,2);
end

X_recon = zeros(m, n);
Y_recon = zeros(m, n);
figure;
for i = 1:10
    % Original positions
    x_og = zeros(m, 1);
    y_og = zeros(m, 1);
    for j = 2:m
        x_og(j) = x_og(j-1) + space_scale * cos(phi_unwrapped(j-1, i));
        y_og(j) = y_og(j-1) + space_scale * sin(phi_unwrapped(j-1, i));
    end
    
    % Reconstructed positions
    x_recon = zeros(m, 1);
    y_recon = zeros(m, 1);
    for j = 2:m
        x_recon(j) = x_recon(j-1) + space_scale * cos(phi_recon(j-1, i));
        y_recon(j) = y_recon(j-1) + space_scale * sin(phi_recon(j-1, i));
    end
    
    plot(x_og, y_og, 'b-');
    hold on;
    plot(x_recon, y_recon, 'r--');
end
hold off;
xlabel('X');
ylabel('Y');
title('Original vs Reconstructed Shapes');
legend('Original', 'Reconstructed');
axis equal;
grid on;