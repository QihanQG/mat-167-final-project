%% Part 2

load('fluidE.mat'); 

[m, n] = size(XX); 

dx = diff(XX, 1, 1);  
dy = diff(YY, 1, 1);  
phi = atan2(dy, dx);
phi_unwrapped = zeros(size(phi));

for i = 1:n
    phi_unwrapped(:,i) = unwrap(phi(:,i));
end

delta_s = space_scale; % Distance between two adjacent points in micrometers 
delta_t = time_scale;  % Time between frames in seconds

s = (0:size(phi_unwrapped,1)-1)' * delta_s; % Position where phi is defined
t = (0:n-1)' * delta_t; % When phi is calculated                   

figure;
pcolor(s, t, phi_unwrapped');
shading flat;
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
for k = 1:4
    plot(U(:,k))
    hold on
end
for k = 1:4
    figure(2)
    plot(V(:,k))
    hold on
end
s = diag(S);
figure(3)
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
for k = 1:size(phi_unwrapped, 2)
    phi_recon(:,k) = phi_mean + V_1_fit(k) * U(:,1) + V_2_fit(k) * U(:,2);
end

m_spatial = size(phi_unwrapped, 1) + 1;  
n_time = size(phi_unwrapped, 2);         
X_recon = zeros(m_spatial, n_time);
Y_recon = zeros(m_spatial, n_time);
period_length = min(30, n_time);

figure;
for k = 1:min(5, period_length)
    X_orig_k = [0; cumsum(space_scale * cos(phi_unwrapped(:,k)))];
    Y_orig_k = [0; cumsum(space_scale * sin(phi_unwrapped(:,k)))];
    
    X_recon_k = [0; cumsum(space_scale * cos(phi_recon(:,k)))];
    Y_recon_k = [0; cumsum(space_scale * sin(phi_recon(:,k)))];
    
    plot(X_orig_k, Y_orig_k, 'b-');  
    hold on;
    plot(X_recon_k, Y_recon_k, 'r--');  
end
xlabel('X  (\mum)');
ylabel('Y  (\mum)');
title('1st Period: Original vs Reconstructed Flagellar Shapes');
legend('Original', 'Reconstructed');
axis equal;
grid on;
hold off;