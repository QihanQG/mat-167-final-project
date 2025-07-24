load('fluidE.mat'); 

[m, n] = size(XX); 

dx = diff(XX, 1, 1);  
dy = diff(YY, 1, 1);  
phi = atan2(dy, dx);
phi_unwrapped = zeros(size(phi));

%loop through each time step ( columns )
for i = 1:n
    phi_unwrapped(:,i) = unwrap(phi(:,i));
end

delta_s = space_scale; %distance between two adjacent points in micrometers 
delta_t = time_scale;  %time between frame in seconds

s = (0:size(phi_unwrapped,1)-1)' * delta_s; %position where phi is defined
t = (0:n-1)' * delta_t; % When phi is calcualted.                    

figure;
pcolor(s, t, phi_unwrapped');
shading flat;
colormap(jet);
colorbar;
xlabel('(\mum): Position along the flagellum where \phi(s,t) is measured ', 'FontSize', 12);
ylabel('Time t (s): Time at which \phi(s,t) is calculated', 'FontSize', 12);
title('Kymograph of Flagellar Tangent Angles \phi(s,t)', 'FontSize', 14);
axis tight;
grid off;