
clear;
source_start = 1;

c0 = 3E8;
eps_0 = 8.8542149e-12;          % vacuum permittivity
mu_0 = 1.2566370614e-6;         % vacuum permeability

%Space
Nz = 100;
zMax  = 2e-6;
dz = zMax/Nz;
z = linspace(0, zMax, Nz);

lambda = 600E-9;
freq = 2*pi*c0/lambda;

%time
Nt = 300;
dt = 1/5*1/freq;
tMax = dt * Nt;
t = linspace(0, tMax, Nt);

eps_r = 2+0.9*sin(freq/2.*t +pi/4);
mu_r = 2+0.9*sin(freq/2.*t);

ep=eps_r*eps_0;
mu= mu_r*mu_0;

E0 = sin(freq.*t);

Hy = zeros(1, Nz);
Ex = zeros(1, Nz);

for n = 1: Nt-1
   %Set time-dependant source
   Ex(source_start) = sin(freq*n*dt); 
    for k = 1: Nz-1
        Hy(k) =mu(n)/mu(n+1)* Hy(k) + dt/mu(n+1)*(Ex(k)-Ex(k+1))/dz;
    end

   for k = 2: Nz-1
        Ex(k) = ep(n)/ep(n+1)*Ex(k) + dt/ep(n+1)*(Hy(k-1)-Hy(k))/dz;
   end

    % Plotting
   plot3(0*z, z, Ex, 'b', 'LineWidth', 2);
   axis([-0.01 0.01 0 zMax -5 5]);
   hold on
   h = quiver3(0*z(1:3:Nz), z(1:3:Nz), 0*z(1:3:Nz), 0*z(1:3:Nz), 0*z(1:3:Nz), Ex(1:3:Nz), 0, 'c', 'LineWidth', 1);
   set (h, "maxheadsize", 0.0);
   plot3(Hy, z, 0*z, 'r', 'LineWidth', 2)
   h = quiver3(0*z(1:3:Nz), z(1:3:Nz), 0*z(1:3:Nz), Hy(1:3:Nz), 0*z(1:3:Nz), 0*z(1:3:Nz), 0, 'm', 'LineWidth', 1);
   set (h, "maxheadsize", 0.0);
   grid on
   myTitle1 =  sprintf("Simulation of Electromagnetic Field in 1 Dimension, Sinusoidal Permeability and Permittivity");
   title(myTitle1)
   xlabel("Magnetic Field (T)");
   zlabel("Electric Field (V/m)");
   ylabel("z (m)");
   pause(0.01)
   hold off

end