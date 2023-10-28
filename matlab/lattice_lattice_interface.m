clear;

left_x = -100;
right_x = 150;
left_y = -100;
right_y = 100;

m_1 = 1;
m_2 = 0.5;
c = 0.75;
omega = 1;
a = 1;

k_1 = asin(omega/2 * sqrt(m_1/c)) * 2 / a;

particles_num_x = left_x:a:right_x;
particles_num_y = left_y:a:right_y;

particles_mass = zeros(int32((-left_y+right_y+1)/a),int32((-left_x+right_x+1)/a));
for j=1:length(particles_num_x)
    for i=1:length(particles_num_y)
%        if j < length(particles_num_x)/2
        if j < -left_x
            particles_mass(i,j) = m_1;
        else
            particles_mass(i,j) = m_2;
        end
    end
end

angle_step = 3;
%for j=int32(length(particles_num_x)/2):length(particles_num_x)
for j=-left_x:length(particles_num_x)
    for i=angle_step:length(particles_num_y)
        particles_mass(i,j) = m_1;
    end
    angle_step = angle_step + 3;
end

particles_disp = zeros(int32((-left_y+right_y+1)/a),int32((-left_x+right_x+1)/a));
particles_vel = zeros(int32((-left_y+right_y+1)/a),int32((-left_x+right_x+1)/a));

beta = 0.1;
n_0 = -35;
u_0 = 1;
g_1 = a * sqrt(c/m_1) * cos(k_1*a/2);

for i=1:length(particles_num_x)
    for j=1:length(particles_num_y)
        if i < length(particles_num_x)
            particles_disp(j, i) = u_0 * exp(-beta^2/2 * (particles_num_x(i) - n_0)^2) * sin(particles_num_x(i) * k_1);
            particles_vel(j, i) = -u_0 * exp(-beta^2/2 * (particles_num_x(i) - n_0)^2) * (omega * cos(k_1*(particles_num_x(i))) - beta^2*g_1/a*(particles_num_x(i)-n_0)*sin(particles_num_x(i) * k_1));
        end
    end
end

figure; hold on
[X,Y] = meshgrid(particles_num_x,particles_num_y);
s1 = surf(X,Y,particles_disp,'FaceAlpha',0.9);
s2 = surf(X,Y,particles_mass,'FaceAlpha',0.5);
title('Цветовая карта перемещений в момент времени t = 0');
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
colorbar;
s1.EdgeColor = 'none';
s2.EdgeColor = 'none';
view(17,22);
hold off


dt = 0.001;
t_max = 150;
times = 0:dt:t_max;

for t=times
    for ind_x=1:length(particles_num_x)
        for ind_y=1:length(particles_num_y)
            if ind_x ~= 1 && ind_y ~= 1 && ind_x ~= length(particles_num_x) && ind_y ~= length(particles_num_y)
                particles_vel(ind_y,ind_x) = particles_vel(ind_y,ind_x) + c/particles_mass(ind_y, ind_x)*(particles_disp(ind_y-1, ind_x)+particles_disp(ind_y+1, ind_x)+particles_disp(ind_y, ind_x-1)+particles_disp(ind_y, ind_x+1)-4*particles_disp(ind_y, ind_x))*dt;
                particles_disp(ind_y, ind_x) = particles_disp(ind_y, ind_x) + particles_vel(ind_y, ind_x) * dt;
            end
            if ind_y == 1
                particles_vel(ind_y,ind_x) = particles_vel(ind_y+1,ind_x);
                particles_disp(ind_y,ind_x) = particles_disp(ind_y+1,ind_x);
            end
            if ind_y == length(particles_num_y)
                particles_vel(ind_y,ind_x) = particles_vel(ind_y-1,ind_x);
                particles_disp(ind_y,ind_x) = particles_disp(ind_y-1,ind_x);
            end
        end
    end
    if t==60
        figure; hold on
        s1 = surf(X,Y,particles_disp,'FaceAlpha',0.9);
        s2 = surf(X,Y,particles_mass,'FaceAlpha',0.5);
        title('Цветовая карта перемещений в момент времени t = 60');
        xlabel('Номер частицы по оси Ox');
        ylabel('Номер частицы по оси Oy');
        colorbar;
        s1.EdgeColor = 'none';
        s2.EdgeColor = 'none';
        view(17,22);
        hold off
    end
    if t==120
        figure; hold on
        s1 = surf(X,Y,particles_disp,'FaceAlpha',0.9);
        s2 = surf(X,Y,particles_mass,'FaceAlpha',0.5);
        title('Цветовая карта перемещений в момент времени t = 120');
        xlabel('Номер частицы по оси Ox');
        ylabel('Номер частицы по оси Oy');
        colorbar;
        s1.EdgeColor = 'none';
        s2.EdgeColor = 'none';
        view(17,22);
        hold off
    end
        
end
figure; hold on
s1 = surf(X,Y,particles_disp,'FaceAlpha',0.9);
s2 = surf(X,Y,particles_mass,'FaceAlpha',0.5);
title('Цветовая карта перемещений в момент времени t = 150');
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
colorbar;
s1.EdgeColor = 'none';
s2.EdgeColor = 'none';
view(17,22);
hold off
