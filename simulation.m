particles = [...
     50 250 0 0 0 0 -3e-3 50 4;...
    100 200 0 0 0 0  3e-3 50 4;...
    150 250 0 0 0 0  3e-3 50 4;...
     50 200 0 0 0 0 -3e-3 50 4;...
    100 250 0 0 0 0  3e-3 50 4;...
    150 200 0 0 0 0 -3e-3 50 4;...
    ];

color_map = {[1 0 0];[0 0 1]};
%color_map = {'red';'blue'};
    
dt = 1.0/60;  % 60 ticks a second
figure(1);
drawnow limitrate


while true
    qs = particles(:,7);
    ms = particles(:,8);
    r2reciplist = pdist(particles(:,1:2)) .^ -3;
    r2recip = squareform(r2reciplist);
    
    combx = nchoosek(particles(:,1), 2);
    dirx = squareform(combx(:,1)-combx(:,2));
    dirx = dirx - 2.*tril(dirx);
    
    comby = nchoosek(particles(:,2), 2);
    diry = squareform(comby(:,1)-comby(:,2));
    diry = diry - 2.*tril(diry);
    
    charge_force_x=diag(qs)*(r2recip.*dirx)*qs.*8.988e9;
    charge_force_y=diag(qs)*(r2recip.*diry)*qs.*8.988e9;
    gravity_force_x=diag(ms)*(r2recip.*dirx)*ms.*6.67e-34;
    gravity_force_y=diag(ms)*(r2recip.*diry)*ms.*6.67e-34;
    
    force = [charge_force_x + gravity_force_x, ...
        charge_force_y + gravity_force_y];
    
%     charge_force=diag(qs)*r2recip*qs.*8.988e9
%     gravity_force=diag(ms)*r2recip*ms.*6.67e-34

    particles(:,5:6) = force ./ repmat(particles(:,8), 1, 2);  % acceleration
    particles(:,3:4) = particles(:,3:4) + particles(:,5:6) .* dt;  % velocity
    particles(:,1:2) = particles(:,1:2) + particles(:,3:4) .* dt;  % position
    
    c = sign(particles(:,7));
    c(c==-1) = 2;
    c = color_map(c);
    c = cell2mat(c);
    scatter(particles(:,1),particles(:,2),pi.*particles(:,9).^2, c,'filled');
    axis([0 400 0 300]);
    
%     disp(particles(:,1:2));
    
    pause(dt);
end