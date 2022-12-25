[time,temp] = meshgrid(-1.414:0.01:+1.414);
Y = 79.94+0.99497.*time+0.5152.*temp+0.25.*time.*temp-1.37625.*time.^2-1*temp.^2;
meshc(time,temp, Y)
colorbar
