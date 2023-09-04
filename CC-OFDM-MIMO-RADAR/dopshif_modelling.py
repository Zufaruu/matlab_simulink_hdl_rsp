vel = 20
N = 32
init_lag = 128
init_distance = 3e8 * (init_lag * 0.1 * 10e-4 /256) / 2
distance_per_samples = 58.594
distance_grid = []
lag_grid = []
cur_vel_distance = init_distance
cur_grid_distance = init_distance
cur_lag = init_lag

distance_grid.append(init_distance)
for i in range(N-1):
    cur_vel_distance += vel
    if(cur_vel_distance < cur_grid_distance + distance_per_samples):
        distance_grid.append(cur_grid_distance)
        lag_grid.append(cur_lag)
    else:
        distance_grid.append(cur_grid_distance + distance_per_samples)
        lag_grid.append(cur_lag)
        cur_lag += 1
        cur_grid_distance += distance_per_samples

print((distance_grid))
print((lag_grid))

