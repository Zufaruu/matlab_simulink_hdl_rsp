vel = 50
N_observation = 32
T_symbol = 0.1 * 10e-3
pri = 0.2e-3
N_sample = 8192

init_lag = 128
init_distance = 3e8 * (init_lag * T_symbol / N_sample) / 2
distance_per_samples = 3e8 * (1 * T_symbol / N_sample) / 2
distance_grid = []
lag_grid = []
cur_vel_distance = init_distance
cur_grid_distance = init_distance
cur_lag = init_lag
lag_grid.append(init_lag)
distance_grid.append(cur_grid_distance)

for i in range(N_observation-1):
    cur_vel_distance += vel * pri
    distance_grid.append(round(cur_vel_distance,4)) 
    if(cur_vel_distance < cur_grid_distance + distance_per_samples):
        lag_grid.append(cur_lag)
    else:
        lag_grid.append(cur_lag)
        cur_lag += 1
        cur_grid_distance += distance_per_samples

print((distance_grid))
print((lag_grid))

