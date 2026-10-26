from VehicleData import falcon9
simulate = True
time = 0.0
time_step = 0.005
simulation_max_time = 60
while simulate and time < simulation_max_time:
    time += time_step
    # falcon9.update(time_step)