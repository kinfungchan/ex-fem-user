import analyser
import database
import benchmarks

def benchmark_taylorImpact():
    # Initialise Input File
    nx = 5
    ny = 50
    input_data = benchmarks.CopperImpactInput(nx=nx, ny=ny)
    bar = analyser.SubdomainSolution(input_data)

    history = database.History(bar.coo, bar.n_nodes, bar.n_elem)
    animate = database.Animation(history, directory="gifs_taylorImpact")
    
    # Apply velocity to all nodes above the rigid wall
    initial_vel_y = -227000.0
    for i in range(bar.n_nodes):
        if bar.coo[i][1] > 1e-6:
            bar.v[i][1] = initial_vel_y
            bar.v_prev[i][1] = initial_vel_y
    
    num_frames = 30
    output_dt = bar.tfinal / (num_frames - 1)
    next_out_t = 0.0
        
    while (bar.t < bar.tfinal):
        if bar.t == 0.0:
            bar.el_state_upd()
        bar.solveq()
        bar.el_state_upd()

        if bar.t >= next_out_t or bar.n == 1:
            history.append_timestep(bar.t, bar.coo, bar.a, bar.v, bar.u, 
                                    bar.res_sxx, bar.res_syy, bar.res_sxy)
            next_out_t += output_dt
        
        print(f"Time step {bar.n:5d}: time = {bar.t*1e6:.2f} us")

    print("Outputting Taylor Impact plots...")

    num_saved_steps = len(history.t)       
    animate.filenames_mesh_displ_pv = animate.save_timesteps(0, num_saved_steps-1, 1, [history.displ], 0, 
                                                                [history.coordinates], [bar.conn], var_name='Displacement (X)', 
                                                                filenames=animate.filenames_mesh_vel)
    animate.create_gif("mesh_ux_pv.gif", animate.filenames_mesh_displ_pv)

    print("Taylor Impact Bar Simulation is done")