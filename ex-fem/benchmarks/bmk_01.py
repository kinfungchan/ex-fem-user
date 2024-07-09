import analyser
import benchmarks
import database

def benchmark_01():
    # Initialise Input File
    input = benchmarks.SteelBarInput()
    bar = analyser.SubdomainSolution(input)

    # Initialise Output 
    history = database.History(bar.coo, bar.n_nodes, bar.n_elem)
    animate = database.Animation(history, directory="mono_gifs")

    while (bar.t < bar.tfinal):
        if bar.t == 0.0:
            bar.el_state_upd()
        bar.solveq()
        bar.el_state_upd()

        print(f"Time step {bar.n}: time = {bar.t:.6g}")

        singlePlot = database.Output(bar.t, bar.coo, bar.a, bar.v, bar.u, bar.res_sxx, bar.res_syy, bar.res_sxy)  
        history.append_timestep(bar.t, bar.coo, bar.a, bar.v, bar.u, bar.res_sxx, bar.res_syy, bar.res_sxy)

        if (bar.n % 20 == 0):
            # just list of x coordinates for y = 0
            coord = bar.coo[:251,0]
            vel = bar.v[:251,0]
            animate.save_single_plot(1, coord, vel, "velocity", "pos", "vel [m/s]", animate.filenames_vel, bar.n, bar.t)
            animate.save_mesh_plot(bar.coo, bar.conn, bar.v, 0, "velocity", "vel [m/s]", animate.filenames_mesh_vel, bar.n, bar.t)

    animate.create_gif("velocity.gif", animate.filenames_vel)
    animate.create_gif("mesh_vel.gif", animate.filenames_mesh_vel)
    node_number = 126
    history.plot_node_dof(node_number - 1, 0)
    history.plot_element_data(125)

    print("Monolithic Simulation 01 is done")