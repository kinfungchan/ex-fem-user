import analyser
import benchmarks
import database

def benchmark_bulkWave():
    # Initialise Input File
    inputWave = benchmarks.Bulk_Wave_Input_Monolithic()
    WaveBlock = analyser.SubdomainSolution(inputWave)

    # Initialise Output 
    history = database.History(WaveBlock.coo, WaveBlock.n_nodes, WaveBlock.n_elem)
    animate = database.Animation(history, directory="gifs_bulkWave")

    while (WaveBlock.t < WaveBlock.tfinal):
        if WaveBlock.t == 0.0:
            WaveBlock.el_state_upd()
        WaveBlock.solveq()
        WaveBlock.el_state_upd()

        print(f"Time step {WaveBlock.n}: time = {WaveBlock.t:.6g}")

        history.append_timestep(WaveBlock.t, WaveBlock.coo, WaveBlock.a, WaveBlock.v, WaveBlock.u, WaveBlock.res_sxx, WaveBlock.res_syy, WaveBlock.res_sxy)

    print("Outputting Bulk Wave plots...")

    animate.filenames_mesh_displ_pv = animate.save_timesteps(0, (WaveBlock.n), 10, [history.displ], 0, [history.coordinates], [WaveBlock.conn], var_name='Displacement (X)', filenames=animate.filenames_mesh_vel, clim=(-1.5e-6, 3.5e-6), zoom=1.1)
    animate.create_gif("mesh_ux_pv.gif", animate.filenames_mesh_displ_pv)

    print("Bulk Wave Simulation is done")