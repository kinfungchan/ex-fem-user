import benchmarks

"""
This main module runs the simulation using the Explicit Finite Element Method
It runs a single domain problem (monolithic) with a single material

"""

def main():
    print("Welcome to ex_fem User!")

    ## Uncommnent Benchmark to Run

    ## 1) Simple 1D Wave Propagation
    # benchmarks.benchmark_01()

    ## 2) Chiappa Bulk Wave Propagation
    # benchmarks.benchmark_bulkWave()

    ## 3) Copper Taylor Impact Bar
    benchmarks.benchmark_taylorImpact()

if __name__ == '__main__':
    main()