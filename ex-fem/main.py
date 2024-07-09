import benchmarks

"""
This main module runs the simulation using the Explicit Finite Element Method
It runs a single domain problem (monolithic) with a single material

"""

def main():
    print("Welcome to ex_fem User!")

    # 1x1 elem cross section Monolithic Wave Propagation (X-Dir)
    benchmarks.benchmark_01()

if __name__ == '__main__':
    main()