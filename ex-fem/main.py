import benchmarks

"""
This main module runs the simulation using the Explicit Finite Element Method
It runs a single domain problem (monolithic) with a single material

"""

def main():
    print("Welcome to ex_fem User!")

    benchmarks.benchmark_01()
    benchmarks.benchmark_bulkWave()

if __name__ == '__main__':
    main()