import benchmarks

"""
This main module runs the simulation using the Explicit Finite Element Method

Helpful References:
 - Belytschko, Ted, et al. Nonlinear finite elements for continua and structures. John Wiley & Sons, 2014.
 - Bombace, Nicola. Dynamic adaptive concurrent multi-scale simulation of wave propagation in 3D media. Diss. University of Oxford, 2018.
 - Chan, Kin Fung. Efficient coupling methods for the dynamic modelling of heterogeneous systems. Diss. University of Oxford, 2025.
"""

def main():
    print("Welcome to ex_fem User!")

    benchmark_choices = {
        "1": ("Simple 1D Wave Propagation", benchmarks.benchmark_01),
        "2": ("Chiappa Bulk Wave Propagation", benchmarks.benchmark_bulkWave),
        "3": ("Copper Taylor Impact Bar", benchmarks.benchmark_taylorImpact),
    }

    print("Choose a benchmark:")
    for number, (name, _) in benchmark_choices.items():
        print(f"{number}) {name}")

    while True:
        choice = input("Enter 1, 2, or 3: ").strip()
        if choice in benchmark_choices:
            benchmark_choices[choice][1]()
            break
        print("Invalid choice. Please enter 1, 2, or 3.")

if __name__ == '__main__':
    main()