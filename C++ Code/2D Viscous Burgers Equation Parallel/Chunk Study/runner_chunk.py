from argparse import ArgumentParser
from csv import writer
from os import sched_getaffinity
from subprocess import run
from sys import stdout
from subprocess import PIPE


def cpu_count():
    return len(sched_getaffinity(0))


def main(N_range, n_range, t_range, c_range, optimizations, repeats_for_average):
    """Run through each set of options and compile/benchmark the test code."""
    # first see if our CSV exists
    try:
        with open('results.csv', 'x') as f:  # open in exclusive mode, fails if already there
            csv = writer(f)
            # write the header row
            csv.writerow(['Optimization', 'N (# of x,y steps)', 'n (# of t steps)', '# of threads', '# of chunk size', '# of samples', 'total time','average time'])
    except FileExistsError:
        pass

    # for each level of optimization
    for opt in optimizations:
        # compile the program, stop if it errors
        run(["g++", "2D_Viscous_Burgers_Arr_Parallel_Chunk.cpp", f"-O{opt}", "-fopenmp"], check=True)
        for c in c_range:
            for N in N_range:
                for n in n_range:
                    for threads in t_range:
                        # print what we're working on
                        print(f"opt={opt}, N={N}, n={n}, threads={threads}, c={c}", end="")
                        # Python tends to buffer until a newline, so we tell it not to
                        stdout.flush()
                        # please make sure the format is correct here
                        test_command = ['./a.out', str(N), str(n), str(threads), str(c)]
                        time_total = 0
                        value_total = 0
                        for samples in range(1, repeats_for_average + 1):
                            print(".", end="")
                            # print that we're working on a new job
                            stdout.flush()
                            # Python tends to buffer until a newline, so we tell it not to
                            output = run(test_command, check=True, stdout=PIPE).stdout.split()
                            # fetch the output of our program and split by whitespace
                            time_total += float(output[0])  # grab the time it took from the first word
                        # print the data
                        print(f"\tTotal Time={time_total}", end="")
                        print(f"\tAverage Time={time_total / samples:.06f}")
                        # then throw it in a CSV
                        with open('results.csv', 'a') as f:
                            csv = writer(f)
                            csv.writerow([opt, N, n, threads, c, samples, time_total, time_total / samples])


if __name__ == '__main__':
    parser = ArgumentParser()
    # make it so we can do multiple optimization levels
    parser.add_argument('-O', action='append', type=int, dest='optimizations', default=[])
    # store the left and right cutoffs
    parser.add_argument('-N_min', action='store', type=int, default=32)
    parser.add_argument('-N_max', action='store', type=int, default=64)
    # the range of segment counts to test
    parser.add_argument('-n_min', action='store', type=int, default=10000)
    parser.add_argument('-n_max', action='store', type=int, default=20000)
    # the range of thread counts to test
    parser.add_argument('-t_min', action='store', type=int, default=1)
    parser.add_argument('-t_max', action='store', type=int, default=cpu_count() * 2)
    # the range of chunk sizes
    parser.add_argument('-c_min', action='store', type=int, default=1)
    parser.add_argument('-c_max', action='store', type=int, default=100)
    # the number of times to repeat a test for average value
    parser.add_argument('-r', action='store', type=int, default=6, dest='repeats')

    # parse arguments
    args = parser.parse_args()
    if not args.optimizations:
        args.optimizations.append(0)
    main(
        range(args.N_min, args.N_max + 1),
        range(args.n_min, args.n_max + 1),
        range(args.t_min, args.t_max + 1),
        range(args.c_min, args.c_max + 1, 1),
        args.optimizations,
        args.repeats
    )
