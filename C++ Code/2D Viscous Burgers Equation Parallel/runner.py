from argparse import ArgumentParser
from csv import writer
from itertools import product
from os import sched_getaffinity
from subprocess import run
from sys import stdout
from subprocess import PIPE


def cpu_count():
    return len(sched_getaffinity(0))


def main(N_range, n_range, t_range, optimizations, repeats_for_average, filename, use_ratio):
    """Run through each set of options and compile/benchmark the test code."""

    # for each level of optimization
    for opt in optimizations:
        # compile the program, stop if it errors
        modes = ["DEFAULT", "STATIC", "DYNAMIC", "GUIDED"]
        for mode in modes:
            csv_name = f'results_{mode.lower()}_{"weak" if use_ratio else "strong"}.csv'
            # first see if our CSV exists
            try:
                with open(csv_name, 'x') as f:  # open in exclusive mode, fails if already there
                    csv = writer(f)
                    # write the header row
                    csv.writerow(['Optimization', 'N (# of x,y steps)', 'n (# of t steps)', '# of threads', '# of samples', 'total time','average time'])
            except FileExistsError:
                pass

            run(["g++", filename, f"-O{opt}", "-fopenmp", "-D" + mode], check=True)
            for N, n in product(N_range, n_range):
                if use_ratio:
                    n = round(N / 32 * 1000)
                for threads in t_range:
                    # print what we're working on
                    print(f"opt={opt}, mode={mode}, N={N}, n={n}, threads={threads}", end="")
                    # Python tends to buffer until a newline, so we tell it not to
                    stdout.flush()
                    # please make sure the format is correct here
                    test_command = ['./a.out', str(N), str(n), str(threads)]
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
                    print(f"\tTotal Time={time_total:.06f}", end="")
                    print(f"\tAverage Time={time_total / samples:.06f}")
                    # then throw it in a CSV
                    with open(csv_name, 'a') as f:
                        csv = writer(f)
                        csv.writerow([opt, N, n, threads, samples, time_total, time_total / samples])


if __name__ == '__main__':
    parser = ArgumentParser()
    # make it so we can do multiple optimization levels
    parser.add_argument('-O', action='append', type=int, dest='optimizations', default=[])
    # store the left and right cutoffs
    parser.add_argument('-N_min', action='store', type=int, default=5)
    parser.add_argument('-N_max', action='store', type=int, default=12)
    # the range of segment counts to test
    parser.add_argument('-n_min', action='store', type=int, default=10000)
    parser.add_argument('-n_max', action='store', type=int, default=20000)
    # the range of thread counts to test
    parser.add_argument('-t_min', action='store', type=int, default=1)
    parser.add_argument('-t_max', action='store', type=int, default=cpu_count() * 2)
    # the number of times to repeat a test for average value
    parser.add_argument('-r', action='store', type=int, default=5, dest='repeats')
    # the number of times to repeat a test for average value
    parser.add_argument('-u', action='store_true', default=False, dest='use_ratio')
    # the file to use
    parser.add_argument('-f', action='store', type=str, default="2D_Viscous_Burgers_Arr_Parallel.cpp", dest='filename')

    # parse arguments
    args = parser.parse_args()
    if not args.optimizations:
        args.optimizations.append(0)
    N_min = args.N_min
    N_max = args.N_max
    n_min = args.n_min
    n_max = args.n_max
    main(
        [1 << x for x in range(N_min, N_max + 1)],
        range(n_min, n_max + 1),
        range(args.t_min, args.t_max + 1),
        args.optimizations,
        args.repeats,
        args.filename,
        args.use_ratio
    )
