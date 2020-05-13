import subprocess

folder = '/home/posypkin/exp/cpp/stts/experiments/'
nexperiments = 20

def run_serial():
    result = subprocess.run(['./kpser', './NxC_test_MPI_100_5000_001'], stdout=subprocess.PIPE, universal_newlines=True)
    return result.stdout

def run_parallel():
    result = subprocess.run(['mpirun', '-np', '4', folder + 'kpdpa', folder + 'NxC_test_MPI_100_5000_001'], stdout=subprocess.PIPE, universal_newlines=True)
    return result.stdout

res_serial = [run_serial() for i in range(1,nexperiments)]
res_parallel = [run_parallel() for i in range(1,nexperiments)]

res_serial = [float(i.strip()) for i in res_serial]
res_parallel = [float(i.strip()) for i in res_parallel]

print("Serial:")
print(res_serial)
print("Serial average = ", sum(res_serial)/nexperiments)

print("\nParallel:")
print(res_parallel)
print("Parallel average = ", sum(res_parallel)/nexperiments)
