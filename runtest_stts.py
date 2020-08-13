import subprocess

folder = '/home/sithuthantsin86/Dropbox/DSc/git/kp-dpa/inputs/NxC_test_MPI/'
nexperiments = 20

def run_serial():
    result = subprocess.run(['./s.out', folder + 'NxC_test_MPI_1600_80000_001'], stdout=subprocess.PIPE, universal_newlines=True)
    return result.stdout

def run_omp():
	result = subprocess.run(['./o.out', folder + 'NxC_test_MPI_1600_80000_001', '4'], stdout=subprocess.PIPE, universal_newlines=True)
	return result.stdout

#def run_parallel():
#    result = subprocess.run(['mpirun', '-np', '4', 'pk.out', folder + 'NxC_test_MPI_6400_320000_001'], stdout=subprocess.PIPE, universal_newlines=True)
#    return result.stdout

def run_parallel_mem_sav():
    result = subprocess.run(['mpirun', '-np', '4', 'pkm.out', folder + 'NxC_test_MPI_1600_80000_001'], stdout=subprocess.PIPE, universal_newlines=True)
    return result.stdout

res_serial = [run_serial() for i in range(1,nexperiments)]
res_omp = [run_omp() for i in range(1,nexperiments)]
#res_parallel = [run_parallel() for i in range(1,nexperiments)]
res_parallel_mem_sav = [run_parallel_mem_sav() for i in range(1,nexperiments)]

res_serial = [float(i.strip()) for i in res_serial]
res_omp = [float(i.strip()) for i in res_omp]
#res_parallel = [float(i.strip()) for i in res_parallel]
res_parallel_mem_sav = [float(i.strip()) for i in res_parallel_mem_sav]

print("Serial:")
print(res_serial)
print("Serial average = ", sum(res_serial)/nexperiments)

print("OpenMP:")
print(res_omp)
print("OpenMP average = ", sum(res_omp)/nexperiments)

#print("\nParallel:")
#print(res_parallel)
#print("Parallel average = ", sum(res_parallel)/nexperiments)

print("\nMPI:")
print(res_parallel_mem_sav)
print("MPI average = ", sum(res_parallel_mem_sav)/nexperiments)