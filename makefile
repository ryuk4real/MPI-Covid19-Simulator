CC = mpiCC
FLAGS = -O3 -std=c++17 -I/usr/include/allegro5 -L/usr/lib -lallegro -lallegro_primitives



all: sequential 1Dparallel 2Dparallel

sequential: sequential.c
	$(CC) src/sequential/main.cpp -O3 -o bin/COVID19-serial.out $(FLAGS)

1Dparallel: 1Dparallel.c
	$(CC) src/1D-parallel-partitioning/main.cpp -O3 -o bin/COVID19-1D-parallel.out $(FLAGS)

2Dparallel: 2Dparallel.c
	$(CC) src/2D-parallel-partitioning/main.cpp -O3 -o bin/COVID19-2D-parallel.out $(FLAGS)