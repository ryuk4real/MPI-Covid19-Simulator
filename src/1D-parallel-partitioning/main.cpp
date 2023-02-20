#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <mpich/mpi.h>
#include "../headers/Settings.hpp"

#define debug
//#define usingGraphics

union Person {
	struct {
		uint8_t isInfected:1;
		uint8_t isImmune:1;
		uint8_t isDead:1;
		uint8_t isVaccinated:1;
		uint8_t daysOfIncubation:2;
		uint8_t daysOfInfection:3;
		uint8_t age:7;
	}values;
	
	unsigned short all;
};

Settings settings = Settings();

#define procs 4
#define root 0

int rows = settings.getMatrixSize();
int cols = settings.getMatrixSize();
int square = settings.getSquareSize();

int numberOfGenerations = settings.getNumberOfGenerations();

int infectionPercentage = settings.getInfectionPercentage();
int immunityPercentage = settings.getImmunityPercentage();
int loseImmunityPercentage = settings.getLoseImmunityPercentage();
int vaccinationPercentage = settings.getVaccinationPercentage();
int deathPercentage = settings.getDeathPercentage();

int millisecondsToWaitForEachGeneration = settings.getMillisecodsToWaitForEachGeneration();


Person * readMatrix = new Person[rows * ((cols / procs) + 2)];
Person * writeMatrix = new Person[rows * ((cols / procs) + 2)];

int rank, left, right, size;


inline void initialize();

#ifdef usingGraphics
    ALLEGRO_DISPLAY * display = NULL;

    // get color from settings
    ALLEGRO_COLOR defaultPersonColor = al_map_rgb(settings.getDefaultPersonColor().r, settings.getDefaultPersonColor().g, settings.getDefaultPersonColor().b);
    ALLEGRO_COLOR infectedColor = al_map_rgb(settings.getInfectedColor().r, settings.getInfectedColor().g, settings.getInfectedColor().b);
    ALLEGRO_COLOR immuneColor = al_map_rgb(settings.getImmuneColor().r, settings.getImmuneColor().g, settings.getImmuneColor().b);
    ALLEGRO_COLOR deadColor = al_map_rgb(settings.getDeadColor().r, settings.getDeadColor().g, settings.getDeadColor().b);
    ALLEGRO_COLOR vaccinatedColor = al_map_rgb(settings.getVaccinatedColor().r, settings.getVaccinatedColor().g, settings.getVaccinatedColor().b);
    ALLEGRO_COLOR incubationColor = al_map_rgb(settings.getIncubationColor().r, settings.getIncubationColor().g, settings.getIncubationColor().b);
#endif //usignGraphics

MPI_Datatype columnType;
MPI_Datatype subMatrixType;
MPI_Comm comm;

inline void sendBorders();
inline void receiveBorders();
inline void update();
inline void updateBorders();
inline void draw(Person * readMatrix);
inline void swap();
inline void finalize();

inline int m(int i, int j) {return j * rows + i;}






int main(int argc, char * argv[])
{
    int source, dest;
    double elapsedTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == root) elapsedTime = MPI_Wtime();

    MPI_Type_contiguous(rows, MPI_UNSIGNED_SHORT, &columnType);
    MPI_Type_commit(&columnType);

    int dims[1] = {size};
    int periods[1] = {1};

    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &comm);
    MPI_Cart_shift(comm, 0, 1, &left, &right);

    #ifdef usingGraphics

        Person * wholeMatrix;

        if (rank == root)
        {
            al_init();
            display = al_create_display(cols * square, rows * square);
            al_init_primitives_addon();
            al_set_app_name("Covid19 Simulation");

            Person * tmp = new Person[rows * cols];
            wholeMatrix = tmp;
        }

    #endif // usingGraphics

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 1; j < cols / procs + 1; ++j)
        {
            readMatrix[m(i, j)].all = 0;
            writeMatrix[m(i, j)].all = 0;
        }
    }

    srand(time(NULL) + rank);
    initialize();


    for (int i = 1; i <= numberOfGenerations; ++i)
    {
        #ifdef debug
            if (rank == root)
                printf("Generation %d\n", i);
        #endif // debug
        
        #ifdef usingGraphics
        MPI_Request request;
        MPI_Isend(&readMatrix[m(0, 2)], cols/procs*2, columnType, root, 0, comm, &request);

            /* Person * buffer = &readMatrix[m(0,1)];
            MPI_Gather(buffer, 1, subMatrixType, wholeMatrix, 1, subMatrixType, root, comm); */

            if (rank == root)
            {
                for (int i = 0; i < size; ++i)
                    MPI_Recv(&wholeMatrix[m(0, i * cols / procs)], cols/procs*2, columnType, i, 0, comm, MPI_STATUS_IGNORE);

                draw(wholeMatrix);
            }

        #endif // usingGraphics

        sendBorders();
        update();
        receiveBorders();
        updateBorders();

        swap();

        MPI_Barrier(comm);

        sleep(millisecondsToWaitForEachGeneration);
    }

    MPI_Barrier(comm);

    if (rank == root)
    {
        elapsedTime = MPI_Wtime() - elapsedTime;
        printf("Elapsed time: %f\n", elapsedTime);
    }

    #ifdef usingGraphics

    if (rank == root)
    {
        delete [] wholeMatrix;
        al_destroy_display(display);
    }

    #endif // usingGraphics

    
    finalize();

    return 0;
}
















inline void initialize()
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 1; j < cols / procs + 3; ++j)
        {
            readMatrix[m(i,j)].values.age = rand() % 100;
        }
    }

    if (size == 1) readMatrix[m(rows/2, cols/2)].values.isInfected = 1;
    else
    {
        int centerRank = size / 2;

        if (rank == centerRank)
        {   
            if (size % 2 == 0)
            {
                readMatrix[m(rows/2, 2)].values.isInfected = 1;
            }
            else
            {
                readMatrix[m(rows/2, cols/procs/2)].values.isInfected = 1;
            }
        }
    }
}

inline void finalize()
{
    MPI_Type_free(&columnType);
    MPI_Finalize();
}

inline void update()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 2; j < cols/procs; j++)
        {
            short infectedNeighbours = 0;
            short vaccinatedNeighbours = 0;

            for(int k = -1; k <= 1; k++)
            {
                for(int l = -1; l <= 1; l++)
                {
                    // if the neighbour is not the person itself
                    if(k != 0 || l != 0)
                    {
                        // if the neighbour is infected
                        if(readMatrix[m((i + k + rows) % rows, ((j + l + cols) % cols))].values.isInfected == true &&
                           readMatrix[m((i + k + rows) % rows, ((j + l + cols) % cols))].values.daysOfIncubation >= 2)
                        {
                            ++infectedNeighbours;
                        }

                        // if the neighbour is vaccinated
                        if(readMatrix[m((i + k + rows) % rows, ((j + l + cols) % cols))].values.isVaccinated == true)
                        {
                            ++vaccinatedNeighbours;
                        }
                    }
                }
            }

            // copy the person to the write matrix
            writeMatrix[m(i,j)] = readMatrix[m(i,j)];

            if (readMatrix[m(i,j)].values.isDead || readMatrix[m(i,j)].values.isVaccinated)
                continue;
            
            // if the person is infected
            if (readMatrix[m(i,j)].values.isInfected && !readMatrix[(i,j)].values.isImmune)
            {

                if (readMatrix[m(i,j)].values.daysOfIncubation < 3)
                {
                    ++writeMatrix[m(i,j)].values.daysOfIncubation;
                    continue;
                }
                
                if (readMatrix[m(i,j)].values.daysOfInfection < 7)
                {
                    writeMatrix[m(i,j)].values.daysOfIncubation = 0;
                    ++writeMatrix[m(i,j)].values.daysOfInfection;
                }
                else if (rand()%100 < immunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isInfected = false;
                    writeMatrix[m(i,j)].values.isImmune = true;
                }
                
                if (readMatrix[m(i,j)].values.age >= 65)
                {
                    if (rand()%100 < deathPercentage)
                    {
                        writeMatrix[m(i,j)].values.isInfected = false;
                        writeMatrix[m(i,j)].values.isDead = true;
                    }
                }
                else if (readMatrix[m(i,j)].values.age > 25 && readMatrix[m(i,j)].values.age < 65)
                {
                    if (rand()%100 < deathPercentage / 2)
                    {
                        writeMatrix[m(i,j)].values.isInfected = false;
                        writeMatrix[m(i,j)].values.isDead = true;
                    }
                }
                else
                {
                    if (rand()%100 < deathPercentage / 4)
                    {
                        writeMatrix[m(i,j)].values.isInfected = false;
                        writeMatrix[m(i,j)].values.isDead = true;
                    }
                }
            }
            else // if the person is not infected
            {
                if (infectedNeighbours > 0 &&
                    rand()%100 < infectionPercentage * infectedNeighbours &&
                    readMatrix[m(i,j)].values.isImmune == false)
                {
                    writeMatrix[m(i,j)].values.isInfected = true;
                    continue;
                }

                if (rand()%100000 < vaccinationPercentage)
                {
                    writeMatrix[m(i,j)].values.isVaccinated = true;
                    continue;
                }

                if (vaccinatedNeighbours > 0 &&
                    rand()%250 < vaccinationPercentage * vaccinatedNeighbours)
                {
                    writeMatrix[m(i,j)].values.isVaccinated = true;
                    continue;
                }

                if(rand()%100 < loseImmunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isImmune = false;
                    continue;
                }
            }
        }
    }
}

inline void updateBorders()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 1; j < cols/procs + 1; j += cols / procs - 1)
        {
            short infectedNeighbours = 0;
            short vaccinatedNeighbours = 0;

            for(int k = -1; k <= 1; ++k)
            {
                for(int l = -1; l <= 1; ++l)
                {
                    // if the neighbour is not the person itself
                    if(k != 0 || l != 0)
                    {
                        // if the neighbour is infected
                        if(readMatrix[m((i + k + rows) % rows, ((j + l + cols) % cols))].values.isInfected == true &&
                           readMatrix[m((i + k + rows) % rows, ((j + l + cols) % cols))].values.daysOfIncubation >= 2)
                        {
                            ++infectedNeighbours;
                        }

                        // if the neighbour is vaccinated
                        if(readMatrix[m((i + k + rows) % rows, ((j + l + cols) % cols))].values.isVaccinated == true)
                        {
                            ++vaccinatedNeighbours;
                        }
                    }
                }
            }

            // copy the person to the write matrix
            writeMatrix[m(i,j)] = readMatrix[m(i,j)];

            if (readMatrix[m(i,j)].values.isDead || readMatrix[m(i,j)].values.isVaccinated)
                continue;
            
            // if the person is infected
            if (readMatrix[m(i,j)].values.isInfected && !readMatrix[(i,j)].values.isImmune)
            {

                if (readMatrix[m(i,j)].values.daysOfIncubation < 3)
                {
                    ++writeMatrix[m(i,j)].values.daysOfIncubation;
                    continue;
                }
                
                if (readMatrix[m(i,j)].values.daysOfInfection < 7)
                {
                    writeMatrix[m(i,j)].values.daysOfIncubation = 0;
                    ++writeMatrix[m(i,j)].values.daysOfInfection;
                }
                else if (rand()%100 < immunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isInfected = false;
                    writeMatrix[m(i,j)].values.isImmune = true;
                }
                
                if (readMatrix[m(i,j)].values.age >= 65)
                {
                    if (rand()%100 < deathPercentage)
                    {
                        writeMatrix[m(i,j)].values.isInfected = false;
                        writeMatrix[m(i,j)].values.isDead = true;
                    }
                }
                else if (readMatrix[m(i,j)].values.age > 25 && readMatrix[m(i,j)].values.age < 65)
                {
                    if (rand()%100 < deathPercentage / 2)
                    {
                        writeMatrix[m(i,j)].values.isInfected = false;
                        writeMatrix[m(i,j)].values.isDead = true;
                    }
                }
                else
                {
                    if (rand()%100 < deathPercentage / 4)
                    {
                        writeMatrix[m(i,j)].values.isInfected = false;
                        writeMatrix[m(i,j)].values.isDead = true;
                    }
                }
            }
            else // if the person is not infected
            {
                if (infectedNeighbours > 0 &&
                    rand()%100 < infectionPercentage * infectedNeighbours &&
                    readMatrix[m(i,j)].values.isImmune == false)
                {
                    writeMatrix[m(i,j)].values.isInfected = true;
                    continue;
                }

                if (rand()%100000 < vaccinationPercentage)
                {
                    writeMatrix[m(i,j)].values.isVaccinated = true;
                    continue;
                }

                if (vaccinatedNeighbours > 0 &&
                    rand()%250 < vaccinationPercentage * vaccinatedNeighbours)
                {
                    writeMatrix[m(i,j)].values.isVaccinated = true;
                    continue;
                }

                if(rand()%100 < loseImmunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isImmune = false;
                    continue;
                }
            }
        }
    }
}

inline void sendBorders()
{
    MPI_Request request;

    // send column to the left
    MPI_Isend(&readMatrix[m(0, 1)], 2, columnType, left, 0, comm, &request);

    // send column to the right
    MPI_Isend(&readMatrix[m(0, cols/procs)], 2, columnType, right, 1, comm, &request);
}

inline void receiveBorders()
{
    // receive column from the left
    MPI_Recv(&readMatrix[m(0, 0)], 2, columnType, left, 1, comm, MPI_STATUS_IGNORE);

    // receive column from the right
    MPI_Recv(&readMatrix[m(0, cols/procs + 1)], 2, columnType, right, 0, comm, MPI_STATUS_IGNORE);
}

#ifdef usingGraphics
inline void draw(Person * readMatrix)
{
    al_clear_to_color(defaultPersonColor);

    for (int j = 0; j < rows; ++j)
    {
        for (int i = 0; i < cols; ++i)
        {
            if (readMatrix[m(i,j)].values.isInfected && readMatrix[m(i,j)].values.daysOfIncubation < 3)
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, incubationColor);
            }
            else if (readMatrix[m(i,j)].values.isInfected)
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, infectedColor);
            }
            else if (readMatrix[m(i,j)].values.isImmune)
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, immuneColor);
            }
            else if (readMatrix[m(i,j)].values.isDead)
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, deadColor);
            }
            else if (readMatrix[m(i,j)].values.isVaccinated)
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, vaccinatedColor);
            }
            else 
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, defaultPersonColor);
            }

            // add black layer that scales with person age
            al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, al_map_rgba(0, 0, 0, readMatrix[m(i,j)].values.age));
        }
    }

    al_flip_display();
}
#endif //usingGraphics

inline void swap()
{
    Person * tmp;
    tmp = readMatrix;
    readMatrix = writeMatrix;
    writeMatrix = tmp;
}