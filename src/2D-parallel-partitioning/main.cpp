#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <mpich/mpi.h>
#include "../headers/Settings.hpp"

#define debug
#define usingGraphics

union Person {
	struct {
		uint8_t isInfected:1;
		uint8_t isImune:1;
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
int imunityPercentage = settings.getImmunityPercentage();
int loseImunityPercentage = settings.getLoseImmunityPercentage();
int vaccinationPercentage = settings.getVaccinationPercentage();
int deathPercentage = settings.getDeathPercentage();

int millisecondsToWaitForEachGeneration = settings.getMillisecodsToWaitForEachGeneration();



int rank, left, right, size;
int innerRows = settings.getMatrixSize() / 2;
int innerCols = settings.getMatrixSize() / 2;
int subRows = innerRows + 2;
int subCols = innerCols + 2;
int inner_grid_size = innerRows * innerCols;

Person * readMatrix = new Person[subRows * subCols];
Person * writeMatrix = new Person[subRows * subCols];

inline void initialize();

#ifdef usingGraphics
    ALLEGRO_DISPLAY * display = NULL;

    // get color from settings
    ALLEGRO_COLOR defaultPersonColor = al_map_rgb(settings.getDefaultPersonColor().r, settings.getDefaultPersonColor().g, settings.getDefaultPersonColor().b);
    ALLEGRO_COLOR infectedColor = al_map_rgb(settings.getInfectedColor().r, settings.getInfectedColor().g, settings.getInfectedColor().b);
    ALLEGRO_COLOR imuneColor = al_map_rgb(settings.getImmuneColor().r, settings.getImmuneColor().g, settings.getImmuneColor().b);
    ALLEGRO_COLOR deadColor = al_map_rgb(settings.getDeadColor().r, settings.getDeadColor().g, settings.getDeadColor().b);
    ALLEGRO_COLOR vaccinatedColor = al_map_rgb(settings.getVaccinatedColor().r, settings.getVaccinatedColor().g, settings.getVaccinatedColor().b);
    ALLEGRO_COLOR incubationColor = al_map_rgb(settings.getIncubationColor().r, settings.getIncubationColor().g, settings.getIncubationColor().b);
#endif //usignGraphics

MPI_Datatype column_t;
MPI_Datatype row_t;
MPI_Datatype corner_t;
MPI_Datatype subMatrixType;

inline void sendRows();
inline void sendCols();
inline void receiveRows();
inline void receiveCols();
inline void sendCorners();
inline void receiveCorners();
inline void update();
inline void updateBorders();
inline void draw(Person * readMatrix);
inline void swap();
inline void finalize();


inline int m(int i, int j) {return j * rows + i;}
inline int mm(int i, int j) {return j * subRows + i;}






int main(int argc, char * argv[])
{
    int source, dest;
    double elapsedTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == root) elapsedTime = MPI_Wtime();

    // create subMatrixType that will be used to send the matrix to the root process being aware of the borders
    MPI_Type_vector(innerRows * 2, innerCols * 2, cols, MPI_UNSIGNED_SHORT, &subMatrixType);
    MPI_Type_commit(&subMatrixType);

    MPI_Type_vector(innerRows, 1, cols, MPI_UNSIGNED_SHORT, &column_t);
    MPI_Type_commit(&column_t);

    MPI_Type_vector(1, innerCols, cols, MPI_UNSIGNED_SHORT, &row_t);
    MPI_Type_commit(&row_t);

    MPI_Type_contiguous(1, MPI_UNSIGNED_SHORT, &corner_t);
    MPI_Type_commit(&corner_t);

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

    for (int i = 1; i < innerRows; ++i)
    {
        for (int j = 1; j < innerCols; ++j)
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
        if (rank == root) MPI_Isend(&readMatrix[m(1,1)], 1, subMatrixType, root, 500, MPI_COMM_WORLD, &request);
        if (rank == 1) MPI_Isend(&readMatrix[m(1,1)], 1, subMatrixType, root, 501, MPI_COMM_WORLD, &request);
        if (rank == 2) MPI_Isend(&readMatrix[m(1,1)], 1, subMatrixType, root, 502, MPI_COMM_WORLD, &request);
        if (rank == 3) MPI_Isend(&readMatrix[m(1,1)], 1, subMatrixType, root, 503, MPI_COMM_WORLD, &request);

        MPI_Barrier(MPI_COMM_WORLD);
        
        if (rank == root)
        {
            MPI_Recv(&wholeMatrix[m(0, 0)], 1, subMatrixType, root, 500, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&wholeMatrix[m(0, innerCols-1)], 1, subMatrixType, 1, 501, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&wholeMatrix[m(innerRows-1, 0)], 1, subMatrixType, 2, 502, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&wholeMatrix[m(innerRows-1, innerCols-1)], 1, subMatrixType, 3, 503, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            
            draw(wholeMatrix);
        }
        #endif // usingGraphics

        sendCols();
        sendRows();
        sendCorners();
        update();
        receiveCols();
        receiveRows();
        receiveCorners();
        updateBorders();

        swap();

        MPI_Barrier(MPI_COMM_WORLD);

        sleep(millisecondsToWaitForEachGeneration);
    }

    MPI_Barrier(MPI_COMM_WORLD);

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
    for (int i = 0; i < subRows; ++i)
    {
        for (int j = 0; j < subCols; ++j)
        {
            readMatrix[m(i,j)].values.age = rand() % 100;
        }
    }

    if ( rank == 0)
        readMatrix[m(innerRows-2, innerCols/2-2)].values.isInfected = 1;
}

inline void finalize()
{
    MPI_Type_free(&column_t);
    MPI_Type_free(&row_t);
    MPI_Type_free(&corner_t);
    MPI_Type_free(&subMatrixType);

    delete [] readMatrix;
    delete [] writeMatrix;

    MPI_Finalize();
}


inline void update()
{
    for(int i = 2; i < innerRows-1; i++)
    {
        for(int j = 2; j < innerCols-1; j++)
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
                        if(readMatrix[m(i + k, j + l)].values.isInfected == true &&
                           readMatrix[m(i + k, j + l)].values.daysOfIncubation >= 2)
                        {
                            ++infectedNeighbours;
                        }

                        // if the neighbour is vaccinated
                        if(readMatrix[m(i + k, j + l)].values.isVaccinated == true)
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
            if (readMatrix[m(i,j)].values.isInfected && !readMatrix[m(i,j)].values.isImune)
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
                else if (rand()%100 < imunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isInfected = false;
                    writeMatrix[m(i,j)].values.isImune = true;
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
                    readMatrix[m(i,j)].values.isImune == false)
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

                if(rand()%100 < loseImunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isImune = false;
                    continue;
                }
            }
        }
    }
}

inline void updateBorders()
{
    // update the top and bottom border
    for(int i = 1; i < innerRows; ++i)
    {
        for(int j = 1; j < innerCols; j += innerCols - 1)
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
                        if(readMatrix[m(i + k, j + l)].values.isInfected == true &&
                           readMatrix[m(i + k, j + l)].values.daysOfIncubation >= 2)
                        {
                            ++infectedNeighbours;
                        }

                        // if the neighbour is vaccinated
                        if(readMatrix[m(i + k, j + l)].values.isVaccinated == true)
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
            if (readMatrix[m(i,j)].values.isInfected && !readMatrix[m(i,j)].values.isImune)
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
                else if (rand()%100 < imunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isInfected = false;
                    writeMatrix[m(i,j)].values.isImune = true;
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
                    readMatrix[m(i,j)].values.isImune == false)
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

                if(rand()%100 < loseImunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isImune = false;
                    continue;
                }
            }
        }
    }

    // update the left and right border
    for(int i = 1; i < innerRows; i += innerRows - 1)
    {
        for(int j = 1; j < innerCols; j++)
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
                        if(readMatrix[m(i + k, j + l)].values.isInfected == true &&
                           readMatrix[m(i + k, j + l)].values.daysOfIncubation >= 2)
                        {
                            ++infectedNeighbours;
                        }

                        // if the neighbour is vaccinated
                        if(readMatrix[m(i + k, j + l)].values.isVaccinated == true)
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
            if (readMatrix[m(i,j)].values.isInfected && !readMatrix[m(i,j)].values.isImune)
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
                else if (rand()%100 < imunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isInfected = false;
                    writeMatrix[m(i,j)].values.isImune = true;
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
                    readMatrix[m(i,j)].values.isImune == false)
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

                if(rand()%100 < loseImunityPercentage)
                {
                    writeMatrix[m(i,j)].values.isImune = false;
                    continue;
                }
            }
        }
    }
}

inline void sendRows()
{
    if (rank == 0)
    {
        MPI_Request requests[2];

        // rank 0 send top row to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, row_t, 2, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 0 send bottom row to rank 1
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, row_t, 2, 0, MPI_COMM_WORLD, &requests[1]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
    }

    if (rank == 1)
    {
        MPI_Request requests[2];

        // rank 0 send top row to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, row_t, 3, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 0 send bottom row to rank 1
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, row_t, 3, 0, MPI_COMM_WORLD, &requests[1]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
    }

    if (rank == 2)
    {
        MPI_Request requests[2];

        // rank 0 send top row to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, row_t, 0, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 0 send bottom row to rank 1
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, row_t, 0, 0, MPI_COMM_WORLD, &requests[1]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
    }

    if (rank == 3)
    {
        MPI_Request requests[2];

        // rank 0 send top row to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, row_t, 1, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 0 send bottom row to rank 1
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, row_t, 1, 0, MPI_COMM_WORLD, &requests[1]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
    }

}

inline void receiveRows()
{
    if (rank == 0)
    {
        // rank 0 recive bottom row from rank 1
        MPI_Recv(&readMatrix[m(innerRows-1,2)], 1, row_t, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // rank 0 recive top row from rank 1
        MPI_Recv(&readMatrix[m(0,2)], 1, row_t, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 1)
    {
        // rank 1 recive bottom row from rank 0
        MPI_Recv(&readMatrix[m(innerRows-1,2)], 1, row_t, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // rank 1 recive top row from rank 0
        MPI_Recv(&readMatrix[m(0,2)], 1, row_t, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 2)
    {
        // rank 2 recive bottom row from rank 0
        MPI_Recv(&readMatrix[m(innerRows-1,2)], 1, row_t, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // rank 2 recive top row from rank 0
        MPI_Recv(&readMatrix[m(0,2)], 1, row_t, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 3)
    {
        // rank 2 recive bottom row from rank 0
        MPI_Recv(&readMatrix[m(innerRows-1,2)], 1, row_t, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // rank 2 recive top row from rank 0
        MPI_Recv(&readMatrix[m(0,2)], 1, row_t, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

inline void sendCols()
{
    if (rank == 0)
    {
        MPI_Request requests[2];

        // rank 0 send left col to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, column_t, 1, 0, MPI_COMM_WORLD, &requests[0]);

        //rank 0 send right col to rank 1
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, column_t, 1, 0, MPI_COMM_WORLD, &requests[1]);
    }

    if (rank == 1)
    {
        MPI_Request requests[2];

        // rank 1 send left col to rank 0
        MPI_Isend(&readMatrix[m(2,2)], 1, column_t, 0, 0, MPI_COMM_WORLD, &requests[0]);

        //rank 1 send right col to rank 0
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, column_t, 0, 0, MPI_COMM_WORLD, &requests[1]);
    }

    if (rank == 2)
    {
        MPI_Request requests[2];

        // rank 2 send left col to rank 3
        MPI_Isend(&readMatrix[m(2,2)], 1, column_t, 3, 0, MPI_COMM_WORLD, &requests[0]);

        //rank 2 send right col to rank 3
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, column_t, 3, 0, MPI_COMM_WORLD, &requests[1]);
    }

    if (rank == 3)
    {
        MPI_Request requests[2];

        // rank 3 send left col to rank 2
        MPI_Isend(&readMatrix[m(2,2)], 1, column_t, 2, 0, MPI_COMM_WORLD, &requests[0]);

        //rank 3 send right col to rank 2
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, column_t, 2, 0, MPI_COMM_WORLD, &requests[1]);
    }
}

inline void receiveCols()
{
    if (rank == 0)
    {
        //rank 1 recive left col from rank 0
        MPI_Recv(&readMatrix[m(2,0)], 1, column_t, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //rank 1 recive right col from rank 0
        MPI_Recv(&readMatrix[m(2,innerCols-1)], 1, column_t, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 1)
    {
        //rank 1 recive left col from rank 0
        MPI_Recv(&readMatrix[m(2,0)], 1, column_t, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //rank 1 recive right col from rank 0
        MPI_Recv(&readMatrix[m(2,innerCols-1)], 1, column_t, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 2)
    {
        //rank 2 recive left col from rank 3
        MPI_Recv(&readMatrix[m(2,0)], 1, column_t, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //rank 2 recive right col from rank 3
        MPI_Recv(&readMatrix[m(2,innerCols-1)], 1, column_t, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 3)
    {
        //rank 3 recive left col from rank 2
        MPI_Recv(&readMatrix[m(2,0)], 1, column_t, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //rank 3 recive right col from rank 2
        MPI_Recv(&readMatrix[m(2,innerCols-1)], 1, column_t, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

}

inline void sendCorners()
{
    if (rank == 0)
    {
        MPI_Request requests[4];

        // rank 0 send top left corner to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 0 send top right corner to rank 1
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 0 send bottom left corner to rank 2
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 0 send bottom right corner to rank 2
        MPI_Isend(&readMatrix[m(innerRows-3,innerCols-3)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }

    if (rank == 1)
    {
        MPI_Request requests[4];

        // rank 1 send top left corner to rank 0
        MPI_Isend(&readMatrix[m(2,2)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 1 send top right corner to rank 0
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 1 send bottom left corner to rank 3
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 1 send bottom right corner to rank 3
        MPI_Isend(&readMatrix[m(innerRows-3,innerCols-3)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }

    if (rank == 2)
    {
        MPI_Request requests[4];

        // rank 2 send top left corner to rank 0
        MPI_Isend(&readMatrix[m(2,2)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 2 send top right corner to rank 0
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 2 send bottom left corner to rank 3
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 2 send bottom right corner to rank 3
        MPI_Isend(&readMatrix[m(innerRows-3,innerCols-3)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }

    if (rank == 3)
    {
        MPI_Request requests[4];

        // rank 3 send top left corner to rank 1
        MPI_Isend(&readMatrix[m(2,2)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 3 send top right corner to rank 1
        MPI_Isend(&readMatrix[m(2,innerCols-3)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 3 send bottom left corner to rank 2
        MPI_Isend(&readMatrix[m(innerRows-3,2)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 3 send bottom right corner to rank 2
        MPI_Isend(&readMatrix[m(innerRows-3,innerCols-3)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }
}

inline void receiveCorners()
{
    if (rank == 0)
    {
        MPI_Request requests[4];

        // rank 0 recive top left corner from rank 1
        MPI_Irecv(&readMatrix[m(1,1)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 0 recive top right corner from rank 1
        MPI_Irecv(&readMatrix[m(1,innerCols-2)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 0 recive bottom left corner from rank 2
        MPI_Irecv(&readMatrix[m(innerRows-2,1)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 0 recive bottom right corner from rank 2
        MPI_Irecv(&readMatrix[m(innerRows-2,innerCols-2)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }

    if (rank == 1)
    {
        MPI_Request requests[4];

        // rank 1 recive top left corner from rank 0
        MPI_Irecv(&readMatrix[m(1,1)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 1 recive top right corner from rank 0
        MPI_Irecv(&readMatrix[m(1,innerCols-2)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 1 recive bottom left corner from rank 3
        MPI_Irecv(&readMatrix[m(innerRows-2,1)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 1 recive bottom right corner from rank 3
        MPI_Irecv(&readMatrix[m(innerRows-2,innerCols-2)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }

    if (rank == 2)
    {
        MPI_Request requests[4];

        // rank 2 recive top left corner from rank 0
        MPI_Irecv(&readMatrix[m(1,1)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 2 recive top right corner from rank 0
        MPI_Irecv(&readMatrix[m(1,innerCols-2)], 1, corner_t, 0, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 2 recive bottom left corner from rank 3
        MPI_Irecv(&readMatrix[m(innerRows-2,1)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 2 recive bottom right corner from rank 3
        MPI_Irecv(&readMatrix[m(innerRows-2,innerCols-2)], 1, corner_t, 3, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }

    if (rank == 3)
    {
        MPI_Request requests[4];

        // rank 3 recive top left corner from rank 1
        MPI_Irecv(&readMatrix[m(1,1)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[0]);

        // rank 3 recive top right corner from rank 1
        MPI_Irecv(&readMatrix[m(1,innerCols-2)], 1, corner_t, 1, 0, MPI_COMM_WORLD, &requests[1]);

        // rank 3 recive bottom left corner from rank 2
        MPI_Irecv(&readMatrix[m(innerRows-2,1)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[2]);

        // rank 3 recive bottom right corner from rank 2
        MPI_Irecv(&readMatrix[m(innerRows-2,innerCols-2)], 1, corner_t, 2, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Request_free(&requests[0]);
        MPI_Request_free(&requests[1]);
        MPI_Request_free(&requests[2]);
        MPI_Request_free(&requests[3]);
    }
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
            else if (readMatrix[m(i,j)].values.isImune)
            {
                al_draw_filled_rectangle(j * square, i * square, (j + 1) * square, (i + 1) * square, imuneColor);
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