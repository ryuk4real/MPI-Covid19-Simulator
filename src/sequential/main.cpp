#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <mpich/mpi.h>
#include "../headers/Settings.hpp"

#define debug
#define usingGraphics


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

Person * readMatrix = new Person[settings.getMatrixSize() * settings.getMatrixSize()];
Person * writeMatrix = new Person[settings.getMatrixSize() * settings.getMatrixSize()];

#ifdef usingGraphics

    ALLEGRO_DISPLAY * display;

    rgb infectedRGBColor = settings.getInfectedColor();
    rgb immuneRGBColor = settings.getImmuneColor();
    rgb incubationRGBColor = settings.getIncubationColor();
    rgb deadRGBColor = settings.getDeadColor();
    rgb vaccinatedRGBColor = settings.getVaccinatedColor();
    rgb defaultPersonRGBColor = settings.getDefaultPersonColor();

    ALLEGRO_COLOR infectedColor = al_map_rgb(infectedRGBColor.r, infectedRGBColor.g, infectedRGBColor.b);
    ALLEGRO_COLOR immuneColor = al_map_rgb(immuneRGBColor.r, immuneRGBColor.g, immuneRGBColor.b);
    ALLEGRO_COLOR incubationColor = al_map_rgb(incubationRGBColor.r, incubationRGBColor.g, incubationRGBColor.b);
    ALLEGRO_COLOR deadColor = al_map_rgb(deadRGBColor.r, deadRGBColor.g, deadRGBColor.b);
    ALLEGRO_COLOR vaccinatedColor = al_map_rgb(vaccinatedRGBColor.r, vaccinatedRGBColor.g, vaccinatedRGBColor.b);
    ALLEGRO_COLOR defaultPersonColor = al_map_rgb(defaultPersonRGBColor.r, defaultPersonRGBColor.g, defaultPersonRGBColor.b);

#endif //usingGraphics

#ifdef debug
    double startTime, endTime, totalTime;
#endif //debug








void initialize();
void update();
inline void swap();
void draw();
void finalize();
inline int m(int i, int j);




int main()
{
    MPI_Init(NULL, NULL);


    #ifdef usingGraphics

        al_init();
        al_init_primitives_addon();
        al_set_app_name("COVID-19 Simulation");

        display = al_create_display(cols * square, rows * square);

    #endif //usingGraphics

    #ifdef debug
        startTime = MPI_Wtime();
    #endif //debug

    initialize();

    for (int i = 0; i < numberOfGenerations; i++)
    {
        update();

        #ifdef usingGraphics
            draw();
        #endif //usingGraphics

        swap();

        sleep(millisecondsToWaitForEachGeneration);
    }

    #ifdef debug
        endTime = MPI_Wtime();
        totalTime = endTime - startTime;
        printf("Total time: %f\n", totalTime);
    #endif //debug

    finalize();

}



void initialize()
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            readMatrix[m(i,j)].values.isInfected = 0;
            readMatrix[m(i,j)].values.isImmune = 0;
            readMatrix[m(i,j)].values.isDead = 0;
            readMatrix[m(i,j)].values.isVaccinated = 0;
            readMatrix[m(i,j)].values.daysOfIncubation = 0;
            readMatrix[m(i,j)].values.daysOfInfection = 0;

            readMatrix[m(i,j)].values.age = rand() % 100;

        }
    }
    
    readMatrix[m(rows/2,cols/2)].values.isInfected = 1;
}

void update()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
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

inline void swap(){
    Person * tmp;
    tmp = readMatrix;
    readMatrix = writeMatrix;
    writeMatrix = tmp;
}

void draw()
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

void finalize()
{
    delete [] readMatrix;
    delete [] writeMatrix;
    MPI_Finalize();
}

inline int m(int i, int j)
{
    return ((i * cols) + j);
}