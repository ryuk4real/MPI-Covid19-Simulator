#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <fstream> //ifstream
#include <iostream> //cout

#include "json.hpp" // parsing json file
using json = nlohmann::json;

struct rgb{int r; int g; int b;};

class Settings
{

    private:

        int numberOfGenerations;

        int matrixSize;

        rgb defaultPersonColor;

        rgb immuneColor;

        rgb infectedColor;

        rgb deadColor;

        rgb vaccinatedColor;

        int millisecondsToWaitForEachGeneration;


    public:

        Settings();

        // Getters -------------------------------------------------------------------------------------

        int getNumberOfGenerations() const {return this->numberOfGenerations;}

        int getMatrixSize() const {return this->matrixSize;}

        rgb getDefaultPersonColor() const {return this->defaultPersonColor;}

        rgb getImmuneColor() const {return this->immuneColor;}

        rgb getInfectedColor() const {return this->infectedColor;}

        rgb getDeadColor() const {return this->deadColor;}

        rgb getVaccinatedColor() const {return this->vaccinatedColor;}

        int getMillisecodsToWaitForEachGeneration() const {return this->millisecondsToWaitForEachGeneration;}

        // ---------------------------------------------------------------------------------------------

        // Utils ---------------------------------------------------------------------------------------
        inline int checkRGBValue(int value) const;
        inline int checkPositive(int value) const;

};

Settings::Settings()
{

    // Read json settings file -------------------------------------------------------------------------
    std::ifstream settings("src/settings/settings.json");

    if(!settings){ throw std::runtime_error("ERROR: Couldn't open/find settings.json file"); }
    
    json jsonSettings;
    settings >> jsonSettings;
    //----------------------------------------------------------------------------------------------------

    numberOfGenerations = checkPositive(jsonSettings["numberOfGenerations"]);

    matrixSize = checkPositive(jsonSettings["matrixSize"]);

    defaultPersonColor.r = checkRGBValue(jsonSettings["defaultPersonColor"][0]);
    defaultPersonColor.g = checkRGBValue(jsonSettings["defaultPersonColor"][1]);
    defaultPersonColor.b = checkRGBValue(jsonSettings["defaultPersonColor"][2]);

    immuneColor.r = checkRGBValue(jsonSettings["immuneColor"][0]);
    immuneColor.g = checkRGBValue(jsonSettings["immuneColor"][1]);
    immuneColor.b = checkRGBValue(jsonSettings["immuneColor"][2]);

    infectedColor.r = checkRGBValue(jsonSettings["infectedColor"][0]);
    infectedColor.g = checkRGBValue(jsonSettings["infectedColor"][1]);
    infectedColor.b = checkRGBValue(jsonSettings["infectedColor"][2]);

    deadColor.r = checkRGBValue(jsonSettings["deadColor"][0]);
    deadColor.g = checkRGBValue(jsonSettings["deadColor"][1]);
    deadColor.b = checkRGBValue(jsonSettings["deadColor"][2]);

    vaccinatedColor.r = checkRGBValue(jsonSettings["vaccinatedColor"][0]);
    vaccinatedColor.g = checkRGBValue(jsonSettings["vaccinatedColor"][1]);
    vaccinatedColor.b = checkRGBValue(jsonSettings["vaccinatedColor"][2]);

    millisecondsToWaitForEachGeneration = checkPositive(jsonSettings["millisecondsToWaitForEachGeneration"]);
}

inline int Settings::checkRGBValue(int value) const
{
    if (value >= 0 && value <= 255) { return value; }

    throw std::range_error("ERROR: RGB values must be 0 <= value <= 255");
}

inline int Settings::checkPositive(int value) const
{
    if (value >= 0) return value;

    throw std::range_error("ERROR: Passed a negative value, positive expected.");
}

#endif