#include "ftl/error_modeling.hh"
#include "util/algorithm.hh"



namespace SimpleSSD {

namespace FTL {

ErrorModeling::ErrorModeling(){}

ErrorModeling::ErrorModeling(float temperature, float activationEnergy,
                             float coeffA, float coeffB,
                             float constA, float constB, 
                             float sigma, uint32_t seed) {

  this->roomTemp = 25 + 273.15;
  this->temperature = temperature + 273.15;
  this->activationEnergy = activationEnergy;

  this->coeffA = coeffA / 1000;
  this->coeffB = coeffB / 1000;
  this->constA = constA / 1000;
  this->constA = constB / 1000;

  this->sigma = sigma;

  this-> generator = std::mt19937(seed);
 
}

ErrorModeling::~ErrorModeling() {}

void ErrorModeling::setTemperature(float newTemp) {
  temperature = newTemp;
}

float ErrorModeling::getAterm(uint32_t peCycle) {
  float result;

  if (peCycle < 1) {
      peCycle = 1;
  }

  result = coeffA * log(peCycle) + constA;

  return result;
}

float ErrorModeling::getBterm(uint32_t peCycle) {
  float result;

  if (peCycle < 1) {
      peCycle = 1;
  }

  result = coeffB * log(peCycle) + constB;

  return result;
}

float ErrorModeling::arrhenius(float t2){
  float kb = 8.62 * pow(10, -5);
  float Ea = activationEnergy;

  float t1;

  t1 = t2 * exp((Ea/kb) * (1/roomTemp - 1/temperature));

  return t1;
}

float ErrorModeling::getRBER(float retentionTime, float peCycle){
  
  float rber;
  float aTerm = getAterm(peCycle);
  float bTerm = getBterm(peCycle);

  retentionTime = arrhenius(retentionTime); // convert to time in room temp

  rber = aTerm * log(retentionTime) + bTerm;

  return rber;
}

uint32_t ErrorModeling::getRandError(float retentionTime, float peCycle,
                                     uint32_t pageSize){
  float rber = getRBER(retentionTime, peCycle);

  

  float errorCount;
  float averageError;

  averageError = rber * pageSize * 8;

  std::normal_distribution<double> normal(averageError, sigma);

  errorCount = normal(generator);

  if (errorCount < 0){
    errorCount = 0;
  }
  

  return static_cast<uint32_t>(errorCount + 0.5);
}

}
}
