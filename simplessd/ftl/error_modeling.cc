#include "ftl/error_modeling.hh"
#include "util/algorithm.hh"

#include <iostream>
#include <fstream>


namespace SimpleSSD {

namespace FTL {

ErrorModeling::ErrorModeling(){}

ErrorModeling::ErrorModeling(float temperature, float activationEnergy,
                             float epsilon, float alpha, float beta, float gamma,
                             float k, float m, float n, 
                             float sigma, uint32_t pageSize, uint32_t seed) {

  this->roomTemp = 25 + 273.15;
  this->temperature = temperature + 273.15;
  this->activationEnergy = activationEnergy;

  this->epsilon = epsilon;
  this->alpha = alpha;
  this->beta = beta;
  this->gamma = gamma;

  this->k = k;
  this->m = m;
  this->n = n;

  this->sigma = sigma;
  this->pageSize = pageSize;

  this-> generator = std::mt19937(seed);

  this->layerFactor = std::vector<float>(64);
  //std::cout << "epsilon " << epsilon << std::endl;
  //std::cout << "alpha " << alpha << std::endl;
  //std::cout << "beta " << beta << std::endl;
  //std::cout << "gamma " << gamma << std::endl;
  //std::cout << "k " << k << std::endl;
  //std::cout << "m " << m << std::endl;
  //std::cout << "n " << n << std::endl;
  
  // read File
  /*
  std::ifstream layerFile;
  layerFile.open("/home/rdolf/SimpleSSD/SimpleSSD-Standalone-base/simplessd/ftl/layer_factor.txt");
  std::string line;
  for(int i=0; i<64; i++){
    std::getline(layerFile, line);
    layerFactor[i] = std::stof(line);
    //std::cout << "layer factor " << i<< " " << layerFactor[i] << std::endl;
  }
  layerFile.close();
  */
  float factor = 1.00;
  for(int i=0; i<21; i++){
    layerFactor[i] = factor;
    factor = factor + 0.01;
  }
  factor = 1.00;
  for(int i=21; i<42; i++){
    layerFactor[i] = factor;
    factor = factor + 0.01;
  }
  factor = 1.00;
  for(int i=42; i<64; i++){
    layerFactor[i] = factor;
    factor = factor + 0.01;
  }
  
  
}

ErrorModeling::~ErrorModeling() {}

void ErrorModeling::setTemperature(float newTemp) {
  temperature = newTemp;
}
/*
float ErrorModeling::getAterm(float peCycle) {
  float result;

  if (peCycle < 1) {
      peCycle = 1;
  }

  //result = coeffA * log(peCycle) + constA;
  result = coeffA * peCycle + constA;

  return result;
}

float ErrorModeling::getBterm(float peCycle) {
  float result;

  if (peCycle < 1) {
      peCycle = 1;
  }

  result = coeffB * peCycle + constB;
  
  //result = coeffB * log(peCycle) + constB;
  //result = coeffB * peCycle + constB;


  return result;
}
*/

float ErrorModeling::getLayerFactor(uint32_t layer){
  return layerFactor[layer];
}

float ErrorModeling::arrhenius(float t2){
  float kb = 8.62 * pow(10, -5);
  float Ea = activationEnergy;

  float t1;

  t1 = t2 * exp((Ea/kb) * (1/roomTemp - 1/temperature));

  return t1;
}

/*
Error model :
  average RBER = epsilon + alpha * PE^k + beta * PE^m * retention time^n
  RBER = avg(RBER) * layer_factor
*/
/*
float ErrorModeling::getRBER(float retentionTime, float peCycle, uint32_t layer){
  
  double rber;

  retentionTime = arrhenius(retentionTime); // convert to time in room temp

  retentionTime = retentionTime / 1000000000000 / 60 / 60 / 24; // time unit : day

  //retentionTime = retentionTime * 23; // Acceleration (about 30 days)
  //retentionTime = retentionTime * 19; // Acceleration (about 25 days)

  //std::cout << "Converted retention " << retentionTime << std::endl;
  
  rber = epsilon + alpha * pow(peCycle, k) + beta * pow(peCycle, m) * pow(retentionTime, n);

  //std::cout << "rber " << rber << std::endl;
  
  rber = rber * layerFactor[layer];  

  //std::cout << "rber " << rber << std::endl;

  return rber;
}
*/
float ErrorModeling::getRBER(uint64_t retentionTime, float peCycle, uint32_t layer){ //Y.Luo 3D NAND version
  
  double rber;
  //std::cout << "retentionTime " << retentionTime << std::endl;

  retentionTime = arrhenius(retentionTime); // convert to time in room temp
  
  retentionTime = retentionTime / 1000000000000; // time unit : sec
  
  rber = (alpha * peCycle + beta) * log(retentionTime) + gamma * peCycle + epsilon;
  rber = exp(rber);

  //std::cout << "peCycle " << peCycle << std::endl;
  //std::cout << "retentionTime " << retentionTime << std::endl;
  //std::cout << "rber " << rber << std::endl;
  //std::cout << "layer " << layer << std::endl;
  //std::cout << "layer factor " << layerFactor[layer] << std::endl;
  
  rber = rber * layerFactor[layer];  

  std::normal_distribution<double> normal(rber, sigma);
  
  double randRber = normal(generator);

  //std::cout << "rand RBER " << randRber << std::endl<< std::endl;

  return randRber;
}


uint32_t ErrorModeling::getRandError(float retentionTime, float peCycle,
                                     uint32_t layer){
                                       
  float rber = getRBER(retentionTime, peCycle, layer);
  //std::cout << "rber " << rber << std::endl;
  //std::cout << "pecycle " << peCycle << std::endl;
  //std::cout << "retentionTime " << retentionTime << std::endl;
  

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

uint32_t ErrorModeling::getAverageError(float retentionTime, float peCycle,
                                      uint32_t layer){
                                       
  float rber = getRBER(retentionTime, peCycle, layer);
  //std::cout << "rber " << rber << std::endl;
  //std::cout << "pecycle " << peCycle << std::endl;
  //std::cout << "retentionTime " << retentionTime << std::endl;
  
  float averageError;

  averageError = rber * pageSize * 8;

  

  return static_cast<uint32_t>(averageError + 0.5);
}

}
}
