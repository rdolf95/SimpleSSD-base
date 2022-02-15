/*
 * Copyright (C) 2017 CAMELab
 *
 * This file is part of SimpleSSD.
 *
 * SimpleSSD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SimpleSSD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SimpleSSD.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __FTL_ERROR_MODELING__
#define __FTL_ERROR_MODELING__

#include <cinttypes>
#include <random>


namespace SimpleSSD {

namespace FTL {

class ErrorModeling {
 private:

  uint32_t pageSize;

  float temperature;
  float roomTemp;
  float activationEnergy;

  /*
  float coeffA; // coeff for A
  float coeffB; // coeff for B

  float constA; // const term for A
  float constB; // const term for B
  */

  float epsilon;
  float gamma;
  float alpha;
  float beta;
  float k;
  float m;
  float n;

  std::vector<float> layerFactor;
  
  float sigma;
  std::mt19937 generator;

  float arrhenius(float);   // convert time in high temp to time in room temp
  float getAterm(float);
  float getBterm(float);

  float getLayerFactor(uint32_t);



 public:
  ErrorModeling();
  ErrorModeling(float, float, float, float, float, float,
                float, float, float, float, uint32_t, uint32_t);
  ~ErrorModeling();
  void setTemperature(float);

  //float getRBER(float, float);
  float getRBER(uint64_t, float, uint32_t);
  uint32_t getRandError(float, float, uint32_t);
  uint32_t getAverageError(float, float, uint32_t);
};

}  // namespace FTL

}  // namespace SimpleSSD

#endif
