#include <iostream>
#include <vector>
#include "FastNoiseLite.h"

int main()
{
    // Create an instance of FastNoiseLite
    FastNoiseLite noise;

    // Set noise type
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    // Gather noise data
    std::vector<float> noiseData(128 * 128);
    int index = 0;

    for (int y = 0; y < 128; y++)
    {
        for (int x = 0; x < 128; x++)
        {
            noiseData[index++] = noise.GetNoise((float)x, (float)y);
        }
    }

    // Output noise data
    for (int i = 0; i < noiseData.size(); ++i)
    {
        std::cout << noiseData[i] << " ";
        if ((i + 1) % 128 == 0)
        {
            std::cout << std::endl;
        }
    }

    return 0;
}