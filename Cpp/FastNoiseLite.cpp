#include "FastNoiseLite.h"


FastNoiseLite::FastNoiseLite(int seed)
{
    mSeed = seed;
    mFrequency = 0.01f;
    mNoiseType = NoiseType_OpenSimplex2;
    mRotationType3D = RotationType3D_None;
    mTransformType3D = TransformType3D_DefaultOpenSimplex2;

    mFractalType = FractalType_None;
    mOctaves = 3;
    mLacunarity = 2.0f;
    mGain = 0.5f;
    mWeightedStrength = 0.0f;
    mPingPongStrength = 2.0f;

    mFractalBounding = 1 / 1.75f;

    mCellularDistanceFunction = CellularDistanceFunction_EuclideanSq;
    mCellularReturnType = CellularReturnType_Distance;
    mCellularJitterModifier = 1.0f;

    mDomainWarpType = DomainWarpType_OpenSimplex2;
    mWarpTransformType3D = TransformType3D_DefaultOpenSimplex2;
    mDomainWarpAmp = 1.0f;
}

void FastNoiseLite::SetSeed(int seed) { mSeed = seed; }

/// <summary>
/// Sets frequency for all noise types
/// </summary>
/// <remarks>
/// Default: 0.01
/// </remarks>

void FastNoiseLite::SetFrequency(float frequency) { mFrequency = frequency; }

/// <summary>
/// Sets noise algorithm used for GetNoise(...)
/// </summary>
/// <remarks>
/// Default: OpenSimplex2
/// </remarks>

void FastNoiseLite::SetNoiseType(NoiseType noiseType)
{
    mNoiseType = noiseType;
    UpdateTransformType3D();
}

/// <summary>
/// Sets domain rotation type for 3D Noise and 3D DomainWarp.
/// Can aid in reducing directional artifacts when sampling a 2D plane in 3D
/// </summary>
/// <remarks>
/// Default: None
/// </remarks>

void FastNoiseLite::SetRotationType3D(RotationType3D rotationType3D)
{
    mRotationType3D = rotationType3D;
    UpdateTransformType3D();
    UpdateWarpTransformType3D();
}

/// <summary>
/// Sets method for combining octaves in all fractal noise types
/// </summary>
/// <remarks>
/// Default: None
/// Note: FractalType_DomainWarp... only affects DomainWarp(...)
/// </remarks>

void FastNoiseLite::SetFractalType(FractalType fractalType) { mFractalType = fractalType; }

/// <summary>
/// Sets octave count for all fractal noise types
/// </summary>
/// <remarks>
/// Default: 3
/// </remarks>

void FastNoiseLite::SetFractalOctaves(int octaves)
{
    mOctaves = octaves;
    CalculateFractalBounding();
}

void FastNoiseLite::SetLacunarity(float lacunarity) { mLacunarity = lacunarity; }

/// <summary>
/// Sets octave lacunarity for all fractal noise types
/// </summary>
/// <remarks>
/// Default: 2.0
/// </remarks>

void FastNoiseLite::SetFractalLacunarity(float lacunarity) { mLacunarity = lacunarity; }

/// <summary>
/// Sets octave gain for all fractal noise types
/// </summary>
/// <remarks>
/// Default: 0.5
/// </remarks>

void FastNoiseLite::SetFractalGain(float gain)
{
    mGain = gain;
    CalculateFractalBounding();
}

/// <summary>
/// Sets octave weighting for all none DomainWarp fratal types
/// </summary>
/// <remarks>
/// Default: 0.0
/// Note: Keep between 0...1 to maintain -1...1 output bounding
/// </remarks>

void FastNoiseLite::SetFractalWeightedStrength(float weightedStrength) { mWeightedStrength = weightedStrength; }

/// <summary>
/// Sets strength of the fractal ping pong effect
/// </summary>
/// <remarks>
/// Default: 2.0
/// </remarks>

void FastNoiseLite::SetFractalPingPongStrength(float pingPongStrength) { mPingPongStrength = pingPongStrength; }

/// <summary>
/// Sets distance function used in cellular noise calculations
/// </summary>
/// <remarks>
/// Default: Distance
/// </remarks>

void FastNoiseLite::SetCellularDistanceFunction(CellularDistanceFunction cellularDistanceFunction) { mCellularDistanceFunction = cellularDistanceFunction; }

/// <summary>
/// Sets return type from cellular noise calculations
/// </summary>
/// <remarks>
/// Default: EuclideanSq
/// </remarks>

void FastNoiseLite::SetCellularReturnType(CellularReturnType cellularReturnType) { mCellularReturnType = cellularReturnType; }

/// <summary>
/// Sets the maximum distance a cellular point can move from it's grid position
/// </summary>
/// <remarks>
/// Default: 1.0
/// Note: Setting this higher than 1 will cause artifacts
/// </remarks>

void FastNoiseLite::SetCellularJitter(float cellularJitter) { mCellularJitterModifier = cellularJitter; }

/// <summary>
/// Sets the warp algorithm when using DomainWarp(...)
/// </summary>
/// <remarks>
/// Default: OpenSimplex2
/// </remarks>

void FastNoiseLite::SetDomainWarpType(DomainWarpType domainWarpType)
{
    mDomainWarpType = domainWarpType;
    UpdateWarpTransformType3D();
}

/// <summary>
/// Sets the maximum warp distance from original position when using DomainWarp(...)
/// </summary>
/// <remarks>
/// Default: 1.0
/// </remarks>

void FastNoiseLite::SetDomainWarpAmp(float domainWarpAmp) { mDomainWarpAmp = domainWarpAmp; }

float FastNoiseLite::FastMin(float a, float b) { return a < b ? a : b; }

float FastNoiseLite::FastMax(float a, float b) { return a > b ? a : b; }

float FastNoiseLite::FastAbs(float f) { return f < 0 ? -f : f; }

float FastNoiseLite::FastSqrt(float f) { return sqrtf(f); }

float FastNoiseLite::Lerp(float a, float b, float t) { return a + t * (b - a); }

float FastNoiseLite::InterpHermite(float t) { return t * t * (3 - 2 * t); }

float FastNoiseLite::InterpQuintic(float t) { return t * t * t * (t * (t * 6 - 15) + 10); }

float FastNoiseLite::CubicLerp(float a, float b, float c, float d, float t)
{
    float p = (d - c) - (a - b);
    return t * t * t * p + t * t * ((a - b) - p) + t * (c - a) + b;
}

float FastNoiseLite::PingPong(float t)
{
    t -= (int)(t * 0.5f) * 2;
    return t < 1 ? t : 2 - t;
}

void FastNoiseLite::CalculateFractalBounding()
{
    float gain = FastAbs(mGain);
    float amp = gain;
    float ampFractal = 1.0f;
    for (int i = 1; i < mOctaves; i++)
    {
        ampFractal += amp;
        amp *= gain;
    }
    mFractalBounding = 1 / ampFractal;
}

int FastNoiseLite::Hash(int seed, int xPrimed, int yPrimed)
{
    int hash = seed ^ xPrimed ^ yPrimed;

    hash *= 0x27d4eb2d;
    return hash;
}

int FastNoiseLite::Hash(int seed, int xPrimed, int yPrimed, int zPrimed)
{
    int hash = seed ^ xPrimed ^ yPrimed ^ zPrimed;

    hash *= 0x27d4eb2d;
    return hash;
}

float FastNoiseLite::ValCoord(int seed, int xPrimed, int yPrimed)
{
    int hash = Hash(seed, xPrimed, yPrimed);

    hash *= hash;
    hash ^= hash << 19;
    return hash * (1 / 2147483648.0f);
}

float FastNoiseLite::ValCoord(int seed, int xPrimed, int yPrimed, int zPrimed)
{
    int hash = Hash(seed, xPrimed, yPrimed, zPrimed);

    hash *= hash;
    hash ^= hash << 19;
    return hash * (1 / 2147483648.0f);
}

float FastNoiseLite::GradCoord(int seed, int xPrimed, int yPrimed, float xd, float yd) const
{
    int hash = Hash(seed, xPrimed, yPrimed);
    hash ^= hash >> 15;
    hash &= 127 << 1;

    float xg = Gradients2D[hash];
    float yg = Gradients2D[hash | 1];

    return xd * xg + yd * yg;
}

float FastNoiseLite::GradCoord(int seed, int xPrimed, int yPrimed, int zPrimed, float xd, float yd, float zd) const
{
    int hash = Hash(seed, xPrimed, yPrimed, zPrimed);
    hash ^= hash >> 15;
    hash &= 63 << 2;

    float xg = Gradients3D[hash];
    float yg = Gradients3D[hash | 1];
    float zg = Gradients3D[hash | 2];

    return xd * xg + yd * yg + zd * zg;
}

void FastNoiseLite::GradCoordOut(int seed, int xPrimed, int yPrimed, float& xo, float& yo) const
{
    int hash = Hash(seed, xPrimed, yPrimed) & (255 << 1);

    xo = RandVecs2D[hash];
    yo = RandVecs2D[hash | 1];
}

void FastNoiseLite::GradCoordOut(int seed, int xPrimed, int yPrimed, int zPrimed, float& xo, float& yo, float& zo) const
{
    int hash = Hash(seed, xPrimed, yPrimed, zPrimed) & (255 << 2);

    xo = RandVecs3D[hash];
    yo = RandVecs3D[hash | 1];
    zo = RandVecs3D[hash | 2];
}

void FastNoiseLite::GradCoordDual(int seed, int xPrimed, int yPrimed, float xd, float yd, float& xo, float& yo) const
{
    int hash = Hash(seed, xPrimed, yPrimed);
    int index1 = hash & (127 << 1);
    int index2 = (hash >> 7) & (255 << 1);

    float xg = Gradients2D[index1];
    float yg = Gradients2D[index1 | 1];
    float value = xd * xg + yd * yg;

    float xgo = RandVecs2D[index2];
    float ygo = RandVecs2D[index2 | 1];

    xo = value * xgo;
    yo = value * ygo;
}

void FastNoiseLite::GradCoordDual(int seed, int xPrimed, int yPrimed, int zPrimed, float xd, float yd, float zd, float& xo, float& yo, float& zo) const
{
    int hash = Hash(seed, xPrimed, yPrimed, zPrimed);
    int index1 = hash & (63 << 2);
    int index2 = (hash >> 6) & (255 << 2);

    float xg = Gradients3D[index1];
    float yg = Gradients3D[index1 | 1];
    float zg = Gradients3D[index1 | 2];
    float value = xd * xg + yd * yg + zd * zg;

    float xgo = RandVecs3D[index2];
    float ygo = RandVecs3D[index2 | 1];
    float zgo = RandVecs3D[index2 | 2];

    xo = value * xgo;
    yo = value * ygo;
    zo = value * zgo;
}

void FastNoiseLite::UpdateTransformType3D()
{
    switch (mRotationType3D)
    {
    case RotationType3D_ImproveXYPlanes:
        mTransformType3D = TransformType3D_ImproveXYPlanes;
        break;
    case RotationType3D_ImproveXZPlanes:
        mTransformType3D = TransformType3D_ImproveXZPlanes;
        break;
    default:
        switch (mNoiseType)
        {
        case NoiseType_OpenSimplex2:
        case NoiseType_OpenSimplex2S:
            mTransformType3D = TransformType3D_DefaultOpenSimplex2;
            break;
        default:
            mTransformType3D = TransformType3D_None;
            break;
        }
        break;
    }
}

void FastNoiseLite::UpdateWarpTransformType3D()
{
    switch (mRotationType3D)
    {
    case RotationType3D_ImproveXYPlanes:
        mWarpTransformType3D = TransformType3D_ImproveXYPlanes;
        break;
    case RotationType3D_ImproveXZPlanes:
        mWarpTransformType3D = TransformType3D_ImproveXZPlanes;
        break;
    default:
        switch (mDomainWarpType)
        {
        case DomainWarpType_OpenSimplex2:
        case DomainWarpType_OpenSimplex2Reduced:
            mWarpTransformType3D = TransformType3D_DefaultOpenSimplex2;
            break;
        default:
            mWarpTransformType3D = TransformType3D_None;
            break;
        }
        break;
    }
}

void FastNoiseLite::SingleDomainWarpBasicGrid(int seed, float warpAmp, float frequency, float x, float y, float& xr, float& yr) const
{
    float xf = x * frequency;
    float yf = y * frequency;

    int x0 = FastFloor(xf);
    int y0 = FastFloor(yf);

    float xs = InterpHermite((float)(xf - x0));
    float ys = InterpHermite((float)(yf - y0));

    x0 *= PrimeX;
    y0 *= PrimeY;
    int x1 = x0 + PrimeX;
    int y1 = y0 + PrimeY;

    int hash0 = Hash(seed, x0, y0) & (255 << 1);
    int hash1 = Hash(seed, x1, y0) & (255 << 1);

    float lx0x = Lerp(RandVecs2D[hash0], RandVecs2D[hash1], xs);
    float ly0x = Lerp(RandVecs2D[hash0 | 1], RandVecs2D[hash1 | 1], xs);

    hash0 = Hash(seed, x0, y1) & (255 << 1);
    hash1 = Hash(seed, x1, y1) & (255 << 1);

    float lx1x = Lerp(RandVecs2D[hash0], RandVecs2D[hash1], xs);
    float ly1x = Lerp(RandVecs2D[hash0 | 1], RandVecs2D[hash1 | 1], xs);

    xr += Lerp(lx0x, lx1x, ys) * warpAmp;
    yr += Lerp(ly0x, ly1x, ys) * warpAmp;
}


/// <summary>
/// 2D noise at given position using current settings
/// </summary>
/// <returns>
/// Noise output bounded between -1...1
/// </returns>


float FastNoiseLite::GetNoise(float x, float y) const
{
    TransformNoiseCoordinate(x, y);

    switch (mFractalType)
    {
    default:
        return GenNoiseSingle(mSeed, x, y);
    case FractalType_FBm:
        return GenFractalFBm(x, y);
    case FractalType_Ridged:
        return GenFractalRidged(x, y);
    case FractalType_PingPong:
        return GenFractalPingPong(x, y);
    }
}

/// <summary>
/// 3D noise at given position using current settings
/// </summary>
/// <returns>
/// Noise output bounded between -1...1
/// </returns>


float FastNoiseLite::GetNoise(float x, float y, float z) const
{
    TransformNoiseCoordinate(x, y, z);

    switch (mFractalType)
    {
    default:
        return GenNoiseSingle(mSeed, x, y, z);
    case FractalType_FBm:
        return GenFractalFBm(x, y, z);
    case FractalType_Ridged:
        return GenFractalRidged(x, y, z);
    case FractalType_PingPong:
        return GenFractalPingPong(x, y, z);
    }
}

/// <summary>
/// 2D warps the input position using current domain warp settings
/// </summary>
/// <example>
/// Example usage with GetNoise
/// <code>DomainWarp(x, y)
/// noise = GetNoise(x, y)</code>
/// </example>


void FastNoiseLite::DomainWarp(float& x, float& y) const
{
    switch (mFractalType)
    {
    default:
        DomainWarpSingle(x, y);
        break;
    case FractalType_DomainWarpProgressive:
        DomainWarpFractalProgressive(x, y);
        break;
    case FractalType_DomainWarpIndependent:
        DomainWarpFractalIndependent(x, y);
        break;
    }
}

/// <summary>
/// 3D warps the input position using current domain warp settings
/// </summary>
/// <example>
/// Example usage with GetNoise
/// <code>DomainWarp(x, y, z)
/// noise = GetNoise(x, y, z)</code>
/// </example>


void FastNoiseLite::DomainWarp(float& x, float& y, float& z) const
{

    switch (mFractalType)
    {
    default:
        DomainWarpSingle(x, y, z);
        break;
    case FractalType_DomainWarpProgressive:
        DomainWarpFractalProgressive(x, y, z);
        break;
    case FractalType_DomainWarpIndependent:
        DomainWarpFractalIndependent(x, y, z);
        break;
    }
}


int FastNoiseLite::FastFloor(float f) { return f >= 0 ? (int)f : (int)f - 1; }


int FastNoiseLite::FastRound(float f) { return f >= 0 ? (int)(f + 0.5f) : (int)(f - 0.5f); }


float FastNoiseLite::GenNoiseSingle(int seed, float x, float y) const
{
    switch (mNoiseType)
    {
    case NoiseType_OpenSimplex2:
        return SingleSimplex(seed, x, y);
    case NoiseType_OpenSimplex2S:
        return SingleOpenSimplex2S(seed, x, y);
    case NoiseType_Cellular:
        return SingleCellular(seed, x, y);
    case NoiseType_Perlin:
        return SinglePerlin(seed, x, y);
    case NoiseType_ValueCubic:
        return SingleValueCubic(seed, x, y);
    case NoiseType_Value:
        return SingleValue(seed, x, y);
    default:
        return 0;
    }
}


float FastNoiseLite::GenNoiseSingle(int seed, float x, float y, float z) const
{
    switch (mNoiseType)
    {
    case NoiseType_OpenSimplex2:
        return SingleOpenSimplex2(seed, x, y, z);
    case NoiseType_OpenSimplex2S:
        return SingleOpenSimplex2S(seed, x, y, z);
    case NoiseType_Cellular:
        return SingleCellular(seed, x, y, z);
    case NoiseType_Perlin:
        return SinglePerlin(seed, x, y, z);
    case NoiseType_ValueCubic:
        return SingleValueCubic(seed, x, y, z);
    case NoiseType_Value:
        return SingleValue(seed, x, y, z);
    default:
        return 0;
    }
}


void FastNoiseLite::TransformNoiseCoordinate(float& x, float& y) const
{
    x *= mFrequency;
    y *= mFrequency;

    switch (mNoiseType)
    {
    case NoiseType_OpenSimplex2:
    case NoiseType_OpenSimplex2S:
    {
        const float SQRT3 = (float)1.7320508075688772935274463415059;
        const float F2 = 0.5f * (SQRT3 - 1);
        float t = (x + y) * F2;
        x += t;
        y += t;
    }
    break;
    default:
        break;
    }
}


void FastNoiseLite::TransformNoiseCoordinate(float& x, float& y, float& z) const
{
    x *= mFrequency;
    y *= mFrequency;
    z *= mFrequency;

    switch (mTransformType3D)
    {
    case TransformType3D_ImproveXYPlanes:
    {
        float xy = x + y;
        float s2 = xy * -(float)0.211324865405187;
        z *= (float)0.577350269189626;
        x += s2 - z;
        y = y + s2 - z;
        z += xy * (float)0.577350269189626;
    }
    break;
    case TransformType3D_ImproveXZPlanes:
    {
        float xz = x + z;
        float s2 = xz * -(float)0.211324865405187;
        y *= (float)0.577350269189626;
        x += s2 - y;
        z += s2 - y;
        y += xz * (float)0.577350269189626;
    }
    break;
    case TransformType3D_DefaultOpenSimplex2:
    {
        const float R3 = (float)(2.0 / 3.0);
        float r = (x + y + z) * R3; // Rotation, not skew
        x = r - x;
        y = r - y;
        z = r - z;
    }
    break;
    default:
        break;
    }
}


void FastNoiseLite::TransformDomainWarpCoordinate(float& x, float& y) const
{
    switch (mDomainWarpType)
    {
    case DomainWarpType_OpenSimplex2:
    case DomainWarpType_OpenSimplex2Reduced:
    {
        const float SQRT3 = (float)1.7320508075688772935274463415059;
        const float F2 = 0.5f * (SQRT3 - 1);
        float t = (x + y) * F2;
        x += t;
        y += t;
    }
    break;
    default:
        break;
    }
}


void FastNoiseLite::TransformDomainWarpCoordinate(float& x, float& y, float& z) const
{
    switch (mWarpTransformType3D)
    {
    case TransformType3D_ImproveXYPlanes:
    {
        float xy = x + y;
        float s2 = xy * -(float)0.211324865405187;
        z *= (float)0.577350269189626;
        x += s2 - z;
        y = y + s2 - z;
        z += xy * (float)0.577350269189626;
    }
    break;
    case TransformType3D_ImproveXZPlanes:
    {
        float xz = x + z;
        float s2 = xz * -(float)0.211324865405187;
        y *= (float)0.577350269189626;
        x += s2 - y;
        z += s2 - y;
        y += xz * (float)0.577350269189626;
    }
    break;
    case TransformType3D_DefaultOpenSimplex2:
    {
        const float R3 = (float)(2.0 / 3.0);
        float r = (x + y + z) * R3; // Rotation, not skew
        x = r - x;
        y = r - y;
        z = r - z;
    }
    break;
    default:
        break;
    }
}


float FastNoiseLite::GenFractalFBm(float x, float y) const
{
    int seed = mSeed;
    float sum = 0;
    float amp = mFractalBounding;

    for (int i = 0; i < mOctaves; i++)
    {
        float noise = GenNoiseSingle(seed++, x, y);
        sum += noise * amp;
        amp *= Lerp(1.0f, FastMin(noise + 1, 2) * 0.5f, mWeightedStrength);

        x *= mLacunarity;
        y *= mLacunarity;
        amp *= mGain;
    }

    return sum;
}


float FastNoiseLite::GenFractalFBm(float x, float y, float z) const
{
    int seed = mSeed;
    float sum = 0;
    float amp = mFractalBounding;

    for (int i = 0; i < mOctaves; i++)
    {
        float noise = GenNoiseSingle(seed++, x, y, z);
        sum += noise * amp;
        amp *= Lerp(1.0f, (noise + 1) * 0.5f, mWeightedStrength);

        x *= mLacunarity;
        y *= mLacunarity;
        z *= mLacunarity;
        amp *= mGain;
    }

    return sum;
}


float FastNoiseLite::GenFractalRidged(float x, float y) const
{
    int seed = mSeed;
    float sum = 0;
    float amp = mFractalBounding;

    for (int i = 0; i < mOctaves; i++)
    {
        float noise = FastAbs(GenNoiseSingle(seed++, x, y));
        sum += (noise * -2 + 1) * amp;
        amp *= Lerp(1.0f, 1 - noise, mWeightedStrength);

        x *= mLacunarity;
        y *= mLacunarity;
        amp *= mGain;
    }

    return sum;
}


float FastNoiseLite::GenFractalRidged(float x, float y, float z) const
{
    int seed = mSeed;
    float sum = 0;
    float amp = mFractalBounding;

    for (int i = 0; i < mOctaves; i++)
    {
        float noise = FastAbs(GenNoiseSingle(seed++, x, y, z));
        sum += (noise * -2 + 1) * amp;
        amp *= Lerp(1.0f, 1 - noise, mWeightedStrength);

        x *= mLacunarity;
        y *= mLacunarity;
        z *= mLacunarity;
        amp *= mGain;
    }

    return sum;
}


float FastNoiseLite::GenFractalPingPong(float x, float y) const
{
    int seed = mSeed;
    float sum = 0;
    float amp = mFractalBounding;

    for (int i = 0; i < mOctaves; i++)
    {
        float noise = PingPong((GenNoiseSingle(seed++, x, y) + 1) * mPingPongStrength);
        sum += (noise - 0.5f) * 2 * amp;
        amp *= Lerp(1.0f, noise, mWeightedStrength);

        x *= mLacunarity;
        y *= mLacunarity;
        amp *= mGain;
    }

    return sum;
}


float FastNoiseLite::GenFractalPingPong(float x, float y, float z) const
{
    int seed = mSeed;
    float sum = 0;
    float amp = mFractalBounding;

    for (int i = 0; i < mOctaves; i++)
    {
        float noise = PingPong((GenNoiseSingle(seed++, x, y, z) + 1) * mPingPongStrength);
        sum += (noise - 0.5f) * 2 * amp;
        amp *= Lerp(1.0f, noise, mWeightedStrength);

        x *= mLacunarity;
        y *= mLacunarity;
        z *= mLacunarity;
        amp *= mGain;
    }

    return sum;
}


float FastNoiseLite::SingleSimplex(int seed, float x, float y) const
{
    // 2D OpenSimplex2 case uses the same algorithm as ordinary Simplex.

    const float SQRT3 = 1.7320508075688772935274463415059f;
    const float G2 = (3 - SQRT3) / 6;

    /*
     * --- Skew moved to TransformNoiseCoordinate method ---
     * const float F2 = 0.5f * (SQRT3 - 1);
     * float s = (x + y) * F2;
     * x += s; y += s;
     */

    int i = FastFloor(x);
    int j = FastFloor(y);
    float xi = (float)(x - i);
    float yi = (float)(y - j);

    float t = (xi + yi) * G2;
    float x0 = (float)(xi - t);
    float y0 = (float)(yi - t);

    i *= PrimeX;
    j *= PrimeY;

    float n0, n1, n2;

    float a = 0.5f - x0 * x0 - y0 * y0;
    if (a <= 0)
        n0 = 0;
    else
    {
        n0 = (a * a) * (a * a) * GradCoord(seed, i, j, x0, y0);
    }

    float c = (float)(2 * (1 - 2 * G2) * (1 / G2 - 2)) * t + ((float)(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a);
    if (c <= 0)
        n2 = 0;
    else
    {
        float x2 = x0 + (2 * (float)G2 - 1);
        float y2 = y0 + (2 * (float)G2 - 1);
        n2 = (c * c) * (c * c) * GradCoord(seed, i + PrimeX, j + PrimeY, x2, y2);
    }

    if (y0 > x0)
    {
        float x1 = x0 + (float)G2;
        float y1 = y0 + ((float)G2 - 1);
        float b = 0.5f - x1 * x1 - y1 * y1;
        if (b <= 0)
            n1 = 0;
        else
        {
            n1 = (b * b) * (b * b) * GradCoord(seed, i, j + PrimeY, x1, y1);
        }
    }
    else
    {
        float x1 = x0 + ((float)G2 - 1);
        float y1 = y0 + (float)G2;
        float b = 0.5f - x1 * x1 - y1 * y1;
        if (b <= 0)
            n1 = 0;
        else
        {
            n1 = (b * b) * (b * b) * GradCoord(seed, i + PrimeX, j, x1, y1);
        }
    }

    return (n0 + n1 + n2) * 99.83685446303647f;
}


float FastNoiseLite::SingleOpenSimplex2(int seed, float x, float y, float z) const
{
    // 3D OpenSimplex2 case uses two offset rotated cube grids.

    /*
     * --- Rotation moved to TransformNoiseCoordinate method ---
     * const float R3 = (float)(2.0 / 3.0);
     * float r = (x + y + z) * R3; // Rotation, not skew
     * x = r - x; y = r - y; z = r - z;
     */

    int i = FastRound(x);
    int j = FastRound(y);
    int k = FastRound(z);
    float x0 = (float)(x - i);
    float y0 = (float)(y - j);
    float z0 = (float)(z - k);

    int xNSign = (int)(-1.0f - x0) | 1;
    int yNSign = (int)(-1.0f - y0) | 1;
    int zNSign = (int)(-1.0f - z0) | 1;

    float ax0 = xNSign * -x0;
    float ay0 = yNSign * -y0;
    float az0 = zNSign * -z0;

    i *= PrimeX;
    j *= PrimeY;
    k *= PrimeZ;

    float value = 0;
    float a = (0.6f - x0 * x0) - (y0 * y0 + z0 * z0);

    for (int l = 0;; l++)
    {
        if (a > 0)
        {
            value += (a * a) * (a * a) * GradCoord(seed, i, j, k, x0, y0, z0);
        }

        float b = a + 1;
        int i1 = i;
        int j1 = j;
        int k1 = k;
        float x1 = x0;
        float y1 = y0;
        float z1 = z0;

        if (ax0 >= ay0 && ax0 >= az0)
        {
            x1 += xNSign;
            b -= xNSign * 2 * x1;
            i1 -= xNSign * PrimeX;
        }
        else if (ay0 > ax0 && ay0 >= az0)
        {
            y1 += yNSign;
            b -= yNSign * 2 * y1;
            j1 -= yNSign * PrimeY;
        }
        else
        {
            z1 += zNSign;
            b -= zNSign * 2 * z1;
            k1 -= zNSign * PrimeZ;
        }

        if (b > 0)
        {
            value += (b * b) * (b * b) * GradCoord(seed, i1, j1, k1, x1, y1, z1);
        }

        if (l == 1)
            break;

        ax0 = 0.5f - ax0;
        ay0 = 0.5f - ay0;
        az0 = 0.5f - az0;

        x0 = xNSign * ax0;
        y0 = yNSign * ay0;
        z0 = zNSign * az0;

        a += (0.75f - ax0) - (ay0 + az0);

        i += (xNSign >> 1) & PrimeX;
        j += (yNSign >> 1) & PrimeY;
        k += (zNSign >> 1) & PrimeZ;

        xNSign = -xNSign;
        yNSign = -yNSign;
        zNSign = -zNSign;

        seed = ~seed;
    }

    return value * 32.69428253173828125f;
}


float FastNoiseLite::SingleOpenSimplex2S(int seed, float x, float y) const
{
    // 2D OpenSimplex2S case is a modified 2D simplex noise.

    const float SQRT3 = (float)1.7320508075688772935274463415059;
    const float G2 = (3 - SQRT3) / 6;

    /*
     * --- Skew moved to TransformNoiseCoordinate method ---
     * const float F2 = 0.5f * (SQRT3 - 1);
     * float s = (x + y) * F2;
     * x += s; y += s;
     */

    int i = FastFloor(x);
    int j = FastFloor(y);
    float xi = (float)(x - i);
    float yi = (float)(y - j);

    i *= PrimeX;
    j *= PrimeY;
    int i1 = i + PrimeX;
    int j1 = j + PrimeY;

    float t = (xi + yi) * (float)G2;
    float x0 = xi - t;
    float y0 = yi - t;

    float a0 = (2.0f / 3.0f) - x0 * x0 - y0 * y0;
    float value = (a0 * a0) * (a0 * a0) * GradCoord(seed, i, j, x0, y0);

    float a1 = (float)(2 * (1 - 2 * G2) * (1 / G2 - 2)) * t + ((float)(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a0);
    float x1 = x0 - (float)(1 - 2 * G2);
    float y1 = y0 - (float)(1 - 2 * G2);
    value += (a1 * a1) * (a1 * a1) * GradCoord(seed, i1, j1, x1, y1);

    // Nested conditionals were faster than compact bit logic/arithmetic.
    float xmyi = xi - yi;
    if (t > G2)
    {
        if (xi + xmyi > 1)
        {
            float x2 = x0 + (float)(3 * G2 - 2);
            float y2 = y0 + (float)(3 * G2 - 1);
            float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i + (PrimeX << 1), j + PrimeY, x2, y2);
            }
        }
        else
        {
            float x2 = x0 + (float)G2;
            float y2 = y0 + (float)(G2 - 1);
            float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i, j + PrimeY, x2, y2);
            }
        }

        if (yi - xmyi > 1)
        {
            float x3 = x0 + (float)(3 * G2 - 1);
            float y3 = y0 + (float)(3 * G2 - 2);
            float a3 = (2.0f / 3.0f) - x3 * x3 - y3 * y3;
            if (a3 > 0)
            {
                value += (a3 * a3) * (a3 * a3) * GradCoord(seed, i + PrimeX, j + (PrimeY << 1), x3, y3);
            }
        }
        else
        {
            float x3 = x0 + (float)(G2 - 1);
            float y3 = y0 + (float)G2;
            float a3 = (2.0f / 3.0f) - x3 * x3 - y3 * y3;
            if (a3 > 0)
            {
                value += (a3 * a3) * (a3 * a3) * GradCoord(seed, i + PrimeX, j, x3, y3);
            }
        }
    }
    else
    {
        if (xi + xmyi < 0)
        {
            float x2 = x0 + (float)(1 - G2);
            float y2 = y0 - (float)G2;
            float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i - PrimeX, j, x2, y2);
            }
        }
        else
        {
            float x2 = x0 + (float)(G2 - 1);
            float y2 = y0 + (float)G2;
            float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i + PrimeX, j, x2, y2);
            }
        }

        if (yi < xmyi)
        {
            float x2 = x0 - (float)G2;
            float y2 = y0 - (float)(G2 - 1);
            float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i, j - PrimeY, x2, y2);
            }
        }
        else
        {
            float x2 = x0 + (float)G2;
            float y2 = y0 + (float)(G2 - 1);
            float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
            if (a2 > 0)
            {
                value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i, j + PrimeY, x2, y2);
            }
        }
    }

    return value * 18.24196194486065f;
}


float FastNoiseLite::SingleOpenSimplex2S(int seed, float x, float y, float z) const
{
    // 3D OpenSimplex2S case uses two offset rotated cube grids.

    /*
     * --- Rotation moved to TransformNoiseCoordinate method ---
     * const float R3 = (float)(2.0 / 3.0);
     * float r = (x + y + z) * R3; // Rotation, not skew
     * x = r - x; y = r - y; z = r - z;
     */

    int i = FastFloor(x);
    int j = FastFloor(y);
    int k = FastFloor(z);
    float xi = (float)(x - i);
    float yi = (float)(y - j);
    float zi = (float)(z - k);

    i *= PrimeX;
    j *= PrimeY;
    k *= PrimeZ;
    int seed2 = seed + 1293373;

    int xNMask = (int)(-0.5f - xi);
    int yNMask = (int)(-0.5f - yi);
    int zNMask = (int)(-0.5f - zi);

    float x0 = xi + xNMask;
    float y0 = yi + yNMask;
    float z0 = zi + zNMask;
    float a0 = 0.75f - x0 * x0 - y0 * y0 - z0 * z0;
    float value = (a0 * a0) * (a0 * a0) * GradCoord(seed, i + (xNMask & PrimeX), j + (yNMask & PrimeY), k + (zNMask & PrimeZ), x0, y0, z0);

    float x1 = xi - 0.5f;
    float y1 = yi - 0.5f;
    float z1 = zi - 0.5f;
    float a1 = 0.75f - x1 * x1 - y1 * y1 - z1 * z1;
    value += (a1 * a1) * (a1 * a1) * GradCoord(seed2, i + PrimeX, j + PrimeY, k + PrimeZ, x1, y1, z1);

    float xAFlipMask0 = ((xNMask | 1) << 1) * x1;
    float yAFlipMask0 = ((yNMask | 1) << 1) * y1;
    float zAFlipMask0 = ((zNMask | 1) << 1) * z1;
    float xAFlipMask1 = (-2 - (xNMask << 2)) * x1 - 1.0f;
    float yAFlipMask1 = (-2 - (yNMask << 2)) * y1 - 1.0f;
    float zAFlipMask1 = (-2 - (zNMask << 2)) * z1 - 1.0f;

    bool skip5 = false;
    float a2 = xAFlipMask0 + a0;
    if (a2 > 0)
    {
        float x2 = x0 - (xNMask | 1);
        float y2 = y0;
        float z2 = z0;
        value += (a2 * a2) * (a2 * a2) * GradCoord(seed, i + (~xNMask & PrimeX), j + (yNMask & PrimeY), k + (zNMask & PrimeZ), x2, y2, z2);
    }
    else
    {
        float a3 = yAFlipMask0 + zAFlipMask0 + a0;
        if (a3 > 0)
        {
            float x3 = x0;
            float y3 = y0 - (yNMask | 1);
            float z3 = z0 - (zNMask | 1);
            value += (a3 * a3) * (a3 * a3) * GradCoord(seed, i + (xNMask & PrimeX), j + (~yNMask & PrimeY), k + (~zNMask & PrimeZ), x3, y3, z3);
        }

        float a4 = xAFlipMask1 + a1;
        if (a4 > 0)
        {
            float x4 = (xNMask | 1) + x1;
            float y4 = y1;
            float z4 = z1;
            value += (a4 * a4) * (a4 * a4) * GradCoord(seed2, i + (xNMask & (PrimeX * 2)), j + PrimeY, k + PrimeZ, x4, y4, z4);
            skip5 = true;
        }
    }

    bool skip9 = false;
    float a6 = yAFlipMask0 + a0;
    if (a6 > 0)
    {
        float x6 = x0;
        float y6 = y0 - (yNMask | 1);
        float z6 = z0;
        value += (a6 * a6) * (a6 * a6) * GradCoord(seed, i + (xNMask & PrimeX), j + (~yNMask & PrimeY), k + (zNMask & PrimeZ), x6, y6, z6);
    }
    else
    {
        float a7 = xAFlipMask0 + zAFlipMask0 + a0;
        if (a7 > 0)
        {
            float x7 = x0 - (xNMask | 1);
            float y7 = y0;
            float z7 = z0 - (zNMask | 1);
            value += (a7 * a7) * (a7 * a7) * GradCoord(seed, i + (~xNMask & PrimeX), j + (yNMask & PrimeY), k + (~zNMask & PrimeZ), x7, y7, z7);
        }

        float a8 = yAFlipMask1 + a1;
        if (a8 > 0)
        {
            float x8 = x1;
            float y8 = (yNMask | 1) + y1;
            float z8 = z1;
            value += (a8 * a8) * (a8 * a8) * GradCoord(seed2, i + PrimeX, j + (yNMask & (PrimeY << 1)), k + PrimeZ, x8, y8, z8);
            skip9 = true;
        }
    }

    bool skipD = false;
    float aA = zAFlipMask0 + a0;
    if (aA > 0)
    {
        float xA = x0;
        float yA = y0;
        float zA = z0 - (zNMask | 1);
        value += (aA * aA) * (aA * aA) * GradCoord(seed, i + (xNMask & PrimeX), j + (yNMask & PrimeY), k + (~zNMask & PrimeZ), xA, yA, zA);
    }
    else
    {
        float aB = xAFlipMask0 + yAFlipMask0 + a0;
        if (aB > 0)
        {
            float xB = x0 - (xNMask | 1);
            float yB = y0 - (yNMask | 1);
            float zB = z0;
            value += (aB * aB) * (aB * aB) * GradCoord(seed, i + (~xNMask & PrimeX), j + (~yNMask & PrimeY), k + (zNMask & PrimeZ), xB, yB, zB);
        }

        float aC = zAFlipMask1 + a1;
        if (aC > 0)
        {
            float xC = x1;
            float yC = y1;
            float zC = (zNMask | 1) + z1;
            value += (aC * aC) * (aC * aC) * GradCoord(seed2, i + PrimeX, j + PrimeY, k + (zNMask & (PrimeZ << 1)), xC, yC, zC);
            skipD = true;
        }
    }

    if (!skip5)
    {
        float a5 = yAFlipMask1 + zAFlipMask1 + a1;
        if (a5 > 0)
        {
            float x5 = x1;
            float y5 = (yNMask | 1) + y1;
            float z5 = (zNMask | 1) + z1;
            value += (a5 * a5) * (a5 * a5) * GradCoord(seed2, i + PrimeX, j + (yNMask & (PrimeY << 1)), k + (zNMask & (PrimeZ << 1)), x5, y5, z5);
        }
    }

    if (!skip9)
    {
        float a9 = xAFlipMask1 + zAFlipMask1 + a1;
        if (a9 > 0)
        {
            float x9 = (xNMask | 1) + x1;
            float y9 = y1;
            float z9 = (zNMask | 1) + z1;
            value += (a9 * a9) * (a9 * a9) * GradCoord(seed2, i + (xNMask & (PrimeX * 2)), j + PrimeY, k + (zNMask & (PrimeZ << 1)), x9, y9, z9);
        }
    }

    if (!skipD)
    {
        float aD = xAFlipMask1 + yAFlipMask1 + a1;
        if (aD > 0)
        {
            float xD = (xNMask | 1) + x1;
            float yD = (yNMask | 1) + y1;
            float zD = z1;
            value += (aD * aD) * (aD * aD) * GradCoord(seed2, i + (xNMask & (PrimeX << 1)), j + (yNMask & (PrimeY << 1)), k + PrimeZ, xD, yD, zD);
        }
    }

    return value * 9.046026385208288f;
}


float FastNoiseLite::SingleCellular(int seed, float x, float y) const
{
    int xr = FastRound(x);
    int yr = FastRound(y);

    float distance0 = 1e10f;
    float distance1 = 1e10f;
    int closestHash = 0;

    float cellularJitter = 0.43701595f * mCellularJitterModifier;

    int xPrimed = (xr - 1) * PrimeX;
    int yPrimedBase = (yr - 1) * PrimeY;

    switch (mCellularDistanceFunction)
    {
    default:
    case CellularDistanceFunction_Euclidean:
    case CellularDistanceFunction_EuclideanSq:
        for (int xi = xr - 1; xi <= xr + 1; xi++)
        {
            int yPrimed = yPrimedBase;

            for (int yi = yr - 1; yi <= yr + 1; yi++)
            {
                int hash = Hash(seed, xPrimed, yPrimed);
                int idx = hash & (255 << 1);

                float vecX = (float)(xi - x) + RandVecs2D[idx] * cellularJitter;
                float vecY = (float)(yi - y) + RandVecs2D[idx | 1] * cellularJitter;

                float newDistance = vecX * vecX + vecY * vecY;

                distance1 = FastMax(FastMin(distance1, newDistance), distance0);
                if (newDistance < distance0)
                {
                    distance0 = newDistance;
                    closestHash = hash;
                }
                yPrimed += PrimeY;
            }
            xPrimed += PrimeX;
        }
        break;
    case CellularDistanceFunction_Manhattan:
        for (int xi = xr - 1; xi <= xr + 1; xi++)
        {
            int yPrimed = yPrimedBase;

            for (int yi = yr - 1; yi <= yr + 1; yi++)
            {
                int hash = Hash(seed, xPrimed, yPrimed);
                int idx = hash & (255 << 1);

                float vecX = (float)(xi - x) + RandVecs2D[idx] * cellularJitter;
                float vecY = (float)(yi - y) + RandVecs2D[idx | 1] * cellularJitter;

                float newDistance = FastAbs(vecX) + FastAbs(vecY);

                distance1 = FastMax(FastMin(distance1, newDistance), distance0);
                if (newDistance < distance0)
                {
                    distance0 = newDistance;
                    closestHash = hash;
                }
                yPrimed += PrimeY;
            }
            xPrimed += PrimeX;
        }
        break;
    case CellularDistanceFunction_Hybrid:
        for (int xi = xr - 1; xi <= xr + 1; xi++)
        {
            int yPrimed = yPrimedBase;

            for (int yi = yr - 1; yi <= yr + 1; yi++)
            {
                int hash = Hash(seed, xPrimed, yPrimed);
                int idx = hash & (255 << 1);

                float vecX = (float)(xi - x) + RandVecs2D[idx] * cellularJitter;
                float vecY = (float)(yi - y) + RandVecs2D[idx | 1] * cellularJitter;

                float newDistance = (FastAbs(vecX) + FastAbs(vecY)) + (vecX * vecX + vecY * vecY);

                distance1 = FastMax(FastMin(distance1, newDistance), distance0);
                if (newDistance < distance0)
                {
                    distance0 = newDistance;
                    closestHash = hash;
                }
                yPrimed += PrimeY;
            }
            xPrimed += PrimeX;
        }
        break;
    }

    if (mCellularDistanceFunction == CellularDistanceFunction_Euclidean && mCellularReturnType >= CellularReturnType_Distance)
    {
        distance0 = FastSqrt(distance0);

        if (mCellularReturnType >= CellularReturnType_Distance2)
        {
            distance1 = FastSqrt(distance1);
        }
    }

    switch (mCellularReturnType)
    {
    case CellularReturnType_CellValue:
        return closestHash * (1 / 2147483648.0f);
    case CellularReturnType_Distance:
        return distance0 - 1;
    case CellularReturnType_Distance2:
        return distance1 - 1;
    case CellularReturnType_Distance2Add:
        return (distance1 + distance0) * 0.5f - 1;
    case CellularReturnType_Distance2Sub:
        return distance1 - distance0 - 1;
    case CellularReturnType_Distance2Mul:
        return distance1 * distance0 * 0.5f - 1;
    case CellularReturnType_Distance2Div:
        return distance0 / distance1 - 1;
    default:
        return 0;
    }
}


float FastNoiseLite::SingleCellular(int seed, float x, float y, float z) const
{
    int xr = FastRound(x);
    int yr = FastRound(y);
    int zr = FastRound(z);

    float distance0 = 1e10f;
    float distance1 = 1e10f;
    int closestHash = 0;

    float cellularJitter = 0.39614353f * mCellularJitterModifier;

    int xPrimed = (xr - 1) * PrimeX;
    int yPrimedBase = (yr - 1) * PrimeY;
    int zPrimedBase = (zr - 1) * PrimeZ;

    switch (mCellularDistanceFunction)
    {
    case CellularDistanceFunction_Euclidean:
    case CellularDistanceFunction_EuclideanSq:
        for (int xi = xr - 1; xi <= xr + 1; xi++)
        {
            int yPrimed = yPrimedBase;

            for (int yi = yr - 1; yi <= yr + 1; yi++)
            {
                int zPrimed = zPrimedBase;

                for (int zi = zr - 1; zi <= zr + 1; zi++)
                {
                    int hash = Hash(seed, xPrimed, yPrimed, zPrimed);
                    int idx = hash & (255 << 2);

                    float vecX = (float)(xi - x) + RandVecs3D[idx] * cellularJitter;
                    float vecY = (float)(yi - y) + RandVecs3D[idx | 1] * cellularJitter;
                    float vecZ = (float)(zi - z) + RandVecs3D[idx | 2] * cellularJitter;

                    float newDistance = vecX * vecX + vecY * vecY + vecZ * vecZ;

                    distance1 = FastMax(FastMin(distance1, newDistance), distance0);
                    if (newDistance < distance0)
                    {
                        distance0 = newDistance;
                        closestHash = hash;
                    }
                    zPrimed += PrimeZ;
                }
                yPrimed += PrimeY;
            }
            xPrimed += PrimeX;
        }
        break;
    case CellularDistanceFunction_Manhattan:
        for (int xi = xr - 1; xi <= xr + 1; xi++)
        {
            int yPrimed = yPrimedBase;

            for (int yi = yr - 1; yi <= yr + 1; yi++)
            {
                int zPrimed = zPrimedBase;

                for (int zi = zr - 1; zi <= zr + 1; zi++)
                {
                    int hash = Hash(seed, xPrimed, yPrimed, zPrimed);
                    int idx = hash & (255 << 2);

                    float vecX = (float)(xi - x) + RandVecs3D[idx] * cellularJitter;
                    float vecY = (float)(yi - y) + RandVecs3D[idx | 1] * cellularJitter;
                    float vecZ = (float)(zi - z) + RandVecs3D[idx | 2] * cellularJitter;

                    float newDistance = FastAbs(vecX) + FastAbs(vecY) + FastAbs(vecZ);

                    distance1 = FastMax(FastMin(distance1, newDistance), distance0);
                    if (newDistance < distance0)
                    {
                        distance0 = newDistance;
                        closestHash = hash;
                    }
                    zPrimed += PrimeZ;
                }
                yPrimed += PrimeY;
            }
            xPrimed += PrimeX;
        }
        break;
    case CellularDistanceFunction_Hybrid:
        for (int xi = xr - 1; xi <= xr + 1; xi++)
        {
            int yPrimed = yPrimedBase;

            for (int yi = yr - 1; yi <= yr + 1; yi++)
            {
                int zPrimed = zPrimedBase;

                for (int zi = zr - 1; zi <= zr + 1; zi++)
                {
                    int hash = Hash(seed, xPrimed, yPrimed, zPrimed);
                    int idx = hash & (255 << 2);

                    float vecX = (float)(xi - x) + RandVecs3D[idx] * cellularJitter;
                    float vecY = (float)(yi - y) + RandVecs3D[idx | 1] * cellularJitter;
                    float vecZ = (float)(zi - z) + RandVecs3D[idx | 2] * cellularJitter;

                    float newDistance = (FastAbs(vecX) + FastAbs(vecY) + FastAbs(vecZ)) + (vecX * vecX + vecY * vecY + vecZ * vecZ);

                    distance1 = FastMax(FastMin(distance1, newDistance), distance0);
                    if (newDistance < distance0)
                    {
                        distance0 = newDistance;
                        closestHash = hash;
                    }
                    zPrimed += PrimeZ;
                }
                yPrimed += PrimeY;
            }
            xPrimed += PrimeX;
        }
        break;
    default:
        break;
    }

    if (mCellularDistanceFunction == CellularDistanceFunction_Euclidean && mCellularReturnType >= CellularReturnType_Distance)
    {
        distance0 = FastSqrt(distance0);

        if (mCellularReturnType >= CellularReturnType_Distance2)
        {
            distance1 = FastSqrt(distance1);
        }
    }

    switch (mCellularReturnType)
    {
    case CellularReturnType_CellValue:
        return closestHash * (1 / 2147483648.0f);
    case CellularReturnType_Distance:
        return distance0 - 1;
    case CellularReturnType_Distance2:
        return distance1 - 1;
    case CellularReturnType_Distance2Add:
        return (distance1 + distance0) * 0.5f - 1;
    case CellularReturnType_Distance2Sub:
        return distance1 - distance0 - 1;
    case CellularReturnType_Distance2Mul:
        return distance1 * distance0 * 0.5f - 1;
    case CellularReturnType_Distance2Div:
        return distance0 / distance1 - 1;
    default:
        return 0;
    }
}


float FastNoiseLite::SinglePerlin(int seed, float x, float y) const
{
    int x0 = FastFloor(x);
    int y0 = FastFloor(y);

    float xd0 = (float)(x - x0);
    float yd0 = (float)(y - y0);
    float xd1 = xd0 - 1;
    float yd1 = yd0 - 1;

    float xs = InterpQuintic(xd0);
    float ys = InterpQuintic(yd0);

    x0 *= PrimeX;
    y0 *= PrimeY;
    int x1 = x0 + PrimeX;
    int y1 = y0 + PrimeY;

    float xf0 = Lerp(GradCoord(seed, x0, y0, xd0, yd0), GradCoord(seed, x1, y0, xd1, yd0), xs);
    float xf1 = Lerp(GradCoord(seed, x0, y1, xd0, yd1), GradCoord(seed, x1, y1, xd1, yd1), xs);

    return Lerp(xf0, xf1, ys) * 1.4247691104677813f;
}


float FastNoiseLite::SinglePerlin(int seed, float x, float y, float z) const
{
    int x0 = FastFloor(x);
    int y0 = FastFloor(y);
    int z0 = FastFloor(z);

    float xd0 = (float)(x - x0);
    float yd0 = (float)(y - y0);
    float zd0 = (float)(z - z0);
    float xd1 = xd0 - 1;
    float yd1 = yd0 - 1;
    float zd1 = zd0 - 1;

    float xs = InterpQuintic(xd0);
    float ys = InterpQuintic(yd0);
    float zs = InterpQuintic(zd0);

    x0 *= PrimeX;
    y0 *= PrimeY;
    z0 *= PrimeZ;
    int x1 = x0 + PrimeX;
    int y1 = y0 + PrimeY;
    int z1 = z0 + PrimeZ;

    float xf00 = Lerp(GradCoord(seed, x0, y0, z0, xd0, yd0, zd0), GradCoord(seed, x1, y0, z0, xd1, yd0, zd0), xs);
    float xf10 = Lerp(GradCoord(seed, x0, y1, z0, xd0, yd1, zd0), GradCoord(seed, x1, y1, z0, xd1, yd1, zd0), xs);
    float xf01 = Lerp(GradCoord(seed, x0, y0, z1, xd0, yd0, zd1), GradCoord(seed, x1, y0, z1, xd1, yd0, zd1), xs);
    float xf11 = Lerp(GradCoord(seed, x0, y1, z1, xd0, yd1, zd1), GradCoord(seed, x1, y1, z1, xd1, yd1, zd1), xs);

    float yf0 = Lerp(xf00, xf10, ys);
    float yf1 = Lerp(xf01, xf11, ys);

    return Lerp(yf0, yf1, zs) * 0.964921414852142333984375f;
}


float FastNoiseLite::SingleValueCubic(int seed, float x, float y) const
{
    int x1 = FastFloor(x);
    int y1 = FastFloor(y);

    float xs = (float)(x - x1);
    float ys = (float)(y - y1);

    x1 *= PrimeX;
    y1 *= PrimeY;
    int x0 = x1 - PrimeX;
    int y0 = y1 - PrimeY;
    int x2 = x1 + PrimeX;
    int y2 = y1 + PrimeY;
    int x3 = x1 + (int)((long)PrimeX << 1);
    int y3 = y1 + (int)((long)PrimeY << 1);

    return CubicLerp(CubicLerp(ValCoord(seed, x0, y0), ValCoord(seed, x1, y0), ValCoord(seed, x2, y0), ValCoord(seed, x3, y0), xs), CubicLerp(ValCoord(seed, x0, y1), ValCoord(seed, x1, y1), ValCoord(seed, x2, y1), ValCoord(seed, x3, y1), xs), CubicLerp(ValCoord(seed, x0, y2), ValCoord(seed, x1, y2), ValCoord(seed, x2, y2), ValCoord(seed, x3, y2), xs), CubicLerp(ValCoord(seed, x0, y3), ValCoord(seed, x1, y3), ValCoord(seed, x2, y3), ValCoord(seed, x3, y3), xs), ys) * (1 / (1.5f * 1.5f));
}


float FastNoiseLite::SingleValueCubic(int seed, float x, float y, float z) const
{
    int x1 = FastFloor(x);
    int y1 = FastFloor(y);
    int z1 = FastFloor(z);

    float xs = (float)(x - x1);
    float ys = (float)(y - y1);
    float zs = (float)(z - z1);

    x1 *= PrimeX;
    y1 *= PrimeY;
    z1 *= PrimeZ;

    int x0 = x1 - PrimeX;
    int y0 = y1 - PrimeY;
    int z0 = z1 - PrimeZ;
    int x2 = x1 + PrimeX;
    int y2 = y1 + PrimeY;
    int z2 = z1 + PrimeZ;
    int x3 = x1 + (int)((long)PrimeX << 1);
    int y3 = y1 + (int)((long)PrimeY << 1);
    int z3 = z1 + (int)((long)PrimeZ << 1);


    return CubicLerp(CubicLerp(CubicLerp(ValCoord(seed, x0, y0, z0), ValCoord(seed, x1, y0, z0), ValCoord(seed, x2, y0, z0), ValCoord(seed, x3, y0, z0), xs), CubicLerp(ValCoord(seed, x0, y1, z0), ValCoord(seed, x1, y1, z0), ValCoord(seed, x2, y1, z0), ValCoord(seed, x3, y1, z0), xs), CubicLerp(ValCoord(seed, x0, y2, z0), ValCoord(seed, x1, y2, z0), ValCoord(seed, x2, y2, z0), ValCoord(seed, x3, y2, z0), xs),
                               CubicLerp(ValCoord(seed, x0, y3, z0), ValCoord(seed, x1, y3, z0), ValCoord(seed, x2, y3, z0), ValCoord(seed, x3, y3, z0), xs), ys),
                     CubicLerp(CubicLerp(ValCoord(seed, x0, y0, z1), ValCoord(seed, x1, y0, z1), ValCoord(seed, x2, y0, z1), ValCoord(seed, x3, y0, z1), xs), CubicLerp(ValCoord(seed, x0, y1, z1), ValCoord(seed, x1, y1, z1), ValCoord(seed, x2, y1, z1), ValCoord(seed, x3, y1, z1), xs), CubicLerp(ValCoord(seed, x0, y2, z1), ValCoord(seed, x1, y2, z1), ValCoord(seed, x2, y2, z1), ValCoord(seed, x3, y2, z1), xs),
                               CubicLerp(ValCoord(seed, x0, y3, z1), ValCoord(seed, x1, y3, z1), ValCoord(seed, x2, y3, z1), ValCoord(seed, x3, y3, z1), xs), ys),
                     CubicLerp(CubicLerp(ValCoord(seed, x0, y0, z2), ValCoord(seed, x1, y0, z2), ValCoord(seed, x2, y0, z2), ValCoord(seed, x3, y0, z2), xs), CubicLerp(ValCoord(seed, x0, y1, z2), ValCoord(seed, x1, y1, z2), ValCoord(seed, x2, y1, z2), ValCoord(seed, x3, y1, z2), xs), CubicLerp(ValCoord(seed, x0, y2, z2), ValCoord(seed, x1, y2, z2), ValCoord(seed, x2, y2, z2), ValCoord(seed, x3, y2, z2), xs),
                               CubicLerp(ValCoord(seed, x0, y3, z2), ValCoord(seed, x1, y3, z2), ValCoord(seed, x2, y3, z2), ValCoord(seed, x3, y3, z2), xs), ys),
                     CubicLerp(CubicLerp(ValCoord(seed, x0, y0, z3), ValCoord(seed, x1, y0, z3), ValCoord(seed, x2, y0, z3), ValCoord(seed, x3, y0, z3), xs), CubicLerp(ValCoord(seed, x0, y1, z3), ValCoord(seed, x1, y1, z3), ValCoord(seed, x2, y1, z3), ValCoord(seed, x3, y1, z3), xs), CubicLerp(ValCoord(seed, x0, y2, z3), ValCoord(seed, x1, y2, z3), ValCoord(seed, x2, y2, z3), ValCoord(seed, x3, y2, z3), xs),
                               CubicLerp(ValCoord(seed, x0, y3, z3), ValCoord(seed, x1, y3, z3), ValCoord(seed, x2, y3, z3), ValCoord(seed, x3, y3, z3), xs), ys),
                     zs) *
        (1 / (1.5f * 1.5f * 1.5f));
}


float FastNoiseLite::SingleValue(int seed, float x, float y) const
{
    int x0 = FastFloor(x);
    int y0 = FastFloor(y);

    float xs = InterpHermite((float)(x - x0));
    float ys = InterpHermite((float)(y - y0));

    x0 *= PrimeX;
    y0 *= PrimeY;
    int x1 = x0 + PrimeX;
    int y1 = y0 + PrimeY;

    float xf0 = Lerp(ValCoord(seed, x0, y0), ValCoord(seed, x1, y0), xs);
    float xf1 = Lerp(ValCoord(seed, x0, y1), ValCoord(seed, x1, y1), xs);

    return Lerp(xf0, xf1, ys);
}


float FastNoiseLite::SingleValue(int seed, float x, float y, float z) const
{
    int x0 = FastFloor(x);
    int y0 = FastFloor(y);
    int z0 = FastFloor(z);

    float xs = InterpHermite((float)(x - x0));
    float ys = InterpHermite((float)(y - y0));
    float zs = InterpHermite((float)(z - z0));

    x0 *= PrimeX;
    y0 *= PrimeY;
    z0 *= PrimeZ;
    int x1 = x0 + PrimeX;
    int y1 = y0 + PrimeY;
    int z1 = z0 + PrimeZ;

    float xf00 = Lerp(ValCoord(seed, x0, y0, z0), ValCoord(seed, x1, y0, z0), xs);
    float xf10 = Lerp(ValCoord(seed, x0, y1, z0), ValCoord(seed, x1, y1, z0), xs);
    float xf01 = Lerp(ValCoord(seed, x0, y0, z1), ValCoord(seed, x1, y0, z1), xs);
    float xf11 = Lerp(ValCoord(seed, x0, y1, z1), ValCoord(seed, x1, y1, z1), xs);

    float yf0 = Lerp(xf00, xf10, ys);
    float yf1 = Lerp(xf01, xf11, ys);

    return Lerp(yf0, yf1, zs);
}


void FastNoiseLite::DoSingleDomainWarp(int seed, float amp, float freq, float x, float y, float& xr, float& yr) const
{
    switch (mDomainWarpType)
    {
    case DomainWarpType_OpenSimplex2:
        SingleDomainWarpSimplexGradient(seed, amp * 38.283687591552734375f, freq, x, y, xr, yr, false);
        break;
    case DomainWarpType_OpenSimplex2Reduced:
        SingleDomainWarpSimplexGradient(seed, amp * 16.0f, freq, x, y, xr, yr, true);
        break;
    case DomainWarpType_BasicGrid:
        SingleDomainWarpBasicGrid(seed, amp, freq, x, y, xr, yr);
        break;
    }
}


void FastNoiseLite::DoSingleDomainWarp(int seed, float amp, float freq, float x, float y, float z, float& xr, float& yr, float& zr) const
{
    switch (mDomainWarpType)
    {
    case DomainWarpType_OpenSimplex2:
        SingleDomainWarpOpenSimplex2Gradient(seed, amp * 32.69428253173828125f, freq, x, y, z, xr, yr, zr, false);
        break;
    case DomainWarpType_OpenSimplex2Reduced:
        SingleDomainWarpOpenSimplex2Gradient(seed, amp * 7.71604938271605f, freq, x, y, z, xr, yr, zr, true);
        break;
    case DomainWarpType_BasicGrid:
        SingleDomainWarpBasicGrid(seed, amp, freq, x, y, z, xr, yr, zr);
        break;
    }
}


void FastNoiseLite::DomainWarpSingle(float& x, float& y) const
{
    int seed = mSeed;
    float amp = mDomainWarpAmp * mFractalBounding;
    float freq = mFrequency;

    float xs = x;
    float ys = y;
    TransformDomainWarpCoordinate(xs, ys);

    DoSingleDomainWarp(seed, amp, freq, xs, ys, x, y);
}


void FastNoiseLite::DomainWarpSingle(float& x, float& y, float& z) const
{
    int seed = mSeed;
    float amp = mDomainWarpAmp * mFractalBounding;
    float freq = mFrequency;

    float xs = x;
    float ys = y;
    float zs = z;
    TransformDomainWarpCoordinate(xs, ys, zs);

    DoSingleDomainWarp(seed, amp, freq, xs, ys, zs, x, y, z);
}


void FastNoiseLite::DomainWarpFractalProgressive(float& x, float& y) const
{
    int seed = mSeed;
    float amp = mDomainWarpAmp * mFractalBounding;
    float freq = mFrequency;

    for (int i = 0; i < mOctaves; i++)
    {
        float xs = x;
        float ys = y;
        TransformDomainWarpCoordinate(xs, ys);

        DoSingleDomainWarp(seed, amp, freq, xs, ys, x, y);

        seed++;
        amp *= mGain;
        freq *= mLacunarity;
    }
}


void FastNoiseLite::DomainWarpFractalProgressive(float& x, float& y, float& z) const
{
    int seed = mSeed;
    float amp = mDomainWarpAmp * mFractalBounding;
    float freq = mFrequency;

    for (int i = 0; i < mOctaves; i++)
    {
        float xs = x;
        float ys = y;
        float zs = z;
        TransformDomainWarpCoordinate(xs, ys, zs);

        DoSingleDomainWarp(seed, amp, freq, xs, ys, zs, x, y, z);

        seed++;
        amp *= mGain;
        freq *= mLacunarity;
    }
}


void FastNoiseLite::DomainWarpFractalIndependent(float& x, float& y) const
{
    float xs = x;
    float ys = y;
    TransformDomainWarpCoordinate(xs, ys);

    int seed = mSeed;
    float amp = mDomainWarpAmp * mFractalBounding;
    float freq = mFrequency;

    for (int i = 0; i < mOctaves; i++)
    {
        DoSingleDomainWarp(seed, amp, freq, xs, ys, x, y);

        seed++;
        amp *= mGain;
        freq *= mLacunarity;
    }
}


void FastNoiseLite::DomainWarpFractalIndependent(float& x, float& y, float& z) const
{
    float xs = x;
    float ys = y;
    float zs = z;
    TransformDomainWarpCoordinate(xs, ys, zs);

    int seed = mSeed;
    float amp = mDomainWarpAmp * mFractalBounding;
    float freq = mFrequency;

    for (int i = 0; i < mOctaves; i++)
    {
        DoSingleDomainWarp(seed, amp, freq, xs, ys, zs, x, y, z);

        seed++;
        amp *= mGain;
        freq *= mLacunarity;
    }
}


void FastNoiseLite::SingleDomainWarpBasicGrid(int seed, float warpAmp, float frequency, float x, float y, float z, float& xr, float& yr, float& zr) const
{
    float xf = x * frequency;
    float yf = y * frequency;
    float zf = z * frequency;

    int x0 = FastFloor(xf);
    int y0 = FastFloor(yf);
    int z0 = FastFloor(zf);

    float xs = InterpHermite((float)(xf - x0));
    float ys = InterpHermite((float)(yf - y0));
    float zs = InterpHermite((float)(zf - z0));

    x0 *= PrimeX;
    y0 *= PrimeY;
    z0 *= PrimeZ;
    int x1 = x0 + PrimeX;
    int y1 = y0 + PrimeY;
    int z1 = z0 + PrimeZ;

    int hash0 = Hash(seed, x0, y0, z0) & (255 << 2);
    int hash1 = Hash(seed, x1, y0, z0) & (255 << 2);

    float lx0x = Lerp(RandVecs3D[hash0], RandVecs3D[hash1], xs);
    float ly0x = Lerp(RandVecs3D[hash0 | 1], RandVecs3D[hash1 | 1], xs);
    float lz0x = Lerp(RandVecs3D[hash0 | 2], RandVecs3D[hash1 | 2], xs);

    hash0 = Hash(seed, x0, y1, z0) & (255 << 2);
    hash1 = Hash(seed, x1, y1, z0) & (255 << 2);

    float lx1x = Lerp(RandVecs3D[hash0], RandVecs3D[hash1], xs);
    float ly1x = Lerp(RandVecs3D[hash0 | 1], RandVecs3D[hash1 | 1], xs);
    float lz1x = Lerp(RandVecs3D[hash0 | 2], RandVecs3D[hash1 | 2], xs);

    float lx0y = Lerp(lx0x, lx1x, ys);
    float ly0y = Lerp(ly0x, ly1x, ys);
    float lz0y = Lerp(lz0x, lz1x, ys);

    hash0 = Hash(seed, x0, y0, z1) & (255 << 2);
    hash1 = Hash(seed, x1, y0, z1) & (255 << 2);

    lx0x = Lerp(RandVecs3D[hash0], RandVecs3D[hash1], xs);
    ly0x = Lerp(RandVecs3D[hash0 | 1], RandVecs3D[hash1 | 1], xs);
    lz0x = Lerp(RandVecs3D[hash0 | 2], RandVecs3D[hash1 | 2], xs);

    hash0 = Hash(seed, x0, y1, z1) & (255 << 2);
    hash1 = Hash(seed, x1, y1, z1) & (255 << 2);

    lx1x = Lerp(RandVecs3D[hash0], RandVecs3D[hash1], xs);
    ly1x = Lerp(RandVecs3D[hash0 | 1], RandVecs3D[hash1 | 1], xs);
    lz1x = Lerp(RandVecs3D[hash0 | 2], RandVecs3D[hash1 | 2], xs);

    xr += Lerp(lx0y, Lerp(lx0x, lx1x, ys), zs) * warpAmp;
    yr += Lerp(ly0y, Lerp(ly0x, ly1x, ys), zs) * warpAmp;
    zr += Lerp(lz0y, Lerp(lz0x, lz1x, ys), zs) * warpAmp;
}


void FastNoiseLite::SingleDomainWarpSimplexGradient(int seed, float warpAmp, float frequency, float x, float y, float& xr, float& yr, bool outGradOnly) const
{
    const float SQRT3 = 1.7320508075688772935274463415059f;
    const float G2 = (3 - SQRT3) / 6;

    x *= frequency;
    y *= frequency;

    /*
     * --- Skew moved to TransformNoiseCoordinate method ---
     * const float F2 = 0.5f * (SQRT3 - 1);
     * float s = (x + y) * F2;
     * x += s; y += s;
     */

    int i = FastFloor(x);
    int j = FastFloor(y);
    float xi = (float)(x - i);
    float yi = (float)(y - j);

    float t = (xi + yi) * G2;
    float x0 = (float)(xi - t);
    float y0 = (float)(yi - t);

    i *= PrimeX;
    j *= PrimeY;

    float vx, vy;
    vx = vy = 0;

    float a = 0.5f - x0 * x0 - y0 * y0;
    if (a > 0)
    {
        float aaaa = (a * a) * (a * a);
        float xo, yo;
        if (outGradOnly)
            GradCoordOut(seed, i, j, xo, yo);
        else
            GradCoordDual(seed, i, j, x0, y0, xo, yo);
        vx += aaaa * xo;
        vy += aaaa * yo;
    }

    float c = (float)(2 * (1 - 2 * G2) * (1 / G2 - 2)) * t + ((float)(-2 * (1 - 2 * G2) * (1 - 2 * G2)) + a);
    if (c > 0)
    {
        float x2 = x0 + (2 * (float)G2 - 1);
        float y2 = y0 + (2 * (float)G2 - 1);
        float cccc = (c * c) * (c * c);
        float xo, yo;
        if (outGradOnly)
            GradCoordOut(seed, i + PrimeX, j + PrimeY, xo, yo);
        else
            GradCoordDual(seed, i + PrimeX, j + PrimeY, x2, y2, xo, yo);
        vx += cccc * xo;
        vy += cccc * yo;
    }

    if (y0 > x0)
    {
        float x1 = x0 + (float)G2;
        float y1 = y0 + ((float)G2 - 1);
        float b = 0.5f - x1 * x1 - y1 * y1;
        if (b > 0)
        {
            float bbbb = (b * b) * (b * b);
            float xo, yo;
            if (outGradOnly)
                GradCoordOut(seed, i, j + PrimeY, xo, yo);
            else
                GradCoordDual(seed, i, j + PrimeY, x1, y1, xo, yo);
            vx += bbbb * xo;
            vy += bbbb * yo;
        }
    }
    else
    {
        float x1 = x0 + ((float)G2 - 1);
        float y1 = y0 + (float)G2;
        float b = 0.5f - x1 * x1 - y1 * y1;
        if (b > 0)
        {
            float bbbb = (b * b) * (b * b);
            float xo, yo;
            if (outGradOnly)
                GradCoordOut(seed, i + PrimeX, j, xo, yo);
            else
                GradCoordDual(seed, i + PrimeX, j, x1, y1, xo, yo);
            vx += bbbb * xo;
            vy += bbbb * yo;
        }
    }

    xr += vx * warpAmp;
    yr += vy * warpAmp;
}


void FastNoiseLite::SingleDomainWarpOpenSimplex2Gradient(int seed, float warpAmp, float frequency, float x, float y, float z, float& xr, float& yr, float& zr, bool outGradOnly) const
{
    x *= frequency;
    y *= frequency;
    z *= frequency;

    /*
     * --- Rotation moved to TransformDomainWarpCoordinate method ---
     * const float R3 = (float)(2.0 / 3.0);
     * float r = (x + y + z) * R3; // Rotation, not skew
     * x = r - x; y = r - y; z = r - z;
     */

    int i = FastRound(x);
    int j = FastRound(y);
    int k = FastRound(z);
    float x0 = (float)x - i;
    float y0 = (float)y - j;
    float z0 = (float)z - k;

    int xNSign = (int)(-x0 - 1.0f) | 1;
    int yNSign = (int)(-y0 - 1.0f) | 1;
    int zNSign = (int)(-z0 - 1.0f) | 1;

    float ax0 = xNSign * -x0;
    float ay0 = yNSign * -y0;
    float az0 = zNSign * -z0;

    i *= PrimeX;
    j *= PrimeY;
    k *= PrimeZ;

    float vx, vy, vz;
    vx = vy = vz = 0;

    float a = (0.6f - x0 * x0) - (y0 * y0 + z0 * z0);
    for (int l = 0; l < 2; l++)
    {
        if (a > 0)
        {
            float aaaa = (a * a) * (a * a);
            float xo, yo, zo;
            if (outGradOnly)
                GradCoordOut(seed, i, j, k, xo, yo, zo);
            else
                GradCoordDual(seed, i, j, k, x0, y0, z0, xo, yo, zo);
            vx += aaaa * xo;
            vy += aaaa * yo;
            vz += aaaa * zo;
        }

        float b = a + 1;
        int i1 = i;
        int j1 = j;
        int k1 = k;
        float x1 = x0;
        float y1 = y0;
        float z1 = z0;

        if (ax0 >= ay0 && ax0 >= az0)
        {
            x1 += xNSign;
            b -= xNSign * 2 * x1;
            i1 -= xNSign * PrimeX;
        }
        else if (ay0 > ax0 && ay0 >= az0)
        {
            y1 += yNSign;
            b -= yNSign * 2 * y1;
            j1 -= yNSign * PrimeY;
        }
        else
        {
            z1 += zNSign;
            b -= zNSign * 2 * z1;
            k1 -= zNSign * PrimeZ;
        }

        if (b > 0)
        {
            float bbbb = (b * b) * (b * b);
            float xo, yo, zo;
            if (outGradOnly)
                GradCoordOut(seed, i1, j1, k1, xo, yo, zo);
            else
                GradCoordDual(seed, i1, j1, k1, x1, y1, z1, xo, yo, zo);
            vx += bbbb * xo;
            vy += bbbb * yo;
            vz += bbbb * zo;
        }

        if (l == 1)
            break;

        ax0 = 0.5f - ax0;
        ay0 = 0.5f - ay0;
        az0 = 0.5f - az0;

        x0 = xNSign * ax0;
        y0 = yNSign * ay0;
        z0 = zNSign * az0;

        a += (0.75f - ax0) - (ay0 + az0);

        i += (xNSign >> 1) & PrimeX;
        j += (yNSign >> 1) & PrimeY;
        k += (zNSign >> 1) & PrimeZ;

        xNSign = -xNSign;
        yNSign = -yNSign;
        zNSign = -zNSign;

        seed += 1293373;
    }

    xr += vx * warpAmp;
    yr += vy * warpAmp;
    zr += vz * warpAmp;
}

