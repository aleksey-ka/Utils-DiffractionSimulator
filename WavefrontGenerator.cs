using System;

namespace DiffractionSimulator
{
    public class WavefrontGenerator
    {
        private const int ARRAY_SIZE = 501;
        
        /// <summary>
        /// Generates a spherical wavefront depth map for the aperture
        /// </summary>
        /// <param name="D">Aperture diameter in mm</param>
        /// <param name="f">Focal length in mm</param>
        /// <returns>2D array of z-positions representing the wavefront surface</returns>
        public static double[,] GenerateSphericalWavefront(double D, double f)
        {
            double[,] wavefrontDepth = new double[ARRAY_SIZE, ARRAY_SIZE];
            
            // Calculate physical pixel size
            double aperturePixelSize = D / ARRAY_SIZE; // mm per pixel in aperture
            
            // Center of curvature is at z = f
            // For a sphere centered at (0, 0, f), the equation is: x² + y² + (z-f)² = f²
            // Solving for z: z = f - sqrt(f² - x² - y²)
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    // Convert pixel coordinates to physical coordinates (mm)
                    double x = (i - ARRAY_SIZE / 2.0) * aperturePixelSize;
                    double y = (j - ARRAY_SIZE / 2.0) * aperturePixelSize;
                    
                    // Calculate distance from center of aperture
                    double rSquared = x * x + y * y;
                    
                    // Calculate z position on sphere surface
                    // z = f - sqrt(f² - r²) where r² = x² + y²
                    double z = f - Math.Sqrt(f * f - rSquared);
                    wavefrontDepth[i, j] = z;
                }
            }
            
            return wavefrontDepth;
        }
    }
}
