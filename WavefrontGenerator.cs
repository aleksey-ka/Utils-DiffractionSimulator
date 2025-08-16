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
        /// <param name="sampling">Sampling factor (0.5 = 2x oversampling, 1.0 = normal, 2.0 = 2x undersampling, 4.0 = 4x undersampling, 8.0 = 8x undersampling, 16.0 = 16x undersampling, 32.0 = 32x undersampling)</param>
        /// <returns>2D array of z-positions representing the wavefront surface</returns>
        public static double[,] GenerateSphericalWavefront(double D, double f, double sampling = 1.0)
        {
            // Calculate effective array size based on sampling
            int effectiveArraySize = (int)(ARRAY_SIZE / sampling);
            double[,] wavefrontDepth = new double[effectiveArraySize, effectiveArraySize];
            
            // Calculate physical pixel size
            double aperturePixelSize = D / effectiveArraySize; // mm per pixel in aperture
            
            // Center of curvature is at z = f
            // For a sphere centered at (0, 0, f), the equation is: x² + y² + (z-f)² = f²
            // Solving for z: z = f - sqrt(f² - x² - y²)
            
            for (int i = 0; i < effectiveArraySize; i++)
            {
                for (int j = 0; j < effectiveArraySize; j++)
                {
                    // Convert pixel coordinates to physical coordinates (mm)
                    double x = (i - effectiveArraySize / 2.0) * aperturePixelSize;
                    double y = (j - effectiveArraySize / 2.0) * aperturePixelSize;
                    
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
