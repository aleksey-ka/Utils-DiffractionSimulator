using System;

namespace DiffractionSimulator
{
    public class ApertureMaskGenerator
    {
        private const int ARRAY_SIZE = 501;
        
        /// <summary>
        /// Generates an aperture mask based on the specified shape and parameters
        /// </summary>
        /// <param name="shape">The shape of the aperture ("Square", "Circular", or "Circular with obstr.")</param>
        /// <param name="D">Aperture diameter in mm</param>
        /// <param name="obstructionRatio">Obstruction ratio as percentage (0-100) for annular apertures</param>
        /// <returns>2D array representing the aperture transparency (1.0 = transparent, 0.0 = opaque)</returns>
        public static float[,] GenerateApertureMask(string shape, double D, double obstructionRatio = 30.0)
        {
            float[,] apertureMask = new float[ARRAY_SIZE, ARRAY_SIZE];
            
            switch (shape)
            {
                case "Circular with obstr.":
                    GenerateAnnularAperture(apertureMask, D, obstructionRatio);
                    break;
                case "Circular":
                    GenerateCircularAperture(apertureMask, D);
                    break;
                case "Square":
                default:
                    GenerateSquareAperture(apertureMask);
                    break;
            }
            
            return apertureMask;
        }
        
        private static void GenerateSquareAperture(float[,] apertureMask)
        {
            // Default square aperture (uniform)
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    apertureMask[i, j] = 1.0f;
                }
            }
        }
        
        private static void GenerateCircularAperture(float[,] apertureMask, double D)
        {
            // Create circular aperture
            double aperturePixelSize = D / ARRAY_SIZE; // mm per pixel in aperture
            double radius = D / 2.0; // radius in mm
            double centerX = ARRAY_SIZE / 2.0;
            double centerY = ARRAY_SIZE / 2.0;
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    // Convert pixel coordinates to physical coordinates (mm)
                    double x = (i - centerX) * aperturePixelSize;
                    double y = (j - centerY) * aperturePixelSize;
                    
                    // Check if point is within circle
                    double distanceFromCenter = Math.Sqrt(x * x + y * y);
                    if (distanceFromCenter <= radius)
                    {
                        apertureMask[i, j] = 1.0f;
                    }
                    else
                    {
                        apertureMask[i, j] = 0.0f;
                    }
                }
            }
        }
        
        private static void GenerateAnnularAperture(float[,] apertureMask, double D, double obstructionRatio)
        {
            // Create circular aperture with central obstruction
            double obstructionRatioDecimal = obstructionRatio / 100.0; // Convert percentage to ratio
            double aperturePixelSize = D / ARRAY_SIZE; // mm per pixel in aperture
            double outerRadius = D / 2.0; // outer radius in mm
            double innerRadius = (D * obstructionRatioDecimal) / 2.0; // inner radius (obstruction) in mm
            double centerX = ARRAY_SIZE / 2.0;
            double centerY = ARRAY_SIZE / 2.0;
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    // Convert pixel coordinates to physical coordinates (mm)
                    double x = (i - centerX) * aperturePixelSize;
                    double y = (j - centerY) * aperturePixelSize;
                    
                    // Check if point is within annular region (between inner and outer radius)
                    double distanceFromCenter = Math.Sqrt(x * x + y * y);
                    if (distanceFromCenter <= outerRadius && distanceFromCenter >= innerRadius)
                    {
                        apertureMask[i, j] = 1.0f;
                    }
                    else
                    {
                        apertureMask[i, j] = 0.0f;
                    }
                }
            }
        }
    }
}
