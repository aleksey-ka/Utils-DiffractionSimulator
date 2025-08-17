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
        /// <param name="sampling">Sampling factor (0.5 = 2x oversampling, 1.0 = normal, 2.0 = 2x undersampling, 4.0 = 4x undersampling, 8.0 = 8x undersampling, 16.0 = 16x undersampling, 32.0 = 32x undersampling)</param>
        /// <returns>2D array representing the aperture transparency (1.0 = transparent, 0.0 = opaque)</returns>
        public static float[,] GenerateApertureMask(string shape, double D, double obstructionRatio = 30.0, double sampling = 1.0)
        {
            // Calculate effective array size based on sampling
            int effectiveArraySize = (int)(ARRAY_SIZE / sampling);
            float[,] apertureMask = new float[effectiveArraySize, effectiveArraySize];
            
            switch (shape)
            {
                case "Circular with obstr.":
                    GenerateAnnularAperture(apertureMask, D, obstructionRatio, sampling);
                    break;
                case "Circular":
                    GenerateCircularAperture(apertureMask, D, sampling);
                    break;
                case "Triangular":
                    GenerateTriangularAperture(apertureMask, D, sampling);
                    break;
                case "Square":
                default:
                    GenerateSquareAperture(apertureMask, sampling);
                    break;
            }
            
            return apertureMask;
        }
        
        private static void GenerateSquareAperture(float[,] apertureMask, double sampling)
        {
            // Default square aperture (uniform)
            int arraySize = apertureMask.GetLength(0);
            for (int i = 0; i < arraySize; i++)
            {
                for (int j = 0; j < arraySize; j++)
                {
                    apertureMask[i, j] = 1.0f;
                }
            }
        }
        
        private static void GenerateCircularAperture(float[,] apertureMask, double D, double sampling)
        {
            // Create circular aperture
            int arraySize = apertureMask.GetLength(0);
            double aperturePixelSize = D / arraySize; // mm per pixel in aperture
            double radius = D / 2.0; // radius in mm
            double centerX = (arraySize - 1) / 2.0; // Center pixel for symmetric aperture
            double centerY = (arraySize - 1) / 2.0;
            
            for (int i = 0; i < arraySize; i++)
            {
                for (int j = 0; j < arraySize; j++)
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
        
        private static void GenerateAnnularAperture(float[,] apertureMask, double D, double obstructionRatio, double sampling)
        {
            // Create circular aperture with central obstruction
            int arraySize = apertureMask.GetLength(0);
            double obstructionRatioDecimal = obstructionRatio / 100.0; // Convert percentage to ratio
            double aperturePixelSize = D / arraySize; // mm per pixel in aperture
            double outerRadius = D / 2.0; // outer radius in mm
            double innerRadius = (D * obstructionRatioDecimal) / 2.0; // inner radius (obstruction) in mm
            double centerX = (arraySize - 1) / 2.0; // Center pixel for symmetric aperture
            double centerY = (arraySize - 1) / 2.0;
            
            for (int i = 0; i < arraySize; i++)
            {
                for (int j = 0; j < arraySize; j++)
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
        
        private static void GenerateTriangularAperture(float[,] apertureMask, double D, double sampling)
        {
            // Create equilateral triangle inscribed in the same circle as the circular aperture
            int arraySize = apertureMask.GetLength(0);
            double aperturePixelSize = D / arraySize; // mm per pixel in aperture
            double radius = D / 2.0; // radius of circumscribing circle in mm
            double centerX = (arraySize - 1) / 2.0; // Center pixel for symmetric aperture
            double centerY = (arraySize - 1) / 2.0;
            
            // For equilateral triangle inscribed in circle:
            // - Vertices are at radius distance from center
            // - First vertex points upward (negative Y direction) - vertical orientation
            // - Vertices are separated by 120 degrees
            
            // Calculate triangle vertices in physical coordinates (mm) - same as circular aperture
            // Equilateral triangle with one vertex pointing up (in physical space)
            // Since x maps to i (vertical) and y maps to j (horizontal), we need to rotate
            // Top vertex (pointing upward in physical space)
            double vertex1X_mm = -radius; // negative X is left
            double vertex1Y_mm = 0;
            
            // Bottom-left vertex 
            double vertex2X_mm = radius / 2.0;
            double vertex2Y_mm = -radius * Math.Sqrt(3) / 2.0;
            
            // Bottom-right vertex
            double vertex3X_mm = radius / 2.0;
            double vertex3Y_mm = radius * Math.Sqrt(3) / 2.0;
            
            for (int i = 0; i < arraySize; i++)
            {
                for (int j = 0; j < arraySize; j++)
                {
                    // Convert pixel coordinates to physical coordinates (mm) - same as circular aperture
                    double x = (i - centerX) * aperturePixelSize;
                    double y = (j - centerY) * aperturePixelSize;
                    
                    // Check if point (x, y) is inside the triangle using barycentric coordinates
                    if (IsPointInTriangle(x, y, vertex1X_mm, vertex1Y_mm, vertex2X_mm, vertex2Y_mm, vertex3X_mm, vertex3Y_mm))
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
        
        private static bool IsPointInTriangle(double px, double py, double ax, double ay, double bx, double by, double cx, double cy)
        {
            // Use barycentric coordinates to determine if point P is inside triangle ABC
            double denominator = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
            
            if (Math.Abs(denominator) < 1e-10) return false; // Degenerate triangle
            
            double alpha = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / denominator;
            double beta = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / denominator;
            double gamma = 1 - alpha - beta;
            
            // Point is inside triangle if all barycentric coordinates are non-negative
            return alpha >= 0 && beta >= 0 && gamma >= 0;
        }
    }
}
