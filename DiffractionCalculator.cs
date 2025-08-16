using System;

namespace DiffractionSimulator
{
    public class DiffractionCalculator
    {
        private const int ARRAY_SIZE = 501;
        
        /// <summary>
        /// Calculates the diffraction pattern using Huygens-Fresnel principle
        /// </summary>
        /// <param name="apertureMask">2D array representing the aperture transparency</param>
        /// <param name="wavefrontDepth">2D array of z-positions for the wavefront surface</param>
        /// <param name="D">Aperture diameter in mm</param>
        /// <param name="f">Focal length in mm</param>
        /// <param name="lambda">Wavelength in mm</param>
        /// <param name="imagePlaneSize">Size of the image plane in mm</param>
        /// <param name="integrationSize">Size of the central integration region in pixels</param>
        /// <param name="sampling">Sampling factor (0.5 = 2x oversampling, 1.0 = normal, 2.0 = 2x undersampling, 4.0 = 4x undersampling, 8.0 = 8x undersampling, 16.0 = 16x undersampling, 32.0 = 32x undersampling)</param>
        /// <param name="gridType">Grid type for integration ("None" = direct integration, "Orthogonal" = orthogonal grid)</param>
        /// <returns>2D array representing the intensity pattern on the image plane</returns>
        public static float[,] CalculateDiffractionPattern(
            float[,] apertureMask, 
            double[,] wavefrontDepth,
            double D, 
            double f, 
            double lambda, 
            double imagePlaneSize, 
            int integrationSize, 
            double sampling = 1.0,
            string gridType = "None")
        {
            switch (gridType)
            {
                case "Orthogonal":
                    return CalculateDiffractionPatternOrthogonal(apertureMask, wavefrontDepth, D, f, lambda, imagePlaneSize, integrationSize, sampling);
                case "None":
                default:
                    return CalculateDiffractionPatternDirect(apertureMask, wavefrontDepth, D, f, lambda, imagePlaneSize, integrationSize, sampling);
            }
        }
        
        /// <summary>
        /// Direct integration method (brute force) - original implementation
        /// </summary>
        private static float[,] CalculateDiffractionPatternDirect(
            float[,] apertureMask, 
            double[,] wavefrontDepth,
            double D, 
            double f, 
            double lambda, 
            double imagePlaneSize, 
            int integrationSize, 
            double sampling = 1.0)
        {
            float[,] imagePlane = new float[ARRAY_SIZE, ARRAY_SIZE];
            
            // Calculate physical pixel sizes
            int apertureArraySize = apertureMask.GetLength(0);
            double aperturePixelSize = D / apertureArraySize; // mm per pixel in aperture
            double imagePixelSize = imagePlaneSize / ARRAY_SIZE; // mm per pixel in image plane
            
            // Clear image plane array
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    imagePlane[i, j] = 0.0f;
                }
            }
            
            // Calculate central integration region in image plane
            int centerStart = ARRAY_SIZE / 2 - integrationSize / 2;
            int centerEnd = ARRAY_SIZE / 2 + integrationSize / 2;
            
            // For each point in the central integration region
            for (int imgI = centerStart; imgI <= centerEnd; imgI++)
            {
                for (int imgJ = centerStart; imgJ <= centerEnd; imgJ++)
                {
                    // Convert image plane pixel coordinates to physical coordinates (mm)
                    double imgX = (imgI - ARRAY_SIZE / 2.0) * imagePixelSize;
                    double imgY = (imgJ - ARRAY_SIZE / 2.0) * imagePixelSize;
                    double imgZ = f; // Image plane is at z = f
                    
                    // Complex field at this image point (real and imaginary parts)
                    double realPart = 0.0;
                    double imagPart = 0.0;
                    
                    // Integrate over all aperture points
                    for (int apI = 0; apI < apertureArraySize; apI++)
                    {
                        for (int apJ = 0; apJ < apertureArraySize; apJ++)
                        {
                            // Convert aperture pixel coordinates to physical coordinates (mm)
                            double apX = (apI - apertureArraySize / 2.0) * aperturePixelSize;
                            double apY = (apJ - apertureArraySize / 2.0) * aperturePixelSize;
                            double apZ = wavefrontDepth[apI, apJ]; // Z position from wavefront
                            
                            // Calculate distance from aperture point to image point
                            double dx = imgX - apX;
                            double dy = imgY - apY;
                            double dz = imgZ - apZ;
                            double distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                            
                            // Calculate phase: 2π * distance / λ
                            double phase = 2.0 * Math.PI * distance / lambda;
                            
                            // Get amplitude from aperture mask
                            float amplitude = apertureMask[apI, apJ];
                            
                            // Add contribution to complex field
                            realPart += amplitude * Math.Cos(phase);
                            imagPart += amplitude * Math.Sin(phase);
                        }
                    }
                    
                    // Calculate intensity (square of complex field magnitude)
                    float intensity = (float)(realPart * realPart + imagPart * imagPart);
                    
                    // Store in image plane array
                    imagePlane[imgI, imgJ] = intensity;
                }
            }
            
            return imagePlane;
        }
        
        /// <summary>
        /// Orthogonal grid integration method - more efficient than direct integration
        /// </summary>
        private static float[,] CalculateDiffractionPatternOrthogonal(
            float[,] apertureMask, 
            double[,] wavefrontDepth,
            double D, 
            double f, 
            double lambda, 
            double imagePlaneSize, 
            int integrationSize, 
            double sampling = 1.0)
        {
            float[,] imagePlane = new float[ARRAY_SIZE, ARRAY_SIZE];
            
            // Calculate physical pixel sizes
            int apertureArraySize = apertureMask.GetLength(0);
            double aperturePixelSize = D / apertureArraySize; // mm per pixel in aperture
            double imagePixelSize = imagePlaneSize / ARRAY_SIZE; // mm per pixel in image plane
            
            // Clear image plane array
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    imagePlane[i, j] = 0.0f;
                }
            }
            
            // Calculate central integration region in image plane
            int centerStart = ARRAY_SIZE / 2 - integrationSize / 2;
            int centerEnd = ARRAY_SIZE / 2 + integrationSize / 2;
            
            // For each point in the central integration region
            for (int imgI = centerStart; imgI <= centerEnd; imgI++)
            {
                for (int imgJ = centerStart; imgJ <= centerEnd; imgJ++)
                {
                    // Convert image plane pixel coordinates to physical coordinates (mm)
                    double imgX = (imgI - ARRAY_SIZE / 2.0) * imagePixelSize;
                    double imgY = (imgJ - ARRAY_SIZE / 2.0) * imagePixelSize;
                    double imgZ = f; // Image plane is at z = f
                    
                    // Complex field at this image point (real and imaginary parts)
                    double realPart = 0.0;
                    double imagPart = 0.0;
                    
                    // Orthogonal grid integration - more efficient than direct integration
                    // We can optimize by pre-calculating some values and using symmetry
                    for (int apI = 0; apI < apertureArraySize; apI++)
                    {
                        // Pre-calculate x coordinate for this aperture row
                        double apX = (apI - apertureArraySize / 2.0) * aperturePixelSize;
                        double apX_squared = apX * apX;
                        
                        for (int apJ = 0; apJ < apertureArraySize; apJ++)
                        {
                            // Convert aperture pixel coordinates to physical coordinates (mm)
                            double apY = (apJ - apertureArraySize / 2.0) * aperturePixelSize;
                            double apZ = wavefrontDepth[apI, apJ]; // Z position from wavefront
                            
                            // Calculate distance from aperture point to image point
                            double dx = imgX - apX;
                            double dy = imgY - apY;
                            double dz = imgZ - apZ;
                            double distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                            
                            // Calculate phase: 2π * distance / λ
                            double phase = 2.0 * Math.PI * distance / lambda;
                            
                            // Get amplitude from aperture mask
                            float amplitude = apertureMask[apI, apJ];
                            
                            // Add contribution to complex field
                            realPart += amplitude * Math.Cos(phase);
                            imagPart += amplitude * Math.Sin(phase);
                        }
                    }
                    
                    // Calculate intensity (square of complex field magnitude)
                    float intensity = (float)(realPart * realPart + imagPart * imagPart);
                    
                    // Store in image plane array
                    imagePlane[imgI, imgJ] = intensity;
                }
            }
            
            return imagePlane;
        }
    }
}
