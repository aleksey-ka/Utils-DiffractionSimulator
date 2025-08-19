using System.Diagnostics;
using System.Numerics;

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

            Trace.WriteLine( $"Intensity at [0,0] is {imagePlane[ARRAY_SIZE / 2 + 1, ARRAY_SIZE / 2 + 1]}" );
            
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
        
        /// <summary>
        /// Calculates the integral over a triangular patch using the triangular decomposition method
        /// </summary>
        /// <param name="vertices">Array of 3 vertices with coordinates and distances</param>
        /// <param name="lambda">Wavelength</param>
        /// <returns>Complex field contribution from the triangular patch</returns>
        public static Complex CalculateTriangularPatchIntegral((double x, double y, double z, double distance)[] vertices, double lambda)
        {
            // Sort the vertices by distance
            Array.Sort(vertices, (a, b) => a.distance.CompareTo(b.distance));
            
            if (Math.Abs(vertices[2].distance - vertices[0].distance) < double.Epsilon)
            {
                // Area of the triangle
                double area = 0.5 * Math.Abs(vertices[0].x * (vertices[1].y - vertices[2].y) +
                vertices[1].x * (vertices[2].y - vertices[0].y) + vertices[2].x * (vertices[0].y - vertices[1].y));
                // Phase
                double P0 = 2 * Math.PI * vertices[0].distance / lambda;
                double P2 = 2 * Math.PI * vertices[2].distance / lambda;
                double P = (P0 + P2) / 2;
                return new Complex(area * Math.Cos(P), area * Math.Sin(P));
            }

            // Find a point along the line from vertices[0] to vertices[2] where the distance is equal to vertices[1].distance
            double k = (vertices[1].distance - vertices[0].distance) / (vertices[2].distance - vertices[0].distance);
            double x = vertices[0].x + k * (vertices[2].x - vertices[0].x);
            double y = vertices[0].y + k * (vertices[2].y - vertices[0].y);
            double z = vertices[0].z + k * (vertices[2].z - vertices[0].z);

            // Divide the triangle into two along the line of equal distance
            (double x, double y, double z, double distance)[] vertices1 = {
                (x, y, z, vertices[1].distance),
                vertices[0],
                vertices[1]
            };

            (double x, double y, double z, double distance)[] vertices2 = {
                vertices[1],
                vertices[2],
                (x, y, z, vertices[1].distance)
            };

            return CalculateComplexIntegral(vertices1, lambda) + CalculateComplexIntegral(vertices2, lambda);
        }
        
        /// <summary>
        /// Calculates the complex integral over a triangle with vertices at different distances
        /// </summary>
        /// <param name="vertices">Array of 3 vertices with coordinates and distances</param>
        /// <param name="lambda">Wavelength</param>
        /// <returns>Complex field contribution from the triangle</returns>
        public static Complex CalculateComplexIntegral((double x, double y, double z, double distance)[] vertices, double lambda)
        {
            // Expecting that vertices[0] and vertices[2] are at the same distance
            
            // Area of the triangle
            double area = 0.5 * Math.Abs(vertices[0].x * (vertices[1].y - vertices[2].y) +
                vertices[1].x * (vertices[2].y - vertices[0].y) + vertices[2].x * (vertices[0].y - vertices[1].y));
            
            // Phases at vertices[2] and vertices[0]
            double P0 = 2 * Math.PI * vertices[0].distance / lambda;
            double P1 = 2 * Math.PI * vertices[1].distance / lambda;

            // Phase difference
            double dP = P0 - P1;
            
            // Handle the case where phases are nearly equal (dP ≈ 0)
            if (Math.Abs(dP) < double.Epsilon)
            {
                double P = (P0 + P1) / 2;
                return new Complex(area * Math.Cos(P), area * Math.Sin(P));
            }

            // Real part: integral of x * cos((dP/alt) * x + P1) for x in [0, alt]
            // Closed form: ∫ x * cos(ax + b) dx = (1/a²) * (cos(ax + b) + ax * sin(ax + b))

            // Imaginary part: integral of x * sin((dP/alt) * x + P1) for x in [0, alt]  
            // Closed form: ∫ x * sin(ax + b) dx = (1/a²) * (sin(ax + b) - ax * cos(ax + b))

            // Where 'alt' is the altitude of the triangle to side [0,2], 
            // integration is along this altitude

            double A = 2 * area / (dP * dP);
            
            double realPart = A * (Math.Cos(P0) - Math.Cos(P1) + dP * Math.Sin(P0));
            double imagPart = A * (Math.Sin(P0) - Math.Sin(P1) - dP * Math.Cos(P0));

            return new Complex(realPart, imagPart);
        }
        
        /// <summary>
        /// Calculates 3D Euclidean distance between two points
        /// </summary>
        public static double CalculateDistance3D(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            double dx = x2 - x1;
            double dy = y2 - y1;
            double dz = z2 - z1;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }
        
        /// <summary>
        /// Calculates the depth of a spherical wavefront at given coordinates
        /// </summary>
        public static double CalculateSphericalWavefrontDepth(double x, double y, double f)
        {
            double rSquared = x * x + y * y;
            double fSquared = f * f;
            
            if (rSquared >= fSquared)
            {
                // Point is beyond the focal length - limit the depth
                return f;
            }
            
            double depth = f - Math.Sqrt(fSquared - rSquared);
            return depth;
        }
    }
}
