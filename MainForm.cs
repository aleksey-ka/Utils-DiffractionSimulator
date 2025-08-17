using System.Drawing;
using System.Drawing.Imaging;
using System.Globalization;
using System.Linq;
using System.Diagnostics;

namespace DiffractionSimulator
{
    public partial class MainForm : Form
    {
        private float[,]? appertureMask;
        private double[,]? appertureDepth;
        private float[,]? imagePlane;
        private const int ARRAY_SIZE = 501;
        
        // Physical parameters are now controlled by UI controls

        public MainForm()
        {
            InitializeComponent();
            InitializeDataArray();
            DisplayArrayAsImage();
            DisplayImagePlane();

            apertureShapeComboBox.SelectedIndex = 0;
            samplingComboBox.SelectedIndex = 1; // Default to 1.0 sampling
            gridComboBox.SelectedIndex = 0; // Default to "None" (direct integration)
            UpdateObstructionControlState();
        }

        private void InitializeDataArray()
        {
            // Initialize arrays will be created by UpdateApertureMask based on sampling
            // Initialize with default circular aperture
            UpdateApertureMask();
            
            // Initialize wavefront and image plane arrays
            if (appertureMask == null)
            {
                throw new InvalidOperationException("Aperture mask was not initialized");
            }
            
            int arraySize = appertureMask.GetLength(0);
            appertureDepth = new double[arraySize, arraySize];
            imagePlane = new float[ARRAY_SIZE, ARRAY_SIZE]; // Image plane always stays at full resolution
            
            for (int i = 0; i < arraySize; i++)
            {
                for (int j = 0; j < arraySize; j++)
                {
                    appertureDepth[i, j] = 0.0f;
                }
            }
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    imagePlane[i, j] = 0.0f;
                }
            }
        }

        private void UpdateApertureMask()
        {
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Circular";
            double D_value = (double)D_numericUpDown.Value;
            double obstructionRatio = (double)obstructionRatioNumericUpDown.Value;
            double sampling = double.Parse(samplingComboBox.SelectedItem?.ToString() ?? "1.0", CultureInfo.InvariantCulture);
            
            // Generate aperture mask using the new class
            appertureMask = ApertureMaskGenerator.GenerateApertureMask(selectedShape, D_value, obstructionRatio, sampling);
            
            // Update obstruction control enabled state
            UpdateObstructionControlState();
        }

        private void UpdateObstructionControlState()
        {
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Circular";
            bool shouldShow = (selectedShape == "Circular with obstr.");
            
            obstructionRatioNumericUpDown.Visible = shouldShow;
            obstructionRatio_label.Visible = shouldShow;
        }

        private void DisplayArrayAsImage()
        {
            if (appertureMask == null)
            {
                return; // Nothing to display if mask is not initialized
            }
            
            int arraySize = appertureMask.GetLength(0);
            Bitmap bitmap = new Bitmap(ARRAY_SIZE, ARRAY_SIZE); // Always 501x501 for display
            
            // Scale the aperture data to fit the display size
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    // Map display coordinates to aperture coordinates
                    int apertureI = (int)((i * arraySize) / (double)ARRAY_SIZE);
                    int apertureJ = (int)((j * arraySize) / (double)ARRAY_SIZE);
                    
                    // Clamp to valid range
                    apertureI = Math.Max(0, Math.Min(arraySize - 1, apertureI));
                    apertureJ = Math.Max(0, Math.Min(arraySize - 1, apertureJ));
                    
                    int grayValue = (int)(appertureMask[apertureI, apertureJ] * 255);
                    Color pixelColor = Color.FromArgb(grayValue, grayValue, grayValue);
                    bitmap.SetPixel(j, i, pixelColor);
                }
            }
            
            // Draw grid dots if orthogonal grid mode is selected
            if (gridComboBox.SelectedItem?.ToString() == "Orthogonal")
            {
                DrawOrthogonalGrid(bitmap, arraySize);
            }
            
            apperturePictureBox.Image = bitmap;
        }
        
        private void DrawOrthogonalGrid(Bitmap bitmap, int arraySize)
        {
            if (appertureMask == null)
            {
                return; // Nothing to draw if mask is not initialized
            }
            
            // Get aperture parameters
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Circular";
            double D_value = (double)D_numericUpDown.Value;
            double obstructionRatio = (double)obstructionRatioNumericUpDown.Value;
            
            // Use the new GridMeshGenerator class
            GridMeshGenerator.GenerateOrthogonalGrid(bitmap, selectedShape, D_value, obstructionRatio);
        }
        

        

        
        private List<(double i, double j, bool isBoundary)> AddBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, (double centerI, double centerJ, double width, double height) bounds)
        {
            List<(double i, double j, bool isBoundary)> boundaryNodes = new List<(double i, double j, bool isBoundary)>();
            
            // Create rays from center at 5-degree intervals
            double angleStep = 5.0; // 5 degrees between rays
            int numRays = (int)(360.0 / angleStep); // 72 rays total
            
            for (int i = 0; i < numRays; i++)
            {
                double angleDegrees = i * angleStep;
                double angleRadians = angleDegrees * Math.PI / 180.0;
                
                // Cast ray from center and find boundary intersection - this never fails now
                var boundaryIntersection = FindRayBoundaryIntersection(centerI, centerJ, angleRadians, arraySize);
                
                // The improved algorithm always returns a valid position
                if (boundaryIntersection.HasValue)
                {
                    DrawDotAtPosition(bitmap, boundaryIntersection.Value.Item1, boundaryIntersection.Value.Item2, arraySize, Color.Red);
                    boundaryNodes.Add((boundaryIntersection.Value.Item1, boundaryIntersection.Value.Item2, true));
                }
            }
            
            return boundaryNodes;
        }
        
        private (double i, double j)? FindRayBoundaryIntersection(double centerI, double centerJ, double angle, int arraySize)
        {
            // Cast a ray from center in the given direction and find where it hits the boundary
            double directionI = Math.Cos(angle);
            double directionJ = Math.Sin(angle);
            
            // Start from center and march outward to find approximate boundary
            double step = 1.0; // Start with larger step for initial search
            double distance = 0.5; // Start very close to center
            double maxDistance = Math.Max(arraySize, 500);
            
            // First, find the approximate boundary with larger steps
            while (distance < maxDistance)
            {
                double testI = centerI + distance * directionI;
                double testJ = centerJ + distance * directionJ;
                
                if (!IsPositionIlluminatedDirect(testI, testJ, arraySize))
                {
                    // Found boundary transition - now use iterative refinement
                    return FindPreciseBoundaryWithIterativeRefinement(centerI, centerJ, directionI, directionJ, distance - step, distance, arraySize);
                }
                
                distance += step;
            }
            
            // If we reach here, the ray went to the array boundary without finding aperture boundary
            // Calculate exactly where this ray intersects the array boundary
            
            // Find intersection with array boundaries (at edges 0 and arraySize-1)
            double tMin = double.MaxValue;
            
            // Check intersection with left edge (i = 0)
            if (directionI < 0)
            {
                double t = -centerI / directionI;
                if (t > 0) tMin = Math.Min(tMin, t);
            }
            
            // Check intersection with right edge (i = arraySize-1)
            if (directionI > 0)
            {
                double t = (arraySize - 1 - centerI) / directionI;
                if (t > 0) tMin = Math.Min(tMin, t);
            }
            
            // Check intersection with top edge (j = 0)
            if (directionJ < 0)
            {
                double t = -centerJ / directionJ;
                if (t > 0) tMin = Math.Min(tMin, t);
            }
            
            // Check intersection with bottom edge (j = arraySize-1)
            if (directionJ > 0)
            {
                double t = (arraySize - 1 - centerJ) / directionJ;
                if (t > 0) tMin = Math.Min(tMin, t);
            }
            
            // Calculate the exact intersection point on the ray
            double edgeDistance = tMin * 0.95; // Step back slightly to be just inside
            double edgeI = centerI + edgeDistance * directionI;
            double edgeJ = centerJ + edgeDistance * directionJ;
            
            return (edgeI, edgeJ);
        }
        
        private (double i, double j)? FindPreciseBoundaryWithIterativeRefinement(double centerI, double centerJ, double directionI, double directionJ, double minDistance, double maxDistance, int arraySize)
        {
            // Iterative refinement: go back and forth with decreasing increments
            double currentDistance = minDistance; // Start from last known good position (inside)
            double increment = (maxDistance - minDistance) / 4.0; // Start with 1/4 of the gap
            bool movingForward = true; // Direction of movement
            int maxIterations = 20; // Prevent infinite loops
            
            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                // Test current position
                double testI = centerI + currentDistance * directionI;
                double testJ = centerJ + currentDistance * directionJ;
                bool isIlluminated = IsPositionIlluminatedDirect(testI, testJ, arraySize);
                
                if (movingForward)
                {
                    if (isIlluminated)
                    {
                        // Still inside, continue forward
                        currentDistance += increment;
                    }
                    else
                    {
                        // Hit boundary, reverse direction and reduce increment
                        currentDistance -= increment;
                        movingForward = false;
                        increment *= 0.5; // Halve the increment for finer precision
                    }
                }
                else
                {
                    if (!isIlluminated)
                    {
                        // Still outside, continue backward
                        currentDistance -= increment;
                    }
                    else
                    {
                        // Back inside, reverse direction and reduce increment
                        currentDistance += increment;
                        movingForward = true;
                        increment *= 0.5; // Halve the increment for finer precision
                    }
                }
                
                // Stop when increment becomes very small (high precision achieved)
                if (increment < 0.005) // 0.005 pixel precision
                {
                    break;
                }
                
                // Ensure we don't go beyond bounds
                currentDistance = Math.Max(0, Math.Min(maxDistance * 1.1, currentDistance));
            }
            
            // Final position should be just inside the illuminated area
            double finalI = centerI + currentDistance * directionI;
            double finalJ = centerJ + currentDistance * directionJ;
            
            // Verify the final position is illuminated
            if (IsPositionIlluminatedDirect(finalI, finalJ, arraySize))
            {
                return (finalI, finalJ);
            }
            
            // If not illuminated, step back slightly
            currentDistance -= 0.01;
            finalI = centerI + currentDistance * directionI;
            finalJ = centerJ + currentDistance * directionJ;
            
            return (finalI, finalJ);
        }
        

        

        
        private (double i, double j)? FindPreciseBoundaryPosition(double centerI, double centerJ, double angle, int arraySize)
        {
            // Use binary search to find the precise boundary position
            double minDistance = 0;
            double maxDistance = 300; // Start with a large search range
            
            // First, find a rough range where the boundary exists
            double testDistance = 0;
            while (testDistance < maxDistance)
            {
                double testI = centerI + testDistance * Math.Cos(angle);
                double testJ = centerJ + testDistance * Math.Sin(angle);
                
                if (!IsPositionIlluminated(testI, testJ, arraySize))
                {
                    maxDistance = testDistance;
                    break;
                }
                testDistance += 10; // Coarse steps first
            }
            
            if (maxDistance >= 300) return null; // No boundary found
            
            // Now use binary search for precision
            double precision = 0.1; // Sub-pixel precision
            
            while (maxDistance - minDistance > precision)
            {
                double midDistance = (minDistance + maxDistance) / 2.0;
                double testI = centerI + midDistance * Math.Cos(angle);
                double testJ = centerJ + midDistance * Math.Sin(angle);
                
                if (IsPositionIlluminated(testI, testJ, arraySize))
                {
                    minDistance = midDistance; // Point is inside, move outward
                }
                else
                {
                    maxDistance = midDistance; // Point is outside, move inward
                }
            }
            
            // Final boundary position (just inside the illuminated area)
            double finalDistance = minDistance + precision * 0.5; // Ensure we're just inside
            double boundaryI = centerI + finalDistance * Math.Cos(angle);
            double boundaryJ = centerJ + finalDistance * Math.Sin(angle);
            
            // Verify the final position is illuminated
            if (IsPositionIlluminated(boundaryI, boundaryJ, arraySize))
            {
                return (boundaryI, boundaryJ);
            }
            
            return null;
        }
        
        private void DrawDotAtPosition(Bitmap bitmap, double apertureI, double apertureJ, int arraySize, Color dotColor)
        {
            // Convert aperture coordinates to display coordinates
            int displayI = (int)((apertureI * ARRAY_SIZE) / (double)arraySize);
            int displayJ = (int)((apertureJ * ARRAY_SIZE) / (double)arraySize);
            
            // Clamp to valid display range
            displayI = Math.Max(0, Math.Min(ARRAY_SIZE - 1, displayI));
            displayJ = Math.Max(0, Math.Min(ARRAY_SIZE - 1, displayJ));
            
            // Draw a colored dot (3x3 pixels for visibility)
            for (int di = -1; di <= 1; di++)
            {
                for (int dj = -1; dj <= 1; dj++)
                {
                    int dotI = displayI + di;
                    int dotJ = displayJ + dj;
                    
                    if (dotI >= 0 && dotI < ARRAY_SIZE && dotJ >= 0 && dotJ < ARRAY_SIZE)
                    {
                        bitmap.SetPixel(dotJ, dotI, dotColor);
                    }
                }
            }
        }
        
        private bool IsPositionIlluminated(double i, double j, int arraySize)
        {
            // Check bounds first - positions outside the array are not illuminated
            if (i < 0 || i >= arraySize || j < 0 || j >= arraySize)
            {
                return false;
            }
            
            // Use bilinear interpolation for sub-pixel precision
            int i0 = (int)Math.Floor(i);
            int i1 = (int)Math.Ceiling(i);
            int j0 = (int)Math.Floor(j);
            int j1 = (int)Math.Ceiling(j);
            
            // Clamp coordinates to valid array range
            i0 = Math.Max(0, Math.Min(arraySize - 1, i0));
            i1 = Math.Max(0, Math.Min(arraySize - 1, i1));
            j0 = Math.Max(0, Math.Min(arraySize - 1, j0));
            j1 = Math.Max(0, Math.Min(arraySize - 1, j1));
            
            // Handle edge case where i0 == i1 or j0 == j1
            if (i0 == i1 && j0 == j1)
            {
                return appertureMask![i0, j0] > 0.001; // Use small threshold for precision
            }
            
            // Calculate interpolation weights
            double wi = (i - i0) / (i1 - i0);
            double wj = (j - j0) / (j1 - j0);
            
            // Handle edge cases
            if (i0 == i1)
            {
                // Linear interpolation in j direction only
                double value = (1 - wj) * appertureMask![i0, j0] + wj * appertureMask[i0, j1];
                return value > 0.001; // Use small threshold for precision
            }
            else if (j0 == j1)
            {
                // Linear interpolation in i direction only
                double value = (1 - wi) * appertureMask![i0, j0] + wi * appertureMask[i1, j0];
                return value > 0.001; // Use small threshold for precision
            }
            else
            {
                // Bilinear interpolation
                double value = (1 - wi) * (1 - wj) * appertureMask![i0, j0] +
                              wi * (1 - wj) * appertureMask[i1, j0] +
                              (1 - wi) * wj * appertureMask[i0, j1] +
                              wi * wj * appertureMask[i1, j1];
                return value > 0.001; // Use small threshold for precision
            }
        }
        
        private bool IsPositionIlluminatedDirect(double i, double j, int arraySize)
        {
            // Check bounds first - positions outside the array are not illuminated
            if (i < 0 || i >= arraySize || j < 0 || j >= arraySize)
            {
                return false;
            }
            
            // Use the same mathematical formulas as the original aperture generation
            // Get current aperture parameters
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Circular";
            double D_value = (double)D_numericUpDown.Value;
            double obstructionRatio = (double)obstructionRatioNumericUpDown.Value;
            double sampling = double.Parse(samplingComboBox.SelectedItem?.ToString() ?? "1.0", CultureInfo.InvariantCulture);
            
            // Calculate physical coordinates using the same logic as ApertureMaskGenerator
            double aperturePixelSize = D_value / arraySize; // mm per pixel in aperture
            double centerX = arraySize / 2.0;
            double centerY = arraySize / 2.0;
            
            // Convert pixel coordinates to physical coordinates (mm)
            double x = (i - centerX) * aperturePixelSize;
            double y = (j - centerY) * aperturePixelSize;
            
            switch (selectedShape)
            {
                case "Circular":
                    {
                        double radius = D_value / 2.0; // radius in mm
                        double distanceFromCenter = Math.Sqrt(x * x + y * y);
                        return distanceFromCenter <= radius;
                    }
                    
                case "Circular with obstr.":
                    {
                        double obstructionRatioDecimal = obstructionRatio / 100.0; // Convert percentage to ratio
                        double outerRadius = D_value / 2.0; // outer radius in mm
                        double innerRadius = (D_value * obstructionRatioDecimal) / 2.0; // inner radius (obstruction) in mm
                        double distanceFromCenter = Math.Sqrt(x * x + y * y);
                        return distanceFromCenter <= outerRadius && distanceFromCenter >= innerRadius;
                    }
                    
                case "Square":
                default:
                    // For square aperture, everything within the array is illuminated
                    return true;
            }
        }
        
        private bool IsPositionForBoundaryDetection(double i, double j, int arraySize)
        {
            // Check bounds first - positions outside the array are not illuminated
            if (i < 0 || i >= arraySize || j < 0 || j >= arraySize)
            {
                return false;
            }
            
            // For boundary detection, use the same logic as IsPositionIlluminated
            // but this method is specifically for finding the true aperture boundary
            if (appertureMask == null) return false;
            
            // Use bilinear interpolation for smooth boundary detection
            int i0 = (int)Math.Floor(i);
            int i1 = (int)Math.Ceiling(i);
            int j0 = (int)Math.Floor(j);
            int j1 = (int)Math.Ceiling(j);
            
            // Clamp coordinates to valid array range
            i0 = Math.Max(0, Math.Min(arraySize - 1, i0));
            i1 = Math.Max(0, Math.Min(arraySize - 1, i1));
            j0 = Math.Max(0, Math.Min(arraySize - 1, j0));
            j1 = Math.Max(0, Math.Min(arraySize - 1, j1));
            
            // Handle edge case where i0 == i1 or j0 == j1
            if (i0 == i1 && j0 == j1)
            {
                return appertureMask[i0, j0] > 0.5; // Use 0.5 as the boundary threshold
            }
            
            // Calculate interpolation weights
            double wi = (i - i0) / (i1 - i0);
            double wj = (j - j0) / (j1 - j0);
            
            // Handle edge cases
            if (i0 == i1)
            {
                // Linear interpolation in j direction only
                double value = (1 - wj) * appertureMask[i0, j0] + wj * appertureMask[i0, j1];
                return value > 0.5; // 0.5 threshold for boundary detection
            }
            else if (j0 == j1)
            {
                // Linear interpolation in i direction only
                double value = (1 - wi) * appertureMask[i0, j0] + wi * appertureMask[i1, j0];
                return value > 0.5; // 0.5 threshold for boundary detection
            }
            else
            {
                // Bilinear interpolation
                double value = (1 - wi) * (1 - wj) * appertureMask[i0, j0] +
                              wi * (1 - wj) * appertureMask[i1, j0] +
                              (1 - wi) * wj * appertureMask[i0, j1] +
                              wi * wj * appertureMask[i1, j1];
                return value > 0.5; // 0.5 threshold for boundary detection
            }
        }
        
        private bool IsPositionOnBoundary(double i, double j, int arraySize)
        {
            if (appertureMask == null) return false;
            
            // Check if this position is illuminated
            if (!IsPositionIlluminated(i, j, arraySize)) return false;
            
            // Check if any of the 8 surrounding positions are not illuminated
            for (int di = -1; di <= 1; di++)
            {
                for (int dj = -1; dj <= 1; dj++)
                {
                    if (di == 0 && dj == 0) continue; // Skip center position
                    
                    double testI = i + di;
                    double testJ = j + dj;
                    
                    if (!IsPositionIlluminated(testI, testJ, arraySize))
                    {
                        return true; // This position is on the boundary
                    }
                }
            }
            
            return false; // Position is surrounded by illuminated pixels (interior)
        }
        
        private void DisplayImagePlane()
        {
            if (imagePlane == null)
            {
                return; // Nothing to display if image plane is not initialized
            }
            
            Bitmap bitmap = new Bitmap(ARRAY_SIZE, ARRAY_SIZE);
            
            // Find maximum intensity for normalization
            float maxIntensity = 0.0f;
            float minIntensity = float.MaxValue;
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    if (imagePlane[i, j] > maxIntensity)
                        maxIntensity = imagePlane[i, j];
                    if (imagePlane[i, j] > 0 && imagePlane[i, j] < minIntensity)
                        minIntensity = imagePlane[i, j];
                }
            }
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    // Normalize intensity to 0-255 range for display
                    int grayValue = 0;
                    if (maxIntensity > 0)
                    {
                        if (gammaCorrectedViewCheckBox.Checked)
                        {
                            // Gamma correction: (intensity/maxIntensity)^gamma
                            // Gamma < 1 expands dark areas, compresses bright areas
                            float gamma = 0.3f; // Good for diffraction patterns
                            float normalizedValue = (float)Math.Pow(imagePlane[i, j] / maxIntensity, gamma);
                            grayValue = (int)(Math.Max(0, Math.Min(255, normalizedValue * 255)));
                        }
                        else
                        {
                            // Linear scaling
                            grayValue = (int)((imagePlane[i, j] / maxIntensity) * 255);
                        }
                    }
                    Color pixelColor = Color.FromArgb(grayValue, grayValue, grayValue);
                    bitmap.SetPixel(j, i, pixelColor);
                }
            }
            
            imagePlanePictureBox.Image = bitmap;
        }
        
        private void CalculateDiffractionPattern()
        {
            if (appertureMask == null)
            {
                return; // Cannot calculate if aperture mask is not initialized
            }
            
            // Get values from controls
            double D_value = (double)D_numericUpDown.Value;
            double f_value = (double)f_numericUpDown.Value;
            double lambda_value = (double)lambda_numericUpDown.Value;
            double imagePlaneSize_value = (double)imagePlaneSize_numericUpDown.Value;
            int integrationSize = (int)integrationSize_numericUpDown.Value;
            double sampling = double.Parse(samplingComboBox.SelectedItem?.ToString() ?? "1.0", CultureInfo.InvariantCulture);
            string gridType = gridComboBox.SelectedItem?.ToString() ?? "None";
            
            // Generate wavefront using the new class
            appertureDepth = WavefrontGenerator.GenerateSphericalWavefront(D_value, f_value, sampling);
            
            // Calculate diffraction pattern using the new class
            imagePlane = DiffractionCalculator.CalculateDiffractionPattern(
                appertureMask, 
                appertureDepth, 
                D_value, 
                f_value, 
                lambda_value, 
                imagePlaneSize_value, 
                integrationSize, 
                sampling,
                gridType);
        }
        
        private void GammaCorrectedViewCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            // Refresh the image plane display when checkbox state changes
            DisplayImagePlane();
        }
        
        private void ApertureShapeComboBox_SelectedIndexChanged(object sender, EventArgs e)
        {
            // Update aperture mask when shape changes
            UpdateApertureMask();
            DisplayArrayAsImage();
        }
        
        private void ObstructionRatioNumericUpDown_ValueChanged(object sender, EventArgs e)
        {
            // Update aperture mask when obstruction ratio changes (only for circular with obstruction)
            if (apertureShapeComboBox.SelectedItem?.ToString() == "Circular with obstr.")
            {
                UpdateApertureMask();
                DisplayArrayAsImage();
            }
        }
        
        private void SamplingComboBox_SelectedIndexChanged(object sender, EventArgs e)
        {
            // Update aperture mask when sampling changes
            UpdateApertureMask();
            DisplayArrayAsImage();
        }
        
        private void GridComboBox_SelectedIndexChanged(object sender, EventArgs e)
        {
            // Grid selection changed - refresh the aperture display to show/hide grid
            DisplayArrayAsImage();
        }
        
        private void CalculateButton_Click(object sender, EventArgs e)
        {
            // Calculate diffraction pattern with current parameters
            CalculateDiffractionPattern();
            DisplayImagePlane();
        }

        private void GenerateAndDrawTriangularMesh(Bitmap bitmap, int arraySize, List<(double i, double j, bool isBoundary)> nodes)
        {
            if (nodes.Count < 3) 
            {
                Trace.WriteLine($"Not enough nodes for triangulation: {nodes.Count}");
                return;
            }
            
            Trace.WriteLine($"Starting mesh generation with {nodes.Count} nodes");
            
            // Generate Delaunay triangulation
            var triangles = GenerateDelaunayTriangulation(nodes);
            
            Trace.WriteLine($"Generated {triangles.Count} triangles");
            
            if (triangles.Count > 0)
            {
                // Draw the triangular mesh
                DrawTriangularMesh(bitmap, arraySize, nodes, triangles);
                Trace.WriteLine($"Mesh drawn successfully with {triangles.Count} triangles");
            }
            else
            {
                Trace.WriteLine("No triangles generated - mesh will not be drawn");
            }
        }
        
        private List<(int node1, int node2, int node3)> GenerateDelaunayTriangulation(List<(double i, double j, bool isBoundary)> nodes)
        {
            if (nodes.Count < 3) return new List<(int, int, int)>();
            
            Trace.WriteLine($"Starting Delaunay triangulation with {nodes.Count} nodes");
            
            // Create a super-triangle that contains all nodes
            double minX = nodes.Min(n => n.i) - 10;
            double maxX = nodes.Max(n => n.i) + 10;
            double minY = nodes.Min(n => n.j) - 10;
            double maxY = nodes.Max(n => n.j) + 10;
            
            double dx = maxX - minX;
            double dy = maxY - minY;
            
            // Create super-triangle vertices (much larger than the point set)
            var superTriangle = new List<(double i, double j)>
            {
                (minX - dx, minY - dy),
                (maxX + dx, minY - dy),
                (minX + dx/2, maxY + dy)
            };
            
            // Initialize with super-triangle
            var triangles = new List<(int v1, int v2, int v3)> { (nodes.Count, nodes.Count + 1, nodes.Count + 2) };
            var allNodes = nodes.Select(n => (n.i, n.j)).ToList();
            allNodes.AddRange(superTriangle);
            
            // Insert each point one by one
            for (int pointIndex = 0; pointIndex < nodes.Count; pointIndex++)
            {
                var point = allNodes[pointIndex];
                var badTriangles = new List<int>();
                
                // Find triangles whose circumcircle contains the point
                for (int t = 0; t < triangles.Count; t++)
                {
                    var triangle = triangles[t];
                    if (IsPointInCircumcircle(point, allNodes[triangle.v1], allNodes[triangle.v2], allNodes[triangle.v3]))
                    {
                        badTriangles.Add(t);
                    }
                }
                
                // Find the boundary of the polygonal hole
                var polygon = new List<(int v1, int v2)>();
                foreach (int badTriangle in badTriangles)
                {
                    var triangle = triangles[badTriangle];
                    var edges = new List<(int v1, int v2)>
                    {
                        (triangle.v1, triangle.v2),
                        (triangle.v2, triangle.v3),
                        (triangle.v3, triangle.v1)
                    };
                    
                    foreach (var edge in edges)
                    {
                        var reverseEdge = (edge.v2, edge.v1);
                        if (polygon.Contains(reverseEdge))
                        {
                            polygon.Remove(reverseEdge);
                        }
                        else
                        {
                            polygon.Add(edge);
                        }
                    }
                }
                
                // Remove bad triangles
                badTriangles.Sort((a, b) => b.CompareTo(a)); // Sort descending
                foreach (int badTriangle in badTriangles)
                {
                    triangles.RemoveAt(badTriangle);
                }
                
                // Add new triangles formed by connecting the point to the polygon boundary
                foreach (var edge in polygon)
                {
                    triangles.Add((pointIndex, edge.v1, edge.v2));
                }
            }
            
            // Remove triangles that contain super-triangle vertices
            var finalTriangles = triangles.Where(t => 
                t.v1 < nodes.Count && t.v2 < nodes.Count && t.v3 < nodes.Count
            ).ToList();
            
            Trace.WriteLine($"Delaunay triangulation created {finalTriangles.Count} triangles from {nodes.Count} nodes");
            return finalTriangles.Select(t => (t.v1, t.v2, t.v3)).ToList();
        }
        
        private bool IsPointInCircumcircle((double i, double j) p, (double i, double j) a, (double i, double j) b, (double i, double j) c)
        {
            // Calculate circumcircle using determinant method
            double ax = a.i - p.i;
            double ay = a.j - p.j;
            double bx = b.i - p.i;
            double by = b.j - p.j;
            double cx = c.i - p.i;
            double cy = c.j - p.j;
            
            double det = (ax * ax + ay * ay) * (bx * cy - by * cx) -
                        (bx * bx + by * by) * (ax * cy - ay * cx) +
                        (cx * cx + cy * cy) * (ax * by - ay * bx);
            
            return det > 0;
        }
        
        private void DrawTriangularMesh(Bitmap bitmap, int arraySize, List<(double i, double j, bool isBoundary)> nodes, List<(int node1, int node2, int node3)> triangles)
        {
            // Draw triangles with visible lines
            Color meshColor = Color.Green;
            
            foreach (var triangle in triangles)
            {
                var node1 = nodes[triangle.node1];
                var node2 = nodes[triangle.node2];
                var node3 = nodes[triangle.node3];
                
                // Draw triangle edges
                DrawLine(bitmap, arraySize, node1.i, node1.j, node2.i, node2.j, meshColor);
                DrawLine(bitmap, arraySize, node2.i, node2.j, node3.i, node3.j, meshColor);
                DrawLine(bitmap, arraySize, node3.i, node3.j, node1.i, node1.j, meshColor);
            }
        }
        
        private void DrawLine(Bitmap bitmap, int arraySize, double x1, double y1, double x2, double y2, Color color)
        {
            // Convert aperture coordinates to display coordinates
            int x0 = (int)((x1 * ARRAY_SIZE) / (double)arraySize);
            int y0 = (int)((y1 * ARRAY_SIZE) / (double)arraySize);
            int x1_int = (int)((x2 * ARRAY_SIZE) / (double)arraySize);
            int y1_int = (int)((y2 * ARRAY_SIZE) / (double)arraySize);
            
            // Clamp coordinates
            x0 = Math.Max(0, Math.Min(ARRAY_SIZE - 1, x0));
            y0 = Math.Max(0, Math.Min(ARRAY_SIZE - 1, y0));
            x1_int = Math.Max(0, Math.Min(ARRAY_SIZE - 1, x1_int));
            y1_int = Math.Max(0, Math.Min(ARRAY_SIZE - 1, y1_int));
            
            // Bresenham's line algorithm
            int dx = Math.Abs(x1_int - x0);
            int dy = Math.Abs(y1_int - y0);
            int sx = x0 < x1_int ? 1 : -1;
            int sy = y0 < y1_int ? 1 : -1;
            int err = dx - dy;
            
            int x = x0, y = y0;
            
            while (true)
            {
                // Draw a single pixel line
                if (x >= 0 && x < ARRAY_SIZE && y >= 0 && y < ARRAY_SIZE)
                {
                    bitmap.SetPixel(y, x, color);
                }
                
                if (x == x1_int && y == y1_int) break;
                
                int e2 = 2 * err;
                if (e2 > -dy)
                {
                    err -= dy;
                    x += sx;
                }
                if (e2 < dx)
                {
                    err += dx;
                    y += sy;
                }
            }
        }
        
        private void TriangularDiffractionButton_Click(object sender, EventArgs e)
        {
            // Get the current triangular aperture parameters
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Circular";
            if (selectedShape != "Triangular")
            {
                MessageBox.Show("Please select 'Triangular' aperture shape first.", "Invalid Aperture", MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return;
            }
            
            // Get physical parameters
            double D_value = (double)D_numericUpDown.Value;
            double f_value = (double)f_numericUpDown.Value;
            double lambda_value = (double)lambda_numericUpDown.Value;
            double imagePlaneSize_value = (double)imagePlaneSize_numericUpDown.Value;
            
            // Calculate triangle vertices (same as in ApertureMaskGenerator)
            double radius = D_value / 2.0; // radius of circumscribing circle in mm
            
            // Left vertex (pointing leftward in physical space)
            double vertex1X_mm = -radius;
            double vertex1Y_mm = 0;
            
            // Bottom-left vertex 
            double vertex2X_mm = radius / 2.0;
            double vertex2Y_mm = -radius * Math.Sqrt(3) / 2.0;
            
            // Bottom-right vertex
            double vertex3X_mm = radius / 2.0;
            double vertex3Y_mm = radius * Math.Sqrt(3) / 2.0;
            
            // Call the triangular diffraction calculation method
            CalculateTriangularDiffraction(vertex1X_mm, vertex1Y_mm, vertex2X_mm, vertex2Y_mm, vertex3X_mm, vertex3Y_mm, 
                D_value, f_value, lambda_value, imagePlaneSize_value);
        }
        
        private void CalculateTriangularDiffraction(double v1x, double v1y, double v2x, double v2y, double v3x, double v3y,
            double D, double f, double lambda, double imagePlaneSize)
        {
            // Calculate z-values (depth) for each vertex using the same spherical wavefront logic
            // as used in WavefrontGenerator.GenerateSphericalWavefront
            
            // Calculate depth for each vertex using the spherical wavefront formula
            double v1z = CalculateSphericalWavefrontDepth(v1x, v1y, f);
            double v2z = CalculateSphericalWavefrontDepth(v2x, v2y, f);
            double v3z = CalculateSphericalWavefrontDepth(v3x, v3y, f);

            // Calculate the area of the triangle (in aperture plane)
            double area = 0.5 * Math.Abs(v1x * (v2y - v3y) + v2x * (v3y - v1y) + v3x * (v1y - v2y));
            
            if (imagePlane != null)
            {
                for (int i = 0; i < ARRAY_SIZE; i++)
                {
                    for (int j = 0; j < ARRAY_SIZE; j++)
                    {
                        // Calculate image plane pixel coordinates to physical coordinates (mm)
                        double imagePixelSize = imagePlaneSize / ARRAY_SIZE;
                        double imageCenterX = (ARRAY_SIZE - 1) / 2.0;
                        double imageCenterY = (ARRAY_SIZE - 1) / 2.0;
                        
                        double imageX = (i - imageCenterX) * imagePixelSize;
                        double imageY = (j - imageCenterY) * imagePixelSize;
                        
                        // Calculate distances from this image plane pixel to each vertex
                        double distance1 = CalculateDistance3D(imageX, imageY, f, v1x, v1y, v1z);
                        double distance2 = CalculateDistance3D(imageX, imageY, f, v2x, v2y, v2z);
                        double distance3 = CalculateDistance3D(imageX, imageY, f, v3x, v3y, v3z);

                        if( distance1 == distance2 && distance2 == distance3 )
                        {
                            // Normalized intensity
                            imagePlane[i, j] = 1.0f;
                            continue;
                        }

                        // Sort the vertices by distance
                        (double x, double y, double z, double distance)[] vertices = { 
                            (v1x, v1y, v1z, distance1), 
                            (v2x, v2y, v2z, distance2), 
                            (v3x, v3y, v3z, distance3) 
                        };
                        Array.Sort(vertices, (a, b) => a.distance.CompareTo(b.distance));

                        if( vertices[2].distance == vertices[0].distance )
                        {
                            double magnitude1 = CalculateIntegral(vertices, lambda, area);
                            imagePlane[i, j] = (float)(magnitude1 * magnitude1 / area / area); 
                            continue;
                        }

                        if( vertices[1].distance == vertices[0].distance )
                        {
                            (double x, double y, double z, double distance)[] vertices0 = {
                                vertices[1],     
                                vertices[2],
                                vertices[0]
                            };

                            double magnitude1 = CalculateIntegral( vertices0, lambda, area);
                            imagePlane[i, j] = (float)(magnitude1 * magnitude1 / area / area); 
                            continue;
                        }

                        if( vertices[2].distance == vertices[1].distance )
                        {
                            (double x, double y, double z, double distance)[] vertices0 = {
                                vertices[2],     
                                vertices[0],
                                vertices[1]
                            };

                            double magnitude1 = CalculateIntegral( vertices0, lambda, area);
                            imagePlane[i, j] = (float)(magnitude1 * magnitude1 / area / area); 
                            continue;
                        }

                        // Find the point along the line from vertices[0] to vertices[2] where the distance is equal to vertices[1].distance
                        double k = ( vertices[1].distance - vertices[0].distance ) / ( vertices[2].distance - vertices[0].distance );
                        double x = vertices[0].x + k * ( vertices[2].x - vertices[0].x );
                        double y = vertices[0].y + k * ( vertices[2].y - vertices[0].y );
                        double z = vertices[0].z + k * ( vertices[2].z - vertices[0].z );

                        (double x, double y, double z, double distance)[] vertices2 = { 
                            (x, y, z, vertices[1].distance), 
                            vertices[0], 
                            vertices[1] 
                        };

                        double magnitude = CalculateIntegral(vertices2, lambda, area);

                        (double x, double y, double z, double distance)[] vertices3 = { 
                            (x, y, z, vertices[1].distance), 
                            vertices[2], 
                            vertices[1] 
                        };
                        magnitude += CalculateIntegral(vertices3, lambda, area);
                        
                        // Store in image plane
                        imagePlane[i, j] = (float)(magnitude * magnitude / area / area );
                    }
                }
                
                // Update the display
                DisplayImagePlane();
            }
        }

        private double CalculateIntegral( (double x, double y, double z, double distance)[] vertices, double lambda, double area )
        {
            System.Diagnostics.Trace.Assert( vertices[2].distance == vertices[0].distance );
            System.Diagnostics.Trace.Assert( vertices[1].distance != vertices[0].distance );
            // Find the length of triangle side from vertices[0] to vertices[2]
            double side02 = Math.Sqrt(
                Math.Pow(vertices[0].x - vertices[2].x, 2) +
                Math.Pow(vertices[0].y - vertices[2].y, 2) );
            // Altitude to side02
            double alt = 2 * area / side02;
            // Phase difference between vertices[2] and vertices[0]
            double phase0 = 2 * Math.PI * vertices[0].distance / lambda;
            double phase1 = 2 * Math.PI * vertices[1].distance / lambda;
            // Integral of x * sin( x * (phase0 - phase1) / alt + phase1 ) for x in [0, alt]
            // Closed form solution: ∫ x * sin(ax + b) dx = (1/a²) * (sin(ax + b) - ax * cos(ax + b))
            double a = (phase0 - phase1) / alt;
            double b = phase1;
            double magnitude = (1.0 / (a * a)) * (Math.Sin(a * alt + b) - a * alt * Math.Cos(a * alt + b));
            // Subtract the value at x=0: (1/a²) * (sin(b) - 0 * cos(b)) = (1/a²) * sin(b)
            magnitude -= (1.0 / (a * a)) * Math.Sin(b);

            return magnitude;
        }
        
        private double CalculateDistance3D(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            // Calculate 3D Euclidean distance between two points
            double dx = x2 - x1;
            double dy = y2 - y1;
            double dz = z2 - z1;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }
        
        private double CalculateSphericalWavefrontDepth(double x, double y, double f)
        {
            // Same logic as in WavefrontGenerator.GenerateSphericalWavefront
            // Calculate the depth (z) for a point (x, y) on the spherical wavefront
            
            // For a spherical wavefront with focal length f:
            // z = f - sqrt(f^2 - (x^2 + y^2))
            // This gives the depth relative to the focal plane
            
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

