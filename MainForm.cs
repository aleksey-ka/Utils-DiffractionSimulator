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
            
            // Calculate the center of the aperture (exact center)
            double centerApertureI = (arraySize - 1) / 2.0;
            double centerApertureJ = (arraySize - 1) / 2.0;
            
            // Find the illuminated area bounds to determine optimal grid spacing
            var illuminatedBounds = FindIlluminatedBounds(arraySize);
            if (illuminatedBounds == null) return;
            
            // Create a grid with a node exactly at the center of the aperture
            int gridSize = 13; // Use odd number to ensure center node (changed from 12 to 13)
            
            // Calculate grid spacing based on the actual illuminated area size
            // Use the smaller dimension to ensure we don't have gaps
            double minDimension = Math.Min(illuminatedBounds.Value.width, illuminatedBounds.Value.height);
            double gridSpacing = minDimension / (gridSize - 1); // Ensure edge-to-edge coverage
            
            // Calculate grid starting position to place center node exactly at aperture center
            // With odd gridSize, the center node will be at index (gridSize-1)/2
            int centerIndex = (gridSize - 1) / 2;
            double gridStartI = centerApertureI - centerIndex * gridSpacing;
            double gridStartJ = centerApertureJ - centerIndex * gridSpacing;
            
            // Collect all nodes for triangulation (including boundary nodes)
            List<(double i, double j, bool isBoundary)> allNodes = new List<(double i, double j, bool isBoundary)>();
            
            // Draw regular grid dots (green) with better coverage
            for (int gridI = 0; gridI < gridSize; gridI++)
            {
                for (int gridJ = 0; gridJ < gridSize; gridJ++)
                {
                    double gridI_pos = gridStartI + gridI * gridSpacing;
                    double gridJ_pos = gridStartJ + gridJ * gridSpacing;
                    
                    // Only draw if this position is illuminated AND not on the boundary
                    if (IsPositionIlluminated(gridI_pos, gridJ_pos, arraySize) && !IsPositionOnBoundary(gridI_pos, gridJ_pos, arraySize))
                    {
                        DrawDotAtPosition(bitmap, gridI_pos, gridJ_pos, arraySize, Color.Green);
                        allNodes.Add((gridI_pos, gridJ_pos, false));
                    }
                    else if (IsPositionIlluminated(gridI_pos, gridJ_pos, arraySize))
                    {
                        // Position is illuminated but on boundary - draw as red
                        DrawDotAtPosition(bitmap, gridI_pos, gridJ_pos, arraySize, Color.Red);
                        allNodes.Add((gridI_pos, gridJ_pos, true)); // Add boundary nodes to triangulation
                    }
                    else
                    {
                        // If the grid position is not illuminated, try to find the nearest illuminated position
                        var nearestPosition = FindNearestIlluminatedPosition(gridI_pos, gridJ_pos, arraySize);
                        if (nearestPosition.HasValue)
                        {
                            // Check if the nearest position is on the boundary
                            if (IsPositionOnBoundary(nearestPosition.Value.Item1, nearestPosition.Value.Item2, arraySize))
                            {
                                DrawDotAtPosition(bitmap, nearestPosition.Value.Item1, nearestPosition.Value.Item2, arraySize, Color.Red);
                                allNodes.Add((nearestPosition.Value.Item1, nearestPosition.Value.Item2, true)); // Add boundary nodes to triangulation
                            }
                            else
                            {
                                DrawDotAtPosition(bitmap, nearestPosition.Value.Item1, nearestPosition.Value.Item2, arraySize, Color.Green);
                                allNodes.Add((nearestPosition.Value.Item1, nearestPosition.Value.Item2, false));
                            }
                        }
                    }
                }
            }
            
            // Add additional boundary nodes around the illuminated area and include them in triangulation
            var additionalBoundaryNodes = AddBoundaryNodes(bitmap, arraySize, centerApertureI, centerApertureJ, illuminatedBounds.Value);
            allNodes.AddRange(additionalBoundaryNodes);
            
            // Generate and draw the triangular mesh using all nodes (including boundary nodes)
            if (allNodes.Count >= 3)
            {
                GenerateAndDrawTriangularMesh(bitmap, arraySize, allNodes);
            }
        }
        
        private (double centerI, double centerJ, double width, double height)? FindIlluminatedBounds(int arraySize)
        {
            if (appertureMask == null) return null;
            
            int minI = arraySize, maxI = 0, minJ = arraySize, maxJ = 0;
            bool foundIlluminated = false;
            
            // Find the bounding box of illuminated area
            for (int i = 0; i < arraySize; i++)
            {
                for (int j = 0; j < arraySize; j++)
                {
                    if (appertureMask[i, j] > 0.001)
                    {
                        minI = Math.Min(minI, i);
                        maxI = Math.Max(maxI, i);
                        minJ = Math.Min(minJ, j);
                        maxJ = Math.Max(maxJ, j);
                        foundIlluminated = true;
                    }
                }
            }
            
            if (!foundIlluminated) return null;
            
            double centerI = (minI + maxI) / 2.0;
            double centerJ = (minJ + maxJ) / 2.0;
            double width = maxI - minI;
            double height = maxJ - minJ;
            
            return (centerI, centerJ, width, height);
        }
        
        private (double, double)? FindNearestIlluminatedPosition(double targetI, double targetJ, int arraySize)
        {
            if (appertureMask == null) return null;
            
            // Search in expanding circles around the target position
            double searchRadius = 0;
            double maxSearchRadius = 20; // Maximum search radius
            
            while (searchRadius <= maxSearchRadius)
            {
                // Check positions at current radius
                for (int angleStep = 0; angleStep < 8; angleStep++) // 8 directions
                {
                    double angle = (2.0 * Math.PI * angleStep) / 8.0;
                    double testI = targetI + searchRadius * Math.Cos(angle);
                    double testJ = targetJ + searchRadius * Math.Sin(angle);
                    
                    if (IsPositionIlluminated(testI, testJ, arraySize))
                    {
                        return (testI, testJ);
                    }
                }
                
                searchRadius += 1.0;
            }
            
            return null; // No illuminated position found within search radius
        }
        
        private List<(double i, double j, bool isBoundary)> AddBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, (double centerI, double centerJ, double width, double height) bounds)
        {
            List<(double i, double j, bool isBoundary)> boundaryNodes = new List<(double i, double j, bool isBoundary)>();
            
            // Add boundary nodes at extremely dense intervals for smooth boundary approximation
            int numBoundaryNodes = 360; // Very dense to eliminate any straight edge artifacts
            
            for (int i = 0; i < numBoundaryNodes; i++)
            {
                double angle = (2.0 * Math.PI * i) / numBoundaryNodes;
                
                // Find the precise boundary position using binary search
                var boundaryPosition = FindPreciseBoundaryPosition(centerI, centerJ, angle, arraySize);
                
                if (boundaryPosition.HasValue)
                {
                    DrawDotAtPosition(bitmap, boundaryPosition.Value.Item1, boundaryPosition.Value.Item2, arraySize, Color.Red);
                    boundaryNodes.Add((boundaryPosition.Value.Item1, boundaryPosition.Value.Item2, true));
                }
            }
            
            // Add adaptive boundary nodes to fill any remaining gaps
            AddImprovedAdaptiveBoundaryNodes(bitmap, arraySize, centerI, centerJ, boundaryNodes);
            
            return boundaryNodes;
        }
        
        private void AddImprovedAdaptiveBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, List<(double i, double j, bool isBoundary)> boundaryNodes)
        {
            // Find gaps in the boundary representation and add nodes to fill them
            double minGapDistance = 1.0; // Very small threshold for extremely dense boundary
            double maxGapDistance = 2.0; // Very strict maximum distance for perfectly smooth curves
            
            var additionalNodes = new List<(double i, double j, bool isBoundary)>();
            
            for (int i = 0; i < boundaryNodes.Count; i++)
            {
                var currentNode = boundaryNodes[i];
                var nextNode = boundaryNodes[(i + 1) % boundaryNodes.Count]; // Wrap around
                
                // Calculate distance between consecutive boundary nodes
                double distance = Math.Sqrt(
                    Math.Pow(nextNode.i - currentNode.i, 2) + 
                    Math.Pow(nextNode.j - currentNode.j, 2)
                );
                
                // Only add nodes if gap is significant
                if (distance > maxGapDistance)
                {
                    // Calculate how many intermediate nodes we need for smooth approximation
                    int numIntermediateNodes = Math.Min(2, (int)(distance / minGapDistance));
                    
                    for (int j = 1; j <= numIntermediateNodes; j++)
                    {
                        double ratio = (double)j / (numIntermediateNodes + 1);
                        
                        // Interpolate position between current and next node
                        double interpI = currentNode.i + ratio * (nextNode.i - currentNode.i);
                        double interpJ = currentNode.j + ratio * (nextNode.j - currentNode.j);
                        
                        // Find the precise boundary position for this interpolated direction
                        double angle = Math.Atan2(interpJ - centerJ, interpI - centerI);
                        var precisePosition = FindPreciseBoundaryPosition(centerI, centerJ, angle, arraySize);
                        
                        if (precisePosition.HasValue)
                        {
                            // Check if this position is sufficiently far from existing nodes
                            bool tooClose = false;
                            foreach (var existing in boundaryNodes)
                            {
                                double distToExisting = Math.Sqrt(
                                    Math.Pow(precisePosition.Value.Item1 - existing.i, 2) + 
                                    Math.Pow(precisePosition.Value.Item2 - existing.j, 2)
                                );
                                if (distToExisting < minGapDistance * 0.5) // Very tight proximity check
                                {
                                    tooClose = true;
                                    break;
                                }
                            }
                            
                            if (!tooClose)
                            {
                                DrawDotAtPosition(bitmap, precisePosition.Value.Item1, precisePosition.Value.Item2, arraySize, Color.Red);
                                additionalNodes.Add((precisePosition.Value.Item1, precisePosition.Value.Item2, true));
                            }
                        }
                    }
                }
            }
            
            // Add the additional nodes to the main list
            boundaryNodes.AddRange(additionalNodes);
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
    }
}
