using System.Drawing;
using System.Drawing.Imaging;
using System.Globalization;

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
            
            // Create a denser grid that covers the entire illuminated area
            int gridSize = 12; // Increased from 8 to 12 for better coverage
            
            // Calculate grid spacing based on the actual illuminated area size
            // Use the smaller dimension to ensure we don't have gaps
            double minDimension = Math.Min(illuminatedBounds.Value.width, illuminatedBounds.Value.height);
            double gridSpacing = minDimension / (gridSize - 1); // Ensure edge-to-edge coverage
            
            // Calculate grid starting position to center it in the illuminated area
            double gridStartI = illuminatedBounds.Value.centerI - ((gridSize - 1) / 2.0) * gridSpacing;
            double gridStartJ = illuminatedBounds.Value.centerJ - ((gridSize - 1) / 2.0) * gridSpacing;
            
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
                    }
                    else if (IsPositionIlluminated(gridI_pos, gridJ_pos, arraySize))
                    {
                        // Position is illuminated but on boundary - draw as red
                        DrawDotAtPosition(bitmap, gridI_pos, gridJ_pos, arraySize, Color.Red);
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
                            }
                            else
                            {
                                DrawDotAtPosition(bitmap, nearestPosition.Value.Item1, nearestPosition.Value.Item2, arraySize, Color.Green);
                            }
                        }
                    }
                }
            }
            
            // Add boundary nodes around the illuminated area
            AddBoundaryNodes(bitmap, arraySize, centerApertureI, centerApertureJ, illuminatedBounds.Value);
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
        
        private void AddBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, (double centerI, double centerJ, double width, double height) bounds)
        {
            // Add boundary nodes at regular intervals around the illuminated area
            int numBoundaryNodes = 32; // Increased from 16 to 32 for better coverage
            
            for (int i = 0; i < numBoundaryNodes; i++)
            {
                double angle = (2.0 * Math.PI * i) / numBoundaryNodes;
                
                // Start from the center and move outward to find boundary
                double distance = 0;
                double maxDistance = Math.Max(bounds.width, bounds.height) / 2.0 + 20; // Increased safety margin
                
                while (distance < maxDistance)
                {
                    double testI = centerI + distance * Math.Cos(angle);
                    double testJ = centerJ + distance * Math.Sin(angle);
                    
                    // Check if we've found the boundary
                    if (!IsPositionIlluminated(testI, testJ, arraySize))
                    {
                        // Move back slightly to stay inside illuminated area
                        distance -= 2.0; // Increased step back for safety
                        double boundaryI = centerI + distance * Math.Cos(angle);
                        double boundaryJ = centerJ + distance * Math.Sin(angle);
                        
                        // ROBUST SAFETY CHECK: Ensure we're definitely inside illuminated area
                        if (IsPositionIlluminated(boundaryI, boundaryJ, arraySize))
                        {
                            // Additional verification: check surrounding positions to ensure we're not on the very edge
                            bool isStable = true;
                            for (int checkAngle = -1; checkAngle <= 1; checkAngle++)
                            {
                                for (int checkDist = -1; checkDist <= 1; checkDist++)
                                {
                                    if (checkAngle == 0 && checkDist == 0) continue; // Skip center position
                                    
                                    double checkI = boundaryI + checkDist * Math.Cos(angle + checkAngle * 0.1);
                                    double checkJ = boundaryJ + checkDist * Math.Sin(angle + checkAngle * 0.1);
                                    
                                    if (!IsPositionIlluminated(checkI, checkJ, arraySize))
                                    {
                                        isStable = false;
                                        break;
                                    }
                                }
                                if (!isStable) break;
                            }
                            
                            // Only draw if the position is stable (surrounded by illuminated pixels)
                            if (isStable)
                            {
                                DrawDotAtPosition(bitmap, boundaryI, boundaryJ, arraySize, Color.Red);
                            }
                        }
                        break;
                    }
                    
                    distance += 1.0;
                }
            }
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
    }
}
