using System.Drawing;
using System.Drawing.Imaging;

namespace DiffractionSimulator
{
    public partial class MainForm : Form
    {
        private float[,] appertureMask;
        private double[,] appertureDepth;
        private float[,] imagePlane;
        private const int ARRAY_SIZE = 501;
        
        // Physical parameters are now controlled by UI controls

        public MainForm()
        {
            InitializeComponent();
            InitializeDataArray();
            DisplayArrayAsImage();
            DisplayImagePlane();
        }

        private void InitializeDataArray()
        {
            appertureMask = new float[ARRAY_SIZE, ARRAY_SIZE];
            appertureDepth = new double[ARRAY_SIZE, ARRAY_SIZE];
            imagePlane = new float[ARRAY_SIZE, ARRAY_SIZE];            
            
            // Ensure the dropdown has a default selection
            if (apertureShapeComboBox.SelectedIndex == -1)
            {
                apertureShapeComboBox.SelectedIndex = 0; // Default to "Circular"
            }
            
            // Initialize with default circular aperture
            UpdateApertureMask();
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    appertureDepth[i, j] = 0.0f;
                    imagePlane[i, j] = 0.0f;
                }
            }
        }

        private void UpdateApertureMask()
        {
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Square";
            
            if (selectedShape == "Circular")
            {
                // Create circular aperture
                double D_value = (double)D_numericUpDown.Value;
                double aperturePixelSize = D_value / ARRAY_SIZE; // mm per pixel in aperture
                double radius = D_value / 2.0; // radius in mm
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
                            appertureMask[i, j] = 1.0f;
                        }
                        else
                        {
                            appertureMask[i, j] = 0.0f;
                        }
                    }
                }
            }
            else
            {
                // Default square aperture (uniform)
                for (int i = 0; i < ARRAY_SIZE; i++)
                {
                    for (int j = 0; j < ARRAY_SIZE; j++)
                    {
                        appertureMask[i, j] = 1.0f;
                    }
                }
            }
        }

        private void DisplayArrayAsImage()
        {
            Bitmap bitmap = new Bitmap(ARRAY_SIZE, ARRAY_SIZE);
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    int grayValue = (int)(appertureMask[i, j] * 255);
                    Color pixelColor = Color.FromArgb(grayValue, grayValue, grayValue);
                    bitmap.SetPixel(j, i, pixelColor);
                }
            }
            
            apperturePictureBox.Image = bitmap;
        }
        
        private void DisplayImagePlane()
        {
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
            // Get values from controls
            double D_value = (double)D_numericUpDown.Value;
            double f_value = (double)f_numericUpDown.Value;
            double lambda_value = (double)lambda_numericUpDown.Value;
            double imagePlaneSize_value = (double)imagePlaneSize_numericUpDown.Value;
            int integrationSize = (int)integrationSize_numericUpDown.Value;
            
            // Calculate physical pixel sizes
            double aperturePixelSize = D_value / ARRAY_SIZE; // mm per pixel in aperture
            double imagePixelSize = imagePlaneSize_value / ARRAY_SIZE; // mm per pixel in image plane

            // Center of curvature is at z = f
            // For a sphere centered at (0, 0, f), the equation is: x² + y² + (z-f)² = f²
            // Solving for z: z = f - sqrt(f² - x² - y²)
            
            // Calculate wavefront depth for all aperture points
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
                    double z = f_value - Math.Sqrt(f_value * f_value - rSquared);
                    appertureDepth[i, j] = z;
                }
            }
            
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
            
            // For each point in the central 11x11 image plane region
            for (int imgI = centerStart; imgI <= centerEnd; imgI++)
            {
                for (int imgJ = centerStart; imgJ <= centerEnd; imgJ++)
                {
                    // Convert image plane pixel coordinates to physical coordinates (mm)
                    double imgX = (imgI - ARRAY_SIZE / 2.0) * imagePixelSize;
                    double imgY = (imgJ - ARRAY_SIZE / 2.0) * imagePixelSize;
                    double imgZ = f_value; // Image plane is at z = f
                    
                    // Complex field at this image point (real and imaginary parts)
                    double realPart = 0.0;
                    double imagPart = 0.0;
                    
                    // Integrate over all aperture points
                    for (int apI = 0; apI < ARRAY_SIZE; apI++)
                    {
                        for (int apJ = 0; apJ < ARRAY_SIZE; apJ++)
                        {
                            // Convert aperture pixel coordinates to physical coordinates (mm)
                            double apX = (apI - ARRAY_SIZE / 2.0) * aperturePixelSize;
                            double apY = (apJ - ARRAY_SIZE / 2.0) * aperturePixelSize;
                            double apZ = appertureDepth[apI, apJ]; // Z position from wavefront
                            
                            // Calculate distance from aperture point to image point
                            double dx = imgX - apX;
                            double dy = imgY - apY;
                            double dz = imgZ - apZ;
                            double distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                            
                            // Calculate phase: 2π * distance / λ
                            double phase = 2.0 * Math.PI * distance / lambda_value;
                            
                            // Get amplitude from aperture mask
                            float amplitude = appertureMask[apI, apJ];
                            
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
        
        private void CalculateButton_Click(object sender, EventArgs e)
        {
            // Calculate diffraction pattern with current parameters
            CalculateDiffractionPattern();
            DisplayImagePlane();
        }
    }
}
