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

            apertureShapeComboBox.SelectedIndex = 0;
            UpdateObstructionControlState();
        }

        private void InitializeDataArray()
        {
            appertureMask = new float[ARRAY_SIZE, ARRAY_SIZE];
            appertureDepth = new double[ARRAY_SIZE, ARRAY_SIZE];
            imagePlane = new float[ARRAY_SIZE, ARRAY_SIZE];
            
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
            string selectedShape = apertureShapeComboBox.SelectedItem?.ToString() ?? "Circular";
            double D_value = (double)D_numericUpDown.Value;
            double obstructionRatio = (double)obstructionRatioNumericUpDown.Value;
            
            // Generate aperture mask using the new class
            appertureMask = ApertureMaskGenerator.GenerateApertureMask(selectedShape, D_value, obstructionRatio);
            
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
            
            // Generate wavefront using the new class
            appertureDepth = WavefrontGenerator.GenerateSphericalWavefront(D_value, f_value);
            
            // Calculate diffraction pattern using the new class
            imagePlane = DiffractionCalculator.CalculateDiffractionPattern(
                appertureMask, 
                appertureDepth, 
                D_value, 
                f_value, 
                lambda_value, 
                imagePlaneSize_value, 
                integrationSize);
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
        
        private void CalculateButton_Click(object sender, EventArgs e)
        {
            // Calculate diffraction pattern with current parameters
            CalculateDiffractionPattern();
            DisplayImagePlane();
        }
    }
}
