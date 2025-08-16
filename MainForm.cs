using System.Drawing;
using System.Drawing.Imaging;

namespace DiffractionSimulator
{
    public partial class MainForm : Form
    {
        private PictureBox apperturePictureBox;
        private PictureBox imagePlanePictureBox;
        private CheckBox gammaCorrectedViewCheckBox;
        private NumericUpDown D_numericUpDown;
        private NumericUpDown f_numericUpDown;
        private NumericUpDown lambda_numericUpDown;
        private NumericUpDown imagePlaneSize_numericUpDown;
        private NumericUpDown integrationSize_numericUpDown;
        private Button calculateButton;
        private Label D_label;
        private Label f_label;
        private Label lambda_label;
        private Label imagePlaneSize_label;
        private Label integrationSize_label;
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

        private void InitializeComponent()
        {
            apperturePictureBox = new PictureBox();
            imagePlanePictureBox = new PictureBox();
            gammaCorrectedViewCheckBox = new CheckBox();
            D_numericUpDown = new NumericUpDown();
            f_numericUpDown = new NumericUpDown();
            lambda_numericUpDown = new NumericUpDown();
            imagePlaneSize_numericUpDown = new NumericUpDown();
            integrationSize_numericUpDown = new NumericUpDown();
            calculateButton = new Button();
            D_label = new Label();
            f_label = new Label();
            lambda_label = new Label();
            imagePlaneSize_label = new Label();
            integrationSize_label = new Label();
            ((System.ComponentModel.ISupportInitialize)apperturePictureBox).BeginInit();
            ((System.ComponentModel.ISupportInitialize)imagePlanePictureBox).BeginInit();
            ((System.ComponentModel.ISupportInitialize)D_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)f_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)lambda_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)imagePlaneSize_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)integrationSize_numericUpDown).BeginInit();
            SuspendLayout();
            // 
            // apperturePictureBox
            // 
            apperturePictureBox.BackColor = Color.Black;
            apperturePictureBox.BorderStyle = BorderStyle.FixedSingle;
            apperturePictureBox.Location = new Point(50, 50);
            apperturePictureBox.Name = "apperturePictureBox";
            apperturePictureBox.Size = new Size(501, 501);
            apperturePictureBox.TabIndex = 0;
            apperturePictureBox.TabStop = false;
            // 
            // imagePlanePictureBox
            // 
            imagePlanePictureBox.BackColor = Color.Black;
            imagePlanePictureBox.BorderStyle = BorderStyle.FixedSingle;
            imagePlanePictureBox.Location = new Point(600, 50);
            imagePlanePictureBox.Name = "imagePlanePictureBox";
            imagePlanePictureBox.Size = new Size(501, 501);
            imagePlanePictureBox.TabIndex = 1;
            imagePlanePictureBox.TabStop = false;
            // 
            // gammaCorrectedViewCheckBox
            // 
            gammaCorrectedViewCheckBox.AutoSize = true;
            gammaCorrectedViewCheckBox.ForeColor = Color.White;
            gammaCorrectedViewCheckBox.Location = new Point(968, 571);
            gammaCorrectedViewCheckBox.Name = "gammaCorrectedViewCheckBox";
            gammaCorrectedViewCheckBox.Size = new Size(123, 19);
            gammaCorrectedViewCheckBox.TabIndex = 2;
            gammaCorrectedViewCheckBox.Text = "Gamma Corrected";
            gammaCorrectedViewCheckBox.UseVisualStyleBackColor = false;
            gammaCorrectedViewCheckBox.CheckedChanged += GammaCorrectedViewCheckBox_CheckedChanged;
            // 
            // D_numericUpDown
            // 
            D_numericUpDown.DecimalPlaces = 1;
            D_numericUpDown.Location = new Point(50, 582);
            D_numericUpDown.Maximum = new decimal(new int[] { 10000, 0, 0, 65536 });
            D_numericUpDown.Minimum = new decimal(new int[] { 10, 0, 0, 65536 });
            D_numericUpDown.Name = "D_numericUpDown";
            D_numericUpDown.Size = new Size(80, 23);
            D_numericUpDown.TabIndex = 3;
            D_numericUpDown.Value = new decimal(new int[] { 1000, 0, 0, 65536 });
            // 
            // f_numericUpDown
            // 
            f_numericUpDown.Increment = new decimal(new int[] { 100, 0, 0, 65536 });
            f_numericUpDown.Location = new Point(148, 582);
            f_numericUpDown.Maximum = new decimal(new int[] { 100000, 0, 0, 65536 });
            f_numericUpDown.Minimum = new decimal(new int[] { 1000, 0, 0, 65536 });
            f_numericUpDown.Name = "f_numericUpDown";
            f_numericUpDown.Size = new Size(80, 23);
            f_numericUpDown.TabIndex = 4;
            f_numericUpDown.Value = new decimal(new int[] { 10000, 0, 0, 65536 });
            // 
            // lambda_numericUpDown
            // 
            lambda_numericUpDown.DecimalPlaces = 4;
            lambda_numericUpDown.Increment = new decimal(new int[] { 1, 0, 0, 262144 });
            lambda_numericUpDown.Location = new Point(247, 582);
            lambda_numericUpDown.Maximum = new decimal(new int[] { 1, 0, 0, 131072 });
            lambda_numericUpDown.Minimum = new decimal(new int[] { 1, 0, 0, 262144 });
            lambda_numericUpDown.Name = "lambda_numericUpDown";
            lambda_numericUpDown.Size = new Size(80, 23);
            lambda_numericUpDown.TabIndex = 5;
            lambda_numericUpDown.Value = new decimal(new int[] { 5, 0, 0, 262144 });
            // 
            // imagePlaneSize_numericUpDown
            // 
            imagePlaneSize_numericUpDown.DecimalPlaces = 2;
            imagePlaneSize_numericUpDown.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            imagePlaneSize_numericUpDown.Location = new Point(601, 582);
            imagePlaneSize_numericUpDown.Maximum = new decimal(new int[] { 100, 0, 0, 65536 });
            imagePlaneSize_numericUpDown.Minimum = new decimal(new int[] { 1, 0, 0, 65536 });
            imagePlaneSize_numericUpDown.Name = "imagePlaneSize_numericUpDown";
            imagePlaneSize_numericUpDown.Size = new Size(80, 23);
            imagePlaneSize_numericUpDown.TabIndex = 6;
            imagePlaneSize_numericUpDown.Value = new decimal(new int[] { 5, 0, 0, 65536 });
            // 
            // integrationSize_numericUpDown
            // 
            integrationSize_numericUpDown.Location = new Point(699, 582);
            integrationSize_numericUpDown.Minimum = new decimal(new int[] { 10, 0, 0, 65536 });
            integrationSize_numericUpDown.Name = "integrationSize_numericUpDown";
            integrationSize_numericUpDown.Size = new Size(80, 23);
            integrationSize_numericUpDown.TabIndex = 7;
            integrationSize_numericUpDown.Value = new decimal(new int[] { 400, 0, 0, 65536 });
            // 
            // calculateButton
            // 
            calculateButton.BackColor = Color.FromArgb(64, 64, 64);
            calculateButton.FlatStyle = FlatStyle.Flat;
            calculateButton.ForeColor = Color.White;
            calculateButton.Location = new Point(1001, 596);
            calculateButton.Name = "calculateButton";
            calculateButton.Size = new Size(100, 39);
            calculateButton.TabIndex = 8;
            calculateButton.Text = "Calculate";
            calculateButton.UseVisualStyleBackColor = false;
            calculateButton.Click += CalculateButton_Click;
            // 
            // D_label
            // 
            D_label.AutoSize = true;
            D_label.ForeColor = Color.White;
            D_label.Location = new Point(50, 562);
            D_label.Name = "D_label";
            D_label.Size = new Size(48, 15);
            D_label.TabIndex = 9;
            D_label.Text = "D (mm)";
            D_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // f_label
            // 
            f_label.AutoSize = true;
            f_label.ForeColor = Color.White;
            f_label.Location = new Point(148, 562);
            f_label.Name = "f_label";
            f_label.Size = new Size(44, 15);
            f_label.TabIndex = 10;
            f_label.Text = "f (mm)";
            f_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // lambda_label
            // 
            lambda_label.AutoSize = true;
            lambda_label.ForeColor = Color.White;
            lambda_label.Location = new Point(247, 562);
            lambda_label.Name = "lambda_label";
            lambda_label.Size = new Size(46, 15);
            lambda_label.TabIndex = 11;
            lambda_label.Text = "λ (mm)";
            lambda_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // imagePlaneSize_label
            // 
            imagePlaneSize_label.AutoSize = true;
            imagePlaneSize_label.ForeColor = Color.White;
            imagePlaneSize_label.Location = new Point(601, 562);
            imagePlaneSize_label.Name = "imagePlaneSize_label";
            imagePlaneSize_label.Size = new Size(73, 15);
            imagePlaneSize_label.TabIndex = 12;
            imagePlaneSize_label.Text = "Image (mm)";
            imagePlaneSize_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // integrationSize_label
            // 
            integrationSize_label.AutoSize = true;
            integrationSize_label.ForeColor = Color.White;
            integrationSize_label.Location = new Point(699, 562);
            integrationSize_label.Name = "integrationSize_label";
            integrationSize_label.Size = new Size(47, 15);
            integrationSize_label.TabIndex = 13;
            integrationSize_label.Text = "Int. Size";
            integrationSize_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // MainForm
            // 
            AutoScaleDimensions = new SizeF(7F, 15F);
            AutoScaleMode = AutoScaleMode.Font;
            BackColor = Color.FromArgb(32, 32, 32);
            ClientSize = new Size(1162, 700);
            Controls.Add(apperturePictureBox);
            Controls.Add(imagePlanePictureBox);
            Controls.Add(gammaCorrectedViewCheckBox);
            Controls.Add(D_numericUpDown);
            Controls.Add(f_numericUpDown);
            Controls.Add(lambda_numericUpDown);
            Controls.Add(imagePlaneSize_numericUpDown);
            Controls.Add(integrationSize_numericUpDown);
            Controls.Add(calculateButton);
            Controls.Add(D_label);
            Controls.Add(f_label);
            Controls.Add(lambda_label);
            Controls.Add(imagePlaneSize_label);
            Controls.Add(integrationSize_label);
            ForeColor = Color.White;
            Name = "MainForm";
            Text = "Diffraction Simulator";
            ((System.ComponentModel.ISupportInitialize)apperturePictureBox).EndInit();
            ((System.ComponentModel.ISupportInitialize)imagePlanePictureBox).EndInit();
            ((System.ComponentModel.ISupportInitialize)D_numericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)f_numericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)lambda_numericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)imagePlaneSize_numericUpDown).EndInit();
            ((System.ComponentModel.ISupportInitialize)integrationSize_numericUpDown).EndInit();
            ResumeLayout(false);
            PerformLayout();
        }

        private void InitializeDataArray()
        {
            appertureMask = new float[ARRAY_SIZE, ARRAY_SIZE];
            appertureDepth = new double[ARRAY_SIZE, ARRAY_SIZE];
            imagePlane = new float[ARRAY_SIZE, ARRAY_SIZE];            
            
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                for (int j = 0; j < ARRAY_SIZE; j++)
                {
                    appertureMask[i, j] = 1.0f;
                    appertureDepth[i, j] = 0.0f;
                    imagePlane[i, j] = 0.0f;
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
        
        private void CalculateButton_Click(object sender, EventArgs e)
        {
            // Calculate diffraction pattern with current parameters
            CalculateDiffractionPattern();
            DisplayImagePlane();
        }
    }
}
