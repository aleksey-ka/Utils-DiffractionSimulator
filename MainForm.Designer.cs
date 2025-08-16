namespace DiffractionSimulator
{
    partial class MainForm
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
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
            apertureShape_label = new Label();
            apertureShapeComboBox = new ComboBox();
            apertureTitle_label = new Label();
            imagePlaneTitle_label = new Label();
            obstructionRatio_label = new Label();
            obstructionRatioNumericUpDown = new NumericUpDown();
            ((System.ComponentModel.ISupportInitialize)apperturePictureBox).BeginInit();
            ((System.ComponentModel.ISupportInitialize)imagePlanePictureBox).BeginInit();
            ((System.ComponentModel.ISupportInitialize)D_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)f_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)lambda_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)imagePlaneSize_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)integrationSize_numericUpDown).BeginInit();
            ((System.ComponentModel.ISupportInitialize)obstructionRatioNumericUpDown).BeginInit();
            SuspendLayout();
            // 
            // apperturePictureBox
            // 
            apperturePictureBox.BackColor = Color.Black;
            apperturePictureBox.BorderStyle = BorderStyle.FixedSingle;
            apperturePictureBox.Location = new Point(50, 50);
            apperturePictureBox.Name = "apperturePictureBox";
            apperturePictureBox.Size = new Size(503, 503);
            apperturePictureBox.TabIndex = 0;
            apperturePictureBox.TabStop = false;
            // 
            // imagePlanePictureBox
            // 
            imagePlanePictureBox.BackColor = Color.Black;
            imagePlanePictureBox.BorderStyle = BorderStyle.FixedSingle;
            imagePlanePictureBox.Location = new Point(600, 50);
            imagePlanePictureBox.Name = "imagePlanePictureBox";
            imagePlanePictureBox.Size = new Size(503, 503);
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
            // apertureShape_label
            // 
            apertureShape_label.AutoSize = true;
            apertureShape_label.ForeColor = Color.White;
            apertureShape_label.Location = new Point(450, 562);
            apertureShape_label.Name = "apertureShape_label";
            apertureShape_label.Size = new Size(88, 15);
            apertureShape_label.TabIndex = 14;
            apertureShape_label.Text = "Aperture Shape";
            apertureShape_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // apertureShapeComboBox
            // 
            apertureShapeComboBox.BackColor = Color.FromArgb(64, 64, 64);
            apertureShapeComboBox.DropDownStyle = ComboBoxStyle.DropDownList;
            apertureShapeComboBox.ForeColor = Color.White;
            apertureShapeComboBox.FormattingEnabled = true;
            apertureShapeComboBox.Items.AddRange(new object[] { "Circular", "Square", "Circular with obstr." });
            apertureShapeComboBox.Location = new Point(450, 582);
            apertureShapeComboBox.Name = "apertureShapeComboBox";
            apertureShapeComboBox.Size = new Size(100, 23);
            apertureShapeComboBox.TabIndex = 15;
            apertureShapeComboBox.SelectedIndexChanged += ApertureShapeComboBox_SelectedIndexChanged;
            // 
            // apertureTitle_label
            // 
            apertureTitle_label.AutoSize = true;
            apertureTitle_label.Font = new Font("Segoe UI", 12F, FontStyle.Bold);
            apertureTitle_label.ForeColor = Color.White;
            apertureTitle_label.Location = new Point(50, 25);
            apertureTitle_label.Name = "apertureTitle_label";
            apertureTitle_label.Size = new Size(77, 21);
            apertureTitle_label.TabIndex = 16;
            apertureTitle_label.Text = "Aperture";
            // 
            // imagePlaneTitle_label
            // 
            imagePlaneTitle_label.AutoSize = true;
            imagePlaneTitle_label.Font = new Font("Segoe UI", 12F, FontStyle.Bold);
            imagePlaneTitle_label.ForeColor = Color.White;
            imagePlaneTitle_label.Location = new Point(600, 25);
            imagePlaneTitle_label.Name = "imagePlaneTitle_label";
            imagePlaneTitle_label.Size = new Size(105, 21);
            imagePlaneTitle_label.TabIndex = 17;
            imagePlaneTitle_label.Text = "Image Plane";
            // 
            // obstructionRatio_label
            // 
            obstructionRatio_label.AutoSize = true;
            obstructionRatio_label.ForeColor = Color.White;
            obstructionRatio_label.Location = new Point(451, 604);
            obstructionRatio_label.Name = "obstructionRatio_label";
            obstructionRatio_label.Size = new Size(83, 15);
            obstructionRatio_label.TabIndex = 18;
            obstructionRatio_label.Text = "Obstruction %";
            obstructionRatio_label.TextAlign = ContentAlignment.MiddleCenter;
            // 
            // obstructionRatioNumericUpDown
            // 
            obstructionRatioNumericUpDown.DecimalPlaces = 1;
            obstructionRatioNumericUpDown.Location = new Point(451, 624);
            obstructionRatioNumericUpDown.Name = "obstructionRatioNumericUpDown";
            obstructionRatioNumericUpDown.Size = new Size(80, 23);
            obstructionRatioNumericUpDown.TabIndex = 19;
            obstructionRatioNumericUpDown.Value = new decimal(new int[] { 300, 0, 0, 65536 });
            obstructionRatioNumericUpDown.ValueChanged += ObstructionRatioNumericUpDown_ValueChanged;
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
            Controls.Add(apertureShape_label);
            Controls.Add(apertureShapeComboBox);
            Controls.Add(apertureTitle_label);
            Controls.Add(imagePlaneTitle_label);
            Controls.Add(obstructionRatio_label);
            Controls.Add(obstructionRatioNumericUpDown);
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
            ((System.ComponentModel.ISupportInitialize)obstructionRatioNumericUpDown).EndInit();
            ResumeLayout(false);
            PerformLayout();
        }

        #endregion

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
        private Label apertureShape_label;
        private ComboBox apertureShapeComboBox;
        private Label apertureTitle_label;
        private Label imagePlaneTitle_label;
        private Label obstructionRatio_label;
        private NumericUpDown obstructionRatioNumericUpDown;
    }
}
