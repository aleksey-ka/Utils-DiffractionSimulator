# Diffraction Simulator

A C# Windows Forms application targeting .NET 8 that displays a 501x501 array of random float values as a grayscale image.

## Features

- Dark theme UI
- 501x501 picture box displaying grayscale image
- Random float values from 0.0 to 1.0
- No image scaling - each float represents a single pixel

## Requirements

- .NET 8 SDK
- Windows operating system

## Building and Running

1. Open a terminal in the project directory
2. Build the project:
   ```
   dotnet build
   ```
3. Run the application:
   ```
   dotnet run
   ```

## Project Structure

- `Program.cs` - Application entry point
- `MainForm.cs` - Main form with the picture box
- `DiffractionSimulator.csproj` - Project file

The application creates a 501x501 array of random float values and displays them as a grayscale image where each float value (0.0-1.0) is converted to a grayscale pixel value (0-255).
