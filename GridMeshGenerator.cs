using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Windows.Forms;

namespace DiffractionSimulator
{
    public class GridMeshGenerator
    {
        private const int ARRAY_SIZE = 501;

        /// <summary>
        /// Generates and draws an orthogonal grid with boundary nodes and triangular mesh
        /// </summary>
        /// <param name="bitmap">Bitmap to draw on</param>
        /// <param name="apertureShape">Shape of the aperture</param>
        /// <param name="D_value">Aperture diameter in mm</param>
        /// <param name="obstructionRatio">Obstruction ratio for annular apertures</param>
        public static void GenerateOrthogonalGrid(Bitmap bitmap, string apertureShape, double D_value, double obstructionRatio)
        {
            // Use fixed size for grid generation - sampling should not affect grid appearance
            int fixedArraySize = ARRAY_SIZE; // Always use the display resolution for grid
            
            // Calculate the center of the aperture (exact center)
            double centerApertureI = (fixedArraySize - 1) / 2.0;
            double centerApertureJ = (fixedArraySize - 1) / 2.0;
            
            // Find the illuminated area bounds to determine optimal grid spacing
            var illuminatedBounds = FindIlluminatedBounds(fixedArraySize, apertureShape, D_value, obstructionRatio);
            if (illuminatedBounds == null) return;
            
            // Create a grid with a node exactly at the center of the aperture
            int gridSize = 37; // Use odd number to ensure center node (changed from 12 to 13)
            
            // Calculate grid spacing based on the actual illuminated area size
            double gridSpacing = Math.Min(illuminatedBounds.Value.width, illuminatedBounds.Value.height) / (gridSize - 1);
            
            // Calculate grid start positions to center the grid on the aperture center
            double gridStartI = centerApertureI - (gridSize - 1) * gridSpacing / 2.0;
            double gridStartJ = centerApertureJ - (gridSize - 1) * gridSpacing / 2.0;
            
            List<(double i, double j, bool isBoundary)> allNodes = new List<(double, double, bool)>();
            
            // Generate regular grid nodes inside the illuminated area
            for (int gridI = 0; gridI < gridSize; gridI++)
            {
                for (int gridJ = 0; gridJ < gridSize; gridJ++)
                {
                    double gridI_pos = gridStartI + gridI * gridSpacing;
                    double gridJ_pos = gridStartJ + gridJ * gridSpacing;
                    
                    // Only draw if this exact grid position is illuminated (no fallback searching)
                    if (IsPositionIlluminatedDirect(gridI_pos, gridJ_pos, fixedArraySize, apertureShape, D_value, obstructionRatio))
                    {
                        DrawDotAtPosition(bitmap, gridI_pos, gridJ_pos, fixedArraySize, Color.Green);
                        allNodes.Add((gridI_pos, gridJ_pos, false));
                    }
                    // If position is not illuminated, skip it entirely - no fallback searching
                }
            }
            
            // Add additional boundary nodes around the illuminated area and include them in triangulation
            var additionalBoundaryNodes = AddBoundaryNodes(bitmap, fixedArraySize, centerApertureI, centerApertureJ, illuminatedBounds.Value, apertureShape, D_value, obstructionRatio);
            allNodes.AddRange(additionalBoundaryNodes);
            
            // Generate and draw the triangular mesh using all nodes (including boundary nodes)
            if (allNodes.Count >= 3)
            {
                GenerateAndDrawTriangularMesh(bitmap, fixedArraySize, allNodes, apertureShape, D_value, obstructionRatio);
            }
        }

        private static (double centerI, double centerJ, double width, double height)? FindIlluminatedBounds(int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            // Use direct mathematical calculation instead of mask scanning
            // This provides exact bounds regardless of sampling
            
            // Calculate physical coordinates
            double aperturePixelSize = D_value / arraySize; // mm per pixel in aperture
            double centerI = arraySize / 2.0;
            double centerJ = arraySize / 2.0;
            
            switch (apertureShape)
            {
                case "Circular":
                case "Circular with obstr.":
                    {
                        // For circular apertures, the bounds are determined by the outer radius
                        double radius = D_value / 2.0; // radius in mm
                        double radiusInPixels = radius / aperturePixelSize;
                        
                        double width = radiusInPixels * 2;
                        double height = radiusInPixels * 2;
                        
                        return (centerI, centerJ, width, height);
                    }
                    
                case "Triangular":
                    {
                        // For triangular aperture, use the circumscribing circle bounds (same as circular)
                        double radius = D_value / 2.0; // radius in mm
                        double radiusInPixels = radius / aperturePixelSize;
                        
                        double width = radiusInPixels * 2;
                        double height = radiusInPixels * 2;
                        
                        return (centerI, centerJ, width, height);
                    }
                    
                case "Square":
                default:
                    {
                        // For square aperture, use the full array size
                        return (centerI, centerJ, arraySize, arraySize);
                    }
            }
        }

        private static List<(double, double, bool)> AddBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, 
            (double centerI, double centerJ, double width, double height) bounds, string apertureShape, double D_value, double obstructionRatio)
        {
            List<(double, double, bool)> boundaryNodes = new List<(double, double, bool)>();
            
            if (apertureShape == "Circular with obstr.")
            {
                // For annular apertures, we need both outer and inner boundary nodes
                AddAnnularBoundaryNodes(bitmap, arraySize, centerI, centerJ, apertureShape, D_value, obstructionRatio, boundaryNodes);
            }
            else
            {
                // For simple circular apertures, only outer boundary
                AddSimpleBoundaryNodes(bitmap, arraySize, centerI, centerJ, apertureShape, D_value, obstructionRatio, boundaryNodes);
            }
            
            return boundaryNodes;
        }
        
        private static void AddSimpleBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, 
            string apertureShape, double D_value, double obstructionRatio, List<(double, double, bool)> boundaryNodes)
        {
            // Create rays from center at 3-degree intervals for outer boundary
            int numRays = 120; // 360 / 3 = 120 rays
            for (int rayIndex = 0; rayIndex < numRays; rayIndex++)
            {
                double angleDegrees = rayIndex * 3.0;
                double angleRadians = angleDegrees * Math.PI / 180.0;
                
                var boundaryIntersection = FindRayBoundaryIntersection(centerI, centerJ, angleRadians, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (boundaryIntersection.HasValue)
                {
                    DrawDotAtPosition(bitmap, boundaryIntersection.Value.Item1, boundaryIntersection.Value.Item2, arraySize, Color.Red);
                    boundaryNodes.Add((boundaryIntersection.Value.Item1, boundaryIntersection.Value.Item2, true));
                }
            }
        }
        
        private static void AddAnnularBoundaryNodes(Bitmap bitmap, int arraySize, double centerI, double centerJ, 
            string apertureShape, double D_value, double obstructionRatio, List<(double, double, bool)> boundaryNodes)
        {
            // Calculate radii
            double outerRadius = D_value / 2.0; // mm
            double innerRadius = (D_value * obstructionRatio / 100.0) / 2.0; // mm
            double aperturePixelSize = D_value / arraySize;
            double outerRadiusPixels = outerRadius / aperturePixelSize;
            double innerRadiusPixels = innerRadius / aperturePixelSize;
            
            // Calculate appropriate number of rays for each boundary to maintain similar node spacing
            int outerRays = 120; // Same density as before for outer boundary
            int innerRays = Math.Max(12, (int)(120 * innerRadiusPixels / outerRadiusPixels)); // Scale by radius ratio, minimum 12 rays
            
            // Add outer boundary nodes (need special method for annular apertures)
            for (int rayIndex = 0; rayIndex < outerRays; rayIndex++)
            {
                double angleDegrees = rayIndex * (360.0 / outerRays);
                double angleRadians = angleDegrees * Math.PI / 180.0;
                
                var outerIntersection = FindOuterBoundaryIntersection(centerI, centerJ, angleRadians, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (outerIntersection.HasValue)
                {
                    DrawDotAtPosition(bitmap, outerIntersection.Value.Item1, outerIntersection.Value.Item2, arraySize, Color.Red);
                    boundaryNodes.Add((outerIntersection.Value.Item1, outerIntersection.Value.Item2, true));
                }
            }
            
            // Add inner boundary nodes (around the central obstruction)
            for (int rayIndex = 0; rayIndex < innerRays; rayIndex++)
            {
                double angleDegrees = rayIndex * (360.0 / innerRays);
                double angleRadians = angleDegrees * Math.PI / 180.0;
                
                var innerIntersection = FindInnerBoundaryIntersection(centerI, centerJ, angleRadians, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (innerIntersection.HasValue)
                {
                    DrawDotAtPosition(bitmap, innerIntersection.Value.Item1, innerIntersection.Value.Item2, arraySize, Color.Red);
                    boundaryNodes.Add((innerIntersection.Value.Item1, innerIntersection.Value.Item2, true));
                }
            }
        }

        private static (double, double)? FindRayBoundaryIntersection(double centerI, double centerJ, double angleRadians, int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            double dirI = Math.Cos(angleRadians);
            double dirJ = Math.Sin(angleRadians);
            
            // Start from center and march outward along the ray
            double step = 0.1;
            double distance = 0;
            
            // Ensure we start from an illuminated position
            double startDistance = 1.0; // Start a bit away from center to ensure illuminated
            
            // March outward to find the boundary
            bool lastWasIlluminated = true;
            double lastIlluminatedDistance = startDistance;
            
            // March outward with coarse steps first
            for (distance = startDistance; distance < arraySize; distance += step)
            {
                double currentI = centerI + distance * dirI;
                double currentJ = centerJ + distance * dirJ;
                
                bool isIlluminated = IsPositionIlluminatedDirect(currentI, currentJ, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (lastWasIlluminated && !isIlluminated)
                {
                    // Found the transition from illuminated to non-illuminated
                    // Use precise boundary refinement
                    var precisePosition = FindPreciseBoundaryWithIterativeRefinement(
                        centerI, centerJ, dirI, dirJ, lastIlluminatedDistance, distance, arraySize, apertureShape, D_value, obstructionRatio);
                    
                    if (precisePosition.HasValue)
                    {
                        return precisePosition;
                    }
                }
                
                if (isIlluminated)
                {
                    lastIlluminatedDistance = distance;
                }
                lastWasIlluminated = isIlluminated;
            }
            
            // If we reach here, we hit array boundary without finding aperture boundary
            // Calculate exact intersection with array boundary
            double maxDistance = Math.Min(
                Math.Abs((arraySize - 1 - centerI) / dirI),
                Math.Abs((arraySize - 1 - centerJ) / dirJ)
            );
            if (dirI < 0) maxDistance = Math.Min(maxDistance, Math.Abs(centerI / dirI));
            if (dirJ < 0) maxDistance = Math.Min(maxDistance, Math.Abs(centerJ / dirJ));
            
            double edgeDistance = Math.Max(0, maxDistance - 1); // Stay within bounds
            double finalI = centerI + edgeDistance * dirI;
            double finalJ = centerJ + edgeDistance * dirJ;
            
            return (finalI, finalJ);
        }

        private static (double, double)? FindInnerBoundaryIntersection(double centerI, double centerJ, double angleRadians, int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            double dirI = Math.Cos(angleRadians);
            double dirJ = Math.Sin(angleRadians);
            
            // For inner boundary, we start from center and move outward until we hit the inner edge of the illuminated area
            // This is the transition from non-illuminated (obstruction) to illuminated (annular region)
            
            double step = 0.1;
            double distance = 0;
            
            // Start from center (which is not illuminated in annular aperture) and march outward
            for (distance = step; distance < arraySize; distance += step)
            {
                double currentI = centerI + distance * dirI;
                double currentJ = centerJ + distance * dirJ;
                
                bool isIlluminated = IsPositionIlluminatedDirect(currentI, currentJ, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (isIlluminated)
                {
                    // Found the transition from non-illuminated (obstruction) to illuminated
                    // Use precise boundary refinement to find exact inner edge
                    var precisePosition = FindPreciseBoundaryWithIterativeRefinement(
                        centerI, centerJ, dirI, dirJ, distance - step, distance, arraySize, apertureShape, D_value, obstructionRatio);
                    
                    if (precisePosition.HasValue)
                    {
                        return precisePosition;
                    }
                    
                    // Fallback: use the current position
                    return (currentI, currentJ);
                }
            }
            
            // Should not reach here for valid annular apertures
            return null;
        }

        private static (double, double)? FindOuterBoundaryIntersection(double centerI, double centerJ, double angleRadians, int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            double dirI = Math.Cos(angleRadians);
            double dirJ = Math.Sin(angleRadians);
            
            // For annular apertures, we need to find the LAST transition from illuminated to non-illuminated
            // as we march outward (this will be the outer boundary, not the inner one)
            
            double step = 0.1;
            double lastIlluminatedDistance = 0;
            bool foundIlluminated = false;
            
            // March outward and track the last illuminated position
            for (double distance = step; distance < arraySize; distance += step)
            {
                double currentI = centerI + distance * dirI;
                double currentJ = centerJ + distance * dirJ;
                
                bool isIlluminated = IsPositionIlluminatedDirect(currentI, currentJ, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (isIlluminated)
                {
                    lastIlluminatedDistance = distance;
                    foundIlluminated = true;
                }
                else if (foundIlluminated)
                {
                    // Found transition from illuminated to non-illuminated (outer boundary)
                    var precisePosition = FindPreciseBoundaryWithIterativeRefinement(
                        centerI, centerJ, dirI, dirJ, lastIlluminatedDistance, distance, arraySize, apertureShape, D_value, obstructionRatio);
                    
                    if (precisePosition.HasValue)
                    {
                        return precisePosition;
                    }
                    
                    // Fallback: use the last known illuminated position
                    double fallbackI = centerI + lastIlluminatedDistance * dirI;
                    double fallbackJ = centerJ + lastIlluminatedDistance * dirJ;
                    return (fallbackI, fallbackJ);
                }
            }
            
            // If we reach here and found illuminated area, use the last illuminated position
            if (foundIlluminated)
            {
                double finalI = centerI + lastIlluminatedDistance * dirI;
                double finalJ = centerJ + lastIlluminatedDistance * dirJ;
                return (finalI, finalJ);
            }
            
            return null;
        }

        private static (double, double)? FindPreciseBoundaryWithIterativeRefinement(double centerI, double centerJ, double dirI, double dirJ, 
            double illuminatedDistance, double nonIlluminatedDistance, int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            double currentDistance = illuminatedDistance;
            double stepSize = (nonIlluminatedDistance - illuminatedDistance) / 2.0;
            bool movingForward = true;
            
            // Iterative refinement with decreasing step size
            for (int iteration = 0; iteration < 20 && stepSize > 0.001; iteration++)
            {
                if (movingForward)
                {
                    currentDistance += stepSize;
                }
                else
                {
                    currentDistance -= stepSize;
                }
                
                double testI = centerI + currentDistance * dirI;
                double testJ = centerJ + currentDistance * dirJ;
                
                bool isIlluminated = IsPositionIlluminatedDirect(testI, testJ, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (isIlluminated != movingForward)
                {
                    // We overshot, reverse direction and reduce step size
                    movingForward = !movingForward;
                    stepSize *= 0.5;
                }
            }
            
            // Ensure final position is just inside the illuminated area
            currentDistance -= 0.1; // Step back slightly to ensure we're inside
            
            double finalI = centerI + currentDistance * dirI;
            double finalJ = centerJ + currentDistance * dirJ;
            
            return (finalI, finalJ);
        }

        private static bool IsPositionIlluminatedDirect(double i, double j, int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            // Check array bounds
            if (i < 0 || i >= arraySize || j < 0 || j >= arraySize)
                return false;
            
            // Calculate physical coordinates - same logic as ApertureMaskGenerator
            double aperturePixelSize = D_value / arraySize; // mm per pixel in aperture
            double centerX = (arraySize - 1) / 2.0; // Same as ApertureMaskGenerator
            double centerY = (arraySize - 1) / 2.0;
            double x = (i - centerX) * aperturePixelSize;
            double y = (j - centerY) * aperturePixelSize;
            double distanceFromCenter = Math.Sqrt(x * x + y * y);
            double radius = D_value / 2.0;
            
            switch (apertureShape)
            {
                case "Circular":
                    return distanceFromCenter <= radius;
                    
                case "Circular with obstr.":
                    double innerRadius = (D_value * obstructionRatio / 100.0) / 2.0;
                    return distanceFromCenter <= radius && distanceFromCenter >= innerRadius;
                    
                case "Triangular":
                    return IsPointInTriangleGrid(x, y, D_value);
                    
                case "Square":
                default:
                    return true; // Square aperture covers the entire array
            }
        }

        private static bool IsPointInTriangleGrid(double x, double y, double D_value)
        {
            // Physical coordinates are already provided, just use them directly
            double radius = D_value / 2.0; // radius of circumscribing circle in mm
            
            // Calculate triangle vertices in physical coordinates (mm) - exactly same as mask generation
            // Left vertex (pointing leftward in physical space)
            double vertex1X_mm = -radius; // negative X is left
            double vertex1Y_mm = 0;
            
            // Bottom-left vertex 
            double vertex2X_mm = radius / 2.0;
            double vertex2Y_mm = -radius * Math.Sqrt(3) / 2.0;
            
            // Bottom-right vertex
            double vertex3X_mm = radius / 2.0;
            double vertex3Y_mm = radius * Math.Sqrt(3) / 2.0;
            
            // Use barycentric coordinates to determine if point is inside triangle
            return IsPointInTriangleBary(x, y, vertex1X_mm, vertex1Y_mm, vertex2X_mm, vertex2Y_mm, vertex3X_mm, vertex3Y_mm);
        }
        
        private static bool IsPointInTriangleBary(double px, double py, double ax, double ay, double bx, double by, double cx, double cy)
        {
            // Use barycentric coordinates to determine if point P is inside triangle ABC
            double denominator = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
            
            if (Math.Abs(denominator) < 1e-10) return false; // Degenerate triangle
            
            double alpha = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / denominator;
            double beta = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / denominator;
            double gamma = 1 - alpha - beta;
            
            // Point is inside triangle if all barycentric coordinates are non-negative
            return alpha >= 0 && beta >= 0 && gamma >= 0;
        }

        private static void DrawDotAtPosition(Bitmap bitmap, double i, double j, int arraySize, Color color)
        {
            // Convert to bitmap coordinates - same as MainForm.cs display
            int pixelX = (int)Math.Round(j); // j maps to x (column) - same as MainForm
            int pixelY = (int)Math.Round(i); // i maps to y (row) - same as MainForm
            
            // Draw a 3x3 dot centered at the position
            for (int di = -1; di <= 1; di++)
            {
                for (int dj = -1; dj <= 1; dj++)
                {
                    int drawX = pixelX + dj;
                    int drawY = pixelY + di;
                    
                    if (drawX >= 0 && drawX < bitmap.Width && drawY >= 0 && drawY < bitmap.Height)
                    {
                        bitmap.SetPixel(drawX, drawY, color);
                    }
                }
            }
        }

        private static void GenerateAndDrawTriangularMesh(Bitmap bitmap, int arraySize, List<(double i, double j, bool isBoundary)> allNodes, string apertureShape, double D_value, double obstructionRatio)
        {
            if (allNodes.Count < 3)
            {
                return;
            }
            
            var triangles = GenerateDelaunayTriangulation(allNodes);
            
            // Filter out triangles that span across non-illuminated areas (for annular apertures)
            var filteredTriangles = FilterValidTriangles(triangles, allNodes, arraySize, apertureShape, D_value, obstructionRatio);
            
            DrawTriangularMesh(bitmap, filteredTriangles, allNodes, arraySize);
        }
        
        private static List<(int node1, int node2, int node3)> FilterValidTriangles(List<(int node1, int node2, int node3)> triangles, 
            List<(double i, double j, bool isBoundary)> nodes, int arraySize, string apertureShape, double D_value, double obstructionRatio)
        {
            var validTriangles = new List<(int node1, int node2, int node3)>();
            
            foreach (var triangle in triangles)
            {
                var node1 = nodes[triangle.node1];
                var node2 = nodes[triangle.node2];
                var node3 = nodes[triangle.node3];
                
                // Calculate triangle centroid
                double centroidI = (node1.i + node2.i + node3.i) / 3.0;
                double centroidJ = (node1.j + node2.j + node3.j) / 3.0;
                
                // Check if triangle centroid is in illuminated area
                bool centroidIlluminated = IsPositionIlluminatedDirect(centroidI, centroidJ, arraySize, apertureShape, D_value, obstructionRatio);
                
                if (centroidIlluminated)
                {
                    validTriangles.Add(triangle);
                }
            }
            
            return validTriangles;
        }

        private static List<(int node1, int node2, int node3)> GenerateDelaunayTriangulation(List<(double i, double j, bool isBoundary)> nodes)
        {
            var triangles = new List<(int v1, int v2, int v3)>();
            
            if (nodes.Count < 3) return new List<(int, int, int)>();
            
            // Create a list of nodes for the algorithm
            var allNodes = nodes.Select(n => (n.i, n.j)).ToList();
            
            // Create a super-triangle that encompasses all points
            double minI = allNodes.Min(p => p.i) - 10;
            double maxI = allNodes.Max(p => p.i) + 10;
            double minJ = allNodes.Min(p => p.j) - 10;
            double maxJ = allNodes.Max(p => p.j) + 10;
            
            double dx = maxI - minI;
            double dy = maxJ - minJ;
            
            // Create super-triangle vertices (much larger than the point set)
            var superTriangle = new List<(double i, double j)>
            {
                (minI - dx, minJ - dy),
                (maxI + dx, minJ - dy),
                (minI + dx/2, maxJ + dy)
            };
            
            // Add super-triangle to nodes
            allNodes.AddRange(superTriangle);
            
            // Initialize with super-triangle
            triangles.Add((nodes.Count, nodes.Count + 1, nodes.Count + 2));
            
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
                
                // Find the boundary of the polygonal hole using the working algorithm from MainForm
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
                
                // Remove bad triangles (sort descending to preserve indices)
                badTriangles.Sort((a, b) => b.CompareTo(a));
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
            
            return finalTriangles;
        }

        private static bool IsPointInCircumcircle((double i, double j) point, (double i, double j) a, (double i, double j) b, (double i, double j) c)
        {
            double ax = a.i, ay = a.j;
            double bx = b.i, by = b.j;
            double cx = c.i, cy = c.j;
            double px = point.i, py = point.j;
            
            double ax_ = ax - px;
            double ay_ = ay - py;
            double bx_ = bx - px;
            double by_ = by - py;
            double cx_ = cx - px;
            double cy_ = cy - py;
            
            double det = (ax_ * ax_ + ay_ * ay_) * (bx_ * cy_ - by_ * cx_) +
                        (bx_ * bx_ + by_ * by_) * (cx_ * ay_ - cy_ * ax_) +
                        (cx_ * cx_ + cy_ * cy_) * (ax_ * by_ - ay_ * bx_);
            
            return det > 0;
        }

        private static void DrawTriangularMesh(Bitmap bitmap, List<(int node1, int node2, int node3)> triangles, List<(double i, double j, bool isBoundary)> nodes, int arraySize)
        {
            Color meshColor = Color.Green;
            
            foreach (var triangle in triangles)
            {
                var node1 = nodes[triangle.node1];
                var node2 = nodes[triangle.node2];
                var node3 = nodes[triangle.node3];
                
                // Draw triangle edges
                DrawLine(bitmap, (node1.i, node1.j), (node2.i, node2.j), meshColor);
                DrawLine(bitmap, (node2.i, node2.j), (node3.i, node3.j), meshColor);
                DrawLine(bitmap, (node3.i, node3.j), (node1.i, node1.j), meshColor);
            }
        }

        private static void DrawLine(Bitmap bitmap, (double i, double j) start, (double i, double j) end, Color color)
        {
            // Use the same coordinate system as MainForm.cs display
            int x1 = (int)Math.Round(start.j); // j maps to x (column) - same as MainForm
            int y1 = (int)Math.Round(start.i); // i maps to y (row) - same as MainForm
            int x2 = (int)Math.Round(end.j);
            int y2 = (int)Math.Round(end.i);
            
            // Bresenham's line algorithm
            int dx = Math.Abs(x2 - x1);
            int dy = Math.Abs(y2 - y1);
            int sx = x1 < x2 ? 1 : -1;
            int sy = y1 < y2 ? 1 : -1;
            int err = dx - dy;
            
            int x = x1, y = y1;
            
            while (true)
            {
                // Draw a single pixel line
                if (x >= 0 && x < bitmap.Width && y >= 0 && y < bitmap.Height)
                {
                    bitmap.SetPixel(x, y, color);
                }
                
                if (x == x2 && y == y2) break;
                
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
