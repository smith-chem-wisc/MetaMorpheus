using System;
using System.Linq;

namespace GuiFunctions;

/// <summary>
/// The Hungarian algorithm, also known as the Kuhn-Munkres algorithm or Munkres assignment algorithm, is a powerful combinatorial optimization algorithm used to solve the assignment problem. The goal is to find the optimal assignment between two sets of equal size(e.g., psms and precursor envelopes) to minimize the total cost or maximize the total value, given a cost matrix representing the cost of assigning each psms to each envelopes where cost is mass error.
/// </summary>
public static class HungarianAlgorithm
{
    /// <summary>
    /// Solves the assignment problem (also known as the minimum cost bipartite matching problem) using the Hungarian algorithm.
    /// Given a cost matrix where each entry represents the cost of assigning a row (e.g., a PSM) to a column (e.g., an isotopic envelope),
    /// this method finds the optimal one-to-one assignment that minimizes the total cost.
    /// If the input matrix is not square, it is padded with zeros to make it square.
    /// Returns an array where each element at index i indicates the column assigned to row i, or -1 if no assignment is possible.
    /// </summary>
    /// <param name="costMatrix">
    /// A 2D array of assignment costs. Rows typically represent items to assign (e.g., PSMs), columns represent targets (e.g., envelopes).
    /// Use a large value (e.g., double.MaxValue) to indicate forbidden assignments.
    /// </param>
    /// <returns>
    /// An array of assignments: result[i] = j means row i is assigned to column j; -1 means row i is unassigned.
    /// </returns>
    public static int[] FindAssignments(double[,] costMatrix)
    {
        int rows = costMatrix.GetLength(0);
        int cols = costMatrix.GetLength(1);
        int dim = Math.Max(rows, cols);

        // Pad to square matrix if needed
        double[,] cost = new double[dim, dim];
        for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            cost[i, j] = (i < rows && j < cols) ? costMatrix[i, j] : 0;

        // Step 1: Subtract row minima
        for (int i = 0; i < dim; i++)
        {
            double min = double.MaxValue;
            for (int j = 0; j < dim; j++)
                if (cost[i, j] < min) min = cost[i, j];
            for (int j = 0; j < dim; j++)
                cost[i, j] -= min;
        }

        // Step 2: Subtract column minima
        for (int j = 0; j < dim; j++)
        {
            double min = double.MaxValue;
            for (int i = 0; i < dim; i++)
                if (cost[i, j] < min) min = cost[i, j];
            for (int i = 0; i < dim; i++)
                cost[i, j] -= min;
        }

        // Marked zeros
        int[] rowCover = new int[dim];
        int[] colCover = new int[dim];
        int[,] mask = new int[dim, dim];

        // Step 3: Cover zeros with minimum number of lines
        int step = 1;
        int[] pathRow = new int[dim * 2];
        int[] pathCol = new int[dim * 2];
        int pathCount = 0;

        while (true)
        {
            switch (step)
            {
                case 1:
                    // Star zeros
                    for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        if (cost[i, j] == 0 && rowCover[i] == 0 && colCover[j] == 0)
                        {
                            mask[i, j] = 1;
                            rowCover[i] = 1;
                            colCover[j] = 1;
                        }
                    Array.Clear(rowCover, 0, dim);
                    Array.Clear(colCover, 0, dim);
                    step = 2;
                    break;

                case 2:
                    // Cover columns with starred zeros
                    for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        if (mask[i, j] == 1)
                            colCover[j] = 1;
                    if (colCover.Sum() >= dim)
                        step = 7;
                    else
                        step = 3;
                    break;

                case 3:
                    // Prime uncovered zeros
                    bool done = false;
                    int zRow = -1, zCol = -1;
                    while (!done)
                    {
                        (zRow, zCol) = FindZero(cost, rowCover, colCover, dim);
                        if (zRow == -1)
                        {
                            step = 5;
                            done = true;
                        }
                        else
                        {
                            mask[zRow, zCol] = 2;
                            // Find the column in this row where mask[zRow, col] == 1
                            int starCol = -1;
                            for (int j = 0; j < dim; j++)
                            {
                                if (mask[zRow, j] == 1)
                                {
                                    starCol = j;
                                    break;
                                }
                            }
                            if (starCol != -1)
                            {
                                rowCover[zRow] = 1;
                                colCover[starCol] = 0;
                            }
                            else
                            {
                                step = 4;
                                pathRow[0] = zRow;
                                pathCol[0] = zCol;
                                pathCount = 1;
                                done = true;
                            }
                        }
                    }
                    break;

                case 4:
                    // Augmenting path
                    bool found;
                    do
                    {
                        int row = FindStarInCol(mask, pathCol[pathCount - 1], dim);
                        if (row != -1)
                        {
                            pathRow[pathCount] = row;
                            pathCol[pathCount] = pathCol[pathCount - 1];
                            pathCount++;
                        }
                        else
                        {
                            found = false;
                            break;
                        }
                        int col = FindPrimeInRow(mask, pathRow[pathCount - 1], dim);
                        pathRow[pathCount] = pathRow[pathCount - 1];
                        pathCol[pathCount] = col;
                        pathCount++;
                        found = true;
                    } while (found);

                    // Convert path
                    for (int i = 0; i < pathCount; i++)
                        mask[pathRow[i], pathCol[i]] = mask[pathRow[i], pathCol[i]] == 1 ? 0 : 1;

                    // Clear covers and primes
                    Array.Clear(rowCover, 0, dim);
                    Array.Clear(colCover, 0, dim);
                    for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        if (mask[i, j] == 2)
                            mask[i, j] = 0;
                    step = 2;
                    break;

                case 5:
                    // Add minimum uncovered value to every element of covered rows, subtract from uncovered columns
                    double minVal = double.MaxValue;
                    for (int i = 0; i < dim; i++)
                        if (rowCover[i] == 0)
                            for (int j = 0; j < dim; j++)
                                if (colCover[j] == 0 && cost[i, j] < minVal)
                                    minVal = cost[i, j];
                    for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                    {
                        if (rowCover[i] == 1)
                            cost[i, j] += minVal;
                        if (colCover[j] == 0)
                            cost[i, j] -= minVal;
                    }
                    step = 3;
                    break;

                case 7:
                    // Done
                    int[] result = new int[rows];
                    for (int i = 0; i < rows; i++)
                    {
                        result[i] = -1;
                        for (int j = 0; j < cols; j++)
                            if (mask[i, j] == 1)
                                result[i] = j;
                    }
                    return result;
            }
        }
    }

    private static (int, int) FindZero(double[,] cost, int[] rowCover, int[] colCover, int dim)
    {
        for (int i = 0; i < dim; i++)
            if (rowCover[i] == 0)
                for (int j = 0; j < dim; j++)
                    if (cost[i, j] == 0 && colCover[j] == 0)
                        return (i, j);
        return (-1, -1);
    }

    private static int FindStarInCol(int[,] mask, int col, int dim)
    {
        for (int i = 0; i < dim; i++)
            if (mask[i, col] == 1)
                return i;
        return -1;
    }

    private static int FindPrimeInRow(int[,] mask, int row, int dim)
    {
        for (int j = 0; j < dim; j++)
            if (mask[row, j] == 2)
                return j;
        return -1;
    }
}