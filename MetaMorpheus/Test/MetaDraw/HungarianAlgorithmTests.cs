using GuiFunctions.MetaDraw.Chimeras;
using NUnit.Framework;

namespace Test.MetaDraw;

[TestFixture]
public class HungarianAlgorithmTests
{
    [Test]
    public void SimpleSquareMatrix_AssignsDiagonal()
    {
        // Arrange
        double[,] cost = new double[,]
        {
            { 1, 2, 3 },
            { 2, 1, 3 },
            { 3, 2, 1 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        Assert.That(assignment, Is.EqualTo(new[] { 0, 1, 2 }));
    }

    [Test]
    public void RectangularMatrix_AssignsOptimal()
    {
        // Arrange: 2x3 matrix, optimal is row 0 to col 1, row 1 to col 2
        double[,] cost = new double[,]
        {
            { 10, 1, 10 },
            { 10, 10, 1 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        Assert.That(assignment[0], Is.EqualTo(1));
        Assert.That(assignment[1], Is.EqualTo(2));
    }

    [Test]
    public void ChargeMismatchPenalty_AvoidsHighCost()
    {
        // Arrange: simulate charge mismatch with high penalty
        double[,] cost = new double[,]
        {
            { 0, 1000 },
            { 1000, 0 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        Assert.That(assignment[0], Is.EqualTo(0));
        Assert.That(assignment[1], Is.EqualTo(1));
    }

    [Test]
    public void Case3_PrimeUncoveredZeros_LeadsToAugmentingPath()
    {
        // Arrange: This matrix will require priming and an augmenting path
        double[,] cost = new double[,]
        {
            { 4, 1, 3 },
            { 2, 0, 5 },
            { 3, 2, 2 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        // The optimal assignment is [1,0,2] (row 0->1, row 1->0, row 2->2)
        Assert.That(assignment, Is.EqualTo(new[] { 1, 0, 2 }));
    }

    [Test]
    public void Case4_AugmentingPath_TogglesAssignments()
    {
        // Arrange: This matrix will require an augmenting path to toggle assignments
        double[,] cost = new double[,]
        {
            { 1, 2, 3 },
            { 2, 4, 6 },
            { 3, 6, 9 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        // The optimal assignment is [2,1,0]
        Assert.That(assignment, Is.EqualTo(new[] { 2, 1, 0 }));
    }

    [Test]
    public void Case5_AddMinimumUncoveredValue_CreatesNewZero()
    {
        // Arrange: This matrix will require adding the minimum uncovered value to create a new zero
        double[,] cost = new double[,]
        {
            { 10, 19, 8 },
            { 10, 18, 7 },
            { 13, 16, 9 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        // The optimal assignment is [0,2,1]
        Assert.That(assignment, Is.EqualTo(new[] { 0, 2, 1 }));
    }

    [Test]
    public void ForbiddenAssignments_AreNotChosen()
    {
        // Arrange: Use double.MaxValue to forbid certain assignments
        double M = double.MaxValue;
        double[,] cost = new double[,]
        {
            { 1, M, 3 },
            { 2, 1, M },
            { M, 2, 1 }
        };

        // Act
        int[] assignment = HungarianAlgorithm.FindAssignments(cost);

        // Assert
        // The optimal assignment is [0,1,2]
        Assert.That(assignment, Is.EqualTo(new[] { 0, 1, 2 }));
    }
}
