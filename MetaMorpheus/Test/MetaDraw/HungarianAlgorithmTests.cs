using GuiFunctions;
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
}