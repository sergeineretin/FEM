using FEM.Points;

namespace FEM.Problems
{
    public sealed class TestBoundaryValueProblem : Problem
    {
        public override double RightPart(Point2D point, int subarea)
        {
            //return Math.Log(1 / Math.Sqrt(point.X * point.X + point.Y * point.Y + 100), 3);
            return -6 * (point.Y + point.X) + 2 * (point.X * point.X * point.X + point.Y * point.Y * point.Y);
            //return -12*(point.X * point.X + point.Y * point.Y) + (point.X * point.X * point.X * point.X + point.Y * point.Y * point.Y * point.Y);
        }
        public override double Theta(Point2D point, int subarea)
        {
            return 0;
        }
        public override double Beta(Point2D point, int subaria)
        {
            return 0;
        }
        public override double UBeta(Point2D point, int subaria)
        {
            return 0;
        }
        public override double Lambda(int subarea)
        {
            return 1;
        }
        public override double Gamma(int subarea)
        {
            return 2;
        }
        public override double Bc1(Point2D point, int subaria)
        {
            //return 1;
            return point.X * point.X * point.X + point.Y * point.Y * point.Y;
            //return point.X * point.X * point.X * point.X + point.Y * point.Y * point.Y * point.Y;
        }

        public override double GetPointConcentratedSource(Point2D point, int subaria)
        {
            throw new NotImplementedException();
        }
    }
}
