using FEM.Points;

namespace FEM.Problems
{
    public abstract class Problem
    {
        public abstract double Gamma(int subarea);
        public abstract double GetPointConcentratedSource(Point2D point, int subaria);
        public abstract double Lambda(int subarea);
        public abstract double RightPart(Point2D point, int subarea);
        public abstract double Theta(Point2D point, int subarea);
        public abstract double Beta(Point2D point, int subaria);
        public abstract double UBeta(Point2D point, int subaria);
        public virtual double Bc1(Point2D point, int subaria) { return 0; }

    }
}
