using FEM.Points;

namespace FEM
{
    struct Subinterval
    {
        public int NumberOfSplits { get; private set; }
        public double Coeff { get; private set; }
        public Subinterval(int numberOfSplits, double coeff)
        {
            NumberOfSplits = numberOfSplits;
            Coeff = coeff;
        }
        public override string ToString()
        {
            return $"{NumberOfSplits} {Coeff}";
        }
    }
    internal static class SegmentSeparator
    {
        private static List<double> GetSublist(double first, double last, Subinterval subinterval)
        {
            List<double> sublist = new List<double>();
            double interval = 0;

            if (subinterval.Coeff == 1)
                interval = (last - first) / subinterval.NumberOfSplits;
            else
                interval = (last - first) * (subinterval.Coeff - 1)
                    / (Math.Pow(subinterval.Coeff, subinterval.NumberOfSplits) - 1);

            sublist.Add(first + interval);
            interval *= subinterval.Coeff;
            for (int i = 1; i < subinterval.NumberOfSplits - 1; i++)
            {
                sublist.Add(sublist.Last() + interval);
                interval *= subinterval.Coeff;
            }
            return sublist;
        }
        public static List<Point2D> GetPointsOfLine(Point2D firstPoint, Point2D secondPoint, Subinterval intervalSplit)
        {
            List<Point2D> points = new List<Point2D>();
            points.Add(firstPoint);

            double xProjection = secondPoint.X - firstPoint.X;
            double yProjection = secondPoint.Y - firstPoint.Y;

            double length = Math.Sqrt(xProjection * xProjection + yProjection * yProjection);
            List<double> midwayIntervals = GetSublist(0, length, intervalSplit);

            foreach (double interval in midwayIntervals)
            {
                points.Add(new(firstPoint.X + interval * xProjection / length, firstPoint.Y + interval * yProjection / length));
            }
            points.Add(secondPoint);

            return points;
        }
    }
}
