namespace FEM.Points
{
    public sealed class Point2D : Point, IComparable<Point2D>
    {
        public Point2D(double x, double y)
        {
            X = x;
            Y = y;
        }
        public double X { get; set; }
        public double Y { get; set; }

        public int CompareTo(Point2D other)
        {
            if (other == null)
            {
                // Define your behavior when 'other' is null, depending on your requirements.
                // For example, you might consider throwing an exception or returning a default value.
                throw new ArgumentNullException(nameof(other));
            }
            // Сначала сравниваем по Y
            int xComparison = Y.CompareTo(other.Y);

            // Если X равны, сравниваем по X
            return xComparison == 0 ? X.CompareTo(other.X) : xComparison;
        }
        public override bool Equals(object? obj)
        {
            return obj is Point2D d &&
                   X == d.X &&
                   Y == d.Y;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(X, Y);
        }

        public override string ToString()
        {
            return $"{X} {Y}";
        }

    }
}
