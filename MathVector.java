/**
 * A MathVector class that allows basic vector operations for both 2D and 3D
 * vectors.
 * 
 * @author franklin
 */
public class MathVector implements Comparable<MathVector>
{
    /**
     * Stores the vector data for this MathVector.
     */
    private double[] vector;

    /**
     * Index position for the x-coordinate.
     */
    private static final int X_SPOT = 0;

    /**
     * Index position for the y-coordinate.
     */
    private static final int Y_SPOT = 1;

    /**
     * Index position for the z-coordinate.
     */
    private static final int Z_SPOT = 2;

    /**
     * Number of elements required for 3D vectors.
     */
    private static final int NUM_FOR_3D_VECTORS = 3;

    /**
     * Constructor for 2D vectors.
     * 
     * @param x the x-coordinate
     * @param y the y-coordinate
     */
    public MathVector(double x, double y)
    {
        vector = new double[2];
        vector[X_SPOT] = x;
        vector[Y_SPOT] = y;
    }

    /**
     * Constructor for 3D vectors.
     * 
     * @param x the x-coordinate
     * @param y the y-coordinate
     * @param z the z-coordinate
     */
    public MathVector(double x, double y, double z)
    {
        vector = new double[NUM_FOR_3D_VECTORS];
        vector[X_SPOT] = x;
        vector[Y_SPOT] = y;
        vector[Z_SPOT] = z;
    }

    /**
     * Private constructor that initializes a vector of a given size.
     * 
     * @param vectorSize the number of dimensions (2 or 3)
     */
    private MathVector(int vectorSize)
    {
        vector = new double[vectorSize];
    }

    /**
     * Returns the number of dimensions of the vector.
     * 
     * @return the vector size
     */
    public int getVectorSize()
    {
        return vector.length;
    }

    /**
     * Returns the x-coordinate value.
     * 
     * @return the x-value
     */
    public double getXValue()
    {
        return vector[X_SPOT];
    }

    /**
     * Returns the y-coordinate value.
     * 
     * @return the y-value
     */
    public double getYValue()
    {
        return vector[Y_SPOT];
    }

    /**
     * Returns the z-coordinate value for 3D vectors.
     * 
     * @return the z-value
     * @throws IllegalArgumentException if called on a 2D vector
     */
    public double getZValue()
    {
        if (getVectorSize() != NUM_FOR_3D_VECTORS)
        {
            throw new IllegalArgumentException("Attempted to access Z value of a 2D vector.");
        }
        return vector[Z_SPOT];
    }

    /**
     * Checks if the vector is a 3D vector.
     * 
     * @return true if the vector is 3D, false otherwise
     */
    public boolean is3DVector()
    {
        return vector.length == NUM_FOR_3D_VECTORS;
    }

    /**
     * Checks if the vector is a 2D vector.
     * 
     * @return true if the vector is 2D, false otherwise
     */
    public boolean is2DVector()
    {
        return vector.length == NUM_FOR_3D_VECTORS - 1;
    }

    /**
     * Calculates and returns the magnitude of the vector.
     * 
     * @return the magnitude of the vector
     */
    public double getMagnitude()
    {
        double xSquare = getXValue() * getXValue();
        double ySquare = getYValue() * getYValue();
        if (is3DVector())
        {
            double zSquare = getZValue() * getZValue();
            return Math.sqrt(xSquare + ySquare + zSquare);
        }
        return Math.sqrt(xSquare + ySquare);
    }

    /**
     * Adds this vector to another vector and returns the result.
     * 
     * @param otherVector the vector to add
     * @return the resulting vector after addition
     * @throws IllegalArgumentException if the vectors are not the same size
     */
    public MathVector add(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException("Cannot add vectors of different dimensions.");
        }
        MathVector next = new MathVector(otherVector.getVectorSize());
        next.vector[X_SPOT] = getXValue() + otherVector.getXValue();
        next.vector[Y_SPOT] = getYValue() + otherVector.getYValue();
        if (is3DVector())
        {
            next.vector[Z_SPOT] = getZValue() + otherVector.getZValue();
        }
        return next;
    }

    /**
     * Subtracts another vector from this vector and returns the result.
     * 
     * @param otherVector the vector to subtract
     * @return the resulting vector after subtraction
     * @throws IllegalArgumentException if the vectors are not the same size
     */
    public MathVector subtract(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException("Cannot subtract vectors of different dimensions.");
        }
        MathVector next = new MathVector(otherVector.getVectorSize());
        next.vector[X_SPOT] = getXValue() - otherVector.getXValue();
        next.vector[Y_SPOT] = getYValue() - otherVector.getYValue();
        if (is3DVector())
        {
            next.vector[Z_SPOT] = getZValue() - otherVector.getZValue();
        }
        return next;
    }

    /**
     * Multiplies this vector by a scalar value and returns the result.
     * 
     * @param value the scalar to multiply by
     * @return the resulting vector after scalar multiplication
     */
    public MathVector scalarMultiplication(double value)
    {
        double nextX = getXValue() * value;
        double nextY = getYValue() * value;
        if (is3DVector())
        {
            double nextZ = getZValue() * value;
            return new MathVector(nextX, nextY, nextZ);
        }
        return new MathVector(nextX, nextY);
    }

    /**
     * Negates this vector (reverses its direction).
     * 
     * @return the negated vector
     */
    public MathVector negate()
    {
        return scalarMultiplication(-1);
    }

    /**
     * Returns the unit vector (a vector with magnitude 1) in the same direction as
     * this vector.
     * 
     * @return the unit vector
     */
    public MathVector getUnitVector()
    {
        if (is3DVector())
        {
            return new MathVector((getXValue() / getMagnitude()), (getYValue() / getMagnitude()),
                    (getZValue() / getMagnitude()));
        }
        return new MathVector((getXValue() / getMagnitude()), (getYValue() / getMagnitude()));
    }

    /**
     * Computes the dot product of this vector and another vector.
     * 
     * @param otherVector the other vector
     * @return the dot product
     * @throws IllegalArgumentException if the vectors are not the same size
     */
    public double dotProduct(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Cannot compute the dot product of vectors with different dimensions.");
        }
        double nextX = getXValue() * otherVector.getXValue();
        double nextY = getYValue() * otherVector.getYValue();
        if (is3DVector())
        {
            double nextZ = getZValue() * otherVector.getZValue();
            return nextX + nextY + nextZ;
        }
        return nextX + nextY;
    }

    /**
     * Computes the cross product of this vector and another 3D vector.
     * 
     * @param otherVector the other 3D vector
     * @return the resulting vector from the cross product
     * @throws IllegalArgumentException if either vector is not 3D
     */
    public MathVector crossProduct3D(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Both vectors must be 3D to compute the cross product.");
        }
        double nextX = (getYValue() * otherVector.getZValue())
                - (getZValue() * otherVector.getYValue());
        double nextY = (getZValue() * otherVector.getXValue())
                - (getXValue() * otherVector.getZValue());
        double nextZ = (getXValue() * otherVector.getYValue())
                - (getYValue() * otherVector.getXValue());
        return new MathVector(nextX, nextY, nextZ);
    }

    /**
     * Computes the cross product of this vector and another 2D vector.
     * 
     * @param otherVector the other 2D vector
     * @return the scalar value resulting from the cross product
     * @throws IllegalArgumentException if either vector is not 2D
     */
    public double crossProduct2D(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Both vectors must be 2D to compute the cross product.");
        }
        double value1 = getXValue() * otherVector.getYValue();
        double value2 = getYValue() * otherVector.getXValue();
        return value1 - value2;
    }

    /**
     * Computes the scalar triple product of this vector with two other vectors. The
     * scalar triple product is calculated as the dot product of this vector and the
     * cross product of vectors a and b.
     * 
     * @param a the first vector for the cross product.
     * @param b the second vector for the cross product.
     * @return the scalar triple product.
     * @throws IllegalArgumentException if the vectors are not of the same size.
     */
    public double scalarTripleProduct(MathVector a, MathVector b)
    {
        if (getVectorSize() != a.getVectorSize() || getVectorSize() != b.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Vectors must have the same dimensions for scalar triple product.");
        }
        return dotProduct(a.crossProduct3D(b));
    }

    /**
     * Computes the scalar projection of vector b onto this vector. The scalar
     * projection is the magnitude of the projection of one vector onto another.
     * 
     * @param b the vector to project.
     * @return the scalar projection of b onto this vector.
     * @throws IllegalArgumentException if the vectors are not of the same size.
     */
    public double scalarProjectionOf(MathVector b)
    {
        if (getVectorSize() != b.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Vectors must have the same dimensions for scalar projection.");
        }
        return dotProduct(b) / getMagnitude();
    }

    /**
     * Computes the vector projection of vector b onto this vector. The vector
     * projection is the projection of one vector onto another.
     * 
     * @param b the vector to project.
     * @return the vector projection of b onto this vector.
     * @throws IllegalArgumentException if the vectors are not of the same size.
     */
    public MathVector vectorProjectionOf(MathVector b)
    {
        if (getVectorSize() != b.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Vectors must have the same dimensions for vector projection.");
        }

        double xPost = scalarProjectionOf(b) * getXValue() / getMagnitude();
        double yPost = scalarProjectionOf(b) * getYValue() / getMagnitude();
        if (is3DVector())
        {
            double zPost = scalarProjectionOf(b) * getZValue() / getMagnitude();
            return new MathVector(xPost, yPost, zPost);
        }
        return new MathVector(xPost, yPost);
    }

    /**
     * Checks if this vector is parallel to another vector. Two vectors are parallel
     * if their dot product equals the product of their magnitudes.
     * 
     * @param otherVector the vector to compare against.
     * @return true if the vectors are parallel, false otherwise.
     * @throws IllegalArgumentException if the vectors are not of the same size.
     */
    public boolean isParallel(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Vectors must have the same dimensions to check parallelism.");
        }
        return Math.abs(dotProduct(otherVector)) == (getMagnitude() * otherVector.getMagnitude());
    }

    /**
     * Checks if this vector is perpendicular to another vector. Two vectors are
     * perpendicular if their dot product is zero.
     * 
     * @param otherVector the vector to compare against.
     * @return true if the vectors are perpendicular, false otherwise.
     * @throws IllegalArgumentException if the vectors are not of the same size.
     */
    public boolean isPerpendicular(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Vectors must have the same dimensions to check perpendicularity.");
        }
        return dotProduct(otherVector) == 0;
    }

    /**
     * Calculates the Euclidean distance between this vector and another vector.
     * 
     * @param otherVector the other vector
     * @return the distance between the two vectors
     * @throws IllegalArgumentException if the vectors are not the same size
     */
    public double distance(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Cannot compute the distance between vectors of different dimensions.");
        }
        double xDiff = getXValue() - otherVector.getXValue();
        double yDiff = getYValue() - otherVector.getYValue();
        if (is3DVector())
        {
            double zDiff = getZValue() - otherVector.getZValue();
            return Math.sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);
        }
        return Math.sqrt(xDiff * xDiff + yDiff * yDiff);
    }

    /**
     * Calculates the angle in radians between this vector and another vector.
     * 
     * @param otherVector the other vector
     * @return the angle in radians between the two vectors
     * @throws IllegalArgumentException if the vectors are not the same size
     */
    public double angleBetween(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            throw new IllegalArgumentException(
                    "Cannot compute the angle between vectors of different dimensions.");
        }
        double dotProd = dotProduct(otherVector);
        double magnitudes = getMagnitude() * otherVector.getMagnitude();
        return Math.acos(dotProd / magnitudes);
    }

    /**
     * Checks if this vector is equal to another vector.
     * 
     * @param otherVector the other vector
     * @return true if the vectors are equal, false otherwise
     */
    public boolean equals(MathVector otherVector)
    {
        if (getVectorSize() != otherVector.getVectorSize())
        {
            return false;
        }
        for (int i = 0; i < getVectorSize(); i++)
        {
            if (vector[i] != otherVector.vector[i])
            {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns a string representation of the vector.
     * 
     * @return a string in the form "<x, y, z>" for 3D vectors or "<x, y>" for 2D
     *         vectors
     */
    public String toString()
    {
        String s = "<" + getXValue() + ", " + getYValue();
        if (is3DVector())
        {
            s += ", " + getZValue();
        }
        s += ">";
        return s;
    }

    /**
     * Compares this vector to another vector first by size (2D vectors before 3D
     * vectors), then by the highest x-value (in descending order).
     * 
     * @param otherVector the other vector to compare with
     * @return a negative integer if this vector comes before the other vector, a
     *         positive integer if this vector comes after the other vector, or 0 if
     *         they are considered equal.
     */
    public int compareTo(MathVector otherVector)
    {
        if (this.is2DVector() && otherVector.is3DVector())
        {
            return -1; // 2D comes before 3D
        }
        else if (this.is3DVector() && otherVector.is2DVector())
        {
            return 1; // 3D comes after 2D
        }
        return Double.compare(otherVector.getXValue(), this.getXValue());
    }

}
