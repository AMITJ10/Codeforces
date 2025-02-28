
#include <bits/stdc++.h>
using namespace std;

// Constants
const int maxPoints = 1e5 + 10;
const double epsilon = 1e-8;

// Global variables
int numPoints, numSegments, functionCount, eventCount;
double segmentLengths[maxPoints], totalPerimeter, aCoeff, bCoeff, cCoeff, dCoeff;

/**
 * @brief Calculates the square of a number.
 * @param value The number to be squared.
 * @return The square of the input value.
 */
inline double square(double value) { return value * value; }

/**
 * @brief Structure to represent a 2D point.
 */
struct Point {
    double x, y;

    /**
     * @brief Constructor for the Point structure.
     * @param xCoord The x-coordinate of the point. Defaults to 0.
     * @param yCoord The y-coordinate of the point. Defaults to 0.
     */
    Point(double xCoord = 0, double yCoord = 0) : x(xCoord), y(yCoord) {}

    /**
     * @brief Overloads the addition operator for Point objects.
     * @param other The other Point object to add.
     * @return A new Point object representing the sum of the two points.
     */
    Point operator+(const Point& other) const { return Point(x + other.x, y + other.y); }

    /**
     * @brief Overloads the subtraction operator for Point objects.
     * @param other The other Point object to subtract.
     * @return A new Point object representing the difference of the two points.
     */
    Point operator-(const Point& other) const { return Point(x - other.x, y - other.y); }

    /**
     * @brief Overloads the division operator for Point objects with a scalar.
     * @param scalar The scalar value to divide by.
     * @return A new Point object representing the division of the point by the scalar.
     */
    Point operator/(double scalar) const { return Point(x / scalar, y / scalar); }

    /**
     * @brief Overloads the multiplication operator for Point objects with a scalar.
     * @param scalar The scalar value to multiply by.
     * @return A new Point object representing the multiplication of the point by the scalar.
     */
    Point operator*(double scalar) const { return Point(x * scalar, y * scalar); }

    /**
     * @brief Calculates the squared magnitude (length) of the Point vector.
     * @return The squared magnitude of the Point.
     */
    double magnitudeSquared() const { return square(x) + square(y); }
} polygonPoints[maxPoints], directionVectors[maxPoints], startPoint, endPoint;

/**
 * @brief Calculates the Euclidean distance between two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @return The Euclidean distance between the two points.
 */
inline double distanceBetween(const Point& p1, const Point& p2) {
    return sqrt((p1 - p2).magnitudeSquared());
}

/**
 * @brief Structure to represent a quadratic function.
 *        Used in the sweep line algorithm to represent distance functions.
 */
struct QuadraticFunction {
    double startParam, endParam, a, b, c;
} quadraticFunctions[maxPoints << 4];

/**
 * @brief Structure to represent an event in the sweep line algorithm.
 *        Events are sorted by position and then by value to ensure stability.
 */
struct Event {
    double position;
    int value;

    /**
     * @brief Overloads the less than operator for Event objects.
     *        Used for sorting events in the sweep line algorithm.
     * @param other The other Event object to compare with.
     * @return True if this event should come before the other event in sorted order, false otherwise.
     */
    bool operator<(const Event& other) const {
        if (fabs(position - other.position) > epsilon) return position < other.position;
        return value < other.value;  // Ensures stability in sorting
    }
} events[maxPoints << 4];

/**
 * @brief Checks if a given limit is achievable using a sweep line algorithm.
 *        This function determines if it's possible to divide the polygon perimeter into segments
 *        such that the squared distance between corresponding points on different segments is no more than 'limit'.
 * @param limit The squared distance limit to check for achievability.
 * @return True if the limit is achievable, meaning a valid segmentation exists; false otherwise.
 */
bool isLimitAchievable(double limit) {
    eventCount = 0;
    double start, end, discriminant, root1, root2;
    for (int i = 1; i <= functionCount; ++i) {
        aCoeff = quadraticFunctions[i].a;
        bCoeff = quadraticFunctions[i].b;
        cCoeff = quadraticFunctions[i].c - limit;
        start = quadraticFunctions[i].startParam;
        end = quadraticFunctions[i].endParam;

        if (aCoeff < -epsilon) return false;  // Safer than assert

        if (fabs(aCoeff) < epsilon) {
            if (fabs(bCoeff) < epsilon) {
                if (cCoeff < epsilon) {
                    events[++eventCount] = Event{start, 1};
                    events[++eventCount] = Event{end, -1};
                }
            } else {
                double z = -cCoeff / bCoeff;
                if (bCoeff > 0 && z > 0) {
                    events[++eventCount] = Event{start, 1};
                    events[++eventCount] = Event{min(std::nextafter(start + z, end), end), -1};
                } else if (bCoeff < 0 && start + z < end) {
                    events[++eventCount] = Event{max(std::nextafter(start + z, 0.0), 0.0), 1};
                    events[++eventCount] = Event{end, -1};
                }
            }
        } else {
            discriminant = bCoeff * bCoeff - 4 * aCoeff * cCoeff;
            if (discriminant < 0) continue;
            discriminant = sqrt(discriminant);
            root1 = (-bCoeff - discriminant) / (2 * aCoeff);
            root2 = (-bCoeff + discriminant) / (2 * aCoeff);
            if (start + root1 > end || root2 < 0) continue;
            events[++eventCount] = Event{std::nextafter(start + max(0.0, root1), end), 1};
            events[++eventCount] = Event{std::nextafter(min(end, start + root2), start), -1};
        }
    }
    sort(events + 1, events + eventCount + 1);
    for (int i = 1, activeCount = 0; i <= eventCount; ++i) {
        activeCount += events[i].value;
        if (activeCount == numSegments) return true;
    }
    return false;
}

/**
 * @brief Main function of the program.
 *        Reads input, sets up the problem, performs binary search to find the minimum limit, and outputs the result.
 * @return 0 if the program executes successfully.
 */
int main() {
    int leftIndex = 1, rightIndex = 1, segmentId = 1;
    Point directionVector, deltaVector, startPoint, endPoint;
    double low = 0, high, mid, segmentLength, currentLength = 0, leftRemaining, rightRemaining, averageLength;

    // Input reading
    scanf("%d%d", &numPoints, &numSegments);
    for (int i = 1; i <= numPoints; ++i) {
        scanf("%lf%lf", &polygonPoints[i].x, &polygonPoints[i].y);
    }
    polygonPoints[numPoints + 1] = polygonPoints[1];

    // Calculate segment lengths and total perimeter
    for (int i = 1; i <= numPoints; ++i) {
        segmentLengths[i] = distanceBetween(polygonPoints[i], polygonPoints[i + 1]);
        directionVectors[i] = (polygonPoints[i + 1] - polygonPoints[i]) / segmentLengths[i];
        totalPerimeter += segmentLengths[i];
    }

    // Initialize variables for the sweep
    averageLength = rightRemaining = totalPerimeter / numSegments;
    while (segmentLengths[rightIndex] < rightRemaining + epsilon) {
        rightRemaining -= segmentLengths[rightIndex++];
    }
    startPoint = polygonPoints[1];
    endPoint = polygonPoints[rightIndex] + directionVectors[rightIndex] * rightRemaining;
    leftRemaining = segmentLengths[1];
    rightRemaining = segmentLengths[rightIndex] - rightRemaining;

    // Generate quadratic functions for each segment
    while (currentLength + (1e-5) < totalPerimeter) {
        segmentLength = min(segmentId * averageLength - currentLength, min(leftRemaining, rightRemaining));
        deltaVector = endPoint - startPoint;
        directionVector = directionVectors[rightIndex] - directionVectors[leftIndex];
        quadraticFunctions[++functionCount] = QuadraticFunction{
            currentLength - (segmentId - 1) * averageLength,
            currentLength - (segmentId - 1) * averageLength + segmentLength,
            directionVector.magnitudeSquared(),
            2 * (directionVector.x * deltaVector.x + directionVector.y * deltaVector.y),
            deltaVector.magnitudeSquared()
        };
        currentLength += segmentLength;
        if (currentLength + (1e-5) > segmentId * averageLength) segmentId++;
        startPoint = startPoint + (directionVectors[leftIndex] * segmentLength);
        endPoint = endPoint + (directionVectors[rightIndex] * segmentLength);
        if (segmentLength + epsilon > leftRemaining) {
            leftIndex = leftIndex % numPoints + 1;
            leftRemaining = segmentLengths[leftIndex];
        } else {
            leftRemaining -= segmentLength;
        }
        if (segmentLength + epsilon > rightRemaining) {
            rightIndex = rightIndex % numPoints + 1;
            rightRemaining = segmentLengths[rightIndex];
        } else {
            rightRemaining -= segmentLength;
        }
    }

    // Binary search to find the minimum achievable limit
    for (high = averageLength; high - low > epsilon;) {
        mid = (low + high) / 2.0;
        isLimitAchievable(mid * mid) ? high = mid : low = mid;
    }

    // Output the result
    printf("%.8lf", low);
    return 0;
}
