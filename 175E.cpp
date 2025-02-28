#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace std;

/**
 * @brief Computes the overlapping length between two intervals.
 *
 * @param start1  Start of the first interval.
 * @param end1    End of the first interval.
 * @param start2  Start of the second interval.
 * @param end2    End of the second interval.
 * @return double Overlapping length between the two intervals.
 */
double computeOverlap(double start1, double end1, double start2, double end2) {
    return max(min(end1, end2) - max(start1, start2), 0.0);
}

int main() {
    // Input variables
    int numFireTowers, numElectricTowers, numSlowTowers, fireTowerDamage, electricTowerDamage;
    double fireTowerRadius, electricTowerRadius, slowTowerRadius;

    // Read input values
    scanf("%d%d%d%lf%lf%lf%d%d", &numFireTowers, &numElectricTowers, &numSlowTowers, &fireTowerRadius, &electricTowerRadius, &slowTowerRadius, &fireTowerDamage, &electricTowerDamage);

    // Adjust radii based on the problem's constraints
    fireTowerRadius = sqrt(fireTowerRadius * fireTowerRadius - 1);
    electricTowerRadius = sqrt(electricTowerRadius * electricTowerRadius - 1);
    slowTowerRadius = sqrt(slowTowerRadius * slowTowerRadius - 1);

    // Compute total elements and mask for subset iteration
    int totalTowers = numFireTowers + numElectricTowers + numSlowTowers;
    int subsetMask = 1 << totalTowers;

    double maxDamage = 0.0;
    double overlapDifferences[120];

    // Iterate over all subsets
    for (int subset = 0; subset < subsetMask; subset++) {
        // Ensure the subset contains exactly `numSlowTowers` elements
        if (__builtin_popcount(subset) != numSlowTowers) continue;

        double totalWeight = 0.0;
        int count = 0;

        // Iterate over remaining elements
        for (int remaining = subsetMask - subset - 1; remaining; remaining -= remaining & -remaining) {
            int towerPosition1 = __builtin_ctz(remaining) >> 1; // Extract position
            double fireTowerWeight = fireTowerDamage * fireTowerRadius * 2;
            double electricTowerWeight = electricTowerDamage * electricTowerRadius * 2;

            // Iterate over towers in the subset
            for (int k = subset; k; k -= k & -k) {
                int towerPosition2 = __builtin_ctz(k) >> 1;
                fireTowerWeight += fireTowerDamage * computeOverlap(towerPosition1 - slowTowerRadius, towerPosition1 + slowTowerRadius, towerPosition2 - fireTowerRadius, towerPosition2 + fireTowerRadius);
                electricTowerWeight += electricTowerDamage * computeOverlap(towerPosition1 - slowTowerRadius, towerPosition1 + slowTowerRadius, towerPosition2 - electricTowerRadius, towerPosition2 + electricTowerRadius);
            }

            totalWeight += fireTowerWeight;
            overlapDifferences[count++] = fireTowerWeight - electricTowerWeight;
        }

        // Sort overlap differences and subtract the smallest `numElectricTowers` elements
        sort(overlapDifferences, overlapDifferences + count);
        for (int i = 0; i < numElectricTowers; i++) {
            totalWeight -= overlapDifferences[i];
        }

        // Update the maximum damage found
        maxDamage = max(maxDamage, totalWeight);
    }

    // Output the final result with precision
    printf("%.10f\n", maxDamage);
    return 0;
}
