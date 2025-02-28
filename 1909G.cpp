# include <bits/stdc++.h>
using namespace std;

/**
 * @brief Main function to solve the problem of counting valid (x, y, z) splits.
 *
 * Reads the input strings and computes the number of valid splits where
 * sourceString = x + y + z and targetString = x + (y repeated k times) + z.
 */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long sourceLength, targetLength;
    cin >> sourceLength >> targetLength;
    string sourceString, targetString;
    cin >> sourceString >> targetString;

    const long long kPrimeBase = 31;
    const long long kModValue[2] = {1000000007, 1000000009};
    vector<array<long long, 2>> power(targetLength + 1, {1, 1});

    for (long long i = 0; i < targetLength; i++) {
        for (long long j = 0; j < 2; j++) {
            power[i + 1][j] = (power[i][j] * kPrimeBase) % kModValue[j];
        }
    }

    auto computeHash = [&](const string& str) {
        long long length = str.size();
        vector<array<long long, 2>> hashValues(length + 1, {0, 0});
        for (long long i = 0; i < length; i++) {
            for (long long j = 0; j < 2; j++) {
                hashValues[i + 1][j] = (hashValues[i][j] + (str[i] - 'a' + 1) * power[i][j]) % kModValue[j];
            }
        }
        return hashValues;
    };

    auto sourceHash = computeHash(sourceString);
    auto targetHash = computeHash(targetString);

    auto calculateHash = [&](long long type, long long left, long long right) {
        array<long long, 2> result = {0, 0};
        vector<array<long long, 2>>* hashPtr = (type == 0) ? &sourceHash : &targetHash;
        for (long long j = 0; j < 2; j++) {
            result[j] = ((*hashPtr)[right + 1][j] - (*hashPtr)[left][j] + kModValue[j]) % kModValue[j];
            result[j] = (power[targetLength - left][j] * result[j]) % kModValue[j];
        }
        return result;
    };

    auto compareSubstrings = [&](long long sourceType, long long sourceLeft, long long sourceRight, long long targetType, long long targetLeft, long long targetRight) {
        assert((sourceRight - sourceLeft) == (targetRight - targetLeft));
        return calculateHash(sourceType, sourceLeft, sourceRight) == calculateHash(targetType, targetLeft, targetRight);
    };

    long long firstMismatch = sourceLength;
    for (long long i = 0; i < sourceLength; i++) {
        if (sourceString[i] != targetString[i]) {
            firstMismatch = i;
            break;
        }
    }

    long long lastMismatch = -1;
    for (long long i = sourceLength - 1; i >= 0; i--) {
        if (sourceString[i] != targetString[i + (targetLength - sourceLength)]) {
            lastMismatch = i;
            break;
        }
    }

    long long validSplitCount = 0;
    for (long long yLength = 1; yLength <= targetLength; yLength++) {
        if ((targetLength - sourceLength) % yLength != 0) continue;

        auto isValid = [&](long long xLength) {
            long long zLength = sourceLength - xLength - yLength;
            if (zLength < 0) return false;
            long long repeatedYLength = targetLength - xLength - zLength;
            if (!compareSubstrings(0, 0, xLength - 1, 1, 0, xLength - 1)) return false;
            if (!compareSubstrings(0, sourceLength - zLength, sourceLength - 1, 1, targetLength - zLength, targetLength - 1)) return false;
            if (!compareSubstrings(0, xLength, xLength + yLength - 1, 1, xLength, xLength + yLength - 1)) return false;
            if (!compareSubstrings(1, xLength, xLength + repeatedYLength - yLength - 1, 1, xLength + yLength, xLength + repeatedYLength - 1)) return false;
            return true;
        };

        if (isValid(lastMismatch + 1)) {
            validSplitCount += max(0LL, firstMismatch - lastMismatch - yLength);
        }
    }

    cout << validSplitCount << "\n";
    return 0;
}
