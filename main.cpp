#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

typedef struct Point {
    double x;
    double y;
} Point;

void calculateDistances(vector<vector<double >> &citiesCoordinatesArray,
                        vector<vector<double >> &distanceArray);

double countDistance(Point, Point);

void searchForPath(vector<vector<double> > &distanceArray);

double distance(vector<int> &permutation, vector<vector<double >> &distanceArray);

void reverseArray(int startIndex, int endIndex, vector<int> &array);


int main() {

    unsigned int numberOfCities;

    string input;

    cin >> numberOfCities;

    vector<vector<double> > citiesCoordinates(numberOfCities, vector<double>(3));

    for (int i = 0; i < numberOfCities; i++) {
        cin >> input;
        citiesCoordinates[i][0] = stod(input);
        cin >> input;
        citiesCoordinates[i][1] = stod(input);
        cin >> input;
        citiesCoordinates[i][2] = stod(input);
    }

    vector<vector<double> > distanceBetweenCities(numberOfCities, vector<double>(numberOfCities));

    calculateDistances(citiesCoordinates, distanceBetweenCities);

    searchForPath(distanceBetweenCities);

    return 0;
}


void calculateDistances(vector<vector<double >> &citiesCoordinatesArray,
                        vector<vector<double >> &distanceArray) {
    for (int columnIndex = 0; columnIndex < citiesCoordinatesArray.size(); columnIndex++) {
        for (int rowIndex = columnIndex + 1; rowIndex < citiesCoordinatesArray.size(); rowIndex++) {
            Point firstPoint, secondPoint;

            firstPoint.x = citiesCoordinatesArray[columnIndex][1];
            firstPoint.y = citiesCoordinatesArray[columnIndex][2];

            secondPoint.x = citiesCoordinatesArray[rowIndex][1];
            secondPoint.y = citiesCoordinatesArray[rowIndex][2];

            distanceArray[columnIndex][rowIndex] = countDistance(firstPoint, secondPoint);
            distanceArray[rowIndex][columnIndex] = distanceArray[columnIndex][rowIndex];
        }
        distanceArray[columnIndex][columnIndex] = 0;
    }
}

double countDistance(Point firstLocation, Point secondLocation) {
    return sqrt(pow((firstLocation.x - secondLocation.x), 2) + pow((firstLocation.y - secondLocation.y), 2));
}

double distance(vector<int> &permutation,
                vector<vector<double >> &distanceArray) {

    double distance = 0;

    for (int i = 1; i < permutation.size(); i++) {
        distance += distanceArray[permutation[i]][permutation[i - 1]];
    }
    distance += distanceArray[permutation[0]][permutation[permutation.size() - 1]];

    return distance;

}

int findGreedyNext(int previousCity,
                   bool *visitedList,
                   vector<vector<double >> &distanceArray) {

    int shortestNext = 0;
    double shortestPath = INFINITY;

    for (int checkedCity = 1; checkedCity < distanceArray.size(); checkedCity++) {
        if (visitedList[checkedCity] == 0) {
            if (distanceArray[previousCity][checkedCity] < shortestPath) {
                shortestPath = distanceArray[previousCity][checkedCity];
                shortestNext = checkedCity;
            }
        }
    }

    visitedList[shortestNext] = 1;

    return shortestNext;

}

vector<int> searchForInitialPath(vector<vector<double >> &distanceArray) {

    vector<int> newPermutation;

    int previousCity = 0;
    bool visitedList[distanceArray.size()];

    memset(visitedList, 0, sizeof(bool) * distanceArray.size());

    visitedList[previousCity] = 1;

    for (int i = 0; i < distanceArray.size(); i++) {
        newPermutation.push_back(findGreedyNext(previousCity, visitedList, distanceArray));
        previousCity = newPermutation[i];
    }
    newPermutation.push_back(0);

    return newPermutation;
}

double opt2(vector<int> &permutation,
            int firstEdgeFirstIndex,
            int secondEdgeFirstIndex,
            vector<vector<double> > &distanceArray) {
    double preOptDistance = distanceArray[permutation[firstEdgeFirstIndex]][permutation[firstEdgeFirstIndex + 1]] +
                            distanceArray[permutation[secondEdgeFirstIndex]][permutation[secondEdgeFirstIndex + 1]];

    double postOptDistance = distanceArray[permutation[firstEdgeFirstIndex]][permutation[secondEdgeFirstIndex]] +
                             distanceArray[permutation[firstEdgeFirstIndex + 1]][permutation[secondEdgeFirstIndex + 1]];

    if (postOptDistance < preOptDistance) {
        reverseArray(firstEdgeFirstIndex + 1, secondEdgeFirstIndex, permutation);
        return preOptDistance - postOptDistance;

    }
    return 0.0;
}


void reverseArray(int startIndex,
                  int endIndex,
                  vector<int> &array) {
    int tmp = 0;
    while (startIndex < endIndex) {
        tmp = array[startIndex];
        array[startIndex++] = array[endIndex];
        array[endIndex--] = tmp;

    }
}


double opt2Loop(int kValue,
                int kMinValue,
                int iAdder,
                vector<int> &optPermutation,
                double &optShortestDistance,
                vector<vector<double> > &distanceArray) {
    double difference = 0;
    double differenceSum = 0;
    if (iAdder == 0) {
        for (int k = kValue; k > kMinValue; k--) {
            for (int i = 0; i < optPermutation.size() - 3; i += k) {
                for (int j = i + 2; j < optPermutation.size() - 1; j++) {
                    difference = opt2(optPermutation, i, j, distanceArray);
                    differenceSum += difference;
                    optShortestDistance -= difference;
                }
            }
        }
    } else {

        for (int k = kValue; k > kMinValue; k--) {
            for (int i = 0; i < optPermutation.size() - 3; i += iAdder) {
                for (int j = i + 2; j < optPermutation.size() - 1; j++) {
                    difference = opt2(optPermutation, i, j, distanceArray);

                    differenceSum += difference;
                    optShortestDistance -= difference;
                }
            }
        }
    }

    return differenceSum;
}

void searchForPath(vector<vector<double> > &distanceArray) {

    vector<int> bestPermutation = searchForInitialPath(distanceArray);

    double shortestDistance = distance(bestPermutation, distanceArray);

    double optShortestDistance;
    vector<int> afterOptBest;
    optShortestDistance = shortestDistance;
    vector<unsigned int> neighbourhood;

    int CARENCY = 15;
    int ITERATIONS_MAX = 10000;
    int MIN_TO_LEAVE_TABU = 1;
    int numberOfIterations = 0;

    double initAverage = 0;

    uniform_int_distribution<int> uni(0, distanceArray.size());

    vector<int> tabuPermutation(bestPermutation);
    int tabuList[distanceArray.size()];

    int currentElement = 0;
    double averageDistance = 0;
    int iter = 0;
    while (numberOfIterations < ITERATIONS_MAX) {
        numberOfIterations++;
        iter++;
        while (tabuList[iter] != 0 && iter < distanceArray.size()) {
            iter++;
        }
        if (iter == distanceArray.size()) break;
        currentElement = iter;
        averageDistance = 0;
        for (unsigned int i = 0; i < distanceArray.size(); i++) {
            averageDistance += distanceArray[i][currentElement];
        }
        averageDistance /= distanceArray.size();

        initAverage = averageDistance;

        for (unsigned int i = 0; i < distanceArray.size(); i++) {
            if (distanceArray[currentElement][i] - averageDistance <
                averageDistance &&
                tabuList[i] < MIN_TO_LEAVE_TABU) {
                neighbourhood.push_back(i);
            }
        }

        for (auto firstNeighbour: neighbourhood) {
            for (auto secondNeighbour: neighbourhood) {
                if (firstNeighbour < secondNeighbour) {
                    opt2(tabuPermutation, firstNeighbour, secondNeighbour, distanceArray);
                }
            }
        }


        for (auto neighbour: neighbourhood) {
            tabuList[neighbour] += CARENCY * (distanceArray[neighbour][currentElement] / initAverage)
                                   - CARENCY * (numberOfIterations / ITERATIONS_MAX);
        }
        neighbourhood.clear();
    }

    vector<int> optPermutation(tabuPermutation);

    double ratio = distanceArray.size() / 14900.0;
    opt2Loop(int(80 * ratio), int(50 * ratio), int(25 * ratio), optPermutation,
             optShortestDistance, distanceArray);

    opt2Loop(int(50 * ratio), int(40 * ratio), int(15 * ratio), optPermutation, optShortestDistance,
             distanceArray);

    opt2Loop(int(50 * ratio), int(10 * ratio), 0, optPermutation, optShortestDistance,
             distanceArray);

    opt2Loop(int(50 * ratio), int(25 * ratio), 0, optPermutation, optShortestDistance,
             distanceArray);

    opt2Loop(10, 1, 1, optPermutation, optShortestDistance, distanceArray);

    optPermutation.insert(optPermutation.begin(), 0);
    optPermutation.pop_back();

    cout << distance(optPermutation, distanceArray) << endl;

    for (auto city: optPermutation) {
        cerr << city << " ";
    }

    return;
}
