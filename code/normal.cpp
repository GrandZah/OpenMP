#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <cstdio>

using namespace std;

double halton(int index, int base) {
    double result = 0;
    double f = 1.0 / base;
    int i = index;

    while (i > 0) {
        result = result + f * (i % base);
        i = i / base;
        f = f / base;
    }

    return result;
}

pair<double, int> calculateWithOpenMP(double diag, int totalPoints, int numThreads) {
    int pointsInsideCircle = 0;
    int size_chunk = max(1, static_cast<int>(pow(2, static_cast<int>(ceil(1.6 * log10(totalPoints)))) / 8 * max(1, numThreads)));

    #pragma omp parallel
    {
        int localCount = 0;

        #pragma omp for schedule(dynamic, size_chunk)
        for (int i = 1; i <= totalPoints; i++) {
            double x = halton(i, 2) * diag;
            double y = halton(i, 3) * diag;
            double z = halton(i, 5) * diag;

            if (abs(x) + abs(y) + abs(z) <= diag) {
                localCount++;
            }
        }

        #pragma omp atomic
        pointsInsideCircle += localCount;
    }

    return make_pair((static_cast<double>(pointsInsideCircle) / totalPoints) *
        (diag * diag * diag) * 8
        , omp_get_num_threads());
}


double calculateSequentially(double diag, int totalPoints) {
    int pointsInsideCircle = 0;

    for (int i = 0; i < totalPoints; i++) {
        double x = halton(i, 2) * diag;
        double y = halton(i, 3) * diag;
        double z = halton(i, 5) * diag;

        if (abs(x) + abs(y) + abs(z) <= diag) {
            pointsInsideCircle++;
        }
    }

    return (static_cast<double>(pointsInsideCircle) / totalPoints) * 
        (diag * diag * diag) * 8;
}


int main(int argc, char* argv[]) {
    setlocale(LC_ALL, "rus");
    setlocale(LC_ALL, "C");

    if (argc < 4) {
        cerr << "There are not enough arguments. Usage: " << argv[0] << " <number of streams> <input filename> <output filename>" << endl;
        return 1;
    }

    int numThreads = atoi(argv[1]);
    string inputFileName = argv[2];
    string outputFileName = argv[3];

    int totalPoints;
    double p1[3], p2[3], p3[3];
    FILE* file;


    file = fopen(inputFileName.c_str(), "r");

    if (file == NULL) {
        printf("Ошибка открытия файла");
        return 1;
    }

    fscanf(file, "%d\n", &totalPoints);
    fscanf(file, "(%lf %lf %lf)\n", &p1[0], &p1[1], &p1[2]);
    fscanf(file, "(%lf %lf %lf)\n", &p2[0], &p2[1], &p2[2]);
    fscanf(file, "(%lf %lf %lf)\n", &p3[0], &p3[1], &p3[2]);

    fclose(file);

    double r12 = sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
    double r13 = sqrt(pow(p1[0] - p3[0], 2) + pow(p1[1] - p3[1], 2) + pow(p1[2] - p3[2], 2));
    double r23 = sqrt(pow(p2[0] - p3[0], 2) + pow(p2[1] - p3[1], 2) + pow(p2[2] - p3[2], 2));

    double diag = 0;
    if (r12 == r13) {
        diag = r23 / 2;
    }
    else if (r12 == r23) {
        diag = r13 / 2;
    }
    else {
        diag = r12 / 2;
    }
    double realVolumeOct = pow(diag, 3) / 2 / 3 * 8;

    double start_time, end_time;
    double volumeOct = 0;

    if (numThreads == -1) {
        start_time = omp_get_wtime();
        volumeOct = calculateSequentially(diag, totalPoints);
        end_time = omp_get_wtime();
        printf("Time (%i thread(s)): %g ms\n", 0, (end_time - start_time) * 1000);

    }
    else if (numThreads >= 0) {
        if (numThreads == 0) {
        }
        else {
            omp_set_num_threads(numThreads);
        }
        int numThreadsToPrint = 0;
        start_time = omp_get_wtime();
        pair<double, int> result = calculateWithOpenMP(diag, totalPoints, numThreads);
        volumeOct = result.first;
        numThreadsToPrint = result.second;
        end_time = omp_get_wtime();

        printf("Time (%i thread(s)): %g ms\n", numThreadsToPrint, (end_time - start_time) * 1000);

    }
    else {
        cerr << "Unsupported number of threads" << endl;
        return 1;
    }


    FILE* file2 = fopen(outputFileName.c_str(), "w");


    try {
        fprintf(file2, "%g %g\n", realVolumeOct, volumeOct);
    }
    catch (const exception& e) {
        cerr << "Error writing the file: " << e.what() << endl;
        return 1;
    }
    fclose(file2);


    return 0;

}