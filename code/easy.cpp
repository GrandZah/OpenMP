#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <fstream>

using namespace std;

double calculateWithOpenMP(double radius, int totalPoints, int numThreads) {
    int pointsInsideCircle = 0;
    int size_chunk = max(1, int(pow(2, (int)(ceil(1.6 * log10(totalPoints)))) / 8 * max(1, numThreads)));

    #pragma omp parallel
    {
        int localCount = 0;
        int thread_id = omp_get_thread_num();
        srand(static_cast<unsigned int>(time(nullptr)) + thread_id);
        numThreads = omp_get_num_threads();

        #pragma omp for schedule(dynamic, size_chunk)
        for (int i = 0; i < totalPoints; i++) {
            double x = (double)rand() / RAND_MAX * radius;
            double y = (double)rand() / RAND_MAX * radius;

            if ((x * x + y * y) <= (radius * radius)) {
                localCount++;
            }
        }

        #pragma omp atomic
        pointsInsideCircle += localCount;
    }

    return (static_cast<double>(pointsInsideCircle) / totalPoints) * (4 * radius * radius), numThreads;
}


double calculateSequentially(double radius, int totalPoints) {
    int pointsInsideCircle = 0;

    for (int i = 0; i < totalPoints; i++) {
        double x = (double)rand() / RAND_MAX * radius;
        double y = (double)rand() / RAND_MAX * radius;

        if ((x * x + y * y) <= (radius * radius)) {
            pointsInsideCircle++;
        }
    }

    return (static_cast<double>(pointsInsideCircle) / totalPoints) * (4 * radius * radius);
}


int main(int argc, char* argv[]) {
    setlocale(LC_ALL, "rus");

    if (argc < 4) {
        cerr << "There are not enough arguments. Usage: " << argv[0] << " <number of streams> <input filename> <output filename>" << endl;
        return 1;
    }

    int numThreads = atoi(argv[1]);
    string inputFileName = argv[2];
    string outputFileName = argv[3];

    double radius;
    int totalPoints;

    ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        cerr << "The file could not be opened " << inputFileName << endl;
        return 1;
    }

    try {
        if (inputFile >> totalPoints >> radius) {
            if (radius > 0 && totalPoints > 0) {
            }
            else {
                std::cerr << "Numbers are not positive" << std::endl;
                return 1;
            }
        }
        else {
            std::cerr << "Error reading numbers from file" << std::endl;
            return 1;
        }
    }
    catch (const exception& e) {
        cerr << "Error reading the file: " << e.what() << endl;
        return 1;
    }
    inputFile.close();



    double start_time, end_time;
    double circleArea = 0;

    if (numThreads == -1) {
        start_time = omp_get_wtime();
        circleArea = calculateSequentially(radius, totalPoints);
        end_time = omp_get_wtime();
        printf("Time (%i thread(s)): %g ms\n", 0, (end_time - start_time) * 1000);

    }
    else if (numThreads >= 0) {
        if (numThreads == 0) {
        }
        else {
            omp_set_num_threads(numThreads);
        }
        int numThreadsToPrint;
        start_time = omp_get_wtime();
        circleArea, numThreadsToPrint = calculateWithOpenMP(radius, totalPoints, numThreads);
        end_time = omp_get_wtime();

        printf("Time (%i thread(s)): %g ms\n", numThreadsToPrint, (end_time - start_time) * 1000);

    }
    else {
        cerr << "Unsupported number of threads" << endl;
        return 1;
    }


    FILE* file = fopen(outputFileName.c_str(), "w");

    if (file == NULL) {
        cerr << "File opening error!" << endl;
        return 1;
    }

    try {
        fprintf(file, "%g %g\n", 3.14159265358979323846 * radius * radius, circleArea);
    }
    catch (const exception& e) {
        cerr << "Error writing the file: " << e.what() << endl;
        return 1;
    }
    fclose(file);


    return 0;
}