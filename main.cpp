#include <iostream>
#include <vector>
#include <fstream>
#include <ctgmath>
#include <algorithm>
#include <omp.h>
#include <chrono>

using ns = std::chrono::nanoseconds;
using get_time = std::chrono::steady_clock ;

//someglobal constants
const int numThreads {4};
const int readLimiter {10};
const int edgeLength {500};

//some mutable global variables, such as results array
float results[edgeLength][edgeLength][edgeLength];
int per;

// structs for our data
struct Atom{
    float xpos;
    float ypos;
    float zpos;
    float energy;
    Atom(float xpos, float ypos, float zpos, float energy) : xpos{xpos}, ypos{ypos}, zpos{zpos}, energy{energy}{}
    Atom() : xpos{0.0f}, ypos{0.0f}, zpos{0.0f}, energy{0.0f} {}

};

//data vecotr for reading input
Atom atoms [readLimiter];



//global functions
void readAtomsFromFile(){
    //reading atoms_vector from atoms_vector set file
    std::ifstream myfile("atoms.txt", std::ios_base::in);

    float a {0};
    float x_in {0};
    float y_in {0};
    float z_in {0};
    float e_in {0};
    int counter{0};

    myfile >> a; // take the 1000 at the start out
    int pos {0};
    while (myfile >> a) {
        if (counter == 0) { x_in = a; }
        if (counter == 1) { y_in = a; }
        if (counter == 2) { z_in = a; }
        if (counter == 3) {
            e_in = a;
            counter = 0;
            atoms[pos] = Atom { x_in, y_in, z_in, e_in };
            //std::cout << x_in <<" " << y_in << " " << z_in << " "<< e_in << " " << pos << std::endl;
            //std::cout << atoms[pos].xpos <<" " << atoms[pos].ypos << " " << atoms[pos].zpos << " "<< atoms[pos].energy << std::endl;
            pos++;
            continue;
        }

        counter++;
        //if (counter%3==0) std::cout << pos << std::endl;
        //if((++limiter/3) > readLimiter) break;
        if(pos > readLimiter) break;
        if (myfile.eof()) break;
    }
}

float euclidianDistance(int x, int y, int z, Atom atom){

    float a = ((float)x) - atom.xpos;
    float b = ((float)y) - atom.ypos;
    float c= ((float)z) - atom.zpos;

    return std::sqrt( (a*a) + (b*b) + (c*c) );

}

float calculate(int x, int y, int z, Atom* localAtoms){
        float result {0};
        //#pragma omp private (result)

        float euclidianDistanceMin{0};
        float euclidianDistanceCurrent{0};

        for (int i {0}; i < 10; i++){
            euclidianDistanceMin = euclidianDistance(x, y, z, localAtoms[i]);
            for (int j = i+1; j < readLimiter; j++){
                euclidianDistanceCurrent = euclidianDistance(x, y, z, localAtoms[j]);
                if (euclidianDistanceCurrent < euclidianDistanceMin){
                    std::swap(localAtoms[i],  localAtoms[j]);
                    euclidianDistanceMin = euclidianDistanceCurrent;
                }
            }
            result += ((localAtoms[i].energy)/euclidianDistanceMin);
        }

        //std::cout << euclidianDistance(x, y, z, localAtoms[0]) << " " << euclidianDistance(x, y, z, localAtoms[1]) << " " << euclidianDistance(x, y, z, localAtoms[2]) << std::endl;
/*
        std::sort(localAtoms, (localAtoms+readLimiter), [x,y,z](Atom a, Atom b)->bool{ return euclidianDistance(x, y, z, a) < euclidianDistance(x, y, z, b);});

        for (int i {0}; i < 10; i++){ //must be smalle rthan readLimiter, otherwise we get undefined behaviour because we read beyond the array
            result += ((localAtoms[i].energy)/(euclidianDistance(x, y, z, localAtoms[i])));
            //std::cout << localAtoms[i].energy << std::endl;
        }
*/
        //std::cout << result << std::endl;
        return result;
}

void calculateLoops(){

    //for all points in the grid do the following
    #pragma omp parallel shared(results, atoms) num_threads(numThreads)
    {
        //std::cout << omp_get_num_threads() << std::endl;
        Atom localAtoms[readLimiter];
        std::copy(atoms, atoms+(sizeof(atoms)/sizeof(*atoms)), localAtoms);

        #pragma omp for schedule(dynamic) ordered
                for (int x = 0; x < edgeLength; x++) {
                    for (int y = 0; y < edgeLength; y++) {
                        for (int z = 0; z < edgeLength; z++) {
                            results[x][y][z] = calculate(x,y,z, localAtoms);
                        }
                    }
                    if (x!=0&&!(x%5))
                        std::cout << ++per << '%'<< std::endl;
                }
    }
}

void logResults(){
    std::ofstream out_file;
    out_file.open ("gridpoints.txt");
    out_file << results [0][0][0] << std::endl;
    out_file << results [23][456][289] << std::endl;
    out_file << results [368][34][189] << std::endl;
    out_file << results [451][216][45] << std::endl;
    out_file << results [499][499][499] << std::endl;
    out_file.close();


}

void logExecutionTime(long execTime){
    std::ofstream out_file;
    out_file.open ("duration.txt");
    out_file << "execution duration: " << ((double)execTime)*(1,66667e-11) << " minutes" << std::endl;
    out_file.close();
}

void logError(std::exception &e){
    std::ofstream out_file;
    out_file.open ("log.txt");
    out_file << "execution halted due to error: " << e.what() << std::endl;
    out_file.close();
}

// do what what a c++ program must do: execute the main function
int main() {
    try{
        // read input
        readAtomsFromFile();

        //(parallel) calculation of all values with time measuremennt of execution time
        auto start = get_time::now();
        calculateLoops();
        auto end = get_time::now();

        //write results to file
        logResults();
        logExecutionTime(std::chrono::duration_cast<ns>((end - start)).count());

    } catch (std::exception &e){
        //if an exception is thrown, log it to a file. program is supposed to run as a daemon due to long execution time so no one will look at cmd line.
        logError(e);
    }


    return 0;
}


/*
    //direct sequential calculation for testing purposes
    std::cout << calculate(0,0,0, atoms) << std::endl;;
    std::cout << calculate(23,456,289, atoms) << std::endl;;
    std::cout << calculate(368,34,189, atoms) << std::endl;;
    std::cout << calculate(499,499,499, atoms)<< std::endl;;
*/

/*
    std::cout << results [0][0][0] << std::endl;
    std::cout << results [23][456][289] << std::endl;
    std::cout << results [368][34][189] << std::endl;
    std::cout << results [451][216][45] << std::endl;
    std::cout << results [499][499][499] << std::endl;
*/